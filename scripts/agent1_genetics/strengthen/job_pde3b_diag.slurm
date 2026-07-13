#!/bin/bash
#SBATCH -J pde3b_diag
#SBATCH -o /path/to/cad-genetics/logs/claude_jobs/pde3b_diag_%j.out
#SBATCH -e /path/to/cad-genetics/logs/claude_jobs/pde3b_diag_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A YOUR_ALLOCATION
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
cd /path/to/cad-genetics/results/agent1_genetics/strengthen
Rscript - <<'RS'
suppressMessages({library(data.table); library(survival)})
OUT <- "."; PH <- "/path/to/cad-genetics/results/agent1_genetics/phenotypes"
invn <- function(x) qnorm((frank(x,ties.method="average")-0.5)/sum(!is.na(x)))
cov <- fread(file.path(PH,"covariates.txt")); qt <- fread(file.path(PH,"pheno_quantitative.txt"))
cmap <- fread("pde3b_conseq.map", header=FALSE, col.names=c("vid","cls"))
setkey(cmap, vid)
foldkeep <- function(raw){                        # fold to minor-allele dosage, keep MAF<0.01
  vc <- names(raw)[7:ncol(raw)]
  af <- sapply(vc, function(v) mean(raw[[v]], na.rm=TRUE)/2)
  M <- as.data.table(lapply(vc, function(v){x<-raw[[v]];a<-af[v]; if(is.na(a)) rep(0,length(x)) else if(a>0.5) 2-x else x}))
  setnames(M, vc); maf <- pmin(af,1-af); maf[is.na(maf)] <- 1
  list(M=M, vc=vc, keep=vc[maf<0.01], maf=maf)
}
vid_of <- function(col) sub("_[^_]*$","",col)     # raw col -> variant id

## ---------- severity classes + P_trend (uses ALL variants) ----------
rawA <- fread("pde3b_all.raw"); fk <- foldkeep(rawA); M <- fk$M
cls  <- cmap[.(vid_of(fk$keep)), cls]              # class per kept rare variant
carr_cls <- function(c){ v <- fk$keep[which(cls==c)]; if(!length(v)) return(rep(0L,nrow(M)));
                         as.integer(rowSums(M[,..v],na.rm=TRUE)>=1) }
plofc <- carr_cls("pLoF"); misc <- carr_cls("missense"); sync <- carr_cls("synonymous")
sev <- fifelse(plofc==1,2L, fifelse(misc==1,1L, 0L))     # 0 none, 1 missense, 2 pLoF (pLoF precedence)
dd <- data.table(IID=as.integer(rawA$IID), sev=sev, syn=sync, plof=plofc, mis=misc)
dd <- merge(dd, cov, by="IID"); dd <- merge(dd, qt[,.(IID,TRIGLY)], by="IID")
dd <- dd[!is.na(TRIGLY)]; dd[, y := invn(TRIGLY)]
ctrl <- paste(c("age","sex",paste0("PC",1:10)), collapse=" + ")
ft <- lm(as.formula(paste("y ~ sev +", ctrl)), data=dd[!(sev==0 & syn==1)])   # trend over none/mis/pLoF
s  <- summary(ft)$coefficients["sev",]
synf <- summary(lm(as.formula(paste("y ~ syn +",ctrl)), data=dd[sev==0]))$coefficients["syn",]
fwrite(data.table(term=c("severity_trend","synonymous_control"),
  beta=c(s[1],synf[1]), se=c(s[2],synf[2]), P=c(s[4],synf[4]),
  n_pLoF=sum(dd$plof), n_missense=sum(dd$mis), n_syn=sum(dd$syn)), "pde3b_ptrend.csv")
cat(sprintf("P_trend=%.2e  syn-control beta=%.3f P=%.2g\n", s[4], synf[1], synf[4]))

## ---------- leave-one-variant-out + per-variant (pLoF mask) ----------
rawP <- fread("pde3b_plof.raw"); fp <- foldkeep(rawP); Mp <- fp$M; keepP <- fp$keep
dp <- data.table(IID=as.integer(rawP$IID))
dp <- merge(dp, cov, by="IID"); dp <- merge(dp, qt[,.(IID,TRIGLY)], by="IID")
ord <- match(dp$IID, as.integer(rawP$IID))
Mp <- Mp[ord]; dp <- dp[!is.na(TRIGLY)]; keeprows <- !is.na(dp$TRIGLY)
Mp <- Mp[which(keeprows)]; dp[, y := invn(TRIGLY)]
full_carr <- as.integer(rowSums(Mp[,..keepP],na.rm=TRUE)>=1); dp[, carr := full_carr]
bfull <- summary(lm(as.formula(paste("y ~ carr +",ctrl)), data=dp))$coefficients["carr",]
lovo <- rbindlist(lapply(keepP, function(v){
  rest <- setdiff(keepP, v)
  cc <- as.integer(rowSums(Mp[,..rest],na.rm=TRUE)>=1)
  f <- summary(lm(as.formula(paste("yy ~ cc +",ctrl)), data=cbind(dp[,.(age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)], yy=dp$y, cc=cc)))$coefficients["cc",]
  sv <- as.integer(Mp[[v]]>=1)
  fs <- tryCatch(summary(lm(as.formula(paste("yy ~ vv +",ctrl)), data=cbind(dp[,.(age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)], yy=dp$y, vv=sv)))$coefficients["vv",], error=function(e) c(NA,NA,NA,NA))
  data.table(variant=vid_of(v), n_carriers=sum(sv,na.rm=TRUE),
    drop_beta=f[1], drop_P=f[4], single_beta=fs[1], single_P=fs[4])}))
setorder(lovo,-n_carriers)
fwrite(data.table(term="full_pLoF_burden", beta=bfull[1], se=bfull[2], P=bfull[4],
                  n_variants=length(keepP), n_carriers=sum(full_carr)), "pde3b_burden_full.csv")
fwrite(lovo, "pde3b_lovo.csv")
cat(sprintf("full pLoF burden beta=%.3f P=%.2e ; LOVO beta range %.3f..%.3f\n",
    bfull[1], bfull[4], min(lovo$drop_beta), max(lovo$drop_beta)))

## ---------- survival diagnostics (KM / Schoenfeld / power) ----------
carr <- data.table(IID=as.integer(rawP$IID), carrier=as.integer(rowSums(fp$M[,..keepP],na.rm=TRUE)>=1))
D <- "/path/to/ukb_phenotypes"
rd <- function(fid,nm){d<-fread(sprintf("%s/fid%s.csv",D,fid)); x<-data.table(IID=as.integer(d[[1]]), v=as.Date(as.character(d[[2]]))); setnames(x,"v",nm); x}
base<-rd("53","base"); cad<-rd("131306","cad"); dth<-rd("40000","death")
m <- Reduce(function(a,b) merge(a,b,by="IID",all.x=TRUE), list(cov, carr, base, cad, dth))
admin <- max(m$cad, na.rm=TRUE); m <- m[!is.na(base)]; m <- m[is.na(cad) | cad > base]
m[, event := as.integer(!is.na(cad) & cad > base)]
m[, cens := pmin(fifelse(is.na(death), admin, death), admin, na.rm=TRUE)]
m[, endd := as.Date(fifelse(event==1, cad, cens))]; m[, time := as.numeric(endd-base)/365.25]
m <- m[time>0 & !is.na(carrier) & !is.na(age) & !is.na(sex)]
fit <- coxph(Surv(time,event) ~ carrier + age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=m)
zph <- cox.zph(fit)
fwrite(data.table(term=rownames(zph$table), chisq=zph$table[,"chisq"], df=zph$table[,"df"], P=zph$table[,"p"]), "pde3b_zph.csv")
sf <- survfit(Surv(time,event) ~ carrier, data=m)
grid <- seq(0, 14, by=0.25)
km <- rbindlist(lapply(c(0,1), function(g){
  ss <- summary(sf[paste0("carrier=",g)], times=grid, extend=TRUE)
  data.table(carrier=g, time=ss$time, cuminc=1-ss$surv, n_risk=ss$n.risk)}))
fwrite(km, "pde3b_km.csv")
p <- mean(m$carrier); d <- sum(m$event); se <- 1/sqrt(d*p*(1-p)); za <- qnorm(0.975)
hr <- seq(0.5,1.6,by=0.02); pow <- pnorm(abs(log(hr))/se - za)
det80 <- exp(-(za+qnorm(0.8))*se)
fwrite(data.table(HR=hr, power=pow), "pde3b_power.csv")
fwrite(data.table(events=d, carriers=sum(m$carrier), carrier_events=sum(m$event[m$carrier==1]),
  prop_carrier=p, detectable_HR_80pct=det80, zph_carrier_P=zph$table["carrier","p"]), "pde3b_surv_summary.csv")
cat(sprintf("events=%d carriers=%d detectable HR@80%%=%.3f (protective side %.3f) ; Schoenfeld carrier P=%.3f\n",
    d, sum(m$carrier), det80, 1/det80, zph$table["carrier","p"]))
cat("PDE3B DIAG DONE\n")
RS
