suppressMessages({library(data.table); library(survival)})
raw <- fread("pde3b_plof.raw")
vc <- names(raw)[7:ncol(raw)]
af <- sapply(vc, function(v) mean(raw[[v]], na.rm=TRUE)/2)
M <- as.data.table(lapply(vc, function(v){x<-raw[[v]];a<-af[v];if(is.na(a))rep(0,length(x)) else if(a>0.5) 2-x else x})); setnames(M,vc)
maf <- pmin(af,1-af); maf[is.na(maf)]<-1; keep <- vc[maf<0.01]
raw[, ncopy := rowSums(M[,..keep], na.rm=TRUE)]
carr <- data.table(IID=as.integer(raw$IID), carrier=as.integer(raw$ncopy>=1))
cat(sprintf("pLoF variants kept: %d ; carriers: %d\n", length(keep), sum(carr$carrier)))
D <- "/path/to/ukb_phenotypes"
rd <- function(fid,nm){d<-fread(sprintf("%s/fid%s.csv",D,fid)); data.table(IID=as.integer(d[[1]]), x=as.Date(as.character(d[[2]])))[,setNames(.SD,c("IID",nm))]}
base<-rd("53","base"); cad<-rd("131306","cad"); dth<-rd("40000","death")
cov <- fread("/path/to/cad-genetics/results/agent1_genetics/phenotypes/covariates.txt")
m <- Reduce(function(a,b) merge(a,b,by="IID",all.x=TRUE), list(cov, carr, base, cad, dth))
admin <- max(m$cad, na.rm=TRUE)
m <- m[!is.na(base)]
prev <- sum(!is.na(m$cad) & m$cad <= m$base)
m <- m[is.na(cad) | cad > base]
m[, event := as.integer(!is.na(cad) & cad > base)]
m[, cens := pmin(fifelse(is.na(death), admin, death), admin, na.rm=TRUE)]
m[, endd := as.Date(fifelse(event==1, cad, cens))]
m[, time := as.numeric(endd - base)/365.25]
m <- m[time>0 & !is.na(carrier) & !is.na(age) & !is.na(sex)]
cat(sprintf("analysis N=%d ; prevalent excluded=%d ; incident CAD=%d ; carriers=%d ; carrier events=%d\n",
    nrow(m), prev, sum(m$event), sum(m$carrier), sum(m$event[m$carrier==1])))
fit <- coxph(Surv(time,event) ~ carrier + age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=m)
ci <- summary(fit)$conf.int["carrier",]; pv <- summary(fit)$coefficients["carrier","Pr(>|z|)"]
cat(sprintf("\n=== PDE3B pLoF carrier -> incident CAD (Cox PH) ===\nHR=%.3f  95%% CI %.3f-%.3f  P=%.3f\n",
    ci["exp(coef)"], ci["lower .95"], ci["upper .95"], pv))
fwrite(data.table(analysis_N=nrow(m), carriers=sum(m$carrier), incident_CAD=sum(m$event),
    carrier_CAD=sum(m$event[m$carrier==1]), HR=round(ci["exp(coef)"],3),
    CI_lo=round(ci["lower .95"],3), CI_hi=round(ci["upper .95"],3), P=signif(pv,3)), "pde3b_cox_cad.csv")
