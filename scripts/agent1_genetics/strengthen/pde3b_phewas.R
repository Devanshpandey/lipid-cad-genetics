suppressMessages(library(data.table))
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen"
PH  <- "/path/to/cad-genetics/results/agent1_genetics/phenotypes"
raw <- fread(file.path(OUT,"pde3b_plof.raw"))
vc  <- names(raw)[7:ncol(raw)]
# counted-allele freq per variant; fold to MINOR-allele dosage; keep rare pLoF (MAF<0.01) = REGENIE mask1
af <- sapply(vc, function(v) mean(raw[[v]], na.rm=TRUE)/2)
M  <- as.data.table(lapply(vc, function(v){ x<-raw[[v]]; a<-af[v]; if(is.na(a)) rep(0,length(x)) else if(a>0.5) 2-x else x })); setnames(M, vc)
maf  <- pmin(af, 1-af); maf[is.na(maf)] <- 1; keep <- vc[maf < 0.01]
cat("pLoF variants kept (MAF<0.01):", length(keep), "of", length(vc), "\n")
raw[, ncopy := rowSums(M[, ..keep], na.rm=TRUE)]
raw[, carrier := as.integer(ncopy >= 1)]
geno <- raw[, .(IID, carrier, ncopy)]
cat("exome samples:", nrow(geno)," pLoF carriers:", sum(geno$carrier),
    " carrier_freq:", signif(mean(geno$carrier),4), "\n")
cov <- fread(file.path(PH,"covariates.txt")); bin <- fread(file.path(PH,"pheno_binary.txt"))
qt  <- fread(file.path(PH,"pheno_quantitative.txt"))
d <- merge(geno, cov, by="IID"); d <- merge(d, bin[,-("FID")], by="IID", all.x=TRUE)
d <- merge(d, qt[,-("FID")], by="IID", all.x=TRUE)
ct <- paste(c("carrier","age","sex",paste0("PC",1:10)), collapse=" + ")
binp <- c("CAD","MI","REVASC","STROKE","HF","CV_DEATH","MACE","AF","PAD","STATIN_USE")
qtp  <- c("LDL_C","HDL_C","TRIGLY","TOT_CHOL","APOA1","APOB","LPA","nonHDL_C","eGFR","CRP","HBA1C")
res <- list()
for (p in binp){ dd<-d[!is.na(get(p))]
  f<-tryCatch(glm(as.formula(paste(p,"~",ct)),family=binomial,data=dd),error=function(e)NULL); if(is.null(f))next
  s<-summary(f)$coefficients; if(!("carrier"%in%rownames(s)))next
  res[[p]]<-data.table(pheno=p,type="binary",N=nrow(dd),n_case=sum(dd[[p]]==1,na.rm=TRUE),
    n_carrier_case=sum(dd$carrier==1&dd[[p]]==1,na.rm=TRUE),effect=exp(s["carrier","Estimate"]),
    effect_lab="OR",se=s["carrier","Std. Error"],P=s["carrier","Pr(>|z|)"]) }
for (p in qtp){ dd<-d[!is.na(get(p))]
  dd[, .y:=qnorm((frank(get(p),ties.method="average")-0.5)/nrow(dd))]
  f<-tryCatch(lm(as.formula(paste(".y ~",ct)),data=dd),error=function(e)NULL); if(is.null(f))next
  s<-summary(f)$coefficients; if(!("carrier"%in%rownames(s)))next
  res[[p]]<-data.table(pheno=p,type="quant(SD)",N=nrow(dd),n_case=NA,n_carrier_case=NA,
    effect=s["carrier","Estimate"],effect_lab="beta_SD",se=s["carrier","Std. Error"],P=s["carrier","Pr(>|t|)"]) }
R<-rbindlist(res); R[,P_bonf:=p.adjust(P,method="bonferroni")]; setorder(R,P)
fwrite(R, file.path(OUT,"pde3b_phewas.csv"))
cat("\n=== PDE3B pLoF-carrier PheWAS (sorted by P) ===\n")
print(R[,.(pheno,type,N,n_carrier_case,effect=round(effect,3),effect_lab,P=signif(P,3),P_bonf=signif(P_bonf,3))])
cat("\nBonferroni 0.05/",nrow(R)," = ",signif(0.05/nrow(R),3),"\n",sep="")
