suppressMessages({library(data.table); library(MendelianRandomization)})
has_presso <- requireNamespace("MRPRESSO", quietly=TRUE)
MERGED <- "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen/mr_robustness"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
EXP <- c(LDL_C="lipids_LDL_C", HDL_C="lipids_HDL_C", TRIGLY="lipids_TRIGLY", LPA="lipids_LPA")
rd <- function(f){ d<-fread(cmd=paste("gunzip -c", file.path(MERGED, paste0(f,".txt.gz"))))
  d[, .(SNP=ID, chr=CHROM, pos=GENPOS, ea=toupper(ALLELE1), oa=toupper(ALLELE0),
        eaf=A1FREQ, beta=BETA, se=SE, p=10^(-LOG10P))] }
cad <- rd("outcomes_CAD")
amb <- function(e,o)(e=="A"&o=="T")|(e=="T"&o=="A")|(e=="C"&o=="G")|(e=="G"&o=="C")
res <- list()
for(nm in names(EXP)){
  ex <- rd(EXP[[nm]]); sig <- ex[p<5e-8][order(p)]; sig<-sig[!amb(ea,oa)]
  kb<-10000L; keep<-c(); uc<-c(); up<-c()
  for(i in seq_len(nrow(sig))){ r<-sig[i]; if(!any(uc==r$chr & abs(up-r$pos)<kb*1000)){keep<-c(keep,r$SNP);uc<-c(uc,r$chr);up<-c(up,r$pos)} }
  iv <- ex[SNP %in% keep]
  m <- merge(iv[,.(SNP,ea,oa,bx=beta,sx=se)], cad[,.(SNP,cea=ea,by=beta,sy=se)], by="SNP")
  m[, by := fifelse(cea==ea, by, -by)]; m<-m[!is.na(by)&!is.na(sy)]
  inp <- mr_input(bx=m$bx, bxse=m$sx, by=m$by, byse=m$sy)
  ivw<-mr_ivw(inp); eg<-mr_egger(inp); wm<-mr_median(inp, weighting="weighted")
  steiger_frac <- mean(abs(m$bx/m$sx) > abs(m$by/m$sy))   # exposure-primary fraction
  loo <- sapply(seq_len(nrow(m)), function(i) mr_ivw(mr_input(bx=m$bx[-i],bxse=m$sx[-i],by=m$by[-i],byse=m$sy[-i]))@Estimate)
  pg<-NA; pcor<-NA
  if(has_presso){ df<-data.frame(bx=m$bx,sx=m$sx,by=m$by,sy=m$sy)
    pr<-tryCatch(MRPRESSO::mr_presso(BetaOutcome="by",BetaExposure="bx",SdOutcome="sy",SdExposure="sx",
        data=df,OUTLIERtest=TRUE,DISTORTIONtest=FALSE,NbDistribution=1000,SignifThreshold=0.05),error=function(e)NULL)
    if(!is.null(pr)){ pg<-pr[["MR-PRESSO results"]][["Global Test"]]$Pvalue
      est<-pr[["Main MR results"]]; pcor<-exp(est[["Causal Estimate"]][nrow(est)]) } }
  res[[nm]]<-data.table(exposure=nm,nSNP=nrow(m),
    IVW_OR=round(exp(ivw@Estimate),3),IVW_P=signif(ivw@Pvalue,3),
    Egger_OR=round(exp(eg@Estimate),3),Egger_intP=signif(eg@Pvalue.Int,3),
    WtdMed_OR=round(exp(wm@Estimate),3),WtdMed_P=signif(wm@Pvalue,3),
    Steiger_frac=round(steiger_frac,3),
    LOO_OR_min=round(exp(min(loo)),3),LOO_OR_max=round(exp(max(loo)),3),
    PRESSO_globalP=as.character(pg),PRESSO_corrOR=round(as.numeric(pcor),3))
  cat(nm,"done\n") }
R<-rbindlist(res,fill=TRUE); fwrite(R, file.path(OUT,"mr_robustness_CAD.csv")); print(R)
cat("MRPRESSO available:", has_presso, "\n")
