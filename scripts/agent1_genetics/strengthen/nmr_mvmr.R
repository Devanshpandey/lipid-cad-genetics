#!/usr/bin/env Rscript
## NMR LDL particle-number (LDL_P) vs LDL cholesterol (LDL_C_NMR): univariable + MVMR on CAD.
## Outcome = CARDIoGRAMplusC4D (external, no UKB overlap). Exposures RINT-scaled => per-SD effects.
## Usage: Rscript nmr_mvmr.R <LDL_P.regenie.gz> <LDL_C_NMR.regenie.gz> <cardiogram.txt> <outdir>
suppressMessages({library(data.table); library(MendelianRandomization)})
a <- commandArgs(TRUE); fp_lp<-a[1]; fp_lc<-a[2]; fp_cad<-a[3]; outdir<-a[4]
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f){ d<-fread(f)
  d[, .(SNP=ID, chr=CHROM, pos=GENPOS, ea=toupper(ALLELE1), oa=toupper(ALLELE0),
        beta=BETA, se=SE, chisq=CHISQ, p=10^(-LOG10P))] }
lp<-rd(fp_lp); lc<-rd(fp_lc)
lam<-function(d) median(d$chisq,na.rm=TRUE)/0.4549
cat(sprintf("[QC] lambda_GC  LDL_P=%.3f  LDL_C_NMR=%.3f\n", lam(lp), lam(lc)))
cat(sprintf("[QC] gw-sig(P<5e-8)  LDL_P=%d  LDL_C_NMR=%d\n", sum(lp$p<5e-8), sum(lc$p<5e-8)))
## instruments: union of gw-sig, greedy 10Mb clump by min-P
sig<-rbind(lp[p<5e-8,.(SNP,chr,pos,p)], lc[p<5e-8,.(SNP,chr,pos,p)])
sig<-sig[order(p)][!duplicated(SNP)]
kb<-10000L; keep<-character(0); uc<-integer(0); up<-numeric(0)
for(i in seq_len(nrow(sig))){ r<-sig[i]
  if(!any(uc==r$chr & abs(up-r$pos)<kb*1000)){ keep<-c(keep,r$SNP); uc<-c(uc,r$chr); up<-c(up,r$pos) } }
ivs<-keep
cat(sprintf("[MVMR] %d union gw-sig -> %d independent instruments\n", nrow(sig), length(ivs)))
ref<-lp[SNP %in% ivs, .(SNP, rea=ea, roa=oa)]
amb<-function(e,o)(e=="A"&o=="T")|(e=="T"&o=="A")|(e=="C"&o=="G")|(e=="G"&o=="C")
ref<-ref[!amb(rea,roa)]  ## drop strand-ambiguous (palindromic) instruments
alx<-function(d){ m<-merge(ref, d[SNP %in% ivs], by="SNP")
  m[, b:=fifelse(ea==rea,beta,fifelse(ea==roa,-beta,NA_real_))]; m[,.(SNP,b,se)] }
xp<-alx(lp); xc<-alx(lc)
cad<-fread(fp_cad)
cad<-cad[, .(SNP=markername, ea=toupper(effect_allele), oa=toupper(noneffect_allele),
             beta=as.numeric(beta), se=as.numeric(se_dgc))]
co<-merge(ref, cad, by="SNP")
co[, yb:=fifelse(ea==rea,beta,fifelse(ea==roa,-beta,NA_real_))]
M<-Reduce(function(x,y) merge(x,y,by="SNP"),
   list(xp[,.(SNP,bxp=b,sxp=se)], xc[,.(SNP,bxc=b,sxc=se)], co[,.(SNP,by=yb,sy=se)]))
M<-M[complete.cases(M)]
cat(sprintf("[MVMR] %d SNPs matched across LDL_P, LDL_C_NMR, CAD\n", nrow(M)))
fwrite(M, file.path(outdir,"nmr_mvmr_snps.csv"))
uni<-function(bx,sx){ o<-mr_ivw(mr_input(bx=bx,bxse=sx,by=M$by,byse=M$sy))
  data.table(OR=exp(o@Estimate),lo=exp(o@CILower),hi=exp(o@CIUpper),P=o@Pvalue,nsnp=length(bx)) }
res<-rbind(cbind(model="univariable",exposure="LDL_P",uni(M$bxp,M$sxp)),
           cbind(model="univariable",exposure="LDL_C_NMR",uni(M$bxc,M$sxc)))
mvi<-mr_mvinput(bx=cbind(M$bxp,M$bxc), bxse=cbind(M$sxp,M$sxc), by=M$by, byse=M$sy)
mv<-mr_mvivw(mvi)
res<-rbind(res, data.table(model="MVMR",exposure=c("LDL_P","LDL_C_NMR"),
   OR=exp(mv@Estimate),lo=exp(mv@CILower),hi=exp(mv@CIUpper),P=mv@Pvalue,nsnp=nrow(M)))
if(requireNamespace("MVMR",quietly=TRUE)){
  fm<-tryCatch(MVMR::format_mvmr(BXGs=cbind(M$bxp,M$bxc),BYG=M$by,seBXGs=cbind(M$sxp,M$sxc),seBYG=M$sy,RSID=M$SNP),error=function(e)NULL)
  if(!is.null(fm)){ sr<-tryCatch(MVMR::strength_mvmr(fm,gencov=0),error=function(e)NULL)
    if(!is.null(sr)) cat("[MVMR] conditional F (LDL_P / LDL_C_NMR):", paste(round(as.numeric(sr),1),collapse=" / "),"\n") } }
fwrite(res, file.path(outdir,"nmr_mvmr_results.csv"))
cat("\n===== RESULTS: OR per SD on CARDIoGRAM CAD =====\n"); print(res)
