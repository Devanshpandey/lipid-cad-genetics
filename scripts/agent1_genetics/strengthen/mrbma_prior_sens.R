#!/usr/bin/env Rscript
# MR-BMA prior-sensitivity: MIPs across prior inclusion pi and prior sd sigma.
# Reuses the saved 345-instrument list (fast PASS-2-only extraction).
suppressPackageStartupMessages({library(data.table)})
set.seed(42)
DIR <- "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen"
EXP <- c("LDL_C","APOB","HDL_C","APOA1","TRIGLY","LPA")
IDF <- file.path(OUT, "mrbma_instruments.txt")
COLS<- c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","BETA","SE","LOG10P")
lf  <- function(tr) file.path(DIR, sprintf("lipids_%s.txt.gz", tr))
inst_read <- function(f) fread(cmd = paste0("gunzip -c ", shQuote(f),
  " | awk 'NR==FNR{a[$1];next} FNR==1||($3 in a)' ", shQuote(IDF), " -"),
  select=COLS, na.strings="NA")

fin <- readLines(IDF); cat(sprintf("[1] %d instruments\n", length(fin)))
# reference alleles from LDL_C
ldl <- inst_read(lf("LDL_C")); setkey(ldl, ID)
refA <- ldl[fin, .(ID, rA1=ALLELE1, rA0=ALLELE0)]; setkey(refA, ID)
getcol <- function(f) {
  d <- inst_read(f); setkey(d, ID); d <- d[refA][, .(ID, ALLELE1, BETA, SE, rA1, rA0)]
  b <- fifelse(d$ALLELE1==d$rA1, d$BETA, fifelse(d$ALLELE1==d$rA0, -d$BETA, NA_real_))
  list(id=d$ID, b=setNames(b,d$ID), se=setNames(d$SE,d$ID))
}
bx <- matrix(NA_real_, length(fin), length(EXP), dimnames=list(fin,EXP)); sx <- bx
for (tr in EXP) { g <- getcol(lf(tr)); bx[g$id,tr] <- g$b[g$id]; sx[g$id,tr] <- g$se[g$id] }
oc <- getcol(file.path(DIR,"outcomes_CAD.txt.gz"))
by <- setNames(rep(NA_real_,length(fin)),fin); sy <- by; by[oc$id]<-oc$b[oc$id]; sy[oc$id]<-oc$se[oc$id]
ok <- complete.cases(bx)&complete.cases(sx)&!is.na(by)&!is.na(sy)
bx<-bx[ok,,drop=FALSE]; sx<-sx[ok,,drop=FALSE]; by<-by[ok]; sy<-sy[ok]; n<-nrow(bx)
saveRDS(list(bx=bx,sx=sx,by=by,sy=sy), file.path(OUT,"mrbma_harmonized.rds"))
cat(sprintf("[2] %d complete instruments; cached RDS\n", n))

k <- ncol(bx); M <- as.matrix(expand.grid(rep(list(0:1),k))); colnames(M)<-EXP
bma <- function(sig2, pin) {
  yw <- by/sy; Xs <- sweep(bx/sy, 2, apply(bx,2,sd), "/")
  lml <- numeric(nrow(M))
  for (m in seq_len(nrow(M))) {
    cols <- which(M[m,]==1)
    if(!length(cols)){lml[m]<- -0.5*(n*log(2*pi)+sum(yw^2));next}
    Xm<-Xs[,cols,drop=FALSE]; Sig<-diag(n)+sig2*(Xm%*%t(Xm)); ch<-chol(Sig)
    lml[m]<- -0.5*(n*log(2*pi)+2*sum(log(diag(ch)))+sum(backsolve(ch,yw,transpose=TRUE)^2))
  }
  lp<-lml+rowSums(M)*log(pin)+(k-rowSums(M))*log(1-pin); pp<-exp(lp-max(lp)); pp<-pp/sum(pp)
  setNames(round(colSums(M*pp),3), EXP)
}
cat("\n===== MR-BMA prior sensitivity (MIP) =====\n")
grid <- expand.grid(sigma=c(0.5,0.25), pi=c(0.5,0.25,0.1))
tab <- rbindlist(lapply(seq_len(nrow(grid)), function(i){
  mip <- bma(grid$sigma[i]^2, grid$pi[i])
  data.table(sigma=grid$sigma[i], pi=grid$pi[i], t(mip))
}))
print(tab)
fwrite(tab, file.path(OUT,"mrbma_prior_sensitivity.csv"))
cat("\nDONE\n")
