#!/usr/bin/env Rscript
# MR-BMA + ApoB-vs-LDL discordant-instrument scan (Section 2.4).
# Memory-safe: awk stream-prefilters each GWAS so R holds only small subsets.
# Effect allele (ALLELE1) is identical across all REGENIE GWAS (same genotype file);
# we still verify allele concordance defensively.
suppressPackageStartupMessages({library(data.table)})
set.seed(42)

DIR <- "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen"
EXP <- c("LDL_C","APOB","HDL_C","APOA1","TRIGLY","LPA")
P_GWS <- 5e-8; KB <- 10000L*1000L
IDF <- file.path(OUT, "mrbma_instruments.txt")
lf  <- function(tr) file.path(DIR, sprintf("lipids_%s.txt.gz", tr))

# stream-read only GW-sig rows (LOG10P col 12 > 7.30103)
gws_read <- function(f) fread(cmd = paste0("gunzip -c ", shQuote(f),
  " | awk 'NR==1 || $12+0 > 7.30103'"), na.strings = "NA")
# stream-read only rows whose ID (col 3) is in idfile
inst_read <- function(f, idfile) fread(cmd = paste0("gunzip -c ", shQuote(f),
  " | awk 'NR==FNR{a[$1];next} FNR==1||($3 in a)' ", shQuote(idfile), " -"), na.strings = "NA")

cat("[1] PASS 1: GW-sig extraction per exposure\n")
gws <- list()
for (tr in EXP) {
  d <- gws_read(lf(tr)); d[, pval := 10^(-LOG10P)]
  gws[[tr]] <- d[!is.na(BETA) & !is.na(SE), .(ID,CHROM,GENPOS,ALLELE0,ALLELE1,BETA,SE,pval)]
  cat(sprintf("    %s: %d GW-sig SNPs\n", tr, nrow(gws[[tr]]))); rm(d); gc()
}

cat("[2] instrument selection (union + 10 Mb greedy clump)\n")
allref <- unique(rbindlist(gws)[, .(ID,CHROM,GENPOS,ALLELE0,ALLELE1,pval)])
allref <- allref[allref[, .I[which.min(pval)], by=ID]$V1][order(pval)]
fin <- character(0); fc <- character(0); fp <- integer(0)
for (i in seq_len(nrow(allref))) {
  ch <- as.character(allref$CHROM[i]); po <- as.integer(allref$GENPOS[i]); s <- fc==ch
  if (!length(fp) || !any(s & abs(fp[s]-po) < KB)) { fin<-c(fin,allref$ID[i]); fc<-c(fc,ch); fp<-c(fp,po) }
}
cat(sprintf("    %d independent instruments\n", length(fin)))
writeLines(fin, IDF)
refA <- allref[ID %in% fin, .(ID, rA1=ALLELE1, rA0=ALLELE0)]; setkey(refA, ID)

cat("[3] PASS 2: extract instrument rows per exposure + CAD\n")
getcol <- function(f) {
  d <- inst_read(f, IDF); setkey(d, ID); d <- d[refA][, .(ID, ALLELE1, BETA, SE, rA1, rA0)]
  b <- fifelse(d$ALLELE1==d$rA1, d$BETA, fifelse(d$ALLELE1==d$rA0, -d$BETA, NA_real_))
  list(id=d$ID, b=setNames(b, d$ID), se=setNames(d$SE, d$ID))
}
bx <- matrix(NA_real_, length(fin), length(EXP), dimnames=list(fin, EXP)); sx <- bx
for (tr in EXP) { g <- getcol(lf(tr)); bx[g$id,tr] <- g$b[g$id]; sx[g$id,tr] <- g$se[g$id]
  cat(sprintf("    %s extracted\n", tr)) }
oc <- getcol(file.path(DIR,"outcomes_CAD.txt.gz"))
by <- setNames(rep(NA_real_, length(fin)), fin); sy <- by
by[oc$id] <- oc$b[oc$id]; sy[oc$id] <- oc$se[oc$id]

ok <- complete.cases(bx) & complete.cases(sx) & !is.na(by) & !is.na(sy)
bx <- bx[ok,,drop=FALSE]; sx <- sx[ok,,drop=FALSE]; by <- by[ok]; sy <- sy[ok]
n <- nrow(bx); cat(sprintf("    %d instruments with complete data\n", n))

# ================== DISCORDANT SCAN (ApoB vs LDL-C) ==================
cat("\n===== DISCORDANT SCAN: ApoB vs LDL-C =====\n")
bl <- bx[,"LDL_C"]; ba <- bx[,"APOB"]
r <- cor(bl, ba); R2 <- summary(lm(ba~bl))$r.squared; rsd <- sd(residuals(lm(ba~0+bl)))
disc <- data.table(ID=rownames(bx), bLDL=bl, bApoB=ba, zLDL=bl/sx[,"LDL_C"], zApoB=ba/sx[,"APOB"],
                   stdres=as.numeric(scale(residuals(lm(ba~bl)))))[order(-abs(stdres))]
n_disc <- sum(abs(disc$stdres) > 2)
cat(sprintf("Pearson r(betaApoB,betaLDL) = %.4f\n", r))
cat(sprintf("R^2 (ApoB~LDL) = %.4f ; independent variance 1-R^2 = %.4f\n", R2, 1-R2))
cat(sprintf("residual SD (ApoB|LDL) = %.4f ; |std resid|>2 discordant = %d / %d\n", rsd, n_disc, n))
print(disc[1:10, .(ID, bLDL=round(bLDL,3), bApoB=round(bApoB,3), zLDL=round(zLDL,1),
                   zApoB=round(zApoB,1), stdres=round(stdres,2))])
fwrite(disc, file.path(OUT,"discordant_apob_ldl.csv"))

# ================== MR-BMA (exact 2^k) ==================
cat("\n===== MR-BMA =====\n")
k <- ncol(bx); sig2 <- 0.25; pin <- 0.5
yw <- by/sy; Xs <- sweep(bx/sy, 2, apply(bx,2,sd), "/")
M <- as.matrix(expand.grid(rep(list(0:1), k))); colnames(M) <- EXP
lml <- numeric(nrow(M)); pmu <- matrix(0, nrow(M), k, dimnames=list(NULL,EXP))
for (m in seq_len(nrow(M))) {
  cols <- which(M[m,]==1)
  if (!length(cols)) { lml[m] <- -0.5*(n*log(2*pi)+sum(yw^2)); next }
  Xm <- Xs[,cols,drop=FALSE]; Sig <- diag(n) + sig2*(Xm%*%t(Xm))
  ch <- chol(Sig); lml[m] <- -0.5*(n*log(2*pi) + 2*sum(log(diag(ch))) + sum(backsolve(ch,yw,transpose=TRUE)^2))
  V <- solve(crossprod(Xm)+diag(1/sig2,length(cols))); pmu[m,cols] <- as.numeric(V%*%crossprod(Xm,yw))
}
lp <- lml + rowSums(M)*log(pin) + (k-rowSums(M))*log(1-pin); pp <- exp(lp-max(lp)); pp <- pp/sum(pp)
res <- data.table(exposure=EXP, MIP=round(colSums(M*pp),3), MA_effect=round(colSums(pmu*pp),4))[order(-MIP)]
cat("MIP + model-averaged effects:\n"); print(res)
cat("\nTop 5 models:\n"); for (i in order(-pp)[1:5]) {
  ex <- EXP[which(M[i,]==1)]; if(!length(ex)) ex<-"(null)"; cat(sprintf("  PP=%.3f : %s\n", pp[i], paste(ex,collapse=" + "))) }
fwrite(res, file.path(OUT,"mrbma_mip.csv"))
cat("\nUnivariable IVW per exposure:\n"); w <- 1/sy^2
for (j in seq_len(k)) { b <- sum(w*bx[,j]*by)/sum(w*bx[,j]^2); se <- sqrt(1/sum(w*bx[,j]^2))
  cat(sprintf("  %-7s beta=%+.3f se=%.3f P=%.1e\n", EXP[j], b, se, 2*pnorm(-abs(b/se)))) }
cat("\nDONE\n")
