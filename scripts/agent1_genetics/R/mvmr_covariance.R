#!/usr/bin/env Rscript
# ============================================================
# Two covariance problems, both raised in peer review.
#
# (1) CONDITIONAL F AND Q_A WERE COMPUTED WITH gencov = 0.
#     The four exposure GWAS were run on the SAME participants and the traits are
#     correlated, so the SNP-exposure effect estimates covary. Setting gencov=0
#     ignores that. MVMR::phenocov_mvmr() builds the required covariance from the
#     phenotypic covariance of the exposures; we recompute conditional F and Q_A with
#     it and report both so the effect of the assumption is visible.
#
# (2) THE CAD-vs-STROKE CONTRAST TREATED THE TWO ESTIMATES AS INDEPENDENT.
#     They share instruments and share controls, so the naive SE of the difference is
#     wrong. We bootstrap over instruments, recomputing BOTH outcome models on each
#     resample, which propagates the shared-instrument covariance directly.
#     This does NOT capture covariance induced by shared controls in the outcome
#     GWAS; that component is positive and would shrink the SE further, so the
#     bootstrap interval remains conservative. Stated, not assumed.
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(MendelianRandomization); library(MVMR); library(optparse)
})
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat_dir", type = "character"),
  make_option("--pheno",       type = "character"),
  make_option("--outcomes",    type = "character", default = "CAD,STROKE,MI,REVASC"),
  make_option("--nboot",       type = "integer",   default = 2000L),
  make_option("--cis_list",    type = "character", default = ""),
  make_option("--out",         type = "character")
)))

TR <- c("LDL_C","HDL_C","TRIGLY","LPA")
LPA_CHR <- "6"; LPA_LO <- 160500000L; LPA_HI <- 161700000L

lip <- lapply(TR, function(t) {
  d <- fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", t)))))
  d[, trait := t][]
}); names(lip) <- TR

# instrument selection: rank on LOG10P (pval underflows), distance-prune, exempt LPA cis
prune <- function(d, kb) {
  g <- d[LOG10P > -log10(5e-8) & !is.na(BETA) & !is.na(SE)][order(-LOG10P)]
  ks <- character(0); kp <- integer(0); kc <- character(0)
  for (i in seq_len(nrow(g))) {
    ch <- as.character(g$CHROM[i]); po <- as.integer(g$GENPOS[i]); sm <- kc == ch
    if (!length(kp) || !any(abs(kp[sm] - po) < kb * 1000L)) {
      ks <- c(ks, g$ID[i]); kp <- c(kp, po); kc <- c(kc, ch) }
  }
  ks
}
# Polygenic exposures: 10 Mb distance prune. Lp(a): the locus-dense cis instrument
# built by PLINK r^2 clumping (r^2<0.01 within 1 Mb), supplied as a file, because a
# 10 Mb rule collapses the ~1.1 Mb LPA locus to a single variant -- the defect this
# project already corrected. Without this the conditional F reported here would not
# match the primary analysis.
ivs <- unique(unlist(lapply(setdiff(TR, "LPA"), function(t) prune(lip[[t]], 10000))))
if (nzchar(opt$cis_list) && file.exists(opt$cis_list)) {
  cl <- fread(opt$cis_list, header = FALSE)[[1]]
  message(sprintf("[cov] cis-LPA instrument list: %d variants from %s", length(cl), opt$cis_list))
  ivs <- unique(c(ivs, cl))
} else stop("--cis_list is required so the Lp(a) instrument matches the primary analysis")
ref <- rbindlist(lapply(TR, function(t) lip[[t]][ID %in% ivs, .(ID, CHROM, GENPOS, LOG10P)]))
ref <- ref[ref[, .I[which.max(LOG10P)], by = ID]$V1][order(-LOG10P)]
incis <- ref[as.character(CHROM) == LPA_CHR & GENPOS >= LPA_LO & GENPOS <= LPA_HI, ID]
poly <- ref[!ID %in% incis]
kept <- character(0); kp <- integer(0); kc <- character(0)
for (i in seq_len(nrow(poly))) {
  ch <- as.character(poly$CHROM[i]); po <- as.integer(poly$GENPOS[i]); sm <- kc == ch
  if (!length(kp) || !any(abs(kp[sm] - po) < 10000 * 1000L)) {
    kept <- c(kept, poly$ID[i]); kp <- c(kp, po); kc <- c(kc, ch) }
}
ivs <- c(kept, incis)
message(sprintf("[cov] %d instruments (%d polygenic + %d cis-LPA)", length(ivs), length(kept), length(incis)))

al <- lip[[1]][ID %in% ivs, .(ID, A1 = ALLELE1, A0 = ALLELE0)]
mat <- function(field) {
  m <- matrix(NA_real_, length(ivs), length(TR), dimnames = list(ivs, TR))
  for (t in TR) {
    d <- merge(lip[[t]][ID %in% ivs], al, by = "ID")
    v <- if (field == "BETA") fifelse(d$ALLELE1 == d$A1, d$BETA, -d$BETA) else d$SE
    m[match(d$ID, ivs), t] <- v
  }
  m
}
bx <- mat("BETA"); bxse <- mat("SE")

# ---- phenotypic covariance of the exposures, for gencov ----
ph <- fread(opt$pheno)
# phenocov_mvmr() expects a phenotypic CORRELATION matrix: it forms the covariance of
# the SNP-exposure effects as Pcor[i,j] * se_i * se_j, so the matrix must be
# dimensionless. Passing a raw covariance matrix is wrong and catastrophically so here,
# because Lp(a) is in nmol/L (variance ~2400) while LDL-C is in mmol/L (variance ~1).
# The GWAS applied a rank-inverse-normal transform, so the relevant correlation is
# between the transformed traits; Spearman on the raw values is invariant to that
# monotone transform and is therefore the right estimate.
Pcov <- cor(as.matrix(ph[, ..TR]), use = "pairwise.complete.obs", method = "spearman")
message("[cov] exposure phenotypic correlation (Spearman, = Pearson on the RINT scale):")
print(round(Pcov, 3))

res <- list()
for (OUT in strsplit(opt$outcomes, ",")[[1]]) {
  od <- fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir, sprintf("outcomes_%s.txt.gz", OUT)))))
  od <- merge(od[ID %in% ivs], al, by = "ID")
  by_ <- rep(NA_real_, length(ivs)); bys <- rep(NA_real_, length(ivs))
  i <- match(od$ID, ivs)
  by_[i] <- fifelse(od$ALLELE1 == od$A1, od$BETA, -od$BETA); bys[i] <- od$SE
  ok <- complete.cases(bx) & complete.cases(bxse) & !is.na(by_) & !is.na(bys)
  BX <- bx[ok,,drop=FALSE]; BXSE <- bxse[ok,,drop=FALSE]; BY <- by_[ok]; BYSE <- bys[ok]

  f0 <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = BXSE, seBYG = BYSE, RSID = rownames(BX))
  s0 <- strength_mvmr(f0, gencov = 0)
  gc <- phenocov_mvmr(Pcov, BXSE)
  s1 <- strength_mvmr(f0, gencov = gc)
  q0 <- tryCatch(pleiotropy_mvmr(f0, gencov = 0)$Qpval,  error = function(e) NA_real_)
  q1 <- tryCatch(pleiotropy_mvmr(f0, gencov = gc)$Qpval, error = function(e) NA_real_)
  fit <- mr_mvivw(mr_mvinput(bx = BX, bxse = BXSE, by = BY, byse = BYSE))
  res[[OUT]] <- data.table(outcome = OUT, exposure = TR, n_iv = nrow(BX),
    OR = exp(fit@Estimate), b = fit@Estimate, se = fit@StdError,
    condF_gencov0 = as.numeric(s0[1, seq_along(TR)]),
    condF_phenocov = as.numeric(s1[1, seq_along(TR)]),
    QA_p_gencov0 = q0, QA_p_phenocov = q1)
  message(sprintf("[cov] %-7s condF gencov0 = %s | phenocov = %s", OUT,
    paste(sprintf("%.1f", as.numeric(s0[1, seq_along(TR)])), collapse="/"),
    paste(sprintf("%.1f", as.numeric(s1[1, seq_along(TR)])), collapse="/")))
}
out <- rbindlist(res)
fwrite(out, opt$out); print(out[exposure == "LPA"])

# ---- instrument bootstrap for the between-outcome contrasts ----
oget <- function(OUT) {
  od <- fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir, sprintf("outcomes_%s.txt.gz", OUT)))))
  od <- merge(od[ID %in% ivs], al, by = "ID")
  b <- rep(NA_real_, length(ivs)); s <- rep(NA_real_, length(ivs)); i <- match(od$ID, ivs)
  b[i] <- fifelse(od$ALLELE1 == od$A1, od$BETA, -od$BETA); s[i] <- od$SE
  list(b = b, s = s)
}
OU <- setNames(lapply(c("CAD","MI","REVASC","STROKE"), oget), c("CAD","MI","REVASC","STROKE"))
base_ok <- complete.cases(bx) & complete.cases(bxse) &
  Reduce(`&`, lapply(OU, function(o) !is.na(o$b) & !is.na(o$s)))
idx <- which(base_ok); k <- which(TR == "LPA")
lpa_b <- function(sel, OUT) {
  fit <- mr_mvivw(mr_mvinput(bx = bx[sel,,drop=FALSE], bxse = bxse[sel,,drop=FALSE],
                             by = OU[[OUT]]$b[sel], byse = OU[[OUT]]$s[sel]))
  fit@Estimate[k]
}
message(sprintf("[cov] bootstrapping %d resamples over %d shared instruments", opt$nboot, length(idx)))
OUTS <- c("CAD","MI","REVASC","STROKE")
B <- matrix(NA_real_, nrow = length(OUTS), ncol = opt$nboot, dimnames = list(OUTS, NULL))
for (bi in seq_len(opt$nboot)) {
  sel <- sample(idx, length(idx), replace = TRUE)
  for (o in OUTS) B[o, bi] <- tryCatch(lpa_b(sel, o), error = function(e) NA_real_)
}
B <- B[, colSums(is.na(B)) == 0, drop = FALSE]
message(sprintf("[cov] %d of %d bootstrap resamples usable", ncol(B), opt$nboot))
lpa_se <- function(OUT) {
  fit <- mr_mvivw(mr_mvinput(bx = bx[idx,,drop=FALSE], bxse = bxse[idx,,drop=FALSE],
                             by = OU[[OUT]]$b[idx], byse = OU[[OUT]]$s[idx]))
  fit@StdError[k]
}
se_stroke <- lpa_se("STROKE")
contr <- rbindlist(lapply(c("CAD","MI","REVASC"), function(o) {
  d <- B[o, ] - B["STROKE", ]
  obs <- lpa_b(idx, o) - lpa_b(idx, "STROKE")
  naive <- sqrt(lpa_se(o)^2 + se_stroke^2)   # the independence assumption under review
  data.table(contrast = paste(o, "vs STROKE"), diff_logOR = obs,
             se_naive_independent = naive, se_bootstrap = sd(d),
             corr_between_estimates = cor(B[o,], B["STROKE",]),
             z_naive = obs/naive, z_bootstrap = obs/sd(d),
             p_bootstrap = 2*pnorm(-abs(obs/sd(d))))
}))
print(contr); fwrite(contr, sub("\\.csv$", "_contrasts.csv", opt$out))
message("[cov] done")
