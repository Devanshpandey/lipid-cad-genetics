#!/usr/bin/env Rscript
# ============================================================================
# MR-BMA: Bayesian model averaging over correlated lipid exposures.
#
# WHY THIS EXISTS. Figure 1d and Supplementary Table S17 report marginal
# inclusion probabilities for six lipid fractions against CAD, but no code for
# that analysis survived anywhere in the project -- the result existed only as a
# CSV. This reimplements it so the panel is reproducible, and re-derives the
# instruments with the corrected distance-pruning rule.
#
# METHOD (Zuber et al., Nat Commun 2020). Instruments are weighted by the
# inverse of the outcome standard error, which puts the regression on a unit
# error-variance scale, and multivariable IVW is a regression through the
# origin. For a model M (a subset of the exposures) with design matrix X_M and
# a normal prior beta ~ N(0, sigma^2 I), the marginal likelihood has the
# closed form
#
#   log ML(M) = -0.5 [ k log sigma^2 + log|X'X + sigma^-2 I| ]
#               -0.5 [ y'y - y'X (X'X + sigma^-2 I)^-1 X'y ]
#
# with k = |M|. Model priors are independent-inclusion: p(M) = p^k (1-p)^(K-k).
# The marginal inclusion probability of an exposure is the summed posterior
# probability of the models containing it. All 2^K - 1 non-empty models are
# enumerated exactly (K = 6, so 63 models); nothing is sampled.
#
# The prior pair (sigma, p) is swept, because MIPs for collinear exposures are
# prior-sensitive and the manuscript reports that sensitivity rather than a
# single grid cell.
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(optparse)})
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat_dir", type = "character"),
  make_option("--outcome",     type = "character", default = "CAD"),
  make_option("--exposures",   type = "character",
              default = "LDL_C,APOB,HDL_C,APOA1,TRIGLY,LPA"),
  make_option("--p_thresh",    type = "double",  default = 5e-8),
  make_option("--kb_clump",    type = "integer", default = 10000L),
  make_option("--out",         type = "character")
)))
TR <- strsplit(opt$exposures, ",")[[1]]; K <- length(TR)

# ---- instrument selection: identical rule to run_mvmr.R, post-bugfix --------
clump <- function(d, kb) {
  g <- d[LOG10P > -log10(opt$p_thresh) & !is.na(BETA) & !is.na(SE)][order(-LOG10P)]
  ks <- character(0); kp <- integer(0); kc <- character(0)
  for (i in seq_len(nrow(g))) {
    ch <- as.character(g$CHROM[i]); po <- as.integer(g$GENPOS[i]); sm <- kc == ch
    if (!length(kp) || !any(abs(kp[sm] - po) < kb * 1000L)) {
      ks <- c(ks, g$ID[i]); kp <- c(kp, po); kc <- c(kc, ch) }
  }
  ks
}
message("[bma] loading exposures")
lip <- lapply(TR, function(t)
  fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir,
        sprintf("lipids_%s.txt.gz", t))))))
names(lip) <- TR
ivs <- unique(unlist(lapply(TR, function(t) clump(lip[[t]], opt$kb_clump))))
ref <- rbindlist(lapply(TR, function(t) lip[[t]][ID %in% ivs, .(ID, CHROM, GENPOS, LOG10P)]))
ref <- ref[ref[, .I[which.max(LOG10P)], by = ID]$V1][order(-LOG10P)]
keep <- character(0); kp <- integer(0); kc <- character(0)
for (i in seq_len(nrow(ref))) {
  ch <- as.character(ref$CHROM[i]); po <- as.integer(ref$GENPOS[i]); sm <- kc == ch
  if (!length(kp) || !any(abs(kp[sm] - po) < opt$kb_clump * 1000L)) {
    keep <- c(keep, ref$ID[i]); kp <- c(kp, po); kc <- c(kc, ch) }
}
message(sprintf("[bma] %d instruments after the union prune", length(keep)))

al <- lip[[1]][ID %in% keep, .(ID, A1 = ALLELE1, A0 = ALLELE0)]
X <- matrix(NA_real_, length(keep), K, dimnames = list(keep, TR))
for (t in TR) {
  d <- merge(lip[[t]][ID %in% keep], al, by = "ID")
  X[match(d$ID, keep), t] <- fifelse(d$ALLELE1 == d$A1, d$BETA,
                              fifelse(d$ALLELE1 == d$A0, -d$BETA, NA_real_))
}
od <- fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir,
      sprintf("outcomes_%s.txt.gz", opt$outcome)))))
od <- merge(od[ID %in% keep], al, by = "ID")
by_ <- rep(NA_real_, length(keep)); bys <- rep(NA_real_, length(keep))
i <- match(od$ID, keep)
by_[i] <- fifelse(od$ALLELE1 == od$A1, od$BETA, fifelse(od$ALLELE1 == od$A0, -od$BETA, NA_real_))
bys[i] <- od$SE
ok <- complete.cases(X) & !is.na(by_) & !is.na(bys)
X <- X[ok, , drop = FALSE]; by_ <- by_[ok]; bys <- bys[ok]
J <- nrow(X)
message(sprintf("[bma] %d instruments complete across %d exposures + outcome", J, K))

# ---- weight to unit error variance -----------------------------------------
Xw <- X / bys; yw <- by_ / bys
yty <- sum(yw^2)

log_ml <- function(idx, s2) {
  k <- length(idx)
  Xm <- Xw[, idx, drop = FALSE]
  A  <- crossprod(Xm) + diag(1 / s2, k)
  ch <- chol(A)
  logdet <- 2 * sum(log(diag(ch)))
  Xty <- crossprod(Xm, yw)
  quad <- sum(backsolve(ch, Xty, transpose = TRUE)^2)
  -0.5 * (k * log(s2) + logdet) - 0.5 * (yty - quad)
}

models <- do.call(c, lapply(1:K, function(k) combn(K, k, simplify = FALSE)))
grid <- expand.grid(sigma = c(0.25, 0.5, 1.0), prior_inclusion = c(0.25, 0.5))
res <- rbindlist(lapply(seq_len(nrow(grid)), function(g) {
  s2 <- grid$sigma[g]^2; p <- grid$prior_inclusion[g]
  lp <- vapply(models, function(m)
        log_ml(m, s2) + length(m) * log(p) + (K - length(m)) * log(1 - p), numeric(1))
  w <- exp(lp - max(lp)); w <- w / sum(w)
  mip <- vapply(seq_len(K), function(j) sum(w[vapply(models, function(m) j %in% m, logical(1))]),
                numeric(1))
  top <- models[[which.max(w)]]
  data.table(prior_sigma = grid$sigma[g], prior_inclusion = p, n_instruments = J,
             t(setNames(round(mip, 3), paste0("MIP_", TR))),
             best_model = paste(TR[top], collapse = "+"), best_model_PP = round(max(w), 3))
}))
print(res)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, opt$out)
message("[bma] wrote ", opt$out)
