#!/usr/bin/env Rscript
# ============================================================
# Tobit check of the censored-block Lp(a) instrument effects, and the one raw-unit result in the manuscript.
#
# Tied blocks lose the variation inside each tail, so simulation leaves ~3% attenuation. A censored-normal likelihood
# (survreg, gaussian, interval2) on the same rank-normal score has no such loss, which sizes the remaining
# bias on real data. Fits per cis variant, with the 16 GWAS covariates and no LOCO term:
#   ols_meas  OLS on LPA_MEAS (measured only)         per residual SD   -> published estimand
#   ols_cens  OLS on LPA_CENS (tied blocks)           per residual SD   -> proposed estimand
#   tobit_z   Tobit on LPA_CENS score, censored at the block quantile bounds, per latent residual SD
#   tobit_nmol, tobit_lognmol  Tobit on raw and log nmol/L, censored at 3.8 and 189
# The nmol/L Gaussian is misspecified (Lp(a) is strongly right-skewed); those rows are exploratory.
#
# Design control: Table S34 reports Lp(a) score x statin -> Lp(a) as +2.64 nmol/L from OLS on measured raw values.
# Refit in the same model with censoring respected.
# Aggregate results only.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(parallel); library(optparse) })
set.seed(42)
opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",   type = "character", help = "lpa_tobit_input.tsv"),
  make_option("--dose",    type = "character", help = "dose_LPA.raw (PLINK --export A)"),
  make_option("--weights", type = "character", help = "lpa_score.txt"),
  make_option("--covar",   type = "character"),
  make_option("--statin",  type = "character", help = "TSV IID statin"),
  make_option("--cores",   type = "integer", default = 32L),
  make_option("--out",     type = "character", help = "Output prefix")
)))
ph <- fread(opt$input, colClasses = c(IID = "character"))
cv <- fread(opt$covar, colClasses = c(FID = "character", IID = "character"))
g  <- fread(opt$dose, colClasses = c(IID = "character"))
w  <- fread(opt$weights, header = FALSE, col.names = c("SNP", "EA", "BETA"))

dose_cols <- setdiff(names(g), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
map <- data.table(col = dose_cols, SNP = sub("_[^_]+$", "", dose_cols), counted = toupper(sub("^.*_", "", dose_cols)))
map <- merge(map, w, by = "SNP"); stopifnot(nrow(map) == nrow(w))

d <- merge(ph, cv[, !"FID"], by = "IID")
d <- d[!is.na(side)]
covs <- setdiff(names(cv), c("FID", "IID"))
message(sprintf("[tobit] N=%d measured=%d low=%d high=%d, %d variants", nrow(d), d[side == 0L, .N], d[side == -1L, .N], d[side == 1L, .N], nrow(map)))

p_lo <- d[, mean(side == -1L)]; p_hi <- d[, mean(side == 1L)]
a <- qnorm(p_lo); b <- qnorm(p_hi, lower.tail = FALSE)
d[, `:=`(zL = fifelse(side == -1L, NA_real_, LPA_CENS), zR = fifelse(side == 1L, NA_real_, LPA_CENS))]
d[side == -1L, zR := a]; d[side == 1L, zL := b]
d[, `:=`(nL = fifelse(side == -1L, NA_real_, fifelse(side == 1L, 189, raw)), nR = fifelse(side == 1L, NA_real_, fifelse(side == -1L, 3.8, raw)))]
d[, `:=`(lL = log(nL), lR = log(nR))]

rhs <- paste(covs, collapse = " + ")
tobit <- function(L, R, x = NULL, dat = d) {
  f <- as.formula(sprintf("Surv(%s, %s, type = 'interval2') ~ %s%s", L, R, if (is.null(x)) "" else paste0(x, " + "), rhs))
  survreg(f, data = dat, dist = "gaussian")
}
# per-SD comparators divide by the covariate-only residual SD, which is how REGENIE scales a quantitative trait
sd_meas <- summary(lm(as.formula(paste("LPA_MEAS ~", rhs)), data = d))$sigma
sd_cens <- summary(lm(as.formula(paste("LPA_CENS ~", rhs)), data = d))$sigma
sc_z <- tobit("zL", "zR")$scale
message(sprintf("[tobit] residual SD: LPA_MEAS %.4f, LPA_CENS %.4f; Tobit latent scale %.4f", sd_meas, sd_cens, sc_z))

g <- g[match(d$IID, IID)]
fit_one <- function(i) {
  x <- g[[map$col[i]]]
  if (map$counted[i] != toupper(map$EA[i])) x <- 2 - x
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  dd <- copy(d[, c("LPA_MEAS", "LPA_CENS", "zL", "zR", "nL", "nR", "lL", "lR", covs), with = FALSE]); dd[, G := x]
  o <- function(y) { cf <- summary(lm(as.formula(paste(y, "~ G +", rhs)), data = dd))$coefficients["G", ]; cf[1:2] }
  t <- function(L, R) { s <- summary(tobit(L, R, "G", dd))$table["G", ]; s[1:2] }
  om <- o("LPA_MEAS"); oc <- o("LPA_CENS"); tz <- t("zL", "zR"); tn <- t("nL", "nR"); tl <- t("lL", "lR")
  data.table(SNP = map$SNP[i], EA = map$EA[i], published_weight = map$BETA[i], eaf = mean(x) / 2,
             ols_meas = om[1] / sd_meas, ols_meas_se = om[2] / sd_meas,
             ols_cens = oc[1] / sd_cens, ols_cens_se = oc[2] / sd_cens,
             tobit_z = tz[1] / sc_z, tobit_z_se = tz[2] / sc_z,
             tobit_nmol = tn[1], tobit_nmol_se = tn[2], tobit_lognmol = tl[1], tobit_lognmol_se = tl[2])
}
res <- rbindlist(mclapply(seq_len(nrow(map)), fit_one, mc.cores = opt$cores, mc.set.seed = FALSE))
res <- res[order(-abs(ols_cens))]
fwrite(res, paste0(opt$out, "_variants.csv"))
slope <- function(y, x) sum(x * y) / sum(x * x)
sm <- c(sprintf("%d cis variants, N=%d", nrow(res), nrow(d)),
        sprintf("slope on ols_meas: ols_cens %.4f, tobit_z %.4f", res[, slope(ols_cens, ols_meas)], res[, slope(tobit_z, ols_meas)]),
        sprintf("slope ols_cens on tobit_z: %.4f (residual attenuation of the tied blocks)", res[, slope(ols_cens, tobit_z)]),
        sprintf("median ratio tobit_z / ols_cens: %.4f", res[, median(tobit_z / ols_cens)]),
        sprintf("lead variant by |ols_cens|: %s ols_meas %.4f ols_cens %.4f tobit_z %.4f tobit_nmol %.2f",
                res$SNP[1], res$ols_meas[1], res$ols_cens[1], res$tobit_z[1], res$tobit_nmol[1]))

# ---- design control: score x statin on the biomarker ----
if (file.exists(opt$statin)) {
  st <- fread(opt$statin, colClasses = c(IID = "character"))
  sc <- numeric(nrow(g))
  for (i in seq_len(nrow(map))) {
    x <- g[[map$col[i]]]
    if (map$counted[i] != toupper(map$EA[i])) x <- 2 - x
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    sc <- sc + x * map$BETA[i]
  }
  e <- cbind(d, score_LPA = as.numeric(scale(sc)))
  e <- merge(e, st[, .(IID, statin)], by = "IID")[!is.na(statin)]
  pcs <- grep("^PC", covs, value = TRUE)
  rh <- paste(c("score_LPA * statin", "age", "sex", pcs), collapse = " + ")
  ols <- function(y) { cf <- summary(lm(as.formula(paste(y, "~", rh)), data = e))$coefficients; cf["score_LPA:statin", ] }
  tb <- function(L, R) { s <- summary(survreg(as.formula(sprintf("Surv(%s, %s, type = 'interval2') ~ %s", L, R, rh)), data = e, dist = "gaussian"))$table; s["score_LPA:statin", ] }
  rows <- list(
    data.table(model = "OLS raw nmol/L, measured only (published form)", t(ols("raw"))),
    data.table(model = "OLS LPA_MEAS", t(ols("LPA_MEAS"))),
    data.table(model = "OLS LPA_CENS", t(ols("LPA_CENS"))),
    data.table(model = "Tobit nmol/L", t(tb("nL", "nR"))),
    data.table(model = "Tobit log nmol/L", t(tb("lL", "lR"))),
    data.table(model = "Tobit z", t(tb("zL", "zR"))))
  ctl <- rbindlist(lapply(rows, function(r) { setnames(r, 2:ncol(r), c("estimate", "se", "stat", "p")[seq_len(ncol(r) - 1)]); r }), fill = TRUE)
  fwrite(ctl, paste0(opt$out, "_statin_control.csv"))
  sm <- c(sm, "score x statin -> biomarker:", capture.output(print(ctl)))
} else sm <- c(sm, sprintf("statin file %s absent; design control skipped", opt$statin))
writeLines(sm, paste0(opt$out, "_summary.txt")); writeLines(sm)
