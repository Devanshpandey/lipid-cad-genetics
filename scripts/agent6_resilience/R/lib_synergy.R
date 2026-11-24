# ============================================================================
# lib_synergy.R - functions for the synergy follow-up (preregistered S1-S3):
# phase split, lipid regions and PRS partition, spline interaction, full-
# quadratic amplification contrast, excess genetic synergy, H3a-style scans.
# Unit-tested in tests/testthat/test_synergy.R. Source lib_resilience.R and
# lib_genetics.R first. Separate file so stage 00-09 hashes are untouched.
# ============================================================================
suppressPackageStartupMessages(library(splines))

# Delta-method SE of a scalar function g(b) (central differences).
delta_se <- function(g, b, V) {
  grad <- vapply(seq_along(b), function(i) {
    h <- 1e-6 * max(1, abs(b[i])); bp <- b; bm <- b
    bp[i] <- b[i] + h; bm[i] <- b[i] - h
    (g(bp) - g(bm)) / (2 * h)
  }, numeric(1))
  sqrt(as.numeric(t(grad) %*% V %*% grad))
}

phase_from_dates <- function(dates, cut) {
  d <- as.IDate(dates)
  fifelse(is.na(d), NA_integer_, fifelse(d < as.IDate(cut), 1L, 2L))
}

# ---- lipid regions and PRS partition ------------------------------------------------
merge_regions <- function(r) {
  r <- copy(r)[order(CHR, start)]
  r[, grp := cumsum(c(1L, as.integer(start[-1] > cummax(end)[-.N]))), by = CHR]
  r[, .(start = min(start), end = max(end)), by = .(CHR, grp)][, grp := NULL][]
}

regions_from_leads <- function(leads, window_bp) {
  merge_regions(leads[, .(CHR = as.integer(CHROM), start = pmax(0, GENPOS - window_bp), end = GENPOS + window_bp)])
}

# TRUE where (chr, pos) lies inside a region (inclusive ends). Regions must be merged.
in_regions <- function(chr, pos, regions) {
  out <- logical(length(chr))
  for (ch in unique(chr)) {
    idx <- which(chr == ch); rr <- regions[CHR == ch][order(start)]
    if (!nrow(rr)) next
    k <- findInterval(pos[idx], rr$start)
    ok <- k > 0
    out[idx[ok]] <- pos[idx[ok]] <= rr$end[k[ok]]
  }
  out
}

# w: CHR SNP BP A1 A2 BETA (PRS-CS output). Returns SNP, A1 and one L<kb>/N<kb> pair per window.
partition_weights <- function(w, regions_by_window) {
  out <- w[, .(SNP, A1)]
  for (nm in names(regions_by_window)) {
    inr <- in_regions(w$CHR, w$BP, regions_by_window[[nm]])
    set(out, j = paste0("L", nm), value = ifelse(inr, w$BETA, 0))
    set(out, j = paste0("N", nm), value = ifelse(inr, 0, w$BETA))
  }
  out
}

# ---- models ---------------------------------------------------------------------
cox_fit <- function(formula, data, where) capture_warnings(coxph(formula, data = data, ties = "efron"), where)
rhs_str <- function(terms, covars, strata_var) paste(c(terms, covars, if (!is.null(strata_var)) sprintf("strata(%s)", strata_var)), collapse = " + ")

# S1a: product M*P beyond 4-df natural-spline main effects of M and P.
fit_spline_interaction <- function(d, covars, df, strata_var = "fold") {
  d <- copy(d); d[, MP := M * P]
  f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s",
                          rhs_str(c(sprintf("ns(M, df = %d)", df), sprintf("ns(P, df = %d)", df), "MP"), covars, strata_var)))
  s <- summary(cox_fit(f, d, "spline interaction"))$coefficients
  c(beta = s["MP", "coef"], se = s["MP", "se(coef)"], p = s["MP", "Pr(>|z|)"])
}

# S2: full quadratic in (M, L, N). a_L = b_ML / b_L, a_N = b_MN / b_N, D = a_N - a_L.
quad_amplification <- function(d, covars, strata_var = "fold") {
  f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s",
                          rhs_str(c("M", "L", "N", "I(M^2)", "I(L^2)", "I(N^2)", "M:L", "M:N", "L:N"), covars, strata_var)))
  fit <- cox_fit(f, d, "quadratic amplification")
  b <- coef(fit); V <- vcov(fit); k <- c("L", "N", "M:L", "M:N")
  check(all(k %in% names(b)), "quadratic model carries the L, N, M:L and M:N terms")
  bb <- b[k]; VV <- V[k, k]
  aL <- function(x) x[[3]] / x[[1]]; aN <- function(x) x[[4]] / x[[2]]; D <- function(x) aN(x) - aL(x)
  est <- D(bb); se <- delta_se(D, bb, VV)
  list(n = nrow(d), events = sum(d$event), beta_M = unname(b["M"]), beta_L = unname(b["L"]), beta_N = unname(b["N"]),
       beta_ML = unname(b["M:L"]), se_ML = sqrt(V["M:L", "M:L"]), beta_MN = unname(b["M:N"]), se_MN = sqrt(V["M:N", "M:N"]),
       a_L = aL(bb), se_aL = delta_se(aL, bb, VV), a_N = aN(bb), se_aN = delta_se(aN, bb, VV),
       D = est, se_D = se, p_D = 2 * pnorm(-abs(est / se)))
}

# S2 sensitivity: products beyond 4-df spline main effects of M, L and N.
spline_amplification <- function(d, covars, df, strata_var = "fold") {
  f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s",
                          rhs_str(c(sprintf("ns(%s, df = %d)", c("M", "L", "N"), df), "M:L", "M:N"), covars, strata_var)))
  s <- summary(cox_fit(f, d, "spline amplification"))$coefficients
  c(beta_ML = s["M:L", "coef"], se_ML = s["M:L", "se(coef)"], beta_MN = s["M:N", "coef"], se_MN = s["M:N", "se(coef)"])
}

# S3: excess synergy E = b_LN - [b_LL (b_N / b_L) + b_NN (b_L / b_N)]; zero under convexity in L + N alone.
excess_synergy <- function(d, covars, strata_var = NULL) {
  f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s",
                          rhs_str(c("L", "N", "I(L^2)", "I(N^2)", "L:N"), covars, strata_var)))
  fit <- cox_fit(f, d, "excess synergy")
  b <- coef(fit); V <- vcov(fit); k <- c("L", "N", "I(L^2)", "I(N^2)", "L:N")
  bb <- b[k]; VV <- V[k, k]
  E <- function(x) x[[5]] - (x[[3]] * x[[2]] / x[[1]] + x[[4]] * x[[1]] / x[[2]])
  e <- E(bb); se <- delta_se(E, bb, VV)
  list(n = nrow(d), events = sum(d$event), beta_L = unname(b["L"]), beta_N = unname(b["N"]), beta_LL = unname(b["I(L^2)"]),
       beta_NN = unname(b["I(N^2)"]), beta_LN = unname(b["L:N"]), se_LN = sqrt(V["L:N", "L:N"]),
       E = e, se_E = se, p_E = 2 * pnorm(-abs(e / se)))
}

# H3a-style scan (mirrors stage 05): PRS x PRS-orthogonal metabolite residual, per metabolite.
h3a_scan <- function(d, Rm, covars, strata_var = "fold") {
  rbindlist(lapply(colnames(Rm), function(cn) {
    x <- Rm[, cn]; ok <- !is.na(x); dd <- d[ok]; dd[, met := x[ok]]
    dd[, met_r := resid(lm(met ~ prs_z))]; dd[, met_r := met_r / sd(met_r)]
    f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s", rhs_str("prs_z * met_r", covars, strata_var)))
    s <- summary(cox_fit(f, dd, "h3a scan"))$coefficients
    data.table(field = cn, n = nrow(dd), beta_int = s["prs_z:met_r", "coef"], se_int = s["prs_z:met_r", "se(coef)"],
               z_int = s["prs_z:met_r", "z"], p_int = s["prs_z:met_r", "Pr(>|z|)"])
  }))
}

# Statin-free re-test of lead variants (DEVIATION 2026-09-29): linear regression of y on each
# genotype column with covariates X (design matrix incl. intercept), by Frisch-Waugh-Lovell:
# y and every genotype column are residualised on X once. Missing dosages are mean-imputed.
retest_leads <- function(y, G, X) {
  G <- as.matrix(G); storage.mode(G) <- "double"
  for (j in seq_len(ncol(G))) if (anyNA(G[, j])) { v <- G[, j]; v[is.na(v)] <- mean(v, na.rm = TRUE); G[, j] <- v }
  q <- qr(X); ry <- qr.resid(q, y); rG <- qr.resid(q, G)
  sxx <- colSums(rG^2); b <- colSums(rG * ry) / sxx
  df <- length(y) - q$rank - 1
  rss <- sum(ry^2) - b^2 * sxx
  se <- sqrt(rss / df / sxx)
  data.table(beta = b, se = se, p = 2 * pt(-abs(b / se), df))
}
