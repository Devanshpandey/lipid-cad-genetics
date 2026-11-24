#!/usr/bin/env Rscript
# ============================================================================
# Does arterial susceptibility amplify the Lp(a) axis as it does the ApoB axis?
# Specification: PREREGISTRATION.md in this directory, committed before this
# is run on real data.
#
# Reuses, unchanged:
#   scripts/agent6_resilience/R/lib_{resilience,genetics,synergy}.R
#       cohort conventions, cox_fit, quad_amplification, excess_synergy
#   scripts/agent1_genetics/lpa_horizon/lpa_horizon.R
#       the Lp(a) exposure model and its gates, loaded into its own environment
#       because both code bases define check() and parse_args()
#
#   --step selftest   synthetic data only
#   --step run        gates (must pass), then the tests
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(optparse) })
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
REPO <- normalizePath(file.path(here, "..", "..", ".."))
A6 <- file.path(REPO, "scripts", "agent6_resilience")
for (f in c("lib_resilience.R", "lib_genetics.R", "lib_synergy.R")) source(file.path(A6, "R", f))
LH <- new.env(); sys.source(file.path(REPO, "scripts", "agent1_genetics", "lpa_horizon", "lpa_horizon.R"), envir = LH)

WINDOW <- 500L; WINDOWS <- c(250L, 500L, 1000L); HALF <- 0.5; SPLINE_DF <- 4L
REF <- list(n_all = 379654L, ev_all = 15978L, n_sf = 324626L, ev_sf = 11740L, s3_beta_LN = 0.0124399151356168,
            n_disc = 181162L, ev_disc = 6540L, s2_beta_MN = 0.0576577396138093, tol = 1e-6)
QUAD <- c("M", "L", "N", "I(M^2)", "I(L^2)", "I(N^2)", "M:L", "M:N", "L:N")
zs <- function(v, ref = v) (v - mean(ref, na.rm = TRUE)) / sd(ref, na.rm = TRUE)

# the full-quadratic model of quad_amplification(), returning every coefficient and the covariance
quad_full <- function(d, covars, strata_var = NULL, terms = QUAD, where = "quadratic") {
  f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s", rhs_str(terms, covars, strata_var)))
  fit <- cox_fit(f, d, where)
  list(b = coef(fit), V = vcov(fit), n = fit$n, events = fit$nevent)
}
# r = b(M:N) / b(M): fractional change of the exposure's own log hazard ratio per unit of N
ratio_of <- function(b, V, num, den) {
  r <- b[[num]] / b[[den]]
  se <- sqrt(V[num, num] / b[[den]]^2 + b[[num]]^2 * V[den, den] / b[[den]]^4 - 2 * b[[num]] * V[num, den] / b[[den]]^3)
  c(r = r, se = se)
}
summ <- function(q, label, m = "M") {
  b <- q$b; V <- q$V; mn <- paste0(m, ":N"); ml <- paste0(m, ":L")
  k <- c("L", "N", ml, mn); bb <- b[k]; VV <- V[k, k]
  aL <- function(x) x[[3]] / x[[1]]; aN <- function(x) x[[4]] / x[[2]]; D <- function(x) aN(x) - aL(x)
  pa <- ratio_of(b, V, mn, m); s <- function(nm) sqrt(V[nm, nm])
  data.table(analysis = label, n = as.integer(q$n), events = as.integer(q$events), beta_M = b[[m]], se_M = s(m), beta_L = b[["L"]], beta_N = b[["N"]],
    beta_MN = b[[mn]], se_MN = s(mn), p_MN = 2 * pnorm(-abs(b[[mn]] / s(mn))), beta_ML = b[[ml]], se_ML = s(ml),
    a_N = aN(bb), se_aN = delta_se(aN, bb, VV), D = D(bb), se_D = delta_se(D, bb, VV),
    r = pa[["r"]], se_r = pa[["se"]], r_lo = pa[["r"]] - 1.96 * pa[["se"]], r_hi = pa[["r"]] + 1.96 * pa[["se"]])
}
# PREREGISTRATION.md, "Primary test and decision rule"
verdict <- function(lp, apob_r, apob_se) {
  lo <- lp$beta_MN - 1.96 * lp$se_MN; hi <- lp$beta_MN + 1.96 * lp$se_MN; thr <- HALF * apob_r
  zc <- (lp$r - apob_r) / sqrt(lp$se_r^2 + apob_se^2); pc <- 2 * pnorm(-abs(zc))
  v <- if (lo > 0 && lp$r >= thr) "AMPLIFIED"
       else if (lo > 0) "AMPLIFIED_LESS_THAN_HALF"
       else if (hi < 0) "REVERSED"
       else if (lp$r_hi < thr && pc < 0.01) "NOT_AMPLIFIED_AXES_DIFFER"
       else "INCONCLUSIVE"
  data.table(verdict = v, lead_grade = v %in% c("AMPLIFIED", "REVERSED", "NOT_AMPLIFIED_AXES_DIFFER"), product_lo = lo, product_hi = hi,
             r = lp$r, r_lo = lp$r_lo, r_hi = lp$r_hi, half_of_apob_r = thr, apob_r = apob_r, apob_se_r = apob_se, contrast = lp$r - apob_r,
             contrast_z = zc, contrast_p = pc)
}
# additive scale: relative excess risk due to interaction for high M with high N
reri <- function(d, covars, label) {
  x <- copy(d); x[, g := factor(paste0(as.integer(hiM), as.integer(hiN)), levels = c("00", "10", "01", "11"))]
  fit <- cox_fit(as.formula(sprintf("Surv(time_analysis, event) ~ %s", rhs_str(c("g", "L"), covars, NULL))), x, "RERI")
  k <- c("g10", "g01", "g11"); b <- coef(fit)[k]; V <- vcov(fit)[k, k]
  g <- function(z) exp(z[[3]]) - exp(z[[1]]) - exp(z[[2]]) + 1; m <- function(z) z[[3]] - z[[1]] - z[[2]]
  data.table(analysis = label, n_both_high = x[g == "11", .N], events_both_high = x[g == "11", sum(event)], hr_M = exp(b[[1]]), hr_N = exp(b[[2]]),
             hr_both = exp(b[[3]]), reri = g(b), se_reri = delta_se(g, b, V), log_multiplicative = m(b), se_log_multiplicative = delta_se(m, b, V))
}
with_LN <- function(d, w, sdT = sd(d$TOTAL)) {
  x <- copy(d); l <- x[[paste0("L", w)]]; n <- x[[paste0("N", w)]]
  x[, `:=`(L = (l - mean(l)) / sdT, N = (n - mean(n)) / sdT)]; x
}

# every test; dall = eligible cohort with exposure columns, ddisc = discovery sample with M_apob and M_lpa
run_tests <- function(dall, ddisc, covars_all, covars_disc, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  base <- with_LN(dall, WINDOW)
  arm <- function(x, col, label, covars = covars_all, terms = QUAD, strata_var = NULL) {
    x <- droplevels(x[!is.na(x[[col]])]); x[, M := get(col)]
    summ(quad_full(x, covars, strata_var, terms, label), label)
  }
  spline_arm <- function(x, label) {
    x <- droplevels(x[!is.na(M_lpa)]); x[, M := M_lpa]
    sp <- spline_amplification(x, covars_all, SPLINE_DF, strata_var = NULL)
    data.table(analysis = label, n = nrow(x), events = as.integer(sum(x$event)), beta_MN = sp[["beta_MN"]], se_MN = sp[["se_MN"]],
               p_MN = 2 * pnorm(-abs(sp[["beta_MN"]] / sp[["se_MN"]])), beta_ML = sp[["beta_ML"]], se_ML = sp[["se_ML"]])
  }
  rows <- list(primary = arm(base, "M_lpa", "predicted Lp(a), all eligible (PRIMARY)"))
  rows$spline <- spline_arm(base, "predicted Lp(a), products beyond spline main effects")
  rows$statin_free <- arm(base[statin == 0], "M_lpa", "predicted Lp(a), statin-free")
  nosex <- setdiff(covars_all, "sex")
  for (sx in levels(base$sex)) rows[[paste0("sex", sx)]] <- arm(base[sex == sx], "M_lpa", paste("predicted Lp(a), sex", sx), covars = nosex)
  for (w in setdiff(WINDOWS, WINDOW)) rows[[paste0("w", w)]] <- arm(with_LN(dall, w), "M_lpa", sprintf("predicted Lp(a), %d kb window", w))
  rows$measured <- arm(base, "M_lpa_measured", "measured Lp(a), assay limits imposed")
  rows$ldl <- arm(base, "M_ldl", "LDL-C score (167 variants)")
  arms <- rbindlist(rows, fill = TRUE)

  # the ApoB reference, re-fitted in the discovery sample to obtain its covariance
  disc <- with_LN(ddisc, WINDOW)
  ap <- copy(disc)[!is.na(M_apob)]; ap[, M := M_apob]
  apob <- summ(quad_full(ap, covars_disc, "fold", QUAD, "ApoB reference"), "measured ApoB, discovery sample (reference)")
  ver <- verdict(arms[1], apob$r, apob$se_r)
  same_side <- function(a) (a$beta_MN - 1.96 * a$se_MN > 0 && ver$product_lo > 0) || (a$beta_MN + 1.96 * a$se_MN < 0 && ver$product_hi < 0)
  ver[, robust := if (verdict %in% c("AMPLIFIED", "REVERSED")) same_side(rows$statin_free) && same_side(rows$spline) else NA]

  # one model with both exposures in the discovery sample
  jd <- copy(disc)[!is.na(M_apob) & !is.na(M_lpa)]; setnames(jd, c("M_apob", "M_lpa"), c("M1", "M2"))
  jt <- c("M1", "M2", "L", "N", "I(M1^2)", "I(M2^2)", "I(L^2)", "I(N^2)", "M1:L", "M1:N", "M2:L", "M2:N", "L:N", "M1:M2")
  jq <- quad_full(jd, covars_disc, "fold", jt, "joint model")
  k <- c("M1", "M2", "M1:N", "M2:N"); jb <- jq$b[k]; jV <- jq$V[k, k]
  dr <- function(z) z[[3]] / z[[1]] - z[[4]] / z[[2]]
  r1 <- ratio_of(jq$b, jq$V, "M1:N", "M1"); r2 <- ratio_of(jq$b, jq$V, "M2:N", "M2"); dse <- delta_se(dr, jb, jV)
  joint <- data.table(n = as.integer(jq$n), events = as.integer(jq$events), beta_apob_N = jq$b[["M1:N"]], se_apob_N = sqrt(jq$V["M1:N", "M1:N"]),
    beta_lpa_N = jq$b[["M2:N"]], se_lpa_N = sqrt(jq$V["M2:N", "M2:N"]), r_apob = r1[["r"]], se_r_apob = r1[["se"]], r_lpa = r2[["r"]], se_r_lpa = r2[["se"]],
    apob_minus_lpa = dr(jb), se_difference = dse, p_difference = 2 * pnorm(-abs(dr(jb) / dse)))

  # additive scale
  rr <- copy(base)[!is.na(M_lpa)]
  rr[, `:=`(hiM = M_lpa >= quantile(M_lpa, 0.90), hiN = N >= quantile(N, 0.80))]
  add <- reri(rr, covars_all, "top decile of predicted Lp(a) with top quintile of the non-lipid score")

  LH$wcsv(rbind(arms, apob, fill = TRUE), file.path(outdir, "product_terms.csv")); LH$wcsv(ver, file.path(outdir, "verdict.csv"))
  LH$wcsv(joint, file.path(outdir, "joint_model.csv")); LH$wcsv(add, file.path(outdir, "additive_scale.csv"))
  list(arms = arms, apob = apob, verdict = ver, joint = joint, add = add)
}

# ---------------------------------------------------------------------------
# Real data, assembled as stages 10 and 11 of the earlier pipeline assemble it
# ---------------------------------------------------------------------------
load_eligible <- function(run_dir, components, cfg, sy) {
  keep <- c("eid", "excl_no_baseline_date", "excl_genetic_qc", "excl_prevalent", "excl_short_followup", "time_years", "event",
            "age", "sex", paste0("PC", 1:10), "genotype_batch", "statin")
  coh <- fread(file.path(run_dir, "01", "cohort.tsv"), select = keep, na.strings = "NA")
  el <- coh[excl_no_baseline_date == 0 & excl_genetic_qc == 0 & excl_prevalent == 0 & excl_short_followup == 0]
  el <- apply_landmark(el, cfg$outcome$landmark_years)
  comp <- fread(components); comp[, eid := as.integer(IID)]
  d <- merge(el, comp[, !c("FID", "IID")], by = "eid")
  d <- d[stats::complete.cases(d[, c("age", "sex", paste0("PC", 1:10), "genotype_batch"), with = FALSE])]
  d[, sex := factor(sex)]
  if (isTRUE(sy$s3$batch_as_factor)) d[, genotype_batch := factor(genotype_batch)]
  d
}
load_discovery <- function(run_dir, components, sy) {
  dt <- fread(file.path(run_dir, "04", "mrs_oof.tsv"), na.strings = "NA"); dt[, sex := factor(sex)]
  f <- paste0("nmr_", sy$primary_field)
  res <- fread(file.path(run_dir, "03", "nmr_technical_resid.tsv"), select = c("eid", f), na.strings = "NA")
  comp <- fread(components); comp[, eid := as.integer(IID)]
  dt <- merge(dt, comp[, !c("FID", "IID")], by = "eid", sort = FALSE)
  dt[, M_apob := res[[f]][match(eid, res$eid)]]
  dt
}
attach_exposures <- function(d, P, ref = NULL) {
  ex <- P[, .(eid = as.integer(IID), pred, lpa, status, score_ldl_raw)]
  x <- merge(d, ex, by = "eid", all.x = TRUE, sort = FALSE)
  meas <- fifelse(x$status == "below_range", LH$LO, fifelse(x$status == "above_range", LH$HI, fifelse(x$status == "measured", x$lpa, NA_real_)))
  if (is.null(ref)) ref <- list(pred = x$pred, meas = meas, ldl = x$score_ldl_raw)
  x[, `:=`(M_lpa = zs(pred, ref$pred), M_lpa_measured = zs(meas, ref$meas), M_ldl = zs(score_ldl_raw, ref$ldl))]
  setattr(x, "ref", ref)
  x
}
real_gates <- function(dall, ddisc, covars_all, covars_disc) {
  s3 <- excess_synergy(with_LN(dall, WINDOW), covars_all)
  ap <- with_LN(ddisc, WINDOW)[!is.na(M_apob)]; ap[, M := M_apob]
  s2 <- quad_amplification(ap, covars_disc)
  mine <- summ(quad_full(ap, covars_disc, "fold"), "check")
  g <- data.table(
    gate = c("eligible cohort", "its events", "statin-free part", "its events", "earlier lipid-region by non-lipid product reproduced",
             "discovery sample", "its events", "earlier ApoB product reproduced", "this script's model equals the earlier function",
             "eligible cohort with predicted Lp(a), %", "|correlation| of predicted Lp(a) with N"),
    value = c(nrow(dall), sum(dall$event), dall[statin == 0, .N], dall[statin == 0, sum(event)], s3$beta_LN, nrow(ap), sum(ap$event), s2$beta_MN,
              abs(mine$beta_MN - s2$beta_MN) + abs(mine$D - s2$D), 100 * mean(!is.na(dall$M_lpa)), abs(cor(dall$M_lpa, with_LN(dall, WINDOW)$N, use = "complete.obs"))),
    target = c(REF$n_all, REF$ev_all, REF$n_sf, REF$ev_sf, REF$s3_beta_LN, REF$n_disc, REF$ev_disc, REF$s2_beta_MN, 0, 99.9, 0.02),
    rule = c(rep("equal", 4), "within 1e-6", "equal", "equal", "within 1e-6", "below 1e-9", "at least", "below"))
  g[, pass := fifelse(rule == "equal", value == target, fifelse(rule == "within 1e-6", abs(value - target) < REF$tol,
              fifelse(rule == "below 1e-9", value < 1e-9, fifelse(rule == "at least", value >= target, value < target))))]
  g
}

# ---------------------------------------------------------------------------
# Self-test (synthetic data only)
# ---------------------------------------------------------------------------
sim_cohort <- function(n, bM, g, seed, bM2 = 0, g2 = 0) {
  set.seed(seed)
  lpa <- pmax(0, rgamma(n, shape = 1.2, scale = 14) + 110 * rbinom(n, 2, 0.08))
  d <- data.table(eid = seq_len(n), M_lpa = zs(lpa), M_apob = rnorm(n), M_ldl = rnorm(n), L500 = rnorm(n, 0, 0.33), N500 = rnorm(n, 0, 0.94),
                  age = runif(n, 40, 70), sex = factor(rbinom(n, 1, 0.5)), statin = rbinom(n, 1, 0.12), genotype_batch = factor(sample(1:3, n, TRUE)),
                  fold = sample(1:5, n, TRUE))
  for (w in c(250, 1000)) { set(d, j = paste0("L", w), value = d$L500 + rnorm(n, 0, 0.05)); set(d, j = paste0("N", w), value = d$N500 + rnorm(n, 0, 0.05)) }
  for (k in 1:10) set(d, j = paste0("PC", k), value = rnorm(n))
  d[, TOTAL := L500 + N500]; d[, M_lpa_measured := M_lpa + rnorm(n, 0, 0.3)]
  sdT <- sd(d$TOTAL); Ls <- d$L500 / sdT; Ns <- d$N500 / sdT
  lin <- bM * d$M_lpa + 0.55 * Ls + 0.55 * Ns + g * d$M_lpa * Ns + bM2 * d$M_apob + g2 * d$M_apob * Ns + 0.04 * (d$age - 55)
  t <- rexp(n, 0.004 * exp(lin)); d[, `:=`(event = as.integer(t <= 11), time_analysis = pmin(t, 11))]
  d
}
selftest <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { if (!isTRUE(cond)) stop("SELFTEST FAILED: ", msg, call. = FALSE); n_ok <<- n_ok + 1L }
  cv <- c("age", "sex", paste0("PC", 1:10), "genotype_batch")

  # a planted product is recovered, and this script's model is the earlier function's model
  d <- sim_cohort(250000, 0.15, 0.045, 41); x <- with_LN(d, WINDOW); x[, M := M_lpa]
  q <- quad_full(x, cv); s <- summ(q, "sim"); old <- quad_amplification(x, cv, strata_var = NULL)
  message(sprintf("[suscept sim] planted product 0.045: %.4f (SE %.4f); r %.3f (SE %.3f), planted 0.30; events %d", s$beta_MN, s$se_MN, s$r, s$se_r, s$events))
  expect(abs(s$beta_MN - 0.045) < 3.5 * s$se_MN && abs(s$r - 0.30) < 3.5 * s$se_r, "a planted product and its ratio are recovered")
  expect(abs(s$beta_MN - old$beta_MN) < 1e-10 && abs(s$a_N - old$a_N) < 1e-10 && abs(s$D - old$D) < 1e-10 && abs(s$se_D - old$se_D) < 1e-8,
         "quad_full and summ reproduce quad_amplification exactly")
  d0 <- sim_cohort(250000, 0.15, 0, 42); x0 <- with_LN(d0, WINDOW); x0[, M := M_lpa]; s0 <- summ(quad_full(x0, cv), "null")
  expect(abs(s0$beta_MN) < 3.5 * s0$se_MN, "no product is found where none is planted")

  # delta-method SE of the ratio agrees with simulation from the coefficient distribution
  set.seed(43); k <- c("M", "M:N"); Z <- matrix(rnorm(2e5 * 2), ncol = 2) %*% chol(q$V[k, k]); B <- sweep(Z, 2, q$b[k], "+")
  expect(abs(sd(B[, 2] / B[, 1]) / s$se_r - 1) < 0.10, "delta-method SE of the ratio matches simulation")

  # every verdict branch
  row <- function(b, se, bM = 0.15, seM = 0.008) { r <- b / bM; sr <- sqrt(se^2 / bM^2 + b^2 * seM^2 / bM^4)
    data.table(beta_MN = b, se_MN = se, r = r, se_r = sr, r_lo = r - 1.96 * sr, r_hi = r + 1.96 * sr) }
  expect(verdict(row(0.045, 0.01), 0.29, 0.075)$verdict == "AMPLIFIED", "verdict: amplified")
  expect(verdict(row(0.018, 0.006), 0.29, 0.075)$verdict == "AMPLIFIED_LESS_THAN_HALF", "verdict: amplified, less than half")
  expect(verdict(row(-0.03, 0.01), 0.29, 0.075)$verdict == "REVERSED", "verdict: reversed")
  expect(verdict(row(0.002, 0.008), 0.29, 0.075)$verdict == "NOT_AMPLIFIED_AXES_DIFFER", "verdict: not amplified, axes differ")
  expect(verdict(row(0.012, 0.010), 0.29, 0.075)$verdict == "INCONCLUSIVE", "verdict: inconclusive")
  v <- verdict(row(0.002, 0.008), 0.29, 0.075)
  expect(v$lead_grade && abs(v$half_of_apob_r - 0.145) < 1e-12 && v$contrast_p < 0.01 && !verdict(row(0.012, 0.010), 0.29, 0.075)$lead_grade, "lead-grade flag follows the rule")

  # additive scale: risks that multiply exactly still show a positive excess on the additive scale
  set.seed(44); n <- 400000
  a <- data.table(hiM = rbinom(n, 1, 0.1) == 1, hiN = rbinom(n, 1, 0.2) == 1, L = rnorm(n), age = runif(n, 40, 70))
  t <- rexp(n, 0.01 * exp(log(1.5) * a$hiM + log(2) * a$hiN)); a[, `:=`(event = as.integer(t <= 10), time_analysis = pmin(t, 10))]
  ra <- reri(a, "age", "sim")
  message(sprintf("[suscept sim] multiplicative truth: RERI %.2f (SE %.2f), expected 0.50; log multiplicative %.3f (SE %.3f)", ra$reri, ra$se_reri, ra$log_multiplicative, ra$se_log_multiplicative))
  expect(abs(ra$reri - 0.5) < 3.5 * ra$se_reri && abs(ra$log_multiplicative) < 3.5 * ra$se_log_multiplicative && abs(ra$hr_both - 3) < 0.4, "additive and multiplicative scales are both reported correctly")

  # the whole analysis on synthetic cohorts: ApoB amplified, Lp(a) not
  dall <- sim_cohort(300000, 0.15, 0, 45); ddisc <- sim_cohort(300000, 0.15, 0, 46, bM2 = 0.20, g2 = 0.06)
  out <- run_tests(dall, ddisc, cv, setdiff(cv, "genotype_batch"), file.path(workdir, "suscept_smoke"))
  files <- c("product_terms.csv", "verdict.csv", "joint_model.csv", "additive_scale.csv")
  expect(all(file.exists(file.path(workdir, "suscept_smoke", files))), "every output file is written")
  expect(nrow(out$arms) == 9 && all(is.finite(out$arms$beta_MN)) && out$arms$analysis[1] == "predicted Lp(a), all eligible (PRIMARY)", "every arm is fitted and the primary is first")
  message(sprintf("[suscept sim] ApoB reference r %.3f (SE %.3f), planted 0.30; Lp(a) r %.3f (SE %.3f), planted 0; verdict %s; joint difference %.3f (SE %.3f)",
                  out$apob$r, out$apob$se_r, out$arms$r[1], out$arms$se_r[1], out$verdict$verdict, out$joint$apob_minus_lpa, out$joint$se_difference))
  expect(abs(out$apob$r - 0.30) < 3.5 * out$apob$se_r && abs(out$arms$r[1]) < 3.5 * out$arms$se_r[1], "the reference and the primary ratio recover their planted values")
  expect(abs(out$joint$apob_minus_lpa - 0.30) < 3.5 * out$joint$se_difference && abs(out$joint$r_apob - out$apob$r) < 0.05, "the joint model recovers the planted difference between the axes")
  again <- verdict(out$arms[1], out$apob$r, out$apob$se_r)
  expect(out$verdict$verdict == again$verdict && out$verdict$verdict %in% c("AMPLIFIED", "AMPLIFIED_LESS_THAN_HALF", "REVERSED", "NOT_AMPLIFIED_AXES_DIFFER", "INCONCLUSIVE") &&
           (is.na(out$verdict$robust) || is.logical(out$verdict$robust)), "the verdict is the rule applied to the primary arm and the reference")
  cat(sprintf("\nSUSCEPTIBILITY SELFTEST PASSED: %d checks\n", n_ok))
}

opt <- optparse::parse_args(OptionParser(option_list = list(
  make_option("--step", type = "character", default = "selftest"),
  make_option("--run_dir", type = "character"), make_option("--components", type = "character"),
  make_option("--raw_dir", type = "character"), make_option("--covar", type = "character"),
  make_option("--dose_lpa", type = "character"), make_option("--w_lpa", type = "character"),
  make_option("--dose_ldl", type = "character"), make_option("--w_ldl", type = "character"),
  make_option("--prs", type = "character"), make_option("--outdir", type = "character"),
  make_option("--workdir", type = "character", default = tempdir()),
  make_option("--censor", type = "character", default = "2022-10-31"),
  make_option("--threads", type = "integer", default = 8L))))
setDTthreads(opt$threads)
if (opt$step == "selftest") { selftest(opt$workdir); quit(save = "no", status = 0L) }
if (opt$step != "run") stop("unknown --step")
cfg <- load_config(file.path(A6, "config.yaml")); sy <- yaml::read_yaml(file.path(A6, "config_synergy.yaml"))
CENS <- as.integer(as.IDate(opt$censor))
R <- LH$load_real(opt, CENS)
if (!LH$gates(R, CENS, file.path(opt$outdir, "gates_recheck"))) { cat("\nGATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
covars_disc <- sy$covariates; covars_all <- c(sy$covariates, "genotype_batch")
dall <- attach_exposures(load_eligible(opt$run_dir, opt$components, cfg, sy), R$P)
ddisc <- attach_exposures(load_discovery(opt$run_dir, opt$components, sy), R$P, ref = attr(dall, "ref"))
g <- real_gates(dall, ddisc, covars_all, covars_disc)
options(width = 250); print(g, digits = 10)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE); fwrite(g, file.path(opt$outdir, "gates_susceptibility.csv"))
if (!isTRUE(all(g$pass))) { cat("\nSUSCEPTIBILITY GATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
cat("\nGATES PASSED\n")
out <- run_tests(dall, ddisc, covars_all, covars_disc, opt$outdir)
cat("\n=== product terms (log hazard per SD of exposure x unit of score) ===\n")
print(rbind(out$arms, out$apob, fill = TRUE)[, .(analysis, n, events, beta_M, se_M, beta_MN, se_MN, p_MN, beta_ML, se_ML, a_N, se_aN, D, se_D, r, r_lo, r_hi)], digits = 3)
cat("\n=== VERDICT ===\n"); print(t(out$verdict))
cat("\n=== one model with both exposures, discovery sample ===\n"); print(out$joint, digits = 3)
cat("\n=== additive scale ===\n"); print(out$add, digits = 3)
report_warnings()
cat("\nSUSCEPTIBILITY_DONE\n")
