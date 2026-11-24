#!/usr/bin/env Rscript
# ============================================================================
# FURTHER POST HOC analyses for the Lp(a)HORIZON genetic emulation (P3-P5).
#
# Written after the pre-registered run (job 3600628) and the first post hoc
# round (job 3600785) had been read. Declared in PREREGISTERED_READOUT.md,
# section "Further post hoc analyses", before being run. Nothing here changes
# the verdict of record.
#
# lpa_horizon.R is sourced unchanged. lpa_horizon_posthoc.R is not sourced
# (its main block is unguarded); its function definitions and two constants
# are evaluated from the parsed file, so the helpers are the ones already run.
#
#   P3  inverse stroke association: selection budget, adjustment, first events
#   P4  calibration against LDL-C under stronger statin corrections
#   P5  linearity of the Lp(a) association into the trial's range
#
#   --step selftest   synthetic data only
#   --step run        gates (must pass), then P3-P5
# ============================================================================
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
source(file.path(here, "lpa_horizon.R"))
for (e in parse(file.path(here, "lpa_horizon_posthoc.R"))) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      (as.character(e[[2]]) %in% c("CTT", "TRIAL_SE") ||
       (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))))) eval(e, globalenv())
}
check(all(vapply(c("add_groups", "fit_strat", "fieller", "ldl_slope", "calibrate", "CTT", "TRIAL_SE"), exists, logical(1))),
      "the first post hoc script did not provide its functions")

STROKE_RISK <- c(RISK, "af_pre", "sbp")
CORONARY <- c("t_mi", "t_rev")
BANDS <- c(-Inf, 50, 100, 150, Inf); BAND_LABELS <- c("<50", "50-99", "100-149", ">=150")
KNOTS <- c(150, 100)

# Cox model on a supplied time/event table, with optional strata on the qualifying
# condition, weights, and extra columns that enter as further terms
fit_tv <- function(D, tv, xs = EXPOSURES, adj = ADJ, strat = FALSE, w = NULL, extra = NULL) {
  cols <- unique(c(xs, adj, if (strat) "qgroup"))
  X <- cbind(D[, cols, with = FALSE], tv)
  if (!is.null(extra)) X <- cbind(X, extra)
  X[, wt := if (is.null(w)) 1 else w]
  X <- X[time > 0 & wt > 0]; X <- X[complete.cases(X)]
  if (sum(X$ev) < 20L) return(NULL)
  keep <- c(xs, names(extra))
  f <- as.formula(paste("Surv(time, ev) ~", paste(c(keep, adj), collapse = " + "), if (strat) "+ strata(qgroup)" else ""))
  m <- if (is.null(w)) coxph(f, data = X) else coxph(f, data = X, weights = wt, robust = TRUE)
  b <- coef(m)[keep]; if (anyNA(b)) return(NULL)
  list(b = b, V = vcov(m)[keep, keep, drop = FALSE], n = m$n, events = sum(X$ev), all = coef(m))
}
# time to `endpoint`, with follow-up also ending at the first of the `censor_at` events
mk_first <- function(D, endpoint, censor_at) {
  cols <- ENDPOINTS[[endpoint]]
  t <- do.call(pmin, c(unname(as.list(D[, cols, with = FALSE])), na.rm = TRUE))
  tc <- do.call(pmin, c(unname(as.list(D[, censor_at, with = FALSE])), na.rm = TRUE))
  lim <- pmin(D$t_end, tc, na.rm = TRUE)
  ev <- as.integer(!is.na(t) & t > D$t0 & t <= lim)
  data.table(time = (ifelse(ev == 1L, t, lim) - D$t0) / 365.25, ev = ev)
}
row_of <- function(f, x, ...) {
  if (is.null(f)) return(NULL)
  b <- f$b[[x]]; se <- sqrt(f$V[x, x])
  data.table(..., exposure = x, n = as.integer(f$n), events = as.integer(f$events), b = b, se = se,
             hr = exp(b), hr_lo = exp(b - 1.96 * se), hr_hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}

# ---------------------------------------------------------------------------
# P3
# ---------------------------------------------------------------------------
# the pre-registered budget, for any endpoint and any set of unmodelled factors; `ref` is
# the cohort in which the exposure is taken to be independent of the factors
budget <- function(A, ref, endpoint, factors) {
  fa <- fit(A, endpoint, adj = c(ADJ, factors))
  if (is.null(fa)) return(NULL)
  rbindlist(lapply(c("x_ldl", "x_prs", "age", "sex", factors), function(k) {
    sA <- coef(lm(A[[k]] ~ A$x_lpa))[[2]]; sC <- coef(lm(ref[[k]] ~ ref$x_lpa))[[2]]
    data.table(endpoint = endpoint, factor = k, in_primary_model = !(k %in% factors), slope_established = sA, slope_primary = sC,
               induced = sA - sC, log_hr_of_factor = fa$all[[k]], implied_bias_log_hr = (sA - sC) * fa$all[[k]])
  }))
}
stroke_checks <- function(A, Cc, w) {
  bud <- rbindlist(lapply(names(ENDPOINTS), function(e) budget(A, Cc, e, STROKE_RISK)))
  Ak <- A[complete.cases(A[, STROKE_RISK, with = FALSE])]
  tvs <- list("all later strokes" = function(D) mk(D, "stroke"),
              "stroke as a first event (censored at infarction or revascularization)" = function(D) mk_first(D, "stroke", CORONARY))
  mods <- list()
  for (fu in names(tvs)) for (st in c(FALSE, TRUE)) {
    lab <- if (st) "stratified on qualifying condition" else "pooled"
    mods[[length(mods) + 1]] <- row_of(fit_tv(A, tvs[[fu]](A), strat = st), "x_lpa", cohort = "A", follow_up = fu, model = lab, adjustment = "primary model")
    mods[[length(mods) + 1]] <- row_of(fit_tv(Ak, tvs[[fu]](Ak), strat = st), "x_lpa", cohort = "A", follow_up = fu, model = lab, adjustment = "primary model, participants with every factor")
    mods[[length(mods) + 1]] <- row_of(fit_tv(Ak, tvs[[fu]](Ak), adj = c(ADJ, STROKE_RISK), strat = st), "x_lpa", cohort = "A", follow_up = fu, model = lab, adjustment = "plus stroke risk factors")
  }
  for (fu in names(tvs))
    mods[[length(mods) + 1]] <- row_of(fit_tv(Cc, tvs[[fu]](Cc), w = w), "x_lpa", cohort = "C standardized", follow_up = fu, model = "pooled", adjustment = "primary model")
  mods <- rbindlist(mods)
  mods[, `:=`(hr_per_144_lower = exp(-b), lower_lo = exp(-b - 1.96 * se), lower_hi = exp(-b + 1.96 * se))]
  factors <- data.table(cohort = c("A", "C"), n = c(nrow(A), nrow(Cc)),
    af_pct = 100 * c(mean(A$af_pre), mean(Cc$af_pre)), sbp_mean = c(mean(A$sbp, na.rm = TRUE), mean(Cc$sbp, na.rm = TRUE)),
    sbp_missing_pct = 100 * c(mean(is.na(A$sbp)), mean(is.na(Cc$sbp))),
    every_factor_pct = 100 * c(mean(complete.cases(A[, STROKE_RISK, with = FALSE])), mean(complete.cases(Cc[, STROKE_RISK, with = FALSE]))))
  list(budget = bud, models = mods, factors = factors)
}

# ---------------------------------------------------------------------------
# P4
# ---------------------------------------------------------------------------
# as ldl_slope(), with any divisor for statin users and an optional subset
ldl_slope_d <- function(D, divisor = 1, w = NULL, keep = NULL) {
  y <- ifelse(D$statin, D$ldl / divisor, D$ldl)
  X <- cbind(D[, c("x_ldl", "x_lpa", ADJ), with = FALSE], y = y, wt = if (is.null(w)) 1 else w)
  if (!is.null(keep)) X <- X[keep]
  X <- X[complete.cases(X) & wt > 0]
  m <- lm(as.formula(paste("y ~", paste(c("x_ldl", "x_lpa", ADJ), collapse = " + "))), data = X, weights = wt)
  coef(m)[["x_ldl"]]
}
statin_calibration <- function(A, Cc, w) {
  nonC <- ldl_slope_d(Cc, w = w, keep = !Cc$statin); nonA <- ldl_slope_d(A, keep = !A$statin)
  sc <- function(D, ww) c("measured LDL-C" = ldl_slope_d(D, 1, ww), "LDL-C / 0.7 in statin users" = ldl_slope_d(D, 0.7, ww),
                          "LDL-C / 0.6 in statin users" = ldl_slope_d(D, 0.6, ww), "LDL-C / 0.5 in statin users" = ldl_slope_d(D, 0.5, ww),
                          "non-users, primary prevention" = nonC)
  slopes <- list(A = c(sc(A, NULL), "non-users, established disease" = nonA), C = sc(Cc, w))
  check(abs(slopes$A[["measured LDL-C"]] - ldl_slope(A)) < 1e-10 && abs(slopes$A[["LDL-C / 0.7 in statin users"]] - ldl_slope(A, TRUE)) < 1e-10,
        "the LDL-C slopes of the first post hoc round are not reproduced")
  cal <- list()
  for (h in c("full", "5y")) {
    hz <- if (h == "5y") H5 else Inf
    fits <- list(list(fit(A, "mace4", hz), "A", "pre-registered (pooled)", "A"),
                 list(fit_strat(A, "mace4", hz), "A", "stratified on qualifying condition", "A"),
                 list(fit(Cc, "mace4", hz, w = w), "C standardized", "pre-registered", "C"))
    for (ft in fits) for (nm in names(slopes[[ft[[4]]]]))
      cal[[length(cal) + 1]] <- calibrate(ft[[1]], slopes[[ft[[4]]]][[nm]], ft[[2]], h, ft[[3]], nm)
  }
  rbindlist(cal)
}

# ---------------------------------------------------------------------------
# P5
# ---------------------------------------------------------------------------
band_of <- function(pred) cut(pred, BANDS, right = FALSE, labels = BAND_LABELS)
range_check <- function(D, label) {
  d <- data.table(band = band_of(D$pred), pred = D$pred, status = D$status, lpa = D$lpa)
  d[, .(cohort = label, n = .N, predicted_mean = mean(pred), measured_n = sum(status != "not_measured"),
        above_assay_range_pct = 100 * sum(status == "above_range") / sum(status != "not_measured"),
        measured_ge150_or_above_pct = 100 * sum(status == "above_range" | (status == "measured" & lpa >= 150)) / sum(status != "not_measured"),
        below_assay_range_pct = 100 * sum(status == "below_range") / sum(status != "not_measured")), keyby = band]
}
linearity <- function(D, cohort, model, w = NULL, strat = FALSE) {
  bands <- list(); shape <- list()
  wt <- if (is.null(w)) rep(1, nrow(D)) else w
  bd <- band_of(D$pred)
  for (e in c("mace4", "mace3")) {
    tv <- mk(D, e)
    dm <- data.table(band2 = as.integer(bd == BAND_LABELS[2]), band3 = as.integer(bd == BAND_LABELS[3]), band4 = as.integer(bd == BAND_LABELS[4]))
    fb <- fit_tv(D, tv, xs = c("x_ldl", "x_prs"), strat = strat, w = w, extra = dm)
    cnt <- data.table(band = bd, ev = tv$ev, pred = D$pred, wt = wt)[tv$time > 0 & wt > 0,
             .(n = .N, events = sum(ev), predicted_mean = weighted.mean(pred, wt)), keyby = band]
    cnt[, band := as.character(band)]
    if (!is.null(fb)) {
      bb <- c(0, fb$b[c("band2", "band3", "band4")]); ss <- c(0, sqrt(diag(fb$V)[c("band2", "band3", "band4")]))
      r <- data.table(cohort = cohort, model = model, endpoint = e, band = BAND_LABELS, b = bb, se = ss)
      r <- merge(r, cnt, by = "band", all.x = TRUE, sort = FALSE)
      r[, `:=`(hr = exp(b), hr_lo = exp(b - 1.96 * se), hr_hi = exp(b + 1.96 * se))]
      gap <- (r$predicted_mean[4] - r$predicted_mean[1]) / DELTA
      r[, `:=`(top_band_slope_per_144 = b[4] / gap, top_band_slope_se = se[4] / gap)]
      bands[[length(bands) + 1]] <- r
    }
    for (k in KNOTS) {
      fs <- fit_tv(D, tv, strat = strat, w = w, extra = data.table(x_above = pmax(D$x_lpa - k / DELTA, 0)))
      if (is.null(fs)) next
      b1 <- fs$b[["x_lpa"]]; b2 <- fs$b[["x_above"]]; V <- fs$V
      ba <- b1 + b2; sa <- sqrt(V["x_lpa", "x_lpa"] + V["x_above", "x_above"] + 2 * V["x_lpa", "x_above"])
      shape[[length(shape) + 1]] <- data.table(cohort = cohort, model = model, endpoint = e, test = sprintf("spline, knot %d nmol/L", k),
        events = as.integer(fs$events), n_above = sum(D$pred >= k), slope_below = b1, se_below = sqrt(V["x_lpa", "x_lpa"]),
        slope_above = ba, se_above = sa, change = b2, se_change = sqrt(V["x_above", "x_above"]),
        p_nonlinear = 2 * pnorm(-abs(b2 / sqrt(V["x_above", "x_above"]))),
        hr_above_per_144_lower = exp(-ba), above_lo = exp(-ba - 1.96 * sa), above_hi = exp(-ba + 1.96 * sa),
        p_above_vs_trial = 2 * pnorm(-abs((-ba - log(TRIAL[["hr"]])) / sqrt(sa^2 + TRIAL_SE^2))))
    }
    fq <- fit_tv(D, tv, strat = strat, w = w, extra = data.table(x_sq = D$x_lpa^2))
    if (!is.null(fq)) shape[[length(shape) + 1]] <- data.table(cohort = cohort, model = model, endpoint = e, test = "quadratic term",
        events = as.integer(fq$events), change = fq$b[["x_sq"]], se_change = sqrt(fq$V["x_sq", "x_sq"]),
        p_nonlinear = 2 * pnorm(-abs(fq$b[["x_sq"]] / sqrt(fq$V["x_sq", "x_sq"]))))
  }
  list(bands = rbindlist(bands), shape = rbindlist(shape, fill = TRUE))
}

run_posthoc2 <- function(P, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  P <- P[complete.cases(P[, c(EXPOSURES, ADJ), with = FALSE])]
  A <- add_groups(P[qual == TRUE]); Cc <- P[broad == FALSE]; w <- std_weights(Cc, A)
  p3 <- stroke_checks(A, Cc, w)
  p4 <- statin_calibration(A, Cc, w)
  rc <- rbind(range_check(A, "A"), range_check(Cc, "C"))
  l1 <- linearity(A, "A", "pre-registered (pooled)"); l2 <- linearity(A, "A", "stratified on qualifying condition", strat = TRUE)
  l3 <- linearity(Cc, "C standardized", "pre-registered", w = w)
  bands <- rbind(l1$bands, l2$bands, l3$bands); shape <- rbind(l1$shape, l2$shape, l3$shape, fill = TRUE)
  wcsv(p3$budget, file.path(outdir, "stroke_budget.csv")); wcsv(p3$models, file.path(outdir, "stroke_models.csv"))
  wcsv(p3$factors, file.path(outdir, "stroke_factors.csv")); wcsv(p4, file.path(outdir, "calibration_statin.csv"))
  wcsv(rc, file.path(outdir, "range_check.csv")); wcsv(bands, file.path(outdir, "linearity_bands.csv"))
  wcsv(shape, file.path(outdir, "linearity_shape.csv"))
  list(budget = p3$budget, models = p3$models, factors = p3$factors, cal = p4, range = rc, bands = bands, shape = shape)
}

# atrial fibrillation before baseline and baseline systolic pressure, from the raw field files
add_stroke_factors <- function(P, raw_dir) {
  rd <- function(fid, cols) {
    f <- file.path(raw_dir, sprintf("fid%s.csv", fid)); h <- names(fread(f, nrows = 0))
    check(all(cols %in% h), sprintf("field %s lacks the expected columns", fid))
    d <- fread(f, select = c(h[1], cols), colClasses = "character"); setnames(d, 1, "IID"); d[IID %in% P$IID]
  }
  af <- rd("131350", "131350-0.0"); af[, d := days(`131350-0.0`)]
  af <- af[!is.na(d) & d < as.integer(as.IDate("2037-01-01"))]
  P[, af_date := af$d[match(IID, af$IID)]]
  P[, af_pre := as.integer(!is.na(af_date) & af_date <= t0)]; P[, af_date := NULL]   # 0/1, so the Cox term is named af_pre
  sb <- rd("4080", c("4080-0.0", "4080-0.1"))
  v <- rowMeans(cbind(num(sb[["4080-0.0"]]), num(sb[["4080-0.1"]])), na.rm = TRUE); v[is.nan(v)] <- NA_real_
  P[, sbp := v[match(IID, sb$IID)]]
  P
}

selftest_posthoc2 <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { check(cond, paste("selftest:", msg)); n_ok <<- n_ok + 1L }
  base <- function(n, freq) {
    D <- data.table(IID = as.character(seq_len(n)), pred = pmax(0, rgamma(n, shape = 1.2, scale = 14) + 110 * rbinom(n, 2, freq)),
                    score_ldl_raw = rnorm(n), prs_n = rnorm(n), score_lpa_raw = rnorm(n), age = round(runif(n, 40, 70)), sex = rbinom(n, 1, 0.5))
    for (k in 1:10) set(D, j = paste0("PC", k), value = rnorm(n))
    D <- finish_exposures(D, NULL)$P
    D[, `:=`(t0 = 14000L, t_end = 18000L, mi_pre = TRUE, is_pre = FALSE, pad_pre = FALSE,
             t_mi = NA_integer_, t_is = NA_integer_, t_cvd = NA_integer_, t_rev = NA_integer_)]
    add_groups(D)
  }
  at <- function(D, rate) { te <- D$t0 + as.integer(ceiling(365.25 * rexp(nrow(D), rate))); fifelse(te <= D$t_end, te, NA_integer_) }

  # P4: a known untreated slope, statin users at 60% of their untreated LDL-C
  set.seed(21); D <- base(200000, 0.08); D[, statin := rbinom(.N, 1, 0.5) == 1]
  D[, ldl := (3.5 + 0.30 * x_ldl + rnorm(.N, 0, 0.8)) * fifelse(statin, 0.6, 1)]
  s <- sapply(c(1, 0.7, 0.6, 0.5), function(d) ldl_slope_d(D, d))
  expect(abs(s[3] - 0.30) < 0.01, "dividing by the true statin effect recovers the untreated LDL-C slope")
  expect(all(diff(s) > 0), "a stronger correction gives a larger slope")
  expect(abs(s[1] - ldl_slope(D)) < 1e-12 && abs(s[2] - ldl_slope(D, TRUE)) < 1e-12, "divisors 1 and 0.7 reproduce the first post hoc round")
  expect(abs(ldl_slope_d(D, keep = !D$statin) - 0.30) < 0.01, "the slope among non-users is the untreated slope when treatment is unrelated to LDL-C")

  # P5: linear truth, then a slope that is flat above the knot
  set.seed(22); D <- base(200000, 0.2); k <- 150 / DELTA
  D[, t_mi := at(D, 0.05 * exp(0.25 * x_lpa + 0.09 * x_ldl))]
  L <- linearity(D, "sim", "linear")
  sp <- L$shape[endpoint == "mace4" & test == "spline, knot 150 nmol/L"]; qd <- L$shape[endpoint == "mace4" & test == "quadratic term"]
  bn <- L$bands[endpoint == "mace4"]
  message(sprintf("[posthoc2 sim] linear truth 0.25: below %.3f, above %.3f (SE %.3f), change %.3f (SE %.3f); top-band slope %.3f (SE %.3f)",
                  sp$slope_below, sp$slope_above, sp$se_above, sp$change, sp$se_change, bn$top_band_slope_per_144[1], bn$top_band_slope_se[1]))
  expect(abs(sp$change) < 3.5 * sp$se_change && abs(qd$change) < 3.5 * qd$se_change, "a linear association shows no change of slope and no quadratic term")
  expect(abs(sp$slope_above - 0.25) < 3.5 * sp$se_above && abs(sp$slope_below - 0.25) < 3.5 * sp$se_below, "both spline slopes recover the planted slope")
  expect(abs(bn$top_band_slope_per_144[1] - 0.25) < 3.5 * bn$top_band_slope_se[1] && nrow(bn) == 4 && bn$b[1] == 0, "the top-band slope recovers the planted slope")
  expect(abs(sp$hr_above_per_144_lower - exp(-sp$slope_above)) < 1e-12 && sp$n_above == sum(D$pred >= 150), "trial-direction hazard ratio and count above the knot are consistent")
  D[, t_mi := at(D, 0.05 * exp(0.25 * x_lpa - 0.25 * pmax(x_lpa - k, 0) + 0.09 * x_ldl))]
  sp2 <- linearity(D, "sim", "flat above")$shape[endpoint == "mace4" & test == "spline, knot 150 nmol/L"]
  message(sprintf("[posthoc2 sim] flat above the knot: above %.3f (SE %.3f), change %.3f (SE %.3f)", sp2$slope_above, sp2$se_above, sp2$change, sp2$se_change))
  expect(sp2$change / sp2$se_change < -3 && abs(sp2$slope_above) < 3.5 * sp2$se_above, "a slope that stops at the knot is detected")
  rc <- range_check(D[, `:=`(status = "measured", lpa = pred)], "sim")
  expect(nrow(rc) == 4 && rc[band == ">=150", measured_ge150_or_above_pct] == 100 && rc[band == "<50", measured_ge150_or_above_pct] == 0, "range check bands are assigned correctly")

  # P3a/b: selection on a stroke risk factor manufactures an inverse association; adjustment removes it
  set.seed(23); Pop <- base(300000, 0.2); Pop[, af_pre := rbinom(.N, 1, 0.3)]
  Pop[, `:=`(bmi = rnorm(.N, 27, 4), smoking_current = rbinom(.N, 1, 0.1), diabetes = rbinom(.N, 1, 0.05), hypertension = rbinom(.N, 1, 0.3), sbp = rnorm(.N, 140, 18))]
  A <- Pop[rbinom(nrow(Pop), 1, plogis(-3 + 2 * x_lpa + 2.5 * af_pre)) == 1]
  A[, t_is := at(A, 0.02 * exp(1.2 * af_pre))]
  un <- fit_tv(A, mk(A, "stroke")); ad <- fit_tv(A, mk(A, "stroke"), adj = c(ADJ, STROKE_RISK))
  bg <- budget(A, Pop, "stroke", STROKE_RISK); imp <- bg[factor == "af_pre", implied_bias_log_hr]; moved <- un$b[["x_lpa"]] - ad$b[["x_lpa"]]
  message(sprintf("[posthoc2 sim] selection: Lp(a) log HR unadjusted %.3f, adjusted %.3f (planted 0); budget %.3f, shift on adjustment %.3f",
                  un$b[["x_lpa"]], ad$b[["x_lpa"]], imp, moved))
  expect(un$b[["x_lpa"]] / sqrt(un$V[1, 1]) < -3, "selection on a stroke risk factor manufactures an inverse association")
  expect(abs(ad$b[["x_lpa"]]) < 3.5 * sqrt(ad$V[1, 1]), "adjusting for the factor removes it")
  expect(imp < 0 && imp / moved > 0.4 && imp / moved < 2.5, "the budget has the sign and the order of magnitude of the bias")
  expect(abs(bg[factor == "af_pre", slope_primary]) < 0.01 && bg[factor == "af_pre", induced] < -0.05, "the induced association is absent in the reference cohort and present after selection")

  # P3c: stroke rate falls after a coronary event that Lp(a) brings forward
  set.seed(24); D <- base(150000, 0.2)
  tc <- D$t0 + as.integer(ceiling(365.25 * rexp(nrow(D), 0.04 * exp(0.6 * D$x_lpa))))
  ts1 <- D$t0 + as.integer(ceiling(365.25 * rexp(nrow(D), 0.02)))
  ts <- ifelse(ts1 <= tc, ts1, tc + as.integer(ceiling(365.25 * rexp(nrow(D), 0.02 * exp(-1.5)))))
  D[, `:=`(t_mi = fifelse(tc <= t_end, tc, NA_integer_), t_is = fifelse(ts <= t_end, as.integer(ts), NA_integer_))]
  al <- fit_tv(D, mk(D, "stroke")); fe <- fit_tv(D, mk_first(D, "stroke", CORONARY))
  message(sprintf("[posthoc2 sim] displacement: Lp(a) log HR on all strokes %.3f (SE %.3f), on first-event strokes %.3f (SE %.3f); planted direct effect 0",
                  al$b[["x_lpa"]], sqrt(al$V[1, 1]), fe$b[["x_lpa"]], sqrt(fe$V[1, 1])))
  expect(al$b[["x_lpa"]] / sqrt(al$V[1, 1]) < -3, "protection after a coronary event manufactures an inverse association with all strokes")
  expect(abs(fe$b[["x_lpa"]]) < 3.5 * sqrt(fe$V[1, 1]) && fe$events < al$events, "censoring at the coronary event removes it")
  tv <- mk_first(D, "stroke", CORONARY)
  expect(all(tv$time <= (pmin(D$t_end, D$t_mi, na.rm = TRUE) - D$t0) / 365.25 + 1e-9) && sum(tv$ev) == sum(!is.na(D$t_is) & (is.na(D$t_mi) | D$t_is <= D$t_mi)),
         "first-event follow-up stops at the coronary event and counts only strokes up to it")

  # smoke run of everything the real run executes
  s <- simulate(300000, 0.25, 0.09, seed = 302); P <- s$P; set.seed(25)
  P[, `:=`(ldl = (3.5 + 0.30 * x_ldl + rnorm(.N, 0, 0.8)) * fifelse(statin, 0.7, 1), af_pre = rbinom(.N, 1, 0.05), sbp = rnorm(.N, 140, 18))]
  P[sample(.N, 5000), sbp := NA_real_]
  out <- run_posthoc2(P, file.path(workdir, "posthoc2_smoke"))
  files <- c("stroke_budget.csv", "stroke_models.csv", "stroke_factors.csv", "calibration_statin.csv", "range_check.csv", "linearity_bands.csv", "linearity_shape.csv")
  expect(all(file.exists(file.path(workdir, "posthoc2_smoke", files))), "every output file is written")
  expect(nrow(out$models) == 14 && nrow(out$cal) == 34 && nrow(out$range) == 8 && nrow(out$bands[cohort == "A"]) == 16, "all models, scales and bands are present")
  ca <- out$cal[cohort == "A" & horizon == "full" & model == "pre-registered (pooled)"]
  expect(abs(ca[ldl_scale == "LDL-C / 0.7 in statin users", ldl_mmol_per_sd] - 0.30) < 0.03 &&
           ca[ldl_scale == "LDL-C / 0.5 in statin users", predicted_trial_hr] < ca[ldl_scale == "measured LDL-C", predicted_trial_hr],
         "the smoke run recovers the planted LDL-C slope and a stronger correction lowers the predicted trial hazard ratio")
  cat(sprintf("\nPOSTHOC2 SELFTEST PASSED: %d checks\n", n_ok))
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--step", type = "character", default = "selftest"),
  make_option("--raw_dir", type = "character"), make_option("--covar", type = "character"),
  make_option("--dose_lpa", type = "character"), make_option("--w_lpa", type = "character"),
  make_option("--dose_ldl", type = "character"), make_option("--w_ldl", type = "character"),
  make_option("--prs", type = "character"), make_option("--outdir", type = "character"),
  make_option("--workdir", type = "character", default = tempdir()),
  make_option("--censor", type = "character", default = "2022-10-31"),
  make_option("--threads", type = "integer", default = 8L))))
setDTthreads(opt$threads)
if (opt$step == "selftest") { selftest_posthoc2(opt$workdir); quit(save = "no", status = 0L) }
check(opt$step == "run", "unknown --step")
CENS <- as.integer(as.IDate(opt$censor))
R <- load_real(opt, CENS)
if (!gates(R, CENS, file.path(opt$outdir, "gates_recheck"))) { cat("\nGATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
cat("\nGATES PASSED\n")
out <- run_posthoc2(add_stroke_factors(R$P, opt$raw_dir), opt$outdir)
options(width = 250)
fmt <- function(b, se) sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
cat("\n=== P3 stroke risk factors at baseline (no outcome) ===\n"); print(out$factors, digits = 3)
cat("\n=== P3a selection budget by endpoint: implied bias in the Lp(a) log hazard ratio from unmodelled factors ===\n")
print(dcast(out$budget[in_primary_model == FALSE], endpoint ~ factor, value.var = "implied_bias_log_hr"), digits = 2)
print(out$budget[in_primary_model == FALSE, .(unmodelled_bias_log_hr = sum(implied_bias_log_hr)), by = endpoint], digits = 3)
cat("\n--- factors for the stroke endpoint ---\n"); print(out$budget[endpoint == "stroke"], digits = 3)
cat("\n=== P3b/c later ischaemic stroke, Lp(a) hazard ratio per 144 nmol/L HIGHER ===\n")
print(out$models[, .(cohort, follow_up, model, adjustment, n, events, hr = fmt(b, se), p)], digits = 3)
cat("\n=== P4 calibration under different statin corrections (trial observed 1.02, 0.90-1.15) ===\n")
print(out$cal[, .(cohort, horizon, model, ldl_scale, ldl_mmol_per_sd, genetic_hr_per_mmol_ldl, genetic_to_trial_ratio, equivalent_ldl_mmol, equiv_lo, equiv_hi,
                  predicted_trial_hr, pred_lo, pred_hi, p_trial_vs_predicted)], digits = 3)
cat("\n=== P5a measured Lp(a) by band of predicted Lp(a) (no outcome) ===\n"); print(out$range, digits = 3)
cat("\n=== P5b hazard ratio by band of predicted Lp(a), against the lowest band ===\n")
print(out$bands[, .(cohort, model, endpoint, band, n, events, predicted_mean, hr = fmt(b, se), top_band_slope_hr_per_144 = fmt(top_band_slope_per_144, top_band_slope_se))], digits = 3)
cat("\n=== P5c/d shape tests (slopes are log hazard ratios per 144 nmol/L higher) ===\n")
print(out$shape, digits = 3)
cat("\nPOSTHOC2_DONE\n")
