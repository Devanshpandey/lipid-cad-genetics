#!/usr/bin/env Rscript
# ============================================================================
# Matched-endpoint analysis for the genetic comparison with Lp(a)HORIZON.
# Specification: PREREGISTERED_READOUT.md, section "Matched-endpoint analysis",
# committed before this is run on real data. Post hoc; it cannot change the
# verdict of record.
#
# The trial counted urgent coronary revascularization; the comparison counts any.
# Record-level episodes are not held, so a revascularization is called
# "acute-linked" when an acute coronary code (I20.0, I21, I22, I24) was first
# recorded from W days before it to one day after it. The matched composite is
# cardiovascular death, infarction, ischaemic stroke or acute-linked
# revascularization.
#
# lpa_horizon.R is sourced unchanged; the function definitions and named
# constants of the post hoc scripts and of lpa_forecast.R are evaluated from
# their parsed files (their main blocks are unguarded, so they are not sourced).
#
#   --step selftest   synthetic data only
#   --step run        gates (must pass), then the estimates
# ============================================================================
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
source(file.path(here, "lpa_horizon.R"))
take <- function(file, consts) for (e in parse(file.path(here, file))) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      (as.character(e[[2]]) %in% consts || (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))))) eval(e, globalenv())
}
take("lpa_horizon_posthoc.R", c("CTT", "TRIAL_SE"))
take("lpa_horizon_posthoc2.R", c("KNOTS"))
take("lpa_forecast.R", c("KNOT", "DIV", "NDRAW", "TRIAL3"))
check(all(vapply(c("fit_tv", "fit_strat", "add_groups", "ldl_slope_d", "forecast", "rdraw", "se_from", "CTT", "KNOT", "DIV", "NDRAW", "TRIAL3"), exists, logical(1))),
      "the earlier scripts did not provide their functions")

ACS <- "^(I200|I21|I22|I24)"          # unstable angina, acute and subsequent infarction, other acute ischaemic heart disease
WINDOWS <- c(7L, 14L, 30L)            # days from the acute code to the procedure; 7 is the definition, the others sensitivity
L0 <- 229.4; L1 <- 65.2               # the trial's time-averaged post-baseline geometric means, nmol/L
TRIAL_UREV <- c(hr = 1.04, lo = 0.80, hi = 1.36)   # urgent coronary revascularization, Cho et al. 2026, Table 2
RECORD <- list(b4 = 0.174533790354549, b3 = 0.101635091628793)   # Lp(a) coefficients of record, four- and three-point, cohort A
COUNTS <- list(urev7 = 786L, elrev = 1180L, m4u7 = 4078L)        # cohort A, counts-only job 3602802
for (W in WINDOWS) {
  ENDPOINTS[[paste0("urev", W)]] <- paste0("t_u", W)
  ENDPOINTS[[paste0("m4u", W)]] <- c("t_mi", "t_is", "t_cvd", paste0("t_u", W))
}
ENDPOINTS$elrev <- "t_el"
MATCHED <- paste0("m4u", WINDOWS); LINKED <- paste0("urev", WINDOWS)

# every revascularization with the acute coronary code first recorded closest before it
# (lag = days from the code to the procedure; -1 allows a code dated the day after)
link_procedures <- function(icd, ops) {
  acs <- unique(icd[grepl(ACS, code), .(IID, acode = code, da = d)])
  rev <- unique(ops[grepl(RX$REV, code), .(IID, pcode = code, dp = d)])
  lk <- merge(rev, acs, by = "IID", allow.cartesian = TRUE)
  lk[, lag := dp - da]; lk <- lk[lag >= -1L]
  setorder(lk, IID, pcode, dp, lag)
  near <- lk[, .SD[1L], by = .(IID, pcode, dp)]
  merge(rev, near[, .(IID, pcode, dp, lag, acode)], by = c("IID", "pcode", "dp"), all.x = TRUE)
}
# dates of the first acute-linked and of the first other revascularization after time zero
add_linked <- function(D, rev, CENS) {
  D <- copy(D)
  r <- merge(rev, D[, .(IID, t0, t_end)], by = "IID")[dp > t0 & dp <= pmin(t_end, CENS)]
  frst <- function(x) {
    if (!nrow(x)) return(rep(NA_integer_, nrow(D)))
    y <- x[, .(t = min(dp)), by = IID]; as.integer(y$t[match(D$IID, y$IID)])
  }
  check(identical(frst(r), as.integer(D$t_rev)), "revascularization dates rebuilt from the codes differ from the endpoint of record")
  for (W in WINDOWS) set(D, j = paste0("t_u", W), value = frst(r[!is.na(lag) & lag <= W]))
  set(D, j = "t_el", value = frst(r[is.na(lag) | lag > max(WINDOWS)]))
  D
}
read_codes <- function(raw_dir, ids) {
  arr0 <- function(fid) {
    f <- file.path(raw_dir, sprintf("fid%s.csv", fid)); h <- names(fread(f, nrows = 0))
    cols <- grep(sprintf("^%s-0\\.[0-9]+$", fid), h, value = TRUE)
    check(length(cols) > 0, sprintf("field %s has no instance-0 columns", fid))
    d <- fread(f, select = c(h[1], cols), colClasses = "character"); setnames(d, 1, "IID"); d <- d[IID %in% ids]
    rbindlist(lapply(cols, function(cn) { v <- d[[cn]]; k <- ok(v)
      data.table(IID = d$IID[k], idx = rep(as.integer(sub("^.*\\.", "", cn)), sum(k)), v = v[k]) }))
  }
  dated0 <- function(fc, fd) {
    a <- arr0(fc); setnames(a, "v", "code"); b <- arr0(fd); setnames(b, "v", "date")
    x <- merge(a, b, by = c("IID", "idx")); x[, d := days(date)]
    x[!is.na(d), .(IID, code, d)]
  }
  list(icd = dated0("41270", "41280"), ops = dated0("41272", "41282"))
}

rows <- function(f, cohort, endpoint, horizon, model) {
  r <- est_rows(f, cohort, endpoint, horizon, model)
  if (is.null(r)) return(NULL)
  r[, `:=`(hr_lower = exp(-b), lower_lo = exp(-b - 1.96 * se), lower_hi = exp(-b + 1.96 * se))]
  r
}
# two event types in the same participants: one Cox model with a stratum per type and every
# covariate free to differ by type, variance clustered on participant
stacked_diff <- function(D, e1, e2, xs = EXPOSURES, adj = ADJ) {
  vars <- c(xs, adj)
  part <- function(e, k) { x <- cbind(D[, c("IID", vars), with = FALSE], mk(D, e)); x[, type := k]; x }
  X <- rbind(part(e1, 1L), part(e2, 2L)); X <- X[time > 0]; X <- X[complete.cases(X)]
  if (X[type == 1L, sum(ev)] < 20L || X[type == 2L, sum(ev)] < 20L) return(NULL)
  for (v in vars) for (k in 1:2) set(X, j = paste0(v, "_", k), value = X[[v]] * (X$type == k))
  terms <- as.vector(outer(vars, 1:2, paste, sep = "_"))
  m <- coxph(as.formula(paste("Surv(time, ev) ~", paste(terms, collapse = " + "), "+ strata(type) + cluster(IID)")), data = X)
  b <- coef(m); V <- vcov(m); if (anyNA(b)) return(NULL)
  rbindlist(lapply(xs, function(x) {
    k1 <- paste0(x, "_1"); k2 <- paste0(x, "_2"); d <- b[[k1]] - b[[k2]]; se <- sqrt(V[k1, k1] + V[k2, k2] - 2 * V[k1, k2])
    data.table(exposure = x, first = e1, second = e2, events_first = X[type == 1L, sum(ev)], events_second = X[type == 2L, sum(ev)],
               b_first = b[[k1]], se_first = sqrt(V[k1, k1]), b_second = b[[k2]], se_second = sqrt(V[k2, k2]),
               diff = d, se_diff = se, p_diff = 2 * pnorm(-abs(d / se)))
  }))
}
# hazard ratio predicted for the trial's fall in Lp(a), with one slope and with a second slope above the knot
predict_trial <- function(D, endpoint, observed, label) {
  tv <- mk(D, endpoint)
  lin <- fit_tv(D, tv); spl <- fit_tv(D, tv, extra = data.table(x_above = pmax(D$x_lpa - KNOT / DELTA, 0)))
  if (is.null(lin) || is.null(spl)) return(NULL)
  s <- ldl_slope_d(D, DIV)
  r <- rbind(forecast(lin, rdraw(lin), s, "linear", L0, L1), forecast(spl, rdraw(spl), s, "spline", L0, L1))
  r[, `:=`(endpoint = endpoint, events = as.integer(lin$events), compared_with = label,
           observed = observed[["hr"]], observed_lo = observed[["lo"]], observed_hi = observed[["hi"]])]
  r[, z := (log(observed) - log(hr)) / sqrt(se_from(observed_lo, observed_hi)^2 + se_from(lo, hi)^2)]
  r[, p_observed_vs_predicted := 2 * pnorm(-abs(z))]
  r[, `:=`(b_lpa = lin$b[["x_lpa"]], se_lpa = sqrt(lin$V["x_lpa", "x_lpa"]),
           p_change_of_slope = 2 * pnorm(-abs(spl$b[["x_above"]] / sqrt(spl$V["x_above", "x_above"]))))]
  r
}

run_matched <- function(P, I, rev, CENS, outdir, real = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE); set.seed(42)
  cc <- function(D) D[complete.cases(D[, c(EXPOSURES, ADJ), with = FALSE])]
  P <- cc(P); I <- cc(I)
  A <- add_groups(add_linked(P[qual == TRUE], rev, CENS)); Cc <- add_linked(P[broad == FALSE], rev, CENS); Dd <- add_linked(I, rev, CENS)
  w <- std_weights(Cc, A)
  ends <- c("mace4", "mace3", MATCHED, "revasc", LINKED, "elrev")
  coh <- list(A = A, C = Cc, D = Dd)
  counts <- rbindlist(lapply(names(coh), function(k) rbindlist(lapply(ends, function(e)
    data.table(cohort = k, endpoint = e, n = nrow(coh[[k]]), events = sum(mk(coh[[k]], e)$ev), events_5y = sum(mk(coh[[k]], e, H5)$ev))))))
  fA4 <- fit(A, "mace4"); fA3 <- fit(A, "mace3")
  if (real) {
    ev <- function(e) counts[cohort == "A" & endpoint == e, events]
    g <- data.table(gate = c("acute-linked revascularizations, cohort A", "other revascularizations, cohort A", "matched composite events, cohort A",
                             "four-point Lp(a) coefficient of record", "three-point Lp(a) coefficient of record"),
                    value = c(ev("urev7"), ev("elrev"), ev("m4u7"), fA4$b[["x_lpa"]], fA3$b[["x_lpa"]]),
                    target = c(COUNTS$urev7, COUNTS$elrev, COUNTS$m4u7, RECORD$b4, RECORD$b3))
    g[, pass := abs(value - target) < 1e-9]
    fwrite(g, file.path(outdir, "gates_matched.csv")); print(g, digits = 15)
    if (!all(g$pass)) return(NULL)
    cat("\nMATCHED GATES PASSED\n")
  }

  est <- list()
  for (hn in c("full", "5y")) {
    h <- if (hn == "full") Inf else H5
    for (e in ends) est <- c(est, list(
      rows(fit(A, e, h), "A", e, hn, "pooled"), rows(fit_strat(A, e, h), "A", e, hn, "stratified on qualifying condition"),
      rows(fit(Cc, e, h, w = w), "C standardized", e, hn, "pooled"), rows(fit(Dd, e, h), "D", e, hn, "pooled")))
  }
  est <- rbindlist(est)

  contrast <- rbindlist(lapply(c("mace4", "mace3", MATCHED), function(e) {
    fA <- fit(A, e); fC <- fit(Cc, e, w = w)
    if (is.null(fA) || is.null(fC)) return(NULL)
    rs <- rbindlist(lapply(EXPOSURES, function(x) ratio_stats(fA, fC, x))); v <- verdict(rs, spec_test(fA, fC))
    cbind(endpoint = e, rs, verdict_by_the_rule = v$verdict, specificity_p = v$specificity_p)
  }))

  pred <- rbindlist(c(
    list(predict_trial(A, "mace4", TRIAL, "trial primary, four-point"), predict_trial(A, "mace3", TRIAL3, "trial three-point"),
         predict_trial(A, "mace3", TRIAL, "trial primary, four-point")),
    lapply(MATCHED, function(e) predict_trial(A, e, TRIAL, "trial primary, four-point")),
    lapply(c(LINKED, "revasc", "elrev"), function(e) predict_trial(A, e, TRIAL_UREV, "trial urgent coronary revascularization"))))

  dd <- function(k, e1, e2) { x <- stacked_diff(coh[[k]], e1, e2); if (is.null(x)) NULL else cbind(cohort = k, x) }
  diffs <- rbindlist(list(dd("A", "urev7", "elrev"), dd("A", "urev30", "elrev"), dd("D", "urev7", "elrev"), dd("C", "urev7", "elrev")))

  wcsv(counts, file.path(outdir, "matched_counts.csv")); wcsv(est, file.path(outdir, "matched_estimates.csv"))
  wcsv(contrast, file.path(outdir, "matched_stage_contrast.csv")); wcsv(pred, file.path(outdir, "matched_predictions.csv"))
  wcsv(diffs, file.path(outdir, "revascularization_difference.csv"))
  list(counts = counts, est = est, contrast = contrast, pred = pred, diffs = diffs, A = A)
}

selftest_matched <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { check(cond, paste("selftest:", msg)); n_ok <<- n_ok + 1L }
  near <- function(a, b, tol = 1e-9) isTRUE(all(abs(a - b) < tol))
  d <- function(s) as.integer(as.IDate(s))

  # 1. linking on a hand-built set of codes
  icd <- data.table(IID = c("1", "2", "3", "4", "6", "7", "8", "11"), code = c("I219", "I214", "I200", "I219", "I214", "I219", "I248", "I219"),
                    d = d(c("2015-03-01", "2015-03-02", "2015-02-10", "2005-01-01", "2016-01-01", "2015-03-06", "2015-03-01", "2015-03-01")))
  ops <- data.table(IID = c("1", "2", "3", "4", "5", "6", "6", "7", "8", "9", "11"),
                    code = c("K751", "K751", "K491", "K751", "K401", "K751", "K401", "K751", "K502", "X999", "K751"),
                    d = d(c("2015-03-03", "2015-03-01", "2015-03-01", "2015-03-01", "2015-03-01", "2012-01-05", "2016-01-05", "2015-03-01", "2015-03-04",
                            "2015-03-01", "2015-03-01")))
  rev <- link_procedures(icd, ops)
  lag_of <- function(id, pc) rev[IID == id & pcode == pc, lag]
  expect(nrow(rev) == 10 && !("9" %in% rev$IID), "only coronary revascularization codes are kept")
  expect(lag_of("1", "K751") == 2L && lag_of("2", "K751") == -1L && lag_of("3", "K491") == 19L && lag_of("8", "K502") == 3L,
         "lags: two days before, the day after, nineteen days before, another acute code")
  expect(lag_of("4", "K751") == d("2015-03-01") - d("2005-01-01") && is.na(lag_of("5", "K401")) && is.na(lag_of("7", "K751")),
         "an old code gives a long lag; no code, or a code more than a day after the procedure, gives none")
  expect(is.na(lag_of("6", "K751")) && lag_of("6", "K401") == 4L, "each procedure of a participant is linked separately")
  D0 <- data.table(IID = c("1", "2", "3", "4", "5", "6", "7", "8", "10", "11"), t0 = d("2008-06-01"), t_end = d("2020-01-01"),
                   t_rev = d(c("2015-03-03", "2015-03-01", "2015-03-01", "2015-03-01", "2015-03-01", "2012-01-05", "2015-03-01", "2015-03-04", NA, NA)))
  D0[IID == "11", t_end := d("2014-01-01")]
  CENS0 <- d("2022-10-31"); L <- add_linked(D0, rev, CENS0)
  v <- function(cn) setNames(L[[cn]], L$IID)
  expect(identical(unname(v("t_u7")), d(c("2015-03-03", "2015-03-01", NA, NA, NA, "2016-01-05", NA, "2015-03-04", NA, NA))),
         "acute-linked within 7 days, including a later procedure of the same participant")
  expect(v("t_u30")[["3"]] == d("2015-03-01") && is.na(v("t_u14")[["3"]]), "a 19-day lag is linked at 30 days and not at 14")
  expect(identical(unname(v("t_el")), d(c(NA, NA, NA, "2015-03-01", "2015-03-01", "2012-01-05", "2015-03-01", NA, NA, NA))),
         "other revascularizations are those with no acute code within 30 days")
  expect(is.na(v("t_u7")[["11"]]) && is.na(v("t_el")[["11"]]), "a procedure after the end of follow-up is not counted")
  D1 <- copy(D0); D1[IID == "5", t_rev := t_rev + 1L]
  expect(inherits(try(add_linked(D1, rev, CENS0), silent = TRUE), "try-error"), "a cohort whose endpoint differs from the codes is refused")

  # 2. the stacked model recovers planted slopes and equals two separate models
  sim2 <- function(n, b1, b2, seed) {
    set.seed(seed)
    D <- data.table(IID = as.character(seq_len(n)), x_lpa = pmax(0, rgamma(n, 1.2, scale = 14) + 110 * rbinom(n, 2, 0.08)) / DELTA,
                    x_ldl = rnorm(n), x_prs = rnorm(n), age = round(runif(n, 40, 70)), sex = rbinom(n, 1, 0.5))
    for (k in 1:10) set(D, j = paste0("PC", k), value = rnorm(n))
    D[, `:=`(t0 = 0L, t_end = as.integer(13 * 365.25), statin = rbinom(n, 1, 0.7) == 1)]
    D[, ldl := (3.5 + 0.30 * x_ldl + rnorm(n, 0, 0.8)) * fifelse(statin, DIV, 1)]
    lin <- 0.08 * D$x_ldl + 0.2 * D$x_prs + 0.3 * D$sex
    t1 <- as.integer(ceiling(365.25 * rexp(n, 0.006 * exp(b1 * D$x_lpa + lin)))); t2 <- as.integer(ceiling(365.25 * rexp(n, 0.008 * exp(b2 * D$x_lpa + lin))))
    D[, `:=`(t_u7 = fifelse(t1 <= t_end, t1, NA_integer_), t_el = fifelse(t2 <= t_end, t2, NA_integer_))]
    D
  }
  S <- sim2(200000, 0.10, 0.40, 11); sd_ <- stacked_diff(S, "urev7", "elrev")[exposure == "x_lpa"]
  f1 <- fit(S, "urev7"); f2 <- fit(S, "elrev")
  message(sprintf("[matched sim] acute-linked %.3f (planted 0.10), other %.3f (planted 0.40), difference %.3f (SE %.3f)", sd_$b_first, sd_$b_second, sd_$diff, sd_$se_diff))
  expect(near(sd_$b_first, f1$b[["x_lpa"]], 1e-6) && near(sd_$b_second, f2$b[["x_lpa"]], 1e-6), "the stacked model gives the same slopes as two separate models")
  expect(abs(sd_$b_first - 0.10) < 3.5 * sd_$se_first && abs(sd_$b_second - 0.40) < 3.5 * sd_$se_second && sd_$p_diff < 1e-6 && sd_$diff < 0,
         "planted slopes are recovered and their difference is detected")
  expect(abs(sd_$se_first / sqrt(f1$V["x_lpa", "x_lpa"]) - 1) < 0.15, "the clustered standard error is close to the model-based one when events rarely share a participant")
  S0 <- sim2(200000, 0.25, 0.25, 12); s0 <- stacked_diff(S0, "urev7", "elrev")[exposure == "x_lpa"]
  expect(abs(s0$diff / s0$se_diff) < 3.5, "equal slopes give no difference")

  # 3. the prediction for the trial and its comparison with an observed result
  set.seed(5); pr <- predict_trial(S, "elrev", TRIAL_UREV, "x")
  b <- fit_tv(S, mk(S, "elrev"))$b[["x_lpa"]]
  h <- pr[model == "linear" & scale == "lifelong"]
  expect(nrow(pr) == 4 && near(h$hr, exp(-b * (L0 - L1) / DELTA)) && h$lo < h$hr && h$hr < h$hi, "the one-slope prediction is exp(-slope x the trial's fall)")
  expect(near(h$z, (log(1.04) - log(h$hr)) / sqrt(se_from(0.80, 1.36)^2 + se_from(h$lo, h$hi)^2)) && near(h$p_observed_vs_predicted, 2 * pnorm(-abs(h$z))),
         "the comparison with the observed result is a z test on the log scale")

  # 4. the whole pipeline on simulated cohorts
  s <- simulate(300000, 0.125, 0.09, seed = 404); P <- s$P; set.seed(41)
  P[, ldl := (3.5 + 0.30 * x_ldl + rnorm(.N, 0, 0.8)) * fifelse(statin, DIV, 1)]
  pr_ <- P[!is.na(t_rev), .(IID, pcode = "K751", dp = t_rev)]
  pr_[, lag := sample(c(1L, 20L, NA_integer_), .N, replace = TRUE, prob = c(0.4, 0.1, 0.5))]; pr_[, acode := fifelse(is.na(lag), NA_character_, "I219")]
  out <- run_matched(P, s$I, pr_, s$CENS, file.path(workdir, "matched_smoke"))
  files <- c("matched_counts.csv", "matched_estimates.csv", "matched_stage_contrast.csv", "matched_predictions.csv", "revascularization_difference.csv")
  expect(all(file.exists(file.path(workdir, "matched_smoke", files))), "every output file is written")
  cn <- function(e) out$counts[cohort == "A" & endpoint == e, events]
  expect(cn("mace3") < cn("m4u7") && cn("m4u7") < cn("m4u30") && cn("m4u30") < cn("mace4") && cn("urev7") + cn("elrev") < cn("revasc") + 1L,
         "the matched composite lies between the three-point and the four-point composites")
  expect(nrow(out$pred) == 4 * (3 + length(MATCHED) + length(LINKED) + 2) && all(c("linear", "spline") %in% out$pred$model), "every prediction is present")
  e4 <- out$est[cohort == "A" & endpoint == "m4u7" & horizon == "full" & model == "pooled" & exposure == "x_lpa"]
  p4 <- out$pred[endpoint == "m4u7" & model == "linear" & scale == "lifelong"]
  expect(nrow(e4) == 1 && near(p4$hr, exp(-e4$b * (L0 - L1) / DELTA)) && near(e4$hr_lower, exp(-e4$b)), "estimates and predictions agree with each other")
  expect(nrow(out$contrast[endpoint == "m4u7"]) == 3 && nrow(out$diffs[cohort == "A" & exposure == "x_lpa"]) == 2, "stage contrast and revascularization difference are present")
  message(sprintf("[matched sim] matched composite, established %.3f (SE %.3f; planted 0.125); events three-point %d, matched %d, four-point %d",
                  e4$b, e4$se, cn("mace3"), cn("m4u7"), cn("mace4")))
  expect(abs(e4$b - 0.125) < 3.5 * e4$se, "the planted slope is recovered on the matched composite")
  cat(sprintf("\nMATCHED SELFTEST PASSED: %d checks\n", n_ok))
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
if (opt$step == "selftest") { selftest_matched(opt$workdir); quit(save = "no", status = 0L) }
check(opt$step == "run", "unknown --step")
CENS <- as.integer(as.IDate(opt$censor))
R <- load_real(opt, CENS)
if (!gates(R, CENS, file.path(opt$outdir, "gates_recheck"))) { cat("\nGATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
cat("\nGATES PASSED\n")
codes <- read_codes(opt$raw_dir, R$P$IID)
out <- run_matched(R$P, R$I, link_procedures(codes$icd, codes$ops), CENS, opt$outdir, real = TRUE)
if (is.null(out)) { cat("\nMATCHED GATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
options(width = 250)
f3 <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
cat("\n=== events ===\n"); print(dcast(out$counts, endpoint ~ cohort, value.var = "events"))
cat("\n=== Lp(a), hazard ratio per 144 nmol/L lower, full follow-up ===\n")
print(dcast(out$est[exposure == "x_lpa" & horizon == "full", .(endpoint, key = paste(cohort, model), v = f3(hr_lower, lower_lo, lower_hi))], endpoint ~ key, value.var = "v"))
cat("\n=== the same, five-year horizon ===\n")
print(dcast(out$est[exposure == "x_lpa" & horizon == "5y", .(endpoint, key = paste(cohort, model), v = f3(hr_lower, lower_lo, lower_hi))], endpoint ~ key, value.var = "v"))
cat("\n=== LDL-C score and non-lipid score, cohort A pooled, per SD lower ===\n")
print(dcast(out$est[exposure != "x_lpa" & horizon == "full" & cohort == "A" & model == "pooled", .(endpoint, exposure, v = f3(hr_lower, lower_lo, lower_hi))], endpoint ~ exposure, value.var = "v"))
cat("\n=== M2: stage contrast (log HR in cohort A as a share of that in cohort C standardized) ===\n")
print(out$contrast[, .(endpoint, exposure, ratio = f3(ratio, ratio_lo, ratio_hi), p_diff, verdict_by_the_rule, specificity_p)], digits = 3)
cat("\n=== M1 and M3: prediction for the trial (229.4 to 65.2 nmol/L) against what it observed ===\n")
print(out$pred[, .(endpoint, events, model, scale, predicted = f3(hr, lo, hi), compared_with, observed = f3(observed, observed_lo, observed_hi), p_observed_vs_predicted, p_change_of_slope)], digits = 3)
cat("\n=== M3: acute-linked against other revascularization (log HR per unit higher exposure) ===\n")
print(out$diffs, digits = 3)
cat("\nMATCHED_DONE\n")
