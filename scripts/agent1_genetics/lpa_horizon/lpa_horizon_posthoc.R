#!/usr/bin/env Rscript
# ============================================================================
# POST HOC analyses for the Lp(a)HORIZON genetic emulation.
#
# Written after the pre-registered run (LS6 job 3600628, verdict
# PARTIAL_PERSISTENCE) had been read. Described in PREREGISTERED_READOUT.md,
# section "Post hoc analyses". Nothing here changes the verdict of record.
#
# This file sources lpa_horizon.R unchanged (its main block does not run when
# sourced) and reuses its data assembly, exposure model and gates, so the
# cohort is the one the pre-registered run analysed.
#
#   P1  qualifying condition: composition, stratified model, components within
#       each condition
#   P2  LDL-C calibration: what trial hazard ratio the genetic estimate predicts
#
#   --step selftest   synthetic data only
#   --step run        gates (must pass), then P1 and P2
# ============================================================================
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
source(file.path(here, "lpa_horizon.R"))

CTT <- 0.78      # rate ratio for major vascular events per 1 mmol/L lower LDL-C (CTT, Lancet 2010)
TRIAL_SE <- (log(TRIAL[["hi"]]) - log(TRIAL[["lo"]])) / (2 * 1.96)

add_groups <- function(D) {
  D[, qgroup := paste0(as.integer(mi_pre), as.integer(is_pre), as.integer(pad_pre))]
  D
}
# the primary model with a separate baseline hazard for each combination of qualifying conditions
fit_strat <- function(D, endpoint, horizon = Inf) {
  X <- cbind(D[, c(EXPOSURES, ADJ, "qgroup"), with = FALSE], mk(D, endpoint, horizon))
  X <- X[time > 0]; X <- X[complete.cases(X)]
  if (sum(X$ev) < 20L) return(NULL)
  f <- as.formula(paste("Surv(time, ev) ~", paste(c(EXPOSURES, ADJ), collapse = " + "), "+ strata(qgroup)"))
  m <- coxph(f, data = X); b <- coef(m)[EXPOSURES]; if (anyNA(b)) return(NULL)
  list(b = b, V = vcov(m)[EXPOSURES, EXPOSURES, drop = FALSE], n = m$n, events = sum(X$ev))
}
# Fieller interval for a / b; NA when the denominator is not clearly away from zero
fieller <- function(a, b, va, vb, cab, t = 1.96) {
  r <- a / b; g <- t^2 * vb / b^2
  if (g >= 1) return(c(est = r, lo = NA_real_, hi = NA_real_))
  centre <- (r - g * cab / vb) / (1 - g)
  half <- t / (abs(b) * (1 - g)) * sqrt(va - 2 * r * cab + r^2 * vb - g * (va - cab^2 / vb))
  c(est = r, lo = centre - half, hi = centre + half)
}
# difference in LDL-C (mmol/L) per SD of the LDL-C score, adjusted as in the outcome model
ldl_slope <- function(D, corrected = FALSE, w = NULL) {
  y <- if (corrected) ifelse(D$statin, D$ldl / 0.7, D$ldl) else D$ldl
  X <- cbind(D[, c("x_ldl", "x_lpa", ADJ), with = FALSE], y = y, wt = if (is.null(w)) 1 else w)
  X <- X[complete.cases(X) & wt > 0]
  m <- lm(as.formula(paste("y ~", paste(c("x_ldl", "x_lpa", ADJ), collapse = " + "))), data = X, weights = wt)
  coef(m)[["x_ldl"]]
}
calibrate <- function(f, s, cohort, horizon, model, scale) {
  if (is.null(f)) return(NULL)
  bL <- f$b[["x_lpa"]]; bD <- f$b[["x_ldl"]]
  fr <- fieller(bL, bD, f$V["x_lpa", "x_lpa"], f$V["x_ldl", "x_ldl"], f$V["x_lpa", "x_ldl"])
  E <- fr * s                                   # mmol/L of LDL-C with the genetic effect of 144 nmol/L of Lp(a)
  seD <- sqrt(f$V["x_ldl", "x_ldl"])
  out <- data.table(cohort = cohort, horizon = horizon, model = model, ldl_scale = scale, events = as.integer(f$events),
    ldl_mmol_per_sd = s,
    genetic_hr_per_mmol_ldl = exp(bD / s), ldl_lo = exp((bD - 1.96 * seD) / s), ldl_hi = exp((bD + 1.96 * seD) / s),
    genetic_to_trial_ratio = (bD / s) / abs(log(CTT)),
    genetic_hr_per_144_lower = exp(-bL),
    equivalent_ldl_mmol = E[["est"]], equiv_lo = E[["lo"]], equiv_hi = E[["hi"]],
    predicted_trial_hr = CTT^E[["est"]], pred_lo = CTT^E[["hi"]], pred_hi = CTT^E[["lo"]])
  out[, z_trial_vs_predicted := (log(TRIAL[["hr"]]) - log(predicted_trial_hr)) /
        sqrt(TRIAL_SE^2 + ((log(pred_hi) - log(pred_lo)) / (2 * 1.96))^2)]
  out[, p_trial_vs_predicted := 2 * pnorm(-abs(z_trial_vs_predicted))]
  out
}

run_posthoc <- function(P, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  P <- P[complete.cases(P[, c(EXPOSURES, ADJ), with = FALSE])]
  A <- add_groups(P[qual == TRUE]); Cc <- P[broad == FALSE]; w <- std_weights(Cc, A)

  # ---- P1a composition (no outcome) ----
  comp <- function(D, label) data.table(group = label, n = nrow(D), predicted_lpa_mean = mean(D$pred),
    predicted_lpa_ge150_pct = 100 * mean(D$pred >= 150), ldl_score_mean = mean(D$x_ldl), nonlipid_score_mean = mean(D$x_prs),
    statin_pct = 100 * mean(D$statin), age_mean = mean(D$age), male_pct = 100 * mean(D$sex == 1))
  subs <- list("prior infarction" = A[mi_pre == TRUE], "prior ischaemic stroke" = A[is_pre == TRUE],
               "peripheral artery disease" = A[pad_pre == TRUE], "infarction, no stroke" = A[mi_pre == TRUE & is_pre == FALSE],
               "stroke, no infarction" = A[is_pre == TRUE & mi_pre == FALSE])
  composition <- rbindlist(c(list(comp(Cc, "primary prevention"), comp(A, "established disease (all)")),
                             lapply(names(subs), function(nm) comp(subs[[nm]], nm)),
                             lapply(sort(unique(A$qgroup)), function(q) comp(A[qgroup == q], paste("infarction/stroke/peripheral =", q)))))

  # ---- P1b stratified model and P1c components within each condition ----
  est <- list()
  for (e in names(ENDPOINTS)) for (h in c("full", "5y")) {
    hz <- if (h == "5y") H5 else Inf
    est[[length(est) + 1]] <- est_rows(fit(A, e, hz), "A", e, h, "pre-registered (pooled)")
    est[[length(est) + 1]] <- est_rows(fit_strat(A, e, hz), "A", e, h, "stratified on qualifying condition")
    for (nm in names(subs)) est[[length(est) + 1]] <- est_rows(fit(subs[[nm]], e, hz), nm, e, h, "within condition")
  }
  est <- rbindlist(est)
  est[, `:=`(hr = exp(b), hr_lo = exp(b - 1.96 * se), hr_hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))]
  est[exposure == "x_lpa", `:=`(hr_per_144_lower = exp(-b), lower_lo = exp(-b - 1.96 * se), lower_hi = exp(-b + 1.96 * se))]

  contrast <- rbindlist(lapply(c("mace4", "mace3"), function(e) {
    fC <- fit(Cc, e, w = w); fS <- fit_strat(A, e); fP <- fit(A, e)
    one <- function(fA, label) {
      rs <- rbindlist(lapply(EXPOSURES, function(x) ratio_stats(fA, fC, x)))
      v <- verdict(rs, spec_test(fA, fC))
      data.table(endpoint = e, established_model = label, rs, posthoc_verdict = v$verdict, specificity_p = v$specificity_p)
    }
    rbind(one(fP, "pre-registered (pooled)"), one(fS, "stratified on qualifying condition"))
  }))

  # ---- P2 calibration against LDL-C ----
  slopes <- list(A = c(measured = ldl_slope(A), corrected = ldl_slope(A, TRUE)),
                 C = c(measured = ldl_slope(Cc, w = w), corrected = ldl_slope(Cc, TRUE, w = w)))
  cal <- list()
  for (h in c("full", "5y")) {
    hz <- if (h == "5y") H5 else Inf
    fits <- list(list(fit(A, "mace4", hz), "A", "pre-registered (pooled)", "A"),
                 list(fit_strat(A, "mace4", hz), "A", "stratified on qualifying condition", "A"),
                 list(fit(Cc, "mace4", hz, w = w), "C standardized", "pre-registered", "C"))
    for (ft in fits) for (sc in c("measured", "corrected"))
      cal[[length(cal) + 1]] <- calibrate(ft[[1]], slopes[[ft[[4]]]][[sc]], ft[[2]], h, ft[[3]],
                                          if (sc == "measured") "measured LDL-C" else "LDL-C / 0.7 in statin users")
  }
  cal <- rbindlist(cal)

  wcsv(composition, file.path(outdir, "composition.csv")); wcsv(est, file.path(outdir, "estimates_posthoc.csv"))
  wcsv(contrast, file.path(outdir, "stage_contrast_posthoc.csv")); wcsv(cal, file.path(outdir, "calibration.csv"))
  list(composition = composition, est = est, contrast = contrast, cal = cal)
}

selftest_posthoc <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { check(cond, paste("selftest:", msg)); n_ok <<- n_ok + 1L }

  # Fieller: agrees with the delta method when the denominator is precise, refuses when it is not
  f1 <- fieller(2, 10, 0.01, 0.0001, 0)
  expect(abs(f1[["est"]] - 0.2) < 1e-12 && abs((f1[["hi"]] - f1[["lo"]]) / 2 - 1.96 * sqrt(0.01 / 100 + 4 * 0.0001 / 1e4)) < 1e-4, "Fieller matches the delta method for a precise denominator")
  expect(all(is.na(fieller(1, 0.1, 0.01, 0.01, 0)[c("lo", "hi")])), "Fieller returns no interval when the denominator may be zero")
  set.seed(11); a <- rnorm(5000, 0.3, 0.05); b <- rnorm(5000, 0.1, 0.02)
  covers <- mapply(function(ai, bi) { z <- fieller(ai, bi, 0.05^2, 0.02^2, 0); is.na(z[["lo"]]) || (z[["lo"]] < 3 && 3 < z[["hi"]]) }, a, b)
  expect(abs(mean(covers) - 0.95) < 0.01, "Fieller interval has 95% coverage in simulation")

  # composition artefact: stroke patients have lower Lp(a) and a higher event rate, true Lp(a) effect zero
  set.seed(12); n <- 40000
  D <- data.table(pred = pmax(0, rgamma(n, shape = 1.2, scale = 14) + 110 * rbinom(n, 2, 0.08)), score_ldl_raw = rnorm(n), prs_n = rnorm(n),
                  score_lpa_raw = rnorm(n), age = round(runif(n, 40, 70)), sex = rbinom(n, 1, 0.5))
  for (k in 1:10) set(D, j = paste0("PC", k), value = rnorm(n))
  D <- finish_exposures(D, NULL)$P
  D[, is_pre := rbinom(n, 1, plogis(0.4 - 3 * x_lpa)) == 1]; D[, `:=`(mi_pre = !is_pre, pad_pre = FALSE)]
  D[, t0 := 14000L]; D[, t_end := t0 + 4000L]
  te <- D$t0 + as.integer(ceiling(365.25 * rexp(n, 0.03 * exp(1.0 * D$is_pre + 0.09 * D$x_ldl))))
  D[, `:=`(t_mi = fifelse(te <= t_end, te, NA_integer_), t_is = NA_integer_, t_cvd = NA_integer_, t_rev = NA_integer_)]
  D <- add_groups(D)
  pooled <- fit(D, "mace4"); strat <- fit_strat(D, "mace4")
  message(sprintf("[posthoc sim] Lp(a) log HR pooled %.3f (SE %.3f), stratified %.3f (SE %.3f); planted 0",
                  pooled$b[["x_lpa"]], sqrt(pooled$V[1, 1]), strat$b[["x_lpa"]], sqrt(strat$V[1, 1])))
  expect(pooled$b[["x_lpa"]] / sqrt(pooled$V[1, 1]) < -3, "pooling conditions manufactures an inverse Lp(a) association (the artefact exists)")
  expect(abs(strat$b[["x_lpa"]]) < 3 * sqrt(strat$V[1, 1]), "stratifying on the condition removes it")
  expect(abs(strat$b[["x_ldl"]] - 0.09) < 4 * sqrt(strat$V[2, 2]), "stratified model recovers the LDL-C effect")

  # calibration: planted LDL-C slope and effects give the planted equivalence
  s <- simulate(300000, 0.25, 0.09, seed = 301)
  P <- s$P; set.seed(13); P[, ldl := 3.5 + 0.30 * x_ldl + rnorm(.N, 0, 0.8)]
  Cc <- P[broad == FALSE]
  sl <- ldl_slope(Cc); expect(abs(sl - 0.30) < 0.01, "LDL-C slope per SD of score is recovered")
  expect(abs(ldl_slope(Cc, corrected = TRUE) - sl) > 0.005, "the statin-corrected scale differs from the measured one")
  cl <- calibrate(fit(Cc, "mace4"), sl, "C", "full", "sim", "measured LDL-C")
  truth <- 0.25 / 0.09 * 0.30
  expect(abs(cl$equivalent_ldl_mmol - truth) < (cl$equiv_hi - cl$equiv_lo) && cl$pred_lo < cl$predicted_trial_hr && cl$predicted_trial_hr < cl$pred_hi &&
           abs(cl$predicted_trial_hr - CTT^cl$equivalent_ldl_mmol) < 1e-12, "equivalence is close to the planted value and the predicted trial hazard ratio follows from it")
  expect(abs(cl$genetic_to_trial_ratio - (0.09 / 0.30) / abs(log(CTT))) < 0.25, "genetic-to-trial ratio is computed on the per-mmol/L scale")

  # smoke run of everything the real run executes
  out <- run_posthoc(P, file.path(workdir, "posthoc_smoke"))
  expect(all(file.exists(file.path(workdir, "posthoc_smoke", c("composition.csv", "estimates_posthoc.csv", "stage_contrast_posthoc.csv", "calibration.csv")))), "every output file is written")
  expect(nrow(out$cal) == 12 && nrow(out$contrast) == 12 && nrow(out$est[model == "stratified on qualifying condition" & exposure == "x_lpa"]) >= 10, "all cohorts, models and horizons are present")
  cat(sprintf("\nPOSTHOC SELFTEST PASSED: %d checks\n", n_ok))
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
if (opt$step == "selftest") { selftest_posthoc(opt$workdir); quit(save = "no", status = 0L) }
check(opt$step == "run", "unknown --step")
CENS <- as.integer(as.IDate(opt$censor))
R <- load_real(opt, CENS)
if (!gates(R, CENS, file.path(opt$outdir, "gates_recheck"))) { cat("\nGATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
cat("\nGATES PASSED\n")
out <- run_posthoc(R$P, opt$outdir)
options(width = 230)
cat("\n=== P1a composition (no outcome) ===\n"); print(out$composition, digits = 3)
cat("\n=== P1b/c Lp(a), hazard ratio per 144 nmol/L HIGHER, full follow-up ===\n")
print(dcast(out$est[exposure == "x_lpa" & horizon == "full", .(cohort, model, endpoint, v = sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi))],
            cohort + model ~ endpoint, value.var = "v"))
cat("\n=== same, five-year horizon ===\n")
print(dcast(out$est[exposure == "x_lpa" & horizon == "5y", .(cohort, model, endpoint, v = sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi))],
            cohort + model ~ endpoint, value.var = "v"))
cat("\n=== LDL-C score, per SD, full follow-up ===\n")
print(dcast(out$est[exposure == "x_ldl" & horizon == "full", .(cohort, model, endpoint, v = sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi))],
            cohort + model ~ endpoint, value.var = "v"))
cat("\n=== events per cell, full follow-up ===\n")
print(dcast(out$est[exposure == "x_lpa" & horizon == "full"], cohort + model ~ endpoint, value.var = "events"))
cat("\n=== stage contrast, post hoc ===\n")
print(out$contrast[, .(endpoint, established_model, exposure, b_established, se_established, b_primary, p_present, p_diff, ratio, ratio_lo, ratio_hi, posthoc_verdict, specificity_p)], digits = 3)
cat("\n=== P2 calibration against LDL-C (trial observed 1.02, 0.90-1.15) ===\n"); print(out$cal, digits = 3)
cat("\nPOSTHOC_DONE\n")
