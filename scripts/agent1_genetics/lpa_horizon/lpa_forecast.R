#!/usr/bin/env Rscript
# ============================================================================
# Forecasts for Lp(a) outcome trials that have not reported
# (OCEAN(a)-Outcomes NCT05581303, ACCLAIM-Lp(a) NCT06292013), with
# Lp(a)HORIZON in retrospect as a check. Specification:
# FORECAST_PREREGISTRATION.md, committed before this is run on real data.
#
# lpa_horizon.R is sourced unchanged; the function definitions and named
# constants of the two post hoc scripts are evaluated from their parsed files
# (their main blocks are unguarded, so they are not sourced).
#
#   --step selftest   synthetic data only
#   --step run        gates (must pass), then the forecasts
# ============================================================================
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
source(file.path(here, "lpa_horizon.R"))
take <- function(file, consts) for (e in parse(file.path(here, file))) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      (as.character(e[[2]]) %in% consts || (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))))) eval(e, globalenv())
}
take("lpa_horizon_posthoc.R", c("CTT", "TRIAL_SE"))
take("lpa_horizon_posthoc2.R", c("STROKE_RISK", "CORONARY", "BANDS", "BAND_LABELS", "KNOTS"))
check(all(vapply(c("fit_tv", "ldl_slope_d", "add_groups", "CTT", "TRIAL_SE"), exists, logical(1))), "the post hoc scripts did not provide their functions")

KNOT <- 150; DIV <- 0.6; NDRAW <- 20000L
ENDPOINTS$cor3 <- c("t_chd", "t_mi", "t_rev"); ENDPOINTS$cor2 <- c("t_chd", "t_mi")
TRIAL3 <- c(hr = 1.01, lo = 0.88, hi = 1.16)     # Lp(a)HORIZON three-component composite (Cho et al. 2026, Table 2)
FC <- list(ocean = 12128L, ocean_cor3 = 3455L, acc_est = 20140L, acc_est_mace4 = 6006L, risk55 = 228783L, risk55_mace4 = 20071L)
yes <- function(x) x %in% 1

grid <- CJ(L0 = c(225, 250, 275, 300), lowering = c(0.80, 0.90, 0.97))
SCEN <- rbind(
  data.table(trial = "Lp(a)HORIZON", L0 = 229.4, L1 = 65.2, lowering = 1 - 65.2 / 229.4, central = TRUE),
  grid[, .(trial = "OCEAN(a)-Outcomes", L0, L1 = L0 * (1 - lowering), lowering, central = L0 == 275 & lowering == 0.90)],
  grid[, .(trial = "ACCLAIM-Lp(a)", L0, L1 = L0 * (1 - lowering), lowering, central = L0 == 250 & lowering == 0.90)])

TARGETS <- data.table(
  trial    = c(rep("Lp(a)HORIZON", 2), rep("OCEAN(a)-Outcomes", 3), rep("ACCLAIM-Lp(a)", 9)),
  stratum  = c(rep("established disease", 2), rep("prior infarction and/or angioplasty", 3), rep("established disease", 2),
               rep("at risk: age 55+, high-risk factors", 2), rep("at risk: age 55+, no vascular disease", 2),
               "at risk: high-risk factors or vascular disease without event", "at risk: all primary prevention", "at risk: all primary prevention"),
  cohort   = c("horizon", "horizon", "ocean", "ocean", "ocean", "acc_est", "acc_est", "acc_risk", "acc_risk", "acc_risk_55", "acc_risk_55",
               "acc_risk_recorded", "primary_prevention", "primary_prevention"),
  endpoint = c("mace4", "mace3", "cor3", "cor2", "mace3", "mace4", "mace3", "mace4", "mace3", "mace4", "mace3", "mace4", "mace4", "mace3"),
  role     = c("check", "check", "primary", "secondary", "secondary", "primary", "secondary", "primary", "secondary",
               rep("sensitivity", 5)))

cohort_flags <- function(P) {
  nrisk <- yes(P$diabetes) + yes(P$hypertension) + yes(P$smoking_current)
  evr <- P$qual | P$pci_pre | P$cabg_pre
  high <- !P$broad & P$age >= 55 & (yes(P$diabetes) | nrisk >= 2)
  list(horizon = P$qual, ocean = P$mi_pre | P$pci_pre, acc_est = evr, acc_risk = high, acc_risk_55 = !P$broad & P$age >= 55,
       acc_risk_recorded = high | (P$broad & !evr & P$age >= 55), primary_prevention = !P$broad)
}

# the parts of a fall from L0 to L1 (nmol/L), in units of DELTA: all of it, and the part above the knot
drops <- function(L0, L1) c(total = (L0 - L1) / DELTA, above = (max(L0, KNOT) - max(L1, KNOT)) / DELTA)
eta_of <- function(B, model, L0, L1) {
  d <- drops(L0, L1)
  if (model == "linear") B[, "x_lpa"] * d[["total"]] else B[, "x_lpa"] * d[["total"]] + B[, "x_above"] * d[["above"]]
}
rdraw <- function(f, n = NDRAW) {
  p <- length(f$b); Z <- matrix(rnorm(n * p), n, p)
  B <- sweep(Z %*% chol(f$V), 2, f$b, "+"); colnames(B) <- names(f$b); B
}
pt_of <- function(f) matrix(f$b, 1, dimnames = list(NULL, names(f$b)))
kappa <- function(B, s) log(1 / CTT) * s / B[, "x_ldl"]      # trial-scale log HR per unit of lifelong log HR, through LDL-C
qq <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))
se_from <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.96)

# hazard ratio for lowering Lp(a) from L0 to L1, both scales
forecast <- function(f, B, s, model, L0, L1) {
  e0 <- eta_of(pt_of(f), model, L0, L1); e <- eta_of(B, model, L0, L1)
  k0 <- kappa(pt_of(f), s); k <- kappa(B, s)
  life <- c(exp(-e0), qq(exp(-e))); sc <- c(exp(-e0 * k0), qq(exp(-e * k)))
  data.table(model = model, scale = c("lifelong", "scaled to trial duration"), L0 = L0, L1 = L1,
             hr = c(life[1], sc[1]), lo = c(life[2], sc[2]), hi = c(life[3], sc[3]),
             ldl_score_nonpositive_draws_pct = 100 * mean(B[, "x_ldl"] <= 0))
}
# log hazard ratio in the established stratum as a fraction of that in the at-risk stratum
stage_ratio <- function(fe, Be, se_, fr, Br, sr, model, L0, L1) {
  e0 <- c(eta_of(pt_of(fe), model, L0, L1), eta_of(pt_of(fr), model, L0, L1))
  ee <- eta_of(Be, model, L0, L1); er <- eta_of(Br, model, L0, L1)
  k0 <- c(kappa(pt_of(fe), se_), kappa(pt_of(fr), sr)); ke <- kappa(Be, se_); kr <- kappa(Br, sr)
  data.table(model = model, scale = c("lifelong", "scaled to trial duration"), L0 = L0, L1 = L1,
             ratio = c(e0[1] / e0[2], e0[1] * k0[1] / (e0[2] * k0[2])),
             lo = c(qq(ee / er)[1], qq(ee * ke / (er * kr))[1]), hi = c(qq(ee / er)[2], qq(ee * ke / (er * kr))[2]),
             at_risk_benefits_more_pct = 100 * c(mean(er > ee), mean(er * kr > ee * ke)))
}

run_forecast <- function(P, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(42)
  P <- P[complete.cases(P[, c(EXPOSURES, ADJ), with = FALSE])]
  COH <- cohort_flags(P)
  fits <- list(); models <- list(); fc <- list()
  for (i in seq_len(nrow(TARGETS))) {
    tg <- TARGETS[i]; key <- paste(tg$cohort, tg$endpoint)
    if (is.null(fits[[key]])) {
      D <- P[COH[[tg$cohort]]]; tv <- mk(D, tg$endpoint)
      lin <- fit_tv(D, tv); spl <- fit_tv(D, tv, extra = data.table(x_above = pmax(D$x_lpa - KNOT / DELTA, 0)))
      check(!is.null(lin) && !is.null(spl), paste("model could not be fitted:", key))
      fits[[key]] <- list(linear = lin, spline = spl, Bl = rdraw(lin), Bs = rdraw(spl), s = ldl_slope_d(D, DIV),
                          n = nrow(D), n_above_knot = sum(D$pred >= KNOT))
    }
    ft <- fits[[key]]
    models[[i]] <- data.table(tg, n = ft$n, events = as.integer(ft$linear$events), n_above_knot = ft$n_above_knot, ldl_mmol_per_sd = ft$s,
      b_lpa = ft$linear$b[["x_lpa"]], se_lpa = sqrt(ft$linear$V["x_lpa", "x_lpa"]), b_ldl = ft$linear$b[["x_ldl"]], se_ldl = sqrt(ft$linear$V["x_ldl", "x_ldl"]),
      b_below = ft$spline$b[["x_lpa"]], b_change_above = ft$spline$b[["x_above"]], se_change_above = sqrt(ft$spline$V["x_above", "x_above"]),
      p_change = 2 * pnorm(-abs(ft$spline$b[["x_above"]] / sqrt(ft$spline$V["x_above", "x_above"]))))
    sc <- SCEN[trial == tg$trial]
    fc[[i]] <- rbindlist(lapply(seq_len(nrow(sc)), function(j) cbind(tg, central = sc$central[j], lowering = sc$lowering[j], rbind(
      forecast(ft$linear, ft$Bl, ft$s, "linear", sc$L0[j], sc$L1[j]), forecast(ft$spline, ft$Bs, ft$s, "spline", sc$L0[j], sc$L1[j])))))
  }
  models <- rbindlist(models); fc <- rbindlist(fc)

  # ACCLAIM-Lp(a): the stage contrast, at the central scenario (the linear ratio does not depend on it)
  cs <- SCEN[trial == "ACCLAIM-Lp(a)" & central == TRUE]
  ratios <- rbindlist(lapply(c("acc_risk", "acc_risk_55", "acc_risk_recorded", "primary_prevention"), function(rk) rbindlist(lapply(c("mace4", "mace3"), function(e) {
    E <- fits[[paste("acc_est", e)]]; Rk <- fits[[paste(rk, e)]]
    if (is.null(E) || is.null(Rk)) return(NULL)
    an <- ratio_stats(E$linear, Rk$linear, "x_lpa")
    cbind(at_risk_cohort = rk, endpoint = e, rbind(
      stage_ratio(E$linear, E$Bl, E$s, Rk$linear, Rk$Bl, Rk$s, "linear", cs$L0, cs$L1),
      stage_ratio(E$spline, E$Bs, E$s, Rk$spline, Rk$Bs, Rk$s, "spline", cs$L0, cs$L1)),
      analytic_linear_ratio = an$ratio, analytic_lo = an$ratio_lo, analytic_hi = an$ratio_hi, p_stage_difference = an$p_diff)
  }))))

  # Lp(a)HORIZON in retrospect
  obs <- list(mace4 = TRIAL, mace3 = TRIAL3)
  chk <- fc[trial == "Lp(a)HORIZON"]
  chk[, `:=`(observed = vapply(endpoint, function(e) obs[[e]][["hr"]], 0), observed_lo = vapply(endpoint, function(e) obs[[e]][["lo"]], 0),
             observed_hi = vapply(endpoint, function(e) obs[[e]][["hi"]], 0))]
  chk[, z := (log(observed) - log(hr)) / sqrt(se_from(observed_lo, observed_hi)^2 + se_from(lo, hi)^2)]
  chk[, p_observed_vs_forecast := 2 * pnorm(-abs(z))]

  wcsv(models, file.path(outdir, "forecast_models.csv")); wcsv(fc, file.path(outdir, "forecasts.csv"))
  wcsv(ratios, file.path(outdir, "acclaim_stage_ratio.csv")); wcsv(chk, file.path(outdir, "horizon_check.csv"))
  list(models = models, fc = fc, ratios = ratios, chk = chk, fits = fits, COH = COH, P = P)
}

# prior coronary angioplasty and bypass grafting, and death from coronary heart disease
add_forecast_columns <- function(P, raw_dir, CENS) {
  arr0 <- function(fid) {
    f <- file.path(raw_dir, sprintf("fid%s.csv", fid)); h <- names(fread(f, nrows = 0))
    cols <- grep(sprintf("^%s-0\\.[0-9]+$", fid), h, value = TRUE)
    check(length(cols) > 0, sprintf("field %s has no instance-0 columns", fid))
    d <- fread(f, select = c(h[1], cols), colClasses = "character"); setnames(d, 1, "IID"); d <- d[IID %in% P$IID]
    rbindlist(lapply(cols, function(cn) { v <- d[[cn]]; k <- ok(v)
      data.table(IID = d$IID[k], idx = rep(as.integer(sub("^.*\\.", "", cn)), sum(k)), v = v[k]) }))
  }
  a <- arr0("41272"); setnames(a, "v", "code"); b <- arr0("41282"); setnames(b, "v", "date")
  ops <- merge(a, b, by = c("IID", "idx")); ops[, d := days(date)]
  ops <- merge(ops[!is.na(d)], P[, .(IID, t0)], by = "IID")[d <= t0]
  sr <- arr0("20004"); sr[, code := as.integer(num(v))]
  pci <- unique(c(ops[grepl("^(K49|K50|K75)", code), IID], sr[code == 1070L, IID]))
  cabg <- unique(c(ops[grepl("^K4[0-6]", code), IID], sr[code == 1095L, IID]))
  P[, `:=`(pci_pre = IID %in% pci, cabg_pre = IID %in% cabg)]
  P[, t_chd := fifelse(!is.na(tdeath) & ok(cause) & grepl("^I2[0-5]", cause) & tdeath > t0 & tdeath <= CENS, tdeath, NA_integer_)]
  P
}
forecast_gates <- function(P) {
  P <- P[complete.cases(P[, c(EXPOSURES, ADJ), with = FALSE])]; COH <- cohort_flags(P)
  ev <- function(k, e) sum(mk(P[COH[[k]]], e)$ev)
  g <- data.table(gate = c("OCEAN(a)-like cohort", "its coronary-composite events", "ACCLAIM established cohort", "its four-point events",
                           "age 55+ without vascular disease", "its four-point events"),
                  value = c(sum(COH$ocean), ev("ocean", "cor3"), sum(COH$acc_est), ev("acc_est", "mace4"), sum(COH$acc_risk_55), ev("acc_risk_55", "mace4")),
                  target = unlist(FC))
  g[, pass := value == target]; print(g); isTRUE(all(g$pass))
}

selftest_forecast <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { check(cond, paste("selftest:", msg)); n_ok <<- n_ok + 1L }
  near <- function(a, b, tol = 1e-9) isTRUE(all(abs(a - b) < tol))

  # the split of a fall around the knot
  expect(near(drops(229.4, 65.2), c(164.2, 79.4) / DELTA) && near(drops(120, 20), c(100, 0) / DELTA) && near(drops(300, 200), c(100, 100) / DELTA),
         "a fall is split correctly into its parts above and below the knot")

  # forecasts from known coefficients
  V0 <- diag(1e-12, 3); nm <- c("x_lpa", "x_ldl", "x_prs"); dimnames(V0) <- list(nm, nm)
  f <- list(b = setNames(c(0.2, 0.1, 0.3), nm), V = V0)
  set.seed(1); r <- forecast(f, rdraw(f), 0.30, "linear", 275, 27.5)
  expect(near(r$hr[1], exp(-0.2 * 247.5 / DELTA), 1e-9) && near(r$lo[1], r$hr[1], 1e-4) && near(r$hi[1], r$hr[1], 1e-4), "linear lifelong forecast is exp(-slope x fall)")
  expect(near(r$hr[2], CTT^((0.2 * 247.5 / DELTA) / 0.1 * 0.30), 1e-9), "scaled forecast is the statin-trial rate ratio raised to the LDL-C equivalent")
  nm4 <- c(nm, "x_above"); V4 <- diag(1e-12, 4); dimnames(V4) <- list(nm4, nm4)
  f4 <- list(b = setNames(c(0.2, 0.1, 0.3, -0.2), nm4), V = V4); B4 <- rdraw(f4)
  expect(near(forecast(f4, B4, 0.30, "spline", 300, 200)$hr[1], 1, 1e-9) &&
           near(forecast(f4, B4, 0.30, "spline", 229.4, 65.2)$hr[1], exp(-0.2 * 84.8 / DELTA), 1e-9), "a slope that is flat above the knot contributes nothing there")

  # simulated intervals agree with the analytic ones where those exist
  V1 <- diag(c(0.03, 0.015, 0.015)^2); dimnames(V1) <- list(nm, nm); f1 <- list(b = setNames(c(0.2, 0.1, 0.3), nm), V = V1)
  set.seed(2); B1 <- rdraw(f1); r1 <- forecast(f1, B1, 0.30, "linear", 250, 25)
  expect(abs(r1$lo[1] / exp(-(0.2 + 1.96 * 0.03) * 225 / DELTA) - 1) < 0.01 && abs(r1$hi[1] / exp(-(0.2 - 1.96 * 0.03) * 225 / DELTA) - 1) < 0.01,
         "the simulated interval matches the analytic one for the linear model")
  expect(near(cov(B1)["x_lpa", "x_ldl"], 0, 2e-5) && abs(sd(B1[, "x_lpa"]) / 0.03 - 1) < 0.02, "draws have the fitted covariance")
  set.seed(3); f2 <- list(b = setNames(c(0.4, 0.1, 0.3), nm), V = V1); B2 <- rdraw(f2)
  sr <- stage_ratio(f1, B1, 0.30, f2, B2, 0.30, "linear", 250, 25)
  expect(near(sr$ratio[1], 0.5) && sr$lo[1] < 0.5 && 0.5 < sr$hi[1] && sr$at_risk_benefits_more_pct[1] > 99, "the stage ratio and its interval are right for known slopes")
  expect(near(sr$ratio[2], 0.5, 1e-9), "with equal LDL-C scaling in both strata the scaled ratio equals the lifelong one")

  # planted stage-specific slopes in simulated cohorts, through the whole pipeline
  s <- simulate(300000, 0.125, 0.09, seed = 303); P <- s$P; set.seed(31)
  P[, `:=`(pci_pre = FALSE, cabg_pre = FALSE, t_chd = NA_integer_, ldl = (3.5 + 0.30 * x_ldl + rnorm(.N, 0, 0.8)) * fifelse(statin, DIV, 1))]
  out <- run_forecast(P, file.path(workdir, "forecast_smoke"))
  files <- c("forecast_models.csv", "forecasts.csv", "acclaim_stage_ratio.csv", "horizon_check.csv")
  expect(all(file.exists(file.path(workdir, "forecast_smoke", files))), "every output file is written")
  expect(nrow(out$models) == nrow(TARGETS) && nrow(out$fc) == 4 * (2 + 12 * 12) && nrow(out$chk) == 8, "every target, scenario, model and scale is present")
  m <- out$models[cohort == "horizon" & endpoint == "mace4"]; mp <- out$models[cohort == "primary_prevention" & endpoint == "mace4"]
  message(sprintf("[forecast sim] slope established %.3f (planted 0.125), primary prevention %.3f (planted 0.25); LDL-C per SD %.3f (planted 0.30)",
                  m$b_lpa, mp$b_lpa, m$ldl_mmol_per_sd))
  expect(abs(m$b_lpa - 0.125) < 3.5 * m$se_lpa && abs(mp$b_lpa - 0.25) < 3.5 * mp$se_lpa && abs(m$ldl_mmol_per_sd - 0.30) < 0.03, "planted slopes and the LDL-C scale are recovered")
  rr <- out$ratios[at_risk_cohort == "primary_prevention" & endpoint == "mace4" & model == "linear" & scale == "lifelong"]
  message(sprintf("[forecast sim] stage ratio %.2f (%.2f to %.2f), analytic %.2f (%.2f to %.2f); planted 0.50", rr$ratio, rr$lo, rr$hi, rr$analytic_linear_ratio, rr$analytic_lo, rr$analytic_hi))
  ase <- (rr$analytic_hi - rr$analytic_lo) / (2 * 1.96)
  expect(near(rr$ratio, rr$analytic_linear_ratio) && abs(rr$ratio - 0.5) < 3.5 * ase && abs((rr$hi - rr$lo) / (rr$analytic_hi - rr$analytic_lo) - 1) < 0.25,
         "the stage ratio recovers the planted value and its interval agrees with the delta method")
  h <- out$chk[endpoint == "mace4" & model == "linear" & scale == "lifelong"]
  expect(near(h$hr, exp(-m$b_lpa * 164.2 / DELTA)) && h$observed == TRIAL[["hr"]] && h$p_observed_vs_forecast > 0 && h$p_observed_vs_forecast <= 1,
         "the retrospective check uses the published levels and the observed result")
  g <- out$fc[trial == "OCEAN(a)-Outcomes" & endpoint == "cor3" & model == "linear" & scale == "lifelong"]
  bo <- out$models[cohort == "ocean" & endpoint == "cor3", b_lpa]
  expect(nrow(g) == 12 && sum(g$central) == 1 && all(diff(g[lowering == 0.90][order(L0), hr]) * sign(bo) < 0), "the scenario grid is complete and a larger fall gives a larger forecast effect")
  cat(sprintf("\nFORECAST SELFTEST PASSED: %d checks\n", n_ok))
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
if (opt$step == "selftest") { selftest_forecast(opt$workdir); quit(save = "no", status = 0L) }
check(opt$step == "run", "unknown --step")
CENS <- as.integer(as.IDate(opt$censor))
R <- load_real(opt, CENS)
if (!gates(R, CENS, file.path(opt$outdir, "gates_recheck"))) { cat("\nGATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
P <- add_forecast_columns(R$P, opt$raw_dir, CENS)
if (!forecast_gates(P)) { cat("\nFORECAST GATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
cat("\nGATES PASSED\n")
out <- run_forecast(P, opt$outdir)
options(width = 250)
f3 <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
cat("\n=== cohorts and fitted slopes (log hazard ratio per 144 nmol/L higher) ===\n")
print(out$models[, .(trial, stratum = substr(stratum, 1, 44), endpoint, role, n, events, n_above_knot, ldl_mmol_per_sd, b_lpa, se_lpa, b_below, b_change_above, se_change_above, p_change)], digits = 3)
cat("\n=== PRIMARY: ACCLAIM-Lp(a) stage contrast, log HR established / log HR at risk ===\n")
print(out$ratios[, .(at_risk_cohort, endpoint, model, scale, ratio = f3(ratio, lo, hi), at_risk_benefits_more_pct, analytic = f3(analytic_linear_ratio, analytic_lo, analytic_hi), p_stage_difference)], digits = 3)
cat("\n=== forecasts at the central scenarios (hazard ratio for the trial's lowering) ===\n")
print(out$fc[central == TRUE, .(trial, stratum = substr(stratum, 1, 44), endpoint, role, model, scale, L0, L1, forecast = f3(hr, lo, hi), ldl_score_nonpositive_draws_pct)], digits = 3)
cat("\n=== OCEAN(a)-Outcomes coronary composite and ACCLAIM-Lp(a) four-point, over the grid (linear, lifelong | spline, lifelong) ===\n")
gg <- out$fc[role == "primary" & scale == "lifelong", .(trial, stratum = substr(stratum, 1, 30), L0, lowering, model, v = f3(hr, lo, hi))]
print(dcast(gg, trial + stratum + L0 + lowering ~ model, value.var = "v"))
cat("\n=== CHECK: Lp(a)HORIZON in retrospect ===\n")
print(out$chk[, .(endpoint, model, scale, forecast = f3(hr, lo, hi), observed = f3(observed, observed_lo, observed_hi), p_observed_vs_forecast)], digits = 3)
cat("\nFORECAST_DONE\n")
