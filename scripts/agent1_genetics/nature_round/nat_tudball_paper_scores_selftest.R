#!/usr/bin/env Rscript
# ============================================================================
# Selftest for nat_tudball_paper_scores.R. SYNTHETIC DATA ONLY, seed 42; no UKB file is touched.
#
# The reproduction check is only worth anything if it can fail, and if "reproduces" means the real code.
# So this test does not trust the copied logic: it RUNS the port's own 03/04 scripts and the paper's own
# statin_interaction_incident.R and statin_interaction_robustness.R on the synthetic files, and requires the
# main script to match their output. Then it breaks things on purpose:
#   world A  planted collider, no true A x G: mechanism must be seen; IPTW must be flagged UNINFORMATIVE.
#            Statin lowering is set so the raw log LDL0 SMD is near the port's -1.31, so its overlap problem is
#            real rather than borderline
#   world B  weak selection, TRUE A x score_LDL, none for score_LPA: balance passes and the effect is found
#   world C  degenerate exposure: NA with a note, no number, no crash
#   world D  everyone has all scores, and baseline LDL-C is measured WITHOUT noise. Two jobs: (1) the sample then
#            equals the port's, so the analysis path (not just the copied positive-control path) can be held to the
#            real scripts up to the rescaling of G; (2) with no noise there is nothing to regress to the mean,
#            so the collider verdict has to be able to come out negative
#   plus     published targets on synthetic data must stop the script; bypass on an unmarked cohort is refused;
#            one shifted target must stop the script
# Row order in the synthetic cohort is shuffled on purpose: a merge() that sorted by eid would change the
# lasso folds and break the IPTW reproduction.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--workdir", type = "character"),
  make_option("--main",    type = "character", help = "path to nat_tudball_paper_scores.R"),
  make_option("--repo_r",  type = "character", help = "scripts/agent1_genetics/R (real paper scripts)"),
  make_option("--port_scripts", type = "character", default = path.expand("/path/to/tudball-tacc-port/scripts")),
  make_option("--n", type = "integer", default = 40000L)
)))
for (a in c("workdir", "main", "repo_r")) if (is.null(opt[[a]])) stop("--", a, " is required", call. = FALSE)
dir.create(opt$workdir, recursive = TRUE, showWarnings = FALSE)
RS <- file.path(R.home("bin"), "Rscript")

nfail <- 0L; npass <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok)
  if (ok) npass <<- npass + 1L else nfail <<- nfail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
}
run <- function(script, args, log, env = character()) {
  st <- system2(RS, c(shQuote(script), args), stdout = log, stderr = log, env = env)
  st
}
rd <- function(d, f) fread(file.path(d, f), skip = 1L)

# ---- synthetic world ----------------------------------------------------------------------------------
gen_world <- function(dir, n, seed = 42L, sel = 1, conf = 1, drop = 0.5, delta_ldl = 0, delta_lpa = 0,
                      write_lpa = TRUE, no_statin0 = FALSE, no_prior = FALSE, e0_sd = .15, full_cover = FALSE) {
  set.seed(seed)
  raw <- file.path(dir, "raw"); dir.create(raw, recursive = TRUE, showWarnings = FALSE)
  eid <- sample(1000000L:2000000L, n)                     # unsorted on purpose
  male <- rbinom(n, 1, .46); age0 <- runif(n, 40, 70)
  wb <- rbinom(n, 1, .9)
  pcs <- matrix(rnorm(n * 10), n, 10, dimnames = list(NULL, paste0("pc", 1:10)))
  bmi <- rnorm(n, 27.4, 4.5); townsend <- rnorm(n, -1.3, 3)
  dm <- rbinom(n, 1, .05)
  diabetes_self <- ifelse(dm == 1, 1L, sample(c(-3L, -1L, 0L), n, TRUE, c(.002, .002, .996)))
  smoking <- sample(c(-3L, 0L, 1L, 2L), n, TRUE, c(.004, .55, .35, .096))
  alcohol <- sample(c(-3L, 1:6), n, TRUE, c(.003, .2, .23, .26, .11, .11, .087))
  degree <- rbinom(n, 1, .33)
  # genotypes: dose_EA drives truth; the raw file counts allele A, so a third of the variants need the flip
  mk_snps <- function(k, prefix) {
    f <- runif(k, .1, .5); b <- runif(k, .02, .08)
    ea <- ifelse(seq_len(k) %% 3 == 0, "C", "A")
    dose_ea <- vapply(f, function(p) rbinom(n, 2, p), numeric(n))
    counted <- dose_ea; counted[, ea != "A"] <- 2 - dose_ea[, ea != "A"]
    colnames(counted) <- sprintf("%s%d_A", prefix, seq_len(k))
    counted[sample(length(counted), round(.005 * length(counted)))] <- NA
    list(counted = counted, weights = data.table(SNP = sprintf("%s%d", prefix, seq_len(k)), EA = ea, BETA = b),
         true = as.numeric(scale(as.numeric(dose_ea %*% b))))
  }
  L <- mk_snps(30, "rsldl"); P <- mk_snps(6, "rslpa")
  zL <- L$true; zP <- P$true
  dose_keep <- rbinom(n, 1, .97) == 1                     # dosage files do not cover the whole cohort
  if (full_cover) dose_keep[] <- TRUE
  wr_dose <- function(x, nm) fwrite(cbind(data.table(FID = eid[dose_keep], IID = eid[dose_keep], PAT = 0L, MAT = 0L, SEX = male[dose_keep] + 1L, PHENOTYPE = -9L),
                                          as.data.table(x$counted[dose_keep, , drop = FALSE])),
                                    file.path(dir, nm), sep = "\t", na = "NA")
  wr_dose(L, "dose_LDL.raw"); wr_dose(P, "dose_LPA.raw")
  fwrite(L$weights, file.path(dir, "ldl_score.txt"), sep = " ", col.names = FALSE)
  fwrite(P$weights, file.path(dir, "lpa_score.txt"), sep = " ", col.names = FALSE)

  # PRS as imperfect measures of the same genetic components
  z <- function(x) as.numeric(scale(x))
  ldl_prs <- .75 * zL + sqrt(1 - .75^2) * rnorm(n)
  cad_prs <- .45 * zL + .35 * zP + sqrt(1 - .45^2 - .35^2) * rnorm(n)
  cvd_prs <- .40 * zL + .30 * zP + sqrt(1 - .40^2 - .30^2) * rnorm(n)
  neg_prs <- rnorm(n)
  prs_na <- rbinom(n, 1, .02) == 1
  if (full_cover) prs_na[] <- FALSE
  ldl_prs[prs_na] <- NA; cad_prs[prs_na] <- NA; cvd_prs[prs_na] <- NA; neg_prs[prs_na] <- NA
  zna <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

  t0 <- as.IDate("2006-03-01") + round(runif(n, 0, 4 * 365))
  prior <- if (no_prior) rep(0L, n) else rbinom(n, 1, plogis(-3.5 + .04 * (age0 - 55) + .6 * male + .4 * zL))
  # the collider: baseline LDL-C is measured with noise, and treatment responds to the noisy value
  u <- rnorm(n); true_log <- log(3.5) + .08 * zL + .012 * zP + .2 * u
  e0 <- rnorm(n, 0, e0_sd); e1 <- rnorm(n, 0, .15)
  e2 <- rnorm(n, 0, .15)   # the measurement a prescriber acts on between visits, not the study's baseline draw
  pre <- true_log + e0
  lp <- -2.0 + sel * z(pre) + conf * (.07 * (age0 - 55) + .5 * male + 1.2 * dm + 3 * prior + .04 * (bmi - 27))
  statin0 <- rbinom(n, 1, plogis(lp))
  if (no_statin0) statin0[] <- 0L
  ldl0_obs <- exp(pre - drop * statin0)
  # NO true A x G unless planted: statin effect is a constant on the log hazard
  hz <- exp(log(.009) + .22 * zL + .20 * zP + .20 * u + .06 * (age0 - 55) + .7 * male + .5 * dm + 1.0 * prior -
            .25 * statin0 + delta_ldl * statin0 * zL + delta_lpa * statin0 * zP)
  tinc <- t0 + round(rexp(n, hz) * 365.25)
  tcad <- tinc
  prev_fo <- (prior == 1 & runif(n) < .9) | (prior == 0 & runif(n) < .02 & !no_prior)
  tcad[prev_fo] <- t0[prev_fo] - round(runif(sum(prev_fo), 1, 8 * 365))
  cens <- as.IDate("2022-10-31")
  tdeath <- t0 + round(rexp(n, .006 * exp(.08 * (age0 - 55))) * 365.25)
  tdeath[tdeath > as.IDate("2022-12-19")] <- NA
  tcad[tcad > cens] <- NA
  # acute events are a labelled subset; people with prior disease can still have one
  acs <- tinc; acs[runif(n) > .45] <- NA; acs[!is.na(acs) & acs > cens] <- NA

  # repeat visit
  attend <- rbinom(n, 1, .6) == 1
  # Initiation follows true LDL-C through e2. If it followed the study's own baseline draw e0 instead, the unadjusted
  # change score would be positive too, because initiators would be selected on a high e0 that then regresses;
  # that is a different mechanism from conditioning on LDL0, and this world plants only the conditioning one.
  a1 <- ifelse(statin0 == 1, 1L, rbinom(n, 1, plogis(-2.6 + .6 * z(true_log + e2) + .05 * (age0 - 55) + .4 * dm)))
  ldl1 <- exp(true_log + e1 - .37 * a1)
  statin1 <- a1; statin1[!attend] <- NA; ldl1[!attend] <- NA
  ldl0 <- ldl0_obs; ldl0[rbinom(n, 1, .05) == 1] <- NA

  # covariates with some missingness, so the complete-case step does real work
  miss <- function(x, p) { x[rbinom(length(x), 1, p) == 1] <- NA; x }
  hdl <- miss(1.45 - .1 * male + rnorm(n, 0, .35) - .01 * (bmi - 27), .08)
  tg <- miss(exp(rnorm(n, .4, .5)), .01); crp <- miss(exp(rnorm(n, .3, 1)), .01)
  a1c <- miss(36 + 12 * dm + rnorm(n, 0, 4), .03)
  sbp <- miss(138 + .5 * (age0 - 55) + rnorm(n, 0, 18), .01); dbp <- miss(82 + rnorm(n, 0, 10), .01)
  bmi_o <- miss(bmi, .01); townsend_o <- miss(townsend, .002)

  end <- pmin(t0 + round(10 * 365.25), fifelse(is.na(tdeath), cens, tdeath), cens)
  CAD1 <- as.integer(!is.na(acs) & acs <= end)
  time_y <- as.numeric(fifelse(CAD1 == 1L, acs, end) - t0) / 365.25

  coh <- data.table(eid = eid, t0 = t0, sex = male, age0 = age0, birth_year = 1950L, birth_month = 6L, bmi = bmi_o, townsend = townsend_o,
                    white_british = wb, as.data.table(pcs), LDL0 = ldl0, HDL0 = hdl, TC0 = ldl0 + 1.5, TG0 = tg, CRP0 = crp, glucose0 = 5.1, HbA1c0 = a1c,
                    LDL1 = ldl1, t1 = fifelse(attend, t0 + 1500L, as.IDate(NA)), diabetes_self = diabetes_self, smoking = smoking, alcohol = alcohol,
                    sbp = sbp, dbp = dbp, degree = degree, statin0 = statin0, statin1 = statin1, chol_med0 = statin0, prior_CAD_self = prior,
                    cvd_prs = cvd_prs, cad_prs = cad_prs, ldl_prs = ldl_prs, neg_prs = neg_prs, geno_batch = 1L, geno_sex = male,
                    sex_aneuploidy = rbinom(n, 1, .001), het_outlier = rbinom(n, 1, .005),
                    cvd_prs_z = zna(cvd_prs), cad_prs_z = zna(cad_prs), ldl_prs_z = zna(ldl_prs), neg_prs_z = zna(neg_prs),
                    death_date = tdeath, prior_CAD_hes = prior, prior_CAD = prior, first_acs_date = acs, first_revasc_date = as.IDate(NA),
                    CAD1 = CAD1, cad_time_years = time_y, CAD_or_revasc1 = CAD1, SYNTHETIC_COHORT = 1L)
  fwrite(coh, file.path(dir, "cohort.csv.gz"))

  # per-field files in the UKB layout: eid,"<ID>-<inst>.<arr>"
  wr_fid <- function(fid, ...) {
    cols <- list(...); dt <- as.data.table(c(list(eid = eid), cols))
    setnames(dt, names(dt), c("eid", sprintf("\"%s-%s\"", fid, names(cols))))
    fwrite(dt, file.path(raw, sprintf("fid%s.csv", fid)), quote = FALSE, na = "")
  }
  wr_fid("53", `0.0` = t0)
  wr_fid("40000", `0.0` = tdeath)
  CADF <- c("131296","131298","131300","131302","131304","131306")
  which_f <- sample(CADF, n, TRUE)
  for (f in CADF) { v <- tcad; v[which_f != f] <- NA; wr_fid(f, `0.0` = v) }
  # paper's exposure: union of a sex-specific question, broader than the 20003 statin list and missing for a few
  chol <- as.integer(statin0 == 1 | rbinom(n, 1, .03) == 1)
  code <- ifelse(chol == 1L, 1L, sample(c(-7L, 2L), n, TRUE, c(.9, .1)))
  q_missing <- rbinom(n, 1, .01) == 1
  f6153 <- ifelse(male == 0 & !q_missing, code, NA_integer_); f6177 <- ifelse(male == 1 & !q_missing, code, NA_integer_)
  wr_fid("6153", `0.0` = f6153); wr_fid("6177", `0.0` = f6177)
  if (write_lpa) {
    t_lpa <- log(20) + 1.0 * zP + .6 * rnorm(n)
    v0 <- exp(t_lpa + rnorm(n, 0, .1)); v1 <- exp(t_lpa + rnorm(n, 0, .1) + .05 * a1)
    v1[!attend] <- NA
    r0 <- ifelse(v0 < 3.8, 4L, ifelse(v0 > 189, 5L, 1L)); r1 <- ifelse(v1 < 3.8, 4L, ifelse(v1 > 189, 5L, 1L))
    v0[r0 %in% 4:5] <- NA; v1[r1 %in% 4:5] <- NA
    wr_fid("30790", `0.0` = v0, `1.0` = v1); wr_fid("30796", `0.0` = r0, `1.0` = r1)
  }
  cv <- data.table(FID = eid, IID = eid, age = age0, sex = male, bmi = miss(bmi, .02), smoking_current = as.integer(smoking == 2L),
                   diabetes = dm, hypertension = rbinom(n, 1, .25), as.data.table(pcs))
  setnames(cv, paste0("pc", 1:10), paste0("PC", 1:10))
  fwrite(cv, file.path(dir, "covariates.txt"), sep = "\t")
  invisible(dir)
}

main_args <- function(dir, out, extra = character()) c(
  "--cohort", shQuote(file.path(dir, "cohort.csv.gz")),
  "--scores", shQuote(paste(file.path(dir, c("dose_LPA.raw", "dose_LDL.raw")), collapse = ",")),
  "--weights", shQuote(paste(file.path(dir, c("lpa_score.txt", "ldl_score.txt")), collapse = ",")),
  "--raw_dir", shQuote(file.path(dir, "raw")), "--covar", shQuote(file.path(dir, "covariates.txt")),
  "--outdir", shQuote(out), extra)

# ---- world A: planted collider, no interaction ------------------------------------------------------------
A <- file.path(opt$workdir, "worldA"); B <- file.path(opt$workdir, "worldB"); C <- file.path(opt$workdir, "worldC")
cat("== generating world A\n"); gen_world(A, opt$n)

cat("== running the REAL port and paper scripts on world A, to build reproduction targets\n")
ran <- file.path(A, "real"); dir.create(ran, showWarnings = FALSE)
coh <- file.path(A, "cohort.csv.gz")
for (g in c("ldl_prs_z", "cad_prs_z")) {
  st <- run(file.path(opt$port_scripts, "03_analysis_cad.R"), c(shQuote(coh), shQuote(file.path(ran, g))),
            file.path(ran, paste0("03_", g, ".log")), env = paste0("PRS_VAR=", g))
  check(sprintf("real port 03 ran on synthetic data (%s)", g), st == 0L)
}
st <- run(file.path(opt$port_scripts, "04_analysis_ldl_repeat_visit.R"), c(shQuote(coh), shQuote(file.path(ran, "rep"))),
          file.path(ran, "04.log"), env = "PRS_VAR=ldl_prs_z")
check("real port 04 ran on synthetic data", st == 0L)
pap <- c("--scores", shQuote(paste(file.path(A, c("dose_LPA.raw", "dose_LDL.raw")), collapse = ",")),
         "--weights", shQuote(paste(file.path(A, c("lpa_score.txt", "ldl_score.txt")), collapse = ",")),
         "--raw_dir", shQuote(file.path(A, "raw")), "--covar", shQuote(file.path(A, "covariates.txt")))
check("real statin_interaction_incident.R ran on synthetic data",
      run(file.path(opt$repo_r, "statin_interaction_incident.R"), c(pap, "--out", shQuote(file.path(ran, "inc.csv"))), file.path(ran, "inc.log")) == 0L)
check("real statin_interaction_robustness.R ran on synthetic data",
      run(file.path(opt$repo_r, "statin_interaction_robustness.R"), c(pap, "--out", shQuote(file.path(ran, "rob.csv"))), file.path(ran, "rob.log")) == 0L)

build_targets <- function() {
  lab <- c("Unadjusted" = "Unadjusted", "Adjusted for baseline LDL + covariates" = "Adjusted", "IPTW / MSM" = "IPTW")
  tg <- list()
  for (g in c("ldl_prs_z", "cad_prs_z")) {
    d <- file.path(ran, g)
    m <- fread(file.path(d, "cad_interaction_models.csv"))[term == "A_x_G"]
    tg[[length(tg) + 1L]] <- m[, .(check_id = sprintf("coronary|%s|logit|%s|A_x_G|all", g, lab[model]), target = estimate, kind = "coef")]
    x <- fread(file.path(d, "cad_cox_sensitivity.csv"))[term == "A_x_G"]
    tg[[length(tg) + 1L]] <- x[, .(check_id = sprintf("coronary|%s|cox|%s|A_x_G|all", g, fifelse(model == "Cox IPTW", "IPTW", "Unadjusted")), target = estimate, kind = "coef")]
    s <- fread(file.path(d, "cad_interaction_by_ldl_tertile.csv"))[term == "A_x_G"]
    tg[[length(tg) + 1L]] <- s[, .(check_id = sprintf("coronary|%s|logit|%s|A_x_G|%s", g, lab[model], stratum), target = estimate, kind = "coef")]
    b <- fread(file.path(d, "iptw_balance.csv"))
    tg[[length(tg) + 1L]] <- b[, .(check_id = sprintf("balance|%s|%s", g, variable), target = smd_weighted, kind = "coef")]
    ss <- fread(file.path(d, "sample_summary.csv")); st <- unique(s[, .(stratum, n, events)])
    tg[[length(tg) + 1L]] <- data.table(check_id = sprintf("count|%s|%s", g, c("N_precc", "statin_users_precc", "events_precc", "n_cc", "events_cc")),
                                        target = c(ss$N, ss$statin_users, ss$events, sum(st$n), sum(st$events)), kind = "count")
  }
  d <- file.path(ran, "rep")
  m <- fread(file.path(d, "ldl_interaction_models.csv"))[term == "A_x_G"]
  tg[[length(tg) + 1L]] <- m[, .(check_id = sprintf("repeat|ldl_prs_z|lm|%s|A_x_G|all", lab[model]), target = estimate, kind = "coef")]
  s <- fread(file.path(d, "ldl_interaction_by_ldl0_tertile.csv"))[term == "A_x_G"]
  tg[[length(tg) + 1L]] <- s[, .(check_id = sprintf("repeat|ldl_prs_z|lm|%s|A_x_G|%s", lab[model], stratum), target = estimate, kind = "coef")]
  x <- fread(file.path(d, "ldl_change_score_initiators.csv"))
  cl <- c("Absolute change in LDL-C" = "Absolute", "Unadjusted relative change in LDL-C" = "Relative", "Adjusted relative change in LDL-C" = "RelativeAdjusted")
  tg[[length(tg) + 1L]] <- x[, .(check_id = sprintf("change|ldl_prs_z|%s", cl[model]), target = estimate, kind = "coef")]
  ss <- fread(file.path(d, "sample_summary.csv"))
  tg[[length(tg) + 1L]] <- data.table(check_id = c("count|repeat|N", "count|repeat|initiators"), target = c(ss$N, ss$initiators), kind = "count")
  i <- fread(file.path(ran, "inc.csv"))
  gi <- function(p) log(i[term == p, HR])
  tg[[length(tg) + 1L]] <- data.table(check_id = c("paper|score_LDL|interaction", "paper|score_LPA|interaction", "paper|score_LDL|main", "paper|score_LPA|main"),
    target = c(gi("LDL-C score x baseline statin (ratio of HRs)"), gi("Lp(a) score x baseline statin (ratio of HRs)"),
               gi("LDL-C score -> incident CAD (HR/SD)"), gi("Lp(a) score -> incident CAD (HR/SD)")), kind = "coef")
  tg[[length(tg) + 1L]] <- data.table(check_id = c("count|paper|n", "count|paper|events", "count|paper|prevalent_excluded"),
                                      target = c(i$n[1], i$events[1], i$prevalent_excluded[1]), kind = "count")
  r <- fread(file.path(ran, "rob.csv"))[model == "primary (as reported)"]
  tg[[length(tg) + 1L]] <- data.table(check_id = c("paper_cc|score_LDL|interaction", "paper_cc|score_LPA|interaction"),
    target = log(c(r[term == "score_LDL:statin", estimate], r[term == "score_LPA:statin", estimate])), kind = "coef")
  tg[[length(tg) + 1L]] <- data.table(check_id = c("count|paper_cc|n", "count|paper_cc|events"), target = c(r$n[1], r$events[1]), kind = "count")
  rbindlist(tg)
}
TG <- build_targets()
tgf <- file.path(A, "targets.csv"); fwrite(TG, tgf)
check(sprintf("reproduction targets harvested from real scripts (%d)", nrow(TG)), nrow(TG) >= 50L && !anyNA(TG$target))

# ---- the hard stop must actually stop ---------------------------------------------------------------------
cat("== hard-stop tests\n")
o <- file.path(A, "out_published"); lg <- file.path(A, "published.log")
st <- run(opt$main, main_args(A, o, "--pc_only"), lg)
check("published targets on synthetic data: script stops with status 3", st == 3L)
check("  and says POSITIVE CONTROL FAILED", any(grepl("POSITIVE CONTROL FAILED", readLines(lg))))

shift <- copy(TG); shift[check_id == "coronary|ldl_prs_z|logit|Unadjusted|A_x_G|all", target := target + 0.01]
sf <- file.path(A, "targets_shift.csv"); fwrite(shift, sf)
st <- run(opt$main, main_args(A, file.path(A, "out_shift"), c("--pc_only", "--pc_targets_synthetic_only", shQuote(sf))), file.path(A, "shift.log"))
check("one target shifted by 0.01: script stops with status 3", st == 3L)

shift2 <- copy(TG); shift2[check_id == "count|paper|n", target := target + 1]
sf2 <- file.path(A, "targets_shift2.csv"); fwrite(shift2, sf2)
st <- run(opt$main, main_args(A, file.path(A, "out_shift2"), c("--pc_only", "--pc_targets_synthetic_only", shQuote(sf2))), file.path(A, "shift2.log"))
check("paper n off by one: script stops with status 3", st == 3L)

# bypass on a cohort without the marker
un <- file.path(opt$workdir, "unmarked"); dir.create(un, showWarnings = FALSE)
for (f in c("dose_LDL.raw", "dose_LPA.raw", "ldl_score.txt", "lpa_score.txt", "covariates.txt", "raw"))
  file.symlink(file.path(A, f), file.path(un, f))
cu <- fread(coh); cu[, SYNTHETIC_COHORT := NULL]; fwrite(cu, file.path(un, "cohort.csv.gz")); rm(cu)
lg <- file.path(un, "bypass.log")
st <- run(opt$main, main_args(un, file.path(un, "out"), "--bypass_positive_control_synthetic_only"), lg)
check("bypass on an unmarked cohort is REFUSED", st != 0L && any(grepl("REFUSED", readLines(lg))))
lg <- file.path(un, "tg.log")
st <- run(opt$main, main_args(un, file.path(un, "out"), c("--pc_targets_synthetic_only", shQuote(tgf))), lg)
check("replacement targets on an unmarked cohort are REFUSED", st != 0L && any(grepl("REFUSED", readLines(lg))))

# ---- world A, full run, real targets, no bypass -----------------------------------------------------------
cat("== world A full run (no bypass; must match the real scripts)\n")
oa <- file.path(A, "out"); lg <- file.path(A, "full.log")
t1 <- Sys.time(); st <- run(opt$main, main_args(A, oa, c("--pc_targets_synthetic_only", shQuote(tgf))), lg)
cat(sprintf("   run time %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
check("world A exit status 0, positive control passes against the real scripts", st == 0L)
if (st != 0L) { cat(tail(readLines(lg), 40), sep = "\n"); quit(save = "no", status = 1) }
pc <- rd(oa, "positive_control_report.csv")
check(sprintf("all %d reproduction rows pass; worst |diff| on coefficients %.2e", nrow(pc), pc[kind == "coef", max(abs(diff))]), all(pc$pass))
check("IPTW rows reproduce too (shuffled row order did not break the lasso folds)", pc[grepl("IPTW", check_id), all(pass)] && pc[grepl("IPTW", check_id), .N] >= 10L)

files <- list.files(oa, "\\.csv$", full.names = TRUE)
check("every CSV opens with the '#' design line naming prevalent-user, no primary care, not a target-trial emulation",
      all(vapply(files, function(f) { l <- readLines(f, 1L); grepl("^#", l) && grepl("Prevalent-user", l) && grepl("no primary-care", l) && grepl("NOT a target-trial", l) }, logical(1))))
tb <- lapply(files, fread, skip = 1L)
check("no individual identifier column in any export", !any(vapply(tb, function(d) any(c("eid", "IID", "FID") %in% names(d)), logical(1))))
ctcols <- c("n", "n_A1", "n_A0", "events", "events_A1", "events_A0", "initiators", "n_repeat_sample", "v0_too_low", "v0_too_high", "v1_too_low",
            "v1_too_high", "v1_too_low_initiators", "v1_too_high_initiators", "events_logistic", "events_logistic_A1", "events_cox", "events_cox_A1", "prior_CAD_retained")
check("no count below 5 in any export",
      !any(vapply(tb, function(d) { cc <- intersect(ctcols, names(d)); length(cc) && any(vapply(cc, function(x) any(d[[x]] < 5, na.rm = TRUE), logical(1))) }, logical(1))))

R <- rd(oa, "coronary_estimators.csv"); D <- rd(oa, "coronary_balance_verdicts.csv"); Bt <- rd(oa, "coronary_balance_tables.csv")
ST <- rd(oa, "coronary_by_ldl_tertile.csv"); LADD <- rd(oa, "ladder_paper_to_port.csv")
Vtxt <- readLines(file.path(oa, "verdicts.txt"))
check("4 specs x (6 scores + mutual) present in the coronary table", uniqueN(R[, .(spec, G)]) == 28L)
check("each run has 3 estimators in both logistic and Cox", R[, uniqueN(paste(family, estimator)), by = .(spec, G)][, all(V1 == 6L)])
check("no NA where the design is fine", R[spec == "C1_port_as_published", !anyNA(estimate)])
check("mutual model reports both interactions and their contrast",
      all(c("A_x_score_LDL", "A_x_score_LPA", "A_x_score_LDL - A_x_score_LPA") %in% R[grepl("^mutual", G), term]))
check("balance scored on every propensity column, not the port's five", Bt[, .N, by = .(spec, G)][, min(N)] > 20L)
check("the verdict is the rule: UNINFORMATIVE exactly when some weighted SMD >= 0.1, in every run",
      D[, all((verdict == "UNINFORMATIVE") == (max_abs_smd >= 0.1))] && Bt[, .(m = max(abs(smd_weighted))), by = .(spec, G)][D, on = .(spec, G)][, all(abs(m - max_abs_smd) < 1e-12)])
worst_raw <- Bt[spec == "C1_port_as_published" & G == "ldl_prs_z" & variable == "log_LDL0", smd_raw]
check(sprintf("world A calibrated to the port's overlap problem: raw log LDL0 SMD %.2f (port -1.31)", worst_raw), worst_raw < -1 && worst_raw > -1.7)
check("planted strong indication gives UNINFORMATIVE balance in C1", D[spec == "C1_port_as_published", all(verdict == "UNINFORMATIVE")])
check("IPTW rows carry the uninformative tag when balance fails, and no row is called readable",
      R[spec == "C1_port_as_published" & estimator == "IPTW_MSM", all(grepl("UNINFORMATIVE", read))])
check("verdict file states UNINFORMATIVE for the coronary arm", any(grepl("^\\[C1_port_as_published\\] UNINFORMATIVE", Vtxt)))
check("verdict file carries the design header and the positive-control status", grepl("NOT a target-trial", Vtxt[1]) && any(grepl("Positive controls: PASSED", Vtxt)))
check("tertile table: 3 strata x 4 estimators per single-score run", ST[, .N, by = .(spec, G)][, all(N == 12L)])
ng <- R[G == "neg_prs_z" & term == "A_x_neg_prs_z" & !is.na(estimate)]
check(sprintf("negative-control score: %d of %d rows exclude the null (expect about 5%%; <= 15%% required)", ng[lower > 0 | upper < 0, .N], nrow(ng)),
      ng[lower > 0 | upper < 0, .N] <= ceiling(.15 * nrow(ng)))
check("specs differ in what they should: C4 has more exposed than C3 (paper's question is broader)",
      { sm <- rd(oa, "coronary_samples.csv"); sm[spec == "C4_paper_design", n_A1] > sm[spec == "C3_paper_outcome", n_A1] })
check("C2, C3, C4 are the same people; C1 is larger because prior CAD is kept",
      { sm <- rd(oa, "coronary_samples.csv"); uniqueN(sm[spec != "C1_port_as_published", n]) == 1L && sm[spec == "C1_port_as_published", n] > sm[spec == "C2_prevalent_removed", n] })

check("ladder has the published-model step 0 and one row per step on the LDL track", identical(sort(LADD[track == "score_LDL", step]), c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L)))
check("ladder step 0 equals the real paper script's output",
      { i <- fread(file.path(ran, "inc.csv")); abs(LADD[track == "score_LDL" & step == 0, HR] - i[term == "LDL-C score x baseline statin (ratio of HRs)", HR]) < 1e-8 })
check("ladder has no NA on world A", !anyNA(LADD$estimate))

CH <- rd(oa, "repeat_ldl_change_scores.csv"); RI <- rd(oa, "repeat_ldl_interaction.csv"); EF <- rd(oa, "repeat_score_effect_on_ldl0.csv")
vl <- function(g) Vtxt[grepl(sprintf("^\\[%s\\] unadj", g), Vtxt)]
chg <- function(g, ch) CH[G == g & change == ch]
check(sprintf("collider planted: LDL0-adjusted change score positive for score_LDL (%.3f), unadjusted near 0 (%.3f)",
              chg("score_LDL", "Relative_log_adjusted_LDL0")$estimate, chg("score_LDL", "Relative_log")$estimate),
      chg("score_LDL", "Relative_log_adjusted_LDL0")$lower > 0 && chg("score_LDL", "Relative_log")$lower < 0 && chg("score_LDL", "Relative_log")$upper > 0)
# Regressing on a noisy baseline: slope on G = beta_G * var(e0) / (var(T | G) + var(e0)) = .08 * .0225 / (.04 + .0225).
check(sprintf("LDL0-adjusted change score has the size the planted noise implies: CI (%.3f, %.3f) covers %.4f",
              chg("score_LDL", "Relative_log_adjusted_LDL0")$lower, chg("score_LDL", "Relative_log_adjusted_LDL0")$upper, .08 * .0225 / .0625),
      chg("score_LDL", "Relative_log_adjusted_LDL0")$lower < .0288 && chg("score_LDL", "Relative_log_adjusted_LDL0")$upper > .0288)
check(sprintf("same for ldl_prs_z (%.3f vs %.3f)", chg("ldl_prs_z", "Relative_log_adjusted_LDL0")$estimate, chg("ldl_prs_z", "Relative_log")$estimate),
      chg("ldl_prs_z", "Relative_log_adjusted_LDL0")$lower > 0 && chg("ldl_prs_z", "Relative_log")$lower < 0 && chg("ldl_prs_z", "Relative_log")$upper > 0)
check("verdict: SUPPORTS_COLLIDER for score_LDL and ldl_prs_z", grepl("SUPPORTS_COLLIDER", vl("score_LDL")) && grepl("SUPPORTS_COLLIDER", vl("ldl_prs_z")))
check(sprintf("off-target contrast: score_LPA adjusted change (%.3f) smaller than score_LDL's (%.3f)",
              chg("score_LPA", "Relative_log_adjusted_LDL0")$estimate, chg("score_LDL", "Relative_log_adjusted_LDL0")$estimate),
      abs(chg("score_LPA", "Relative_log_adjusted_LDL0")$estimate) < abs(chg("score_LDL", "Relative_log_adjusted_LDL0")$estimate))
check(sprintf("score_LPA's LDL-C effect is stated (%.3f per SD, planted about 0.012), and smaller than score_LDL's (%.3f)",
              EF[G == "score_LPA" & effect == "log_LDL0_per_SD", estimate], EF[G == "score_LDL" & effect == "log_LDL0_per_SD", estimate]),
      abs(EF[G == "score_LPA" & effect == "log_LDL0_per_SD", estimate]) < abs(EF[G == "score_LDL" & effect == "log_LDL0_per_SD", estimate]))
check("verdict: neg_prs_z does not support the collider", !grepl("SUPPORTS_COLLIDER", vl("neg_prs_z")))
check("repeat arm has a balance table now, unlike the port", nrow(rd(oa, "repeat_balance_tables.csv")) > 0L)

LP <- rd(oa, "repeat_lpa_optional.csv"); LC <- rd(oa, "repeat_lpa_counts.csv")
check("Lp(a) arm: both handling rules, both scores, no NA", LP[, uniqueN(paste(handling, G))] == 4L && !anyNA(LP$estimate))
check("Lp(a) reportability counted explicitly, in both tails", LC[!is.na(v1_too_low), v1_too_low] > 0 && LC[!is.na(v1_too_high), v1_too_high] > 0)
check("bound keeps more people than drop", LC[handling == "bound", n] > LC[handling == "drop", n])

# ---- world B: selection weak enough for balance, TRUE interaction planted ------------------------------------
cat("== world B: weak selection, true A x score_LDL = -0.35, A x score_LPA = 0\n")
gen_world(B, opt$n, seed = 42L, sel = .1, conf = .15, drop = .02, delta_ldl = -.35, delta_lpa = 0)
ob <- file.path(B, "out"); lg <- file.path(B, "full.log")
st <- run(opt$main, main_args(B, ob, "--bypass_positive_control_synthetic_only"), lg)
check("world B exit status 0 with the synthetic-only bypass", st == 0L)
if (st != 0L) { cat(tail(readLines(lg), 40), sep = "\n"); quit(save = "no", status = 1) }
RB <- rd(ob, "coronary_estimators.csv"); DB <- rd(ob, "coronary_balance_verdicts.csv")
# Recovery is read in C3, whose exposure IS the variable carrying the planted effect. C4's 6153|6177 exposure adds
# people the effect was not planted in, so C4 should come out attenuated, and that is checked as well.
ip <- function(g, tm = paste0("A_x_", g), sp = "C3_paper_outcome") RB[spec == sp & G == g & family == "cox" & estimator == "IPTW_MSM" & term == tm]
check(sprintf("balance passes when selection is weak (max |SMD| %.3f)", DB[spec == "C4_paper_design", max(max_abs_smd)]), DB[spec == "C4_paper_design", all(verdict == "BALANCED")])
check(sprintf("planted LDL interaction recovered under IPTW: %.3f (%.3f, %.3f), truth -0.35", ip("score_LDL")$estimate, ip("score_LDL")$lower, ip("score_LDL")$upper),
      ip("score_LDL")$upper < 0 && ip("score_LDL")$lower < -.35 && ip("score_LDL")$upper > -.35)
check(sprintf("no planted Lp(a) interaction, none found: %.3f (%.3f, %.3f)", ip("score_LPA")$estimate, ip("score_LPA")$lower, ip("score_LPA")$upper),
      ip("score_LPA")$lower < 0 && ip("score_LPA")$upper > 0)
check(sprintf("exposure misclassification attenuates: C4 %.3f is closer to 0 than C3 %.3f", ip("score_LDL", sp = "C4_paper_design")$estimate, ip("score_LDL")$estimate),
      abs(ip("score_LDL", sp = "C4_paper_design")$estimate) < abs(ip("score_LDL")$estimate))
ctr <- RB[spec == "C4_paper_design" & grepl("^mutual", G) & family == "cox" & estimator == "IPTW_MSM" & term == "A_x_score_LDL - A_x_score_LPA"]
check(sprintf("asymmetry contrast excludes 0: %.3f (%.3f, %.3f)", ctr$estimate, ctr$lower, ctr$upper), ctr$upper < 0)
check("readable tag only when balanced", RB[spec == "C4_paper_design" & estimator == "IPTW_MSM", all(grepl("^readable", read))])
Vb <- readLines(file.path(ob, "verdicts.txt"))
check("verdict line reads BALANCED and names the interaction", any(grepl("^\\[C4_paper_design\\] BALANCED; IPTW interaction present for:.*score_LDL", Vb)))
check("bypass is recorded in the verdict file", any(grepl("BYPASSED", Vb)))

# ---- world C: degenerate exposure ------------------------------------------------------------------------------
cat("== world C: nobody on statin at baseline per 20003; nobody with prior CAD; no Lp(a) fields\n")
gen_world(C, 8000L, seed = 42L, no_statin0 = TRUE, no_prior = TRUE, write_lpa = FALSE)
cc <- fread(file.path(C, "cohort.csv.gz")); cc[, statin1 := fifelse(is.na(statin1), NA_integer_, 0L)]; fwrite(cc, file.path(C, "cohort.csv.gz")); rm(cc)
oc <- file.path(C, "out"); lg <- file.path(C, "full.log")
st <- run(opt$main, main_args(C, oc, "--bypass_positive_control_synthetic_only"), lg)
check("world C exit status 0: degenerate design does not crash", st == 0L)
if (st != 0L) { cat(tail(readLines(lg), 40), sep = "\n"); quit(save = "no", status = 1) }
RC <- rd(oc, "coronary_estimators.csv"); Vc <- readLines(file.path(oc, "verdicts.txt"))
check("no exposure variance under 20003: C1-C3 estimates are all NA, each with a note",
      RC[spec != "C4_paper_design", all(is.na(estimate)) && all(nzchar(note))])
check("the paper's own exposure still has variance: C4 is estimated", RC[spec == "C4_paper_design" & estimator == "Unadjusted", !anyNA(estimate)])
check("verdict says NOT_EVALUABLE, not a confident line", any(grepl("^\\[C1_port_as_published\\] NOT_EVALUABLE", Vc)))
check("no initiators: repeat rows NA with a note", { ri <- rd(oc, "repeat_ldl_interaction.csv"); all(is.na(ri$estimate)) && all(nzchar(ri$note)) })
check("repeat verdict NOT_EVALUABLE", any(grepl("^\\[score_LDL\\] unadj.*NOT_EVALUABLE", Vc)))
check("Lp(a) arm reports not_run with a note when fields are absent",
      { l <- rd(oc, "repeat_lpa_optional.csv"); all(l$handling == "not_run") && all(is.na(l$estimate)) && all(nzchar(l$note)) })

# ---- world D: analysis path against the real scripts; no baseline noise ---------------------------------------
cat("== world D: full score coverage, baseline LDL-C measured without noise\n")
DW <- file.path(opt$workdir, "worldD")
gen_world(DW, 24000L, seed = 42L, e0_sd = 0, full_cover = TRUE)
cd <- file.path(DW, "cohort.csv.gz"); rd03 <- file.path(DW, "real03"); rd04 <- file.path(DW, "real04")
check("real port 03 and 04 ran on world D",
      run(file.path(opt$port_scripts, "03_analysis_cad.R"), c(shQuote(cd), shQuote(rd03)), file.path(DW, "03.log"), env = "PRS_VAR=ldl_prs_z") == 0L &&
      run(file.path(opt$port_scripts, "04_analysis_ldl_repeat_visit.R"), c(shQuote(cd), shQuote(rd04)), file.path(DW, "04.log"), env = "PRS_VAR=ldl_prs_z") == 0L)
od <- file.path(DW, "out"); lg <- file.path(DW, "full.log")
st <- run(opt$main, main_args(DW, od, "--bypass_positive_control_synthetic_only"), lg)
check("world D exit status 0", st == 0L)
if (st != 0L) { cat(tail(readLines(lg), 40), sep = "\n"); quit(save = "no", status = 1) }
RD <- rd(od, "coronary_estimators.csv")[spec == "C1_port_as_published" & G == "ldl_prs_z" & term == "A_x_ldl_prs_z"]
pm <- fread(file.path(rd03, "cad_interaction_models.csv"))[term == "A_x_G"]
px <- fread(file.path(rd03, "cad_cox_sensitivity.csv"))[term == "A_x_G"]
ps <- unique(fread(file.path(rd03, "cad_interaction_by_ldl_tertile.csv"))[, .(stratum, n)])
check(sprintf("same sample as the real port (%d = %d)", RD$n[1], sum(ps$n)), RD$n[1] == sum(ps$n))
# The analysis path standardises G in its own sample, so its interaction is the port's times that sample's SD of G.
cw <- fread(cd)[white_british == 1 & (is.na(sex_aneuploidy) | sex_aneuploidy != 1) & (is.na(het_outlier) | het_outlier != 1) & !is.na(statin0) & !is.na(CAD1) & !is.na(LDL0) & LDL0 > .5 & LDL0 < 10]
cvs <- c("bmi", "townsend", "sbp", "dbp", "HbA1c0", "CRP0", "HDL0", "TG0", "smoking", "alcohol", "diabetes_self", "degree", "prior_CAD", paste0("pc", 1:10), "age0")
cw <- cw[complete.cases(cw[, ..cvs])]; sg <- sd(cw$ldl_prs_z)
br <- function(fam, e) RD[family == fam & estimator == e, estimate]
check(sprintf("analysis path = real port x SD(G), logistic unadjusted (%.6f vs %.6f)", br("logistic", "Unadjusted"), pm[model == "Unadjusted", estimate] * sg),
      abs(br("logistic", "Unadjusted") - pm[model == "Unadjusted", estimate] * sg) < 1e-6)
check("  logistic adjusted", abs(br("logistic", "Adjusted_LDL0_covariates") - pm[model == "Adjusted for baseline LDL + covariates", estimate] * sg) < 1e-6)
check("  Cox unadjusted", abs(br("cox", "Unadjusted") - px[model == "Cox unadjusted", estimate] * sg) < 1e-6)
# glmnet standardises internally, so rescaling G should leave the fit, and hence the weights, equal up to its tolerance.
check(sprintf("  logistic IPTW, lasso refit on rescaled G (|diff| %.1e)", abs(br("logistic", "IPTW_MSM") - pm[model == "IPTW / MSM", estimate] * sg)),
      abs(br("logistic", "IPTW_MSM") - pm[model == "IPTW / MSM", estimate] * sg) < 2e-3)
check(sprintf("  Cox IPTW, lasso refit on rescaled G (|diff| %.1e)", abs(br("cox", "IPTW_MSM") - px[model == "Cox IPTW", estimate] * sg)),
      abs(br("cox", "IPTW_MSM") - px[model == "Cox IPTW", estimate] * sg) < 2e-3)

rw <- fread(cd)[!is.na(LDL0) & !is.na(LDL1) & LDL0 > .5 & LDL0 < 10 & LDL1 > .5 & LDL1 < 10 & white_british == 1 & !is.na(statin0) & !is.na(statin1) & statin0 == 0]
rcv <- c("bmi", "townsend", "diabetes_self", "smoking", "prior_CAD", paste0("pc", 1:10), "age0")
rw <- rw[complete.cases(rw[, ..rcv])]; sr <- sd(rw$ldl_prs_z)
pc4 <- fread(file.path(rd04, "ldl_change_score_initiators.csv")); pi4 <- fread(file.path(rd04, "ldl_interaction_models.csv"))[term == "A_x_G"]
CD <- rd(od, "repeat_ldl_change_scores.csv")[G == "ldl_prs_z"]; ID <- rd(od, "repeat_ldl_interaction.csv")[G == "ldl_prs_z" & term == "A_x_G"]
check("repeat arm, analysis path = real port x SD(G): change scores",
      abs(CD[change == "Relative_log", estimate] - pc4[model == "Unadjusted relative change in LDL-C", estimate] * sr) < 1e-8 &&
      abs(CD[change == "Relative_log_adjusted_LDL0", estimate] - pc4[model == "Adjusted relative change in LDL-C", estimate] * sr) < 1e-8)
check("repeat arm, analysis path = real port x SD(G): full-sample unadjusted interaction",
      abs(ID[estimator == "Unadjusted", estimate] - pi4[model == "Unadjusted", estimate] * sr) < 1e-8)

Vd <- readLines(file.path(od, "verdicts.txt"))
vd <- function(g) Vd[grepl(sprintf("^\\[%s\\] unadj", g), Vd)]
adjD <- rd(od, "repeat_ldl_change_scores.csv")[G == "score_LDL" & change == "Relative_log_adjusted_LDL0"]
check(sprintf("no baseline noise: adjusted change score collapses to %.4f (%.4f, %.4f), planted exactly 0", adjD$estimate, adjD$lower, adjD$upper), adjD$lower < 0 && adjD$upper > 0)
check("no baseline noise: the verdict can say no, score_LDL is NOT SUPPORTS_COLLIDER", !grepl("SUPPORTS_COLLIDER", vd("score_LDL")) && !grepl("NOT_EVALUABLE", vd("score_LDL")))

cat(sprintf("\nSELFTEST: %d passed, %d failed\n", npass, nfail))
quit(save = "no", status = if (nfail) 1L else 0L)
