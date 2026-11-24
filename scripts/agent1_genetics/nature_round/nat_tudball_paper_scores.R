#!/usr/bin/env Rscript
# ============================================================================
# Nature round, analysis D: Tudball's design on the paper's scores.
# How results are read is fixed in PREREGISTERED_READOUT.md section D (do not edit it).
#
# ATTRIBUTION
#   Design and estimators: Tudball, Cirulli, Holmes, Davey Smith, Richards. medRxiv 2026,
#   doi 10.64898/2026.10.18.26365857; repository prs-drug-response, MIT licence. Adapted
#   (copied, not edited) from the UKB port /path/to/tudball-tacc-port/scripts/
#   03_analysis_cad.R and 04_analysis_ldl_repeat_visit.R. That port departs from Tudball's
#   landmark design because there are no primary-care records on Lonestar6.
#
# WHAT THIS IS NOT
#   Prevalent-user exposure, no primary-care records, not a target-trial emulation. The
#   coefficient on A is not a causal statin effect.
#
# WHY THERE ARE FOUR COHORT SPECIFICATIONS, NOT ONE
#   The port's published rows differ from the paper's analysis in four independent ways, and
#   each one could move the LDL score x statin estimate on its own:
#     (a) exposure: port = field 20003 statin codes; paper = fields 6153|6177 code 1
#     (b) outcome:  port = HES primary or death I21-I24 within 10 y; paper = first-occurrence
#                   I20-I25 over all follow-up
#     (c) sample:   port KEEPS prior CAD. RESULTS_SUMMARY.md section 3 row D says "no prior CAD",
#                   but 03_analysis_cad.R never excludes it: prior_CAD is a covariate, its
#                   coefficient (OR 3.08) is in cad_interaction_models.csv, and only the adjusted
#                   model conditions on it. The published HR 0.935 and 0.989 therefore include
#                   people with prevalent disease, which is the exact reverse-causation failure
#                   that statin_interaction_incident.R was written to remove.
#     (d) score:    port = UKB standard PRS; paper = its own weighted allelic scores
#   C1..C4 change one thing at a time, and ladder.csv walks from the published paper estimate to
#   the published port estimate one step at a time, so a similarity between 0.93 and 0.935 cannot
#   be mistaken for replication and a difference can be attributed.
#     C1 acute coronary | statin0 20003    | port sample, prior CAD kept  (the port's own rows)
#     C2 acute coronary | statin0 20003    | prevalent first-occurrence I20-I25 removed
#     C3 paper I20-I25  | statin0 20003    | same people as C2
#     C4 paper I20-I25  | paper 6153|6177  | same people as C2 (closest to the paper)
#
# POSITIVE CONTROLS (hard stop; any failure ends the script before any new estimate is written)
#   PC1: port's coronary and repeat-visit rows for ldl_prs_z and cad_prs_z, recomputed by a verbatim
#        copy of the port's logic, must equal its published CSVs to 3 decimals.
#   PC2: the paper's interaction, built by a verbatim copy of statin_interaction_incident.R's prep,
#        must equal Table S35 (n 379,896; 32,053 events) and the complete-covariate rows in Fig. 2d.
#   Seed: 20261024, not the project's 42, because the lasso folds come from the RNG and PC1 cannot
#   reproduce under any other seed. The seed is reset immediately before each cv.glmnet call, which
#   is what each of the port's separate per-PRS processes effectively did.
#
# WHAT A "BALANCED" VERDICT DOES NOT MEAN
#   Baseline LDL-C in a prevalent user is an on-treatment value, so the propensity model and the balance
#   table work on the lowered number. Weighted SMD < 0.1 therefore shows balance on measured, partly
#   post-treatment values, not exchangeability on the LDL-C that led to the prescription. The selftest's
#   world A has no true A x G anywhere, and still gives a negative A x score_LDL under weights that pass
#   the rule. So BALANCED is required before an IPTW row can be read; it does not by itself make the row
#   causal, and it cannot certify the absence of selection. Readout section B is the test for that.
#
# DEGENERATE DESIGNS emit NA with a note, never a number. Output is summary-level only and any cell
# below 5 is suppressed. Every CSV begins with one '#' design line: read with fread(f, skip = 1).
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(optparse); library(survival); library(glmnet); library(sandwich)
})

PORT_SEED <- 20261024L
SMD_LIMIT <- 0.1                       # ANALYSIS_PLAN.md section 6, restated in the readout
PLAN_MIN_INITIATORS <- 800L            # ANALYSIS_PLAN.md section 6: below this the LDL arm is uninformative
CELL_MIN <- 5L
TOL_COEF <- 5e-4                       # "to 3 decimals"
ALL_G <- c("score_LDL", "score_LPA", "ldl_prs_z", "cad_prs_z", "cvd_prs_z", "neg_prs_z")
REPEAT_G <- c("ldl_prs_z", "score_LDL", "score_LPA", "neg_prs_z")
HDR <- paste0("# Tudball et al. 2026 design adapted to the paper's scores (MIT-licensed code, adaptation disclosed). ",
              "Prevalent-user exposure; no primary-care records; NOT a target-trial emulation; ",
              "coefficient on A is not a causal statin effect. Cells < 5 suppressed.")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cohort",  type = "character", help = "port cohort_with_outcomes.csv.gz"),
  make_option("--scores",  type = "character", help = "dose_LPA.raw,dose_LDL.raw (plink --export A)"),
  make_option("--weights", type = "character", help = "lpa_score.txt,ldl_score.txt"),
  make_option("--raw_dir", type = "character", help = "UKB per-field directory with fid<ID>.csv"),
  make_option("--covar",   type = "character", help = "covariates.txt, used only for positive control 2"),
  make_option("--outdir",  type = "character"),
  make_option("--censor",  type = "character", default = "2022-10-31"),
  make_option("--nfolds",  type = "integer",   default = 10L),
  make_option("--min_initiators_na", type = "integer", default = 100L,
              help = "fewer initiators than this: repeat-visit rows are NA [default %default]"),
  make_option("--pc_only", action = "store_true", default = FALSE,
              help = "stop after the positive controls"),
  make_option("--bypass_positive_control_synthetic_only", action = "store_true", default = FALSE,
              help = "refused unless every cohort row carries SYNTHETIC_COHORT == 1"),
  make_option("--pc_targets_synthetic_only", type = "character", default = NULL,
              help = "replace the published targets; refused unless every row carries SYNTHETIC_COHORT == 1")
)))
for (a in c("cohort", "scores", "weights", "raw_dir", "covar", "outdir"))
  if (is.null(opt[[a]])) stop("missing required argument --", a)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

write_hdr <- function(dt, name) {
  # A cohort identifier in an export would break the summary-only rule, so refuse rather than strip.
  if (any(c("eid", "IID", "FID") %in% names(dt))) stop("refusing to export individual identifiers in ", name)
  f <- file.path(opt$outdir, name)
  writeLines(HDR, f); fwrite(dt, f, append = TRUE, col.names = TRUE)
}

# ---- raw UKB field reader; instance 0 is the same regex as the paper's rdcsv -----------------------
rdcsv <- function(fid, inst = 0L) {
  f <- file.path(opt$raw_dir, sprintf("fid%s.csv", fid))
  d <- fread(f, na.strings = c("", "NA"))
  col <- grep(sprintf("^\"?%s-%d\\.0\"?$", fid, inst), names(d), value = TRUE)[1]
  if (is.na(col)) stop(sprintf("field %s instance %d not in %s", fid, inst, f))
  setnames(d, c("eid", col), c("IID", "v"))
  d[, .(IID = as.character(IID), v)]
}
# Scores, explicit allele orientation: verbatim from statin_interaction_incident.R.
build_score <- function(raw_file, weight_file, label) {
  g <- fread(raw_file); w <- fread(weight_file, header = FALSE, col.names = c("SNP","EA","BETA"))
  dc <- setdiff(names(g), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
  map <- merge(data.table(col = dc, SNP = sub("_[^_]+$", "", dc), counted = sub("^.*_", "", dc)),
               w, by = "SNP")
  sc <- numeric(nrow(g)); flip <- 0L
  for (i in seq_len(nrow(map))) {
    d <- g[[map$col[i]]]
    if (toupper(map$counted[i]) != toupper(map$EA[i])) { d <- 2 - d; flip <- flip + 1L }
    d[is.na(d)] <- mean(d, na.rm = TRUE); sc <- sc + d * map$BETA[i]
  }
  message(sprintf("  [%s] %d variants (%d flipped)", label, nrow(map), flip))
  data.table(IID = as.character(g$IID), s = sc)
}

# ---- published targets (from the port's CSVs in results_from_tacc and from Table S35) --------------
# Stored as coefficients (log scale), so a Cox HR of 0.989 is compared as log(0.989).
PUBLISHED <- fread(text = "check_id,target,kind
coronary|ldl_prs_z|logit|Unadjusted|A_x_G|all,-0.0640217732831038,coef
coronary|ldl_prs_z|logit|Adjusted|A_x_G|all,-0.0186344799730344,coef
coronary|ldl_prs_z|logit|IPTW|A_x_G|all,-0.00896376160291473,coef
coronary|ldl_prs_z|cox|Unadjusted|A_x_G|all,-0.0675166497391548,coef
coronary|ldl_prs_z|cox|IPTW|A_x_G|all,-0.0111496282415541,coef
coronary|ldl_prs_z|logit|Unadjusted|A_x_G|T1,-0.0990950141723426,coef
coronary|ldl_prs_z|logit|IPTW|A_x_G|T1,-0.194766305589468,coef
coronary|ldl_prs_z|logit|Unadjusted|A_x_G|T2,-0.0381487537530905,coef
coronary|ldl_prs_z|logit|IPTW|A_x_G|T2,0.0027135410021015,coef
coronary|ldl_prs_z|logit|Unadjusted|A_x_G|T3,0.115984396898606,coef
coronary|ldl_prs_z|logit|IPTW|A_x_G|T3,0.177996311750585,coef
balance|ldl_prs_z|age0,0.26041673615795,coef
balance|ldl_prs_z|log_LDL0,-0.241502399706206,coef
balance|ldl_prs_z|G,0.0889078727845493,coef
balance|ldl_prs_z|bmi,0.197068178262105,coef
balance|ldl_prs_z|townsend,0.0873461984790714,coef
count|ldl_prs_z|N_precc,388122,count
count|ldl_prs_z|statin_users_precc,63588,count
count|ldl_prs_z|events_precc,7453,count
count|ldl_prs_z|n_cc,317098,count
count|ldl_prs_z|events_cc,6090,count
coronary|cad_prs_z|logit|Unadjusted|A_x_G|all,-0.0685170111832526,coef
coronary|cad_prs_z|logit|Adjusted|A_x_G|all,-0.102267884560991,coef
coronary|cad_prs_z|logit|IPTW|A_x_G|all,-0.0583724465800701,coef
coronary|cad_prs_z|cox|Unadjusted|A_x_G|all,-0.076227266243318,coef
coronary|cad_prs_z|cox|IPTW|A_x_G|all,-0.0617537380104055,coef
coronary|cad_prs_z|logit|Unadjusted|A_x_G|T1,-0.0487307402585414,coef
coronary|cad_prs_z|logit|IPTW|A_x_G|T1,-0.00580527759378214,coef
coronary|cad_prs_z|logit|Unadjusted|A_x_G|T2,0.0313409827737192,coef
coronary|cad_prs_z|logit|IPTW|A_x_G|T2,0.0225522471642595,coef
coronary|cad_prs_z|logit|Unadjusted|A_x_G|T3,-0.202040481455531,coef
coronary|cad_prs_z|logit|IPTW|A_x_G|T3,-0.198712700865019,coef
balance|cad_prs_z|age0,0.236105993906821,coef
balance|cad_prs_z|log_LDL0,-0.239104874323789,coef
balance|cad_prs_z|G,0.0791958430418706,coef
balance|cad_prs_z|bmi,0.178716541152494,coef
balance|cad_prs_z|townsend,0.0739954021035439,coef
count|cad_prs_z|N_precc,388122,count
count|cad_prs_z|statin_users_precc,63588,count
count|cad_prs_z|events_precc,7453,count
count|cad_prs_z|n_cc,317098,count
count|cad_prs_z|events_cc,6090,count
repeat|ldl_prs_z|lm|Unadjusted|A_x_G|all,-0.00734272597799308,coef
repeat|ldl_prs_z|lm|Adjusted|A_x_G|all,-0.00159845731826125,coef
repeat|ldl_prs_z|lm|IPTW|A_x_G|all,-0.00904275760272749,coef
repeat|ldl_prs_z|lm|Unadjusted|A_x_G|T1,0.0206522903960157,coef
repeat|ldl_prs_z|lm|IPTW|A_x_G|T1,0.022276327030468,coef
repeat|ldl_prs_z|lm|Unadjusted|A_x_G|T2,0.0282765184844502,coef
repeat|ldl_prs_z|lm|IPTW|A_x_G|T2,0.0265262201478709,coef
repeat|ldl_prs_z|lm|Unadjusted|A_x_G|T3,-0.000511626434232384,coef
repeat|ldl_prs_z|lm|IPTW|A_x_G|T3,-0.00585901757266657,coef
change|ldl_prs_z|Absolute,-0.0790127220432361,coef
change|ldl_prs_z|Relative,0.00460723791466372,coef
change|ldl_prs_z|RelativeAdjusted,0.0431033983725438,coef
count|repeat|N,12235,count
count|repeat|initiators,1250,count
count|paper|n,379896,count
count|paper|events,32053,count
count|paper|prevalent_excluded,21847,count
count|paper_cc|n,354185,count
count|paper_cc|events,29939,count
")
# Paper's rows from Table S35 and exports/agent1/revision/statin_robustness.csv, as published HRs.
PUBLISHED <- rbind(PUBLISHED, data.table(
  check_id = c("paper|score_LDL|interaction", "paper|score_LPA|interaction", "paper|score_LDL|main", "paper|score_LPA|main",
               "paper_cc|score_LDL|interaction", "paper_cc|score_LPA|interaction"),
  target = log(c(0.926553374411201, 0.991254541855647, 1.0919400155506, 1.09620739463031,
                 0.920167656327393, 0.987031077972047)),
  kind = "coef"))

# ---- verbatim port logic, kept apart from the analysis path below so that nothing added for the new
#      rows (common sample, re-standardisation, extra columns) can leak into the reproduction ----------
port_tidy <- function(fit, label) {
  V <- sandwich::vcovHC(fit, type = "HC0")
  est <- coef(fit); se <- sqrt(diag(V))
  data.table(model = label, term = names(est), estimate = est, se = se)
}
port_coronary_exact <- function(cohort, G_var) {
  d <- copy(cohort)
  d <- d[!is.na(t0) & !is.na(age0) & !is.na(sex)]
  if ("white_british" %in% names(d)) d <- d[white_british == 1]
  if ("sex_aneuploidy" %in% names(d)) d <- d[is.na(sex_aneuploidy) | sex_aneuploidy != 1]
  if ("het_outlier" %in% names(d)) d <- d[is.na(het_outlier) | het_outlier != 1]
  d <- d[!is.na(statin0) & !is.na(CAD1)]
  d <- d[!is.na(LDL0) & LDL0 > 0.5 & LDL0 < 10]
  d <- d[!is.na(get(G_var))]
  d[, A := as.integer(statin0 == 1)]
  d[, G := get(G_var)]
  d[, A_x_G := A * G]
  d[, log_LDL0 := log(LDL0)]
  d[, ldl_tertile := cut(LDL0, quantile(LDL0, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("T1", "T2", "T3"))]
  d[, sex := factor(sex)]
  for (v in c("smoking", "alcohol", "diabetes_self", "degree", "prior_CAD")) if (v %in% names(d)) d[, (v) := factor(get(v))]
  out <- list(
    data.table(check_id = sprintf("count|%s|N_precc", G_var), observed = nrow(d)),
    data.table(check_id = sprintf("count|%s|statin_users_precc", G_var), observed = sum(d$A)),
    data.table(check_id = sprintf("count|%s|events_precc", G_var), observed = sum(d$CAD1)))
  covars <- intersect(c("bmi", "townsend", "sbp", "dbp", "HbA1c0", "CRP0", "HDL0", "TG0",
                        "smoking", "alcohol", "diabetes_self", "degree", "prior_CAD",
                        paste0("pc", 1:10)), names(d))
  d <- d[complete.cases(d[, c(covars, "age0"), with = FALSE])]
  out[[4]] <- data.table(check_id = sprintf("count|%s|n_cc", G_var), observed = nrow(d))
  out[[5]] <- data.table(check_id = sprintf("count|%s|events_cc", G_var), observed = sum(d$CAD1))

  f1 <- glm(CAD1 ~ age0 + sex + A + G + A_x_G, family = binomial, data = d)
  rhs <- paste(c("age0", "sex", "log_LDL0", covars, "A + G + A_x_G"), collapse = " + ")
  f2 <- glm(as.formula(paste("CAD1 ~", rhs)), family = binomial, data = d)
  X <- model.matrix(as.formula(paste("~", paste(c("age0", "sex", "log_LDL0", covars, "G"), collapse = " + "))), d)[, -1]
  set.seed(PORT_SEED)
  cvfit <- cv.glmnet(X, d$A, family = "binomial", alpha = 1, nfolds = opt$nfolds)
  pA <- as.numeric(predict(cvfit, X, s = "lambda.min", type = "response"))
  pA_marg <- mean(d$A)
  d[, w := ifelse(A == 1, pA_marg / pA, (1 - pA_marg) / (1 - pA))]
  qs <- quantile(d$w, c(0.01, 0.99)); d[, w := pmin(pmax(w, qs[1]), qs[2])]
  f3 <- glm(CAD1 ~ age0 + sex + A + G + A_x_G, family = quasibinomial, data = d, weights = w)
  r <- rbindlist(list(port_tidy(f1, "Unadjusted"), port_tidy(f2, "Adjusted"), port_tidy(f3, "IPTW")))[term == "A_x_G"]
  out[[6]] <- r[, .(check_id = sprintf("coronary|%s|logit|%s|A_x_G|all", G_var, model), observed = estimate)]
  rm(f1, f2, f3, X, cvfit); gc(FALSE)

  cx  <- coxph(Surv(cad_time_years, CAD1) ~ age0 + sex + A + G + A_x_G, data = d)
  cxw <- coxph(Surv(cad_time_years, CAD1) ~ age0 + sex + A + G + A_x_G, data = d, weights = w, robust = TRUE)
  out[[7]] <- data.table(check_id = sprintf("coronary|%s|cox|%s|A_x_G|all", G_var, c("Unadjusted", "IPTW")),
                         observed = c(coef(cx)["A_x_G"], coef(cxw)["A_x_G"]))
  # Port's own unadjusted-vs-robust SE rule is kept here only; the analysis path uses robust everywhere.
  out[[8]] <- rbindlist(lapply(levels(d$ldl_tertile), function(tt) {
    s <- d[ldl_tertile == tt]
    a <- coef(glm(CAD1 ~ age0 + sex + A + G + A_x_G, family = binomial, data = s))["A_x_G"]
    b <- coef(glm(CAD1 ~ age0 + sex + A + G + A_x_G, family = quasibinomial, data = s, weights = w))["A_x_G"]
    data.table(check_id = sprintf("coronary|%s|logit|%s|A_x_G|%s", G_var, c("Unadjusted", "IPTW"), tt),
               observed = c(a, b))
  }))
  out[[9]] <- rbindlist(lapply(c("age0", "log_LDL0", "G", "bmi", "townsend"), function(v) if (v %in% names(d)) {
    m1 <- weighted.mean(d[A == 1][[v]], d[A == 1]$w); m0 <- weighted.mean(d[A == 0][[v]], d[A == 0]$w)
    s <- sqrt((var(d[A == 1][[v]]) + var(d[A == 0][[v]])) / 2)
    data.table(check_id = sprintf("balance|%s|%s", G_var, v), observed = (m1 - m0) / s)
  }))
  rbindlist(out)
}
port_repeat_exact <- function(cohort, G_var) {
  d <- copy(cohort)
  d <- d[!is.na(LDL0) & !is.na(LDL1) & LDL0 > 0.5 & LDL0 < 10 & LDL1 > 0.5 & LDL1 < 10]
  if ("white_british" %in% names(d)) d <- d[white_british == 1]
  d <- d[!is.na(statin0) & !is.na(statin1) & statin0 == 0]
  d <- d[!is.na(get(G_var))]
  d[, A := as.integer(statin1 == 1)][, G := get(G_var)][, A_x_G := A * G]
  d[, log_LDL0 := log(LDL0)][, log_LDL1 := log(LDL1)]
  d[, ldl_tertile := cut(LDL0, quantile(LDL0, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("T1", "T2", "T3"))]
  d[, sex := factor(sex)]
  covars <- intersect(c("bmi", "townsend", "diabetes_self", "smoking", "prior_CAD", paste0("pc", 1:10)), names(d))
  for (v in intersect(c("smoking", "diabetes_self", "prior_CAD"), covars)) d[, (v) := factor(get(v))]
  d <- d[complete.cases(d[, c(covars, "age0"), with = FALSE])]
  out <- list(data.table(check_id = c("count|repeat|N", "count|repeat|initiators"), observed = c(nrow(d), sum(d$A))))
  f1 <- lm(log_LDL1 ~ A + G + A_x_G, data = d)
  f2 <- lm(as.formula(paste("log_LDL1 ~", paste(c("age0", "sex", "log_LDL0", covars, "A + G + A_x_G"), collapse = " + "))), data = d)
  X <- model.matrix(as.formula(paste("~", paste(c("age0", "sex", "log_LDL0", covars, "G"), collapse = " + "))), d)[, -1]
  set.seed(PORT_SEED)
  cvfit <- cv.glmnet(X, d$A, family = "binomial", alpha = 1, nfolds = opt$nfolds)
  pA <- as.numeric(predict(cvfit, X, s = "lambda.min", type = "response")); pm <- mean(d$A)
  d[, w := ifelse(A == 1, pm / pA, (1 - pm) / (1 - pA))]; qs <- quantile(d$w, c(0.01, 0.99)); d[, w := pmin(pmax(w, qs[1]), qs[2])]
  f3 <- lm(log_LDL1 ~ age0 + sex + A + G + A_x_G, data = d, weights = w)
  r <- rbindlist(list(port_tidy(f1, "Unadjusted"), port_tidy(f2, "Adjusted"), port_tidy(f3, "IPTW")))[term == "A_x_G"]
  out[[2]] <- r[, .(check_id = sprintf("repeat|%s|lm|%s|A_x_G|all", G_var, model), observed = estimate)]
  out[[3]] <- rbindlist(lapply(levels(d$ldl_tertile), function(tt) {
    s <- d[ldl_tertile == tt]
    data.table(check_id = sprintf("repeat|%s|lm|%s|A_x_G|%s", G_var, c("Unadjusted", "IPTW"), tt),
               observed = c(coef(lm(log_LDL1 ~ A + G + A_x_G, data = s))["A_x_G"],
                            coef(lm(log_LDL1 ~ age0 + sex + A + G + A_x_G, data = s, weights = w))["A_x_G"]))
  }))
  ini <- d[A == 1]
  ini[, abs_change := LDL1 - LDL0][, rel_change := log_LDL1 - log_LDL0]
  out[[4]] <- data.table(check_id = sprintf("change|%s|%s", G_var, c("Absolute", "Relative", "RelativeAdjusted")),
                         observed = c(coef(lm(abs_change ~ G, data = ini))["G"], coef(lm(rel_change ~ G, data = ini))["G"],
                                      coef(lm(rel_change ~ G + log_LDL0, data = ini))["G"]))
  rbindlist(out)
}

# ---- paper's data, verbatim from statin_interaction_incident.R --------------------------------------
build_paper_tables <- function() {
  raws <- strsplit(opt$scores, ",")[[1]]; wts <- strsplit(opt$weights, ",")[[1]]
  sc <- NULL
  for (k in seq_along(raws)) {
    # basename, because the paper's grepl on the whole path would label both files score_LPA under any
    # directory whose name contains LPA. Same labels as the paper for the paths it actually ran on.
    nm <- if (grepl("LPA", basename(raws[k]))) "score_LPA" else "score_LDL"
    d <- build_score(raws[k], wts[k], nm); setnames(d, "s", nm)
    sc <- if (is.null(sc)) d else merge(sc, d, by = "IID")
  }
  if (!all(c("score_LDL", "score_LPA") %in% names(sc))) stop("--scores must name one LPA and one LDL dosage file")
  base <- rdcsv("53");    setnames(base, "v", "t0")
  dth  <- rdcsv("40000"); setnames(dth,  "v", "tdeath")
  CADF <- c("131296","131298","131300","131302","131304","131306")  # I20 I21 I22 I23 I24 I25
  ev <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
               lapply(CADF, function(f) { d <- rdcsv(f); setnames(d,"v",paste0("d",f)); d }))
  ev[, tcad := do.call(pmin, c(lapply(paste0("d",CADF), function(c) as.IDate(get(c))), na.rm = TRUE))]
  st <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
               list(rdcsv("6153")[, .(IID, f = v)], rdcsv("6177")[, .(IID, m = v)]))
  st[, statin := as.integer((!is.na(f) & f == 1) | (!is.na(m) & m == 1))]
  st[is.na(f) & is.na(m), statin := NA_integer_]
  list(sc = sc, base = base[, .(IID, t0)], ev = ev[, .(IID, tcad)], st = st[, .(IID, statin)],
       dth = dth[, .(IID, tdeath)])
}
paper_positive_control <- function(pt) {
  cv <- fread(opt$covar); cv[, IID := as.character(IID)]
  d <- Reduce(function(a,b) merge(a,b,by="IID"), list(pt$sc, pt$base, pt$ev, pt$st, cv))
  d <- merge(d, pt$dth, by = "IID", all.x = TRUE)
  d[, `:=`(t0 = as.IDate(t0), tdeath = as.IDate(tdeath))]
  cens <- as.IDate(opt$censor)
  d <- d[!is.na(t0) & !is.na(statin)]
  prev <- d[!is.na(tcad) & tcad <= t0, .N]
  d <- d[is.na(tcad) | tcad > t0]
  d[, tend  := pmin(fifelse(is.na(tcad), cens, tcad), fifelse(is.na(tdeath), cens, tdeath), cens, na.rm = TRUE)]
  d[, event := as.integer(!is.na(tcad) & tcad <= tend)]
  d[, fu    := as.numeric(tend - t0) / 365.25]
  d <- d[fu > 0]
  pcs <- grep("^PC", names(d), value = TRUE)[1:10]
  PCS <- paste(pcs, collapse = " + ")
  fitone <- function(dd, S) {
    dd <- copy(dd); dd[, (S) := scale(get(S))[, 1]]
    m <- coxph(as.formula(sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", S, PCS)), data = dd)
    list(coef = coef(m), V = vcov(m), n = m$n, e = m$nevent)
  }
  LP <- fitone(d, "score_LPA"); LD <- fitone(d, "score_LDL")
  IND <- c("bmi","smoking_current","diabetes","hypertension")
  dc <- d[complete.cases(d[, ..IND])]
  LPc <- fitone(dc, "score_LPA"); LDc <- fitone(dc, "score_LDL")
  obs <- rbind(
    data.table(check_id = c("paper|score_LDL|interaction", "paper|score_LPA|interaction",
                            "paper|score_LDL|main", "paper|score_LPA|main"),
               observed = c(LD$coef["score_LDL:statin"], LP$coef["score_LPA:statin"],
                            LD$coef["score_LDL"], LP$coef["score_LPA"])),
    data.table(check_id = c("count|paper|n", "count|paper|events", "count|paper|prevalent_excluded"),
               observed = c(LD$n, LD$e, prev)),
    data.table(check_id = c("paper_cc|score_LDL|interaction", "paper_cc|score_LPA|interaction"),
               observed = c(LDc$coef["score_LDL:statin"], LPc$coef["score_LPA:statin"])),
    data.table(check_id = c("count|paper_cc|n", "count|paper_cc|events"), observed = c(LDc$n, LDc$e)))
  L0 <- rbind(
    data.table(G = "score_LDL", estimate = LD$coef["score_LDL:statin"], se = sqrt(LD$V["score_LDL:statin", "score_LDL:statin"])),
    data.table(G = "score_LPA", estimate = LP$coef["score_LPA:statin"], se = sqrt(LP$V["score_LPA:statin", "score_LPA:statin"])))
  L0[, `:=`(n = LD$n, events = LD$e)]
  list(obs = obs, ladder0 = L0)
}

# ---- load cohort, decide the positive-control mode ----------------------------------------------------
stamp("reading cohort")
cohort <- fread(opt$cohort)
is_synth <- "SYNTHETIC_COHORT" %in% names(cohort) && nrow(cohort) > 0 && all(cohort$SYNTHETIC_COHORT %in% 1L)
if ((opt$bypass_positive_control_synthetic_only || !is.null(opt$pc_targets_synthetic_only)) && !is_synth)
  stop("REFUSED: positive-control bypass or replacement targets requested on a cohort that is not marked synthetic. ",
       "These flags exist only for the selftest.")

stamp("building the paper's scores and fields")
pt <- build_paper_tables()

# ---- POSITIVE CONTROLS: hard stop ---------------------------------------------------------------------
if (opt$bypass_positive_control_synthetic_only) {
  message("!! POSITIVE CONTROL BYPASSED: synthetic cohort, selftest only.")
  pc_status <- "BYPASSED (synthetic cohort)"
  pc2 <- NULL
} else {
  targets <- if (is.null(opt$pc_targets_synthetic_only)) PUBLISHED else fread(opt$pc_targets_synthetic_only)
  stamp("positive control 1: port rows (ldl_prs_z, cad_prs_z, repeat visit)")
  obs <- rbindlist(list(port_coronary_exact(cohort, "ldl_prs_z"), port_coronary_exact(cohort, "cad_prs_z"),
                        port_repeat_exact(cohort, "ldl_prs_z")))
  stamp("positive control 2: paper's published interaction")
  pc2 <- paper_positive_control(pt)
  obs <- rbind(obs, pc2$obs)
  pcr <- merge(targets, obs, by = "check_id", all.x = TRUE)
  pcr[, diff := observed - target]
  pcr[, pass := !is.na(observed) & ifelse(kind == "count", observed == target, abs(diff) < TOL_COEF)]
  write_hdr(pcr, "positive_control_report.csv")
  if (!all(pcr$pass)) {
    bad <- pcr[pass == FALSE]
    message("POSITIVE CONTROL FAILED: ", nrow(bad), " of ", nrow(pcr), " checks. No estimate from this run may be read.")
    print(bad[, .(check_id, target, observed, diff)], nrows = 200)
    message("Deterministic rows (no lasso) failing means the cohort or score inputs changed; only IPTW/balance rows ",
            "failing points to the glmnet or BLAS build rather than the data. See positive_control_report.csv.")
    quit(save = "no", status = 3)
  }
  pc_status <- sprintf("PASSED: %d of %d checks within %g (coefficients) or exact (counts)", sum(pcr$pass), nrow(pcr), TOL_COEF)
  message("POSITIVE CONTROL ", pc_status)
}
if (opt$pc_only) { stamp("--pc_only: stopping after positive controls"); quit(save = "no", status = 0) }

# ---- attach the paper's fields to the port cohort WITHOUT reordering rows ------------------------------
# match(), not merge(): merge sorts by key, and a different row order is a different set of lasso folds.
iid <- as.character(cohort$eid)
idx_sc <- match(iid, pt$sc$IID)
cohort[, `:=`(score_LDL = pt$sc$score_LDL[idx_sc], score_LPA = pt$sc$score_LPA[idx_sc])]
cohort[, t0p := as.IDate(pt$base$t0)[match(iid, pt$base$IID)]]
cohort[, tcad := pt$ev$tcad[match(iid, pt$ev$IID)]]
cohort[, tdeath := as.IDate(pt$dth$tdeath)[match(iid, pt$dth$IID)]]
cohort[, statin_paper := pt$st$statin[match(iid, pt$st$IID)]]
rm(pt, idx_sc, iid); gc(FALSE)

# t0 appears in both sources and the port's cad_time_years is built from its own t0. If they disagree the
# C2-vs-C3 comparison would be comparing different time origins.
t0chk <- cohort[!is.na(t0) & !is.na(t0p), mean(as.IDate(t0) != t0p)]
if (is.na(t0chk) || t0chk > 1e-3) stop(sprintf("port t0 and field 53 disagree in a fraction %s of participants", format(t0chk)))

# ---- analysis base sample: the port's restrictions, then everything complete on all six scores ----------
base <- cohort[!is.na(t0) & !is.na(age0) & !is.na(sex)]
if ("white_british" %in% names(base)) base <- base[white_british == 1]
if ("sex_aneuploidy" %in% names(base)) base <- base[is.na(sex_aneuploidy) | sex_aneuploidy != 1]
if ("het_outlier" %in% names(base)) base <- base[is.na(het_outlier) | het_outlier != 1]
base <- base[!is.na(statin0) & !is.na(CAD1)]
base <- base[!is.na(LDL0) & LDL0 > 0.5 & LDL0 < 10]
n_port <- nrow(base)
cover <- mean(!is.na(base$score_LDL) & !is.na(base$score_LPA))
message(sprintf("paper's scores cover %.1f%% of the port sample", 100 * cover))
if (is.na(cover) || cover < 0.5) stop("paper's scores cover under half the port sample: eid/IID mismatch or wrong dosage file")
for (g in ALL_G) if (!g %in% names(base)) stop("score column missing: ", g)
base <- base[complete.cases(base[, ..ALL_G])]

base[, log_LDL0 := log(LDL0)]
base[, sex := factor(sex)]
for (v in c("smoking", "alcohol", "diabetes_self", "degree", "prior_CAD")) if (v %in% names(base)) base[, (v) := factor(get(v))]
COVARS <- intersect(c("bmi", "townsend", "sbp", "dbp", "HbA1c0", "CRP0", "HDL0", "TG0",
                      "smoking", "alcohol", "diabetes_self", "degree", "prior_CAD",
                      paste0("pc", 1:10)), names(base))
base <- base[complete.cases(base[, c(COVARS, "age0"), with = FALSE])]

# Paper's outcome, built as in statin_interaction_incident.R, on the port sample. event10/fu10 exist because
# the port's outcome stops at 10 years and the paper's does not; ladder.csv separates that from the code list.
cens <- as.IDate(opt$censor)
base[, tend := pmin(fifelse(is.na(tcad), cens, tcad), fifelse(is.na(tdeath), cens, tdeath), cens, na.rm = TRUE)]
base[, event := as.integer(!is.na(tcad) & tcad <= tend)]
base[, fu := as.numeric(tend - t0p) / 365.25]
base[, event10 := as.integer(event == 1 & fu <= 10)]
base[, fu10 := pmin(fu, 10)]
base[, incident_ok := !is.na(t0p) & (is.na(tcad) | tcad > t0p) & !is.na(fu) & fu > 0 & !is.na(statin_paper)]
n_prior_kept <- base[incident_ok == TRUE & prior_CAD == 1, .N]
stamp(sprintf("port sample %d -> complete on all scores and covariates %d -> incident-eligible %d (of which %s still flagged prior_CAD by HES/self-report)",
              n_port, nrow(base), sum(base$incident_ok), if (n_prior_kept < CELL_MIN) "<5" else n_prior_kept))

SPECS <- list(
  list(id = "C1_port_as_published", outcome = "acute_coronary_I21_I24_10y", exposure = "statin0_field20003", sample = "port_prior_CAD_kept"),
  list(id = "C2_prevalent_removed", outcome = "acute_coronary_I21_I24_10y", exposure = "statin0_field20003", sample = "incident_eligible"),
  list(id = "C3_paper_outcome",     outcome = "paper_firstocc_I20_I25",     exposure = "statin0_field20003", sample = "incident_eligible"),
  list(id = "C4_paper_design",      outcome = "paper_firstocc_I20_I25",     exposure = "statin_6153_6177",   sample = "incident_eligible"))

spec_frame <- function(sp) {
  s <- if (sp$sample == "port_prior_CAD_kept") copy(base) else copy(base[incident_ok == TRUE])
  s[, A := as.integer(if (sp$exposure == "statin0_field20003") statin0 == 1 else statin_paper == 1)]
  if (sp$outcome == "acute_coronary_I21_I24_10y") {
    s[, `:=`(Ybin = as.integer(CAD1), Ytime = cad_time_years, Yevent = as.integer(CAD1))]
  } else {
    # Logistic reads the 10-year event, as the port's CAD1 does; Cox reads all follow-up, as the paper does.
    s[, `:=`(Ybin = event10, Ytime = fu, Yevent = event)]
  }
  s <- droplevels(s)
  for (g in ALL_G) s[, (g) := as.numeric(scale(get(g)))]     # standardised in THIS analysis sample
  for (g in ALL_G) s[, (paste0("A_x_", g)) := A * get(g)]
  s[, ldl_tertile := tryCatch(cut(LDL0, quantile(LDL0, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("T1", "T2", "T3")),
                              error = function(e) factor(rep(NA_character_, .N), levels = c("T1", "T2", "T3")))]
  s
}
active_covars <- function(s, covars) {
  # A factor with one level (prior_CAD after an exclusion, say) would make glm error on contrasts.
  keep <- vapply(covars, function(v) uniqueN(s[[v]]) > 1L, logical(1))
  if (any(!keep)) message("  dropped single-valued covariates: ", paste(covars[!keep], collapse = ", "))
  covars[keep]
}

# ---- estimation helpers. Every failure path returns NA rows with a note, never a number ----------------
safely <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) { w <<- c(w, paste0("ERROR: ", conditionMessage(e))); NULL }),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  list(val = v, warn = w)
}
na_rows <- function(terms, contrast, note) {
  tt <- c(terms, if (!is.null(contrast)) paste(contrast, collapse = " - "))
  data.table(term = tt, estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_, note = note)
}
coef_rows <- function(est, V, terms, contrast = NULL, note = "") {
  if (!all(terms %in% names(est)) || anyNA(est[terms])) return(na_rows(terms, contrast, "term aliased or absent"))
  se <- sqrt(diag(V))[terms]
  if (any(!is.finite(se))) return(na_rows(terms, contrast, "non-finite SE"))
  r <- data.table(term = terms, estimate = unname(est[terms]), se = unname(se))
  if (!is.null(contrast)) {
    a <- contrast[1]; b <- contrast[2]
    r <- rbind(r, data.table(term = paste(contrast, collapse = " - "), estimate = est[[a]] - est[[b]],
                             se = sqrt(V[a, a] + V[b, b] - 2 * V[a, b])))
  }
  r[, `:=`(lower = estimate - 1.96 * se, upper = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)), note = note)]
  r[]
}
suppress_pair <- function(a, b) if (is.na(a) || is.na(b) || a < CELL_MIN || b < CELL_MIN) c(NA_integer_, NA_integer_) else c(a, b)
cell_info <- function(s, ev) {
  y <- s[[ev]]; A <- s$A
  y_ok <- !is.na(y)
  na1 <- sum(A == 1 & y_ok); na0 <- sum(A == 0 & y_ok)
  e1 <- sum(y[A == 1 & y_ok] == 1L); e0 <- sum(y[A == 0 & y_ok] == 1L)
  ok <- all(c(na1 - e1, na0 - e0, e1, e0) >= CELL_MIN)
  np <- suppress_pair(na1, na0)
  ep <- if (ok) c(e1, e0) else c(NA_integer_, NA_integer_)
  list(ok = ok, n = if (sum(y_ok) < CELL_MIN) NA_integer_ else sum(y_ok),
       n_A1 = np[1], n_A0 = np[2], events = if (e1 + e0 < CELL_MIN) NA_integer_ else e1 + e0,
       events_A1 = ep[1], events_A0 = ep[2])
}
fit_logit <- function(f, s, weighted, terms, contrast = NULL) {
  r <- safely(if (weighted) glm(f, family = quasibinomial, data = s, weights = w_iptw)
              else glm(f, family = binomial, data = s))
  if (is.null(r$val)) return(na_rows(terms, contrast, paste(r$warn, collapse = "; ")))
  if (!isTRUE(r$val$converged)) return(na_rows(terms, contrast, "glm did not converge"))
  V <- sandwich::vcovHC(r$val, type = "HC0")
  coef_rows(coef(r$val), V, terms, contrast, note = paste(unique(r$warn), collapse = "; "))
}
fit_cox <- function(f, s, weighted, terms, contrast = NULL) {
  r <- safely(if (weighted) coxph(f, data = s, weights = w_iptw, robust = TRUE) else coxph(f, data = s, robust = TRUE))
  if (is.null(r$val)) return(na_rows(terms, contrast, paste(r$warn, collapse = "; ")))
  if (any(grepl("converge|infinite|Ran out", r$warn))) return(na_rows(terms, contrast, paste(unique(r$warn), collapse = "; ")))
  coef_rows(coef(r$val), vcov(r$val), terms, contrast, note = paste(unique(r$warn), collapse = "; "))
}

# Lasso propensity as in the port. Balance is scored on EVERY column of the propensity model, not the port's
# five: "every weighted SMD < 0.1" is the rule, and checking five of forty could not establish it.
fit_iptw <- function(s, gs, covars, nfolds) {
  empty <- list(ok = FALSE, w = rep(NA_real_, nrow(s)), bal = data.table(), verdict = "NOT_EVALUABLE", worst = NA_character_,
                max_smd = NA_real_, diag = data.table())
  if (sum(s$A == 1) < 2 * CELL_MIN || sum(s$A == 0) < 2 * CELL_MIN) return(empty)
  X <- model.matrix(as.formula(paste("~", paste(c("age0", "sex", "log_LDL0", covars, gs), collapse = " + "))), s)[, -1, drop = FALSE]
  set.seed(PORT_SEED)
  cvr <- safely(cv.glmnet(X, s$A, family = "binomial", alpha = 1, nfolds = nfolds))
  if (is.null(cvr$val)) { message("  lasso failed: ", paste(cvr$warn, collapse = "; ")); return(empty) }
  pA <- as.numeric(predict(cvr$val, X, s = "lambda.min", type = "response"))
  pm <- mean(s$A)
  w <- ifelse(s$A == 1, pm / pA, (1 - pm) / (1 - pA))
  if (any(!is.finite(w))) { message("  non-finite weights"); return(empty) }
  qs <- quantile(w, c(0.01, 0.99)); w <- pmin(pmax(w, qs[1]), qs[2])
  i1 <- s$A == 1
  bal <- rbindlist(lapply(colnames(X), function(v) {
    x <- X[, v]
    v1 <- var(x[i1]); v0 <- var(x[!i1]); sdp <- sqrt((v1 + v0) / 2)
    if (!is.finite(sdp) || sdp == 0) return(NULL)      # a constant column cannot be imbalanced
    data.table(variable = v,
               smd_raw = (mean(x[i1]) - mean(x[!i1])) / sdp,
               smd_weighted = (weighted.mean(x[i1], w[i1]) - weighted.mean(x[!i1], w[!i1])) / sdp)
  }))
  ms <- max(abs(bal$smd_weighted)); worst <- bal$variable[which.max(abs(bal$smd_weighted))]
  verdict <- if (anyNA(bal$smd_weighted)) "NOT_EVALUABLE" else if (ms >= SMD_LIMIT) "UNINFORMATIVE" else "BALANCED"
  big <- sum(i1) >= 100 && sum(!i1) >= 100     # quantiles of very small arms are effectively individual values
  dg <- data.table(
    n_over_limit = sum(abs(bal$smd_weighted) >= SMD_LIMIT), n_columns = nrow(bal),
    w_q01 = qs[[1]], w_median = median(w), w_mean = mean(w), w_q99 = qs[[2]],
    pA_A1_q01 = if (big) quantile(pA[i1], .01)[[1]] else NA_real_, pA_A1_q50 = if (big) median(pA[i1]) else NA_real_,
    pA_A1_q99 = if (big) quantile(pA[i1], .99)[[1]] else NA_real_, pA_A0_q01 = if (big) quantile(pA[!i1], .01)[[1]] else NA_real_,
    pA_A0_q50 = if (big) median(pA[!i1]) else NA_real_, pA_A0_q99 = if (big) quantile(pA[!i1], .99)[[1]] else NA_real_,
    frac_A1_above_A0_q99 = if (big) mean(pA[i1] > quantile(pA[!i1], .99)) else NA_real_)
  list(ok = TRUE, w = w, bal = bal, verdict = verdict, worst = worst, max_smd = ms, diag = dg)
}

run_one <- function(s, sp, gs, covars, label) {
  terms <- c("A", gs, paste0("A_x_", gs))
  contrast <- if (length(gs) == 2L) paste0("A_x_", gs) else NULL
  base_rhs <- paste(c("age0", "sex", "A", gs, paste0("A_x_", gs)), collapse = " + ")
  adj_rhs  <- paste(c("age0", "sex", "log_LDL0", covars, "A", gs, paste0("A_x_", gs)), collapse = " + ")
  iptw <- fit_iptw(s, gs, covars, opt$nfolds)
  s[, w_iptw := iptw$w]
  ci_b <- cell_info(s, "Ybin"); ci_t <- cell_info(s, "Yevent")
  deg <- function(ci, tt) if (!ci$ok) na_rows(tt, contrast, "degenerate design: an exposure-by-outcome cell < 5") else NULL
  noW <- function(tt) na_rows(tt, contrast, "IPTW not estimable: propensity model failed or an arm too small")
  mk <- function(rows, fam, est, ci, read) cbind(
    data.table(spec = sp$id, outcome = sp$outcome, exposure = sp$exposure, sample = sp$sample, G = label,
               family = fam, estimator = est), rows,
    data.table(n = ci$n, n_A1 = ci$n_A1, n_A0 = ci$n_A0, events = ci$events, events_A1 = ci$events_A1,
               events_A0 = ci$events_A0, iptw_balance = iptw$verdict, read = read))
  out <- list()
  yb <- function(r) as.formula(paste("Ybin ~", r)); yt <- function(r) as.formula(paste("Surv(Ytime, Yevent) ~", r))
  add <- function(rows) out[[length(out) + 1L]] <<- rows
  add(mk(if (is.null(d <- deg(ci_b, terms))) fit_logit(yb(base_rhs), s, FALSE, terms, contrast) else d, "logistic", "Unadjusted", ci_b, "conditional on age and sex only"))
  add(mk(if (is.null(d <- deg(ci_b, terms))) fit_logit(yb(adj_rhs),  s, FALSE, terms, contrast) else d, "logistic", "Adjusted_LDL0_covariates", ci_b, "conditions on suspected colliders"))
  iread <- if (iptw$verdict == "BALANCED") "readable: every weighted SMD < 0.1" else "UNINFORMATIVE by pre-specified rule; supports neither reading"
  add(mk(if (!is.null(d <- deg(ci_b, terms))) d else if (!iptw$ok) noW(terms) else fit_logit(yb(base_rhs), s, TRUE, terms, contrast), "logistic", "IPTW_MSM", ci_b, iread))
  add(mk(if (is.null(d <- deg(ci_t, terms))) fit_cox(yt(base_rhs), s, FALSE, terms, contrast) else d, "cox", "Unadjusted", ci_t, "conditional on age and sex only"))
  add(mk(if (is.null(d <- deg(ci_t, terms))) fit_cox(yt(adj_rhs),  s, FALSE, terms, contrast) else d, "cox", "Adjusted_LDL0_covariates", ci_t, "conditions on suspected colliders"))
  add(mk(if (!is.null(d <- deg(ci_t, terms))) d else if (!iptw$ok) noW(terms) else fit_cox(yt(base_rhs), s, TRUE, terms, contrast), "cox", "IPTW_MSM", ci_t, iread))
  res <- rbindlist(out)
  res[, `:=`(ratio = exp(estimate), ratio_lower = exp(lower), ratio_upper = exp(upper))]

  strat <- NULL
  if (length(gs) == 1L && !anyNA(s$ldl_tertile)) {
    aterm <- paste0("A_x_", gs)
    strat <- rbindlist(lapply(levels(s$ldl_tertile), function(tt) {
      st <- s[ldl_tertile == tt]
      cb <- cell_info(st, "Ybin"); ct <- cell_info(st, "Yevent")
      one <- function(ci, fn, f, weighted, fam, est) {
        rows <- if (!ci$ok) na_rows(aterm, NULL, "degenerate design in stratum: a cell < 5")
                else if (weighted && !iptw$ok) noW(aterm) else fn(f, st, weighted, aterm)
        cbind(data.table(spec = sp$id, G = label, stratum = tt, family = fam, estimator = est), rows,
              data.table(n = ci$n, n_A1 = ci$n_A1, n_A0 = ci$n_A0, events = ci$events,
                         read = if (weighted) "descriptive: global weights, stratum balance not evaluated" else "descriptive"))
      }
      rbind(one(cb, fit_logit, yb(base_rhs), FALSE, "logistic", "Unadjusted"),
            one(cb, fit_logit, yb(base_rhs), TRUE,  "logistic", "IPTW_MSM"),
            one(ct, fit_cox,   yt(base_rhs), FALSE, "cox", "Unadjusted"),
            one(ct, fit_cox,   yt(base_rhs), TRUE,  "cox", "IPTW_MSM"))
    }))
    strat[, `:=`(ratio = exp(estimate), ratio_lower = exp(lower), ratio_upper = exp(upper))]
  }
  bal <- if (nrow(iptw$bal)) cbind(data.table(spec = sp$id, G = label), iptw$bal)[, exceeds_limit := abs(smd_weighted) >= SMD_LIMIT] else NULL
  dg <- cbind(data.table(spec = sp$id, G = label, verdict = iptw$verdict, max_abs_smd = iptw$max_smd, worst_variable = iptw$worst),
              if (nrow(iptw$diag)) iptw$diag else NULL)
  s[, w_iptw := NULL]
  gc(FALSE)
  list(res = res, strat = strat, bal = bal, diag = dg)
}

# ---- coronary arm --------------------------------------------------------------------------------------
all_res <- list(); all_strat <- list(); all_bal <- list(); all_diag <- list(); sample_rows <- list()
for (sp in SPECS) {
  stamp("spec ", sp$id)
  s <- spec_frame(sp)
  covars <- active_covars(s, COVARS)
  ci_b <- cell_info(s, "Ybin"); ci_t <- cell_info(s, "Yevent")
  sample_rows[[sp$id]] <- data.table(spec = sp$id, outcome = sp$outcome, exposure = sp$exposure, sample = sp$sample,
    n = ci_b$n, n_A1 = ci_b$n_A1, n_A0 = ci_b$n_A0, events_logistic = ci_b$events, events_logistic_A1 = ci_b$events_A1,
    events_cox = ci_t$events, events_cox_A1 = ci_t$events_A1,
    prior_CAD_retained = { x <- sum(s$prior_CAD == "1"); if (x < CELL_MIN) NA_integer_ else x })
  runs <- c(lapply(ALL_G, function(g) list(gs = g, label = g)),
            list(list(gs = c("score_LDL", "score_LPA"), label = "mutual:score_LDL+score_LPA")))
  for (r in runs) {
    stamp("  ", r$label)
    o <- run_one(s, sp, r$gs, covars, r$label)
    all_res[[length(all_res) + 1L]] <- o$res
    if (!is.null(o$strat)) all_strat[[length(all_strat) + 1L]] <- o$strat
    if (!is.null(o$bal)) all_bal[[length(all_bal) + 1L]] <- o$bal
    all_diag[[length(all_diag) + 1L]] <- o$diag
  }
  rm(s); gc(FALSE)
}
RES <- rbindlist(all_res, fill = TRUE); STRAT <- rbindlist(all_strat, fill = TRUE)
BAL <- rbindlist(all_bal, fill = TRUE); DIAG <- rbindlist(all_diag, fill = TRUE)
write_hdr(RES, "coronary_estimators.csv")
write_hdr(STRAT, "coronary_by_ldl_tertile.csv")
write_hdr(BAL, "coronary_balance_tables.csv")
write_hdr(DIAG, "coronary_balance_verdicts.csv")
write_hdr(rbindlist(sample_rows), "coronary_samples.csv")

# ---- ladder from the published paper estimate to the published port estimate -----------------------------
# Same model family throughout (Cox), one change per step, so each step's shift is attributable.
ladder_rows <- list()
add_ladder <- function(step, what, G, s, time, ev, rhs) {
  ci <- cell_info(s, ev)
  aterm <- paste0("A_x_", G)
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time, ev, rhs))
  r <- if (!ci$ok) na_rows(aterm, NULL, "degenerate design: a cell < 5") else {
    v <- safely(coxph(f, data = s))
    if (is.null(v$val) || any(grepl("converge|infinite|Ran out", v$warn))) na_rows(aterm, NULL, paste(unique(v$warn), collapse = "; "))
    else coef_rows(coef(v$val), vcov(v$val), aterm)
  }
  ladder_rows[[length(ladder_rows) + 1L]] <<- cbind(data.table(track = NA_character_, step = step, what_changed = what, G = G), r,
                                                    data.table(n = ci$n, events = ci$events))
}
PCS <- paste(paste0("pc", 1:10), collapse = " + ")
s_c4 <- spec_frame(SPECS[[4]]); s_c3 <- spec_frame(SPECS[[3]]); s_c2 <- spec_frame(SPECS[[2]]); s_c1 <- spec_frame(SPECS[[1]])
for (trk in c("score_LDL", "score_LPA")) {
  n0 <- length(ladder_rows)
  pm <- function(g) paste(c("A", g, paste0("A_x_", g), "age0", "sex", PCS), collapse = " + ")
  add_ladder(1, "paper model, exposure, outcome, score on the port's restricted sample (white British, LDL0 range, complete covariates)", trk, s_c4, "fu", "event", pm(trk))
  add_ladder(2, "+ follow-up truncated at 10 y (port's horizon)", trk, s_c4, "fu10", "event10", pm(trk))
  add_ladder(3, "+ exposure -> field 20003 statin codes", trk, s_c3, "fu10", "event10", pm(trk))
  add_ladder(4, "+ outcome -> acute coronary I21-I24", trk, s_c2, "cad_time_years", "CAD1", pm(trk))
  add_ladder(5, "+ sample -> prior CAD kept (port as published)", trk, s_c1, "cad_time_years", "CAD1", pm(trk))
  g_end <- if (trk == "score_LDL") "ldl_prs_z" else trk
  if (trk == "score_LDL") add_ladder(6, "+ score -> UKB standard LDL PRS 26250", g_end, s_c1, "cad_time_years", "CAD1", pm(g_end))
  add_ladder(7, "+ model -> age and sex only (port's unadjusted; 0.935 for ldl_prs_z before re-standardisation)", g_end, s_c1, "cad_time_years", "CAD1",
             paste(c("age0", "sex", "A", g_end, paste0("A_x_", g_end)), collapse = " + "))
  for (i in seq(n0 + 1L, length(ladder_rows))) ladder_rows[[i]][, track := trk]
}
rm(s_c1, s_c2, s_c3, s_c4); gc(FALSE)
LAD <- rbindlist(ladder_rows, fill = TRUE)
if (!is.null(pc2)) {
  L0 <- pc2$ladder0[, .(track = G, step = 0, what_changed = "published paper model, own sample (Table S35)", G, term = paste0(G, ":statin"),
                        estimate, se, lower = estimate - 1.96 * se, upper = estimate + 1.96 * se,
                        p = 2 * pnorm(-abs(estimate / se)), note = "", n, events)]
  LAD <- rbind(L0, LAD, fill = TRUE)
}
LAD[, `:=`(HR = exp(estimate), HR_lower = exp(lower), HR_upper = exp(upper))]
setorder(LAD, track, step)
write_hdr(LAD, "ladder_paper_to_port.csv")

# ---- repeat-visit arm ------------------------------------------------------------------------------------
stamp("repeat-visit arm")
rv <- cohort[!is.na(LDL0) & !is.na(LDL1) & LDL0 > 0.5 & LDL0 < 10 & LDL1 > 0.5 & LDL1 < 10]
if ("white_british" %in% names(rv)) rv <- rv[white_british == 1]
rv <- rv[!is.na(statin0) & !is.na(statin1) & statin0 == 0]
rv <- rv[complete.cases(rv[, ..REPEAT_G])]
rv[, A := as.integer(statin1 == 1)]
rv[, log_LDL0 := log(LDL0)][, log_LDL1 := log(LDL1)]
rv[, sex := factor(sex)]
RV_COV <- intersect(c("bmi", "townsend", "diabetes_self", "smoking", "prior_CAD", paste0("pc", 1:10)), names(rv))
for (v in intersect(c("smoking", "diabetes_self", "prior_CAD"), RV_COV)) rv[, (v) := factor(get(v))]
rv <- rv[complete.cases(rv[, c(RV_COV, "age0"), with = FALSE])]
if (nrow(rv)) {
  rv <- droplevels(rv)
  for (g in REPEAT_G) if (isTRUE(sd(rv[[g]]) > 0)) rv[, (g) := as.numeric(scale(get(g)))]
  rv[, ldl_tertile := tryCatch(cut(LDL0, quantile(LDL0, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("T1", "T2", "T3")),
                               error = function(e) factor(rep(NA_character_, .N), levels = c("T1", "T2", "T3")))]
}
n_ini <- sum(rv$A == 1L); n_non <- sum(rv$A == 0L)
ini_ok <- n_ini >= opt$min_initiators_na && n_non >= opt$min_initiators_na
message(sprintf("repeat-visit sample %d; initiators %s", nrow(rv), if (n_ini < CELL_MIN) "<5" else n_ini))
RV_COV_A <- if (nrow(rv)) active_covars(rv, RV_COV) else character()

fit_lm <- function(f, s, weighted, terms) {
  r <- safely(if (weighted) lm(f, data = s, weights = w_iptw) else lm(f, data = s))
  if (is.null(r$val)) return(na_rows(terms, NULL, paste(r$warn, collapse = "; ")))
  coef_rows(coef(r$val), sandwich::vcovHC(r$val, type = "HC0"), terms, NULL, paste(unique(r$warn), collapse = "; "))
}
rv_int <- list(); rv_str <- list(); rv_chg <- list(); rv_bal <- list(); rv_dia <- list(); rv_eff <- list()
for (g in REPEAT_G) {
  note0 <- sprintf("too few initiators (< %d) or non-initiators", opt$min_initiators_na)
  if (!ini_ok) {
    rv_int[[g]] <- cbind(data.table(G = g, estimator = c("Unadjusted", "Adjusted_LDL0_covariates", "IPTW_MSM")),
                         rbindlist(replicate(3, na_rows("A_x_G", NULL, note0), simplify = FALSE)))
    rv_chg[[g]] <- cbind(data.table(G = g, change = c("Absolute_mmol_L", "Relative_log", "Relative_log_adjusted_LDL0")),
                         rbindlist(replicate(3, na_rows("G", NULL, note0), simplify = FALSE)))
    rv_dia[[g]] <- data.table(G = g, verdict = "NOT_EVALUABLE", max_abs_smd = NA_real_, worst_variable = NA_character_)
    rv_eff[[g]] <- cbind(data.table(G = g, effect = c("log_LDL0_per_SD", "LDL0_mmol_L_per_SD")),
                         rbindlist(replicate(2, na_rows("G", NULL, note0), simplify = FALSE)))
    next
  }
  s <- copy(rv)
  s[, G := get(g)][, A_x_G := A * G]
  iptw <- fit_iptw(s, "G", RV_COV_A, opt$nfolds)
  s[, w_iptw := iptw$w]
  adj <- as.formula(paste("log_LDL1 ~", paste(c("age0", "sex", "log_LDL0", RV_COV_A, "A + G + A_x_G"), collapse = " + ")))
  msm <- log_LDL1 ~ age0 + sex + A + G + A_x_G
  tt <- c("A", "G", "A_x_G")
  rv_int[[g]] <- rbind(
    cbind(data.table(G = g, estimator = "Unadjusted"), fit_lm(log_LDL1 ~ A + G + A_x_G, s, FALSE, tt)),
    cbind(data.table(G = g, estimator = "Adjusted_LDL0_covariates"), fit_lm(adj, s, FALSE, tt)),
    cbind(data.table(G = g, estimator = "IPTW_MSM"),
          if (iptw$ok) fit_lm(msm, s, TRUE, tt) else na_rows(tt, NULL, "IPTW not estimable")))[, iptw_balance := iptw$verdict]
  if (!anyNA(s$ldl_tertile)) rv_str[[g]] <- rbindlist(lapply(levels(s$ldl_tertile), function(tl) {
    st <- s[ldl_tertile == tl]
    a1 <- sum(st$A == 1L); a0 <- sum(st$A == 0L)
    good <- a1 >= CELL_MIN && a0 >= CELL_MIN
    na_n <- na_rows("A_x_G", NULL, "a stratum arm < 5")
    rbind(cbind(data.table(G = g, stratum = tl, estimator = "Unadjusted"), if (good) fit_lm(log_LDL1 ~ A + G + A_x_G, st, FALSE, "A_x_G") else na_n),
          cbind(data.table(G = g, stratum = tl, estimator = "IPTW_MSM"), if (good && iptw$ok) fit_lm(msm, st, TRUE, "A_x_G") else na_n))[
            , `:=`(n = if (good) nrow(st) else NA_integer_, initiators = if (good) a1 else NA_integer_)]
  }))
  ini <- s[A == 1L]
  ini[, abs_change := LDL1 - LDL0][, rel_change := log_LDL1 - log_LDL0]
  rv_chg[[g]] <- rbind(
    cbind(data.table(G = g, change = "Absolute_mmol_L"), fit_lm(abs_change ~ G, ini, FALSE, "G")),
    cbind(data.table(G = g, change = "Relative_log"), fit_lm(rel_change ~ G, ini, FALSE, "G")),
    cbind(data.table(G = g, change = "Relative_log_adjusted_LDL0"), fit_lm(rel_change ~ G + log_LDL0, ini, FALSE, "G")))[, n := nrow(ini)]
  # The score's own effect on LDL-C sizes the expected collider signal: regression to the mean scales with it.
  # Measured LDL-C carries Lp(a)-cholesterol, so an Lp(a) score is not literally off-target, and this row
  # states by how much rather than leaving "off-target" as an assertion.
  rv_eff[[g]] <- rbind(
    cbind(data.table(G = g, effect = "log_LDL0_per_SD"), fit_lm(log_LDL0 ~ G, s, FALSE, "G")),
    cbind(data.table(G = g, effect = "LDL0_mmol_L_per_SD"), fit_lm(LDL0 ~ G, s, FALSE, "G")))
  if (nrow(iptw$bal)) rv_bal[[g]] <- cbind(data.table(G = g), iptw$bal)[, exceeds_limit := abs(smd_weighted) >= SMD_LIMIT]
  rv_dia[[g]] <- cbind(data.table(G = g, verdict = iptw$verdict, max_abs_smd = iptw$max_smd, worst_variable = iptw$worst),
                       if (nrow(iptw$diag)) iptw$diag else NULL)
}
RV_INT <- rbindlist(rv_int, fill = TRUE); RV_CHG <- rbindlist(rv_chg, fill = TRUE)
RV_EFF <- rbindlist(rv_eff, fill = TRUE); RV_DIA <- rbindlist(rv_dia, fill = TRUE)
RV_STR <- rbindlist(rv_str, fill = TRUE); RV_BAL <- rbindlist(rv_bal, fill = TRUE)
write_hdr(RV_INT, "repeat_ldl_interaction.csv")
write_hdr(RV_STR, "repeat_ldl_by_ldl0_tertile.csv")
write_hdr(RV_CHG, "repeat_ldl_change_scores.csv")
write_hdr(RV_EFF, "repeat_score_effect_on_ldl0.csv")
write_hdr(RV_BAL, "repeat_balance_tables.csv")
write_hdr(RV_DIA, "repeat_balance_verdicts.csv")
write_hdr(data.table(n = if (nrow(rv) < CELL_MIN) NA_integer_ else nrow(rv), initiators = if (n_ini < CELL_MIN) NA_integer_ else n_ini,
                     floor_NA = opt$min_initiators_na, floor_uninformative = PLAN_MIN_INITIATORS),
          "repeat_sample.csv")

# ---- optional arm: measured repeat Lp(a) (30790-1.0), reportability explicit -------------------------------
# UKB coding 4917: 4 = too low, 5 = too high (verified against LPA genotypes, not the 3/4 one might assume).
# Lp(a) is missing about 17% of the time in both tails, and that missingness is a function of LPA genotype.
# Dropping is therefore selection on the outcome and on G, which is exactly the structure this analysis is
# about. Both handling rules are reported: "drop" and "bound" (limit values set to the observed extreme of that
# visit, keeping the people but flattening their values). Neither is a censored-data model; no readout rule
# depends on this arm, and it is descriptive.
lpa_rows <- list(); lpa_cnt <- list()
lpa_files <- file.path(opt$raw_dir, c("fid30790.csv", "fid30796.csv"))
lpa_NA <- function(mode, g, note) {
  rbind(cbind(data.table(handling = mode, G = g, model = c("Unadjusted", "Adjusted_LPA0_covariates")),
              rbindlist(replicate(2, na_rows("A_x_G", NULL, note), simplify = FALSE))),
        cbind(data.table(handling = mode, G = g, model = c("Change_relative", "Change_relative_adjusted_LPA0")),
              rbindlist(replicate(2, na_rows("G", NULL, note), simplify = FALSE))))
}
lpa_raw <- if (all(file.exists(lpa_files)) && nrow(rv))
  tryCatch(list(v0 = rdcsv("30790", 0L), v1 = rdcsv("30790", 1L), r0 = rdcsv("30796", 0L), r1 = rdcsv("30796", 1L)),
           error = function(e) { message("Lp(a) arm: ", conditionMessage(e)); NULL })
if (is.null(lpa_raw)) {
  note <- "Lp(a) or reportability field absent, or empty repeat sample; arm not run"
  message("Lp(a) arm: ", note)
  lpa_rows <- lapply(c("score_LPA", "score_LDL"), function(g) lpa_NA("not_run", g, note))
} else {
  v0 <- lpa_raw$v0; v1 <- lpa_raw$v1; r0 <- lpa_raw$r0; r1 <- lpa_raw$r1
  id <- as.character(rv$eid)
  rv[, `:=`(lpa0 = as.numeric(v0$v[match(id, v0$IID)]), lpa1 = as.numeric(v1$v[match(id, v1$IID)]),
            rep0 = as.integer(r0$v[match(id, r0$IID)]), rep1 = as.integer(r1$v[match(id, r1$IID)]))]
  rm(v0, v1, r0, r1, id)
  sup <- function(x) if (x < CELL_MIN) NA_integer_ else as.integer(x)
  lo0 <- min(rv$lpa0, na.rm = TRUE); hi0 <- max(rv$lpa0, na.rm = TRUE)
  lo1 <- min(rv$lpa1, na.rm = TRUE); hi1 <- max(rv$lpa1, na.rm = TRUE)
  lpa_cnt[[1]] <- data.table(
    n_repeat_sample = sup(nrow(rv)),
    # UKB data-coding 4917: 2 and 4 too low, 3 and 5 too high.
    v0_too_low = sup(rv[rep0 %in% c(2L, 4L), .N]), v0_too_high = sup(rv[rep0 %in% c(3L, 5L), .N]),
    v1_too_low = sup(rv[rep1 %in% c(2L, 4L), .N]), v1_too_high = sup(rv[rep1 %in% c(3L, 5L), .N]),
    v1_too_low_initiators = sup(rv[rep1 %in% c(2L, 4L) & A == 1L, .N]), v1_too_high_initiators = sup(rv[rep1 %in% c(3L, 5L) & A == 1L, .N]))
  for (mode in c("drop", "bound")) {
    s <- copy(rv)
    if (mode == "bound") {
      s[is.na(lpa0) & rep0 %in% c(2L, 4L), lpa0 := lo0]; s[is.na(lpa0) & rep0 %in% c(3L, 5L), lpa0 := hi0]
      s[is.na(lpa1) & rep1 %in% c(2L, 4L), lpa1 := lo1]; s[is.na(lpa1) & rep1 %in% c(3L, 5L), lpa1 := hi1]
    } else {
      # Reportability is also checked directly in case a limit value was stored rather than left blank.
      s <- s[(is.na(rep0) | !rep0 %in% 2:5) & (is.na(rep1) | !rep1 %in% 2:5)]
    }
    s <- s[!is.na(lpa0) & !is.na(lpa1) & lpa0 > 0 & lpa1 > 0]
    s[, `:=`(log_lpa0 = log(lpa0), log_lpa1 = log(lpa1))]
    ai <- sum(s$A == 1L); an <- sum(s$A == 0L)
    lpa_cnt[[length(lpa_cnt) + 1L]] <- data.table(handling = mode, n = sup(nrow(s)), initiators = sup(ai))
    for (g in c("score_LPA", "score_LDL")) {
      if (ai < opt$min_initiators_na || an < opt$min_initiators_na) {
        lpa_rows[[length(lpa_rows) + 1L]] <- lpa_NA(mode, g, sprintf("too few initiators (< %d) with repeat Lp(a)", opt$min_initiators_na))
        next
      }
      st <- copy(s)
      st <- droplevels(st)
      cv_l <- active_covars(st, RV_COV)
      st[, G := as.numeric(scale(get(g)))][, A_x_G := A * G]
      fa <- as.formula(paste("log_lpa1 ~", paste(c("age0", "sex", "log_lpa0", cv_l, "A + G + A_x_G"), collapse = " + ")))
      ini <- st[A == 1L]; ini[, rel := log_lpa1 - log_lpa0]
      lpa_rows[[length(lpa_rows) + 1L]] <- rbind(
        cbind(data.table(handling = mode, G = g, model = "Unadjusted"), fit_lm(log_lpa1 ~ A + G + A_x_G, st, FALSE, "A_x_G")),
        cbind(data.table(handling = mode, G = g, model = "Adjusted_LPA0_covariates"), fit_lm(fa, st, FALSE, "A_x_G")),
        cbind(data.table(handling = mode, G = g, model = "Change_relative"), fit_lm(rel ~ G, ini, FALSE, "G")),
        cbind(data.table(handling = mode, G = g, model = "Change_relative_adjusted_LPA0"), fit_lm(rel ~ G + log_lpa0, ini, FALSE, "G")))
    }
  }
}
write_hdr(rbindlist(lpa_rows, fill = TRUE), "repeat_lpa_optional.csv")
if (length(lpa_cnt)) write_hdr(rbindlist(lpa_cnt, fill = TRUE), "repeat_lpa_counts.csv")

# ---- verdict lines ---------------------------------------------------------------------------------------
fmt <- function(est, lo, hi, dp = 3) if (is.na(est)) "NA" else sprintf(paste0("%.", dp, "f (%.", dp, "f, %.", dp, "f)"), est, lo, hi)
excl <- function(lo, hi) !is.na(lo) && !is.na(hi) && (lo > 0 || hi < 0)
V <- c(HDR, "", paste("Positive controls:", pc_status),
       sprintf("Port sample %d -> analysis sample %d (complete on all six scores). The port's published coronary rows kept prior CAD; C2-C4 remove it.", n_port, nrow(base)),
       "IPTW rule (readout D; ANALYSIS_PLAN.md section 6): any weighted SMD >= 0.1 across all propensity columns means UNINFORMATIVE, not support for either reading.",
       "BALANCED is necessary, not sufficient: baseline LDL-C in prevalent users is on-treatment, so balance on it is not exchangeability on the pre-treatment LDL-C. A readable IPTW row is still an association.",
       "", "== Coronary arm: A x G per SD, ratio scale, by spec ==")
for (sp in SPECS) for (g in c(ALL_G)) for (fam in c("logistic", "cox")) {
  pick <- function(e) RES[spec == sp$id & G == g & family == fam & estimator == e & term == paste0("A_x_", g)]
  u <- pick("Unadjusted"); a <- pick("Adjusted_LDL0_covariates"); i <- pick("IPTW_MSM")
  dv <- DIAG[spec == sp$id & G == g]
  iv <- if (dv$verdict == "BALANCED")
    sprintf("BALANCED (max |SMD| %.3f): IPTW CI %s 1", dv$max_abs_smd, if (excl(i$lower, i$upper)) "EXCLUDES" else "includes")
  else if (dv$verdict == "UNINFORMATIVE")
    sprintf("UNINFORMATIVE (max |SMD| %.3f on %s; %d of %d columns over limit)", dv$max_abs_smd, dv$worst_variable, dv$n_over_limit, dv$n_columns)
  else "NOT_EVALUABLE"
  V <- c(V, sprintf("[%s][%s][%s] unadj %s | adj %s | IPTW %s -> %s", sp$id, g, fam,
                    fmt(u$ratio, u$ratio_lower, u$ratio_upper), fmt(a$ratio, a$ratio_lower, a$ratio_upper),
                    fmt(i$ratio, i$ratio_lower, i$ratio_upper), iv))
}
V <- c(V, "", "== Asymmetry: A x score_LDL minus A x score_LPA, mutually adjusted, log scale ==")
for (sp in SPECS) for (fam in c("logistic", "cox")) for (e in c("Unadjusted", "Adjusted_LDL0_covariates", "IPTW_MSM")) {
  r <- RES[spec == sp$id & G == "mutual:score_LDL+score_LPA" & family == fam & estimator == e & term == "A_x_score_LDL - A_x_score_LPA"]
  dv <- DIAG[spec == sp$id & G == "mutual:score_LDL+score_LPA"]
  tag <- if (is.na(r$estimate)) "NOT_EVALUABLE" else
    if (e == "IPTW_MSM" && dv$verdict != "BALANCED") sprintf("UNINFORMATIVE (max |SMD| %.3f)", dv$max_abs_smd) else
    if (excl(r$lower, r$upper)) "CI excludes 0" else "CI includes 0"
  V <- c(V, sprintf("[%s][%s][%s] %s -> %s", sp$id, fam, e, fmt(r$estimate, r$lower, r$upper), tag))
}
V <- c(V, "", "== Pre-specified coronary reading per spec, port's two PRS variants plus the paper's scores (ANALYSIS_PLAN.md section 6) ==",
       "The plan's sign language was written for the LDL arm. Coronary interactions in the port run negative, so 'IPTW interaction present' is read two-sided.")
for (sp in SPECS) {
  gg <- c("ldl_prs_z", "cad_prs_z", "score_LDL", "score_LPA")
  dv <- DIAG[spec == sp$id & G %in% gg]
  if (any(dv$verdict == "NOT_EVALUABLE")) { V <- c(V, sprintf("[%s] NOT_EVALUABLE", sp$id)); next }
  ng <- RES[spec == sp$id & G == "neg_prs_z" & term == "A_x_neg_prs_z"]
  ng_sig <- ng[!is.na(estimate) & (estimator != "IPTW_MSM" | iptw_balance == "BALANCED") & (lower > 0 | upper < 0)]
  ng_txt <- if (nrow(ng_sig)) sprintf("NEGATIVE CONTROL SIGNAL in %d row(s): undermines any reading of this spec", nrow(ng_sig)) else "negative control null"
  if (any(dv$verdict != "BALANCED")) {
    V <- c(V, sprintf("[%s] UNINFORMATIVE: weighted SMD >= %.1f in %s; %s", sp$id, SMD_LIMIT,
                      paste(dv[verdict != "BALANCED", G], collapse = ", "), ng_txt))
  } else {
    ip <- RES[spec == sp$id & G %in% gg & estimator == "IPTW_MSM" & grepl("^A_x_", term) & !is.na(estimate) & (lower > 0 | upper < 0)]
    # ANALYSIS_PLAN.md section 6 also calls poor overlap uninformative but sets no number, so it is shown, not ruled on.
    V <- c(V, sprintf("[%s] BALANCED; IPTW interaction present for: %s; %s; overlap, largest share of users above the non-users' 99th propensity percentile: %.3f", sp$id,
                      if (nrow(ip)) paste(unique(ip$G), collapse = ", ") else "none", ng_txt, max(dv$frac_A1_above_A0_q99, na.rm = TRUE)))
  }
}
V <- c(V, "", "== Repeat-visit arm: rule = null full-sample interaction + positive LDL0-adjusted change score supports the collider mechanism ==",
       "Full-sample null means the unadjusted CI includes 0 and, when balance allows reading it, the IPTW CI includes 0. A CI containing 0 is a null only to the width shown.")
for (g in REPEAT_G) {
  u <- RV_INT[G == g & estimator == "Unadjusted" & term == "A_x_G"]; a <- RV_INT[G == g & estimator == "Adjusted_LDL0_covariates" & term == "A_x_G"]
  i <- RV_INT[G == g & estimator == "IPTW_MSM" & term == "A_x_G"]
  c0 <- RV_CHG[G == g & change == "Relative_log"]; c1 <- RV_CHG[G == g & change == "Relative_log_adjusted_LDL0"]
  ef <- RV_EFF[G == g & effect == "log_LDL0_per_SD"]; dv <- RV_DIA[G == g]
  v <- if (!ini_ok || is.na(u$estimate) || is.na(c1$estimate)) "NOT_EVALUABLE (NA)"
       else if (n_ini < PLAN_MIN_INITIATORS) sprintf("UNINFORMATIVE: %d initiators < %d", n_ini, PLAN_MIN_INITIATORS)
       else {
         iptw_null <- if (dv$verdict == "BALANCED") !excl(i$lower, i$upper) else NA
         null_full <- !excl(u$lower, u$upper) && !isFALSE(iptw_null)
         pos <- c1$lower > 0
         paste0(if (null_full && pos) "SUPPORTS_COLLIDER" else if (!null_full) "NOT_SUPPORTED: full-sample interaction not null" else "NOT_SUPPORTED: adjusted change score not positive",
                if (is.na(iptw_null)) sprintf(" [IPTW limb unreadable, max |SMD| %.3f]", dv$max_abs_smd) else "")
       }
  V <- c(V, sprintf("[%s] unadj %s | adj %s | IPTW %s | change %s -> LDL0-adj %s | score on log LDL0 %s -> %s", g,
                    fmt(u$estimate, u$lower, u$upper), fmt(a$estimate, a$lower, a$upper), fmt(i$estimate, i$lower, i$upper),
                    fmt(c0$estimate, c0$lower, c0$upper), fmt(c1$estimate, c1$lower, c1$upper),
                    fmt(ef$estimate, ef$lower, ef$upper), v))
}
V <- c(V, "", "score_LPA is the off-target contrast. Its adjusted change score is only comparable with score_LDL's relative to its effect on log LDL0, stated above.",
       "Repeat Lp(a) arm: descriptive; no readout rule. Lp(a) reportability drops depend on genotype; see repeat_lpa_counts.csv.")
writeLines(V, file.path(opt$outdir, "verdicts.txt"))
cat(V, sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))
stamp("done -> ", opt$outdir)
