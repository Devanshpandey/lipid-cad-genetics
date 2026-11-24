#!/usr/bin/env Rscript
# ============================================================================
# COUNTS ONLY: feasibility of forecasting OCEAN(a)-Outcomes (NCT05581303) and
# ACCLAIM-Lp(a) (NCT06292013) from the UK Biobank emulation cohorts.
#
# Scoping step. No model is fitted and no outcome is tabulated against Lp(a),
# genotype or any score: each row gives a cohort's size, its event counts and,
# separately, how many participants are predicted above the trials' thresholds.
#
# lpa_horizon.R is sourced unchanged for the data assembly and exposure model.
# ============================================================================
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
source(file.path(here, "lpa_horizon.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--raw_dir", type = "character"), make_option("--covar", type = "character"),
  make_option("--dose_lpa", type = "character"), make_option("--w_lpa", type = "character"),
  make_option("--dose_ldl", type = "character"), make_option("--w_ldl", type = "character"),
  make_option("--prs", type = "character"), make_option("--outdir", type = "character"),
  make_option("--censor", type = "character", default = "2022-10-31"),
  make_option("--threads", type = "integer", default = 8L))))
setDTthreads(opt$threads)
CENS <- as.integer(as.IDate(opt$censor))
R <- load_real(opt, CENS)
P <- R$P[complete.cases(R$P[, c(EXPOSURES, ADJ), with = FALSE])]
check(P[qual == TRUE, .N] == FEAS$n_A && P[broad == FALSE, .N] == FEAS$n_C, "cohorts differ from the emulation's")

# prior coronary revascularization by type: OPCS-4 with dates, and self-reported operations
arr0 <- function(fid) {
  f <- file.path(opt$raw_dir, sprintf("fid%s.csv", fid)); h <- names(fread(f, nrows = 0))
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
pci <- unique(c(ops[grepl("^(K49|K50|K75)", code), IID], sr[code == 1070L, IID]))    # 1070 coronary angioplasty +/- stent
cabg <- unique(c(ops[grepl("^K4[0-6]", code), IID], sr[code == 1095L, IID]))         # 1095 coronary artery bypass grafts
P[, `:=`(pci_pre = IID %in% pci, cabg_pre = IID %in% cabg)]

# coronary composite of OCEAN(a)-Outcomes: death from coronary heart disease, infarction, coronary revascularization
P[, t_chd := fifelse(!is.na(tdeath) & ok(cause) & grepl("^I2[0-5]", cause) & tdeath > t0 & tdeath <= CENS, tdeath, NA_integer_)]
ENDPOINTS$cor3 <- c("t_chd", "t_mi", "t_rev")

yes <- function(x) x %in% 1
P[, nrisk := yes(diabetes) + yes(hypertension) + yes(smoking_current)]
cohorts <- list(
  "emulation: established disease (infarction, ischaemic stroke or peripheral disease)" = P$qual,
  "emulation: primary prevention" = !P$broad,
  "OCEAN(a)-like: prior infarction and/or coronary angioplasty" = P$mi_pre | P$pci_pre,
  "OCEAN(a)-like, prior infarction only" = P$mi_pre,
  "OCEAN(a)-like, angioplasty without infarction" = P$pci_pre & !P$mi_pre,
  "ACCLAIM established: prior event or revascularization" = P$qual | P$pci_pre | P$cabg_pre,
  "ACCLAIM at risk: vascular disease recorded, no event or revascularization, age 55+" = P$broad & !(P$qual | P$pci_pre | P$cabg_pre) & P$age >= 55,
  "ACCLAIM at risk: no vascular disease, age 55+" = !P$broad & P$age >= 55,
  "ACCLAIM at risk: no vascular disease, age 55+, diabetes" = !P$broad & P$age >= 55 & yes(P$diabetes),
  "ACCLAIM at risk: no vascular disease, age 55+, two or more of diabetes, hypertension, smoking" = !P$broad & P$age >= 55 & P$nrisk >= 2)

one <- function(nm) {
  D <- P[cohorts[[nm]]]
  ev <- function(e, h = Inf) sum(mk(D, e, h)$ev)
  m <- D$status != "not_measured"
  data.table(cohort = nm, n = nrow(D), age_mean = mean(D$age), male_pct = 100 * mean(D$sex == 1), statin_pct = 100 * mean(D$statin),
    followup_median_y = median(mk(D, "mace4")$time),
    mace4_5y = ev("mace4", H5), mace4_full = ev("mace4"), mace3_5y = ev("mace3", H5), mace3_full = ev("mace3"),
    coronary3_5y = ev("cor3", H5), coronary3_full = ev("cor3"),
    predicted_mean = mean(D$pred), predicted_ge150 = sum(D$pred >= 150), predicted_ge175 = sum(D$pred >= 175), predicted_ge200 = sum(D$pred >= 200),
    measured_above_assay_range_pct = 100 * sum(D$status[m] == "above_range") / sum(m))
}
out <- rbindlist(lapply(names(cohorts), one))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
wcsv(out, file.path(opt$outdir, "forecast_cohort_counts.csv"))
options(width = 250)
cat("\n=== cohort sizes and event counts (no exposure involved) ===\n")
print(out[, .(cohort = substr(cohort, 1, 70), n, age_mean, male_pct, statin_pct, followup_median_y, mace4_5y, mace4_full, mace3_5y, mace3_full, coronary3_5y, coronary3_full)], digits = 3)
cat("\n=== predicted Lp(a) above the trials' thresholds (no outcome involved) ===\n")
print(out[, .(cohort = substr(cohort, 1, 70), n, predicted_mean, predicted_ge150, predicted_ge175, predicted_ge200, measured_above_assay_range_pct)], digits = 3)
cat(sprintf("\nprior coronary angioplasty %d, prior bypass grafting %d, either %d\n", sum(P$pci_pre), sum(P$cabg_pre), sum(P$pci_pre | P$cabg_pre)))
cat("\nFEAS_DONE\n")
