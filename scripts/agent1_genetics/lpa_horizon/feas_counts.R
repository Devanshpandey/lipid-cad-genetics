#!/usr/bin/env Rscript
# ============================================================================
# Lp(a)HORIZON genetic emulation: FEASIBILITY COUNTS ONLY.
#
# Question this answers: does UK Biobank on LS6 hold enough participants with
# established atherosclerotic disease, and enough later events, to test whether
# genetically lower Lp(a) changes the risk of a further event (the trial's
# question; Cho et al., NEJM 2026, HR 1.02 for a 144 nmol/L time-averaged
# difference)?
#
# WHAT THIS SCRIPT MUST NEVER DO
#   Tabulate or model any outcome against LPA genotype, the LPA score or measured
#   Lp(a). It reports (i) cohort sizes and event counts ignoring exposure and
#   (ii) the exposure distribution ignoring outcome. The power table combines the
#   two with a formula. Keep it that way: the analysis plan is to be fixed before
#   any exposure-outcome estimate exists.
#
# DATA LIMIT THAT SHAPES THE DESIGN
#   LS6 has no record-level hospital data. Fields 41270/41280 (ICD-10) and
#   41272/41282 (OPCS-4) hold the FIRST date of each distinct code. A second
#   infarction is therefore visible only if it was coded with a different
#   4-character code (I21.4 after I21.9, or I22.x). Events of a new type, new
#   revascularization codes and all deaths are complete.
#
# Coding sources checked on 2026-11-20 against the UK Biobank showcase
# (codown.cgi): data-coding 6 (field 20002) 1075 heart attack, 1074 angina,
# 1081 stroke, 1082 TIA, 1583 ischaemic stroke, 1067 peripheral vascular disease,
# 1087 leg claudication; data-coding 5 (field 20004) 1070 coronary angioplasty,
# 1095 CABG; data-coding 100605 (field 6150) 1 heart attack, 2 angina, 3 stroke.
# Lp(a) reportability (field 30796, data-coding 4917): 2 and 4 too low, 3 and 5
# too high; reportable range 3.8-189 nmol/L.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse); library(survival) })
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--raw_dir", type = "character"),
  make_option("--covar",   type = "character"),
  make_option("--dose",    type = "character", help = "PLINK --export A file of the LPA cis variants"),
  make_option("--outdir",  type = "character"),
  make_option("--censor",  type = "character", default = "2022-10-31"),
  make_option("--delta",   type = "double",    default = 144, help = "nmol/L, the trial's time-averaged difference"),
  make_option("--threads", type = "integer",   default = 8L)
)))
setDTthreads(opt$threads)
RAW <- opt$raw_dir
CENS <- as.integer(as.IDate(opt$censor))
MIN_CELL <- 5L
H5 <- as.integer(round(5 * 365.25))
LO <- 3.8; HI <- 189
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf(...))
check <- function(cond, msg) if (!isTRUE(cond)) stop("CHECK FAILED: ", msg, call. = FALSE)
ok <- function(x) !is.na(x) & nzchar(x)
days <- function(x) { x[!ok(x)] <- NA_character_; as.integer(as.IDate(x)) }
num <- function(x) { x[!ok(x)] <- NA_character_; suppressWarnings(as.numeric(x)) }

# ---- cohort: the manuscript's QC'd European set with LPA dosages ----------
cv <- fread(opt$covar); cv[, IID := as.character(IID)]
g  <- fread(opt$dose);  g[, IID := as.character(IID)]
COHORT <- intersect(cv$IID, g$IID)
say("[cohort] covariates %d, LPA dosages %d, both %d", nrow(cv), nrow(g), length(COHORT))
check(length(COHORT) > 380000 && length(COHORT) <= nrow(cv), "cohort size outside the expected range")

hdr_of <- function(fid) names(fread(file.path(RAW, sprintf("fid%s.csv", fid)), nrows = 0))
one <- function(fid) {                       # baseline value of a single-valued field
  h <- hdr_of(fid); col <- sprintf("%s-0.0", fid)
  check(col %in% h, sprintf("field %s has no column %s", fid, col))
  d <- fread(file.path(RAW, sprintf("fid%s.csv", fid)), select = c(h[1], col), colClasses = "character")
  setnames(d, c("IID", "v")); d[IID %in% COHORT]
}
arr <- function(fid) {                       # baseline-instance array of a field, long, blanks dropped
  h <- hdr_of(fid); cols <- grep(sprintf("^%s-0\\.[0-9]+$", fid), h, value = TRUE)
  check(length(cols) > 0, sprintf("field %s has no instance-0 columns", fid))
  d <- fread(file.path(RAW, sprintf("fid%s.csv", fid)), select = c(h[1], cols), colClasses = "character")
  setnames(d, 1, "IID"); d <- d[IID %in% COHORT]
  rbindlist(lapply(cols, function(cn) { v <- d[[cn]]; k <- ok(v)
    data.table(IID = d$IID[k], idx = as.integer(sub("^.*\\.", "", cn)), v = v[k]) }))
}
dated <- function(fid_code, fid_date, label) {   # code array joined to its date array
  a <- arr(fid_code); setnames(a, "v", "code")
  b <- arr(fid_date); setnames(b, "v", "date")
  x <- merge(a, b, by = c("IID", "idx"), all.x = TRUE); x[, d := days(date)]
  say("[%s] %d rows in %d people; %.3f%% without a date", label, nrow(x), uniqueN(x$IID), 100 * mean(is.na(x$d)))
  check(mean(is.na(x$d)) < 0.01, paste(label, "codes and dates are misaligned"))
  x[!is.na(d), .(IID, code, d)]
}
single_date <- function(fid) {               # first-occurrence / algorithmic date fields
  x <- one(fid); x[, d := days(v)]
  x[!is.na(d) & d < as.integer(as.IDate("2037-01-01")), .(IID, d)]   # 2037-07-07 = date in the future
}

P <- merge(cv[IID %in% COHORT, .(IID, age, sex)], one("53")[, .(IID, t0 = days(v))], by = "IID")
P <- P[!is.na(t0)]
P <- merge(P, one("40000")[, .(IID, tdeath = days(v))], by = "IID", all.x = TRUE)
P <- merge(P, one("40001")[, .(IID, cause = v)], by = "IID", all.x = TRUE)
P <- merge(P, one("191")[, .(IID, tlost = days(v))], by = "IID", all.x = TRUE)
P[, t_end := pmin(tdeath, tlost, CENS, na.rm = TRUE)]
say("[cohort] with a baseline date %d; deaths before censoring %d (cause recorded %d)",
    nrow(P), P[!is.na(tdeath) & tdeath <= CENS, .N], P[!is.na(tdeath) & tdeath <= CENS & ok(cause), .N])
check(P[!is.na(tdeath) & tdeath <= CENS, .N] > 20000, "implausibly few deaths")
T0 <- P[, .(IID, t0)]

icd <- merge(dated("41270", "41280", "ICD-10"), T0, by = "IID")
ops <- merge(dated("41272", "41282", "OPCS-4"), T0, by = "IID")
fo  <- function(fids) merge(rbindlist(lapply(fids, single_date)), T0, by = "IID")
sr6150 <- arr("6150"); sr20002 <- arr("20002"); sr20004 <- arr("20004")
sr <- function(tab, codes) unique(tab[as.integer(num(v)) %in% codes, IID])

RX_MI  <- "^(I21|I22|I23|I252)"            # I25.2 = old myocardial infarction
RX_MIE <- "^(I21|I22)"                     # outcome: acute or subsequent infarction
RX_IS  <- "^(I63|I64)"
RX_PAD <- "^(I702|I738|I739|I743|I744|I745)"
RX_CHD <- "^I2[0-5]"
RX_TIA <- "^G45"
RX_REV <- "^(K4[0-6]|K49|K50|K75)"
pre  <- function(tab, rx) unique(tab[d <= t0 & grepl(rx, code), IID])
post <- function(tab, rx) tab[d > t0 & d <= CENS & grepl(rx, code), .(t = min(d)), by = IID]

fo_mi  <- fo(c("131298", "131300", "131302")); fo_mie <- fo(c("131298", "131300"))
fo_is  <- fo(c("131366", "131368"));           fo_chd <- fo(c("131296", "131298", "131300", "131302", "131304", "131306"))
al_mi  <- fo("42000"); al_is <- fo("42008"); al_st <- fo("42006")

# ---- disease at baseline -------------------------------------------------
mi_pre  <- unique(c(pre(icd, RX_MI), fo_mi[d <= t0, IID], al_mi[d <= t0, IID], sr(sr6150, 1L), sr(sr20002, 1075L)))
is_pre  <- unique(c(pre(icd, RX_IS), fo_is[d <= t0, IID], al_is[d <= t0, IID], sr(sr20002, 1583L)))
pad_pre <- unique(c(pre(icd, RX_PAD), sr(sr20002, c(1067L, 1087L))))
qual    <- unique(c(mi_pre, is_pre, pad_pre))                    # the trial's three qualifying conditions
broad   <- unique(c(qual, pre(icd, RX_CHD), fo_chd[d <= t0, IID], pre(ops, RX_REV), pre(icd, RX_TIA),
                    al_st[d <= t0, IID], sr(sr6150, c(2L, 3L)), sr(sr20002, c(1074L, 1081L, 1082L)),
                    sr(sr20004, c(1070L, 1095L))))
P[, `:=`(mi_pre = IID %in% mi_pre, is_pre = IID %in% is_pre, pad_pre = IID %in% pad_pre,
         qual = IID %in% qual, broad = IID %in% broad)]
say("[baseline] heart-attack self-report (6150=1) %d; MI %d, ischaemic stroke %d, PAD %d; qualifying %d; any vascular disease %d",
    length(intersect(sr(sr6150, 1L), P$IID)), sum(P$mi_pre), sum(P$is_pre), sum(P$pad_pre), sum(P$qual), sum(P$broad))
check(sum(P$mi_pre) > 5000 && sum(P$mi_pre) < 30000, "prior-MI count implausible")

# ---- events after baseline (first of each type) ---------------------------
first_of <- function(...) rbindlist(list(...))[, .(t = min(t)), by = IID]
add <- function(x, nm) { setnames(x, "t", nm); P <<- merge(P, x, by = "IID", all.x = TRUE); invisible() }
add(first_of(post(icd, RX_MIE), fo_mie[d > t0 & d <= CENS, .(t = min(d)), by = IID],
             al_mi[d > t0 & d <= CENS, .(IID, t = d)]), "t_mi")
add(first_of(post(icd, RX_IS), fo_is[d > t0 & d <= CENS, .(t = min(d)), by = IID],
             al_is[d > t0 & d <= CENS, .(IID, t = d)]), "t_is")
add(post(ops, RX_REV), "t_rev")
P[, t_cvd := fifelse(!is.na(tdeath) & ok(cause) & grepl("^I", cause) & tdeath > t0 & tdeath <= CENS, tdeath, NA_integer_)]
# Events dated after death or loss to follow-up are dropped (follow-up ends there).
for (cn in c("t_mi", "t_is", "t_rev")) P[!is.na(get(cn)) & get(cn) > t_end, (cn) := NA_integer_]

flag <- function(D, cols, horizon = Inf) {
  t <- do.call(pmin, c(as.list(D[, ..cols]), na.rm = TRUE))
  as.integer(!is.na(t) & t > D$t0 & t <= pmin(D$t_end, D$t0 + horizon))
}

# ---- descriptive covariates (not exposures) ------------------------------
st <- unique(rbind(arr("6153"), arr("6177"))[as.integer(num(v)) == 1L, IID])
P[, statin := IID %in% st]
P <- merge(P, one("30780")[, .(IID, ldl = num(v))], by = "IID", all.x = TRUE)

describe <- function(D, label) {
  fu <- (pmin(D$t_end, CENS) - D$t0) / 365.25
  out <- data.table(cohort = label, n = nrow(D), age_median = median(D$age), male_pct = 100 * mean(D$sex == 1),
    statin_pct = 100 * mean(D$statin), ldl_median_mmol = median(D$ldl, na.rm = TRUE),
    followup_median_y = median(fu), person_years = sum(pmax(fu, 0)))
  for (h in c("5y", "full")) {
    hz <- if (h == "5y") H5 else Inf
    set(out, j = paste0("mace3_", h),   value = sum(flag(D, c("t_mi", "t_is", "t_cvd"), hz)))
    set(out, j = paste0("mace4_", h),   value = sum(flag(D, c("t_mi", "t_is", "t_cvd", "t_rev"), hz)))
    set(out, j = paste0("mi_", h),      value = sum(flag(D, "t_mi", hz)))
    set(out, j = paste0("stroke_", h),  value = sum(flag(D, "t_is", hz)))
    set(out, j = paste0("cvdeath_", h), value = sum(flag(D, "t_cvd", hz)))
    set(out, j = paste0("revasc_", h),  value = sum(flag(D, "t_rev", hz)))
  }
  out
}

# ---- second design: people free of disease at baseline who have a first
#      qualifying event in follow-up, then followed from 30 days after it ----
I <- P[broad == FALSE]
I[, index := pmin(t_mi, t_is, na.rm = TRUE)]
I <- I[!is.na(index) & index <= t_end & t_end > index + 30L]
I[, index_type := fifelse(!is.na(t_mi) & t_mi == index, "MI", "IS")]
later <- function(tab, rx, nm) {
  x <- merge(tab[grepl(rx, code)], I[, .(IID, index)], by = "IID")[d > index + 30L & d <= CENS, .(t = min(d)), by = IID]
  setnames(x, "t", nm); x
}
I[, c("t_mi", "t_is", "t_rev", "t0") := NULL]
I <- merge(I, later(icd, RX_MIE, "t_mi"), by = "IID", all.x = TRUE)
I <- merge(I, later(icd, RX_IS, "t_is"), by = "IID", all.x = TRUE)
I <- merge(I, later(ops, RX_REV, "t_rev"), by = "IID", all.x = TRUE)
# first-occurrence fields add the other event type (never seen before the index event)
I <- merge(I, fo_mie[, .(fm = min(d)), by = IID], by = "IID", all.x = TRUE)
I <- merge(I, fo_is[, .(fs = min(d)), by = IID], by = "IID", all.x = TRUE)
I[index_type == "IS" & !is.na(fm) & fm > index + 30L & fm <= CENS, t_mi := pmin(t_mi, fm, na.rm = TRUE)]
I[index_type == "MI" & !is.na(fs) & fs > index + 30L & fs <= CENS, t_is := pmin(t_is, fs, na.rm = TRUE)]
I[, t0 := index + 30L]
I[, t_cvd := fifelse(!is.na(t_cvd) & t_cvd > t0, t_cvd, NA_integer_)]
for (cn in c("t_mi", "t_is", "t_rev")) I[!is.na(get(cn)) & get(cn) > t_end, (cn) := NA_integer_]
I[, age := age + (t0 - P$t0[match(IID, P$IID)]) / 365.25]

counts <- rbindlist(list(
  describe(P[qual == TRUE], "A. established disease at baseline (prior MI, ischaemic stroke or PAD)"),
  describe(P[mi_pre == TRUE], "A1. prior MI"),
  describe(P[is_pre == TRUE], "A2. prior ischaemic stroke"),
  describe(P[pad_pre == TRUE], "A3. peripheral artery disease"),
  describe(P[broad == TRUE], "B. any vascular disease at baseline (adds angina, revascularization, TIA, unspecified stroke)"),
  describe(P[broad == FALSE], "C. primary prevention (none of the above)"),
  describe(I, "D. first MI or ischaemic stroke during follow-up, followed from day 30")), fill = TRUE)

# how much same-type recurrence the first-date-per-code data can see
A <- P[mi_pre == TRUE]
say("[ascertainment] prior MI: %d; of them a later infarction is visible in %d (%.1f%%), a later stroke in %d, CV death in %d",
    nrow(A), sum(flag(A, "t_mi")), 100 * mean(flag(A, "t_mi")), sum(flag(A, "t_is")), sum(flag(A, "t_cvd")))
say("[ascertainment] cohort A crude 5-year 4-point event proportion %.1f%% (trial placebo arm: 12.1%% over a median 4.7 years)",
    100 * counts[1, mace4_5y] / counts[1, n])

# ---- exposure: genotype-predicted Lp(a) in nmol/L, censored-normal model ----
E <- merge(one("30790")[, .(IID, lpa = num(v))], one("30796")[, .(IID, rep = as.integer(num(v)))], by = "IID", all = TRUE)
E[, status := fifelse(rep %in% c(2L, 4L), "below_range", fifelse(rep %in% c(3L, 5L), "above_range",
                 fifelse(!is.na(lpa), "measured", "not_measured")), na = "not_measured")]
E[is.na(status), status := "not_measured"]
E <- E[IID %in% P$IID]
cat("\n=== Lp(a) reportability code by value present (exposure only) ===\n")
print(E[, .N, by = .(reportability = rep, value_present = !is.na(lpa))][order(reportability)])
gcols <- setdiff(names(g), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G <- as.matrix(g[match(P$IID, IID), ..gcols])
for (j in seq_len(ncol(G))) { m <- is.na(G[, j]); if (any(m)) G[m, j] <- mean(G[, j], na.rm = TRUE) }
rownames(G) <- P$IID
fitrows <- E[status != "not_measured"]
lower <- fifelse(fitrows$status == "below_range", NA_real_, fifelse(fitrows$status == "above_range", HI, fitrows$lpa))
upper <- fifelse(fitrows$status == "below_range", LO, fifelse(fitrows$status == "above_range", NA_real_, fitrows$lpa))
Gf <- G[fitrows$IID, , drop = FALSE]
tob <- survreg(Surv(lower, upper, type = "interval2") ~ Gf, dist = "gaussian")
b <- coef(tob); nb <- sum(is.na(b)); b[is.na(b)] <- 0
P[, pred := as.numeric(b[1] + G %*% b[-1])]
P <- merge(P, E[, .(IID, lpa, status)], by = "IID", all.x = TRUE); P[is.na(status), status := "not_measured"]
i72 <- grep("^rs10455872_", gcols)
say("[exposure] censored-normal model on %d people (%d below range, %d above range), %d variants, %d dropped as collinear; residual SD %.1f nmol/L",
    nrow(fitrows), sum(fitrows$status == "below_range"), sum(fitrows$status == "above_range"), ncol(G), nb, tob$scale)
if (length(i72) == 1) say("[exposure] joint coefficient for %s: %.1f nmol/L per counted allele (marginal estimate on file: 109.9 per G allele)",
                          gcols[i72], b[1 + i72])
say("[exposure] sanity (exposure only): mean predicted Lp(a) %.0f below range, %.0f measured, %.0f above range; r2 with measured values %.3f",
    P[status == "below_range", mean(pred)], P[status == "measured", mean(pred)], P[status == "above_range", mean(pred)],
    P[status == "measured", cor(pred, lpa)^2])

expo <- function(D, label) data.table(cohort = label, n = nrow(D),
  measured = sum(D$status == "measured"), below_range = sum(D$status == "below_range"),
  above_range_gt189 = sum(D$status == "above_range"), not_measured = sum(D$status == "not_measured"),
  pred_sd = sd(D$pred), pred_p05 = quantile(D$pred, .05), pred_p50 = quantile(D$pred, .5),
  pred_p95 = quantile(D$pred, .95), pred_p99 = quantile(D$pred, .99), pred_ge150 = sum(D$pred >= 150))
I <- merge(I, P[, .(IID, pred, status)], by = "IID")
exposure <- rbindlist(list(expo(P[qual == TRUE], "A"), expo(P[broad == TRUE], "B"), expo(P[broad == FALSE], "C"), expo(I, "D")))

# ---- power: SE of the log hazard ratio per `delta` nmol/L ~ delta / (SD * sqrt(events)) ----
pw <- rbindlist(lapply(seq_len(4), function(k) {
  cc <- counts[c(1, 5, 6, 7)][k]; sdp <- exposure$pred_sd[k]
  rbindlist(lapply(c("mace3_5y", "mace4_5y", "mace3_full", "mace4_full"), function(e) {
    D <- cc[[e]]; se <- opt$delta / (sdp * sqrt(D))
    data.table(cohort = exposure$cohort[k], endpoint = e, events = D, pred_sd = sdp, se_log_hr = se,
      ci_if_null_lo = exp(-1.96 * se), ci_if_null_hi = exp(1.96 * se), hr_detectable_80pct = exp(-2.80 * se),
      power_vs_0.85 = pnorm(abs(log(0.85)) / se - 1.96), power_vs_0.67 = pnorm(abs(log(0.67)) / se - 1.96))
  }))
}))

mask <- function(x) { for (cn in names(x)) if (is.integer(x[[cn]])) set(x, which(x[[cn]] < MIN_CELL), cn, NA_integer_); x }
fwrite(mask(counts), file.path(opt$outdir, "cohort_counts.csv"))
fwrite(mask(exposure), file.path(opt$outdir, "exposure_summary.csv"))
fwrite(pw, file.path(opt$outdir, "power.csv"))
options(width = 250)
cat("\n=== cohort counts ===\n"); print(t(counts), quote = FALSE)
cat("\n=== exposure summary ===\n"); print(exposure)
cat("\n=== power (hazard ratio per", opt$delta, "nmol/L lower genetically predicted Lp(a)) ===\n"); print(pw, digits = 3)
cat("\nFEASIBILITY_DONE\n")
