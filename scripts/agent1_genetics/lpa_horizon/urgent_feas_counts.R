#!/usr/bin/env Rscript
# ============================================================================
# COUNTS ONLY: can urgent coronary revascularization be approximated from the
# hospital data held on LS6?
#
# Scoping step for a reviewer's request. Lp(a)HORIZON counted urgent coronary
# revascularization leading to hospitalization; the genetic comparison counts any
# coronary revascularization. Record-level episodes with the method of admission
# are not held. What is held is the FIRST date of each distinct ICD-10 and OPCS-4
# code. Proxy examined here: a revascularization is "acute-linked" if an acute
# coronary code (unstable angina I20.0, infarction I21 or I22, other acute
# ischaemic heart disease I24) was first recorded from W days before the procedure
# to one day after it. A repeat acute event under a code already on record is
# invisible, so the proxy undercounts.
#
# No model is fitted and nothing is tabulated against Lp(a), genotype or any
# score. Output: event counts by cohort, the lag from acute code to first
# procedure, and how the four-point composite divides by type of first event.
# lpa_horizon.R is sourced unchanged for the data assembly.
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
cc <- function(D) D[complete.cases(D[, c(EXPOSURES, ADJ), with = FALSE])]
P <- cc(R$P); I <- cc(R$I)
check(P[qual == TRUE, .N] == FEAS$n_A && P[broad == FALSE, .N] == FEAS$n_C, "cohorts differ from the comparison's")

# every dated ICD-10 and OPCS-4 code (first date of each distinct code)
arr0 <- function(fid) {
  f <- file.path(opt$raw_dir, sprintf("fid%s.csv", fid)); h <- names(fread(f, nrows = 0))
  cols <- grep(sprintf("^%s-0\\.[0-9]+$", fid), h, value = TRUE)
  check(length(cols) > 0, sprintf("field %s has no instance-0 columns", fid))
  d <- fread(f, select = c(h[1], cols), colClasses = "character"); setnames(d, 1, "IID"); d <- d[IID %in% P$IID]
  rbindlist(lapply(cols, function(cn) { v <- d[[cn]]; k <- ok(v)
    data.table(IID = d$IID[k], idx = rep(as.integer(sub("^.*\\.", "", cn)), sum(k)), v = v[k]) }))
}
dated0 <- function(fc, fd) {
  a <- arr0(fc); setnames(a, "v", "code"); b <- arr0(fd); setnames(b, "v", "date")
  x <- merge(a, b, by = c("IID", "idx")); x[, d := days(date)]
  x[!is.na(d), .(IID, code, d)]
}
icd <- dated0("41270", "41280"); ops <- dated0("41272", "41282")
ACS <- "^(I200|I21|I22|I24)"
acs <- unique(icd[grepl(ACS, code), .(IID, acode = code, da = d)])
rev <- unique(ops[grepl(RX$REV, code), .(IID, pcode = code, dp = d)])
say("dated codes: %d acute coronary codes in %d participants; %d revascularization codes in %d participants",
    nrow(acs), uniqueN(acs$IID), nrow(rev), uniqueN(rev$IID))

# for each procedure, the acute code first recorded closest before it (lag = days from code to procedure; -1 allows the day after)
lk <- merge(rev, acs, by = "IID", allow.cartesian = TRUE)
lk[, lag := dp - da]; lk <- lk[lag >= -1L]
setorder(lk, IID, pcode, dp, lag)
near <- lk[, .SD[1L], by = .(IID, pcode, dp)]
rev <- merge(rev, near[, .(IID, pcode, dp, lag, acode)], by = c("IID", "pcode", "dp"), all.x = TRUE)
rev[, family := fifelse(is.na(acode), NA_character_, fifelse(grepl("^(I21|I22)", acode), "infarction",
                 fifelse(grepl("^I200", acode), "unstable angina", "other acute")))]

WINDOWS <- c(7L, 14L, 30L)
ENDPOINTS$urev7 <- "t_u7"; ENDPOINTS$urev14 <- "t_u14"; ENDPOINTS$urev30 <- "t_u30"; ENDPOINTS$elrev <- "t_el"
ENDPOINTS$mace4u7 <- c("t_mi", "t_is", "t_cvd", "t_u7"); ENDPOINTS$mace4u30 <- c("t_mi", "t_is", "t_cvd", "t_u30")
cohorts <- list("established disease at baseline (cohort A)" = P[qual == TRUE],
                "before a first event (cohort C)" = P[broad == FALSE],
                "from 30 days after a first event in follow-up (cohort D)" = I)

build <- function(D) {
  D <- copy(D)
  r <- merge(rev, D[, .(IID, t0, t_end)], by = "IID")[dp > t0 & dp <= pmin(t_end, CENS)]
  frst <- function(x) { y <- x[, .(t = min(dp)), by = IID]; y$t[match(D$IID, y$IID)] }
  # the first revascularization rebuilt here must be the endpoint the comparison used
  check(identical(as.integer(frst(r)), as.integer(D$t_rev)), "revascularization dates differ from the comparison's endpoint")
  for (W in WINDOWS) set(D, j = paste0("t_u", W), value = as.integer(frst(r[!is.na(lag) & lag <= W])))
  set(D, j = "t_el", value = as.integer(frst(r[is.na(lag) | lag > max(WINDOWS)])))
  for (fam in c("infarction", "unstable angina", "other acute"))
    set(D, j = paste0("t_u7_", gsub(" ", "_", fam)), value = as.integer(frst(r[!is.na(lag) & lag <= 7L & family == fam])))
  list(D = D, r = r)
}
B <- lapply(cohorts, build)

counts <- rbindlist(lapply(names(B), function(nm) {
  D <- B[[nm]]$D
  ev <- function(e, h = Inf) sum(mk(D, e, h)$ev)
  inwin <- function(t) sum(!is.na(t) & t > D$t0 & t <= D$t_end)
  pre <- function(tab, col) 100 * mean(D$IID %in% unique(merge(tab, D[, .(IID, t0)], by = "IID")[get(col) <= t0, IID]))
  data.table(cohort = nm, n = nrow(D),
    revasc_any = ev("revasc"), revasc_any_5y = ev("revasc", H5),
    acute_linked_7d = ev("urev7"), acute_linked_14d = ev("urev14"), acute_linked_30d = ev("urev30"), acute_linked_7d_5y = ev("urev7", H5),
    acute7_infarction = inwin(D$t_u7_infarction), acute7_unstable_angina = inwin(D$t_u7_unstable_angina), acute7_other = inwin(D$t_u7_other_acute),
    not_acute_linked_30d = ev("elrev"),
    three_point = ev("mace3"), four_point_any_revasc = ev("mace4"), four_point_acute7 = ev("mace4u7"), four_point_acute30 = ev("mace4u30"),
    three_point_5y = ev("mace3", H5), four_point_any_revasc_5y = ev("mace4", H5), four_point_acute7_5y = ev("mace4u7", H5),
    prior_revasc_code_pct = pre(rev, "dp"), prior_acute_code_pct = pre(acs, "da"))
}))

lags <- rbindlist(lapply(names(B), function(nm) {
  D <- B[[nm]]$D; r <- B[[nm]]$r
  f <- merge(r, D[!is.na(t_rev), .(IID, t_rev)], by = "IID")[dp == t_rev]
  setorder(f, IID, lag, na.last = TRUE); f <- f[, .SD[1L], by = IID]
  f[, bin := fifelse(is.na(lag), "8 no acute code on record before it", fifelse(lag <= 0L, "1 same day or the day after",
            fifelse(lag <= 2L, "2 one or two days before", fifelse(lag <= 7L, "3 three to 7 days", fifelse(lag <= 14L, "4 eight to 14 days",
            fifelse(lag <= 30L, "5 fifteen to 30 days", fifelse(lag <= 365L, "6 thirty-one to 365 days", "7 more than a year")))))))]
  x <- f[, .(first_procedures = .N), by = bin][order(bin)]
  x[, `:=`(cohort = nm, pct = 100 * first_procedures / sum(first_procedures))]
  x[, .(cohort, lag_from_nearest_acute_code = bin, first_procedures, pct)]
}))

types <- rbindlist(lapply(names(B), function(nm) {
  D <- B[[nm]]$D; r <- B[[nm]]$r
  inw <- function(t) fifelse(!is.na(t) & t > D$t0 & t <= D$t_end, t, NA_integer_)
  tm <- inw(D$t_mi); ti <- inw(D$t_is); tc <- inw(D$t_cvd); tr <- inw(D$t_rev)
  first <- pmin(tm, ti, tc, tr, na.rm = TRUE); ev <- !is.na(first)
  check(sum(ev) == sum(mk(D, "mace4")$ev), "four-point events differ from the comparison's")
  is_ <- function(t) ev & !is.na(t) & t == first
  death <- is_(tc); mi <- is_(tm) & !death; st <- is_(ti) & !death & !mi; alone <- ev & !(death | mi | st)
  fr <- merge(r, data.table(IID = D$IID[alone], first = first[alone]), by = "IID")[dp == first]
  setorder(fr, IID, lag, na.last = TRUE); fr <- fr[, .SD[1L], by = IID]
  check(nrow(fr) == sum(alone), "a first-event revascularization could not be matched to its procedure")
  data.table(cohort = nm, four_point_events = sum(ev), first_cv_death = sum(death), first_infarction = sum(mi), first_ischaemic_stroke = sum(st),
    first_revascularization_alone = sum(alone), alone_acute_linked_7d = fr[!is.na(lag) & lag <= 7L, .N],
    alone_acute_linked_8_to_30d = fr[!is.na(lag) & lag > 7L & lag <= 30L, .N], alone_not_acute_linked = fr[is.na(lag) | lag > 30L, .N],
    alone_pct_of_four_point = 100 * sum(alone) / sum(ev))
}))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
wcsv(counts, file.path(opt$outdir, "urgent_proxy_counts.csv"))
wcsv(lags, file.path(opt$outdir, "lag_distribution.csv"))
wcsv(types, file.path(opt$outdir, "first_event_types.csv"))
options(width = 250)
cat("\n=== revascularization events, any and acute-linked (no exposure involved) ===\n")
print(counts[, .(cohort = substr(cohort, 1, 46), n, revasc_any, acute_linked_7d, acute_linked_14d, acute_linked_30d, not_acute_linked_30d, revasc_any_5y, acute_linked_7d_5y)])
cat("\n=== acute-linked within 7 days, by the acute code; and how complete the link can be ===\n")
print(counts[, .(cohort = substr(cohort, 1, 46), acute7_infarction, acute7_unstable_angina, acute7_other, prior_revasc_code_pct, prior_acute_code_pct)], digits = 3)
cat("\n=== composites ===\n")
print(counts[, .(cohort = substr(cohort, 1, 46), three_point, four_point_any_revasc, four_point_acute7, four_point_acute30, three_point_5y, four_point_any_revasc_5y, four_point_acute7_5y)])
cat("\n=== lag from the nearest earlier acute code to the first revascularization ===\n")
print(lags, digits = 3)
cat("\n=== first event of the four-point composite, by type ===\n")
print(types, digits = 3)
cat("\nURGENT_FEAS_DONE\n")
