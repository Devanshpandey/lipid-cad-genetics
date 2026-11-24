#!/usr/bin/env Rscript
# ============================================================================
# Genetic emulation of Lp(a)HORIZON in UK Biobank.
#
# Question: genetically higher Lp(a) raises the risk of a FIRST vascular event.
# Does it still raise the risk of a further event once disease is established,
# which is the population the trial treated (Cho et al., NEJM 2026: HR 1.02,
# 0.90-1.15, for a 144 nmol/L time-averaged reduction)?
#
# The rules for reading the output are in PREREGISTERED_READOUT.md in this folder
# and were committed before any exposure-outcome estimate from this script
# existed. The verdict() function below is that document, nothing looser.
#
# Steps
#   --step selftest  synthetic data only: unit tests of the data assembly on a
#                    fixture, censored-model and score checks, planted-effect
#                    recovery, every verdict branch, and a full smoke run
#   --step gates     real data, no new exposure-outcome estimate: assembly counts,
#                    exposure model, reproduction of the published score result
#   --step test      the gates again, then every estimate and the verdict
#
# Failure behaviour: a failed gate exits with status 2 and writes gates.csv only.
# Outputs are summary level; counts below 5 are suppressed.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse); library(survival) })
set.seed(42)

DELTA <- 144                 # nmol/L, the trial's time-averaged difference
LO <- 3.8; HI <- 189         # reportable range of the UK Biobank Lp(a) assay
H5 <- as.integer(round(5 * 365.25))
MIN_CELL <- 5L
PUB <- list(n = 379896L, events = 32053L, hr = c(1.096, 1.083, 1.110))  # manuscript Table oS34, Lp(a) score
TRIAL <- c(hr = 1.02, lo = 0.90, hi = 1.15)
FEAS <- list(n_A = 16955L, ev_A = 4790L, n_C = 371735L)   # job 3600440, same definitions
EXPOSURES <- c("x_lpa", "x_ldl", "x_prs")
ADJ <- c("age", "sex", paste0("PC", 1:10))
RISK <- c("bmi", "smoking_current", "diabetes", "hypertension")
ENDPOINTS <- list(mace4 = c("t_mi", "t_is", "t_cvd", "t_rev"), mace3 = c("t_mi", "t_is", "t_cvd"),
                  mi = "t_mi", stroke = "t_is", cvdeath = "t_cvd", revasc = "t_rev")
RX <- list(MI = "^(I21|I22|I23|I252)", MIE = "^(I21|I22)", IS = "^(I63|I64)",
           PAD = "^(I702|I738|I739|I743|I744|I745)", CHD = "^I2[0-5]", TIA = "^G45",
           REV = "^(K4[0-6]|K49|K50|K75)")
META <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")

say <- function(...) message(sprintf(...))
check <- function(cond, msg) if (!isTRUE(cond)) stop("CHECK FAILED: ", msg, call. = FALSE)
ok <- function(x) !is.na(x) & nzchar(x)
days <- function(x) { x[!ok(x)] <- NA_character_; as.integer(as.IDate(x)) }
num <- function(x) { x[!ok(x)] <- NA_character_; suppressWarnings(as.numeric(x)) }
zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)

# ---------------------------------------------------------------------------
# Data assembly. LS6 holds the FIRST date of each distinct hospital code only
# (41270/41280, 41272/41282), so a later event of the same type is visible only
# under a different 4-character code. New event types and deaths are complete.
# ---------------------------------------------------------------------------
build_cohort <- function(RAW, cv, ids, CENS) {
  hdr_of <- function(fid) names(fread(file.path(RAW, sprintf("fid%s.csv", fid)), nrows = 0))
  one <- function(fid) {
    h <- hdr_of(fid); col <- sprintf("%s-0.0", fid)
    check(col %in% h, sprintf("field %s has no column %s", fid, col))
    d <- fread(file.path(RAW, sprintf("fid%s.csv", fid)), select = c(h[1], col), colClasses = "character")
    setnames(d, c("IID", "v")); d[IID %in% ids]
  }
  arr <- function(fid) {
    h <- hdr_of(fid); cols <- grep(sprintf("^%s-0\\.[0-9]+$", fid), h, value = TRUE)
    check(length(cols) > 0, sprintf("field %s has no instance-0 columns", fid))
    d <- fread(file.path(RAW, sprintf("fid%s.csv", fid)), select = c(h[1], cols), colClasses = "character")
    setnames(d, 1, "IID"); d <- d[IID %in% ids]
    rbindlist(lapply(cols, function(cn) { v <- d[[cn]]; k <- ok(v)
      data.table(IID = d$IID[k], idx = rep(as.integer(sub("^.*\\.", "", cn)), sum(k)), v = v[k]) }))
  }
  dated <- function(fc, fd, label) {
    a <- arr(fc); setnames(a, "v", "code"); b <- arr(fd); setnames(b, "v", "date")
    x <- merge(a, b, by = c("IID", "idx"), all.x = TRUE); x[, d := days(date)]
    check(nrow(x) == 0 || mean(is.na(x$d)) < 0.01, paste(label, "codes and dates are misaligned"))
    x[!is.na(d), .(IID, code, d)]
  }
  single_date <- function(fid) {
    x <- one(fid); x[, d := days(v)]
    x[!is.na(d) & d < as.integer(as.IDate("2037-01-01")), .(IID, d)]   # 2037-07-07 codes a future date
  }

  P <- merge(cv[IID %in% ids], one("53")[, .(IID, t0 = days(v))], by = "IID")
  n_in <- nrow(P); P <- P[!is.na(t0)]
  P <- merge(P, one("40000")[, .(IID, tdeath = days(v))], by = "IID", all.x = TRUE)
  P <- merge(P, one("40001")[, .(IID, cause = v)], by = "IID", all.x = TRUE)
  P <- merge(P, one("191")[, .(IID, tlost = days(v))], by = "IID", all.x = TRUE)
  P[, t_end := pmin(tdeath, tlost, CENS, na.rm = TRUE)]
  T0 <- P[, .(IID, t0)]
  icd <- merge(dated("41270", "41280", "ICD-10"), T0, by = "IID")
  ops <- merge(dated("41272", "41282", "OPCS-4"), T0, by = "IID")
  fo <- function(fids) merge(rbindlist(lapply(fids, single_date)), T0, by = "IID")
  sr6150 <- arr("6150"); sr20002 <- arr("20002"); sr20004 <- arr("20004")
  sr <- function(tab, codes) unique(tab[as.integer(num(v)) %in% codes, IID])
  pre  <- function(tab, rx) unique(tab[d <= t0 & grepl(rx, code), IID])
  post <- function(tab, rx) tab[d > t0 & d <= CENS & grepl(rx, code), .(t = min(d)), by = IID]
  fo_mi <- fo(c("131298", "131300", "131302")); fo_mie <- fo(c("131298", "131300"))
  fo_is <- fo(c("131366", "131368"))
  fo_chd <- fo(c("131296", "131298", "131300", "131302", "131304", "131306"))
  al_mi <- fo("42000"); al_is <- fo("42008"); al_st <- fo("42006")

  mi_pre  <- unique(c(pre(icd, RX$MI), fo_mi[d <= t0, IID], al_mi[d <= t0, IID], sr(sr6150, 1L), sr(sr20002, 1075L)))
  is_pre  <- unique(c(pre(icd, RX$IS), fo_is[d <= t0, IID], al_is[d <= t0, IID], sr(sr20002, 1583L)))
  pad_pre <- unique(c(pre(icd, RX$PAD), sr(sr20002, c(1067L, 1087L))))
  qual    <- unique(c(mi_pre, is_pre, pad_pre))
  broad   <- unique(c(qual, pre(icd, RX$CHD), fo_chd[d <= t0, IID], pre(ops, RX$REV), pre(icd, RX$TIA),
                      al_st[d <= t0, IID], sr(sr6150, c(2L, 3L)), sr(sr20002, c(1074L, 1081L, 1082L)),
                      sr(sr20004, c(1070L, 1095L))))
  P[, `:=`(mi_pre = IID %in% mi_pre, is_pre = IID %in% is_pre, pad_pre = IID %in% pad_pre,
           qual = IID %in% qual, broad = IID %in% broad)]

  first_of <- function(...) rbindlist(list(...))[, .(t = min(t)), by = IID]
  add <- function(P, x, nm) { setnames(x, "t", nm); merge(P, x, by = "IID", all.x = TRUE) }
  P <- add(P, first_of(post(icd, RX$MIE), fo_mie[d > t0 & d <= CENS, .(t = min(d)), by = IID],
                       al_mi[d > t0 & d <= CENS, .(IID, t = d)]), "t_mi")
  P <- add(P, first_of(post(icd, RX$IS), fo_is[d > t0 & d <= CENS, .(t = min(d)), by = IID],
                       al_is[d > t0 & d <= CENS, .(IID, t = d)]), "t_is")
  P <- add(P, post(ops, RX$REV), "t_rev")
  P[, t_cvd := fifelse(!is.na(tdeath) & ok(cause) & grepl("^I", cause) & tdeath > t0 & tdeath <= CENS, tdeath, NA_integer_)]
  for (cn in c("t_mi", "t_is", "t_rev")) P[!is.na(get(cn)) & get(cn) > t_end, (cn) := NA_integer_]

  # baseline descriptors; `statin_first` and `tcad` exist only to reproduce the published model
  s <- rbind(arr("6153"), arr("6177")); s0 <- s[idx == 0L]
  P[, statin := IID %in% unique(s[as.integer(num(v)) == 1L, IID])]
  answered <- unique(s0$IID); first_statin <- unique(s0[as.integer(num(v)) == 1L, IID])
  P[, statin_first := fifelse(IID %in% answered, as.integer(IID %in% first_statin), NA_integer_)]
  P <- merge(P, fo_chd[, .(tcad = min(d)), by = IID], by = "IID", all.x = TRUE)
  P <- merge(P, one("30780")[, .(IID, ldl = num(v))], by = "IID", all.x = TRUE)
  E <- merge(one("30790")[, .(IID, lpa = num(v))], one("30796")[, .(IID, rep = as.integer(num(v)))], by = "IID", all = TRUE)
  E[, status := fifelse(rep %in% c(2L, 4L), "below_range", fifelse(rep %in% c(3L, 5L), "above_range",
                   fifelse(!is.na(lpa), "measured", "not_measured")))]
  P <- merge(P, E[, .(IID, lpa, status)], by = "IID", all.x = TRUE); P[is.na(status), status := "not_measured"]

  # second established-disease cohort: free of vascular disease at baseline, first MI or
  # ischaemic stroke during follow-up, followed from 30 days after it
  I <- P[broad == FALSE]
  I[, index := pmin(t_mi, t_is, na.rm = TRUE)]
  I <- I[!is.na(index) & index <= t_end & t_end > index + 30L]
  I[, index_type := fifelse(!is.na(t_mi) & t_mi == index, "MI", "IS")]
  ix <- I[, .(IID, index)]
  later <- function(tab, rx, nm) {
    x <- merge(tab[grepl(rx, code)], ix, by = "IID")[d > index + 30L & d <= CENS, .(t = min(d)), by = IID]
    setnames(x, "t", nm); x
  }
  I[, c("t_mi", "t_is", "t_rev", "t0") := NULL]
  I <- merge(I, later(icd, RX$MIE, "t_mi"), by = "IID", all.x = TRUE)
  I <- merge(I, later(icd, RX$IS, "t_is"), by = "IID", all.x = TRUE)
  I <- merge(I, later(ops, RX$REV, "t_rev"), by = "IID", all.x = TRUE)
  I <- merge(I, fo_mie[, .(fm = min(d)), by = IID], by = "IID", all.x = TRUE)
  I <- merge(I, fo_is[, .(fs = min(d)), by = IID], by = "IID", all.x = TRUE)
  I[index_type == "IS" & !is.na(fm) & fm > index + 30L & fm <= CENS, t_mi := pmin(t_mi, fm, na.rm = TRUE)]
  I[index_type == "MI" & !is.na(fs) & fs > index + 30L & fs <= CENS, t_is := pmin(t_is, fs, na.rm = TRUE)]
  I[, t0 := index + 30L]
  I[, t_cvd := fifelse(!is.na(t_cvd) & t_cvd > t0, t_cvd, NA_integer_)]
  for (cn in c("t_mi", "t_is", "t_rev")) I[!is.na(get(cn)) & get(cn) > t_end, (cn) := NA_integer_]
  I[, age := age + (t0 - P$t0[match(IID, P$IID)]) / 365.25]
  I[, c("fm", "fs") := NULL]
  list(P = P, I = I, n_in = n_in)
}

# Paper's score: explicit allele orientation, mean imputation (statin_interaction_incident.R)
build_score <- function(g, w) {
  dc <- setdiff(names(g), META)
  map <- merge(data.table(col = dc, SNP = sub("_[^_]+$", "", dc), counted = sub("^.*_", "", dc)), w, by = "SNP")
  sc <- numeric(nrow(g))
  for (i in seq_len(nrow(map))) {
    d <- g[[map$col[i]]]
    if (toupper(map$counted[i]) != toupper(map$EA[i])) d <- 2 - d
    d[is.na(d)] <- mean(d, na.rm = TRUE); sc <- sc + d * map$BETA[i]
  }
  list(s = sc, n_variants = nrow(map))
}

# Genotype-predicted Lp(a) in nmol/L: censored-normal (Tobit) model, left-censored at 3.8
# and right-censored at 189, fitted jointly on the cis variants
fit_exposure <- function(lpa, status, G) {
  k <- status != "not_measured"
  lower <- ifelse(status[k] == "below_range", NA, ifelse(status[k] == "above_range", HI, lpa[k]))
  upper <- ifelse(status[k] == "below_range", LO, ifelse(status[k] == "above_range", NA, lpa[k]))
  Gf <- G[k, , drop = FALSE]
  m <- survreg(Surv(lower, upper, type = "interval2") ~ Gf, dist = "gaussian")
  b <- coef(m); dropped <- sum(is.na(b)); b[is.na(b)] <- 0
  pred <- as.numeric(b[1] + G %*% b[-1])
  ms <- status == "measured"
  list(pred = pred, coef = setNames(b[-1], colnames(G)), scale = m$scale, n = sum(k), dropped = dropped,
       r2 = cor(pred[ms], lpa[ms])^2)
}

finish_exposures <- function(P, I) {
  P[, x_lpa := pred / DELTA]            # per 144 nmol/L HIGHER predicted Lp(a)
  P[, x_ldl := zs(score_ldl_raw)]       # per SD of the paper's LDL-C score
  P[, x_prs := zs(prs_n)]               # per SD of the UK Biobank-free non-lipid CAD score
  P[, x_lpa_paper := zs(score_lpa_raw)] # per SD of the paper's Lp(a) score
  if (!is.null(I)) {
    keep <- c("IID", "pred", "x_lpa", "x_ldl", "x_prs", "x_lpa_paper")
    I <- merge(I[, setdiff(names(I), keep[-1]), with = FALSE], P[, ..keep], by = "IID")
  }
  list(P = P, I = I)
}

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
mk <- function(D, endpoint, horizon = Inf) {
  cols <- ENDPOINTS[[endpoint]]
  t <- do.call(pmin, c(unname(as.list(D[, ..cols])), na.rm = TRUE))
  lim <- pmin(D$t_end, D$t0 + horizon)
  ev <- as.integer(!is.na(t) & t > D$t0 & t <= lim)
  data.table(time = (ifelse(ev == 1L, t, lim) - D$t0) / 365.25, ev = ev)
}
fit <- function(D, endpoint, horizon = Inf, xs = EXPOSURES, adj = ADJ, w = NULL, zph = FALSE) {
  X <- cbind(D[, c(xs, adj), with = FALSE], mk(D, endpoint, horizon))
  X[, wt := if (is.null(w)) 1 else w]
  X <- X[time > 0 & wt > 0]; X <- X[complete.cases(X)]
  if (sum(X$ev) < 20L) return(NULL)
  f <- as.formula(paste("Surv(time, ev) ~", paste(c(xs, adj), collapse = " + ")))
  m <- if (is.null(w)) coxph(f, data = X) else coxph(f, data = X, weights = wt, robust = TRUE)
  b <- coef(m)[xs]; if (anyNA(b)) return(NULL)
  out <- list(b = b, V = vcov(m)[xs, xs, drop = FALSE], n = m$n, events = sum(X$ev), all = coef(m))
  if (zph) out$zph <- tryCatch(cox.zph(m, terms = FALSE)$table[xs, "p"], error = function(e) rep(NA_real_, length(xs)))
  out
}
# weights that give the primary-prevention cohort the age-by-sex mix of the established-disease cohort
std_weights <- function(Cd, Ad) {
  cell <- function(D) paste(cut(D$age, c(-Inf, 44, 49, 54, 59, 64, Inf)), D$sex)
  share <- function(D) { t <- prop.table(table(cell(D))); setNames(as.numeric(t), names(t)) }
  pa <- share(Ad); pc <- share(Cd)
  w <- unname(pa[cell(Cd)] / pc[cell(Cd)]); w[is.na(w)] <- 0
  w / mean(w)
}
ratio_stats <- function(A, C, x) {
  bA <- A$b[[x]]; bC <- C$b[[x]]; vA <- A$V[x, x]; vC <- C$V[x, x]
  d <- bA - bC; sed <- sqrt(vA + vC); rho <- bA / bC; ser <- sqrt(vA / bC^2 + bA^2 * vC / bC^4)
  data.table(exposure = x, b_established = bA, se_established = sqrt(vA), b_primary = bC, se_primary = sqrt(vC),
             p_present = pnorm(bA / sqrt(vA), lower.tail = FALSE), p_anchor = pnorm(bC / sqrt(vC), lower.tail = FALSE),
             diff = d, se_diff = sed, p_diff = 2 * pnorm(-abs(d / sed)),
             ratio = rho, ratio_lo = rho - 1.96 * ser, ratio_hi = rho + 1.96 * ser)
}
# one-sided test that the Lp(a) ratio is smaller than the LDL-C ratio (delta method; the two
# cohorts are disjoint, so their covariance blocks are independent)
spec_test <- function(A, C, x1 = "x_lpa", x2 = "x_ldl") {
  th <- c(A$b[[x1]], A$b[[x2]], C$b[[x1]], C$b[[x2]])
  V <- matrix(0, 4, 4); V[1:2, 1:2] <- A$V[c(x1, x2), c(x1, x2)]; V[3:4, 3:4] <- C$V[c(x1, x2), c(x1, x2)]
  g <- th[1] / th[3] - th[2] / th[4]
  gr <- c(1 / th[3], -1 / th[4], -th[1] / th[3]^2, th[2] / th[4]^2)
  se <- sqrt(drop(t(gr) %*% V %*% gr))
  list(diff = g, se = se, p_one = pnorm(g / se))
}
# PREREGISTERED_READOUT.md, section "Verdict"
verdict <- function(rs, sp) {
  L <- rs[exposure == "x_lpa"]; D <- rs[exposure == "x_ldl"]
  anchors <- L$p_anchor < 0.001 && D$p_anchor < 0.001
  P   <- L$p_present < 0.025
  Att <- L$p_diff < 0.05 && L$b_established < L$b_primary
  C1  <- D$p_present < 0.025
  S   <- sp$p_one < 0.05
  v <- if (!anchors) "NO_ANCHOR"
       else if (P && !Att) "GENETICS_DOES_NOT_REPRODUCE_TRIAL"
       else if (P && Att)  "PARTIAL_PERSISTENCE"
       else if (!Att)      "INDETERMINATE"
       else if (!C1)       "UNINFORMATIVE_DESIGN"
       else if (S)         "STAGE_SPECIFIC_FOR_LPA"
       else                "ATTENUATED_NOT_SPECIFIC"
  data.table(verdict = v, anchors_ok = anchors, lpa_present = P, lpa_attenuated = Att, ldl_control_present = C1,
             lpa_more_attenuated_than_ldl = S, specificity_diff = sp$diff, specificity_se = sp$se, specificity_p = sp$p_one)
}
primary_analysis <- function(P) {
  A <- P[qual == TRUE]; Cc <- P[broad == FALSE]
  w <- std_weights(Cc, A)
  fA <- fit(A, "mace4", zph = TRUE); fC <- fit(Cc, "mace4", w = w)
  check(!is.null(fA) && !is.null(fC), "primary models could not be fitted")
  rs <- rbindlist(lapply(EXPOSURES, function(x) ratio_stats(fA, fC, x)))
  sp <- spec_test(fA, fC)
  list(rs = rs, verdict = verdict(rs, sp), fA = fA, fC = fC, w = w)
}

est_rows <- function(f, cohort, endpoint, horizon, model) {
  if (is.null(f)) return(NULL)
  data.table(cohort = cohort, endpoint = endpoint, horizon = horizon, model = model, exposure = names(f$b),
             n = as.integer(f$n), events = as.integer(f$events), b = unname(f$b), se = unname(sqrt(diag(f$V))))
}
bias_budget <- function(A, Cc) {
  fa <- fit(A, "mace4", adj = c(ADJ, RISK))
  if (is.null(fa)) return(NULL)
  rbindlist(lapply(c("x_ldl", "x_prs", "age", "sex", RISK), function(k) {
    sA <- coef(lm(A[[k]] ~ A$x_lpa))[[2]]; sC <- coef(lm(Cc[[k]] ~ Cc$x_lpa))[[2]]
    data.table(factor = k, in_primary_model = !(k %in% RISK), slope_established = sA, slope_primary = sC,
               induced = sA - sC, log_hr_of_factor = fa$all[[k]], implied_bias_log_hr = (sA - sC) * fa$all[[k]])
  }))
}
published_control <- function(P, CENS) {
  d <- P[!is.na(statin_first) & !is.na(score_lpa_raw) & !is.na(score_ldl_raw)]
  n_start <- nrow(d); prev <- d[!is.na(tcad) & tcad <= t0, .N]
  d <- d[is.na(tcad) | tcad > t0]
  d[, tend := pmin(fifelse(is.na(tcad), CENS, tcad), fifelse(is.na(tdeath), CENS, tdeath), CENS)]
  d[, event := as.integer(!is.na(tcad) & tcad <= tend)]; d[, fu := (tend - t0) / 365.25]; d <- d[fu > 0]
  d[, sL := zs(score_lpa_raw)]
  f <- as.formula(paste("Surv(fu, event) ~ sL * statin_first + age + sex +", paste(paste0("PC", 1:10), collapse = " + ")))
  m <- coxph(f, data = d); cf <- summary(m)$coefficients["sL", ]
  list(n_start = n_start, prevalent = prev, n = as.integer(m$n), events = as.integer(m$nevent),
       hr = exp(cf[[1]]), lo = exp(cf[[1]] - 1.96 * cf[[3]]), hi = exp(cf[[1]] + 1.96 * cf[[3]]))
}
wcsv <- function(x, f) {
  x <- copy(x)
  for (cn in names(x)) if (is.integer(x[[cn]])) set(x, which(x[[cn]] < MIN_CELL), cn, NA_integer_)
  fwrite(x, f)
}

run_analysis <- function(P, I, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  n0 <- nrow(P); P <- P[complete.cases(P[, c(EXPOSURES, "x_lpa_paper", ADJ), with = FALSE])]
  say("[analysis] %d of %d participants have every exposure and covariate", nrow(P), n0)
  if (!is.null(I)) I <- I[IID %in% P$IID]
  pa <- primary_analysis(P)
  A <- P[qual == TRUE]; Cc <- P[broad == FALSE]; w <- pa$w

  cohorts <- list(A = list(A, NULL), A1_prior_MI = list(P[mi_pre == TRUE], NULL), A2_prior_stroke = list(P[is_pre == TRUE], NULL),
                  A3_PAD = list(P[pad_pre == TRUE], NULL), B = list(P[broad == TRUE], NULL),
                  C = list(Cc, NULL), C_std = list(Cc, w))
  if (!is.null(I) && nrow(I) > 0) cohorts$D <- list(I, NULL)
  est <- list()
  for (cn in names(cohorts)) for (e in names(ENDPOINTS)) for (h in c("full", "5y")) {
    if (grepl("^A[123]", cn) && e != "mace4") next
    est[[length(est) + 1]] <- est_rows(fit(cohorts[[cn]][[1]], e, if (h == "5y") H5 else Inf, w = cohorts[[cn]][[2]]), cn, e, h, "joint")
  }
  for (cn in c("A", "C_std", "D")) {
    if (is.null(cohorts[[cn]])) next
    D <- cohorts[[cn]][[1]]; ww <- cohorts[[cn]][[2]]
    for (x in EXPOSURES) est[[length(est) + 1]] <- est_rows(fit(D, "mace4", xs = x, w = ww), cn, "mace4", "full", "one exposure")
    est[[length(est) + 1]] <- est_rows(fit(D, "mace4", xs = c("x_lpa_paper", "x_ldl", "x_prs"), w = ww), cn, "mace4", "full", "paper Lp(a) score")
    est[[length(est) + 1]] <- est_rows(fit(D, "mace4", adj = c(ADJ, RISK), w = ww), cn, "mace4", "full", "risk-factor adjusted")
  }
  est[[length(est) + 1]] <- est_rows(fit(A[statin == TRUE], "mace4"), "A statin users", "mace4", "full", "joint")
  est <- rbindlist(est)
  est[, `:=`(hr = exp(b), hr_lo = exp(b - 1.96 * se), hr_hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))]
  # for Lp(a): the trial's direction, per 144 nmol/L LOWER
  est[exposure == "x_lpa", `:=`(hr_per_144_lower = exp(-b), lower_lo = exp(-b - 1.96 * se), lower_hi = exp(-b + 1.96 * se))]

  # pooled established-disease estimate (A and D are disjoint) and the trial-matched comparison
  pick <- function(cn, h) est[cohort == cn & endpoint == "mace4" & horizon == h & model == "joint" & exposure == "x_lpa"]
  trial_se <- (log(TRIAL[["hi"]]) - log(TRIAL[["lo"]])) / (2 * 1.96)
  extra <- rbindlist(lapply(c("full", "5y"), function(h) {
    a <- pick("A", h); d <- pick("D", h); if (nrow(a) == 0) return(NULL)
    rows <- list(data.table(quantity = "A", horizon = h, b = a$b, se = a$se))
    if (nrow(d) == 1) { wA <- 1 / a$se^2; wD <- 1 / d$se^2
      rows[[2]] <- data.table(quantity = "A + D pooled", horizon = h, b = (wA * a$b + wD * d$b) / (wA + wD), se = sqrt(1 / (wA + wD))) }
    rbindlist(rows)
  }))
  extra[, `:=`(hr_per_144_lower = exp(-b), lower_lo = exp(-b - 1.96 * se), lower_hi = exp(-b + 1.96 * se),
               z_vs_trial = (-b - log(TRIAL[["hr"]])) / sqrt(se^2 + trial_se^2))]
  extra[, p_vs_trial := 2 * pnorm(-abs(z_vs_trial))]

  bb <- bias_budget(A, Cc)
  ver <- copy(pa$verdict)
  ver[, `:=`(zph_p_lpa_established = pa$fA$zph[[1]],
             unmodelled_bias_log_hr = if (is.null(bb)) NA_real_ else bb[in_primary_model == FALSE, sum(implied_bias_log_hr)],
             bias_material = if (is.null(bb)) NA else abs(bb[in_primary_model == FALSE, sum(implied_bias_log_hr)]) >
                               0.25 * abs(pa$rs[exposure == "x_lpa", b_primary]))]
  wcsv(pa$rs, file.path(outdir, "stage_contrast.csv")); wcsv(ver, file.path(outdir, "verdict.csv"))
  wcsv(est, file.path(outdir, "estimates.csv")); wcsv(extra, file.path(outdir, "trial_matched.csv"))
  if (!is.null(bb)) wcsv(bb, file.path(outdir, "bias_budget.csv"))
  list(rs = pa$rs, verdict = ver, est = est, extra = extra, bias = bb)
}

# ---------------------------------------------------------------------------
# Real data
# ---------------------------------------------------------------------------
load_real <- function(opt, CENS) {
  cv <- fread(opt$covar); cv[, IID := as.character(IID)]; if ("FID" %in% names(cv)) cv[, FID := NULL]
  gl <- fread(opt$dose_lpa); gl[, IID := as.character(IID)]
  co <- build_cohort(opt$raw_dir, cv, intersect(cv$IID, gl$IID), CENS)
  P <- co$P
  gcols <- setdiff(names(gl), META)
  G <- as.matrix(gl[match(P$IID, IID), ..gcols])
  for (j in seq_len(ncol(G))) { m <- is.na(G[, j]); if (any(m)) G[m, j] <- mean(G[, j], na.rm = TRUE) }
  ex <- fit_exposure(P$lpa, P$status, G)
  P[, pred := ex$pred]
  rdw <- function(f) fread(f, header = FALSE, col.names = c("SNP", "EA", "BETA"))
  sl <- build_score(gl, rdw(opt$w_lpa)); P[, score_lpa_raw := sl$s[match(P$IID, gl$IID)]]
  gd <- fread(opt$dose_ldl); gd[, IID := as.character(IID)]
  sd_ <- build_score(gd, rdw(opt$w_ldl)); P[, score_ldl_raw := sd_$s[match(P$IID, gd$IID)]]
  prs <- fread(opt$prs); prs[, IID := as.character(IID)]
  check("N500" %in% names(prs), "polygenic score file lacks the non-lipid component N500")
  P[, prs_n := prs$N500[match(P$IID, prs$IID)]]
  fe <- finish_exposures(P, co$I)
  list(P = fe$P, I = fe$I, ex = ex, n_in = co$n_in, nv = c(lpa = sl$n_variants, ldl = sd_$n_variants), gcols = gcols)
}
gates <- function(R, CENS, outdir) {
  P <- R$P; A <- P[qual == TRUE]
  evA <- sum(mk(A, "mace4")$ev)
  i72 <- grep("^rs10455872_", R$gcols)
  c72 <- if (length(i72) == 1) abs(R$ex$coef[[i72]]) else NA_real_
  pc <- published_control(P, CENS)
  g <- rbindlist(list(
    data.table(gate = "cohort with a baseline date", value = nrow(P), target = "395,000-402,200", pass = nrow(P) > 395000 & nrow(P) <= 402200),
    data.table(gate = "established-disease cohort size equals the feasibility count", value = nrow(A), target = as.character(FEAS$n_A), pass = nrow(A) == FEAS$n_A),
    data.table(gate = "its 4-point events equal the feasibility count", value = evA, target = as.character(FEAS$ev_A), pass = evA == FEAS$ev_A),
    data.table(gate = "primary-prevention cohort size equals the feasibility count", value = P[broad == FALSE, .N], target = as.character(FEAS$n_C), pass = P[broad == FALSE, .N] == FEAS$n_C),
    data.table(gate = "exposure model: variants dropped as collinear", value = R$ex$dropped, target = "0", pass = R$ex$dropped == 0),
    data.table(gate = "exposure model: |rs10455872| nmol/L per allele", value = round(c72, 1), target = "100-120", pass = isTRUE(c72 > 100 && c72 < 120)),
    data.table(gate = "exposure model: r2 with in-range measurements", value = round(R$ex$r2, 3), target = "> 0.45", pass = R$ex$r2 > 0.45),
    data.table(gate = "share with all three exposures", value = round(mean(complete.cases(P[, ..EXPOSURES])), 4), target = "> 0.99", pass = mean(complete.cases(P[, ..EXPOSURES])) > 0.99),
    data.table(gate = "published model: participants", value = pc$n, target = as.character(PUB$n), pass = pc$n == PUB$n),
    data.table(gate = "published model: events", value = pc$events, target = as.character(PUB$events), pass = pc$events == PUB$events),
    data.table(gate = "published model: Lp(a) score HR", value = round(pc$hr, 3), target = as.character(PUB$hr[1]), pass = abs(pc$hr - PUB$hr[1]) < 0.0015),
    data.table(gate = "published model: lower limit", value = round(pc$lo, 3), target = as.character(PUB$hr[2]), pass = abs(pc$lo - PUB$hr[2]) < 0.0015),
    data.table(gate = "published model: upper limit", value = round(pc$hi, 3), target = as.character(PUB$hr[3]), pass = abs(pc$hi - PUB$hr[3]) < 0.0015)))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE); fwrite(g, file.path(outdir, "gates.csv"))
  say("[gates] variants in scores: Lp(a) %d, LDL-C %d; exposure correlations (no outcome): Lp(a)~LDL-C score %.3f, Lp(a)~non-lipid score %.3f, LDL-C~non-lipid %.3f",
      R$nv[["lpa"]], R$nv[["ldl"]], cor(P$x_lpa, P$x_ldl, use = "complete.obs"), cor(P$x_lpa, P$x_prs, use = "complete.obs"),
      cor(P$x_ldl, P$x_prs, use = "complete.obs"))
  options(width = 200); print(g)
  isTRUE(all(g$pass))
}

# ---------------------------------------------------------------------------
# Self-test (synthetic data only)
# ---------------------------------------------------------------------------
simulate <- function(n, cL, cD, seed, cP = 0.10) {
  set.seed(seed)
  bL <- 0.25; bD <- 0.09; bP <- 0.30; bU <- 0.5
  pred <- pmax(0, rgamma(n, shape = 1.2, scale = 14) + 110 * rbinom(n, 2, 0.08) + rnorm(n, 0, 5))
  P <- data.table(IID = as.character(seq_len(n)), pred = pred, score_ldl_raw = rnorm(n), prs_n = rnorm(n),
                  age = round(runif(n, 40, 70)), sex = rbinom(n, 1, 0.46))
  for (k in 1:10) set(P, j = paste0("PC", k), value = rnorm(n))
  P[, `:=`(bmi = rnorm(n, 27, 4), smoking_current = rbinom(n, 1, 0.1), diabetes = rbinom(n, 1, 0.05), hypertension = rbinom(n, 1, 0.3))]
  P[, score_lpa_raw := pred + rnorm(n, 0, 3)]
  U <- rnorm(n)
  P <- finish_exposures(P, NULL)$P
  lin1 <- bL * P$x_lpa + bD * P$x_ldl + bP * P$x_prs + bU * U + 0.5 * P$sex
  lin2 <- cL * P$x_lpa + cD * P$x_ldl + cP * P$x_prs + bU * U + 0.3 * P$sex
  prev <- rexp(n, 0.0009 * exp(lin1)) < (P$age - 30)
  t0 <- 14000L; CENS <- t0 + as.integer(13 * 365.25)
  te <- t0 + as.integer(ceiling(365.25 * ifelse(prev, rexp(n, 0.025 * exp(lin2)), rexp(n, 0.004 * exp(lin1)))))
  td <- t0 + as.integer(ceiling(365.25 * rexp(n, 0.006)))
  P[, `:=`(t0 = t0, tdeath = fifelse(td <= CENS, td, NA_integer_), t_end = pmin(td, CENS))]
  type <- sample(c("t_mi", "t_is", "t_rev"), n, replace = TRUE, prob = c(0.45, 0.25, 0.30))
  for (cn in c("t_mi", "t_is", "t_rev")) set(P, j = cn, value = ifelse(type == cn & te <= P$t_end, te, NA_integer_))
  P[, t_cvd := NA_integer_]
  P[, `:=`(qual = prev, broad = prev, mi_pre = prev & type == "t_mi", is_pre = prev & type == "t_is", pad_pre = prev & type == "t_rev",
           statin = rbinom(n, 1, ifelse(prev, 0.8, 0.12)) == 1, status = "measured", lpa = pred)]
  P[, `:=`(statin_first = as.integer(statin), tcad = fifelse(prev, t0 - 100L, fifelse(te <= t_end, te, NA_integer_)))]
  # second cohort: first MI or stroke after baseline, then a further event
  I <- P[broad == FALSE & (!is.na(t_mi) | !is.na(t_is))]
  I[, index := pmin(t_mi, t_is, na.rm = TRUE)]; I <- I[t_end > index + 30L]
  k <- match(I$IID, P$IID)
  nxt <- I$index + 30L + as.integer(ceiling(365.25 * rexp(nrow(I), 0.03 * exp(lin2[k]))))
  I[, `:=`(t0 = index + 30L, t_mi = fifelse(nxt <= t_end, nxt, NA_integer_), t_is = NA_integer_, t_rev = NA_integer_)]
  I[, age := age + (t0 - 14000L) / 365.25]
  list(P = P, I = I, CENS = CENS, planted = c(bL = bL, bD = bD, cL = cL, cD = cD))
}

fixture <- function(dir) {
  T0 <- "2008-06-01"
  ppl <- list(
    healthy = list(), no_date = list(t0 = ""),
    primary_mi = list(hes = list(c("I219", "2015-03-01")), fo = c("131298" = "2015-03-01")),
    prev_mi_hes = list(hes = list(c("I214", "2005-01-01"), c("I219", "2012-05-05")), fo = c("131298" = "2005-01-01")),
    prev_mi_self = list(s6150 = "1", hes = list(c("I219", "2011-01-01")), fo = c("131298" = "1900-01-01")),
    prev_mi_same = list(hes = list(c("I219", "2004-02-02")), fo = c("131298" = "2004-02-02"), death = c("2016-02-02", "I251")),
    old_mi = list(hes = list(c("I252", "2007-01-01"))),
    prev_stroke = list(hes = list(c("I639", "2006-01-01"), c("I210", "2014-04-04")), fo = c("131366" = "2006-01-01", "131298" = "2014-04-04")),
    pad_self = list(s20002 = "1087"),
    angina_only = list(s6150 = "2"),
    revasc_pre = list(ops = list(c("K751", "2007-07-07"), c("K491", "2010-10-10"))),
    future_fo = list(fo = c("131298" = "2037-07-07")),
    lost_then_event = list(lost = "2012-01-01", hes = list(c("I639", "2013-01-01")), fo = c("131366" = "2013-01-01")),
    noncv_death = list(death = c("2015-05-05", "C349")),
    stroke_then_mi = list(hes = list(c("I639", "2012-01-01"), c("I214", "2014-06-01")), fo = c("131366" = "2012-01-01", "131298" = "2014-06-01"),
                          ops = list(c("K751", "2012-01-05"), c("K401", "2015-01-01"))),
    mi_dies_fast = list(hes = list(c("I219", "2013-01-01")), fo = c("131298" = "2013-01-01"), death = c("2013-01-10", "I219")),
    tia_pre = list(hes = list(c("G459", "2006-06-06"))),
    statin_slot2 = list(m6153 = c("2", "1")),
    lpa_high = list(rep = "5"), lpa_low = list(rep = "4"), lpa_meas = list(lpa = "50.5", rep = "1"))
  ids <- as.character(1000 + seq_along(ppl)); names(ids) <- names(ppl)
  g <- function(p, k, default = "") if (is.null(p[[k]])) default else p[[k]]
  wide <- function(fid, ncol, getter) {
    d <- data.table(eid = unname(ids))
    for (j in seq_len(ncol)) set(d, j = sprintf("%s-0.%d", fid, j - 1L), value = vapply(ppl, function(p) { v <- getter(p); if (length(v) >= j) v[[j]] else "" }, ""))
    fwrite(d, file.path(dir, sprintf("fid%s.csv", fid)), quote = TRUE)
  }
  wide("53", 1, function(p) g(p, "t0", T0))
  wide("40000", 1, function(p) g(p, "death", c("", ""))[1]); wide("40001", 1, function(p) g(p, "death", c("", ""))[2])
  wide("191", 1, function(p) g(p, "lost"))
  wide("41270", 3, function(p) vapply(g(p, "hes", list()), `[`, "", 1)); wide("41280", 3, function(p) vapply(g(p, "hes", list()), `[`, "", 2))
  wide("41272", 2, function(p) vapply(g(p, "ops", list()), `[`, "", 1)); wide("41282", 2, function(p) vapply(g(p, "ops", list()), `[`, "", 2))
  for (f in c("131296", "131298", "131300", "131302", "131304", "131306", "131366", "131368", "42000", "42006", "42008"))
    wide(f, 1, function(p) { v <- g(p, "fo", character()); if (f %in% names(v)) v[[f]] else "" })
  wide("6150", 2, function(p) g(p, "s6150", character())); wide("20002", 2, function(p) g(p, "s20002", character()))
  wide("20004", 1, function(p) g(p, "s20004", character()))
  wide("6153", 2, function(p) g(p, "m6153", "-7")); wide("6177", 2, function(p) character())
  wide("30780", 1, function(p) "3.5"); wide("30790", 1, function(p) g(p, "lpa")); wide("30796", 1, function(p) g(p, "rep"))
  cv <- data.table(IID = unname(ids), age = 60, sex = 1)
  for (k in 1:10) set(cv, j = paste0("PC", k), value = 0)
  cv[, `:=`(bmi = 27, smoking_current = 0, diabetes = 0, hypertension = 0)]
  list(cv = cv, ids = ids)
}

selftest <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { check(cond, paste("selftest:", msg)); n_ok <<- n_ok + 1L }
  D <- function(s) as.integer(as.IDate(s))

  # 1. data assembly on a fixture
  fx_dir <- file.path(workdir, "fixture"); dir.create(fx_dir, recursive = TRUE, showWarnings = FALSE)
  fx <- fixture(fx_dir); CENS <- D("2022-10-31")
  co <- build_cohort(fx_dir, fx$cv, unname(fx$ids), CENS); P <- co$P; I <- co$I
  row <- function(nm) P[IID == fx$ids[[nm]]]
  expect(nrow(P) == length(fx$ids) - 1L && nrow(row("no_date")) == 0, "a participant without a baseline date is dropped")
  expect(!row("healthy")$broad && row("healthy")$t_end == CENS && is.na(row("healthy")$t_mi), "healthy participant has no disease and no event")
  expect(!row("primary_mi")$broad && row("primary_mi")$t_mi == D("2015-03-01"), "first MI after baseline is an event in primary prevention")
  expect(row("prev_mi_hes")$mi_pre && row("prev_mi_hes")$qual && row("prev_mi_hes")$t_mi == D("2012-05-05"), "later MI under a new 4-character code is a recurrent event")
  expect(row("prev_mi_self")$mi_pre && row("prev_mi_self")$t_mi == D("2011-01-01"), "self-reported MI qualifies; first hospital MI afterwards counts")
  expect(row("prev_mi_same")$mi_pre && is.na(row("prev_mi_same")$t_mi) && row("prev_mi_same")$t_cvd == D("2016-02-02"), "no new code means no visible MI; cardiovascular death counts")
  expect(row("old_mi")$mi_pre, "I25.2 (old MI) qualifies")
  expect(row("prev_stroke")$is_pre && !row("prev_stroke")$mi_pre && row("prev_stroke")$t_mi == D("2014-04-04"), "MI after a prior stroke is an event")
  expect(row("pad_self")$pad_pre && row("pad_self")$qual, "self-reported claudication qualifies")
  expect(row("angina_only")$broad && !row("angina_only")$qual, "angina alone is vascular disease but not a qualifying condition")
  expect(row("revasc_pre")$broad && !row("revasc_pre")$qual && row("revasc_pre")$t_rev == D("2010-10-10"), "revascularization before baseline excludes from primary prevention; a new procedure code is an event")
  expect(!row("future_fo")$broad && is.na(row("future_fo")$t_mi), "the 2037 placeholder date is ignored")
  expect(row("lost_then_event")$t_end == D("2012-01-01") && is.na(row("lost_then_event")$t_is), "events after loss to follow-up are dropped")
  expect(is.na(row("noncv_death")$t_cvd) && row("noncv_death")$t_end == D("2015-05-05"), "non-cardiovascular death censors without an event")
  expect(row("tia_pre")$broad && !row("tia_pre")$qual, "TIA is vascular disease but not qualifying")
  expect(row("statin_slot2")$statin && row("statin_slot2")$statin_first == 0L, "medication in the second slot is found; the paper's first-slot rule is kept separately")
  expect(row("lpa_high")$status == "above_range" && row("lpa_low")$status == "below_range" && row("lpa_meas")$status == "measured" && row("healthy")$status == "not_measured", "Lp(a) reportability codes map to the right status")
  expect(row("prev_mi_same")$tcad == D("2004-02-02"), "the paper's CAD date is the earliest first-occurrence date")
  s <- I[IID == fx$ids[["stroke_then_mi"]]]
  expect(nrow(s) == 1 && s$index_type == "IS" && s$t0 == D("2012-01-01") + 30L && s$t_mi == D("2014-06-01") && s$t_rev == D("2015-01-01"), "second cohort: follow-up starts 30 days after the index stroke; a procedure inside 30 days is not an event")
  expect(nrow(I[IID == fx$ids[["mi_dies_fast"]]]) == 0 && fx$ids[["primary_mi"]] %in% I$IID, "second cohort requires survival beyond 30 days")
  expect(sum(mk(P[qual == TRUE], "mace4")$ev) == 4L && sum(mk(P[qual == TRUE], "mace4", H5)$ev) == 2L, "event flags and the five-year horizon count correctly")

  # 2. paper score orientation and the censored exposure model
  g <- data.table(IID = c("a", "b"), rs1_A = c(0, 2), rs2_G = c(1, NA))
  sc <- build_score(g, data.table(SNP = c("rs1", "rs2"), EA = c("A", "C"), BETA = c(0.5, 2)))
  expect(isTRUE(all.equal(sc$s, c(0 + (2 - 1) * 2, 1 + (2 - 1) * 2))) && sc$n_variants == 2, "score flips to the effect allele and mean-imputes")
  set.seed(7); n <- 60000
  G <- cbind(v1 = rbinom(n, 2, 0.08), v2 = rbinom(n, 2, 0.3), v3 = rbinom(n, 2, 0.5))
  y <- 30 + 100 * G[, 1] + 40 * G[, 2] - 20 * G[, 3] + rnorm(n, 0, 40)
  st <- ifelse(y < LO, "below_range", ifelse(y > HI, "above_range", "measured")); st[sample(n, 3000)] <- "not_measured"
  ex <- fit_exposure(ifelse(st == "measured", y, NA), st, G)
  expect(max(abs(ex$coef - c(100, 40, -20))) < 2.5 && abs(ex$scale - 40) < 1.5 && mean(st == "above_range") > 0.01 && mean(st == "below_range") > 0.05, "censored model recovers planted effects despite both tails being cut")
  naive <- coef(lm(y ~ G, subset = st == "measured"))[-1]
  expect(abs(naive[[1]] - 100) > 10, "the same data analysed without censoring is biased (the check has teeth)")

  # 3. every verdict branch from hand-made inputs
  mkrs <- function(bA, bC, dA, dC, seA = 0.04, seC = 0.02) rbind(
    data.table(exposure = "x_lpa", b_established = bA, b_primary = bC, p_present = pnorm(bA / seA, lower.tail = FALSE),
               p_anchor = pnorm(bC / seC, lower.tail = FALSE), p_diff = 2 * pnorm(-abs((bA - bC) / sqrt(seA^2 + seC^2)))),
    data.table(exposure = "x_ldl", b_established = dA, b_primary = dC, p_present = pnorm(dA / 0.015, lower.tail = FALSE),
               p_anchor = pnorm(dC / 0.006, lower.tail = FALSE), p_diff = 1))
  v <- function(rs, p) verdict(rs, list(diff = 0, se = 1, p_one = p))$verdict
  expect(v(mkrs(0.25, 0.25, 0.09, 0.09), 0.5) == "GENETICS_DOES_NOT_REPRODUCE_TRIAL", "verdict: persists")
  expect(v(mkrs(0.12, 0.30, 0.09, 0.09), 0.5) == "PARTIAL_PERSISTENCE", "verdict: partial")
  expect(v(mkrs(0.00, 0.25, 0.09, 0.09), 0.01) == "STAGE_SPECIFIC_FOR_LPA", "verdict: stage specific")
  expect(v(mkrs(0.00, 0.25, 0.09, 0.09), 0.40) == "ATTENUATED_NOT_SPECIFIC", "verdict: attenuated, not specific")
  expect(v(mkrs(0.00, 0.25, 0.00, 0.09), 0.40) == "UNINFORMATIVE_DESIGN", "verdict: control absent")
  expect(v(mkrs(0.05, 0.10, 0.09, 0.09, seA = 0.06), 0.40) == "INDETERMINATE", "verdict: imprecise")
  expect(v(mkrs(0.00, 0.02, 0.09, 0.09), 0.01) == "NO_ANCHOR", "verdict: no primary-prevention effect to compare with")

  # 4. planted effects: three scenarios through the real primary analysis
  scen <- list(c(cL = 0, cD = 0.09, want = "STAGE_SPECIFIC_FOR_LPA"), c(cL = 0.25, cD = 0.09, want = "GENETICS_DOES_NOT_REPRODUCE_TRIAL"),
               c(cL = 0, cD = 0, want = "UNINFORMATIVE_DESIGN"))
  for (i in seq_along(scen)) {
    s <- simulate(400000, as.numeric(scen[[i]][["cL"]]), as.numeric(scen[[i]][["cD"]]), seed = 100 + i)
    pa <- primary_analysis(s$P); L <- pa$rs[exposure == "x_lpa"]; Dl <- pa$rs[exposure == "x_ldl"]
    say("[sim %d] Lp(a) established %.3f (SE %.3f) planted %.2f | primary %.3f (SE %.3f) planted %.2f | LDL-C established %.3f planted %.2f | verdict %s",
        i, L$b_established, L$se_established, s$planted[["cL"]], L$b_primary, L$se_primary, s$planted[["bL"]], Dl$b_established, s$planted[["cD"]], pa$verdict$verdict)
    expect(abs(L$b_primary - s$planted[["bL"]]) < 4 * L$se_primary, sprintf("scenario %d recovers the primary-prevention Lp(a) effect", i))
    expect(abs(L$b_established - s$planted[["cL"]]) < 4 * L$se_established, sprintf("scenario %d recovers the established-disease Lp(a) effect", i))
    expect(abs(Dl$b_established - s$planted[["cD"]]) < 4 * Dl$se_established, sprintf("scenario %d recovers the LDL-C effect", i))
    expect(pa$verdict$verdict == scen[[i]][["want"]], sprintf("scenario %d gives verdict %s", i, scen[[i]][["want"]]))
  }

  # 5. smoke run of everything the real run executes after the gates
  s <- simulate(120000, 0, 0.09, seed = 200)
  out <- run_analysis(s$P, s$I, file.path(workdir, "smoke"))
  pc <- published_control(s$P, s$CENS)
  expect(all(file.exists(file.path(workdir, "smoke", c("stage_contrast.csv", "verdict.csv", "estimates.csv", "trial_matched.csv", "bias_budget.csv")))), "every output file is written")
  expect(nrow(out$est[cohort == "D" & exposure == "x_lpa"]) > 0 && nrow(out$extra[quantity == "A + D pooled"]) == 2 && pc$n > 0, "second cohort, pooling and the published-model code all run")
  expect(all(is.finite(out$bias$implied_bias_log_hr)) && is.finite(out$verdict$unmodelled_bias_log_hr), "bias budget is finite")
  cat(sprintf("\nSELFTEST PASSED: %d checks\n", n_ok))
}

# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
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
  if (opt$step == "selftest") { selftest(opt$workdir); quit(save = "no", status = 0L) }
  check(opt$step %in% c("gates", "test"), "unknown --step")
  CENS <- as.integer(as.IDate(opt$censor))
  R <- load_real(opt, CENS)
  passed <- gates(R, CENS, opt$outdir)
  if (!passed) { cat("\nGATES FAILED: no estimate is to be read.\n"); quit(save = "no", status = 2L) }
  cat("\nGATES PASSED\n")
  if (opt$step == "gates") quit(save = "no", status = 0L)
  out <- run_analysis(R$P, R$I, opt$outdir)
  options(width = 220)
  cat("\n=== stage contrast (log hazard ratios; Lp(a) per 144 nmol/L higher, scores per SD) ===\n"); print(out$rs, digits = 3)
  cat("\n=== verdict ===\n"); print(t(out$verdict), quote = FALSE)
  cat("\n=== trial-matched (hazard ratio per 144 nmol/L LOWER predicted Lp(a); trial 1.02, 0.90-1.15) ===\n"); print(out$extra, digits = 3)
  cat("\n=== 4-point composite, joint model ===\n")
  print(out$est[endpoint == "mace4" & model == "joint", .(cohort, horizon, exposure, n, events, hr, hr_lo, hr_hi, p, hr_per_144_lower)], digits = 3)
  cat("\n=== selection (index-event) bias budget ===\n"); print(out$bias, digits = 3)
  cat("\nANALYSIS_DONE\n")
}
