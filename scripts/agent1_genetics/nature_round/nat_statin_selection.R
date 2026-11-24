#!/usr/bin/env Rscript
# ============================================================================
# Nature round, analysis B: is the statin x score asymmetry explained by who
# gets statins?  Read nature_round/PREREGISTERED_READOUT.md section B first; the
# verdict rules at the bottom of this script are that section, nothing looser.
#
# WHAT IS BEING TESTED
#   Published (Cox on incident CAD, score*statin + age + sex + PC1-10):
#     LDL-C score x statin ratio of HRs 0.93 (0.90-0.95)
#     Lp(a) score x statin ratio of HRs 0.99 (0.97-1.01)
#   Statin use is a common effect of LDL-C and of cardiovascular risk, so
#   comparing score slopes across statin strata conditions on a collider even
#   when fitted on the whole cohort (Tudball 2026). None of the published models
#   contained measured LDL-C, the main prescribing indication.
#
# WHAT THE MODELS BELOW CAN AND CANNOT DO
#   Measured LDL-C is itself downstream of both the score and the statin, so the
#   LDL-C-adjusted Cox models and the propensity models that contain it condition
#   on a descendant of the collider. They are diagnostics: if the interaction moves
#   when LDL-C enters, the interaction is sensitive to how selection is handled.
#   They do not recover an unbiased effect, and /0.7 is one constant applied to
#   everyone, not each person's untreated LDL-C. That is why the readout requires
#   the simulation as well as IPTW before the word ROBUST may be used.
#
# FAILURE BEHAVIOUR
#   - The published estimates must reproduce first (hard stop, exit status 2,
#     only positive_control.csv is written).
#   - Degenerate designs return NA, never a number.
#   - Outputs are summary level and pass through wcsv(), which refuses anything
#     that looks individual level; counts below 5 are suppressed.
#
# Run: see slurm/nat_B_selection.slurm. Selftest (synthetic data only):
#   Rscript nat_statin_selection.R --selftest --selftest_dir <scratch dir>
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(optparse); library(survival); library(splines); library(parallel)
})
set.seed(42)

# Published values the data prep must reproduce. Manuscript Results and Table oS34:
# LDL-C 0.927 (0.903-0.951), Lp(a) 0.991 (0.968-1.015), main effects 1.092 and 1.096,
# non-users 327,254 / 23,535 events, statin users 52,642 / 8,518 events.
# Compared at two decimals, as written in the paper. The Lp(a) upper limit is the catch:
# the table prints 1.015 and the text 1.01, so the true value sits on the rounding edge
# and 1.0151 must not stop the run. A value therefore passes if it matches at 2 decimals
# OR equals the table at 3 decimals; the 3-decimal targets are what the *_3dp fields hold.
# The cohort flow is printed, but only the prevalent exclusion is checked. The "start" of
# 402,119 this file used to enforce is not a published number: it traces to an unrelated Lp(a)
# diagnostic (logs/claude_jobs/lpa_diag_3542525.out, "merged rows 402200 in-field 402119"), and
# the 376 "dropped" was the arithmetic residual 402,119 - 21,847 - 379,896 rather than a count
# anyone reported. Jobs 3545604 and 3549233 both stopped on that pair while every genuine target
# matched. The real start of this cohort is 401,743, which the flow now prints without a target.
# Of the enforced rows, the four estimates are in the manuscript; the counts (n, events, statin
# users, events in statin users, prevalent) come from a prior run of the paper's own script, not
# from the text, so they verify reproducibility rather than agreement with a published table.
PUBLISHED <- list(
  ldl_int  = c(0.93, 0.90, 0.95), lpa_int  = c(0.99, 0.97, 1.01),
  ldl_main = c(1.09, 1.08, 1.11), lpa_main = c(1.10, 1.08, 1.11),
  ldl_int_3dp  = c(0.927, 0.903, 0.951), lpa_int_3dp  = c(0.991, 0.968, 1.015),
  ldl_main_3dp = c(1.092, 1.078, 1.106), lpa_main_3dp = c(1.096, 1.083, 1.110),
  n = 379896L, events = 32053L, n_statin = 52642L, events_statin = 8518L,
  prevalent = 21847L,
  # Full precision from the same published run (nat_tudball_paper_scores.R PC2). The 2-decimal test
  # above passes a rebuilt score giving 0.934; these do not. Absent from the selftest's targets.
  ldl_int_full = 0.92655, lpa_int_full = 0.99125)

SMD_LIMIT <- 0.1       # PREREGISTERED_READOUT.md section B; ANALYSIS_PLAN.md section 6
MIN_EVENTS <- 10L      # fewer events than this in a statin stratum and the model is not fitted
MIN_CELL <- 5L         # export rule

# ---- small utilities --------------------------------------------------------
sup <- function(x) { x <- as.integer(x); x[!is.na(x) & x < MIN_CELL] <- NA_integer_; x }

# One gate for every file this script writes, so that a coding slip cannot put a
# row per participant on disk. No analysis output here has anywhere near 5,000 rows.
wcsv <- function(dt, outdir, name) {
  if (nrow(dt) > 5000L) stop("refusing to write ", name, ": ", nrow(dt), " rows looks individual level")
  if (any(tolower(names(dt)) %in% c("iid", "fid", "eid", "id")))
    stop("refusing to write ", name, ": identifier column present")
  fwrite(dt, file.path(outdir, name))
  invisible(dt)
}

# Runs a model and keeps its warnings: a fit that did not converge has to become
# NA rather than a number with a tidy-looking confidence interval.
guarded <- function(expr) {
  msgs <- character()
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL }),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = val, msgs = msgs)
}
BAD_FIT <- "converg|iteration|singular|numerically 0 or 1|infinite"

na_row <- function(model, term, note, n = NA_integer_, events = NA_integer_)
  data.table(model = model, term = term, estimate = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
             P = NA_real_, n = sup(n), events = sup(events), note = note)

# R names an interaction by variable-introduction order, not by how it is written:
# in ~ score_LPA*statin + score_LDL*statin the second is "statin:score_LDL". The
# published script silently lost a coefficient to this (commit 04c6c6f).
term_index <- function(cf, term) {
  i <- which(rownames(cf) == term)
  if (!length(i) && grepl(":", term, fixed = TRUE))
    i <- which(rownames(cf) == paste(rev(strsplit(term, ":", fixed = TRUE)[[1]]), collapse = ":"))
  i
}

design_problem <- function(dd, score) {
  if (!nrow(dd)) return("empty design")
  if (length(unique(dd$statin)) < 2L) return("only one statin stratum")
  ev <- dd[, .(e = sum(event), v = var(get(score))), by = statin]
  if (any(ev$e < MIN_EVENTS)) return(sprintf("fewer than %d events in a statin stratum", MIN_EVENTS))
  if (any(is.na(ev$v) | ev$v <= 0)) return("no score variation in a statin stratum")
  NA_character_
}

# Cox wrapper. Uses 1.96, not qnorm(.975), because the published script did and the
# Lp(a) upper limit sits at a rounding boundary.
fit_cox <- function(formula_txt, dd, model, term, score, wt = NULL, note = "") {
  prob <- design_problem(dd, score)
  if (!is.na(prob)) return(na_row(model, term, paste("NA:", prob), nrow(dd), sum(dd$event)))
  f <- as.formula(formula_txt)
  g <- guarded(if (is.null(wt)) coxph(f, data = dd)
               else { dd2 <- copy(dd); dd2[, .w := wt]; coxph(f, data = dd2, weights = .w, robust = TRUE) })
  if (is.null(g$value) || any(grepl(BAD_FIT, g$msgs, ignore.case = TRUE)))
    return(na_row(model, term, paste("NA: fit failed:", paste(g$msgs, collapse = "; ")), nrow(dd), sum(dd$event)))
  cf <- summary(g$value)$coefficients
  i <- term_index(cf, term)
  if (!length(i) || is.na(cf[i, "coef"]))
    return(na_row(model, term, "NA: term not estimated", nrow(dd), sum(dd$event)))
  se <- if ("robust se" %in% colnames(cf)) cf[i, "robust se"] else cf[i, "se(coef)"]
  data.table(model = model, term = term, estimate = exp(cf[i, "coef"]),
             lo95 = exp(cf[i, "coef"] - 1.96 * se), hi95 = exp(cf[i, "coef"] + 1.96 * se),
             P = 2 * pnorm(-abs(cf[i, "coef"] / se)), n = sup(g$value$n), events = sup(g$value$nevent),
             note = note)
}

# ---- data prep, copied from R/statin_interaction_incident.R ------------------
# Deliberately the same statements in the same order; the selftest runs the published
# script on the same synthetic files and requires identical estimates.
rdcsv <- function(fid, raw_dir) {
  f <- file.path(raw_dir, sprintf("fid%s.csv", fid))
  d <- fread(f, na.strings = c("", "NA"))
  col <- grep(sprintf("^\"?%s-0\\.0\"?$", fid), names(d), value = TRUE)[1]
  setnames(d, c("eid", col), c("IID", "v"))
  d[, .(IID = as.character(IID), v)]
}

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

prep_cohort <- function(opt) {
  raws <- strsplit(opt$scores, ",")[[1]]; wts <- strsplit(opt$weights, ",")[[1]]
  sc <- NULL
  for (k in seq_along(raws)) {
    nm <- if (grepl("LPA", raws[k])) "score_LPA" else "score_LDL"
    x <- build_score(raws[k], wts[k], nm); setnames(x, "s", nm)
    sc <- if (is.null(sc)) x else merge(sc, x, by = "IID")
  }
  base <- rdcsv("53", opt$raw_dir);    setnames(base, "v", "t0")
  dth  <- rdcsv("40000", opt$raw_dir); setnames(dth,  "v", "tdeath")
  CADF <- c("131296","131298","131300","131302","131304","131306")  # I20 I21 I22 I23 I24 I25
  ev <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
               lapply(CADF, function(f) { x <- rdcsv(f, opt$raw_dir); setnames(x,"v",paste0("d",f)); x }))
  ev[, tcad := do.call(pmin, c(lapply(paste0("d",CADF), function(c) as.IDate(get(c))), na.rm = TRUE))]
  st <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
               list(rdcsv("6153", opt$raw_dir)[, .(IID, f = v)], rdcsv("6177", opt$raw_dir)[, .(IID, m = v)]))
  st[, statin := as.integer((!is.na(f) & f == 1) | (!is.na(m) & m == 1))]
  st[is.na(f) & is.na(m), statin := NA_integer_]
  cv <- fread(opt$covar); cv[, IID := as.character(IID)]
  # The published flow starts before the covariate file is joined (job 3545604: 401,743 here vs
  # 402,119 published, and 0 vs 376 lost later). Count the start there; whoever the inner join
  # drops is counted with the 376. Every estimate and the final cohort are unaffected.
  d0 <- Reduce(function(a,b) merge(a,b,by="IID"), list(sc, base[, .(IID,t0)], ev[, .(IID,tcad)], st[, .(IID,statin)]))
  d0 <- d0[!is.na(t0) & !is.na(statin)]
  n_start_precov <- nrow(d0)
  lost_cov <- d0[!IID %in% cv$IID]
  message(sprintf("  [flow] %d before the covariate join; %d not in the covariate file, %d of them with CAD before baseline",
                  n_start_precov, nrow(lost_cov), lost_cov[!is.na(tcad) & as.IDate(tcad) <= as.IDate(t0), .N]))
  rm(d0)
  d <- Reduce(function(a,b) merge(a,b,by="IID"), list(sc, base[, .(IID,t0)], ev[, .(IID,tcad)],
                                                      st[, .(IID,statin)], cv))
  d <- merge(d, dth[, .(IID,tdeath)], by = "IID", all.x = TRUE)
  d[, `:=`(t0 = as.IDate(t0), tdeath = as.IDate(tdeath))]
  cens <- as.IDate(opt$censor)
  d <- d[!is.na(t0) & !is.na(statin)]
  n_start <- nrow(d)
  prev <- d[!is.na(tcad) & tcad <= t0, .N]
  d <- d[is.na(tcad) | tcad > t0]
  d[, tend  := pmin(fifelse(is.na(tcad), cens, tcad), fifelse(is.na(tdeath), cens, tdeath), cens, na.rm = TRUE)]
  d[, event := as.integer(!is.na(tcad) & tcad <= tend)]
  d[, fu    := as.numeric(tend - t0) / 365.25]
  n_prezero <- nrow(d)
  d <- d[fu > 0]
  d[, score_LPA := scale(score_LPA)[,1]]; d[, score_LDL := scale(score_LDL)[,1]]
  # ---- end of the published prep. ----
  # The published n of 379,896 is coxph's m$n, after its listwise deletion on age, sex and
  # PCs, while the scores above were standardized BEFORE that deletion. Doing the deletion
  # here, after scaling, leaves every estimate unchanged and makes nrow(d) mean the published
  # n, which the simulation and the descriptive tables rely on. Below this, columns only.
  pcs0 <- grep("^PC", names(d), value = TRUE)[1:10]
  n_scaled <- nrow(d)
  d <- d[complete.cases(d[, c("age", "sex", pcs0), with = FALSE])]
  n_cov_missing <- n_scaled - nrow(d)

  add_field <- function(fid, nm, required = TRUE) {
    f <- file.path(opt$raw_dir, sprintf("fid%s.csv", fid))
    if (!file.exists(f)) {
      if (required) stop("required field file missing: ", f)
      warning("optional field file missing: ", f, call. = FALSE); d[, (nm) := NA_real_]; return(invisible())
    }
    x <- rdcsv(fid, opt$raw_dir); setnames(x, "v", nm)
    d <<- merge(d, x, by = "IID", all.x = TRUE)
  }
  add_field("30780", "LDL"); add_field("30790", "LPA"); add_field("26252", "score_MEL")
  add_field("30786", "LDL_rep", required = FALSE); add_field("30796", "LPA_rep", required = FALSE)
  d[, LDL_corr := fifelse(statin == 1L, LDL / 0.7, LDL)]         # same 0.7 the LDL GWAS phenotype used
  d[, score_MEL := scale(score_MEL)[,1]]
  # Lp(a) is heavily right skewed; its log is the scale the adjustment is plausibly linear on
  d[, LPA_log := log(LPA)]; d[, LDL_log := log(LDL)]
  d[!is.finite(LPA_log), LPA_log := NA_real_]; d[!is.finite(LDL_log), LDL_log := NA_real_]
  # Any two-level numeric coding gives the same fit; a factor would break the simulation matrix
  if (!is.numeric(d$sex)) stop("covariate sex must be numeric: the simulation and propensity models use it as a number")
  list(d = d, counts = list(start = n_start_precov, prevalent = prev, no_followup = n_prezero - n_scaled,
                            covariate_missing = n_cov_missing + (n_start_precov - n_start)))
}

# 6153 and 6177 are multi-select. The published statin flag reads array 0 only. If a
# code 1 sits in a later array the "non-user" arm contains users, which matters for a
# paper about selection, so count it. A diagnostic only: the flag is not changed.
statin_array_diagnostic <- function(d, raw_dir) {
  ids <- d[statin == 0L, IID]; hit <- character()
  for (fid in c("6153", "6177")) {
    x <- fread(file.path(raw_dir, sprintf("fid%s.csv", fid)), na.strings = c("", "NA"))
    cols <- grep(sprintf("^\"?%s-0\\.[1-9]\"?$", fid), names(x), value = TRUE)
    if (!length(cols)) next
    any1 <- rowSums(as.matrix(x[, ..cols]) == 1, na.rm = TRUE) > 0
    hit <- union(hit, as.character(x$eid[any1]))
  }
  length(intersect(ids, hit))
}

# ---- positive control: hard stop --------------------------------------------
positive_control <- function(d, counts, opt) {
  pcs <- grep("^PC", names(d), value = TRUE)[1:10]; PCS <- paste(pcs, collapse = " + ")
  tg <- opt$published
  one <- function(S) {
    m <- coxph(as.formula(sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", S, PCS)), data = d)
    cf <- summary(m)$coefficients
    r <- function(i) c(exp(cf[i, 1]), exp(cf[i, 1] - 1.96 * cf[i, 3]), exp(cf[i, 1] + 1.96 * cf[i, 3]))
    list(main = r(S), int = r(term_index(cf, paste0(S, ":statin"))), statin = r("statin"), P = cf[term_index(cf, paste0(S, ":statin")), 5],
         n = m$n, events = m$nevent)
  }
  L <- one("score_LDL"); P <- one("score_LPA")
  est <- function(lbl, got, want, want3 = NULL) data.table(check = paste0(lbl, c(" estimate"," lo95"," hi95")),
    observed = got, published = want,
    pass = !is.na(got) & (abs(got - want) <= 0.005 + 1e-9 | (if (is.null(want3)) FALSE else abs(got - want3) <= 0.0005 + 1e-9)))
  cnt <- function(lbl, got, want) data.table(check = lbl, observed = as.numeric(sup(got)),
    published = as.numeric(want), pass = isTRUE(as.integer(got) == as.integer(want)))
  # Flow steps with no published target: printed and exported, never a pass or a fail.
  obs <- function(lbl, got) data.table(check = paste(lbl, "[no published target: reported, not checked]"),
    observed = as.numeric(sup(got)), published = NA_real_, pass = NA)
  tab <- rbindlist(list(
    est("LDL-C x statin ratio of HRs", L$int, tg$ldl_int, tg$ldl_int_3dp), est("Lp(a) x statin ratio of HRs", P$int, tg$lpa_int, tg$lpa_int_3dp),
    est("LDL-C score HR/SD", L$main, tg$ldl_main, tg$ldl_main_3dp),        est("Lp(a) score HR/SD", P$main, tg$lpa_main, tg$lpa_main_3dp),
    if (!is.null(tg$ldl_int_full)) data.table(check = c("LDL-C ratio of HRs, full precision", "Lp(a) ratio of HRs, full precision"),
      observed = c(L$int[1], P$int[1]), published = c(tg$ldl_int_full, tg$lpa_int_full),
      pass = abs(c(L$int[1], P$int[1]) - c(tg$ldl_int_full, tg$lpa_int_full)) <= 5e-4),
    # the model's own n, as the published table reports; nrow(d) must equal it after prep's deletion
    obs("start: baseline date, medication report, scores", counts$start),
    cnt("prevalent excluded", counts$prevalent, tg$prevalent),
    obs("no follow-up or missing covariate", counts$no_followup + counts$covariate_missing),
    cnt("n", L$n, tg$n), cnt("events", L$events, tg$events),
    cnt("statin users", sum(d$statin == 1L), tg$n_statin),
    cnt("events in statin users", sum(d$event[d$statin == 1L]), tg$events_statin)))
  if (L$n != nrow(d)) stop("internal: model n differs from nrow(d) after listwise deletion")
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  wcsv(tab, opt$outdir, "positive_control.csv")
  print(tab)
  if (!isTRUE(all(tab$pass[!is.na(tab$pass)]))) {
    bad <- tab[!is.na(pass) & pass == FALSE, check]
    stop(structure(class = c("positive_control_failure", "error", "condition"), list(
      message = paste0("POSITIVE CONTROL FAILED, no estimate is read until this is fixed: ",
                       paste(bad, collapse = "; ")), call = NULL)))
  }
  message("[B] positive control passed: published interactions reproduce to 2 decimals")
  list(ldl = L, lpa = P, table = tab)
}

# ---- (1) who gets statins ----------------------------------------------------
fit_logit <- function(rhs, dd, model, terms) {
  if (length(unique(na.omit(dd$statin))) < 2L || min(table(dd$statin)) < MIN_EVENTS)
    return(rbindlist(lapply(terms, function(t) na_row(model, t, "NA: statin use has no contrast", nrow(dd)))))
  g <- guarded(glm(as.formula(paste("statin ~", rhs)), data = dd, family = binomial()))
  if (is.null(g$value) || any(grepl(BAD_FIT, g$msgs, ignore.case = TRUE)))
    return(rbindlist(lapply(terms, function(t) na_row(model, t, paste("NA: fit failed:", paste(g$msgs, collapse = "; ")), nrow(dd)))))
  cf <- summary(g$value)$coefficients
  rbindlist(lapply(terms, function(t) {
    if (!t %in% rownames(cf)) return(na_row(model, t, "NA: term not estimated", nobs(g$value)))
    data.table(model = model, term = t, estimate = exp(cf[t, 1]), lo95 = exp(cf[t, 1] - 1.96 * cf[t, 2]),
               hi95 = exp(cf[t, 1] + 1.96 * cf[t, 2]), P = cf[t, 4], n = sup(nobs(g$value)),
               events = sup(sum(g$value$y)), note = "OR per SD; events column = statin users")
  }))
}

selection_tables <- function(d, PCS) {
  adj <- paste("age + sex +", PCS)
  logit <- rbindlist(list(
    fit_logit(paste("score_LDL +", adj), d, "LDL-C score alone", "score_LDL"),
    fit_logit(paste("score_LPA +", adj), d, "Lp(a) score alone", "score_LPA"),
    fit_logit(paste("score_LDL + score_LPA +", adj), d, "both scores", c("score_LDL", "score_LPA")),
    fit_logit(paste("score_MEL +", adj), d[!is.na(score_MEL)], "melanoma score alone (negative control)", "score_MEL")))
  lin_one <- function(y, S, k, dd) {
    lbl <- sprintf("%s ~ %s, %s", y, S, if (k == 1L) "statin users" else "non-users")
    s <- dd[statin == k & !is.na(get(y))]
    if (nrow(s) < MIN_EVENTS || var(s[[S]]) <= 0)
      return(data.table(model = lbl, beta = NA_real_, lo95 = NA_real_, hi95 = NA_real_, P = NA_real_, n = sup(nrow(s)), note = "NA: too few rows"))
    m <- lm(as.formula(sprintf("%s ~ %s + %s", y, S, adj)), data = s); cf <- summary(m)$coefficients
    data.table(model = lbl, beta = cf[S, 1], lo95 = cf[S, 1] - 1.96 * cf[S, 2], hi95 = cf[S, 1] + 1.96 * cf[S, 2],
               P = cf[S, 4], n = sup(nobs(m)), note = "per SD of score, measured baseline value, no statin correction")
  }
  lin_int <- function(y, S, dd) {
    s <- dd[!is.na(get(y))]; lbl <- sprintf("%s ~ %s x statin (pooled)", y, S)
    if (length(unique(s$statin)) < 2L)
      return(data.table(model = lbl, beta = NA_real_, lo95 = NA_real_, hi95 = NA_real_, P = NA_real_, n = sup(nrow(s)), note = "NA: one stratum"))
    m <- lm(as.formula(sprintf("%s ~ %s * statin + %s", y, S, adj)), data = s); cf <- summary(m)$coefficients
    i <- term_index(cf, paste0(S, ":statin"))
    data.table(model = lbl, beta = cf[i, 1], lo95 = cf[i, 1] - 1.96 * cf[i, 2], hi95 = cf[i, 1] + 1.96 * cf[i, 2],
               P = cf[i, 4], n = sup(nobs(m)), note = "interaction term")
  }
  grid <- list(c("LDL", "score_LDL"), c("LDL_log", "score_LDL"), c("LPA", "score_LPA"), c("LPA_log", "score_LPA"))
  lin <- rbindlist(lapply(grid, function(g) rbind(lin_one(g[1], g[2], 0L, d), lin_one(g[1], g[2], 1L, d), lin_int(g[1], g[2], d))))
  list(logit = logit, linear = lin)
}

# ---- (2) Cox with the measured biomarker -------------------------------------
measured_cox <- function(d, PCS) {
  out <- list()
  base <- function(S) sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", S, PCS)
  add <- function(S, dd, extra, lbl, note = "") {
    out[[length(out) + 1L]] <<- fit_cox(paste(base(S), if (nzchar(extra)) paste("+", extra) else ""),
                                        dd, lbl, paste0(S, ":statin"), S, note = note)
  }
  dL <- d[!is.na(LDL)]; dP <- d[!is.na(LPA)]
  # Same-subset refit first: otherwise the change from adding LDL-C mixes with the
  # change from dropping the people without an LDL-C result.
  for (S in c("score_LDL", "score_LPA")) {
    add(S, dL, "", "base, LDL-C-measured subset")
    add(S, dL, "LDL", "+ measured LDL-C (raw)", "on-treatment value in users")
    add(S, dL, "LDL_corr", "+ LDL-C /0.7 in users", "single constant correction")
    add(S, dL, "LDL + LDL:statin", "+ LDL-C (raw) + LDL-C x statin")
    add(S, dL, "LDL_corr + LDL_corr:statin", "+ LDL-C /0.7 + LDL-C/0.7 x statin")
  }
  # The Lp(a) reportability tails are missing because of the value itself (17% in the
  # two tails together): dropping them truncates on the biomarker the score predicts.
  # So the base refit is shown in the same subset and the n is stated in every row.
  for (S in c("score_LPA", "score_LDL")) {
    add(S, dP, "", "base, Lp(a)-measured subset", "reportability tails excluded")
    add(S, dP, "LPA", "+ measured Lp(a) (raw)", "reportability tails excluded")
    add(S, dP, "LPA_log", "+ measured Lp(a) (log)", "reportability tails excluded")
    add(S, dP, "LPA + LPA:statin", "+ Lp(a) (raw) + Lp(a) x statin", "reportability tails excluded")
    add(S, dP, "LPA_log + LPA_log:statin", "+ Lp(a) (log) + log Lp(a) x statin", "reportability tails excluded")
  }
  rbindlist(out)
}

# ---- (3) IPTW ----------------------------------------------------------------
smd_one <- function(x, A, w) {
  ok <- !is.na(x); x <- x[ok]; A <- A[ok]; w <- w[ok]
  if (sum(A == 1L) < 2L || sum(A == 0L) < 2L) return(c(NA_real_, NA_real_))
  s <- sqrt((var(x[A == 1L]) + var(x[A == 0L])) / 2)
  if (!is.finite(s) || s <= 0) return(c(NA_real_, NA_real_))
  c((mean(x[A == 1L]) - mean(x[A == 0L])) / s,
    (weighted.mean(x[A == 1L], w[A == 1L]) - weighted.mean(x[A == 0L], w[A == 0L])) / s)
}

# counts_for_rule: the readout says EVERY weighted SMD, so everything in the table counts,
# including what is left out of the propensity (melanoma score, log Lp(a)), with one
# structural exception.
# Raw LDL-C and /0.7 LDL-C differ by the statin factor itself, so weights that balance one
# leave the other imbalanced by construction, and requiring both on one set of weights would
# make the rule impossible rather than strict. The other form is reported and not counted;
# each form gets its own arm. An NA SMD fails.
balance_table <- function(dd, covs, in_ps, sw, sw_untrunc, not_counted = character()) {
  rbindlist(lapply(covs, function(v) {
    a <- smd_one(dd[[v]], dd$statin, sw); b <- smd_one(dd[[v]], dd$statin, sw_untrunc)
    data.table(covariate = v, in_propensity = v %in% in_ps, counts_for_rule = !(v %in% not_counted),
               smd_raw = a[1], smd_weighted = a[2], smd_weighted_untruncated = b[2],
               pass = !is.na(a[2]) && abs(a[2]) < SMD_LIMIT)
  }))
}

# The published truncation rule (top 1% only) so the estimate is comparable with the
# published IPTW row. Truncation trades balance for stability, hence both SMD columns.
stabilized <- function(e, A) {
  pt <- mean(A); sw <- ifelse(A == 1L, pt / e, (1 - pt) / (1 - e)); list(raw = sw, trunc = pmin(sw, quantile(sw, 0.99)))
}
cstat <- function(e, A) { r <- rank(e); n1 <- as.numeric(sum(A == 1L)); n0 <- as.numeric(sum(A == 0L)); (sum(r[A == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
ess <- function(w) sum(w)^2 / sum(w^2)

iptw_block <- function(d, pcs, PCS, use_glmnet = TRUE) {
  IND <- c("bmi", "smoking_current", "diabetes", "hypertension")
  cc <- d[complete.cases(d[, c(IND, "LDL", "age", "sex", pcs), with = FALSE])]
  base_cov <- c("score_LDL", "score_LPA", "age", "sex", IND, pcs)
  bal <- list(); est <- list(); diag <- list(); chosen <- list()
  outcome <- function(S) sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", S, PCS)
  run_outcome <- function(dd, sw, lbl, note) {
    for (S in c("score_LDL", "score_LPA", "score_MEL")) {
      keep <- !is.na(dd[[S]])
      est[[length(est) + 1L]] <<- fit_cox(outcome(S), dd[keep], lbl, paste0(S, ":statin"), S, wt = sw[keep], note = note)
    }
  }
  # Unweighted in the IPTW subset, so that weighting and sample restriction can be told apart
  for (S in c("score_LDL", "score_LPA", "score_MEL"))
    est[[length(est) + 1L]] <- fit_cox(outcome(S), cc[!is.na(get(S))], "unweighted, IPTW subset", paste0(S, ":statin"), S)

  if (length(unique(cc$statin)) < 2L) {
    return(list(balance = data.table(), est = rbindlist(est), diag = data.table(),
                chosen = data.table(form = c("raw", "corrected"), spec = NA_character_, balanced = FALSE, max_abs_smd = NA_real_)))
  }

  # DESIGN STAGE, blind to the outcome: no line in this loop touches event or fu. For
  # each LDL-C form, step up a fixed ladder (linear, then splines) and stop at the first
  # specification with every weighted SMD below the limit. Balance is what a propensity
  # model is for, so refining it before anyone sees an outcome is not a forking path. If
  # nothing passes, the linear model is kept and the arm is labelled uninformative.
  # Why two forms: raw LDL-C in users is on treatment, so balancing it does NOT balance
  # pre-treatment LDL-C, which is the actual indication; /0.7 targets that but with one
  # constant. Neither is clearly right. See iptw_arm_state for how the verdict reads them.
  for (form in c("raw", "corrected")) {
    lv <- if (form == "raw") "LDL" else "LDL_corr"; other <- if (form == "raw") "LDL_corr" else "LDL"
    specs <- list(
      linear = paste("statin ~", paste(c(base_cov, lv), collapse = " + ")),
      spline = paste("statin ~", paste(c("score_LDL", "score_LPA", "sex", "smoking_current", "diabetes", "hypertension", pcs,
                                          "ns(age, 4)", "ns(bmi, 4)", sprintf("ns(%s, 4)", lv)), collapse = " + ")))
    picked <- NULL
    for (sp in names(specs)) {
      g <- guarded(glm(as.formula(specs[[sp]]), data = cc, family = binomial()))
      if (is.null(g$value) || any(grepl("converg", g$msgs, ignore.case = TRUE))) { diag[[length(diag) + 1L]] <- data.table(form = form, spec = sp, note = paste("NA: propensity failed:", paste(g$msgs, collapse = "; "))); next }
      e <- fitted(g$value); W <- stabilized(e, cc$statin)
      # log Lp(a) is not in the propensity and is missing in the assay tails, but the rule says
      # every SMD and the score it tracks is half of this paper's contrast. Counting it can only
      # block a ROBUST call, which is the direction an interpretive choice here should err.
      b <- balance_table(cc, c(base_cov, lv, other, "score_MEL", "LPA_log"), c(base_cov, lv), W$trunc, W$raw,
                         not_counted = other)
      b[, `:=`(form = form, spec = sp)]
      bal[[length(bal) + 1L]] <- b
      mx <- b[counts_for_rule == TRUE, max(abs(smd_weighted))]
      ok <- b[counts_for_rule == TRUE, all(pass)]
      diag[[length(diag) + 1L]] <- data.table(form = form, spec = sp, c_statistic = cstat(e, cc$statin),
        max_abs_smd = mx, weight_max = max(W$trunc), ess_users = ess(W$trunc[cc$statin == 1L]),
        ess_nonusers = ess(W$trunc[cc$statin == 0L]), n = sup(nrow(cc)), users = sup(sum(cc$statin == 1L)),
        note = paste(setdiff(g$msgs, ""), collapse = "; "))
      if (is.null(picked) || (ok && !picked$ok)) picked <- list(spec = sp, ok = ok, mx = mx, W = W)
      if (ok) break
    }
    if (is.null(picked)) { chosen[[form]] <- data.table(form = form, spec = NA_character_, balanced = FALSE, max_abs_smd = NA_real_); next }
    chosen[[form]] <- data.table(form = form, spec = picked$spec, balanced = picked$ok, max_abs_smd = picked$mx)
    run_outcome(cc, picked$W$trunc, sprintf("IPTW, %s LDL-C", form),
                sprintf("%s propensity; weights truncated at the 99th centile; balanced=%s", picked$spec, picked$ok))
  }

  # Sensitivity only, never read by the verdict: the Tudball-style lasso propensity
  if (use_glmnet && requireNamespace("glmnet", quietly = TRUE)) {
    X <- as.matrix(cc[, c(base_cov, "LDL"), with = FALSE])
    g <- guarded(glmnet::cv.glmnet(X, cc$statin, family = "binomial", alpha = 1, nfolds = 10))
    if (!is.null(g$value)) {
      e <- as.numeric(predict(g$value, X, s = "lambda.min", type = "response")); W <- stabilized(e, cc$statin)
      b <- balance_table(cc, c(base_cov, "LDL", "LDL_corr", "score_MEL", "LPA_log"), c(base_cov, "LDL"), W$trunc, W$raw,
                         not_counted = "LDL_corr")
      b[, `:=`(form = "raw", spec = "lasso_sensitivity")]; bal[[length(bal) + 1L]] <- b
      run_outcome(cc, W$trunc, "IPTW lasso, raw LDL-C (sensitivity)", "not used by the verdict")
    }
  } else message("[B] glmnet unavailable: lasso sensitivity skipped")

  list(balance = rbindlist(bal, fill = TRUE), est = rbindlist(est), diag = rbindlist(diag, fill = TRUE),
       chosen = rbindlist(chosen))
}

# ---- (5) selection-only simulation -------------------------------------------
# WHAT THIS SIMULATION IS
#   A bound on a mechanism. It asks: if statin selection were the only thing going on,
#   and the drug multiplied everyone's hazard by the same number whatever their score,
#   what ratio of HRs would the published model report? It does NOT estimate a causal
#   effect and it cannot show that selection IS the explanation; reproducing 0.93 only
#   shows that 0.93 does not require a true interaction.
#
# ASSUMPTIONS, each of which could be wrong
#   - Scores, age and sex: resampled jointly from the real cohort, so the skewed Lp(a)
#     score keeps its shape. Sampling from survivors and then depleting again double
#     counts slightly; the effect on a standardized score is negligible.
#   - Pre-treatment LDL-C z = r_L G_L + r_P G_P + sigma e, with r and sigma taken from
#     /0.7-corrected LDL-C, the only observable proxy for untreated LDL-C.
#   - U is a single unmeasured N(0,1) risk factor, independent of the scores, acting
#     with the same log coefficient s on statin log-odds and on log hazard. Smoking,
#     blood pressure, diabetes and family history are all folded into U.
#   - Selection is on LDL-C, U, age, sex, and a direct Lp(a) score term. The G_L route is
#     wholly through LDL-C. Coefficients are fitted to observed OR(statin|score), the
#     statin prevalence, and the age and sex coefficients, separately in each grid cell.
#   - Hazard: exponential, log-linear in G_L, G_P, LDL-C z, age, sex and U, fitted so the
#     simulated non-user Cox coefficients match the observed ones and total events match.
#     Exponential is harmless since Cox ignores the baseline shape, except through the
#     interplay with censoring.
#   - Measured LDL-C in simulated users = untreated x 0.7 x exp(N(0, tau)), tau from
#     --statin_response_sd (default 0.2: an assumption, about the size of published
#     between-person variation in statin LDL-C response, not fitted here). Without this,
#     dividing by 0.7 would hand the simulated IPTW the true untreated value, which no real
#     analysis has, and make it look cleaner than the real one.
#   - U is wholly absent from the simulated propensity, whereas the real propensity contains
#     BMI, smoking, diabetes and hypertension. Simulated IPTW therefore leaves at least as
#     much residual selection as a real one with the same U; for a ROBUST claim that errs
#     on the side of not claiming it. Measurement error in LDL-C is not modelled.
#   - Statin: fixed from baseline, multiplicative, identical at every score. There is no
#     initiation or stopping during follow-up. Real non-users start statins later, which
#     would change the non-user slope, and that is not modelled.
#   - Prevalent disease: the pre-baseline hazard shares the same linear predictor and
#     removes the observed fraction before statin assignment. This is the index-event part
#     of selection. Statin effects before baseline are not modelled.
#   - Censoring: resampled from follow-up of real non-cases, independent of everything.
#     PCs are left out of the simulated outcome model since nothing there depends on them.
#   The grid over U strength is the honest way to show how much depends on the unknown.
sim_targets <- function(d, primary) {
  cc <- d[!is.na(LDL_corr)]
  if (nrow(cc) < 1000L || sum(d$statin == 1L) < MIN_EVENTS || sum(d$event) < MIN_EVENTS) return(NULL)
  mu <- mean(cc$LDL_corr); sdl <- sd(cc$LDL_corr)
  cc[, Lz := (LDL_corr - mu) / sdl]; cc[, Lz_nc := (LDL - mu) / sdl]
  r <- summary(lm(Lz ~ score_LDL + score_LPA, data = cc))
  amean <- mean(d$age)
  lg <- glm(statin ~ score_LDL + score_LPA + I(age - amean) + sex, data = d, family = binomial())
  ch <- coxph(Surv(fu, event) ~ score_LDL + score_LPA + Lz_nc + I(age - amean) + sex, data = cc[statin == 0L])
  list(N = nrow(d), events = sum(d$event), p_statin = mean(d$statin), prev_frac = primary$prev_frac,
       mu = mu, sdl = sdl, r_L = coef(r)["score_LDL", 1], r_P = coef(r)["score_LPA", 1], sigma_e = r$sigma,
       t_logit = unname(coef(lg)), t_cox = unname(coef(ch)),
       pool = data.frame(G_L = d$score_LDL, G_P = d$score_LPA, age = d$age - amean, sex = d$sex),
       fu_pool = d$fu[d$event == 0L], obs = primary)
}

draw_base <- function(cfg, seed) {
  set.seed(seed)
  npop <- ceiling(cfg$N / (1 - cfg$prev_frac) * 1.02)
  # column by column: `[.data.frame` with repeated rows builds unique row names, seconds per draw
  idx <- sample.int(nrow(cfg$pool), npop, TRUE); p <- lapply(cfg$pool, function(v) v[idx])
  list(G_L = p$G_L, G_P = p$G_P, age = p$age, sex = p$sex, U = rnorm(npop), eps = rnorm(npop),
       uPrev = runif(npop), uA = runif(npop), uT = runif(npop), uC = runif(npop))
}

# Uniforms are fixed in draw_base, so during calibration a parameter change moves the outcome
# smoothly (common random numbers), and fixed-point updates converge rather than chase noise.
realize <- function(b, par, cfg, s, log_hr, solve_kappa = FALSE) {
  Lz <- cfg$r_L * b$G_L + cfg$r_P * b$G_P + cfg$sigma_e * b$eps
  lp0 <- par$b_GL * b$G_L + par$b_GP * b$G_P + par$b_L * Lz + par$b_age * b$age + par$b_sex * b$sex + s * b$U
  if (solve_kappa) {
    if (cfg$prev_frac <= 0) par$lkappa <- -Inf
    else par$lkappa <- uniroot(function(lk) mean(-expm1(-exp(lk + lp0))) - cfg$prev_frac, c(-40, 10))$root
  }
  keep <- b$uPrev > -expm1(-exp(par$lkappa + lp0))
  ii <- which(keep); if (length(ii) > cfg$N) ii <- ii[seq_len(cfg$N)]
  eta <- par$g0 + par$gL * Lz[ii] + par$gP * b$G_P[ii] + par$g_age * b$age[ii] + par$g_sex * b$sex[ii] + s * b$U[ii]
  A <- as.integer(b$uA[ii] < plogis(eta))
  Ti <- -log(b$uT[ii]) / exp(par$lh + lp0[ii] + log_hr * A)
  C <- cfg$fu_pool[pmax(1L, ceiling(b$uC[ii] * length(cfg$fu_pool)))]
  x <- data.frame(G_L = as.numeric(scale(b$G_L[ii])), G_P = as.numeric(scale(b$G_P[ii])), Lz = Lz[ii],
                  age = b$age[ii], sex = b$sex[ii], A = A, time = pmin(Ti, C), E = as.integer(Ti <= C))
  list(x = x, par = par)
}

calibrate_cell <- function(cfg, s, log_hr, seed, max_iter = 25L, tol = 0.005) {
  tl <- cfg$t_logit; tc <- cfg$t_cox
  par <- list(g0 = tl[1], gL = tl[2] / cfg$r_L, gP = tl[3], g_age = tl[4], g_sex = tl[5],
              b_GL = tc[1], b_GP = tc[2], b_L = tc[3], b_age = tc[4], b_sex = tc[5],
              lh = log(cfg$events / (cfg$N * mean(cfg$fu_pool))), lkappa = NA_real_)
  b <- draw_base(cfg, seed); conv <- FALSE; res <- rep(NA_real_, 11L); it <- 0L
  for (it in seq_len(max_iter)) {
    z <- guarded({
      r <- realize(b, par, cfg, s, log_hr, solve_kappa = TRUE); x <- r$x
      lg <- unname(glm.fit(cbind(1, x$G_L, x$G_P, x$age, x$sex), x$A, family = binomial())$coefficients)
      nu <- x$A == 0L
      cx <- coef(coxph(Surv(time, E) ~ G_L + G_P + Lz + age + sex, data = x[nu, ]))
      list(r = r, x = x, lg = lg, cx = cx)
    })
    if (is.null(z$value) || any(is.na(z$value$lg)) || any(is.na(z$value$cx))) { conv <- FALSE; break }
    par <- z$value$r$par; lg <- z$value$lg; cx <- z$value$cx; x <- z$value$x
    dl <- tl - lg; dc <- tc - unname(cx); de <- log(cfg$events / max(1, sum(x$E)))
    res <- c(dl, dc, de)
    if (all(abs(res) < tol)) { conv <- TRUE; break }
    par$g0 <- par$g0 + dl[1]; par$gL <- par$gL + dl[2] / cfg$r_L; par$gP <- par$gP + dl[3]
    par$g_age <- par$g_age + dl[4]; par$g_sex <- par$g_sex + dl[5]
    par$b_GL <- par$b_GL + dc[1]; par$b_GP <- par$b_GP + dc[2]; par$b_L <- par$b_L + dc[3]
    par$b_age <- par$b_age + dc[4]; par$b_sex <- par$b_sex + dc[5]; par$lh <- par$lh + de
  }
  list(converged = conv, iterations = it, par = par, max_residual = max(abs(res)), s = s, log_hr = log_hr)
}

# One replicate: the published model, on a cohort with no true interaction by construction.
sim_rep <- function(cell, cfg, seed, secondary) {
  x <- realize(draw_base(cfg, seed), cell$par, cfg, cell$s, cell$log_hr)$x
  out <- c(events = sum(x$E), user_events = sum(x$E[x$A == 1L]), users = sum(x$A))
  grab <- function(f, dd, w = NULL) {
    m <- if (is.null(w)) coxph(f, data = dd) else coxph(f, data = cbind(dd, .w = w), weights = .w)
    cf <- coef(m); if (anyNA(cf)) stop("non-estimable"); cf
  }
  for (S in c("L", "P")) {
    x$G <- if (S == "L") x$G_L else x$G_P
    cf <- grab(Surv(time, E) ~ G * A + age + sex, x)
    out[paste0("primary_", S)] <- cf[["G:A"]]; out[paste0("primary_main_", S)] <- cf[["G"]]; out[paste0("primary_statin_", S)] <- cf[["A"]]
  }
  if (secondary) {
    L <- cfg$mu + cfg$sdl * x$Lz
    tau <- if (is.null(cfg$resp_sd)) 0 else cfg$resp_sd
    raw <- ifelse(x$A == 1L, 0.7 * exp(rnorm(nrow(x), 0, tau)), 1) * L
    meas <- list(raw = raw, corrected = ifelse(x$A == 1L, raw / 0.7, raw))
    for (form in names(meas)) {
      x$M <- meas[[form]]
      ps <- glm.fit(cbind(1, x$G_L, x$G_P, x$age, x$sex, x$M), x$A, family = binomial())$fitted.values
      w <- stabilized(ps, x$A)$trunc
      for (S in c("L", "P")) {
        x$G <- if (S == "L") x$G_L else x$G_P
        out[sprintf("ldladj_%s_%s", form, S)] <- grab(Surv(time, E) ~ G * A + age + sex + M, x)[["G:A"]]
        out[sprintf("iptw_%s_%s", form, S)]   <- grab(Surv(time, E) ~ G * A + age + sex, x, w)[["G:A"]]
      }
    }
  }
  out
}

summ_roh <- function(v) {
  v <- v[is.finite(v)]
  list(n_ok = length(v), mean_roh = exp(mean(v)), q025 = exp(quantile(v, 0.025, names = FALSE)),
       q975 = exp(quantile(v, 0.975, names = FALSE)))
}

run_sim <- function(cfg, opt) {
  hrs <- as.numeric(strsplit(as.character(opt$grid_hr), ",")[[1]]); us <- as.numeric(strsplit(as.character(opt$grid_u), ",")[[1]])
  cells <- CJ(statin_hr = hrs, u_strength = us, sorted = FALSE); cells[, cell := .I]
  ckdir <- file.path(opt$outdir, "sim_cells"); dir.create(ckdir, recursive = TRUE, showWarnings = FALSE)
  key <- list(N = cfg$N, events = cfg$events, t_logit = cfg$t_logit, t_cox = cfg$t_cox, r = c(cfg$r_L, cfg$r_P, cfg$sigma_e),
              n_reps = opt$n_reps, n_sec = opt$n_reps_secondary, prev = cfg$prev_frac, resp_sd = cfg$resp_sd,
              code = unname(tools::md5sum(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
  setDTthreads(1L)
  run_one_cell <- function(i) {
    ck <- file.path(ckdir, sprintf("cell_%03d.rds", i))
    if (file.exists(ck)) {
      sv <- readRDS(ck)
      # A cell short of replicates (a worker killed for memory, say) is rerun, not trusted
      n_good <- sum(vapply(sv$res$reps, function(v) length(v) > 1L, TRUE))
      # cal NULL is an error (memory, say), not a cell that failed to calibrate, so it reruns.
      complete <- !is.null(sv$res$cal) && (is.null(sv$res$reps) || n_good >= ceiling(0.9 * opt$n_reps))
      if (identical(sv$key, key) && isTRUE(all.equal(sv$cell_id, as.list(cells[i]))) && complete) return(sv$res)
    }
    t0 <- Sys.time()
    cal <- guarded(calibrate_cell(cfg, cells$u_strength[i], log(cells$statin_hr[i]), seed = 42L + 1000000L * i))$value
    if (is.null(cal) || !isTRUE(cal$converged)) {
      res <- list(cal = cal, reps = NULL)
    } else {
      reps <- mclapply(seq_len(opt$n_reps), function(r) {
        v <- guarded(sim_rep(cal, cfg, 42L + 1000000L * i + r, secondary = r <= opt$n_reps_secondary))$value
        if (is.null(v)) NA_real_ else v
      }, mc.cores = opt$cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
      res <- list(cal = cal, reps = reps)
    }
    saveRDS(list(key = key, cell_id = as.list(cells[i]), res = res), ck)     # replicate-level summaries only
    message(sprintf("[sim] cell %d/%d HR=%.2f U=%.2f converged=%s (%.1f min)", i, nrow(cells), cells$statin_hr[i],
                    cells$u_strength[i], isTRUE(cal$converged), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    res
  }
  raw <- lapply(seq_len(nrow(cells)), run_one_cell)

  grid <- list(); calib <- list()
  need <- ceiling(0.9 * opt$n_reps); need2 <- ceiling(0.9 * opt$n_reps_secondary)
  for (i in seq_len(nrow(cells))) {
    cal <- raw[[i]]$cal; reps <- raw[[i]]$reps
    calib[[i]] <- cbind(cells[i], data.table(converged = isTRUE(cal$converged), iterations = if (is.null(cal)) NA_integer_ else cal$iterations,
      max_residual = if (is.null(cal)) NA_real_ else cal$max_residual), if (!is.null(cal)) as.data.table(cal$par[c("g0","gL","gP","g_age","g_sex","b_GL","b_GP","b_L","b_age","b_sex","lh","lkappa")]))
    m <- NULL
    if (!is.null(reps)) { ok <- reps[vapply(reps, function(v) length(v) > 1L, TRUE)]; if (length(ok)) m <- rbindlist(lapply(ok, function(v) as.list(v)), fill = TRUE) }
    rows <- list()
    add_row <- function(est, sc, v, thresh, sh = NA_real_, ev = NA_real_, uev = NA_real_, main = NA_real_) {
      s <- summ_roh(v); bad <- s$n_ok < thresh
      rows[[length(rows) + 1L]] <<- cbind(cells[i, .(statin_hr, u_strength)], data.table(estimator = est, score = sc, n_reps_ok = s$n_ok,
        mean_roh = if (bad) NA_real_ else s$mean_roh, q025 = if (bad) NA_real_ else s$q025, q975 = if (bad) NA_real_ else s$q975,
        mean_statin_main_hr = sh, mean_nonuser_score_hr = main, mean_events = ev, mean_user_events = uev))
    }
    if (is.null(m) || nrow(m) < need) {
      for (sc in c("LDL", "LPA", "LPA_minus_LDL")) rows[[length(rows) + 1L]] <- cbind(cells[i, .(statin_hr, u_strength)], data.table(estimator = "primary", score = sc, n_reps_ok = if (is.null(m)) 0L else nrow(m), mean_roh = NA_real_, q025 = NA_real_, q975 = NA_real_, mean_statin_main_hr = NA_real_, mean_nonuser_score_hr = NA_real_, mean_events = NA_real_, mean_user_events = NA_real_))
    } else {
      add_row("primary", "LDL", m$primary_L, need, exp(mean(m$primary_statin_L, na.rm = TRUE)), mean(m$events), mean(m$user_events), exp(mean(m$primary_main_L, na.rm = TRUE)))
      add_row("primary", "LPA", m$primary_P, need, exp(mean(m$primary_statin_P, na.rm = TRUE)), mean(m$events), mean(m$user_events), exp(mean(m$primary_main_P, na.rm = TRUE)))
      add_row("primary", "LPA_minus_LDL", m$primary_P - m$primary_L, need)
      if (opt$n_reps_secondary > 0L) for (form in c("raw", "corrected")) for (S in c("L", "P")) {
        sc <- if (S == "L") "LDL" else "LPA"
        cn1 <- sprintf("ldladj_%s_%s", form, S); cn2 <- sprintf("iptw_%s_%s", form, S)
        add_row(sprintf("ldl_adjusted_%s", form), sc, if (cn1 %in% names(m)) m[[cn1]] else NA_real_, need2)
        add_row(sprintf("iptw_%s", form), sc, if (cn2 %in% names(m)) m[[cn2]] else NA_real_, need2)
      }
    }
    grid[[i]] <- rbindlist(rows, fill = TRUE)
  }
  list(grid = rbindlist(grid, fill = TRUE), calibration = rbindlist(calib, fill = TRUE))
}

# Does selection alone reach the observed attenuation?
#   in_range   the observed value lies inside the cell's 2.5-97.5% Monte Carlo range
#   reproduced the cell gets to an attenuation at least that large (q025 <= observed),
#              or its mean falls inside the observed CI. Overshooting counts: a mechanism
#              able to make 0.85 is able to make 0.93.
# Any cell in the prespecified grid is enough. A ROBUST claim is the strong one, so it
# carries the burden; cells whose statin main effect misses the observed one are
# flagged in the output, not dropped from the rule.
sim_compare <- function(grid, est_name, score_name, obs) {
  g <- grid[estimator == est_name & score == score_name]
  if (!nrow(g) || any(is.na(obs))) return(list(reproduced = NA, in_range = NA, incomplete = TRUE, k = NA_integer_, n_cells = nrow(g), rows = g))
  g <- copy(g)
  g[, in_range := q025 <= obs[1] & obs[1] <= q975]
  g[, reproduced := q025 <= obs[1] | mean_roh <= obs[3]]
  incomplete <- anyNA(g$mean_roh)
  any_rep <- any(g$reproduced, na.rm = TRUE)
  list(reproduced = if (any_rep) TRUE else if (incomplete) NA else FALSE,
       in_range = if (any(g$in_range, na.rm = TRUE)) TRUE else if (incomplete) NA else FALSE,
       incomplete = incomplete, k = sum(g$reproduced, na.rm = TRUE), n_cells = nrow(g), rows = g)
}

# The asymmetry is the paper's claim, so the informative question is whether ONE selection
# scenario gives both observed values at once. A grid that reaches 0.93 in one cell and 0.99
# in another shows less. Descriptive; the preregistered rule reads the LDL-C row alone.
sim_joint <- function(grid, ob_L, ob_P, ob_main_L, ob_main_P, ob_statin) {
  w <- function(sc, pre) { g <- grid[estimator == "primary" & score == sc,
                                     .(statin_hr, u_strength, m = mean_roh, lo = q025, hi = q975, sh = mean_statin_main_hr, mh = mean_nonuser_score_hr)]
    setnames(g, c("m", "lo", "hi", "mh"), paste0(pre, c("_mean_roh", "_q025", "_q975", "_nonuser_score_hr"))); g }
  j <- merge(w("LDL", "ldl"), w("LPA", "lpa")[, !"sh"], by = c("statin_hr", "u_strength"))
  j[, `:=`(ldl_in_range = ldl_q025 <= ob_L & ob_L <= ldl_q975, lpa_in_range = lpa_q025 <= ob_P & ob_P <= lpa_q975)]
  j[, both_in_range := ldl_in_range & lpa_in_range]
  # Would a reader accept this cell? Same observable the design has: the statin main effect.
  j[, statin_hr_within_10pct_of_observed := abs(log(sh / ob_statin)) < log(1.1)]
  j[, `:=`(observed_ldl_roh = ob_L, observed_lpa_roh = ob_P, observed_statin_main_hr = ob_statin,
           observed_ldl_nonuser_hr = ob_main_L, observed_lpa_nonuser_hr = ob_main_P)]
  setnames(j, "sh", "sim_statin_main_hr"); j[]
}

# For each true statin HR, the U strength at which the simulated statin main effect
# equals the observed one (log-linear interpolation) and the ratio of HRs there.
# Descriptive: narrows the grid using the one extra observable the design has.
sim_interpolate <- function(grid, obs_statin) {
  rbindlist(lapply(unique(grid$statin_hr), function(h) rbindlist(lapply(c("LDL", "LPA"), function(sc) {
    g <- grid[statin_hr == h & estimator == "primary" & score == sc & !is.na(mean_roh)][order(u_strength)]
    star <- NA_real_; roh <- NA_real_
    if (nrow(g) >= 2L && min(g$mean_statin_main_hr) <= obs_statin && obs_statin <= max(g$mean_statin_main_hr)) {
      lg <- log(g$mean_statin_main_hr)
      for (k in seq_len(nrow(g) - 1L)) if ((lg[k] - log(obs_statin)) * (lg[k + 1L] - log(obs_statin)) <= 0) {
        fr <- (log(obs_statin) - lg[k]) / (lg[k + 1L] - lg[k]); if (!is.finite(fr)) fr <- 0
        star <- g$u_strength[k] + fr * (g$u_strength[k + 1L] - g$u_strength[k])
        roh <- exp(log(g$mean_roh[k]) + fr * (log(g$mean_roh[k + 1L]) - log(g$mean_roh[k]))); break
      }
    }
    data.table(statin_hr = h, score = sc, u_strength_matching_observed_statin_hr = star, interpolated_mean_roh = roh,
               observed_statin_hr = obs_statin, note = if (is.na(star)) "observed statin HR outside this row of the grid" else "")
  }))))
}

# ---- verdict, PREREGISTERED_READOUT.md section B ---------------------------------
# Pure function so every branch can be tested with a table of inputs.
# WHICH IPTW ARMS COUNT. The readout says "IPTW including measured LDL-C" and "every weighted
# SMD < 0.1", and does not say which form of LDL-C. Raw LDL-C in a user is on treatment: the
# drug itself moves it by about 1 SD, so it separates users from non-users for a reason that is
# not indication, and top-1% truncated weights cannot close that gap (SMD -0.54 in the selftest,
# whose only selection is on untreated LDL-C and age; the Tudball port's arm was uninformative
# too). Demanding balance in that arm would make ROBUST unreachable whatever the truth, a verdict
# fixed in advance rather than a strict one. So, per arm: an unbalanced arm is uninformative and,
# as the readout says, supports neither reading; every balanced arm must lie below 1; at least one
# must be balanced. Balance is judged before any outcome is touched (iptw_block), so this is not a
# choice of the arm with the better answer. require_both = TRUE gives the stricter reading and
# is written to verdict.csv alongside, so that any dependence on this choice is visible.
iptw_arm_state <- function(r) {
  if (is.na(r$estimate) || is.na(r$lo95) || is.na(r$hi95)) return("UNAVAILABLE")
  if (!isTRUE(r$balanced)) return("UNINFORMATIVE")
  # The claim is an ATTENUATION, so the whole interval has to lie below 1. A ratio above 1 with
  # its CI excluding 1 would be a different finding and must not read as this one.
  if (r$hi95 < 1) "CI_EXCLUDES_1" else if (r$lo95 > 1) "CI_ABOVE_1" else "CI_INCLUDES_1"
}
# Recorded in PREREGISTERED_READOUT.md, Amendments. The /0.7 arm has to be balanced: it is the one
# balanced on the indication. A raw arm that balances can block ROBUST and can never license it,
# because weights balancing on-treatment LDL-C leave the /0.7 value about 1.3 SD apart.
combine_arm_states <- function(st, forms, require_both) {
  if (length(st) < 2L || "UNAVAILABLE" %in% st) return("UNAVAILABLE")
  if (require_both && "UNINFORMATIVE" %in% st) return("UNINFORMATIVE")
  if (st[forms == "corrected"] == "UNINFORMATIVE") return("UNINFORMATIVE")
  use <- setdiff(st, "UNINFORMATIVE")
  if (!length(use)) return("UNINFORMATIVE")
  if (all(use == "CI_EXCLUDES_1")) "CI_EXCLUDES_1" else if ("CI_ABOVE_1" %in% use) "CI_ABOVE_1" else "CI_INCLUDES_1"
}

decide_verdict <- function(primary_ldl, primary_lpa, primary_mel, iptw, sim_primary, sim_iptw, require_both = FALSE) {
  # iptw: data.table(form, score, estimate, lo95, hi95, balanced, max_abs_smd); one row per form x score
  excl1 <- function(lo, hi) !is.na(lo) && !is.na(hi) && (lo > 1 || hi < 1)
  fmt <- function(v) if (anyNA(v[1:3])) "NA" else sprintf("%.3f (%.3f-%.3f)", v[1], v[2], v[3])
  # --- negative control. A balanced IPTW arm counts too; an unbalanced one is uninformative.
  neg_rows <- iptw[score == "MEL" & balanced == TRUE]
  neg_bad <- excl1(primary_mel[2], primary_mel[3]) || (nrow(neg_rows) > 0L && any(mapply(excl1, neg_rows$lo95, neg_rows$hi95)))
  neg_avail <- !anyNA(primary_mel[1:3])
  neg <- if (!neg_avail) "UNAVAILABLE" else if (neg_bad) "CONFOUNDED" else "PASS"
  # --- IPTW, LDL-C, arm by arm
  L <- iptw[score == "LDL"]
  forms_ok <- nrow(L) >= 2L && all(c("raw", "corrected") %in% L$form)
  arm_st <- if (forms_ok) vapply(seq_len(nrow(L)), function(i) iptw_arm_state(L[i]), "") else character()
  iptw_state <- if (forms_ok) combine_arm_states(arm_st, L$form, require_both) else "UNAVAILABLE"
  primary_below <- !anyNA(primary_ldl[1:3]) && primary_ldl[3] < 1
  # --- simulation. NA means incomplete; incomplete can never license ROBUST.
  sim_ok <- isFALSE(sim_primary$reproduced) && isFALSE(sim_iptw$reproduced)
  reasons <- character()
  if (neg == "CONFOUNDED") reasons <- c(reasons, "negative-control interaction CI excludes 1: every interaction in this section is treated as confounded")
  if (neg == "UNAVAILABLE") reasons <- c(reasons, "negative control unavailable")
  if (iptw_state == "UNINFORMATIVE") reasons <- c(reasons, sprintf("weighted SMD >= %.1f (%s): IPTW is uninformative, not support for either reading", SMD_LIMIT, paste(sprintf("%s max %.3f", L$form, L$max_abs_smd), collapse = ", ")))
  if (iptw_state == "UNAVAILABLE") reasons <- c(reasons, "IPTW estimate unavailable")
  if (iptw_state == "CI_INCLUDES_1") reasons <- c(reasons, "IPTW ratio of HRs CI includes 1")
  if (iptw_state == "CI_ABOVE_1") reasons <- c(reasons, "IPTW ratio of HRs lies above 1: not an attenuation")
  if (!primary_below) reasons <- c(reasons, "primary ratio of HRs CI does not lie below 1")
  if (isTRUE(sim_primary$reproduced)) reasons <- c(reasons, "selection-only simulation reproduces the observed attenuation")
  else if (is.na(sim_primary$reproduced)) reasons <- c(reasons, "selection-only simulation incomplete")
  if (isTRUE(sim_iptw$reproduced)) reasons <- c(reasons, "selection-only simulation reproduces the attenuation left after IPTW")
  else if (is.na(sim_iptw$reproduced) && iptw_state == "CI_EXCLUDES_1") reasons <- c(reasons, "simulated IPTW comparison incomplete")
  robust <- neg == "PASS" && primary_below && iptw_state == "CI_EXCLUDES_1" && sim_ok
  ldl <- if (robust) "ROBUST" else "NOT SEPARABLE FROM SELECTION"
  # --- Lp(a): only ever a bounded null, with its interval. Confounded negative control voids it.
  lpa_null <- !anyNA(primary_lpa[1:3]) && !excl1(primary_lpa[2], primary_lpa[3]) && neg != "CONFOUNDED"
  lpa_iptw <- iptw[score == "LPA" & balanced == TRUE]
  lpa_iptw_ok <- nrow(lpa_iptw) == 0L || !any(mapply(excl1, lpa_iptw$lo95, lpa_iptw$hi95))
  lpa <- if (lpa_null && lpa_iptw_ok) sprintf("BOUNDED NULL: ratio of HRs %s; excludes attenuation beyond %.1f%%", fmt(primary_lpa), 100 * (1 - primary_lpa[2]))
         else "NO NULL CLAIM"
  if (lpa_null && lpa_iptw_ok && nrow(lpa_iptw) == 0L) lpa <- paste0(lpa, "; IPTW arm unbalanced so it adds no support")
  data.table(
    item = c("LDL_C_ATTENUATION", "LDL_C_REASONS", "IPTW_STATE", "IPTW_ARMS", "NEGATIVE_CONTROL", "LPA_NO_MODIFICATION", "CONSEQUENCE"),
    verdict = c(ldl, if (length(reasons)) paste(reasons, collapse = " | ") else "all preregistered conditions met",
                iptw_state,
                if (forms_ok) paste(sprintf("%s: %s, max SMD %.3f, ratio of HRs %s", L$form, arm_st, L$max_abs_smd, vapply(seq_len(nrow(L)), function(i) fmt(c(L$estimate[i], L$lo95[i], L$hi95[i])), "")), collapse = " | ") else "NA",
                neg, lpa,
                if (robust) "LDL-C attenuation may be called robust to measured selection; it remains an association with baseline statin use"
                else "The asymmetry leaves the abstract as a finding and stays as a described association"))
}

# ---- orchestration -----------------------------------------------------------
run_analysis <- function(opt) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  setDTthreads(opt$cores)
  pr <- prep_cohort(opt); d <- pr$d
  # Positive control is first: nothing but its own table is written if it fails
  pc <- positive_control(d, pr$counts, opt)
  pcs <- grep("^PC", names(d), value = TRUE)[1:10]; PCS <- paste(pcs, collapse = " + ")

  cor_or_na <- function(a, b) { ok <- !is.na(a) & !is.na(b); if (sum(ok) < 10L) NA_real_ else cor(a[ok], b[ok]) }
  tails <- function(v, code) if (all(is.na(v))) NA_integer_ else sup(sum(v %in% code))
  n_arr <- guarded(statin_array_diagnostic(d, opt$raw_dir))$value
  summary_tab <- data.table(
    quantity = c("n", "events", "statin users", "events in statin users", "median follow-up (y)",
                 "prevalent excluded", "no follow-up excluded", "missing age/sex/PCs excluded (coxph listwise)", "measured LDL-C n", "LDL-C missing n",
                 "measured Lp(a) n", "Lp(a) missing n", "Lp(a) too low (30796=2,4)", "Lp(a) too high (30796=3,5)",
                 "LDL-C too low (30786=2,4)", "LDL-C too high (30786=3,5)", "melanoma score n",
                 "r(LDL-C /0.7, LDL-C score), published 0.36", "r(Lp(a), Lp(a) score), published 0.68",
                 "non-users with statin code 1 in a later 6153/6177 array"),
    value = c(sup(nrow(d)), sup(sum(d$event)), sup(sum(d$statin == 1L)), sup(sum(d$event[d$statin == 1L])), median(d$fu),
              sup(pr$counts$prevalent), sup(pr$counts$no_followup), sup(pr$counts$covariate_missing), sup(sum(!is.na(d$LDL))), sup(sum(is.na(d$LDL))),
              sup(sum(!is.na(d$LPA))), sup(sum(is.na(d$LPA))), tails(d$LPA_rep, c(2L, 4L)), tails(d$LPA_rep, c(3L, 5L)),
              tails(d$LDL_rep, c(2L, 4L)), tails(d$LDL_rep, c(3L, 5L)), sup(sum(!is.na(d$score_MEL))),
              cor_or_na(d$LDL_corr, d$score_LDL), cor_or_na(d$LPA, d$score_LPA), if (is.null(n_arr)) NA_real_ else sup(n_arr)))
  wcsv(summary_tab, opt$outdir, "cohort_summary.csv"); print(summary_tab)

  message("[B] (1) selection into statin use")
  sel <- selection_tables(d, PCS)
  wcsv(sel$logit, opt$outdir, "selection_logistic.csv"); wcsv(sel$linear, opt$outdir, "measured_by_statin_stratum.csv")
  print(sel$logit[, .(model, term, OR = round(estimate, 3), lo = round(lo95, 3), hi = round(hi95, 3), P = signif(P, 3), n)])

  message("[B] (2) Cox with measured biomarkers")
  cm <- measured_cox(d, PCS)
  wcsv(cm, opt$outdir, "cox_measured_adjusted.csv")
  print(cm[, .(model, term, HR = round(estimate, 3), lo = round(lo95, 3), hi = round(hi95, 3), n)])

  message("[B] (3) stabilized truncated IPTW")
  ip <- iptw_block(d, pcs, PCS, use_glmnet = opt$use_glmnet)
  wcsv(ip$balance, opt$outdir, "iptw_balance.csv"); wcsv(ip$est, opt$outdir, "iptw_estimates.csv"); wcsv(ip$diag, opt$outdir, "iptw_diagnostics.csv")
  wcsv(ip$chosen, opt$outdir, "iptw_chosen_propensity.csv")
  print(ip$chosen)
  print(ip$est[, .(model, term, HR = round(estimate, 3), lo = round(lo95, 3), hi = round(hi95, 3), n, note)])

  message("[B] (4) negative control: melanoma score x statin")
  neg <- fit_cox(sprintf("Surv(fu, event) ~ score_MEL * statin + age + sex + %s", PCS), d[!is.na(score_MEL)],
                 "melanoma, primary model", "score_MEL:statin", "score_MEL")
  wcsv(rbind(neg, ip$est[term == "score_MEL:statin"]), opt$outdir, "negative_control.csv"); print(neg)

  iptw_long <- function() {
    if (!nrow(ip$est)) return(data.table(form = character(), score = character(), estimate = numeric(), lo95 = numeric(), hi95 = numeric(), balanced = logical(), max_abs_smd = numeric()))
    rbindlist(lapply(c("raw", "corrected"), function(fm) {
      lbl <- sprintf("IPTW, %s LDL-C", fm); ch <- ip$chosen[form == fm]
      rbindlist(lapply(c("LDL", "LPA", "MEL"), function(sc) {
        r <- ip$est[model == lbl & term == sprintf("score_%s:statin", sc)]
        if (!nrow(r)) r <- data.table(estimate = NA_real_, lo95 = NA_real_, hi95 = NA_real_)
        data.table(form = fm, score = sc, estimate = r$estimate, lo95 = r$lo95, hi95 = r$hi95,
                   balanced = isTRUE(ch$balanced), max_abs_smd = if (nrow(ch)) ch$max_abs_smd else NA_real_)
      }))
    }))
  }
  IP <- iptw_long()
  wcsv(IP, opt$outdir, "iptw_verdict_inputs.csv")

  prim_L <- c(pc$ldl$int); prim_P <- c(pc$lpa$int); prim_M <- c(neg$estimate, neg$lo95, neg$hi95)
  grid <- NULL; sim_L <- list(reproduced = NA, in_range = NA, incomplete = TRUE, k = NA, n_cells = 0L); sim_P <- sim_L; sim_I <- sim_L; sim_I_all <- sim_L
  joint <- NULL; interp <- NULL; sim_L_match <- NULL
  if (!isTRUE(opt$skip_sim)) {
    message("[B] (5) selection-only simulation")
    primary <- list(ldl = prim_L, lpa = prim_P, statin = c(pc$ldl$statin), prev_frac = pr$counts$prevalent / (pr$counts$prevalent + nrow(d)))
    cfg <- sim_targets(d, primary)
    if (!is.null(cfg)) cfg$resp_sd <- opt$statin_response_sd
    if (is.null(cfg) || !is.finite(cfg$r_L) || abs(cfg$r_L) < 0.02) {
      warning("selection-only simulation not run: targets degenerate (LDL-C score unrelated to LDL-C, or too few events)", call. = FALSE)
    } else {
      pr$d <- NULL; rm(d); gc(verbose = FALSE)   # the workers only need cfg; don't fork the whole cohort
      s <- run_sim(cfg, opt); grid <- s$grid
      wcsv(s$calibration, opt$outdir, "sim_calibration.csv")
      grid[, `:=`(observed_ldl_roh = prim_L[1], observed_lpa_roh = prim_P[1], observed_statin_hr = pc$ldl$statin[1])]
      sim_L <- sim_compare(grid, "primary", "LDL", prim_L); sim_P <- sim_compare(grid, "primary", "LPA", prim_P)
      # The readout says calibrated to the observed values. The rule stays any-cell, the reading
      # that makes ROBUST harder; these are the cells that also match the observed statin main
      # effect, and if there are none the grid has not tested the observed scenario.
      mcell <- grid[estimator == "primary" & score == "LDL" & !is.na(mean_statin_main_hr) & abs(log(mean_statin_main_hr / pc$ldl$statin[1])) < log(1.1), .(statin_hr, u_strength)]
      sim_L_match <- sim_compare(merge(grid, mcell, by = c("statin_hr", "u_strength")), "primary", "LDL", prim_L)
      if (!nrow(mcell) && isFALSE(sim_L$reproduced)) sim_L$reproduced <- NA
      # Like with like: the IPTW estimate is compared with the same estimator on simulated cohorts
      I <- IP[score == "LDL"]
      ires <- lapply(c("raw", "corrected"), function(fm) { r <- I[form == fm]; if (!nrow(r) || is.na(r$estimate)) list(reproduced = NA, incomplete = TRUE) else sim_compare(grid, paste0("iptw_", fm), "LDL", c(r$estimate, r$lo95, r$hi95)) })
      pool_rep <- function(z) list(reproduced = if (!length(z)) NA else if (any(vapply(z, function(r) isTRUE(r$reproduced), TRUE))) TRUE else if (any(vapply(z, function(r) is.na(r$reproduced), TRUE))) NA else FALSE)
      sim_I_all <- pool_rep(ires)
      # An unbalanced real arm says nothing, so neither does its simulated counterpart
      sim_I <- pool_rep(ires[vapply(c("raw", "corrected"), function(fm) isTRUE(I[form == fm, balanced]), TRUE)])
      grid <- merge(grid, sim_L$rows[, .(statin_hr, u_strength, estimator, score, ldl_in_range = in_range, ldl_reproduced = reproduced)],
                    by = c("statin_hr", "u_strength", "estimator", "score"), all.x = TRUE)
      wcsv(grid, opt$outdir, "sim_grid.csv")
      interp <- sim_interpolate(grid, pc$ldl$statin[1]); wcsv(interp, opt$outdir, "sim_interpolated_to_observed_statin_hr.csv")
      joint <- sim_joint(grid, prim_L[1], prim_P[1], pc$ldl$main[1], pc$lpa$main[1], pc$ldl$statin[1])
      wcsv(joint, opt$outdir, "sim_joint_cells.csv")
      g <- grid[estimator == "primary"]
      rng <- rbindlist(lapply(c("LDL", "LPA"), function(sc) {
        x <- g[score == sc]; ob <- if (sc == "LDL") prim_L else prim_P; cmp <- if (sc == "LDL") sim_L else sim_P
        data.table(score = sc, observed = ob[1], observed_lo95 = ob[2], observed_hi95 = ob[3], cells = nrow(x), cells_failed = sum(is.na(x$mean_roh)),
                   sim_min_mean = suppressWarnings(min(x$mean_roh, na.rm = TRUE)), sim_max_mean = suppressWarnings(max(x$mean_roh, na.rm = TRUE)),
                   sim_lowest_q025 = suppressWarnings(min(x$q025, na.rm = TRUE)), sim_highest_q975 = suppressWarnings(max(x$q975, na.rm = TRUE)),
                   observed_in_range_any_cell = cmp$in_range, reproduced_any_cell = cmp$reproduced, cells_reproducing = cmp$k)
      }))
      # The asymmetry itself, per cell and replicate: Lp(a) ratio over LDL-C ratio. Descriptive.
      x <- g[score == "LPA_minus_LDL"]; ob <- prim_P[1] / prim_L[1]
      rng <- rbind(rng, data.table(score = "LPA_minus_LDL (ratio of the two ratios)", observed = ob, cells = nrow(x), cells_failed = sum(is.na(x$mean_roh)),
        sim_min_mean = suppressWarnings(min(x$mean_roh, na.rm = TRUE)), sim_max_mean = suppressWarnings(max(x$mean_roh, na.rm = TRUE)),
        sim_lowest_q025 = suppressWarnings(min(x$q025, na.rm = TRUE)), sim_highest_q975 = suppressWarnings(max(x$q975, na.rm = TRUE)),
        observed_in_range_any_cell = any(x$q025 <= ob & ob <= x$q975, na.rm = TRUE)), fill = TRUE)
      wcsv(rng, opt$outdir, "sim_range_summary.csv"); print(rng)
    }
  } else message("[B] simulation skipped by option: the verdict cannot be ROBUST")

  v <- decide_verdict(prim_L, prim_P, prim_M, IP, sim_L, sim_I)
  strict <- decide_verdict(prim_L, prim_P, prim_M, IP, sim_L, sim_I_all, require_both = TRUE)
  v <- rbind(v, data.table(item = "LDL_C_ATTENUATION_IF_BOTH_ARMS_MUST_BALANCE", verdict = strict[item == "LDL_C_ATTENUATION", verdict]))
  simline <- function(cmp, lbl, ob) data.table(item = lbl, verdict = sprintf(
    "observed %.3f: in a Monte Carlo range = %s; attenuation reproduced by selection alone = %s (%s of %s cells)",
    ob[1], as.character(cmp$in_range), as.character(cmp$reproduced), as.character(cmp$k), as.character(cmp$n_cells)))
  info <- rbind(simline(sim_L, "SIM_LDL_C", prim_L), simline(sim_P, "SIM_LPA", prim_P))
  if (!is.null(sim_L_match)) info <- rbind(info, simline(sim_L_match, "SIM_LDL_C_CELLS_MATCHING_OBSERVED_STATIN_HR", prim_L))
  if (!is.null(joint)) {
    jn <- sum(joint$both_in_range, na.rm = TRUE); jp <- sum(joint$both_in_range & joint$statin_hr_within_10pct_of_observed, na.rm = TRUE)
    it <- interp[score == "LDL" & !is.na(interpolated_mean_roh)]
    info <- rbind(info, data.table(item = c("SIM_JOINT", "SIM_AT_OBSERVED_STATIN_HR"), verdict = c(
      sprintf("descriptive: both observed ratios fall in the Monte Carlo range of the same cell in %d of %d cells; %d of those also match the observed statin main effect within 10%%", jn, nrow(joint), jp),
      if (nrow(it)) sprintf("descriptive: at U strength matching the observed statin HR %.3f, selection alone gives LDL-C ratio of HRs %.3f to %.3f across true statin HR %s", pc$ldl$statin[1], min(it$interpolated_mean_roh), max(it$interpolated_mean_roh), paste(it$statin_hr, collapse = "/"))
      else sprintf("descriptive: the grid does not reach the observed statin HR %.3f", pc$ldl$statin[1]))))
  }
  v <- rbind(v, info,
    data.table(item = "SCOPE", verdict = "Simulation bounds a mechanism; it does not estimate a causal effect. Baseline statin use is a single prevalent-user measurement. The simulated IPTW is linear, without PCs or the BMI, smoking, diabetes and hypertension proxies, on the full N; the real arm can use the spline propensity on complete cases, so that comparison is approximate near the boundary."))
  wcsv(v, opt$outdir, "verdict.csv")
  writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))
  for (i in seq_len(nrow(v))) cat(sprintf("VERDICT %s: %s\n", v$item[i], v$verdict[i]))
  invisible(list(verdict = v, pc = pc, selection = sel, cox = cm, iptw = ip, IP = IP, neg = neg, grid = grid,
                 sim_L = sim_L, sim_P = sim_P, sim_I = sim_I, sim_L_match = sim_L_match, joint = joint, interp = interp))
}

# ---- selftest: synthetic data only ------------------------------------------
# Writes the raw files in the real on-disk formats, so that prep_cohort is tested too
# (quoted headers, sex-specific statin fields, blank dates, flipped alleles, no-follow-up).
st_make_files <- function(dir, scn) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE); raw_dir <- file.path(dir, "raw"); dir.create(raw_dir, showWarnings = FALSE)
  set.seed(scn$seed); n <- scn$N; iid <- 1000000L + seq_len(n)
  mk_score <- function(J, beta_fun, maf_fun, tag) {
    maf <- maf_fun(J); beta <- beta_fun(J)
    effd <- matrix(rbinom(n * J, 2L, rep(maf, each = n)), n, J)
    S <- as.numeric(effd %*% beta)
    flipped <- rep(c(FALSE, FALSE, TRUE), length.out = J)
    cnt <- effd; cnt[, flipped] <- 2L - cnt[, flipped]
    cnt[sample.int(length(cnt), ceiling(0.002 * length(cnt)))] <- NA
    snp <- sprintf("rs%s%03d", tag, seq_len(J)); counted <- ifelse(flipped, "G", "A")
    dt <- data.table(FID = iid, IID = iid, PAT = 0L, MAT = 0L, SEX = 1L, PHENOTYPE = -9L)
    colnames(cnt) <- paste0(snp, "_", counted); dt <- cbind(dt, as.data.table(cnt))
    fwrite(dt, file.path(dir, sprintf("dose_%s.raw", tag)), sep = " ", na = "NA", quote = FALSE)
    fwrite(data.table(snp, "A", beta), file.path(dir, sprintf("%s_score.txt", tolower(tag))), sep = " ", col.names = FALSE, quote = FALSE)
    as.numeric(scale(S))
  }
  G_L <- mk_score(24L, function(J) rnorm(J, 0.1, 0.03), function(J) runif(J, 0.1, 0.5), "LDL")
  G_P <- mk_score(8L, function(J) c(1.0, rep(0.15, J - 1L)), function(J) c(0.12, runif(J - 1L, 0.1, 0.5)), "LPA")
  G_M <- rnorm(n)
  U <- rnorm(n); Lz <- 0.36 * G_L + sqrt(1 - 0.36^2) * rnorm(n)
  age <- runif(n, 40, 70); sex <- rbinom(n, 1L, 0.46)
  bmi <- 27 + 0.8 * U + rnorm(n, 0, 4); smoke <- rbinom(n, 1L, plogis(-2 + 0.3 * U)); dia <- rbinom(n, 1L, plogis(-3 + 0.5 * U)); hyp <- rbinom(n, 1L, plogis(-1 + 0.3 * U))
  lpa_hz <- 0.09; ldl_hz <- 0.05
  eta <- scn$gL * Lz + scn$gU * U + 0.07 * (age - 55) + 0.4 * sex + 0.5 * dia + 0.4 * hyp + 0.03 * (bmi - 27)
  g0 <- uniroot(function(g) mean(plogis(g + eta)) - 0.14, c(-15, 10))$root
  A <- rbinom(n, 1L, plogis(g0 + eta))
  lp <- ldl_hz * G_L + lpa_hz * G_P + scn$bL * Lz + scn$bU * U + 0.07 * (age - 55) + 0.5 * sex + 0.2 * dia + 0.1 * smoke
  lp_t <- lp + log(scn$hr) * A + log(scn$roh_ldl) * G_L * A + log(scn$roh_lpa) * G_P * A + log(scn$roh_mel) * G_M * A
  h <- exp(uniroot(function(lh) mean(-expm1(-exp(lh + lp_t) * 13)) - 0.09, c(-15, 5))$root)
  Ti <- -log(runif(n)) / (h * exp(lp_t))
  # prevalent CAD: same risk, dated before baseline
  prev <- runif(n) < -expm1(-exp(uniroot(function(lk) mean(-expm1(-exp(lk + lp))) - 0.05, c(-30, 10))$root + lp))
  t0 <- as.IDate("2006-03-13") + sample.int(1400L, n, TRUE)
  tcad <- t0 + pmin(ceiling(Ti * 365.25), 40000L); tcad[prev] <- t0[prev] - sample.int(5000L, sum(prev), TRUE)
  tdeath <- t0 + pmin(ceiling(-log(runif(n)) / 0.006 * 365.25), 40000L)
  # some events far after the admin censor stay in the file, as first-occurrence fields do
  tcad[tcad > as.IDate("2024-01-01")] <- NA; tdeath[tdeath > as.IDate("2022-12-19")] <- NA
  zf <- sample.int(n, 60L); tdeath[zf] <- t0[zf]; tcad[zf] <- NA   # no follow-up time
  w1 <- function(fid, v) {
    con <- file.path(raw_dir, sprintf("fid%s.csv", fid))
    writeLines(sprintf('eid,"%s-0.0","%s-1.0"', fid, fid), con)
    fwrite(data.table(eid = iid, a = v, b = NA), con, append = TRUE, col.names = FALSE, na = "", quote = FALSE)
  }
  t0c <- as.character(t0); t0c[sample.int(n, 25L)] <- NA; w1("53", t0c)
  w1("40000", as.character(tdeath))
  # the event date goes in one of six first-occurrence fields; a second date in a later code tests pmin
  which_f <- sample.int(6L, n, TRUE); fids <- c("131296","131298","131300","131302","131304","131306")
  later <- sample.int(n, floor(n / 10))
  for (k in 1:6) {
    x <- rep(NA_character_, n); x[which_f == k] <- as.character(tcad[which_f == k])
    x[later][which_f[later] != k & k == 6L & !is.na(tcad[later])] <- as.character(tcad[later][which_f[later] != k & k == 6L & !is.na(tcad[later])] + 300L)
    w1(fids[k], x)
  }
  # 6153 women, 6177 men: the published fix depends on reading both
  f6153 <- ifelse(sex == 0L, ifelse(A == 1L, 1L, -7L), NA); f6177 <- ifelse(sex == 1L, ifelse(A == 1L, 1L, -7L), NA)
  nm <- sample.int(n, 40L); f6153[nm] <- NA; f6177[nm] <- NA
  w1("6153", f6153); w1("6177", f6177)
  ldl <- (3.7 + 0.9 * Lz) * ifelse(A == 1L, 0.7, 1) * exp(rnorm(n, 0, 0.05)); ldl[ldl <= 0.3] <- 0.3
  rep_l <- rep(NA_integer_, n); miss <- sample.int(n, floor(0.06 * n)); ldl[miss] <- NA; rep_l[miss[seq_len(floor(length(miss) / 2))]] <- 5L
  w1("30780", ldl); w1("30786", rep_l)
  lpa <- exp(3.2 + 1.2 * G_P + rnorm(n, 0, 0.5)); rep_p <- rep(NA_integer_, n)
  rep_p[lpa < 3.8] <- 4L; rep_p[lpa > 189] <- 5L; lpa[!is.na(rep_p)] <- NA
  w1("30790", lpa); w1("30796", rep_p)
  mel <- 1 + 0.3 * G_M; mel[sample.int(n, floor(0.01 * n))] <- NA; w1("26252", mel)
  pcsm <- matrix(rnorm(n * 10L), n, 10L, dimnames = list(NULL, paste0("PC", 1:10)))
  cv <- data.table(FID = iid, IID = iid, age = round(age, 1), sex = sex, bmi = bmi, smoking_current = smoke, diabetes = dia, hypertension = hyp)
  cv[sample.int(n, floor(0.05 * n)), bmi := NA]
  cv[sample.int(n, 37L), age := NA]    # coxph drops these after scaling; the prep has to drop the same rows
  cv <- cbind(cv, as.data.table(pcsm)); cv[sample.int(n, 23L), PC3 := NA]
  fwrite(cv, file.path(dir, "covariates.txt"), sep = " ", na = "NA", quote = FALSE)
  list(raw_dir = raw_dir,
       scores_arg = paste(file.path(dir, "dose_LPA.raw"), file.path(dir, "dose_LDL.raw"), sep = ","),
       weights_arg = paste(file.path(dir, "lpa_score.txt"), file.path(dir, "ldl_score.txt"), sep = ","),
       covar = file.path(dir, "covariates.txt"))
}

st_opt <- function(files, outdir, published, extra = list()) {
  o <- list(scores = files$scores_arg, weights = files$weights_arg, raw_dir = files$raw_dir, covar = files$covar,
            outdir = outdir, censor = "2022-10-31", published = published, cores = max(1L, min(8L, detectCores() - 1L)),
            n_reps = 40L, n_reps_secondary = 10L, grid_hr = "0.8,1", grid_u = "0,0.6,1.2", skip_sim = FALSE, use_glmnet = FALSE, statin_response_sd = 0.2)
  modifyList(o, extra)
}

# Runs the published script itself on the synthetic files and turns its output into the
# targets. A pass therefore means the copied prep is equivalent, not merely plausible.
st_published_targets <- function(files, dir, script) {
  out <- file.path(dir, "published_incident.csv")
  st <- system2("Rscript", c(shQuote(script), "--scores", shQuote(files$scores_arg), "--weights", shQuote(files$weights_arg),
                             "--raw_dir", shQuote(files$raw_dir), "--covar", shQuote(files$covar), "--out", shQuote(out)),
                stdout = FALSE, stderr = FALSE)
  if (st != 0L) stop("selftest: published script failed on synthetic files, exit ", st)
  r <- fread(out)
  pick <- function(p) { x <- r[grepl(p, term, fixed = TRUE)]; c(x$HR, x$lo95, x$hi95) }
  list(values = r, targets = list(
    ldl_int = pick("LDL-C score x"), lpa_int = pick("Lp(a) score x"),
    ldl_main = pick("LDL-C score ->"), lpa_main = pick("Lp(a) score ->"),
    n = r$n[1], events = r$events[1], prevalent = r$prevalent_excluded[1]))
}

run_selftest <- function(opt) {
  base <- opt$selftest_dir; if (is.null(base) || !nzchar(base)) stop("--selftest_dir is required")
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  me <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  pub <- if (nzchar(opt$published_script)) opt$published_script else file.path(dirname(me), "..", "R", "statin_interaction_incident.R")
  if (!file.exists(pub)) stop("selftest: published script not found at ", pub)
  unlink(file.path(base, c("pc_fail", "pc_nudge", "pc_nudge2", "pc_nudge3", "out_A", "out_B", "out_D")), recursive = TRUE)
  n_fail <- 0L
  check <- function(lbl, cond) {
    ok <- isTRUE(cond); cat(sprintf("%-4s %s\n", if (ok) "PASS" else "FAIL", lbl)); if (!ok) n_fail <<- n_fail + 1L; invisible(ok)
  }
  fill_counts <- function(tg, files) {
    # statin counts are not in the published CSV: take them from the prep, then pin every other target
    p0 <- prep_cohort(list(scores = files$scores_arg, weights = files$weights_arg, raw_dir = files$raw_dir, covar = files$covar, censor = "2022-10-31")); p <- p0$d
    tg$n_statin <- sum(p$statin == 1L); tg$events_statin <- sum(p$event[p$statin == 1L]); tg
  }

  # ---- T1 unit tests that need no cohort ----
  cat("--- T1 helpers\n")
  check("sup() suppresses counts below 5", identical(sup(c(0L, 4L, 5L, 100L)), c(NA_integer_, NA_integer_, 5L, 100L)))
  check("wcsv refuses identifier columns", inherits(try(wcsv(data.table(IID = 1), base, "x.csv"), silent = TRUE), "try-error"))
  check("wcsv refuses individual-level row counts", inherits(try(wcsv(data.table(a = 1:6000), base, "x.csv"), silent = TRUE), "try-error"))
  s <- smd_one(c(1, 2, 3, 4), c(1L, 1L, 0L, 0L), rep(1, 4))
  check("smd_one matches the hand calculation", isTRUE(all.equal(s[1], (1.5 - 3.5) / sqrt(0.5))) && isTRUE(all.equal(s[2], s[1])))
  check("smd_one is NA with no variation", all(is.na(smd_one(c(1, 1, 1, 1), c(1L, 1L, 0L, 0L), rep(1, 4)))))
  mkI <- function(est, lo, hi, bal, mel = c(NA, NA, NA)) rbindlist(lapply(c("raw", "corrected"), function(f)
    data.table(form = f, score = c("LDL", "LPA", "MEL"), estimate = c(est, 1, mel[1]), lo95 = c(lo, .95, mel[2]), hi95 = c(hi, 1.05, mel[3]), balanced = bal, max_abs_smd = if (bal) .05 else .24)))
  no <- list(reproduced = FALSE); yes <- list(reproduced = TRUE); na <- list(reproduced = NA)
  good <- c(.93, .90, .95); lpa <- c(.99, .97, 1.01); mel <- c(1, .97, 1.03)
  V <- function(...) { v <- decide_verdict(...); setNames(v$verdict, v$item) }
  a <- V(good, lpa, mel, mkI(.93, .90, .96, TRUE), no, no)
  check("verdict: all conditions met -> ROBUST", a["LDL_C_ATTENUATION"] == "ROBUST")
  check("verdict: Lp(a) bounded null states the interval", grepl("BOUNDED NULL.*0\\.970-1\\.010.*3\\.0%", a["LPA_NO_MODIFICATION"]))
  a <- V(good, lpa, mel, mkI(.93, .90, .96, FALSE), no, no)
  check("verdict: SMD >= 0.1 -> UNINFORMATIVE and not robust", a["LDL_C_ATTENUATION"] == "NOT SEPARABLE FROM SELECTION" && a["IPTW_STATE"] == "UNINFORMATIVE")
  a <- V(good, lpa, mel, mkI(.97, .93, 1.01, TRUE), no, no)
  check("verdict: IPTW CI includes 1 -> not robust", a["LDL_C_ATTENUATION"] == "NOT SEPARABLE FROM SELECTION" && a["IPTW_STATE"] == "CI_INCLUDES_1")
  a <- V(good, lpa, mel, mkI(.93, .90, .96, TRUE), yes, no)
  check("verdict: simulation reproduces -> not robust", a["LDL_C_ATTENUATION"] == "NOT SEPARABLE FROM SELECTION")
  a <- V(good, lpa, mel, mkI(.93, .90, .96, TRUE), na, na)
  check("verdict: simulation incomplete -> not robust", a["LDL_C_ATTENUATION"] == "NOT SEPARABLE FROM SELECTION" && grepl("incomplete", a["LDL_C_REASONS"]))
  a <- V(good, lpa, c(.95, .92, .98), mkI(.93, .90, .96, TRUE), no, no)
  check("verdict: negative control fires -> CONFOUNDED, no Lp(a) null", a["NEGATIVE_CONTROL"] == "CONFOUNDED" && a["LDL_C_ATTENUATION"] != "ROBUST" && a["LPA_NO_MODIFICATION"] == "NO NULL CLAIM")
  a <- V(good, lpa, mel, mkI(1.05, 1.02, 1.08, TRUE), no, no)
  check("verdict: IPTW ratio above 1 is not an attenuation", a["LDL_C_ATTENUATION"] != "ROBUST" && a["IPTW_STATE"] == "CI_ABOVE_1")
  a <- V(c(.98, .95, 1.01), lpa, mel, mkI(.93, .90, .96, TRUE), no, no)
  check("verdict: primary CI including 1 cannot be ROBUST", a["LDL_C_ATTENUATION"] != "ROBUST")
  mixed <- mkI(.93, .90, .96, TRUE); mixed[form == "corrected" & score == "LDL", `:=`(estimate = .97, lo95 = .93, hi95 = 1.01)]
  a <- V(good, lpa, mel, mixed, no, no)
  check("verdict: one LDL-C form with CI including 1 is not enough", a["LDL_C_ATTENUATION"] != "ROBUST" && a["IPTW_STATE"] == "CI_INCLUDES_1")
  one_arm <- function(bal_raw, bal_cor, cor = c(.93, .90, .96)) { x <- mkI(.93, .90, .96, TRUE); x[form == "raw", `:=`(balanced = bal_raw, max_abs_smd = if (bal_raw) .05 else .54)]; x[form == "corrected", `:=`(balanced = bal_cor, max_abs_smd = if (bal_cor) .05 else .3)]; x[form == "corrected" & score == "LDL", `:=`(estimate = cor[1], lo95 = cor[2], hi95 = cor[3])]; x }
  a <- V(good, lpa, mel, one_arm(FALSE, TRUE), no, no)
  check("verdict: unbalanced raw arm + balanced corrected arm below 1 -> ROBUST", a["LDL_C_ATTENUATION"] == "ROBUST" && grepl("raw: UNINFORMATIVE", a["IPTW_ARMS"]))
  a <- V(good, lpa, mel, one_arm(FALSE, TRUE), no, no, require_both = TRUE)
  check("verdict: the same with both arms required -> UNINFORMATIVE", a["LDL_C_ATTENUATION"] != "ROBUST" && a["IPTW_STATE"] == "UNINFORMATIVE")
  a <- V(good, lpa, mel, one_arm(FALSE, TRUE, c(.97, .93, 1.01)), no, no)
  check("verdict: an unbalanced arm cannot rescue a balanced arm that includes 1", a["LDL_C_ATTENUATION"] != "ROBUST" && a["IPTW_STATE"] == "CI_INCLUDES_1")
  a <- V(good, lpa, mel, one_arm(TRUE, FALSE, c(.80, .75, .85)), no, no)
  check("verdict: a balanced raw arm alone cannot license ROBUST", a["IPTW_STATE"] == "UNINFORMATIVE" && a["LDL_C_ATTENUATION"] != "ROBUST")
  x <- mkI(.93, .90, .96, TRUE); x[form == "raw" & score == "LDL", `:=`(estimate = .97, lo95 = .93, hi95 = 1.01)]
  a <- V(good, lpa, mel, x, no, no)
  check("verdict: a balanced raw arm including 1 blocks ROBUST", a["LDL_C_ATTENUATION"] != "ROBUST")
  x <- one_arm(FALSE, TRUE); x[form == "raw" & score == "LDL", `:=`(estimate = NA_real_, lo95 = NA_real_, hi95 = NA_real_)]
  a <- V(good, lpa, mel, x, no, no)
  check("verdict: a failed fit in either arm blocks the call", a["LDL_C_ATTENUATION"] != "ROBUST" && a["IPTW_STATE"] == "UNAVAILABLE")
  bt <- balance_table(data.table(statin = rep(c(0L, 1L), each = 50), a = rnorm(100) + rep(c(0, 1), each = 50), b = rnorm(100) + rep(c(0, 1), each = 50)),
                      c("a", "b"), "a", rep(1, 100), rep(1, 100), not_counted = "b")
  check("balance table separates structural exceptions from counted rows", identical(bt$counts_for_rule, c(TRUE, FALSE)) && !bt$pass[1])
  tt <- PUBLISHED; d_edge <- abs(1.0151 - tt$lpa_int_3dp[3]) <= 0.0005 + 1e-9
  check("rounding-edge Lp(a) upper limit 1.0151 is not a false stop", d_edge && abs(1.0151 - tt$lpa_int[3]) > 0.005)
  a <- V(good, lpa, c(NA, NA, NA), mkI(.93, .90, .96, TRUE), no, no)
  check("verdict: negative control missing -> not robust", a["NEGATIVE_CONTROL"] == "UNAVAILABLE" && a["LDL_C_ATTENUATION"] != "ROBUST")
  a <- V(good, c(.90, .85, .95), mel, mkI(.93, .90, .96, TRUE), no, no)
  check("verdict: Lp(a) CI excludes 1 -> no null claim", a["LPA_NO_MODIFICATION"] == "NO NULL CLAIM")
  a <- V(c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA), mkI(NA, NA, NA, FALSE), na, na)
  check("verdict: all NA stays NA and not robust", a["LDL_C_ATTENUATION"] == "NOT SEPARABLE FROM SELECTION" && a["LPA_NO_MODIFICATION"] == "NO NULL CLAIM")
  gg <- data.table(estimator = "primary", score = "LDL", statin_hr = c(.7, 1), u_strength = 0, mean_roh = c(.99, .985), q025 = c(.97, .965), q975 = c(1.01, 1.0))
  check("sim_compare: not reached", isFALSE(sim_compare(gg, "primary", "LDL", c(.93, .90, .95))$reproduced))
  check("sim_compare: reached", isTRUE(sim_compare(gg, "primary", "LDL", c(.968, .95, .96))$reproduced))
  gg2 <- copy(gg); gg2[1, `:=`(mean_roh = NA_real_, q025 = NA_real_, q975 = NA_real_)]
  check("sim_compare: failed cell gives NA rather than FALSE", is.na(sim_compare(gg2, "primary", "LDL", c(.93, .90, .95))$reproduced))
  check("sim_compare: a reaching cell still counts despite a failed cell", isTRUE(sim_compare(gg2, "primary", "LDL", c(.968, .95, .96))$reproduced))

  # ---- T2 degenerate designs must give NA ----
  cat("--- T2 degenerate designs\n")
  set.seed(42)
  dg <- data.table(fu = rexp(2000, .1) + .01, event = rbinom(2000, 1, .2), score_X = rnorm(2000), statin = 0L, age = rnorm(2000, 55, 8), sex = rbinom(2000, 1, .5))
  r <- fit_cox("Surv(fu, event) ~ score_X * statin + age + sex", dg, "m", "score_X:statin", "score_X")
  check("no statin users -> NA", is.na(r$estimate) && grepl("^NA:", r$note))
  dg2 <- copy(dg); dg2[, statin := rbinom(.N, 1, .3)]; dg2[, event := 0L]; dg2[1:3, event := 1L]
  r <- fit_cox("Surv(fu, event) ~ score_X * statin + age + sex", dg2, "m", "score_X:statin", "score_X")
  check("too few events -> NA, event count suppressed", is.na(r$estimate) && is.na(r$events))
  dg3 <- copy(dg); dg3[, statin := rbinom(.N, 1, .3)]; dg3[, score_X := 1]
  r <- fit_cox("Surv(fu, event) ~ score_X * statin + age + sex", dg3, "m", "score_X:statin", "score_X")
  check("constant score -> NA", is.na(r$estimate))
  dg4 <- copy(dg); dg4[, statin := rbinom(.N, 1, .3)]; dg4[, z := score_X * 2]
  r <- fit_cox("Surv(fu, event) ~ score_X * statin + z:statin + age + sex", dg4, "m", "z:statin", "score_X")
  check("collinear term -> NA, not a number", is.na(r$estimate))
  r <- fit_logit("score_X + age", dg, "m", "score_X")
  check("logistic with one class -> NA", is.na(r$estimate))
  check("IPTW with one statin class returns no estimate",
        { z <- copy(dg); for (v in c("bmi","smoking_current","diabetes","hypertension","LDL","LDL_corr",paste0("PC",1:10),"score_LDL","score_LPA","score_MEL")) z[, (v) := rnorm(.N)]
          o <- iptw_block(z, paste0("PC",1:10), paste(paste0("PC",1:10), collapse=" + "), FALSE); all(is.na(o$est$estimate)) && !any(o$chosen$balanced) })

  # ---- T3 positive control and equivalence with the published script ----
  cat("--- T3 positive control\n")
  scA <- list(N = 150000L, seed = 42L, gL = 1.6, gU = 1.3, bL = 0.35, bU = 1.1, hr = 0.8, roh_ldl = 1, roh_lpa = 1, roh_mel = 1)
  dA <- file.path(base, "scenA"); fA <- st_make_files(dA, scA)
  pA <- st_published_targets(fA, dA, pub)
  tgA <- fill_counts(pA$targets, fA)
  # Real exit code through the command line, with the real (wrong for this data) published targets
  outbad <- file.path(base, "pc_fail"); unlink(outbad, recursive = TRUE)
  st <- system2("Rscript", c(shQuote(me), "--scores", shQuote(fA$scores_arg), "--weights", shQuote(fA$weights_arg), "--raw_dir", shQuote(fA$raw_dir),
                             "--covar", shQuote(fA$covar), "--outdir", shQuote(outbad), "--cores", "2", "--no_glmnet"), stdout = FALSE, stderr = FALSE)
  check("positive-control failure exits with status 2", st == 2L)
  check("positive-control failure writes only its own table", identical(list.files(outbad), "positive_control.csv"))
  check("positive-control table records the failure", any(fread(file.path(outbad, "positive_control.csv"))$pass == FALSE))
  # A 0.01 miss on the interaction alone must trip the stop
  nudged <- tgA; nudged$ldl_int <- tgA$ldl_int + c(0.0101, 0, 0)
  e <- tryCatch(run_analysis(st_opt(fA, file.path(base, "pc_nudge"), nudged, list(skip_sim = TRUE))), positive_control_failure = function(c) "stopped")
  check("a 0.01 shift in the LDL-C target stops the run", identical(e, "stopped"))
  nudged <- tgA; nudged$events <- tgA$events + 1L
  e <- tryCatch(run_analysis(st_opt(fA, file.path(base, "pc_nudge2"), nudged, list(skip_sim = TRUE))), positive_control_failure = function(c) "stopped")
  check("a one-event difference stops the run", identical(e, "stopped"))
  nudged <- tgA; nudged$prevalent <- tgA$prevalent + 1L
  e <- tryCatch(run_analysis(st_opt(fA, file.path(base, "pc_nudge3"), nudged, list(skip_sim = TRUE))), positive_control_failure = function(c) "stopped")
  check("a one-person difference in the prevalent exclusion stops the run", identical(e, "stopped"))

  # ---- T4 planted selection, no true interaction ----
  cat("--- T4 scenario A: planted selection, no true interaction\n")
  ra <- run_analysis(st_opt(fA, file.path(base, "out_A"), tgA))
  vA <- setNames(ra$verdict$verdict, ra$verdict$item)
  check("copied prep matches the published script exactly (LDL-C)", isTRUE(all.equal(ra$pc$ldl$int, tgA$ldl_int, tolerance = 1e-10)))
  check("copied prep matches the published script exactly (Lp(a))", isTRUE(all.equal(ra$pc$lpa$int, tgA$lpa_int, tolerance = 1e-10)))
  cat(sprintf("     observed LDL-C ratio of HRs %.3f (%.3f-%.3f)\n", ra$pc$ldl$int[1], ra$pc$ldl$int[2], ra$pc$ldl$int[3]))
  check("the planted mechanism yields a CI-excluding-1 attenuation (test is not vacuous)", ra$pc$ldl$int[3] < 1)
  check("selection alone reproduces the attenuation", isTRUE(ra$sim_L$reproduced))
  check("verdict is NOT SEPARABLE FROM SELECTION", vA["LDL_C_ATTENUATION"] == "NOT SEPARABLE FROM SELECTION")
  # The engine against the truth: the planted data used statin HR 0.8 and unmeasured strength
  # about 1.2, so the cell nearest that truth must place the observed ratios in its own range.
  # A simulation that reached 0.93 only in some unrelated corner would pass the line above.
  jt <- ra$joint[statin_hr == 0.8 & u_strength == 1.2]
  cat(sprintf("     truth cell: LDL-C %.3f (%.3f-%.3f), Lp(a) %.3f (%.3f-%.3f); statin main HR simulated %.2f, observed %.2f\n",
              jt$ldl_mean_roh, jt$ldl_q025, jt$ldl_q975, jt$lpa_mean_roh, jt$lpa_q025, jt$lpa_q975, jt$sim_statin_main_hr, ra$pc$ldl$statin[1]))
  check("the cell nearest the planted truth covers observed LDL-C and Lp(a) together", isTRUE(jt$both_in_range))
  check("the cell nearest the planted truth matches the observed statin main effect within 10%", isTRUE(jt$statin_hr_within_10pct_of_observed))
  check("a cell with no unmeasured selection falls short of the same attenuation", ra$joint[statin_hr == 0.8 & u_strength == 0, ldl_mean_roh] > ra$pc$ldl$int[1])
  check("selection on statin use is visible in the logistic table", ra$selection$logit[model == "LDL-C score alone", estimate] > 1.1)
  check("melanoma score does not predict statin use", { z <- ra$selection$logit[model == "melanoma score alone (negative control)"]; z$lo95 < 1 && z$hi95 > 1 })
  check("all outputs are summary level", { fs <- list.files(file.path(base, "out_A"), pattern = "csv$", full.names = TRUE); all(vapply(fs, function(f) { x <- fread(f); nrow(x) < 5000L && !any(tolower(names(x)) %in% c("iid", "eid")) }, TRUE)) })

  # ---- T5 planted true interaction ----
  cat("--- T5 scenario B: planted true interaction, weak measured selection\n")
  scB <- list(N = 150000L, seed = 42L, gL = 0.3, gU = 0, bL = 0.1, bU = 0, hr = 0.8, roh_ldl = 0.78, roh_lpa = 1, roh_mel = 1)
  dB <- file.path(base, "scenB"); fB <- st_make_files(dB, scB)
  tgB <- fill_counts(st_published_targets(fB, dB, pub)$targets, fB)
  rb <- run_analysis(st_opt(fB, file.path(base, "out_B"), tgB))
  vB <- setNames(rb$verdict$verdict, rb$verdict$item)
  cat(sprintf("     observed LDL-C ratio of HRs %.3f (%.3f-%.3f)\n", rb$pc$ldl$int[1], rb$pc$ldl$int[2], rb$pc$ldl$int[3]))
  print(rb$IP)
  check("/0.7 arm is balanced", isTRUE(rb$iptw$chosen[form == "corrected", balanced]))
  check("raw arm is not balanced: the drug moves raw LDL-C, as it will in real data", isFALSE(rb$iptw$chosen[form == "raw", balanced]))
  check("IPTW CI excludes 1 in the balanced arm", rb$IP[score == "LDL" & form == "corrected", hi95] < 1)
  check("selection alone does NOT reproduce the attenuation", isFALSE(rb$sim_L$reproduced))
  check("verdict is ROBUST", vB["LDL_C_ATTENUATION"] == "ROBUST")
  check("the stricter both-arms reading is reported and differs", vB["LDL_C_ATTENUATION_IF_BOTH_ARMS_MUST_BALANCE"] == "NOT SEPARABLE FROM SELECTION")
  check("Lp(a) stays a bounded null with its interval", grepl("^BOUNDED NULL", vB["LPA_NO_MODIFICATION"]))
  check("A and B are separated", vA["LDL_C_ATTENUATION"] != vB["LDL_C_ATTENUATION"])

  # ---- T6 planted negative control and Lp(a) interaction ----
  cat("--- T6 scenario D: the same as B, plus melanoma x statin and Lp(a) x statin\n")
  scD <- modifyList(scB, list(roh_mel = 0.85, roh_lpa = 0.85))
  dD <- file.path(base, "scenD"); fD <- st_make_files(dD, scD)
  tgD <- fill_counts(st_published_targets(fD, dD, pub)$targets, fD)
  rd <- run_analysis(st_opt(fD, file.path(base, "out_D"), tgD, list(skip_sim = TRUE)))
  vD <- setNames(rd$verdict$verdict, rd$verdict$item)
  check("melanoma interaction detected", rd$neg$hi95 < 1)
  check("negative control fires -> CONFOUNDED", vD["NEGATIVE_CONTROL"] == "CONFOUNDED")
  check("no robust call and no Lp(a) null claim", vD["LDL_C_ATTENUATION"] != "ROBUST" && vD["LPA_NO_MODIFICATION"] == "NO NULL CLAIM")
  check("skipped simulation cannot yield ROBUST", vD["LDL_C_ATTENUATION"] != "ROBUST")

  # ---- T7 the simulation engine itself ----
  cat("--- T7 simulation calibration and null behaviour\n")
  # No selection whatsoever: statin independent of everything. The engine has to return a null
  # interaction and hand back the statin HR it was given.
  set.seed(42); N <- 150000L
  pool <- data.frame(G_L = rnorm(N), G_P = rnorm(N), age = runif(N, -15, 15), sex = rbinom(N, 1, .5))
  cfg0 <- list(N = N, events = 13000L, p_statin = .14, prev_frac = .05, mu = 3.7, sdl = 0.9, r_L = .36, r_P = 0.02, sigma_e = sqrt(1 - .36^2),
               t_logit = c(qlogis(.14), 0, 0, 0, 0), t_cox = c(.05, .09, .1, .07, .5), pool = pool, fu_pool = runif(50000, 3, 16))
  cl <- calibrate_cell(cfg0, s = 0, log_hr = log(0.75), seed = 42L)
  check("calibration converges", cl$converged)
  check("calibration matches the requested moments", cl$max_residual < 0.005)
  check("no-selection statin selection coefficient is near zero", abs(cl$par$gL) < 0.05)
  reps <- rbindlist(lapply(mclapply(1:30, function(r) as.list(sim_rep(cl, cfg0, 42L + 5000L + r, secondary = FALSE)), mc.cores = max(1L, min(8L, detectCores() - 1L)), mc.set.seed = FALSE), as.data.table))
  check("no selection -> interaction null (LDL-C)", abs(mean(reps$primary_L)) < 2.5 * sd(reps$primary_L) / sqrt(30) + 0.005)
  check("no selection -> interaction null (Lp(a))", abs(mean(reps$primary_P)) < 2.5 * sd(reps$primary_P) / sqrt(30) + 0.005)
  check("no selection -> statin HR recovered", abs(mean(reps$primary_statin_L) - log(.75)) < 0.02)
  # Selection on LDL-C and on U must attenuate the LDL-C slope more than the Lp(a) slope
  # 48 replicates at .5 / 1.5: a 48-replicate grid gave ratios of 0.957 and 0.978, SE 0.003 each.
  # The bias is small at weaker settings (0.980 at .25 / 1), which 12 replicates cannot resolve.
  cfg1 <- modifyList(cfg0, list(t_logit = c(qlogis(.14), .5, .02, .07, .4)))
  cl1 <- calibrate_cell(cfg1, s = 1.5, log_hr = log(0.75), seed = 42L + 1L)
  rp1 <- rbindlist(lapply(mclapply(1:48, function(r) as.list(sim_rep(cl1, cfg1, 42L + 6000L + r, secondary = FALSE)), mc.cores = max(1L, min(8L, detectCores() - 1L)), mc.set.seed = FALSE), as.data.table))
  se_L <- sd(rp1$primary_L) / sqrt(48); se_D <- sd(rp1$primary_L - rp1$primary_P) / sqrt(48)
  cat(sprintf("     selection on: LDL-C ratio %.3f, Lp(a) ratio %.3f\n", exp(mean(rp1$primary_L)), exp(mean(rp1$primary_P))))
  check("selection alone attenuates the LDL-C slope", mean(rp1$primary_L) + 2.5 * se_L < 0)
  check("selection attenuates LDL-C more than Lp(a)", mean(rp1$primary_L - rp1$primary_P) + 2.5 * se_D < 0)
  # Resubmission after the 12 h limit: finished cells are reused, a short cell is rerun
  op_c <- list(outdir = file.path(base, "resume"), grid_hr = "0.8", grid_u = "0,1", n_reps = 6L, n_reps_secondary = 0L, cores = max(1L, min(8L, detectCores() - 1L)))
  unlink(op_c$outdir, recursive = TRUE)
  small <- modifyList(cfg0, list(t_logit = c(qlogis(.14), .25, .02, .07, .4), N = 30000L, events = 2600L, fu_pool = cfg1$fu_pool[1:20000]))
  r1 <- run_sim(small, op_c)
  cks <- file.path(op_c$outdir, "sim_cells", "cell_001.rds"); stamp <- file.mtime(file.path(op_c$outdir, "sim_cells", "cell_002.rds"))
  sv <- readRDS(cks); sv$res$reps <- rep(list(NA_real_), 6L); saveRDS(sv, cks)
  Sys.sleep(1.1); r2 <- run_sim(small, op_c)
  check("resume: a finished cell is not recomputed", file.mtime(file.path(op_c$outdir, "sim_cells", "cell_002.rds")) == stamp)
  check("resume: a cell short of replicates is rerun, to the same numbers", isTRUE(all.equal(r1$grid, r2$grid)) && !anyNA(r2$grid[estimator == "primary", mean_roh]))
  # Calibration that cannot be met must fail rather than report a number
  cf_bad <- modifyList(cfg0, list(t_logit = c(qlogis(.14), 50, 0, 0, 0)))
  cb <- guarded(calibrate_cell(cf_bad, s = 0, log_hr = 0, seed = 42L, max_iter = 4L))$value
  check("impossible calibration target is reported as not converged", is.null(cb) || !cb$converged)

  cat("--- T8 outputs\n")
  check("no individual-level files were created by the analysis", { fs <- list.files(file.path(base, "out_A"), recursive = TRUE); !any(grepl("\\.raw$|covariates", fs)) })
  check("verdict lines were written", all(c("LDL_C_ATTENUATION", "SIM_LDL_C", "SIM_JOINT", "SIM_AT_OBSERVED_STATIN_HR") %in% fread(file.path(base, "out_A", "verdict.csv"))$item))
  check("the balance table reports every covariate with its rule status", { b <- fread(file.path(base, "out_A", "iptw_balance.csv")); all(c("score_LDL", "score_MEL", "LDL", "LDL_corr", "LPA_log") %in% b$covariate) && "counts_for_rule" %in% names(b) })
  check("the balance table carries the row the verdict used", file.exists(file.path(base, "out_A", "iptw_verdict_inputs.csv")) && file.exists(file.path(base, "out_A", "iptw_chosen_propensity.csv")))

  if (n_fail) { cat(sprintf("SELFTEST FAILED: %d check(s)\n", n_fail)); quit(save = "no", status = 1L) }
  cat("SELFTEST PASSED\n")
}

# ---- command line -----------------------------------------------------------
main <- function() {
  o <- parse_args(OptionParser(option_list = list(
    make_option("--scores",  type = "character", help = "dose_LPA.raw,dose_LDL.raw"),
    make_option("--weights", type = "character", help = "lpa_score.txt,ldl_score.txt"),
    make_option("--raw_dir", type = "character"),
    make_option("--covar",   type = "character"),
    make_option("--outdir",  type = "character"),
    make_option("--censor",  type = "character", default = "2022-10-31"),
    make_option("--cores",   type = "integer",   default = 16L),
    make_option("--n_reps",  type = "integer",   default = 200L, help = "replicates per grid cell, at least 200"),
    make_option("--n_reps_secondary", type = "integer", default = 50L, help = "of those, how many also run LDL-adjusted and IPTW estimators"),
    make_option("--grid_hr", type = "character", default = "0.6,0.7,0.8,0.9,1"),
    make_option("--grid_u",  type = "character", default = "0,0.25,0.5,0.75,1,1.25"),
    make_option("--statin_response_sd", type = "double", default = 0.2, help = "SD of log statin LDL-C response in simulated users; an assumption"),
    make_option("--skip_sim", action = "store_true", default = FALSE),
    make_option("--no_glmnet", action = "store_true", default = FALSE),
    make_option("--glmnet_lib", type = "character", default = "", help = "library searched LAST, for glmnet only"),
    make_option("--selftest", action = "store_true", default = FALSE),
    make_option("--selftest_dir", type = "character", default = ""),
    make_option("--published_script", type = "character", default = "")
  )))
  # Appended, not R_LIBS_USER: that would put this library ahead of the conda env and could
  # swap the survival or data.table underneath a positive control tied to their versions.
  if (nzchar(o$glmnet_lib)) .libPaths(c(.libPaths(), o$glmnet_lib))
  if (o$selftest) { run_selftest(o); return(invisible()) }
  for (a in c("scores", "weights", "raw_dir", "covar", "outdir")) if (is.null(o[[a]])) stop("--", a, " is required")
  if (o$n_reps < 200L && !o$skip_sim) stop("--n_reps must be at least 200 (preregistered analysis plan)")
  o$published <- PUBLISHED; o$use_glmnet <- !o$no_glmnet
  r <- tryCatch(run_analysis(o), positive_control_failure = function(c) { message(conditionMessage(c)); quit(save = "no", status = 2L) })
  message("[B] wrote ", o$outdir)
}
if (sys.nframe() == 0L) main()
