# ============================================================================
# lib_resilience.R - shared, pure functions for the metabolome-resilience
# pipeline. Every function here is unit-tested in tests/testthat.
# No function reads global state except .checks (assertion log for manifests).
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(jsonlite)
  library(yaml)
  library(digest)
})

# ---------------------------------------------------------------------------
# Logging, assertions, CLI
# ---------------------------------------------------------------------------
.checks <- new.env()
.checks$log <- list()

log_msg <- function(fmt, ...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(fmt, ...)))
}

# A check is a hard invariant. Failure aborts the stage; every check (pass or
# fail) is recorded in the stage manifest so a reader can audit what was tested.
check <- function(cond, fmt, ...) {
  msg <- sprintf(fmt, ...)
  ok  <- isTRUE(cond)
  .checks$log[[length(.checks$log) + 1]] <- list(check = msg, pass = ok)
  if (!ok) stop("CHECK FAILED: ", msg, call. = FALSE)
  log_msg("check ok: %s", msg)
  invisible(TRUE)
}

reset_checks <- function() .checks$log <- list()

# Warnings are never silently dropped: every distinct message is counted, logged and
# written to the stage manifest.
.checks$warnings <- list()
capture_warnings <- function(expr, where = "") {
  withCallingHandlers(expr, warning = function(w) {
    key <- paste0(where, ": ", conditionMessage(w))
    .checks$warnings[[key]] <- (if (is.null(.checks$warnings[[key]])) 0L else .checks$warnings[[key]]) + 1L
    invokeRestart("muffleWarning")
  })
}
report_warnings <- function() {
  if (!length(.checks$warnings)) { log_msg("warnings: none"); return(invisible(NULL)) }
  for (k in names(.checks$warnings)) log_msg("WARNING x%d  %s", .checks$warnings[[k]], k)
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE), required = character(0)) {
  out <- list(); i <- 1L
  while (i <= length(args)) {
    k <- args[i]
    if (!startsWith(k, "--")) stop("Unexpected argument: ", k)
    k <- sub("^--", "", k)
    if (i + 1L <= length(args) && !startsWith(args[i + 1L], "--")) {
      out[[k]] <- args[i + 1L]; i <- i + 2L
    } else {
      out[[k]] <- TRUE; i <- i + 1L
    }
  }
  missing <- setdiff(required, names(out))
  if (length(missing)) stop("Missing required arguments: ", paste0("--", missing, collapse = " "))
  out
}

load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$.path <- normalizePath(path)
  cfg$.md5  <- digest::digest(file = path, algo = "md5")
  cfg
}

# ---------------------------------------------------------------------------
# Manifest: provenance record written by every stage
# ---------------------------------------------------------------------------
file_fingerprint <- function(path) {
  if (!file.exists(path)) return(list(path = path, exists = FALSE))
  sz <- file.info(path)$size
  list(path = normalizePath(path), exists = TRUE, size_bytes = sz,
       mtime = format(file.info(path)$mtime, "%Y-%m-%dT%H:%M:%S"),
       md5 = if (sz < 2e9) digest::digest(file = path, algo = "md5") else "skipped_gt_2GB")
}

git_commit <- function() {
  x <- tryCatch(system("git rev-parse HEAD 2>/dev/null", intern = TRUE), error = function(e) character(0))
  dirty <- tryCatch(length(system("git status --porcelain 2>/dev/null", intern = TRUE)) > 0, error = function(e) NA)
  list(commit = if (length(x)) x[1] else NA_character_, dirty = dirty)
}

# md5 of the stage script and of this library, so a manifest is tied to the code that made it
code_fingerprint <- function() {
  sf <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
  lib <- file.path(dirname(sf), "lib_resilience.R")
  list(stage_script = basename(sf), stage_md5 = if (!is.na(sf) && file.exists(sf)) digest::digest(file = sf, algo = "md5") else NA,
       lib_md5 = if (!is.na(sf) && file.exists(lib)) digest::digest(file = lib, algo = "md5") else NA)
}

write_manifest <- function(stage, out_dir, cfg, inputs = character(0), outputs = character(0), extra = list()) {
  pk <- c("data.table", "survival", "glmnet", "jsonlite", "yaml", "digest")
  m <- list(
    stage = stage,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    git = git_commit(),
    config = list(path = cfg$.path, md5 = cfg$.md5, seed = cfg$seed),
    code = code_fingerprint(),
    r_version = R.version.string,
    packages = setNames(lapply(pk, function(p) as.character(packageVersion(p))), pk),
    inputs  = lapply(inputs, file_fingerprint),
    outputs = lapply(outputs, file_fingerprint),
    checks  = .checks$log,
    warnings = lapply(names(.checks$warnings), function(k) list(message = k, count = .checks$warnings[[k]])),
    extra   = extra
  )
  path <- file.path(out_dir, sprintf("%s.manifest.json", stage))
  jsonlite::write_json(m, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  report_warnings()
  log_msg("manifest written: %s (%d checks, all passed; %d distinct warnings)", path, length(.checks$log), length(.checks$warnings))
  invisible(path)
}

# ---------------------------------------------------------------------------
# UKB split-CSV phenotype store readers
# ---------------------------------------------------------------------------
field_path <- function(pheno_dir, field) file.path(pheno_dir, sprintf("fid%s.csv", field))

# Read a whole field file. eid is coerced to integer, values to character
# unless numeric=TRUE. Missing = "" or "NA".
read_field <- function(pheno_dir, field, numeric = FALSE) {
  f <- field_path(pheno_dir, field)
  check(file.exists(f), "field file exists: fid%s.csv", field)
  dt <- fread(f, na.strings = c("", "NA"), colClasses = "character")
  check("eid" %in% names(dt), "fid%s.csv has an eid column", field)
  dt[, eid := as.integer(eid)]
  if (numeric) for (cn in setdiff(names(dt), "eid")) set(dt, j = cn, value = suppressWarnings(as.numeric(dt[[cn]])))
  dt
}

# Column names for a given instance, e.g. "41270-0.0" ... "41270-0.258"
field_cols <- function(dt, field, instance = 0) {
  pat <- sprintf("^%s-%d\\.\\d+$", field, instance)
  cols <- grep(pat, names(dt), value = TRUE)
  cols[order(as.integer(sub(".*\\.", "", cols)))]
}

# Read one instance of one field -> data.table(eid, <name>) using array 0 only.
read_scalar <- function(pheno_dir, field, name, instance = 0, numeric = TRUE) {
  dt <- read_field(pheno_dir, field, numeric = numeric)
  cols <- field_cols(dt, field, instance)
  check(length(cols) >= 1, "fid%s has instance %d columns", field, instance)
  out <- dt[, .(eid, value = get(cols[1]))]
  setnames(out, "value", name)
  out
}

# Sample-level QC flag (coded 1 or missing). If the field has no column for the
# requested instance (nobody flagged at that visit), everyone gets 0.
read_flag <- function(pheno_dir, field, name, instance = 0) {
  dt <- read_field(pheno_dir, field, numeric = TRUE)
  cols <- field_cols(dt, field, instance)
  v <- if (length(cols)) as.integer(!is.na(dt[[cols[1]]]) & dt[[cols[1]]] == 1) else rep(0L, nrow(dt))
  out <- data.table(eid = dt$eid, value = v); setnames(out, "value", name); out
}

# Any array element of instance `instance` equal to one of `values` -> 1/0
read_any_flag <- function(pheno_dir, field, name, values, instance = 0) {
  dt <- read_field(pheno_dir, field)
  cols <- field_cols(dt, field, instance)
  check(length(cols) >= 1, "fid%s has instance %d columns", field, instance)
  vals <- as.character(values)
  flag <- rowSums(sapply(cols, function(cn) dt[[cn]] %in% vals)) > 0
  out <- data.table(eid = dt$eid, value = as.integer(flag))
  setnames(out, "value", name)
  out
}

# Convert paired wide code/date arrays (41270/41280 or 41272/41282) to long.
wide_pairs_to_long <- function(codes_dt, dates_dt, code_field, date_field, source, instance = 0) {
  cc <- field_cols(codes_dt, code_field, instance)
  dc <- field_cols(dates_dt, date_field, instance)
  check(length(cc) == length(dc), "%s/%s arrays have equal length (%d vs %d)", code_field, date_field, length(cc), length(dc))
  m <- merge(codes_dt[, c("eid", cc), with = FALSE], dates_dt[, c("eid", dc), with = FALSE], by = "eid")
  long_c <- melt(m[, c("eid", cc), with = FALSE], id.vars = "eid", variable.name = "k", value.name = "code")
  long_d <- melt(m[, c("eid", dc), with = FALSE], id.vars = "eid", variable.name = "k", value.name = "date")
  long_c[, k := as.integer(sub(".*\\.", "", k))]
  long_d[, k := as.integer(sub(".*\\.", "", k))]
  long <- merge(long_c, long_d, by = c("eid", "k"))
  long <- long[!is.na(code) & code != ""]
  long[, date := as.IDate(date)]
  long[, source := source]
  long[, .(eid, source, code, date)]
}

# ---------------------------------------------------------------------------
# Outcome construction (pure; unit-tested with hand-built tables)
# ---------------------------------------------------------------------------
code_matches <- function(code, prefixes) {
  code <- gsub("\\.", "", as.character(code))
  if (length(prefixes) == 0) return(rep(FALSE, length(code)))
  Reduce(`|`, lapply(prefixes, function(p) startsWith(code, p)))
}

# events: data.table(eid, source in {"icd10","opcs4","death"}, code, date)
# base:   data.table(eid, baseline_date, death_date, lost_date, self_report_prevalent)
# Returns one row per eid in base with: prevalent (0/1), event (0/1), time_years,
# end_date, first_event_date, reason (event/death/lost/admin).
build_incident_outcome <- function(base, events, def, admin_censor_date) {
  admin <- as.IDate(admin_censor_date)
  ev <- copy(events)
  ev[, hit := (source == "icd10" & code_matches(code, def$icd10_prefixes)) |
              (source == "opcs4" & code_matches(code, def$opcs4_prefixes)) |
              (source == "death" & code_matches(code, def$death_icd10_prefixes))]
  ev <- ev[hit == TRUE & !is.na(date)]
  ev <- merge(ev, base[, .(eid, baseline_date)], by = "eid")
  prev_eids  <- unique(ev[date <= baseline_date, eid])
  first_post <- ev[date > baseline_date, .(first_event_date = min(date)), by = eid]

  out <- copy(base)
  out[, prevalent := as.integer(eid %in% prev_eids | (!is.na(self_report_prevalent) & self_report_prevalent == 1))]
  out <- merge(out, first_post, by = "eid", all.x = TRUE)
  for (cn in c("first_event_date", "death_date", "lost_date")) set(out, j = cn, value = as.IDate(out[[cn]]))
  out[, end_date := as.IDate(pmin(as.integer(first_event_date), as.integer(death_date), as.integer(lost_date), as.integer(admin), na.rm = TRUE))]
  out[, event := as.integer(!is.na(first_event_date) & first_event_date <= end_date)]
  out[, reason := fifelse(event == 1, "event",
                  fifelse(!is.na(death_date) & death_date == end_date, "death",
                  fifelse(!is.na(lost_date) & lost_date == end_date, "lost", "admin")))]
  out[, time_years := as.numeric(end_date - baseline_date) / 365.25]
  out
}

# Landmark: drop anyone whose follow-up ended before baseline + landmark_years,
# restart the clock at the landmark.
apply_landmark <- function(dt, landmark_years) {
  x <- dt[time_years > landmark_years]
  x[, time_analysis := time_years - landmark_years]
  x
}

# ---------------------------------------------------------------------------
# NMR processing
# ---------------------------------------------------------------------------
# X: numeric matrix (n x p) of raw concentrations. Returns log-transformed,
# winsorised matrix. Zeros -> half the minimum positive value per column.
prep_nmr_matrix <- function(X, log_transform = TRUE, winsor_sd = 5) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (log_transform) {
      pos <- v[!is.na(v) & v > 0]
      if (length(pos) == 0) { X[, j] <- NA_real_; next }
      v[!is.na(v) & v <= 0] <- min(pos) / 2
      v <- log(v)
    }
    if (is.finite(winsor_sd)) {
      mu <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
      if (is.finite(s) && s > 0) v <- pmin(pmax(v, mu - winsor_sd * s), mu + winsor_sd * s)
    }
    X[, j] <- v
  }
  X
}

# Residualise every column of Y on design matrix X (complete cases per column),
# then z-score. NA rows stay NA. Returns matrix same shape as Y.
residualize_columns <- function(Y, X) {
  Y <- as.matrix(Y); X <- as.matrix(X)
  check(nrow(Y) == nrow(X), "residualize_columns: Y and X have same rows")
  ok_x <- stats::complete.cases(X)
  out <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  qrx <- qr(X[ok_x, , drop = FALSE])
  for (j in seq_len(ncol(Y))) {
    idx <- which(ok_x & !is.na(Y[, j]))
    if (length(idx) < ncol(X) + 2) next
    if (length(idx) == sum(ok_x)) {
      r <- qr.resid(qrx, Y[idx, j])
    } else {
      r <- qr.resid(qr(X[idx, , drop = FALSE]), Y[idx, j])
    }
    s <- sd(r)
    out[idx, j] <- if (is.finite(s) && s > 0) (r - mean(r)) / s else NA_real_
  }
  out
}

# Build technical design matrix from a covariate data.table (factors expanded).
technical_design <- function(cov, vars) {
  parts <- list(intercept = rep(1, nrow(cov)))
  for (v in vars) {
    if (v == "age2")    { parts[[v]] <- cov$age^2; next }
    if (v == "age_sex") { parts[[v]] <- cov$age * (cov$sex == 2); next }
    if (v %in% c("spectrometer", "genotype_batch", "shipment_plate")) {
      f <- factor(cov[[v]])
      if (nlevels(f) > 1) {
        mm <- model.matrix(~ f)[, -1, drop = FALSE]
        # model.matrix drops NA rows; re-expand
        full <- matrix(NA_real_, nrow(cov), ncol(mm)); full[!is.na(f), ] <- mm
        colnames(full) <- paste0(v, "_", seq_len(ncol(mm)))
        for (k in seq_len(ncol(full))) parts[[colnames(full)[k]]] <- full[, k]
      }
      next
    }
    parts[[v]] <- as.numeric(cov[[v]])
  }
  D <- do.call(cbind, parts)
  # drop non-intercept columns with zero variance (e.g. a QC flag nobody carries in this subset)
  sds <- apply(D, 2, sd, na.rm = TRUE)
  D[, colnames(D) == "intercept" | (is.finite(sds) & sds > 0), drop = FALSE]
}

# ---------------------------------------------------------------------------
# Decomposition on the PRS (vectorised; per-column OLS with intercept)
# R: n x p matrix (technically-residualised, z-scored metabolites); prs: length n.
# Returns list(beta, r2, fitted, resid, mean_r, mean_prs) so that the same
# training-set parameters can be applied to a test set.
# ---------------------------------------------------------------------------
decompose_on_prs <- function(R, prs) {
  R <- as.matrix(R)
  check(length(prs) == nrow(R), "decompose_on_prs: prs length matches rows")
  mean_prs <- mean(prs); pc <- prs - mean_prs
  mean_r <- colMeans(R, na.rm = TRUE)
  Rc <- sweep(R, 2, mean_r)
  Rc0 <- Rc; Rc0[is.na(Rc0)] <- 0
  beta <- as.vector(crossprod(Rc0, pc)) / sum(pc^2)
  fitted <- outer(pc, beta)
  resid  <- Rc - fitted
  r2 <- beta^2 * var(prs) / apply(R, 2, var, na.rm = TRUE)
  list(beta = setNames(beta, colnames(R)), r2 = setNames(r2, colnames(R)),
       fitted = fitted, resid = resid, mean_r = mean_r, mean_prs = mean_prs)
}

apply_decomposition <- function(dec, R_new, prs_new) {
  Rc <- sweep(as.matrix(R_new), 2, dec$mean_r)
  Rc - outer(prs_new - dec$mean_prs, dec$beta)
}

# ---------------------------------------------------------------------------
# Fold assignment: stratified on event status, deterministic given seed
# ---------------------------------------------------------------------------
make_folds <- function(event, k, seed) {
  set.seed(seed)
  f <- integer(length(event))
  for (e in unique(event)) {
    idx <- which(event == e)
    f[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  f
}

# ---------------------------------------------------------------------------
# Nested-CV residual metabolome score.
# R: n x p matrix (complete cases only), prs, time, event, fold (outer).
# Within each outer training set: PRS-decomposition -> z-score -> cv.glmnet Cox.
# Out-of-fold linear predictor returned; nothing in the test fold touches fitting.
# ---------------------------------------------------------------------------
fit_residual_score_cv <- function(R, prs, time, event, fold, alpha = 0.5,
                                  inner_folds = 5, lambda_rule = "lambda.min", seed = 42,
                                  standardize_within_test_fold = TRUE) {
  R <- as.matrix(R)
  n <- nrow(R)
  check(all(stats::complete.cases(R)), "residual score input matrix is complete")
  check(length(prs) == n && length(time) == n && length(event) == n && length(fold) == n, "residual score inputs aligned")
  check(all(time > 0), "all follow-up times positive")
  folds <- sort(unique(fold))
  mrs <- rep(NA_real_, n)
  coefs <- matrix(NA_real_, ncol(R), length(folds), dimnames = list(colnames(R), paste0("fold", folds)))
  lambdas <- setNames(numeric(length(folds)), paste0("fold", folds))
  n_nonzero <- lambdas
  for (k in folds) {
    tr <- which(fold != k); te <- which(fold == k)
    check(length(intersect(tr, te)) == 0, "fold %d: train/test disjoint", k)
    dec  <- decompose_on_prs(R[tr, , drop = FALSE], prs[tr])
    Xtr  <- dec$resid
    sd_tr <- apply(Xtr, 2, sd)
    keep <- which(is.finite(sd_tr) & sd_tr > 0)
    Xtr  <- sweep(Xtr[, keep, drop = FALSE], 2, sd_tr[keep], "/")
    Xte  <- apply_decomposition(dec, R[te, , drop = FALSE], prs[te])
    Xte  <- sweep(Xte[, keep, drop = FALSE], 2, sd_tr[keep], "/")
    set.seed(seed + k)
    inner_id <- make_folds(event[tr], inner_folds, seed + 100 * k)
    cvfit <- capture_warnings(cv.glmnet(Xtr, Surv(time[tr], event[tr]), family = "cox", alpha = alpha,
                       foldid = inner_id, standardize = FALSE), where = sprintf("cv.glmnet outer fold %d", k))
    lam_pos <- match(cvfit[[lambda_rule]], cvfit$lambda)
    log_msg("  outer fold %d: lambda path length %d, selected index %d (edge of path = %s)", k, length(cvfit$lambda), lam_pos, lam_pos == length(cvfit$lambda))
    lam <- cvfit[[lambda_rule]]
    lp  <- as.vector(predict(cvfit, newx = Xte, s = lam, type = "link"))
    cf  <- as.vector(coef(cvfit, s = lam))
    coefs[keep, paste0("fold", k)] <- cf
    lambdas[paste0("fold", k)] <- lam
    n_nonzero[paste0("fold", k)] <- sum(cf != 0)
    if (standardize_within_test_fold) {
      s <- sd(lp)
      lp <- if (is.finite(s) && s > 0) (lp - mean(lp)) / s else lp - mean(lp)
    }
    mrs[te] <- lp
    log_msg("  outer fold %d: n_train=%d n_test=%d lambda=%.4g nonzero=%d", k, length(tr), length(te), lam, sum(cf != 0))
  }
  check(all(!is.na(mrs)), "every individual received an out-of-fold score")
  check(sd(mrs) > 0, "out-of-fold score has non-zero variance")
  list(mrs = mrs, coefs = coefs, lambdas = lambdas, n_nonzero = n_nonzero)

}

# ---------------------------------------------------------------------------
# Cox tests
# ---------------------------------------------------------------------------
cox_table <- function(fit, terms = NULL) {
  s <- summary(fit)$coefficients
  ci <- exp(confint(fit))
  dt <- data.table(term = rownames(s), beta = s[, "coef"], se = s[, "se(coef)"],
                   hr = s[, "exp(coef)"], hr_lo = ci[, 1], hr_hi = ci[, 2],
                   z = s[, "z"], p = s[, "Pr(>|z|)"])
  if (!is.null(terms)) dt <- dt[term %in% terms]
  dt
}

# Primary multiplicative interaction model, stratified on the outer CV fold.
fit_interaction_cox <- function(dt, covars, strata_var = "fold", exposure = "mrs_z") {
  rhs <- paste(c(sprintf("prs_z * %s", exposure), covars, sprintf("strata(%s)", strata_var)), collapse = " + ")
  f <- as.formula(sprintf("Surv(time_analysis, event) ~ %s", rhs))
  coxph(f, data = dt, ties = "efron")
}

# Li & Ji (2005) effective number of independent tests from a correlation matrix
meff_li_ji <- function(C) {
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  ev <- abs(ev)
  sum(as.integer(ev >= 1) + (ev - floor(ev)))
}

# ---------------------------------------------------------------------------
# CHIP (Vlasschaert 2023 Blood curated calls, UKB return format)
#   ids_file:     one column of exome-screened participant IDs (the denominator)
#   variant_file: Participant ID, n variants, variant | Array 0..6, VAF | Array 0..6
# Returns one row per screened eid; carriers flagged, everyone else screened
# and negative. Non-screened people are NOT in the table (unknown, not 0).
# ---------------------------------------------------------------------------
CHIP_GENE_GROUPS <- list(
  DNMT3A = "DNMT3A", TET2 = "TET2", ASXL1 = "ASXL1",
  spliceosome = c("SF3B1", "SRSF2", "U2AF1", "ZRSR2"),
  JAK2 = "JAK2", PPM1D = "PPM1D", TP53 = "TP53"
)

parse_chip <- function(ids_file, variant_file, large_vaf = 0.10, min_screened = 1000) {
  ids <- fread(ids_file, header = TRUE, colClasses = "character")
  eid_scr <- unique(as.integer(ids[[1]]))
  eid_scr <- eid_scr[!is.na(eid_scr)]
  check(length(eid_scr) >= min_screened, "CHIP screened-ID list has >= %d IDs (%d)", min_screened, length(eid_scr))
  v <- fread(variant_file, colClasses = "character")
  setnames(v, 1, "eid"); v[, eid := as.integer(eid)]
  var_cols <- grep("variant  \\| Array|variant \\| Array", names(v), value = TRUE)
  var_cols <- var_cols[!grepl("VAF", var_cols)]
  vaf_cols <- grep("VAF", names(v), value = TRUE)
  check(length(var_cols) == length(vaf_cols) && length(var_cols) >= 1, "CHIP variant and VAF array columns pair up (%d/%d)", length(var_cols), length(vaf_cols))
  long <- rbindlist(lapply(seq_along(var_cols), function(k) {
    data.table(eid = v$eid, variant = v[[var_cols[k]]], vaf = suppressWarnings(as.numeric(v[[vaf_cols[k]]])))
  }))
  long <- long[!is.na(variant) & variant != ""]
  long[, gene := sub(":.*", "", variant)]
  check(all(long$eid %in% eid_scr), "every CHIP carrier is in the screened-ID list")
  check(all(long$vaf >= 0.01 & long$vaf <= 1, na.rm = TRUE), "CHIP VAFs lie in [0.01, 1]")
  per <- long[, .(chip_n_variants = .N, chip_max_vaf = max(vaf, na.rm = TRUE), genes = list(unique(gene))), by = eid]
  out <- data.table(eid = eid_scr, wes_screened = 1L)
  out <- merge(out, per, by = "eid", all.x = TRUE)
  out[is.na(chip_n_variants), chip_n_variants := 0L]
  out[, chip_any := as.integer(chip_n_variants > 0)]
  out[, chip_large := as.integer(chip_any == 1 & !is.na(chip_max_vaf) & chip_max_vaf >= large_vaf)]
  for (g in names(CHIP_GENE_GROUPS)) {
    gs <- CHIP_GENE_GROUPS[[g]]
    set(out, j = paste0("chip_", g), value = as.integer(vapply(out$genes, function(x) !is.null(x) && any(x %in% gs), logical(1))))
  }
  known <- unlist(CHIP_GENE_GROUPS)
  out[, chip_other := as.integer(chip_any == 1 & vapply(genes, function(x) !is.null(x) && !any(x %in% known), logical(1)))]
  out[, genes := NULL]
  check(sum(out$chip_any) == length(unique(long$eid)), "carrier count equals distinct carrier IDs")
  out
}
