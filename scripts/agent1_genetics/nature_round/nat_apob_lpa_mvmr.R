#!/usr/bin/env Rscript
# ============================================================
# Nature round, analysis A: does the Lp(a) effect survive conditioning on ApoB?
#
# Why this script exists
# ----------------------
# Every published model conditions Lp(a) on LDL-C. Lp(a) particles carry apoB-100, so
# an LPA-raising allele also raises measured ApoB. A reviewer's first question is what
# is left of the Lp(a) coefficient once ApoB, not LDL-C, is the competing exposure.
# PREREGISTERED_READOUT.md section A fixes how the answer is read; this script applies
# those thresholds literally and does not get to choose a favourable instrument set.
#
# Order of operations, and why it is in this order
# ------------------------------------------------
#  1. POSITIVE CONTROL FIRST. The LDL-C analogue must reproduce the published CAD
#     estimates (LDL-C OR 1.25, Lp(a) OR 1.15, +/-0.01) on the published instrument
#     rule. If not, the run exits non-zero having computed no ApoB estimate at all: the
#     readout says no estimate is read until the control is fixed, and an estimate that
#     was never produced cannot be read.
#  2. Models x instrument variants x outcomes, written long.
#  3. Bootstrap of the Lp(a) log-OR difference between the ApoB and LDL-C models on a
#     MATCHED instrument set, so the percentage attenuation has an interval.
#  4. Readout verdict.
#
# Instrument variants
#   published        run_mvmr.R: per-trait 10 Mb greedy, union, second-pass 10 Mb on the
#                    max LOG10P across the model's exposures. Leaves the LPA locus as a
#                    single variant, so the Lp(a) coefficient rests on rs10455872.
#   cis_dense        mvmr_covariance.R: polygenic exposures as above, Lp(a) replaced by
#                    the PLINK-clumped cis list, window variants exempt from pruning.
#   noLPAwin_nonLPA  ApoB/LDL-C/HDL-C/TG may not select an instrument from the LPA
#                    window; Lp(a) keeps its own. The narrow reading of the readout.
#   noLPAwin_all     Nothing from the window at all. Lp(a) rests on trans variants only,
#                    so this asks a different question and must not be read as the
#                    primary effect. Emits NA if Lp(a) has no GWS trans instrument.
#   cis_noLPAwin_nonLPA  cis_dense, but ApoB/LDL-C/HDL-C/TG may not add a window variant.
#                    This is where the exclusion can bite: under the published rule any
#                    window variant ApoB selects is pruned by rs10455872 at the second
#                    pass anyway, whereas cis_dense exempts the window from pruning, so a
#                    trait's own window lead (possibly in LD with the cis list) stays in.
#
# Usage: see slurm/nat_A_apob_lpa.slurm. Summary-level output only: the phenotype matrix
# is read solely to form an exposure correlation matrix and is never written out.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(MVMR)
  library(optparse)
  library(parallel)
})
set.seed(42)

# ---- constants ------------------------------------------------------------
LOG10P_GWS <- -log10(5e-8)
KB_CLUMP   <- 10000L
# GRCh37/hg19, matching the REGENIE sumstats (rs10455872 = chr6:161,010,118).
LPA_CHR <- "6"; LPA_LO <- 160500000L; LPA_HI <- 161700000L
# What the noLPAwin_* variants exclude from selection. Wider than the window, the cis region of
# run_mr.R: rs531596603 (chr6:160,179,747, Lp(a) LOG10P 438) lies outside 160.5-161.7 Mb and is
# plainly locus-driven, so the narrow window would label a cis-driven estimate "trans only".
LPA_EXCL_LO <- 157000000L; LPA_EXCL_HI <- 163000000L
LPA_LEAD <- "rs10455872"
LPA_REPORT <- c("rs10455872", "rs3798220")
EXPOSURES_ALL <- c("LDL_C", "APOB", "HDL_C", "TRIGLY", "LPA")

# Exposure order matters only for ties and for the Egger orientation (first exposure),
# so LDL_C / APOB go first, as in run_mvmr.R.
MODELS <- list(
  LDL4  = c("LDL_C", "HDL_C", "TRIGLY", "LPA"),
  APOB4 = c("APOB",  "HDL_C", "TRIGLY", "LPA"),
  LDL3  = c("LDL_C", "TRIGLY", "LPA"),
  APOB3 = c("APOB",  "TRIGLY", "LPA"),
  LDL2  = c("LDL_C", "LPA"),
  APOB2 = c("APOB",  "LPA")
)
# Each ApoB model is compared with its same-structure LDL-C analogue. A 2-exposure ApoB
# model compared with the 4-exposure LDL-C value would mix up "ApoB vs LDL-C" with
# "fewer covariates".
REF_OF <- c(APOB4 = "LDL4", APOB3 = "LDL3", APOB2 = "LDL2")
VARIANTS <- c("published", "cis_dense", "noLPAwin_nonLPA", "noLPAwin_all", "cis_noLPAwin_nonLPA")
BOOT_VARIANTS <- c("published", "cis_dense")

# ---- small helpers --------------------------------------------------------
in_lpa_window <- function(chrom, pos) {
  as.character(chrom) == LPA_CHR & pos >= LPA_LO & pos <= LPA_HI
}
in_lpa_exclusion <- function(chrom, pos) {
  as.character(chrom) == LPA_CHR & pos >= LPA_EXCL_LO & pos <= LPA_EXCL_HI
}

log_msg <- function(...) message(sprintf("[natA %s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# Exact copy of the corrected rule in run_mvmr.R, including both bug fixes: rank on
# LOG10P (10^-LOG10P underflows to 0 above ~323 and would pick by file order), and the
# window test on same-chromosome positions only (the old guard recycled two vectors of
# different lengths and let correlated variants through). Plain vectors instead of
# dt$col[i] only because LDL-C has hundreds of thousands of GWS rows.
distance_prune <- function(ids, chrs, poss, kb_window = KB_CLUMP) {
  keep_snps <- character(0); keep_pos <- integer(0); keep_chr <- character(0)
  win <- as.integer(kb_window) * 1000L
  for (i in seq_along(ids)) {
    same_chr <- keep_chr == chrs[i]
    if (length(keep_pos) == 0L || !any(abs(keep_pos[same_chr] - poss[i]) < win)) {
      keep_snps <- c(keep_snps, ids[i]); keep_pos <- c(keep_pos, poss[i]); keep_chr <- c(keep_chr, chrs[i])
    }
  }
  keep_snps
}

greedy_clump <- function(dt, kb_window = KB_CLUMP) {
  gws <- dt[LOG10P > LOG10P_GWS & !is.na(BETA) & !is.na(SE)]
  if (nrow(gws) == 0L) return(character(0))
  gws <- gws[order(-LOG10P)]
  distance_prune(gws$ID, as.character(gws$CHROM), as.integer(gws$GENPOS), kb_window)
}

# One reference row per ID carrying its strongest LOG10P across the traits. Same
# construction and same tie order as run_mvmr.R (trait order, then file order).
ref_rows <- function(ids, traits, rows) {
  r <- rbindlist(lapply(traits, function(tr) rows[[tr]][ID %in% ids, .(ID, CHROM, GENPOS, LOG10P)]))
  if (!nrow(r)) return(r)
  r <- r[r[, .I[which.max(LOG10P)], by = ID]$V1]
  r[order(-LOG10P)]
}

second_pass <- function(ref, kb_window = KB_CLUMP) {
  distance_prune(ref$ID, as.character(ref$CHROM), as.integer(ref$GENPOS), kb_window)
}

# ---- instrument rules -----------------------------------------------------
select_instruments <- function(variant, traits, rows, ivs_full, ivs_nowin, cis_ids) {
  non_lpa <- setdiff(traits, "LPA")
  if (variant == "published") {
    ids <- unique(unlist(ivs_full[traits]))
    return(second_pass(ref_rows(ids, traits, rows)))
  }
  if (variant == "noLPAwin_nonLPA") {
    ids <- unique(c(unlist(ivs_nowin[non_lpa]), ivs_full[["LPA"]]))
    return(second_pass(ref_rows(ids, traits, rows)))
  }
  if (variant == "noLPAwin_all") {
    ids <- unique(unlist(ivs_nowin[traits]))
    return(second_pass(ref_rows(ids, traits, rows)))
  }
  if (variant %in% c("cis_dense", "cis_noLPAwin_nonLPA")) {
    # mvmr_covariance.R: a 10 Mb rule collapses the ~1.1 Mb LPA locus to one variant,
    # so window variants are exempt from the second pass. Lp(a)'s own trans hits are
    # not added, which is why this set is not a superset of "published".
    src <- if (variant == "cis_dense") ivs_full else ivs_nowin
    ids <- unique(c(unlist(src[non_lpa]), cis_ids))
    r <- ref_rows(ids, traits, rows)
    incis <- r[in_lpa_window(CHROM, GENPOS), ID]
    kept <- second_pass(r[!ID %in% incis])
    return(c(kept, incis))
  }
  stop("unknown instrument variant: ", variant)
}

# ---- extraction (never load a whole GWAS into R) ---------------------------
AWK_PROGRAM <- '
BEGIN {
  OFS = "\\t"; nw = split(want, w, ",")
  if (idfile != "") { while ((getline line < idfile) > 0) ids[line] = 1; close(idfile) }
}
NR == 1 {
  for (i = 1; i <= NF; i++) c[$i] = i
  for (j = 1; j <= nw; j++) {
    if (!(w[j] in c)) { print "missing column: " w[j] > "/dev/stderr"; exit 3 }
    ci[j] = c[w[j]]
  }
  if (thrcol != "") { if (!(thrcol in c)) { print "missing column: " thrcol > "/dev/stderr"; exit 3 }; tc = c[thrcol] }
  if (idcol  != "") { if (!(idcol  in c)) { print "missing column: " idcol  > "/dev/stderr"; exit 3 }; ic = c[idcol] }
  out = w[1]; for (j = 2; j <= nw; j++) out = out OFS w[j]
  print out; next
}
{
  if (thrcol != "" && ($tc + 0) <= thr) next
  v = ""
  if (idcol != "") {
    v = $ic
    if (first == "1") { split(v, a, ","); v = a[1] }
    if (!(v in ids)) next
  }
  out = ""
  for (j = 1; j <= nw; j++) {
    x = $ci[j]
    if (idcol != "" && ci[j] == ic) x = v
    out = (j == 1) ? x : out OFS x
  }
  print out
}'

run_shell <- function(cmd) {
  # bash + pipefail: a truncated .gz must be an error, not a silently short table.
  st <- system2("bash", c("-c", shQuote(paste("set -o pipefail;", cmd))))
  if (!identical(as.integer(st), 0L)) stop("shell pipeline failed (status ", st, "): ", cmd)
  invisible(TRUE)
}

reader_cmd <- function(f) if (grepl("\\.gz$", f)) paste("gunzip -c", shQuote(f)) else paste("cat", shQuote(f))

# Cache key = input identity + extraction parameters. Stale build artefacts are this
# project's recurring failure mode, so a changed input can never reuse an old table.
cache_file <- function(cache_dir, tag, ...) {
  tmp <- tempfile(); on.exit(unlink(tmp))
  writeLines(paste(..., sep = "|"), tmp)
  file.path(cache_dir, sprintf("%s.%s.tsv", tag, unname(tools::md5sum(tmp))))
}

extract_table <- function(in_file, out_file, want, fs = "", thrcol = "", thr = 0,
                          idfile = "", idcol = "", first = "0") {
  if (file.exists(out_file) && file.exists(paste0(out_file, ".ok"))) return(out_file)
  tmp_prog <- tempfile(fileext = ".awk"); writeLines(AWK_PROGRAM, tmp_prog)
  tmp_out <- paste0(out_file, ".part.", Sys.getpid())
  on.exit(unlink(c(tmp_prog, tmp_out)))
  awk <- paste0("awk ", if (nzchar(fs)) paste0("-F", shQuote(fs), " ") else "",
                "-v want=", shQuote(paste(want, collapse = ",")),
                " -v thrcol=", shQuote(thrcol), " -v thr=", format(thr, digits = 15),
                " -v idfile=", shQuote(idfile), " -v idcol=", shQuote(idcol),
                " -v first=", shQuote(first), " -f ", shQuote(tmp_prog))
  run_shell(paste(reader_cmd(in_file), "|", awk, ">", shQuote(tmp_out)))
  file.rename(tmp_out, out_file); file.create(paste0(out_file, ".ok"))
  out_file
}

parallel_stop_on_error <- function(jobs, f, threads) {
  # mclapply returns a try-error (or NULL if a worker died) rather than raising: both
  # must stop the run, or a failed extraction silently becomes a missing column.
  res <- mclapply(jobs, f, mc.cores = max(1L, threads), mc.preschedule = FALSE)
  bad <- vapply(res, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
  if (any(bad)) stop("extraction failed: ", paste(names(jobs)[bad], collapse = ", "), "\n",
                     paste(vapply(res[bad], function(x) as.character(x)[1], ""), collapse = "\n"))
  res
}

# Allele columns are forced to character: in a short extract a column holding only "T"
# would otherwise be read as logical TRUE and no allele would ever match.
read_tsv <- function(f, chr_cols) {
  fread(f, sep = "\t", quote = "", na.strings = c("NA", ".", ""), showProgress = FALSE,
        colClasses = list(character = chr_cols))
}

REGENIE_COLS <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")

load_regenie <- function(f) {
  d <- read_tsv(f, c("CHROM", "ID", "ALLELE0", "ALLELE1"))
  d <- d[!is.na(ID)]
  d[, CHROM := as.character(CHROM)]
  d[, GENPOS := as.integer(GENPOS)]
  d[]
}

EXTERNAL_FORMATS <- list(
  finngen    = list(want = c("rsids", "alt", "ref", "beta", "sebeta", "af_alt"), fs = "\t",
                    idcol = "rsids", first = "1"),
  cardiogram = list(want = c("markername", "effect_allele", "noneffect_allele", "beta", "se_dgc", "effect_allele_freq"),
                    fs = "", idcol = "markername", first = "0")
)

# ---- alignment -------------------------------------------------------------
# Same strand within UKB, so the allele check is exact, as in run_mvmr.R. Duplicate rsIDs
# (multi-allelic sites) keep the row whose alleles match the reference, so a merge can
# never expand and count one variant twice.
align_ukb <- function(ref, d) {
  m <- merge(ref[, .(ID, rA1, rA0)], d[, .(ID, ALLELE1, ALLELE0, BETA, SE, LOG10P)], by = "ID")
  m[, sgn := fifelse(ALLELE1 == rA1 & ALLELE0 == rA0, 1, fifelse(ALLELE1 == rA0 & ALLELE0 == rA1, -1, NA_real_))]
  m <- m[!is.na(sgn)]
  m <- m[!duplicated(ID)]
  m[, b := BETA * sgn]
  m[, .(ID, b, se = SE, LOG10P)]
}

is_palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
flip_base <- function(x) chartr("ACGT", "TGCA", x)

# Harmonization logic reused from run_external_mr.R: exact / flipped / strand-flipped;
# palindromes resolved by frequency, dropped when ambiguous (EAF 0.42-0.58); one row per
# rsID after alignment because FinnGen lists multi-allelic rsIDs on several rows.
harmonize_external <- function(ref, ext) {
  m <- merge(ref[, .(ID, rA1 = toupper(rA1), rA0 = toupper(rA0), rEAF)],
             ext[, .(ID, oEA = toupper(oEA), oOA = toupper(oOA), oBETA, oSE, oEAF)], by = "ID")
  if (!nrow(m)) return(data.table(ID = character(0), b = numeric(0), se = numeric(0)))
  m[, palin := is_palindromic(rA1, rA0)]
  m[, aligned := fifelse(rA1 == oEA & rA0 == oOA, 1,
                 fifelse(rA1 == oOA & rA0 == oEA, -1,
                 fifelse(rA1 == flip_base(oEA) & rA0 == flip_base(oOA), 1,
                 fifelse(rA1 == flip_base(oOA) & rA0 == flip_base(oEA), -1, NA_real_))))]
  # Ambiguous in EITHER source is dropped: UKB 0.40 with FinnGen 0.53 would otherwise be flipped
  # although its alleles match. run_external_mr.R tests the UKB frequency only.
  amb_o <- m[palin == TRUE & !is.na(aligned) & !is.na(oEAF) & oEAF > 0.42 & oEAF < 0.58 & !is.na(rEAF) & !(rEAF > 0.42 & rEAF < 0.58), uniqueN(ID)]
  if (amb_o > 0L) message(sprintf("[harmonise] %d palindromes ambiguous only in the outcome, dropped; run_external_mr.R keeps them", amb_o))
  m[palin == TRUE & !is.na(aligned), aligned := fifelse(
    is.na(rEAF) | is.na(oEAF) | (rEAF > 0.42 & rEAF < 0.58) | (oEAF > 0.42 & oEAF < 0.58), NA_real_,
    fifelse((rEAF - 0.5) * (oEAF - 0.5) > 0, 1, -1))]
  m <- m[!is.na(aligned)]
  m <- m[!duplicated(ID)]
  m[, .(ID, b = oBETA * aligned, se = oSE)]
}

# Assemble one model's design. Rows dropped for missing data in any exposure or the
# outcome are dropped for all columns, so the matrices can never be misaligned.
build_design <- function(ids, traits, rows, outcome) {
  ref0 <- rows[[traits[1]]][ID %in% ids & !duplicated(ID)]
  ref <- ref0[, .(ID, rA1 = ALLELE1, rA0 = ALLELE0, rEAF = A1FREQ, CHROM, GENPOS)]
  ids <- ids[ids %in% ref$ID]
  mk <- function() matrix(NA_real_, length(ids), length(traits), dimnames = list(ids, traits))
  BX <- mk(); BXSE <- mk(); LP <- mk()
  for (tr in traits) {
    a <- align_ukb(ref[ID %in% ids], rows[[tr]][ID %in% ids])
    i <- match(a$ID, ids)
    BX[i, tr] <- a$b; BXSE[i, tr] <- a$se; LP[i, tr] <- a$LOG10P
  }
  o <- if (outcome$type == "ukb") align_ukb(ref[ID %in% ids], outcome$dt)
       else harmonize_external(ref[ID %in% ids], outcome$dt)
  BY <- BYSE <- rep(NA_real_, length(ids))
  i <- match(o$ID, ids)
  keep <- !is.na(i)
  BY[i[keep]] <- o$b[keep]; BYSE[i[keep]] <- o$se[keep]
  ok <- complete.cases(BX) & complete.cases(BXSE) & !is.na(BY) & !is.na(BYSE) & BYSE > 0
  r <- ref[match(ids[ok], ID)]
  list(BX = BX[ok, , drop = FALSE], BXSE = BXSE[ok, , drop = FALSE], LOG10P = LP[ok, , drop = FALSE],
       BY = BY[ok], BYSE = BYSE[ok], ids = ids[ok],
       in_window = in_lpa_window(r$CHROM, r$GENPOS))
}

subset_design <- function(D, traits) {
  D$BX <- D$BX[, traits, drop = FALSE]; D$BXSE <- D$BXSE[, traits, drop = FALSE]
  D$LOG10P <- D$LOG10P[, traits, drop = FALSE]
  D
}

# ---- estimation -------------------------------------------------------------
is_correlation_matrix <- function(P) {
  is.matrix(P) && nrow(P) == ncol(P) && !anyNA(P) && isSymmetric(unname(P), tol = 1e-8) &&
    all(abs(diag(P) - 1) < 1e-8) && all(abs(P) <= 1 + 1e-8)
}

# Identical to mr_mvivw(model = "random"): weighted LS, SE scaled by max(sigma, 1). A
# separate copy because the bootstrap needs it thousands of times; fit_model() checks
# the two agree every time, so it cannot silently drift from the published estimator.
wls_ivw <- function(X, y, se_y) {
  sw <- 1 / se_y
  Xw <- X * sw; yw <- y * sw
  f <- .lm.fit(Xw, yw)
  K <- ncol(X); n <- nrow(X)
  if (f$rank < K || n <= K) return(list(b = rep(NA_real_, K), se = rep(NA_real_, K)))
  b <- f$coefficients
  sigma <- sqrt(sum(f$residuals^2) / (n - K))
  XtXi <- chol2inv(chol(crossprod(Xw)))
  list(b = b, se = sqrt(diag(XtXi)) * max(sigma, 1))
}

model_na <- function(D, traits, reason) {
  data.table(exposure = traits, n_instruments = nrow(D$BX),
             n_lpa_window = sum(D$in_window), lead_in_model = LPA_LEAD %in% D$ids,
             b = NA_real_, se = NA_real_, OR = NA_real_, OR_lo95 = NA_real_, OR_hi95 = NA_real_, p = NA_real_,
             condF_phenocov = NA_real_, condF_gencov0 = NA_real_,
             Q_A_stat_phenocov = NA_real_, Q_A_p_phenocov = NA_real_,
             egger_b = NA_real_, egger_se = NA_real_, egger_intercept = NA_real_, egger_intercept_p = NA_real_,
             median_b = NA_real_, median_se = NA_real_, median_note = "not_fitted",
             status = reason)
}

# Degenerate designs return NA with the reason named, never a number:
#  - fewer than K+2 variants (no residual degrees of freedom worth the name; Egger impossible)
#  - an exposure with no genome-wide-significant instrument in the set
#  - rank deficiency or extreme collinearity: lm() silently drops an aliased column and
#    mr_mvivw() then returns a shorter coefficient vector than there are exposures
#  - the LPA locus absent from the outcome: the Lp(a) coefficient would then rest on
#    trans variants while being reported as if it were the locus effect
# Weighted median (project convention: IVW + median + Egger) is fitted on request only:
# it bootstraps internally, so it runs on the CAD rows, not on every endpoint. Its Lp(a)
# coefficient is emitted as NA with a handful of window variants: the median regression
# ignores a single row among 161, so under the published rule that number would come from
# trans variants alone while being labelled the Lp(a) effect.
MEDIAN_MIN_WINDOW <- 5L
fit_model <- function(D, traits, Pcor, require_lpa_locus = TRUE, median = FALSE) {
  K <- length(traits); n <- nrow(D$BX)
  if (n < K + 2L) return(model_na(D, traits, sprintf("too_few_instruments(%d<%d)", n, K + 2L)))
  hits <- colSums(D$LOG10P > LOG10P_GWS, na.rm = TRUE)
  if (any(hits == 0L)) return(model_na(D, traits, paste0("not_instrumented:", paste(names(hits)[hits == 0L], collapse = "+"))))
  if (require_lpa_locus && "LPA" %in% traits && sum(D$in_window) == 0L)
    return(model_na(D, traits, "lpa_locus_absent"))
  Xw <- D$BX / D$BYSE
  if (qr(Xw)$rank < K) return(model_na(D, traits, "rank_deficient"))
  kp <- kappa(sweep(Xw, 2, sqrt(colSums(Xw^2)), "/"), exact = TRUE)
  if (!is.finite(kp) || kp > 1e6) return(model_na(D, traits, sprintf("collinear(kappa=%.3g)", kp)))
  if (!is_correlation_matrix(Pcor)) stop("Pcor must be a correlation matrix: phenocov_mvmr() scales it by the SEs")
  P <- Pcor[traits, traits, drop = FALSE]

  inp <- mr_mvinput(bx = D$BX, bxse = D$BXSE, by = D$BY, byse = D$BYSE, exposure = traits, outcome = "outcome")
  ivw <- mr_mvivw(inp)
  if (length(ivw@Estimate) != K) return(model_na(D, traits, "mr_mvivw_dropped_a_coefficient"))
  chk <- wls_ivw(D$BX, D$BY, D$BYSE)
  if (max(abs(chk$b - ivw@Estimate)) > 1e-8 || max(abs(chk$se - ivw@StdError)) > 1e-8)
    stop("internal IVW disagrees with mr_mvivw(); the bootstrap would not be the published estimator")
  egg <- tryCatch(mr_mvegger(inp), error = function(e) NULL)
  med_b <- med_se <- rep(NA_real_, K); med_note <- rep("not_fitted", K)
  if (median) {
    med <- tryCatch(mr_mvmedian(inp, iterations = 2000L), error = function(e) NULL)
    if (is.null(med) || length(med@Estimate) != K) {
      med_note[] <- "mr_mvmedian_failed"
    } else {
      med_b <- med@Estimate; med_se <- med@StdError; med_note[] <- "ok"
      if ("LPA" %in% traits && sum(D$in_window) < MEDIAN_MIN_WINDOW) {
        i <- match("LPA", traits)
        med_b[i] <- NA_real_; med_se[i] <- NA_real_
        med_note[i] <- sprintf("lpa_window_variants_lt_%d", MEDIAN_MIN_WINDOW)
      }
    }
  }

  # strength_mvmr()/pleiotropy_mvmr() print to stdout and warn when gencov = 0, by design.
  f0 <- format_mvmr(BXGs = D$BX, BYG = D$BY, seBXGs = D$BXSE, seBYG = D$BYSE, RSID = D$ids)
  gc <- phenocov_mvmr(P, D$BXSE)
  s1 <- NULL; s0 <- NULL; qa <- NULL
  invisible(capture.output({
    s1 <- strength_mvmr(f0, gencov = gc)
    s0 <- suppressWarnings(strength_mvmr(f0, gencov = 0))
    qa <- pleiotropy_mvmr(f0, gencov = gc)
  }))
  data.table(exposure = traits, n_instruments = n, n_lpa_window = sum(D$in_window),
             lead_in_model = LPA_LEAD %in% D$ids,
             b = ivw@Estimate, se = ivw@StdError,
             OR = exp(ivw@Estimate),
             OR_lo95 = exp(ivw@Estimate - 1.96 * ivw@StdError),
             OR_hi95 = exp(ivw@Estimate + 1.96 * ivw@StdError),
             p = ivw@Pvalue,
             condF_phenocov = as.numeric(s1[1, seq_len(K)]),
             condF_gencov0 = as.numeric(s0[1, seq_len(K)]),
             Q_A_stat_phenocov = qa$Qstat, Q_A_p_phenocov = qa$Qpval,
             egger_b = if (is.null(egg)) NA_real_ else egg@Estimate,
             egger_se = if (is.null(egg)) NA_real_ else egg@StdError.Est,
             egger_intercept = if (is.null(egg)) NA_real_ else egg@Intercept,
             egger_intercept_p = if (is.null(egg)) NA_real_ else egg@Pvalue.Int,
             median_b = med_b, median_se = med_se, median_note = med_note,
             status = "ok")
}

# ---- bootstrap of the attenuation -------------------------------------------
# Three flavours, because the obvious one fails here.
#  naive       resample all instruments (the project's precedent, and what was asked).
#              Under the published rule the Lp(a) coefficient rests on ONE variant, which
#              is absent from ~37% of resamples, and those estimates are identified by
#              trans variants only. The interval then describes that mixture, not the
#              attenuation. Reported so the failure is visible, not used for the read.
#  stratified  resample inside and outside the LPA window separately, so the locus is
#              always represented. Holding one variant fixed removes its sampling error,
#              so with a single LPA variant this undercovers (selftest: about 80%).
#  strat_param the interval the readout uses. Outside the window: resample instruments and
#              add nothing, because a pairs bootstrap already carries their sampling error
#              and heterogeneity. Inside the window: keep every variant and draw its
#              sampling error, by ~ N(by, (sigma * byse)^2) and bx with the exposure
#              correlation between the estimates, identically for both models. sigma is
#              max(1, residual SD) from the reference fit, the same overdispersion the
#              IVW standard error carries. An earlier version also added noise outside the
#              window and resampled inside it, which counted the same error twice and
#              covered 100% in the selftest.
# Ignored: covariance between by and bx from UKB sample overlap. FinnGen has none.
boot_attenuation <- function(D, t_ref, t_alt, Pcor, nboot, seed = 42L) {
  traits <- colnames(D$BX)
  stopifnot(all(c(t_ref, t_alt) %in% traits), "LPA" %in% t_ref, "LPA" %in% t_alt)
  ir <- match(t_ref, traits); ia <- match(t_alt, traits)
  kr <- match("LPA", t_ref); ka <- match("LPA", t_alt)
  BX <- D$BX; BXSE <- D$BXSE; BY <- D$BY; BYSE <- D$BYSE
  n <- nrow(BX); nx <- ncol(BX)
  if (n < length(t_ref) + 2L) return(NULL)
  if (!is_correlation_matrix(Pcor)) stop("Pcor must be a correlation matrix")
  U <- chol(Pcor[traits, traits, drop = FALSE])
  w_in <- which(D$in_window); w_out <- which(!D$in_window)
  fit_two <- function(sel, bx, by) {
    sw <- 1 / BYSE[sel]; yw <- by * sw; Xw <- bx * sw
    f1 <- .lm.fit(Xw[, ir, drop = FALSE], yw); f2 <- .lm.fit(Xw[, ia, drop = FALSE], yw)
    if (f1$rank < length(ir) || f2$rank < length(ia)) return(c(NA_real_, NA_real_))
    c(f1$coefficients[kr], f2$coefficients[ka])
  }
  obs <- fit_two(seq_len(n), BX, BY)
  se_ref <- wls_ivw(BX[, ir, drop = FALSE], BY, BYSE)$se[kr]
  se_alt <- wls_ivw(BX[, ia, drop = FALSE], BY, BYSE)$se[ka]
  f_ref <- .lm.fit(BX[, ir, drop = FALSE] / BYSE, BY / BYSE)
  sigma_ref <- max(1, sqrt(sum(f_ref$residuals^2) / (n - length(ir))))
  by_noise <- BYSE * sigma_ref
  n_in <- length(w_in)
  flavours <- c("naive", "stratified", "strat_param")
  out <- vector("list", length(flavours)); names(out) <- flavours
  for (fl in flavours) {
    set.seed(seed)
    B <- matrix(NA_real_, nboot, 2L); nowin <- 0L
    for (bi in seq_len(nboot)) {
      re_out <- if (length(w_out)) w_out[sample.int(length(w_out), length(w_out), replace = TRUE)]
      sel <- switch(fl,
        naive = sample.int(n, n, replace = TRUE),
        stratified = c(if (n_in) w_in[sample.int(n_in, n_in, replace = TRUE)], re_out),
        strat_param = c(w_in, re_out))
      bx <- BX[sel, , drop = FALSE]; by <- BY[sel]
      if (fl == "strat_param" && n_in) {
        # The first n_in rows of sel are the window variants, each exactly once.
        Z <- matrix(rnorm(n_in * nx), n_in, nx) %*% U
        bx[seq_len(n_in), ] <- bx[seq_len(n_in), , drop = FALSE] + BXSE[w_in, , drop = FALSE] * Z
        by[seq_len(n_in)] <- by[seq_len(n_in)] + by_noise[w_in] * rnorm(n_in)
      }
      if (!any(D$in_window[sel])) nowin <- nowin + 1L
      B[bi, ] <- fit_two(sel, bx, by)
    }
    d <- B[, 2] - B[, 1]
    pct <- 100 * (B[, 1] - B[, 2]) / B[, 1]
    use <- is.finite(d) & is.finite(pct)
    ref_nonpos <- mean(B[use, 1] <= 0)
    usable <- sum(use)
    # A percentile interval on a ratio whose denominator crosses zero is meaningless:
    # it would look precise while being driven by near-zero divisions. NA instead.
    pct_ok <- usable >= 0.9 * nboot && ref_nonpos <= 0.025
    qd <- if (usable > 0L) quantile(d[use], c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
    qp <- if (pct_ok) quantile(pct[use], c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
    out[[fl]] <- data.table(flavour = fl, nboot = nboot, n_usable = usable,
                            frac_no_window_variant = nowin / nboot, frac_ref_nonpositive = ref_nonpos,
                            diff_se = if (usable > 1L) sd(d[use]) else NA_real_,
                            diff_lo95 = qd[1], diff_hi95 = qd[2],
                            atten_pct_lo95 = qp[1], atten_pct_hi95 = qp[2])
  }
  res <- rbindlist(out)
  res[, `:=`(n_matched = n, n_matched_window = length(w_in),
             b_ref_matched = obs[1], se_ref_matched = se_ref,
             b_alt_matched = obs[2], se_alt_matched = se_alt,
             sigma_ref = sigma_ref, diff_obs = obs[2] - obs[1],
             atten_pct_obs = 100 * (obs[1] - obs[2]) / obs[1])]
  res[]
}

# ---- the readout, literally --------------------------------------------------
# Thresholds from PREREGISTERED_READOUT.md section A. "Within 25%" is two-sided in the
# text, and the rules say nothing about the log OR RISING by more than 25%, so that case
# is named rather than folded into SUPPORTS.
cohort_category <- function(b_ref, b, lo, hi) {
  if (anyNA(c(b_ref, b, lo, hi))) return("NOT_READABLE")
  if (b_ref <= 0) return("NOT_READABLE")   # "falls by" has no meaning without a positive reference
  if (lo <= 1 && hi >= 1) return("WEAKENED")
  a <- 100 * (b_ref - b) / b_ref
  if (a > 50) return("WEAKENED")
  if (a > 25) return("ATTENUATED")
  if (a >= -25) return("SUPPORTS")
  "OUTSIDE_RULES_INCREASE"
}

combine_categories <- function(cats) {
  # "WEAKENED: ... in either cohort", so it outranks an unreadable cohort.
  for (k in c("WEAKENED", "NOT_READABLE", "ATTENUATED", "OUTSIDE_RULES_INCREASE", "SUPPORTS"))
    if (k %in% cats) return(k)
  stop("no categories to combine")
}

atten_crosses_threshold <- function(lo, hi) {
  if (anyNA(c(lo, hi))) return(NA)
  any(vapply(c(-25, 25, 50), function(t) lo < t && hi > t, logical(1)))
}

# ---- per-allele effects of the Lp(a) instruments ---------------------------------
# The whole question exists because an Lp(a)-raising allele raises ApoB too: this table
# shows by how much, in SD units per allele, oriented to the Lp(a)-raising allele.
lpa_allele_effects <- function(ids, rows) {
  ref <- rows[["LPA"]][ID %in% ids & !duplicated(ID)]
  if (!nrow(ref)) return(data.table())
  ref <- ref[, .(ID, CHROM, GENPOS, rA1 = ALLELE1, rA0 = ALLELE0, rEAF = A1FREQ)]
  got <- lapply(c("LPA", "APOB", "LDL_C", "HDL_C", "TRIGLY"), function(tr) {
    a <- align_ukb(ref, rows[[tr]])[, .(ID, b, se)]
    setnames(a, c("b", "se"), paste0(c("beta_", "se_"), tr)); a
  })
  res <- Reduce(function(x, y) merge(x, y, by = "ID", all.x = TRUE), got, init = ref)
  res[, sgn := sign(beta_LPA)]
  for (tr in c("LPA", "APOB", "LDL_C", "HDL_C", "TRIGLY")) {
    col <- paste0("beta_", tr); res[, (col) := get(col) * sgn]
  }
  res[, `:=`(raising_allele = fifelse(sgn > 0, rA1, rA0), other_allele = fifelse(sgn > 0, rA0, rA1),
             raising_eaf = fifelse(sgn > 0, rEAF, 1 - rEAF))]
  # Approximate SE of a ratio from the delta method, ignoring the covariance of the two
  # betas. At the lead variant the denominator is so precise that this hardly matters.
  ratio_se <- function(num, se_num, den, se_den) abs(num / den) * sqrt((se_num / num)^2 + (se_den / den)^2)
  res[, `:=`(apob_sd_per_lpa_sd = beta_APOB / beta_LPA,
             apob_ratio_se = ratio_se(beta_APOB, se_APOB, beta_LPA, se_LPA),
             ldl_sd_per_lpa_sd = beta_LDL_C / beta_LPA,
             ldl_ratio_se = ratio_se(beta_LDL_C, se_LDL_C, beta_LPA, se_LPA))]
  res[, in_lpa_window := in_lpa_window(CHROM, GENPOS)]
  res[, c("sgn", "rA1", "rA0", "rEAF") := NULL]
  setcolorder(res, c("ID", "CHROM", "GENPOS", "raising_allele", "other_allele", "raising_eaf"))
  res[order(-abs(beta_LPA / se_LPA))]
}

# ---- main -------------------------------------------------------------------
parse_opts <- function() {
  parse_args(OptionParser(option_list = list(
    make_option("--sumstat_dir", type = "character", help = "Arm's merged sumstats (lipids_*.txt.gz, outcomes_*.txt.gz)"),
    make_option("--fallback_dir", type = "character", default = "",
                help = "Used per file when absent from --sumstat_dir; every use is flagged in the output"),
    make_option("--pheno", type = "character", help = "Phenotype matrix with exposure columns named by trait ID"),
    make_option("--arm", type = "character", default = "statinfix"),
    make_option("--arm_type", type = "character", default = "primary", help = "primary | sensitivity"),
    make_option("--control_stamp", type = "character", default = "",
                help = "sensitivity arms: positive_control.PASSED written by the primary arm"),
    make_option("--outcomes", type = "character", default = "CAD,MI,REVASC,STROKE,HF,CV_DEATH,MACE,PAD",
                help = "UKB endpoints; default = the paper's eight"),
    make_option("--finngen", type = "character", default = ""),
    make_option("--cardiogram", type = "character", default = ""),
    make_option("--cis_list", type = "character", help = "cis-LPA rsIDs, one per line"),
    make_option("--keep", type = "character", default = "", help = "FID IID list the GWAS of this arm was restricted to; the exposure correlation is computed on it"),
    make_option("--nboot", type = "integer", default = 2000L),
    make_option("--threads", type = "integer", default = 4L),
    make_option("--cache_dir", type = "character"),
    make_option("--outdir", type = "character"),
    make_option("--pc_ldl_or", type = "double", default = 1.25),
    make_option("--pc_lpa_or", type = "double", default = 1.15),
    make_option("--pc_tol", type = "double", default = 0.01),
    make_option("--pc_n", type = "integer", default = 161L, help = "Published instruments, reported not enforced; -1 = no comparison"),
    make_option("--pc_cis_n", type = "integer", default = 205L, help = "Published cis-dense instruments, reported not enforced; -1 = no comparison"),
    make_option("--pc_instruments", type = "character", default = "",
                help = "Published .instruments.tsv, compared with the rebuilt set; reported not enforced")
  )))
}

resolve_file <- function(name, primary, fallback) {
  p <- file.path(primary, name)
  if (file.exists(p)) return(list(path = p, fallback = FALSE))
  if (nzchar(fallback)) {
    q <- file.path(fallback, name)
    if (file.exists(q)) return(list(path = q, fallback = TRUE))
  }
  stop("missing input: ", name, " in ", primary, if (nzchar(fallback)) paste(" or", fallback))
}

main <- function() {
  opt <- parse_opts()
  for (a in c("sumstat_dir", "pheno", "cis_list", "cache_dir", "outdir"))
    if (is.null(opt[[a]])) stop("--", a, " is required")
  if (!opt$arm_type %in% c("primary", "sensitivity")) stop("--arm_type must be primary or sensitivity")
  if (opt$nboot < 1000L) stop("--nboot must be >= 1000: the percentage attenuation is read off this interval")
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(opt$cache_dir, recursive = TRUE, showWarnings = FALSE)
  setDTthreads(min(4L, opt$threads))
  notes <- character(0)
  add_note <- function(...) { s <- sprintf(...); notes <<- c(notes, s); log_msg("NOTE %s", s) }
  verdict_path <- file.path(opt$outdir, "VERDICT.txt")
  # Everything this script writes goes first. Otherwise a rerun that stops at the positive
  # control leaves the last run's estimates next to its INVALID verdict, and stale build
  # artefacts are this project's recurring failure mode.
  unlink(file.path(opt$outdir, c("VERDICT.txt", "positive_control.PASSED", "positive_control.FAILED", "positive_control_rows.csv",
                                 "exposure_correlation.csv", "instrument_sets.tsv", "results_long.csv", "bootstrap_attenuation.csv",
                                 "lpa_instrument_allele_effects.csv", "lpa_attenuation_all_outcomes.csv", "readout_compact.csv", "run_manifest.txt")))

  if (opt$arm_type == "sensitivity") {
    # A sensitivity arm has no control of its own to pass (its exposures differ by
    # construction), so it is only allowed to run after the primary arm's control did.
    if (!nzchar(opt$control_stamp) || !file.exists(opt$control_stamp))
      stop("sensitivity arm needs --control_stamp from a primary arm that passed; none at '", opt$control_stamp, "'")
    if (!startsWith(readLines(opt$control_stamp, n = 1L), "PASSED"))
      stop("--control_stamp does not record a pass: ", opt$control_stamp)
  }

  # ---- inputs ----
  exposure_files <- lapply(EXPOSURES_ALL, function(tr) resolve_file(sprintf("lipids_%s.txt.gz", tr), opt$sumstat_dir, opt$fallback_dir))
  names(exposure_files) <- EXPOSURES_ALL
  outcomes <- trimws(strsplit(opt$outcomes, ",")[[1]])
  if (!"CAD" %in% outcomes) stop("CAD must be among --outcomes: it is the primary endpoint")
  outcome_files <- lapply(outcomes, function(o) resolve_file(sprintf("outcomes_%s.txt.gz", o), opt$sumstat_dir, opt$fallback_dir))
  names(outcome_files) <- outcomes
  src_note <- function(lst) vapply(names(lst), function(n) sprintf("%s=%s%s", n, lst[[n]]$path, if (lst[[n]]$fallback) " [FALLBACK]" else ""), "")
  for (tr in EXPOSURES_ALL) if (exposure_files[[tr]]$fallback)
    add_note("exposure %s is NOT from this arm's directory: it was taken from %s", tr, opt$fallback_dir)
  for (o in outcomes) if (outcome_files[[o]]$fallback)
    add_note("outcome %s is NOT from this arm's directory: it was taken from %s", o, opt$fallback_dir)
  add_note("exposure correlation for conditional F and Q_A comes from %s", opt$pheno)

  # ApoB in merged_statinfix is a symlink into merged/, i.e. not statin-corrected while
  # LDL-C is. Comparing the two models then also compares two statin corrections.
  apob_path <- exposure_files[["APOB"]]$path
  apob_link <- Sys.readlink(apob_path)
  add_note("APOB file %s%s; statin handling of this GWAS follows its directory: merged_statinfix/merged APOB is NOT statin-corrected while LDL-C there is, so the statin-naive arm is the like-for-like comparison of the two models",
           apob_path, if (!is.na(apob_link) && nzchar(apob_link)) sprintf(" (symlink -> %s)", apob_link) else "")
  if (opt$arm_type == "primary" && !nzchar(opt$finngen))
    stop("primary arm needs --finngen: the readout is defined on UK Biobank AND FinnGen")
  ext_inputs <- list()
  if (nzchar(opt$finngen)) {
    if (!file.exists(opt$finngen)) stop("missing --finngen: ", opt$finngen)
    ext_inputs[["FinnGen_R11"]] <- list(path = opt$finngen, format = "finngen")
  }
  if (nzchar(opt$cardiogram)) {
    if (!file.exists(opt$cardiogram)) stop("missing --cardiogram: ", opt$cardiogram)
    ext_inputs[["CARDIoGRAMplusC4D"]] <- list(path = opt$cardiogram, format = "cardiogram")
  } else add_note("CARDIoGRAM not supplied; that cohort is absent from this run")
  if (!file.exists(opt$cis_list)) stop("missing --cis_list: ", opt$cis_list)
  cis_ids <- unique(trimws(readLines(opt$cis_list)))
  cis_ids <- cis_ids[nzchar(cis_ids)]
  if (length(cis_ids) < 2L) stop("--cis_list has ", length(cis_ids), " variants: not a locus-dense instrument")
  log_msg("cis-LPA list: %d variants", length(cis_ids))

  # ---- exposure correlation: Spearman, on all five traits (see mvmr_covariance.R) ----
  # fread(select=) only warns about an absent column, so check the header first.
  ph_cols <- names(fread(opt$pheno, nrows = 0L))
  if (length(setdiff(EXPOSURES_ALL, ph_cols)))
    stop("--pheno lacks column(s): ", paste(setdiff(EXPOSURES_ALL, ph_cols), collapse = ", "))
  ph <- fread(opt$pheno, select = c(if (nzchar(opt$keep)) "IID", EXPOSURES_ALL), showProgress = FALSE)
  if (nzchar(opt$keep)) {
    kp <- fread(opt$keep, header = FALSE, col.names = c("FID", "IID"))
    ph <- ph[as.character(IID) %chin% as.character(kp$IID)]
    if (!nrow(ph)) stop("--keep matches no row of --pheno")
    ph[, IID := NULL]
  }
  Pcor <- cor(as.matrix(ph), use = "pairwise.complete.obs", method = "spearman")
  n_ph <- nrow(ph); rm(ph); invisible(gc())
  if (!is_correlation_matrix(Pcor)) stop("exposure correlation matrix is not valid (NA or degenerate column in --pheno)")
  if (Pcor["LDL_C", "APOB"] < 0.5) add_note("Spearman LDL-C/ApoB = %.3f is implausibly low: check the --pheno columns", Pcor["LDL_C", "APOB"])
  log_msg("exposure Spearman correlation from %d rows:", n_ph)
  print(round(Pcor, 3))
  fwrite(as.data.table(Pcor, keep.rownames = "trait"), file.path(opt$outdir, "exposure_correlation.csv"))

  # ---- pass 1: GWS rows per exposure ----
  log_msg("pass 1: GWS rows from %d exposure files on %d threads", length(exposure_files), opt$threads)
  file_id <- function(p) paste(normalizePath(p), file.size(p), as.numeric(file.mtime(p)))
  gws_paths <- parallel_stop_on_error(exposure_files, function(x) {
    out <- cache_file(opt$cache_dir, "gws", file_id(x$path), "7.3", paste(REGENIE_COLS, collapse = ","))
    extract_table(x$path, out, REGENIE_COLS, thrcol = "LOG10P", thr = 7.3)
  }, opt$threads)
  gws <- lapply(gws_paths, load_regenie)
  for (tr in EXPOSURES_ALL) {
    log_msg("  %-7s %d rows above the pre-filter", tr, nrow(gws[[tr]]))
    if (!nrow(gws[[tr]])) stop("no GWS rows in ", tr, ": a lipid GWAS cannot have none; check the file")
  }

  ivs_full <- lapply(gws, greedy_clump)
  ivs_nowin <- lapply(gws, function(d) greedy_clump(d[!in_lpa_exclusion(CHROM, GENPOS)]))
  for (tr in EXPOSURES_ALL) log_msg("  %-7s %d instruments (%d with the LPA window excluded)", tr, length(ivs_full[[tr]]), length(ivs_nowin[[tr]]))
  rm(gws); invisible(gc())

  # ---- pass 2: every candidate instrument in every file ----
  # Selection needs each candidate's row in every trait, in file order, to reproduce
  # run_mvmr.R's tie order. Two passes keep this at a few thousand rows per file.
  all_ids <- sort(unique(c(unlist(ivs_full), unlist(ivs_nowin), cis_ids, LPA_REPORT)))
  idfile <- file.path(opt$cache_dir, sprintf("ids.%s.txt", opt$arm))
  writeLines(all_ids, idfile)
  ids_md5 <- unname(tools::md5sum(idfile))
  log_msg("pass 2: %d candidate variants, %d UKB files + %d external", length(all_ids), length(exposure_files) + length(outcome_files), length(ext_inputs))
  ukb_jobs <- c(lapply(exposure_files, function(x) list(path = x$path, tag = "lip")),
                setNames(lapply(outcome_files, function(x) list(path = x$path, tag = "out")), paste0("OUT_", outcomes)))
  ukb_paths <- parallel_stop_on_error(ukb_jobs, function(x) {
    out <- cache_file(opt$cache_dir, x$tag, file_id(x$path), ids_md5)
    extract_table(x$path, out, REGENIE_COLS, idfile = idfile, idcol = "ID")
  }, opt$threads)
  rows <- lapply(ukb_paths[EXPOSURES_ALL], load_regenie)
  # The cis list sits on scratch with no script in the repo that builds it, so check what it holds.
  cr <- rows[["LPA"]][ID %in% cis_ids & !duplicated(ID)]
  if (nrow(cr) != length(cis_ids)) stop("--cis_list: ", length(cis_ids) - nrow(cr), " of ", length(cis_ids), " IDs absent from the Lp(a) GWAS")
  if (!all(in_lpa_window(cr$CHROM, cr$GENPOS))) stop("--cis_list: ", sum(!in_lpa_window(cr$CHROM, cr$GENPOS)), " IDs outside the LPA window")
  log_msg("cis-LPA list on the Lp(a) GWAS: %d IDs, all in the window, lowest LOG10P %.1f", nrow(cr), min(cr$LOG10P))
  OUTCOMES <- list()
  for (o in outcomes) {
    d <- load_regenie(ukb_paths[[paste0("OUT_", o)]])
    if (!nrow(d)) stop("no instruments found in UKB outcome ", o)
    OUTCOMES[[paste("UKB", o)]] <- list(cohort = "UKB", outcome = o, type = "ukb", dt = d)
  }
  if (length(ext_inputs)) {
    ext_paths <- parallel_stop_on_error(ext_inputs, function(x) {
      fm <- EXTERNAL_FORMATS[[x$format]]
      out <- cache_file(opt$cache_dir, "ext", file_id(x$path), ids_md5, x$format)
      extract_table(x$path, out, fm$want, fs = fm$fs, idfile = idfile, idcol = fm$idcol, first = fm$first)
    }, opt$threads)
    for (nm in names(ext_inputs)) {
      want <- EXTERNAL_FORMATS[[ext_inputs[[nm]]$format]]$want
      d <- read_tsv(ext_paths[[nm]], want[1:3])
      setnames(d, want, c("ID", "oEA", "oOA", "oBETA", "oSE", "oEAF"))
      d <- d[!is.na(ID) & !is.na(oBETA) & !is.na(oSE)]
      if (!nrow(d)) stop("no instruments found in ", nm, ": check rsID matching in ", ext_inputs[[nm]]$path)
      OUTCOMES[[paste(nm, "CAD")]] <- list(cohort = nm, outcome = "CAD", type = "external", dt = d)
    }
  }
  n_dup <- sum(vapply(rows, function(d) sum(duplicated(d$ID)), integer(1)))
  if (n_dup > 0L) add_note("%d duplicate rsID rows among candidate variants; the allele-matching row is kept", n_dup)

  # Instrument sets per variant and model, frozen before any outcome is looked at.
  SETS <- lapply(setNames(VARIANTS, VARIANTS), function(v)
    lapply(MODELS, function(tr) select_instruments(v, tr, rows, ivs_full, ivs_nowin, cis_ids)))
  fwrite(rbindlist(lapply(VARIANTS, function(v) rbindlist(lapply(names(MODELS), function(m) {
    r <- ref_rows(SETS[[v]][[m]], MODELS[[m]], rows)
    r[, .(variant = v, model = m, ID, CHROM, GENPOS, max_LOG10P = LOG10P, in_lpa_window = in_lpa_window(CHROM, GENPOS))]
  })))), file.path(opt$outdir, "instrument_sets.tsv"), sep = "\t")

  # ---- positive control: nothing with ApoB in it has been fitted yet ----
  ukb_cad <- OUTCOMES[["UKB CAD"]]
  ctrl <- fit_model(build_design(SETS$published$LDL4, MODELS$LDL4, rows, ukb_cad), MODELS$LDL4, Pcor)
  ctrl_ldl <- ctrl[exposure == "LDL_C"]; ctrl_lpa <- ctrl[exposure == "LPA"]
  pc_line <- sprintf("LDL_C OR %.4f (want %.2f), LPA OR %.4f (want %.2f), tolerance %.2f, %d instruments, condF gencov0 LPA %.1f",
                     ctrl_ldl$OR, opt$pc_ldl_or, ctrl_lpa$OR, opt$pc_lpa_or, opt$pc_tol, ctrl_ldl$n_instruments, ctrl_lpa$condF_gencov0)
  log_msg("positive control (LDL4, published instruments, UKB CAD): %s", pc_line)
  if (opt$arm_type == "primary") {
    pass <- !anyNA(c(ctrl_ldl$OR, ctrl_lpa$OR)) &&
      abs(ctrl_ldl$OR - opt$pc_ldl_or) <= opt$pc_tol && abs(ctrl_lpa$OR - opt$pc_lpa_or) <= opt$pc_tol
    lead_ok <- LPA_LEAD %in% SETS$published$LDL4
    if (!pass || !lead_ok) {
      why <- paste(c(if (!pass) "control ORs out of tolerance",
                     if (!lead_ok) sprintf("%s not chosen as the LPA instrument (LOG10P ranking?)", LPA_LEAD)), collapse = "; ")
      txt <- sprintf("VERDICT_LITERAL: INVALID | positive control failed: %s | %s | status=%s", why, pc_line, ctrl_ldl$status)
      writeLines(txt, verdict_path)
      writeLines(c(txt, sprintf("written %s", format(Sys.time()))), file.path(opt$outdir, "positive_control.FAILED"))
      fwrite(ctrl, file.path(opt$outdir, "positive_control_rows.csv"))
      message(txt)
      message("[natA] STOP: no ApoB estimate was computed; the readout forbids reading any until this is fixed.")
      quit(save = "no", status = 2L)
    }
    writeLines(c(paste("PASSED", pc_line), sprintf("written %s", format(Sys.time()))), file.path(opt$outdir, "positive_control.PASSED"))
    if (opt$pc_n >= 0L && ctrl_ldl$n_instruments != opt$pc_n)
      add_note("control ORs pass but the instrument count is %d, published %d: the set is not identical to the published one", ctrl_ldl$n_instruments, opt$pc_n)
    if (nzchar(opt$pc_instruments)) {
      # Not enforced, so an unreadable comparator must say the comparison did not happen
      # rather than leave a reader to assume the set was checked.
      if (file.exists(opt$pc_instruments)) {
        pub <- fread(opt$pc_instruments, select = "ID")$ID
        got <- build_design(SETS$published$LDL4, MODELS$LDL4, rows, ukb_cad)$ids
        add_note("published instrument file %s: %d in common, %d only published, %d only rebuilt", opt$pc_instruments,
                 length(intersect(pub, got)), length(setdiff(pub, got)), length(setdiff(got, pub)))
      } else add_note("--pc_instruments %s does not exist: the rebuilt set was NOT compared with the published one", opt$pc_instruments)
    }
  } else {
    add_note("sensitivity arm: LDL4 control reads %s; not enforced because this arm's exposures differ", pc_line)
  }

  # ---- all models ----
  log_msg("fitting %d models x %d variants x %d outcomes", length(MODELS), length(VARIANTS), length(OUTCOMES))
  fits <- list()
  for (key in names(OUTCOMES)) {
    O <- OUTCOMES[[key]]
    for (v in VARIANTS) for (m in names(MODELS)) {
      D <- build_design(SETS[[v]][[m]], MODELS[[m]], rows, O)
      r <- fit_model(D, MODELS[[m]], Pcor, require_lpa_locus = v != "noLPAwin_all",
                     median = O$outcome == "CAD" && v %in% BOOT_VARIANTS)
      r[, `:=`(arm = opt$arm, cohort = O$cohort, outcome = O$outcome, variant = v, model = m,
               sample_overlap_with_exposures = O$cohort == "UKB")]
      fits[[length(fits) + 1L]] <- r
    }
  }
  long <- rbindlist(fits)
  setcolorder(long, c("arm", "cohort", "outcome", "variant", "model", "exposure"))
  # Attenuation against the same-structure LDL-C analogue on its own rule-built set.
  lpa <- long[exposure == "LPA"]
  ref_lpa <- lpa[model %in% unique(REF_OF), .(cohort, outcome, variant, ref_model = model, b_ref = b, se_ref = se, OR_ref = OR, OR_ref_lo95 = OR_lo95, OR_ref_hi95 = OR_hi95, ref_status = status)]
  comp <- lpa[model %in% names(REF_OF)]
  comp[, ref_model := REF_OF[model]]
  comp <- merge(comp, ref_lpa, by = c("cohort", "outcome", "variant", "ref_model"), all.x = TRUE)
  comp[, atten_pct := fifelse(!is.na(b_ref) & b_ref > 0, 100 * (b_ref - b) / b_ref, NA_real_)]
  comp[, category := mapply(cohort_category, b_ref, b, OR_lo95, OR_hi95)]
  comp[, arm := opt$arm]
  fwrite(long, file.path(opt$outdir, "results_long.csv"))
  fwrite(comp[, .(arm, cohort, outcome, variant, model, ref_model, n_instruments, n_lpa_window, lead_in_model,
                  OR, OR_lo95, OR_hi95, p, OR_ref, OR_ref_lo95, OR_ref_hi95, atten_pct, category, status, ref_status)],
         file.path(opt$outdir, "lpa_attenuation_all_outcomes.csv"))

  # ---- bootstrap on matched instruments ----
  log_msg("bootstrap: %d resamples, %d variants x %d model pairs x %d outcomes", opt$nboot, length(BOOT_VARIANTS), length(REF_OF), length(OUTCOMES))
  boots <- list()
  for (v in BOOT_VARIANTS) for (alt in names(REF_OF)) {
    ref_m <- REF_OF[[alt]]
    traits <- c(MODELS[[ref_m]][1], MODELS[[alt]][1], setdiff(MODELS[[ref_m]], MODELS[[ref_m]][1]))
    ids <- select_instruments(v, traits, rows, ivs_full, ivs_nowin, cis_ids)
    for (key in names(OUTCOMES)) {
      O <- OUTCOMES[[key]]
      D <- build_design(ids, traits, rows, O)
      lpa_ok <- sum(D$in_window) > 0L
      if (!lpa_ok || nrow(D$BX) < length(MODELS[[alt]]) + 2L) next
      b <- boot_attenuation(D, MODELS[[ref_m]], MODELS[[alt]], Pcor, opt$nboot)
      if (is.null(b)) next
      b[, `:=`(arm = opt$arm, cohort = O$cohort, outcome = O$outcome, variant = v, model = alt, ref_model = ref_m)]
      boots[[length(boots) + 1L]] <- b
    }
  }
  boot <- rbindlist(boots)
  if (!nrow(boot)) stop("bootstrap produced no rows: no model pair was estimable in any cohort")
  setcolorder(boot, c("arm", "cohort", "outcome", "variant", "model", "ref_model", "flavour"))
  fwrite(boot, file.path(opt$outdir, "bootstrap_attenuation.csv"))

  # ---- allele effects ----
  # All of Lp(a)'s own instruments, trans ones included (in_lpa_window says which is which):
  # a trans variant that moves ApoB as well is exactly what this table is there to show.
  win_ids <- unique(c(LPA_REPORT, cis_ids, ivs_full[["LPA"]],
                      unlist(lapply(SETS$published, function(s) intersect(s, rows$LPA[in_lpa_window(CHROM, GENPOS), ID])))))
  allele <- lpa_allele_effects(win_ids, rows)
  fwrite(allele, file.path(opt$outdir, "lpa_instrument_allele_effects.csv"))
  lead_row <- if (nrow(allele)) allele[ID == LPA_LEAD] else data.table()
  if (!nrow(lead_row) || anyNA(lead_row[, .(beta_LPA, beta_APOB, beta_LDL_C)])) {
    if (opt$arm_type == "primary") stop(LPA_LEAD, " is missing from LPA, APOB or LDL_C: the per-allele report cannot be produced")
    add_note("%s is missing from at least one of LPA, APOB, LDL_C in this arm", LPA_LEAD)
  }

  # ---- compact readout ----
  bp <- boot[flavour == "strat_param", .(cohort, outcome, variant, model, n_matched, atten_pct_matched = atten_pct_obs,
                                         atten_pct_lo95, atten_pct_hi95, diff_obs, diff_lo95, diff_hi95,
                                         b_ref_matched, se_ref_matched, b_alt_matched, se_alt_matched)]
  compact <- merge(comp, bp, by = c("cohort", "outcome", "variant", "model"), all.x = TRUE)
  # NA, not NOT_READABLE, for variants the bootstrap does not run on: nothing was attempted.
  compact[, `:=`(category_matched = mapply(function(n, r, a, s) if (is.na(n)) NA_character_ else cohort_category(r, a, exp(a - 1.96 * s), exp(a + 1.96 * s)),
                                           n_matched, b_ref_matched, b_alt_matched, se_alt_matched),
                 atten_ci_crosses_threshold = mapply(atten_crosses_threshold, atten_pct_lo95, atten_pct_hi95))]
  compact[, `:=`(median_OR = exp(median_b), median_OR_lo95 = exp(median_b - 1.96 * median_se), median_OR_hi95 = exp(median_b + 1.96 * median_se))]
  compact <- compact[outcome == "CAD"][order(match(cohort, c("UKB", "FinnGen_R11", "CARDIoGRAMplusC4D")), variant, model)]
  keep <- c("arm", "cohort", "outcome", "variant", "model", "ref_model", "n_instruments", "n_lpa_window", "lead_in_model",
            "OR", "OR_lo95", "OR_hi95", "p", "OR_ref", "OR_ref_lo95", "OR_ref_hi95", "atten_pct", "category",
            "n_matched", "atten_pct_matched", "atten_pct_lo95", "atten_pct_hi95", "atten_ci_crosses_threshold", "category_matched",
            "diff_obs", "diff_lo95", "diff_hi95", "condF_phenocov", "condF_gencov0", "Q_A_p_phenocov", "egger_b", "egger_intercept", "egger_intercept_p",
            "median_OR", "median_OR_lo95", "median_OR_hi95", "median_note", "status", "ref_status")
  fwrite(compact[, ..keep], file.path(opt$outdir, "readout_compact.csv"))

  # ---- verdict ----
  fmt_cell <- function(r, cohort) {
    if (!nrow(r)) return(sprintf("%s absent", cohort))
    sprintf("%s Lp(a) OR %.3f (%.3f-%.3f; %d instruments) vs LDL-C model %.3f, attenuation %s [%s]", cohort, r$OR, r$OR_lo95, r$OR_hi95,
            r$n_instruments, r$OR_ref, if (is.na(r$atten_pct)) "NA" else sprintf("%.1f%%", r$atten_pct), r$status)
  }
  catg <- function(r, use_matched = FALSE) {
    x <- if (!nrow(r)) NA_character_ else if (use_matched) r$category_matched else r$category
    if (is.na(x)) "NOT_READABLE" else x
  }
  verdict_for <- function(v, use_matched = FALSE, m = "APOB4") {
    u <- compact[cohort == "UKB" & variant == v & model == m]
    f <- compact[cohort == "FinnGen_R11" & variant == v & model == m]
    list(verdict = combine_categories(c(catg(u, use_matched), catg(f, use_matched))), u = u, f = f)
  }
  lines <- character(0)
  head_tag <- if (opt$arm_type == "primary") "" else "SENSITIVITY_ARM_"
  V <- verdict_for("published")
  lines <- c(lines, sprintf("%sVERDICT_LITERAL: %s | arm=%s, model APOB4 vs LDL4, published instrument rule | %s | %s",
                            head_tag, V$verdict, opt$arm, fmt_cell(V$u, "UKB"), fmt_cell(V$f, "FinnGen")))
  VM <- verdict_for("published", use_matched = TRUE)
  fmt_m <- function(r, cohort) if (!nrow(r) || is.na(r$atten_pct_matched)) sprintf("%s NA", cohort) else
    sprintf("%s %.1f%% (95%% CI %s; %d matched instruments)", cohort, r$atten_pct_matched,
            if (is.na(r$atten_pct_lo95)) "NA, denominator unstable" else sprintf("%.1f to %.1f", r$atten_pct_lo95, r$atten_pct_hi95), r$n_matched)
  lines <- c(lines, sprintf("%sVERDICT_MATCHED: %s | attenuation on instruments matched between the two models | %s | %s",
                            head_tag, VM$verdict, fmt_m(VM$u, "UKB"), fmt_m(VM$f, "FinnGen")))
  VC <- verdict_for("cis_dense")
  lines <- c(lines, sprintf("%sVERDICT_CIS_DENSE: %s | %s | %s", head_tag, VC$verdict, fmt_cell(VC$u, "UKB"), fmt_cell(VC$f, "FinnGen")))
  # The paper does not use one rule in both cohorts: UK Biobank rests on rs10455872 alone,
  # FinnGen on the locus-dense set. The readout's "same instrument rules" most naturally
  # means the line above, which applies the distance rule to both, so this pairing does not
  # replace it; it is printed because it is the comparison a reader of the paper will make.
  VP <- combine_categories(c(catg(V$u), catg(VC$f)))
  lines <- c(lines, sprintf("%sVERDICT_PAPER_PAIRING: %s | UKB on the published rule, FinnGen on cis-dense, as in the paper | %s | %s",
                            head_tag, VP, fmt_cell(V$u, "UKB"), fmt_cell(VC$f, "FinnGen")))
  if (length(unique(c(V$verdict, VM$verdict, VC$verdict, VP))) > 1L)
    lines <- c(lines, "VERDICTS_DISAGREE: the literal read depends on the instrument set; report all three and do not cite the most favourable.")
  crosses <- c(if (nrow(VM$u)) VM$u$atten_ci_crosses_threshold, if (nrow(VM$f)) VM$f$atten_ci_crosses_threshold)
  if (any(crosses %in% TRUE)) lines <- c(lines, "ATTENUATION_CI_CROSSES_A_THRESHOLD: the point estimate falls in one readout band but the interval spans a 25% or 50% boundary.")
  # The readout's sensitivity list, read by the same rule. Never used to replace the primary.
  for (v in setdiff(VARIANTS, c("noLPAwin_all"))) for (m in names(REF_OF)) {
    if (v == "published" && m == "APOB4") next
    if (v == "cis_dense" && m == "APOB4") next
    S <- verdict_for(v, m = m)
    lines <- c(lines, sprintf("SENSITIVITY %s %s: %s | %s | %s", v, m, S$verdict, fmt_cell(S$u, "UKB"), fmt_cell(S$f, "FinnGen")))
  }
  S <- verdict_for("noLPAwin_all")
  lines <- c(lines, sprintf("TRANS_ONLY noLPAwin_all APOB4, a different question and not part of the read: %s | %s | %s", S$verdict, fmt_cell(S$u, "UKB"), fmt_cell(S$f, "FinnGen")))
  pm <- long[cohort == "UKB" & outcome == "CAD" & variant == "published" & model == "APOB4"]
  lines <- c(lines, sprintf("PRIMARY_MODEL_UKB_CAD: %s", paste(sprintf("%s OR %.3f (%.3f-%.3f) condF %.1f", pm$exposure, pm$OR, pm$OR_lo95, pm$OR_hi95, pm$condF_phenocov), collapse = "; ")))
  if (nrow(lead_row)) lines <- c(lines, sprintf("LEAD_%s: per Lp(a)-raising allele %s: Lp(a) %+.3f SD, ApoB %+.3f SD (ApoB/Lp(a) %.3f), LDL-C %+.3f SD (LDL-C/Lp(a) %.3f)",
                                                LPA_LEAD, lead_row$raising_allele, lead_row$beta_LPA, lead_row$beta_APOB, lead_row$apob_sd_per_lpa_sd, lead_row$beta_LDL_C, lead_row$ldl_sd_per_lpa_sd))
  if (nrow(V$u)) lines <- c(lines, sprintf("STRENGTH: UKB APOB4 published Lp(a) condF phenocov %.1f, Q_A P %.3g, Egger intercept P %s",
                                           V$u$condF_phenocov, V$u$Q_A_p_phenocov, format(V$u$egger_intercept_p, digits = 3)))
  ref_cis <- ref_lpa[cohort == "UKB" & outcome == "CAD" & variant == "cis_dense" & ref_model == "LDL4"]
  if (opt$arm_type == "primary" && nrow(ref_cis)) {
    n_cis <- long[cohort == "UKB" & outcome == "CAD" & variant == "cis_dense" & model == "LDL4" & exposure == "LPA", n_instruments]
    lines <- c(lines, sprintf("CHECK_NOT_ENFORCED: cis-dense LDL4 UKB CAD Lp(a) OR %.3f on %d instruments (published 1.149 on %d)", ref_cis$OR_ref, n_cis, opt$pc_cis_n))
    ref_fin <- ref_lpa[cohort == "FinnGen_R11" & variant == "cis_dense" & ref_model == "LDL4"]
    if (nrow(ref_fin)) lines <- c(lines, sprintf("CHECK_NOT_ENFORCED: cis-dense LDL4 FinnGen Lp(a) OR %.3f (published 1.18 on 185 instruments from PLINK r2 clumping; this script's polygenic rule is distance-only, so the set differs)", ref_fin$OR_ref))
  }
  if (length(notes)) lines <- c(lines, paste("NOTE:", notes))
  writeLines(lines, verdict_path)
  writeLines(c(sprintf("arm=%s type=%s", opt$arm, opt$arm_type), sprintf("MVMR %s, MendelianRandomization %s, data.table %s, R %s",
                       packageVersion("MVMR"), packageVersion("MendelianRandomization"), packageVersion("data.table"), getRversion()),
               sprintf("nboot=%d seed=42 Spearman rows=%d", opt$nboot, n_ph), src_note(exposure_files), src_note(outcome_files),
               vapply(names(ext_inputs), function(n) sprintf("%s=%s", n, ext_inputs[[n]]$path), ""), sprintf("cis_list=%s (%d)", opt$cis_list, length(cis_ids))),
             file.path(opt$outdir, "run_manifest.txt"))
  message(paste(lines, collapse = "\n"))
  log_msg("done; outputs in %s", opt$outdir)
  invisible(NULL)
}

if (sys.nframe() == 0L) main()
