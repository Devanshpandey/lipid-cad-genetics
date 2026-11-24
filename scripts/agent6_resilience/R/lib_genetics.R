# ============================================================================
# lib_genetics.R - pure functions for stage 09 (genetics of the residual
# metabolome score, preregistered H4c): REGENIE output QC, lambda, distance-
# based loci, allele harmonisation, two-sample MR (IVW, weighted median,
# MR-Egger), LDSC and munge log parsing. Unit-tested in
# tests/testthat/test_genetics.R. Source lib_resilience.R first (check()).
# Kept separate from lib_resilience.R so that editing it does not invalidate
# the hash-gated manifests of stages 00-07. Since 2026-09-29 the pipeline uses
# only the QC, loci and LDSC parsing functions; the harmonisation and MR
# functions are kept (tested, and equal to MendelianRandomization) for reuse.
# ============================================================================

REGENIE_COLS <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST",
                  "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

# ---------------------------------------------------------------------------
# REGENIE output
# ---------------------------------------------------------------------------
read_regenie_files <- function(files, chroms) {
  check(length(files) == length(chroms), "one REGENIE file per chromosome")
  parts <- lapply(seq_along(files), function(i) {
    f <- files[i]
    check(file.exists(f), "REGENIE output present: %s", basename(f))
    hdr <- names(fread(f, nrows = 0))
    check(identical(hdr, REGENIE_COLS), "REGENIE header as expected in %s", basename(f))
    d <- fread(f, select = REGENIE_COLS[1:12],
               colClasses = list(character = c("ID", "ALLELE0", "ALLELE1", "TEST")))
    check(nrow(d) > 0, "%s has rows", basename(f))
    check(all(d$CHROM == chroms[i]), "%s holds only chromosome %s", basename(f), chroms[i])
    d
  })
  rbindlist(parts)
}

# -log10 P(chi2_1 > x) without underflow
log10p_from_chisq <- function(chisq) -pchisq(chisq, 1, lower.tail = FALSE, log.p = TRUE) / log(10)

# -log10 two-sided normal p from z, without underflow
log10p_from_z <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

qc_regenie <- function(d, n_expected, min_call_rate, min_maf) {
  counts <- c(input = nrow(d))
  d <- d[TEST == "ADD"];                                              counts["after_test_add"] <- nrow(d)
  d <- d[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(LOG10P) &
         is.finite(A1FREQ) & is.finite(N) & is.finite(CHISQ)];          counts["after_finite"] <- nrow(d)
  d <- d[N >= min_call_rate * n_expected];                            counts["after_call_rate"] <- nrow(d)
  d <- d[pmin(A1FREQ, 1 - A1FREQ) >= min_maf];                        counts["after_maf"] <- nrow(d)
  dup <- unique(d$ID[duplicated(d$ID)])
  d <- d[!ID %chin% dup];                                             counts["after_dedup"] <- nrow(d)
  list(d = d, counts = counts)
}

lambda_gc <- function(chisq) median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

# Greedy distance-based loci: the most significant remaining variant becomes a
# lead; any variant within window_bp of an existing lead on the same chromosome
# belongs to that lead.
define_loci <- function(d, thr_log10p, window_bp) {
  sig <- d[LOG10P >= thr_log10p][order(-LOG10P)]
  if (!nrow(sig)) return(sig[0][, n_sig_in_window := integer(0)][])
  is_lead <- logical(nrow(sig)); lc <- integer(0); lp <- numeric(0)
  for (i in seq_len(nrow(sig))) {
    if (!any(lc == sig$CHROM[i] & abs(lp - sig$GENPOS[i]) <= window_bp)) {
      is_lead[i] <- TRUE; lc <- c(lc, sig$CHROM[i]); lp <- c(lp, sig$GENPOS[i])
    }
  }
  leads <- sig[is_lead]
  leads[, n_sig_in_window := vapply(seq_len(.N), function(j)
    sum(sig$CHROM == CHROM[j] & abs(sig$GENPOS - GENPOS[j]) <= window_bp), integer(1))]
  leads[]
}

# ---------------------------------------------------------------------------
# Allele harmonisation. exp and out: data.tables with SNP, A1 (effect), A2,
# BETA, SE plus optional extra columns. Returns out BETA expressed per copy of
# the exposure A1 allele, a 'flip' flag (TRUE when out alleles were swapped,
# so any out-side allele frequency must become 1 - f), and counts per rule.
# Palindromic exposure SNPs and non-SNVs are dropped.
# ---------------------------------------------------------------------------
COMPLEMENT <- c(A = "T", C = "G", G = "C", T = "A")
is_snv <- function(a1, a2) nchar(a1) == 1L & nchar(a2) == 1L & a1 %chin% names(COMPLEMENT) &
  a2 %chin% names(COMPLEMENT) & a1 != a2
is_palindromic <- function(a1, a2) (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
  (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")

harmonize_alleles <- function(exp, out) {
  e <- copy(exp); o <- copy(out)
  e[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  o[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  m <- merge(e, o, by = "SNP", suffixes = c(".exp", ".out"))
  counts <- c(exposure = nrow(e), outcome_matched = nrow(m))
  m <- m[is_snv(A1.exp, A2.exp) & is_snv(A1.out, A2.out)];  counts["snv"] <- nrow(m)
  m <- m[!is_palindromic(A1.exp, A2.exp)];                  counts["non_palindromic"] <- nrow(m)
  cA1 <- unname(COMPLEMENT[m$A1.out]); cA2 <- unname(COMPLEMENT[m$A2.out])
  direct    <- m$A1.exp == m$A1.out & m$A2.exp == m$A2.out
  swapped   <- m$A1.exp == m$A2.out & m$A2.exp == m$A1.out
  strand    <- !direct & !swapped & m$A1.exp == cA1 & m$A2.exp == cA2
  strand_sw <- !direct & !swapped & !strand & m$A1.exp == cA2 & m$A2.exp == cA1
  m[, status := fifelse(direct, "direct", fifelse(swapped, "swapped",
                fifelse(strand, "strand_flip", fifelse(strand_sw, "strand_flip_swapped", "mismatch"))))]
  m[, flip := status %chin% c("swapped", "strand_flip_swapped")]
  m[flip == TRUE, BETA.out := -BETA.out]
  for (s in c("direct", "swapped", "strand_flip", "strand_flip_swapped", "mismatch")) counts[s] <- sum(m$status == s)
  m <- m[status != "mismatch"];                             counts["harmonised"] <- nrow(m)
  list(dat = m, counts = counts)
}

# ---------------------------------------------------------------------------
# Two-sample MR. Formulas follow MendelianRandomization defaults (simple
# first-order weights, normal-distribution p values) so the stage can
# cross-check them against that package at run time.
# ---------------------------------------------------------------------------
mr_ivw_fn <- function(bx, by, sy) {
  k <- length(bx); w <- sy^-2
  est <- sum(w * bx * by) / sum(w * bx^2)
  se_fixed <- 1 / sqrt(sum(w * bx^2))
  Q <- sum(w * (by - est * bx)^2)
  sigma <- sqrt(Q / (k - 1))
  model <- if (k > 3) "random" else "fixed"
  se <- if (model == "random") se_fixed * max(1, sigma) else se_fixed
  list(method = "IVW", k = k, estimate = est, se = se, p = 2 * pnorm(-abs(est / se)),
       Q = Q, Q_df = k - 1, Q_p = pchisq(Q, k - 1, lower.tail = FALSE), model = model)
}

i2_gx <- function(bx_abs, sx) {
  w <- 1 / sx^2; mu <- sum(bx_abs * w) / sum(w)
  Q <- sum(w * (bx_abs - mu)^2); k <- length(bx_abs)
  max(0, (Q - (k - 1)) / Q)
}

mr_egger_fn <- function(bx, by, sy, sx = NULL) {
  k <- length(bx); s <- sign(bx); bxo <- abs(bx); byo <- by * s
  fit <- lm(byo ~ bxo, weights = sy^-2); sm <- summary(fit); sig <- sm$sigma
  est <- unname(coef(fit)[2]); se <- sm$coefficients[2, 2] / min(sig, 1)
  int <- unname(coef(fit)[1]); int_se <- sm$coefficients[1, 2] / min(sig, 1)
  list(method = "MR-Egger", k = k, estimate = est, se = se, p = 2 * pnorm(-abs(est / se)),
       intercept = int, intercept_se = int_se, intercept_p = 2 * pnorm(-abs(int / int_se)),
       # Bowden 2016 (IJE): with inverse-variance-weighted Egger, I2_GX is computed on
       # bx/sy with standard errors sx/sy (matches MendelianRandomization)
       I2GX = if (is.null(sx)) NA_real_ else i2_gx(bxo / sy, sx / sy))
}

weighted_median <- function(b, w) {
  o <- order(b); b <- b[o]; w <- w[o]
  ws <- (cumsum(w) - 0.5 * w) / sum(w)
  below <- max(which(ws < 0.5))
  b[below] + (b[below + 1] - b[below]) * (0.5 - ws[below]) / (ws[below + 1] - ws[below])
}

mr_wm_fn <- function(bx, by, sx, sy, n_boot = 10000, seed = 42) {
  check(length(bx) >= 3, "weighted median needs at least 3 instruments")
  ratio <- by / bx; w <- (sy / abs(bx))^-2
  est <- weighted_median(ratio, w)
  set.seed(seed)
  boot <- vapply(seq_len(n_boot), function(i) {
    weighted_median(rnorm(length(by), by, sy) / rnorm(length(bx), bx, sx), w)
  }, numeric(1))
  se <- sd(boot)
  list(method = "Weighted median", k = length(bx), estimate = est, se = se, p = 2 * pnorm(-abs(est / se)))
}

# ---------------------------------------------------------------------------
# LDSC / munge log parsing. Labels are located with startsWith() so digits in
# the labels themselves ("h2", "Chi^2") can never be read as values.
# ---------------------------------------------------------------------------
.num_rx <- "([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
.line_rest <- function(lines, label, which = 1L) {
  hit <- which(startsWith(lines, label))
  if (length(hit) < which) return(NA_character_)
  substring(lines[hit[which]], nchar(label) + 1L)
}
.pair <- function(rest) {
  if (is.na(rest)) return(c(NA_real_, NA_real_))
  m <- regmatches(rest, regexec(paste0("^\\s*", .num_rx, "\\s*\\(\\s*", .num_rx, "\\s*\\)"), rest, perl = TRUE))[[1]]
  if (length(m) < 3) c(NA_real_, NA_real_) else as.numeric(m[2:3])
}
.single <- function(rest) {
  if (is.na(rest)) return(NA_real_)
  m <- regmatches(rest, regexec(paste0("^\\s*", .num_rx, "\\s*$"), rest, perl = TRUE))[[1]]
  if (length(m) < 2) NA_real_ else as.numeric(m[2])
}
.grab_int <- function(lines, rx) {
  m <- regmatches(lines, regexec(rx, lines))
  v <- unlist(lapply(m, function(x) if (length(x) >= 2) x[2] else NULL))
  if (length(v)) as.integer(v[1]) else NA_integer_
}

parse_ldsc_h2_log <- function(lines) {
  h2 <- .pair(.line_rest(lines, "Total Observed scale h2:"))
  ic <- .pair(.line_rest(lines, "Intercept:"))
  ra <- .pair(.line_rest(lines, "Ratio:"))
  list(h2 = h2[1], h2_se = h2[2],
       lambda_gc = .single(.line_rest(lines, "Lambda GC:")),
       mean_chi2 = .single(.line_rest(lines, "Mean Chi^2:")),
       intercept = ic[1], intercept_se = ic[2], ratio = ra[1], ratio_se = ra[2],
       n_snps = .grab_int(lines, "After merging with regression SNP LD, ([0-9]+) SNPs remain"))
}

parse_ldsc_rg_log <- function(lines) {
  i <- which(startsWith(lines, "Summary of Genetic Correlation Results"))
  check(length(i) == 1, "LDSC rg log has exactly one summary table")
  hdr <- strsplit(trimws(lines[i + 1]), "\\s+")[[1]]
  val <- strsplit(trimws(lines[i + 2]), "\\s+")[[1]]
  check(length(hdr) == length(val) && all(c("rg", "se", "z", "p", "gcov_int", "gcov_int_se") %in% hdr),
        "LDSC rg summary table is well formed")
  s <- setNames(as.list(val), hdr)
  num <- function(k) suppressWarnings(as.numeric(s[[k]]))
  h1 <- .pair(.line_rest(lines, "Total Observed scale h2:", 1L))
  h2 <- .pair(.line_rest(lines, "Total Observed scale h2:", 2L))
  gc <- .pair(.line_rest(lines, "Total Observed scale gencov:"))
  list(rg = num("rg"), rg_se = num("se"), z = num("z"), p = num("p"),
       gcov_int = num("gcov_int"), gcov_int_se = num("gcov_int_se"),
       h2_p1 = h1[1], h2_p1_se = h1[2], h2_p2 = h2[1], h2_p2_se = h2[2],
       gencov = gc[1], gencov_se = gc[2],
       n_snps = .grab_int(lines, "([0-9]+) SNPs with valid alleles"))
}

parse_munge_log <- function(lines) .grab_int(lines, "Writing summary statistics for ([0-9]+) SNPs")

parse_md5_file <- function(path) {
  x <- readLines(path); x <- x[nzchar(x)]
  data.table(md5 = sub("\\s.*$", "", x), path = sub("^[0-9a-fA-F]+\\s+\\*?", "", x))
}
