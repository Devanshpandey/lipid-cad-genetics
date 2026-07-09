#!/usr/bin/env Rscript
# ============================================================
# Multivariable Mendelian Randomization (MVMR)
# All 7 lipid traits as joint exposures → one CAD/MACE outcome
#
# Scientific question:
#   Which lipid fractions independently predispose to CAD
#   after accounting for the others?
#
# Expected findings from literature:
#   ApoB  → CAD  independent of LDL-C (particle number > cholesterol content)
#   LDL-C → CAD  attenuates when ApoB is in model
#   Lp(a) → CAD  independent of LDL-C and ApoB
#   HDL-C → null after adjusting for TG (HDL hypothesis largely refuted)
#   TG    → attenuates after LDL-C/ApoB adjustment
#
# Methods:
#   MVMR-IVW  (MendelianRandomization::mr_mvivw)
#   Conditional F-statistics per exposure (MVMR package)
#   Q_A pleiotropy test (MVMR package)
#   Sensitivity: MVMR-Egger, MVMR-Lasso (if available)
#
# Usage:
#   Rscript run_mvmr.R \
#     --outcome     CAD \
#     --sumstat_dir /path/to/merged \
#     --out         /path/to/mvmr_CAD.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(optparse)
})
set.seed(42)

# Optional MVMR package for conditional F-stats and Q_A
has_mvmr <- requireNamespace("MVMR", quietly = TRUE)
if (has_mvmr) {
  library(MVMR)
  message("[MVMR] MVMR package available — conditional F-stats and Q_A enabled")
} else {
  message("[MVMR] MVMR package not found — install with install.packages('MVMR')")
  message("[MVMR] Proceeding with MendelianRandomization::mr_mvivw only")
}

opt_list <- list(
  make_option("--outcome",      type = "character", help = "Outcome trait name (e.g. CAD)"),
  make_option("--sumstat_dir",  type = "character", help = "Dir with merged REGENIE sumstats"),
  make_option("--p_thresh",     type = "double",  default = 5e-8,   help = "GWS threshold"),
  make_option("--r2_clump",     type = "double",  default = 0.001,  help = "LD r2 for clumping"),
  make_option("--kb_clump",     type = "integer", default = 10000L, help = "Window kb for clumping"),
  make_option("--min_ivs",      type = "integer", default = 10L,    help = "Min instruments for MVMR"),
  make_option("--out",          type = "character", help = "Output CSV path")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Based on LDSC rg results (step 6):
# LDL_C ↔ ApoB rg = 1.074 (near-identical — collinear in MVMR)
# HDL_C ↔ ApoA1 rg = 1.013 (near-identical — collinear in MVMR)
# TOT_CHOL is a linear combination of LDL_C + HDL_C + TG/5 — redundant
# MVMR exposure set: {LDL_C, HDL_C, TRIGLY, LPA} — 4 independent traits
# LPA included but expect wide CIs (single-locus instrument set)
LIPID_TRAITS <- c("LDL_C", "HDL_C", "TRIGLY", "LPA")

message(sprintf("[MVMR] Outcome: %s", opt$outcome))
message(sprintf("[MVMR] Joint exposures: %s", paste(LIPID_TRAITS, collapse = ", ")))

# ── 1. Load all lipid sumstats ──────────────────────────────
message("[MVMR] Loading lipid summary statistics...")
lipid_data <- list()
for (trait in LIPID_TRAITS) {
  f <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", trait))
  if (!file.exists(f)) {
    message(sprintf("  WARNING: %s not found — excluding from MVMR", f))
    next
  }
  dt <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA")
  dt[, pval  := 10^(-LOG10P)]
  dt[, trait := trait]
  lipid_data[[trait]] <- dt
  message(sprintf("  Loaded %s: %d variants", trait, nrow(dt)))
}

active_traits <- names(lipid_data)
if (length(active_traits) < 2) stop("[MVMR] Need at least 2 exposure traits.")

# ── 2. Select instruments: union of GWS SNPs across all lipids ──
message("[MVMR] Selecting instruments (GWS SNPs across all lipid traits)...")

# Greedy LD-distance clumping per trait, then take union
greedy_clump <- function(dt, p_thresh, kb_window) {
  gws <- dt[pval < p_thresh & !is.na(BETA) & !is.na(SE)]
  if (nrow(gws) == 0) return(character(0))
  gws <- gws[order(pval)]
  keep_snps <- character(0)
  keep_pos  <- integer(0)
  keep_chr  <- character(0)
  for (i in seq_len(nrow(gws))) {
    chr <- as.character(gws$CHROM[i])
    pos <- as.integer(gws$GENPOS[i])
    same_chr <- keep_chr == chr
    if (length(keep_pos) == 0 ||
        !any(same_chr & abs(keep_pos[same_chr] - pos) < kb_window * 1000L)) {
      keep_snps <- c(keep_snps, gws$ID[i])
      keep_pos  <- c(keep_pos,  pos)
      keep_chr  <- c(keep_chr,  chr)
    }
  }
  keep_snps
}

all_instruments <- character(0)
instruments_per_trait <- list()
for (trait in active_traits) {
  ivs <- greedy_clump(lipid_data[[trait]], opt$p_thresh, opt$kb_clump)
  instruments_per_trait[[trait]] <- ivs
  message(sprintf("  %s: %d independent GWS instruments", trait, length(ivs)))
  all_instruments <- union(all_instruments, ivs)
}

# Second-pass clump on the union (remove duplicates within 1Mb across traits)
# Build a reference row per SNP using the trait with lowest p-value
snp_ref <- rbindlist(lapply(active_traits, function(tr) {
  lipid_data[[tr]][ID %in% all_instruments,
                   .(ID, CHROM, GENPOS, pval)]
}))
snp_ref <- snp_ref[snp_ref[, .I[which.min(pval)], by = ID]$V1]
snp_ref <- snp_ref[order(pval)]

final_instruments <- character(0)
final_pos  <- integer(0)
final_chr  <- character(0)
for (i in seq_len(nrow(snp_ref))) {
  chr <- as.character(snp_ref$CHROM[i])
  pos <- as.integer(snp_ref$GENPOS[i])
  same_chr <- final_chr == chr
  if (length(final_pos) == 0 ||
      !any(same_chr & abs(final_pos[same_chr] - pos) < opt$kb_clump * 1000L)) {
    final_instruments <- c(final_instruments, snp_ref$ID[i])
    final_pos <- c(final_pos, pos)
    final_chr <- c(final_chr, chr)
  }
}

n_ivs <- length(final_instruments)
message(sprintf("[MVMR] Final instrument set: %d independent SNPs", n_ivs))

if (n_ivs < opt$min_ivs) {
  message(sprintf("[MVMR] Too few instruments (%d < %d) — skipping", n_ivs, opt$min_ivs))
  quit(status = 0)
}

# ── 3. Build beta/SE matrix for all exposures at instrument SNPs ──
message("[MVMR] Building exposure beta matrix...")

# Reference alleles from the first available trait
ref_alleles <- lipid_data[[active_traits[1]]][ID %in% final_instruments,
                                               .(ID, ALLELE1, ALLELE0, CHROM, GENPOS)]

# For each trait, extract beta and SE at instrument SNPs
beta_mat <- matrix(NA_real_, nrow = n_ivs, ncol = length(active_traits),
                   dimnames = list(final_instruments, active_traits))
se_mat   <- matrix(NA_real_, nrow = n_ivs, ncol = length(active_traits),
                   dimnames = list(final_instruments, active_traits))

for (trait in active_traits) {
  dt <- lipid_data[[trait]][ID %in% final_instruments]
  # Align alleles to reference: flip beta if effect allele differs
  dt <- merge(dt, ref_alleles[, .(ID, ref_A1 = ALLELE1, ref_A0 = ALLELE0)],
              by = "ID", all.x = TRUE)
  dt[, beta_aligned := fifelse(
    ALLELE1 == ref_A1, BETA,
    fifelse(ALLELE1 == ref_A0, -BETA, NA_real_)
  )]
  # Fill matrix
  idx <- match(dt$ID, final_instruments)
  beta_mat[idx, trait] <- dt$beta_aligned
  se_mat[idx,   trait] <- dt$SE
}

# Drop SNPs missing data in any exposure
complete_rows <- complete.cases(beta_mat) & complete.cases(se_mat)
beta_mat <- beta_mat[complete_rows, , drop = FALSE]
se_mat   <- se_mat[complete_rows,   , drop = FALSE]
n_complete <- nrow(beta_mat)
message(sprintf("[MVMR] %d SNPs with complete data across all exposures", n_complete))

if (n_complete < opt$min_ivs) {
  message(sprintf("[MVMR] Too few complete SNPs (%d) — skipping", n_complete))
  quit(status = 0)
}

# ── 4. Load outcome data ─────────────────────────────────────
message(sprintf("[MVMR] Loading outcome: %s", opt$outcome))
out_file <- file.path(opt$sumstat_dir, sprintf("outcomes_%s.txt.gz", opt$outcome))
if (!file.exists(out_file)) stop("Outcome file missing: ", out_file)

out_dt <- fread(cmd = paste("gunzip -c", shQuote(out_file)), na.strings = "NA")
out_dt[, pval := 10^(-LOG10P)]
out_sub <- out_dt[ID %in% rownames(beta_mat)]

if (nrow(out_sub) == 0) stop("[MVMR] No instrument SNPs found in outcome data")

# Align outcome betas to reference alleles
out_sub <- merge(out_sub, ref_alleles[ID %in% rownames(beta_mat),
                                       .(ID, ref_A1 = ALLELE1, ref_A0 = ALLELE0)],
                 by = "ID", all.x = TRUE)
out_sub[, beta_aligned := fifelse(
  ALLELE1 == ref_A1, BETA,
  fifelse(ALLELE1 == ref_A0, -BETA, NA_real_)
)]

# Order to match beta_mat rows
out_sub <- out_sub[match(rownames(beta_mat), ID)]
valid   <- !is.na(out_sub$beta_aligned) & !is.na(out_sub$SE)
beta_mat <- beta_mat[valid, , drop = FALSE]
se_mat   <- se_mat[valid,   , drop = FALSE]
beta_out <- out_sub$beta_aligned[valid]
se_out   <- out_sub$SE[valid]

n_final <- nrow(beta_mat)
message(sprintf("[MVMR] Final: %d SNPs with complete exposure + outcome data", n_final))
if (n_final < opt$min_ivs) {
  message(sprintf("[MVMR] Too few final SNPs (%d) — skipping", n_final))
  quit(status = 0)
}

# ── 5. Run MVMR-IVW (MendelianRandomization package) ────────
message("[MVMR] Running MVMR-IVW...")
mr_input <- mr_mvinput(
  bx      = beta_mat,
  bxse    = se_mat,
  by      = beta_out,
  byse    = se_out,
  exposure = active_traits,
  outcome  = opt$outcome
)

mvivw_res <- mr_mvivw(mr_input)

# MVMR-Egger sensitivity (tests directional pleiotropy)
mvegger_res <- tryCatch(
  mr_mvegger(mr_input),
  error = function(e) { message("  MVMR-Egger failed: ", e$message); NULL }
)

# ── 6. MVMR package: conditional F-stats + Q_A ──────────────
cond_f   <- rep(NA_real_, length(active_traits))
names(cond_f) <- active_traits
Q_A_stat <- NA_real_
Q_A_pval <- NA_real_

if (has_mvmr) {
  message("[MVMR] Computing conditional F-statistics (MVMR package)...")
  tryCatch({
    mvmr_input <- format_mvmr(
      BXGs  = beta_mat,
      BYG   = beta_out,
      seBXGs = se_mat,
      seBYG  = se_out,
      RSID   = rownames(beta_mat)
    )
    fstats <- strength_mvmr(mvmr_input, gencov = 0)  # gencov=0 = no genetic covariance
    cond_f[active_traits] <- fstats$exposure[match(active_traits, rownames(fstats))]

    # Q_A: test for pleiotropy
    qa_res <- pleiotropy_mvmr(mvmr_input)
    Q_A_stat <- qa_res$Qstat
    Q_A_pval <- qa_res$Qpval
    message(sprintf("  Q_A pleiotropy: stat=%.2f, p=%.3f", Q_A_stat, Q_A_pval))
  }, error = function(e) {
    message("  MVMR conditional F-stats failed: ", e$message)
  })
}

# ── 7. Compile results ───────────────────────────────────────
message("[MVMR] Compiling results...")

n_exp <- length(active_traits)
results <- data.table(
  outcome          = opt$outcome,
  exposure         = active_traits,
  n_instruments    = n_final,
  # MVMR-IVW
  mvivw_b          = mvivw_res@Estimate,
  mvivw_se         = mvivw_res@StdError,
  mvivw_pval       = mvivw_res@Pvalue,
  mvivw_OR         = exp(mvivw_res@Estimate),
  mvivw_OR_lo95    = exp(mvivw_res@Estimate - 1.96 * mvivw_res@StdError),
  mvivw_OR_hi95    = exp(mvivw_res@Estimate + 1.96 * mvivw_res@StdError),
  # MVMR-Egger
  mvegger_b        = if (!is.null(mvegger_res)) mvegger_res@Estimate        else rep(NA_real_, n_exp),
  mvegger_se       = if (!is.null(mvegger_res)) mvegger_res@StdError.Est    else rep(NA_real_, n_exp),
  mvegger_pval     = if (!is.null(mvegger_res)) mvegger_res@Pvalue.Est      else rep(NA_real_, n_exp),
  mvegger_intercept   = if (!is.null(mvegger_res)) rep(mvegger_res@Intercept,   n_exp) else rep(NA_real_, n_exp),
  mvegger_intercept_p = if (!is.null(mvegger_res)) rep(mvegger_res@Pvalue.Int,  n_exp) else rep(NA_real_, n_exp),
  # Instrument strength
  cond_F           = cond_f[active_traits],
  weak_instrument  = cond_f[active_traits] < 10,
  # Pleiotropy
  Q_A_stat         = Q_A_stat,
  Q_A_pval         = Q_A_pval,
  pleiotropy_flag  = !is.na(Q_A_pval) & Q_A_pval < 0.05
)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(results, opt$out, sep = ",")
message(sprintf("[MVMR] Results written: %s", opt$out))

# ── 8. Print interpretation ──────────────────────────────────
message(sprintf("\n[MVMR] === %s — MVMR-IVW results (%d SNPs) ===", opt$outcome, n_final))
for (i in seq_len(nrow(results))) {
  sig <- if (!is.na(results$mvivw_pval[i]) && results$mvivw_pval[i] < 0.05) "*" else ""
  ind <- if (!is.na(results$mvivw_pval[i]) && results$mvivw_pval[i] < 0.05/n_exp) "**BONF**" else ""
  message(sprintf("  %-10s  OR=%.3f [%.3f-%.3f]  p=%.2e  condF=%.1f  %s%s",
    results$exposure[i],
    results$mvivw_OR[i], results$mvivw_OR_lo95[i], results$mvivw_OR_hi95[i],
    results$mvivw_pval[i],
    results$cond_F[i],
    sig, ind))
}
if (!is.na(Q_A_pval) && Q_A_pval < 0.05) {
  message(sprintf("  *** Q_A pleiotropy flagged: p=%.3f ***", Q_A_pval))
}
message("\n[MVMR] Interpretation guide:")
message("  OR > 1 = exposure independently INCREASES outcome risk")
message("  OR < 1 = exposure independently DECREASES outcome risk")
message("  condF < 10 = weak instrument warning for that exposure")
message("  Bonferroni threshold: p < ", round(0.05/n_exp, 4), " (", n_exp, " exposures)")
message("[MVMR] Complete.")
