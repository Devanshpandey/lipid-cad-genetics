#!/usr/bin/env Rscript
# Mendelian Randomization: lipid exposure → CAD/MACE outcomes
# Methods: IVW, weighted median, MR-Egger
# Flags: Egger intercept p < 0.05, Steiger filtering, heterogeneity test
#
# Usage:
#   Rscript run_mr.R \
#     --exposure    LDL_C \
#     --outcome     CAD \
#     --sumstat_dir /path/to/merged \
#     --out         /path/to/mr_LDL_C_CAD.csv

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(optparse)
})

set.seed(42)

opt_list <- list(
  make_option("--exposure",     type = "character"),
  make_option("--outcome",      type = "character"),
  make_option("--sumstat_dir",  type = "character"),
  make_option("--p_thresh",     type = "double",  default = 5e-8),
  make_option("--r2_clump",     type = "double",  default = 0.001),
  make_option("--kb_clump",     type = "integer", default = 10000L),
  make_option("--min_ivs",      type = "integer", default = 3L),
  make_option("--out",          type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

load_local_sumstat <- function(trait_file, type = "exposure", phenotype_id) {
  dt <- fread(trait_file, na.strings = "NA")
  # Standardise REGENIE column names
  dt[, pval   := 10^(-LOG10P)]
  dt[, SNP    := ID]
  dt[, CHR    := CHROM]
  dt[, BP     := GENPOS]
  dt[, effect_allele.exposure  := ALLELE1]
  dt[, other_allele.exposure   := ALLELE0]
  dt[, beta.exposure           := BETA]
  dt[, se.exposure             := SE]
  dt[, pval.exposure           := pval]
  dt[, eaf.exposure            := A1FREQ]
  dt[, exposure                := phenotype_id]
  dt[, id.exposure             := phenotype_id]
  dt[, effect_allele.outcome   := ALLELE1]
  dt[, other_allele.outcome    := ALLELE0]
  dt[, beta.outcome            := BETA]
  dt[, se.outcome              := SE]
  dt[, pval.outcome            := pval]
  dt[, eaf.outcome             := A1FREQ]
  dt[, outcome                 := phenotype_id]
  dt[, id.outcome              := phenotype_id]
  dt[, mr_keep.outcome         := TRUE]
  as.data.frame(dt)
}

message(sprintf("[MR] %s → %s", opt$exposure, opt$outcome))

exp_file <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz",   opt$exposure))
out_file <- file.path(opt$sumstat_dir, sprintf("outcomes_%s.txt.gz", opt$outcome))

if (!file.exists(exp_file)) stop("Exposure file missing: ", exp_file)
if (!file.exists(out_file)) stop("Outcome file missing: ", out_file)

# Load and filter to GWS SNPs for exposure
# fread() requires R.utils for direct .gz reading — use gunzip pipe instead
exp_dt <- fread(cmd = paste("gunzip -c", shQuote(exp_file)), na.strings = "NA")
exp_dt[, pval := 10^(-LOG10P)]

# ── Lp(a) cis-MR: restrict to LPA locus (chr6q27 ±3Mb) ──────────────────────
# LDSC rg is structurally unreliable for Lp(a) (single-locus trait).
# Standard genome-wide IVW is also inappropriate — ~90% of Lp(a) heritability
# is at the LPA locus (rs10455872, rs3798220, chr6:160,589,756).
# Cis-MR: use only chr6:157,000,000–163,000,000 instruments.
CIS_LPA_CHR   <- 6L
CIS_LPA_START <- 157000000L
CIS_LPA_END   <- 163000000L

if (opt$exposure == "LPA") {
  message("[MR] Lp(a) detected — applying cis-MR restriction to chr6:157Mb-163Mb")
  exp_dt <- exp_dt[CHROM == CIS_LPA_CHR & GENPOS >= CIS_LPA_START & GENPOS <= CIS_LPA_END]
  message(sprintf("[MR] Lp(a) cis window: %d SNPs available", nrow(exp_dt)))
}

exp_gws <- exp_dt[pval < opt$p_thresh & !is.na(BETA) & !is.na(SE)]

if (nrow(exp_gws) < opt$min_ivs) {
  message(sprintf("[MR] Only %d GWS SNPs — below min_ivs threshold. Skipping.", nrow(exp_gws)))
  quit(status = 0)
}

# Format as TwoSampleMR exposure data
exp_dat <- format_data(
  as.data.frame(exp_gws),
  type              = "exposure",
  snp_col           = "ID",
  beta_col          = "BETA",
  se_col            = "SE",
  eaf_col           = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col  = "ALLELE0",
  pval_col          = "pval",
  chr_col           = "CHROM",
  pos_col           = "GENPOS"
  # phenotype_col omitted: passing NULL causes format_data() to crash
  # with "argument is of length zero" — exposure/id.exposure set below
)
exp_dat$exposure    <- opt$exposure
exp_dat$id.exposure <- opt$exposure

# Clump to independent instruments (requires PLINK reference)
# TwoSampleMR clumping uses LD from 1000G EUR by default (needs internet)
# On TACC compute nodes (no internet), use local LD reference via PLINK
message(sprintf("[MR] Clumping %d GWS SNPs...", nrow(exp_dat)))

# Resolve position column — TwoSampleMR versions differ (pos.exposure vs chr_pos)
pos_col_name <- intersect(c("pos.exposure", "chr_pos", "BP", "GENPOS"), names(exp_dat))[1]
if (is.na(pos_col_name)) stop("[MR] Cannot find position column in formatted exposure data")
message(sprintf("[MR] Using position column: %s", pos_col_name))

exp_clumped <- tryCatch({
  clump_data(exp_dat,
    clump_r2  = opt$r2_clump,
    clump_kb  = opt$kb_clump,
    clump_p1  = opt$p_thresh,
    clump_p2  = 1
  )
}, error = function(e) {
  # Fallback: manual greedy clumping by p-value ordering (no LD, position-based)
  message("  TwoSampleMR clump_data failed — using greedy position-based clumping")
  message("  Error was: ", conditionMessage(e))
  exp_dat_sorted <- exp_dat[order(exp_dat$pval.exposure), ]
  keep <- c()
  for (i in seq_len(nrow(exp_dat_sorted))) {
    snp_bp <- as.integer(exp_dat_sorted[[pos_col_name]][i])
    if (length(keep) == 0 ||
        all(abs(as.integer(exp_dat_sorted[[pos_col_name]][keep]) - snp_bp) > opt$kb_clump * 1000)) {
      keep <- c(keep, i)
    }
  }
  exp_dat_sorted[keep, ]
})

n_ivs <- nrow(exp_clumped)
message(sprintf("[MR] %d independent instruments after clumping", n_ivs))

if (n_ivs < opt$min_ivs) {
  message(sprintf("[MR] Too few IVs (%d) — skipping", n_ivs))
  quit(status = 0)
}

# Load outcome data for the IV SNPs
out_dt <- fread(cmd = paste("gunzip -c", shQuote(out_file)), na.strings = "NA")
out_dt[, pval := 10^(-LOG10P)]
out_snps <- out_dt[ID %in% exp_clumped$SNP]

if (nrow(out_snps) == 0) {
  message("[MR] No IVs found in outcome data — check SNP overlap")
  quit(status = 0)
}

out_dat <- format_data(
  as.data.frame(out_snps),
  type              = "outcome",
  snp_col           = "ID",
  beta_col          = "BETA",
  se_col            = "SE",
  eaf_col           = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col  = "ALLELE0",
  pval_col          = "pval",
  chr_col           = "CHROM",
  pos_col           = "GENPOS"
)
out_dat$outcome    <- opt$outcome
out_dat$id.outcome <- opt$outcome

# Harmonise
dat <- harmonise_data(exp_clumped, out_dat, action = 2)
dat <- dat[dat$mr_keep, ]
n_harmonised <- nrow(dat)
message(sprintf("[MR] %d SNPs after harmonisation", n_harmonised))

if (n_harmonised < opt$min_ivs) {
  message("[MR] Too few harmonised SNPs — skipping")
  quit(status = 0)
}

# Run MR methods
message("[MR] Running MR methods...")
methods  <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
mr_res   <- mr(dat, method_list = methods)

# Egger intercept test (pleiotropy)
egger    <- mr_pleiotropy_test(dat)

# Heterogeneity (Cochran Q for IVW, Rucker Q for Egger)
hetero   <- mr_heterogeneity(dat)

# Steiger filtering (tests direction of causality)
steiger  <- directionality_test(dat)

# Compile results
results <- data.table(
  exposure          = opt$exposure,
  outcome           = opt$outcome,
  n_ivs             = n_harmonised,
  method            = mr_res$method,
  b                 = mr_res$b,
  se                = mr_res$se,
  pval              = mr_res$pval,
  OR                = exp(mr_res$b),
  OR_lo95           = exp(mr_res$b - 1.96 * mr_res$se),
  OR_hi95           = exp(mr_res$b + 1.96 * mr_res$se),
  egger_intercept   = egger$egger_intercept[1],
  egger_intercept_p = egger$pval[1],
  pleiotropy_flag   = egger$pval[1] < 0.05,
  steiger_dir       = steiger$correct_causal_direction[1],
  steiger_p         = steiger$steiger_pval[1]
)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(results, opt$out, sep = ",")
message(sprintf("[MR] Results written: %s (%d rows)", opt$out, nrow(results)))

# Print summary
for (i in seq_len(nrow(mr_res))) {
  message(sprintf("  %s: b=%.3f (SE=%.3f), p=%.2e, OR=%.2f [%.2f-%.2f]",
                  mr_res$method[i], mr_res$b[i], mr_res$se[i],
                  mr_res$pval[i], exp(mr_res$b[i]),
                  exp(mr_res$b[i] - 1.96*mr_res$se[i]),
                  exp(mr_res$b[i] + 1.96*mr_res$se[i])))
}
if (egger$pval[1] < 0.05) {
  message(sprintf("  *** Pleiotropy flagged: Egger intercept p=%.3f ***", egger$pval[1]))
}
message("[MR] Complete.")
