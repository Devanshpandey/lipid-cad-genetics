#!/usr/bin/env Rscript
# Fine-mapping with SuSiE at colocalized loci
# Inputs: merged sumstats + coloc results (PP.H4 > 0.8 loci)
# Outputs: credible sets CSV with top variants per locus
#
# Usage:
#   Rscript run_finemapping.R \
#     --coloc_results /path/to/coloc_results_combined.csv \
#     --sumstat_dir   /path/to/merged \
#     --ld_plink_dir  /path/to/ukb_plink \
#     --out           /path/to/finemapping_credible_sets.csv

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(optparse)
  library(Matrix)
})

set.seed(42)

opt_list <- list(
  make_option("--coloc_results", type = "character"),
  make_option("--sumstat_dir",   type = "character"),
  make_option("--ld_plink_dir",  type = "character"),
  make_option("--window_kb",     type = "integer", default = 500L),
  make_option("--pp4_thresh",    type = "double",  default = 0.8),
  make_option("--out",           type = "character", default = "finemapping_credible_sets.csv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

message("[finemapping] Loading colocalized loci...")
coloc_res <- fread(opt$coloc_results)
coloc_sig <- coloc_res[PP.H4 >= opt$pp4_thresh]
message(sprintf("[finemapping] %d colocalized loci (PP.H4 >= %.1f)", nrow(coloc_sig), opt$pp4_thresh))

if (nrow(coloc_sig) == 0) {
  message("[finemapping] No colocalized loci — exiting")
  quit(status = 0)
}

# Prioritise loci: deduplicate by sentinel SNP (same locus may appear for multiple
# exposure×outcome pairs — run SuSiE once per unique locus×exposure, highest PP.H4 first)
setorder(coloc_sig, -PP.H4)
coloc_sig[, locus_key := paste(exposure, chr, round(bp / 1e6), sep = "_")]
coloc_sig <- coloc_sig[!duplicated(locus_key)]
message(sprintf("[finemapping] %d unique exposure×locus combinations after deduplication", nrow(coloc_sig)))

compute_ld_matrix <- function(snp_list, chr, geno_prefix, tmp_dir) {
  snp_file <- file.path(tmp_dir, "snps.txt")
  writeLines(snp_list, snp_file)
  ld_prefix <- file.path(tmp_dir, "ld")
  # Use plink1.9 for LD: --r square produces correlation matrix (not r²)
  # plink2 --r-unphased flag not available in v2.00a5LM (June 2023 build)
  plink1_bin <- "/path/to/software/plink"
  if (!file.exists(plink1_bin)) {
    plink2_bin <- Sys.getenv("PLINK2", unset = "/path/to/software/plink2")
    plink1_bin <- sub("plink2$", "plink", plink2_bin)
  }

  # Remove stale ld.* files from previous runs
  file.remove(list.files(tmp_dir, pattern = "^ld\\.", full.names = TRUE))

  # plink1.9 --r square → <prefix>.ld (space-separated square correlation matrix)
  cmd <- sprintf(
    "%s --bfile %s --chr %s --extract %s --r square --out %s --threads 4 2>&1 | tail -3",
    plink1_bin, geno_prefix, chr, snp_file, ld_prefix
  )
  system(cmd)

  ld_file <- paste0(ld_prefix, ".ld")
  if (!file.exists(ld_file)) {
    message("  plink1.9 --r square produced no .ld file — skipping locus")
    return(NULL)
  }
  as.matrix(fread(ld_file, header = FALSE))
}

tmp_dir <- tempdir()
dir.create(tmp_dir, showWarnings = FALSE)
results <- list()

for (i in seq_len(nrow(coloc_sig))) {
  row       <- coloc_sig[i]
  exposure  <- row$exposure
  outcome   <- row$outcome
  sentinel  <- row$sentinel_snp
  chr       <- row$chr
  bp        <- row$bp
  w         <- opt$window_kb * 1000L

  message(sprintf("[finemapping] Locus %d/%d: %s vs %s at %s (chr%s:%s)",
                  i, nrow(coloc_sig), exposure, outcome, sentinel, chr, bp))

  # Load sumstats in window (fread requires gunzip pipe — R.utils not available on TACC)
  exp_file <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", exposure))
  exp_dt   <- fread(cmd = paste("gunzip -c", shQuote(exp_file)), na.strings = "NA")
  exp_dt[, pval := 10^(-LOG10P)]
  exp_win  <- exp_dt[CHROM == chr & GENPOS >= bp - w & GENPOS <= bp + w & !is.na(BETA) & !is.na(SE)]

  if (nrow(exp_win) < 10) next

  # Compute LD matrix using plink1.9
  ld_mat <- compute_ld_matrix(exp_win$ID, chr, opt$ld_plink_dir, tmp_dir)
  if (is.null(ld_mat)) {
    message("  LD matrix computation failed — skipping locus")
    next
  }

  # Match SNP order to LD matrix
  n_snps <- min(nrow(exp_win), nrow(ld_mat))
  exp_win <- exp_win[1:n_snps]
  ld_mat  <- ld_mat[1:n_snps, 1:n_snps]

  # Remove SNPs with any NA in their LD row (monomorphic or missing genotypes)
  # These arise in high-LD regions (APOE locus) and cause XtX to have NAs
  na_rows <- apply(is.na(ld_mat), 1, any)
  if (any(na_rows)) {
    n_na <- sum(na_rows)
    message(sprintf("  Removing %d SNPs with NA in LD matrix", n_na))
    exp_win <- exp_win[!na_rows]
    ld_mat  <- ld_mat[!na_rows, !na_rows]
  }
  if (nrow(exp_win) < 10) {
    message("  Too few SNPs after NA removal — skipping locus")
    next
  }

  # Force symmetry (floating-point rounding from plink1.9 --r square)
  ld_mat <- (ld_mat + t(ld_mat)) / 2
  diag(ld_mat) <- 1.0

  # Run SuSiE with z-scores
  z_scores <- exp_win$BETA / exp_win$SE
  n_samples <- median(exp_win$N, na.rm = TRUE)

  susie_res <- tryCatch({
    susie_rss(
      z          = z_scores,
      R          = ld_mat,
      n          = n_samples,
      L          = 10,
      coverage   = 0.95,
      min_abs_corr = 0.5
    )
  }, error = function(e) {
    message("  SuSiE error: ", e$message)
    NULL
  })

  if (is.null(susie_res)) next

  # Extract credible sets
  cs <- susie_res$sets$cs
  if (length(cs) == 0) {
    message("  No credible sets found")
    next
  }

  for (cs_idx in seq_along(cs)) {
    cs_snps <- exp_win$ID[cs[[cs_idx]]]
    cs_pip  <- susie_res$pip[cs[[cs_idx]]]
    top_snp <- cs_snps[which.max(cs_pip)]
    results[[length(results) + 1]] <- data.table(
      exposure    = exposure,
      outcome     = outcome,
      sentinel    = sentinel,
      chr         = chr,
      locus_bp    = bp,
      cs_id       = cs_idx,
      cs_size     = length(cs_snps),
      top_snp     = top_snp,
      top_pip     = max(cs_pip),
      cs_snps     = paste(cs_snps, collapse = ";"),
      cs_pips     = paste(round(cs_pip, 4), collapse = ";"),
      PP.H4       = row$PP.H4
    )
  }
}

if (length(results) > 0) {
  out_dt <- rbindlist(results)
  fwrite(out_dt, opt$out, sep = ",")
  message(sprintf("[finemapping] %d credible sets written: %s", nrow(out_dt), opt$out))
} else {
  message("[finemapping] No credible sets identified")
}
message("[finemapping] Complete.")
