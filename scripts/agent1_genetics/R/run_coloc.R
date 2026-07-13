#!/usr/bin/env Rscript
# Cross-trait colocalization: lipid loci × CAD/MACE loci
# Uses coloc.abf (Bayesian coloc) with PP.H4 > 0.8 threshold
# Window: ±500kb around each sentinel lipid SNP
#
# Usage:
#   Rscript run_coloc.R \
#     --lipid_trait  LDL_C \
#     --outcome      CAD \
#     --sumstat_dir  /path/to/gwas_summary_stats/merged \
#     --out          /path/to/coloc_results.csv

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(optparse)
})

set.seed(42)

opt_list <- list(
  make_option("--lipid_trait",  type = "character"),
  make_option("--outcome",      type = "character"),
  make_option("--sumstat_dir",  type = "character"),
  make_option("--window_kb",    type = "integer", default = 500L),
  make_option("--pp4_thresh",   type = "double",  default = 0.8),
  make_option("--gwas_p_thresh",type = "double",  default = 5e-8),
  make_option("--out",          type = "character", default = "coloc_results.csv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

load_sumstat <- function(trait, dir) {
  f <- file.path(dir, paste0(trait, ".txt.gz"))
  if (!file.exists(f)) stop(sprintf("File not found: %s", f))
  # fread() requires R.utils for direct .gz reading — use gunzip pipe instead
  dt <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA")
  # Standardise column names from REGENIE output
  setnames(dt,
    old = intersect(c("ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1",
                      "A1FREQ", "N", "BETA", "SE", "LOG10P"),
                    names(dt)),
    new = intersect(c("SNP", "CHR", "BP", "A2", "A1",
                      "freq", "N", "beta", "se", "log10p"),
                    c("SNP", "CHR", "BP", "A2", "A1", "freq", "N", "beta", "se", "log10p"))
  )
  dt[, pval := 10^(-log10p)]
  dt
}

run_coloc_window <- function(exposure_dt, outcome_dt, sentinel_snp, window_kb, pp4_thresh) {
  snp_row <- exposure_dt[SNP == sentinel_snp]
  if (nrow(snp_row) == 0) return(NULL)
  chr <- snp_row$CHR[1]
  bp  <- snp_row$BP[1]
  w   <- window_kb * 1000L

  exp_win <- exposure_dt[CHR == chr & BP >= bp - w & BP <= bp + w]
  out_win <- outcome_dt[CHR == chr  & BP >= bp - w & BP <= bp + w]

  shared <- intersect(exp_win$SNP, out_win$SNP)
  if (length(shared) < 50) return(NULL)

  exp_win <- exp_win[SNP %in% shared][order(BP)]
  out_win <- out_win[SNP %in% shared][order(BP)]
  setkey(exp_win, SNP); setkey(out_win, SNP)
  out_win <- out_win[exp_win$SNP]

  D1 <- list(
    beta   = exp_win$beta,
    varbeta= exp_win$se^2,
    N      = median(exp_win$N, na.rm = TRUE),
    type   = "quant",
    snp    = exp_win$SNP,
    MAF    = exp_win$freq
  )
  D2 <- list(
    beta   = out_win$beta,
    varbeta= out_win$se^2,
    N      = median(out_win$N, na.rm = TRUE),
    type   = "cc",
    snp    = out_win$SNP,
    MAF    = out_win$freq
  )

  tryCatch({
    res <- coloc.abf(D1, D2)$summary
    data.table(
      sentinel_snp = sentinel_snp,
      chr          = chr,
      bp           = bp,
      nsnps        = res["nsnps"],
      PP.H0        = res["PP.H0.abf"],
      PP.H1        = res["PP.H1.abf"],
      PP.H2        = res["PP.H2.abf"],
      PP.H3        = res["PP.H3.abf"],
      PP.H4        = res["PP.H4.abf"],
      coloc_sig    = res["PP.H4.abf"] >= pp4_thresh
    )
  }, error = function(e) {
    message("  coloc error at ", sentinel_snp, ": ", e$message)
    NULL
  })
}

message(sprintf("[coloc] %s vs %s", opt$lipid_trait, opt$outcome))

# ── Lp(a) skip: coloc is inappropriate for single-locus traits ───────────────
# ~90% of Lp(a) heritability is at LPA (chr6q27). coloc.abf assumes two
# independent association signals can be tested for shared causal variants
# across a genomic window — meaningless when the exposure has exactly one locus.
# Lp(a) → outcome causal estimates come from cis-MR (run_mr.R), not coloc.
if (opt$lipid_trait == "LPA") {
  message("[coloc] Lp(a) is a single-locus trait — coloc.abf not appropriate. Skipping.")
  message("[coloc] Use cis-MR results for Lp(a) → outcome evidence.")
  quit(status = 0)
}

message("[coloc] Loading summary statistics...")

exp_dt  <- load_sumstat(paste0("lipids_", opt$lipid_trait), opt$sumstat_dir)
out_dt  <- load_sumstat(paste0("outcomes_", opt$outcome),   opt$sumstat_dir)

# Identify genome-wide significant sentinel SNPs in the exposure
sig_snps <- exp_dt[pval < opt$gwas_p_thresh]
if (nrow(sig_snps) == 0) {
  message("[coloc] No GWS SNPs — exiting")
  quit(status = 0)
}

# LD clumping to get independent sentinels (p < 5e-8, r2 < 0.001)
# Requires PLINK — write temp file and clump
tmp_dir  <- tempdir()
tmp_pval <- file.path(tmp_dir, "pval.txt")
fwrite(sig_snps[, .(SNP, pval)], tmp_pval, sep = "\t")

# UKB_GENO_STEP2 exported by SLURM script via config/tacc_paths.sh
# NOTE: plink2 does not support --clump; use plink1.9 for clumping.
# plink1.9 is available at the same software directory.
geno <- Sys.getenv("UKB_GENO_STEP2")
plink1_bin <- "/path/to/software/plink"  # plink1.9
if (!file.exists(plink1_bin)) {
  plink1_bin <- sub("plink2$", "plink", Sys.getenv("PLINK2",
    unset = "/path/to/software/plink2"))
}
clump_cmd <- sprintf(
  "%s --bfile %s --clump %s --clump-snp-field SNP --clump-field pval --clump-p1 5e-8 --clump-p2 1e-4 --clump-r2 0.001 --clump-kb 10000 --out %s 2>/dev/null",
  plink1_bin, geno, tmp_pval, file.path(tmp_dir, "clump")
)
system(clump_cmd)
clump_out <- file.path(tmp_dir, "clump.clumped")
if (file.exists(clump_out)) {
  sentinels <- fread(clump_out)$SNP
} else {
  # Fall back: just top SNP per 1Mb window
  setkey(sig_snps, CHR, BP)
  sentinels <- sig_snps[pval == min(pval), SNP]
}
message(sprintf("[coloc] Found %d independent sentinel SNPs", length(sentinels)))

# Run coloc at each locus
results <- rbindlist(lapply(sentinels, function(s) {
  run_coloc_window(exp_dt, out_dt, s, opt$window_kb, opt$pp4_thresh)
}), fill = TRUE)

if (nrow(results) > 0) {
  results[, exposure := opt$lipid_trait]
  results[, outcome  := opt$outcome]
  setcolorder(results, c("exposure", "outcome", "sentinel_snp", "chr", "bp",
                         "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "coloc_sig"))
  fwrite(results, opt$out, sep = ",")
  n_coloc <- sum(results$coloc_sig, na.rm = TRUE)
  message(sprintf("[coloc] %d loci tested, %d colocalized (PP.H4 >= %.1f)",
                  nrow(results), n_coloc, opt$pp4_thresh))
} else {
  message("[coloc] No loci passed coloc analysis")
}
message("[coloc] Complete.")
