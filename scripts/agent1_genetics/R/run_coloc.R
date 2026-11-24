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
#     --sentinels    /path/to/sentinels_LDL_C.tsv \
#     --out          /path/to/coloc_results.csv
#
# Sentinels come from make_coloc_sentinels.R — run that first, once per
# exposure. This script no longer derives them (see --sentinels below).

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
  # Sentinels are defined ONCE per exposure by make_coloc_sentinels.R and
  # cached. They are never derived here: doing so re-ran an expensive clump
  # for each of the 9 outcomes and, worse, hid a silent fallback that
  # collapsed the sentinel set to a handful of loci (see that script's header).
  make_option("--sentinels",    type = "character"),
  make_option("--window_kb",    type = "integer", default = 500L),
  make_option("--pp4_thresh",   type = "double",  default = 0.8),
  make_option("--out",          type = "character", default = "coloc_results.csv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

load_sumstat <- function(trait, dir) {
  f <- file.path(dir, paste0(trait, ".txt.gz"))
  if (!file.exists(f)) stop(sprintf("File not found: %s", f))
  # fread() requires R.utils for direct .gz reading — use gunzip pipe instead
  dt <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA")
  # NOTE: do not add a `pval := 10^(-log10p)` column. It underflows to exactly
  # 0 for -log10P > 323.3, and any ranking or `== min()` comparison on it then
  # silently selects every variant below the double-precision floor. Threshold
  # and sort on log10p instead.
  # Standardise column names from REGENIE output
  setnames(dt,
    old = intersect(c("ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1",
                      "A1FREQ", "N", "BETA", "SE", "LOG10P"),
                    names(dt)),
    new = intersect(c("SNP", "CHR", "BP", "A2", "A1",
                      "freq", "N", "beta", "se", "log10p"),
                    c("SNP", "CHR", "BP", "A2", "A1", "freq", "N", "beta", "se", "log10p"))
  )
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

# ── Sentinels: read the cached, pre-clumped list ──────────────────────────
# Produced by make_coloc_sentinels.R. Required, and validated here: an empty
# or missing file is a hard error, never a reason to guess a sentinel set.
if (is.null(opt$sentinels) || !nzchar(opt$sentinels)) {
  stop("--sentinels is required (run make_coloc_sentinels.R first)")
}
if (!file.exists(opt$sentinels)) {
  stop(sprintf("Sentinel file not found: %s", opt$sentinels))
}
sent_dt <- fread(opt$sentinels)
if (nrow(sent_dt) == 0) {
  message("[coloc] Sentinel file is empty — no GWS loci for this exposure. Exiting.")
  quit(status = 0)
}
if (!"SNP" %in% names(sent_dt)) stop("Sentinel file has no SNP column: ", opt$sentinels)
sentinels <- sent_dt$SNP

message(sprintf("[coloc] %d sentinel loci on %d chromosomes (from %s)",
                length(sentinels), uniqueN(sent_dt$CHR), basename(opt$sentinels)))

# Run coloc at each locus.
# Split by chromosome first. run_coloc_window filters with CHR == chr & BP in
# range, which is a full vector scan of the entire table for EVERY sentinel;
# with a few hundred sentinels over millions of rows that dominates runtime and
# is pure waste, since a sentinel can only match rows on its own chromosome.
# Splitting once makes each window scan one chromosome instead of the genome.
# Performance only -- the selected rows are identical.
exp_by_chr <- split(exp_dt, exp_dt$CHR)
out_by_chr <- split(out_dt, out_dt$CHR)

results <- rbindlist(lapply(seq_len(nrow(sent_dt)), function(i) {
  s  <- sent_dt$SNP[i]
  ch <- as.character(sent_dt$CHR[i])
  e  <- exp_by_chr[[ch]]; o <- out_by_chr[[ch]]
  if (is.null(e) || is.null(o)) return(NULL)
  run_coloc_window(e, o, s, opt$window_kb, opt$pp4_thresh)
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
