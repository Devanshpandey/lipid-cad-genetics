#!/usr/bin/env Rscript
# Define independent lipid sentinels for colocalization — ONCE per exposure.
#
# Why this script exists
# ----------------------
# Sentinel selection used to live inside run_coloc.R, where it (a) re-ran for
# every one of the 9 outcomes and (b) read the LD panel from an environment
# variable that config/tacc_paths.sh never exported. PLINK therefore got an
# empty --bfile, failed silently behind 2>/dev/null, and a fallback branch
# picked `sig_snps[pval == min(pval)]`. Because pval was computed as
# 10^(-LOG10P), everything with -log10P > 323.3 underflowed to exactly 0, so
# "the most significant SNP" resolved to *every variant below the double
# precision floor* — i.e. only the 2-5 most extreme loci per trait.
#
# Two rules follow, and both are enforced below:
#   1. The LD panel is an explicit argument. Never read it from the
#      environment, and stop() if it is missing.
#   2. Never rank or compare on 10^(-LOG10P). LOG10P is the primitive; keep
#      every threshold and every sort on that scale.
#
# Usage:
#   Rscript make_coloc_sentinels.R \
#     --trait LDL_C --sumstat_dir <dir> --ld_bfile <plink prefix> \
#     --plink <plink1.9 binary> --out <sentinels.tsv>

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

set.seed(42)

opt_list <- list(
  make_option("--trait",       type = "character"),
  make_option("--sumstat_dir", type = "character"),
  make_option("--ld_bfile",    type = "character"),
  make_option("--plink",       type = "character",
              default = "/path/to/software/plink"),
  make_option("--gwas_p_thresh", type = "double", default = 5e-8),
  make_option("--clump_r2",    type = "double",  default = 0.001),
  make_option("--clump_kb",    type = "integer", default = 1000L),
  # Coloc uses +/-500 kb windows, so sentinels closer than this produce
  # near-duplicate windows and near-duplicate coloc rows. Collapse them.
  make_option("--min_sep_kb",  type = "integer", default = 500L),
  make_option("--out",         type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (req in c("trait", "sumstat_dir", "ld_bfile", "out")) {
  if (is.null(opt[[req]]) || !nzchar(opt[[req]])) {
    stop(sprintf("--%s is required and must be non-empty", req))
  }
}
if (!file.exists(paste0(opt$ld_bfile, ".bed"))) {
  stop(sprintf("LD panel not found: %s.bed", opt$ld_bfile))
}
if (!file.exists(opt$plink)) {
  stop(sprintf("PLINK 1.9 binary not found: %s  (plink2 does not support --clump)",
               opt$plink))
}

LOG10P_THRESH <- -log10(opt$gwas_p_thresh)

f <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", opt$trait))
if (!file.exists(f)) stop(sprintf("Summary statistics not found: %s", f))

message(sprintf("[sentinels] %s  (GWS: -log10P > %.5f)", opt$trait, LOG10P_THRESH))
dt <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA",
            select = c("ID", "CHROM", "GENPOS", "LOG10P"))
setnames(dt, c("ID", "CHROM", "GENPOS", "LOG10P"), c("SNP", "CHR", "BP", "LOG10P"))

sig <- dt[!is.na(LOG10P) & LOG10P > LOG10P_THRESH]
message(sprintf("[sentinels]   %d genome-wide-significant variants on %d chromosomes",
                nrow(sig), uniqueN(sig$CHR)))
if (nrow(sig) == 0) {
  fwrite(data.table(SNP = character(), CHR = integer(), BP = integer(),
                    LOG10P = numeric()), opt$out, sep = "\t")
  quit(status = 0)
}

# ---- PLINK clump -----------------------------------------------------------
# PLINK parses the P column into a double, so 1e-660 would underflow to 0 the
# same way R's did. Floor at 1e-300 (inside the normal double range) and sort
# by LOG10P so the strongest variant at a locus is offered first; ties at the
# floor only affect WHICH variant indexes an already-discovered locus, not
# which loci are discovered, and the +/-500 kb window is insensitive to it.
tmp_dir  <- tempdir()
tmp_pval <- file.path(tmp_dir, sprintf("clumpin_%s.txt", opt$trait))
setorder(sig, -LOG10P)
fwrite(sig[, .(SNP, P = pmax(10^(-LOG10P), 1e-300))], tmp_pval, sep = "\t")

clump_prefix <- file.path(tmp_dir, sprintf("clump_%s", opt$trait))
clump_cmd <- sprintf(
  "%s --bfile %s --clump %s --clump-snp-field SNP --clump-field P --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --out %s",
  shQuote(opt$plink), shQuote(opt$ld_bfile), shQuote(tmp_pval),
  opt$gwas_p_thresh, opt$gwas_p_thresh, opt$clump_r2, opt$clump_kb,
  shQuote(clump_prefix)
)
message("[sentinels]   ", clump_cmd)
# No 2>/dev/null: a failure here must be visible, not silently absorbed.
rc <- system(clump_cmd)

clump_out <- paste0(clump_prefix, ".clumped")
if (rc != 0 || !file.exists(clump_out)) {
  log_file <- paste0(clump_prefix, ".log")
  if (file.exists(log_file)) {
    message("---- PLINK log ----"); writeLines(readLines(log_file)); message("-------------------")
  }
  stop(sprintf("PLINK --clump failed (exit %d, .clumped %s). Refusing to fall back to a heuristic sentinel set.",
               rc, if (file.exists(clump_out)) "present" else "missing"))
}

clumped <- fread(clump_out)
message(sprintf("[sentinels]   PLINK returned %d index SNPs on %d chromosomes",
                nrow(clumped), uniqueN(clumped$CHR)))

sent <- sig[SNP %in% clumped$SNP][order(CHR, BP)]
if (nrow(sent) == 0) stop("No clumped index SNPs matched the summary statistics.")

# ---- collapse sentinels closer than min_sep_kb ------------------------------
# r2 < 0.001 leaves many conditionally-independent index SNPs inside a single
# heterogeneous locus (LDLR returned dozens). Their +/-500 kb coloc windows are
# effectively the same window, so keep the strongest per locus.
sep <- opt$min_sep_kb * 1000L
setorder(sent, -LOG10P)
keep <- logical(nrow(sent))
taken_chr <- integer(0); taken_bp <- integer(0)
for (i in seq_len(nrow(sent))) {
  ci <- sent$CHR[i]; bi <- sent$BP[i]
  if (any(taken_chr == ci & abs(taken_bp - bi) <= sep)) next
  keep[i] <- TRUE
  taken_chr <- c(taken_chr, ci); taken_bp <- c(taken_bp, bi)
}
sent <- sent[keep][order(CHR, BP)]

message(sprintf("[sentinels]   %d independent loci after >=%d kb separation",
                nrow(sent), opt$min_sep_kb))
fwrite(sent[, .(SNP, CHR, BP, LOG10P)], opt$out, sep = "\t")
message("[sentinels] wrote ", opt$out)
