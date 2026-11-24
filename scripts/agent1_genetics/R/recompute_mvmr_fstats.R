#!/usr/bin/env Rscript
# ============================================================
# Recompute MVMR instrument-strength statistics from data
#
# Why this script exists:
#   The four-exposure MVMR reports two different F-statistics that are
#   easily conflated:
#     first-stage F   — strength of each exposure's own instruments
#     conditional F   — Sanderson-Windmeijer strength of each exposure
#                       AFTER projecting out the other exposures, in the
#                       shared (union) instrument set actually used
#   The manuscript quotes first-stage F = 164 and conditional F = 13 for
#   Lp(a). This script regenerates both from the summary statistics so the
#   quoted values are traceable rather than hardcoded.
#
#   It also emits the diagnostic that explains the gap: how many of the
#   union instruments actually carry each exposure, and specifically how
#   many independent LPA-region variants survive the 10 Mb distance rule.
#   Lp(a) is close to monogenic, so a genome-wide selection rule shared
#   with polygenic exposures retains almost no LPA signal.
#
# Usage:
#   Rscript recompute_mvmr_fstats.R \
#     --outcome     CAD \
#     --sumstat_dir /path/to/merged \
#     --out         exports/agent1/mvmr_fstats_CAD.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})
set.seed(42)

has_mvmr <- requireNamespace("MVMR", quietly = TRUE)
if (has_mvmr) library(MVMR) else
  message("[FSTAT] MVMR package not found — conditional F will be NA. install.packages('MVMR')")

opt_list <- list(
  make_option("--outcome",     type = "character", help = "Outcome trait name (e.g. CAD)"),
  make_option("--sumstat_dir", type = "character", help = "Dir with merged REGENIE sumstats"),
  make_option("--p_thresh",    type = "double",  default = 5e-8),
  make_option("--kb_clump",    type = "integer", default = 10000L,
              help = "Distance-pruning window in kb (matches run_mvmr.R)"),
  make_option("--out",         type = "character", help = "Output CSV path")
)
opt <- parse_args(OptionParser(option_list = opt_list))

LIPID_TRAITS <- c("LDL_C", "HDL_C", "TRIGLY", "LPA")
# LPA locus. NOTE: these REGENIE sumstats are GRCh37/hg19 -- rs10455872 is at
# chr6:161,010,118 here, matching CARDIoGRAMplusC4D's bp_hg19 column (the hg38
# position is 160,589,086, ~420 kb lower). hg19 LPA spans 160,952,515-161,087,407;
# the window below adds flanks for LPAL2/PLG/SLC22A3 and contains both lead
# variants (rs10455872 161,010,118; rs3798220 160,961,137).
LPA_CHR <- "6"; LPA_START <- 160500000L; LPA_END <- 161700000L

# ── identical pruning rule to run_mvmr.R: distance-only, no r2 ──
# NOTE: order by LOG10P descending, NOT by pval ascending. pval = 10^(-LOG10P)
# underflows to exactly 0 for any LOG10P > ~323, which for Lp(a) is 1,085 variants.
# order(pval) then ties them all and data.table's stable sort keeps whichever comes
# first in FILE order -- rs531596603 (chr6:160,179,747, LOG10P 438) rather than
# rs10455872 (161,010,118, LOG10P 14,632). This is the bug in run_mvmr.R.
greedy_clump <- function(dt, p_thresh, kb_window) {
  gws <- dt[LOG10P > -log10(p_thresh) & !is.na(BETA) & !is.na(SE)][order(-LOG10P)]
  if (!nrow(gws)) return(character(0))
  keep_snps <- character(0); keep_pos <- integer(0); keep_chr <- character(0)
  for (i in seq_len(nrow(gws))) {
    chr <- as.character(gws$CHROM[i]); pos <- as.integer(gws$GENPOS[i])
    same <- keep_chr == chr
    if (!length(keep_pos) ||
        !any(same & abs(keep_pos[same] - pos) < kb_window * 1000L)) {
      keep_snps <- c(keep_snps, gws$ID[i])
      keep_pos  <- c(keep_pos, pos); keep_chr <- c(keep_chr, chr)
    }
  }
  keep_snps
}

# ── 1. Load exposures ──
lipid_data <- list()
for (tr in LIPID_TRAITS) {
  f <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", tr))
  if (!file.exists(f)) { message("  WARNING: missing ", f); next }
  dt <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA")
  dt[, pval := 10^(-LOG10P)]
  lipid_data[[tr]] <- dt
}
active <- names(lipid_data)

# ── 2. Per-trait instruments + first-stage F ──
# first-stage F per variant = (beta/se)^2; per trait we report the mean over
# that trait's own independent instruments (the conventional summary).
per_trait <- list(); ivs_by_trait <- list()
for (tr in active) {
  ivs <- greedy_clump(lipid_data[[tr]], opt$p_thresh, opt$kb_clump)
  ivs_by_trait[[tr]] <- ivs
  d <- lipid_data[[tr]][ID %in% ivs]
  Fv <- (d$BETA / d$SE)^2
  n_lpa <- lipid_data[[tr]][ID %in% ivs & as.character(CHROM) == LPA_CHR &
                            GENPOS >= LPA_START & GENPOS <= LPA_END, .N]
  # how many independent LPA-region signals exist BEFORE the 10 Mb rule
  gws_locus <- lipid_data[[tr]][LOG10P > -log10(opt$p_thresh) & as.character(CHROM) == LPA_CHR &
                                GENPOS >= LPA_START & GENPOS <= LPA_END]
  per_trait[[tr]] <- data.table(
    exposure               = tr,
    n_own_instruments      = length(ivs),
    first_stage_F_mean     = mean(Fv, na.rm = TRUE),
    first_stage_F_median   = median(Fv, na.rm = TRUE),
    n_instruments_at_LPA   = n_lpa,
    n_gws_variants_at_LPA  = nrow(gws_locus)
  )
  message(sprintf("[FSTAT] %-7s own IVs %4d | first-stage F %6.1f | at LPA locus %d (of %d GWS)",
                  tr, length(ivs), mean(Fv, na.rm = TRUE), n_lpa, nrow(gws_locus)))
}
per_trait <- rbindlist(per_trait)

# ── 3. Union instrument set, pruned again at the same window ──
all_ivs <- Reduce(union, ivs_by_trait)
snp_ref <- rbindlist(lapply(active, function(tr)
  lipid_data[[tr]][ID %in% all_ivs, .(ID, CHROM, GENPOS, LOG10P)]))
snp_ref <- snp_ref[snp_ref[, .I[which.max(LOG10P)], by = ID]$V1][order(-LOG10P)]
final_ivs <- character(0); fpos <- integer(0); fchr <- character(0)
for (i in seq_len(nrow(snp_ref))) {
  chr <- as.character(snp_ref$CHROM[i]); pos <- as.integer(snp_ref$GENPOS[i])
  same <- fchr == chr
  if (!length(fpos) || !any(same & abs(fpos[same] - pos) < opt$kb_clump * 1000L)) {
    final_ivs <- c(final_ivs, snp_ref$ID[i]); fpos <- c(fpos, pos); fchr <- c(fchr, chr)
  }
}
message(sprintf("[FSTAT] Union instrument set: %d SNPs", length(final_ivs)))

# ── 4. Aligned exposure matrices ──
ref <- lipid_data[[active[1]]][ID %in% final_ivs, .(ID, ALLELE1, ALLELE0)]
bx <- bxse <- matrix(NA_real_, length(final_ivs), length(active),
                     dimnames = list(final_ivs, active))
for (tr in active) {
  d <- merge(lipid_data[[tr]][ID %in% final_ivs], ref[, .(ID, rA1 = ALLELE1, rA0 = ALLELE0)],
             by = "ID", all.x = TRUE)
  d[, ba := fifelse(ALLELE1 == rA1, BETA, fifelse(ALLELE1 == rA0, -BETA, NA_real_))]
  idx <- match(d$ID, final_ivs)
  bx[idx, tr] <- d$ba; bxse[idx, tr] <- d$SE
}

# ── 5. Outcome ──
of <- file.path(opt$sumstat_dir, sprintf("outcomes_%s.txt.gz", opt$outcome))
if (!file.exists(of)) stop("Outcome file missing: ", of)
od <- fread(cmd = paste("gunzip -c", shQuote(of)), na.strings = "NA")
od <- merge(od[ID %in% final_ivs], ref[, .(ID, rA1 = ALLELE1, rA0 = ALLELE0)],
            by = "ID", all.x = TRUE)
od[, ba := fifelse(ALLELE1 == rA1, BETA, fifelse(ALLELE1 == rA0, -BETA, NA_real_))]
od <- od[match(final_ivs, ID)]

keep <- complete.cases(bx) & complete.cases(bxse) & !is.na(od$ba) & !is.na(od$SE)
bx <- bx[keep, , drop = FALSE]; bxse <- bxse[keep, , drop = FALSE]
by <- od$ba[keep]; byse <- od$SE[keep]
message(sprintf("[FSTAT] %d SNPs with complete exposure + outcome data", nrow(bx)))

# ── 6. Sanderson-Windmeijer conditional F ──
cond_F <- setNames(rep(NA_real_, length(active)), active)
if (has_mvmr) {
  fmt <- format_mvmr(BXGs = bx, BYG = by, seBXGs = bxse, seBYG = byse,
                     RSID = rownames(bx))
  st  <- strength_mvmr(r_input = fmt, gencov = 0)
  cond_F[seq_along(active)] <- as.numeric(st[1, seq_along(active)])
}

# ── 7. Write ──
res <- merge(per_trait,
             data.table(exposure = active,
                        conditional_F = as.numeric(cond_F[active]),
                        n_union_instruments = nrow(bx)),
             by = "exposure", sort = FALSE)
res[, outcome := opt$outcome]
res[, F_ratio_firststage_over_conditional := first_stage_F_mean / conditional_F]
setcolorder(res, c("outcome", "exposure", "n_own_instruments", "n_union_instruments",
                   "first_stage_F_mean", "first_stage_F_median", "conditional_F",
                   "F_ratio_firststage_over_conditional",
                   "n_instruments_at_LPA", "n_gws_variants_at_LPA"))
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, opt$out)
print(res)
message("[FSTAT] Wrote ", opt$out)
