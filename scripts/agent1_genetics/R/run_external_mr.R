#!/usr/bin/env Rscript
# ============================================================
# Two-sample MR: UK Biobank lipid instruments -> EXTERNAL CAD outcome
#
# Why this script exists
# ----------------------
# The manuscript reports external replication in FinnGen (Table S7) and
# CARDIoGRAMplusC4D (Table S24), but no script in the repository produced
# them. Table S24 also omits Lp(a) entirely and the text states "Lp(a) not
# available" for CARDIoGRAMplusC4D -- which is wrong: cad.add.160614 carries
# rs10455872 (p=5.7e-39), rs55730499 (p=5.4e-39) and rs3798220 (p=4.7e-09).
# Lp(a) was never submitted to that analysis.
#
# This script makes both replications reproducible and includes Lp(a).
#
# Methodological fixes relative to run_mvmr.R
# -------------------------------------------
#  1. REAL r^2-based LD clumping via PLINK (run_mvmr.R declared --r2_clump but
#     never used it; its greedy_clump() prunes on distance alone, so a 10 Mb
#     window collapses the whole ~1.1 Mb LPA locus to a single variant --
#     which is why Lp(a) conditional F was 13 against 40-52 for the polygenic
#     exposures). Falls back to distance pruning only if no LD panel is given,
#     and says so in the output.
#  2. --cis_dense lets a near-monogenic exposure keep multiple conditionally
#     independent signals at its own locus (the design Wang et al. used for
#     Lp(a): r^2<=0.4 within 660 kb, and r^2<=0.01 within 1 Mb for their pQTL
#     MR). Without this, Lp(a) is a one-SNP instrument.
#  3. Palindromic variants resolved by allele frequency where the outcome
#     reports it, dropped only when genuinely ambiguous (MAF 0.42-0.58).
#
# Usage
#   Rscript run_external_mr.R \
#     --sumstat_dir   $AGENT1_OUT/gwas_summary_stats/merged \
#     --outcome_file  data/external/cardiogram/cad.add.160614.website.txt \
#     --outcome_format cardiogram \
#     --exposures     LDL_C,HDL_C,TRIGLY,LPA \
#     --mode          univariable \
#     --ld_ref        data/reference/1000G_EUR \
#     --cis_dense     "LPA:6:160500000:161700000:0.01:1000" \
#   NOTE: cis_dense coordinates are GRCh37/hg19 to match the REGENIE sumstats
#   (rs10455872 = chr6:161,010,118). run_mr.R uses a wider 157-163 Mb window;
#   the narrow window here is the deliberate choice for a cis instrument.
#     --out           exports/agent1/mr_cardiogram_with_LPA.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(optparse)
})
set.seed(42)

opt_list <- list(
  make_option("--sumstat_dir",    type = "character"),
  make_option("--outcome_file",   type = "character"),
  make_option("--outcome_format", type = "character", default = "cardiogram",
              help = "cardiogram | finngen | mvp | regenie"),
  make_option("--outcome_label",  type = "character", default = ""),
  make_option("--exposures",      type = "character", default = "LDL_C,HDL_C,TRIGLY,LPA"),
  make_option("--mode",           type = "character", default = "univariable",
              help = "univariable | mvmr"),
  make_option("--p_thresh",       type = "double",    default = 5e-8),
  make_option("--r2_clump",       type = "double",    default = 0.001),
  make_option("--kb_clump",       type = "integer",   default = 10000L),
  make_option("--ld_ref",         type = "character", default = "",
              help = "PLINK bfile prefix for r2 clumping. If empty, falls back to distance pruning."),
  make_option("--plink",          type = "character", default = "plink"),
  make_option("--cis_dense",      type = "character", default = "",
              help = "Semicolon-separated TRAIT:CHR:START:END:R2:KB entries"),
  make_option("--min_ivs",        type = "integer",   default = 3L),
  make_option("--out",            type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

EXPOSURES <- trimws(strsplit(opt$exposures, ",")[[1]])
OUTLAB <- if (nzchar(opt$outcome_label)) opt$outcome_label else
  toupper(opt$outcome_format)

# ── outcome column maps ─────────────────────────────────────
# Each entry maps the external file's columns onto a common schema.
FORMATS <- list(
  cardiogram = list(snp = "markername", ea = "effect_allele", oa = "noneffect_allele",
                    beta = "beta", se = "se_dgc", eaf = "effect_allele_freq", p = "p_dgc"),
  finngen    = list(snp = "rsids", ea = "alt", oa = "ref",
                    beta = "beta", se = "sebeta", eaf = "af_alt", p = "pval"),
  mvp        = list(snp = "SNP", ea = "A1", oa = "A2",
                    beta = "BETA", se = "SE", eaf = "A1FREQ", p = "P"),
  regenie    = list(snp = "ID", ea = "ALLELE1", oa = "ALLELE0",
                    beta = "BETA", se = "SE", eaf = "A1FREQ", p = NA),
  # MEGASTROKE (Malik 2018) METAL output. Allele1 is the effect allele and alleles are
  # lower-case; harmonize() upper-cases both sides before matching.
  megastroke = list(snp = "MarkerName", ea = "Allele1", oa = "Allele2",
                    beta = "Effect", se = "StdErr", eaf = "Freq1", p = "P-value")
)
if (!opt$outcome_format %in% names(FORMATS))
  stop("Unknown --outcome_format: ", opt$outcome_format)
FM <- FORMATS[[opt$outcome_format]]

is_palindromic <- function(a1, a2) {
  p <- paste0(toupper(a1), toupper(a2))
  p %in% c("AT", "TA", "CG", "GC")
}

# ── 1. Load UKB exposure sumstats ───────────────────────────
message("[EXT] Loading exposure summary statistics...")
lipid <- list()
for (tr in EXPOSURES) {
  f <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", tr))
  if (!file.exists(f)) { message("  WARNING: missing ", f, " -- skipping ", tr); next }
  d <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA")
  d[, pval := 10^(-LOG10P)]
  lipid[[tr]] <- d
  message(sprintf("  %s: %d variants", tr, nrow(d)))
}
EXPOSURES <- names(lipid)
if (!length(EXPOSURES)) stop("No exposure summary statistics found in ", opt$sumstat_dir)

# ── 2. Instrument selection ─────────────────────────────────
# r2-based clumping via PLINK when an LD panel is supplied; distance pruning
# otherwise. The method actually used is recorded in the output.
# pval = 10^(-LOG10P) underflows to exactly 0 above LOG10P ~323 (1,085 Lp(a)
# variants do). Selecting on it ties them all and picks by file order, missing the
# true lead. Everywhere below we rank on LOG10P instead. PLINK still needs a P
# column, so extreme values are rescaled monotonically onto [1e-300, p_thresh]:
# strictly order-preserving, never underflows, always below the clump threshold.
plink_safe_p <- function(log10p, p_thresh) {
  lo <- -log10(p_thresh); hi <- max(log10p, na.rm = TRUE)
  if (hi <= 300) return(10^(-log10p))
  10^(-(lo + (300 - lo) * (log10p - lo) / (hi - lo)))
}

clump_r2 <- function(d, r2, kb, tag) {
  gws <- d[LOG10P > -log10(opt$p_thresh) & !is.na(BETA) & !is.na(SE)]
  if (!nrow(gws)) return(character(0))
  tmp <- tempfile(); on.exit(unlink(paste0(tmp, c("", ".clumped", ".log", ".nosex")), force = TRUE))
  fwrite(gws[, .(SNP = ID, P = plink_safe_p(LOG10P, opt$p_thresh))], tmp, sep = " ")
  cmd <- sprintf("%s --bfile %s --clump %s --clump-p1 %g --clump-p2 1 --clump-r2 %g --clump-kb %d --out %s",
                 shQuote(opt$plink), shQuote(opt$ld_ref), shQuote(tmp),
                 opt$p_thresh, r2, as.integer(kb), shQuote(tmp))
  st <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
  cf <- paste0(tmp, ".clumped")
  if (st != 0 || !file.exists(cf)) {
    message(sprintf("  [%s] PLINK clump failed -- falling back to distance pruning", tag))
    return(NULL)
  }
  cl <- fread(cf)
  as.character(cl$SNP)
}

clump_distance <- function(d, kb) {
  gws <- d[LOG10P > -log10(opt$p_thresh) & !is.na(BETA) & !is.na(SE)][order(-LOG10P)]
  if (!nrow(gws)) return(character(0))
  ks <- character(0); kp <- integer(0); kc <- character(0)
  for (i in seq_len(nrow(gws))) {
    ch <- as.character(gws$CHROM[i]); po <- as.integer(gws$GENPOS[i])
    sm <- kc == ch
    if (!length(kp) || !any(abs(kp[sm] - po) < kb * 1000L)) {
      ks <- c(ks, gws$ID[i]); kp <- c(kp, po); kc <- c(kc, ch)
    }
  }
  ks
}

# cis-dense overrides, e.g. "LPA:6:160500000:161700000:0.01:1000"
cis_spec <- list()
if (nzchar(opt$cis_dense)) {
  for (e in strsplit(opt$cis_dense, ";")[[1]]) {
    p <- strsplit(trimws(e), ":")[[1]]
    if (length(p) != 6) stop("Bad --cis_dense entry (need TRAIT:CHR:START:END:R2:KB): ", e)
    cis_spec[[p[1]]] <- list(chr = p[2], start = as.numeric(p[3]), end = as.numeric(p[4]),
                             r2 = as.numeric(p[5]), kb = as.numeric(p[6]))
  }
}

ivs <- list(); how <- list()
for (tr in EXPOSURES) {
  d <- lipid[[tr]]
  sp <- cis_spec[[tr]]
  if (!is.null(sp)) {
    d <- d[as.character(CHROM) == sp$chr & GENPOS >= sp$start & GENPOS <= sp$end]
    r2 <- sp$r2; kb <- sp$kb; lab <- sprintf("cis-dense r2<%g/%gkb", r2, kb)
  } else {
    r2 <- opt$r2_clump; kb <- opt$kb_clump; lab <- sprintf("genome-wide r2<%g/%gkb", r2, kb)
  }
  sel <- NULL
  if (nzchar(opt$ld_ref)) sel <- clump_r2(d, r2, kb, tr)
  if (is.null(sel)) { sel <- clump_distance(d, kb); lab <- paste0(lab, " [DISTANCE FALLBACK]") }
  ivs[[tr]] <- sel; how[[tr]] <- lab
  message(sprintf("[EXT] %-7s %4d instruments (%s)", tr, length(sel), lab))
}

# ── 3. Load + harmonize outcome ─────────────────────────────
message("[EXT] Loading outcome: ", opt$outcome_file)
need <- unlist(FM[c("snp", "ea", "oa", "beta", "se")], use.names = FALSE)
rdr <- if (grepl("\\.gz$", opt$outcome_file))
  paste("gunzip -c", shQuote(opt$outcome_file)) else NULL
out <- if (is.null(rdr)) fread(opt$outcome_file, na.strings = c("NA", ".")) else
  fread(cmd = rdr, na.strings = c("NA", "."))
miss <- setdiff(need, names(out))
if (length(miss)) stop("Outcome file lacks column(s): ", paste(miss, collapse = ", "))

setnames(out, FM$snp, "SNP"); setnames(out, FM$ea, "oEA"); setnames(out, FM$oa, "oOA")
setnames(out, FM$beta, "oBETA"); setnames(out, FM$se, "oSE")
if (!is.na(FM$eaf) && FM$eaf %in% names(out)) setnames(out, FM$eaf, "oEAF") else out[, oEAF := NA_real_]
# FinnGen packs multiple rsIDs per row
if (opt$outcome_format == "finngen") out[, SNP := tstrsplit(SNP, ",", fixed = TRUE)[[1]]]
out <- out[!is.na(SNP) & SNP != "" & !is.na(oBETA) & !is.na(oSE)]
setkey(out, SNP)
message(sprintf("  outcome variants: %d", nrow(out)))

harmonize <- function(tr, snps) {
  ex <- lipid[[tr]][ID %in% snps,
                    .(SNP = ID, eEA = ALLELE1, eOA = ALLELE0,
                      eBETA = BETA, eSE = SE, eEAF = A1FREQ)]
  m <- merge(ex, out[, .(SNP, oEA, oOA, oBETA, oSE, oEAF)], by = "SNP")
  if (!nrow(m)) return(m)
  m[, `:=`(eEA = toupper(eEA), eOA = toupper(eOA),
           oEA = toupper(oEA), oOA = toupper(oOA))]
  m[, palin := is_palindromic(eEA, eOA)]
  # same orientation / flipped / strand-flipped
  flip_base <- function(x) chartr("ACGT", "TGCA", x)
  m[, aligned := fifelse(eEA == oEA & eOA == oOA,  1,
              fifelse(eEA == oOA & eOA == oEA, -1,
              fifelse(eEA == flip_base(oEA) & eOA == flip_base(oOA),  1,
              fifelse(eEA == flip_base(oOA) & eOA == flip_base(oEA), -1, NA_real_))))]
  # palindromes: resolve by frequency, drop if genuinely ambiguous
  m[palin == TRUE, aligned := fifelse(
      is.na(eEAF) | is.na(oEAF) | (eEAF > 0.42 & eEAF < 0.58), NA_real_,
      fifelse((eEAF - 0.5) * (oEAF - 0.5) > 0, 1, -1))]
  m <- m[!is.na(aligned)]
  # FinnGen (and other outcome files) list the same rsID on several rows at
  # multi-allelic sites -- 141,873 such rsIDs in R11. merge() on a non-unique key
  # expands cartesian-style, inflating the instrument count and narrowing intervals
  # because duplicated rows count as independent observations. Allele alignment
  # removes most, but keep exactly one row per rsID to be certain.
  nd <- sum(duplicated(m$SNP))
  if (nd > 0) {
    m <- m[!duplicated(SNP)]
    message(sprintf("  [%s] dropped %d duplicate rsID matches (multi-allelic outcome rows)", tr, nd))
  }
  m[, oBETA := oBETA * aligned]
  m[]
}

# ── 4. Estimate ─────────────────────────────────────────────
rows <- list()
if (opt$mode == "univariable") {
  for (tr in EXPOSURES) {
    h <- harmonize(tr, ivs[[tr]])
    message(sprintf("[EXT] %-7s %d instruments harmonized against outcome", tr, nrow(h)))
    if (nrow(h) < opt$min_ivs) {
      message(sprintf("  too few instruments (%d) -- recording as not estimable", nrow(h)))
      rows[[tr]] <- data.table(exposure = tr, outcome = OUTLAB, n_instruments = nrow(h),
                               method = NA_character_, OR_perSD = NA_real_,
                               CI_lo95 = NA_real_, CI_hi95 = NA_real_, pval = NA_real_,
                               egger_intercept = NA_real_, egger_intercept_p = NA_real_,
                               instrument_selection = how[[tr]])
      next
    }
    inp <- mr_input(bx = h$eBETA, bxse = h$eSE, by = h$oBETA, byse = h$oSE, snps = h$SNP)
    ivw <- mr_ivw(inp); wme <- mr_median(inp, weighting = "weighted")
    egg <- tryCatch(mr_egger(inp), error = function(e) NULL)
    mk <- function(nm, est, lo, hi, p, ei = NA_real_, eip = NA_real_)
      data.table(exposure = tr, outcome = OUTLAB, n_instruments = nrow(h), method = nm,
                 OR_perSD = exp(est), CI_lo95 = exp(lo), CI_hi95 = exp(hi), pval = p,
                 egger_intercept = ei, egger_intercept_p = eip,
                 instrument_selection = how[[tr]])
    r <- list(mk("IVW", ivw@Estimate, ivw@CILower, ivw@CIUpper, ivw@Pvalue),
              mk("weighted median", wme@Estimate, wme@CILower, wme@CIUpper, wme@Pvalue))
    if (!is.null(egg)) r <- c(r, list(mk("MR-Egger", egg@Estimate, egg@CILower.Est,
                                         egg@CIUpper.Est, egg@Pvalue.Est,
                                         egg@Intercept, egg@Pvalue.Int)))
    rows[[tr]] <- rbindlist(r)
  }
} else if (opt$mode == "mvmr") {
  union_ivs <- Reduce(union, ivs)
  # Second-pass pruning of the union. run_mvmr.R does this; without it the union
  # keeps variants within the clump window of one another ACROSS traits, so the
  # instruments are correlated and the intervals come out too narrow.
  # Critically, variants inside a --cis_dense window are exempt: they were already
  # r^2-pruned at the locus and applying the genome-wide distance rule here would
  # collapse the locus to a single variant -- the very defect this script fixes.
  in_cis <- function(chr, pos) {
    for (sp in cis_spec)
      if (as.character(chr) == sp$chr && pos >= sp$start && pos <= sp$end) return(TRUE)
    FALSE
  }
  pos_ref <- rbindlist(lapply(EXPOSURES, function(tr)
    lipid[[tr]][ID %in% union_ivs, .(ID, CHROM, GENPOS, LOG10P)]))
  pos_ref <- pos_ref[pos_ref[, .I[which.max(LOG10P)], by = ID]$V1][order(-LOG10P)]
  cis_keep <- pos_ref[mapply(in_cis, CHROM, GENPOS), ID]
  poly     <- pos_ref[!mapply(in_cis, CHROM, GENPOS)]
  kept <- character(0); kp <- integer(0); kc <- character(0)
  for (i in seq_len(nrow(poly))) {
    ch <- as.character(poly$CHROM[i]); po <- as.integer(poly$GENPOS[i])
    sm <- kc == ch
    if (!length(kp) || !any(abs(kp[sm] - po) < opt$kb_clump * 1000L)) {
      kept <- c(kept, poly$ID[i]); kp <- c(kp, po); kc <- c(kc, ch)
    }
  }
  message(sprintf("[EXT] union %d -> %d (%d polygenic pruned at %g kb + %d cis-dense retained)",
                  length(union_ivs), length(kept) + length(cis_keep),
                  length(kept), opt$kb_clump, length(cis_keep)))
  union_ivs <- c(kept, cis_keep)
  h <- harmonize(EXPOSURES[1], union_ivs)[, .(SNP, oBETA, oSE)]
  bx <- bxse <- matrix(NA_real_, nrow(h), length(EXPOSURES),
                       dimnames = list(h$SNP, EXPOSURES))
  for (tr in EXPOSURES) {
    d <- lipid[[tr]][ID %in% h$SNP]
    i <- match(d$ID, h$SNP)
    bx[i, tr] <- d$BETA; bxse[i, tr] <- d$SE
  }
  ok <- complete.cases(bx) & complete.cases(bxse)
  bx <- bx[ok, , drop = FALSE]; bxse <- bxse[ok, , drop = FALSE]
  h <- h[ok]
  message(sprintf("[EXT] MVMR: %d instruments across %d exposures", nrow(bx), length(EXPOSURES)))
  fit <- mr_mvivw(mr_mvinput(bx = bx, bxse = bxse, by = h$oBETA, byse = h$oSE))
  rows[["mvmr"]] <- data.table(
    exposure = EXPOSURES, outcome = OUTLAB, n_instruments = nrow(bx), method = "MVMR-IVW",
    OR_perSD = exp(fit@Estimate), CI_lo95 = exp(fit@CILower), CI_hi95 = exp(fit@CIUpper),
    pval = fit@Pvalue, egger_intercept = NA_real_, egger_intercept_p = NA_real_,
    instrument_selection = vapply(EXPOSURES, function(t) how[[t]], character(1)))
} else stop("--mode must be 'univariable' or 'mvmr'")

res <- rbindlist(rows, fill = TRUE)
res[, design := sprintf("two-sample %s (UKB instruments -> %s)", opt$mode, OUTLAB)]
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, opt$out)
print(res)
message("[EXT] Wrote ", opt$out)
