#!/usr/bin/env Rscript
# Nature round r1, Fig. 4a: multi-signal colocalization (coloc.susie) at all five loci of panel a.
#
# Why: the published Supplementary Table (tab:oS16, coloc_susie) covers PCSK9 and LDLR only, so panel a
# printed "n/a" at SORT1/CELSR2, APOE/APOC1 and TRIB1. Those three were never actually tested. The first
# five-locus attempt (strengthen/coloc_susie.R, job 3288056, 9 Jul 2026) returned no credible sets anywhere
# because it called runsusie(repeat_until_convergence=FALSE) inside tryCatch, so every fit that needed more
# than 100 iterations became "no signal"; it was cancelled and the working rerun (strengthen/coloc_lean.R,
# 10 Jul) was restricted to two loci. This script is coloc_lean.R's protocol, run at all five.
#
# Protocol, fixed before any new locus was run (8 Nov 2026):
#  * UK Biobank LDL-C (merged/, the lineage behind the published table) against UK Biobank CAD;
#    window +-120 kb; in-sample LD from the 402,199-person UK Biobank bed (PLINK 1.9 --r square);
#    coloc::runsusie defaults (L = 10, repeat until convergence); CAD case fraction s = 0.123;
#    coloc.susie default priors. Per locus: credible sets per trait, every credible-set pair, the maximum
#    PP.H4 over pairs and the number of pairs above 0.8.
#  * Window centres: PCSK9 and LDLR exactly as published (55,505,647 and 11,197,598). SORT1/CELSR2,
#    APOE/APOC1 and TRIB1 are centred on the fine-mapped lead variant that Fig. 4a prints for the locus
#    (rs140584594, rs429358, rs28601761), looked up in the LDL-C statistics.
#  * Reproduction control: PCSK9 and LDLR run first and must return the published table exactly (SNPs
#    1,411 / 1,486; LDL-C signals 9 / 10; CAD signals 1 / 1; pairs 1 / 2; max PP.H4 1.000 / 0.936).
#    Otherwise exit status 2 and no new locus is run.
#  * Alignment, which coloc_lean.R assumed without checking: each window SNP must occur once on the
#    chromosome in the bim, the LD matrix is labelled from the bim rows PLINK read, and ALLELE1/ALLELE0
#    must equal bim A1/A2 (REGENIE's effect allele on bed input). A failed check is corrected
#    (reorder; flip the sign of beta) and counted in the output, never silently truncated.
#  * Sensitivity, reported in the table but not plotted: TRIB1 with triglycerides as the exposure (the
#    lipid FinnGen's TRIB1 replication used).
#  * Figure rule: if every primary locus returns a coloc.susie result (at least one credible set for each
#    trait), panel a plots max PP.H4 at all five; otherwise the column is removed from panel a and the
#    table reports every locus.
#  * Diagnostic: susieR::estimate_s_rss for each trait (consistency of z-scores with the LD matrix;
#    0 is perfect, larger is worse). Reported, not used to exclude.
#
# Amendment 1 (8 Nov 2026, after run 1 = job 3565988, whose outputs are kept in run1_prespecified/):
#  the centre rule assumed the lead variant printed in Fig. 4a lies at its locus. For SORT1/CELSR2 it does not:
#  rs140584594 is at chr1:110,232,983 (GRCh37, in GSTM1), 415 kb from rs12740374, so run 1's window held no
#  SORT1 signal (max LDL-C -log10P 32, CAD 3.0, no CAD credible set) at a locus whose UKB LDL-C-CAD coloc.abf
#  PP.H4 is 0.92. Rule added, for every locus: the window must contain the strongest LDL-C variant within
#  +-500 kb of the panel's UK Biobank sentinel (the Table S4 LDL-C -> CAD row that coloc.abf tested); a new
#  locus whose window does not is re-centred on that variant. Controls are only checked (their centres are the
#  published ones). Written after run 1 showed SORT1 with no CAD credible set; no SORT1 PP.H4 had been seen.
#
# Usage: Rscript nat_coloc_susie_5loci.R OUTDIR [THREADS]
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc)})
set.seed(42)
a <- commandArgs(TRUE)
OUT <- if (length(a) >= 1) a[1] else stop("usage: nat_coloc_susie_5loci.R OUTDIR [THREADS]")
THREADS <- if (length(a) >= 2) as.integer(a[2]) else 8L
REPO <- "/path/to/cad-genetics"
M <- file.path(REPO, "results/agent1_genetics/gwas_summary_stats/merged")
LDBF <- "/path/to/ukb_genotypes/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"
PLINK1 <- "/path/to/software/plink"
S_CAD <- 0.123; W <- 1.2e5
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("R", as.character(getRversion()), "susieR", as.character(packageVersion("susieR")),
    "coloc", as.character(packageVersion("coloc")), "data.table", as.character(packageVersion("data.table")), "\n")

loci <- data.table(
  locus      = c("PCSK9", "LDLR", "SORT1/CELSR2", "APOE/APOC1", "TRIB1", "TRIB1"),
  exposure   = c("LDL_C", "LDL_C", "LDL_C", "LDL_C", "LDL_C", "TRIGLY"),
  role       = c("control", "control", "primary", "primary", "primary", "sensitivity"),
  chr        = c(1L, 19L, 1L, 19L, 8L, 8L),
  centre_snp = c(NA, NA, "rs140584594", "rs429358", "rs28601761", "rs28601761"),
  centre_bp  = c(55505647, 11197598, NA, NA, NA, NA),
  # the Table S4 LDL-C -> CAD sentinel behind the panel's UK Biobank column (Amendment 1)
  sentinel_snp = c("rs11591147", "rs10423733", "rs4970834", "rs531660643", "rs28601761", "rs28601761"))
PUBLISHED <- list(PCSK9 = c(nsnp = 1411, n1 = 9, n2 = 1, npair = 1, maxpp = 1.000),
                  LDLR  = c(nsnp = 1486, n1 = 10, n2 = 1, npair = 2, maxpp = 0.936))

# ---------------------------------------------------------------- statistics for the three chromosomes
rd_chr <- function(trait) {
  f <- file.path(M, sprintf("%s.txt.gz", trait))
  stopifnot(file.exists(f))
  dt <- fread(cmd = sprintf("gunzip -c %s | awk 'NR==1 || $1==1 || $1==8 || $1==19'", shQuote(f)), na.strings = "NA")
  dt[!is.na(BETA) & !is.na(SE)]
}
ST <- list(LDL_C = rd_chr("lipids_LDL_C"), TRIGLY = rd_chr("lipids_TRIGLY"), CAD = rd_chr("outcomes_CAD"))
for (k in names(ST)) cat(k, nrow(ST[[k]]), "rows on chr 1/8/19\n")
for (i in which(is.na(loci$centre_bp))) {
  hit <- ST$LDL_C[ID == loci$centre_snp[i]]
  if (nrow(hit) != 1 || hit$CHROM != loci$chr[i]) stop("centre SNP ", loci$centre_snp[i], " not found once on chr ", loci$chr[i])
  set(loci, i, "centre_bp", as.numeric(hit$GENPOS))
}
# Amendment 1: the window must contain the strongest LDL-C variant within +-500 kb of the sentinel
loci[, `:=`(centre_snp_prespecified = centre_snp, centre_bp_prespecified = centre_bp, lead_ldl_snp = NA_character_,
            lead_ldl_bp = NA_real_, window_contains_lead = NA, recentred = FALSE)]
for (i in seq_len(nrow(loci))) {
  sb <- ST$LDL_C[ID == loci$sentinel_snp[i] & CHROM == loci$chr[i]]$GENPOS
  if (length(sb) != 1) stop("sentinel ", loci$sentinel_snp[i], " not found once")
  reg <- ST$LDL_C[CHROM == loci$chr[i] & abs(GENPOS - sb) <= 5e5]
  ld1 <- reg[which.max(LOG10P)]
  inside <- abs(ld1$GENPOS - loci$centre_bp[i]) <= W
  set(loci, i, c("lead_ldl_snp", "lead_ldl_bp", "window_contains_lead"), list(ld1$ID, as.numeric(ld1$GENPOS), inside))
  if (!inside && loci$role[i] != "control")
    set(loci, i, c("centre_snp", "centre_bp", "recentred"), list(ld1$ID, as.numeric(ld1$GENPOS), TRUE))
}
print(loci)

win <- function(dt, chr, bp) dt[CHROM == chr & GENPOS >= bp - W & GENPOS <= bp + W][order(GENPOS)]

run_locus <- function(L) {
  tag <- gsub("[^A-Za-z0-9]", "_", paste(L$locus, L$exposure, sep = "_"))
  message("== ", L$locus, " (", L$exposure, ", ", L$role, ") ==")
  e <- win(ST[[L$exposure]], L$chr, L$centre_bp); o <- win(ST$CAD, L$chr, L$centre_bp)
  stopifnot(!anyDuplicated(e$ID), !anyDuplicated(o$ID))
  com <- intersect(e$ID, o$ID); e <- e[ID %in% com][order(GENPOS)]; o <- o[ID %in% com][order(GENPOS)]
  stopifnot(identical(e$ID, o$ID), all(e$ALLELE0 == o$ALLELE0 & e$ALLELE1 == o$ALLELE1))
  td <- tempfile(); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  writeLines(e$ID, file.path(td, "s.txt"))
  # the bim rows PLINK will read, in its order (all rows on the chromosome with a requested ID)
  bim <- fread(cmd = sprintf("awk 'NR==FNR{a[$1];next} ($1==%d) && ($2 in a)' %s %s", L$chr,
                             shQuote(file.path(td, "s.txt")), shQuote(paste0(LDBF, ".bim"))),
               header = FALSE, col.names = c("chr", "ID", "cm", "bp", "A1", "A2"))
  system(sprintf("%s --bfile %s --chr %d --extract %s --r square --out %s --threads %d >/dev/null 2>&1",
                 PLINK1, LDBF, L$chr, file.path(td, "s.txt"), file.path(td, "ld"), THREADS))
  if (!file.exists(file.path(td, "ld.ld"))) stop("PLINK wrote no LD matrix for ", L$locus)
  ld <- as.matrix(fread(file.path(td, "ld.ld"), header = FALSE))
  stopifnot(nrow(ld) == nrow(bim), ncol(ld) == nrow(bim))
  n_missing_bim <- sum(!e$ID %in% bim$ID)
  dup <- bim$ID[duplicated(bim$ID)]
  lean_exact <- n_missing_bim == 0 && !length(dup) && identical(bim$ID, e$ID)
  keep_b <- !bim$ID %in% dup
  bim <- bim[keep_b]; ld <- ld[keep_b, keep_b, drop = FALSE]
  e <- e[match(bim$ID, ID)]; o <- o[match(bim$ID, ID)]
  stopifnot(identical(e$ID, bim$ID), identical(o$ID, bim$ID))
  same <- e$ALLELE1 == bim$A1 & e$ALLELE0 == bim$A2
  swap <- e$ALLELE1 == bim$A2 & e$ALLELE0 == bim$A1
  n_flip <- sum(swap & !same); n_allele_drop <- sum(!same & !swap)
  e[swap & !same, BETA := -BETA]; o[swap & !same, BETA := -BETA]
  ok <- same | swap
  e <- e[ok]; o <- o[ok]; ld <- ld[ok, ok, drop = FALSE]
  bad <- apply(is.na(ld), 1, any)
  if (any(bad)) { e <- e[!bad]; o <- o[!bad]; ld <- ld[!bad, !bad, drop = FALSE] }
  lean_exact <- lean_exact && n_flip == 0 && n_allele_drop == 0
  ld <- (ld + t(ld)) / 2; diag(ld) <- 1; rownames(ld) <- colnames(ld) <- e$ID
  D1 <- list(beta = e$BETA, varbeta = e$SE^2, snp = e$ID, position = e$GENPOS, type = "quant",
             N = as.integer(median(e$N)), sdY = 1, LD = ld)
  D2 <- list(beta = o$BETA, varbeta = o$SE^2, snp = o$ID, position = o$GENPOS, type = "cc",
             N = as.integer(median(o$N)), s = S_CAD, LD = ld)
  fit <- function(D) tryCatch(runsusie(D), error = function(x) { message("  runsusie ERROR: ", conditionMessage(x)); NULL })
  s1 <- fit(D1); s2 <- fit(D2)
  if (is.null(s1) || is.null(s2)) stop("runsusie failed at ", L$locus, "; refusing to report it as 'no signal'")
  n1 <- length(s1$sets$cs); n2 <- length(s2$sets$cs)
  srss <- function(D) tryCatch(estimate_s_rss(D$beta / sqrt(D$varbeta), D$LD, D$N), error = function(x) NA_real_)
  cs_tab <- function(s, trait) {
    if (!length(s$sets$cs)) return(NULL)
    rbindlist(lapply(seq_along(s$sets$cs), function(k) {
      idx <- s$sets$cs[[k]]; top <- idx[which.max(s$pip[idx])]
      data.table(locus = L$locus, exposure = L$exposure, trait = trait, cs = names(s$sets$cs)[k],
                 size = length(idx), top_snp = e$ID[top], top_pip = round(s$pip[top], 4),
                 min_abs_corr = round(s$sets$purity[k, "min.abs.corr"], 4))
    }))
  }
  cst <- rbind(cs_tab(s1, L$exposure), cs_tab(s2, "CAD"))
  maxpp <- NA_real_; npair <- 0L; hit1 <- hit2 <- NA_character_; pairs <- NULL
  if (n1 > 0 && n2 > 0) {
    cs <- coloc.susie(s1, s2)
    if (!is.null(cs$summary) && nrow(cs$summary) > 0) {
      pairs <- data.table(locus = L$locus, exposure = L$exposure, as.data.table(cs$summary))
      j <- which.max(pairs$PP.H4.abf)
      maxpp <- pairs$PP.H4.abf[j]; npair <- sum(pairs$PP.H4.abf > 0.8)
      hit1 <- pairs$hit1[j]; hit2 <- pairs$hit2[j]
    }
  }
  res <- data.table(locus = L$locus, exposure = L$exposure, role = L$role, chr = L$chr,
                    centre_snp = ifelse(is.na(L$centre_snp), "", L$centre_snp), centre_bp = L$centre_bp,
                    sentinel_snp = L$sentinel_snp, lead_ldl_snp = L$lead_ldl_snp,
                    window_contains_lead = L$window_contains_lead, recentred = L$recentred,
                    window_kb = W / 1000, nsnp = nrow(e), n_signals_exposure = n1, n_signals_CAD = n2,
                    n_pairs_tested = if (is.null(pairs)) 0L else nrow(pairs),
                    n_coloc_pairs_PP4gt0.8 = npair, max_PP.H4 = round(maxpp, 3),
                    hit_exposure_at_max = hit1, hit_CAD_at_max = hit2,
                    max_log10P_exposure = round(max(e$LOG10P), 1), max_log10P_CAD = round(max(o$LOG10P), 1),
                    s_rss_exposure = round(srss(D1), 4), s_rss_CAD = round(srss(D2), 4),
                    n_missing_from_bim = n_missing_bim, n_duplicate_bim_ids = length(dup),
                    n_allele_flipped = n_flip, n_allele_dropped = n_allele_drop, n_dropped_ld_na = sum(bad),
                    alignment_as_coloc_lean = lean_exact)
  print(res)
  list(res = res, pairs = pairs, cs = cst)
}

out <- list()
for (i in which(loci$role == "control")) out[[length(out) + 1]] <- run_locus(loci[i])
R0 <- rbindlist(lapply(out, `[[`, "res"))
ctrl_ok <- all(vapply(names(PUBLISHED), function(g) {
  r <- R0[locus == g]; p <- PUBLISHED[[g]]
  nrow(r) == 1 && r$nsnp == p[["nsnp"]] && r$n_signals_exposure == p[["n1"]] && r$n_signals_CAD == p[["n2"]] &&
    r$n_coloc_pairs_PP4gt0.8 == p[["npair"]] && isTRUE(all.equal(r$max_PP.H4, p[["maxpp"]], tolerance = 1e-9))
}, logical(1)))
verdict <- c(sprintf("reproduction control (PCSK9, LDLR vs published table): %s", if (ctrl_ok) "PASS" else "FAIL"))
if (!ctrl_ok) {
  fwrite(R0, file.path(OUT, "coloc_susie5_CONTROL_FAILED.csv"))
  writeLines(verdict, file.path(OUT, "VERDICT.txt")); cat(verdict, sep = "\n"); quit(status = 2)
}
for (i in which(loci$role != "control")) out[[length(out) + 1]] <- run_locus(loci[i])
R <- rbindlist(lapply(out, `[[`, "res"))
P <- rbindlist(lapply(out, `[[`, "pairs"), fill = TRUE)
C <- rbindlist(lapply(out, `[[`, "cs"), fill = TRUE)
prim <- R[role %in% c("control", "primary")]
all_tested <- all(prim$n_signals_exposure > 0 & prim$n_signals_CAD > 0 & !is.na(prim$max_PP.H4))
verdict <- c(verdict,
  sprintf("alignment identical to coloc_lean.R at every locus: %s", all(R$alignment_as_coloc_lean)),
  sprintf("Amendment 1: windows re-centred on the lead LDL-C variant: %s",
          if (any(R$recentred)) paste(unique(R[recentred == TRUE, paste0(locus, " -> ", centre_snp)]), collapse = "; ") else "none"),
  sprintf("every window contains its locus's lead LDL-C variant (controls included): %s", all(R$window_contains_lead | R$recentred)),
  sprintf("every primary locus returned a coloc.susie result: %s -> %s", all_tested,
          if (all_tested) "panel a plots max PP.H4 at all five loci" else "remove the coloc.susie column from panel a"),
  sprintf("  %-13s %-6s max PP.H4 %s (%d of %d pairs >0.8; %d %s / %d CAD sets)", R$locus, R$exposure,
          format(R$max_PP.H4, nsmall = 3), R$n_coloc_pairs_PP4gt0.8,
          R$n_pairs_tested, R$n_signals_exposure, R$exposure, R$n_signals_CAD))
fwrite(R, file.path(OUT, "coloc_susie5_summary.csv"))
fwrite(P, file.path(OUT, "coloc_susie5_pairs.csv"))
fwrite(C, file.path(OUT, "coloc_susie5_credible_sets.csv"))
writeLines(verdict, file.path(OUT, "VERDICT.txt"))
cat("\n=== VERDICT ===\n"); cat(verdict, sep = "\n")
cat("ALL_DONE\n")
