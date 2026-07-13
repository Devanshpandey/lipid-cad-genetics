#!/usr/bin/env Rscript
# L=10 vs L=20 SuSiE sensitivity at the LDLR / LDL-C locus, with proper allele
# harmonization between REGENIE z-scores and plink LD, and a stabilized fit.
# Reviewer concern: is the "ten credible sets" result a hard L=10 ceiling artifact?
suppressMessages({library(data.table); library(susieR)})
set.seed(42)

BFILE  <- "/corral/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE_QC_400k/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"
PLINK  <- "/path/to/software/plink"
SUM    <- "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged/lipids_LDL_C.txt.gz"
CHR    <- 19L; BP <- 11202306L; W <- 500000L
OUT    <- "/path/to/cad-genetics/results/agent1_genetics/strengthen/susie_L_sensitivity"
tmp <- tempfile(); dir.create(tmp)

cat("[1] loading sumstats window\n")
d <- fread(cmd = paste("gunzip -c", shQuote(SUM)), na.strings = "NA")
w <- d[CHROM == CHR & GENPOS >= BP - W & GENPOS <= BP + W & !is.na(BETA) & !is.na(SE)]
cat(sprintf("    %d SNPs in +/-%dkb window\n", nrow(w), W/1000L))

snpf <- file.path(tmp, "snps.txt"); writeLines(w$ID, snpf)
ldp  <- file.path(tmp, "ld")
cat("[2] computing square LD (plink1.9)\n")
system(sprintf("%s --bfile %s --chr %d --extract %s --r square --write-snplist --out %s --threads 8 2>&1 | tail -2",
               PLINK, BFILE, CHR, snpf, ldp))
ld  <- as.matrix(fread(paste0(ldp, ".ld")))
ids <- readLines(paste0(ldp, ".snplist"))          # LD row/col order (bfile order)

cat("[3] allele harmonization (REGENIE ALLELE1 vs plink A1)\n")
bim <- fread(BFILE_BIM <- paste0(BFILE, ".bim"), header = FALSE)
setnames(bim, c("chr","ID","cm","bp","A1","A2"))
bim <- bim[match(ids, ID)]                          # A1/A2 in LD order
w   <- w[match(ids, ID)]                            # sumstats in LD order
stopifnot(nrow(w) == length(ids), nrow(bim) == length(ids))
z   <- w$BETA / w$SE
same <- w$ALLELE1 == bim$A1 & w$ALLELE0 == bim$A2   # z already w.r.t. A1
flip <- w$ALLELE1 == bim$A2 & w$ALLELE0 == bim$A1   # flip z sign
keep <- same | flip
z[flip] <- -z[flip]
cat(sprintf("    %d same-orient, %d flipped, %d dropped (allele mismatch)\n",
            sum(same), sum(flip), sum(!keep)))
w <- w[keep]; z <- z[keep]; ld <- ld[keep, keep]

# drop any SNP with NA in its LD row
na_rows <- apply(is.na(ld), 1, any)
if (any(na_rows)) { w <- w[!na_rows]; z <- z[!na_rows]; ld <- ld[!na_rows, !na_rows] }
ld <- (ld + t(ld)) / 2; diag(ld) <- 1
N  <- as.integer(median(w$N, na.rm = TRUE))
cat(sprintf("    %d SNPs into SuSiE; N=%d\n", nrow(w), N))

run <- function(L) {
  cat(sprintf("[4] SuSiE L=%d\n", L))
  fit <- tryCatch(
    susie_rss(z = z, R = ld, n = N, L = L, coverage = 0.95, min_abs_corr = 0.5,
              estimate_residual_variance = FALSE, check_prior = FALSE, max_iter = 500),
    error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(data.table(L=L, n_cs=NA, cs_id=NA, size=NA, top_snp=NA, top_pip=NA, purity=NA))
  cs <- fit$sets$cs
  if (is.null(cs) || length(cs) == 0)
    return(data.table(L=L, n_cs=0L, cs_id=NA, size=NA, top_snp=NA, top_pip=NA, purity=NA))
  pur <- fit$sets$purity
  rbindlist(lapply(seq_along(cs), function(k) {
    idx <- cs[[k]]; pip <- fit$pip[idx]; j <- idx[which.max(pip)]
    data.table(L=L, n_cs=length(cs), cs_id=names(cs)[k], size=length(idx),
               top_snp=w$ID[j], top_pip=round(max(pip),4),
               purity=round(pur[k,"min.abs.corr"],4))
  }))
}

res <- list()
for (L in c(10L, 15L, 20L)) {
  r <- run(L); res[[as.character(L)]] <- r
  fwrite(rbindlist(res, fill=TRUE), paste0(OUT, ".csv"))   # incremental
}
res <- rbindlist(res, fill = TRUE)
cat("\n===== RESULT =====\n"); print(res)
cat("\n--- credible sets per L (n_cs, distinct top variants) ---\n")
print(res[!is.na(top_snp), .(n_cs = .N, n_distinct_topsnp = uniqueN(top_snp)), by = L])
cat("DONE\n")
