#!/usr/bin/env Rscript
# Agent 2 — Step 2: GTEx eQTL colocalization
# For each colocalized locus from Agent 1, runs coloc between GWAS
# summary stats and GTEx v8 eQTL data for key tissues.
#
# Tissues: Liver, Artery_Coronary, Artery_Aorta, Heart_Left_Ventricle,
#          Whole_Blood (proxy for macrophage)
#
# Usage:
#   Rscript 02_gtex_eqtl_coloc.R \
#     --credible_sets /path/to/finemapping_credible_sets.csv \
#     --sumstat_dir   /path/to/merged \
#     --eqtl_dir      /path/to/reference/eqtl_catalogue \
#     --out           /path/to/agent2_genes/gtex_eqtl_coloc.csv

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(optparse)
})

set.seed(42)

opt_list <- list(
  make_option("--credible_sets", type = "character"),
  make_option("--sumstat_dir",   type = "character"),
  make_option("--eqtl_dir",      type = "character"),
  make_option("--window_kb",     type = "integer", default = 500L),
  make_option("--pp4_thresh",    type = "double",  default = 0.5),
  make_option("--out",           type = "character", default = "gtex_eqtl_coloc.csv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

TISSUES <- c(
  "Liver",
  "Artery_Coronary",
  "Artery_Aorta",
  "Heart_Left_Ventricle",
  "Whole_Blood"
)

message("[eQTL coloc] Loading credible sets...")
cs <- fread(opt$credible_sets)
# Deduplicate: one unique locus per exposure
cs[, locus_key := paste(exposure, chr, round(locus_bp / 1e6), sep = "_")]
cs_unique <- cs[!duplicated(locus_key)]
message(sprintf("  %d unique exposure x locus combinations", nrow(cs_unique)))

results <- list()

for (i in seq_len(nrow(cs_unique))) {
  row      <- cs_unique[i]
  exposure <- row$exposure
  chr      <- row$chr
  bp       <- row$locus_bp
  sentinel <- row$top_snp
  w        <- opt$window_kb * 1000L

  message(sprintf("\n[eQTL coloc] Locus %d/%d: %s chr%s:%s",
                  i, nrow(cs_unique), exposure, chr, bp))

  # Load GWAS summary stats for this locus
  exp_file <- file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", exposure))
  if (!file.exists(exp_file)) {
    message("  GWAS file missing — skipping")
    next
  }

  gwas <- fread(cmd = paste("gunzip -c", shQuote(exp_file)), na.strings = "NA")
  gwas[, pval := 10^(-LOG10P)]
  gwas_win <- gwas[CHROM == chr & GENPOS >= bp - w & GENPOS <= bp + w &
                   !is.na(BETA) & !is.na(SE) & SE > 0]

  if (nrow(gwas_win) < 50) {
    message("  Too few GWAS SNPs in window — skipping")
    next
  }

  gwas_win[, varbeta := SE^2]
  gwas_n <- median(gwas_win$N, na.rm = TRUE)

  for (tissue in TISSUES) {
    eqtl_file <- file.path(opt$eqtl_dir, paste0(tissue, ".tsv.gz"))
    eqtl_tbi  <- paste0(eqtl_file, ".tbi")

    if (!file.exists(eqtl_file)) {
      message(sprintf("  %s eQTL file not found — skipping tissue", tissue))
      next
    }

    # Query eQTL file with tabix for the locus window
    region <- sprintf("%s:%d-%d", chr, max(1L, bp - w), bp + w)
    eqtl_raw <- tryCatch({
      fread(cmd = sprintf("tabix %s %s", shQuote(eqtl_file), region),
            header = FALSE, na.strings = "NA")
    }, error = function(e) {
      message(sprintf("  tabix failed for %s: %s", tissue, e$message))
      NULL
    })

    if (is.null(eqtl_raw) || nrow(eqtl_raw) == 0) {
      message(sprintf("  No eQTL data for %s in this region", tissue))
      next
    }

    # eQTL Catalogue column order (GTEx V8 harmonised):
    # molecular_trait_id, chromosome, position, ref, alt, variant,
    # ma_samples, maf, pvalue, beta, se, type, ac, an, r2, molecular_trait_object_id,
    # gene_id, median_tpm, rsid
    col_names <- c("gene_id","chromosome","position","ref","alt","variant",
                   "ma_samples","maf","pvalue","beta","se","type","ac","an",
                   "r2","molecular_trait_object_id","gene_name","median_tpm","rsid")
    if (ncol(eqtl_raw) >= length(col_names)) {
      setnames(eqtl_raw, seq_along(col_names), col_names)
    } else {
      # Fallback: minimal columns
      setnames(eqtl_raw, 1:min(ncol(eqtl_raw), 6), c("gene_id","chromosome","position","ref","alt","variant")[1:min(ncol(eqtl_raw),6)])
    }

    eqtl_raw <- eqtl_raw[!is.na(beta) & !is.na(se) & se > 0]
    if (nrow(eqtl_raw) == 0) next

    # Get unique genes tested at this locus
    genes_at_locus <- unique(eqtl_raw$gene_id)
    message(sprintf("  %s: %d genes, %d eQTL SNPs",
                    tissue, length(genes_at_locus), nrow(eqtl_raw)))

    for (gene in genes_at_locus) {
      eqtl_gene <- eqtl_raw[gene_id == gene]
      if (nrow(eqtl_gene) < 20) next

      # Match SNPs between GWAS and eQTL by position
      eqtl_gene[, pos_key := as.character(position)]
      gwas_win[,  pos_key := as.character(GENPOS)]

      common_pos <- intersect(eqtl_gene$pos_key, gwas_win$pos_key)
      if (length(common_pos) < 20) next

      eqtl_sub <- eqtl_gene[pos_key %in% common_pos]
      gwas_sub  <- gwas_win[pos_key  %in% common_pos]

      # Align by position (take first match per position)
      eqtl_sub <- eqtl_sub[!duplicated(pos_key)]
      gwas_sub  <- gwas_sub[!duplicated(pos_key)]
      setkey(eqtl_sub, pos_key)
      setkey(gwas_sub,  pos_key)
      merged <- merge(gwas_sub[, .(pos_key, BETA, SE, pval, varbeta)],
                      eqtl_sub[, .(pos_key, beta, se, pvalue, gene_id, gene_name)],
                      by = "pos_key")
      if (nrow(merged) < 20) next

      # Run coloc
      coloc_res <- tryCatch({
        coloc.abf(
          dataset1 = list(beta   = merged$BETA,
                          varbeta = merged$varbeta,
                          N      = gwas_n,
                          type   = "quant",
                          snp    = merged$pos_key),
          dataset2 = list(beta   = merged$beta,
                          varbeta = merged$se^2,
                          N      = nrow(merged),
                          type   = "quant",
                          snp    = merged$pos_key)
        )
      }, error = function(e) NULL)

      if (is.null(coloc_res)) next

      pp <- coloc_res$summary
      results[[length(results) + 1]] <- data.table(
        exposure   = exposure,
        locus_chr  = chr,
        locus_bp   = bp,
        sentinel   = sentinel,
        tissue     = tissue,
        gene_id    = gene,
        gene_name  = eqtl_gene$gene_name[1],
        n_snps     = nrow(merged),
        PP.H0      = pp["PP.H0.abf"],
        PP.H1      = pp["PP.H1.abf"],
        PP.H2      = pp["PP.H2.abf"],
        PP.H3      = pp["PP.H3.abf"],
        PP.H4      = pp["PP.H4.abf"],
        eqtl_coloc_sig = pp["PP.H4.abf"] >= opt$pp4_thresh
      )
    }
  }
}

if (length(results) > 0) {
  out_dt <- rbindlist(results)
  setorder(out_dt, -PP.H4)
  fwrite(out_dt, opt$out, sep = ",")
  n_sig <- sum(out_dt$eqtl_coloc_sig, na.rm = TRUE)
  message(sprintf("\n[eQTL coloc] %d gene-tissue pairs tested, %d with PP.H4 >= %.1f",
                  nrow(out_dt), n_sig, opt$pp4_thresh))
  message(sprintf("[eQTL coloc] Results written: %s", opt$out))

  # Print top hits
  top <- head(out_dt[eqtl_coloc_sig == TRUE][order(-PP.H4)], 20)
  if (nrow(top) > 0) {
    message("\nTop eQTL colocalized genes:")
    for (j in seq_len(nrow(top))) {
      message(sprintf("  %s | %s | %s | PP.H4=%.3f",
                      top$gene_name[j], top$tissue[j],
                      top$exposure[j],  top$PP.H4[j]))
    }
  }
} else {
  message("[eQTL coloc] No results produced")
}

message("[eQTL coloc] Complete.")
