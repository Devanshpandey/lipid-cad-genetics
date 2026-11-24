#!/usr/bin/env Rscript
# cis-QTL colocalization: gene expression (or plasma protein) vs GWAS traits.
#
# Replaces the uncommitted TACC heredoc that produced Table S27 (which covered
# 10 genes x 3 tissues x LDL-C only) and supersedes
# scripts/agent2_genes/R/02_gtex_eqtl_coloc.R, whose hardcoded 19-name column
# vector does not match the 18-column eQTL Catalogue GTEx_V8/ge format: the
# length test fails, it drops into a 6-column fallback, and then references
# beta/se columns that do not exist. Column names here come from the file's own
# header, never from a positional vector.
#
# Two additions that matter for interpretation:
#   * every gene is tested against the OUTCOME (CAD) as well as the lipid
#     exposures, so "no colocalization" can be attributed to an absent outcome
#     signal (H1) rather than to distinct causal variants (H3);
#   * PP.H0-H4 are all reported, not just H3/H4.
#
# Usage:
#   Rscript run_qtl_coloc.R --genes genes.tsv --qtl_dir <extracts> \
#     --tissues Liver,Artery_Aorta --sumstat_dir <merged> \
#     --traits lipids_LDL_C,outcomes_CAD --out qtl_coloc.csv

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(optparse)
})
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--genes",       type = "character"),
  make_option("--qtl_dir",     type = "character"),
  make_option("--tissues",     type = "character"),
  make_option("--sumstat_dir", type = "character"),
  make_option("--traits",      type = "character"),
  make_option("--window_kb",   type = "integer", default = 500L),
  make_option("--min_snps",    type = "integer", default = 30L),
  make_option("--out",         type = "character")
)))
for (r in c("genes","qtl_dir","tissues","sumstat_dir","traits","out"))
  if (is.null(opt[[r]]) || !nzchar(opt[[r]])) stop(sprintf("--%s is required", r))

genes   <- fread(opt$genes, header = FALSE, col.names = c("sym","ens","chr"))
tissues <- strsplit(opt$tissues, ",")[[1]]
traits  <- strsplit(opt$traits,  ",")[[1]]
W       <- opt$window_kb * 1000L

message(sprintf("[qtl-coloc] %d genes x %d tissues x %d traits",
                nrow(genes), length(tissues), length(traits)))

# ---- load QTL extracts (header-driven column names) -------------------------
qtl <- list()
for (t in tissues) {
  f <- file.path(opt$qtl_dir, sprintf("%s_qtl.tsv.gz", t))
  if (!file.exists(f)) { message("  missing extract: ", f); next }
  d <- fread(cmd = sprintf("zcat %s", shQuote(f)))
  need <- c("gene_id","rsid","beta","se","maf","an","position")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(sprintf("%s extract lacks columns: %s", t, paste(miss, collapse=", ")))
  qtl[[t]] <- d[!is.na(beta) & !is.na(se) & se > 0 & maf > 0]
  message(sprintf("  %-22s %8d rows, %d genes", t, nrow(qtl[[t]]), uniqueN(qtl[[t]]$gene_id)))
}

res <- list()

for (tr in traits) {
  f <- file.path(opt$sumstat_dir, sprintf("%s.txt.gz", tr))
  if (!file.exists(f)) { message("  missing sumstat: ", f); next }
  is_cc <- grepl("^outcomes_", tr)
  message(sprintf("[qtl-coloc] trait %s (type=%s)", tr, if (is_cc) "cc" else "quant"))
  G <- fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA",
             select = c("CHROM","GENPOS","ID","A1FREQ","N","BETA","SE"))
  setnames(G, c("CHR","BP","rsid","freq","N","gb","gse"))
  G <- G[!is.na(gb) & !is.na(gse) & gse > 0]

  for (i in seq_len(nrow(genes))) {
    g <- genes[i]
    for (t in names(qtl)) {
      e <- qtl[[t]][gene_id == g$ens]
      if (nrow(e) < opt$min_snps) next
      e <- e[!duplicated(rsid)]
      # Restrict the GWAS by CHROMOSOME only, then merge on rsID.
      # Do NOT filter on eQTL `position`: the eQTL Catalogue is GRCh38 and these
      # GWAS are GRCh37 (~466 kb apart at PCSK9), so a position window taken from
      # one build silently drops real variants in the other. rsID is build-agnostic,
      # and the extract already contains only this gene's cis variants, so the
      # merge itself is what bounds the window.
      gw <- G[CHR == as.integer(g$chr)]
      if (nrow(gw) < opt$min_snps) next
      gw <- gw[!duplicated(rsid)]

      m <- merge(e[, .(rsid, eb = beta, ese = se, emaf = pmin(maf, 1 - maf), ean = an)],
                 gw[, .(rsid, gb, gse, gmaf = pmin(freq, 1 - freq), gn = N)], by = "rsid")
      m <- m[emaf > 0 & gmaf > 0]
      if (nrow(m) < opt$min_snps) next

      d_gwas <- list(beta = m$gb, varbeta = m$gse^2, snp = m$rsid,
                     MAF = m$gmaf, N = max(m$gn),
                     type = if (is_cc) "cc" else "quant")
      d_qtl  <- list(beta = m$eb, varbeta = m$ese^2, snp = m$rsid,
                     MAF = m$emaf, N = round(max(m$ean) / 2), type = "quant")

      cr <- tryCatch(suppressWarnings(coloc.abf(d_gwas, d_qtl)),
                     error = function(er) NULL)
      if (is.null(cr)) next
      s <- cr$summary
      res[[paste(tr, g$sym, t)]] <- data.table(
        trait = tr, gene = g$sym, ensembl = g$ens, tissue = t, n_snp = nrow(m),
        PP.H0 = round(s["PP.H0.abf"],3), PP.H1 = round(s["PP.H1.abf"],3),
        PP.H2 = round(s["PP.H2.abf"],3), PP.H3 = round(s["PP.H3.abf"],3),
        PP.H4 = round(s["PP.H4.abf"],3))
    }
  }
  rm(G); invisible(gc())
}

out <- rbindlist(res)
if (nrow(out) == 0) stop("No colocalization tests completed — check gene IDs and extracts.")
setorder(out, -PP.H4)
fwrite(out, opt$out)
message(sprintf("[qtl-coloc] %d tests written to %s", nrow(out), opt$out))
message(sprintf("[qtl-coloc] PP.H4 > 0.8: %d", nrow(out[PP.H4 > 0.8])))
print(head(out, 25))
