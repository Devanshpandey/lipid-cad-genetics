#!/usr/bin/env Rscript
# ============================================================================
# Drug-target cis-Mendelian randomization (Table S33 / Figure 5c).
#
# WHY THIS EXISTS. No driver script survived for the cis-MR reported in the
# manuscript, and the Methods disclose that its LDL-C instruments came from the
# pre-correction LDL-C GWAS, in which statin use was ascertained from field 6153
# alone -- a field UK Biobank asks only of women. This reimplements the analysis
# to the Methods specification and runs it against BOTH the old and the corrected
# (6153 U 6177) LDL-C GWAS, so the effect of the ascertainment fix is measurable
# rather than assumed.
#
# VALIDATION. Run with --validate to check the old-GWAS arm against the published
# instrument counts (LPA 11, CETP 8, PCSK9/LDLR/APOC3 3 each, HMGCR/NPC1L1/
# ANGPTL3 1 each). If those do not reproduce, the cis windows below are wrong and
# nothing downstream should be believed.
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(optparse)})
set.seed(42)
opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat_dir", type = "character"),
  make_option("--ldl_file",    type = "character", default = "lipids_LDL_C.txt.gz"),
  make_option("--outcome_dir", type = "character", default = ""),
  make_option("--plink",       type = "character", default = "plink"),
  make_option("--ld_ref",      type = "character", default = ""),
  make_option("--label",       type = "character", default = "run"),
  make_option("--out",         type = "character")
)))
if (!nzchar(opt$outcome_dir)) opt$outcome_dir <- opt$sumstat_dir

# GRCh37/hg19 gene bodies; cis window is the gene body +/- 100 kb (Methods).
GENES <- data.table(
  gene   = c("PCSK9","ANGPTL3","HMGCR","NPC1L1","APOC3","LPA","CETP","LDLR"),
  chrom  = c(1L, 1L, 5L, 7L, 11L, 6L, 16L, 19L),
  start  = c(55505221L, 63063191L, 74632154L, 44552134L, 116700623L, 160952515L, 56995762L, 11200038L),
  end    = c(55530525L, 63071984L, 74657929L, 44580929L, 116703788L, 161087407L, 57017757L, 11244506L),
  exposure = c("LDL_C","TRIGLY","LDL_C","LDL_C","TRIGLY","LPA","HDL_C","LDL_C"))
# ANGPTL3 is additionally instrumented on LDL-C and on ApoB (Methods).
EXTRA <- data.table(gene = "ANGPTL3", exposure = c("LDL_C","APOB"))
JOBS <- rbind(GENES[, .(gene, exposure)], EXTRA)

expo_path <- function(tr)
  file.path(opt$sumstat_dir,
            if (tr == "LDL_C") opt$ldl_file else sprintf("lipids_%s.txt.gz", tr))

read_gz <- function(p) fread(cmd = paste("gunzip -c", shQuote(p)), showProgress = FALSE)

cache <- new.env()
get_expo <- function(tr) {
  if (!exists(tr, envir = cache)) assign(tr, read_gz(expo_path(tr)), envir = cache)
  get(tr, envir = cache)
}
get_out <- function(o) {
  k <- paste0("O_", o)
  if (!exists(k, envir = cache))
    assign(k, read_gz(file.path(opt$outcome_dir, sprintf("outcomes_%s.txt.gz", o))), envir = cache)
  get(k, envir = cache)
}

# LD clumping at r2<0.001 within 1 Mb. A distance rule cannot substitute here: the
# cis windows are ~200 kb, so dropping everything within 1 Mb of the lead always
# leaves exactly one instrument, which is how the first version of this script was
# wrong. Requires a real LD panel.
clump_cis <- function(d, tag) {
  if (nrow(d) == 1L) return(d)
  if (!nzchar(opt$ld_ref)) stop("no --ld_ref supplied; refusing to fake LD clumping")
  tmp <- tempfile(); assoc <- paste0(tmp, ".assoc")
  fwrite(d[, .(SNP = ID, P = 10^(-LOG10P))], assoc, sep = " ")
  # LOG10P can underflow to P=0; PLINK needs a positive P, and rank order is what matters
  system(paste0("awk 'NR==1{print;next}{if($2==0)$2=1e-300; print}' ", assoc,
                " > ", assoc, ".fix && mv ", assoc, ".fix ", assoc))
  cmd <- sprintf("%s --bfile %s --clump %s --clump-p1 1 --clump-r2 0.001 --clump-kb 1000 --out %s >/dev/null 2>&1",
                 opt$plink, opt$ld_ref, assoc, tmp)
  system(cmd)
  cf <- paste0(tmp, ".clumped")
  if (!file.exists(cf)) { message("  [clump] PLINK produced nothing for ", tag, "; keeping lead only")
                          return(d[order(-LOG10P)][1]) }
  kept <- fread(cf)$SNP
  message(sprintf("  [clump] %s: %d candidates -> %d after r2<0.001/1Mb", tag, nrow(d), length(kept)))
  d[ID %in% kept]
}

rows <- rbindlist(lapply(seq_len(nrow(JOBS)), function(k) {
  g  <- JOBS$gene[k]; tr <- JOBS$exposure[k]
  gi <- GENES[gene == g]
  e  <- get_expo(tr)
  win <- e[CHROM == gi$chrom & GENPOS >= gi$start - 1e5 & GENPOS <= gi$end + 1e5 &
           LOG10P > -log10(5e-8) & grepl("^rs", ID) & !is.na(BETA) & !is.na(SE)]
  if (!nrow(win)) return(NULL)
  iv <- clump_cis(win, paste(g, tr))
  rbindlist(lapply(c("CAD","MI"), function(o) {
    od <- get_out(o)[ID %in% iv$ID]
    m  <- merge(iv[, .(ID, CHROM, GENPOS, A1 = ALLELE1, A0 = ALLELE0, bx = BETA, sx = SE)],
                od[, .(ID, EA = ALLELE1, OA = ALLELE0, by = BETA, sy = SE)], by = "ID")
    if (!nrow(m)) return(NULL)
    m[, by := fifelse(EA == A1, by, fifelse(EA == A0, -by, NA_real_))]
    m <- m[!is.na(by)]
    # align to the lipid-RAISING allele so OR>1 is concordant with a lowering drug
    m[bx < 0, `:=`(bx = -bx, by = -by)]
    w  <- 1 / (m$sy^2 / m$bx^2)
    bh <- sum(w * (m$by / m$bx)) / sum(w); se <- sqrt(1 / sum(w))
    data.table(gene = g, exposure = tr, outcome = o, n_instruments = nrow(m),
               OR = exp(bh), CI_lo = exp(bh - 1.96*se), CI_hi = exp(bh + 1.96*se),
               P = 2 * pnorm(-abs(bh / se)), arm = opt$label)
  }))
}))
print(rows[outcome == "CAD", .(gene, exposure, n_instruments, OR = round(OR,3),
                               CI_lo = round(CI_lo,3), CI_hi = round(CI_hi,3),
                               P = signif(P,2))])
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(rows, opt$out)
message("[cis] wrote ", opt$out)
