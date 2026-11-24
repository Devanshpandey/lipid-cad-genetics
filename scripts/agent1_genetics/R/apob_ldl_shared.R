#!/usr/bin/env Rscript
# ============================================================================
# ApoB vs LDL-C per-variant genetic effects: shared variance and discordant loci.
#
# WHY. Supplementary Fig. S5b and Table S18 report "94% shared variance (r=0.97),
# 19 discordant loci" across 345 instruments. 345 is the PRE-bugfix instrument
# count: the distance-pruning guard compared misaligned vectors (a recycling bug)
# and let correlated variants through. This recomputes the statistic on the
# corrected instrument sets, using the SAME selection path as run_mrbma.R /
# run_mvmr.R post-fix: per-trait clump, union, prune the union again, then
# require completeness across the model's exposures and the outcome.
#
#   six  = the six MR-BMA exposures        (direct successor to the 345 set)
#   apob = LDL_C, APOB, TRIGLY             (the LDL-C/ApoB model)
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(optparse)})
set.seed(42)
opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat_dir", type = "character"),
  make_option("--outcome",     type = "character", default = "CAD"),
  make_option("--kb_clump",    type = "integer",   default = 10000L),
  make_option("--p_thresh",    type = "double",    default = 5e-8),
  make_option("--out",         type = "character")
)))

clump <- function(d, kb) {
  g <- d[LOG10P > -log10(opt$p_thresh) & !is.na(BETA) & !is.na(SE)][order(-LOG10P)]
  ks <- character(0); kp <- integer(0); kc <- character(0)
  for (i in seq_len(nrow(g))) {
    ch <- as.character(g$CHROM[i]); po <- as.integer(g$GENPOS[i]); sm <- kc == ch
    if (!length(kp) || !any(abs(kp[sm] - po) < kb * 1000L)) {
      ks <- c(ks, g$ID[i]); kp <- c(kp, po); kc <- c(kc, ch) }
  }
  ks
}
prune_union <- function(ivs, lip, TR) {
  ref <- rbindlist(lapply(TR, function(t) lip[[t]][ID %in% ivs, .(ID, CHROM, GENPOS, LOG10P)]))
  ref <- ref[ref[, .I[which.max(LOG10P)], by = ID]$V1][order(-LOG10P)]
  keep <- character(0); kp <- integer(0); kc <- character(0)
  for (i in seq_len(nrow(ref))) {
    ch <- as.character(ref$CHROM[i]); po <- as.integer(ref$GENPOS[i]); sm <- kc == ch
    if (!length(kp) || !any(abs(kp[sm] - po) < opt$kb_clump * 1000L)) {
      keep <- c(keep, ref$ID[i]); kp <- c(kp, po); kc <- c(kc, ch) }
  }
  keep
}

ALL <- c("LDL_C","APOB","HDL_C","APOA1","TRIGLY","LPA")
lip <- lapply(ALL, function(t) fread(cmd = paste("gunzip -c",
        shQuote(file.path(opt$sumstat_dir, sprintf("lipids_%s.txt.gz", t))))))
names(lip) <- ALL
od0 <- fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir,
        sprintf("outcomes_%s.txt.gz", opt$outcome)))))

report <- function(tag, TR) {
  ivs  <- unique(unlist(lapply(TR, function(t) clump(lip[[t]], opt$kb_clump))))
  keep <- prune_union(ivs, lip, TR)
  al <- lip[[TR[1]]][ID %in% keep, .(ID, A1 = ALLELE1, A0 = ALLELE0)]
  X <- matrix(NA_real_, length(keep), length(TR), dimnames = list(keep, TR))
  for (t in TR) {
    d <- merge(lip[[t]][ID %in% keep], al, by = "ID")
    X[match(d$ID, keep), t] <- fifelse(d$ALLELE1 == d$A1, d$BETA,
                                fifelse(d$ALLELE1 == d$A0, -d$BETA, NA_real_))
  }
  od <- merge(od0[ID %in% keep], al, by = "ID")
  by_ <- rep(NA_real_, length(keep)); i <- match(od$ID, keep)
  by_[i] <- fifelse(od$ALLELE1 == od$A1, od$BETA, fifelse(od$ALLELE1 == od$A0, -od$BETA, NA_real_))
  ok <- complete.cases(X) & !is.na(by_)
  X <- X[ok, , drop = FALSE]
  fit <- lm(X[, "APOB"] ~ X[, "LDL_C"]); r2 <- summary(fit)$r.squared
  sr <- rstandard(fit); disc <- sum(abs(sr) > 2)
  cat(sprintf("%-5s n=%d  r=%.3f  r2=%.4f  shared=%.0f%%  independent=%.0f%%  discordant(|z|>2)=%d\n",
              tag, nrow(X), sqrt(r2), r2, 100*r2, 100*(1-r2), disc))
  all <- data.table(ID = rownames(X), beta_LDL = X[,"LDL_C"], beta_APOB = X[,"APOB"],
                    std_resid = round(sr, 3))[order(-abs(std_resid))]
  list(sum = data.table(set = tag, n_instruments = nrow(X), r = round(sqrt(r2), 3),
                        r2 = round(r2, 4), shared_pct = round(100*r2, 1),
                        independent_pct = round(100*(1-r2), 1), n_discordant = disc),
       top = all[, set := tag])
}
a <- report("six",  ALL); b <- report("apob", c("LDL_C","APOB","TRIGLY"))
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(rbind(a$sum, b$sum), opt$out)
fwrite(rbind(a$top, b$top), sub("\\.csv$", "_discordant.csv", opt$out))
message("[shared] wrote ", opt$out)
