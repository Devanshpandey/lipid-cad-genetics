#!/usr/bin/env Rscript
# ============================================================
# Compare the censored-block Lp(a) GWAS with the published one.
#
# Control gate: LPA_MEAS (our ranks, no --apply-rint, step 1 on scores rather than raw nmol/L) must reproduce the
# published BETAs at LPA-locus GWS variants. If it does not, any difference in the censored arm cannot be
# attributed to the censored participants, so the script stops non-zero and the chain stops with it.
# Writes summary statistics only.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--archived", type = "character", help = "merged/lipids_LPA.txt.gz"),
  make_option("--meas",     type = "character", help = "LPA_MEAS sumstats"),
  make_option("--cens",     type = "character", help = "LPA_CENS sumstats"),
  make_option("--min_r",    type = "double", default = 0.995),
  make_option("--slope_tol", type = "double", default = 0.02),
  make_option("--out",      type = "character", help = "Output prefix")
)))
cols <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P")
rd <- function(f) fread(cmd = paste("gunzip -c", shQuote(f)), select = cols, na.strings = "NA")
meas <- rd(opt$meas); cens <- rd(opt$cens); arch <- rd(opt$archived)
arch <- arch[CHROM %in% unique(meas$CHROM)]
message(sprintf("[cmp] variants: archived %d, meas %d, cens %d, chromosomes %s",
                nrow(arch), nrow(meas), nrow(cens), paste(sort(unique(meas$CHROM)), collapse = ",")))

sfx <- function(d, s) { d <- copy(d); setnames(d, c("A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P"), paste0(c("A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P"), s)); d }
m <- merge(sfx(arch, "_arch"), sfx(meas, "_meas")[, !c("CHROM", "GENPOS")], by = c("ID", "ALLELE0", "ALLELE1"))
m <- merge(m, sfx(cens, "_cens")[, !c("CHROM", "GENPOS")], by = c("ID", "ALLELE0", "ALLELE1"))
if (nrow(m) < 0.99 * nrow(arch)) stop("GATE: only ", nrow(m), " of ", nrow(arch), " archived variants aligned")

cis <- m[CHROM == 6 & GENPOS >= 160500000 & GENPOS <= 161700000]
cis_gws <- cis[LOG10P_arch > -log10(5e-8) & !is.na(BETA_meas) & !is.na(BETA_cens)]
slope <- function(y, x) sum(x * y) / sum(x * x)

ctl <- cis_gws[LOG10P_arch > 10]
ctl_r <- cor(ctl$BETA_meas, ctl$BETA_arch); ctl_s <- slope(ctl$BETA_meas, ctl$BETA_arch)
ctl_sign <- mean(sign(ctl$BETA_meas) == sign(ctl$BETA_arch))
lines <- c(
  sprintf("LPA-locus GWS variants in published GWAS: %d (control set, LOG10P>10: %d)", nrow(cis_gws), nrow(ctl)),
  sprintf("CONTROL  LPA_MEAS vs published: r=%.5f slope=%.4f sign=%.4f", ctl_r, ctl_s, ctl_sign),
  sprintf("N at lead: published %d, LPA_MEAS %d, LPA_CENS %d",
          m[ID == "rs10455872", N_arch], m[ID == "rs10455872", N_meas], m[ID == "rs10455872", N_cens]),
  sprintf("CENSORED LPA_CENS vs published, all locus GWS: slope=%.4f median beta ratio=%.4f median chi2 ratio=%.3f",
          cis_gws[, slope(BETA_cens, BETA_arch)], cis_gws[, median(BETA_cens / BETA_arch)], cis_gws[, median(CHISQ_cens / CHISQ_arch)]),
  sprintf("         sign agreement %.4f, %d variants above LOG10P 7.3 in LPA_CENS vs %d published",
          cis_gws[, mean(sign(BETA_cens) == sign(BETA_arch))], cis[LOG10P_cens > -log10(5e-8), .N], cis[LOG10P_arch > -log10(5e-8), .N]))

# by effect size: the simulation predicts the largest-effect variants move differently
cis_gws[, bin := cut(abs(BETA_arch), c(0, 0.1, 0.25, 0.5, 1, Inf), right = FALSE)]
by_size <- cis_gws[, .(n = .N, slope_cens = slope(BETA_cens, BETA_arch), slope_meas = slope(BETA_meas, BETA_arch),
                       med_chi2_ratio = median(CHISQ_cens / CHISQ_arch)), keyby = bin]
by_maf <- cis_gws[, maf := pmin(A1FREQ_arch, 1 - A1FREQ_arch)][
  , .(n = .N, slope_cens = slope(BETA_cens, BETA_arch)), keyby = .(maf_bin = cut(maf, c(0, 0.01, 0.05, 0.5)))]

# calibration away from the locus: the blocks must not inflate the null
far <- m[!(CHROM == 6 & GENPOS >= 155000000 & GENPOS <= 166000000) & A1FREQ_arch > 0.01 & A1FREQ_arch < 0.99]
lam <- function(x) median(x, na.rm = TRUE) / qchisq(0.5, 1)
lines <- c(lines, sprintf("lambdaGC, MAF>1%%, outside chr6:155-166Mb, %d variants: published %.3f, LPA_MEAS %.3f, LPA_CENS %.3f",
                          nrow(far), lam(far$CHISQ_arch), lam(far$CHISQ_meas), lam(far$CHISQ_cens)))
gws_far <- far[LOG10P_arch > -log10(5e-8)]
if (nrow(gws_far))
  lines <- c(lines, sprintf("Off-locus published GWS (%d): LPA_CENS/published slope=%.4f, median chi2 ratio=%.3f; off-locus GWS count published %d, LPA_CENS %d",
                            nrow(gws_far), gws_far[, slope(BETA_cens, BETA_arch)], gws_far[, median(CHISQ_cens / CHISQ_arch)],
                            far[LOG10P_arch > -log10(5e-8), .N], far[LOG10P_cens > -log10(5e-8), .N]))

known <- c("rs10455872", "rs3798220", "rs41272114", "rs74617384", "rs143431368", "rs140570886", "rs118039278", "rs1800769", "rs6415084", "rs73596816")
kv <- cis[ID %in% known][order(-LOG10P_arch)]
kv[, ratio_cens := BETA_cens / BETA_arch]
kv[, ratio_meas := BETA_meas / BETA_arch]
keep_cols <- c("ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ_arch", "N_arch", "N_meas", "N_cens", "BETA_arch", "BETA_meas", "BETA_cens",
               "SE_arch", "SE_cens", "LOG10P_arch", "LOG10P_meas", "LOG10P_cens", "ratio_meas", "ratio_cens")
fwrite(kv[, ..keep_cols], paste0(opt$out, "_known_variants.csv"))
fwrite(cis_gws[, .(ID, GENPOS, ALLELE0, ALLELE1, A1FREQ_arch, BETA_arch, BETA_meas, BETA_cens, SE_arch, SE_cens, LOG10P_arch, LOG10P_cens)],
       paste0(opt$out, "_cis_gws.csv"))
writeLines(c(lines, "by |published beta|:", capture.output(print(by_size)), "by MAF:", capture.output(print(by_maf)),
             "known variants:", capture.output(print(kv[, .(ID, A1FREQ_arch, BETA_arch, BETA_meas, BETA_cens, ratio_meas, ratio_cens, LOG10P_arch, LOG10P_cens)]))),
           paste0(opt$out, "_summary.txt"))
writeLines(readLines(paste0(opt$out, "_summary.txt")))

if (is.na(ctl_r) || ctl_r < opt$min_r) stop(sprintf("GATE: control r=%.5f < %.3f", ctl_r, opt$min_r))
if (abs(ctl_s - 1) > opt$slope_tol) stop(sprintf("GATE: control slope=%.4f outside 1 +/- %.2f", ctl_s, opt$slope_tol))
if (ctl_sign < 0.999) stop(sprintf("GATE: control sign agreement %.4f", ctl_sign))
message("[cmp] control gate passed")
