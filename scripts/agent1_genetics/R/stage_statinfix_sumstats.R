#!/usr/bin/env Rscript
# ============================================================
# Stage the statin-corrected LDL-C GWAS into merged/ format so the MVMR can
# consume it alongside the archived HDL_C / TRIGLY / LPA / outcome sumstats.
#
# run_gwas.sh emits:  CHR BP SNP EA NEA EAF N BETA SE CHISQ LOG10P P
# merged/ format is:  CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
# In REGENIE, ALLELE1 is the EFFECT allele (A1FREQ is its frequency), so
# ALLELE1<-EA and ALLELE0<-NEA. A silent flip here would invert every LDL-C
# estimate downstream, so the mapping is VERIFIED against the archived file
# before anything is written.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
set.seed(42)
opt <- parse_args(OptionParser(option_list=list(
  make_option("--new",     type="character"),
  make_option("--archive", type="character"),
  make_option("--out",     type="character"),
  make_option("--min_r",   type="double", default=0.98)
)))

new <- fread(cmd=paste("gunzip -c", shQuote(opt$new)), na.strings="NA")
arc <- fread(cmd=paste("gunzip -c", shQuote(opt$archive)), na.strings="NA")
message(sprintf("[stage] new=%d variants  archive=%d variants", nrow(new), nrow(arc)))

out <- data.table(CHROM=new$CHR, GENPOS=new$BP, ID=new$SNP,
                  ALLELE0=new$NEA, ALLELE1=new$EA, A1FREQ=new$EAF,
                  N=new$N, TEST="ADD", BETA=new$BETA, SE=new$SE,
                  CHISQ=new$CHISQ, LOG10P=new$LOG10P, EXTRA=NA)

# ---- verification: same-allele variants must agree in sign and magnitude ----
m <- merge(out[, .(ID, ALLELE0, ALLELE1, BETA)],
           arc[, .(ID, aA0=ALLELE0, aA1=ALLELE1, aBETA=BETA)], by="ID")
same <- m[ALLELE0==aA0 & ALLELE1==aA1]
flip <- m[ALLELE0==aA1 & ALLELE1==aA0]
message(sprintf("[stage] overlap=%d  same-orientation=%d  flipped=%d", nrow(m), nrow(same), nrow(flip)))
if (nrow(same) < 1000) stop("too few same-orientation variants to verify mapping")
# Restrict the concordance check to genome-wide-significant variants. Across ALL
# ~14.8M variants the statistic is dominated by null SNPs whose BETA is noise, and
# the corrected phenotype genuinely perturbs that noise (40,637 people had LDL-C
# raised by ~1.18 mmol/L) -- an all-variant check gives r~0.92 / sign~0.88 even
# when the mapping is perfect. Among GWS variants, where signal dominates, the
# same comparison gives r=0.998 and sign agreement 1.000.
gws <- same[!is.na(aBETA) & ID %in% arc[LOG10P > -log10(5e-8), ID]]
message(sprintf("[stage] concordance evaluated on %d GWS variants", nrow(gws)))
if (nrow(gws) < 500) stop("too few GWS variants to verify mapping")
r   <- cor(gws$BETA, gws$aBETA)
sgn <- mean(sign(gws$BETA) == sign(gws$aBETA))
message(sprintf("[stage] correlation with archive BETA = %.4f ; sign agreement = %.3f", r, sgn))
if (is.na(r) || r < opt$min_r) stop(sprintf("VERIFY FAILED: r=%.4f < %.2f -- allele mapping is wrong, refusing to write", r, opt$min_r))
if (sgn < 0.95) stop(sprintf("VERIFY FAILED: sign agreement %.3f -- effect alleles are inverted", sgn))
message("[stage] allele mapping verified (ALLELE1 = EA)")
# effect sizes should be LARGER after correction, not merely correlated
message(sprintf("[stage] mean |BETA| (GWS) archive=%.5f  corrected=%.5f  ratio=%.3f",
        mean(abs(gws$aBETA)), mean(abs(gws$BETA)), mean(abs(gws$BETA))/mean(abs(gws$aBETA))))

fwrite(out, opt$out, sep="\t", na="NA", quote=FALSE, compress="gzip")
message("[stage] wrote ", opt$out)
