#!/usr/bin/env bash
# ============================================================
# Focused cis-LPA-region instrument MR: Lp(a) -> FinnGen CAD.
# Uses genome-wide-significant, distance-clumped LPA-region variants
# (led by rs10455872) as instruments; IVW through the origin.
#
# Inputs: $MERGED/lipids_LPA.txt.gz (b37), $FINNGEN_CHD (finngen_R11_I9_CHD.gz)
# Output: prints per-variant Wald ratios + IVW estimate.
# Usage: MERGED=... FINNGEN_CHD=... bash cis_lpa_mr.sh
# ============================================================
set -euo pipefail
W=${W:-./lpamr}; mkdir -p "$W"
# LPA region (chr6 b37 ~160.5-161.7 Mb), GWS variants, greedy 250 kb clump
zcat "$MERGED/lipids_LPA.txt.gz" | \
  awk '$1==6 && $2>=160500000 && $2<=161700000 && $12>7.30103{print $1,$2,$3,$4,$5,$9,$10,$12}' | \
  sort -k8 -grk8 > "$W/lpa_gws.txt"
awk '{keep=1; for(i in kept) if(($2-kept[i])<250000 && ($2-kept[i])>-250000) keep=0;
      if(keep){print; kept[NR]=$2}}' "$W/lpa_gws.txt" > "$W/lpa_instr.txt"
awk '{print $3}' "$W/lpa_instr.txt" > "$W/lpa_rs.txt"
{ zcat "$MERGED/lipids_LPA.txt.gz" | head -1; zcat "$MERGED/lipids_LPA.txt.gz" | grep -Fwf "$W/lpa_rs.txt"; } > "$W/lpa_ukb.txt"
{ zcat "$FINNGEN_CHD" | head -1; zcat "$FINNGEN_CHD" | grep -Fwf "$W/lpa_rs.txt"; } > "$W/lpa_fg.txt"

Rscript - "$W" <<'RS'
W <- commandArgs(trailingOnly=TRUE)[1]
u <- read.table(file.path(W,"lpa_ukb.txt"), header=TRUE)[,c("ID","ALLELE0","ALLELE1","BETA","SE")]
colnames(u) <- c("rsid","A0","A1","bx","sx")
f <- read.table(file.path(W,"lpa_fg.txt"), header=TRUE, comment.char="")
colnames(f) <- c("chrom","pos","ref","alt","rsid","genes","pval","mlogp","beta","sebeta","af","afc","afco")
m <- merge(u, f[,c("rsid","ref","alt","beta","sebeta")], by="rsid")
m$flip <- ifelse(toupper(m$alt)==toupper(m$A1),1,ifelse(toupper(m$ref)==toupper(m$A1),-1,NA))
m <- m[!is.na(m$flip),]; m$by <- m$beta*m$flip
for (i in seq_len(nrow(m)))
  cat(sprintf("  %-12s Lp(a)beta=%.3f CADbeta=%+.4f WaldOR=%.3f\n", m$rsid[i], m$bx[i], m$by[i], exp(m$by[i]/m$bx[i])))
w <- 1/(m$sebeta^2/m$bx^2); ratio <- m$by/m$bx
b <- sum(w*ratio)/sum(w); se <- sqrt(1/sum(w))
cat(sprintf("\nIVW Lp(a)->CAD OR/SD = %.3f (95%% CI %.3f-%.3f) P=%.2e (%d cis-LPA instruments)\n",
  exp(b), exp(b-1.96*se), exp(b+1.96*se), 2*pnorm(-abs(b/se)), nrow(m)))
RS
