#!/usr/bin/env bash
# ============================================================
# External colocalization replication: UK Biobank lipids x FinnGen CAD
# at the five key lipid-CAD loci (PCSK9, SORT1/CELSR2, LDLR, APOE, TRIB1).
#
# UK Biobank sumstats are GRCh37, FinnGen R11 is GRCh38 -> match by rsID.
# Runs coloc.abf under the default and a stringent (p12=1e-6) prior.
#
# Inputs:
#   $MERGED/lipids_<TRAIT>.txt.gz   REGENIE lipid sumstats (b37; cols incl. ID, ALLELE0/1, A1FREQ, N, BETA, SE)
#   $FINNGEN_CHD                    finngen_R11_I9_CHD.gz (public; cols chrom pos ref alt rsids ... beta sebeta af_alt)
# Output: finngen_coloc_results.csv
# Usage: MERGED=... FINNGEN_CHD=... bash finngen_colocalization.sh
# ============================================================
set -euo pipefail
W=${W:-./fgcoloc}; mkdir -p "$W"
cat > "$W/loci.txt" <<'EOF'
PCSK9 1 55000000 56000000 LDL_C
SORT1_CELSR2 1 109300000 110300000 LDL_C
LDLR 19 10700000 11700000 LDL_C
APOE_APOC1 19 44900000 45900000 LDL_C
TRIB1 8 126000000 127000000 TRIGLY
EOF
while read -r name chr lo hi lipid; do
  zcat "$MERGED/lipids_${lipid}.txt.gz" | awk -v c="$chr" -v lo="$lo" -v hi="$hi" \
    'NR==1||($1==c&&$2>=lo&&$2<=hi){print $3,$4,$5,$6,$7,$9,$10}' > "$W/${name}_ukb.txt"
  awk 'NR>1{print $1}' "$W/${name}_ukb.txt" | grep '^rs' > "$W/${name}_rs.txt"
  { zcat "$FINNGEN_CHD" | head -1; zcat "$FINNGEN_CHD" | grep -Fwf "$W/${name}_rs.txt"; } > "$W/${name}_fg.txt"
done < "$W/loci.txt"

Rscript - "$W" <<'RS'
suppressPackageStartupMessages(library(coloc))
W <- commandArgs(trailingOnly=TRUE)[1]
loci <- read.table(file.path(W,"loci.txt"), col.names=c("name","chr","lo","hi","lipid"))
res <- data.frame()
for (i in seq_len(nrow(loci))) {
  nm <- loci$name[i]
  u <- read.table(file.path(W,paste0(nm,"_ukb.txt")), header=TRUE)
  colnames(u) <- c("rsid","A0","A1","af","N","beta","se")
  f <- read.table(file.path(W,paste0(nm,"_fg.txt")), header=TRUE, comment.char="")
  colnames(f) <- c("chrom","pos","ref","alt","rsid","genes","pval","mlogp","beta","sebeta","af","afc","afco")
  m <- merge(u, f, by="rsid")
  m <- m[!is.na(m$beta.x)&!is.na(m$se)&!is.na(m$beta.y)&!is.na(m$sebeta)&m$se>0&m$sebeta>0,]
  m$flip <- ifelse(toupper(m$alt)==toupper(m$A1),1,ifelse(toupper(m$ref)==toupper(m$A1),-1,NA))
  m <- m[!is.na(m$flip),]; m$bf <- m$beta.y*m$flip
  m$maf <- pmin(m$af.x,1-m$af.x); m <- m[m$maf>0&m$maf<1 & !duplicated(m$rsid),]
  D1 <- list(beta=m$beta.x, varbeta=m$se^2, N=round(mean(m$N)), MAF=m$maf, type="quant", snp=m$rsid)
  D2 <- list(beta=m$bf, varbeta=m$sebeta^2, type="cc", s=0.13, N=394000, MAF=m$maf, snp=m$rsid)
  c1 <- coloc.abf(D1,D2); c2 <- coloc.abf(D1,D2,p12=1e-6)
  res <- rbind(res, data.frame(locus=nm, nsnp=nrow(m),
    PPH4_default=round(c1$summary["PP.H4.abf"],3),
    PPH4_p12_1e6=round(c2$summary["PP.H4.abf"],3)))
}
write.csv(res, file.path(W,"finngen_coloc_results.csv"), row.names=FALSE)
print(res, row.names=FALSE)
RS
