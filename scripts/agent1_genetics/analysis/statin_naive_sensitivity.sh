#!/usr/bin/env bash
# ============================================================
# Medication-correction sensitivity: compare LDL-C instrument effects in
# statin-naive participants (uncorrected LDL-C) vs the full corrected analysis
# (LDL-C /0.7 for lipid-lowering users). Also compares CAD outcome effects
# under minimal vs full covariate sets (collider/overadjustment check).
#
# Inputs: PLINK BED/BIM/FAM genotypes ($BFILE), phenotype tables, covariate file,
#         instrument rsID list ($IVS). Requires PLINK 2.0.
# Usage: BFILE=... PB=pheno_binary.txt PQ=pheno_quantitative.txt COV=covariates.txt IVS=instruments.txt \
#        bash statin_naive_sensitivity.sh
# ============================================================
set -euo pipefail
PLINK2=${PLINK2:-plink2}; W=${W:-./medtest}; mkdir -p "$W"
MINCOV="age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
FULLCOV="age,sex,bmi,smoking_current,diabetes,hypertension,${MINCOV}"

# --- (5) statin-naive LDL-C instrument effects vs full ---
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="STATIN_USE")c=i; next} $c==0{print $1,$2}' "$PB" > "$W/statin_naive.txt"
for tag in full naive; do
  KEEP=""; [ "$tag" = naive ] && KEEP="--keep $W/statin_naive.txt"
  $PLINK2 --bfile "$BFILE" $KEEP --extract "$IVS" --pheno "$PQ" --pheno-name LDL_C \
    --covar "$COV" --covar-name "$MINCOV" --glm hide-covar --covar-variance-standardize \
    --out "$W/ldl_${tag}"
done

# --- (4) CAD outcome effects: minimal vs full covariates ---
for tag in minimal fullc; do
  CN=$MINCOV; [ "$tag" = fullc ] && CN=$FULLCOV
  $PLINK2 --bfile "$BFILE" --1 --extract "$IVS" --pheno "$PB" --pheno-name CAD \
    --covar "$COV" --covar-name "$CN" --glm firth-fallback hide-covar --covar-variance-standardize \
    --out "$W/cad_${tag}"
done
echo "Compare betas across $W/ldl_full vs ldl_naive (LDL correction) and cad_minimal vs cad_fullc (covariates)."
