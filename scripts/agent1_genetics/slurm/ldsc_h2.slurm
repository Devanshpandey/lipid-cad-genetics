#!/bin/bash
#SBATCH -J ldsch2
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/ldsch2_%j.out
#SBATCH -e /path/to/cad-genetics/ldsch2_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_LDSC}"; set -u
W="${LDSC_WEIGHTS_DIR}/weights.hm3_noMHC."
O=$_REPO/results/agent1_genetics/ldsc_statinfix; mkdir -p $O
# LDSC intercept separates polygenicity from confounding; lambda_GC alone cannot.
# Run BOTH the corrected and control LDL-C so the intercepts are directly comparable.
for TAG in control statinfix; do
  G=$_REPO/results/gwas/LDL_C_${TAG}/LDL_C.gwas.tsv.gz
  echo "=== munge $TAG ==="
  python "${MUNGE_PY}" --sumstats $G --out $O/LDL_C_${TAG} \
    --snp SNP --a1 EA --a2 NEA --p P --signed-sumstats BETA,0 --N-col N
  echo "=== h2 $TAG ==="
  python "${LDSC_PY}" --h2 $O/LDL_C_${TAG}.sumstats.gz \
    --ref-ld-chr "${LDSC_REF_PREFIX}" --w-ld-chr "$W" --out $O/h2_LDL_C_${TAG}
done
echo "=== INTERCEPTS ==="
grep -H -E "Intercept|Total Observed scale h2|Lambda GC|Ratio" $O/h2_LDL_C_*.log
echo ALL_DONE
