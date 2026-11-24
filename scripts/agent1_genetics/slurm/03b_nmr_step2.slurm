#!/bin/bash
#SBATCH -J nmr_gwas
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL
set -euo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
CHR="${SLURM_ARRAY_TASK_ID}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"
OUT_DIR="${AGENT1_OUT}/gwas_summary_stats/nmr"
mkdir -p "${OUT_DIR}"
echo "[$(date)] REGENIE Step 2 — NMR GWAS — Chr ${CHR}"
regenie \
  --step 2 \
  --bed "${UKB_GENO_STEP2}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_nmr.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
  --pred "${AGENT1_OUT}/regenie_step1/step1_nmr_pred.list" \
  --apply-rint \
  --bsize 400 \
  --chr "${CHR}" \
  --minMAC 5 \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${OUT_DIR}/nmr_chr${CHR}"
echo "[$(date)] Chr ${CHR} NMR GWAS complete"
