#!/bin/bash
#------------------------------------------------------------
# Statin-naive LDL-C / ApoB / triglyceride GWAS — REGENIE Step 2.
# Submit: sbatch --array=1-22 --dependency=afterok:<step1id> 08g_statinnaive_step2.slurm
#------------------------------------------------------------
#SBATCH -J snaive2
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/snaive2_%A_%a.out
#SBATCH -e /path/to/cad-genetics/snaive2_%A_%a.err
set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"
CHR="${SLURM_ARRAY_TASK_ID}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"
PH="${AGENT1_OUT}/phenotypes"
OUT_DIR="${AGENT1_OUT}/gwas_summary_stats/lipids_statinnaive"
mkdir -p "${OUT_DIR}"
regenie \
  --step 2 \
  --bed "${UKB_GENO_STEP2}" \
  --keep "${PH}/regenie_keep_statinnaive.txt" \
  --phenoFile "${PH}/pheno_quantitative.txt" \
  --phenoColList LDL_C,APOB,TRIGLY \
  --covarFile "${PH}/covariates.txt" \
  --pred "${AGENT1_OUT}/regenie_step1_statinnaive/step1_naive_pred.list" \
  --apply-rint \
  --bsize 400 \
  --chr "${CHR}" \
  --minMAC 5 \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${OUT_DIR}/lipids_chr${CHR}"
echo "[$(date)] chr ${CHR} complete"
