#!/bin/bash
#------------------------------------------------------------
# SLURM array: REGENIE Step 2 — GWAS for CAD/MACE outcomes
# Submit: sbatch --array=1-22 04_regenie_step2_outcomes.slurm
#------------------------------------------------------------
#SBATCH -J cad_gwas_outcomes
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -A OTH26002
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=devanshpandey@utexas.edu

set -euo pipefail
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

CHR="${SLURM_ARRAY_TASK_ID}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"

OUT_DIR="${AGENT1_OUT}/gwas_summary_stats/outcomes"
mkdir -p "${OUT_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_ARRAY_JOB_ID:-0}_${CHR}] REGENIE Step 2 — Outcomes GWAS — Chr ${CHR}"

regenie \
  --step 2 \
  --bed "${UKB_GENO_STEP2}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_binary.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
  --pred "${AGENT1_OUT}/regenie_step1/step1_bt_pred.list" \
  --bt \
  --firth --approx \
  --pThresh 0.01 \
  --bsize 400 \
  --chr "${CHR}" \
  --minMAC 5 \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${OUT_DIR}/outcomes_chr${CHR}"

echo "[$(date)] Chr ${CHR} outcomes GWAS complete"
