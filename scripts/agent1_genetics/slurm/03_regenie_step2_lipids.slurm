#!/bin/bash
#------------------------------------------------------------
# SLURM array: REGENIE Step 2 — GWAS for lipid traits
# One job per chromosome (array index = chromosome number)
# Submits 22 jobs automatically via SLURM_ARRAY_TASK_ID
#
# Submit: sbatch --array=1-22 03_regenie_step2_lipids.slurm
#------------------------------------------------------------
#SBATCH -J cad_gwas_lipids
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
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

CHR="${SLURM_ARRAY_TASK_ID}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"

OUT_DIR="${AGENT1_OUT}/gwas_summary_stats/lipids"
mkdir -p "${OUT_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_ARRAY_JOB_ID:-0}_${CHR}] REGENIE Step 2 — Lipid GWAS — Chr ${CHR}"

# UKB_GENO_STEP2 = merged BED across all chromosomes (MAF > 0.001).
# REGENIE accepts --chr to restrict testing to one chromosome from a merged file.
regenie \
  --step 2 \
  --bed "${UKB_GENO_STEP2}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_quantitative.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
  --pred "${AGENT1_OUT}/regenie_step1/step1_qt_pred.list" \
  --apply-rint \
  --bsize 400 \
  --chr "${CHR}" \
  --minMAC 5 \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${OUT_DIR}/lipids_chr${CHR}"

echo "[$(date)] Chr ${CHR} lipid GWAS complete"
