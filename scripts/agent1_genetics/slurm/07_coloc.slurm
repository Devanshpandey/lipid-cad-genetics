#!/bin/bash
#------------------------------------------------------------
# SLURM array: Colocalization — all lipid × outcome pairs
# Array index maps to a lipid×outcome combination
# 7 lipids × 9 outcomes = 63 pairs (0-62)
# New outcomes vs v1: added AF (atrial fibrillation), PAD (peripheral artery disease)
# Submit: sbatch --array=0-62 07_coloc.slurm
#------------------------------------------------------------
#SBATCH -J cad_coloc
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 08:00:00
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

set +u  # conda activate scripts reference unbound vars (binutils ADDR2LINE)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u

LIPID_TRAITS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA)
OUTCOME_TRAITS=(CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD)
N_OUTCOMES=${#OUTCOME_TRAITS[@]}

# Map array index to lipid×outcome pair
IDX="${SLURM_ARRAY_TASK_ID}"
lipid_idx=$(( IDX / N_OUTCOMES ))
outcome_idx=$(( IDX % N_OUTCOMES ))
LIPID="${LIPID_TRAITS[$lipid_idx]}"
OUTCOME="${OUTCOME_TRAITS[$outcome_idx]}"

OUT_DIR="${AGENT1_OUT}/coloc"
mkdir -p "${OUT_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_ARRAY_JOB_ID:-0}_${IDX}] Coloc: ${LIPID} vs ${OUTCOME}"

Rscript "${SCRIPTS_DIR}/agent1_genetics/R/run_coloc.R" \
  --lipid_trait  "${LIPID}" \
  --outcome      "${OUTCOME}" \
  --sumstat_dir  "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --out          "${OUT_DIR}/${LIPID}_vs_${OUTCOME}_coloc.csv"

echo "[$(date)] Done: ${LIPID} vs ${OUTCOME}"
