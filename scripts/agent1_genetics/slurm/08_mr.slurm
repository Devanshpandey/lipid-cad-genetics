#!/bin/bash
#------------------------------------------------------------
# SLURM array: Mendelian Randomization — all lipid → outcome pairs
# 7 lipids × 9 outcomes = 63 pairs
# New outcomes vs v1: added AF (atrial fibrillation), PAD (peripheral artery disease)
# LPA uses cis-MR mode automatically (chr6:157-163Mb instruments only)
# Submit: sbatch --array=0-62 08_mr.slurm
#------------------------------------------------------------
#SBATCH -J cad_mr
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
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

# 7 lipid × 9 outcome = 63 pairs
LIPID_TRAITS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA)
OUTCOME_TRAITS=(CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD)
N_OUTCOMES=${#OUTCOME_TRAITS[@]}

IDX="${SLURM_ARRAY_TASK_ID}"
LIPID="${LIPID_TRAITS[$((IDX / N_OUTCOMES))]}"
OUTCOME="${OUTCOME_TRAITS[$((IDX % N_OUTCOMES))]}"

OUT_DIR="${AGENT1_OUT}/mr"
mkdir -p "${OUT_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_ARRAY_JOB_ID:-0}_${IDX}] MR: ${LIPID} → ${OUTCOME}"

Rscript "${SCRIPTS_DIR}/agent1_genetics/R/run_mr.R" \
  --exposure    "${LIPID}" \
  --outcome     "${OUTCOME}" \
  --sumstat_dir "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --out         "${OUT_DIR}/${LIPID}_to_${OUTCOME}_mr.csv"

echo "[$(date)] Done: ${LIPID} → ${OUTCOME}"
