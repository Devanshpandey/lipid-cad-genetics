#!/bin/bash
#------------------------------------------------------------
# SLURM: Multivariable Mendelian Randomization (MVMR)
#   All 7 lipid traits as joint exposures → one outcome at a time
#   7 jobs, one per outcome (submit as array or loop)
#
#   Submit as array: sbatch --array=0-6 08b_mvmr.slurm
#   Or looped via run_jupyter.sh (no sbatch needed)
#------------------------------------------------------------
#SBATCH -J cad_mvmr
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 02:00:00
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

OUTCOME_TRAITS=(CAD MI REVASC STROKE HF CV_DEATH MACE)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
OUTCOME="${OUTCOME_TRAITS[$IDX]}"

OUT_DIR="${AGENT1_OUT}/mvmr"
mkdir -p "${OUT_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_ARRAY_JOB_ID:-0}_${IDX}] MVMR: all lipids → ${OUTCOME}"

Rscript "${SCRIPTS_DIR}/agent1_genetics/R/run_mvmr.R" \
  --outcome     "${OUTCOME}" \
  --sumstat_dir "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --out         "${OUT_DIR}/mvmr_${OUTCOME}.csv"

echo "[$(date)] Done: MVMR → ${OUTCOME}"
