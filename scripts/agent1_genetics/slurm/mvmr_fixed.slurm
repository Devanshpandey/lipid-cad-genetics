#!/bin/bash
#SBATCH -J mvmrfix
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/mvmrfix_%A_%a.out
#SBATCH -e /path/to/cad-genetics/mvmrfix_%A_%a.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u
O=(CAD MI REVASC STROKE HF CV_DEATH MACE)
T=${O[${SLURM_ARRAY_TASK_ID:-0}]}
mkdir -p $_REPO/exports/agent1/mvmr_full
Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $_REPO/results/agent1_genetics/gwas_summary_stats/merged --out $_REPO/exports/agent1/mvmr_full/mvmr_$T.csv
