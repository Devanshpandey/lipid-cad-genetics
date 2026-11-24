#!/bin/bash
#SBATCH -J fstat
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/fstat_%j.out
#SBATCH -e /path/to/cad-genetics/fstat_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u
Rscript $_REPO/scripts/agent1_genetics/R/recompute_mvmr_fstats.R --outcome CAD --sumstat_dir $_REPO/results/agent1_genetics/gwas_summary_stats/merged --out $_REPO/exports/agent1/mvmr_fstats_CAD_fixed.csv
