#!/bin/bash
#SBATCH -J mvmrall
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/mvmrall_%j.out
#SBATCH -e /path/to/cad-genetics/mvmrall_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u
mkdir -p $_REPO/exports/agent1/mvmr_full
for T in CAD MI REVASC STROKE HF CV_DEATH MACE; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $_REPO/results/agent1_genetics/gwas_summary_stats/merged --out $_REPO/exports/agent1/mvmr_full/mvmr_$T.csv
done
echo ALL_DONE
