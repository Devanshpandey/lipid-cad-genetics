#!/bin/bash
#SBATCH -J minadj2
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/minadj2_%j.out
#SBATCH -e /path/to/cad-genetics/minadj2_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# Minimally-adjusted sensitivity arm (Table S32), re-run after the distance-pruning fix
# so that it is comparable with the corrected 161-instrument primary model.
# This arm's LDL-C is NOT statin-ascertainment-corrected.
M=$_REPO/results/agent1_genetics/gwas_summary_stats/minadj_merged
mkdir -p $_REPO/exports/agent1/mvmr_minadj_prunefix
for T in CAD MI MACE; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $M \
    --out $_REPO/exports/agent1/mvmr_minadj_prunefix/mvmr_$T.csv
done
echo ALL_DONE
