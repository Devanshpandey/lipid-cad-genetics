#!/bin/bash
#SBATCH -J mvmrmin
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/mvmrmin_%j.out
#SBATCH -e /path/to/cad-genetics/mvmrmin_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# Minimally-adjusted sensitivity arm, re-run with the LOG10P-ranked instrument
# selection. Its published numbers (236 instruments; conditional F 52.9/51.9/39.6/13.0)
# came from the p-value-underflow selection and are not comparable to the corrected
# primary model. NOTE: this arm's LDL-C is NOT statin-ascertainment-corrected.
M=$_REPO/results/agent1_genetics/gwas_summary_stats/minadj_merged
mkdir -p $_REPO/exports/agent1/mvmr_minadj
for T in CAD MI MACE; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $M \
    --out $_REPO/exports/agent1/mvmr_minadj/mvmr_$T.csv
done
echo ALL_DONE
