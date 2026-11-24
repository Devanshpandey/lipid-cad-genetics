#!/bin/bash
#SBATCH -J statrob
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/statrob_%j.out
#SBATCH -e /path/to/cad-genetics/statrob_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
W=$SCRATCH/statint
Rscript $_REPO/scripts/agent1_genetics/R/statin_interaction_robustness.R \
  --scores  $W/dose_LPA.raw,$W/dose_LDL.raw \
  --weights $W/lpa_score.txt,$W/ldl_score.txt \
  --raw_dir /path/to/ukb_phenotypes \
  --covar   $_REPO/results/agent1_genetics/phenotypes/covariates.txt \
  --out     $_REPO/exports/agent1/revision/statin_robustness.csv
echo ALL_DONE
