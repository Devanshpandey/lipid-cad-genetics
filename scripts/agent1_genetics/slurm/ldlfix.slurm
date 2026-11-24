#!/bin/bash
#SBATCH -J ldlfix
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:45:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/ldlfix_%j.out
#SBATCH -e /path/to/cad-genetics/ldlfix_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u
Rscript $_REPO/scripts/agent1_genetics/R/fix_statin_ldl.R \
  --pheno   $_REPO/results/agent1_genetics/phenotypes/pheno_quantitative.txt \
  --raw_dir /path/to/ukb_phenotypes \
  --out     $_REPO/results/agent1_genetics/phenotypes/pheno_quantitative_statinfix.txt
