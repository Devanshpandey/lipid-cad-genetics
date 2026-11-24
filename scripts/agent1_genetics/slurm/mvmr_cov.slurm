#!/bin/bash
#SBATCH -J mvmrcov
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/mvmrcov_%j.out
#SBATCH -e /path/to/cad-genetics/mvmrcov_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
Rscript $_REPO/scripts/agent1_genetics/R/mvmr_covariance.R \
  --sumstat_dir $_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix \
  --pheno       $_REPO/results/agent1_genetics/phenotypes/pheno_quantitative_statinfix.txt \
  --outcomes    CAD,STROKE,MI,REVASC \
  --nboot       1000 \
  --cis_list    /path/to/scratch/statint/lpa_cis_ids.txt \
  --out         $_REPO/exports/agent1/mvmr_covariance.csv
echo ALL_DONE
