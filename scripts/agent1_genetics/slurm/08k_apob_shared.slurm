#!/bin/bash
#SBATCH -J apobshare
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/apobshare_%j.out
#SBATCH -e /path/to/cad-genetics/apobshare_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
ln -sf $_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_APOB.txt.gz  $S/lipids_APOB.txt.gz
ln -sf $_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_APOA1.txt.gz $S/lipids_APOA1.txt.gz
Rscript $_REPO/scripts/agent1_genetics/R/apob_ldl_shared.R \
  --sumstat_dir $S \
  --out $_REPO/exports/agent1/revision/apob_ldl_shared.csv
echo ALL_DONE
