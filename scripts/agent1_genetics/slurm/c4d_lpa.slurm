#!/bin/bash
#SBATCH -J c4dlpa
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/c4dlpa_%j.out
#SBATCH -e /path/to/cad-genetics/c4dlpa_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u
Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
  --sumstat_dir $_REPO/results/agent1_genetics/gwas_summary_stats/merged \
  --outcome_file $_REPO/data/external/cardiogram/cad.add.160614.website.txt \
  --outcome_format cardiogram \
  --outcome_label CARDIoGRAMplusC4D \
  --exposures LDL_C,HDL_C,TRIGLY,LPA \
  --mode univariable \
  --ld_ref "${UKB_GENO_STEP1}" \
  --plink /path/to/software/plink \
  --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
  --out $_REPO/exports/agent1/mr_cardiogram_with_LPA.csv
