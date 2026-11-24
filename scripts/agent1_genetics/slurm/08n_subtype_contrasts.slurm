#!/bin/bash
#SBATCH -J subcon
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 03:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/subcon_%j.out
#SBATCH -e /path/to/cad-genetics/subcon_%j.err
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
OUT=$_REPO/exports/agent1/revision
Rscript $_REPO/scripts/agent1_genetics/R/stroke_subtype_contrasts.R \
  --sumstat_dir $S \
  --megastroke_dir $_REPO/data/external/megastroke \
  --ld_ref /path/to/scratch/lpa_ldref \
  --plink /path/to/software/plink \
  --nboot 1000 \
  --out_prefix $OUT/stroke_subtype_contrasts
echo ALL_DONE
