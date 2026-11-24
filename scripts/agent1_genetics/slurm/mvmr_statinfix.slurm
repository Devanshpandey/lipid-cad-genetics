#!/bin/bash
#SBATCH -J mvmrfix2
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/mvmrfix2_%j.out
#SBATCH -e /path/to/cad-genetics/mvmrfix2_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
M=$_REPO/results/agent1_genetics/gwas_summary_stats/merged
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
mkdir -p $S
# everything except LDL_C is unchanged -> symlink, do not duplicate 7 GB
for f in $M/*.txt.gz; do b=$(basename $f); [ "$b" = "lipids_LDL_C.txt.gz" ] || ln -sf $f $S/$b; done
echo "[1] staging corrected LDL_C (with allele-mapping verification)"
Rscript $_REPO/scripts/agent1_genetics/R/stage_statinfix_sumstats.R \
  --new     $_REPO/results/gwas/LDL_C_statinfix/LDL_C.gwas.tsv.gz \
  --archive $M/lipids_LDL_C.txt.gz \
  --out     $S/lipids_LDL_C.txt.gz
echo "[2] MVMR across all outcomes"
mkdir -p $_REPO/exports/agent1/mvmr_statinfix
for T in CAD MI REVASC STROKE HF CV_DEATH MACE; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $S \
    --out $_REPO/exports/agent1/mvmr_statinfix/mvmr_$T.csv
done
echo ALL_DONE
