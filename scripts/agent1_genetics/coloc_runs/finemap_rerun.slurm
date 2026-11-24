#!/bin/bash
#------------------------------------------------------------
# SuSiE fine-mapping at the CORRECTED colocalized loci.
#
# The deposited TableS5 covers 5 loci, because fine-mapping was run against the
# buggy coloc output. The corrected LD-clumped scan gives 70 colocalized
# sentinels / 108 unique exposure x locus combinations; this re-runs SuSiE over
# all of them. Array = one task per exposure so no single task carries 108 loci.
#
# Submit: sbatch --array=0-5 finemap_rerun.slurm
#------------------------------------------------------------
#SBATCH -J finemap_rerun
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o finemap_%A_%a.out
#SBATCH -e finemap_%A_%a.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
_REPO=/path/to/cad-genetics
source "${_REPO}/config/tacc_paths.sh"
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u

EXPS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB)
EX=${EXPS[$SLURM_ARRAY_TASK_ID]}
W=$SCRATCH/finemap_rerun
cd $W

[ -s "coloc_${EX}.csv" ] || { echo "FATAL: missing coloc_${EX}.csv"; exit 1; }
echo "=== SuSiE fine-mapping: ${EX} ==="
Rscript "${_REPO}/scripts/agent1_genetics/R/run_finemapping.R" \
  --coloc_results "coloc_${EX}.csv" \
  --sumstat_dir   "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --ld_plink_dir  "${UKB_GENO_STEP2}" \
  --window_kb     500 \
  --pp4_thresh    0.8 \
  --out           "out/finemap_${EX}.csv"
echo "TASK_DONE ${EX}"
