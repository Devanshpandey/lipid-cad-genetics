#!/bin/bash
#SBATCH -J mvmr_apob
#SBATCH -o logs/mvmr_apob_%j.out
#SBATCH -e logs/mvmr_apob_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:40:00
#SBATCH -A YOUR_ALLOCATION
cd /path/to/cad-genetics
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env
mkdir -p results/agent1_genetics/mvmr_apob
for O in CAD MI MACE; do
  echo "=== MVMR-ApoB outcome $O ==="
  Rscript scripts/agent1_genetics/R/run_mvmr_apob.R \
    --outcome $O \
    --sumstat_dir results/agent1_genetics/gwas_summary_stats/merged \
    --out results/agent1_genetics/mvmr_apob/mvmr_apob_${O}.csv
done
echo "ALL DONE"
