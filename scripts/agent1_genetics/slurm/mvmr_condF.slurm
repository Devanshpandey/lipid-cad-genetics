#!/bin/bash
#SBATCH -J mvmr_cF
#SBATCH -o logs/mvmr_condF_%j.out
#SBATCH -e logs/mvmr_condF_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A YOUR_ALLOCATION
cd /path/to/cad-genetics
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env
mkdir -p results/agent1_genetics/mvmr_condF
echo "=== 4-exposure model (CAD) ==="
Rscript scripts/agent1_genetics/R/run_mvmr.R --outcome CAD \
  --sumstat_dir results/agent1_genetics/gwas_summary_stats/merged \
  --out results/agent1_genetics/mvmr_condF/mvmr4_CAD.csv
echo "=== LDL/ApoB/TG model (CAD) ==="
Rscript scripts/agent1_genetics/R/run_mvmr_apob.R --outcome CAD \
  --sumstat_dir results/agent1_genetics/gwas_summary_stats/merged \
  --out results/agent1_genetics/mvmr_condF/mvmrApoB_CAD.csv
echo ALL DONE
