#!/bin/bash
#SBATCH -J apobgcn
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/apobgcn_%j.out
#SBATCH -e /path/to/cad-genetics/apobgcn_%j.err
# Same SNP-exposure covariance correction (job 3448312) applied to the STATIN-NAIVE arm.
# Correcting only the primary would leave Table S38 internally inconsistent: primary
# conditional F would read 48/49 beside a statin-naive 4.8/5.1 computed with gencov=0.
# The phenotypic correlation must be estimated IN the statin-naive subsample, not the
# full cohort, because that is the sample the arm's GWAS were run in.
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
P=$_REPO/results/agent1_genetics/phenotypes
SN=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinnaive
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
W=$SCRATCH/apobgcn; mkdir -p $W
OUT=$_REPO/exports/agent1/mvmr_apob_gencov; mkdir -p $OUT

for O in CAD MI MACE; do ln -sf $S/outcomes_${O}.txt.gz $SN/outcomes_${O}.txt.gz 2>/dev/null; done

# restrict the phenotype table to the statin-naive analysis sample
awk 'NR==FNR{k[$1];k[$2];next} FNR==1{print;next} ($1 in k)||($2 in k)' \
  $P/regenie_keep_statinnaive.txt $P/pheno_quantitative.txt > $W/pheno_naive.txt
echo "statin-naive phenotype rows: $(( $(wc -l < $W/pheno_naive.txt) - 1 ))"

for T in CAD MI MACE; do
  echo "########## statin-naive $T ##########"
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr_apob.R \
    --sumstat_dir $SN --outcome $T --pheno $W/pheno_naive.txt \
    --out $OUT/mvmr_apob_statinnaive_${T}.csv || echo "  $T FAILED"
done
echo ALL_DONE
