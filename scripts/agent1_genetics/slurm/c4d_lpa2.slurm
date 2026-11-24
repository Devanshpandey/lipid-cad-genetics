#!/bin/bash
#SBATCH -J c4dlpa2
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/c4dlpa2_%j.out
#SBATCH -e /path/to/cad-genetics/c4dlpa2_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
M=$_REPO/results/agent1_genetics/gwas_summary_stats/merged
REF=$SCRATCH/lpa_ldref
# Stage 1: union of genome-wide-significant variant IDs across the four exposures.
# The previous run scanned the 785 GB panel once PER TRAIT and timed out; extracting
# the needed variants ONCE turns four full scans into one.
echo "[1] collecting GWS variants"
for T in LDL_C HDL_C TRIGLY LPA; do zcat $M/lipids_$T.txt.gz | awk 'NR>1 && $12>7.30103 {print $3}'; done | sort -u > $SCRATCH/gws_union.txt
wc -l < $SCRATCH/gws_union.txt
echo "[2] extracting LD subset (single scan)"
"${PLINK2}" --bfile "${UKB_GENO_STEP1}" --extract $SCRATCH/gws_union.txt --make-bed --out $REF --threads ${SLURM_NTASKS} --memory 100000
echo "[3] MR"
Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
  --sumstat_dir $M \
  --outcome_file $_REPO/data/external/cardiogram/cad.add.160614.website.txt \
  --outcome_format cardiogram --outcome_label CARDIoGRAMplusC4D \
  --exposures LDL_C,HDL_C,TRIGLY,LPA --mode univariable \
  --ld_ref $REF --plink /path/to/software/plink \
  --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
  --out $_REPO/exports/agent1/mr_cardiogram_with_LPA.csv
