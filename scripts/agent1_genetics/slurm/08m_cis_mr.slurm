#!/bin/bash
#SBATCH -J cismr
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/cismr_%j.out
#SBATCH -e /path/to/cad-genetics/cismr_%j.err
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
G=$_REPO/results/agent1_genetics/gwas_summary_stats
S=$G/merged_statinfix
OUT=$_REPO/exports/agent1/revision
mkdir -p $OUT
# non-LDL exposures live in the original merged dir; symlink so one dir has them all
for TR in TRIGLY HDL_C LPA APOB; do
  [ -e $S/lipids_${TR}.txt.gz ] || ln -sf $G/merged/lipids_${TR}.txt.gz $S/lipids_${TR}.txt.gz
done

echo "############ OLD arm: pre-correction LDL-C (field 6153 only) -- validation"
Rscript $_REPO/scripts/agent1_genetics/R/run_cis_mr.R \
  --sumstat_dir $S --outcome_dir $S \
  --ldl_file ../merged/lipids_LDL_C.txt.gz \
  --ld_ref /path/to/scratch/lpa_ldref --plink /path/to/software/plink --label old_6153only \
  --out $OUT/cis_mr_old6153.csv || echo "OLD ARM FAILED"

echo "############ NEW arm: corrected LDL-C (6153 U 6177)"
Rscript $_REPO/scripts/agent1_genetics/R/run_cis_mr.R \
  --sumstat_dir $S --outcome_dir $S \
  --ldl_file lipids_LDL_C.txt.gz \
  --ld_ref /path/to/scratch/lpa_ldref --plink /path/to/software/plink --label corrected_union \
  --out $OUT/cis_mr_corrected.csv || echo "NEW ARM FAILED"
echo ALL_DONE
