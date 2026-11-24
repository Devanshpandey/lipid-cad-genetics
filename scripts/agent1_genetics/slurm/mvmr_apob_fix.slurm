#!/bin/bash
#SBATCH -J apobfix
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/apobfix_%j.out
#SBATCH -e /path/to/cad-genetics/apobfix_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# Table S3b: focused LDL-C / ApoB / triglyceride model, the basis for the
# identifiability claim in Section 2.4. Re-run with LOG10P-ranked selection and the
# statin-corrected LDL-C. ApoB is deliberately NOT medication-corrected, so LDL-C and
# ApoB are not differentially rescaled in treated participants.
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
ln -sf $_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_APOB.txt.gz $S/lipids_APOB.txt.gz
mkdir -p $_REPO/exports/agent1/mvmr_apob_fixed
for T in CAD MI MACE; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr_apob.R --outcome $T --sumstat_dir $S \
    --out $_REPO/exports/agent1/mvmr_apob_fixed/mvmr_apob_$T.csv
done
echo ALL_DONE
