#!/bin/bash
#SBATCH -J fgrerun
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/fgrerun_%j.out
#SBATCH -e /path/to/cad-genetics/fgrerun_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
D=$_REPO/data/external/finngen; mkdir -p $D
# Table S7 currently reports a two-sample MVMR run with the ORIGINAL instruments --
# before the p-value-underflow fix and before the LDL-C statin correction. It is the
# one replication column not yet recomputed, so the paper's three columns are
# currently produced three different ways. R11 is the release cited in the manuscript.
if [ ! -s $D/finngen_R11_I9_CHD.gz ]; then
  echo "[dl] downloading FinnGen R11 I9_CHD"
  curl -sSL --retry 3 -o $D/finngen_R11_I9_CHD.gz \
    https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_I9_CHD.gz
fi
ls -lh $D/finngen_R11_I9_CHD.gz
echo "[dl] download present"
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
  --sumstat_dir $S \
  --outcome_file $D/finngen_R11_I9_CHD.gz \
  --outcome_format finngen --outcome_label FinnGen_R11_CAD \
  --exposures LDL_C,HDL_C,TRIGLY,LPA --mode mvmr \
  --ld_ref $SCRATCH/lpa_ldref --plink /path/to/software/plink \
  --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
  --out $_REPO/exports/agent1/mvmr_finngen_R11_corrected.csv
echo ALL_DONE
