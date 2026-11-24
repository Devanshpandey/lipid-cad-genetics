#!/bin/bash
#SBATCH -J lpaHfeas2
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:40:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaHfeas2_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaHfeas2_%j.err
# COUNTS ONLY: sizes and event counts of the cohorts that would be used to forecast
# OCEAN(a)-Outcomes and ACCLAIM-Lp(a). Fits no outcome model; the pre-registered script is
# sourced, not edited, and must stay 6e1c268df17113e7959d4078ee830b64.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
DIR=$_REPO/scripts/agent1_genetics/lpa_horizon
SCRIPT=$DIR/forecast_feas_counts.R
BASE_MD5=$(md5sum "$DIR/lpa_horizon.R" | cut -d' ' -f1)
echo "start $(date) on $(hostname); counts md5 $(md5sum "$SCRIPT" | cut -d' ' -f1); pre-registered md5 $BASE_MD5"
[ "$BASE_MD5" = "6e1c268df17113e7959d4078ee830b64" ] || { echo "pre-registered script has changed"; exit 1; }
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
PRS=$_REPO/results/agent6_resilience/run_cad_ukbfree_v2/10_synergy/prs_components.tsv
OUT=$_REPO/exports/agent1/lpa_horizon/forecast_feasibility
for f in "$COV" "$PRS" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
status=0
Rscript "$SCRIPT" --raw_dir "$RAW" --covar "$COV" \
  --dose_lpa "$W/dose_LPA.raw" --w_lpa "$W/lpa_score.txt" --dose_ldl "$W/dose_LDL.raw" --w_ldl "$W/ldl_score.txt" \
  --prs "$PRS" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}" || status=$?
echo "exit status $status; end $(date)"
exit "$status"
