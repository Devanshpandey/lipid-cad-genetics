#!/bin/bash
#SBATCH -J lpaHp3
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaHp3_%x_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaHp3_%x_%j.err
# MATCHED-ENDPOINT analysis for the genetic comparison with Lp(a)HORIZON (PREREGISTERED_READOUT.md,
# section "Matched-endpoint analysis"). STEP=selftest (default; synthetic only) or run. No earlier
# script is edited: lpa_horizon.R must stay 6e1c268df17113e7959d4078ee830b64,
# lpa_horizon_posthoc.R 03a7ee5e873ab0800c34cac935f92a0e, lpa_horizon_posthoc2.R
# 408990dd06369cb0490f552f327b7b81 and lpa_forecast.R af96702bf5b9dcbba0bea8a6adb527a7.
set -euo pipefail
umask 077
STEP=${STEP:-selftest}
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
DIR=$_REPO/scripts/agent1_genetics/lpa_horizon
SCRIPT=$DIR/lpa_horizon_posthoc3.R
echo "step $STEP; start $(date) on $(hostname); script md5 $(md5sum "$SCRIPT" | cut -d' ' -f1)"
want() { [ "$(md5sum "$DIR/$1" | cut -d' ' -f1)" = "$2" ] || { echo "$1 has changed"; exit 1; }; }
want lpa_horizon.R 6e1c268df17113e7959d4078ee830b64
want lpa_horizon_posthoc.R 03a7ee5e873ab0800c34cac935f92a0e
want lpa_horizon_posthoc2.R 408990dd06369cb0490f552f327b7b81
want lpa_forecast.R af96702bf5b9dcbba0bea8a6adb527a7
Rscript -e 'invisible(parse(commandArgs(TRUE)[1])); cat("parse ok\n")' "$SCRIPT"

if [ "$STEP" = "selftest" ]; then
  WORK=$(mktemp -d "${SCRATCH}/lpaH_matched.XXXXXX")
  Rscript "$SCRIPT" --step selftest --workdir "$WORK" --threads "${SLURM_NTASKS:-8}"
  rm -rf "$WORK"
  echo "end $(date)"; exit 0
fi

: "${COMMIT:?pass COMMIT=<hash of the commit holding the declaration in PREREGISTERED_READOUT.md and this script>}"
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
PRS=$_REPO/results/agent6_resilience/run_cad_ukbfree_v2/10_synergy/prs_components.tsv
OUT=$_REPO/exports/agent1/lpa_horizon/matched_endpoint
for f in "$COV" "$PRS" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
{ echo "commit $COMMIT"; echo "job ${SLURM_JOB_ID:-none}"; echo "date $(date)"; md5sum "$SCRIPT" "$DIR/lpa_forecast.R" "$DIR/lpa_horizon_posthoc2.R" "$DIR/lpa_horizon_posthoc.R" "$DIR/lpa_horizon.R"; } > "$OUT/PROVENANCE.txt"
status=0
Rscript "$SCRIPT" --step run --raw_dir "$RAW" --covar "$COV" \
  --dose_lpa "$W/dose_LPA.raw" --w_lpa "$W/lpa_score.txt" --dose_ldl "$W/dose_LDL.raw" --w_ldl "$W/ldl_score.txt" \
  --prs "$PRS" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}" || status=$?
echo "exit status $status; end $(date)"
exit "$status"
