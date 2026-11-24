#!/bin/bash
#SBATCH -J lpaSus
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaSus_%x_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaSus_%x_%j.err
# Lp(a) by arterial susceptibility (see PREREGISTRATION.md). STEP=selftest (default; synthetic
# only) or run. Reused code must be unchanged: lpa_horizon.R 6e1c268df17113e7959d4078ee830b64,
# agent6 lib_synergy.R bfad4277288f93c091769bfb7c447b61, lib_resilience.R
# e3859ff7aae063e10a57f402d7b1bb52, config_synergy.yaml a80fe645b0891b221b1cca4b60cfeb44.
set -euo pipefail
umask 077
STEP=${STEP:-selftest}
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
DIR=$_REPO/scripts/agent1_genetics/lpa_susceptibility
A6=$_REPO/scripts/agent6_resilience
SCRIPT=$DIR/lpa_susceptibility.R
chk() { local got; got=$(md5sum "$1" | cut -d' ' -f1); [ "$got" = "$2" ] || { echo "CHANGED: $1 ($got)"; exit 1; }; }
chk "$_REPO/scripts/agent1_genetics/lpa_horizon/lpa_horizon.R" 6e1c268df17113e7959d4078ee830b64
chk "$A6/R/lib_synergy.R" bfad4277288f93c091769bfb7c447b61
chk "$A6/R/lib_resilience.R" e3859ff7aae063e10a57f402d7b1bb52
chk "$A6/config_synergy.yaml" a80fe645b0891b221b1cca4b60cfeb44
echo "step $STEP; start $(date) on $(hostname); script md5 $(md5sum "$SCRIPT" | cut -d' ' -f1)"
Rscript -e 'invisible(parse(commandArgs(TRUE)[1])); cat("parse ok\n")' "$SCRIPT"

if [ "$STEP" = "selftest" ]; then
  WORK=$(mktemp -d "${SCRATCH}/lpaSus.XXXXXX")
  Rscript "$SCRIPT" --step selftest --workdir "$WORK" --threads "${SLURM_NTASKS:-8}"
  rm -rf "$WORK"
  echo "end $(date)"; exit 0
fi

: "${COMMIT:?pass COMMIT=<hash of the commit holding PREREGISTRATION.md and this script>}"
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
RUN=$_REPO/results/agent6_resilience/run_cad_ukbfree_v2
PRS=$RUN/10_synergy/prs_components.tsv
OUT=$_REPO/exports/agent1/lpa_susceptibility
for f in "$COV" "$PRS" "$RUN/01/cohort.tsv" "$RUN/04/mrs_oof.tsv" "$RUN/03/nmr_technical_resid.tsv" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
{ echo "commit $COMMIT"; echo "job ${SLURM_JOB_ID:-none}"; echo "date $(date)"; md5sum "$SCRIPT" "$A6/R/lib_synergy.R" "$A6/R/lib_resilience.R" "$_REPO/scripts/agent1_genetics/lpa_horizon/lpa_horizon.R"; } > "$OUT/PROVENANCE.txt"
status=0
Rscript "$SCRIPT" --step run --run_dir "$RUN" --components "$PRS" --raw_dir "$RAW" --covar "$COV" \
  --dose_lpa "$W/dose_LPA.raw" --w_lpa "$W/lpa_score.txt" --dose_ldl "$W/dose_LDL.raw" --w_ldl "$W/ldl_score.txt" \
  --prs "$PRS" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}" || status=$?
echo "exit status $status; end $(date)"
exit "$status"
