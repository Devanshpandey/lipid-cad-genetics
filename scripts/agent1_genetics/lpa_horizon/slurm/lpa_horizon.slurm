#!/bin/bash
#SBATCH -J lpaH
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaH_%x_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaH_%x_%j.err
# Genetic emulation of Lp(a)HORIZON. STEP=selftest (default; synthetic data only), gates or test.
# Rules: scripts/agent1_genetics/lpa_horizon/PREREGISTERED_READOUT.md. The real-data steps
# refuse to run unless that file is committed in the worktree this code was copied from
# (COMMIT must be passed and is written to the output folder).
set -euo pipefail
umask 077
STEP=${STEP:-selftest}
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SCRIPT=$_REPO/scripts/agent1_genetics/lpa_horizon/lpa_horizon.R
echo "step $STEP; start $(date) on $(hostname); script md5 $(md5sum "$SCRIPT" | cut -d' ' -f1)"
Rscript -e 'invisible(parse(commandArgs(TRUE)[1])); cat("parse ok\n")' "$SCRIPT"

if [ "$STEP" = "selftest" ]; then
  WORK=$(mktemp -d "${SCRATCH}/lpaH_selftest.XXXXXX")
  Rscript "$SCRIPT" --step selftest --workdir "$WORK" --threads "${SLURM_NTASKS:-8}"
  rm -rf "$WORK"
  echo "end $(date)"; exit 0
fi

: "${COMMIT:?pass COMMIT=<hash of the commit holding the readout rules and this script>}"
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
PRS=$_REPO/results/agent6_resilience/run_cad_ukbfree_v2/10_synergy/prs_components.tsv
OUT=$_REPO/exports/agent1/lpa_horizon/emulation
for f in "$COV" "$PRS" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
{ echo "commit $COMMIT"; echo "job ${SLURM_JOB_ID:-none}"; echo "step $STEP"; echo "date $(date)"; md5sum "$SCRIPT"; } > "$OUT/PROVENANCE_${STEP}.txt"
status=0
Rscript "$SCRIPT" --step "$STEP" --raw_dir "$RAW" --covar "$COV" \
  --dose_lpa "$W/dose_LPA.raw" --w_lpa "$W/lpa_score.txt" --dose_ldl "$W/dose_LDL.raw" --w_ldl "$W/ldl_score.txt" \
  --prs "$PRS" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}" || status=$?
echo "exit status $status; end $(date)"
exit "$status"
