#!/bin/bash
#SBATCH -J natB_sel
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/natB_selection_%j.out
#SBATCH -e /path/to/cad-genetics/logs/natB_selection_%j.err
# Nature round, analysis B: selection into statin use. Readout rules are in
# scripts/agent1_genetics/nature_round/PREREGISTERED_READOUT.md section B.
#
# The R script stops with exit status 2 if the published interactions do not reproduce.
# A resubmission resumes the simulation from sim_cells/ when nothing has changed, so if the
# 12 h limit is hit, sbatch this same file again; verdict.csv only appears when all cells exist.
# vm-small has about 25 GB across 16 cores; see --cores below before raising --n_reps.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# glmnet lives here and drives only a sensitivity row; the script logs it if absent. It goes in
# through --glmnet_lib, searched last. R_LIBS_USER would search it first, and any survival or
# data.table in there would then replace the conda env's under the positive control.
GLMNET_LIB=/path/to/R/library
# Parallelism is by fork, one fit per process. Threaded BLAS would oversubscribe the 16 cores.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
OUT=$_REPO/exports/agent1/nature_round/statin_selection
SCRIPT=$_REPO/scripts/agent1_genetics/nature_round/nat_statin_selection.R

missing=0
for f in "$SCRIPT" "$COV" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; missing=1; }
done
for fid in 53 40000 131296 131298 131300 131302 131304 131306 6153 6177 30780 30790 30786 30796 26252; do
  [ -s "$RAW/fid${fid}.csv" ] || { echo "MISSING INPUT: $RAW/fid${fid}.csv"; missing=1; }
done
[ "$missing" -eq 0 ] || { echo "Stopping before R: inputs missing (statint scores: rerun 08e's score build)"; exit 1; }
# Fail here, not hours into the queue. glmnet is optional by design; the rest are not.
Rscript -e 'for (p in c("data.table","optparse","survival","splines","parallel")) if (!requireNamespace(p, quietly = TRUE)) { message("MISSING R PACKAGE: ", p); quit(save = "no", status = 1L) }; cat("R", as.character(getRversion()), "survival", as.character(packageVersion("survival")), "from", find.package("survival"), "\n")'
mkdir -p "$OUT" "$_REPO/logs"

echo "start $(date) on $(hostname)"
status=0
Rscript "$SCRIPT" \
  --scores  "$W/dose_LPA.raw,$W/dose_LDL.raw" \
  --weights "$W/lpa_score.txt,$W/ldl_score.txt" \
  --raw_dir "$RAW" \
  --covar   "$COV" \
  --outdir  "$OUT" \
  --glmnet_lib "$GLMNET_LIB" \
  --cores   "${SLURM_NTASKS:-16}" \
  --n_reps  200 \
  --n_reps_secondary 50 || status=$?
if [ "$status" -eq 2 ]; then
  echo "POSITIVE CONTROL FAILED: see $OUT/positive_control.csv. No estimate is to be read."
elif [ "$status" -ne 0 ]; then
  echo "FAILED with status $status"
else
  echo "ALL_DONE $(date)"
  cat "$OUT/verdict.csv"
fi
exit "$status"
