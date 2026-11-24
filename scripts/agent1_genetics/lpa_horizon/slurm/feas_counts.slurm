#!/bin/bash
#SBATCH -J lpaH_feas
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaH_feas_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaH_feas_%j.err
# Feasibility counts for the Lp(a)HORIZON genetic emulation. Counts and the exposure
# distribution only: no outcome is tabulated against Lp(a) or LPA genotype.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
DOSE=/path/to/scratch/statint/dose_LPA.raw
OUT=$_REPO/exports/agent1/lpa_horizon/feasibility
for f in "$COV" "$DOSE"; do [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }; done
mkdir -p "$OUT" "$_REPO/logs"
echo "start $(date) on $(hostname)"
Rscript "$_REPO/scripts/agent1_genetics/lpa_horizon/feas_counts.R" \
  --raw_dir "$RAW" --covar "$COV" --dose "$DOSE" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}"
echo "end $(date)"
