#!/bin/bash
#SBATCH -J lpaHurg
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:40:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaHurg_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaHurg_%j.err
# COUNTS ONLY: how many coronary revascularizations can be linked to a newly recorded acute
# coronary code, as a stand-in for the urgent procedures Lp(a)HORIZON counted. Fits no outcome
# model and tabulates nothing against Lp(a) or genotype. The script specified in advance is
# sourced, not edited, and must stay 6e1c268df17113e7959d4078ee830b64.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
DIR=$_REPO/scripts/agent1_genetics/lpa_horizon
SCRIPT=$DIR/urgent_feas_counts.R
BASE_MD5=$(md5sum "$DIR/lpa_horizon.R" | cut -d' ' -f1)
echo "start $(date) on $(hostname); counts md5 $(md5sum "$SCRIPT" | cut -d' ' -f1); base md5 $BASE_MD5"
[ "$BASE_MD5" = "6e1c268df17113e7959d4078ee830b64" ] || { echo "the base script has changed"; exit 1; }
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
PRS=$_REPO/results/agent6_resilience/run_cad_ukbfree_v2/10_synergy/prs_components.tsv
OUT=$_REPO/exports/agent1/lpa_horizon/urgent_feasibility
for f in "$COV" "$PRS" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
status=0
Rscript "$SCRIPT" --raw_dir "$RAW" --covar "$COV" \
  --dose_lpa "$W/dose_LPA.raw" --w_lpa "$W/lpa_score.txt" --dose_ldl "$W/dose_LDL.raw" --w_ldl "$W/ldl_score.txt" \
  --prs "$PRS" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}" || status=$?
printf 'job %s\ndate %s\ncounts_md5 %s\nbase_md5 %s\nstatus %s\n' "${SLURM_JOB_ID:-none}" "$(date +%F)" "$(md5sum "$SCRIPT" | cut -d' ' -f1)" "$BASE_MD5" "$status" > "$OUT/PROVENANCE.txt"
echo "exit status $status; end $(date)"
exit "$status"
