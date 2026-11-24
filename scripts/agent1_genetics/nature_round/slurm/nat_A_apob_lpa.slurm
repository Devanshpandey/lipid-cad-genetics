#!/bin/bash
#SBATCH -J natA_apobl
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/natA_apob_lpa_%j.out
#SBATCH -e /path/to/cad-genetics/logs/natA_apob_lpa_%j.err
# Nature round, analysis A: ApoB + Lp(a) multivariable MR. Readout rules are in
# scripts/agent1_genetics/nature_round/PREREGISTERED_READOUT.md section A.
#
# Arms, in order:
#   statinfix    primary. Enforces the positive control: exit status 2 and no ApoB estimate
#                if the LDL-C analogue does not reproduce OR 1.25 / 1.15 on CAD.
#   statinnaive  LDL_C, APOB, TRIGLY from merged_statinnaive. LPA and HDL_C are not there and
#                come from merged_statinfix; the script flags that in VERDICT.txt.
#   lpacens      censored Lp(a) GWAS (merged_lpacens holds only lipids_LPA). Skipped, loudly,
#                until both the GWAS and its phenotype file exist. When they do, run only it:
#                  ARMS=lpacens sbatch --export=ALL scripts/agent1_genetics/nature_round/slurm/nat_A_apob_lpa.slurm
#
# Sensitivity arms refuse to start without the stamp written by a passing statinfix arm.
# The selftest runs first on this R installation and a failure stops the job. Skip it only
# when resubmitting an arm: SKIP_SELFTEST=1.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# Parallelism is by fork, one gunzip|awk per process. Threaded BLAS would oversubscribe.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ARMS=${ARMS:-statinfix statinnaive lpacens}
SKIP_SELFTEST=${SKIP_SELFTEST:-0}
CORES=${SLURM_NTASKS:-16}

ND=$_REPO/scripts/agent1_genetics/nature_round
SCRIPT=$ND/nat_apob_lpa_mvmr.R
GW=$_REPO/results/agent1_genetics/gwas_summary_stats
FIX=$GW/merged_statinfix
NAIVE=$GW/merged_statinnaive
CENS=$GW/merged_lpacens
PHENO=$_REPO/results/agent1_genetics/phenotypes/pheno_quantitative_statinfix.txt
# The statin-naive GWAS used the unrecoded phenotype on a keep list (08g), so its correlation must too.
PHENO_NAIVE=$_REPO/results/agent1_genetics/phenotypes/pheno_quantitative.txt
KEEP_NAIVE=$_REPO/results/agent1_genetics/phenotypes/regenie_keep_statinnaive.txt
PHENO_CENS=$_REPO/results/agent1_genetics/phenotypes/lpa_censored/pheno_quantitative_statinfix_lpacens.txt
FINNGEN=$_REPO/data/external/finngen/finngen_R11_I9_CHD.gz
CARDIOGRAM=$_REPO/data/external/cardiogram/cad.add.160614.website.txt
CIS=/path/to/scratch/statint/lpa_cis_ids.txt
PUB_IV=$_REPO/exports/agent1/mvmr_prunefix/mvmr_CAD.instruments.tsv
OUT=$_REPO/exports/agent1/nature_round/apob_lpa
CACHE=$SCRATCH/natA_apob_lpa_cache
STAMP=$OUT/statinfix/positive_control.PASSED

missing=0
for f in "$SCRIPT" "$ND/nat_apob_lpa_mvmr_selftest.R" "$PHENO" "$PHENO_NAIVE" "$KEEP_NAIVE" "$FINNGEN" "$CARDIOGRAM" "$CIS" \
         "$_REPO/scripts/agent1_genetics/R/run_mvmr.R" "$_REPO/scripts/agent1_genetics/R/mvmr_covariance.R" \
         "$_REPO/scripts/agent1_genetics/R/run_external_mr.R"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; missing=1; }
done
for t in lipids_LDL_C lipids_APOB lipids_HDL_C lipids_TRIGLY lipids_LPA \
         outcomes_CAD outcomes_MI outcomes_REVASC outcomes_STROKE outcomes_HF outcomes_CV_DEATH outcomes_MACE outcomes_PAD; do
  [ -s "$FIX/$t.txt.gz" ] || { echo "MISSING INPUT: $FIX/$t.txt.gz"; missing=1; }
done
for t in lipids_LDL_C lipids_APOB lipids_TRIGLY; do
  [ -s "$NAIVE/$t.txt.gz" ] || { echo "MISSING INPUT: $NAIVE/$t.txt.gz"; missing=1; }
done
[ "$missing" -eq 0 ] || { echo "Stopping before R: inputs missing"; exit 1; }
mkdir -p "$OUT" "$CACHE" "$_REPO/logs"

echo "start $(date) on $(hostname); arms: $ARMS"
Rscript -e 'cat("R", as.character(getRversion()), "MVMR", as.character(packageVersion("MVMR")), "MendelianRandomization", as.character(packageVersion("MendelianRandomization")), "data.table", as.character(packageVersion("data.table")), "\n")'

if [ "$SKIP_SELFTEST" != "1" ]; then
  echo "=== selftest (synthetic data only) ==="
  # The selftest empties its --workdir, so this path is its own and nothing else lives there.
  if ! Rscript "$ND/nat_apob_lpa_mvmr_selftest.R" --workdir "$SCRATCH/natA_apob_lpa_selftest" --threads "$CORES"; then
    echo "SELFTEST FAILED: nothing was run on real data."
    exit 1
  fi
fi

run_arm() {
  local arm=$1; shift
  echo "=== arm $arm $(date) ==="
  Rscript "$SCRIPT" --arm "$arm" --threads "$CORES" --cache_dir "$CACHE" \
    --outdir "$OUT/$arm" --cis_list "$CIS" --finngen "$FINNGEN" --cardiogram "$CARDIOGRAM" --nboot 2000 "$@"
}

overall=0
for arm in $ARMS; do
  status=0
  case "$arm" in
    statinfix)
      run_arm statinfix --arm_type primary --sumstat_dir "$FIX" --pheno "$PHENO" --pc_instruments "$PUB_IV" || status=$?
      if [ "$status" -eq 2 ]; then
        echo "POSITIVE CONTROL FAILED: see $OUT/statinfix/VERDICT.txt. No estimate is to be read; stopping."
        exit 2
      elif [ "$status" -ne 0 ]; then
        echo "statinfix arm FAILED with status $status; sensitivity arms cannot run without it."
        exit "$status"
      fi
      ;;
    statinnaive)
      run_arm statinnaive --arm_type sensitivity --sumstat_dir "$NAIVE" --fallback_dir "$FIX" \
        --pheno "$PHENO_NAIVE" --keep "$KEEP_NAIVE" --control_stamp "$STAMP" || status=$?
      ;;
    lpacens)
      if [ ! -s "$CENS/lipids_LPA.txt.gz" ] || [ ! -s "$PHENO_CENS" ]; then
        echo "SKIPPED lpacens: $CENS/lipids_LPA.txt.gz or $PHENO_CENS does not exist yet. The readout lists this arm as 'when available', and it has not been read."
        continue
      fi
      run_arm lpacens --arm_type sensitivity --sumstat_dir "$CENS" --fallback_dir "$FIX" \
        --pheno "$PHENO_CENS" --control_stamp "$STAMP" || status=$?
      ;;
    *) echo "Unknown arm: $arm"; exit 1 ;;
  esac
  if [ "$status" -ne 0 ]; then
    echo "arm $arm FAILED with status $status"
    overall=$status
  fi
done

echo "=== verdicts ==="
for arm in $ARMS; do
  if [ -s "$OUT/$arm/VERDICT.txt" ]; then
    echo "--- $arm"; cat "$OUT/$arm/VERDICT.txt"
  fi
done
[ "$overall" -eq 0 ] && echo "ALL_DONE $(date)"
exit "$overall"
