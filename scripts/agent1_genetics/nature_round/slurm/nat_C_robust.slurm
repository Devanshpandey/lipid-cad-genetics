#!/bin/bash
#SBATCH -J natC_robust
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/natC_robust_%j.out
#SBATCH -e /path/to/cad-genetics/logs/natC_robust_%j.err
# Nature round, analysis C: weak-instrument-robust inference on the LDL-C minus ApoB contrast.
# Readout rules are in scripts/agent1_genetics/nature_round/PREREGISTERED_READOUT.md section C.
#
# STEPS picks what runs; the default is everything, in the only order that can give a result:
#   selftest    synthetic data only; a failure stops the job
#   calibrate   simulated collinear exposures; writes calibration_status.txt. A gated scenario out
#               of 93-97% coverage stops the job, and no real number is produced.
#   real        refuses to start unless calibration_status.txt says PASSED for this exact script
#               (md5), so any edit to the R script means calibrating again.
# Resubmit only the real step, calibration kept, after adding an arm:
#   STEPS=real sbatch --export=ALL scripts/agent1_genetics/nature_round/slurm/nat_C_robust.slurm
# Secondary arms whose files are absent are skipped with a message and named in the verdict.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# Parallelism is by fork; threaded BLAS would oversubscribe the 16 cores.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

STEPS=${STEPS:-selftest calibrate real}
# One core left for the parent, which does the extracts and collects the forks.
CORES=$(( ${SLURM_NTASKS:-16} - 1 ))

ND=$_REPO/scripts/agent1_genetics/nature_round
SCRIPT=$ND/nat_ldl_apob_robust.R
PH=$_REPO/results/agent1_genetics/phenotypes
GW=$_REPO/results/agent1_genetics/gwas_summary_stats
OUT=$_REPO/exports/agent1/nature_round/robust_contrast
# Instrument rows are summary-level but large; they stay out of the export folder.
CACHE=$SCRATCH/natC_robust_cache
ST=$SCRATCH/natC_robust_selftest
mkdir -p "$OUT" "$CACHE"

# What the primary analysis cannot run without. Statin-naive and NMR inputs may be absent.
missing=0
for f in "$SCRIPT" "$GW/merged_statinfix/lipids_LDL_C.txt.gz" "$GW/merged_statinfix/lipids_APOB.txt.gz" \
         "$GW/merged_statinfix/lipids_TRIGLY.txt.gz" "$GW/merged_statinfix/outcomes_CAD.txt.gz" \
         "$PH/pheno_quantitative_statinfix.txt" "$_REPO/data/external/finngen/finngen_R11_I9_CHD.gz"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; missing=1; }
done
[ "$missing" = 0 ] || { echo "Primary inputs missing; stopping before any step."; exit 1; }
Rscript -e 'for (p in c("data.table", "optparse", "MendelianRandomization", "MVMR", "boot")) { if (!requireNamespace(p, quietly = TRUE)) stop("R package missing: ", p); cat(p, as.character(packageVersion(p)), "\n") }'

echo "[$(date)] script md5 $(md5sum "$SCRIPT" | cut -d' ' -f1), steps: $STEPS, cores: $CORES"
for step in $STEPS; do
  echo "[$(date)] ===== $step ====="
  case "$step" in
    selftest)
      rm -rf "$ST"
      # 1000 replicates give MC SE 0.007, and 0.97 is then crossed by noise (AR 0.972, then 0.948).
      Rscript "$SCRIPT" --selftest --selftest_dir "$ST" --ncores "$CORES" --selftest_reps 10000
      ;;
    calibrate)
      Rscript "$SCRIPT" --calibrate --calibration_dir "$OUT" --ncores "$CORES" --n_reps 10000
      ;;
    real)
      Rscript "$SCRIPT" --repo "$_REPO" --out_dir "$OUT" --calibration_dir "$OUT" --cache_dir "$CACHE" --ncores "$CORES" \
        --nmr_dir /path/to/archive/results/agent1_genetics/gwas_summary_stats/nmr
      ;;
    *)
      echo "unknown step: $step"; exit 1
      ;;
  esac
done
echo "[$(date)] done; outputs in $OUT"
case " $STEPS " in *" real "*) cat "$OUT/robust_contrast_verdict.txt" ;; esac
exit 0
