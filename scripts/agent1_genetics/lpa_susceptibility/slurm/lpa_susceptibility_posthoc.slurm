#!/bin/bash
#SBATCH -J lpaSusPH
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpaSusPH_%x_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpaSusPH_%x_%j.err
# POST HOC rerun with the LPA locus removed from the lipid-region score (PREREGISTRATION.md,
# "Post hoc analysis"). STEP=selftest (default; synthetic only) or run. Reused code must be
# unchanged; the pre-registered script must stay bdd4cd26b171756c65b104a53d063ec2.
set -euo pipefail
umask 077
STEP=${STEP:-selftest}
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
DIR=$_REPO/scripts/agent1_genetics/lpa_susceptibility
A6=$_REPO/scripts/agent6_resilience
SCRIPT=$DIR/lpa_susceptibility_posthoc.R
chk() { local got; got=$(md5sum "$1" | cut -d' ' -f1); [ "$got" = "$2" ] || { echo "CHANGED: $1 ($got)"; exit 1; }; }
chk "$DIR/lpa_susceptibility.R" bdd4cd26b171756c65b104a53d063ec2
chk "$_REPO/scripts/agent1_genetics/lpa_horizon/lpa_horizon.R" 6e1c268df17113e7959d4078ee830b64
chk "$A6/R/lib_synergy.R" bfad4277288f93c091769bfb7c447b61
chk "$A6/R/lib_resilience.R" e3859ff7aae063e10a57f402d7b1bb52
chk "$A6/config_synergy.yaml" a80fe645b0891b221b1cca4b60cfeb44
echo "step $STEP; start $(date) on $(hostname); script md5 $(md5sum "$SCRIPT" | cut -d' ' -f1)"
Rscript -e 'invisible(parse(commandArgs(TRUE)[1])); cat("parse ok\n")' "$SCRIPT"
WORK=$(mktemp -d "${SCRATCH}/lpaSusPH.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

if [ "$STEP" = "selftest" ]; then
  Rscript "$SCRIPT" --step selftest --workdir "$WORK" --threads "${SLURM_NTASKS:-8}"
  echo "end $(date)"; exit 0
fi

: "${COMMIT:?pass COMMIT=<hash of the commit holding the post hoc declaration and this script>}"
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
RUN=$_REPO/results/agent6_resilience/run_cad_ukbfree_v2
PRS=$RUN/10_synergy/prs_components.tsv
ALLW=$_REPO/results/agent6_resilience/prs/cad_ukbfree_all_chr.txt
SUBSET=${TACC_SCRATCH}/agent6_prs_subset/ukb400k_prs_snps
REFCSV=$_REPO/exports/agent1/lpa_susceptibility/product_terms.csv
OUT=$_REPO/exports/agent1/lpa_susceptibility/posthoc_lpa_removed
for f in "$COV" "$PRS" "$ALLW" "$REFCSV" "$RUN/01/cohort.tsv" "$RUN/10_synergy/prs_component_weights.tsv" "$RUN/10_synergy/lipid_regions.tsv" \
         "$SUBSET.bed" "$UKB_SAMPLE_EID" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
{ echo "commit $COMMIT"; echo "job ${SLURM_JOB_ID:-none}"; echo "date $(date)"; md5sum "$SCRIPT" "$DIR/lpa_susceptibility.R"; } > "$OUT/PROVENANCE.txt"

echo "[$(date)] LPA-region weights"
Rscript "$SCRIPT" --step weights --weights_all "$ALLW" --component_weights "$RUN/10_synergy/prs_component_weights.tsv" \
  --regions "$RUN/10_synergy/lipid_regions.tsv" --outdir "$WORK"
cp "$WORK/lpa_region_weights.tsv" "$OUT/lpa_region_weights.tsv"     # variant weights only, no individual data
echo "[$(date)] score the LPA region (same options as the original components)"
"${PLINK2}" --bfile "$SUBSET" --keep "$UKB_SAMPLE_EID" --extract "$WORK/lpa_region.snps" \
  --score "$WORK/lpa_region_weights.tsv" 1 2 header-read cols=+scoresums --score-col-nums 3 \
  --threads 8 --memory 20000 --out "$WORK/lpa_region" > "$WORK/plink.out" 2>&1 || { tail -20 "$WORK/plink.out"; exit 1; }
grep -E "variants|samples|--score" "$WORK/plink.out" | tail -6
status=0
Rscript "$SCRIPT" --step run --lpa_sscore "$WORK/lpa_region.sscore" --reference "$REFCSV" --run_dir "$RUN" --components "$PRS" \
  --raw_dir "$RAW" --covar "$COV" --dose_lpa "$W/dose_LPA.raw" --w_lpa "$W/lpa_score.txt" --dose_ldl "$W/dose_LDL.raw" --w_ldl "$W/ldl_score.txt" \
  --prs "$PRS" --outdir "$OUT" --threads "${SLURM_NTASKS:-8}" || status=$?
echo "exit status $status; end $(date)"
exit "$status"
