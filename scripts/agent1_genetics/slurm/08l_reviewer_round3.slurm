#!/bin/bash
#SBATCH -J rev3
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 04:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/rev3_%j.out
#SBATCH -e /path/to/cad-genetics/rev3_%j.err
# ============================================================================
# Round-3 reviewer items. Blocks are independent; a failing block warns and the
# job continues, so one bad input cannot cost the whole run.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
G=$_REPO/results/agent1_genetics/gwas_summary_stats
S=$G/merged_statinfix
SN=$G/merged_statinnaive
OUT=$_REPO/exports/agent1/revision
mkdir -p $OUT $SN

echo "############ [0] merge the statin-naive per-chromosome GWAS"
for TR in LDL_C APOB TRIGLY; do
  DST=$SN/lipids_${TR}.txt.gz
  if [ -s "$DST" ]; then echo "  $TR already merged"; continue; fi
  FIRST=$G/lipids_statinnaive/lipids_chr1_${TR}.regenie.gz
  if [ ! -s "$FIRST" ]; then echo "  MISSING $FIRST"; continue; fi
  ( zcat $FIRST | head -1
    for C in $(seq 1 22); do
      F=$G/lipids_statinnaive/lipids_chr${C}_${TR}.regenie.gz
      [ -s "$F" ] && zcat $F | tail -n +2
    done ) | gzip -c > $DST
  echo "  $TR -> $(zcat $DST | wc -l) rows"
done
# outcomes come from the corrected set; only the EXPOSURE differs in this arm
for O in CAD MI MACE; do ln -sf $S/outcomes_${O}.txt.gz $SN/outcomes_${O}.txt.gz; done

echo "############ [1] statin-naive LDL-C / ApoB joint MVMR + conditional F"
for T in CAD MI MACE; do
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr_apob.R \
    --sumstat_dir $SN --outcome $T \
    --out $OUT/mvmr_apob_statinnaive_${T}.csv || echo "  BLOCK1 $T FAILED"
done

echo "############ [2] between-subtype heterogeneity: omnibus Q and LAS contrasts"
Rscript $_REPO/scripts/agent1_genetics/R/stroke_subtype_contrasts.R \
  --sumstat_dir $S \
  --megastroke_dir $_REPO/data/external/megastroke \
  --nboot 1000 \
  --out_prefix $OUT/stroke_subtype_contrasts || echo "  BLOCK2 FAILED"

echo ALL_DONE
