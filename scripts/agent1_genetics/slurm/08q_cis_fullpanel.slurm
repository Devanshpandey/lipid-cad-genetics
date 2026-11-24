#!/bin/bash
#SBATCH -J cisfull
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 03:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/cisfull_%j.out
#SBATCH -e /path/to/cad-genetics/cisfull_%j.err
# ============================================================================
# WHY THIS RERUN EXISTS.
#   08o built the cis LD panel with `--thin-indiv-count 20000`. UKB_GENO_STEP2 is
#   filtered at MAF>=0.001, so at 20,000 individuals the rarest retained variants
#   carry ~40 minor alleles. r2 estimated on ~40 carriers is noise, and noise is
#   almost always below the r2<0.001 clumping threshold -- so rare variants survive
#   clumping as "independent" instruments. That is why LPA returned 80 instruments
#   from a 335 kb window at r2<0.001 while the subtype-contrast analysis returned
#   only 45 from a 1.2 Mb window at the LOOSER r2<0.01. The difference is the panel,
#   not the biology, and 80 is an artifact.
#
#   Fix: drop the thinning. Clump against all 402,199 participants -- the same panel
#   the cis-LPA subtype-contrast instrument already uses, so the two analyses in the
#   manuscript stop describing two different reference panels for the same locus.
#
# VALIDATION. The old-GWAS arm should approach the published counts (LPA 11, CETP 8).
#   If LPA stays near 80 with the full panel, the rare-variant explanation is wrong
#   and nothing here should be written up.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
PLINK2=/path/to/software/plink2
PLINK=/path/to/software/plink
W=$SCRATCH/cis_ldref_full
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
OUT=$_REPO/exports/agent1/revision
mkdir -p $OUT

cat > $SCRATCH/cis_windows.txt <<'BED'
1	55405221	55630525	PCSK9
1	62963191	63171984	ANGPTL3
5	74532154	74757929	HMGCR
7	44452134	44680929	NPC1L1
11	116600623	116803788	APOC3
6	160852515	161187407	LPA
16	56895762	57117757	CETP
19	11100038	11344506	LDLR
BED

echo "=== extracting cis panel at FULL sample size (no --thin-indiv-count) ==="
$PLINK2 --bfile $UKB_GENO_STEP2 \
  --extract bed1 $SCRATCH/cis_windows.txt \
  --max-alleles 2 --snps-only \
  --make-bed --out $W --threads 32 --memory 120000 || { echo "EXTRACT FAILED"; exit 1; }
echo "  panel: $(wc -l < $W.bim) variants, $(wc -l < $W.fam) participants"

# per-variant MAF, so the rare-variant explanation is checkable rather than asserted
$PLINK2 --bfile $W --freq --out $W.freq --threads 32 >/dev/null 2>&1

for ARM in old new; do
  if [ "$ARM" = old ]; then LDL=../merged/lipids_LDL_C.txt.gz; LAB=old_6153only; else LDL=lipids_LDL_C.txt.gz; LAB=corrected_union; fi
  echo "=== ARM $LAB (full panel) ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_cis_mr.R \
    --sumstat_dir $S --outcome_dir $S --ldl_file $LDL \
    --ld_ref $W --plink $PLINK --label $LAB \
    --out $OUT/cis_mr_fullpanel_${LAB}.csv || echo "  ARM $LAB FAILED"
done

# Re-derive the cis-LPA subtype-contrast instrument against the SAME panel object so
# the manuscript can name one panel. 45 at r2<0.01 is the published value; if this
# reproduces it, the two analyses are finally on identical footing.
echo "=== cis-LPA contrast instrument, same panel ==="
Rscript $_REPO/scripts/agent1_genetics/R/stroke_subtype_contrasts.R \
  --sumstat_dir $S \
  --megastroke_dir $_REPO/data/external/megastroke \
  --ld_ref $W --plink $PLINK \
  --nboot 1000 \
  --out_prefix $OUT/stroke_subtype_contrasts_fullpanel || echo "CONTRAST FAILED"
echo ALL_DONE
