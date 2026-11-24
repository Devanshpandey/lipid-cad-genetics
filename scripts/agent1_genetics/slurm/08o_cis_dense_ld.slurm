#!/bin/bash
#SBATCH -J cisdense
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 03:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/cisdense_%j.out
#SBATCH -e /path/to/cad-genetics/cisdense_%j.err
# Build a DENSE cis-window LD panel and rerun the drug-target cis-MR against it.
# lpa_ldref is a genome-wide-significant UNION subset and is ~2x too sparse where
# instruments cluster (749 vs 1,640 variants in the LPA cis window), which is why
# it recovered only 4 of the published 11 LPA instruments.
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
PLINK2=/path/to/software/plink2
PLINK=/path/to/software/plink
[ -x "$PLINK2" ] || PLINK2=$(which plink2 2>/dev/null || echo "$PLINK")
W=$SCRATCH/cis_ldref
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
OUT=$_REPO/exports/agent1/revision

# cis windows: gene body +/- 100 kb, GRCh37
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
echo "=== extracting dense cis panel from the full UKB genotypes ==="
$PLINK2 --bfile $UKB_GENO_STEP2 \
  --extract bed1 $SCRATCH/cis_windows.txt \
  --max-alleles 2 --snps-only \
  --thin-indiv-count 20000 \
  --make-bed --out $W --threads 32 --memory 120000 || { echo "EXTRACT FAILED"; exit 1; }
echo "  panel: $(wc -l < $W.bim) variants, $(wc -l < $W.fam) samples"
for G in PCSK9 LPA CETP APOC3; do :; done

for ARM in old new; do
  if [ "$ARM" = old ]; then LDL=../merged/lipids_LDL_C.txt.gz; LAB=old_6153only; else LDL=lipids_LDL_C.txt.gz; LAB=corrected_union; fi
  echo "=== ARM $LAB ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_cis_mr.R \
    --sumstat_dir $S --outcome_dir $S --ldl_file $LDL \
    --ld_ref $W --plink $PLINK --label $LAB \
    --out $OUT/cis_mr_dense_${LAB}.csv || echo "  ARM $LAB FAILED"
done
echo ALL_DONE
