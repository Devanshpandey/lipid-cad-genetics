#!/bin/bash
# ============================================================================
# DRIVER — Multi-ancestry causal MR: GLGC-{ANC} lipid -> MVP-{ANC} CAD
# GLGC lipid exposures (AFR/HIS) are already on TACC in data/external/glgc/.
# The MVP ancestry-stratified CAD outcomes require dbGaP/EGA phs002453 access
# (GCST90132304 = African, GCST90132303 = Hispanic). Download them, drop the
# files at the paths below, edit the column names to match, then run this once.
# ============================================================================
set -e
BASE=/path/to/cad-genetics
GLGC=$BASE/data/external/glgc
ENG=$BASE/results/agent1_genetics/strengthen/multiancestry_mr.R
OUT=$BASE/results/agent1_genetics/strengthen/multiancestry_mr
mkdir -p $OUT; rm -f $OUT/multiancestry_mr_results.csv
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env

# >>>>>>>>>>>>>>>>>>>> FILL IN AFTER dbGaP phs002453 DOWNLOAD <<<<<<<<<<<<<<<<<<<
declare -A CAD
CAD[AFR]="$BASE/data/external/mvp_cad/GCST90132304.tsv"   # MVP African American CAD
CAD[HIS]="$BASE/data/external/mvp_cad/GCST90132303.tsv"   # MVP Hispanic/Latin American CAD
# Edit these to match the downloaded CAD file's header (GWAS-Catalog harmonised defaults shown):
CR=rsid; CEA=effect_allele; COA=other_allele; CB=beta; CSE=standard_error; CP=p_value
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

glgc_file(){ local t=$1 a=$2
  if [ "$a" = HIS ]; then echo "$GLGC/${t}_INV_${a}_1KGP3_ALL.meta.singlevar.results.gz";
  else echo "$GLGC/${t}_INV_${a}_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"; fi; }

for ANC in AFR HIS; do
  CADF=${CAD[$ANC]}
  if [ ! -s "$CADF" ]; then echo "SKIP $ANC — CAD file missing: $CADF (download from dbGaP phs002453)"; continue; fi
  for T in LDL HDL logTG; do
    EXPF=$(glgc_file $T $ANC)
    echo "======== $T x $ANC ========"
    Rscript $ENG --exp_full "$EXPF" --cad_full "$CADF" --anc $ANC --trait $T --out $OUT \
      --cad_rsid $CR --cad_ea $CEA --cad_oa $COA --cad_beta $CB --cad_se $CSE --cad_p $CP
  done
done
echo; echo "=================== RESULTS ==================="
[ -f $OUT/multiancestry_mr_results.csv ] && column -t -s, $OUT/multiancestry_mr_results.csv || echo "no results (CAD files not present yet)"
echo
echo "OVERLAP-ROBUST (optional): GLGC and MVP share MVP samples, so also run CAUSE"
echo "  (needs: install.packages('cause'); genome-wide MVP CAD + GLGC sumstats; LD pruning)."
echo "  Or MRlap (cross-trait-LDSC overlap correction). See README_multiancestry_mr.md."
