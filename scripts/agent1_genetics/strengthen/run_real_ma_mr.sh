set -e
GLGC=/path/to/cad-genetics/data/external/glgc
MVP=/path/to/cad-genetics/data/external/mvp_cad
OUT=/path/to/cad-genetics/results/agent1_genetics/strengthen/multiancestry_mr
ENG=/path/to/cad-genetics/results/agent1_genetics/strengthen/multiancestry_mr.R
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
COLS="--cad_rsid SNP_ID --cad_ea ea --cad_ref ref --cad_alt alt --cad_or or --cad_ci ci --cad_p pval"
run(){ local ANC=$1 GA=$2 CAD=$3
  for T in LDL HDL logTG; do
    if [ "$GA" = HIS ]; then EXP=$GLGC/${T}_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz
    else EXP=$GLGC/${T}_INV_${GA}_HRC_1KGP3_others_ALL.meta.singlevar.results.gz; fi
    echo "==== $T  ${ANC} (GLGC-$GA lipid -> MVP-$ANC CAD) ===="
    Rscript $ENG --exp_full "$EXP" --cad_full "$MVP/$CAD" --anc $ANC --trait $T --out $OUT $COLS 2>&1 | grep -vE "built under|Warning|^Analysing|Loading|API|Harmonising|Removing|MRPRESSO|outlier|^rs"
  done; }
run AFR AFR mvp_cad_AFR.txt.gz
run HIS HIS mvp_cad_AMR.txt.gz
echo "REAL MULTI-ANCESTRY MR DONE"
