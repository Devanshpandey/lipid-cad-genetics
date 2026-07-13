#!/bin/bash
#SBATCH -J ma_sens
#SBATCH -o /path/to/cad-genetics/logs/claude_jobs/ma_sens_%j.out
#SBATCH -e /path/to/cad-genetics/logs/claude_jobs/ma_sens_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:30:00
#SBATCH -A YOUR_ALLOCATION
set -e
GLGC=/path/to/cad-genetics/data/external/glgc
MVP=/path/to/cad-genetics/data/external/mvp_cad
OUT=/path/to/cad-genetics/results/agent1_genetics/strengthen/multiancestry_mr
ENG=/path/to/cad-genetics/results/agent1_genetics/strengthen/multiancestry_mr_sens.R
rm -f $OUT/multiancestry_mr_sensitivity.csv
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
COLS="--cad_rsid SNP_ID --cad_ea ea --cad_ref ref --cad_alt alt --cad_or or --cad_ci ci --cad_p pval"
run(){ local ANC=$1 GA=$2 CAD=$3
  for T in LDL HDL logTG; do
    if [ "$GA" = HIS ]; then EXP=$GLGC/${T}_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz
    else EXP=$GLGC/${T}_INV_${GA}_HRC_1KGP3_others_ALL.meta.singlevar.results.gz; fi
    echo "==== $T $ANC sensitivity ===="
    Rscript $ENG --exp_full "$EXP" --cad_full "$MVP/$CAD" --anc $ANC --trait $T --out $OUT $COLS 2>&1 | grep -vE "built under|^Analysing|Loading|API|Harmonising|Removing|^rs"
  done; }
run AFR AFR mvp_cad_AFR.txt.gz
run HIS HIS mvp_cad_AMR.txt.gz
echo "MA SENSITIVITY DONE"
