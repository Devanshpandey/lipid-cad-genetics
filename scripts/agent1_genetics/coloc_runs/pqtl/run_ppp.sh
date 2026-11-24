#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate pyt_env
cd $SCRATCH/pqtl2/ukbppp
export PPPDIR=$SCRATCH/pqtl2/ukbppp OUTDIR=$SCRATCH/pqtl2/ukbppp
python3 /path/to/cad-genetics/scripts/agent1_genetics/ukbppp_coloc_v2.py > $SCRATCH/pqtl2/ppp_coloc.log 2>&1
echo PPP_COLOC_DONE >> $SCRATCH/pqtl2/ppp_coloc.log
