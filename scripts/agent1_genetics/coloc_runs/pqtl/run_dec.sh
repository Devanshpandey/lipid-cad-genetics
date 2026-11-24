#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate pyt_env
python3 /path/to/cad-genetics/scripts/agent1_genetics/decode_coloc_v2.py > $SCRATCH/pqtl2/decode_coloc.log 2>&1
echo DECODE_COLOC_DONE >> $SCRATCH/pqtl2/decode_coloc.log
