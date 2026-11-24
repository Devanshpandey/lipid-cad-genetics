#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate pyt_env
cd $SCRATCH/coloc_dist && python3 sentinels.py > sent.log 2>&1
echo SENT_DONE >> sent.log
