#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env
W=$SCRATCH/qtl_coloc
Rscript /path/to/cad-genetics/scripts/agent1_genetics/R/run_qtl_coloc.R \
  --genes "$W/genes.tsv" \
  --qtl_dir "$W" \
  --tissues "Liver,Artery_Aorta,Artery_Coronary,Heart_Left_Ventricle,Whole_Blood" \
  --sumstat_dir "$W/gwas_subset" \
  --traits "lipids_LDL_C,lipids_TRIGLY,lipids_HDL_C,lipids_APOB,outcomes_CAD" \
  --window_kb 500 \
  --out "$W/cis_eqtl_coloc_v2.csv" > $W/eqtl_coloc.log 2>&1
echo EQTL_COLOC_DONE >> $W/eqtl_coloc.log
