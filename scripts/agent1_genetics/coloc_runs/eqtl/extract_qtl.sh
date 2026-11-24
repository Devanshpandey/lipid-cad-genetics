#!/bin/bash
set -u
REF=/path/to/cad-genetics/data/reference/eqtl_catalogue
GW=/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged
W=$SCRATCH/qtl_coloc
ENS=$(awk '{print $2}' $W/genes.tsv | paste -sd'|')
for TIS in Liver Artery_Aorta Artery_Coronary Heart_Left_Ventricle Whole_Blood; do
  [ -s "$W/${TIS}_qtl.tsv.gz" ] && { echo "$TIS cached" >> $W/ex.log; continue; }
  # col 7 = gene_id, col 18 = rsid (header-defined layout, 18 cols)
  zcat $REF/$TIS.tsv.gz | awk -F'\t' -v ids="$ENS" \
    'BEGIN{n=split(ids,a,"|");for(i=1;i<=n;i++)w[a[i]]=1} NR==1||($7 in w)' \
    | gzip > $W/${TIS}_qtl.tsv.gz
  echo "$TIS rows=$(zcat $W/${TIS}_qtl.tsv.gz|wc -l) genes=$(zcat $W/${TIS}_qtl.tsv.gz|tail -n+2|cut -f7|sort -u|wc -l)" >> $W/ex.log
done
# rsID union across tissues -> build-agnostic key for subsetting the GWAS
zcat $W/*_qtl.tsv.gz | tail -n +2 | cut -f18 | grep '^rs' | sort -u > $W/qtl_rsids.txt
echo "qtl rsids: $(wc -l < $W/qtl_rsids.txt)" >> $W/ex.log
for TR in lipids_LDL_C lipids_TRIGLY lipids_HDL_C lipids_APOB outcomes_CAD; do
  [ -s "$W/gwas_subset/${TR}.txt.gz" ] && continue
  zcat $GW/${TR}.txt.gz | awk -v f=$W/qtl_rsids.txt \
    'BEGIN{while((getline l<f)>0) k[l]=1} NR==1||($3 in k)' | gzip > $W/gwas_subset/${TR}.txt.gz
  echo "$TR subset=$(zcat $W/gwas_subset/${TR}.txt.gz|wc -l)" >> $W/ex.log
done
echo EXTRACT_DONE >> $W/ex.log
