#!/bin/bash
# Extract UKB lipid instruments from CARDIoGRAM CAD + GLGC ancestry-stratified lipids.
set -e
BASE=/path/to/cad-genetics
MERGED=$BASE/results/agent1_genetics/gwas_summary_stats/merged
EXT=$BASE/data/external
OUT=$BASE/results/agent1_genetics/strengthen/external_mr
mkdir -p $OUT
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env

echo "[1] UKB instruments (GW-sig, 10 Mb greedy clump) for LDL_C, HDL_C, TRIGLY"
Rscript -e '
suppressMessages(library(data.table))
MERGED <- Sys.getenv("MERGED"); OUT <- Sys.getenv("OUT")
clump <- function(tr){
  d <- fread(cmd=paste0("gunzip -c ",MERGED,"/lipids_",tr,".txt.gz | awk \x27NR==1||$12>7.30103\x27"))
  d <- d[!is.na(BETA)&!is.na(SE)][order(-LOG10P)]
  keep<-logical(nrow(d)); kc<-character(0); kp<-integer(0)
  for(i in seq_len(nrow(d))){ch<-as.character(d$CHROM[i]);po<-as.integer(d$GENPOS[i]);s<-kc==ch
    if(!length(kp)||!any(s&abs(kp[s]-po)<1e7)){keep[i]<-TRUE;kc<-c(kc,ch);kp<-c(kp,po)}}
  d[keep, .(rsID=ID, effA=ALLELE1, othA=ALLELE0, beta=BETA, se=SE, trait=tr)]
}
ins <- rbindlist(lapply(c("LDL_C","HDL_C","TRIGLY"), clump))
fwrite(ins, paste0(OUT,"/ukb_instruments.csv"))
writeLines(unique(ins$rsID), paste0(OUT,"/instrument_rsids.txt"))
cat(sprintf("   instruments: %d rows, %d unique SNPs (LDL %d, HDL %d, TG %d)\n",
  nrow(ins), uniqueN(ins$rsID), sum(ins$trait=="LDL_C"), sum(ins$trait=="HDL_C"), sum(ins$trait=="TRIGLY")))
'
IDS=$OUT/instrument_rsids.txt
echo "[2] extract instruments from CARDIoGRAM CAD"
awk 'NR==FNR{a[$1];next} FNR==1||($1 in a)' $IDS $EXT/cardiogram/cad.add.160614.website.txt > $OUT/cardiogram_ins.txt
echo "   CARDIoGRAM instrument rows: $(($(wc -l < $OUT/cardiogram_ins.txt)-1))"
echo "[3] extract instruments from GLGC ancestry-specific lipids"
for f in $EXT/glgc/*.gz; do
  b=$(basename $f .meta.singlevar.results.gz)
  zcat $f | awk 'NR==FNR{a[$1];next} FNR==1||($1 in a)' $IDS - > $OUT/glgc_${b}.txt
  echo "   $b: $(($(wc -l < $OUT/glgc_${b}.txt)-1)) rows"
done
echo "EXTRACTION DONE"
