#!/bin/bash
#SBATCH -J xanc_scatter
#SBATCH -o /path/to/cad-genetics/logs/claude_jobs/xanc_scatter_%j.out
#SBATCH -e /path/to/cad-genetics/logs/claude_jobs/xanc_scatter_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:30:00
#SBATCH -A YOUR_ALLOCATION
set -e
G=/path/to/cad-genetics/data/external/glgc
OUT=/path/to/cad-genetics/results/agent1_genetics/strengthen/xancestry
mkdir -p $OUT
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
for T in LDL logTG; do
  EUR=$G/${T}_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz
  AFR=$G/${T}_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz
  HIS=$G/${T}_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz
  echo "[$T] selecting EUR gwsig, 1Mb-pruned instruments"
  zcat $EUR | awk 'NR>1 && $12!="" && ($12+0)<5e-8 {
      b=int($3/1000000); key=$2"_"b; p=$12+0;
      if(!(key in bp) || p<bp[key]){bp[key]=p; line[key]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9"\t"$10"\t"$12}
    } END{for(k in line) print line[k]}' > $OUT/${T}_eur_instruments.tsv
  cut -f1 $OUT/${T}_eur_instruments.tsv | sort -u > $OUT/${T}_snps.txt
  echo "[$T] EUR instruments: $(wc -l < $OUT/${T}_eur_instruments.tsv)"
  for A in AFR HIS; do
    F=$AFR; [ "$A" = HIS ] && F=$HIS
    zcat $F | awk -v S=$OUT/${T}_snps.txt 'BEGIN{while((getline l<S)>0) want[l]=1}
      NR>1 && ($1 in want){print $1"\t"$4"\t"$5"\t"$9"\t"$10}' > $OUT/${T}_${A}.tsv
    echo "[$T] $A extracted: $(wc -l < $OUT/${T}_${A}.tsv)"
  done
done
Rscript - <<'RS'
library(data.table)
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen/xancestry"
merge_anc <- function(trait){
  eur <- fread(file.path(OUT,sprintf("%s_eur_instruments.tsv",trait)), header=FALSE,
               col.names=c("rsID","chr","pos","REF","ALT","eur_beta","eur_se","eur_p"))
  out <- eur[, .(rsID,chr,pos,REF,ALT,eur_beta,eur_se,eur_p)]
  for(A in c("AFR","HIS")){
    a <- fread(file.path(OUT,sprintf("%s_%s.tsv",trait,A)), header=FALSE,
               col.names=c("rsID","REF2","ALT2","beta2","se2"))
    a <- a[!duplicated(rsID)]
    m <- merge(out[, .(rsID,REF,ALT)], a, by="rsID")
    m[, align := fifelse(REF==REF2 & ALT==ALT2, 1, fifelse(REF==ALT2 & ALT==REF2, -1, NA_real_))]
    m <- m[!is.na(align)]
    m[, bb := beta2*align]
    lo <- tolower(A)
    add <- m[, .(rsID, bb, se2)]; setnames(add, c("rsID", paste0(lo,"_beta"), paste0(lo,"_se")))
    out <- merge(out, add, by="rsID", all.x=TRUE)
  }
  fwrite(out, file.path(OUT, sprintf("xanc_%s.csv", trait)))
  cat(sprintf("%s: n=%d  EUR-AFR=%d  EUR-HIS=%d\n", trait, nrow(out),
      sum(!is.na(out$afr_beta)), sum(!is.na(out$his_beta))))
  for(A in c("afr","his")){
    ok <- out[!is.na(get(paste0(A,"_beta")))]
    cat(sprintf("   EUR vs %s: r=%.3f  sign-conc=%.3f  n=%d\n", toupper(A),
        cor(ok$eur_beta, ok[[paste0(A,"_beta")]]),
        mean(sign(ok$eur_beta)==sign(ok[[paste0(A,"_beta")]])), nrow(ok)))
  }
}
merge_anc("LDL"); merge_anc("logTG")
RS
echo "XANC SCATTER DONE"
