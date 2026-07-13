#!/bin/bash
#SBATCH -J eqtl_coloc
#SBATCH -o /path/to/cad-genetics/logs/claude_jobs/eqtl_coloc_%j.out
#SBATCH -e /path/to/cad-genetics/logs/claude_jobs/eqtl_coloc_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:30:00
#SBATCH -A YOUR_ALLOCATION
set -e
REF=/path/to/cad-genetics/data/reference/eqtl_catalogue
GWAS=/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/lipids
OUT=/path/to/cad-genetics/results/agent1_genetics/strengthen/eqtl_coloc
mkdir -p $OUT
cat > $OUT/genes.txt <<G
PCSK9 ENSG00000169174 1
SORT1 ENSG00000134243 1
CELSR2 ENSG00000143126 1
PSRC1 ENSG00000134222 1
LDLR ENSG00000130164 19
APOE ENSG00000130203 19
APOC1 ENSG00000130208 19
TRIB1 ENSG00000173334 8
APOB ENSG00000084674 2
ANGPTL3 ENSG00000132855 1
G
ENS=$(awk '{print $2}' $OUT/genes.txt | paste -sd'|')
for TIS in Liver Artery_Aorta Artery_Coronary; do
  echo "extracting $TIS eQTL for effector genes..."
  zcat $REF/$TIS.tsv.gz | awk -F'\t' -v ids="$ENS" 'BEGIN{n=split(ids,a,"|");for(i=1;i<=n;i++)want[a[i]]=1} NR==1||($7 in want)' | gzip > $OUT/${TIS}_eqtl.tsv.gz
  echo "  rows: $(zcat $OUT/${TIS}_eqtl.tsv.gz | wc -l)"
done
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
Rscript - <<'RS'
suppressMessages({library(data.table); library(coloc)})
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen/eqtl_coloc"
GW  <- "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/lipids"
genes <- fread(file.path(OUT,"genes.txt"), header=FALSE, col.names=c("sym","ens","chr"))
tissues <- c("Liver","Artery_Aorta","Artery_Coronary")
eqtl <- lapply(tissues, function(t) fread(cmd=sprintf("zcat %s/%s_eqtl.tsv.gz",OUT,t)))
names(eqtl) <- tissues
gw <- list()
for(c in unique(genes$chr)) gw[[as.character(c)]] <- fread(cmd=sprintf("zcat %s/lipids_chr%s_LDL_C.regenie.gz",GW,c))
res <- list()
for(i in 1:nrow(genes)){
  g <- genes[i]; G <- gw[[as.character(g$chr)]]
  for(t in tissues){
    e <- eqtl[[t]][gene_id==g$ens & maf>0 & !is.na(beta) & !is.na(se) & se>0]
    if(nrow(e)<30) next
    e <- e[!duplicated(rsid)]
    m <- merge(e[,.(rsid, eb=beta, ese=se, emaf=pmin(maf,1-maf), ean=an)],
               G[,.(rsid=ID, gb=BETA, gse=SE, gmaf=A1FREQ, gn=N)], by="rsid")
    m <- m[ese>0 & gse>0 & emaf>0 & pmin(gmaf,1-gmaf)>0]
    m[, gmaf2:=pmin(gmaf,1-gmaf)]
    if(nrow(m)<30) next
    d1 <- list(beta=m$gb, varbeta=m$gse^2, snp=m$rsid, MAF=m$gmaf2, N=max(m$gn), type="quant")
    d2 <- list(beta=m$eb, varbeta=m$ese^2, snp=m$rsid, MAF=m$emaf, N=round(max(m$ean)/2), type="quant")
    cr <- tryCatch(suppressWarnings(coloc.abf(d1,d2)), error=function(er) NULL)
    if(is.null(cr)) next
    s <- cr$summary
    res[[paste(g$sym,t)]] <- data.table(gene=g$sym, tissue=t, n_snp=nrow(m),
      PP.H3=round(s["PP.H3.abf"],3), PP.H4=round(s["PP.H4.abf"],3))
  }
}
out <- rbindlist(res); setorder(out, -PP.H4)
fwrite(out, file.path(OUT,"cis_eqtl_coloc.csv")); print(out)
cat("\nColoc PP.H4>0.8:", nrow(out[PP.H4>0.8]), "gene-tissue pairs\n")
RS
echo "EQTL COLOC DONE"
