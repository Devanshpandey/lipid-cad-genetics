suppressPackageStartupMessages({library(data.table);library(susieR);library(coloc)})
set.seed(42)
M<-"/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
LDBF<-"/corral/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE_QC_400k/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"
PLINK1<-"/path/to/software/plink"; S_CAD<-0.123; W<-1.2e5
loci<-data.table(name=c("PCSK9","LDLR"),chr=c(1,19),bp=c(55505647,11197598))
rd<-function(f,chr,bp){dt<-fread(cmd=paste("gunzip -c",shQuote(f)),na.strings="NA")
  dt[CHROM==chr & GENPOS>=bp-W & GENPOS<=bp+W & !is.na(BETA)&!is.na(SE)][order(GENPOS)]}
res<-list()
for(i in 1:nrow(loci)){L<-loci[i]; cat("==",L$name,"==\n")
  e<-rd(file.path(M,"lipids_LDL_C.txt.gz"),L$chr,L$bp); o<-rd(file.path(M,"outcomes_CAD.txt.gz"),L$chr,L$bp)
  com<-intersect(e$ID,o$ID); e<-e[ID%in%com][order(GENPOS)]; o<-o[ID%in%com][order(GENPOS)]
  td<-tempfile();dir.create(td);writeLines(e$ID,file.path(td,"s.txt"))
  system(sprintf("%s --bfile %s --chr %s --extract %s --r square --out %s --threads 8 >/dev/null 2>&1",PLINK1,LDBF,L$chr,file.path(td,"s.txt"),file.path(td,"ld")))
  ld<-as.matrix(fread(file.path(td,"ld.ld"),header=FALSE))
  n<-min(nrow(e),nrow(ld));e<-e[1:n];o<-o[1:n];ld<-ld[1:n,1:n]
  bad<-apply(is.na(ld),1,any);if(any(bad)){e<-e[!bad];o<-o[!bad];ld<-ld[!bad,!bad]}
  ld<-(ld+t(ld))/2;diag(ld)<-1; rownames(ld)<-colnames(ld)<-e$ID
  D1<-list(beta=e$BETA,varbeta=e$SE^2,snp=e$ID,position=e$GENPOS,type="quant",N=as.integer(median(e$N)),sdY=1,LD=ld)
  D2<-list(beta=o$BETA,varbeta=o$SE^2,snp=o$ID,position=o$GENPOS,type="cc",N=as.integer(median(o$N)),s=S_CAD,LD=ld)
  s1<-runsusie(D1); s2<-runsusie(D2)
  n1<-length(s1$sets$cs); n2<-length(s2$sets$cs); mx<-NA; npr<-0
  if(n1>0&&n2>0){cs<-coloc.susie(s1,s2)
    if(!is.null(cs$summary)&&nrow(cs$summary)>0){mx<-max(cs$summary$PP.H4.abf);npr<-sum(cs$summary$PP.H4.abf>0.8)}}
  res[[L$name]]<-data.table(locus=L$name,nsnp=nrow(e),n_signals_LDL=n1,n_signals_CAD=n2,n_coloc_pairs_PP4gt0.8=npr,max_PP.H4=round(mx,3))
  print(res[[L$name]])}
R<-rbindlist(res); fwrite(R,file.path("/path/to/cad-genetics/results/agent1_genetics/strengthen","coloc_susie_lean.csv"))
cat("\n=== DONE ===\n"); print(R)
