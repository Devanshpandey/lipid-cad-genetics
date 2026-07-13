suppressPackageStartupMessages({library(data.table);library(susieR);library(coloc);library(Matrix)})
set.seed(42)
MERGED <- "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
LDBF   <- "/corral/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE_QC_400k/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"
PLINK1 <- "/path/to/software/plink"
OUT    <- "/path/to/cad-genetics/results/agent1_genetics/strengthen"
S_CAD  <- 0.123; W <- 5e5
loci <- data.table(name=c("PCSK9","SORT1_CELSR2","TRIB1","LDLR","APOE"),
                   exp=c("LDL_C","LDL_C","TRIGLY","LDL_C","LDL_C"),
                   chr=c(1,1,8,19,19), bp=c(55505647,109814880,126475770,11185919,45231821))
rd <- function(f,chr,bp){dt<-fread(cmd=paste("gunzip -c",shQuote(f)),na.strings="NA")
  dt[CHROM==chr & GENPOS>=bp-W & GENPOS<=bp+W & !is.na(BETA) & !is.na(SE)][order(GENPOS)]}
ld_of <- function(ids,chr,td){writeLines(ids,file.path(td,"s.txt"))
  file.remove(list.files(td,"^ld\\.",full.names=TRUE))
  system(sprintf("%s --bfile %s --chr %s --extract %s --r square --out %s --threads 4 >/dev/null 2>&1",
    PLINK1,LDBF,chr,file.path(td,"s.txt"),file.path(td,"ld")))
  f<-file.path(td,"ld.ld"); if(!file.exists(f))return(NULL); as.matrix(fread(f,header=FALSE))}
res<-list()
for(i in seq_len(nrow(loci))){
  L<-loci[i]; message("== ",L$name," ==")
  e<-rd(file.path(MERGED,sprintf("lipids_%s.txt.gz",L$exp)),L$chr,L$bp)
  o<-rd(file.path(MERGED,"outcomes_CAD.txt.gz"),L$chr,L$bp)
  com<-intersect(e$ID,o$ID); e<-e[ID %in% com][order(GENPOS)]; o<-o[ID %in% com][order(GENPOS)]
  stopifnot(all(e$ID==o$ID))
  td<-tempfile(); dir.create(td); ld<-ld_of(e$ID,L$chr,td)
  if(is.null(ld)){message("  no LD");next}
  n<-min(nrow(e),nrow(ld)); e<-e[1:n];o<-o[1:n];ld<-ld[1:n,1:n]
  bad<-apply(is.na(ld),1,any); if(any(bad)){e<-e[!bad];o<-o[!bad];ld<-ld[!bad,!bad]}
  ld<-(ld+t(ld))/2; diag(ld)<-1
  rownames(ld)<-colnames(ld)<-e$ID
  D1<-list(beta=e$BETA,varbeta=e$SE^2,snp=e$ID,position=e$GENPOS,type="quant",N=as.integer(median(e$N)),sdY=1,LD=ld)
  D2<-list(beta=o$BETA,varbeta=o$SE^2,snp=o$ID,position=o$GENPOS,type="cc",N=as.integer(median(o$N)),s=S_CAD,LD=ld)
  s1<-tryCatch(runsusie(D1,repeat_until_convergence=FALSE),error=function(e)NULL)
  s2<-tryCatch(runsusie(D2,repeat_until_convergence=FALSE),error=function(e)NULL)
  n1<-if(!is.null(s1$sets$cs))length(s1$sets$cs) else 0
  n2<-if(!is.null(s2$sets$cs))length(s2$sets$cs) else 0
  maxpp<-NA; npair<-0
  if(!is.null(s1)&&!is.null(s2)&&n1>0&&n2>0){
    cs<-tryCatch(coloc.susie(s1,s2),error=function(e)NULL)
    if(!is.null(cs)&&!is.null(cs$summary)&&nrow(cs$summary)>0){
      maxpp<-max(cs$summary$PP.H4.abf); npair<-sum(cs$summary$PP.H4.abf>0.8)}
  }
  res[[L$name]]<-data.table(locus=L$name,exposure=L$exp,n_sig_lipid=n1,n_sig_CAD=n2,
    n_coloc_pairs_PP4gt0.8=npair,max_PP.H4_susie=round(maxpp,3))
  print(res[[L$name]])
}
R<-rbindlist(res); fwrite(R,file.path(OUT,"coloc_susie_results.csv"))
cat("\n=== coloc.susie multi-signal results ===\n"); print(R)
