suppressPackageStartupMessages({library(data.table); library(MVMR)})
W <- Sys.getenv("W"); M <- Sys.getenv("M")
traits <- c("LDL_C","APOB","TRIGLY")
# greedy distance clumping (10 Mb) per trait on GWS SNPs
clump <- function(dt){
  dt <- dt[order(-LOG10P)]; keep_id<-character(); keep_chr<-integer(); keep_pos<-numeric()
  for(i in seq_len(nrow(dt))){
    c<-dt$CHROM[i]; p<-dt$GENPOS[i]
    if(!any(keep_chr==c & abs(keep_pos-p)<1e7)){ keep_id<-c(keep_id,dt$ID[i]); keep_chr<-c(keep_chr,c); keep_pos<-c(keep_pos,p)}
  }; keep_id
}
ivs <- character()
for(t in traits){
  g <- fread(file.path(W,paste0("gws_",t,".txt")))
  ivs <- union(ivs, clump(g))
}
cat("union instruments:", length(ivs), "\n")
idf <- file.path(W,"ivids.txt"); writeLines(ivs, idf)
# grep instrument rows from full sumstats (memory-light)
load_iv <- function(f){
  dt <- fread(cmd=paste0("{ zcat ",f," | head -1; zcat ",f," | grep -Fwf ",idf,"; }"))
  dt[ID %in% ivs]
}
E <- lapply(traits, function(t) load_iv(file.path(M,paste0("lipids_",t,".txt.gz"))))
names(E) <- traits
O <- load_iv(file.path(M,"outcomes_CAD.txt.gz"))
# common IDs across all
common <- Reduce(intersect, c(lapply(E,function(x)x$ID), list(O$ID)))
cat("instruments with complete data:", length(common), "\n")
ref <- E[[1]][match(common,ID)]  # reference alleles from trait1
BX <- se <- matrix(NA, length(common), length(traits), dimnames=list(common,traits))
for(t in traits){ d<-E[[t]][match(common,ID)]
  flip <- ifelse(d$ALLELE1==ref$ALLELE1, 1, ifelse(d$ALLELE1==ref$ALLELE0, -1, NA))
  BX[,t]<-d$BETA*flip; se[,t]<-d$SE }
od <- O[match(common,ID)]
of <- ifelse(od$ALLELE1==ref$ALLELE1,1,ifelse(od$ALLELE1==ref$ALLELE0,-1,NA))
BY <- od$BETA*of; seY <- od$SE
ok <- complete.cases(BX,se,BY,seY)
BX<-BX[ok,,drop=F]; se<-se[ok,,drop=F]; BY<-BY[ok]; seY<-seY[ok]
cat("final SNPs:", nrow(BX), "\n")
fm <- format_mvmr(BXGs=BX, BYG=BY, seBXGs=se, seBYG=seY, RSID=rownames(BX))
fs <- strength_mvmr(fm, gencov=0)
cat("=== Conditional F-statistics (LDL_C / APOB / TRIGLY) ===\n")
print(fs)
