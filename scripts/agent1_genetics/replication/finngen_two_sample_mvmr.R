suppressPackageStartupMessages({library(data.table); library(MendelianRandomization)})
W<-Sys.getenv("W"); M<-Sys.getenv("M"); FG<-Sys.getenv("FG")
traits<-c("LDL_C","HDL_C","TRIGLY","LPA")
clump<-function(dt){dt<-dt[order(-LOG10P)];kid<-character();kc<-integer();kp<-numeric()
 for(i in seq_len(nrow(dt))){c<-dt$CHROM[i];p<-dt$GENPOS[i]
  if(!any(kc==c&abs(kp-p)<1e7)){kid<-c(kid,dt$ID[i]);kc<-c(kc,c);kp<-c(kp,p)}};kid}
ivs<-character(); for(t in traits){g<-fread(file.path(W,paste0("gws_",t,".txt")));ivs<-union(ivs,clump(g))}
cat("instruments:",length(ivs),"\n")
idf<-file.path(W,"ivids_fg.txt"); writeLines(ivs,idf)
gi<-function(f,cmd_id="ID"){fread(cmd=paste0("{ zcat ",f," | head -1; zcat ",f," | grep -Fwf ",idf,"; }"))}
E<-lapply(traits,function(t){d<-gi(file.path(M,paste0("lipids_",t,".txt.gz")));d[ID%in%ivs]}); names(E)<-traits
# FinnGen: match by rsids (col5), effect allele = alt (col4), ref=col3
fg<-gi(FG); setnames(fg,c("chrom","pos","ref","alt","rsids","genes","pval","mlogp","beta","sebeta","af","afc","afco"))
fg<-fg[rsids%in%ivs]
common<-Reduce(intersect,c(lapply(E,function(x)x$ID),list(fg$rsids)))
cat("common instruments (UKB∩FinnGen):",length(common),"\n")
ref<-E[[1]][match(common,ID)]
BX<-se<-matrix(NA,length(common),length(traits),dimnames=list(common,traits))
for(t in traits){d<-E[[t]][match(common,ID)]
 fl<-ifelse(d$ALLELE1==ref$ALLELE1,1,ifelse(d$ALLELE1==ref$ALLELE0,-1,NA));BX[,t]<-d$BETA*fl;se[,t]<-d$SE}
f<-fg[match(common,rsids)]
of<-ifelse(f$alt==ref$ALLELE1,1,ifelse(f$ref==ref$ALLELE1,-1,NA))
BY<-f$beta*of; seY<-f$sebeta
ok<-complete.cases(BX,se,BY,seY); BX<-BX[ok,,drop=F];se<-se[ok,,drop=F];BY<-BY[ok];seY<-seY[ok]
cat("final SNPs for two-sample MVMR:",nrow(BX),"\n")
inp<-mr_mvinput(bx=BX,bxse=se,by=BY,byse=seY,exposure=traits,outcome="FinnGen_CAD")
res<-mr_mvivw(inp)
eg<-tryCatch(mr_mvegger(inp),error=function(e)NULL)
cat("=== Two-sample MVMR: UKB lipids -> FinnGen CAD (I9_CHD) ===\n")
for(i in seq_along(traits)) cat(sprintf("%-7s OR=%.3f (%.3f-%.3f) p=%.2e\n",traits[i],
 exp(res@Estimate[i]),exp(res@Estimate[i]-1.96*res@StdError[i]),exp(res@Estimate[i]+1.96*res@StdError[i]),res@Pvalue[i]))
if(!is.null(eg)) cat(sprintf("MR-Egger intercept p=%.3f\n",eg@Pvalue.Int))
