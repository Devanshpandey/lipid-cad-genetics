suppressPackageStartupMessages(library(coloc))
W<-Sys.getenv("W")
loci<-read.table(file.path(W,"loci.txt"),col.names=c("name","chr","lo","hi","lipid"),stringsAsFactors=FALSE)
res<-data.frame()
for(i in 1:nrow(loci)){
  nm<-loci$name[i]
  u<-read.table(file.path(W,paste0(nm,"_ukb.txt")),header=TRUE,stringsAsFactors=FALSE)
  colnames(u)<-c("rsid","A0","A1","af","N","beta","se")
  f<-read.table(file.path(W,paste0(nm,"_fg.txt")),header=TRUE,comment.char="",stringsAsFactors=FALSE)
  colnames(f)<-c("chrom","pos","ref","alt","rsid","genes","pval","mlogp","beta","sebeta","af","afc","afco")
  m<-merge(u,f,by="rsid",suffixes=c(".u",".f"))
  m<-m[!is.na(m$beta.u)&!is.na(m$se)&!is.na(m$beta.f)&!is.na(m$sebeta)&m$se>0&m$sebeta>0,]
  # harmonize FinnGen beta to UKB effect allele (A1)
  m$flip<-ifelse(toupper(m$alt)==toupper(m$A1),1,ifelse(toupper(m$ref)==toupper(m$A1),-1,NA))
  m<-m[!is.na(m$flip),]
  m$bf<-m$beta.f*m$flip
  m$maf<-pmin(m$af.u,1-m$af.u); m<-m[m$maf>0&m$maf<1,]
  m<-m[!duplicated(m$rsid),]
  D1<-list(beta=m$beta.u,varbeta=m$se^2,N=round(mean(m$N)),MAF=m$maf,type="quant",snp=m$rsid)
  D2<-list(beta=m$bf,varbeta=m$sebeta^2,type="cc",s=0.13,N=394000,MAF=m$maf,snp=m$rsid)
  c1<-coloc.abf(D1,D2)  # default priors
  c2<-coloc.abf(D1,D2,p12=1e-6)  # stringent prior
  res<-rbind(res,data.frame(locus=nm,lipid=loci$lipid[i],nsnp=nrow(m),
    PPH4_default=round(c1$summary["PP.H4.abf"],3),
    PPH4_p12_1e6=round(c2$summary["PP.H4.abf"],3),
    PPH3_default=round(c1$summary["PP.H3.abf"],3)))
}
write.csv(res,file.path(W,"finngen_coloc_results.csv"),row.names=FALSE)
cat("=== External colocalization: UKB lipids x FinnGen CAD ===\n")
print(res,row.names=FALSE)
