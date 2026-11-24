suppressMessages(library(data.table))
set.seed(42)
g <- fread("gene_ranking_scores_corrected.csv")
ap <- function(s,y){ o<-order(-s); y<-y[o]; tp<-cumsum(y); prec<-tp/seq_along(y)
                     if(sum(y)==0) return(NA_real_); sum(prec*y)/sum(y) }
cat("=== average precision (canonical recovery) ===\n")
cat(sprintf("  burden only                         AP=%.3f\n", ap(g$burden_only, g$canonical)))
for (v in c("integrated","integrated_w500","integrated_w100","integrated_body")) {
  cat(sprintf("  %-35s AP=%.3f\n", v, ap(g[[v]], g$canonical)))
}
cat("\n=== paired bootstrap, integrated vs burden-only (2000 reps) ===\n")
for (v in c("integrated","integrated_w500","integrated_w100")) {
  d <- replicate(2000, { i<-sample(nrow(g), replace=TRUE)
                         ap(g[[v]][i], g$canonical[i]) - ap(g$burden_only[i], g$canonical[i]) })
  ci <- quantile(d, c(.025,.975), na.rm=TRUE)
  cat(sprintf("  %-18s dAP=%+.3f  95%% CI [%+.3f, %+.3f]  P=%.3f\n",
      v, mean(d,na.rm=TRUE), ci[1], ci[2], mean(d<=0,na.rm=TRUE)*2))
}

cat("\n=== the claim the Results paragraph rests on ===\n")
# 'seven of the 27 canonical genes carry no exome-wide burden signal'
g[, burden_neg := lip_lp < 5.6]
bn <- g[canonical==1 & burden_neg==TRUE]
cat(sprintf("  canonical genes with no exome-wide lipid burden signal (-log10P<5.6): %d\n", nrow(bn)))
cat(sprintf("    %s\n", paste(sort(bn$gene), collapse=", ")))
for (v in c("coloc","coloc_w500","coloc_w100")) {
  rec <- bn[get(v)==1]$gene; mis <- bn[get(v)==0]$gene
  cat(sprintf("  %-11s recovers %d/%d: %s\n", v, length(rec), nrow(bn), paste(sort(rec),collapse=", ")))
  cat(sprintf("  %-11s misses      : %s\n", "", paste(sort(mis),collapse=", ")))
}
