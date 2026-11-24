suppressMessages(library(data.table))
g <- fread("gene_ranking_scores.csv")
# reverse-engineer the score: integrated = 6*min(lip_lp,60)/60 + 2*coloc + 2*I_out
g[, bt := 6*pmin(lip_lp,60)/60]
g[, I_out := as.integer(round((integrated - bt - 2*coloc)/2))]
g[, recon := bt + 2*coloc + 2*I_out]
cat(sprintf("formula check: max |recon - integrated| = %.4g (0 confirms 6:2:2, cap 60)\n", max(abs(g$recon-g$integrated))))
auc <- function(s,y){r<-rank(s);np<-sum(y==1);nn<-sum(y==0);if(np==0||nn==0)return(NA);(sum(r[y==1])-np*(np+1)/2)/(np*nn)}
score <- function(wB,wC,wO,cap) wB*pmin(g$lip_lp,cap)/cap + wC*g$coloc + wO*g$I_out
# ---- WEIGHT SENSITIVITY ----
grid <- CJ(wB=c(4,6,8), wC=c(1,2,4), wO=c(0,2,3), cap=c(20,40,60,Inf))
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i){
  s <- score(grid$wB[i],grid$wC[i],grid$wO[i],grid$cap[i])
  rk <- frank(-s, ties.method="min")
  data.table(grid[i], auc=auc(s,g$canonical), pde3b_rank=rk[g$gene=="PDE3B"],
             pde3b_pct=100*rk[g$gene=="PDE3B"]/nrow(g),
             top10_canon=sum(g$canonical[order(-s)][1:10]))
}))
cat(sprintf("\n=== WEIGHT SENSITIVITY (%d specifications) ===\n", nrow(res)))
cat(sprintf("canonical-recovery AUROC: median %.3f  range %.3f-%.3f\n", median(res$auc), min(res$auc), max(res$auc)))
cat(sprintf("PDE3B rank: median %d  range %d-%d  (always top %.1f%%)\n", as.integer(median(res$pde3b_rank)),
    min(res$pde3b_rank), max(res$pde3b_rank), max(res$pde3b_pct)))
cat(sprintf("PDE3B in top decile (<=1800) in %d/%d specs; in top 1%% (<=180) in %d/%d\n",
    sum(res$pde3b_rank<=1800), nrow(res), sum(res$pde3b_rank<=180), nrow(res)))
cat(sprintf("top-10 genes canonical count: median %d range %d-%d\n", as.integer(median(res$top10_canon)), min(res$top10_canon), max(res$top10_canon)))
fwrite(res, "score_weight_sensitivity.csv")
# ---- STRATIFIED AUROC ----
s14 <- fread("TableS14.csv")
drug <- s14[grepl("drug target|Approved", target_category, ignore.case=TRUE), gene]
mend <- s14[grepl("Mendelian", target_category, ignore.case=TRUE), gene]
cat(sprintf("\n=== STRATIFIED AUROC (integrated score; negatives = non-canonical background) ===\n"))
cat(sprintf("approved/trial drug targets (n=%d): AUROC=%.3f\n", length(drug), auc(g$integrated, as.integer(g$gene %in% drug))))
cat(sprintf("Mendelian dyslipidemia genes (n=%d): AUROC=%.3f\n", length(mend), auc(g$integrated, as.integer(g$gene %in% mend))))
cat(sprintf("all canonical targets (n=%d): AUROC=%.3f\n", sum(g$canonical), auc(g$integrated, g$canonical)))
fwrite(data.table(stratum=c("approved/trial drug target","Mendelian dyslipidemia","all canonical"),
  n=c(length(drug),length(mend),sum(g$canonical)),
  AUROC=c(auc(g$integrated,as.integer(g$gene%in%drug)), auc(g$integrated,as.integer(g$gene%in%mend)), auc(g$integrated,g$canonical))),
  "stratified_auroc.csv")
