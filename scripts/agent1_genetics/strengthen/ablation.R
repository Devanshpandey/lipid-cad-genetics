suppressMessages(library(data.table))
g <- fread("gene_ranking_scores.csv")
g[, `:=`(zlip=as.numeric(scale(lip_lp)), zout=as.numeric(scale(out_lp)))]
auc <- function(p,y){r<-rank(p);np<-sum(y==1);nn<-sum(y==0);(sum(r[y==1])-np*(np+1)/2)/(np*nn)}
mods <- list(
 "burden only (lipid)"        = canonical ~ zlip,
 "+ colocalization"           = canonical ~ zlip + coloc,
 "+ outcome burden"           = canonical ~ zlip + zout,
 "full (burden+coloc+outcome)"= canonical ~ zlip + zout + coloc)
cat("=== LAYER ABLATION: AUROC for recovering canonical targets (n=",sum(g$canonical),"of",nrow(g),") ===\n")
fits <- list()
for (nm in names(mods)) {
  f <- glm(mods[[nm]], family=binomial, data=g); fits[[nm]] <- f
  cat(sprintf("  %-30s AUROC=%.3f\n", nm, auc(predict(f,type="response"), g$canonical)))
}
cat(sprintf("  %-30s AUROC=%.3f\n", "integrated score (as reported)", auc(g$integrated, g$canonical)))
cat("\n=== nested LRT vs burden-only ===\n")
cat(sprintf("  +coloc            : LRT P=%.2e\n", anova(fits[[1]],fits[[2]],test='Chisq')$`Pr(>Chi)`[2]))
cat(sprintf("  +outcome          : LRT P=%.2e\n", anova(fits[[1]],fits[[3]],test='Chisq')$`Pr(>Chi)`[2]))
cat(sprintf("  full vs burden    : LRT P=%.2e\n", anova(fits[[1]],fits[[4]],test='Chisq')$`Pr(>Chi)`[2]))
