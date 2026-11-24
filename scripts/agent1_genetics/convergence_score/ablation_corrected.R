suppressMessages(library(data.table))
g <- fread("gene_ranking_scores_corrected.csv")
g[, `:=`(zlip=as.numeric(scale(lip_lp)), zout=as.numeric(scale(out_lp)))]
auc <- function(p,y){r<-rank(p);np<-sum(y==1);nn<-sum(y==0);(sum(r[y==1])-np*(np+1)/2)/(np*nn)}

run <- function(colocvar, label, scorevar) {
  g[, coloc_use := get(colocvar)]
  cat(sprintf("\n=== %s  (coloc=1 for %d of %d genes; %d of %d canonical) ===\n",
      label, sum(g$coloc_use), nrow(g), sum(g$coloc_use & g$canonical), sum(g$canonical)))
  mods <- list("burden only (lipid)"=canonical ~ zlip,
               "+ colocalization"   =canonical ~ zlip + coloc_use,
               "+ outcome burden"   =canonical ~ zlip + zout,
               "full"               =canonical ~ zlip + zout + coloc_use)
  fits <- list()
  for (nm in names(mods)) {
    f <- glm(mods[[nm]], family=binomial, data=g); fits[[nm]] <- f
    cat(sprintf("  %-28s AUROC=%.3f\n", nm, auc(predict(f,type="response"), g$canonical)))
  }
  cat(sprintf("  %-28s AUROC=%.3f\n", "integrated score", auc(g[[scorevar]], g$canonical)))
  cat(sprintf("  LRT +coloc vs burden-only : P=%.3g\n",
      anova(fits[[1]],fits[[2]],test='Chisq')$`Pr(>Chi)`[2]))
  cat(sprintf("  LRT full  vs burden-only : P=%.3g\n",
      anova(fits[[1]],fits[[4]],test='Chisq')$`Pr(>Chi)`[2]))
}

cat("################ PUBLISHED (buggy coloc annotation) ################\n")
run("coloc","published","integrated")
cat("\n################ CORRECTED ################\n")
run("coloc_w500","corrected, +/-500kb = the Methods definition","integrated_w500")
run("coloc_w100","corrected, +/-100kb cis convention","integrated_w100")
run("coloc_body","corrected, strict gene-body containment","integrated_body")
