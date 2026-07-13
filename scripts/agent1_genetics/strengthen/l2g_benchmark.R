suppressMessages(library(data.table))
# ENSG->symbol map from GTEx egenes
map <- unique(fread(cmd="zcat /path/to/cad-genetics/data/reference/gtex_v8/Liver.v8.egenes.txt.gz | cut -f1,2", header=TRUE))
map[, ensg := sub("\\..*","",gene_id)]
m2 <- unique(map[, .(ensg, gene=gene_name)])
# L2G: max score per ENSG
l2g <- fread("/path/to/cad-genetics/results/agent2_genes/opentargets_l2g.csv")
l2g[, ensg := sub("\\..*","",gene_name)]
lmax <- l2g[, .(l2g=max(l2g_score,na.rm=TRUE)), by=ensg]
lmax <- merge(lmax, m2, by="ensg")[, .(gene, l2g)]
lmax <- lmax[lmax[, .I[which.max(l2g)], by=gene]$V1]   # one row per symbol
# convergence
g <- fread("/path/to/cad-genetics/results/agent1_genetics/strengthen/gene_ranking_scores.csv")
m <- merge(g, lmax, by="gene")
auc <- function(s,y){r<-rank(s);np<-sum(y==1);nn<-sum(y==0);if(np==0||nn==0)return(NA);(sum(r[y==1])-np*(np+1)/2)/(np*nn)}
cat("=== L2G BENCHMARK (genes scored by BOTH convergence and Open Targets L2G) ===\n")
cat(sprintf("matched genes: %d ; canonical positives in set: %d\n", nrow(m), sum(m$canonical)))
cat(sprintf("AUROC integrated (convergence): %.3f\n", auc(m$integrated, m$canonical)))
cat(sprintf("AUROC L2G (Open Targets)      : %.3f\n", auc(m$l2g, m$canonical)))
cat(sprintf("AUROC burden_only             : %.3f\n", auc(m$burden_only, m$canonical)))
cat(sprintf("Spearman cor(integrated, L2G) : %.3f\n", cor(m$integrated, m$l2g, method="spearman", use="complete.obs")))
canon <- g[canonical==1]; cl <- merge(canon, lmax, by="gene", all.x=TRUE)
cat(sprintf("\n=== COMPLEMENTARITY ===\ncanonical total: %d ; with an L2G score: %d ; rare-variant-only (no L2G): %d\n",
    nrow(canon), sum(!is.na(cl$l2g)), sum(is.na(cl$l2g))))
cat("canonical invisible to L2G:", paste(sort(cl[is.na(l2g),gene]),collapse=", "), "\n")
cat(sprintf("\nPDE3B in L2G set? %s\n", ifelse("PDE3B"%in%lmax$gene, sprintf("yes, L2G=%.3f", lmax[gene=="PDE3B",l2g]), "NO (no GWAS-coloc locus -> unscored by L2G)")))
fwrite(m[,.(gene,canonical,integrated,burden_only,l2g,coloc)], "l2g_benchmark.csv")
cat("wrote l2g_benchmark.csv\n")
