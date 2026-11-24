#!/usr/bin/env Rscript
# ============================================================
# Does statin therapy modify the Lp(a) -> CAD effect?
#
# DESIGN. The obvious analysis -- restrict to statin users and re-estimate -- is
# biased. Statin prescription is caused by LDL-C and by cardiovascular risk, so it
# is a collider on SNP -> LDL-C -> statin <- risk factors -> CAD. Conditioning on it
# induces association between the instrument and confounders and breaks MR
# exchangeability (index-event bias). Instead we fit an INTERACTION on the full
# cohort, which never conditions on statin use:
#
#     CAD ~ score * statin + age + sex + PC1..PC10
#
# Lp(a) is the exposure least exposed to the residual concern, because statins do
# not lower Lp(a) (they raise it slightly), so Lp(a) is close to unrelated to the
# prescribing decision.
#
# POSITIVE CONTROLS on the biomarker itself, to show the test can detect an
# interaction when one truly exists:
#   LDL score x statin -> LDL-C level : interaction EXPECTED (statins compress the
#                                       genotype-phenotype gradient)
#   LPA score x statin -> Lp(a) level : interaction NOT expected (statins do not
#                                       lower Lp(a))
# If the LDL control is null the machinery is not sensitive and the primary null
# means nothing.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
set.seed(42)

opt <- parse_args(OptionParser(option_list=list(
  make_option("--scores",  type="character", help="PLINK --export A .raw files, LPA then LDL"),
  make_option("--weights", type="character", help="matching weight files (SNP, effect allele, beta)"),
  make_option("--pheno_q", type="character"),
  make_option("--pheno_b", type="character"),
  make_option("--covar",   type="character"),
  make_option("--statin",  type="character", help="TSV: IID, statin (0/1), union of 6153 and 6177"),
  make_option("--out",     type="character")
)))

# Scores are built HERE from exported dosages rather than with PLINK2 --score.
# --score did not orient to our effect allele: a 45-variant Lp(a) score came out at
# r = -0.004 with measured Lp(a) while rs10455872 ALONE gives r = -0.64. PLINK
# --export A names each column <SNP>_<counted allele>, so the counted allele is
# recoverable from the header and we flip explicitly where it is not the effect
# allele. Orientation is then ours rather than inferred.
build_score <- function(raw_file, weight_file, label) {
  g <- fread(raw_file)
  w <- fread(weight_file, header = FALSE, col.names = c("SNP","EA","BETA"))
  dose_cols <- setdiff(names(g), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
  map <- data.table(col = dose_cols,
                    SNP = sub("_[^_]+$", "", dose_cols),
                    counted = sub("^.*_", "", dose_cols))
  map <- merge(map, w, by = "SNP")
  if (!nrow(map)) stop("no overlap between dosages and weights for ", label)
  sc <- numeric(nrow(g)); flipped <- 0L
  for (i in seq_len(nrow(map))) {
    d <- g[[map$col[i]]]
    if (toupper(map$counted[i]) != toupper(map$EA[i])) { d <- 2 - d; flipped <- flipped + 1L }
    d[is.na(d)] <- mean(d, na.rm = TRUE)
    sc <- sc + d * map$BETA[i]
  }
  message(sprintf("  [%s] %d variants scored (%d flipped to the effect allele)",
                  label, nrow(map), flipped))
  data.table(IID = as.character(g$IID), s = sc)
}
raws <- strsplit(opt$scores, ",")[[1]]
wts  <- strsplit(opt$weights, ",")[[1]]
sc <- NULL
for (k in seq_along(raws)) {
  nm <- if (grepl("LPA", raws[k])) "score_LPA" else "score_LDL"
  d  <- build_score(raws[k], wts[k], nm); setnames(d, "s", nm)
  sc <- if (is.null(sc)) d else merge(sc, d, by="IID")
}
q  <- fread(opt$pheno_q);  b <- fread(opt$pheno_b)
cv <- fread(opt$covar);    st <- fread(opt$statin)
for (D in list(sc,q,b,cv,st)) D[, IID := as.character(IID)]

d <- Reduce(function(x,y) merge(x,y,by="IID"), list(sc, q[,.(IID,LDL_C,LPA)], b[,.(IID,CAD)], cv, st))
d <- d[!is.na(statin)]
pcs <- grep("^PC", names(d), value=TRUE)[1:10]
d[, score_LPA := scale(score_LPA)[,1]]
d[, score_LDL := scale(score_LDL)[,1]]
message(sprintf("[int] N=%d  CAD cases=%d  statin users=%d (%.1f%%)",
        nrow(d), sum(d$CAD==1,na.rm=TRUE), sum(d$statin==1), 100*mean(d$statin==1)))
# A score that does not predict its own biomarker cannot support any statement about
# effect modification, so this is a hard gate rather than a warning.
for (chk in list(c("score_LPA","LPA"), c("score_LDL","LDL_C"))) {
  r <- cor(d[[chk[1]]], d[[chk[2]]], use="complete.obs")
  message(sprintf("[int] %s vs measured %s: r = %+.3f", chk[1], chk[2], r))
  if (is.na(r) || r < 0.10)
    stop(sprintf("%s does not predict %s (r=%.3f) -- refusing to model", chk[1], chk[2], r))
}

fit <- function(y, score, fam, dat) {
  f <- as.formula(sprintf("%s ~ %s * statin + age + sex + %s", y, score, paste(pcs, collapse=" + ")))
  m <- glm(f, data=dat, family=fam)
  cf <- summary(m)$coefficients
  r  <- cf[grep(":statin", rownames(cf)), ]
  s  <- cf[score, ]
  list(main=s, inter=r, n=nobs(m))
}
row <- function(lbl, r, expo=TRUE) {
  est <- if (expo) exp(r["Estimate"]) else r["Estimate"]
  lo  <- if (expo) exp(r["Estimate"]-1.96*r["Std. Error"]) else r["Estimate"]-1.96*r["Std. Error"]
  hi  <- if (expo) exp(r["Estimate"]+1.96*r["Std. Error"]) else r["Estimate"]+1.96*r["Std. Error"]
  data.table(term=lbl, estimate=est, lo95=lo, hi95=hi, P=r[4])
}

out <- list()
# --- positive / negative controls on the biomarker ---
c1 <- fit("LDL_C", "score_LDL", gaussian(), d)
out[["ctrl_LDL_biomarker_main"]]  <- row("LDL score -> LDL-C (beta/SD)", c1$main, FALSE)
out[["ctrl_LDL_biomarker_inter"]] <- row("LDL score x statin -> LDL-C  [control: expect interaction]", c1$inter, FALSE)
c2 <- fit("LPA", "score_LPA", gaussian(), d)
out[["ctrl_LPA_biomarker_main"]]  <- row("Lp(a) score -> Lp(a) (beta/SD)", c2$main, FALSE)
out[["ctrl_LPA_biomarker_inter"]] <- row("Lp(a) score x statin -> Lp(a)  [control: expect none]", c2$inter, FALSE)
# --- primary: does statin therapy modify the Lp(a) -> CAD effect? ---
p1 <- fit("CAD", "score_LPA", binomial(), d)
out[["primary_LPA_CAD_main"]]  <- row("Lp(a) score -> CAD (OR/SD)", p1$main)
out[["primary_LPA_CAD_inter"]] <- row("Lp(a) score x statin -> CAD (ratio of ORs)", p1$inter)
p2 <- fit("CAD", "score_LDL", binomial(), d)
out[["LDL_CAD_main"]]  <- row("LDL score -> CAD (OR/SD)", p2$main)
out[["LDL_CAD_inter"]] <- row("LDL score x statin -> CAD (ratio of ORs)", p2$inter)

res <- rbindlist(out, idcol="analysis")
res[, N := nrow(d)]
dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
fwrite(res, opt$out)
print(res)

# Stratified estimates for interpretation only -- these condition on the collider
# and are reported alongside the interaction, not in place of it.
strat <- rbindlist(lapply(c(0,1), function(k) {
  m <- glm(as.formula(sprintf("CAD ~ score_LPA + age + sex + %s", paste(pcs, collapse=" + "))),
           data=d[statin==k], family=binomial())
  cf <- summary(m)$coefficients["score_LPA", ]
  data.table(stratum=ifelse(k==1,"statin users","non-users"), n=nobs(m),
             OR=exp(cf[1]), lo95=exp(cf[1]-1.96*cf[2]), hi95=exp(cf[1]+1.96*cf[2]), P=cf[4])
}))
print(strat)
fwrite(strat, sub("\\.csv$", "_stratified.csv", opt$out))
message("[int] wrote ", opt$out)
