#!/usr/bin/env Rscript
# ============================================================================
# Multi-ancestry causal MR: GLGC-{ANC} lipid  ->  MVP-{ANC} CAD
# Methods: IVW, MR-Egger (+intercept), weighted median, MR-RAPS, MR-PRESSO.
# Overlap-robust CAUSE is gated (needs genome-wide sumstats + install) — see README.
#
# Self-contained: streams the huge GLGC exposure gz and the MVP CAD outcome via
# awk prefilters, so only small instrument subsets ever enter memory.
#
# Usage (called by run_multiancestry_mr.sh):
#   Rscript multiancestry_mr.R --exp_full <GLGC gz> --cad_full <MVP CAD tsv> \
#       --anc AFR --trait LDL --out <dir> \
#       --cad_rsid rsid --cad_ea effect_allele --cad_oa other_allele \
#       --cad_beta beta --cad_se standard_error --cad_p p_value
# ============================================================================
suppressMessages({library(data.table); library(optparse)})
set.seed(42)
op <- OptionParser()
op <- add_option(op, "--exp_full", type="character")   # GLGC ancestry lipid .gz (exposure)
op <- add_option(op, "--cad_full", type="character")   # MVP ancestry CAD sumstats (outcome)
op <- add_option(op, "--anc", type="character")
op <- add_option(op, "--trait", type="character")      # LDL / HDL / logTG
op <- add_option(op, "--out", type="character")
op <- add_option(op, "--p_iv", type="double", default=5e-8)
op <- add_option(op, "--cad_rsid", type="character", default="rsid")
op <- add_option(op, "--cad_ea",   type="character", default="effect_allele")
op <- add_option(op, "--cad_oa",   type="character", default="other_allele")
op <- add_option(op, "--cad_beta", type="character", default="beta")
op <- add_option(op, "--cad_se",   type="character", default="standard_error")
op <- add_option(op, "--cad_p",    type="character", default="p_value")
op <- add_option(op, "--cad_eaf",  type="character", default="")   # optional CAD allele-freq col
# OR/CI mode (MVP SAIGE: gives OR + "lo,hi" CI instead of beta/se; ref/alt + effect-allele col)
op <- add_option(op, "--cad_or",   type="character", default="")   # if set -> OR/CI mode
op <- add_option(op, "--cad_ci",   type="character", default="ci")
op <- add_option(op, "--cad_ref",  type="character", default="ref")
op <- add_option(op, "--cad_alt",  type="character", default="alt")
o <- parse_args(op)
tmp <- tempfile()

## 1) exposure instruments: GW-sig in this ancestry, greedy 10 Mb clump
##    GLGC cols: rsID CHROM POS_b37 REF ALT N N_studies POOLED_ALT_AF EFFECT_SIZE SE pvalue_neg_log10 pvalue ...
gws <- fread(cmd=sprintf("zcat %s | awk 'NR==1 || $12<%g'", o$exp_full, o$p_iv))
gws <- gws[!is.na(EFFECT_SIZE)&!is.na(SE)][order(pvalue)]
if (nrow(gws) < 5) { gws <- fread(cmd=sprintf("zcat %s | awk 'NR==1 || $12<1e-6'", o$exp_full))[!is.na(EFFECT_SIZE)&!is.na(SE)][order(pvalue)]
                     cat(sprintf("[%s %s] <5 instruments at 5e-8; relaxed to 1e-6\n", o$trait, o$anc)) }
keep<-logical(nrow(gws)); kc<-character(0); kp<-integer(0)
for(i in seq_len(nrow(gws))){ch<-as.character(gws$CHROM[i]);po<-as.integer(gws$POS_b37[i]);s<-kc==ch
  if(!length(kp)||!any(s&abs(kp[s]-po)<1e7)){keep[i]<-TRUE;kc<-c(kc,ch);kp<-c(kp,po)}}
ivs <- gws[keep, .(SNP=rsID, ea=toupper(ALT), oa=toupper(REF), beta=EFFECT_SIZE, se=SE, eaf=POOLED_ALT_AF, p=pvalue)]
writeLines(ivs$SNP, tmp)
cat(sprintf("[%s %s] %d independent instruments\n", o$trait, o$anc, nrow(ivs)))

## 2) outcome effects at instruments (stream-extract by rsID)
catcmd <- if (grepl("\\.(gz|bgz)$", o$cad_full)) paste("zcat", shQuote(o$cad_full)) else paste("cat", shQuote(o$cad_full))
cad <- fread(cmd=sprintf("%s | awk 'NR==FNR{a[$1];next} FNR==1||($1 in a)' %s -", catcmd, tmp))
cn <- names(cad); pick <- function(x) cn[tolower(cn)==tolower(x)][1]
if (nzchar(o$cad_or)) {                       # ---- OR/CI mode (MVP SAIGE) ----
  d <- data.table(SNP=cad[[pick(o$cad_rsid)]], ea=toupper(cad[[pick(o$cad_ea)]]),
                  ref=toupper(cad[[pick(o$cad_ref)]]), alt=toupper(cad[[pick(o$cad_alt)]]),
                  orr=as.numeric(cad[[pick(o$cad_or)]]), ci=as.character(cad[[pick(o$cad_ci)]]))
  cip <- tstrsplit(d$ci, ",", fixed=TRUE)
  d[, obeta := log(orr)]
  d[, ose := (log(as.numeric(cip[[2]])) - log(as.numeric(cip[[1]]))) / (2*1.96)]
  cad <- d[, .(SNP, oea=ea, ooa=fifelse(ea==alt, ref, alt), obeta, ose, oeaf=NA_real_)]
} else {                                      # ---- beta/se mode ----
  cad[, oeaf := if (nzchar(o$cad_eaf) && !is.na(pick(o$cad_eaf))) as.numeric(get(pick(o$cad_eaf))) else NA_real_]
  cad <- cad[, .(SNP=get(pick(o$cad_rsid)), oea=toupper(get(pick(o$cad_ea))), ooa=toupper(get(pick(o$cad_oa))),
                 obeta=as.numeric(get(pick(o$cad_beta))), ose=as.numeric(get(pick(o$cad_se))), oeaf)]
}
cad <- cad[is.finite(obeta) & is.finite(ose)]

## 3) harmonize + format for TwoSampleMR
suppressMessages(library(TwoSampleMR))
exp_dat <- data.frame(SNP=ivs$SNP, beta.exposure=ivs$beta, se.exposure=ivs$se,
  effect_allele.exposure=ivs$ea, other_allele.exposure=ivs$oa, eaf.exposure=ivs$eaf,
  pval.exposure=ivs$p, exposure=o$trait, id.exposure=o$trait)
out_dat <- data.frame(SNP=cad$SNP, beta.outcome=cad$obeta, se.outcome=cad$ose,
  effect_allele.outcome=cad$oea, other_allele.outcome=cad$ooa, eaf.outcome=cad$oeaf,
  outcome="CAD", id.outcome="CAD")
dat <- harmonise_data(exp_dat, out_dat, action=2)
dat <- dat[dat$mr_keep, ]
cat(sprintf("[%s %s] %d SNPs after harmonization\n", o$trait, o$anc, nrow(dat)))

## 4) SENSITIVITY: IVW / weighted median / weighted mode / MR-Egger + MR-PRESSO (outlier-corrected)
suppressMessages({library(TwoSampleMR); library(MRPRESSO)})
mm <- mr(dat, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
mm$OR <- exp(mm$b); mm$lo <- exp(mm$b-1.96*mm$se); mm$hi <- exp(mm$b+1.96*mm$se)
row <- function(method,nsnp,OR=NA,lo=NA,hi=NA,P=NA) data.table(ancestry=o$anc,trait=o$trait,method=method,n_snp=nsnp,OR=OR,CI_lo=lo,CI_hi=hi,P=P)
rows <- list(row(mm$method, mm$nsnp, round(mm$OR,3), round(mm$lo,3), round(mm$hi,3), signif(mm$pval,3)))
pin <- data.frame(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
pr <- tryCatch(mr_presso(BetaOutcome="b_out",BetaExposure="b_exp",SdOutcome="se_out",SdExposure="se_exp",
      OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pin,NbDistribution=2000,SignifThreshold=0.05),
      error=function(e){cat("PRESSO err:",conditionMessage(e),"\n");NULL})
if(!is.null(pr)){
  gp <- pr[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  mrr <- pr[["Main MR results"]]
  nout <- tryCatch(length(pr[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]),error=function(e)NA)
  rows <- c(rows, list(row("MR-PRESSO global", nrow(dat), P=signif(as.numeric(gp),3))))
  cr <- mrr[nrow(mrr),]
  cb <- as.numeric(cr[["Causal Estimate"]]); cs <- as.numeric(cr[["Sd"]]); cp <- as.numeric(cr[["P-value"]])
  rows <- c(rows, list(row(sprintf("MR-PRESSO corrected (%s outliers)", ifelse(is.na(nout),"?",nout)),
       nrow(dat), round(exp(cb),3), round(exp(cb-1.96*cs),3), round(exp(cb+1.96*cs),3), signif(cp,3))))
}
outdt <- rbindlist(rows, fill=TRUE); print(outdt)
fwrite(outdt, file.path(o$out,"multiancestry_mr_sensitivity.csv"), append=file.exists(file.path(o$out,"multiancestry_mr_sensitivity.csv")))
