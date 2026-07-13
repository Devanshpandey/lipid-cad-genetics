#!/usr/bin/env Rscript
# (A) Multi-cohort MR: UKB lipid instruments -> CARDIoGRAM CAD (3rd cohort after FinnGen)
# (B) Trans-ancestry transferability: UKB lipid instruments in GLGC EUR/AFR/HIS
suppressMessages(library(data.table))
OUT <- "/path/to/cad-genetics/results/agent1_genetics/strengthen/external_mr"
ins <- fread(file.path(OUT,"ukb_instruments.csv"))   # rsID, effA, othA, beta, se, trait

## ---------- (A) UKB instruments -> CARDIoGRAM CAD ----------
cad <- fread(file.path(OUT,"cardiogram_ins.txt"))
setnames(cad, c("markername","effect_allele","noneffect_allele","beta","se_dgc"),
              c("rsID","cad_eff","cad_oth","cad_beta","cad_se"), skip_absent=TRUE)
cad <- cad[, .(rsID, cad_eff=toupper(cad_eff), cad_oth=toupper(cad_oth),
               cad_beta=as.numeric(cad_beta), cad_se=as.numeric(cad_se))]
ivw <- function(tr){
  d <- merge(ins[trait==tr], cad, by="rsID")
  d[, ea:=toupper(effA)][, oa:=toupper(othA)]
  d[, bcad := fifelse(cad_eff==ea, cad_beta, fifelse(cad_eff==oa, -cad_beta, NA_real_))]
  d <- d[!is.na(bcad) & !is.na(cad_se)]
  w <- 1/d$cad_se^2; b <- sum(w*d$beta*d$bcad)/sum(w*d$beta^2); se <- sqrt(1/sum(w*d$beta^2))
  data.table(trait=tr, n_iv=nrow(d), OR=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), P=2*pnorm(-abs(b/se)))
}
cat("=== (A) MULTI-COHORT MR: UKB lipid instruments -> CARDIoGRAM CAD (per SD) ===\n")
A <- rbindlist(lapply(c("LDL_C","HDL_C","TRIGLY"), ivw))
A[, `:=`(OR=round(OR,3), lo=round(lo,3), hi=round(hi,3), P=signif(P,3))]
print(A); fwrite(A, file.path(OUT,"mr_cardiogram.csv"))

## ---------- (B) Trans-ancestry transferability (GLGC) ----------
glgc <- function(trait_file){
  g <- fread(file.path(OUT, trait_file))
  g[, .(rsID, ALT=toupper(ALT), REF=toupper(REF), eff=as.numeric(EFFECT_SIZE), se=as.numeric(SE))]
}
align_to_ukb <- function(u, g){  # align GLGC effect to UKB effect allele
  m <- merge(u, g, by="rsID")
  m[, ea:=toupper(effA)][, oa:=toupper(othA)]
  m[, g_aln := fifelse(ALT==ea, eff, fifelse(ALT==oa, -eff, NA_real_))]
  m[!is.na(g_aln)]
}
cat("\n=== (B) TRANS-ANCESTRY LIPID-EFFECT TRANSFERABILITY (GLGC; incl. MVP) ===\n")
map <- list(LDL_C=c(EUR="LDL_INV_EUR_HRC_1KGP3_others_ALL", AFR="LDL_INV_AFR_HRC_1KGP3_others_ALL", HIS="LDL_INV_HIS_1KGP3_ALL"),
            TRIGLY=c(EUR="logTG_INV_EUR_HRC_1KGP3_others_ALL", AFR="logTG_INV_AFR_HRC_1KGP3_others_ALL", HIS="logTG_INV_HIS_1KGP3_ALL"))
Bres <- list()
for (tr in names(map)) {
  u <- ins[trait==tr]
  for (anc in names(map[[tr]])) {
    g <- glgc(paste0("glgc_", map[[tr]][anc], ".txt"))
    m <- align_to_ukb(u, g)
    r <- cor(m$beta, m$g_aln); sc <- mean(sign(m$beta)==sign(m$g_aln))
    Bres[[paste(tr,anc)]] <- data.table(trait=tr, ancestry=anc, n_iv=nrow(m),
        effect_r=round(r,3), sign_concordance=round(sc,3))
    cat(sprintf("  %-7s %s: n=%d  effect r=%.3f  sign-concordance=%.1f%%\n", tr, anc, nrow(m), r, 100*sc))
  }
}
fwrite(rbindlist(Bres), file.path(OUT,"transferability_glgc.csv"))

## prioritized loci lead variants across ancestries (LDL)
cat("\n=== prioritized LDL loci: effect across ancestries (aligned to UKB LDL-lowering) ===\n")
leads <- c(PCSK9="rs11591147", SORT1="rs12740374", LDLR="rs6511720", APOB="rs1367117", APOE="rs7412")
uL <- ins[trait=="LDL_C"]
tab <- data.table(gene=names(leads), rsID=leads)
for (anc in c("EUR","AFR","HIS")) {
  g <- glgc(paste0("glgc_", map[["LDL_C"]][anc], ".txt")); m <- align_to_ukb(uL, g)
  tab[[anc]] <- round(m$g_aln[match(leads, m$rsID)], 3)
}
tab[["UKB"]] <- round(uL$beta[match(leads, uL$rsID)], 3)
print(tab); fwrite(tab, file.path(OUT,"prioritized_loci_ancestry.csv"))
cat("\nDONE\n")
