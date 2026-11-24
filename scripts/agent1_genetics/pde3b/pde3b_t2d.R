suppressPackageStartupMessages({library(data.table); library(coloc)})
W <- Sys.getenv("SCRATCH"); W <- file.path(W, "pde3b")

tg  <- fread(file.path(W,"ukb_lipids_TRIGLY_pde3b.tsv"))
hdl <- fread(file.path(W,"ukb_lipids_HDL_C_pde3b.tsv"))
fg  <- fread(file.path(W,"fg_t2d_pde3b.tsv"))
setnames(fg, gsub("^#","",names(fg)))

cat(sprintf("UKB TG window: %d variants | FinnGen T2D window: %d variants\n", nrow(tg), nrow(fg)))
cat(sprintf("FinnGen T2D max -log10P in window: %.2f\n", max(fg$mlogp, na.rm=TRUE)))
cat(sprintf("UKB TG      max -log10P in window: %.2f  (%s)\n",
            max(tg$LOG10P,na.rm=TRUE), tg[which.max(LOG10P)]$ID))

# FinnGen can list several rsids per row; expand and keep the first usable one
fg <- fg[!is.na(beta) & !is.na(sebeta) & sebeta > 0 & rsids != ""]
fg[, rsid := tstrsplit(rsids, ",", fixed=TRUE)[[1]]]
fg <- fg[!duplicated(rsid)]

tg <- tg[!is.na(BETA) & !is.na(SE) & SE > 0]
m <- merge(tg[, .(rsid=ID, a0=toupper(ALLELE0), a1=toupper(ALLELE1),
                  b_tg=BETA, se_tg=SE, f_tg=A1FREQ, n_tg=N, lp_tg=LOG10P)],
           fg[, .(rsid, ref=toupper(ref), alt=toupper(alt),
                  b_t2d=beta, se_t2d=sebeta, f_t2d=af_alt, lp_t2d=mlogp)], by="rsid")
cat(sprintf("merged on rsID: %d variants\n", nrow(m)))

# Harmonise to the UKB effect allele (REGENIE ALLELE1); FinnGen beta is for alt
m[, align := fifelse(a1==alt & a0==ref, 1, fifelse(a1==ref & a0==alt, -1, NA_real_))]
cat(sprintf("  aligned %d, flipped %d, dropped %d\n",
            sum(m$align==1,na.rm=TRUE), sum(m$align==-1,na.rm=TRUE), sum(is.na(m$align))))
m <- m[!is.na(align)]
m[, b_t2d := b_t2d*align]
m[, f_t2d := fifelse(align==-1, 1-f_t2d, f_t2d)]
m[, maf := pmin(f_tg, 1-f_tg)]
m <- m[maf > 0.001]

NC <- 71728; NK <- 369007; NT <- NC+NK
D1 <- list(beta=m$b_tg,  varbeta=m$se_tg^2,  snp=m$rsid, MAF=m$maf, N=max(m$n_tg), type="quant")
D2 <- list(beta=m$b_t2d, varbeta=m$se_t2d^2, snp=m$rsid, MAF=m$maf, N=NT, s=NC/NT, type="cc")
r <- suppressWarnings(coloc.abf(D1, D2))
s <- r$summary
cat(sprintf("\n=== coloc: UKB triglycerides vs FinnGen T2D at PDE3B (%d SNPs) ===\n", s["nsnps"]))
cat(sprintf("  PP.H0=%.3f  PP.H1(TG only)=%.3f  PP.H2(T2D only)=%.3f  PP.H3(distinct)=%.3f  PP.H4(shared)=%.3f\n",
            s["PP.H0.abf"],s["PP.H1.abf"],s["PP.H2.abf"],s["PP.H3.abf"],s["PP.H4.abf"]))

# cis-MR: Wald ratio at the lead TG variant (single instrument, no LD panel needed)
lead <- m[which.max(lp_tg)]
wald <- lead$b_t2d/lead$b_tg
se_w <- abs(lead$se_t2d/lead$b_tg)
cat(sprintf("\n=== cis-MR (Wald ratio, lead TG variant %s) ===\n", lead$rsid))
cat(sprintf("  TG   beta=%+.4f SD (se %.4f, -log10P %.1f)\n", lead$b_tg, lead$se_tg, lead$lp_tg))
cat(sprintf("  T2D  beta=%+.4f log-OR (se %.4f, -log10P %.2f)\n", lead$b_t2d, lead$se_t2d, lead$lp_t2d))
cat(sprintf("  T2D OR per 1 SD LOWER TG via PDE3B = %.3f (95%% CI %.3f-%.3f), P=%.3g\n",
            exp(-wald), exp(-wald-1.96*se_w), exp(-wald+1.96*se_w),
            2*pnorm(-abs(wald/se_w))))
bl <- 0.4227  # TG SD lowering achieved by PDE3B LoF burden
cat(sprintf("  At the TG reduction PDE3B LoF achieves (%.3f SD): OR = %.3f (95%% CI %.3f-%.3f)\n",
            bl, exp(-wald*bl), exp((-wald-1.96*se_w)*bl), exp((-wald+1.96*se_w)*bl)))
fwrite(m, file.path(W,"pde3b_tg_vs_t2d_harmonised.csv"))
