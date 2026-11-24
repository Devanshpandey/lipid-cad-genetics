#!/usr/bin/env Rscript
# Stage 10a - lipid / non-lipid partition of the UKB-free CAD PRS (preregistered S2).
#   --step leads     : ApoB GWAS -> QC -> distance-based lead variants (apob_leads.tsv, apob_leads.snps)
#   --step retest    : DEVIATION 2026-09-29 - re-test every lead for ApoB in statin-free participants
#                      (removes leads created by statin prescription, e.g. 9p21)
#   --step partition : lipid regions from replicated leads per window; positive/negative control checks;
#                      split PRS-CS weights (L + N = total)
#   --step assemble  : read the PLINK2 multi-column score, check L + N equals the total PRS for every
#                      person and window, write prs_components.tsv
# Usage:
#   Rscript R/10a_prs_partition.R --step leads     --config_syn FILE --apob_gwas FILE --out DIR
#   Rscript R/10a_prs_partition.R --step retest    --config_syn FILE --raw FILE --cohort FILE --pheno_dir DIR --out DIR
#   Rscript R/10a_prs_partition.R --step partition --config_syn FILE --weights FILE --out DIR
#   Rscript R/10a_prs_partition.R --step assemble  --config_syn FILE --sscore FILE --total_prs FILE --out DIR
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "lib_resilience.R")); source(file.path(here, "lib_genetics.R")); source(file.path(here, "lib_synergy.R"))
a <- parse_args(required = c("step", "config_syn", "out"))
sy <- yaml::read_yaml(a$config_syn); S2 <- sy$s2
cfg <- list(.path = normalizePath(a$config_syn), .md5 = digest::digest(file = a$config_syn, algo = "md5"), seed = sy$seed)
dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
wins <- as.character(S2$windows_kb); prim <- as.character(S2$primary_window_kb)

if (a$step == "leads") {
  check(!is.null(a$apob_gwas) && file.exists(a$apob_gwas), "--apob_gwas supplied and present")
  g <- fread(a$apob_gwas, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P"),
             colClasses = list(character = c("ID", "ALLELE0", "ALLELE1", "TEST")))
  log_msg("ApoB GWAS: %d variants, max N %d", nrow(g), max(g$N))
  q <- qc_regenie(g, max(g$N), S2$gwas_min_call_rate, S2$gwas_min_maf); rm(g)
  log_msg("ApoB GWAS QC: %s", paste(sprintf("%s=%d", names(q$counts), q$counts), collapse = " | "))
  leads <- define_loci(q$d, S2$gw_log10p, S2$locus_merge_kb * 1000)
  check(nrow(leads) >= 20, "ApoB GWAS yields at least 20 genome-wide significant loci (%d)", nrow(leads))
  fwrite(leads[, .(CHROM, GENPOS, ID, ALLELE0, ALLELE1, BETA, SE, LOG10P, n_sig_in_window)], file.path(a$out, "apob_leads.tsv"), sep = "\t")
  writeLines(leads$ID, file.path(a$out, "apob_leads.snps"))
  write_manifest("10a_leads", a$out, cfg, inputs = a$apob_gwas, outputs = file.path(a$out, c("apob_leads.tsv", "apob_leads.snps")),
                 extra = list(n_leads = nrow(leads), qc = as.list(q$counts)))
} else if (a$step == "retest") {
  R <- S2$retest
  for (k in c("raw", "cohort", "pheno_dir")) check(!is.null(a[[k]]) && file.exists(a[[k]]), "--%s supplied and present", k)
  leads <- fread(file.path(a$out, "apob_leads.tsv"), colClasses = list(character = c("ID", "ALLELE0", "ALLELE1")))
  cv <- c("age", "sex", paste0("PC", 1:10), "genotype_batch")
  coh <- fread(a$cohort, select = c("eid", "excl_genetic_qc", "statin", cv), na.strings = "NA")
  apob <- read_scalar(a$pheno_dir, R$apob_field, "apob")
  d <- merge(coh[excl_genetic_qc == 0 & statin == 0], apob, by = "eid")
  d <- d[!is.na(apob) & stats::complete.cases(d[, cv, with = FALSE])]
  d[, y := qnorm((frank(apob, ties.method = "average") - 0.5) / .N)]
  raw <- fread(a$raw); raw[, eid := as.integer(IID)]
  gcols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "eid"))
  ids <- sub("_[^_]+$", "", gcols); counted <- sub("^.*_", "", gcols)
  check(setequal(ids, leads$ID), "genotypes exported for exactly the %d leads (%d columns)", nrow(leads), length(gcols))
  L0 <- leads[match(ids, ID)]
  sgn <- ifelse(counted == L0$ALLELE1, 1, ifelse(counted == L0$ALLELE0, -1, NA_real_))
  check(!anyNA(sgn), "every exported counted allele is one of the GWAS alleles")
  miss <- sum(!d$eid %in% raw$eid)
  check(miss / nrow(d) < 0.001, "statin-free participants without genotypes are negligible (%d of %d; the sample list and .fam differ slightly)", miss, nrow(d))
  d <- d[eid %in% raw$eid]
  G <- as.matrix(raw[match(d$eid, eid), gcols, with = FALSE])
  X <- model.matrix(~ age + factor(sex) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + factor(genotype_batch), data = d)
  r <- retest_leads(d$y, G, X)
  res <- data.table(ID = ids, counted = counted, beta_statinfree = r$beta * sgn, se_statinfree = r$se, p_statinfree = r$p)
  res <- merge(leads[, .(CHROM, GENPOS, ID, gwas_beta = BETA, gwas_log10p = LOG10P)], res, by = "ID")
  res[, same_sign := sign(beta_statinfree) == sign(gwas_beta)]
  res[, keep := p_statinfree < R$p_max & (!isTRUE(R$require_same_sign) | same_sign)]
  check(mean(res$same_sign) >= R$min_frac_same_sign, "most ApoB leads replicate in sign among statin-free participants (%.1f%%)", 100 * mean(res$same_sign))
  log_msg("statin-free re-test (n = %d): %d of %d leads kept; dropped: %s", nrow(d), sum(res$keep), nrow(res),
          paste(res[keep == FALSE, sprintf("%s chr%d:%d (GWAS -log10p %.1f, statin-free p %.2g)", ID, CHROM, GENPOS, gwas_log10p, p_statinfree)], collapse = "; "))
  fwrite(res[order(CHROM, GENPOS)], file.path(a$out, "apob_leads_retest.tsv"), sep = "\t")
  write_manifest("10a_retest", a$out, cfg, inputs = c(a$raw, a$cohort), outputs = file.path(a$out, "apob_leads_retest.tsv"),
                 extra = list(n_statinfree = nrow(d), n_leads = nrow(res), n_kept = sum(res$keep), dropped = res[keep == FALSE, ID]))
} else if (a$step == "partition") {
  check(!is.null(a$weights) && file.exists(a$weights), "--weights supplied and present")
  rt <- fread(file.path(a$out, "apob_leads_retest.tsv"))
  leads <- rt[keep == TRUE]
  check(nrow(leads) >= 20, "at least 20 ApoB leads replicate among statin-free participants (%d)", nrow(leads))
  regs <- setNames(lapply(S2$windows_kb, function(kb) regions_from_leads(leads, kb * 1000)), wins)
  ctl <- rbind(data.table(name = names(S2$positive_controls), type = "positive",
                          CHR = sapply(S2$positive_controls, `[[`, 1), BP = sapply(S2$positive_controls, `[[`, 2)),
               data.table(name = names(S2$negative_controls), type = "negative",
                          CHR = sapply(S2$negative_controls, `[[`, 1), BP = sapply(S2$negative_controls, `[[`, 2)))
  for (w in wins) set(ctl, j = paste0("inside_", w), value = in_regions(ctl$CHR, ctl$BP, regs[[w]]))
  for (w in wins) check(all(ctl[type == "positive"][[paste0("inside_", w)]]), "positive controls (LDLR, PCSK9, APOB, APOE) inside lipid regions at +/- %s kb", w)
  check(!any(ctl[type == "negative"][[paste0("inside_", prim)]]), "negative controls (9p21, PHACTR1) outside lipid regions at the primary +/- %s kb", prim)
  w <- fread(a$weights, header = FALSE, col.names = c("CHR", "SNP", "BP", "A1", "A2", "BETA"))
  check(!anyDuplicated(w$SNP), "PRS-CS weights have unique SNP ids (%d)", nrow(w))
  pw <- partition_weights(w, regs)
  for (x in wins) check(isTRUE(all.equal(pw[[paste0("L", x)]] + pw[[paste0("N", x)]], w$BETA, tolerance = 0)), "L%s + N%s equals the total weight for every SNP", x, x)
  share <- rbindlist(lapply(wins, function(x) data.table(window_kb = as.integer(x), n_regions = nrow(regs[[x]]),
    region_mb = sum(regs[[x]]$end - regs[[x]]$start) / 1e6, snps_lipid = sum(pw[[paste0("L", x)]] != 0),
    share_abs_weight_lipid = sum(abs(pw[[paste0("L", x)]])) / sum(abs(w$BETA)))))
  log_msg("partition: %s", paste(sprintf("+/-%skb: %d regions, %.0f Mb, %d SNPs, %.1f%% of |weight|", share$window_kb, share$n_regions,
                                         share$region_mb, share$snps_lipid, 100 * share$share_abs_weight_lipid), collapse = " | "))
  fwrite(pw, file.path(a$out, "prs_component_weights.tsv"), sep = "\t")
  fwrite(rbindlist(regs, idcol = "window_kb"), file.path(a$out, "lipid_regions.tsv"), sep = "\t")
  fwrite(ctl, file.path(a$out, "partition_controls.tsv"), sep = "\t"); fwrite(share, file.path(a$out, "partition_summary.csv"))
  write_manifest("10a_partition", a$out, cfg, inputs = c(file.path(a$out, "apob_leads_retest.tsv"), a$weights),
                 outputs = file.path(a$out, c("prs_component_weights.tsv", "lipid_regions.tsv", "partition_controls.tsv")),
                 extra = list(n_leads_used = nrow(leads), share = share))
} else if (a$step == "assemble") {
  for (k in c("sscore", "total_prs")) check(!is.null(a[[k]]) && file.exists(a[[k]]), "--%s supplied and present", k)
  s <- fread(a$sscore); setnames(s, 1, "FID")
  tot <- fread(a$total_prs)
  check(identical(names(tot), c("FID", "IID", "PRS")), "total PRS file has FID IID PRS")
  cols <- c(paste0("L", wins, "_SUM"), paste0("N", wins, "_SUM"))
  check(all(cols %in% names(s)), "PLINK2 wrote every component sum (%s)", paste(cols, collapse = ","))
  comp <- s[, c("IID", cols), with = FALSE]; setnames(comp, cols, sub("_SUM$", "", cols))
  comp <- merge(comp, tot[, .(IID, TOTAL = PRS)], by = "IID")
  check(nrow(comp) == nrow(tot) && nrow(comp) == nrow(s), "components and total PRS cover the same %d people", nrow(comp))
  # PLINK2 prints scores to 6 significant digits, so L + N and the printed total can differ by the
  # rounding of three numbers. Bound each person's deviation by exactly that rounding.
  half6 <- function(v) 0.5 * 10^(floor(log10(pmax(abs(v), 1e-12))) - 5)
  for (x in wins) {
    Lx <- comp[[paste0("L", x)]]; Nx <- comp[[paste0("N", x)]]
    dev <- abs(Lx + Nx - comp$TOTAL); bound <- half6(Lx) + half6(Nx) + half6(comp$TOTAL) + 1e-12
    check(all(dev <= bound), "L%s + N%s reproduces the total PRS for every person within 6-significant-digit rounding (max deviation %.3g; %d exceed the bound)",
          x, x, max(dev), sum(dev > bound))
    check(max(dev) < 1e-4 * sd(comp$TOTAL), "largest L%s + N%s deviation is negligible against the PRS SD (%.3g)", x, x, max(dev) / sd(comp$TOTAL))
  }
  comp[, FID := IID]; setcolorder(comp, c("FID", "IID"))
  fwrite(comp, file.path(a$out, "prs_components.tsv"), sep = "\t")
  cr <- sapply(wins, function(x) cor(comp[[paste0("L", x)]], comp[[paste0("N", x)]]))
  log_msg("assembled %d people; cor(L, N) by window: %s", nrow(comp), paste(sprintf("%s=%.3f", wins, cr), collapse = " "))
  write_manifest("10a_assemble", a$out, cfg, inputs = c(a$sscore, a$total_prs), outputs = file.path(a$out, "prs_components.tsv"),
                 extra = list(n = nrow(comp), cor_L_N = as.list(cr)))
} else stop("unknown --step: ", a$step)
