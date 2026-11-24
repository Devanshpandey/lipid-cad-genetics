#!/usr/bin/env Rscript
# ============================================================================
# Stroke-subtype contrasts and Lp(a) instrument sensitivity (MEGASTROKE)
#
# Reviewer point: "significant in one subtype and non-significant in another is
# not, by itself, evidence that their effects differ." This computes the formal
# between-subtype contrasts the claim needs.
#
#   H0: beta_LAS == beta_CES        H0: beta_LAS == beta_SVS
#
# The subtype GWAS share controls, so their MR estimates are positively
# correlated and an independence-assuming SE of the difference is too large.
# We bootstrap over INSTRUMENTS, refitting every subtype model on each resample,
# which propagates the covariance induced by the shared instrument set. It does
# not capture the additional positive correlation induced by shared controls,
# which would shrink the SE of the difference further -- so the test reported
# here is conservative in that direction.
#
# Also runs the Lp(a) instrument sensitivities the reviewer asked for:
#   (i)  single variant rs10455872 alone
#   (ii) the cis instrument pruned at the stricter r2 < 0.001
#   (iii) maximum observed pairwise r2 within the cis instrument set
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(optparse)})
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat_dir", type = "character"),
  make_option("--megastroke_dir", type = "character"),
  make_option("--ld_ref",  type = "character", default = ""),
  make_option("--plink",   type = "character", default = "plink"),
  make_option("--nboot",   type = "integer",   default = 1000L),
  make_option("--out_prefix", type = "character")
)))

LPA_CHR <- "6"; LPA_LO <- 160500000L; LPA_HI <- 161700000L
SUB <- c(LAS = "MEGASTROKE.3.LAS.EUR.out", CES = "MEGASTROKE.4.CES.EUR.out",
         SVS = "MEGASTROKE.5.SVS.EUR.out", AIS = "MEGASTROKE.2.AIS.EUR.out",
         AS  = "MEGASTROKE.1.AS.EUR.out")

message("[contrast] loading UKB Lp(a) sumstats")
lpa <- fread(cmd = paste("gunzip -c", shQuote(file.path(opt$sumstat_dir, "lipids_LPA.txt.gz"))))
lpa <- lpa[as.character(CHROM) == LPA_CHR & GENPOS >= LPA_LO & GENPOS <= LPA_HI &
           LOG10P > -log10(5e-8) & !is.na(BETA) & !is.na(SE)]
setorder(lpa, -LOG10P)
message(sprintf("[contrast] %d genome-wide-significant variants in the LPA cis window", nrow(lpa)))

# ---- LD-aware pruning with PLINK, ranked on LOG10P (pval underflows) --------
prune_r2 <- function(dt, r2) {
  if (!nzchar(opt$ld_ref)) stop("--ld_ref required")
  tmp <- tempfile()
  # monotone rescale of LOG10P onto a representable p-value range; PLINK ranks on
  # P and exact zeros would tie, which is the bug that lost rs10455872 before.
  z <- dt$LOG10P
  p <- 10^(-(1 + (z - min(z)) / (max(z) - min(z)) * 299))
  fwrite(data.table(SNP = dt$ID, P = p), paste0(tmp, ".assoc"), sep = "\t")
  system2(opt$plink, c("--bfile", opt$ld_ref, "--clump", paste0(tmp, ".assoc"),
                       "--clump-p1", "1", "--clump-r2", as.character(r2),
                       "--clump-kb", "1000", "--out", tmp), stdout = NULL, stderr = NULL)
  cl <- paste0(tmp, ".clumped")
  if (!file.exists(cl)) return(character(0))
  fread(cl)$SNP
}

message("[contrast] pruning cis instrument at r2<0.01 and r2<0.001")
ivs_01  <- prune_r2(lpa, 0.01)
ivs_001 <- prune_r2(lpa, 0.001)
message(sprintf("[contrast] r2<0.01: %d variants | r2<0.001: %d variants",
                length(ivs_01), length(ivs_001)))

# ---- maximum observed pairwise r2 inside the primary cis instrument set -----
max_r2 <- NA_real_
if (length(ivs_01) > 1) {
  tmp <- tempfile(); writeLines(ivs_01, paste0(tmp, ".snps"))
  system2(opt$plink, c("--bfile", opt$ld_ref, "--extract", paste0(tmp, ".snps"),
                       "--r2", "--ld-window-r2", "0", "--ld-window", "99999",
                       "--ld-window-kb", "2000", "--out", tmp), stdout = NULL, stderr = NULL)
  if (file.exists(paste0(tmp, ".ld"))) {
    ld <- fread(paste0(tmp, ".ld"))
    max_r2 <- max(ld$R2, na.rm = TRUE)
    message(sprintf("[contrast] max pairwise r2 among the %d cis instruments: %.4f",
                    length(ivs_01), max_r2))
  }
}

# ---- harmonize one MEGASTROKE file to a given instrument set ----------------
load_out <- function(f, ids) {
  d <- fread(file.path(opt$megastroke_dir, f), showProgress = FALSE)
  setnames(d, c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value"),
           c("SNP", "EA", "OA", "BETA_O", "SE_O", "P_O"), skip_absent = TRUE)
  d <- d[SNP %in% ids]
  d[, `:=`(EA = toupper(EA), OA = toupper(OA))]
  d[!duplicated(SNP)]
}

ivw <- function(bx, bxse, by, byse) {
  w <- 1 / byse^2
  b <- sum(w * bx * by) / sum(w * bx^2)
  se <- sqrt(1 / sum(w * bx^2))
  list(b = b, se = se, p = 2 * pnorm(-abs(b / se)), n = length(bx))
}

build <- function(ids) {
  ex <- lpa[ID %in% ids, .(SNP = ID, EA = toupper(ALLELE1), OA = toupper(ALLELE0),
                           BETA_X = BETA, SE_X = SE)]
  outs <- lapply(SUB, load_out, ids = ids)
  common <- Reduce(intersect, c(list(ex$SNP), lapply(outs, `[[`, "SNP")))
  ex <- ex[SNP %in% common][order(SNP)]
  aligned <- lapply(outs, function(o) {
    o <- o[SNP %in% common][order(SNP)]
    m <- merge(o, ex[, .(SNP, EA_X = EA, OA_X = OA)], by = "SNP")
    setorder(m, SNP)
    flip <- fifelse(m$EA == m$EA_X, 1, fifelse(m$EA == m$OA_X, -1, NA_real_))
    data.table(SNP = m$SNP, BETA_O = m$BETA_O * flip, SE_O = m$SE_O)
  })
  ok <- Reduce(`&`, lapply(aligned, function(a) !is.na(a$BETA_O)))
  ex <- ex[ok]; aligned <- lapply(aligned, function(a) a[ok])
  list(ex = ex, out = aligned)
}

report <- function(tag, ids) {
  d <- build(ids)
  n <- nrow(d$ex)
  message(sprintf("\n[contrast] === %s : %d instruments matched across all 5 outcomes ===", tag, n))
  est <- rbindlist(lapply(names(SUB), function(k) {
    r <- ivw(d$ex$BETA_X, d$ex$SE_X, d$out[[k]]$BETA_O, d$out[[k]]$SE_O)
    data.table(instrument = tag, subtype = k, n_iv = n, OR = exp(r$b),
               lo = exp(r$b - 1.96 * r$se), hi = exp(r$b + 1.96 * r$se),
               b = r$b, se = r$se, P = r$p)
  }))
  print(est[, .(subtype, n_iv, OR = round(OR, 3), P = signif(P, 3))])

  # bootstrap over instruments, refitting every subtype on each resample
  B <- opt$nboot
  bmat <- matrix(NA_real_, B, length(SUB), dimnames = list(NULL, names(SUB)))
  for (i in seq_len(B)) {
    s <- sample.int(n, n, replace = TRUE)
    for (k in names(SUB))
      bmat[i, k] <- ivw(d$ex$BETA_X[s], d$ex$SE_X[s],
                        d$out[[k]]$BETA_O[s], d$out[[k]]$SE_O[s])$b
  }
  contr <- rbindlist(lapply(c("CES", "SVS", "AIS", "AS"), function(k) {
    diff  <- est[subtype == "LAS", b] - est[subtype == k, b]
    se_i  <- sqrt(est[subtype == "LAS", se]^2 + est[subtype == k, se]^2)
    se_b  <- sd(bmat[, "LAS"] - bmat[, k])
    r     <- cor(bmat[, "LAS"], bmat[, k])
    # se_b == 0 means every resample returned the same variant (single-variant
    # design). The contrast is undefined, not infinitely significant: emit NA so a
    # degenerate design can never be mistaken for P = 0 downstream.
    ok <- is.finite(se_b) && se_b > 0
    data.table(instrument = tag, contrast = paste0("LAS vs ", k), delta_logOR = diff,
               se_independent = se_i, se_bootstrap = se_b, r_bootstrap = r,
               z = if (ok) diff / se_b else NA_real_,
               P = if (ok) 2 * pnorm(-abs(diff / se_b)) else NA_real_)
  }))
  print(contr[, .(contrast, delta = round(delta_logOR, 3), se_ind = round(se_independent, 4),
                  se_boot = round(se_bootstrap, 4), r = round(r_bootstrap, 2),
                  z = round(z, 2), P = signif(P, 3))])

  # Cochran's Q across the three TOAST subtypes, bootstrap covariance
  idx <- c("LAS", "CES", "SVS")
  V <- cov(bmat[, idx]); bvec <- est[match(idx, subtype), b]
  # A single-variant design gives every resample the same variant, so V is singular
  # and Q is undefined. That is a property of the design, not an error: report it as
  # undefined rather than letting solve() abort the script before anything is written.
  qq <- tryCatch({
    Vi <- solve(V); w <- rowSums(Vi)
    bbar <- sum(w * bvec) / sum(w)
    Q <- as.numeric(t(bvec - bbar) %*% Vi %*% (bvec - bbar))
    list(Q = Q, P = pchisq(Q, df = 2, lower.tail = FALSE))
  }, error = function(e) {
    message(sprintf("[contrast] %s: Q undefined (%s)", tag, conditionMessage(e)))
    list(Q = NA_real_, P = NA_real_)
  })
  Q <- qq$Q; Qp <- qq$P
  if (!is.na(Q))
    message(sprintf("[contrast] %s: Cochran Q across LAS/CES/SVS = %.2f, df=2, P=%.3g", tag, Q, Qp))
  list(est = est, contr = contr,
       q = data.table(instrument = tag, Q = Q, df = 2L, P = Qp, max_pairwise_r2 = max_r2))
}

res_main  <- report("cis r2<0.01", ivs_01)
res_strict<- report("cis r2<0.001", ivs_001)
res_single<- report("rs10455872 only", "rs10455872")

fwrite(rbindlist(list(res_main$est, res_strict$est, res_single$est)),
       paste0(opt$out_prefix, "_estimates.csv"))
fwrite(rbindlist(list(res_main$contr, res_strict$contr, res_single$contr)),
       paste0(opt$out_prefix, "_contrasts.csv"))
fwrite(rbindlist(list(res_main$q, res_strict$q, res_single$q)),
       paste0(opt$out_prefix, "_heterogeneity.csv"))
message("\n[contrast] done")
