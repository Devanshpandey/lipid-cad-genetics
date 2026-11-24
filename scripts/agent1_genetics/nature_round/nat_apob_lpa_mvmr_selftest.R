#!/usr/bin/env Rscript
# ============================================================
# Selftest for nat_apob_lpa_mvmr.R. SYNTHETIC DATA ONLY, seed 42.
#
# The test that matters: the same script, on synthetic GWAS files on disk, must reach
# DIFFERENT verdicts when a planted Lp(a) effect survives conditioning on ApoB and when
# Lp(a) acts only through its ApoB content, with the LDL-C analogue passing the positive
# control in both. A script that returned SUPPORTS for both would be useless however
# clean its output looked. Also checked: the published implementations (run_mvmr.R,
# mvmr_covariance.R, run_external_mr.R) are RUN on the same files and must return the
# same instruments and estimates, so "reproduces the published logic" is tested, not
# asserted; and the attenuation interval is tested for coverage.
#
# Usage: Rscript nat_apob_lpa_mvmr_selftest.R --workdir /path/to/scratch
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--workdir", type = "character"),
  make_option("--threads", type = "integer", default = 2L),
  make_option("--nboot", type = "integer", default = 1000L),
  make_option("--cov_reps", type = "integer", default = 200L, help = "Monte Carlo replicates for the coverage check")
)))
stopifnot(!is.null(opt$workdir))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
HERE <- dirname(this_file)
MAIN <- file.path(HERE, "nat_apob_lpa_mvmr.R")
R_DIR <- normalizePath(file.path(HERE, "..", "R"))
WORK <- opt$workdir
unlink(WORK, recursive = TRUE); dir.create(WORK, recursive = TRUE)

source(MAIN)          # library mode: sys.nframe() > 0, so main() does not run
set.seed(42)

N_PASS <- 0L; N_FAIL <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok)
  cat(sprintf("%-4s %s\n", if (ok) "PASS" else "FAIL", label))
  if (ok) N_PASS <<- N_PASS + 1L else N_FAIL <<- N_FAIL + 1L
  invisible(ok)
}
near <- function(a, b, tol) isTRUE(!is.na(a) && !is.na(b) && abs(a - b) <= tol)

# ---- 1. unit tests of the pieces that have bitten this project before ------------
cat("\n== unit tests ==\n")

# LOG10P ranking: the underflow bug picked a low-position variant over the real lead.
d <- data.table(ID = c("rs531596603", "rs10455872", "rsX"), CHROM = "6", GENPOS = c(160179747L, 161010118L, 100000000L),
                BETA = c(.1, 1, .2), SE = .01, LOG10P = c(438, 14632, 50))
check("rank on LOG10P keeps the true lead, not the first underflowing row in file order",
      identical(greedy_clump(d), c("rs10455872", "rsX")))

# The recycling bug: same-chromosome kept variants not a prefix of the kept list.
d <- data.table(ID = paste0("v", 1:6), CHROM = c("1", "2", "1", "2", "1", "1"),
                GENPOS = c(50e6, 50e6, 80e6, 90e6, 50e6 + 26, 80e6 + 5e6), BETA = 1, SE = 1, LOG10P = 100:95)
got <- greedy_clump(d)
check("distance prune has no within-window survivor on interleaved chromosomes",
      identical(got, c("v1", "v2", "v3", "v4")))

# Harmonization: exact, flipped, strand-flipped, palindrome by frequency, ambiguous
# palindrome dropped, multi-allelic rows reduced to the matching one.
ref <- data.table(ID = paste0("s", 1:7), rA1 = c("A", "A", "A", "A", "A", "C", "A"), rA0 = c("G", "G", "G", "T", "T", "T", "C"),
                  rEAF = c(.2, .2, .2, .2, .5, .3, .3))
ext <- data.table(ID = c("s1", "s2", "s3", "s4", "s5", "s6", "s6", "s7"),
                  oEA = c("A", "G", "T", "A", "A", "T", "C", "G"), oOA = c("G", "A", "C", "T", "T", "G", "T", "T"),
                  oBETA = c(1, 1, 1, 1, 1, 9, 1, 1), oSE = 1, oEAF = c(.2, .8, .2, .8, .5, .1, .3, .7))
h <- harmonize_external(ref, ext); setkey(h, ID)
check("external harmonization: flip, strand flip, palindrome sign, ambiguous drop, multi-allelic row",
      identical(h$ID, c("s1", "s2", "s3", "s4", "s6", "s7")) && all(h[.("s1", "s3", "s6"), b] == 1) && all(h[.("s2", "s4", "s7"), b] == -1))

# Estimator and guards.
set.seed(42)
gX <- function(n, K) { X <- matrix(rnorm(n * K, 0, .05), n, K, dimnames = list(paste0("r", 1:n), paste0("T", 1:K))); X }
mkD <- function(X, win = 1L, y = NULL) {
  n <- nrow(X); K <- ncol(X)
  list(BX = X, BXSE = matrix(.003, n, K, dimnames = dimnames(X)), LOG10P = matrix(50, n, K, dimnames = dimnames(X)),
       BY = if (is.null(y)) as.numeric(X %*% rep(.2, K) + rnorm(n, 0, .003)) else y, BYSE = rep(.003, n),
       ids = rownames(X), in_window = seq_len(n) <= win)
}
X <- gX(60, 4); colnames(X) <- c("APOB", "HDL_C", "TRIGLY", "LPA")
P <- diag(4); dimnames(P) <- list(colnames(X), colnames(X)); P["APOB", "HDL_C"] <- P["HDL_C", "APOB"] <- -.3
D <- mkD(X)
r <- fit_model(D, colnames(X), P)
check("fit_model returns status ok on a sound design", all(r$status == "ok") && !anyNA(r$b))
cat("  (mr_mvivw agreement with the internal estimator is enforced inside fit_model on every fit)\n")
K2 <- c("APOB", "LPA")
r2 <- fit_model(subset_design(D, K2), K2, P[K2, K2])
check("two-exposure model gives conditional F and Q_A under phenocov", all(r2$status == "ok") && !anyNA(r2$condF_phenocov) && !anyNA(r2$Q_A_p_phenocov))
Xd <- X; Xd[, "HDL_C"] <- Xd[, "APOB"] * 2
check("rank-deficient design gives NA, not a number", all(is.na(fit_model(mkD(Xd), colnames(X), P)$OR)) &&
        all(fit_model(mkD(Xd), colnames(X), P)$status == "rank_deficient"))
check("too few instruments gives NA", all(grepl("too_few", fit_model(mkD(X[1:5, ]), colnames(X), P)$status)))
Dn <- D; Dn$LOG10P[, "LPA"] <- 1
check("exposure with no GWS instrument gives NA", all(grepl("not_instrumented:LPA", fit_model(Dn, colnames(X), P)$status)))
Dw <- D; Dw$in_window[] <- FALSE
check("LPA locus absent from the outcome gives NA", all(fit_model(Dw, colnames(X), P)$status == "lpa_locus_absent"))
check("a covariance matrix is refused where phenocov needs a correlation matrix",
      inherits(try(fit_model(D, colnames(X), P * 2400), silent = TRUE), "try-error"))

# The readout, as a truth table. b_ref = log(1.15).
br <- log(1.15)
tt <- function(frac, lo = 1.02, hi = 1.5) cohort_category(br, br * (1 - frac), lo, hi)
check("readout: 0% fall = SUPPORTS", tt(0) == "SUPPORTS")
check("readout: exactly 25% fall = SUPPORTS (within 25%)", cohort_category(1, .75, 1.1, 1.5) == "SUPPORTS")
check("readout: 30% fall = ATTENUATED", tt(.30) == "ATTENUATED")
check("readout: exactly 50% fall = ATTENUATED", cohort_category(1, .5, 1.1, 1.5) == "ATTENUATED")
check("readout: 60% fall = WEAKENED", tt(.60) == "WEAKENED")
check("readout: CI includes 1 = WEAKENED whatever the fall", tt(.05, .99, 1.3) == "WEAKENED")
check("readout: sign flip with CI excluding 1 = WEAKENED", tt(2, .7, .9) == "WEAKENED")
check("readout: rise over 25% is not called SUPPORTS", tt(-.40) == "OUTSIDE_RULES_INCREASE")
check("readout: NA estimate = NOT_READABLE", cohort_category(br, NA_real_, NA_real_, NA_real_) == "NOT_READABLE")
check("readout: non-positive reference = NOT_READABLE", cohort_category(-.1, .1, 1.05, 1.3) == "NOT_READABLE")
check("readout: WEAKENED in either cohort wins", combine_categories(c("SUPPORTS", "WEAKENED")) == "WEAKENED")
check("readout: an unreadable cohort blocks SUPPORTS", combine_categories(c("SUPPORTS", "NOT_READABLE")) == "NOT_READABLE")
check("readout: ATTENUATED beats SUPPORTS", combine_categories(c("ATTENUATED", "SUPPORTS")) == "ATTENUATED")

# The window exclusion: under cis_dense a trait's own window lead outside the cis list would
# otherwise stay in, exempt from pruning. This is the only variant where the exclusion bites.
mk_rows <- function(ids, pos, lp) data.table(ID = ids, CHROM = "6", GENPOS = pos, LOG10P = lp)
rw <- list(APOB = mk_rows(c("rsWIN", "rsCL", "rsFAR"), c(161000000L, 160900000L, 100000000L), c(300, 200, 100)),
           LPA  = mk_rows(c("rsWIN", "rsCL", "rsFAR"), c(161000000L, 160900000L, 100000000L), c(50, 900, 1)))
fl <- list(APOB = c("rsWIN", "rsFAR"), LPA = "rsCL"); nw <- list(APOB = "rsFAR", LPA = character(0))
check("cis_dense keeps ApoB's own window lead; cis_noLPAwin_nonLPA removes it and keeps the cis list",
      setequal(select_instruments("cis_dense", c("APOB", "LPA"), rw, fl, nw, "rsCL"), c("rsWIN", "rsCL", "rsFAR")) &&
        setequal(select_instruments("cis_noLPAwin_nonLPA", c("APOB", "LPA"), rw, fl, nw, "rsCL"), c("rsCL", "rsFAR")))

# ---- 2. synthetic world on disk --------------------------------------------------
# Effects are in SD units per allele, SEs from N = 350,000, so the scale matches the
# real GWAS. LDL-C and ApoB are nearly collinear across the polygenic loci, as in the
# real data. At the LPA locus every variant moves ApoB by c_apob and LDL-C by c_ldl per
# unit of Lp(a); that ratio is what separates the scenarios.
SCENARIOS <- list(
  # Lp(a) has its own effect: ApoB conditioning should change nothing.
  SURVIVES = list(driver = "LDL_C", theta = log(1.25), direct = log(1.15), c_apob = .10, c_ldl = .10),
  # Lp(a) acts ONLY through apoB-100 on the particle. The LDL-C model still shows OR 1.15,
  # because LDL-C captures less of that ApoB content than ApoB itself.
  APOB_ONLY = list(driver = "APOB", theta = log(1.25), direct = 0, c_apob = .70, c_ldl = .70 - log(1.15) / log(1.25)),
  # 35% of the LDL-model Lp(a) log OR is ApoB content.
  PARTIAL = list(driver = "APOB", theta = log(1.25), direct = .65 * log(1.15), c_apob = .30, c_ldl = .30 - .35 * log(1.15) / log(1.25)),
  # Does not match the published control values; must stop the run.
  BAD_CONTROL = list(driver = "LDL_C", theta = log(1.40), direct = log(1.30), c_apob = .10, c_ldl = .10)
)
TRAITS <- c("LDL_C", "APOB", "HDL_C", "TRIGLY", "LPA")
THETA_HDL <- -.07; THETA_TG <- .10; N_GWAS <- 350000

make_true_effects <- function() {
  set.seed(42)
  nloc <- 13L
  loci <- rbindlist(lapply(1:22, function(ch) data.table(CHROM = ch, locus_pos = as.integer(5e6 + 12e6 * (seq_len(nloc) - 1L)))))
  L <- nrow(loci)
  g_ldl <- rnorm(L, 0, .05)
  # A dozen genome-wide-significant trans Lp(a) variants, so noLPAwin_all is estimable.
  # The rest are tiny: the locus must dominate the Lp(a) coefficient, as it does in reality.
  g_lpa <- rnorm(L, 0, .003); g_lpa[sample(L, 12)] <- .018 * sample(c(-1, 1), 12, TRUE)
  lead <- data.table(loci, GENPOS = loci$locus_pos, role = "lead", LDL_C = g_ldl, APOB = g_ldl + rnorm(L, 0, .004),
                     HDL_C = -.3 * g_ldl + rnorm(L, 0, .03), TRIGLY = .2 * g_ldl + rnorm(L, 0, .035), LPA = g_lpa)
  fl1 <- copy(lead)[, `:=`(GENPOS = locus_pos + 60000L, role = "flank1")]
  # On a fifth of loci the ApoB lead is a different variant from the LDL-C lead, so the
  # two instrument sets differ and "matched" means something.
  sw <- rbinom(L, 1, .2) == 1
  fl1[, `:=`(LDL_C = .7 * LDL_C, HDL_C = .7 * HDL_C, TRIGLY = .7 * TRIGLY, LPA = .7 * LPA,
             APOB = fifelse(sw, 1.2 * APOB, .7 * APOB))]
  fl2 <- copy(lead)[, `:=`(GENPOS = locus_pos - 90000L, role = "flank2", LDL_C = .4 * LDL_C, APOB = .4 * APOB,
                           HDL_C = .4 * HDL_C, TRIGLY = .4 * TRIGLY, LPA = .4 * LPA)]
  nul <- data.table(CHROM = sample(1:22, 1500, TRUE), locus_pos = 0L, GENPOS = sample(1e6:2e8, 1500), role = "null",
                    LDL_C = 0, APOB = 0, HDL_C = 0, TRIGLY = 0, LPA = 0)
  poly <- rbindlist(list(lead, fl1, fl2, nul))
  poly[, `:=`(ID = sprintf("rs%d", 1000000L + seq_len(.N)), is_lpa_locus = FALSE)]
  w <- data.table(CHROM = 6L, locus_pos = 0L,
                  GENPOS = c(160179747L, 160961137L, 161010118L, 160700000L, 160850000L, 161150000L, 161300000L, 161500000L),
                  role = "lpa", ID = c("rs531596603", "rs3798220", "rs10455872", "rsCIS1", "rsCIS2", "rsCIS3", "rsCIS4", "rsCIS5"),
                  g = c(.10, .85, 1.0, .35, .25, -.20, .15, .30),
                  is_lpa_locus = TRUE)
  list(poly = poly, lpa = w)
}
TRUE_FX <- make_true_effects()
CIS_IDS <- c("rs10455872", "rs3798220", "rsCIS1", "rsCIS2", "rsCIS3", "rsCIS4", "rsCIS5")

alleles_for <- function(n) {
  set.seed(42)
  pal <- rbinom(n, 1, .08) == 1
  pairs <- list(c("A", "G"), c("C", "T"), c("A", "C"), c("G", "T"), c("G", "A"), c("T", "C"))
  a <- t(vapply(seq_len(n), function(i) if (pal[i]) sample(list(c("A", "T"), c("T", "A"), c("C", "G"), c("G", "C")), 1)[[1]] else pairs[[sample(6, 1)]], c("", "")))
  f <- fifelse(pal & rbinom(n, 1, .5) == 1, runif(n, .45, .55), runif(n, .1, .4))
  list(A1 = a[, 1], A0 = a[, 2], f = f, pal = pal)
}

write_world <- function(dir, sc, seed = 42L, outcomes = c("CAD", "MI", "REVASC", "STROKE"), with_externals = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  poly <- copy(TRUE_FX$poly)
  lp <- copy(TRUE_FX$lpa)
  lp[, `:=`(LPA = g, APOB = sc$c_apob * g, LDL_C = sc$c_ldl * g, HDL_C = 0, TRIGLY = 0)]
  a <- rbindlist(list(poly[, .(CHROM, GENPOS, ID, role, LDL_C, APOB, HDL_C, TRIGLY, LPA)],
                      lp[, .(CHROM, GENPOS, ID, role, LDL_C, APOB, HDL_C, TRIGLY, LPA)]))
  setorder(a, CHROM, GENPOS)
  al <- alleles_for(nrow(a))
  a[, `:=`(A1 = al$A1, A0 = al$A0, f = al$f, pal = al$pal)]
  a[ID == "rs10455872", `:=`(A1 = "G", A0 = "A", f = .25, pal = FALSE)]
  a[ID == "rs3798220", `:=`(A1 = "C", A0 = "T", f = .02, pal = FALSE)]
  a[ID == "rs531596603", `:=`(A1 = "G", A0 = "T", f = .30, pal = FALSE)]
  a[, se := 1 / sqrt(2 * f * (1 - f) * N_GWAS)]
  log10p <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
  to_regenie <- function(b, s) data.table(CHROM = a$CHROM, GENPOS = a$GENPOS, ID = a$ID, ALLELE0 = a$A0, ALLELE1 = a$A1, A1FREQ = a$f,
                                          N = N_GWAS, TEST = "ADD", BETA = b, SE = s, CHISQ = (b / s)^2, LOG10P = log10p(b / s), EXTRA = "NA")
  for (tr in TRAITS) {
    b <- a[[tr]] + a$se * rnorm(nrow(a))
    fwrite(to_regenie(b, a$se), file.path(dir, sprintf("lipids_%s.txt.gz", tr)), sep = " ")
  }
  true_by <- sc$theta * a[[sc$driver]] + THETA_HDL * a$HDL_C + THETA_TG * a$TRIGLY + sc$direct * a$LPA
  # Small balanced pleiotropy so the model sits in the random-effects regime as the real one does.
  pleio <- rnorm(nrow(a), 0, .001) * (a$role != "lpa")
  for (o in outcomes) {
    s <- 1 / sqrt(2 * a$f * (1 - a$f) * 2e6)
    fwrite(to_regenie(true_by + pleio + s * rnorm(nrow(a)), s), file.path(dir, sprintf("outcomes_%s.txt.gz", o)), sep = " ")
  }
  # Phenotype matrix: only a correlation matrix is ever taken from this.
  Pm <- matrix(c(1, .9, -.2, .3, .05, .9, 1, -.3, .4, .08, -.2, -.3, 1, -.4, .01, .3, .4, -.4, 1, -.02, .05, .08, .01, -.02, 1), 5, 5, dimnames = list(TRAITS, TRAITS))
  ph <- as.data.table(matrix(rnorm(5000 * 5), 5000, 5) %*% chol(Pm)); setnames(ph, TRAITS)
  ph[, `:=`(FID = .I, IID = .I)]
  fwrite(ph[, c("FID", "IID", TRAITS), with = FALSE], file.path(dir, "pheno.txt"), sep = " ")
  writeLines(CIS_IDS, file.path(dir, "lpa_cis_ids.txt"))
  if (with_externals) {
    # FinnGen: rsID match only, with everything the real file does to you.
    sf <- 1 / sqrt(2 * a$f * (1 - a$f) * 1.2e6)
    bf <- true_by + pleio + sf * rnorm(nrow(a))
    flip <- rbinom(nrow(a), 1, .5) == 1
    ref <- fifelse(flip, a$A1, a$A0); alt <- fifelse(flip, a$A0, a$A1)
    bo <- fifelse(flip, -bf, bf); af <- fifelse(flip, 1 - a$f, a$f)
    # Strand flip on a tenth of non-palindromes.
    st <- !a$pal & rbinom(nrow(a), 1, .1) == 1
    ref[st] <- chartr("ACGT", "TGCA", ref[st]); alt[st] <- chartr("ACGT", "TGCA", alt[st])
    # Palindromes reported on the other strand: labels look aligned, the sign is not.
    # Resolved by frequency when clear, and must be dropped when the EAF is near 0.5.
    pal_opp <- a$pal
    ref[pal_opp] <- a$A0[pal_opp]; alt[pal_opp] <- a$A1[pal_opp]
    bo[pal_opp] <- -bf[pal_opp]; af[pal_opp] <- 1 - a$f[pal_opp]
    fg <- data.table(`#chrom` = a$CHROM, pos = a$GENPOS + 1000L, ref = ref, alt = alt, rsids = a$ID, nearest_genes = "GENE",
                     pval = 2 * pnorm(-abs(bo / sf)), mlogp = log10p(bo / sf), beta = bo, sebeta = sf, af_alt = af, af_alt_cases = af, af_alt_controls = af)
    # Packed rsIDs, and multi-allelic rows that a cartesian merge would double-count.
    fg[seq(10, .N, by = 40), rsids := paste0(rsids, ",rsALIAS", .I)]
    dup <- fg[seq(7, .N, by = 25)]
    dup[, `:=`(alt = fifelse(alt == "A", "C", "A"), ref = fifelse(ref == "G", "T", "G"), beta = 0.5)]
    dup <- dup[!(paste0(alt, ref) %in% c("AT", "TA", "CG", "GC"))]
    fg <- rbind(fg, dup); setorder(fg, `#chrom`, pos)
    fwrite(fg, file.path(dir, "finngen_R11_I9_CHD.gz"), sep = "\t")
    # CARDIoGRAM: lacks the LPA window here, so Lp(a) must come back NA, not a trans-only number.
    sc2 <- 1 / sqrt(2 * a$f * (1 - a$f) * 8e5)
    cg <- data.table(markername = a$ID, chr = a$CHROM, bp_hg19 = a$GENPOS, effect_allele = a$A1, noneffect_allele = a$A0, effect_allele_freq = a$f,
                     median_info = .99, model = "FIXED", beta = true_by + pleio + sc2 * rnorm(nrow(a)), se_dgc = sc2, p_dgc = .5, het_pvalue = .5, n_studies = 48L)
    cg <- cg[!in_lpa_window(chr, bp_hg19)]
    fwrite(cg, file.path(dir, "cad.add.160614.website.txt"), sep = "\t")
  }
  invisible(a)
}

run_main <- function(world, out, extra = character(0), arm = "statinfix", cache = file.path(WORK, "cache"),
                     outcomes = "CAD,MI", sumstat_dir = world) {
  args <- c(MAIN, "--sumstat_dir", sumstat_dir, "--pheno", file.path(world, "pheno.txt"), "--arm", arm,
            "--outcomes", outcomes, "--finngen", file.path(world, "finngen_R11_I9_CHD.gz"),
            "--cardiogram", file.path(world, "cad.add.160614.website.txt"), "--cis_list", file.path(world, "lpa_cis_ids.txt"),
            "--nboot", opt$nboot, "--threads", opt$threads, "--cache_dir", cache, "--outdir", out,
            "--pc_n", "-1", "--pc_cis_n", "-1", extra)
  log_f <- paste0(out, ".log")
  st <- system2("Rscript", args, stdout = log_f, stderr = log_f)
  list(status = as.integer(st), log = log_f, out = out)
}
verdict_of <- function(out, tag = "VERDICT_LITERAL") {
  f <- file.path(out, "VERDICT.txt")
  if (!file.exists(f)) return(NA_character_)
  l <- grep(paste0("^(SENSITIVITY_ARM_)?", tag, ":"), readLines(f), value = TRUE)
  if (!length(l)) return(NA_character_)
  sub(" .*$", "", sub(paste0("^(SENSITIVITY_ARM_)?", tag, ": "), "", l[1]))
}

# ---- 3. the scenarios ------------------------------------------------------------
cat("\n== SURVIVES: planted Lp(a) effect independent of ApoB ==\n")
w1 <- file.path(WORK, "w_SURVIVES"); A1 <- write_world(w1, SCENARIOS$SURVIVES)
r1 <- run_main(w1, file.path(WORK, "o_SURVIVES"))
check("exits 0", r1$status == 0L)
if (r1$status != 0L) { cat(tail(readLines(r1$log), 30), sep = "\n"); stop("cannot continue") }
check("verdict SUPPORTS", identical(verdict_of(r1$out), "SUPPORTS"))
check("positive_control.PASSED written", file.exists(file.path(r1$out, "positive_control.PASSED")))
L1 <- fread(file.path(r1$out, "results_long.csv"))
isets <- fread(file.path(r1$out, "instrument_sets.tsv"))
pub <- isets[variant == "published" & model == "LDL4"]
check("published rule keeps rs10455872 and drops the underflow decoy rs531596603",
      "rs10455872" %in% pub$ID && !"rs531596603" %in% pub$ID)
check("published rule leaves the LPA locus as one variant", sum(pub$in_lpa_window) == 1L)
check("no two published instruments within 10 Mb", pub[order(CHROM, GENPOS), all(diff(GENPOS)[diff(CHROM) == 0] >= 1e7)])
cis <- isets[variant == "cis_dense" & model == "LDL4"]
check("cis-dense retains the full cis list", all(CIS_IDS %in% cis$ID))
check("excluding the window from all selection leaves no window variant", isets[variant == "noLPAwin_all", sum(in_lpa_window)] == 0L)
check("the exclusion also removes the locus-driven decoy just outside the window", !"rs531596603" %in% isets[variant == "noLPAwin_all", ID])
check("LPA window excluded from non-LPA selection keeps Lp(a)'s own lead", "rs10455872" %in% isets[variant == "noLPAwin_nonLPA" & model == "APOB4", ID])
g1 <- function(L, co, v, m, e, o = "CAD") L[cohort == co & outcome == o & variant == v & model == m & exposure == e]
check("LDL4 UKB control near the planted 1.25 / 1.15",
      near(g1(L1, "UKB", "published", "LDL4", "LDL_C")$OR, 1.25, .01) && near(g1(L1, "UKB", "published", "LDL4", "LPA")$OR, 1.15, .01))
check("APOB4 UKB Lp(a) recovers the planted 1.15", near(g1(L1, "UKB", "published", "APOB4", "LPA")$OR, 1.15, .02))
check("FinnGen harmonization through the traps recovers the planted Lp(a) and driver",
      near(g1(L1, "FinnGen_R11", "published", "APOB4", "LPA")$OR, 1.15, .02) && near(g1(L1, "FinnGen_R11", "published", "LDL4", "LDL_C")$OR, 1.25, .02))
check("FinnGen dropped ambiguous palindromes that UKB kept", g1(L1, "FinnGen_R11", "published", "LDL4", "LPA")$n_instruments < g1(L1, "UKB", "published", "LDL4", "LPA")$n_instruments)
check("two- and three-exposure models and all variants are present and fitted on UKB CAD",
      L1[cohort == "UKB" & outcome == "CAD" & exposure == "LPA", .N] == length(MODELS) * length(VARIANTS) &&
        all(L1[cohort == "UKB" & outcome == "CAD" & exposure == "LPA", status] == "ok"))
cg_lpa <- L1[cohort == "CARDIoGRAMplusC4D" & exposure == "LPA" & variant %in% c("published", "cis_dense")]
check("CARDIoGRAM without the LPA window: NA with a named reason, not a trans-only number",
      all(is.na(cg_lpa$OR)) && all(cg_lpa$status == "lpa_locus_absent"))
cg_nw <- L1[cohort == "CARDIoGRAMplusC4D" & exposure == "LPA" & variant == "noLPAwin_all"]
check("CARDIoGRAM noLPAwin_all (window excluded by design) is still estimable where trans instruments exist", any(cg_nw$status == "ok"))
check("condF phenocov differs from gencov0 and both are finite on the primary model",
      with(g1(L1, "UKB", "published", "APOB4", "LPA"), is.finite(condF_phenocov) && is.finite(condF_gencov0) && abs(condF_phenocov - condF_gencov0) > 1e-6))
check("Egger intercept and Q_A present on the primary model",
      with(g1(L1, "UKB", "published", "APOB4", "LPA"), !is.na(egger_intercept_p) && !is.na(Q_A_p_phenocov)))
mp <- g1(L1, "UKB", "published", "APOB4", "LPA"); mc <- g1(L1, "UKB", "cis_dense", "APOB4", "LPA")
check("weighted median: NA with a reason for Lp(a) on one window variant, a number for ApoB beside it",
      is.na(mp$median_b) && grepl("lpa_window_variants_lt", mp$median_note) && !is.na(g1(L1, "UKB", "published", "APOB4", "APOB")$median_b))
check("weighted median: Lp(a) recovers the planted 1.15 on the cis-dense set", near(exp(mc$median_b), 1.15, .03) && mc$median_note == "ok")
check("weighted median runs on CAD only", all(L1[outcome != "CAD", median_note] == "not_fitted"))
AO <- fread(file.path(r1$out, "lpa_attenuation_all_outcomes.csv"))
check("attenuation table covers the other endpoints too", setequal(unique(AO[cohort == "UKB", outcome]), c("CAD", "MI")))
vt <- readLines(file.path(r1$out, "VERDICT.txt"))
check("verdict file reads the sensitivity models and the window-excluded sets",
      all(vapply(c("SENSITIVITY published APOB3:", "SENSITIVITY cis_dense APOB2:", "SENSITIVITY noLPAwin_nonLPA APOB4:", "SENSITIVITY cis_noLPAwin_nonLPA APOB4:", "TRANS_ONLY", "PRIMARY_MODEL_UKB_CAD:"),
                 function(k) any(grepl(k, vt, fixed = TRUE)), logical(1))))
AE <- fread(file.path(r1$out, "lpa_instrument_allele_effects.csv"))
ae <- AE[ID == "rs10455872"]
check("per-allele table: rs10455872 raising allele oriented, Lp(a) about +1 SD, ApoB/Lp(a) near the planted 0.10",
      nrow(ae) == 1L && ae$raising_allele == "G" && near(ae$beta_LPA, 1, .05) && near(ae$apob_sd_per_lpa_sd, .10, .03))
RC <- fread(file.path(r1$out, "readout_compact.csv"))
check("per-allele table includes Lp(a)'s own trans instruments as well as the locus", AE[in_lpa_window == FALSE, .N] > 0L && AE[in_lpa_window == TRUE, .N] >= length(CIS_IDS))
check("paper-pairing verdict (UKB published rule, FinnGen cis-dense) is printed and matches its two cells",
      identical(verdict_of(r1$out, "VERDICT_PAPER_PAIRING"), combine_categories(c(fread(file.path(r1$out, "readout_compact.csv"))[cohort == "UKB" & variant == "published" & model == "APOB4", category],
                                                                                 fread(file.path(r1$out, "readout_compact.csv"))[cohort == "FinnGen_R11" & variant == "cis_dense" & model == "APOB4", category]))))
check("compact readout covers UKB, FinnGen and CARDIoGRAM on CAD", setequal(unique(RC$cohort), c("UKB", "FinnGen_R11", "CARDIoGRAMplusC4D")))
BA <- fread(file.path(r1$out, "bootstrap_attenuation.csv"))
check("bootstrap ran the three flavours with >=1000 resamples", setequal(unique(BA$flavour), c("naive", "stratified", "strat_param")) && all(BA$nboot >= 1000L))
ba1 <- BA[cohort == "UKB" & outcome == "CAD" & variant == "published" & model == "APOB4"]
check("bootstrap naive flavour exposes the lost lead variant (~37% of resamples)", near(ba1[flavour == "naive", frac_no_window_variant], exp(-1), .05))
check("bootstrap stratified flavours always contain the LPA locus", all(ba1[flavour != "naive", frac_no_window_variant] == 0))
ba1p <- ba1[flavour == "strat_param"]
# Not "covers exactly 0": ApoB is an imperfect stand-in for the LDL-C driver here, so the
# ApoB coefficient differs and the true attenuation is slightly off zero by construction.
check("attenuation CI present, narrow and nowhere near the 25% boundary",
      !is.na(ba1p$atten_pct_lo95) && ba1p$atten_pct_lo95 > -10 && ba1p$atten_pct_hi95 < 10)

cat("\n== published implementations run on the same files must agree ==\n")
cad_only <- file.path(WORK, "o_pub_run_mvmr.csv")
st <- system2("Rscript", c(file.path(R_DIR, "run_mvmr.R"), "--outcome", "CAD", "--sumstat_dir", w1, "--out", cad_only),
              stdout = paste0(cad_only, ".log"), stderr = paste0(cad_only, ".log"))
check("run_mvmr.R ran", as.integer(st) == 0L)
PUB <- fread(cad_only); PIV <- fread(sub("\\.csv$", ".instruments.tsv", cad_only))
mine <- L1[cohort == "UKB" & outcome == "CAD" & variant == "published" & model == "LDL4"]
check("published variant: same instrument set as run_mvmr.R", setequal(PIV$ID, build_ids <- isets[variant == "published" & model == "LDL4", ID]))
check("published variant: same estimates as run_mvmr.R (1e-8)", all(abs(mine[match(PUB$exposure, exposure), b] - PUB$mvivw_b) < 1e-8))
check("published variant: same SE as run_mvmr.R (1e-8)", all(abs(mine[match(PUB$exposure, exposure), se] - PUB$mvivw_se) < 1e-8))
check("published variant: same Egger intercept as run_mvmr.R", near(mine$egger_intercept[1], PUB$mvegger_intercept[1], 1e-8))
cov_out <- file.path(WORK, "o_pub_cov.csv")
st <- system2("Rscript", c(file.path(R_DIR, "mvmr_covariance.R"), "--sumstat_dir", w1, "--pheno", file.path(w1, "pheno.txt"),
                           "--outcomes", "CAD", "--nboot", "5", "--cis_list", file.path(w1, "lpa_cis_ids.txt"), "--out", cov_out),
              stdout = paste0(cov_out, ".log"), stderr = paste0(cov_out, ".log"))
check("mvmr_covariance.R ran", as.integer(st) == 0L)
COV <- fread(cov_out)[outcome == "CAD"]
micd <- L1[cohort == "UKB" & outcome == "CAD" & variant == "cis_dense" & model == "LDL4"]
check("cis_dense variant: same instrument count and estimates as mvmr_covariance.R",
      micd$n_instruments[1] == COV$n_iv[1] && all(abs(micd[match(COV$exposure, exposure), b] - COV$b) < 1e-8))
check("cis_dense variant: same phenocov conditional F as mvmr_covariance.R",
      all(abs(micd[match(COV$exposure, exposure), condF_phenocov] - COV$condF_phenocov) < 1e-6))
ext_out <- file.path(WORK, "o_pub_ext.csv")
st <- system2("Rscript", c(file.path(R_DIR, "run_external_mr.R"), "--sumstat_dir", w1, "--outcome_file", file.path(w1, "finngen_R11_I9_CHD.gz"),
                           "--outcome_format", "finngen", "--outcome_label", "FG", "--exposures", "LDL_C,HDL_C,TRIGLY,LPA",
                           "--mode", "mvmr", "--out", ext_out),
              stdout = paste0(ext_out, ".log"), stderr = paste0(ext_out, ".log"))
check("run_external_mr.R ran", as.integer(st) == 0L)
EXT <- fread(ext_out)
mfg <- L1[cohort == "FinnGen_R11" & outcome == "CAD" & variant == "published" & model == "LDL4"]
check("FinnGen harmonization: same instrument count and OR as run_external_mr.R on the distance rule",
      mfg$n_instruments[1] == EXT$n_instruments[1] && all(abs(mfg[match(EXT$exposure, exposure), OR] - EXT$OR_perSD) < 1e-8))

cat("\n== APOB_ONLY: Lp(a) acts only through its ApoB content ==\n")
w2 <- file.path(WORK, "w_APOB_ONLY"); write_world(w2, SCENARIOS$APOB_ONLY)
r2 <- run_main(w2, file.path(WORK, "o_APOB_ONLY"))
check("exits 0", r2$status == 0L)
if (r2$status != 0L) { cat(tail(readLines(r2$log), 30), sep = "\n"); stop("cannot continue") }
L2 <- fread(file.path(r2$out, "results_long.csv"))
check("control still passes: the LDL-C model alone cannot tell the scenarios apart",
      near(g1(L2, "UKB", "published", "LDL4", "LPA")$OR, 1.15, .01) && near(g1(L2, "UKB", "published", "LDL4", "LDL_C")$OR, 1.25, .01))
check("verdict WEAKENED", identical(verdict_of(r2$out), "WEAKENED"))
check("THE TEST: verdict differs from the scenario where the effect survives", !identical(verdict_of(r1$out), verdict_of(r2$out)))
check("ApoB model Lp(a) OR is near null while the LDL-C model's is 1.15", near(g1(L2, "UKB", "published", "APOB4", "LPA")$OR, 1, .02))
RC2 <- fread(file.path(r2$out, "readout_compact.csv"))
check("cis-dense reaches the same verdict", identical(verdict_of(r2$out, "VERDICT_CIS_DENSE"), "WEAKENED"))
check("two-exposure ApoB model also near null", near(g1(L2, "UKB", "published", "APOB2", "LPA")$OR, 1, .03))
AE2 <- fread(file.path(r2$out, "lpa_instrument_allele_effects.csv"))[ID == "rs10455872"]
check("per-allele table reports the planted ApoB/Lp(a) 0.70", near(AE2$apob_sd_per_lpa_sd, .70, .03))
ba2 <- fread(file.path(r2$out, "bootstrap_attenuation.csv"))[cohort == "UKB" & outcome == "CAD" & variant == "published" & model == "APOB4" & flavour == "strat_param"]
check("matched attenuation near 100%, CI excludes the 50% boundary", near(ba2$atten_pct_obs, 100, 15) && ba2$atten_pct_lo95 > 50)

cat("\n== PARTIAL: 35% of the LDL-model estimate is ApoB content ==\n")
w3 <- file.path(WORK, "w_PARTIAL"); write_world(w3, SCENARIOS$PARTIAL)
r3 <- run_main(w3, file.path(WORK, "o_PARTIAL"))
check("exits 0", r3$status == 0L)
if (r3$status != 0L) { cat(tail(readLines(r3$log), 30), sep = "\n"); stop("cannot continue") }
check("verdict ATTENUATED", identical(verdict_of(r3$out), "ATTENUATED"))
ba3 <- fread(file.path(r3$out, "bootstrap_attenuation.csv"))[cohort == "UKB" & outcome == "CAD" & variant == "published" & model == "APOB4" & flavour == "strat_param"]

cat("\n== BAD_CONTROL: positive control must stop the run ==\n")
w4 <- file.path(WORK, "w_BAD_CONTROL"); write_world(w4, SCENARIOS$BAD_CONTROL)
r4 <- run_main(w4, file.path(WORK, "o_BAD_CONTROL"))
check("exits with status 2", r4$status == 2L)
check("verdict INVALID", identical(verdict_of(r4$out), "INVALID"))
check("no estimate files written", !any(file.exists(file.path(r4$out, c("results_long.csv", "readout_compact.csv", "bootstrap_attenuation.csv", "positive_control.PASSED")))))
check("no ApoB model appears in anything written", !any(grepl("APOB4|APOB3|APOB2", unlist(lapply(setdiff(list.files(r4$out, full.names = TRUE), file.path(r4$out, "instrument_sets.tsv")), readLines)))))
# Run with the comparison options on: nothing else exercises them, and they only report.
r4b <- run_main(w4, file.path(WORK, "o_BAD_CONTROL_lead"), extra = c("--pc_ldl_or", "1.40", "--pc_lpa_or", "1.30", "--pc_n", "161",
                                                                     "--pc_instruments", sub("\\.csv$", ".instruments.tsv", cad_only)))
check("the same data pass when the control values are the planted ones, so the stop was the control and not a crash", r4b$status == 0L)
v4b <- readLines(file.path(r4b$out, "VERDICT.txt"))
check("an instrument count unlike the published one is reported, not enforced", any(grepl("the instrument count is [0-9]+, published 161", v4b)))
check("the published instrument file is compared and the overlap reported", any(grepl("published instrument file .*: [0-9]+ in common, [0-9]+ only published, [0-9]+ only rebuilt", v4b)))

stale <- file.path(WORK, "o_stale"); dir.create(stale)
for (f in c("results_long.csv", "readout_compact.csv", "bootstrap_attenuation.csv", "positive_control.PASSED", "VERDICT.txt"))
  file.copy(file.path(r1$out, f), file.path(stale, f))
r4s <- run_main(w4, stale)
check("a failed control removes the previous run's estimates from the same directory",
      r4s$status == 2L && identical(verdict_of(stale), "INVALID") &&
        !any(file.exists(file.path(stale, c("results_long.csv", "readout_compact.csv", "bootstrap_attenuation.csv", "positive_control.PASSED")))))

cat("\n== FinnGen without the LPA locus: the read must not fall back on UKB alone ==\n")
nf <- file.path(WORK, "w_noFGlocus"); dir.create(nf)
for (f in setdiff(list.files(w1), "finngen_R11_I9_CHD.gz")) file.symlink(file.path(w1, f), file.path(nf, f))
fg0 <- fread(file.path(w1, "finngen_R11_I9_CHD.gz"))
fwrite(fg0[!sub(",.*$", "", rsids) %in% c(CIS_IDS, "rs531596603")], file.path(nf, "finngen_R11_I9_CHD.gz"), sep = "\t")
rnf <- run_main(nf, file.path(WORK, "o_noFGlocus"), extra = c("--pc_instruments", file.path(WORK, "no_such_file.tsv")))
check("exits 0 with verdict NOT_READABLE, not SUPPORTS", rnf$status == 0L && identical(verdict_of(rnf$out), "NOT_READABLE"))
check("a missing published instrument file is reported as a comparison not made", any(grepl("was NOT compared with the published one", readLines(file.path(rnf$out, "VERDICT.txt")))))

cat("\n== sensitivity arm: statin-naive layout, fallback files, stamp gate ==\n")
nv <- file.path(WORK, "w_naive_dir"); dir.create(nv)
for (f in c("lipids_LDL_C", "lipids_APOB", "lipids_TRIGLY", "outcomes_CAD", "outcomes_MI"))
  file.symlink(file.path(w1, paste0(f, ".txt.gz")), file.path(nv, paste0(f, ".txt.gz")))
rn0 <- run_main(w1, file.path(WORK, "o_naive_nostamp"), sumstat_dir = nv, arm = "statinnaive",
                extra = c("--arm_type", "sensitivity", "--fallback_dir", w1))
check("refused without the primary arm's control stamp", rn0$status != 0L)
rn <- run_main(w1, file.path(WORK, "o_naive"), sumstat_dir = nv, arm = "statinnaive",
               extra = c("--arm_type", "sensitivity", "--fallback_dir", w1, "--control_stamp", file.path(r1$out, "positive_control.PASSED")))
check("runs with the stamp", rn$status == 0L)
if (rn$status != 0L) cat(tail(readLines(rn$log), 30), sep = "\n")
vn <- readLines(file.path(rn$out, "VERDICT.txt"))
check("LPA and HDL_C taken from the fallback directory are flagged in the verdict file",
      any(grepl("exposure LPA is NOT from this arm", vn)) && any(grepl("exposure HDL_C is NOT from this arm", vn)))
check("sensitivity verdict line is labelled as such", any(grepl("^SENSITIVITY_ARM_VERDICT_LITERAL:", vn)))
check("APOB symlink carried into the note", any(grepl("symlink ->", vn)))

# Censored Lp(a) GWAS: a directory holding ONLY lipids_LPA, as merged_lpacens will. Planted
# here as half the per-allele effect, so the file actually used can be read off the result.
cz <- file.path(WORK, "w_lpacens_dir"); dir.create(cz)
lp0 <- fread(file.path(w1, "lipids_LPA.txt.gz"))
lp0[, `:=`(BETA = BETA / 2, CHISQ = (BETA / 2 / SE)^2)]
lp0[, LOG10P := -(pnorm(-abs(BETA / SE), log.p = TRUE) + log(2)) / log(10)]
fwrite(lp0, file.path(cz, "lipids_LPA.txt.gz"), sep = " ")
rc <- run_main(w1, file.path(WORK, "o_lpacens"), sumstat_dir = cz, arm = "lpacens",
               extra = c("--arm_type", "sensitivity", "--fallback_dir", w1, "--control_stamp", file.path(r1$out, "positive_control.PASSED")))
check("censored-Lp(a) layout runs through --sumstat_dir with no code path of its own", rc$status == 0L)
if (rc$status != 0L) cat(tail(readLines(rc$log), 30), sep = "\n")
aec <- fread(file.path(rc$out, "lpa_instrument_allele_effects.csv"))[ID == "rs10455872"]
check("the censored file is the one used (planted per-allele Lp(a) 0.5 SD), and other exposures are flagged as fallback",
      near(aec$beta_LPA, .5, .03) && any(grepl("exposure APOB is NOT from this arm", readLines(file.path(rc$out, "VERDICT.txt")))) &&
        !any(grepl("exposure LPA is NOT from this arm", readLines(file.path(rc$out, "VERDICT.txt")))))
wstamp <- file.path(WORK, "stamp_failed.txt"); writeLines("FAILED nothing", wstamp)
rw_ <- run_main(w1, file.path(WORK, "o_badstamp"), sumstat_dir = cz, arm = "lpacens",
                extra = c("--arm_type", "sensitivity", "--fallback_dir", w1, "--control_stamp", wstamp))
check("a stamp that does not record a pass is refused", rw_$status != 0L)

cat("\n== failure modes ==\n")
bad <- file.path(WORK, "w_bad_gz"); dir.create(bad)
for (f in list.files(w1)) file.symlink(file.path(w1, f), file.path(bad, f))
unlink(file.path(bad, "lipids_APOB.txt.gz"))
raw <- readBin(file.path(w1, "lipids_APOB.txt.gz"), "raw", file.size(file.path(w1, "lipids_APOB.txt.gz")))
writeBin(raw[seq_len(length(raw) %/% 2L)], file.path(bad, "lipids_APOB.txt.gz"))
rb <- run_main(bad, file.path(WORK, "o_bad_gz"), cache = file.path(WORK, "cache_bad"))
check("a truncated .gz stops the run", rb$status != 0L && !file.exists(file.path(rb$out, "results_long.csv")))
miss <- file.path(WORK, "w_missing"); dir.create(miss)
for (f in setdiff(list.files(w1), "lipids_LPA.txt.gz")) file.symlink(file.path(w1, f), file.path(miss, f))
rm_ <- run_main(miss, file.path(WORK, "o_missing"))
check("a missing exposure with no fallback stops the run", rm_$status != 0L)
rs <- system2("Rscript", c(MAIN, "--sumstat_dir", w1, "--pheno", file.path(w1, "pheno.txt"), "--finngen", file.path(w1, "finngen_R11_I9_CHD.gz"),
                           "--cis_list", file.path(w1, "lpa_cis_ids.txt"), "--nboot", "100", "--cache_dir", file.path(WORK, "cache"), "--outdir", file.path(WORK, "o_smallboot")),
              stdout = FALSE, stderr = FALSE)
check("--nboot below 1000 is refused", as.integer(rs) != 0L)
nophen <- file.path(WORK, "pheno_noAPOB.txt")
fwrite(fread(file.path(w1, "pheno.txt"))[, !"APOB"], nophen, sep = " ")
rp <- system2("Rscript", c(MAIN, "--sumstat_dir", w1, "--pheno", nophen, "--finngen", file.path(w1, "finngen_R11_I9_CHD.gz"),
                           "--cis_list", file.path(w1, "lpa_cis_ids.txt"), "--nboot", "1000", "--cache_dir", file.path(WORK, "cache"), "--outdir", file.path(WORK, "o_nophen")),
              stdout = FALSE, stderr = FALSE)
check("a phenotype matrix lacking an exposure column is refused", as.integer(rp) != 0L)
tmpf <- tempfile(); writeLines("x", tmpf); k1 <- cache_file(WORK, "t", "a", 1); k2 <- cache_file(WORK, "t", "a", 2)
check("cache key changes with input identity", k1 != k2)

# ---- 4. does the attenuation interval cover? --------------------------------------
# In memory only. The truth is the noise-free WLS on the true effects under the same
# weights, per scenario. A percentile interval that looks tidy but undercovers would make
# the readout's band comparison mean nothing. Run at two scales: the precise outcome the
# worlds above use, and the real one, where rs10455872 has SE(log OR) about 0.015 and
# Q_A says sigma is about 2, so the denominator of the percentage is genuinely noisy.
cat("\n== coverage of the attenuation interval ==\n")
coverage_run <- function(sc, nrep, nboot, single_lpa, n_out, overdisp, truth_only = FALSE) {
  a <- copy(TRUE_FX$poly)[role == "lead"]
  lp <- copy(TRUE_FX$lpa)[ID %in% (if (single_lpa) "rs10455872" else CIS_IDS)]
  lp[, `:=`(LPA = g, APOB = sc$c_apob * g, LDL_C = sc$c_ldl * g, HDL_C = 0, TRIGLY = 0)]
  G <- as.matrix(rbind(a[, ..TRAITS], lp[, ..TRAITS]))
  n <- nrow(G); win <- c(rep(FALSE, nrow(a)), rep(TRUE, nrow(lp)))
  f <- c(rep(.3, nrow(a)), rep(.07, nrow(lp)))
  sx <- 1 / sqrt(2 * f * (1 - f) * N_GWAS); sy <- 1 / sqrt(2 * f * (1 - f) * n_out)
  by_true <- sc$theta * G[, sc$driver] + THETA_HDL * G[, "HDL_C"] + THETA_TG * G[, "TRIGLY"] + sc$direct * G[, "LPA"]
  tr_ref <- c("LDL_C", "HDL_C", "TRIGLY", "LPA"); tr_alt <- c("APOB", "HDL_C", "TRIGLY", "LPA")
  wfit <- function(tr) .lm.fit(G[, tr] / sy, by_true / sy)$coefficients[4]
  truth <- 100 * (wfit(tr_ref) - wfit(tr_alt)) / wfit(tr_ref)
  if (truth_only) return(unname(truth))
  Pm <- diag(5); dimnames(Pm) <- list(TRAITS, TRAITS)
  Pm["LDL_C", "APOB"] <- Pm["APOB", "LDL_C"] <- .9
  U <- chol(Pm)
  res <- matrix(NA_real_, nrep, 3, dimnames = list(NULL, c("naive", "stratified", "strat_param")))
  # One seed per replicate and per bootstrap, all drawn from the fixed seed 42.
  set.seed(42); rep_seeds <- sample.int(1e8, 2L * nrep)
  est <- wid <- numeric(nrep)
  for (i in seq_len(nrep)) {
    set.seed(rep_seeds[i])
    BX <- G + sx * (matrix(rnorm(n * 5), n, 5) %*% U)
    # overdisp multiplies every variant's SE: the random-effects model the IVW SE assumes.
    D <- list(BX = BX, BXSE = matrix(sx, n, 5, dimnames = dimnames(G)), BY = by_true + sy * overdisp * rnorm(n), BYSE = sy, in_window = win,
              ids = paste0("r", 1:n), LOG10P = matrix(50, n, 5, dimnames = dimnames(G)))
    b <- boot_attenuation(D, tr_ref, tr_alt, Pm, nboot, seed = rep_seeds[nrep + i])
    res[i, ] <- vapply(colnames(res), function(fl) {
      x <- b[flavour == fl]
      if (is.na(x$atten_pct_lo95)) NA_real_ else as.numeric(x$atten_pct_lo95 <= truth && x$atten_pct_hi95 >= truth)
    }, numeric(1))
    x <- b[flavour == "strat_param"]; est[i] <- x$atten_pct_obs; wid[i] <- x$atten_pct_hi95 - x$atten_pct_lo95
  }
  list(truth = truth, cover = colMeans(res, na.rm = TRUE), defined = colMeans(!is.na(res)), sd_est = sd(est), se_boot = mean(wid, na.rm = TRUE) / 3.92)
}
# The PARTIAL world's estimand is the noise-free WLS, not the 35 used to build it: 35% is of
# log(1.15), and the WLS weights move it slightly. Its "published" arm has one LPA variant.
partial_truth <- coverage_run(SCENARIOS$PARTIAL, 0L, 0L, TRUE, 2e6, 1.3, truth_only = TRUE)
cat(sprintf("  PARTIAL world: truth %.2f%%, CI %.2f-%.2f\n", partial_truth, ba3$atten_pct_lo95, ba3$atten_pct_hi95))
check("matched attenuation CI covers the PARTIAL world's noise-free truth", ba3$atten_pct_lo95 < partial_truth && ba3$atten_pct_hi95 > partial_truth)

for (cs in list(list(lab = "one LPA variant, precise outcome", single = TRUE, n_out = 2e6, od = 1.3),
                list(lab = "cis-dense LPA,   precise outcome", single = FALSE, n_out = 2e6, od = 1.3),
                list(lab = "one LPA variant, real-scale outcome", single = TRUE, n_out = 140000, od = 2),
                list(lab = "cis-dense LPA,   real-scale outcome", single = FALSE, n_out = 140000, od = 2))) {
  cv <- coverage_run(SCENARIOS$PARTIAL, opt$cov_reps, 300L, cs$single, cs$n_out, cs$od)
  cat(sprintf("  %s, truth %.1f%%: coverage naive %.3f (defined %.2f), stratified %.3f (%.2f), strat_param %.3f (%.2f); empirical SD %.2f vs bootstrap SE %.2f\n",
              cs$lab, cv$truth, cv$cover["naive"], cv$defined["naive"], cv$cover["stratified"], cv$defined["stratified"],
              cv$cover["strat_param"], cv$defined["strat_param"], cv$sd_est, cv$se_boot))
  check(sprintf("strat_param interval covers 90-99%%: %s", cs$lab),
        cv$cover["strat_param"] >= .90 && cv$cover["strat_param"] <= .99 && cv$defined["strat_param"] > .99)
  check(sprintf("bootstrap SE within 25%% of the empirical SD: %s", cs$lab), abs(cv$se_boot / cv$sd_est - 1) < .25)
}

cat(sprintf("\n==== selftest: %d passed, %d failed ====\n", N_PASS, N_FAIL))
if (N_FAIL > 0L) quit(save = "no", status = 1L)
cat("VERDICTS: SURVIVES=", verdict_of(r1$out), " APOB_ONLY=", verdict_of(r2$out), " PARTIAL=", verdict_of(r3$out), " BAD_CONTROL=", verdict_of(r4$out), "\n", sep = "")
