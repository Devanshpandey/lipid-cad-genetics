#!/usr/bin/env Rscript
# ============================================================================
# POST HOC: the Lp(a) by lipid-region product with the LPA locus removed from
# the lipid-region score. Declared in PREREGISTRATION.md, "Post hoc analysis",
# after the pre-registered run (job 3602201) had been read. It cannot change
# the verdict of that run.
#
# lpa_susceptibility.R is not sourced (its main block is unguarded); its
# function definitions and named constants are evaluated from the parsed file.
#
#   --step selftest   synthetic data only
#   --step weights    write the LPA-region weights for PLINK2 (no outcome data)
#   --step run        gates (must pass), then the three models
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(optparse) })
args0 <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])))
REPO <- normalizePath(file.path(here, "..", "..", ".."))
A6 <- file.path(REPO, "scripts", "agent6_resilience")
for (f in c("lib_resilience.R", "lib_genetics.R", "lib_synergy.R")) source(file.path(A6, "R", f))
LH <- new.env(); sys.source(file.path(REPO, "scripts", "agent1_genetics", "lpa_horizon", "lpa_horizon.R"), envir = LH)
take <- function(file, consts) for (e in parse(file)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      (as.character(e[[2]]) %in% consts || (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))))) eval(e, globalenv())
}
take(file.path(here, "lpa_susceptibility.R"), c("WINDOW", "WINDOWS", "HALF", "SPLINE_DF", "REF", "QUAD"))
stopifnot(all(vapply(c("quad_full", "summ", "load_eligible", "attach_exposures", "excess_synergy", "QUAD", "REF"), exists, logical(1))))

LPA_GENE <- c(160952515, 161087407)     # GRCh37, chromosome 6
RS10455872 <- 161010118
QUAD3 <- c("M", "L", "S", "N", "I(M^2)", "I(L^2)", "I(S^2)", "I(N^2)", "M:L", "M:S", "M:N", "L:N")

# weights of the lipid-region component that lie in the one lipid region overlapping the LPA gene
lpa_region_weights <- function(w, cw, regions, window = WINDOW) {
  reg <- regions[window_kb == window & CHR == 6 & start <= LPA_GENE[2] & end >= LPA_GENE[1]]
  if (nrow(reg) != 1L) stop("expected exactly one lipid region overlapping the LPA gene, found ", nrow(reg))
  if (!(reg$start <= RS10455872 && RS10455872 <= reg$end)) stop("the region does not contain rs10455872")
  inside <- w[CHR == 6 & BP >= reg$start & BP <= reg$end, SNP]
  lc <- paste0("L", window); nc <- paste0("N", window)
  out <- cw[SNP %in% inside & get(lc) != 0, .(SNP, A1, S_LPA = get(lc))]
  if (cw[SNP %in% out$SNP, any(get(nc) != 0)]) stop("a variant of the region carries a non-lipid weight")
  list(weights = out, region = reg, n_inside = length(inside))
}
# centred components in units of the SD of the total score; L = Lr + S exactly
components <- function(d, window = WINDOW) {
  x <- copy(d); sdT <- sd(x$TOTAL); l <- x[[paste0("L", window)]]; n <- x[[paste0("N", window)]]; s <- x$S_LPA
  x[, `:=`(L_orig = (l - mean(l)) / sdT, S = (s - mean(s)) / sdT, L_rest = ((l - s) - mean(l - s)) / sdT, N = (n - mean(n)) / sdT)]
  x
}
fit_models <- function(d, covars, col, label) {
  x <- droplevels(d[!is.na(d[[col]])]); x[, M := get(col)]
  one <- function(lcol, terms, model) {
    y <- copy(x); y[, L := get(lcol)]
    q <- quad_full(y, covars, NULL, terms, model); r <- summ(q, label)
    r[, model := model]
    if ("M:S" %in% names(q$b)) r[, `:=`(beta_MS = q$b[["M:S"]], se_MS = sqrt(q$V["M:S", "M:S"]), p_MS = 2 * pnorm(-abs(q$b[["M:S"]] / sqrt(q$V["M:S", "M:S"]))))]
    r[, p_ML := 2 * pnorm(-abs(beta_ML / se_ML))]
    r
  }
  rbindlist(list(one("L_orig", QUAD, "1 original lipid-region score"), one("L_rest", QUAD, "2 LPA region removed"),
                 one("L_rest", QUAD3, "3 LPA region as its own term")), fill = TRUE)
}
describe <- function(d) {
  data.table(n = nrow(d), cor_M_L_original = cor(d$M_lpa, d$L_orig, use = "complete.obs"), cor_M_LPA_region = cor(d$M_lpa, d$S, use = "complete.obs"),
             cor_M_L_rest = cor(d$M_lpa, d$L_rest, use = "complete.obs"), cor_LPA_region_L_original = cor(d$S, d$L_orig),
             lpa_region_share_of_L_variance = var(d$S) / var(d$L_orig), max_identity_error = max(abs(d$L_orig - d$L_rest - d$S)))
}
run_posthoc <- function(dall, covars, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  d <- components(dall)
  res <- rbind(fit_models(d, covars, "M_lpa", "predicted Lp(a), all eligible"),
               fit_models(d[statin == 0], covars, "M_lpa", "predicted Lp(a), statin-free"),
               fit_models(d, covars, "M_lpa_measured", "measured Lp(a), assay limits imposed"), fill = TRUE)
  ds <- describe(d)
  LH$wcsv(res, file.path(outdir, "product_terms_lpa_removed.csv")); LH$wcsv(ds, file.path(outdir, "components.csv"))
  list(res = res, describe = ds, d = d)
}

selftest <- function(workdir) {
  n_ok <- 0L
  expect <- function(cond, msg) { if (!isTRUE(cond)) stop("SELFTEST FAILED: ", msg, call. = FALSE); n_ok <<- n_ok + 1L }
  fails <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

  # weights builder on a fixture
  w <- data.table(CHR = c(6, 6, 6, 6, 1), SNP = paste0("rs", 1:5), BP = c(159991741, 161010118, 162021039, 162021040, 161000000), A1 = "A", A2 = "G", BETA = 0.01)
  cw <- data.table(SNP = paste0("rs", 1:5), A1 = "A", L500 = c(0.01, 0.02, 0, 0, 0.03), N500 = c(0, 0, 0, 0.04, 0))
  reg <- data.table(window_kb = c(500, 500, 250), CHR = c(6, 6, 6), start = c(159991741, 1000, 160500000), end = c(162021039, 2000, 161500000))
  lw <- lpa_region_weights(w, cw, reg)
  expect(identical(lw$weights$SNP, c("rs1", "rs2")) && all(lw$weights$S_LPA == c(0.01, 0.02)) && lw$n_inside == 3L, "only lipid-weighted variants inside the LPA region are kept, ends inclusive")
  expect(fails(lpa_region_weights(w, cw, reg[window_kb == 250])) && fails(lpa_region_weights(w, cw, rbind(reg, reg[1]))), "no region or two regions is an error")
  expect(fails(lpa_region_weights(w, copy(cw)[SNP == "rs2", N500 := 0.5], reg)), "a variant with both weights is an error")

  # same locus on both sides: a product planted with the LPA-region component only
  set.seed(51); n <- 400000
  lpa <- pmax(0, rgamma(n, shape = 1.2, scale = 14) + 110 * rbinom(n, 2, 0.08)); M <- zs(lpa)
  d <- data.table(eid = seq_len(n), M_lpa = M, M_lpa_measured = M + rnorm(n, 0, 0.3), statin = rbinom(n, 1, 0.12), age = runif(n, 40, 70),
                  sex = factor(rbinom(n, 1, 0.5)), genotype_batch = factor(sample(1:3, n, TRUE)))
  for (k in 1:10) set(d, j = paste0("PC", k), value = rnorm(n))
  S <- 0.15 * M + rnorm(n, 0, 0.15); Lo <- rnorm(n, 0, 0.30); N <- rnorm(n, 0, 0.94)
  d[, `:=`(S_LPA = S + 0.2, L500 = Lo + S + 0.2, N500 = N)]; d[, TOTAL := L500 + N500]
  sdT <- sd(d$TOTAL)
  lin <- 0.15 * M + 0.55 * Lo / sdT + 0.55 * N / sdT + 0.40 * S / sdT - 0.80 * M * S / sdT + 0.04 * (d$age - 55)
  t <- rexp(n, 0.004 * exp(lin)); d[, `:=`(event = as.integer(t <= 11), time_analysis = pmin(t, 11))]
  cv <- c("age", "sex", paste0("PC", 1:10), "genotype_batch")
  out <- run_posthoc(d, cv, file.path(workdir, "posthoc_smoke"))
  r <- out$res[analysis == "predicted Lp(a), all eligible"]; m1 <- r[model %like% "^1"]; m2 <- r[model %like% "^2"]; m3 <- r[model %like% "^3"]
  message(sprintf("[posthoc sim] Lp(a) x lipid-region: original %.3f (SE %.3f), LPA region removed %.3f (SE %.3f); Lp(a) x LPA region %.3f (SE %.3f), planted -0.80",
                  m1$beta_ML, m1$se_ML, m2$beta_ML, m2$se_ML, m3$beta_MS, m3$se_MS))
  expect(out$describe$max_identity_error < 1e-9 && abs(out$describe$cor_M_L_rest) < 0.01 && out$describe$cor_M_LPA_region > 0.5, "the two parts sum to the original score and only the LPA part tracks Lp(a)")
  expect(m1$beta_ML / m1$se_ML < -2.5, "with the locus inside the score the product is negative (the artefact exists)")
  expect(abs(m2$beta_ML) < 3.5 * m2$se_ML, "with the locus removed the product is not different from zero")
  expect(abs(m3$beta_MS + 0.80) < 3.5 * m3$se_MS && abs(m3$beta_ML) < 3.5 * m3$se_ML, "entered on its own, the LPA-region component carries the planted product")
  expect(nrow(out$res) == 9 && all(file.exists(file.path(workdir, "posthoc_smoke", c("product_terms_lpa_removed.csv", "components.csv")))), "every model and output is present")
  cat(sprintf("\nPOSTHOC SELFTEST PASSED: %d checks\n", n_ok))
}

opt <- optparse::parse_args(OptionParser(option_list = list(
  make_option("--step", type = "character", default = "selftest"),
  make_option("--weights_all", type = "character"), make_option("--component_weights", type = "character"),
  make_option("--regions", type = "character"), make_option("--lpa_sscore", type = "character"),
  make_option("--reference", type = "character"),
  make_option("--run_dir", type = "character"), make_option("--components", type = "character"),
  make_option("--raw_dir", type = "character"), make_option("--covar", type = "character"),
  make_option("--dose_lpa", type = "character"), make_option("--w_lpa", type = "character"),
  make_option("--dose_ldl", type = "character"), make_option("--w_ldl", type = "character"),
  make_option("--prs", type = "character"), make_option("--outdir", type = "character"),
  make_option("--workdir", type = "character", default = tempdir()),
  make_option("--censor", type = "character", default = "2022-10-31"),
  make_option("--threads", type = "integer", default = 8L))))
setDTthreads(opt$threads)
if (opt$step == "selftest") { selftest(opt$workdir); quit(save = "no", status = 0L) }
if (opt$step == "weights") {
  w <- fread(opt$weights_all, header = FALSE, col.names = c("CHR", "SNP", "BP", "A1", "A2", "BETA"))
  lw <- lpa_region_weights(w, fread(opt$component_weights), fread(opt$regions))
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(lw$weights, file.path(opt$outdir, "lpa_region_weights.tsv"), sep = "\t"); writeLines(lw$weights$SNP, file.path(opt$outdir, "lpa_region.snps"))
  cat(sprintf("LPA lipid region chr6:%d-%d; %d score variants inside, %d with a lipid-region weight\n", lw$region$start, lw$region$end, lw$n_inside, nrow(lw$weights)))
  quit(save = "no", status = 0L)
}
if (opt$step != "run") stop("unknown --step")
cfg <- load_config(file.path(A6, "config.yaml")); sy <- yaml::read_yaml(file.path(A6, "config_synergy.yaml"))
CENS <- as.integer(as.IDate(opt$censor))
R <- LH$load_real(opt, CENS)
if (!LH$gates(R, CENS, file.path(opt$outdir, "gates_recheck"))) { cat("\nGATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
covars_all <- c(sy$covariates, "genotype_batch")
dall <- attach_exposures(load_eligible(opt$run_dir, opt$components, cfg, sy), R$P)
s <- fread(opt$lpa_sscore); setnames(s, 1, "FID")
if (!("S_LPA_SUM" %in% names(s))) stop("the PLINK2 output has no S_LPA_SUM column")
dall[, S_LPA := s$S_LPA_SUM[match(eid, as.integer(s$IID))]]
have <- mean(!is.na(dall$S_LPA)); n_all <- nrow(dall); ev_all <- sum(dall$event)
dall <- dall[!is.na(S_LPA)]
d0 <- components(dall); d0[, `:=`(M = M_lpa, L = L_orig)]
orig <- summ(quad_full(d0[!is.na(M)], covars_all), "check")
ref <- fread(opt$reference)[analysis == "predicted Lp(a), all eligible (PRIMARY)"]
g <- data.table(
  gate = c("eligible cohort", "its events", "LPA-region score present, %", "parts sum to the original lipid-region score (max error)",
           "pre-registered Lp(a) x lipid-region product reproduced", "pre-registered Lp(a) x non-lipid product reproduced"),
  value = c(n_all, ev_all, 100 * have, max(abs(d0$L_orig - d0$L_rest - d0$S)), orig$beta_ML, orig$beta_MN),
  target = c(REF$n_all, REF$ev_all, 99.9, 1e-9, ref$beta_ML, ref$beta_MN),
  rule = c("equal", "equal", "at least", "below", "within 1e-6", "within 1e-6"))
g[, pass := fifelse(rule == "equal", value == target, fifelse(rule == "at least", value >= target, fifelse(rule == "below", value < target, abs(value - target) < 1e-6)))]
options(width = 250); print(g, digits = 10)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE); fwrite(g, file.path(opt$outdir, "gates_posthoc.csv"))
if (!isTRUE(all(g$pass))) { cat("\nPOST HOC GATES FAILED: nothing is to be read.\n"); quit(save = "no", status = 2L) }
cat("\nGATES PASSED\n")
out <- run_posthoc(dall, covars_all, opt$outdir)
cat("\n=== components of the lipid-region score (no outcome) ===\n"); print(out$describe, digits = 3)
cat("\n=== products with predicted or measured Lp(a) ===\n")
print(out$res[, .(analysis, model, n, events, beta_ML, se_ML, p_ML, beta_MS, se_MS, p_MS, beta_MN, se_MN, a_N, D, se_D)], digits = 3)
report_warnings()
cat("\nPOSTHOC_DONE\n")
