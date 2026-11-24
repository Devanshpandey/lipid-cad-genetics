#!/usr/bin/env Rscript
# ============================================================
# Lp(a) phenotype that keeps participants outside the assay's reportable range.
#
# UKB field 30790 has no value outside 3.8-189 nmol/L. The reason is in field 30796, data-coding 4917:
#   1 reportable; 2 / 3 reportable at assay but not after corrections (too low / too high);
#   4 / 5 not reportable at assay (too low / too high).
# prep_phenotypes.R reads 30790 only, so the published Lp(a) GWAS dropped both tails. This script writes two
# columns over the SAME analysis sample as that GWAS:
#   LPA_MEAS  RINT over measured values only. Same ranks as REGENIE --apply-rint: the reproduction control.
#   LPA_CENS  rint_censored(): the too-low / too-high participants kept as tied blocks.
#
# Scores are computed here, so REGENIE must run WITHOUT --apply-rint: it would re-rank the tied blocks.
# REGENIE ranks after intersecting fam, keep and complete covariates, so we do the same; --expect_n and
# --expect_measured are the counts in the published REGENIE log and are hard gates.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--raw_dir", type = "character", help = "Directory of fid*.csv"),
  make_option("--keep",    type = "character", help = "regenie_keep.txt of the published GWAS"),
  make_option("--fam",     type = "character", help = "fam of the step 2 genotypes"),
  make_option("--covar",   type = "character", help = "covariates.txt of the published GWAS"),
  make_option("--pheno_q", type = "character", help = "pheno_quantitative_statinfix.txt"),
  make_option("--expect_n",        type = "integer", default = 373164L),
  make_option("--expect_measured", type = "integer", default = 283714L),
  make_option("--lib",     type = "character", help = "lib_rint_censored.R"),
  make_option("--out",     type = "character", help = "Output directory")
)))
source(opt$lib)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
gate <- function(ok, ...) if (!isTRUE(ok)) stop("GATE: ", sprintf(...), call. = FALSE)

read_field <- function(fid) {
  d <- fread(file.path(opt$raw_dir, sprintf("fid%s.csv", fid)), na.strings = c("", "NA"))
  col <- grep(sprintf("^\"?%s-0\\.0\"?$", fid), names(d), value = TRUE)
  gate(length(col) == 1L, "fid%s has %d instance-0 columns", fid, length(col))
  d[, .(IID = as.character(eid), v = as.numeric(get(col)))]
}

keep <- fread(opt$keep, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
fam  <- fread(opt$fam,  header = FALSE, select = 1:2, col.names = c("FID", "IID"), colClasses = "character")
cv   <- fread(opt$covar, colClasses = c(FID = "character", IID = "character"))
cc   <- cv[complete.cases(cv), IID]
analysis <- keep[IID %in% fam$IID & IID %in% cc]
q <- fread(opt$pheno_q, colClasses = c(FID = "character", IID = "character"))
# REGENIE's published 373,164 also needs at least one non-missing trait in the 11-trait file. That cannot change
# the Lp(a) ranks (anyone measured has a trait), so ranks are over `analysis` and this only checks the sample.
tr <- setdiff(names(q), c("FID", "IID"))
any_trait <- q[rowSums(!is.na(q[, ..tr])) > 0L, IID]
n_pub <- analysis[IID %in% any_trait, .N]
gate(n_pub == opt$expect_n, "published-definition sample %d, published REGENIE log says %d", n_pub, opt$expect_n)

val  <- read_field("30790"); setnames(val, "v", "raw")
code <- read_field("30796"); setnames(code, "v", "code")
# prep_phenotypes.R wrote IIDs with fwrite, which puts a round EID in scientific notation (3e+06), so field joins
# use a numeric key. IID stays as written: REGENIE matches those strings against the fam, and so must `analysis`.
ekey <- function(x) sprintf("%.0f", as.numeric(x))
d <- copy(keep)[, key := ekey(IID)]
d <- merge(d, val[, .(key = ekey(IID), raw)], by = "key", all.x = TRUE)
d <- merge(d, code[, .(key = ekey(IID), code)], by = "key", all.x = TRUE)
gate(nrow(d) == nrow(keep), "EID join changed the row count")
n_sci <- sum(grepl("e+", keep$IID, fixed = TRUE))
n_nofield <- sum(!d$key %in% ekey(val$IID))

# ---- gates, as Causal Atlas censored_masks() ----
bad <- setdiff(unique(na.omit(d$code)), 1:5)
gate(!length(bad), "codes outside coding 4917: %s", paste(head(bad), collapse = ","))
d[, low := code %in% c(2, 4)]
d[, high := code %in% c(3, 5)]
gate(any(d$low) && any(d$high), "nobody coded too low / too high: field or format problem")
gate(d[(low | high) & !is.na(raw), .N] == 0L, "%d censored participants have a value", d[(low | high) & !is.na(raw), .N])
gate(d[!is.na(raw) & is.na(code), .N] <= 0.01 * d[!is.na(raw), .N], "more than 1%% of measured have no 30796 code")
eps <- 1e-9 * 189
gate(d[raw < 3.8 - eps | raw > 189 + eps, .N] == 0L, "measured values outside 3.8-189")

# the phenotype file behind the published GWAS must be exactly this field
m <- merge(q[, .(IID, LPA)], d[, .(IID, raw)], by = "IID")
gate(nrow(m) == nrow(q) && m[is.na(LPA) != is.na(raw), .N] == 0L && m[!is.na(LPA), max(abs(LPA - raw))] < 1e-6,
     "pheno_q LPA is not field 30790")

# ---- scores, over the analysis sample ----
a <- d[IID %in% analysis$IID]
r_meas <- rint_censored(a$raw)
r_cens <- rint_censored(a$raw, a$low, a$high)
a[, `:=`(LPA_MEAS = r_meas$z, LPA_CENS = r_cens$z)]
gate(a[!is.na(LPA_MEAS), .N] == opt$expect_measured, "%d measured in sample, published REGENIE log says %d",
     a[!is.na(LPA_MEAS), .N], opt$expect_measured)
gate(abs(a[, mean(LPA_CENS, na.rm = TRUE)]) < 1e-3, "LPA_CENS mean is not 0")
gate(abs(a[, sd(LPA_CENS, na.rm = TRUE)] - r_cens$expected_sd) < 2e-3, "LPA_CENS SD differs from expected")
gate(a[low == TRUE, max(LPA_CENS)] < a[is.na(raw) == FALSE, min(LPA_CENS)] &&
     a[high == TRUE, min(LPA_CENS)] > a[is.na(raw) == FALSE, max(LPA_CENS)], "blocks do not bracket measured scores")

# REGENIE keeps fam order irrelevant; FID/IID from keep
fwrite(a[, .(FID, IID, LPA_MEAS, LPA_CENS)], file.path(opt$out, "pheno_lpa_censored.txt"), sep = "\t", na = "NA", quote = FALSE)
a[, side := fifelse(low, -1L, fifelse(high, 1L, fifelse(!is.na(raw), 0L, NA_integer_)))]
fwrite(a[, .(IID, raw, side, LPA_MEAS, LPA_CENS)], file.path(opt$out, "lpa_tobit_input.tsv"), sep = "\t", na = "NA", quote = FALSE)

# mvmr_covariance.R reads its phenotypic correlation from a table with an LPA column
q2 <- copy(q)
q2[, LPA := a$LPA_CENS[match(IID, a$IID)]]
fwrite(q2, file.path(opt$out, "pheno_quantitative_statinfix_lpacens.txt"), sep = "\t", na = "NA", quote = FALSE)

# ---- aggregates only ----
n_all <- a[!is.na(LPA_CENS), .N]
p_lo <- a[low == TRUE, .N] / n_all; p_hi <- a[high == TRUE, .N] / n_all
s <- c(
  sprintf("keep file:        %d (%d IIDs in scientific notation; %d absent from the phenotype store)", nrow(keep), n_sci, n_nofield),
  sprintf("analysis sample:  %d (fam, keep, complete covariates); %d with any of the 11 traits (published)", nrow(a), n_pub),
  "field 30796 code counts, keep file:",
  capture.output(print(d[, .N, keyby = code])),
  "field 30796 code counts, analysis sample:",
  capture.output(print(a[, .N, keyby = code])),
  sprintf("measured:         %d", a[side == 0L, .N]),
  sprintf("too low:          %d (%.2f%% of analysable)", a[low == TRUE, .N], 100 * p_lo),
  sprintf("too high:         %d (%.2f%% of analysable)", a[high == TRUE, .N], 100 * p_hi),
  sprintf("never measured:   %d", a[is.na(LPA_CENS), .N]),
  sprintf("LPA_CENS N:       %d", n_all),
  sprintf("low block score:  %.4f ; quantile bound %.4f", a[low == TRUE, LPA_CENS[1]], qnorm(p_lo)),
  sprintf("high block score: %.4f ; quantile bound %.4f", a[high == TRUE, LPA_CENS[1]], qnorm(p_hi, lower.tail = FALSE)),
  sprintf("LPA_CENS mean %.5f SD %.5f (expected %.5f)", a[, mean(LPA_CENS, na.rm = TRUE)], a[, sd(LPA_CENS, na.rm = TRUE)], r_cens$expected_sd),
  sprintf("LPA_MEAS mean %.5f SD %.5f", a[, mean(LPA_MEAS, na.rm = TRUE)], a[, sd(LPA_MEAS, na.rm = TRUE)]),
  sprintf("measured in-range SD on the censored score: %.4f", a[side == 0L, sd(LPA_CENS)])
)
writeLines(s, file.path(opt$out, "lpa_censored_summary.txt"))
writeLines(s)
message("[prep_lpa] done")
