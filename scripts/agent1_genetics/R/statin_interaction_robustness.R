#!/usr/bin/env Rscript
# ============================================================================
# Robustness of the baseline-statin x genotype interaction on INCIDENT CAD.
#
# Reviewer points addressed here, in order:
#  (1) proportional hazards for score, statin and score x statin (cox.zph)
#  (2) adjustment for baseline treatment-indication covariates
#  (3) covariate x score and covariate x statin terms, because a gene-environment
#      interaction can be manufactured by unmodelled covariate interactions
#  (4) propensity weighting (stabilized IPTW) for non-random statin assignment
#  (5) both scores in one model, mutually adjusted
#
# Data prep is identical to statin_interaction_incident.R (prevalent disease
# excluded, events strictly after baseline, baseline statin report).
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse); library(survival) })
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--scores",  type = "character"),
  make_option("--weights", type = "character"),
  make_option("--raw_dir", type = "character"),
  make_option("--covar",   type = "character"),
  make_option("--out",     type = "character"),
  make_option("--censor",  type = "character", default = "2022-10-31")
)))

rdcsv <- function(fid) {
  f <- file.path(opt$raw_dir, sprintf("fid%s.csv", fid))
  d <- fread(f, na.strings = c("", "NA"))
  col <- grep(sprintf("^\"?%s-0\\.0\"?$", fid), names(d), value = TRUE)[1]
  setnames(d, c("eid", col), c("IID", "v")); d[, .(IID = as.character(IID), v)]
}
build_score <- function(raw_file, weight_file, label) {
  g <- fread(raw_file); w <- fread(weight_file, header = FALSE, col.names = c("SNP","EA","BETA"))
  dc <- setdiff(names(g), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
  map <- merge(data.table(col = dc, SNP = sub("_[^_]+$", "", dc), counted = sub("^.*_", "", dc)),
               w, by = "SNP")
  sc <- numeric(nrow(g)); flip <- 0L
  for (i in seq_len(nrow(map))) {
    d <- g[[map$col[i]]]
    if (toupper(map$counted[i]) != toupper(map$EA[i])) { d <- 2 - d; flip <- flip + 1L }
    d[is.na(d)] <- mean(d, na.rm = TRUE); sc <- sc + d * map$BETA[i]
  }
  message(sprintf("  [%s] %d variants (%d flipped)", label, nrow(map), flip))
  data.table(IID = as.character(g$IID), s = sc)
}

raws <- strsplit(opt$scores, ",")[[1]]; wts <- strsplit(opt$weights, ",")[[1]]
sc <- NULL
for (k in seq_along(raws)) {
  nm <- if (grepl("LPA", raws[k])) "score_LPA" else "score_LDL"
  d <- build_score(raws[k], wts[k], nm); setnames(d, "s", nm)
  sc <- if (is.null(sc)) d else merge(sc, d, by = "IID")
}
base <- rdcsv("53");    setnames(base, "v", "t0")
dth  <- rdcsv("40000"); setnames(dth,  "v", "tdeath")
CADF <- c("131296","131298","131300","131302","131304","131306")
ev <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
             lapply(CADF, function(f) { d <- rdcsv(f); setnames(d,"v",paste0("d",f)); d }))
ev[, tcad := do.call(pmin, c(lapply(paste0("d",CADF), function(c) as.IDate(get(c))), na.rm = TRUE))]
st <- Reduce(function(a,b) merge(a,b,by="IID",all=TRUE),
             list(rdcsv("6153")[, .(IID, f = v)], rdcsv("6177")[, .(IID, m = v)]))
st[, statin := as.integer((!is.na(f) & f == 1) | (!is.na(m) & m == 1))]
st[is.na(f) & is.na(m), statin := NA_integer_]
cv <- fread(opt$covar); cv[, IID := as.character(IID)]
d <- Reduce(function(a,b) merge(a,b,by="IID"), list(sc, base[, .(IID,t0)], ev[, .(IID,tcad)],
                                                    st[, .(IID,statin)], cv))
d <- merge(d, dth[, .(IID,tdeath)], by = "IID", all.x = TRUE)
d[, `:=`(t0 = as.IDate(t0), tdeath = as.IDate(tdeath))]
cens <- as.IDate(opt$censor)
d <- d[!is.na(t0) & !is.na(statin)]
prev <- d[!is.na(tcad) & tcad <= t0, .N]
d <- d[is.na(tcad) | tcad > t0]
d[, tend  := pmin(fifelse(is.na(tcad), cens, tcad), fifelse(is.na(tdeath), cens, tdeath), cens, na.rm = TRUE)]
d[, event := as.integer(!is.na(tcad) & tcad <= tend)]
d[, fu    := as.numeric(tend - t0) / 365.25]
d <- d[fu > 0]
pcs <- grep("^PC", names(d), value = TRUE)[1:10]
IND <- c("bmi","smoking_current","diabetes","hypertension")
d <- d[complete.cases(d[, ..IND])]
d[, score_LPA := scale(score_LPA)[,1]]; d[, score_LDL := scale(score_LDL)[,1]]
message(sprintf("[rob] analysed %d, events %d, statin users %d (%.1f%%), prevalent excluded %d",
                nrow(d), sum(d$event), sum(d$statin == 1), 100*mean(d$statin == 1), prev))

PCS <- paste(pcs, collapse = " + ")
grab <- function(m, term, model, note = "") {
  cf <- summary(m)$coefficients
  i <- which(rownames(cf) == term)
  if (!length(i) && grepl(":", term, fixed = TRUE)) {
    # R labels an interaction by the order its variables appear in the formula's
    # variable list, not the order written. In
    #   ~ score_LPA * statin + score_LDL * statin
    # `statin` is introduced by the FIRST product, so the second interaction is
    # named "statin:score_LDL", not "score_LDL:statin". Exact matching therefore
    # silently dropped that coefficient from the mutually-adjusted model, and the
    # figure quietly plotted the row with one arm missing. Try the reversal.
    i <- which(rownames(cf) == paste(rev(strsplit(term, ":", fixed = TRUE)[[1]]),
                                     collapse = ":"))
  }
  if (!length(i)) {
    warning(sprintf("grab(): term '%s' not found in model '%s'; available: %s",
                    term, model, paste(rownames(cf), collapse = ", ")),
            call. = FALSE)
    return(NULL)
  }
  se <- if ("robust se" %in% colnames(cf)) cf[i, "robust se"] else cf[i, "se(coef)"]
  z  <- cf[i, "coef"] / se
  data.table(model = model, term = term, estimate = exp(cf[i, "coef"]),
             lo95 = exp(cf[i, "coef"] - 1.96*se), hi95 = exp(cf[i, "coef"] + 1.96*se),
             P = 2*pnorm(-abs(z)), n = m$n, events = m$nevent, note = note)
}
res <- list(); ph <- list()

for (S in c("score_LPA", "score_LDL")) {
  # (0) primary, as published
  f0 <- as.formula(sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", S, PCS))
  m0 <- coxph(f0, data = d)
  res[[length(res)+1]] <- grab(m0, paste0(S, ":statin"), "primary (as reported)")
  res[[length(res)+1]] <- grab(m0, S, "primary (as reported)")

  # (1) proportional hazards
  z <- cox.zph(m0)
  tb <- as.data.table(z$table, keep.rownames = "term")
  tb[, `:=`(score = S, model = "primary")]
  ph[[length(ph)+1]] <- tb

  # (2) + baseline indication covariates
  f1 <- as.formula(sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s + %s",
                           S, paste(IND, collapse = " + "), PCS))
  m1 <- coxph(f1, data = d)
  res[[length(res)+1]] <- grab(m1, paste0(S, ":statin"), "+ indication covariates",
                               "bmi, smoking, diabetes, hypertension as main effects")

  # (3) + covariate x score and covariate x statin interactions
  f2 <- as.formula(sprintf(
    "Surv(fu, event) ~ %s * statin + age + sex + %s + %s + %s + %s",
    S, paste(IND, collapse = " + "),
    paste(sprintf("%s:%s", S, IND), collapse = " + "),
    paste(sprintf("statin:%s", IND), collapse = " + "), PCS))
  m2 <- coxph(f2, data = d)
  res[[length(res)+1]] <- grab(m2, paste0(S, ":statin"), "+ covariate interactions",
                               "score x covariate and statin x covariate terms included")

  # (4) stabilized IPTW for statin assignment
  ps <- glm(as.formula(sprintf("statin ~ age + sex + %s + %s", paste(IND, collapse = " + "), PCS)),
            data = d, family = binomial())
  e  <- fitted(ps); pt <- mean(d$statin)
  d[, sw := fifelse(statin == 1, pt/e, (1-pt)/(1-e))]
  d[, sw := pmin(sw, quantile(sw, 0.99))]   # truncate the top 1% of weights
  cstat <- mean(sample(e[d$statin == 1], 2e4, TRUE) > sample(e[d$statin == 0], 2e4, TRUE))
  m3 <- coxph(as.formula(sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", S, PCS)),
              data = d, weights = sw, robust = TRUE)
  res[[length(res)+1]] <- grab(m3, paste0(S, ":statin"), "stabilized IPTW",
                               sprintf("propensity c-statistic %.3f; weights truncated at the 99th centile", cstat))
}

# (5) both scores mutually adjusted, each with its own statin interaction
m4 <- coxph(as.formula(sprintf(
  "Surv(fu, event) ~ score_LPA * statin + score_LDL * statin + age + sex + %s", PCS)), data = d)
for (t in c("score_LPA:statin", "score_LDL:statin", "score_LPA", "score_LDL"))
  res[[length(res)+1]] <- grab(m4, t, "both scores, mutually adjusted")

out <- rbindlist(Filter(Negate(is.null), res))
print(out[, .(model, term, est = round(estimate, 4), lo = round(lo95, 4),
              hi = round(hi95, 4), P = signif(P, 3))])
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(out, opt$out)
phdt <- rbindlist(ph, fill = TRUE)
print(phdt)
fwrite(phdt, sub("\\.csv$", "_ph.csv", opt$out))
message("[rob] wrote ", opt$out)
