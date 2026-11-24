#!/usr/bin/env Rscript
# ============================================================
# Statin x genotype interaction on INCIDENT CAD.
#
# WHY THIS SUPERSEDES statin_interaction.R
#   The earlier analysis used a prevalent+incident CAD binary with statin exposure
#   measured at baseline. If an event preceded the prescription, a genotype x statin
#   interaction can reflect treatment selection or reverse causation rather than any
#   modification of genetic risk: people with prior CAD are precisely the people put on
#   statins. Restricting to events AFTER baseline removes that specific failure mode.
#
# DESIGN
#   - exclude prevalent CAD/MI/revascularization (first-occurrence date <= baseline)
#   - outcome = first CAD event strictly after the baseline assessment date
#   - statin exposure = baseline medication report, which therefore precedes follow-up
#   - Cox proportional hazards, time from baseline to first event / death / censor
#   - the interaction is still fitted on the FULL eligible cohort (no restriction to
#     statin users), so it does not condition on the collider
#
# WHAT THIS STILL DOES NOT DO
#   Baseline statin status is a single point measurement: it cannot capture initiation,
#   discontinuation, adherence or intensity over follow-up. The estimand is therefore
#   effect modification associated with BASELINE statin use, not a causal treatment
#   effect. Report it as such.
# ============================================================
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
  setnames(d, c("eid", col), c("IID", "v"))
  d[, .(IID = as.character(IID), v)]
}

# ---- scores, built with explicit allele orientation (see statin_interaction.R) ----
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

# ---- dates ----
base <- rdcsv("53");    setnames(base, "v", "t0")
dth  <- rdcsv("40000"); setnames(dth,  "v", "tdeath")
CADF <- c("131296","131298","131300","131302","131304","131306")  # I20 I21 I22 I23 I24 I25
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

n0 <- nrow(d)
d <- d[!is.na(t0) & !is.na(statin)]
# EXCLUDE PREVALENT DISEASE — the whole point of this reanalysis
prev <- d[!is.na(tcad) & tcad <= t0, .N]
d <- d[is.na(tcad) | tcad > t0]
d[, tend  := pmin(fifelse(is.na(tcad), cens, tcad), fifelse(is.na(tdeath), cens, tdeath), cens, na.rm = TRUE)]
d[, event := as.integer(!is.na(tcad) & tcad <= tend)]
d[, fu    := as.numeric(tend - t0) / 365.25]
d <- d[fu > 0]
pcs <- grep("^PC", names(d), value = TRUE)[1:10]
d[, score_LPA := scale(score_LPA)[,1]]; d[, score_LDL := scale(score_LDL)[,1]]

message(sprintf("[inc] start %d -> excluded %d prevalent CAD -> analysed %d", n0, prev, nrow(d)))
message(sprintf("[inc] incident events %d | statin users %d (%.1f%%) | median follow-up %.1f y",
        sum(d$event), sum(d$statin == 1), 100*mean(d$statin == 1), median(d$fu)))

fit <- function(score) {
  f <- as.formula(sprintf("Surv(fu, event) ~ %s * statin + age + sex + %s", score, paste(pcs, collapse = " + ")))
  m <- coxph(f, data = d); cf <- summary(m)$coefficients
  list(main = cf[score, ], inter = cf[grep(":statin", rownames(cf)), ], n = m$n, e = m$nevent)
}
row <- function(lbl, r) data.table(term = lbl, HR = exp(r[1]),
  lo95 = exp(r[1] - 1.96*r[3]), hi95 = exp(r[1] + 1.96*r[3]), P = r[5])

L <- fit("score_LPA"); D <- fit("score_LDL")
res <- rbindlist(list(
  row("Lp(a) score -> incident CAD (HR/SD)",            L$main),
  row("Lp(a) score x baseline statin (ratio of HRs)",   L$inter),
  row("LDL-C score -> incident CAD (HR/SD)",            D$main),
  row("LDL-C score x baseline statin (ratio of HRs)",   D$inter)))
res[, `:=`(n = L$n, events = L$e, prevalent_excluded = prev)]
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, opt$out); print(res)

strat <- rbindlist(lapply(c(0,1), function(k) {
  m <- coxph(as.formula(sprintf("Surv(fu, event) ~ score_LPA + age + sex + %s", paste(pcs, collapse=" + "))),
             data = d[statin == k]); cf <- summary(m)$coefficients["score_LPA", ]
  data.table(stratum = ifelse(k==1,"baseline statin users","non-users"), n = m$n, events = m$nevent,
             HR = exp(cf[1]), lo95 = exp(cf[1]-1.96*cf[3]), hi95 = exp(cf[1]+1.96*cf[3]), P = cf[5])
}))
print(strat); fwrite(strat, sub("\\.csv$", "_stratified.csv", opt$out))
message("[inc] wrote ", opt$out)
