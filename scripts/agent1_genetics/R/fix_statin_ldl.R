#!/usr/bin/env Rscript
# ============================================================
# Rebuild LDL-C with a correctly ascertained statin correction.
#
# The bug: prep_phenotypes.R read UK Biobank field 6153 only. UKB splits the
# medication touchscreen question by sex -- 6153 is asked of women, 6177 of men
# -- so statin use was ascertained in women ONLY and every male was NA. The
# LDL-C statin correction (LDL/0.7) therefore reached 27,549 women and none of
# the 42,678 male users. True cohort statin use is 70,227 (17.5%), matching the
# published UKB baseline rate.
#
# This script verifies the existing LDL_C construction, then rebuilds it from
# the raw field with the union-based correction. Everything else is untouched.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
set.seed(42)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--pheno",     type="character"),
  make_option("--raw_dir",   type="character"),
  make_option("--out",       type="character"),
  make_option("--factor",    type="double", default=0.7)
)))

rd <- function(fid, dir) {
  f <- file.path(dir, sprintf("fid%s.csv", fid))
  if (!file.exists(f)) stop("missing ", f)
  d <- fread(f, na.strings=c("","NA"))
  col <- grep(sprintf("^\"?%s-0\\.0\"?$", fid), names(d), value=TRUE)[1]
  if (is.na(col)) stop("no -0.0 column in ", f)
  setnames(d, c("eid", col), c("eid","v"))
  d[, .(eid=as.character(eid), v)]
}

ph <- fread(opt$pheno, na.strings="NA")
ph[, eid := as.character(IID)]
ldl_raw <- rd("30780", opt$raw_dir); setnames(ldl_raw, "v", "ldl_raw")
f6153   <- rd("6153",  opt$raw_dir); setnames(f6153,   "v", "med_f")
f6177   <- rd("6177",  opt$raw_dir); setnames(f6177,   "v", "med_m")

d <- Reduce(function(a,b) merge(a,b,by="eid",all.x=TRUE),
            list(ph[, .(eid, LDL_C)], ldl_raw, f6153, f6177))
d[, statin_f   := !is.na(med_f) & med_f == 1]
d[, statin_any := (!is.na(med_f) & med_f==1) | (!is.na(med_m) & med_m==1)]

message(sprintf("[fix] statin users: women-only(old)=%d  union(new)=%d  newly corrected=%d",
        sum(d$statin_f), sum(d$statin_any), sum(d$statin_any & !d$statin_f)))

# --- verify the existing column really is raw, with /0.7 applied to women users only ---
d[, expected_old := fifelse(statin_f, ldl_raw/opt$factor, ldl_raw)]
chk <- d[!is.na(LDL_C) & !is.na(expected_old)]
mad <- max(abs(chk$LDL_C - chk$expected_old), na.rm=TRUE)
message(sprintf("[fix] reconstruction check on %d rows: max|observed-expected| = %.3g", nrow(chk), mad))
if (mad > 1e-6) {
  message("[fix] FAILED: current LDL_C is NOT raw-with-female-only-correction.")
  message("[fix] Refusing to rebuild -- the construction differs from the assumed model.")
  quit(status = 1)
}
message("[fix] verified: current LDL_C = raw, /0.7 applied to female statin users only")

# --- rebuild with correct ascertainment ---
d[, LDL_C_fixed := fifelse(statin_any, ldl_raw/opt$factor, ldl_raw)]
ch <- d[!is.na(LDL_C) & !is.na(LDL_C_fixed) & abs(LDL_C - LDL_C_fixed) > 1e-9]
message(sprintf("[fix] values changed: %d (all male statin users); mean shift %+.3f mmol/L",
        nrow(ch), mean(ch$LDL_C_fixed - ch$LDL_C)))
message(sprintf("[fix] cohort mean LDL_C: %.4f -> %.4f mmol/L",
        mean(d$LDL_C, na.rm=TRUE), mean(d$LDL_C_fixed, na.rm=TRUE)))

out <- copy(ph)
out[, LDL_C := d$LDL_C_fixed[match(out$eid, d$eid)]]
out[, eid := NULL]
fwrite(out, opt$out, sep="\t", na="NA", quote=FALSE)
message("[fix] wrote ", opt$out)
