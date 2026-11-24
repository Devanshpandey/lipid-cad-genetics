#!/usr/bin/env Rscript
# Build a type-2 diabetes case/control phenotype for the PDE3B target-safety test.
#
# Design notes that matter:
#  * Cases   = ICD-10 E11, or self-reported doctor-diagnosed diabetes (field 2443)
#              in someone with no E10 code.
#  * Excluded (NA, not control) = E10 (type 1) and E12/E13/E14 (malnutrition-
#              related / other / unspecified). Calling ambiguous diabetes a
#              CONTROL would bias the case-control contrast toward the null,
#              which is the direction that would falsely reassure us here.
#  * Controls = no E10-E14 code and self-report == 0.
#  * Covariates DROP `diabetes` (it defines the outcome) and DROP `bmi`
#              (a mediator of the PDE3B->T2D path; adjusting for it would
#              attenuate exactly the effect we are testing). BMI-adjusted is
#              written separately as a sensitivity covariate file.

suppressPackageStartupMessages({library(data.table); library(optparse)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--icd_file",   type="character"),
  make_option("--covar_file", type="character"),
  make_option("--out_pheno",  type="character"),
  make_option("--out_covar",  type="character"),
  make_option("--out_covar_bmi", type="character")
)))

covar <- fread(opt$covar_file)
message(sprintf("[t2d] covariates: %d samples", nrow(covar)))

keep <- c("Patient EID","E10","E11","E12","E13","E14")
icd  <- fread(opt$icd_file, select = keep)
setnames(icd, "Patient EID", "IID")
icd[, IID := as.character(IID)]
covar[, IID := as.character(IID)]
message(sprintf("[t2d] ICD rows: %d", nrow(icd)))

m <- merge(covar[, .(FID, IID, diabetes)], icd, by = "IID", all.x = TRUE)
for (cc in c("E10","E11","E12","E13","E14")) m[is.na(get(cc)), (cc) := 0L]

m[, T2D := NA_integer_]
m[E11 == 1 | (diabetes == 1 & E10 == 0), T2D := 1L]
m[is.na(T2D) & E10 == 0 & E12 == 0 & E13 == 0 & E14 == 0 &
    !is.na(diabetes) & diabetes == 0, T2D := 0L]
# anyone left (E10 / E12-E14 only, or missing self-report) stays NA
m[E10 == 1, T2D := NA_integer_]

message(sprintf("[t2d] cases=%d  controls=%d  excluded/NA=%d",
                sum(m$T2D == 1, na.rm=TRUE), sum(m$T2D == 0, na.rm=TRUE), sum(is.na(m$T2D))))
message(sprintf("[t2d]   E11 ICD cases: %d | self-report-only cases: %d | T1D excluded: %d",
                sum(m$E11 == 1), sum(m$E11 == 0 & m$diabetes == 1 & m$E10 == 0, na.rm=TRUE),
                sum(m$E10 == 1)))

fwrite(m[!is.na(T2D), .(FID, IID, T2D)], opt$out_pheno, sep = "\t", na = "NA", quote = FALSE)

pcs <- grep("^PC", names(covar), value = TRUE)
# na = "NA": fwrite writes NA as an EMPTY field by default, which REGENIE
# rejects with "incorrectly formatted covariate file". The upstream
# covariates.txt carries literal NA (29,012 of them) and REGENIE simply
# drops those samples.
# quote = FALSE: setting na to a non-empty string flips fwrite's quote="auto"
# into quoting mode, which quotes the HEADER too ("FID" "IID" ...). REGENIE
# then rejects the file with 'header ... must start with: FID IID'. Both
# arguments are required together.
fwrite(covar[, c("FID","IID","age","sex", pcs), with = FALSE], opt$out_covar, sep = "\t", na = "NA", quote = FALSE)
fwrite(covar[, c("FID","IID","age","sex","bmi", pcs), with = FALSE], opt$out_covar_bmi, sep = "\t", na = "NA", quote = FALSE)
message("[t2d] wrote ", opt$out_pheno, " / ", opt$out_covar, " / ", opt$out_covar_bmi)
