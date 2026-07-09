#!/usr/bin/env Rscript
# Prepare REGENIE-format phenotype and covariate files from UKB split-CSV phenotype store.
#
# Phenotype store structure:
#   Each field is its own file: fid{field_id}.csv
#   Columns: eid, {field_id}-0.0, {field_id}-1.0, ...  (one column per UKB visit)
#   We use visit 0 (baseline, suffix -0.0) for all traits.
#
# Outputs (all in --out directory):
#   regenie_keep.txt            FID IID list (REGENIE --keep format)
#   pheno_quantitative.txt      FID IID LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA
#                                           nonHDL_C eGFR CRP HBA1C
#   pheno_binary.txt            FID IID CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD
#                                           STATIN_USE
#   covariates.txt              FID IID age sex bmi smoking_current diabetes hypertension PC1..PC10
#   phenotype_summary.txt       case/control counts + missing rates per trait
#
# Paper-driven additions (v2 — 2026-06):
#   - New quantitative: nonHDL_C (=TOT_CHOL−HDL_C), eGFR (CKD-EPI 2021),
#                       CRP (field 30710), HBA1C (field 30750)
#   - New binary: AF (I48), PAD (I70/I73/I74), STATIN_USE
#   - Fixed MACE: now includes CV_DEATH (previously only MI+STROKE+HF)
#   - REVASC: supplemented with OPCS4 PCI/CABG codes from field 41200
#
# Usage:
#   Rscript prep_phenotypes.R \
#     --pheno_dir /corral/.../pheno_split_into_files_011924 \
#     --icd_file  /corral/.../binary_ICD_011924.txt \
#     --sample    /corral/.../geno_qc_eids_400k_white_british_2023-11-29.txt \
#     --out       /path/to/output_dir

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})
set.seed(42)

opt_list <- list(
  make_option("--pheno_dir", type = "character", help = "Directory of fid*.csv files"),
  make_option("--icd_file",  type = "character", help = "binary_ICD_011924.txt"),
  make_option("--sample",    type = "character", help = "White British QC-passed EID list"),
  make_option("--out",       type = "character", default = ".", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: load a single fid*.csv and return baseline visit column ----
load_field <- function(field_id, pheno_dir, visit = "0.0") {
  f <- file.path(pheno_dir, sprintf("fid%s.csv", field_id))
  if (!file.exists(f)) {
    message(sprintf("  WARNING: fid%s.csv not found — skipping", field_id))
    return(NULL)
  }
  dt <- fread(f, na.strings = c("", "NA"))
  col <- sprintf("%s-%s", field_id, visit)
  if (!col %in% names(dt)) {
    col <- grep(sprintf("^.?%s-%s.?$", field_id, visit), names(dt), value = TRUE)[1]
  }
  if (is.na(col) || !col %in% names(dt)) {
    message(sprintf("  WARNING: column %s-%s not found in fid%s.csv", field_id, visit, field_id))
    return(NULL)
  }
  result <- dt[, .(eid = get("eid"), value = get(col))]
  setnames(result, "value", as.character(field_id))
  result
}

# ---- Load sample list ----
message("[prep] Loading sample list...")
sample_raw <- fread(opt$sample, header = FALSE)
if (ncol(sample_raw) == 1) {
  samples <- data.table(FID = sample_raw[[1]], IID = sample_raw[[1]])
} else {
  samples <- data.table(FID = sample_raw[[1]], IID = sample_raw[[2]])
}
message(sprintf("[prep] %d QC-passed samples", nrow(samples)))

# ---- Field → column name mappings ----
LIPID_FIELDS <- list(
  "30780" = "LDL_C",      # LDL direct (mmol/L)
  "30760" = "HDL_C",      # HDL cholesterol (mmol/L)
  "30870" = "TRIGLY",     # Triglycerides (mmol/L)
  "30690" = "TOT_CHOL",   # Total cholesterol (mmol/L)
  "30630" = "APOA1",      # Apolipoprotein A1 (g/L)
  "30640" = "APOB",       # Apolipoprotein B (g/L)
  "30790" = "LPA"         # Lipoprotein(a) nmol/L (field 30790 preferred over NMR 23400)
)

BIOMARKER_FIELDS <- list(
  "30700" = "CREATININE",  # Creatinine μmol/L — used for eGFR calculation
  "30710" = "CRP",         # C-reactive protein mg/L
  "30750" = "HBA1C"        # HbA1c mmol/mol
)

COVAR_FIELDS <- list(
  "21003" = "age",
  "31"    = "sex",         # 1=male, 2=female
  "21001" = "bmi",
  "20116" = "smoking",     # 0=never, 1=previous, 2=current
  "2443"  = "diabetes",    # 1=yes, 0/−7=no
  "4080"  = "sbp",         # Systolic BP (for hypertension flag)
  "6153"  = "statin_med"   # Cholesterol-lowering medication (1=yes)
)

PC_FIELD <- "22009"

# ---- Load fields ----
message("[prep] Loading lipid fields...")
pheno_list <- lapply(names(LIPID_FIELDS), function(fid) {
  dt <- load_field(fid, opt$pheno_dir)
  if (!is.null(dt)) setnames(dt, fid, LIPID_FIELDS[[fid]])
  dt
})
pheno_list <- Filter(Negate(is.null), pheno_list)

message("[prep] Loading biomarker fields...")
biomarker_list <- lapply(names(BIOMARKER_FIELDS), function(fid) {
  dt <- load_field(fid, opt$pheno_dir)
  if (!is.null(dt)) setnames(dt, fid, BIOMARKER_FIELDS[[fid]])
  dt
})
biomarker_list <- Filter(Negate(is.null), biomarker_list)

message("[prep] Loading covariate fields...")
covar_list <- lapply(names(COVAR_FIELDS), function(fid) {
  dt <- load_field(fid, opt$pheno_dir)
  if (!is.null(dt)) setnames(dt, fid, COVAR_FIELDS[[fid]])
  dt
})
covar_list <- Filter(Negate(is.null), covar_list)

# ---- Load ancestry PCs ----
message("[prep] Loading ancestry PCs...")
pc_file <- file.path(opt$pheno_dir, sprintf("fid%s.csv", PC_FIELD))
pc_list <- list()
if (file.exists(pc_file)) {
  pc_raw <- fread(pc_file, na.strings = c("", "NA"))
  for (i in 1:10) {
    col <- sprintf("%s-0.%d", PC_FIELD, i)
    if (col %in% names(pc_raw)) {
      pc_list[[i]] <- pc_raw[, .(eid, value = get(col))]
      setnames(pc_list[[i]], "value", sprintf("PC%d", i))
    }
  }
  if (length(pc_list) > 0) {
    pc_dt <- Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE), pc_list)
    message(sprintf("[prep] Loaded %d PC columns", length(pc_list)))
  } else {
    message("[prep] WARNING: No PC columns found in fid22009.csv")
    pc_dt <- NULL
  }
} else {
  message("[prep] WARNING: fid22009.csv not found — PCs must be added manually")
  pc_dt <- NULL
}

# ---- Merge everything on EID ----
message("[prep] Merging fields on EID...")
base <- data.table(eid = samples$IID, FID = samples$FID, IID = samples$IID)

merge_all <- function(base, dt_list) {
  for (dt in dt_list) {
    if (!is.null(dt)) base <- merge(base, dt, by = "eid", all.x = TRUE)
  }
  base
}

pheno_dt  <- merge_all(base, pheno_list)
biomark_dt <- merge_all(base, biomarker_list)
covar_dt  <- merge_all(base, covar_list)
if (!is.null(pc_dt)) covar_dt <- merge(covar_dt, pc_dt, by = "eid", all.x = TRUE)
# Bring biomarkers into pheno_dt — drop FID/IID from biomark_dt to avoid x/y suffix collision
biomark_extra <- setdiff(names(biomark_dt), c("eid", "FID", "IID"))
pheno_dt <- merge(pheno_dt, biomark_dt[, c("eid", biomark_extra), with = FALSE], by = "eid", all.x = TRUE)

message(sprintf("[prep] Merged: %d samples", nrow(pheno_dt)))

# ---- Statin correction for LDL-C ----
if ("statin_med" %in% names(covar_dt) && "LDL_C" %in% names(pheno_dt)) {
  statin_eids <- covar_dt[statin_med == 1, eid]
  n_statin <- length(statin_eids)
  pheno_dt[eid %in% statin_eids & !is.na(LDL_C), LDL_C := LDL_C / 0.7]
  message(sprintf("[prep] Statin correction applied: %d individuals (LDL_C / 0.7)", n_statin))
} else {
  message("[prep] WARNING: statin_med or LDL_C not found — no statin correction applied")
}

# ---- Lp(a) units check ----
if ("LPA" %in% names(pheno_dt)) {
  lpa_median <- median(pheno_dt$LPA, na.rm = TRUE)
  lpa_max    <- max(pheno_dt$LPA, na.rm = TRUE)
  message(sprintf("[prep] Lp(a) fid30790: median=%.1f max=%.1f nmol/L", lpa_median, lpa_max))
  if (lpa_max < 20) message("[prep] WARNING: Lp(a) values suspiciously low — check fid30790")
}

# ---- Derived quantitative traits ----

# nonHDL cholesterol = Total cholesterol − HDL-C  (reflects ApoB-containing particles)
if (all(c("TOT_CHOL", "HDL_C") %in% names(pheno_dt))) {
  pheno_dt[, nonHDL_C := TOT_CHOL - HDL_C]
  message(sprintf("[prep] nonHDL_C derived: median=%.2f mmol/L (n_valid=%d)",
    median(pheno_dt$nonHDL_C, na.rm = TRUE), sum(!is.na(pheno_dt$nonHDL_C))))
} else {
  message("[prep] WARNING: TOT_CHOL or HDL_C missing — nonHDL_C not derived")
}

# eGFR — CKD-EPI 2021 (race-free, Levey et al. NEJM 2021)
# Input: CREATININE (μmol/L) — convert to mg/dL (/88.4), age (years), sex (1=M,2=F)
if ("CREATININE" %in% names(pheno_dt) && "age" %in% names(covar_dt) && "sex" %in% names(covar_dt)) {
  pheno_dt <- merge(pheno_dt,
    covar_dt[, .(eid, age_for_egfr = age, sex_for_egfr = sex)], by = "eid", all.x = TRUE)
  pheno_dt[, Scr_mgdl := CREATININE / 88.4]

  # CKD-EPI 2021: κ and α differ by sex
  pheno_dt[, kappa := ifelse(sex_for_egfr == 2, 0.7, 0.9)]
  pheno_dt[, alpha := ifelse(sex_for_egfr == 2, -0.241, -0.302)]
  pheno_dt[, sex_factor := ifelse(sex_for_egfr == 2, 1.012, 1.0)]
  pheno_dt[, ratio := Scr_mgdl / kappa]
  pheno_dt[, eGFR := 142 *
    pmin(ratio, 1)^alpha *
    pmax(ratio, 1)^(-1.200) *
    (0.9938^age_for_egfr) *
    sex_factor]
  # Set to NA where inputs are missing
  pheno_dt[is.na(CREATININE) | is.na(age_for_egfr) | is.na(sex_for_egfr), eGFR := NA_real_]
  pheno_dt[, c("Scr_mgdl","kappa","alpha","sex_factor","ratio","age_for_egfr","sex_for_egfr") := NULL]
  message(sprintf("[prep] eGFR (CKD-EPI 2021) derived: median=%.1f mL/min/1.73m² (n_valid=%d)",
    median(pheno_dt$eGFR, na.rm = TRUE), sum(!is.na(pheno_dt$eGFR))))
} else {
  message("[prep] WARNING: CREATININE, age, or sex missing — eGFR not derived")
}

# ---- Build covariate file ----
if ("smoking" %in% names(covar_dt)) {
  covar_dt[, smoking_current := as.integer(smoking == 2)]
  covar_dt[, smoking := NULL]
}
if ("diabetes" %in% names(covar_dt)) {
  covar_dt[diabetes == -7, diabetes := NA]
  covar_dt[, diabetes := as.integer(diabetes == 1)]
}
if ("sbp" %in% names(covar_dt)) {
  covar_dt[, hypertension := as.integer(sbp >= 140)]
  covar_dt[, sbp := NULL]
}
covar_dt[, statin_med := NULL]  # used only for LDL correction; exposed below as phenotype

pc_names <- paste0("PC", 1:10)
pc_names_present <- pc_names[pc_names %in% names(covar_dt)]
cov_cols <- c("FID", "IID", "age", "sex", "bmi", "smoking_current",
              "diabetes", "hypertension", pc_names_present)
cov_cols <- cov_cols[cov_cols %in% names(covar_dt)]
out_covar <- covar_dt[, ..cov_cols]

# ---- Helper: derive CV_DEATH from death registry fields 40001 / 40002 ----
derive_cv_death <- function(base_eids, pheno_dir) {
  eids_char <- as.character(base_eids)
  cv_case   <- rep(0L, length(eids_char))
  found_any <- FALSE
  for (fid in c("40001", "40002")) {
    f <- file.path(pheno_dir, sprintf("fid%s.csv", fid))
    if (!file.exists(f)) { message(sprintf("  [CV_DEATH] WARNING: fid%s.csv not found", fid)); next }
    dt <- fread(f, na.strings = c("", "NA"))
    if (!"eid" %in% names(dt)) setnames(dt, names(dt)[1], "eid")
    dt[, eid := as.character(eid)]
    code_cols <- setdiff(names(dt), "eid")
    dt_sub <- dt[eid %in% eids_char]
    if (nrow(dt_sub) == 0) next
    found_any <- TRUE
    has_cv <- dt_sub[, .(eid = eid, cv = apply(.SD, 1, function(row) {
      vals <- as.character(row); any(startsWith(vals[!is.na(vals) & vals != "NA"], "I"))
    })), .SDcols = code_cols]
    cv_eids <- has_cv[cv == TRUE, eid]
    cv_case[eids_char %in% cv_eids] <- 1L
    message(sprintf("  [CV_DEATH] fid%s: %d with circulatory cause (I*)", fid, length(cv_eids)))
  }
  if (!found_any) { message("  [CV_DEATH] No death registry files — CV_DEATH set to NA"); return(rep(NA_integer_, length(eids_char))) }
  cv_case
}

# ---- Helper: derive binary outcome from ICD-10 prefixes ----
derive_outcome <- function(icd_dt, icd_prefixes, n_rows = NULL) {
  if (is.null(icd_dt)) return(rep(NA_integer_, n_rows %||% 0))
  cols <- names(icd_dt)
  matched <- unlist(lapply(icd_prefixes, function(p) grep(sprintf("^%s", p), cols, value = TRUE)))
  if (length(matched) == 0) return(rep(NA_integer_, nrow(icd_dt)))
  as.integer(rowSums(icd_dt[, ..matched], na.rm = TRUE) > 0)
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Load ICD-10 binary outcomes ----
message("[prep] Loading ICD binary outcomes...")
icd_dt <- NULL
if (file.exists(opt$icd_file)) {
  peek <- fread(opt$icd_file, nrows = 3)
  message(sprintf("[prep]   ICD file columns (first 8): %s", paste(names(peek)[1:min(8, ncol(peek))], collapse = ", ")))
  icd_dt <- fread(opt$icd_file, na.strings = c("", "NA"))
  if ("f.eid" %in% names(icd_dt)) setnames(icd_dt, "f.eid", "eid") else if (!"eid" %in% names(icd_dt)) setnames(icd_dt, names(icd_dt)[1], "eid")
  icd_dt[, eid := as.character(eid)]
  message(sprintf("[prep]   ICD file: %d samples, %d columns", nrow(icd_dt), ncol(icd_dt)))
} else {
  message(sprintf("[prep] WARNING: ICD file not found: %s", opt$icd_file))
}

# ---- OPCS4 procedure codes for REVASC ----
# Field 41200: OPCS operations (main); 41201: OPCS operations (secondary)
# PCI: K49.x (percutaneous transluminal balloon angioplasty of coronary artery)
#       K50.x (other percutaneous transluminal operations on coronary artery — stenting)
# CABG: K40.x–K46.x (saphenous vein CABG, internal mammary CABG, etc.)
derive_opcs_revasc <- function(base_eids, pheno_dir) {
  eids_char <- as.character(base_eids)
  revasc    <- rep(0L, length(eids_char))
  found_any <- FALSE
  revasc_prefixes <- c("K40","K41","K42","K43","K44","K45","K46","K49","K50")
  for (fid in c("41200","41201")) {
    f <- file.path(pheno_dir, sprintf("fid%s.csv", fid))
    if (!file.exists(f)) next
    dt <- fread(f, na.strings = c("", "NA"))
    if (!"eid" %in% names(dt)) setnames(dt, names(dt)[1], "eid")
    dt[, eid := as.character(eid)]
    code_cols <- setdiff(names(dt), "eid")
    dt_sub <- dt[eid %in% eids_char]
    if (nrow(dt_sub) == 0) next
    found_any <- TRUE
    has_revasc <- dt_sub[, .(eid = eid, rv = apply(.SD, 1, function(row) {
      vals <- as.character(row)
      any(sapply(revasc_prefixes, function(p) any(startsWith(vals[!is.na(vals) & vals != "NA"], p))))
    })), .SDcols = code_cols]
    rv_eids <- has_revasc[rv == TRUE, eid]
    revasc[eids_char %in% rv_eids] <- 1L
    message(sprintf("  [REVASC_OPCS] fid%s: %d individuals with PCI/CABG codes", fid, length(rv_eids)))
  }
  if (!found_any) { message("  [REVASC_OPCS] No OPCS files found — OPCS REVASC supplement skipped"); return(rep(NA_integer_, length(eids_char))) }
  revasc
}

# ---- Build binary phenotype table ----
bin_pheno <- data.table(FID = base$FID, IID = base$IID, eid = as.character(base$eid))

if (!is.null(icd_dt)) {
  icd_merged <- merge(bin_pheno, icd_dt, by = "eid", all.x = TRUE)

  bin_pheno[, CAD    := derive_outcome(icd_merged, c("I20","I21","I22","I23","I24","I25"))]
  bin_pheno[, MI     := derive_outcome(icd_merged, c("I21","I22"))]
  bin_pheno[, STROKE := derive_outcome(icd_merged, c("I63","I64"))]
  bin_pheno[, HF     := derive_outcome(icd_merged, c("I50"))]
  bin_pheno[, AF     := derive_outcome(icd_merged, c("I48"))]
  bin_pheno[, PAD    := derive_outcome(icd_merged, c("I70","I73","I74"))]

  # REVASC: ICD Z95/Z98 plus OPCS PCI/CABG codes
  revasc_icd  <- derive_outcome(icd_merged, c("Z95","Z98"))
  revasc_opcs <- derive_opcs_revasc(base$IID, opt$pheno_dir)
  bin_pheno[, REVASC := as.integer(
    (!is.na(revasc_icd)  & revasc_icd  == 1) |
    (!is.na(revasc_opcs) & revasc_opcs == 1)
  )]
  message(sprintf("[prep]   REVASC (ICD+OPCS combined): %d cases", sum(bin_pheno$REVASC == 1, na.rm = TRUE)))

  # CV_DEATH from death registry
  bin_pheno[, CV_DEATH := derive_cv_death(base$IID, opt$pheno_dir)]

  # MACE = MI | STROKE | HF | CV_DEATH  (v2: added CV_DEATH)
  bin_pheno[, MACE := as.integer(
    (!is.na(MI)       & MI       == 1) |
    (!is.na(STROKE)   & STROKE   == 1) |
    (!is.na(HF)       & HF       == 1) |
    (!is.na(CV_DEATH) & CV_DEATH == 1)
  )]

  # STATIN_USE binary phenotype (from field 6153 — separate from LDL correction)
  statin_dt <- load_field("6153", opt$pheno_dir)
  if (!is.null(statin_dt)) {
    statin_dt[, eid := as.character(eid)]
    bin_pheno <- merge(bin_pheno, statin_dt, by = "eid", all.x = TRUE)
    setnames(bin_pheno, "6153", "STATIN_USE")
    bin_pheno[, STATIN_USE := as.integer(STATIN_USE == 1)]
  } else {
    bin_pheno[, STATIN_USE := NA_integer_]
    message("[prep] WARNING: fid6153.csv not found — STATIN_USE set to NA")
  }

  # Report case counts
  for (outcome in c("CAD","MI","REVASC","STROKE","HF","CV_DEATH","MACE","AF","PAD","STATIN_USE")) {
    if (!outcome %in% names(bin_pheno)) next
    n_case <- sum(bin_pheno[[outcome]] == 1, na.rm = TRUE)
    n_ctrl <- sum(bin_pheno[[outcome]] == 0, na.rm = TRUE)
    n_na   <- sum(is.na(bin_pheno[[outcome]]))
    message(sprintf("[prep]   %-12s  cases=%6d  controls=%6d  NA=%d", outcome, n_case, n_ctrl, n_na))
  }

} else {
  for (col in c("CAD","MI","REVASC","STROKE","HF","CV_DEATH","MACE","AF","PAD","STATIN_USE")) {
    bin_pheno[, (col) := NA_integer_]
  }
  message("[prep] WARNING: Binary outcomes set to NA — populate from ICD data")
}

# ---- Write REGENIE --keep file ----
keep_file <- file.path(opt$out, "regenie_keep.txt")
fwrite(base[, .(FID, IID)], keep_file, sep = "\t", col.names = FALSE)
message(sprintf("[prep] Keep file: %s (%d samples)", keep_file, nrow(base)))

# ---- Write quantitative phenotype file ----
qt_core_cols   <- unname(unlist(LIPID_FIELDS))          # LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA
qt_derived_cols <- c("nonHDL_C", "eGFR", "CRP", "HBA1C")
qt_all_cols    <- c(qt_core_cols, qt_derived_cols)
qt_out_cols    <- c("FID", "IID", qt_all_cols[qt_all_cols %in% names(pheno_dt)])
out_qt         <- file.path(opt$out, "pheno_quantitative.txt")
fwrite(pheno_dt[, ..qt_out_cols], out_qt, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("[prep] Quantitative phenotypes: %s  (%d traits)", out_qt, length(qt_out_cols) - 2))

# ---- Write binary phenotype file ----
bin_out_cols <- c("FID","IID","CAD","MI","REVASC","STROKE","HF","CV_DEATH","MACE","AF","PAD","STATIN_USE")
bin_out_cols <- bin_out_cols[bin_out_cols %in% names(bin_pheno)]
out_bin <- file.path(opt$out, "pheno_binary.txt")
fwrite(bin_pheno[, ..bin_out_cols], out_bin, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("[prep] Binary phenotypes: %s  (%d traits)", out_bin, length(bin_out_cols) - 2))

# ---- Write covariate file ----
out_cov <- file.path(opt$out, "covariates.txt")
fwrite(out_covar, out_cov, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("[prep] Covariates: %s", out_cov))

# ---- Write summary file ----
summary_lines <- c(
  sprintf("Prepared: %s", Sys.time()),
  sprintf("Samples:           %d", nrow(base)),
  sprintf("Quantitative:      %s", paste(qt_out_cols[-c(1,2)], collapse=", ")),
  sprintf("Binary outcomes:   %s", paste(bin_out_cols[-c(1,2)], collapse=", ")),
  sprintf("Covariates:        %s", paste(cov_cols[-c(1,2)], collapse=", ")),
  sprintf("Ancestry PCs:      %d loaded", length(pc_names_present)),
  "",
  "== Missing rates (quantitative) =="
)
for (trait in qt_out_cols[-c(1,2)]) {
  if (trait %in% names(pheno_dt)) {
    n_miss <- sum(is.na(pheno_dt[[trait]]))
    pct    <- round(100 * n_miss / nrow(pheno_dt), 1)
    summary_lines <- c(summary_lines, sprintf("  %-12s  missing=%d (%.1f%%)", trait, n_miss, pct))
  }
}
writeLines(summary_lines, file.path(opt$out, "phenotype_summary.txt"))
message("[prep] Complete.")
