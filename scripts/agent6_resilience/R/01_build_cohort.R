#!/usr/bin/env Rscript
# Stage 01 - build the analysis cohort table.
# Output: cohort.tsv (one row per eid: covariates, statin flag, incident CAD survival,
#         raw baseline NMR), nmr_repeat.tsv (instance-1 NMR + covariates), cohort_flow.json.
# Usage: Rscript R/01_build_cohort.R --config config.yaml --pheno_dir DIR --sample FILE --out DIR
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "lib_resilience.R"))
a   <- parse_args(required = c("config", "pheno_dir", "sample", "out"))
cfg <- load_config(a$config); F <- cfg$fields
dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
set.seed(cfg$seed)
flow <- list()

# ---- 1. sample list --------------------------------------------------------
samp <- fread(a$sample, header = FALSE)
base <- data.table(eid = as.integer(samp[[1]]))
base <- unique(base)
flow$sample_list <- nrow(base)
log_msg("sample list: %d", nrow(base))

# ---- 2. scalar covariates --------------------------------------------------
add <- function(dt, field, name, numeric = TRUE, instance = 0) merge(dt, read_scalar(a$pheno_dir, field, name, instance, numeric), by = "eid", all.x = TRUE)
base <- add(base, F$assessment_date, "baseline_date", numeric = FALSE)
base[, baseline_date := as.IDate(baseline_date)]
base <- add(base, F$age, "age")
base <- add(base, F$sex, "sex")
base <- add(base, F$genetic_sex, "genetic_sex")
base <- add(base, F$genetic_ethnic_grouping, "genetic_ethnic")
base <- add(base, F$sex_aneuploidy, "sex_aneuploidy")
base <- add(base, F$het_missing_outlier, "het_missing_outlier")
base <- add(base, F$genotype_batch, "genotype_batch")
base <- add(base, F$fasting_time, "fasting_time")
base <- add(base, F$bmi, "bmi")
base <- add(base, F$smoking_status, "smoking_status")
base <- add(base, F$alcohol_freq, "alcohol_freq")
base <- add(base, F$physical_activity_ipaq, "physical_activity")
base <- add(base, F$sleep_duration, "sleep_duration")
base <- add(base, F$townsend, "townsend")
# UKB special codes are not quantities: -1 = do not know, -3 = prefer not to answer
for (v in c("smoking_status", "alcohol_freq", "physical_activity", "sleep_duration", "bmi")) {
  n_neg <- sum(base[[v]] < 0, na.rm = TRUE)
  if (n_neg > 0) { log_msg("recoding %d negative special codes to NA in %s", n_neg, v); set(base, i = which(base[[v]] < 0), j = v, value = NA) }
}
check(all(base$smoking_status %in% c(0, 1, 2, NA)), "smoking_status takes only values 0/1/2/NA")
base <- add(base, F$nmr_spectrometer, "spectrometer")
base <- add(base, F$nmr_shipment_plate, "shipment_plate", numeric = FALSE)
for (q in c("low_protein", "high_lactate", "high_pyruvate", "low_glucose"))
  base <- merge(base, read_flag(a$pheno_dir, F[[paste0("nmr_", q)]], paste0("qc_", q)), by = "eid", all.x = TRUE)
log_msg("NMR sample flags at baseline: low_protein=%d high_lactate=%d high_pyruvate=%d low_glucose=%d",
        sum(base$qc_low_protein, na.rm = TRUE), sum(base$qc_high_lactate, na.rm = TRUE), sum(base$qc_high_pyruvate, na.rm = TRUE), sum(base$qc_low_glucose, na.rm = TRUE))
base <- add(base, F$death_date, "death_date", numeric = FALSE)
base <- add(base, F$death_cause_primary, "death_cause", numeric = FALSE)
base <- add(base, F$lost_to_followup_date, "lost_date", numeric = FALSE)
base[, death_date := as.IDate(death_date)]
base[, lost_date := as.IDate(lost_date)]

pcs <- read_field(a$pheno_dir, F$ancestry_pcs, numeric = TRUE)
pc_cols <- field_cols(pcs, F$ancestry_pcs, 0)[seq_len(F$n_pcs)]
check(length(pc_cols) == F$n_pcs && all(!is.na(pc_cols)), "%d ancestry PC columns found", F$n_pcs)
pcs <- pcs[, c("eid", pc_cols), with = FALSE]; setnames(pcs, pc_cols, paste0("PC", seq_len(F$n_pcs)))
base <- merge(base, pcs, by = "eid", all.x = TRUE)

# ---- 3. statin at baseline (three sources, either sex) --------------------
st_f <- read_any_flag(a$pheno_dir, F$medication_female, "statin_6153", cfg$statin$medication_value)
st_m <- read_any_flag(a$pheno_dir, F$medication_male,   "statin_6177", cfg$statin$medication_value)
st_t <- read_any_flag(a$pheno_dir, F$treatment_codes,   "statin_20003", cfg$statin$treatment_codes)
base <- Reduce(function(x, y) merge(x, y, by = "eid", all.x = TRUE), list(base, st_f, st_m, st_t))
base[, statin := as.integer((statin_6153 %in% 1) | (statin_6177 %in% 1) | (statin_20003 %in% 1))]
check(sum(base$statin_6177 %in% 1) > 0 || nrow(base) < 5000, "male medication field contributes statin users")
log_msg("statin at baseline: %d (%.1f%%) [6153=%d, 6177=%d, 20003=%d]", sum(base$statin), 100 * mean(base$statin),
        sum(base$statin_6153 %in% 1), sum(base$statin_6177 %in% 1), sum(base$statin_20003 %in% 1))

# ---- 4. self-reported prevalent CAD ---------------------------------------
sr1 <- read_any_flag(a$pheno_dir, F$self_report_illness, "sr_20002", cfg$outcome$self_report_20002)
sr2 <- read_any_flag(a$pheno_dir, F$vascular_heart_problems, "sr_6150", cfg$outcome$self_report_6150)
base <- merge(merge(base, sr1, by = "eid", all.x = TRUE), sr2, by = "eid", all.x = TRUE)
base[, self_report_prevalent := as.integer((sr_20002 %in% 1) | (sr_6150 %in% 1))]

# ---- 5. HES + death events -> long ----------------------------------------
icd  <- wide_pairs_to_long(read_field(a$pheno_dir, F$hes_icd10_codes), read_field(a$pheno_dir, F$hes_icd10_dates),
                           F$hes_icd10_codes, F$hes_icd10_dates, "icd10")
opcs <- wide_pairs_to_long(read_field(a$pheno_dir, F$hes_opcs4_codes), read_field(a$pheno_dir, F$hes_opcs4_dates),
                           F$hes_opcs4_codes, F$hes_opcs4_dates, "opcs4")
death <- base[!is.na(death_date) & !is.na(death_cause), .(eid, source = "death", code = death_cause, date = death_date)]
events <- rbindlist(list(icd, opcs, death), use.names = TRUE)
events <- events[eid %in% base$eid]
log_msg("events long table: %d rows (icd10=%d opcs4=%d death=%d)", nrow(events), nrow(icd), nrow(opcs), nrow(death))
check(all(events$date <= as.IDate(cfg$outcome$admin_censor_date) | events$source == "death", na.rm = TRUE) ||
      mean(events$date > as.IDate(cfg$outcome$admin_censor_date), na.rm = TRUE) < 0.001,
      "essentially no HES events after the admin censor date")

surv <- build_incident_outcome(base[, .(eid, baseline_date, death_date, lost_date, self_report_prevalent)],
                               events, cfg$outcome, cfg$outcome$admin_censor_date)
base <- merge(base, surv[, .(eid, prevalent, first_event_date, end_date, event, reason, time_years)], by = "eid", all.x = TRUE)

# ---- 6. NMR baseline + repeat ----------------------------------------------
nmr_tab <- fread(file.path(dirname(cfg$.path), F$nmr_field_table))
nmr_ids <- nmr_tab[class %in% cfg$nmr$use_class, field_id]
read_nmr_instance <- function(inst) {
  cols <- lapply(nmr_ids, function(fid) {
    dt <- read_field(a$pheno_dir, fid, numeric = TRUE)
    cc <- field_cols(dt, fid, inst)
    if (!length(cc)) return(NULL)
    out <- dt[, .(eid, v = get(cc[1]))]; setnames(out, "v", paste0("nmr_", fid)); out
  })
  cols <- Filter(Negate(is.null), cols)
  Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), cols)
}
nmr0 <- read_nmr_instance(0)
check(ncol(nmr0) - 1 == length(nmr_ids), "all %d primary NMR fields read at baseline", length(nmr_ids))
nmr_cols <- setdiff(names(nmr0), "eid")
nmr0[, nmr_n_missing := rowSums(is.na(.SD)), .SDcols = nmr_cols]
base <- merge(base, nmr0, by = "eid", all.x = TRUE)
base[, has_nmr := as.integer(!is.na(nmr_n_missing) & nmr_n_missing <= cfg$nmr$max_missing_frac_per_person * length(nmr_cols))]

nmr1 <- read_nmr_instance(1)
if (!is.null(nmr1) && ncol(nmr1) > 1) {
  nmr1[, n_missing := rowSums(is.na(.SD)), .SDcols = setdiff(names(nmr1), "eid")]
  nmr1 <- nmr1[n_missing <= cfg$nmr$max_missing_frac_per_person * length(nmr_cols)]
  rep_cov <- merge(read_scalar(a$pheno_dir, F$age, "age", instance = 1), read_scalar(a$pheno_dir, F$fasting_time, "fasting_time", instance = 1), by = "eid", all = TRUE)
  rep_cov <- merge(rep_cov, read_scalar(a$pheno_dir, F$nmr_spectrometer, "spectrometer", instance = 1), by = "eid", all = TRUE)
  for (q in c("low_protein", "high_lactate", "high_pyruvate", "low_glucose"))
    rep_cov <- merge(rep_cov, read_flag(a$pheno_dir, F[[paste0("nmr_", q)]], paste0("qc_", q), instance = 1), by = "eid", all.x = TRUE)
  nmr1 <- merge(nmr1, rep_cov, by = "eid", all.x = TRUE)
  nmr1 <- nmr1[qc_low_protein %in% 0]
  nmr1 <- merge(nmr1, base[, .(eid, sex, genotype_batch)], by = "eid")
  fwrite(nmr1, file.path(a$out, "nmr_repeat.tsv"), sep = "\t", na = "NA")
  log_msg("repeat NMR: %d individuals", nrow(nmr1))
}

# ---- 7. exclusion flow (flags only; the table keeps everyone) --------------
base[, excl_no_baseline_date := as.integer(is.na(baseline_date))]
base[, excl_genetic_qc := as.integer((!is.na(sex_aneuploidy) & sex_aneuploidy == 1) |
                                     (!is.na(het_missing_outlier) & het_missing_outlier == 1) |
                                     (!is.na(genetic_sex) & !is.na(sex) & genetic_sex != sex))]
base[, excl_no_nmr := as.integer(has_nmr != 1)]
base[, excl_nmr_quality := as.integer(isTRUE(cfg$nmr$drop_low_protein) & !is.na(qc_low_protein) & qc_low_protein == 1)]
base[, excl_prevalent := as.integer(prevalent == 1)]
base[, excl_short_followup := as.integer(is.na(time_years) | time_years <= cfg$outcome$landmark_years)]
base[, eligible := as.integer(excl_no_baseline_date == 0 & excl_genetic_qc == 0 & excl_no_nmr == 0 &
                              excl_nmr_quality == 0 & excl_prevalent == 0 & excl_short_followup == 0)]

step <- function(name, keep) { flow[[name]] <<- sum(keep); keep }
k <- rep(TRUE, nrow(base))
k <- step("has_baseline_date", k & base$excl_no_baseline_date == 0)
k <- step("passes_genetic_qc", k & base$excl_genetic_qc == 0)
k <- step("has_nmr_baseline",  k & base$excl_no_nmr == 0)
k <- step("nmr_not_low_protein", k & base$excl_nmr_quality == 0)
k <- step("no_prevalent_cad",  k & base$excl_prevalent == 0)
k <- step("followup_beyond_landmark", k & base$excl_short_followup == 0)
flow$eligible_all <- sum(k)
flow$eligible_statin_free <- sum(k & base$statin == 0)
flow$events_eligible_all <- sum(base$event[k] == 1)
flow$events_eligible_statin_free <- sum(base$event[k & base$statin == 0] == 1)
check(sum(base$eligible) == flow$eligible_all, "eligible flag agrees with flow count")
check(flow$eligible_all > 0, "non-empty eligible cohort")
check(all(base[eligible == 1, time_years] > cfg$outcome$landmark_years), "eligible follow-up exceeds landmark")
check(all(base[eligible == 1 & event == 1, first_event_date > baseline_date]), "incident events strictly after baseline")
check(all(base[eligible == 1, prevalent] == 0), "no prevalent cases among eligible")
log_msg("cohort flow: %s", paste(sprintf("%s=%d", names(flow), unlist(flow)), collapse = " | "))

fwrite(base, file.path(a$out, "cohort.tsv"), sep = "\t", na = "NA")
jsonlite::write_json(flow, file.path(a$out, "cohort_flow.json"), auto_unbox = TRUE, pretty = TRUE)
write_manifest("01_build_cohort", a$out, cfg, inputs = a$sample,
               outputs = file.path(a$out, c("cohort.tsv", "cohort_flow.json")),
               extra = list(flow = flow, n_nmr_fields = length(nmr_cols)))
log_msg("Stage 01 complete.")
