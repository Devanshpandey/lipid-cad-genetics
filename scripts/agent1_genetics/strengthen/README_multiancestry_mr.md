# Multi-ancestry causal MR — ready-to-run pipeline

**Goal:** within-ancestry two-sample MR of lipids → CAD in African and Hispanic
ancestries: `GLGC-{ANC} lipid (exposure) → MVP-{ANC} CAD (outcome)`.

## Status
- **Exposure (lipids by ancestry): already on TACC** — `data/external/glgc/`
  (GLGC 2021; LDL/HDL/logTG for AFR + HIS; includes MVP). Downloaded, verified.
- **Outcome (CAD by ancestry): NOT yet present** — needs one access step.

## The one blocking step: dbGaP/EGA phs002453
The MVP ancestry-stratified CAD summary statistics are registered on the GWAS
Catalog but the files are gated:
- `GCST90132304` = MVP **African American** CAD
- `GCST90132303` = MVP **Hispanic/Latin American** CAD
- Distributed via **dbGaP/EGA phs002453** ("MVP Summary Results from Non-Sensitive
  Omics Studies") — summary-level, **no VA affiliation needed**, but a data-access
  request (institutional + DUC signoff) is required. It is *not* an open `wget`.

## Run in one pass once the CAD files land
1. Download `GCST90132304` (AFR) and `GCST90132303` (HIS) via phs002453 to
   `data/external/mvp_cad/`.
2. Open `run_multiancestry_mr.sh`, set the two `CAD[...]` paths and the six
   `C*` column names to match the downloaded headers (GWAS-Catalog *harmonised*
   files use `hm_rsid / hm_effect_allele / hm_other_allele / hm_beta /
   standard_error / p_value`; raw files differ — check the header).
3. `bash run_multiancestry_mr.sh`

Output: `multiancestry_mr/multiancestry_mr_results.csv` with, for each
ancestry × trait: **IVW, MR-Egger (+intercept), weighted median, MR-RAPS,
MR-PRESSO global test** (OR per SD lipid, 95% CI, P).

## Instruments & effect scale
- Instruments = GW-significant (P<5e-8; auto-relaxed to 1e-6 if <5) in the
  **target-ancestry** GLGC lipid GWAS, greedy 10 Mb clumping (matches the rest of
  the pipeline). GLGC effects are per-SD (inverse-normal), CAD beta is log-OR →
  IVW slope exponentiates to **OR per SD**.

## Sample overlap (important)
GLGC-{ANC} lipids and MVP-{ANC} CAD **share the MVP participants**, so the two
samples overlap. IVW/median/Egger/RAPS are reported as primary; for a formal
overlap correction also run one of:
- **CAUSE** (`install.packages("cause")`): models shared-sample correlation +
  correlated pleiotropy. Needs **genome-wide** GLGC + MVP-CAD sumstats (not just
  instruments) and LD pruning.
- **MRlap**: corrects overlap + weak-instrument + winner's-curse via cross-trait
  LDSC. Needs genome-wide sumstats + LD scores.
The overlap biases a strong instrument (LDL) modestly toward the observational
estimate; report the caveat explicitly.

## Files
- `multiancestry_mr.R` — the MR engine (self-contained; streams the big files).
- `run_multiancestry_mr.sh` — driver looping AFR/HIS × LDL/HDL/logTG.
