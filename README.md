# Convergent common- and rare-variant genetics of the lipid–CAD axis (UK Biobank)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21285448.svg)](https://doi.org/10.5281/zenodo.21285448)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Analysis code for **"Convergent common- and rare-variant genetics maps lipid pathways
and candidate targets for coronary artery disease"** (Pandey & Narasimhan).

Archived on Zenodo: **https://doi.org/10.5281/zenodo.21285448** (concept DOI, always
resolves to the latest version).

This repository contains **code only**. No UK Biobank or individual-level data are
included or tracked (see `.gitignore`). UK Biobank data are available to approved
researchers via the [UK Biobank Access Management System](https://www.ukbiobank.ac.uk);
this work was performed under **UK Biobank application 65439**.

## Overview

Six complementary human-genetic analyses of the lipid–coronary artery disease (CAD)
axis, integrated in a single harmonized cohort of 402,200 UK Biobank participants and
469,835 sequenced exomes, with every conclusion tested outside the discovery sample.
Multivariable MR isolates LDL-C and triglycerides as the independent causal exposures
(HDL-C lacking robust target-specific support; LDL-C and ApoB genetically
inseparable); colocalization and cis-eQTL analysis separate regulatory from
coding-driven effector genes; rare-variant burden recovers canonical lipid genes; and
a discovery-only convergence score recovers established drug targets and nominates
candidates (e.g. *PDE3B*). Replication spans FinnGen, CARDIoGRAMplusC4D and, within
ancestry, African- and Hispanic/admixed-American participants of the Million Veteran
Program.

## Pipeline

Run on the TACC Lonestar6 cluster.

| Analysis | Method | Tool | Script |
|---|---|---|---|
| GWAS | Whole-genome regression | REGENIE v4 | `agent1_genetics/slurm/02–04_*` |
| Rare-variant burden | Gene-based collapsing (exome) | REGENIE | `agent1_genetics/slurm/03b_burden_test.slurm`, `prepare_burden_masks.*` |
| Genetic correlation | LD-score regression | LDSC | `agent1_genetics/slurm/06_ldsc_rg.slurm` |
| Causal inference | Multivariable MR (+ conditional F) | MendelianRandomization, MVMR | `agent1_genetics/R/run_mvmr.R`, `replication/conditional_f_mvmr.R` |
| Colocalization | Bayesian coloc | coloc | `agent1_genetics/R/run_coloc.R` |
| Fine-mapping | Sum of Single Effects | susieR | `agent1_genetics/R/run_finemapping.R` |
| Gene prioritization | Convergence score (common+rare) | in-house | `agent1_genetics/analysis/integrated_convergence_score.py` |
| Allelic series / burden QC | Consequence-stratified burden | in-house | `agent1_genetics/analysis/allelic_series.py`, `burden_diagnostics.py` |

### External replication and sensitivity

| Analysis | Outcome dataset | Script |
|---|---|---|
| Two-sample MVMR | FinnGen r11 CAD (`I9_CHD`) | `replication/finngen_two_sample_mvmr.R` |
| Colocalization replication | FinnGen CAD | `replication/finngen_colocalization.sh` |
| *cis*-LPA restricted MR | FinnGen CAD | `replication/cis_lpa_mr.sh` |
| Statin-naive / covariate sensitivity | UK Biobank | `analysis/statin_naive_sensitivity.sh` |

> The CARDIoGRAMplusC4D two-sample MR, the cross-ancestry GLGC→Million Veteran Program
> MR (with MR-PRESSO outlier correction), the MR-BMA analysis, and the cis-eQTL
> (GTEx / eQTL Catalogue) colocalization are described in the paper; their scripts are
> being added to `replication/` and `analysis/` in the next release.

## Figures

`scripts/figures/` regenerates every main, Central Illustration, and supplementary
display item from the paper's aggregate summary tables (`figstyle.py` holds the shared
style). See `scripts/figures/README.md`. Code only; figure inputs are the
Supplementary Tables archived on Zenodo.

## External datasets (all public)

| Dataset | Use | Accession / source |
|---|---|---|
| FinnGen release 11 (`I9_CHD`) | CAD replication (MR, coloc) | https://www.finngen.fi/en |
| CARDIoGRAMplusC4D | CAD replication (MR) | GWAS Catalog `GCST003116` |
| Global Lipids Genetics Consortium 2021 | Ancestry-stratified lipid instruments | https://csg.sph.umich.edu/willer/public/glgc-lipids2021/ |
| Million Veteran Program | Ancestry-stratified CAD outcome | dbGaP `phs002453` |
| Multi-ancestry exome (Koyama et al.) | Rare-variant concordance | dbGaP `phs001672` |
| GTEx v8 / eQTL Catalogue | cis-eQTL colocalization | https://www.gtexportal.org, https://www.ebi.ac.uk/eqtl/ |
| Open Targets / ChEMBL | Drug-target annotation | https://www.opentargets.org, https://www.ebi.ac.uk/chembl/ |

## Reproducing

1. Obtain UK Biobank access (application) and stage genotype/exome/phenotype data.
2. Edit paths in `config/tacc_paths.sh` to your environment.
3. Create the environment: `conda env create -f envs/environment.yml`.
4. Run: `bash scripts/agent1_genetics/run_all.sh` (submits the SLURM chain).
5. Regenerate figures from summary tables: see `scripts/figures/README.md`.

## Key software
REGENIE v4 · PLINK 2.0 · LDSC · coloc · susieR · MendelianRandomization · MVMR ·
MR-PRESSO · MR-BMA · matplotlib

## Citation
If you use this code, please cite the paper (link to be added upon publication),
UK Biobank (Bycroft et al., *Nature* 2018), FinnGen (Kurki et al., *Nature* 2023),
and the other external resources listed above.

## License
MIT (see `LICENSE`).
