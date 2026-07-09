# Lipid–CAD genetic architecture (UK Biobank)

Analysis code for **"Convergent common-variant, causal, and rare-variant genetic
evidence delineates the lipid–coronary artery disease axis in 402,200 UK Biobank
participants"** (Pandey & Narasimhan).

This repository contains **code only**. No UK Biobank or individual-level data are
included or tracked (see `.gitignore`). UK Biobank data are available to approved
researchers via the [UK Biobank Access Management System](https://www.ukbiobank.ac.uk);
this work was performed under **UK Biobank application 65439**.

## Pipeline overview

Five complementary genetic analyses of the lipid–CAD axis, run on the TACC
Lonestar6 cluster:

| Step | Method | Tool | Script |
|---|---|---|---|
| GWAS | Whole-genome regression | REGENIE v4 | `scripts/agent1_genetics/slurm/02–04_*` |
| Rare-variant burden | Gene-based collapsing (exome) | REGENIE | `scripts/agent1_genetics/slurm/03b_burden_test.slurm`, `prepare_burden_masks.*` |
| Genetic correlation | LD-score regression | LDSC | `scripts/agent1_genetics/slurm/06_ldsc_rg.slurm` |
| Causal inference | Multivariable MR (+ conditional F) | MendelianRandomization, MVMR | `R/run_mvmr.R`, `replication/conditional_f_mvmr.R` |
| **External replication** | **Two-sample MVMR vs FinnGen CAD** | MendelianRandomization | `replication/finngen_two_sample_mvmr.R` |
| Colocalization | Bayesian coloc | coloc | `R/run_coloc.R` |
| Fine-mapping | Sum of Single Effects | susieR | `R/run_finemapping.R` |

Downstream agents (gene prioritization, networks, subtyping) are scaffolded under
`scripts/agent2_genes/`, `agent3_networks/`, `agent4_subtypes/`.

## Reproducing

1. Obtain UK Biobank access (application) and stage genotype/exome/phenotype data.
2. Edit paths in `config/tacc_paths.sh` to your environment.
3. Create the environment: `conda env create -f envs/environment.yml`.
4. Run: `bash scripts/agent1_genetics/run_all.sh` (submits the SLURM chain).

External replication additionally uses public **FinnGen** release-11 summary
statistics (endpoint `I9_CHD`), downloaded from the FinnGen public bucket.

## Key software
REGENIE v4 · PLINK 2.0 · LDSC · coloc · susieR · MendelianRandomization · MVMR

## Citation
If you use this code, please cite the paper (link to be added upon publication)
and UK Biobank (Bycroft et al., *Nature* 2018) and FinnGen (Kurki et al.,
*Nature* 2023).

## License
MIT (see `LICENSE`).
