# Changelog

## v2.0.0 (24 November 2026)

Brings the repository up to the current paper. Earlier files are kept unless named here.

**Added**
- Genetic comparison with Lp(a)HORIZON, its post hoc checks, the trial-proxy endpoint, the trial forecasts and the amplification
  test (`A1/lpa_horizon`, `A1/lpa_susceptibility`), with their protocol documents and history (`protocols/`) and summary-level
  outputs (`results_summary/`).
- The polygenic score the comparison adjusts for, and the libraries the amplification test uses (`scripts/agent6_resilience`).
- Analyses added in review: ApoB with Lp(a), weak-instrument-robust sets, selection into statin use, the design of Tudball et al.,
  multi-signal colocalization (`A1/nature_round`).
- Statin ascertainment correction and statin-interaction analyses; multivariable MR after the distance-pruning and covariance
  corrections; external, stroke-subtype and drug-target MR; colocalization and fine-mapping runs (`A1/coloc_runs`); *cis*-eQTL and
  *cis*-pQTL colocalization; the Lp(a) GWAS that keeps
  out-of-range values; type 2 diabetes and TG:HDL analyses of *PDE3B*; the recomputed gene score.
- `TABLE_SOURCES.tsv` (which output each published table carries) and `CODE_MD5.tsv`.
- Figure code of the current paper (`scripts/figures_nature`): every main and Supplementary Figure, the map of its inputs
  (`INPUTS.tsv`) and `make_layout.py`, which redraws them from this repository and the paper's Supplementary Data.
- `scripts/run_selftests.sh`, `THIRD_PARTY_NOTICES.md`, and `docs/study_design.png` (panel a of Fig. 1, shown in the README;
  `scripts/figures_nature/make_readme_figure.py` draws it).

**Changed**
- `A1/R/prep_phenotypes.R` (statin use read from 6153 and 6177), `A1/R/run_mvmr.R` (instrument ranking, pruning, covariance),
  `A1/R/run_coloc.R` and `A1/slurm/07_coloc.slurm` (one set of sentinels per exposure), `config/tacc_paths.sh` (genotype paths
  exported to child processes), `README.md`, `CITATION.cff`, `.gitignore`.

`A1` is `scripts/agent1_genetics`.

## v1.2.0 (9 July 2026) and earlier

Releases for an earlier version of the paper. Commits after v1.2.0 and before v2.0.0 (12 to 20 July 2026) were not archived.
