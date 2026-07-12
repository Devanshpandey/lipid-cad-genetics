# Figure generation

Scripts that regenerate every main, Central Illustration, and supplementary display
item in the paper. Shared house style lives in `figstyle.py` (colors, shapes,
`panel()`/`save()`/`despine()` helpers).

These scripts are **code only**. They read small, aggregate summary tables (MR
estimates, colocalization posteriors, gene scores, KM curves, per-SNP instrument
effects) from a `data/` subdirectory, which is intentionally **not** tracked here.
No UK Biobank or individual-level data are used. To reproduce the figures, place the
paper's Supplementary Tables (archived on Zenodo, DOI
[10.5281/zenodo.21285448](https://doi.org/10.5281/zenodo.21285448)) into
`scripts/figures/data/` under the expected filenames, then run e.g.:

```bash
cd scripts/figures
python make_fig5_replication.py        # cross-cohort / cross-ancestry MR forest
python make_central_illustration.py    # 4-stage evidence funnel (Figure 1)
```

Requirements: Python 3, `matplotlib`, `numpy` (see `../../envs/environment.yml`).
Figures are written as paired `.pdf` (vector) and `.png` (300 dpi).
