# Figure code of the current paper

Every main and Supplementary Figure of the current version of the paper (revision of 24 November 2026) was drawn by
the scripts in `figure_src/`, copied here unchanged. All figures share one theme, `figure_src/r4/theme.py`. The scripts
read summary tables only and fit nothing.

| Figure | Script |
|---|---|
| Fig. 1, genetically predicted Lp(a) before and after disease onset | `figure_src/r5/fig1.py` |
| Fig. 2, translation into a trial effect, and forecasts | `figure_src/r5/fig2.py` |
| Fig. 3, multivariable MR and replication | `figure_src/r5/fig3.py` |
| Fig. 4, vascular territory and statin use | `figure_src/r4/fig3.py` |
| LDL-C and ApoB as one axis; within-ancestry replication (Supplementary Note 6) | `figure_src/r5/si_moved.py` |
| Dose-response across predicted Lp(a) | `figure_src/r4/fig2.py` (it also draws a figure of an earlier version) |
| eQTL colocalization of LDL-C and ApoB | `figure_src/r4/figS_eqtl_coloc.py` |
| Other Supplementary Figures | `figure_src/r4/si/`, one script per figure |

The r5 scripts import the r4 panel functions (`figure_src/r4/fig1.py`, `fig2.py`, `fig4.py`), so both folders are
needed. `figure_src/r5/build_figures.sh` runs everything in the manuscript's own folder.

## Running it from this repository

The scripts were run inside the manuscript's folder and read their inputs by relative path. `INPUTS.tsv` lists every
input with its md5, the scripts that read it, and where a reader gets it:
- **this repository:** the comparison with Lp(a)HORIZON, in `results_summary/lpa_horizon/`; in `figure_data/`, the
  summary tables that are not among the Supplementary Data; and, for Fig. 4, the three printed lines of the Supplementary
  Information that the script reads by pattern (`figure_data/supplement_lines_for_fig4.tex`);
- **the paper's Supplementary Data**, under its published file name.

One input is not included: `figure_data/si/pde3b_lovo.csv`, the leave-one-variant-out burden estimates behind the
*PDE3B* variant Supplementary Figure, which list per-variant carrier counts below five.

`make_layout.py` builds that layout in a folder of your choice and checks every md5. With `--draw` it runs the scripts:

```
python scripts/figures_nature/make_layout.py figures_work --supplementary-data path/to/supplementary_data --draw
```

With this repository and the paper's Supplementary Data, the four main figures and 19 of the 20 Supplementary Figures
are drawn; the script names the one it cannot draw. Drawn this way, each image and each source-data file was
byte-identical to the one in the paper's version of 24 November 2026.

`make_readme_figure.py` draws `docs/study_design.png`, the README's figure, from a laid-out folder with the paper's own panel
code. `regen_cismr_artifacts.py` writes the files of the drug-target *cis*-MR table; it draws nothing. `scripts/figures/` holds
the figure code of earlier versions of the paper.
