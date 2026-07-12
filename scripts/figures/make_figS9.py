#!/usr/bin/env python3
"""
Supplementary Figure S9 — Cross-ancestry instrument transferability.
Per-SNP UK Biobank lipid-instrument effects vs GLGC African and Hispanic/admixed
effects, for LDL-C and triglycerides (aligned to the UKB effect allele). High
correlation and sign concordance show the instruments transfer across ancestries;
consistent with Supplementary Table S25.

Outputs figS9_xancestry_scatter.{pdf,png}. Data: xanc_LDL.csv, xanc_TG.csv
(UKB instruments -> GLGC AFR/HIS, reproducing the Table S25 computation).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")


def load(f):
    return list(csv.DictReader(open(os.path.join(DATA, f))))


PANELS = [("xanc_LDL.csv", "afr", "LDL-C", "African", LIPID["LDL-C"]),
          ("xanc_LDL.csv", "his", "LDL-C", "Hispanic/admixed", LIPID["LDL-C"]),
          ("xanc_TG.csv", "afr", "Triglycerides", "African", LIPID["TG"]),
          ("xanc_TG.csv", "his", "Triglycerides", "Hispanic/admixed", LIPID["TG"])]

fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.6))
for ax, (f, anc, tlab, anclab, col) in zip(axes.ravel(), PANELS):
    rows = load(f)
    x = np.array([float(r["ukb_beta"]) for r in rows])
    y = np.array([float(r[anc + "_beta"]) if r[anc + "_beta"] not in ("", "NA") else np.nan for r in rows])
    ok = ~np.isnan(y); x, y = x[ok], y[ok]
    r = np.corrcoef(x, y)[0, 1]; sc = np.mean(np.sign(x) == np.sign(y))
    lim = max(np.abs(x).max(), np.abs(y).max()) * 1.08
    ax.axline((0, 0), slope=1, color=GREY, ls=(0, (4, 3)), lw=0.9, zorder=1)
    ax.axhline(0, color=GREY, lw=0.4, zorder=0); ax.axvline(0, color=GREY, lw=0.4, zorder=0)
    ax.scatter(x, y, s=13, color=col, alpha=0.55, edgecolor="none", zorder=2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel(f"UK Biobank {tlab} effect", fontsize=7.6)
    ax.set_ylabel(f"{anclab} effect", fontsize=7.6)
    ax.text(0.05, 0.95, f"$r$ = {r:.2f}\nsign conc. = {sc:.0%}\n$n$ = {len(x)}",
            transform=ax.transAxes, va="top", ha="left", fontsize=7.4, color=INK)
    despine(ax); ax.tick_params(labelsize=6.6)
    ax.set_title(f"{tlab}  ·  {anclab}", fontsize=8.6, fontweight="bold", pad=4)
fig.text(0.5, 1.005, "Lipid instruments transfer across ancestries (per-SNP effects)",
         ha="center", fontsize=10.6, fontweight="bold")
fig.text(0.5, 0.972, "UK Biobank instrument effects vs GLGC African and Hispanic/admixed effects "
         "(aligned to the UKB effect allele)", ha="center", fontsize=7.4, color=SUBTLE)
fig.text(0.015, 0.985, "S9", fontsize=11, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.955])
save(fig, "figS9_xancestry_scatter")
