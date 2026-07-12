#!/usr/bin/env python3
"""
Supplementary Figure S5 — Colocalization and fine-mapping robustness.
  (a) FinnGen external PP.H4 under the default and a stringent (p12=1e-6) prior,
      as paired dots per locus; classified robust / prior-sensitive / fails.
  (b) LDLR credible-set count under increasing SuSiE L (10/15/20): the sets keep
      splitting, quantifying allelic heterogeneity (no convergence at L=20).

Outputs figS5_coloc_finemap.{pdf,png}. Data: TableS9 (FinnGen coloc priors),
susie_L_sensitivity.csv (LDLR L-dependence).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
TEAL, AMBER, RED = LIPID["rare"], "#E08A2B", LIPID["LDL-C"]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.6), gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.33})

# ---------------------------------------------------- (a) paired prior dots
loci = ["PCSK9", "SORT1/CELSR2", "APOE/APOC1", "LDLR", "TRIB1"]
default = [1.000, 0.998, 1.000, 0.957, 0.612]
stringent = [1.000, 0.976, 0.999, 0.691, 0.136]
cls = ["robust", "robust", "robust", "prior-sensitive", "fails"]
clscol = {"robust": TEAL, "prior-sensitive": AMBER, "fails": GREY}
for i, loc in enumerate(loci):
    y = len(loci) - 1 - i; c = clscol[cls[i]]
    axA.plot([default[i], stringent[i]], [y, y], color=c, lw=1.6, zorder=2, alpha=0.7)
    axA.scatter(default[i], y, s=60, facecolor="white", edgecolor=c, lw=1.6, marker="o", zorder=3)
    axA.scatter(stringent[i], y, s=52, color=c, edgecolor="white", lw=0.7, marker="D", zorder=3)
    axA.text(1.06, y, cls[i], ha="left", va="center", fontsize=6.8, color=c, fontweight="bold")
axA.axvline(0.8, color=RED, ls=(0, (3, 3)), lw=1.1)
axA.text(0.8, len(loci) - 0.4, "PP.H4 = 0.8", fontsize=6.2, color=RED, ha="center", va="bottom")
axA.set_yticks(range(len(loci))); axA.set_yticklabels(loci[::-1], fontsize=8, fontstyle="italic")
axA.set_xlim(0, 1.32); axA.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); axA.set_ylim(-0.6, len(loci) - 0.15)
axA.set_xlabel("FinnGen external PP.H4"); despine(axA)
axA.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="w", markeredgecolor=INK, markersize=7, label="default prior"),
                    Line2D([0], [0], marker="D", color="w", markerfacecolor=INK, markersize=7, label="stringent prior")],
           loc="lower left", fontsize=6.8)
axA.text(-0.02, 1.13, "a", transform=axA.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=INK)
axA.text(0.5, 1.125, "External replication vs shared-variant prior", transform=axA.transAxes,
         fontsize=9.4, fontweight="bold", va="top", ha="center", color=INK)
axA.text(0.5, 1.045, "three of five loci robust to both priors", transform=axA.transAxes,
         fontsize=7.2, va="top", ha="center", color=SUBTLE)

# ---------------------------------------------------- (b) LDLR CS count vs L
rows = list(csv.DictReader(open(os.path.join(DATA, "susie_L_sensitivity.csv"))))
Ls = [10, 15, 20]
ncs = [max(int(r["n_cs"]) for r in rows if int(r["L"]) == L) for L in Ls]
axB.bar(range(3), ncs, color=[CB["sky"], CB["blue"], LIPID["outcome"]], width=0.6, zorder=2)
for xi, (L, n) in enumerate(zip(Ls, ncs)):
    axB.text(xi, n + 0.3, str(n), ha="center", fontsize=8, color=INK, fontweight="bold")
axB.plot(range(3), Ls, color=GREY, ls=(0, (3, 3)), lw=1.0, marker="o", ms=4, zorder=3)
axB.text(0.55, 16.2, "L cap", fontsize=6.4, color=GREY, ha="left", va="center")
axB.set_xticks(range(3)); axB.set_xticklabels([f"L={L}" for L in Ls], fontsize=8)
axB.set_ylim(0, 23); axB.set_ylabel("LDLR credible sets recovered")
despine(axB)
axB.text(-0.02, 1.13, "b", transform=axB.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=INK)
axB.text(0.5, 1.125, "LDLR allelic heterogeneity", transform=axB.transAxes,
         fontsize=9.4, fontweight="bold", va="top", ha="center", color=INK)
axB.text(0.5, 1.045, "credible sets keep splitting; no convergence at L=20",
         transform=axB.transAxes, fontsize=7.2, va="top", ha="center", color=SUBTLE)

save(fig, "figS5_coloc_finemap")
