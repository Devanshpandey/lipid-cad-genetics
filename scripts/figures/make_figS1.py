#!/usr/bin/env python3
"""
Supplementary Figure S1 — Exome-wide rare-variant burden recovers trait-specific
gene architecture. One lollipop panel per trait (LDL-C, triglycerides, HDL-C,
Lp(a), CAD, MI): top burden genes by -log10P, coloured by effect direction, with
the exome-wide-significance threshold. Values are true -log10P (per-panel axis,
capped where a single gene dominates, with the true value labelled).

Outputs figS1_burden_validation.{pdf,png}. Data: TableS6 (burden per trait).
"""
import os, sys, csv
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); SUP = os.path.join(HERE, "..", "supplementary")
rows = list(csv.DictReader(open(os.path.join(SUP, "TableS6_rvas_siggenes.csv"))))
bytrait = defaultdict(list)
for r in rows:
    bytrait[r["trait"]].append((r["gene"], float(r["log10P"]), float(r["beta"])))
PANELS = [("LDL_C", "LDL-C"), ("TRIGLY", "Triglycerides"), ("HDL_C", "HDL-C"),
          ("LPA", "Lp(a)"), ("CAD", "CAD"), ("MI", "MI")]
DOWN, UP = LIPID["HDL-C"], LIPID["LDL-C"]   # lowers = blue, raises = red
THR = 5.6   # -log10(2.5e-6)
CAP = 90.0  # display cap; true value labelled when exceeded

fig, axes = plt.subplots(2, 3, figsize=(12.4, 6.6))
for ax, (tcode, tlab) in zip(axes.ravel(), PANELS):
    genes = sorted(bytrait.get(tcode, []), key=lambda g: -g[1])[:8]
    genes = genes[::-1]  # smallest at bottom
    for i, (g, lp, b) in enumerate(genes):
        disp = min(lp, CAP); col = UP if b > 0 else DOWN
        ax.plot([0, disp], [i, i], color=col, lw=1.6, alpha=0.6, zorder=2)
        ax.scatter(disp, i, s=34, color=col, edgecolor="white", lw=0.5, zorder=3)
        lab = f"{lp:.0f}" + ("*" if lp > CAP else "")
        ax.text(disp + CAP * 0.015, i, lab, va="center", ha="left", fontsize=5.6, color=SUBTLE)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([g[0] for g in genes], fontsize=6.6, fontstyle="italic")
    xmax = max(CAP * 1.12, 12)
    ax.axvline(THR, color=GREY, ls=(0, (3, 3)), lw=0.9)
    ax.set_xlim(0, xmax); ax.set_ylim(-0.6, max(len(genes) - 0.4, 1))
    ax.set_xlabel("$-\\log_{10}P$", fontsize=7.2); despine(ax); ax.tick_params(labelsize=6.4)
    ax.set_title(tlab, fontsize=8.6, fontweight="bold", pad=4)
    if tcode in ("CAD", "MI"):
        ax.text(0.96, 0.1, "outcome burden\n(few genes)", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5.8, color=SUBTLE)
fig.text(0.5, 1.005, "Exome-wide rare-variant burden recovers trait-specific gene architecture",
         ha="center", fontsize=10.5, fontweight="bold")
fig.text(0.5, 0.965, "top genes per trait; blue lowers, red raises the trait; dashed = exome-wide "
         "significance (P<2.5×10$^{-6}$); * = value exceeds display cap",
         ha="center", fontsize=7.2, color=SUBTLE)
fig.text(0.012, 0.985, "S1", fontsize=11, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.95])
save(fig, "figS1_burden_validation")
