#!/usr/bin/env python3
"""
Supplementary Figure S8 — Robustness and complementarity of the convergence score.
  (a) Coverage vs Open Targets L2G: the score ranks all 27 canonical targets
      genome-wide; L2G scores only the 5 with common-variant signals (22 act
      through rare coding). Inset: integrated vs L2G among jointly-scored genes.
  (b) Layer ablation: AUROC for canonical recovery, with nested-LRT P vs burden.
  (c) PDE3B genome-wide rank across 108 weight/cap specifications (all top 0.1%).
  (d) Recovery AUROC by positive-label class.

Outputs figS8_score_robustness.{pdf,png}. Data: TableS19 (L2G), S20 (ablation),
S22 (weight sensitivity), S23 (stratified AUROC).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); SUP = os.path.join(HERE, "..", "supplementary")
TEAL, BURG, AMBER = LIPID["rare"], LIPID["outcome"], "#E08A2B"


def hd(ax, letter, title, sub=None):
    ax.text(-0.02, 1.17, letter, transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="top", ha="right", color=INK)
    ax.text(0.5, 1.165, title, transform=ax.transAxes, fontsize=9.3, fontweight="bold",
            va="top", ha="center", color=INK)
    if sub:
        ax.text(0.5, 1.075, sub, transform=ax.transAxes, fontsize=7.0, va="top",
                ha="center", color=SUBTLE)


fig = plt.figure(figsize=(11.2, 8.4))
gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.34)
axA, axB, axC, axD = [fig.add_subplot(gs[i]) for i in range(4)]

# ------------------------------------------------------- (a) L2G coverage
d = list(csv.DictReader(open(os.path.join(SUP, "TableS19_l2g_benchmark.csv"))))
integ = np.array([float(r["integrated"]) for r in d]); l2g = np.array([float(r["l2g"]) for r in d])
canon = np.array([r["canonical"] == "1" for r in d])
rho = spearmanr(integ, l2g).correlation
axA.bar([0, 1], [27, 5], color=[TEAL, GREY], width=0.58, zorder=2)
axA.text(0, 27.6, "27", ha="center", fontsize=8.5, fontweight="bold", color=INK)
axA.text(1, 5.6, "5", ha="center", fontsize=8.5, fontweight="bold", color=INK)
axA.set_xticks([0, 1]); axA.set_xticklabels(["convergence\nscore", "Open Targets\nL2G"], fontsize=7.6)
axA.set_ylim(0, 31); axA.set_xlim(-0.6, 1.6); axA.set_ylabel("canonical targets scored (of 27)")
axA.text(0.02, 0.52, "22 / 27 canonical targets\nact through rare coding\nvariation, outside the\nL2G-scored set",
         transform=axA.transAxes, ha="left", va="top", fontsize=6.6, color=BURG, fontweight="bold")
despine(axA)
# inset scatter among jointly-scored genes (empty upper-right, above the short L2G bar)
ins = axA.inset_axes([0.56, 0.50, 0.42, 0.44])
ins.scatter(l2g[~canon], integ[~canon], s=7, color=GREY, alpha=0.5, edgecolor="none")
ins.scatter(l2g[canon], integ[canon], s=20, color=AMBER, edgecolor="black", lw=0.3)
ins.set_xlabel("L2G score", fontsize=6); ins.set_ylabel("conv. score", fontsize=6)
ins.tick_params(labelsize=5.2); ins.set_title(f"jointly scored: $\\rho$={rho:.2f}", fontsize=6, pad=2)
despine(ins)
hd(axA, "a", "Orthogonal to common-variant L2G", "captures rare-variant genes L2G cannot score")

# ------------------------------------------------------- (b) layer ablation
ab = list(csv.DictReader(open(os.path.join(SUP, "TableS20_layer_ablation.csv"))))
labs = ["burden", "+coloc", "+outcome", "full"]
au = [float(r["AUROC_canonical_recovery"]) for r in ab]
lrt = [r["nested_LRT_P_vs_burden_only"] for r in ab]
cols = [GREY, TEAL, GREY, TEAL]
axB.bar(range(4), au, color=cols, width=0.62, zorder=2)
for xi, (v, p) in enumerate(zip(au, lrt)):
    axB.text(xi, v + 0.004, f"{v:.2f}", ha="center", fontsize=7, color=INK)
    if p not in ("NA", ""):
        axB.text(xi, 0.862, f"P={float(p):.0e}" if float(p) < 0.01 else f"P={float(p):.2f}",
                 ha="center", fontsize=5.6, color=SUBTLE)
axB.set_xticks(range(4)); axB.set_xticklabels(labs, fontsize=7.4)
axB.set_ylim(0.85, 1.0); axB.set_ylabel("AUROC (canonical recovery)")
axB.text(0.03, 0.97, "gain from colocalization;\noutcome burden adds nothing",
         transform=axB.transAxes, ha="left", va="top", fontsize=6.4, color=TEAL, fontweight="bold")
despine(axB)
hd(axB, "b", "Which layer earns the gain", "nested ablation (LRT P vs burden-only)")

# ------------------------------------------------------- (c) PDE3B rank strip
ws = list(csv.DictReader(open(os.path.join(SUP, "TableS22_score_weight_sensitivity.csv"))))
ranks = np.array([int(r["pde3b_rank"]) for r in ws])
rng = np.random.default_rng(42)
axC.scatter(ranks, rng.uniform(-0.3, 0.3, len(ranks)), s=16, color=CB["purple"], alpha=0.5, edgecolor="none")
axC.boxplot([ranks], vert=False, positions=[0], widths=0.5, showfliers=False,
            medianprops=dict(color=INK, lw=1.4), boxprops=dict(color=INK),
            whiskerprops=dict(color=INK), capprops=dict(color=INK))
axC.axvline(18, color=GREY, ls=(0, (3, 3)), lw=0.9)
axC.text(18, 0.55, "top 0.1%", fontsize=6.2, color=GREY, ha="center")
axC.set_yticks([]); axC.set_xlim(0, max(26, ranks.max() + 3)); axC.set_ylim(-0.7, 0.7)
axC.set_xlabel("PDE3B genome-wide rank")
axC.text(0.97, 0.9, f"top 0.1% in all\n{len(ranks)} specifications", transform=axC.transAxes,
         ha="right", va="top", fontsize=6.8, color=INK)
despine(axC, left=False)
hd(axC, "c", "Nomination robust to weighting", "PDE3B rank across 108 weight/cap specs")

# ------------------------------------------------------- (d) positive classes
st = list(csv.DictReader(open(os.path.join(SUP, "TableS23_stratified_auroc.csv"))))
labs = [s["stratum"].replace("approved/trial drug target", "approved/\ntrial")
        .replace("Mendelian dyslipidemia", "Mendelian").replace("all canonical", "all\ncanonical") for s in st]
vals = [float(s["AUROC"]) for s in st]; ns = [s["n"] for s in st]
axD.bar(range(len(st)), vals, color=[AMBER, CB["blue"], TEAL], width=0.6, zorder=2)
for xi, (v, n) in enumerate(zip(vals, ns)):
    axD.text(xi, v + 0.004, f"{v:.2f}\n(n={n})", ha="center", fontsize=6.6, color=INK)
axD.set_xticks(range(len(st))); axD.set_xticklabels(labs, fontsize=7.0)
axD.set_ylim(0.85, 1.02); axD.set_ylabel("AUROC")
despine(axD)
hd(axD, "d", "Recovery holds across target classes", "positive set analysed by evidence type")

save(fig, "figS8_score_robustness")
