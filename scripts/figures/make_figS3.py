#!/usr/bin/env python3
"""
Supplementary Figure S3 — MVMR-IVW vs MVMR-Egger effect concordance.
Per outcome x exposure, the multivariable-IVW causal estimate (x) versus the
pleiotropy-robust MVMR-Egger estimate (y), coloured by exposure, with the identity
line and a fitted line. High concordance (r, sign agreement) indicates the causal
estimates are not driven by the no-pleiotropy assumption.

Outputs figS3_mvmr_robustness.{pdf,png}. Data: TableS3 (mvivw_b, mvegger_b).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); SUP = os.path.join(HERE, "..", "supplementary")
r = list(csv.DictReader(open(os.path.join(SUP, "TableS3_mvmr.csv"))))
EXP = {"LDL_C": ("LDL-C", LIPID["LDL-C"]), "TRIGLY": ("Triglycerides", LIPID["TG"]),
       "HDL_C": ("HDL-C", LIPID["HDL-C"]), "LPA": ("Lp(a)", LIPID["Lp(a)"])}
ivw = np.array([float(x["mvivw_b"]) for x in r]); egg = np.array([float(x["mvegger_b"]) for x in r])
cor = np.corrcoef(ivw, egg)[0, 1]
conc = np.mean(np.sign(ivw) == np.sign(egg))
slope, intercept = np.polyfit(ivw, egg, 1)

fig, ax = plt.subplots(figsize=(5.9, 5.4))
lim = max(np.abs(ivw).max(), np.abs(egg).max()) * 1.15
ax.axline((0, 0), slope=1, color=GREY, ls=(0, (4, 3)), lw=1.0, zorder=1, label="identity")
xs = np.array([-lim, lim])
ax.plot(xs, slope * xs + intercept, color=INK, lw=1.0, zorder=2)
for x in r:
    lab, col = EXP[x["exposure"]]
    ax.scatter(float(x["mvivw_b"]), float(x["mvegger_b"]), s=34, color=col,
               edgecolor="white", lw=0.5, zorder=3)
ax.axhline(0, color=GREY, lw=0.5, zorder=0); ax.axvline(0, color=GREY, lw=0.5, zorder=0)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel("MVMR-IVW effect (log OR/SD)"); ax.set_ylabel("MVMR-Egger effect (log OR/SD)")
ax.text(0.04, 0.96, f"$r$ = {cor:.3f}\nslope = {slope:.2f}\nsign concordance = {conc:.0%}\n$n$ = {len(r)} outcome$\\times$exposure",
        transform=ax.transAxes, va="top", ha="left", fontsize=7.4, color=INK)
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=7, label=l)
           for l, c in EXP.values()]
ax.legend(handles=handles, loc="lower right", fontsize=7, title="exposure", title_fontsize=7)
despine(ax)
ax.text(-0.12, 1.08, "S3", transform=ax.transAxes, fontsize=11, fontweight="bold", ha="left", color=INK)
ax.text(0.5, 1.09, "Effect estimates are concordant between MVMR-IVW and MVMR-Egger",
        transform=ax.transAxes, fontsize=8.8, fontweight="bold", va="top", ha="center", color=INK)
save(fig, "figS3_mvmr_robustness")
