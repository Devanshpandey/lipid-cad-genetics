#!/usr/bin/env python3
"""
Supplementary Figure S2 — MVMR-Egger intercepts.
Directional-pleiotropy intercept (point + 95% CI) for each outcome, ordered by
coronary relevance; intervals overlap zero for all outcomes except a nominal MI
signal (P=0.048, not multiplicity-corrected). SE derived from the reported P.

Outputs figS2_mvmr_pleiotropy.{pdf,png}. Data: TableS3 (mvegger_intercept, P).
"""
import os, sys, csv
from scipy.stats import norm
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); SUP = os.path.join(HERE, "..", "supplementary")
r = list(csv.DictReader(open(os.path.join(SUP, "TableS3_mvmr.csv"))))
seen = {}
for x in r:
    seen.setdefault(x["outcome"], (float(x["mvegger_intercept"]), float(x["mvegger_intercept_p"])))
NICE = {"CAD": "CAD", "MI": "MI", "REVASC": "Revascularization", "MACE": "MACE",
        "STROKE": "Stroke", "HF": "Heart failure", "CV_DEATH": "CV death"}
order = ["CAD", "MI", "REVASC", "MACE", "STROKE", "HF", "CV_DEATH"]

fig, ax = plt.subplots(figsize=(6.6, 4.4))
for i, o in enumerate(order):
    ic, p = seen[o]; y = len(order) - 1 - i
    se = abs(ic) / norm.ppf(1 - p / 2) if p < 1 else abs(ic)
    lo, hi = ic - 1.96 * se, ic + 1.96 * se
    nominal = p < 0.05
    col = CB["red"] if nominal else (LIPID["outcome"] if o in ("CAD", "MI", "REVASC", "MACE") else GREY)
    ax.plot([lo, hi], [y, y], color=col, lw=2.0, solid_capstyle="round", zorder=2)
    ax.scatter(ic, y, s=44, color=col, edgecolor="white", lw=0.7, zorder=3)
    ax.text(0.0135, y, f"P={p:.3f}" + ("  †" if nominal else ""), transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=6.6, color=col if nominal else SUBTLE)
ax.axvline(0, color=INK, ls=(0, (3, 3)), lw=0.9)
ax.set_yticks(range(len(order))); ax.set_yticklabels([NICE[o] for o in order[::-1]], fontsize=8)
ax.set_xlabel("MVMR-Egger intercept (95% CI)")
ax.set_xlim(-0.009, 0.013)
despine(ax)
ax.text(-0.10, 1.10, "S2", transform=ax.transAxes, fontsize=11, fontweight="bold", ha="left", color=INK)
ax.text(0.5, 1.115, "MVMR-Egger intercepts show little directional pleiotropy",
        transform=ax.transAxes, fontsize=9.6, fontweight="bold", va="top", ha="center", color=INK)
ax.text(0.5, 1.04, "† only MI reaches nominal (uncorrected) significance",
        transform=ax.transAxes, fontsize=7.2, va="top", ha="center", color=SUBTLE)
save(fig, "figS2_mvmr_pleiotropy")
