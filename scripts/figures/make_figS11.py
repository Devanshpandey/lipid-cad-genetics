#!/usr/bin/env python3
"""
Supplementary Figure S11 — PDE3B variant-level burden architecture.
  (a) Leave-one-variant-out: the pLoF-burden effect on triglycerides is stable
      when any single qualifying variant is removed (range -0.34 to -0.39 vs the
      full -0.378), so it is not driven by one recurrent variant.
  (b) Per-variant single-variant triglyceride effects for the qualifying variants
      with >=10 carriers (95% CI from P), by carrier count.

Outputs figS11_pde3b_variant.{pdf,png}. Data: pde3b_lovo.csv, pde3b_burden_full.csv.
"""
import os, sys, csv, math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
TEAL = LIPID["rare"]
rows = list(csv.DictReader(open(os.path.join(DATA, "pde3b_lovo.csv"))))
FULL = -0.378
def f(v):
    return float(v) if v not in ("", "NA") else float("nan")


for r in rows:
    r["n_carriers"] = int(r["n_carriers"]); r["drop_beta"] = f(r["drop_beta"])
    r["single_beta"] = f(r["single_beta"]); r["single_P"] = f(r["single_P"])
rows.sort(key=lambda r: -r["n_carriers"])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.7), gridspec_kw={"width_ratios": [1.1, 1.0], "wspace": 0.34})

# ---------------------------------------------------- (a) leave-one-out
db = [r["drop_beta"] for r in rows]
axA.axhspan(min(db), max(db), color="#EAF3F1", zorder=0)
axA.axhline(FULL, color=INK, ls=(0, (4, 3)), lw=1.1, zorder=2, label=f"full burden ({FULL:.2f})")
axA.scatter(range(len(rows)), db, s=18, color=TEAL, alpha=0.7, edgecolor="none", zorder=3)
axA.scatter(0, db[0], s=40, color=CB["red"], edgecolor="white", lw=0.5, zorder=4)
axA.annotate(f"drop most common variant\n(780 carriers): {db[0]:.2f}", (0, db[0]),
             textcoords="offset points", xytext=(16, -2), fontsize=6.4, color=CB["red"], va="center")
axA.set_xlabel("qualifying pLoF variant removed (ranked by carrier count)", fontsize=7.6)
axA.set_ylabel("triglyceride burden effect (SD)")
axA.set_ylim(-0.42, -0.30); axA.legend(loc="lower right", fontsize=7)
despine(axA)
axA.text(-0.02, 1.13, "a", transform=axA.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=INK)
axA.text(0.5, 1.125, "Leave-one-variant-out is stable", transform=axA.transAxes,
         fontsize=9.4, fontweight="bold", va="top", ha="center", color=INK)
axA.text(0.5, 1.045, "no single variant drives the burden signal", transform=axA.transAxes,
         fontsize=7.2, va="top", ha="center", color=SUBTLE)

# ---------------------------------------------------- (b) per-variant effects
big = [r for r in rows if r["n_carriers"] >= 10 and not math.isnan(r["single_beta"])
       and not math.isnan(r["single_P"])]
big.sort(key=lambda r: r["n_carriers"])
for i, r in enumerate(big):
    b, p, n = r["single_beta"], r["single_P"], r["n_carriers"]
    se = abs(b) / norm.ppf(1 - p / 2) if p < 1 else abs(b)
    lo, hi = b - 1.96 * se, b + 1.96 * se
    sig = p < 0.05
    axB.plot([lo, hi], [i, i], color=TEAL if b < 0 else CB["red"], lw=1.8, solid_capstyle="round", zorder=2)
    axB.scatter(b, i, s=34, facecolor=(TEAL if b < 0 else CB["red"]) if sig else "white",
                edgecolor=TEAL if b < 0 else CB["red"], lw=1.2, zorder=3)
    axB.text(1.01, i, f"n={n}", transform=axB.get_yaxis_transform(), va="center", ha="left",
             fontsize=6.0, color=SUBTLE)
axB.axvline(0, color=INK, ls=(0, (3, 3)), lw=0.9)
axB.axvline(FULL, color=GREY, ls=(0, (2, 2)), lw=0.9)
axB.text(FULL, len(big) - 0.3, "full burden", fontsize=6.0, color=GREY, ha="center", rotation=90, va="top")
axB.set_yticks(range(len(big))); axB.set_yticklabels([r["variant"].split(":")[1] for r in big], fontsize=6.2)
axB.set_ylabel("variant (chr11 position)", fontsize=7.6)
axB.set_xlabel("single-variant triglyceride effect (SD, 95% CI)", fontsize=7.6)
axB.set_ylim(-0.7, len(big) - 0.3); despine(axB)
axB.text(-0.02, 1.13, "b", transform=axB.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=INK)
axB.text(0.5, 1.125, "Individual qualifying variants", transform=axB.transAxes,
         fontsize=9.4, fontweight="bold", va="top", ha="center", color=INK)
axB.text(0.5, 1.045, "variants with $\\geq$10 carriers; filled = P<0.05", transform=axB.transAxes,
         fontsize=7.2, va="top", ha="center", color=SUBTLE)
save(fig, "figS11_pde3b_variant")
