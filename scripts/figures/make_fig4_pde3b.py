#!/usr/bin/env python3
"""
Figure 5 — PDE3B rare coding variation modifies lipid levels but does not establish
CAD protection.
  (a) effect-based forests: quantitative phenotypes (SD) and binary outcomes (OR),
      each with 95% CI; Bonferroni-significant filled.
  (b) allelic-severity series (pLoF > missense > synonymous negative control).
  (c) CAD association designs, separated: gene-burden proxy vs prospective carrier Cox.

Outputs fig4_pde3b.{pdf,png} (the manuscript's Figure 5 include).
Data: TableS13 (PheWAS effect/se), pde3b_consequence_class, TableS21 (Cox).
"""
import os, sys, csv, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); SUP = os.path.join(HERE, "..", "supplementary")
RED, ORNG, BLUE, PURP = LIPID["LDL-C"], LIPID["TG"], LIPID["HDL-C"], LIPID["Lp(a)"]
TEAL, BURG, NS = LIPID["rare"], LIPID["outcome"], "#9AA3B0"
AXIS = {"TRIGLY": ORNG, "HDL_C": BLUE, "ApoA1": BLUE, "APOA1": BLUE, "ApoB": RED, "APOB": RED,
        "non-HDL-C": RED, "nonHDL_C": RED, "LDL_C": RED, "TOT_CHOL": RED, "LPA": PURP}
NICE = {"TRIGLY": "Triglycerides", "HDL_C": "HDL-C", "APOA1": "ApoA1", "APOB": "ApoB",
        "nonHDL_C": "non-HDL-C", "LDL_C": "LDL-C", "CRP": "CRP", "eGFR": "eGFR",
        "HBA1C": "HbA1c", "TOT_CHOL": "Total chol.", "LPA": "Lp(a)", "STATIN_USE": "Statin use",
        "CAD": "CAD", "MI": "MI", "MACE": "MACE", "REVASC": "Revascularization", "PAD": "PAD",
        "STROKE": "Stroke", "HF": "Heart failure", "AF": "Atrial fibrillation", "CV_DEATH": "CV death"}

ph = {r["pheno"]: r for r in csv.DictReader(open(os.path.join(SUP, "TableS13_pde3b_phewas.csv")))}


def ci(r):
    eff, se = float(r["effect"]), float(r["se"])
    if r["effect_lab"] == "OR":
        return eff, math.exp(math.log(eff) - 1.96 * se), math.exp(math.log(eff) + 1.96 * se)
    return eff, eff - 1.96 * se, eff + 1.96 * se


fig = plt.figure(figsize=(13.0, 7.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.30)
gA = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.62)
axQ = fig.add_subplot(gA[0]); axB = fig.add_subplot(gA[1])
gR = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.52, height_ratios=[1, 1])
axS = fig.add_subplot(gR[0]); axD = fig.add_subplot(gR[1])

# ------------------------------------------------ (a-left) quantitative forest
quant = ["TRIGLY", "HDL_C", "APOA1", "APOB", "nonHDL_C", "LDL_C", "CRP", "eGFR", "HBA1C"]
for i, p in enumerate(quant):
    key = {"APOA1": "APOA1", "nonHDL_C": "nonHDL_C"}.get(p, p)
    r = ph[key]; y = len(quant) - 1 - i
    e, lo, hi = ci(r); sig = float(r["P_bonf"]) < 0.05
    col = AXIS.get(key, NS)
    axQ.plot([lo, hi], [y, y], color=col, lw=1.9, solid_capstyle="round", zorder=2)
    axQ.scatter(e, y, s=42, facecolor=col if sig else "white", edgecolor=col, lw=1.4, zorder=3)
    axQ.set_yticks(range(len(quant)))
axQ.axvline(0, color=INK, ls=(0, (3, 3)), lw=0.9)
axQ.set_yticks(range(len(quant)))
axQ.set_yticklabels([NICE[k] for k in quant[::-1]], fontsize=7.4)
axQ.set_xlabel("effect (SD, 95% CI)", fontsize=7.8); despine(axQ); axQ.tick_params(labelsize=7)
axQ.text(-0.02, 1.17, "a", transform=axQ.transAxes, fontsize=13.5, fontweight="bold",
         va="top", ha="right", color=INK)
axQ.text(0.5, 1.16, "Quantitative phenotypes", transform=axQ.transAxes, fontsize=8.8,
         fontweight="bold", va="top", ha="center", color=INK)
axQ.text(0.5, 1.08, "1,251 pLoF carriers · filled = Bonferroni-significant",
         transform=axQ.transAxes, fontsize=6.6, va="top", ha="center", color=SUBTLE)

# ------------------------------------------------ (a-right) binary forest
binm = ["STATIN_USE", "CAD", "MI", "MACE", "REVASC", "PAD", "STROKE", "HF", "AF", "CV_DEATH"]
for i, p in enumerate(binm):
    r = ph[p]; y = len(binm) - 1 - i
    e, lo, hi = ci(r); sig = float(r["P_bonf"]) < 0.05
    col = GREY if p == "STATIN_USE" else BURG
    axB.plot([lo, hi], [y, y], color=col, lw=1.9, solid_capstyle="round", zorder=2)
    axB.scatter(e, y, s=42, facecolor=col if sig else "white", edgecolor=col, lw=1.4, zorder=3)
axB.axvline(1, color=INK, ls=(0, (3, 3)), lw=0.9)
axB.set_xscale("log"); axB.set_xticks([0.5, 1, 2]); axB.set_xticklabels(["0.5", "1", "2"])
axB.minorticks_off()
axB.set_yticks(range(len(binm)))
axB.set_yticklabels([NICE[k] for k in binm[::-1]], fontsize=7.4)
axB.set_xlabel("odds ratio (95% CI)", fontsize=7.8); despine(axB); axB.tick_params(labelsize=7)
axB.text(0.5, 1.16, "Binary outcomes", transform=axB.transAxes, fontsize=8.8,
         fontweight="bold", va="top", ha="center", color=INK)
axB.text(0.5, 1.08, "no adverse outcome reaches significance (wide CIs)",
         transform=axB.transAxes, fontsize=6.6, va="top", ha="center", color=SUBTLE)

# ------------------------------------------------ (b) severity series
cc = list(csv.DictReader(open(os.path.join(SUP, "pde3b_consequence_class.csv"))))
labs = {"pLoF": "pLoF", "missense": "missense", "synonymous": "synonymous\n(neg. control)"}
colc = {"pLoF": TEAL, "missense": CB["sky"], "synonymous": NS}
for i, r in enumerate(cc):
    y = len(cc) - 1 - i; cl = r["class"]; b, lo, hi = float(r["beta"]), float(r["lo"]), float(r["hi"])
    axS.plot([lo, hi], [y, y], color=colc[cl], lw=3.0, solid_capstyle="round", zorder=2)
    axS.scatter(b, y, s=70, color=colc[cl], edgecolor="black", lw=0.5, zorder=3)
    axS.text(0.055, y, f"n={int(r['n_carriers']):,}", ha="left", va="center", fontsize=6.0,
             color=SUBTLE)
axS.axvline(0, color=INK, ls=(0, (3, 3)), lw=0.9)
axS.set_xlim(-0.47, 0.16)
axS.set_yticks(range(len(cc))); axS.set_yticklabels([labs[r["class"]] for r in cc[::-1]], fontsize=7.2)
axS.set_xlabel("triglyceride effect (SD, 95% CI)", fontsize=7.6); despine(axS); axS.tick_params(labelsize=7)
axS.text(-0.02, 1.21, "b", transform=axS.transAxes, fontsize=13.5, fontweight="bold",
         va="top", ha="right", color=INK)
axS.text(0.5, 1.20, "Effect scales with variant severity", transform=axS.transAxes,
         fontsize=8.8, fontweight="bold", va="top", ha="center", color=INK)
axS.text(0.5, 1.09, "carrier sets analysed separately; synonymous = negative control",
         transform=axS.transAxes, fontsize=6.5, va="top", ha="center", color=SUBTLE)

# ------------------------------------------------ (c) CAD association designs
designs = [("pLoF-burden proxy", 0.861, 0.713, 1.038, "OR", "o", "burden"),
           ("broad-mask proxy", 0.844, 0.704, 1.012, "OR", "o", "burden"),
           ("carrier incident CAD (Cox)", 1.028, 0.838, 1.261, "HR", "s", "cox")]
axD.axvspan(0.75, 1.25, color="#EEF1F5", zorder=0)
for i, (nm, e, lo, hi, lab, mk, grp) in enumerate(designs):
    y = len(designs) - 1 - i
    col = CB["purple"] if grp == "burden" else BURG
    axD.plot([lo, hi], [y, y], color=col, lw=2.4, solid_capstyle="round", zorder=2)
    axD.scatter(e, y, s=62, marker=mk, color=col, edgecolor="black", lw=0.5, zorder=3)
    axD.text(e, y + 0.24, f"{lab} {e:.2f} ({lo:.2f}-{hi:.2f})", va="bottom", ha="center",
             fontsize=6.2, color=INK)
axD.axvline(1, color=INK, ls=(0, (3, 3)), lw=0.9)
axD.set_yticks(range(len(designs)))
axD.set_yticklabels([d[0] for d in designs[::-1]], fontsize=6.9)
axD.set_xlim(0.60, 1.55); axD.set_ylim(-0.75, 2.75)
axD.set_xlabel("effect on CAD (95% CI)", fontsize=7.6)
despine(axD); axD.tick_params(labelsize=7)
axD.text(0.625, 2.62, "92 events / 1,103 carriers · shaded = effect-ratio band 0.75–1.25 (OR/HR)",
         fontsize=5.6, color=SUBTLE, ha="left", va="top")
axD.text(-0.02, 1.21, "c", transform=axD.transAxes, fontsize=13.5, fontweight="bold",
         va="top", ha="right", color=INK)
axD.text(0.5, 1.20, "CAD protection not established", transform=axD.transAxes,
         fontsize=8.8, fontweight="bold", va="top", ha="center", color=BURG)
axD.text(0.5, 1.09, "burden proxy protective but underpowered; carrier Cox null",
         transform=axD.transAxes, fontsize=6.5, va="top", ha="center", color=SUBTLE)

save(fig, "fig4_pde3b")
