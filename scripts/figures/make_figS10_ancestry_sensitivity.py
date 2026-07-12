#!/usr/bin/env python3
"""
Supplementary Figure S10 — Ancestry-specific MR sensitivity.
For each ancestry (African, Hispanic/admixed) and exposure (LDL-C, triglycerides,
HDL-C), the IVW / weighted-median / weighted-mode / MR-Egger / MR-PRESSO
outlier-corrected estimates on Million Veteran Program CAD (OR per SD, 95% CI),
with the MR-PRESSO global pleiotropy P. LDL-C and triglycerides remain causal under
every pleiotropy-robust estimator in both ancestries, including after MR-PRESSO
outlier removal at loci flagged for pleiotropy.

Data: strengthen/multiancestry_mr/multiancestry_mr_sensitivity.csv (fetched locally).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
mm = list(csv.DictReader(open(os.path.join(DATA, "multiancestry_mr_sensitivity.csv"))))
TC = {"LDL": LIPID["LDL-C"], "logTG": LIPID["TG"], "HDL": LIPID["HDL-C"]}
NICE = {"LDL": "LDL-C", "logTG": "Triglycerides", "HDL": "HDL-C"}
# (method key match, short label, marker)
METH = [("Inverse variance weighted", "IVW", "o"),
        ("Weighted median", "wtd median", "s"),
        ("Weighted mode", "wtd mode", "v"),
        ("MR Egger", "MR-Egger", "^"),
        ("MR-PRESSO corrected", "MR-PRESSO corr.", "D")]


def get(anc, tr, mprefix):
    for r in mm:
        if r["ancestry"] == anc and r["trait"] == tr and r["method"].startswith(mprefix):
            return r
    return None


fig, axes = plt.subplots(1, 2, figsize=(11.4, 6.0), sharex=True)
TR = ["LDL", "logTG", "HDL"]
for ax, (anc, anclab) in zip(axes, [("AFR", "African"), ("HIS", "Hispanic/admixed-American")]):
    yt, yl = [], []; y = 0
    for tr in TR:
        c = TC[tr]
        ax.text(0.6, y + 1.0, NICE[tr], fontsize=8.2, fontweight="bold", color=INK, ha="left")
        gp = get(anc, tr, "MR-PRESSO global")
        for mname, mlab, mk in METH:
            r = get(anc, tr, mname)
            if not r or r["OR"] in ("", "NA"):
                continue
            e, lo, hi = float(r["OR"]), float(r["CI_lo"]), float(r["CI_hi"])
            filled = mname.startswith("MR-PRESSO")
            ax.plot([lo, hi], [y, y], color=c, lw=1.8, solid_capstyle="round", zorder=2)
            ax.scatter(e, y, s=42 if filled else 34, color=c if filled else "white",
                       edgecolor=c, lw=1.2, marker=mk, zorder=3)
            yt.append(y); yl.append(mlab); y -= 1
        if gp and gp["P"] not in ("", "NA"):
            p = float(gp["P"])
            ax.text(2.18, y + 2.5, f"PRESSO global P={p:.3g}" + ("  ‡" if p < 0.05 else ""),
                    fontsize=6.0, color=CB["red"] if p < 0.05 else SUBTLE, ha="right", va="center")
        y -= 0.9
    ax.axvline(1, color=INK, ls=(0, (4, 3)), lw=0.9)
    ax.set_xscale("log"); ax.minorticks_off(); ax.set_xlim(0.55, 2.3)
    ax.set_xticks([0.7, 1.0, 1.5, 2.0]); ax.set_xticklabels(["0.7", "1.0", "1.5", "2.0"])
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=6.8); ax.set_ylim(y + 0.6, 1.7)
    ax.set_xlabel("CAD odds ratio per SD"); despine(ax)
    ax.set_title(anclab, fontsize=9.4, fontweight="bold", pad=6)
axes[0].text(-0.20, 1.05, "S10", transform=axes[0].transAxes, fontsize=11, fontweight="bold",
             color=INK, ha="left", va="bottom")
fig.suptitle("Ancestry-specific MR is robust to pleiotropy-tolerant estimators and outlier removal",
             fontsize=10.4, fontweight="bold", y=1.01)
fig.text(0.5, 0.965, "GLGC ancestry-stratified lipids to Million Veteran Program CAD; "
         "filled diamond = MR-PRESSO outlier-corrected; ‡ = MR-PRESSO global P<0.05",
         ha="center", fontsize=7.2, color=SUBTLE)
save(fig, "figS10_ancestry_sensitivity")
