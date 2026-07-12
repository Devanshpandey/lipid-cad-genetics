#!/usr/bin/env python3
"""
Figure 2 — The LDL/ApoB- and triglyceride-centred architecture replicates across
cohorts and ancestry groups.
  (a) causal estimates by estimand: joint MVMR (UKB + FinnGen) vs univariable IVW
      (CARDIoGRAMplusC4D), kept in separate facets so they are not conflated.
  (b) within-ancestry MR (African, Hispanic/admixed): weighted-median preferred
      estimate + IVW sensitivity; African LDL-C flagged for Egger-intercept.
  (c) cross-ancestry instrument transferability: per-SNP EUR vs AFR/HIS LDL-C
      effects if available (xanc job), else effect-size correlation summary.

Outputs fig5_replication.{pdf,png} (the manuscript's Figure 2 include).
Data: TableS3 (UKB MVMR), S7 (FinnGen MVMR), S24 (CARDIoGRAM), multiancestry_mr_results,
S25 (transferability), data/xanc_ldl.csv (when the cross-ancestry job completes).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
TC = {"LDL-C": LIPID["LDL-C"], "Triglycerides": LIPID["TG"], "HDL-C": LIPID["HDL-C"], "Lp(a)": LIPID["Lp(a)"]}


def hd(ax, letter, title, sub=None, lx=-0.02, cx=0.5):
    ax.text(lx, 1.19, letter, transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="top", ha="right", color=INK)
    ax.text(cx, 1.185, title, transform=ax.transAxes, fontsize=9.2, fontweight="bold",
            va="top", ha="center", color=INK)
    if sub:
        ax.text(cx, 1.09, sub, transform=ax.transAxes, fontsize=7.0, va="top",
                ha="center", color=SUBTLE)


def logx(ax, lo, hi, ticks):
    ax.set_xscale("log"); ax.minorticks_off(); ax.set_xlim(lo, hi)
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks])
    ax.axvline(1, color=INK, ls=(0, (4, 3)), lw=0.9)


fig = plt.figure(figsize=(13.4, 5.0))
outer = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.35, 0.95], wspace=0.42)

# ============================================================= (a) estimands
gA = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=0.10, width_ratios=[1, 0.62])
axMV = fig.add_subplot(gA[0]); axUV = fig.add_subplot(gA[1])
traits = ["LDL-C", "Triglycerides", "HDL-C", "Lp(a)"]
# (trait -> (ukb OR,lo,hi), (finngen OR,lo,hi), (cardiogram OR,lo,hi or None))
mvmr = {
 "LDL-C":        ((1.304, 1.211, 1.405), (1.425, 1.356, 1.498), (1.613, 1.542, 1.688)),
 "Triglycerides":((1.116, 1.021, 1.218), (1.149, 1.081, 1.221), (1.307, 1.253, 1.363)),
 "HDL-C":        ((0.913, 0.847, 0.984), (0.870, 0.828, 0.914), (0.766, 0.733, 0.800)),
 "Lp(a)":        ((1.116, 0.964, 1.291), (1.190, 1.143, 1.239), None),
}
ymap = {t: len(traits) - 1 - i for i, t in enumerate(traits)}
for t in traits:
    y = ymap[t]; c = TC[t]; (u, ul, uh), (f, fl, fh), cg = mvmr[t]
    axMV.plot([ul, uh], [y + 0.16, y + 0.16], color=c, lw=1.9, solid_capstyle="round", zorder=2)
    axMV.scatter(u, y + 0.16, s=34, facecolor="white", edgecolor=c, lw=1.4, marker="o", zorder=3)
    axMV.plot([fl, fh], [y - 0.16, y - 0.16], color=c, lw=1.9, solid_capstyle="round", zorder=2)
    axMV.scatter(f, y - 0.16, s=44, color=c, edgecolor="white", lw=0.7, marker="D", zorder=3)
    if cg:
        g, gl, gh = cg
        axUV.plot([gl, gh], [y, y], color=c, lw=1.9, solid_capstyle="round", zorder=2)
        axUV.scatter(g, y, s=40, color=c, edgecolor="white", lw=0.7, marker="s", zorder=3)
    else:
        axUV.text(1.0, y, "n/a", ha="center", va="center", fontsize=6.2, color=GREY)
logx(axMV, 0.63, 1.85, [0.7, 1.0, 1.4]); logx(axUV, 0.63, 1.85, [0.8, 1.2, 1.6])
axMV.set_yticks(list(ymap.values())); axMV.set_yticklabels(traits, fontsize=7.8)
axUV.set_yticks(list(ymap.values())); axUV.set_yticklabels([])
axMV.set_ylim(-0.6, len(traits) - 0.4); axUV.set_ylim(-0.6, len(traits) - 0.4)
axMV.set_xlabel("CAD OR per SD", fontsize=7.4); axUV.set_xlabel("CAD OR per SD", fontsize=7.4)
despine(axMV); despine(axUV); axMV.tick_params(labelsize=6.8); axUV.tick_params(labelsize=6.8)
axMV.set_title("Joint MVMR", fontsize=7.8, fontweight="bold", pad=3)
axUV.set_title("Univariable", fontsize=7.8, fontweight="bold", pad=3)
axMV.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="w",
                            markeredgecolor=INK, markersize=6, label="UK Biobank"),
                     Line2D([0], [0], marker="D", color="w", markerfacecolor=INK, markersize=6,
                            label="FinnGen")],
            loc="lower right", fontsize=6.2, handletextpad=0.2, borderpad=0.2)
hd(axMV, "a", "Causal estimates by estimand",
   "MVMR (UKB, FinnGen) vs univariable IVW (CARDIoGRAMplusC4D)", cx=0.87)

# ============================================================= (b) ancestry MR
gB = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], wspace=0.12)
mm = list(csv.DictReader(open(os.path.join(DATA, "multiancestry_mr_results.csv"))))


def getv(anc, trait, method):
    for r in mm:
        if r["ancestry"] == anc and r["trait"] == trait and r["method"] == method:
            return r
    return None


TR = [("LDL", "LDL-C"), ("logTG", "Triglycerides"), ("HDL", "HDL-C")]
for k, (anc, anclab) in enumerate([("AFR", "African"), ("HIS", "Hispanic/admixed")]):
    ax = fig.add_subplot(gB[k])
    for i, (tcode, tlab) in enumerate(TR):
        y = len(TR) - 1 - i; c = TC[tlab]
        ivw = getv(anc, tcode, "Inverse variance weighted")
        wm = getv(anc, tcode, "Weighted median")
        icept = getv(anc, tcode, "MR-Egger intercept")
        nsnp = ivw["n_snp"]
        # IVW sensitivity (small open) + weighted-median preferred (large filled)
        ax.plot([float(ivw["CI_lo"]), float(ivw["CI_hi"])], [y + 0.13, y + 0.13],
                color=c, lw=1.4, alpha=0.55, solid_capstyle="round", zorder=2)
        ax.scatter(float(ivw["OR"]), y + 0.13, s=26, facecolor="white", edgecolor=c,
                   lw=1.1, marker="o", zorder=3)
        ax.plot([float(wm["CI_lo"]), float(wm["CI_hi"])], [y - 0.13, y - 0.13],
                color=c, lw=2.4, solid_capstyle="round", zorder=2)
        ax.scatter(float(wm["OR"]), y - 0.13, s=52, color=c, edgecolor="white", lw=0.7,
                   marker="s", zorder=3)
        if icept and icept["P"] and float(icept["P"]) < 0.05:
            ax.text(float(ivw["CI_hi"]) + 0.02, y + 0.13, "‡", fontsize=8, color=CB["red"],
                    va="center", ha="left", fontweight="bold")
        ax.text(0.015, y, f"{tlab}", transform=ax.get_yaxis_transform(), fontsize=6.6,
                va="bottom", ha="left", color=INK)
        ax.text(0.015, y - 0.34, f"n={nsnp}", transform=ax.get_yaxis_transform(),
                fontsize=5.6, va="center", ha="left", color=SUBTLE)
    logx(ax, 0.6, 2.15, [0.7, 1.0, 1.5, 2.0])
    ax.set_yticks([]); ax.set_ylim(-0.7, len(TR) - 0.4)
    ax.set_xlabel("CAD OR per SD", fontsize=7.4); despine(ax); ax.tick_params(labelsize=6.6)
    ax.set_title(anclab, fontsize=7.8, fontweight="bold", pad=3)
    if k == 0:
        hd(ax, "b", "Within-ancestry MR (MVP CAD)",
           "weighted median (filled) + IVW (open); ‡ = Egger-intercept P<0.05", cx=1.06)
        ax.legend(handles=[Line2D([0], [0], marker="s", color="w", markerfacecolor=INK,
                                  markersize=6, label="weighted median"),
                           Line2D([0], [0], marker="o", color="w", markerfacecolor="w",
                                  markeredgecolor=INK, markersize=6, label="IVW")],
                  loc="lower right", fontsize=6.0, handletextpad=0.2, borderpad=0.2)

# ============================================================= (c) transferability
axC = fig.add_subplot(outer[2])
xanc = os.path.join(DATA, "xanc_LDL.csv")
if os.path.exists(xanc):
    rows = list(csv.DictReader(open(xanc)))
    eur = np.array([float(r["ukb_beta"]) for r in rows])
    for A, col, mk in (("afr", LIPID["outcome"], "s"), ("his", CB["purple"], "^")):
        b = np.array([float(r[A + "_beta"]) if r[A + "_beta"] not in ("", "NA") else np.nan for r in rows])
        ok = ~np.isnan(b)
        r = np.corrcoef(eur[ok], b[ok])[0, 1]
        axC.scatter(eur[ok], b[ok], s=11, color=col, alpha=0.6, edgecolor="none", marker=mk,
                    label=f"AFR $r$={r:.2f}" if A == "afr" else f"HIS $r$={r:.2f}")
    lim = np.nanmax(np.abs(eur)) * 1.08
    axC.axline((0, 0), slope=1, color=GREY, ls=(0, (4, 3)), lw=0.9)
    axC.set_xlim(-lim, lim); axC.set_ylim(-lim, lim)
    axC.set_xlabel("UK Biobank LDL-C effect", fontsize=7.4)
    axC.set_ylabel("non-European effect", fontsize=7.4)
    axC.legend(loc="upper left", fontsize=6.6)
    hd(axC, "c", "Instrument transferability", "per-SNP LDL-C effects, UKB vs non-European")
else:
    anc = ["EUR", "AFR", "HIS"]; x = np.arange(3); w = 0.36
    ldl = [0.991, 0.933, 0.904]; tg = [0.982, 0.712, 0.786]
    axC.bar(x - w / 2, ldl, w, color=LIPID["LDL-C"], label="LDL-C")
    axC.bar(x + w / 2, tg, w, color=LIPID["TG"], label="Triglycerides")
    for xi, v in zip(x - w / 2, ldl):
        axC.text(xi, v + 0.01, f"{v:.2f}", ha="center", fontsize=6.0, color=INK)
    for xi, v in zip(x + w / 2, tg):
        axC.text(xi, v + 0.01, f"{v:.2f}", ha="center", fontsize=6.0, color=INK)
    axC.set_xticks(x); axC.set_xticklabels(["European", "African", "Hispanic"], fontsize=7.2)
    axC.set_ylim(0, 1.1); axC.set_ylabel("effect-size correlation ($r$)", fontsize=7.4)
    axC.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=2, fontsize=6.6,
               handlelength=1.0, columnspacing=1.0)
    hd(axC, "c", "Effect transferability", "lipid instruments across ancestries (GLGC)")
despine(axC); axC.tick_params(labelsize=6.8)

save(fig, "fig5_replication")
