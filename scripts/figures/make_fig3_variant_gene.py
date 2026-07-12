#!/usr/bin/env python3
"""
Figure 3 — Common- and rare-variant evidence converges on lipid-CAD loci and genes.
  (a) locus convergence matrix: internal / FinnGen(default,stringent) / coloc.susie / burden
  (b) PCSK9 regional colocalization (LDL-C and CAD association tracks)
  (c) LDLR  regional colocalization (LDL-C and CAD association tracks)
  (d) gene-level convergence matrix: burden -log10P x trait/outcome, + coloc & fine-map

Outputs fig2_variant_gene.{pdf,png} (the manuscript's Figure 3 include).
All values real: coloc from TableS4/S9/S16, fine-map from TableS5, burden from TableS6,
regional tracks from strengthen/regional/*.csv.
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RED, BLUE, TEAL, BURG = LIPID["LDL-C"], LIPID["HDL-C"], LIPID["rare"], LIPID["outcome"]
PPCMAP = LinearSegmentedColormap.from_list("pp", ["#EAF0F6", "#9CC3DE", "#2E6E9E", "#123B5E"])


def hd(ax, letter, title, sub=None):
    ax.text(-0.02, 1.16, letter, transform=ax.transAxes, fontsize=13.5,
            fontweight="bold", va="top", ha="right", color=INK)
    ax.text(0.5, 1.155, title, transform=ax.transAxes, fontsize=9.6,
            fontweight="bold", va="top", ha="center", color=INK)
    if sub:
        ax.text(0.5, 1.075, sub, transform=ax.transAxes, fontsize=7.4,
                va="top", ha="center", color=SUBTLE)


fig = plt.figure(figsize=(13.2, 10.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30, height_ratios=[1.0, 1.0],
                      width_ratios=[1.05, 1.0])

# =============================================== (a) locus convergence matrix
axA = fig.add_subplot(gs[0, 0])
loci = ["PCSK9", "SORT1/CELSR2", "APOE/APOC1", "LDLR", "TRIB1"]
cols = ["UKB\ninternal", "FinnGen\ndefault", "FinnGen\nstringent", "coloc\n.susie", "rare\nburden"]
# PP.H4 values (NaN = not applicable); burden col: 1 = exome-wide burden gene at locus
pp = {
    "PCSK9":       [1.00, 1.00, 1.00, 1.00, 1],
    "SORT1/CELSR2":[0.97, 0.998, 0.976, np.nan, 0],
    "APOE/APOC1":  [1.00, 1.00, 0.999, np.nan, 0],
    "LDLR":        [0.967, 0.957, 0.691, 0.936, 1],
    "TRIB1":       [0.993, 0.612, 0.136, np.nan, 0],
}
leadv = {"PCSK9":"rs11591147 (R46L)", "SORT1/CELSR2":"rs12740374",
         "APOE/APOC1":"rs7412 region", "LDLR":"rs72658867 (+rs5930)", "TRIB1":"rs28601761"}
nrow, ncol = len(loci), 5
for i, loc in enumerate(loci):
    y = nrow - 1 - i
    vals = pp[loc]
    for j in range(4):
        v = vals[j]
        x = j
        if np.isnan(v):
            axA.text(x, y, "n/a", ha="center", va="center", fontsize=6.4, color=GREY)
            continue
        passed = v >= 0.8
        s = 120 + 460 * v
        axA.scatter(x, y, s=s, c=[PPCMAP(v)] if passed else ["white"],
                    edgecolor=PPCMAP(v) if passed else GREY, lw=1.4, zorder=3)
        axA.text(x, y, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=6.3,
                 color="white" if (passed and v > 0.6) else INK, zorder=4, fontweight="bold")
    # rare-burden column: triangle if exome-wide burden gene
    if vals[4]:
        axA.scatter(4, y, s=190, marker="^", c=[TEAL], edgecolor="white", lw=0.8, zorder=3)
    else:
        axA.text(4, y, "–", ha="center", va="center", fontsize=8, color=GREY)
axA.axvline(-0.5, color="none")
for j, c in enumerate(cols):
    axA.text(j, nrow - 0.32, c, ha="center", va="bottom", fontsize=6.9, fontweight="bold",
             color=INK, linespacing=0.95)
axA.set_yticks(range(nrow)); axA.set_yticklabels(loci[::-1], fontsize=8.2, fontstyle="italic")
# fine-mapped lead annotations to the right
for i, loc in enumerate(loci):
    axA.text(4.75, nrow - 1 - i, leadv[loc], ha="left", va="center", fontsize=6.0, color=SUBTLE)
axA.set_xlim(-0.6, 7.6); axA.set_ylim(-0.6, nrow - 0.05)
axA.set_xticks([]); axA.tick_params(left=False)
for sp in axA.spines.values():
    sp.set_visible(False)
hd(axA, "a", "Locus convergence matrix",
   "posterior PP.H4 (circle size/shade); triangle marks an exome-wide burden gene")

# =============================================== (b,c) regional colocalization
def regional(subgs, tag, csvname, lead_label, susie_txt):
    inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=subgs, hspace=0.12, height_ratios=[1, 1])
    axL = fig.add_subplot(inner[0]); axC = fig.add_subplot(inner[1])
    rows = list(csv.DictReader(open(os.path.join(DATA, csvname))))
    pos = np.array([int(r["pos"]) for r in rows]) / 1e6
    ldl = np.array([float(r["ldl_log10P"]) for r in rows])
    cad = np.array([float(r["cad_log10P"]) for r in rows])
    r2 = np.array([float(r["r2"]) for r in rows])
    lead_i = int(np.argmax([r["is_lead"] == "1" for r in rows])) if any(r["is_lead"] == "1" for r in rows) else int(np.argmax(ldl))
    order = np.argsort(r2)
    scL = None
    for ax, yv, lab, col in ((axL, ldl, "LDL-C", RED), (axC, cad, "CAD", BURG)):
        sc = ax.scatter(pos[order], yv[order], c=r2[order], cmap="viridis",
                        s=14, edgecolor="none", vmin=0, vmax=1, zorder=2)
        if ax is axL:
            scL = sc
        ax.scatter(pos[lead_i], yv[lead_i], s=70, marker="D", facecolor=col,
                   edgecolor="white", lw=0.8, zorder=4)
        ax.set_ylabel(f"{lab}\n$-\\log_{{10}}P$", fontsize=7.2, linespacing=0.95)
        ax.set_ylim(-0.05 * yv.max(), yv.max() * 1.18)
        despine(ax)
        ax.tick_params(labelsize=6.6)
    axL.set_xticklabels([])
    axC.set_xlabel(f"chr position (Mb)", fontsize=7.2)
    axL.text(0.98, 0.9, f"{lead_label}", transform=axL.transAxes, ha="right", va="top",
             fontsize=6.6, color=INK, fontstyle="italic")
    axC.text(0.98, 0.9, susie_txt, transform=axC.transAxes, ha="right", va="top",
             fontsize=6.4, color=TEAL, fontweight="bold")
    axL.text(-0.02, 1.30, "b" if tag == "PCSK9" else "c", transform=axL.transAxes,
             fontsize=13.5, fontweight="bold", va="top", ha="right", color=INK)
    axL.text(0.5, 1.30, f"{tag} regional colocalization", transform=axL.transAxes,
             fontsize=9.6, fontweight="bold", va="top", ha="center", color=INK)
    axL.text(0.5, 1.13, "LDL-C and CAD peaks coincide (points shaded by LD $r^2$ to lead)",
             transform=axL.transAxes, fontsize=7.0, va="top", ha="center", color=SUBTLE)
    cax = axL.inset_axes([0.04, 0.62, 0.20, 0.05])
    cb = fig.colorbar(scL, cax=cax, orientation="horizontal")
    cb.set_ticks([0, 1]); cb.ax.tick_params(labelsize=5.6, length=2)
    cb.set_label("LD $r^2$", fontsize=6.0, labelpad=1)
    return axL


regional(gs[0, 1], "PCSK9", "pcsk9_region.csv", "rs11591147 (p.R46L), PIP=1.0",
         "coloc.susie PP.H4 = 1.00")
regional(gs[1, 0], "LDLR", "ldlr_region.csv", "rs72658867 lead; rs5930 coding secondary",
         "coloc.susie PP.H4 = 0.94")

# =============================================== (d) gene convergence matrix
axD = fig.add_subplot(gs[1, 1])
genes = ["LDLR", "PCSK9", "APOB", "APOC3", "ANGPTL3", "LPA", "ABCA1", "CETP", "LCAT", "LPL"]
tcols = ["LDL-C", "ApoB", "TG", "HDL-C", "Lp(a)", "CAD"]
# (log10P, sign): sign +1 raises trait, -1 lowers
burd = {
 "LDLR":   {"LDL-C":(42.4,+1), "CAD":(13.1,+1)},
 "PCSK9":  {"LDL-C":(267,-1), "ApoB":(221,-1)},
 "APOB":   {"LDL-C":(293,-1), "ApoB":(145,-1), "TG":(144,-1)},
 "APOC3":  {"TG":(587,-1), "HDL-C":(324,+1), "LDL-C":(16.9,-1), "ApoB":(16.4,-1)},
 "ANGPTL3":{"TG":(129,-1), "LDL-C":(42.3,-1), "HDL-C":(29,-1)},
 "LPA":    {"Lp(a)":(463,-1)},
 "ABCA1":  {"HDL-C":(247,-1), "LDL-C":(15.9,-1)},
 "CETP":   {"HDL-C":(86.5,+1)},
 "LCAT":   {"HDL-C":(81.5,-1)},
 "LPL":    {"TG":(29.8,+1), "HDL-C":(37,-1)},
}
coloc_g = {"PCSK9", "LDLR"}
drug_g = {"LDLR", "PCSK9", "APOB", "APOC3", "ANGPTL3", "LPA", "CETP"}
CAP = 60.0
for i, g in enumerate(genes):
    y = len(genes) - 1 - i
    for j, t in enumerate(tcols):
        if t in burd.get(g, {}):
            lp, sgn = burd[g][t]
            s = 40 + 300 * min(lp, CAP) / CAP
            col = RED if sgn > 0 else BLUE
            axD.scatter(j, y, s=s, c=[col], edgecolor="black", lw=0.5, zorder=3)
    # coloc + fine-map columns
    axD.scatter(len(tcols), y, s=110, marker="D",
                c=[TEAL] if g in coloc_g else ["white"],
                edgecolor=TEAL if g in coloc_g else GREY, lw=1.1, zorder=3)
    axD.scatter(len(tcols) + 1, y, s=110, marker="D",
                c=[TEAL] if g in coloc_g else ["white"],
                edgecolor=TEAL if g in coloc_g else GREY, lw=1.1, zorder=3)
    axD.text(-0.7, y, g, ha="right", va="center", fontsize=7.6,
             fontstyle="italic", color=INK)
    if g in drug_g:
        axD.scatter(-2.75, y, s=52, marker="*", c=[INK], edgecolor="none", zorder=4)
allc = tcols + ["coloc", "fine-map"]
axD.set_xticks(range(len(allc)))
axD.set_xticklabels(allc, fontsize=6.9, rotation=40, ha="right")
axD.set_yticks([]); axD.set_xlim(-3.2, len(allc) - 0.3); axD.set_ylim(-0.6, len(genes) - 0.4)
for sp in axD.spines.values():
    sp.set_visible(False)
axD.tick_params(bottom=False, left=False)
leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE, markersize=8, label="lowers trait"),
       Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markersize=8, label="raises trait"),
       Line2D([0], [0], marker="D", color="w", markerfacecolor=TEAL, markersize=8, label="coloc / fine-map"),
       Line2D([0], [0], marker="*", color="w", markerfacecolor=INK, markersize=11, label="established gene", linestyle="None")]
axD.legend(handles=leg, loc="lower left", bbox_to_anchor=(-0.02, -0.30), ncol=2,
           fontsize=6.6, handletextpad=0.3, columnspacing=1.0, frameon=False)
hd(axD, "d", "Gene-level convergence matrix",
   "burden $-\\log_{10}P$ (dot size, capped 60); direction by colour")

save(fig, "fig2_variant_gene")
