#!/usr/bin/env python3
"""
make_main_figures.py — consolidated multipanel main figures for Circ GPM.

Folds the six original panels into three comprehensive figures:
  Figure 1  Shared genetic architecture and causal inference
            (a) genetic correlation  (b) multivariable MR
            (c) LDL-C vs ApoB        (d) MR-BMA exposure selection
  Figure 2  Variant- and gene-level resolution
            (a) colocalization  (b) fine-mapping PCSK9
            (c) fine-mapping LDLR  (d) rare-variant burden
  Figure 3  Convergence on drug targets
            (a) rare-variant allelic series  (b) discovery-only convergence score

All values are real (identical to the prior per-panel figures).
Outputs: fig1_architecture, fig2_variant_gene, fig3_convergence  ({png,pdf})
"""
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import NullLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from figstyle import CB, INK, GREY, FAINT, SUBTLE, LIPID, save

np.random.seed(42)
DIV = LinearSegmentedColormap.from_list("div", [CB["blue"], "#DCE9F2", "white",
                                                "#F7DEE3", CB["red"]])
SEQ = LinearSegmentedColormap.from_list("bu", ["#f0f4fa", "#3b6fb0", "#0b3d91"])


def head(ax, letter, title, sub=None, lx=-0.11):
    """Panel letter (left) + centered bold title (+ optional centered subtitle).
    Matches the shared centered header style used across all main figures."""
    ax.text(lx, 1.16, letter, transform=ax.transAxes, fontsize=13.5,
            fontweight="bold", va="top", ha="left", color=INK)
    ax.text(0.5, 1.155, title, transform=ax.transAxes, fontsize=9.6,
            fontweight="bold", va="top", ha="center", color=INK)
    if sub:
        ax.text(0.5, 1.075, sub, transform=ax.transAxes, fontsize=7.4,
                va="top", ha="center", color=SUBTLE)


# ------------------------------------------------------------------ panels
def draw_rg(ax, fig):
    lipids = ["LDL-C", "ApoB", "Triglycerides", "Total chol.", "Lp(a)", "ApoA1", "HDL-C"]
    outc = ["CAD", "MI", "REVASC", "MACE", "Stroke", "HF", "CV death", "AF"]
    rg = np.array([
        [ 0.113, 0.193, 0.137, 0.020, 0.079,-0.062, 0.091,-0.091],
        [ 0.077, 0.169, 0.119, 0.014, 0.052,-0.139, 0.060,-0.115],
        [ 0.189, 0.232, 0.197, 0.026, 0.099, 0.057,-0.003,-0.097],
        [-0.042, 0.051,-0.001,-0.004, 0.001,-0.188, 0.012,-0.115],
        [ 0.228,-0.032, 0.318,-0.078,-0.504,-0.031,-0.129, 0.242],
        [-0.141,-0.174,-0.176,-0.035,-0.076,-0.149,-0.048, 0.010],
        [-0.190,-0.230,-0.218,-0.040,-0.099,-0.111,-0.047, 0.007]])
    p = np.array([
        [6.6e-3,4e-4,1.9e-3,.067,.304,.470,.093,.030],
        [.064,2e-3,8.1e-3,.240,.488,.106,.254,7.5e-3],
        [1.5e-4,2.4e-4,6.3e-4,.022,.171,.419,.976,7.7e-3],
        [.528,.438,.983,.784,.990,.109,.839,6.7e-3],
        [.341,.912,.319,.410,.387,.948,.670,.328],
        [7.4e-4,1.5e-3,1.2e-4,5.9e-3,.361,.058,.342,.763],
        [5.6e-7,2.0e-5,2.4e-7,2.4e-3,.231,.153,.378,.835]])
    norm = TwoSlopeNorm(vmin=-0.25, vcenter=0.0, vmax=0.25)
    im = ax.imshow(rg, cmap=DIV, norm=norm, aspect="auto")
    ax.set_xticks(np.arange(-.5, len(outc), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(lipids), 1), minor=True)
    ax.grid(which="minor", color="white", lw=2.0); ax.tick_params(which="minor", length=0)
    ax.set_xticks(range(len(outc))); ax.set_xticklabels(outc, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(lipids))); ax.set_yticklabels(lipids, fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.axhline(4.5, color="white", lw=3.5)
    for i in range(rg.shape[0]):
        for j in range(rg.shape[1]):
            v, pv = rg[i, j], p[i, j]
            ax.text(j, i, f"{v:+.2f}".replace("+", " "), ha="center", va="center",
                    fontsize=6.2, color=INK if abs(v) < 0.28 else "white",
                    fontweight="bold" if pv < 1e-3 else "regular")
            if pv < 0.05:
                ax.add_patch(plt.Circle((j+0.37, i-0.37), 0.05, color=INK, zorder=5))
    div = make_axes_locatable(ax); cax = div.append_axes("right", size="3.5%", pad=0.08)
    cb = fig.colorbar(im, cax=cax); cb.set_label("$r_g$", fontsize=8.5)
    cb.outline.set_linewidth(0.6); cb.ax.tick_params(length=2, labelsize=7)
    head(ax, "a", "Genetic correlation", "LD-score regression · dot = P<0.05 · bold = P<0.001")


def _forest_row(ax, y, OR, lo, hi, color, sig, size=46):
    ax.plot([lo, hi], [y, y], color=color, lw=2.4 if sig else 1.4,
            alpha=1.0 if sig else 0.55, solid_capstyle="round", zorder=2)
    ax.scatter([OR], [y], s=size if sig else 26, color=color if sig else "white",
               edgecolor=color if not sig else INK, lw=1.2 if not sig else 0.6,
               marker="s" if sig else "o", zorder=3)


def draw_mvmr(ax):
    main = {
        "LDL-C":  {"CAD":(1.304,1.211,1.405,2.6e-12),"MI":(1.360,1.212,1.527,1.7e-7),
                   "MACE":(1.231,1.142,1.327,5.9e-8),"STROKE":(1.212,1.080,1.361,1.1e-3),
                   "HF":(1.134,1.027,1.251,1.2e-2)},
        "Trigly": {"CAD":(1.116,1.021,1.218,1.5e-2),"MI":(1.150,1.003,1.319,4.6e-2),
                   "MACE":(1.031,0.943,1.127,0.505),"STROKE":(0.883,0.769,1.014,7.8e-2),
                   "HF":(1.011,0.899,1.137,0.850)},
        "HDL-C":  {"CAD":(0.913,0.847,0.984,1.7e-2),"MI":(0.895,0.796,1.006,6.2e-2),
                   "MACE":(0.951,0.882,1.026,0.195),"STROKE":(0.951,0.845,1.069,0.400),
                   "HF":(0.967,0.875,1.069,0.516)},
        "Lp(a)":  {"CAD":(1.116,0.964,1.291,0.141),"MI":(1.139,0.918,1.412,0.237),
                   "MACE":(1.050,0.901,1.223,0.531),"STROKE":(0.900,0.711,1.140,0.382),
                   "HF":(1.099,0.899,1.344,0.358)}}
    outcomes = ["CAD","MI","MACE","STROKE","HF"]; exps = ["LDL-C","Trigly","HDL-C","Lp(a)"]
    col = {"LDL-C":LIPID["LDL-C"],"Trigly":LIPID["TG"],"HDL-C":LIPID["HDL-C"],"Lp(a)":LIPID["Lp(a)"]}
    gap, band = 0.19, 1.15; y = 0; yt, yl = [], []
    for out in outcomes:
        for ei, ex in enumerate(exps):
            OR, lo, hi, pv = main[ex][out]
            _forest_row(ax, y+ei*gap, OR, lo, hi, col[ex], pv < 0.05)
        yt.append(y+gap*1.5); yl.append(out); y += band
    ax.axvline(1.0, color=INK, ls=(0,(3,3)), lw=1)
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=9.5, fontweight="bold")
    ax.set_xscale("log"); ax.set_xlim(0.72, 1.62)
    ax.set_xticks([0.8,1.0,1.2,1.4,1.6]); ax.set_xticklabels(["0.8","1.0","1.2","1.4","1.6"])
    ax.xaxis.set_minor_locator(NullLocator()); ax.invert_yaxis()
    ax.set_xlabel("odds ratio per SD (95% CI)")
    ax.legend(handles=[Patch(color=col[e], label=e) for e in exps],
              loc="lower right", fontsize=8, title="exposure", title_fontsize=8, handlelength=1.0)
    head(ax, "b", "Multivariable MR", "247 shared instruments · filled = P<0.05")


def draw_ldl_apob(ax):
    apob = {"LDL-C":(1.432,0.958,2.141,0.08),"ApoB":(0.979,0.647,1.481,0.92),
            "Trigly":(1.102,0.996,1.219,0.059)}
    bcol = {"LDL-C":LIPID["LDL-C"],"ApoB":"#9B2D3F","Trigly":LIPID["TG"]}  # ApoB darker red (same axis)
    for i, ex in enumerate(["LDL-C","ApoB","Trigly"]):
        OR, lo, hi, pv = apob[ex]
        _forest_row(ax, i, OR, lo, hi, bcol[ex], pv < 0.05, size=52)
        ax.text(hi+0.05, i, f"{OR:.2f}\nP={pv:.2g}", va="center", fontsize=7.4, color=GREY)
    ax.axvline(1.0, color=INK, ls=(0,(3,3)), lw=1)
    ax.set_yticks(range(3)); ax.set_yticklabels(["LDL-C","ApoB","Trigly"], fontsize=9.5, fontweight="bold")
    ax.set_xscale("log"); ax.set_xlim(0.55, 2.7); ax.set_xticks([0.6,1.0,1.6,2.4])
    ax.set_xticklabels(["0.6","1.0","1.6","2.4"]); ax.xaxis.set_minor_locator(NullLocator())
    ax.invert_yaxis(); ax.set_xlabel("odds ratio per SD, CAD")
    head(ax, "c", "LDL-C vs ApoB", "$r_g\\approx1.07$: collinear, unresolved")


def draw_mrbma(ax):
    """Panel d: MR-BMA marginal inclusion probabilities against CAD.
    LDL-C and ApoB are near-interchangeable substitutes (MIP sums to ~1)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "supplementary", "TableS17_mrbma_mip.csv")
    rows = list(csv.DictReader(open(path)))
    cols = {"MIP_Lpa":"Lp(a)", "MIP_TG":"Triglycerides", "MIP_LDL_C":"LDL-C",
            "MIP_ApoA1":"ApoA1", "MIP_HDL_C":"HDL-C", "MIP_ApoB":"ApoB"}
    dflt = next(r for r in rows if r["prior_sigma"] == "0.50"
                and r["prior_inclusion"] == "0.50")
    data = []
    for c, lab in cols.items():
        vals = [float(r[c]) for r in rows]
        data.append((lab, float(dflt[c]), min(vals), max(vals)))
    data.sort(key=lambda r: r[1])
    axiscol = {"Lp(a)": LIPID["Lp(a)"], "Triglycerides": LIPID["TG"], "LDL-C": LIPID["LDL-C"],
               "ApoA1": LIPID["HDL-C"], "HDL-C": LIPID["HDL-C"], "ApoB": LIPID["LDL-C"]}
    for i, (lab, v, lo, hi) in enumerate(data):
        c = axiscol[lab]
        ax.barh(i, v, color=c, alpha=0.9, height=0.62, zorder=2)
        ax.plot([lo, hi], [i, i], color=INK, lw=1.2, zorder=3)
        ax.text(v + 0.02 if v < 0.9 else v - 0.02, i, f"{v:.2f}", va="center",
                ha="left" if v < 0.9 else "right", fontsize=7.2,
                color=INK if v < 0.9 else "white")
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels([r[0] for r in data], fontsize=9.5, fontweight="bold")
    ax.set_xlim(0, 1.12); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.axvline(0.5, color=GREY, ls=(0, (2, 2)), lw=0.9)
    ax.set_xlabel("marginal inclusion probability, CAD")
    ax.text(0.72, 1.15, "LDL-C + ApoB $\\approx$ 1.0\n(substitutes)",
            fontsize=6.8, color=LIPID["LDL-C"], va="center", ha="center", fontweight="bold")
    head(ax, "d", "MR-BMA exposure selection",
         "TG and Lp(a) robust · LDL-C and ApoB interchangeable")


def draw_coloc(ax):
    loci = ["APOE / APOC1", "SORT1 / CELSR2", "LDLR", "TRIB1", "PCSK9"]
    band = ["19q13","1p13","19p13","8q24","1p32"]
    pph4 = [0.998, 0.968, 0.961, 0.993, 0.865]
    fg   = [1.000, 0.998, 0.957, 0.612, 1.000]
    order = np.argsort(pph4)
    for k, idx in enumerate(order):
        ax.plot([0, max(pph4[idx], fg[idx])], [k, k], color="#D5DCE4", lw=2.4, zorder=1,
                solid_capstyle="round")
        ax.scatter(pph4[idx], k, s=88, color=CB["green"], edgecolor="white", lw=1.2,
                   zorder=3, label="UKB (internal)" if k == 0 else "")
        ax.scatter(fg[idx], k, s=66, marker="D", facecolor="white", edgecolor=CB["blue"],
                   lw=1.8, zorder=4, label="FinnGen CAD (external)" if k == 0 else "")
        ax.text(1.05, k, f"{fg[idx]:.2f}", va="center", ha="left", fontsize=7,
                color=CB["blue"], fontweight="bold")
    ax.axvline(0.8, color=CB["red"], ls=(0,(3,3)), lw=1.2)
    labs = [f"{loci[i]}\n$_{{\\mathrm{{{band[i]}}}}}$" for i in order]
    ax.set_yticks(range(len(loci))); ax.set_yticklabels(labs, fontsize=8.2)
    ax.set_xlim(0, 1.3); ax.set_ylim(-0.9, len(loci)-0.4); ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_xlabel("PP.H4 (shared causal variant)")
    ax.legend(loc="lower left", fontsize=6.8)
    head(ax, "a", "Colocalization", "internal + external; TRIB1 fails externally")


def draw_finemap(ax, snps, pip, gene, band_, letter, ylab=False, legend=False):
    x = np.arange(len(snps))
    ax.bar(x, pip, width=0.66, color=[CB["red"] if v>=0.95 else CB["sky"] for v in pip],
           edgecolor="white", lw=0.8, zorder=2)
    ax.axhline(0.95, color=GREY, ls=(0,(2,2)), lw=1)
    ax.set_xticks(x); ax.set_xticklabels(snps, rotation=42, ha="right", fontsize=6.4)
    ax.set_ylim(0, 1.12)
    if ylab: ax.set_ylabel("posterior inclusion prob. (PIP)")
    for xi, v in zip(x, pip):
        ax.text(xi, v+0.02, f"{v:.2f}", ha="center", fontsize=6.0, color=INK)
    if legend:
        ax.legend(handles=[Patch(color=CB["red"], label=r"PIP $\geq$ 0.95"),
                           Patch(color=CB["sky"], label="PIP < 0.95")],
                  loc="upper right", fontsize=6.8)
    head(ax, letter, f"Fine-mapping: {gene}", f"{gene} ({band_})")


_R2BINS = [(0.8, "#E4211C"), (0.6, "#F58221"), (0.4, "#4DAF4A"),
           (0.2, "#377EB8"), (0.0, "#2C3E76")]


def _r2color(r):
    for thr, c in _R2BINS:
        if r >= thr:
            return c
    return "#2C3E76"


def draw_region(axL, axC, csv_name, gene, band, chrom, letter, lead_label,
                ylab=False, legend=False):
    """Stacked regional colocalization (mini-LocusZoom): LDL-C over CAD."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), csv_name)
    P, L, C, R = [], [], [], []
    li = None
    for row in csv.DictReader(open(path)):
        P.append(float(row["pos"]) / 1e6); L.append(float(row["ldl_log10P"]))
        C.append(float(row["cad_log10P"])); R.append(float(row["r2"]))
        if row["is_lead"] == "1":
            li = len(P) - 1
    order = sorted(range(len(P)), key=lambda k: R[k])   # low-LD first (background)
    for ax, Y, ttl, tcol in [(axL, L, "LDL-C", CB["red"]), (axC, C, "CAD", INK)]:
        ax.axvline(P[li], color=GREY, ls=(0, (2, 2)), lw=0.7, zorder=1)
        ax.scatter([P[k] for k in order], [Y[k] for k in order],
                   c=[_r2color(R[k]) for k in order], s=9, edgecolor="none", zorder=2)
        ax.scatter([P[li]], [Y[li]], marker="D", s=48, facecolor="#7A3FA0",
                   edgecolor="black", lw=0.6, zorder=5)
        ax.set_ylim(0, max(Y) * 1.18); ax.margins(x=0.02)
        ax.text(0.03, 0.86, ttl, transform=ax.transAxes, fontsize=8, fontweight="bold", color=tcol)
        ax.tick_params(labelsize=7)
        if ylab:
            ax.set_ylabel(r"$-\log_{10}P$", fontsize=8)
    axL.set_xticklabels([])
    axL.text(P[li], max(L)*1.02, lead_label, ha="center", va="bottom", fontsize=6.4, color="#7A3FA0")
    axC.set_xlabel(f"chr{chrom} position (Mb)", fontsize=8)
    if legend:
        from matplotlib.lines import Line2D
        h = [Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=5,
                    label=lab) for (thr, c), lab in
             zip(_R2BINS, ["0.8–1.0", "0.6–0.8", "0.4–0.6", "0.2–0.4", "<0.2"])]
        axC.legend(handles=h, loc="upper right", fontsize=5.6, title="$r^2$",
                   title_fontsize=6, handletextpad=0.2, borderpad=0.3, labelspacing=0.2)
    head(axL, letter, f"Regional colocalization: {gene} ({band})",
         lx=-0.2 if ylab else -0.1)


def draw_burden(ax):
    entries = [
        ("APOB","LDL-C",293,-1),("PCSK9","LDL-C",267,-1),("LDLR","LDL-C",42,+1),
        ("PCSK9","ApoB",221,-1),("APOB","ApoB",146,-1),
        ("APOC3","Trigly",587,-1),("ANGPTL3","Trigly",129,-1),("APOA5","Trigly",67,+1),
        ("APOC3","HDL-C",324,+1),("ABCA1","HDL-C",248,-1),("CETP","HDL-C",86,+1),
        ("LCAT","HDL-C",81,-1),("LPL","HDL-C",37,-1),
        ("LPA","Lp(a)",463,-1),("PLG","Lp(a)",65,-1),
        ("LDLR","CAD",13.1,+1),("LDLR","MI",6.4,+1)]
    order = ["LDL-C","ApoB","Trigly","HDL-C","Lp(a)","CAD","MI"]
    up, dn = CB["red"], CB["blue"]
    x = 0; xt, xl = [], []
    for tr in order:
        genes = [e for e in entries if e[1] == tr]; xs = [x+k for k in range(len(genes))]
        for (g,_,lp,sgn), xx in zip(genes, xs):
            ax.bar(xx, lp, width=0.82, color=up if sgn>0 else dn, edgecolor="white", lw=0.8, zorder=2)
            ax.text(xx, lp*1.12, g, ha="center", va="bottom", fontsize=6.4,
                    fontweight="bold" if lp > 100 else "regular", color=INK)
        xt.append(np.mean(xs)); xl.append(tr); x = xs[-1]+1.9
    ax.axhline(5.6, color=CB["red"], ls=(0,(3,3)), lw=1)
    ax.set_yscale("log"); ax.set_ylim(3, 1500); ax.set_xlim(-1, x-1)
    ax.text(0.015, 0.97, "red dashed = exome-wide significant (P < $2.5\\times10^{-6}$); log scale",
            transform=ax.transAxes, fontsize=6.6, color=INK, ha="left", va="top")
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=9.5, fontweight="bold")
    ax.set_ylabel("gene-based burden  $-\\log_{10}P$  (log scale)")
    ax.legend(handles=[Patch(color=up, label="burden raises trait"),
                       Patch(color=dn, label="burden lowers trait")],
              loc="upper right", fontsize=8)
    head(ax, "d", "Rare-variant burden", "N $\\approx$ 342–360k · disease-outcome signals (CAD, MI) now visible", lx=-0.045)


def draw_allelic(ax):
    G = [
     ("LDLR","LDL-C", 1.447,0.1050,"Statins (upregulate LDLR) lower LDL"),
     ("CETP","HDL-C", 0.914,0.0462,"CETP inhibitors raise HDL"),
     ("LPL","TG",     0.363,0.0316,"LPL pathway (APOC3/ANGPTL3i) lowers TG"),
     ("ABCA1","HDL-C",-0.343,0.0102,"Tangier gene (no approved drug)"),
     ("LPA","Lp(a)", -0.356,0.0077,"Pelacarsen (LPA ASO) lowers Lp(a)"),
     ("ANGPTL3","TG",-0.545,0.0225,"Evinacumab (ANGPTL3 mAb) lowers TG/LDL"),
     ("PCSK9","LDL-C",-0.608,0.0174,"Evolocumab (PCSK9i) lowers LDL"),
     ("LCAT","HDL-C",-1.028,0.0535,"LCAT deficiency (no approved drug)"),
     ("APOB","LDL-C",-1.118,0.0305,"Mipomersen (APOB ASO) lowers LDL"),
     ("APOC3","TG",  -1.178,0.0227,"Volanesorsen (APOC3 ASO) lowers TG")]
    dm = {"LDLR":1,"CETP":1,"LPL":0,"ABCA1":0,"LPA":1,"ANGPTL3":1,"PCSK9":1,"LCAT":0,"APOB":1,"APOC3":1}
    G = sorted(G, key=lambda r: r[2])
    for i,(g,lip,b,se,drug) in enumerate(G):
        lo,hi = b-1.96*se, b+1.96*se; c = CB["blue"] if b<0 else CB["red"]
        ax.plot([lo,hi],[i,i],color=c,lw=2.4,solid_capstyle="round",zorder=2)
        ax.scatter([b],[i],s=52,color=c,edgecolor="black",lw=0.5,zorder=3)
        ax.text(-2.55, i, g, ha="left", va="center", fontsize=8.6, fontweight="bold", fontstyle="italic")
        ax.text(-2.02, i, f"({lip})", ha="left", va="center", fontsize=6.6, color=GREY)
        tick = r"  $\checkmark$" if dm[g] else ""
        ax.text(1.62, i, drug+tick, ha="left", va="center", fontsize=7.2,
                color=CB["green"] if dm[g] else GREY)
    ax.axvline(0, color=INK, ls=(0,(4,3)), lw=1)
    ax.set_yticks([]); ax.set_ylim(-0.7, len(G)-0.3); ax.set_xlim(-2.6, 1.55)
    ax.set_xlabel("rare loss-of-function burden effect on target lipid (SD, 95% CI)")
    ax.text(0.72, len(G)-0.4, r"raises lipid $\rightarrow$", fontsize=7.2, color=CB["red"], ha="center")
    ax.text(-0.92, len(G)-0.4, r"$\leftarrow$ lowers lipid", fontsize=7.2, color=CB["blue"], ha="center")
    ax.legend(handles=[Patch(color=CB["blue"],label="burden lowers lipid"),
                       Patch(color=CB["red"],label="burden raises lipid"),
                       Patch(color=CB["green"],label=r"$\checkmark$ matches drug direction")],
              loc="lower right", fontsize=7.2, bbox_to_anchor=(1.0,-0.02))
    for s in ("left",): ax.spines[s].set_visible(False)
    head(ax, "a", "Rare-variant allelic series recapitulates lipid pharmacology",
         "loss-of-function direction matches the approved/trial therapeutic", lx=-0.02)


def _load_roc():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "roc_curve_points.csv")
    d = {"integrated": ([], []), "burden": ([], [])}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d[r["curve"]][0].append(float(r["fpr"])); d[r["curve"]][1].append(float(r["tpr"]))
    return d


def draw_roc(ax):
    d = _load_roc()
    ax.plot([0, 1], [0, 1], ls=(0, (3, 3)), color=GREY, lw=1, zorder=1)
    ax.plot(d["integrated"][0], d["integrated"][1], color=CB["red"], lw=2.6, zorder=3,
            label="integrated score  (AUROC 0.97, 95% CI 0.93–1.00)")
    ax.plot(d["burden"][0], d["burden"][1], color=CB["blue"], lw=2.1, zorder=2,
            label="rare-variant burden only  (AUROC 0.90, 0.80–0.99)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xlabel("false-positive rate"); ax.set_ylabel("true-positive rate")
    ax.legend(loc="lower right", fontsize=7.0, handlelength=1.4)
    ax.text(0.30, 0.49, r"$\Delta$AUROC $=+0.07$" "\n" r"(bootstrap $P=0.04$)" "\n"
            r"nested LRT $P=2{\times}10^{-12}$", transform=ax.transAxes,
            fontsize=6.8, color=INK, va="center", linespacing=1.4)
    ax.text(0.30, 0.30, "27 positive genes vs\n17,975 background;\n3,000-sample bootstrap",
            transform=ax.transAxes, fontsize=6.2, color=GREY, va="center", linespacing=1.35)
    head(ax, "b", "Convergence score beats a burden-only baseline",
         "recovery of established lipid targets (genome-wide ROC)", lx=-0.16)


def _load_recovery():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "recovery_curve.csv")
    d = {"integrated": ([], []), "burden": ([], [])}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d[r["curve"]][0].append(int(r["k"])); d[r["curve"]][1].append(float(r["recall"]))
    return d


def draw_recovery(ax):
    d = _load_recovery()
    ax.plot(d["integrated"][0], d["integrated"][1], color=CB["red"], lw=2.6, zorder=3,
            label="integrated score")
    ax.plot(d["burden"][0], d["burden"][1], color=CB["blue"], lw=2.1, zorder=2,
            label="rare-variant burden only")
    for k, rb, ri in [(50, 0.63, 0.74), (100, 0.78, 0.89)]:
        ax.plot([k, k], [rb, ri], color=GREY, lw=0.9, ls=(0, (2, 2)), zorder=1)
        ax.scatter([k, k], [rb, ri], s=18, color=[CB["blue"], CB["red"]], zorder=4)
    ax.set_xscale("log"); ax.set_xlim(1, 18002); ax.set_ylim(0, 1.02)
    ax.set_xlabel("top-ranked genes screened"); ax.set_ylabel("fraction of known targets recovered")
    ax.legend(loc="lower right", fontsize=7.4)
    ax.text(0.03, 0.90, "recall@50:  0.74 vs 0.63\nrecall@100:  0.89 vs 0.78",
            transform=ax.transAxes, fontsize=6.8, color=INK, va="top", linespacing=1.5)
    head(ax, "c", "Faster recovery of known targets",
         "integrated score reaches higher recall at each screening depth", lx=-0.16)


# ------------------------------------------------------------------ composites
def figure1():
    fig = plt.figure(figsize=(13.2, 10.4))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34,
                          height_ratios=[1.0, 0.82], width_ratios=[1.05, 1.0])
    draw_rg(fig.add_subplot(gs[0,0]), fig)
    draw_mvmr(fig.add_subplot(gs[0,1]))
    draw_ldl_apob(fig.add_subplot(gs[1,0]))
    draw_mrbma(fig.add_subplot(gs[1,1]))
    save(fig, "fig1_architecture"); plt.close(fig)


def figure2():
    fig = plt.figure(figsize=(13.2, 9.2))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42,
                          height_ratios=[1.0, 0.92], width_ratios=[1.25, 1.0, 1.0])
    draw_coloc(fig.add_subplot(gs[0,0]))
    gb = gs[0,1].subgridspec(2, 1, hspace=0.14)
    draw_region(fig.add_subplot(gb[0]), fig.add_subplot(gb[1]), "pcsk9_region.csv",
                "PCSK9", "1p32", 1, "b", "rs11591147 (R46L)", ylab=True)
    gc = gs[0,2].subgridspec(2, 1, hspace=0.14)
    draw_region(fig.add_subplot(gc[0]), fig.add_subplot(gc[1]), "ldlr_region.csv",
                "LDLR", "19p13", 19, "c", "rs73015024", legend=True)
    draw_burden(fig.add_subplot(gs[1,:]))
    save(fig, "fig2_variant_gene"); plt.close(fig)


def figure3():
    fig = plt.figure(figsize=(11.6, 10.2))
    gs = fig.add_gridspec(2, 1, hspace=0.34, height_ratios=[0.50, 0.50])
    draw_allelic(fig.add_subplot(gs[0]))
    gsb = gs[1].subgridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.32)
    draw_roc(fig.add_subplot(gsb[0]))
    draw_recovery(fig.add_subplot(gsb[1]))
    save(fig, "fig3_convergence"); plt.close(fig)


if __name__ == "__main__":
    # figure2()/figure3() are superseded by make_fig3_variant_gene.py and
    # make_fig4_convergence.py; this script now owns only Figure 1.
    figure1()
    print("Figure 1 written: fig1_architecture.")
