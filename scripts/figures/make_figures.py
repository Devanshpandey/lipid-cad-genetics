#!/usr/bin/env python3
"""Main-text figures 1-4 (striking house style; real values embedded)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Patch, FancyBboxPatch
from matplotlib.ticker import NullLocator
from figstyle import CB, INK, GREY, FAINT, panel, title_block, save

np.random.seed(42)
DIV = LinearSegmentedColormap.from_list("div", [CB["blue"], "#DCE9F2", "white",
                                                "#F7DEE3", CB["red"]])

# ============================================================
# FIGURE 1 — LDSC genetic-correlation heatmap
# ============================================================
def fig1():
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
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    norm = TwoSlopeNorm(vmin=-0.25, vcenter=0.0, vmax=0.25)
    im = ax.imshow(rg, cmap=DIV, norm=norm, aspect="auto")
    # white gridlines between cells
    ax.set_xticks(np.arange(-.5, len(outc), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(lipids), 1), minor=True)
    ax.grid(which="minor", color="white", lw=2.4); ax.tick_params(which="minor", length=0)
    ax.set_xticks(range(len(outc))); ax.set_xticklabels(outc, rotation=35, ha="right", fontsize=9.5)
    ax.set_yticks(range(len(lipids))); ax.set_yticklabels(lipids, fontsize=9.5)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    # thick group divider
    ax.axhline(4.5, color="white", lw=4)
    ax.text(-1.75, 2.0, "ATHEROGENIC", rotation=90, va="center", ha="center",
            fontsize=8, fontweight="bold", color=CB["red"], alpha=0.9)
    ax.text(-1.75, 5.5, "PROTECTIVE", rotation=90, va="center", ha="center",
            fontsize=8, fontweight="bold", color=CB["blue"], alpha=0.9)
    for i in range(rg.shape[0]):
        for j in range(rg.shape[1]):
            v, pv = rg[i, j], p[i, j]
            strong = pv < 1e-3
            ax.text(j, i, f"{v:+.2f}".replace("+", " "),
                    ha="center", va="center", fontsize=7.2,
                    color=INK if abs(v) < 0.28 else "white",
                    fontweight="bold" if strong else "regular")
            if pv < 0.05:  # significance ring
                ax.add_patch(plt.Circle((j+0.37, i-0.37), 0.055, color=INK, zorder=5))
    cbar = fig.colorbar(im, ax=ax, shrink=0.68, pad=0.02, aspect=22)
    cbar.set_label("genetic correlation  $r_g$", fontsize=9)
    cbar.outline.set_linewidth(0.6); cbar.ax.tick_params(length=2)
    title_block(ax, "Genetic correlation of lipid traits with coronary outcomes",
                "LD-score regression · dot = P<0.05 · bold = P<0.001", y=1.11, sy=1.055)
    save(fig, "fig1_ldsc_rg_heatmap"); plt.close(fig)


# ============================================================
# FIGURE 2 — MVMR forest (3 panels)
# ============================================================
def _forest_row(ax, y, OR, lo, hi, color, sig, size=46):
    ax.plot([lo, hi], [y, y], color=color, lw=2.4 if sig else 1.4,
            alpha=1.0 if sig else 0.55, solid_capstyle="round", zorder=2)
    ax.scatter([OR], [y], s=size if sig else 26, color=color if sig else "white",
               edgecolor=color if not sig else INK, lw=1.2 if not sig else 0.6,
               marker="s" if sig else "o", zorder=3)

def fig2():
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
    col = {"LDL-C":CB["red"],"Trigly":CB["amber"],"HDL-C":CB["blue"],"Lp(a)":CB["purple"]}
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.6, 6.0),
                                        gridspec_kw={"width_ratios":[1.7,0.85,1.05]})
    gap, band = 0.19, 1.15; y = 0; yt, yl = [], []
    for oi, out in enumerate(outcomes):
        for ei, ex in enumerate(exps):
            OR, lo, hi, pv = main[ex][out]
            _forest_row(axA, y+ei*gap, OR, lo, hi, col[ex], pv < 0.05)
        yt.append(y+gap*1.5); yl.append(out); y += band
    axA.axvline(1.0, color=INK, ls=(0,(3,3)), lw=1)
    axA.set_yticks(yt); axA.set_yticklabels(yl, fontsize=10.5, fontweight="bold")
    axA.set_xscale("log"); axA.set_xlim(0.72, 1.62)
    axA.set_xticks([0.8,1.0,1.2,1.4,1.6]); axA.set_xticklabels(["0.8","1.0","1.2","1.4","1.6"])
    axA.xaxis.set_minor_locator(NullLocator()); axA.invert_yaxis()
    axA.set_xlabel("odds ratio per SD (95% CI)")
    axA.legend(handles=[Patch(color=col[e], label=e) for e in exps],
               loc="lower right", fontsize=8.5, title="exposure", title_fontsize=8.5,
               handlelength=1.0)
    axA.text(0.0, 1.10, "Four-lipid multivariable MR", transform=axA.transAxes,
             fontsize=10.5, fontweight="bold")
    axA.text(0.0, 1.045, "247 shared instruments · filled = P<0.05", transform=axA.transAxes,
             fontsize=8.2, color=GREY)
    panel(axA, "a", dx=-0.14)
    # Panel B
    apob = {"LDL-C":(1.432,0.958,2.141,0.08),"ApoB":(0.979,0.647,1.481,0.92),
            "Trigly":(1.102,0.996,1.219,0.059)}
    bcol = {"LDL-C":CB["red"],"ApoB":CB["green"],"Trigly":CB["amber"]}
    for i, ex in enumerate(["LDL-C","ApoB","Trigly"]):
        OR, lo, hi, pv = apob[ex]
        _forest_row(axB, i, OR, lo, hi, bcol[ex], pv < 0.05, size=52)
        axB.text(hi+0.05, i, f"{OR:.2f}\nP={pv:.2g}", va="center", fontsize=7.4, color=GREY)
    axB.axvline(1.0, color=INK, ls=(0,(3,3)), lw=1)
    axB.set_yticks(range(3)); axB.set_yticklabels(["LDL-C","ApoB","Trigly"], fontsize=10.5, fontweight="bold")
    axB.set_xscale("log"); axB.set_xlim(0.55, 2.7); axB.set_xticks([0.6,1.0,1.6,2.4])
    axB.set_xticklabels(["0.6","1.0","1.6","2.4"]); axB.xaxis.set_minor_locator(NullLocator())
    axB.invert_yaxis(); axB.set_xlabel("odds ratio per SD, CAD")
    axB.text(0.0, 1.10, "LDL-C vs ApoB", transform=axB.transAxes, fontsize=10.5, fontweight="bold")
    axB.text(0.0, 1.045, "$r_g\\approx1.07$: collinear, unresolved", transform=axB.transAxes,
             fontsize=8.2, color=GREY)
    panel(axB, "b", dx=-0.16)
    # Panel C — replication
    rep = {"LDL-C":((1.304,1.211,1.405),(1.425,1.356,1.498)),
           "Trigly":((1.116,1.021,1.218),(1.149,1.081,1.221)),
           "HDL-C":((0.913,0.847,0.984),(0.870,0.828,0.914)),
           "Lp(a)":((1.116,0.964,1.291),(1.190,1.143,1.239))}
    for i, ex in enumerate(["LDL-C","Trigly","HDL-C","Lp(a)"]):
        (u,ul,uh),(f,fl,fh) = rep[ex]
        axC.plot([ul,uh],[i+0.15,i+0.15], color=GREY, lw=1.8, solid_capstyle="round")
        axC.scatter([u],[i+0.15], s=30, facecolor="white", edgecolor=GREY, lw=1.4, zorder=3,
                    label="UKB in-sample" if i==0 else "")
        axC.plot([fl,fh],[i-0.15,i-0.15], color=col[ex], lw=2.6, solid_capstyle="round")
        axC.scatter([f],[i-0.15], s=50, color=col[ex], edgecolor=INK, lw=0.6, marker="s", zorder=3,
                    label="FinnGen two-sample" if i==0 else "")
    axC.axvline(1.0, color=INK, ls=(0,(3,3)), lw=1)
    axC.set_yticks(range(4)); axC.set_yticklabels(["LDL-C","Trigly","HDL-C","Lp(a)"], fontsize=10.5, fontweight="bold")
    axC.set_xscale("log"); axC.set_xlim(0.8, 1.6); axC.set_xticks([0.8,1.0,1.2,1.4])
    axC.set_xticklabels(["0.8","1.0","1.2","1.4"]); axC.xaxis.set_minor_locator(NullLocator())
    axC.invert_yaxis(); axC.set_xlabel("odds ratio per SD, CAD")
    axC.legend(loc="upper left", fontsize=7.6, bbox_to_anchor=(0.01, 0.99),
               borderaxespad=0.2, handletextpad=0.5)
    axC.text(0.0, 1.10, "External replication", transform=axC.transAxes, fontsize=10.5, fontweight="bold")
    axC.text(0.0, 1.045, "UKB vs FinnGen CAD · Lp(a) emerges", transform=axC.transAxes,
             fontsize=8.2, color=GREY)
    panel(axC, "c", dx=-0.14)
    save(fig, "fig2_mvmr_forest"); plt.close(fig)


# ============================================================
# FIGURE 3 — colocalization (internal+external) + fine-mapping
# ============================================================
def fig3():
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.4, 4.5),
                                        gridspec_kw={"width_ratios":[1.5,1,1]})
    loci = ["APOE / APOC1", "SORT1 / CELSR2", "LDLR", "TRIB1", "PCSK9"]
    band = ["19q13","1p13","19p13","8q24","1p32"]
    pph4 = [0.998, 0.968, 0.961, 0.993, 0.865]
    fg   = [1.000, 0.998, 0.957, 0.612, 1.000]
    order = np.argsort(pph4); yy = np.arange(len(loci))
    for k, idx in enumerate(order):
        axA.plot([0, max(pph4[idx], fg[idx])], [k, k], color="#D5DCE4", lw=2.4, zorder=1,
                 solid_capstyle="round")
        axA.scatter(pph4[idx], k, s=92, color=CB["green"], edgecolor="white", lw=1.2,
                    zorder=3, label="UKB (internal)" if k == 0 else "")
        axA.scatter(fg[idx], k, s=70, marker="D", facecolor="white",
                    edgecolor=CB["blue"], lw=1.8, zorder=4,
                    label="FinnGen CAD (external)" if k == 0 else "")
        axA.text(1.05, k, f"{fg[idx]:.2f}", va="center", ha="left", fontsize=7.2,
                 color=CB["blue"], fontweight="bold")
    axA.axvline(0.8, color=CB["red"], ls=(0,(3,3)), lw=1.2)
    labs = [f"{loci[i]}\n$_{{\\mathrm{{{band[i]}}}}}$" for i in order]
    axA.set_yticks(yy); axA.set_yticklabels(labs, fontsize=8.6)
    axA.set_xlim(0, 1.3); axA.set_ylim(-0.9, len(loci)-0.4); axA.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
    axA.set_xlabel("PP.H4 (shared causal variant)")
    axA.legend(loc="lower left", fontsize=7)
    axA.text(0.0, 1.10, "Colocalization: internal + external", transform=axA.transAxes,
             fontsize=10.2, fontweight="bold")
    axA.text(0.0, 1.045, "TRIB1 does not replicate externally", transform=axA.transAxes,
             fontsize=8.0, color=GREY)
    panel(axA, "a", dx=-0.13)

    def finemap(ax, snps, pip, gene, band_, letter, ylab=False, legend=False):
        x = np.arange(len(snps))
        ax.bar(x, pip, width=0.66, color=[CB["red"] if v>=0.95 else CB["sky"] for v in pip],
               edgecolor="white", lw=0.8, zorder=2)
        ax.axhline(0.95, color=GREY, ls=(0,(2,2)), lw=1)
        ax.set_xticks(x); ax.set_xticklabels(snps, rotation=42, ha="right", fontsize=6.6)
        ax.set_ylim(0, 1.12)
        if ylab: ax.set_ylabel("posterior inclusion prob. (PIP)")
        for xi, v in zip(x, pip):
            ax.text(xi, v+0.02, f"{v:.2f}", ha="center", fontsize=6.2, color=INK)
        if legend:
            ax.legend(handles=[Patch(color=CB["red"], label=r"PIP $\geq$ 0.95"),
                               Patch(color=CB["sky"], label="PIP < 0.95")],
                      loc="upper right", fontsize=7)
        ax.text(0.0, 1.10, f"SuSiE fine-mapping", transform=ax.transAxes,
                fontsize=10.2, fontweight="bold")
        ax.text(0.0, 1.045, f"{gene} ({band_})", transform=ax.transAxes, fontsize=8.0, color=GREY)
        panel(ax, letter, dx=-0.12)
    finemap(axB, ["rs11591147\n(R46L)","rs2479415","rs559531539","rs746223925","rs55862049","rs472495"],
            [1.00,1.00,1.00,0.985,0.725,0.668], "PCSK9", "1p32", "b", ylab=True, legend=True)
    finemap(axC, ["rs72658867","rs145960625","rs118068660","rs5930","rs2738464","rs112107114"],
            [1.00,1.00,1.00,0.965,0.713,0.538], "LDLR", "19p13", "c")
    save(fig, "fig3_coloc_finemap"); plt.close(fig)


# ============================================================
# FIGURE 4 — rare-variant burden (direction-coded)
# ============================================================
def fig4():
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
    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    x = 0; xt, xl = [], []
    for gi, tr in enumerate(order):
        genes = [e for e in entries if e[1] == tr]; xs = [x+k for k in range(len(genes))]
        for (g,_,lp,sgn), xx in zip(genes, xs):
            ax.bar(xx, lp, width=0.82, color=up if sgn>0 else dn, edgecolor="white", lw=0.8, zorder=2)
            ax.text(xx, lp+10, g, ha="center", va="bottom", fontsize=6.8, rotation=0,
                    fontweight="bold" if lp > 100 else "regular", color=INK)
        xt.append(np.mean(xs)); xl.append(tr); x = xs[-1]+1.9
    ax.axhline(5.6, color=INK, ls=(0,(3,3)), lw=1)
    ax.text(-0.9, 632, "–– exome-wide significant (P < $2.5\\times10^{-6}$)",
            fontsize=7, color=INK, ha="left", va="top")
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=10.5, fontweight="bold")
    ax.set_ylabel("gene-based burden  $-\\log_{10}P$"); ax.set_ylim(0, 650); ax.set_xlim(-1, x-1)
    ax.legend(handles=[Patch(color=up, label="burden raises trait"),
                       Patch(color=dn, label="burden lowers trait")],
              loc="upper right", fontsize=8.5)
    ax.text(0.0, 1.06, "Exome-wide rare-variant burden recovers canonical lipid genes",
            transform=ax.transAxes, fontsize=10.8, fontweight="bold")
    ax.text(0.0, 1.015, "N $\\approx$ 342–360k · effect direction colour-coded",
            transform=ax.transAxes, fontsize=8.2, color=GREY)
    save(fig, "fig4_rvas_burden"); plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("Figures 1-4 regenerated (house style).")
