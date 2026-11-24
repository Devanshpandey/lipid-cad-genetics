#!/usr/bin/env python3
"""Supplementary Figure (loci and genes; output figures/fig2_variant_gene): common- and rare-variant convergence.
r4 theme redraw of make_fig4_variant_gene_nature.py; same panels, rows and values.

a  locus convergence matrix: PP.H4 for the UK Biobank LDL-C -> CAD pair of Table S4 at each locus, FinnGen under the
   default and the stringent prior (Table S9) and coloc.susie (LDL-C rows of Table S16); circle area = PP.H4, filled
   at PP.H4 >= 0.8; a triangle marks an exome-wide burden gene at the locus (Table S6); the fine-mapped lead at right.
   Ink and greys only: the FinnGen analysis at TRIB1 is the triglyceride one (Table S9), so no exposure hue is implied.
b  gene-level convergence matrix: every Table S6 burden cell of ten canonical genes, marker area = -log10P (capped 60),
   hue = trait, upward triangle = burden raises the trait, downward = lowers; squares = colocalization and fine-mapping
   support (filled) or its absence (open); star = drug target (Table S14).
c  therapeutic direction: cis-MR OR for CAD per SD higher cognate lipid (Table S32; ANGPTL3 shown for LDL-C), hue = the
   cognate lipid, filled = CI excludes 1; rare-burden -log10P at left (Table S15; Table S6 LDL-C row for LDLR), filled =
   exome-wide; the direction class of Table S32 printed at right.

Every printed value is read from the deposited tables and asserted against the expected print. The SORT1/CELSR2 lead is
rs12740374. Reads source/supplementary and figure_data/si/ed4_gene_support.csv; fits nothing.
"""
import os
import sys

import matplotlib.mathtext as _mt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

# mathtext shrinks subscripts to 0.7, which would set the "10" of log10 below 5 pt; 0.84 x 6 pt = 5.04 pt
_mt.SHRINK_FACTOR = 0.84

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
FD = os.path.join(SRC_DIR, "figure_data", "si")
R = lambda f: T.rows(os.path.join(SUP, f))
LOG10P = r"$-$log$_{10}$ $P$"

# ------------------------------------------------------------------ (a) data, as in the r3 script
LOCI = ["PCSK9", "SORT1/CELSR2", "APOE/APOC1", "LDLR", "TRIB1"]
COLS = ["UKB\nLDL-C–CAD", "FinnGen\ndefault", "FinnGen\nstringent", "SuSiE\ncoloc", "rare\nburden"]
UKB_ROW = {"PCSK9": ("LDL_C", "rs11591147"), "SORT1/CELSR2": ("LDL_C", "rs4970834"),
           "APOE/APOC1": ("LDL_C", "rs531660643"), "LDLR": ("LDL_C", "rs10423733"), "TRIB1": ("LDL_C", "rs28601761")}
FG_KEY = {"PCSK9": "PCSK9", "SORT1/CELSR2": "SORT1_CELSR2", "APOE/APOC1": "APOE_APOC1", "LDLR": "LDLR", "TRIB1": "TRIB1"}
SHIPPED = {"PCSK9": ["1.00", "1.00", "1.00", "1.00"], "SORT1/CELSR2": ["0.92", "1.00", "0.98", "0.86"],
           "APOE/APOC1": ["1.00", "1.00", "1.00", "1.00"], "LDLR": ["0.97", "0.96", "0.69", "0.94"],
           "TRIB1": ["0.99", "0.61", "0.14", "0.99"]}
BURDEN_LOCUS = {"PCSK9": "PCSK9", "LDLR": "LDLR"}
LEADS = {"PCSK9": ("rs11591147", "rs11591147 (R46L)"), "SORT1/CELSR2": ("rs12740374", "rs12740374"),
         "APOE/APOC1": ("rs429358", "rs429358"), "LDLR": ("rs72658867", "rs72658867"), "TRIB1": ("rs28601761", "rs28601761")}


def data_a():
    s4, s5, s6 = R("TableS4_coloc_loci.csv"), R("TableS5_finemapping.csv"), R("TableS6_rvas_siggenes.csv")
    s9 = {r["locus"]: r for r in R("TableS9_finngen_coloc.csv")}
    s16 = {r["locus"]: r for r in R("TableS16_coloc_susie.csv") if r["exposure"] == "LDL_C"}
    pp = {}
    for loc in LOCI:
        exp, snp = UKB_ROW[loc]
        hit = [r for r in s4 if r["exposure"] == exp and r["outcome"] == "CAD" and r["sentinel_snp"] == snp]
        assert len(hit) == 1, (loc, hit)
        fg = s9[FG_KEY[loc]]
        pp[loc] = [float(hit[0]["PP.H4"]), float(fg["PPH4_default"]), float(fg["PPH4_p12_1e6"]), float(s16[loc]["max_PP.H4"])]
        assert not any(np.isnan(v) for v in pp[loc]), loc
        assert ["%.2f" % v for v in pp[loc]] == SHIPPED[loc], (loc, pp[loc])
        assert any(r["top_snp"] == LEADS[loc][0] and float(r["top_pip"]) >= 0.99 for r in s5), loc
        if loc in BURDEN_LOCUS:
            assert any(r["gene"] == BURDEN_LOCUS[loc] for r in s6), loc
    return pp


# ------------------------------------------------------------------ (b) data
GENES = ["LDLR", "PCSK9", "APOB", "APOC3", "ANGPTL3", "LPA", "ABCA1", "CETP", "LCAT", "LPL"]
TCOLS = [("LDL-C", "LDL_C", "LDL-C"), ("ApoB", "APOB", "ApoB"), ("TG", "TRIGLY", "TG"), ("HDL-C", "HDL_C", "HDL-C"),
         ("Lp(a)", "LPA", "Lp(a)"), ("CAD", "CAD", "non-lipid")]
CELLS = {"LDLR": {"LDL-C": +1, "ApoB": +1, "CAD": +1}, "PCSK9": {"LDL-C": -1, "ApoB": -1},
         "APOB": {"LDL-C": -1, "ApoB": -1, "TG": -1, "HDL-C": +1},
         "APOC3": {"TG": -1, "HDL-C": +1, "LDL-C": -1, "ApoB": -1},
         "ANGPTL3": {"TG": -1, "LDL-C": -1, "ApoB": -1, "HDL-C": -1}, "LPA": {"Lp(a)": -1},
         "ABCA1": {"HDL-C": -1, "LDL-C": -1, "TG": -1}, "CETP": {"HDL-C": +1}, "LCAT": {"HDL-C": -1, "TG": +1},
         "LPL": {"TG": +1, "HDL-C": -1}}
CAP = 60.0


def data_b():
    s6 = R("TableS6_rvas_siggenes.csv")
    burd, omitted = {}, []
    for g in GENES:
        burd[g] = {}
        for t, code, _ in TCOLS:
            hit = [r for r in s6 if r["gene"] == g and r["trait"] == code]
            assert len(hit) <= 1, (g, t, hit)
            if not hit:
                assert t not in CELLS[g], (g, t, "drawn but not in Table S6")
                continue
            lp, sgn = float(hit[0]["log10P"]), 1 if float(hit[0]["beta"]) > 0 else -1
            if t in CELLS[g]:
                assert sgn == CELLS[g][t], (g, t, sgn)
                burd[g][t] = (lp, sgn)
            else:
                omitted.append((g, t, lp, sgn))
    assert not omitted, ("Table S6 cells not drawn in panel b", omitted)
    sup = {r["gene"]: r for r in T.rows(os.path.join(FD, "ed4_gene_support.csv"))}
    coloc_g = {g for g, r in sup.items() if r["coloc"] == "1"}
    finemap_g = {g for g, r in sup.items() if r["finemap"] == "1"}
    s14 = {r["gene"]: r["target_category"] for r in R("TableS14_canonical_positive_set.csv")}
    drug_g = {g for g in GENES if "rug target" in s14[g] and "not a direct" not in s14[g]}
    assert drug_g == {"PCSK9", "APOB", "APOC3", "ANGPTL3", "LPA", "CETP"}, drug_g
    return burd, coloc_g, finemap_g, drug_g


# ------------------------------------------------------------------ (c) data
ORDER = ["PCSK9", "APOC3", "LPA", "NPC1L1", "LDLR", "HMGCR", "ANGPTL3", "CETP"]
PRIMARY = {"ANGPTL3": "LDL-C"}
HUE = {"LDL-C": "LDL-C", "Triglycerides": "TG", "Lp(a)": "Lp(a)", "HDL-C": "HDL-C"}


def data_c():
    tgt = {r["gene"]: r for r in R("TableS15_lipid_lowering_targets.csv")}
    cis, s6 = R("TableS32_drug_target_cis_MR.csv"), R("TableS6_rvas_siggenes.csv")
    out = []
    for g in ORDER:
        rows_g = [r for r in cis if r["gene"] == g]
        r = next(x for x in rows_g if x["exposure"] == PRIMARY.get(g, x["exposure"]))
        est, lo, hi = float(r["CAD_OR"]), float(r["CAD_CI_low"]), float(r["CAD_CI_high"])
        conc = r["direction_concordant"]
        cls = "matches the drug" if conc == "yes" else ("matches, not significant" if conc.startswith("yes") else "discordant on common variants")
        b = tgt.get(g)
        if b is not None:
            bp, bfill = float(b["rare_burden_neglog10P"]), b["recovered_exomewide"] == "yes"
        else:       # LDLR is not a Table S15 target: cognate-trait row of Table S6
            hit = [x for x in s6 if x["gene"] == g and x["trait"] == "LDL_C"]
            assert len(hit) == 1 and g == "LDLR", g
            bp, bfill = float(hit[0]["log10P"]), True
        out.append((g, r["exposure"], est, lo, hi, cls, bp, bfill))
    return out


# ------------------------------------------------------------------ drawing
AREA = lambda v: 160 + 300 * v      # circle area in pt^2: 460 at PP.H4 = 1, large enough at 0.14 to hold its value
def draw_a(ov, P, pp):
    xl, xs, xb, xlead = 62, [84, 114, 144, 174], 202, 216
    y0, dy = 64, 27
    for j, c in enumerate(COLS):
        ov.text(xs[j] if j < 4 else xb, 50, c, ha="center", va="bottom", fontsize=T.FS["note"], color=T.MUTED, linespacing=1.05)
    ov.text(xlead, 50, "fine-mapped\nlead variant", ha="left", va="bottom", fontsize=T.FS["note"], color=T.MUTED, linespacing=1.05)
    for i, loc in enumerate(LOCI):
        y = y0 + i * dy
        ov.text(xl, y, loc, ha="right", va="center", fontsize=T.FS["label"], fontstyle="italic")
        for j, v in enumerate(pp[loc]):
            passed = v >= 0.8
            ov.scatter([xs[j]], [y], s=AREA(v), marker="o", zorder=3, clip_on=False,
                       c=[T.tint(T.INK, 0.28) if passed else "white"], edgecolors="none" if passed else T.NULL, linewidths=0.8)
            ov.text(xs[j], y, "%.2f" % v, ha="center", va="center", fontsize=T.FS["small"], fontweight="bold", zorder=4)
            T.src("fig2_variant_gene", "a", loc, COLS[j].replace("\n", " "), v, v, v)
        if loc in BURDEN_LOCUS:
            ov.scatter([xb], [y], s=70, marker="^", c=[T.INK], edgecolors="white", linewidths=0.4, zorder=3, clip_on=False)
        else:
            ov.text(xb, y, "–", ha="center", va="center", fontsize=T.FS["label"], color=T.MUTED)
        ov.text(xlead, y, LEADS[loc][1], ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)


def draw_b(ov, P, burd, coloc_g, finemap_g, drug_g):
    xstar, xl = 278, 326
    xt = [338 + k * 21 for k in range(len(TCOLS))]
    xc, xf = 465, 491
    y0, dy = 64, 13.6
    for (t, _, _), x in zip(TCOLS, xt):
        ov.text(x, 55, t, ha="center", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    for x, t in ((xc, "coloc"), (xf, "fine-map")):
        ov.text(x, 55, t, ha="center", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    for i, g in enumerate(GENES):
        y = y0 + i * dy
        ov.text(xl, y, g, ha="right", va="center", fontsize=T.FS["label"], fontstyle="italic")
        if g in drug_g:
            ov.scatter([xstar], [y], s=22, marker="*", c=[T.INK], edgecolors="none", zorder=4, clip_on=False)
        for (t, _, hue), x in zip(TCOLS, xt):
            if t in burd[g]:
                lp, sgn = burd[g][t]
                ov.scatter([x], [y], s=1.2 * (13 + 97 * min(lp, CAP) / CAP), marker="^" if sgn > 0 else "v",
                           c=[T.TRAIT[hue]], edgecolors="white", linewidths=0.3, zorder=3, clip_on=False)
                T.src("fig2_variant_gene", "b", g, t, lp, lp, lp, sign=sgn)
        for x, member in ((xc, g in coloc_g), (xf, g in finemap_g)):
            ov.scatter([x], [y], s=24, marker="s", c=[T.INK if member else "white"], edgecolors=T.INK, linewidths=0.6,
                       zorder=3, clip_on=False)
    # size key, on the panel's own scale
    yk = 208
    ov.text(280, yk, "marker size, %s:" % LOG10P, ha="left", va="center", fontsize=T.FS["key"])
    for xk, v_, lab_ in ((360, 15, "15"), (386, 30, "30"), (412, CAP, "≥60")):
        ov.scatter([xk], [yk], s=1.2 * (13 + 97 * v_ / CAP), marker="^", c=["white"], edgecolors=T.INK, linewidths=0.5, clip_on=False)
        ov.text(xk + 7, yk, lab_, ha="left", va="center", fontsize=T.FS["key"])


def draw_c(ax, P, rows, xlab, xbur):
    T.null(ax)
    tr = blended_transform_factory(P.fig.transFigure, ax.transData)
    for i, (g, exp, est, lo, hi, cls, bp, bfill) in enumerate(rows):
        c = T.TRAIT[HUE[exp]]
        T.point(ax, i, est, lo, hi, c, "s", filled=T.excludes(lo, hi))
        ax.text(xlab / P.w, i, g, transform=tr, ha="right", va="center", fontsize=T.FS["label"], fontstyle="italic")
        ax.plot([xbur / P.w], [i], transform=tr, marker="o", ms=T.MS - 0.4, ls="none", clip_on=False,
                mfc=T.INK if bfill else "white", mec=T.INK if not bfill else "white", mew=0.8 if not bfill else 0.4)
        ax.text((xbur + 6) / P.w, i, "%.0f" % bp if bp >= 10 else "%.1f" % bp, transform=tr, ha="left", va="center",
                fontsize=T.FS["note"])
        ax.text(1.02, i, cls, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"],
                color=T.MUTED if cls != "matches the drug" else T.INK)
        T.src("fig2_variant_gene", "c", g, exp + " cis-MR", est, lo, hi, burden=bp, class_=cls)
    ax.text(xbur / P.w - 0.004, -0.95, "rare burden\n" + LOG10P, transform=tr, ha="left", va="bottom", fontsize=T.FS["note"],
            color=T.MUTED, linespacing=1.05)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([])
    ax.set_ylim(len(rows) - 0.5, -0.6)
    T.logx(ax, 0.45, 3.0, [0.5, 0.7, 1, 1.5, 2, 3])
    ax.set_xlabel(r"CAD odds ratio per SD higher cognate lipid ($\mathit{cis}$-MR)")
    T.clean(ax)


def build():
    pp = data_a()
    b = data_b()
    c = data_c()
    P = T.Page(392.0)
    ov = P.ax(0, 0, P.w, P.h)
    ov.set_xlim(0, P.w)
    ov.set_ylim(P.h, 0)
    ov.axis("off")
    P.panel(2, 10, "a", "Locus convergence matrix")
    P.key(12, 22, [("o", T.tint(T.INK, 0.28), T.tint(T.INK, 0.28), "PP.H4 ≥ 0.8"), ("o", "white", T.NULL, "below 0.8"),
                   ("^", T.INK, "white", "exome-wide burden gene")], gap=8)
    P.text(12, 32, "circle area scales with PP.H4", fontsize=T.FS["key"], color=T.MUTED)
    draw_a(ov, P, pp)
    P.panel(270, 10, "b", "Gene-level convergence matrix")
    P.key(280, 22, [("^", T.INK, "white", "burden raises trait"), ("v", T.INK, "white", "lowers trait"),
                    ("*", T.INK, T.INK, "drug target")], gap=8)
    P.key(280, 32, [("s", T.INK, "white", "colocalization or fine-map support"), ("s", "white", T.INK, "none")], gap=8)
    draw_b(ov, P, *b)
    P.panel(2, 232, "c", "Therapeutic direction, both allele-frequency classes")
    P.key(12, 244, [("s", T.TRAIT["LDL-C"], "white", "LDL-C"), ("s", T.TRAIT["TG"], "white", "Triglycerides"),
                    ("s", T.TRAIT["Lp(a)"], "white", "Lp(a)"), ("s", T.TRAIT["HDL-C"], "white", "HDL-C"),
                    ("s", "white", T.INK, "open: CI includes 1"), ("o", T.INK, "white", "burden exome-wide"),
                    ("o", "white", T.INK, "not exome-wide")], gap=8)
    axc = P.ax(110, 270, 236, 92)
    draw_c(axc, P, c, xlab=60, xbur=72)
    P.save(os.path.join(SRC_DIR, "figures", "fig2_variant_gene"))


if __name__ == "__main__":
    build()
    print("wrote figures/fig2_variant_gene.pdf")
