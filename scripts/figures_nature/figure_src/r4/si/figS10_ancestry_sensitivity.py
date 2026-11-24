#!/usr/bin/env python3
"""Supplementary Figure: ancestry-specific MR sensitivity (r4 theme redraw of make_figS10_ancestry_sensitivity.py).

For each ancestry (African; Hispanic/admixed-American) and exposure (LDL-C, triglycerides, HDL-C): the IVW,
weighted-median, weighted-mode, MR-Egger and MR-PRESSO outlier-corrected estimates on Million Veteran Program CAD
(OR per SD, 95% CI), with the MR-PRESSO global pleiotropy P where MR-PRESSO reported one (double dagger when < 0.05).

Reads figure_data/si/multiancestry_mr_sensitivity.csv (see figure_data/si/PROVENANCE_B.md); fits nothing.
"""
import os
import sys

from matplotlib.transforms import blended_transform_factory

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DATA = os.path.join(SRC_DIR, "figure_data", "si", "multiancestry_mr_sensitivity.csv")

TRAITS = [("LDL", "LDL-C", "LDL-C"), ("logTG", "Triglycerides", "TG"), ("HDL", "HDL-C", "HDL-C")]
# (method prefix in the table, row label, marker): IVW circle and weighted median square as in Fig. 4b
METH = [("Inverse variance weighted", "IVW", "o"), ("Weighted median", "Weighted median", "s"),
        ("Weighted mode", "Weighted mode", "v"), ("MR Egger", "MR-Egger", "^"),
        ("MR-PRESSO corrected", "MR-PRESSO corrected", "p")]
ANC = [("AFR", "African"), ("HIS", "Hispanic/admixed-American")]
XLO, XHI = 0.4, 3.0
STEP, GAP = 1.0, 0.7          # row pitch and extra gap before each trait header, in data units


def get(rows, anc, tr, prefix):
    hit = [r for r in rows if r["ancestry"] == anc and r["trait"] == tr and r["method"].startswith(prefix)]
    assert len(hit) == 1, (anc, tr, prefix, len(hit))
    return hit[0]


def layout():
    """y of each trait header and each method row (top to bottom, increasing y)."""
    y, heads, slots = 0.0, [], []
    for k, _ in enumerate(TRAITS):
        if k:
            y += GAP
        heads.append(y)
        y += STEP
        slots.append([y + i * STEP for i in range(len(METH))])
        y += len(METH) * STEP
    return heads, slots, y


def draw(ax, rows, anc, heads, slots, ybot, P, show_labels):
    T.null(ax)
    T.band(ax, heads[1] - 0.5, slots[1][-1] + 0.5)
    tr_fig = blended_transform_factory(P.fig.transFigure, ax.transData)
    for (code, name, hue), yh, ys in zip(TRAITS, heads, slots):
        c = T.TRAIT[hue]
        if show_labels:
            ax.text(4 / P.w, yh, name, transform=tr_fig, ha="left", va="center", fontsize=T.FS["label"], fontweight="bold")
        g = get(rows, anc, code, "MR-PRESSO global")
        if g["P"] not in ("", "NA"):
            p = float(g["P"])
            ax.text(0.995, yh, "MR-PRESSO global $P$ = %.3g%s" % (p, "  ‡" if p < 0.05 else ""), transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=T.FS["note"], color=T.INK if p < 0.05 else T.MUTED)
        for (prefix, lab, mk), y in zip(METH, ys):
            r = get(rows, anc, code, prefix)
            e, lo, hi = float(r["OR"]), float(r["CI_lo"]), float(r["CI_hi"])
            T.point(ax, y, e, max(lo, XLO), min(hi, XHI), c, mk, filled=T.excludes(lo, hi),
                    ms=T.MS + 1.0 if mk in ("p", "^", "v") else T.MS)
            for edge, clipped, m in ((XLO, lo < XLO, "<"), (XHI, hi > XHI, ">")):
                if clipped:      # the interval runs past the axis: an arrowhead at the edge
                    ax.plot([edge], [y], marker=m, ms=2.6, mfc=c, mec=c, mew=0.3, ls="none", zorder=3, clip_on=False)
            T.src("S10", anc, "%s, %s" % (name, lab), "OR per SD", e, lo, hi, p=r["P"])
    ys_all = [y for ys in slots for y in ys]
    ax.set_yticks(ys_all)
    ax.set_yticklabels([m[1] for _ in TRAITS for m in METH] if show_labels else [])
    ax.set_ylim(ybot - 0.4, -0.6)
    ax.tick_params(axis="y", pad=4)
    T.logx(ax, XLO, XHI, [0.5, 0.7, 1, 1.5, 2, 3])
    ax.set_xlabel("CAD odds ratio per SD, Million Veteran Program")
    T.clean(ax)


def build():
    rows = T.rows(DATA)
    heads, slots, ybot = layout()
    P = T.Page(292.0)
    P.key(4, 12, [("o", T.INK, "white", "IVW"), ("s", T.INK, "white", "weighted median"), ("v", T.INK, "white", "weighted mode"),
                  ("^", T.INK, "white", "MR-Egger"), ("p", T.INK, "white", "MR-PRESSO corrected")], gap=9)
    x = P.key(4, 23, [("s", "white", T.INK, "open: CI includes 1")], gap=9)
    P.text(x, 23, "‡ MR-PRESSO global $P$ < 0.05; GLGC ancestry-stratified lipid instruments", fontsize=T.FS["key"])
    top, h, w = 52, 206, 190
    for k, ((anc, lab), x0) in enumerate(zip(ANC, (106, 316))):
        P.text(x0, top - 10, lab, fontsize=T.FS["title"], fontweight="bold")
        ax = P.ax(x0, top, w, h)
        draw(ax, rows, anc, heads, slots, ybot, P, k == 0)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS10_ancestry_sensitivity"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS10_ancestry_sensitivity.pdf")
