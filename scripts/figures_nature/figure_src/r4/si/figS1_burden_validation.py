#!/usr/bin/env python3
"""Supplementary Figure: exome-wide rare-variant burden recovers trait-specific gene architecture (r4 theme redraw of
make_figS1.py).

One lollipop panel per trait (LDL-C, triglycerides, HDL-C, Lp(a), CAD, MI): the top eight burden genes by -log10P from
Supplementary Table S6, with the exome-wide threshold (P < 2.5e-6). Hue is the trait (CAD and MI non-lipid grey);
the direction of the burden effect is the marker: an upward triangle raises the trait, a downward triangle lowers it.
The axis is capped at 90; a gene beyond the cap is drawn at the cap and labelled with its true -log10P and an asterisk.

Reads supplementary/TableS6_rvas_siggenes.csv; fits nothing.
"""
import os
import sys
from collections import defaultdict

import matplotlib.mathtext as _mt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

# mathtext shrinks subscripts to 0.7, which would set the "10" of log10 below 5 pt; 0.84 x 6 pt = 5.04 pt
_mt.SHRINK_FACTOR = 0.84

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
S6 = os.path.join(SRC_DIR, "supplementary", "TableS6_rvas_siggenes.csv")
PANELS = [("LDL_C", "LDL-C", "LDL-C"), ("TRIGLY", "Triglycerides", "TG"), ("HDL_C", "HDL-C", "HDL-C"),
          ("LPA", "Lp(a)", "Lp(a)"), ("CAD", "CAD", "non-lipid"), ("MI", "MI", "non-lipid")]
THR = 5.6      # -log10(2.5e-6)
CAP = 90.0     # display cap; the true value is labelled when exceeded
NTOP = 8


def draw(ax, genes, colour, outcome):
    T.null(ax, THR)
    for i, (g, lp, b) in enumerate(genes):
        disp = min(lp, CAP)
        ax.plot([0, disp], [i, i], color=colour, lw=0.8, solid_capstyle="butt", zorder=2)
        ax.plot([disp], [i], marker="^" if b > 0 else "v", ms=T.MS + 0.8, mfc=colour, mec="white", mew=0.4, ls="none", zorder=3,
                clip_on=False)
        ax.text(disp + 3.5, i, "%.0f%s" % (lp, "*" if lp > CAP else ""), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
        T.src("S1", outcome, g, "-log10P", lp, lp, lp, beta=b)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([g[0] for g in genes], fontstyle="italic")
    ax.set_ylim(NTOP - 0.4, -0.6)
    ax.set_xlim(0, max(CAP * 1.12, 12))
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel(r"$-$log$_{10}$ $P$")
    T.clean(ax)


def build():
    by = defaultdict(list)
    for r in T.rows(S6):
        by[r["trait"]].append((r["gene"], float(r["log10P"]), float(r["beta"])))
    P = T.Page(274.0)
    P.key(4, 10, [("^", T.INK, "white", "burden raises the trait"), ("v", T.INK, "white", "burden lowers the trait"),
                  ("line", None, T.NULL, "exome-wide threshold, $P$ = 2.5$\\times$10$^{-6}$")], gap=10)
    P.text(330, 10, "* true value beyond the axis cap of 90", fontsize=T.FS["key"], color=T.MUTED)
    xs, ys, w, h = (52, 222, 392), (34, 152), 108, 82
    for k, (code, lab, hue) in enumerate(PANELS):
        x0, y0 = xs[k % 3], ys[k // 3]
        genes = sorted(by.get(code, []), key=lambda g: -g[1])[:NTOP]
        assert genes, code
        P.text(x0 - 46, y0, lab, fontsize=T.FS["title"], fontweight="bold")
        ax = P.ax(x0, y0 + 8, w, h)
        draw(ax, genes, T.TRAIT[hue], lab)
        if code in ("CAD", "MI"):
            ax.text(0.98, 0.5, "outcome burden\n(few genes)", transform=ax.transAxes, ha="right", va="center",
                    fontsize=T.FS["note"], color=T.MUTED, linespacing=1.1)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS1_burden_validation"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS1_burden_validation.pdf")
