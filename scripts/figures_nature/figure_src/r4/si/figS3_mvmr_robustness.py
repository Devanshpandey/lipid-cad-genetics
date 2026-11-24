#!/usr/bin/env python3
"""Supplementary Figure (figS3_mvmr_robustness): multivariable IVW against multivariable MR-Egger estimates
(log OR per SD) for all 32 exposure x outcome pairs, coloured by exposure, with the identity line (grey) and the
least-squares fitted line (black).

Redraw of figure_src/make_figS3.py in the r4 theme; same points, same statistics.
Reads source/supplementary/TableS3_mvmr.csv (mvivw_b, mvegger_b); single column (88 mm).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
EXP = {"LDL_C": ("LDL-C", "LDL-C"), "TRIGLY": ("Triglycerides", "TG"), "HDL_C": ("HDL-C", "HDL-C"), "LPA": ("Lp(a)", "Lp(a)")}


def build():
    R = T.rows(os.path.join(SUP, "TableS3_mvmr.csv"))
    ivw = np.array([float(x["mvivw_b"]) for x in R])
    egg = np.array([float(x["mvegger_b"]) for x in R])
    cor = np.corrcoef(ivw, egg)[0, 1]
    conc = np.mean(np.sign(ivw) == np.sign(egg))
    slope, intercept = np.polyfit(ivw, egg, 1)
    assert len(R) == 32

    P = T.Page(254.0, w=88 * T.MM)
    P.key(40, 10, [("o", T.TRAIT[h], "white", lab) for lab, h in EXP.values()], gap=9)
    P.key(40, 20, [("line", None, T.NULL, "identity"), ("line", None, T.INK, "fitted line")], gap=9)
    ax = P.ax(40, 30, 200, 200)
    lim = max(np.abs(ivw).max(), np.abs(egg).max()) * 1.15
    ax.axhline(0, color=T.TINT, lw=0.8, zorder=0); ax.axvline(0, color=T.TINT, lw=0.8, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=T.NULL, lw=T.LW["ref"], zorder=1)
    xs = np.array([-lim, lim])
    ax.plot(xs, slope * xs + intercept, color=T.INK, lw=0.8, zorder=2)
    for x, a, b in zip(R, ivw, egg):
        lab, h = EXP[x["exposure"]]
        ax.plot([a], [b], ls="none", marker="o", ms=T.MS, mfc=T.TRAIT[h], mec="white", mew=0.5, zorder=3)
        T.src("S3", "", "%s, %s" % (x["outcome"], lab), "IVW vs Egger", a, b, b)
    ax.text(0.03, 0.97, "r = %.3f\nslope = %.2f\nsign concordance = %.0f%%\nn = %d outcome × exposure" % (cor, slope, 100 * conc, len(R)),
            transform=ax.transAxes, ha="left", va="top", fontsize=T.FS["note"], color=T.MUTED, linespacing=1.3)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    t = [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3]
    ax.set_xticks(t); ax.set_yticks(t)
    ax.set_xlabel("MVMR-IVW effect (log OR per SD)"); ax.set_ylabel("MVMR-Egger effect (log OR per SD)")
    T.clean(ax, left=True)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS3_mvmr_robustness"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS3_mvmr_robustness.pdf")
