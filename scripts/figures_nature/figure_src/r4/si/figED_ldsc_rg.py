#!/usr/bin/env python3
"""Supplementary Figure (figED_ldsc_rg): genome-wide genetic correlation (LD-score regression) between seven lipid
traits and eight cardiovascular endpoints.

Redraw of figure_ed_rg()/draw_rg() in figure_src/make_main_figures.py in the r4 theme. The old function typed the
7 x 8 matrix; this one reads it from source/supplementary/TableS2_ldsc_rg.csv (the same values, unrounded).
Cell shade and printed value: r_g (RdBu_r, -1 to 1); dot, P < 0.05; bold value, P < 0.001 (as before).
"""
import os
import sys

import matplotlib.mathtext as mt
import numpy as np
from matplotlib import cm
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

mt.SHRINK_FACTOR, mt.GROW_FACTOR = 0.84, 1 / 0.84    # keeps the subscript of r_g at or above 5 pt
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
DOT = 2.4                                        # significance dot, pt

LIPIDS = [("LDL_C", "LDL-C"), ("APOB", "ApoB"), ("TRIGLY", "Triglycerides"), ("TOT_CHOL", "Total cholesterol"),
          ("LPA", "Lp(a)"), ("APOA1", "ApoA1"), ("HDL_C", "HDL-C")]
OUTCOMES = [("CAD", "CAD"), ("MI", "MI"), ("REVASC", "Revascularization"), ("MACE", "MACE"),
            ("STROKE", "Ischaemic stroke"), ("HF", "Heart failure"), ("CV_DEATH", "CV death"), ("AF", "AF")]


def matrix():
    t = {(r["exposure"].replace("lipids_", ""), r["outcome"].replace("outcomes_", "")): r for r in T.rows(os.path.join(SUP, "TableS2_ldsc_rg.csv"))}
    rg = np.array([[float(t[(e, o)]["rg"]) for o, _ in OUTCOMES] for e, _ in LIPIDS])
    p = np.array([[float(t[(e, o)]["p"]) for o, _ in OUTCOMES] for e, _ in LIPIDS])
    return rg, p


def build():
    rg, p = matrix()
    cmap, norm = cm.get_cmap("RdBu_r"), Normalize(vmin=-1.0, vmax=1.0)
    nr, nc = rg.shape
    gap = 0.3                                    # the old white rule between Lp(a) and ApoA1
    ys = [i + (gap if i >= 5 else 0.0) for i in range(nr)]

    P = T.Page(222.0)
    P.panel(2, 10, "a", "Genetic correlation")
    P.fig.add_artist(Line2D([81 / P.w], [1 - 21 / P.h], transform=P.fig.transFigure, ls="none", marker="o", ms=DOT, mfc=T.INK, mec="none"))
    t = P.text(86, 23, "P < 0.05", fontsize=T.FS["key"])
    P.text(86 + t.get_window_extent(P.fig.canvas.get_renderer()).width * 72 / P.fig.dpi + 12, 23, "bold value: P < 0.001", fontsize=T.FS["key"])
    ax = P.ax(78, 32, 376, 7 * 21 + 21 * gap)
    for i, y in enumerate(ys):
        for j in range(nc):
            v, pv = rg[i, j], p[i, j]
            ax.add_patch(Rectangle((j - 0.5 + 0.02, y - 0.5 + 0.04), 0.96, 0.92, facecolor=cmap(norm(v)), edgecolor="none", zorder=1))
            ax.text(j, y, ("%+.2f" % v).replace("+", ""), ha="center", va="center", fontsize=T.FS["small"], color=T.INK,
                    fontweight="bold" if pv < 1e-3 else "normal", zorder=3)
            if pv < 0.05:
                ax.plot([j + 0.38], [y - 0.28], ls="none", marker="o", ms=DOT, mfc=T.INK, mec="none", zorder=4)
            T.src("ED", "a", LIPIDS[i][1], OUTCOMES[j][1], v, v, v, p="%.3g" % pv)
    ax.set_xlim(-0.5, nc - 0.5); ax.set_ylim(ys[-1] + 0.5, -0.5)
    ax.set_yticks(ys); ax.set_yticklabels([l for _, l in LIPIDS])
    ax.set_xticks(range(nc)); ax.set_xticklabels([l for _, l in OUTCOMES], rotation=30, ha="right", rotation_mode="anchor")
    ax.tick_params(length=0, pad=3)
    for s in ax.spines.values():
        s.set_visible(False)
    cax = P.ax(463, 32, 2.5 * T.MM, 7 * 21 + 21 * gap)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm, orientation="vertical", ticks=[-1, -0.5, 0, 0.5, 1])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, width=T.LW["axis"], labelsize=T.FS["tick"])
    cb.set_label("$r_g$", fontsize=T.FS["label"], labelpad=4, rotation=0, va="center")
    P.save(os.path.join(SRC_DIR, "figures", "figED_ldsc_rg"))
    return rg, p


if __name__ == "__main__":
    build()
    print("wrote figures/figED_ldsc_rg.pdf")
