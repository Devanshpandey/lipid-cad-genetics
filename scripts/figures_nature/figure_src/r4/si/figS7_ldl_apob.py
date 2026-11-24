#!/usr/bin/env python3
"""Supplementary Figure (figS7_ldl_apob): ApoB against LDL-C per-allele effects at the 151 instruments of the
LDL-C/ApoB model; loci with |standardized residual| > 2 are discordant.

Redraw of figure_src/make_figS7_ldl_apob.py in the r4 theme; same points, same statistics.
Reads source/supplementary/TableS17_apob_ldl_residuals.csv (the deposited copy of the old input,
results/final/manuscript/figures/data_discordant_apob_ldl.csv: identical rows, renamed columns); single column (88 mm).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")


def build():
    R = T.rows(os.path.join(SUP, "TableS17_apob_ldl_residuals.csv"))
    x = np.array([float(r["beta_LDL_C"]) for r in R])
    y = np.array([float(r["beta_ApoB"]) for r in R])
    disc = np.abs(np.array([float(r["std_residual"]) for r in R])) > 2
    r = np.corrcoef(x, y)[0, 1]
    assert len(R) == 151 and disc.sum() == 8 and "%.1f" % (100 * r ** 2) == "94.7", (len(R), disc.sum(), r)

    P = T.Page(244.0, w=88 * T.MM)
    P.key(36, 10, [("o", T.NS, "white", "concordant"), ("o", T.INK, "white", "discordant (n = %d)" % disc.sum())], gap=10)
    ax = P.ax(36, 20, 202, 196)
    ax.plot([-0.7, 0.3], [-0.7, 0.3], color=T.NULL, lw=T.LW["ref"], zorder=1)
    ax.plot(x[~disc], y[~disc], ls="none", marker="o", ms=3.0, mfc=T.NS, mec="white", mew=0.4, zorder=2)
    ax.plot(x[disc], y[disc], ls="none", marker="o", ms=3.6, mfc=T.INK, mec="white", mew=0.5, zorder=3)
    ax.text(0.03, 0.97, "r = %.2f, R² = %.2f (%.1f%% shared)\n%d/%d discordant loci" % (r, r ** 2, 100 * r ** 2, disc.sum(), len(R)),
            transform=ax.transAxes, ha="left", va="top", fontsize=T.FS["note"], color=T.MUTED, linespacing=1.3)
    ax.set_xlim(-0.63, 0.21); ax.set_ylim(-0.68, 0.16)
    ax.set_xticks([-0.6, -0.4, -0.2, 0, 0.2]); ax.set_yticks([-0.6, -0.4, -0.2, 0])
    ax.set_xlabel("LDL-C effect per allele (SD)"); ax.set_ylabel("ApoB effect per allele (SD)")
    T.clean(ax, left=True)
    P.save(os.path.join(SRC_DIR, "figures", "figS7_ldl_apob"))
    return x, y, disc


if __name__ == "__main__":
    build()
    print("wrote figures/figS7_ldl_apob.pdf")
