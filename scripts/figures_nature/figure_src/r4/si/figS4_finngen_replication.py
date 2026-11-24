#!/usr/bin/env python3
"""Supplementary Figure (figS4_finngen_replication): two-sample multivariable MR of the four lipids against FinnGen R11
coronary heart disease; odds ratio per 1-SD higher lipid with 95% CI; open marker = CI includes 1.

Redraw of figS4() in results/final/manuscript/supplementary/make_supp_figures.py in the r4 theme; same rows in the same
order (by odds ratio), same values. Reads source/supplementary/TableS7_finngen_mvmr.csv (md5-identical to the file the
old function read); single column (88 mm).
"""
import os
import sys

import matplotlib.mathtext as mt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

mt.SHRINK_FACTOR, mt.GROW_FACTOR = 0.84, 1 / 0.84    # keeps P-value exponents at or above 5 pt
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
LAB = {"LDL_C": ("LDL-C", "LDL-C"), "HDL_C": ("HDL-C", "HDL-C"), "TRIGLY": ("Triglycerides", "TG"), "LPA": ("Lp(a)", "Lp(a)")}


def build():
    R = T.rows(os.path.join(SUP, "TableS7_finngen_mvmr.csv"))
    n = {r["n_instruments"] for r in R}
    assert len(R) == 4 and len(n) == 1
    rows = sorted(((LAB[r["exposure"]], float(r["OR_perSD"]), float(r["CI_lo95"]), float(r["CI_hi95"]), float(r["pval"])) for r in R),
                  key=lambda z: -z[1])

    P = T.Page(138.0, w=88 * T.MM)
    P.key(48, 10, [("o", T.INK, "white", "95% CI excludes 1"), ("o", "white", T.INK, "includes 1")], gap=9)
    P.text(48, 20, "Two-sample MR, FinnGen R11; %s instruments" % n.pop(), color=T.MUTED)
    ax = P.ax(48, 34, 104, 78)
    T.null(ax)
    for i, ((name, hue), o, lo, hi, p) in enumerate(rows):
        T.point(ax, i, o, lo, hi, T.TRAIT[hue], "o", filled=T.excludes(lo, hi))
        ax.text(1.05, i + 0.11, T.ci_text(o, lo, hi), transform=ax.get_yaxis_transform(), ha="left", va="baseline", fontsize=T.FS["note"])
        ax.text(1.64, i + 0.11, T.p_text(p), transform=ax.get_yaxis_transform(), ha="left", va="baseline", fontsize=T.FS["note"], color=T.MUTED)
        T.src("S4", "", name, "FinnGen MVMR", o, lo, hi, p="%.3g" % p)
    ax.text(1.05, -0.75, "OR (95% CI)", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.text(1.64, -0.75, "P", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0][0] for r in rows])
    ax.set_ylim(len(rows) - 0.5, -0.5)
    T.logx(ax, 0.75, 1.6, [0.8, 1, 1.2, 1.4])
    ax.set_xlabel("Odds ratio for CAD per 1-SD higher lipid")
    T.clean(ax)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS4_finngen_replication"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS4_finngen_replication.pdf")
