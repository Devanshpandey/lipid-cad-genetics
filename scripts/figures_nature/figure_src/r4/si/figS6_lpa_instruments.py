#!/usr/bin/env python3
"""Supplementary Figure (figS6_lpa_instruments): Lp(a) -> CAD under three instrument designs; odds ratio per 1-SD
higher Lp(a) with 95% CI.

Redraw of figure_src/make_figS6_lpa_instruments_nature.py in the r4 theme; same three rows, values and notes.
Rows: Table S3 (primary multivariable model, CAD), Table S7 (FinnGen), Table S23 (CARDIoGRAMplusC4D, IVW), all read
from source/supplementary. No number is typed.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")


def one(f, **kw):
    hit = [r for r in T.rows(os.path.join(SUP, f)) if all(r[k] == v for k, v in kw.items())]
    assert len(hit) == 1, (f, kw, len(hit))
    return hit[0]


def designs():
    s3 = one("TableS3_mvmr.csv", outcome="CAD", exposure="LPA")
    s7 = one("TableS7_finngen_mvmr.csv", exposure="LPA")
    s23 = one("TableS23_mr_cardiogram.csv", exposure="LPA", method="IVW")
    return [("Internal multivariable MR (UK Biobank; %s instruments)" % s3["n_instruments"],
             *(float(s3[k]) for k in ("mvivw_OR", "mvivw_OR_lo95", "mvivw_OR_hi95")),
             "genome-wide rule (rs10455872 at the locus); conditional F = %.0f" % float(s3["cond_F"])),
            ("Two-sample multivariable MR (FinnGen; %s instruments)" % s7["n_instruments"],
             *(float(s7[k]) for k in ("OR_perSD", "CI_lo95", "CI_hi95")),
             r"$\it{cis}$-$\it{LPA}$ locus exempted from pruning; non-overlapping outcome"),
            ("Two-sample univariable MR (CARDIoGRAMplusC4D; %s $\\it{cis}$ variants)" % s23["n_instruments"],
             *(float(s23[k]) for k in ("OR_perSD", "CI_lo95", "CI_hi95")),
             r"locus-dense $\it{cis}$ instrument; non-overlapping outcome")]


def build():
    D = designs()
    P = T.Page(110.0)
    ax = P.ax(196, 8, 244, 72)
    T.null(ax)
    c = T.TRAIT["Lp(a)"]
    for i, (name, o, lo, hi, note) in enumerate(D):
        T.point(ax, i, o, lo, hi, c, "o", filled=T.excludes(lo, hi))
        ax.text(-0.03, i - 0.13, name, transform=ax.get_yaxis_transform(), ha="right", va="baseline", fontsize=T.FS["note"])
        ax.text(-0.03, i + 0.27, note, transform=ax.get_yaxis_transform(), ha="right", va="baseline", fontsize=T.FS["small"], color=T.MUTED)
        ax.text(1.03, i + 0.08, "OR " + T.ci_text(o, lo, hi), transform=ax.get_yaxis_transform(), ha="left", va="baseline", fontsize=T.FS["note"])
        T.src("S6", "", name, "Lp(a) -> CAD", o, lo, hi)
    ax.set_yticks([])
    ax.set_ylim(len(D) - 0.45, -0.55)
    T.logx(ax, 0.97, 1.33, [1, 1.1, 1.2, 1.3])
    ax.set_xlabel("Odds ratio for CAD per SD higher Lp(a) (95% CI)")
    T.clean(ax)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS6_lpa_instruments"))
    return D


if __name__ == "__main__":
    D = build()
    print("wrote supplementary/figS6_lpa_instruments.pdf", [d[1:4] for d in D])
