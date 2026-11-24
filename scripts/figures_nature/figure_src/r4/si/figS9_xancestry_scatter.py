#!/usr/bin/env python3
"""Supplementary Figure: cross-ancestry instrument transferability (r4 theme redraw of make_figS9.py).

Per-SNP UK Biobank lipid-instrument effects against GLGC African and Hispanic/admixed effects, for LDL-C and
triglycerides (aligned to the UK Biobank effect allele), with Pearson r, sign concordance and the number of SNPs.
Markers follow Fig. 4c: open, LDL-C blue or triglyceride orange; square African, triangle Hispanic/admixed.

Reads figure_data/si/xanc_LDL.csv and xanc_TG.csv (see figure_data/si/PROVENANCE_B.md) and checks r against
Supplementary Table S24 (supplementary/TableS24_transancestry_transferability.csv); fits nothing.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FD = os.path.join(SRC_DIR, "figure_data", "si")
S24 = T.rows(os.path.join(SRC_DIR, "supplementary", "TableS24_transancestry_transferability.csv"))

PANELS = [("xanc_LDL.csv", "afr", "LDL-C", "African", "LDL-C", "LDL_C", "African", "s"),
          ("xanc_LDL.csv", "his", "LDL-C", "Hispanic/admixed", "LDL-C", "LDL_C", "Hispanic/admixed American", "^"),
          ("xanc_TG.csv", "afr", "Triglycerides", "African", "TG", "TRIGLY", "African", "s"),
          ("xanc_TG.csv", "his", "Triglycerides", "Hispanic/admixed", "TG", "TRIGLY", "Hispanic/admixed American", "^")]


def draw(ax, f, anc, tlab, anclab, hue, s24_trait, s24_anc, mk):
    rows = T.rows(os.path.join(FD, f))
    x = np.array([float(r["ukb_beta"]) for r in rows])
    y = np.array([float(r[anc + "_beta"]) if r[anc + "_beta"] not in ("", "NA") else np.nan for r in rows])
    ok = ~np.isnan(y)
    x, y = x[ok], y[ok]
    r = np.corrcoef(x, y)[0, 1]
    sc = np.mean(np.sign(x) == np.sign(y))
    want = [w for w in S24 if w["trait"] == s24_trait and w["ancestry"] == s24_anc]
    assert len(want) == 1 and abs(r - float(want[0]["effect_size_correlation_r"])) < 0.00051, (f, anc, r)
    assert int(want[0]["n_instruments"]) == len(x), (f, anc, len(x))
    lim = max(np.abs(x).max(), np.abs(y).max()) * 1.08
    ax.axhline(0, color=T.BAND, lw=0.8, zorder=0)
    ax.axvline(0, color=T.BAND, lw=0.8, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=T.NULL, lw=T.LW["ref"], zorder=1)
    ax.plot(x, y, ls="none", marker=mk, ms=2.6, mfc="white", mec=T.TRAIT[hue], mew=0.5, alpha=0.9, zorder=3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    step = 0.2 if lim > 0.4 else 0.1
    ticks = [k * step for k in range(-int(lim / step), int(lim / step) + 1)]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    labels = [("%g" % round(t, 2)).replace("-", "\u2212") for t in ticks]
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.set_xlabel("UK Biobank %s effect (SD per allele)" % tlab)
    ax.set_ylabel("%s effect (SD per allele)" % anclab)
    ax.text(0.04, 0.97, "$r$ = %.2f\nsign concordance %.0f%%\n$n$ = %d" % (r, sc * 100, len(x)), transform=ax.transAxes,
            ha="left", va="top", fontsize=T.FS["note"], linespacing=1.3)
    T.clean(ax, left=True)
    T.src("S9", tlab, anclab, "Pearson r", r, r, r, sign_concordance="%.3f" % sc, n_snp=len(x))


def build():
    P = T.Page(452.0)
    side = 168
    xs, ys = (70, 330), (20, 244)
    for k, panel in enumerate(PANELS):
        x0, y0 = xs[k % 2], ys[k // 2]
        P.text(x0 - 40, y0, "%s, %s" % (panel[2], panel[3]), fontsize=T.FS["title"], fontweight="bold")
        ax = P.ax(x0, y0 + 10, side, side)
        draw(ax, *panel)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS9_xancestry_scatter"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS9_xancestry_scatter.pdf")
