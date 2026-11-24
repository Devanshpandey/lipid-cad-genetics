#!/usr/bin/env python3
"""Fig. 4: the causal estimates replicate across cohorts and, for LDL-C and triglycerides, across ancestry groups.

a  joint multivariable MR on CAD in UK Biobank and FinnGen; univariable MR against CARDIoGRAMplusC4D
b  within-ancestry MR against MVP CAD with ancestry-stratified GLGC exposures
c  per-SNP LDL-C instrument effects, UK Biobank against GLGC African and Hispanic/admixed-American

Reads the deposited tables (source/supplementary) and figure_data/ (multiancestry_mr_results.csv, xanc_LDL.csv).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theme as T  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
FD = os.path.join(SRC_DIR, "figure_data")
TRAITS = [("LDL_C", "LDL-C", "LDL-C"), ("TRIGLY", "Triglycerides", "TG"), ("HDL_C", "HDL-C", "HDL-C"), ("LPA", "Lp(a)", "Lp(a)")]


def mvmr():
    s3, s7, s23 = (T.rows(os.path.join(SUP, f)) for f in ("TableS3_mvmr.csv", "TableS7_finngen_mvmr.csv", "TableS23_mr_cardiogram.csv"))
    out = {}
    for c, name, _ in TRAITS:
        u = [r for r in s3 if r["outcome"] == "CAD" and r["exposure"] == c]
        f = [r for r in s7 if r["outcome"] == "FinnGen_CAD" and r["exposure"] == c]
        g = [r for r in s23 if r["exposure"] == c and r["method"] == "IVW"]
        assert len(u) == len(f) == len(g) == 1, name
        out[c] = (tuple(float(u[0][k]) for k in ("mvivw_OR", "mvivw_OR_lo95", "mvivw_OR_hi95")),
                  tuple(float(f[0][k]) for k in ("OR_perSD", "CI_lo95", "CI_hi95")),
                  tuple(float(g[0][k]) for k in ("OR_perSD", "CI_lo95", "CI_hi95")))
    return out


def rowsetup(ax, labels, show):
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels if show else [])
    ax.set_ylim(len(labels) - 0.5, -0.5)
    T.clean(ax)


def draw_a(axm, axu):
    M = mvmr()
    for ax in (axm, axu):
        T.null(ax)
    for i, (c, name, hue) in enumerate(TRAITS):
        (u, f, g) = M[c]
        T.point(axm, i - 0.17, *u, T.TRAIT[hue], "s", filled=T.excludes(u[1], u[2]))
        T.point(axm, i + 0.17, *f, T.TRAIT[hue], "^", filled=T.excludes(f[1], f[2]))
        T.point(axu, i, *g, T.TRAIT[hue], "o", filled=T.excludes(g[1], g[2]))
        T.src("4", "a", name, "UK Biobank MVMR", *u); T.src("4", "a", name, "FinnGen MVMR", *f); T.src("4", "a", name, "CARDIoGRAMplusC4D univariable", *g)
    for ax, lab, show in ((axm, "Multivariable", True), (axu, "Univariable", False)):
        rowsetup(ax, [t[1] for t in TRAITS], show)
        T.logx(ax, 0.63, 1.9, [0.8, 1, 1.4, 1.8])
        ax.set_xlabel("CAD OR per SD")
        ax.set_title(lab, fontsize=T.FS["note"], color=T.MUTED, loc="center", pad=3)


def draw_b(axes):
    MM = T.rows(os.path.join(FD, "multiancestry_mr_results.csv"))
    getv = lambda anc, tr, m: next((r for r in MM if r["ancestry"] == anc and r["trait"] == tr and r["method"] == m), None)
    tr = [("LDL", "LDL-C", "LDL-C"), ("logTG", "Triglycerides", "TG"), ("HDL", "HDL-C", "HDL-C")]
    for k, ((anc, lab), ax) in enumerate(zip((("AFR", "African"), ("HIS", "Hispanic/admixed")), axes)):
        T.null(ax)
        for i, (code, name, hue) in enumerate(tr):
            ivw, wm, ic = getv(anc, code, "Inverse variance weighted"), getv(anc, code, "Weighted median"), getv(anc, code, "MR-Egger intercept")
            for r, dy, mk in ((wm, -0.17, "s"), (ivw, 0.17, "o")):
                o, lo, hi = float(r["OR"]), float(r["CI_lo"]), float(r["CI_hi"])
                T.point(ax, i + dy, o, lo, hi, T.TRAIT[hue], mk, filled=T.excludes(lo, hi))
                T.src("4", "b", "%s, %s" % (lab, name), r["method"], o, lo, hi, n_snp=r["n_snp"])
            if ic and ic["P"] and float(ic["P"]) < 0.05:
                ax.text(float(ivw["CI_hi"]) * 1.03, i + 0.17, "‡", fontsize=T.FS["note"], va="center", ha="left")
            ax.text(0.99, i + 0.42, "$n$ = %s" % ivw["n_snp"], transform=ax.get_yaxis_transform(), fontsize=T.FS["small"], va="center", ha="right", color=T.MUTED)
        rowsetup(ax, [t[1] for t in tr], k == 0)
        T.logx(ax, 0.6, 2.2, [0.7, 1, 1.5, 2])
        ax.set_xlabel("CAD OR per SD")
        ax.set_title(lab, fontsize=T.FS["note"], color=T.MUTED, loc="center", pad=3)


def draw_c(ax):
    X = T.rows(os.path.join(FD, "xanc_LDL.csv"))
    S24 = T.rows(os.path.join(SUP, "TableS24_transancestry_transferability.csv"))
    eur = np.array([float(r["ukb_beta"]) for r in X])
    lim = np.nanmax(np.abs(eur)) * 1.08
    ax.plot([-lim, lim], [-lim, lim], color=T.NULL, lw=T.LW["ref"], zorder=1)
    out = []
    for A, name, mk in (("afr", "African", "s"), ("his", "Hispanic/admixed American", "^")):
        b = np.array([float(r[A + "_beta"]) if r[A + "_beta"] not in ("", "NA") else np.nan for r in X])
        ok = ~np.isnan(b)
        r = np.corrcoef(eur[ok], b[ok])[0, 1]
        want = [x for x in S24 if x["trait"] == "LDL_C" and x["ancestry"] == name]
        assert len(want) == 1 and abs(r - float(want[0]["effect_size_correlation_r"])) < 0.00051, (A, r)
        assert not np.any(np.abs(b[ok]) > lim)
        ax.plot(eur[ok], b[ok], ls="none", marker=mk, ms=2.6, mfc="white", mec=T.TRAIT["LDL-C"], mew=0.5, alpha=0.9, zorder=3)
        out.append((A.upper(), mk, r, int(ok.sum())))
        T.src("4", "c", "LDL-C, " + name, "Pearson r", r, r, r, n_snp=int(ok.sum()))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel("UK Biobank, SD per allele"); ax.set_ylabel("GLGC, SD per allele")
    T.clean(ax, left=True)
    return out


def build():
    P = T.Page(214.0)
    P.panel(2, 10, "a", "Across cohorts")
    P.key(12, 22, [("s", T.INK, "white", "UK Biobank"), ("^", T.INK, "white", "FinnGen"), ("o", T.INK, "white", "CARDIoGRAMplusC4D")], gap=7)
    axm = P.ax(58, 46, 62, 132); axu = P.ax(128, 46, 62, 132); draw_a(axm, axu)
    P.panel(204, 10, "b", "Within ancestry (MVP CAD)")
    P.key(214, 22, [("s", T.INK, "white", "weighted median"), ("o", T.INK, "white", "IVW"), ("s", "white", T.INK, "open: CI includes 1")], gap=7)
    axb = [P.ax(254, 46, 56, 132), P.ax(318, 46, 56, 132)]; draw_b(axb)
    P.panel(390, 10, "c", "Instrument transferability")
    axc = P.ax(416, 32, 88, 88)
    out = draw_c(axc)
    y = 150
    for anc, mk, r, n in out:
        P.key(420, y, [(mk, "white", T.TRAIT["LDL-C"], "%s, $r$ = %.2f" % (anc, r))], gap=4)
        y += 10
    P.save(os.path.join(SRC_DIR, "figures", "fig5_replication"))


if __name__ == "__main__":
    build()
    print("wrote Fig. 4")
