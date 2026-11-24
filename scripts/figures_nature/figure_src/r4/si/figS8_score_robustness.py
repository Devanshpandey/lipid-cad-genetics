#!/usr/bin/env python3
"""Supplementary Figure (score robustness), redrawn in the r4 theme.

a  convergence score against Open Targets L2G (Spearman rho)
b  nested layer ablation of the score (AUROC for canonical recovery)
c  PDE3B genome-wide rank across 108 weight/cap specifications
d  recovery (AUROC) by target class

Same content as figure_src/make_figS8_score_robustness_nature.py; reads the deposited tables
source/supplementary/TableS18, S19, S21 and S22 (the same-named files in submission/supplementary_data are a
superseded run). Fits nothing beyond the Spearman correlation the old script printed.
"""
import os
import sys

import numpy as np
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")


def rich(P, x, y, parts, **kw):
    """One line of text in runs, each (string, italic)."""
    rend = P.fig.canvas.get_renderer()
    for s, it in parts:
        t = P.text(x, y, s, fontstyle="italic" if it else "normal", **kw)
        x += t.get_window_extent(rend).width * 72 / P.fig.dpi
    return x


def title(P, x, y, letter, parts, sub=None):
    P.panel(x, y, letter)
    rich(P, x + 9, y, parts, fontsize=T.FS["title"], fontweight="bold")
    if sub:
        P.text(x + 9, y + 9.5, sub, color=T.MUTED)


def draw_a(ax):
    d = T.rows(os.path.join(SUP, "TableS18_l2g_benchmark.csv"))
    l2g = np.array([float(r["l2g"]) for r in d]); integ = np.array([float(r["integrated"]) for r in d])
    can = np.array([r["canonical"] == "1" for r in d])
    assert can.sum() == 5
    rho = spearmanr(integ, l2g).correlation
    ax.plot(l2g[~can], integ[~can], ls="none", marker="o", ms=2.4, mfc=T.NS, mec="none", alpha=0.8, zorder=2)
    ax.plot(l2g[can], integ[can], ls="none", marker="o", ms=T.MS, mfc=T.INK, mec="white", mew=0.5, zorder=3)
    ax.text(0.03, 0.97, r"Spearman $\rho$ = %.2f" % rho, transform=ax.transAxes, ha="left", va="top", fontsize=T.FS["note"])
    ax.set_xlim(-0.03, 0.9); ax.set_ylim(-0.4, 8.8)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8]); ax.set_yticks([0, 2, 4, 6, 8])
    ax.set_xlabel("Open Targets L2G score"); ax.set_ylabel("Convergence (integrated) score")
    T.clean(ax, left=True)
    return rho, len(d)


def dotrows(ax, labels, vals, faces, texts):
    for i, (v, fc, s) in enumerate(zip(vals, faces, texts)):
        ax.plot([v], [i], ls="none", marker="o", ms=T.MS, mfc=fc, mec="white", mew=0.5, zorder=3)
        ax.text(v + 0.004, i, s, ha="left", va="center", fontsize=T.FS["note"])
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.set_xlim(0.85, 1.0); ax.set_xticks([0.85, 0.9, 0.95, 1.0])
    T.clean(ax)


def draw_b(ax):
    ab = T.rows(os.path.join(SUP, "TableS19_layer_ablation.csv"))
    assert len(ab) == 5
    ab = ab[:4]                      # the fifth row (integrated score) is the score itself
    labs = ["Burden only", "+ coloc", "+ outcome burden", "Full"]
    auc = [float(r["AUROC_canonical_recovery"]) for r in ab]
    faces = [T.NS, T.INK, T.NS, T.INK]
    dotrows(ax, labs, auc, faces, ["%.2f" % v for v in auc])
    ax.set_xlabel("AUROC (canonical recovery)")
    return auc


def draw_c(ax):
    ws = T.rows(os.path.join(SUP, "TableS21_score_weight_sensitivity.csv"))
    assert len(ws) == 108
    r = np.array([int(float(x["pde3b_rank"])) for x in ws])
    bins = np.logspace(np.log10(max(r.min(), 1)), np.log10(r.max() * 1.05), 20)
    ax.hist(r, bins=bins, color=T.NS, edgecolor="white", lw=0.3, zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(13, 700); ax.set_ylim(0, 40)
    ax.axvline(180, color=T.INK, ls=(0, (3, 2)), lw=T.LW["ref"], zorder=3)
    ax.text(190, 40, "top 1%\n(rank 180)", ha="left", va="top", fontsize=T.FS["note"], linespacing=1.0)
    ax.xaxis.set_major_locator(FixedLocator([20, 50, 100, 200, 500]))
    ax.xaxis.set_major_formatter(FixedFormatter(["20", "50", "100", "200", "500"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_xlabel("$\\it{PDE3B}$ genome-wide rank"); ax.set_ylabel("Number of weight\nspecifications", linespacing=1.0)
    good = np.array([x["cap"] not in ("", "NA", "inf", "nan") and np.isfinite(float(x["cap"])) and float(x["wB"]) > float(x["wC"]) for x in ws])
    rg, rb = r[good], r[~good]
    ax.text(33, 30, "rank %d–%d in %d/%d\nspecifications (lipid term\ncapped and weighted\nabove coloc)" % (rg.min(), rg.max(), len(rg), len(ws)),
            ha="left", va="top", fontsize=T.FS["note"], linespacing=1.1)
    ax.text(455, 24, "rank %d–%d\nin the other %d" % (rb.min(), rb.max(), len(rb)), ha="right", va="top", fontsize=T.FS["note"], linespacing=1.1)
    T.clean(ax, left=True)
    return rg.min(), rg.max(), len(rg), rb.min(), rb.max(), len(rb)


def draw_d(ax):
    st = T.rows(os.path.join(SUP, "TableS22_stratified_auroc.csv"))
    labs = [s["stratum"].replace("approved/trial drug target", "Approved/trial drug targets")
            .replace("Mendelian dyslipidemia", "Mendelian dyslipidaemia").replace("all canonical", "All canonical") for s in st]
    auc = [float(s["AUROC"]) for s in st]
    dotrows(ax, labs, auc, [T.INK] * len(st), ["%.2f (n = %s)" % (v, s["n"]) for v, s in zip(auc, st)])
    ax.set_xlabel("AUROC")
    return list(zip(labs, auc))


def build():
    P = T.Page(300.0)
    title(P, 2, 10, "a", [("Weakly correlated with common-variant L2G", False)],
          "rare-variant genes fall outside the L2G-scored set")
    P.key(12, 32, [("o", T.INK, "white", "canonical target (n = 5)"), ("o", T.NS, "white", "other genes")], gap=8)
    axa = P.ax(42, 42, 196, 92); rho, n = draw_a(axa)
    title(P, 262, 10, "b", [("Which layer earns the gain", False)], "nested ablation of the score (points, AUROC)")
    P.key(272, 32, [("o", T.INK, "white", "with the coloc layer"), ("o", T.NS, "white", "without it")], gap=8)
    axb = P.ax(340, 44, 140, 84); auc = draw_b(axb)
    title(P, 2, 162, "c", [("PDE3B", True), (" rank depends on weighting", False)], "108 burden:coloc:outcome / cap specifications")
    axc = P.ax(46, 184, 200, 86); ranks = draw_c(axc)
    title(P, 262, 162, "d", [("Recovery holds across target classes", False)], "positive set analysed by evidence type (points, AUROC)")
    axd = P.ax(356, 192, 110, 78); strata = draw_d(axd)
    P.save(os.path.join(SRC_DIR, "figures", "figS8_score_robustness"))
    print("rho %.3f over %d genes; ablation AUROC %s; PDE3B ranks %s; strata %s" % (rho, n, auc, ranks, strata))


if __name__ == "__main__":
    build()
