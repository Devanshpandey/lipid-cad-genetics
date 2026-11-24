#!/usr/bin/env python3
"""Supplementary Figure: colocalization and fine-mapping robustness (r4 theme redraw of make_figS5_coloc.py).

a  FinnGen external PP.H4 at the five replication loci under the default and a stringent (p12 = 1e-6) prior, paired per
   locus and classified robust / prior-sensitive / fails. Circle = default prior, square = stringent prior;
   filled = PP.H4 > 0.8. Ink only: the FinnGen analysis is LDL-C for four loci and triglycerides for TRIB1 (Table S9).
b  LDLR credible sets recovered under an increasing SuSiE signal cap L (10, 15, 20); the fit did not converge at L = 20,
   so that cap has no bar and no count.

Reads supplementary/TableS9_finngen_coloc.csv and figure_data/si/susie_L_sensitivity.csv; fits nothing.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402
import matplotlib.mathtext as _mt  # noqa: E402

# mathtext shrinks sub- and superscripts to 0.7 (4.2 pt at 6 pt); 0.84 x 6 pt = 5.04 pt
_mt.SHRINK_FACTOR = 0.84

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
S9 = os.path.join(SRC_DIR, "supplementary", "TableS9_finngen_coloc.csv")
SUSIE = os.path.join(SRC_DIR, "figure_data", "si", "susie_L_sensitivity.csv")

LOCI = ["PCSK9", "SORT1/CELSR2", "APOE/APOC1", "LDLR", "TRIB1"]
KEY = {"PCSK9": "PCSK9", "SORT1/CELSR2": "SORT1_CELSR2", "APOE/APOC1": "APOE_APOC1", "LDLR": "LDLR", "TRIB1": "TRIB1"}
CLS = ["robust", "robust", "robust", "prior-sensitive", "fails"]
DY = 0.15


def draw_a(ax):
    s9 = {r["locus"]: r for r in T.rows(S9)}
    ax.axvline(0.8, color=T.NULL, lw=T.LW["ref"], zorder=1)
    ax.text(0.79, -0.75, "PP.H4 = 0.8", ha="right", va="center", fontsize=T.FS["note"], color=T.MUTED)
    for i, loc in enumerate(LOCI):
        d, s = float(s9[KEY[loc]]["PPH4_default"]), float(s9[KEY[loc]]["PPH4_p12_1e6"])
        # the two priors sit just above and below the locus row, joined, so coinciding values stay visible
        ax.plot([d, s], [i - DY, i + DY], color=T.NS, lw=T.LW["link"], zorder=2)
        for v, y, mk in ((d, i - DY, "o"), (s, i + DY, "s")):
            if v > 0.8:
                ax.plot([v], [y], marker=mk, ms=T.MS + 0.4, mfc=T.INK, mec="white", mew=0.5, ls="none", zorder=4)
            else:
                ax.plot([v], [y], marker=mk, ms=T.MS + 0.1, mfc="white", mec=T.INK, mew=0.9, ls="none", zorder=4)
        ax.text(1.04, i, CLS[i], ha="left", va="center", fontsize=T.FS["note"], color=T.INK if CLS[i] == "robust" else T.MUTED)
        T.src("S5", "a", loc, "PP.H4 default; stringent", d, s, s)
    ax.text(1.04, -0.75, "3 of 5 loci robust to both priors", ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_yticks(range(len(LOCI)))
    ax.set_yticklabels(LOCI, fontstyle="italic")
    ax.set_ylim(len(LOCI) - 0.45, -1.2)
    ax.set_xlim(0, 1.32)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])
    ax.spines["bottom"].set_bounds(0, 1.0)
    ax.set_xlabel("FinnGen external PP.H4", x=1.0 / 1.32 / 2)
    T.clean(ax)


def draw_b(ax):
    rows = T.rows(SUSIE)
    Ls = [10, 15, 20]
    ncs = [max(int(r["n_cs"]) for r in rows if int(r["L"]) == L) for L in Ls]
    assert ncs[:2] == [10, 14], ncs
    # L = 20 over-parameterized the region and did not converge: no bar and no count, only the words
    for i in range(2):
        ax.bar(i, ncs[i], color=T.NS, width=0.56, lw=0, zorder=2)
        ax.text(i, ncs[i] / 2, str(ncs[i]), ha="center", va="center", fontsize=T.FS["label"], fontweight="bold", zorder=4)
        T.src("S5", "b", "L=%d" % Ls[i], "credible sets", ncs[i], ncs[i], ncs[i])
    ax.text(2, 4.0, "did not\nconverge", ha="center", va="center", fontsize=T.FS["note"], color=T.MUTED, linespacing=1.1)
    ax.plot(range(3), Ls, color=T.NULL, lw=T.LW["ref"], marker="o", ms=2.4, mfc=T.NULL, mec=T.NULL, zorder=3)
    ax.text(1.45, 19.3, "L cap", ha="right", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["L = %d" % L for L in Ls])
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(0, 22)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_ylabel(r"$\it{LDLR}$ credible sets recovered")
    T.clean(ax, left=True)
    ax.tick_params(axis="x", length=0)


def build():
    P = T.Page(168.0)
    P.panel(2, 10, "a", "External replication versus shared-variant prior")
    P.key(12, 22, [("o", "white", T.INK, "default prior"), ("s", "white", T.INK, "stringent prior ($p_{12}$ = 10$^{-6}$)"),
                   ("s", T.INK, "white", "filled: PP.H4 > 0.8")], gap=9)
    axa = P.ax(64, 34, 230, 108)
    draw_a(axa)
    x0 = 340
    P.text(x0, 10, "b", fontsize=T.FS["letter"], fontweight="bold")
    t = P.text(x0 + 9, 10, "LDLR", fontsize=T.FS["title"], fontweight="bold", fontstyle="italic")
    w = t.get_window_extent(P.fig.canvas.get_renderer()).width * 72 / P.fig.dpi
    P.text(x0 + 9 + w, 10, " allelic heterogeneity", fontsize=T.FS["title"], fontweight="bold")
    axb = P.ax(378, 30, 126, 112)
    draw_b(axb)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS5_coloc_finemap"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS5_coloc_finemap.pdf")
