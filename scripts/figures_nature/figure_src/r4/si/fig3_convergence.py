#!/usr/bin/env python3
"""Supplementary Figure (integrated evidence and recovery of established lipid genes), redrawn in the r4 theme.

a  evidence architecture of the 16 top-ranked genes (burden, colocalization, fine-mapping, drug/Mendelian status)
b  precision-recall of the full score against burden alone, bootstrap 95% bands
c  rank of each of the 27 canonical genes under burden alone and under the integrated score
d  AUROC against a label-permutation null

Same content and the same random streams as figure_src/make_fig4_convergence.py (default_rng(1) for the PR bands and
then the permutation null; default_rng(42) for the paired AP bootstrap), so every printed number is unchanged.
Reads source/supplementary (TableS8, TableS5, gene_ranking_scores_corrected.csv) and figure_data/si
(gene_ranking_scores.csv, ed7_gene_coords_GRCh37.csv; see PROVENANCE_C.md).
"""
import os
import sys

import numpy as np
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
from matplotlib.transforms import offset_copy
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
FD = os.path.join(SRC_DIR, "figure_data", "si")
BURDEN = T.NULL          # burden-only score: mid grey; the full / integrated score is ink

CATEG = {"LDLR": "established", "PCSK9": "established", "ANGPTL3": "established", "APOB": "established",
         "LPA": "established", "ABCA1": "established", "APOA5": "established", "APOC3": "established",
         "CETP": "established", "LCAT": "established", "APOA1": "established", "LPL": "established",
         "PLG": "tagging", "SLC22A2": "tagging", "ZNF229": "tagging", "ABO": "candidate", "PDE3B": "candidate"}
DRUGMEND = {"LDLR", "PCSK9", "ANGPTL3", "APOB", "LPA", "ABCA1", "APOA5", "APOC3", "CETP", "LCAT", "APOA1", "LPL"}


def sub(P, x, y, s):
    P.text(x + 9, y + 9.5, s, color=T.MUTED)


def rich(P, x, y, parts, **kw):
    """One line of text in runs, each (string, italic)."""
    rend = P.fig.canvas.get_renderer()
    for s, it in parts:
        t = P.text(x, y, s, fontstyle="italic" if it else "normal", **kw)
        x += t.get_window_extent(rend).width * 72 / P.fig.dpi
    return x


# ------------------------------------------------------------------ a, evidence matrix
def finemap_genes():
    """A resolved credible set (max PIP >= 0.9) within +/-500 kb of the GRCh37 gene body (as ED Fig 5b)."""
    coord = {r["gene"]: (int(r["chr"]), int(r["start"]), int(r["end"])) for r in T.rows(os.path.join(FD, "ed7_gene_coords_GRCh37.csv"))}
    cs = T.rows(os.path.join(SUP, "TableS5_finemapping.csv"))
    w = 500_000
    return {g for g, (c, st, en) in coord.items()
            if any(int(r["chr"]) == c and float(r["top_pip"]) >= 0.9 and (int(float(r["locus_bp"])) - w) <= en
                   and (int(float(r["locus_bp"])) + w) >= st for r in cs)}


def draw_a(ax):
    rows = T.rows(os.path.join(SUP, "TableS8_convergence_scores.csv"))[:16]
    fm = finemap_genes()
    cap = 60.0
    xs = dict(lip=0, cv=1, coloc=2, fm=3, drug=4, cat=5.0, rank=7.6)
    ng = len(rows)
    out = []
    for i, r in enumerate(rows):
        g, y = r["gene"], i
        if i % 2 == 0:
            T.band(ax, y - 0.5, y + 0.5)
        lb = min(float(r["lipid_burden_log10P"]), cap); ob = float(r["outcome_burden_log10P"])
        ax.scatter([xs["lip"]], [y], s=2 + 34 * lb / cap, c=T.INK, edgecolors="none", zorder=3)
        ax.scatter([xs["cv"]], [y], s=2 + 34 * min(ob, 13) / 13, c=T.INK, edgecolors="none", zorder=3)
        if r["coloc_locus"] == "1":
            ax.plot([xs["coloc"]], [y], ls="none", marker="s", ms=T.MS, mfc=T.INK, mec="none", zorder=3)
        if g in fm:
            ax.plot([xs["fm"]], [y], ls="none", marker="s", ms=T.MS, mfc=T.INK, mec="none", zorder=3)
        if g in DRUGMEND:
            ax.plot([xs["drug"]], [y], ls="none", marker="*", ms=T.MS + 1.2, mfc=T.INK, mec="none", zorder=3)
        cat = CATEG.get(g, "candidate")
        ax.text(xs["cat"], y, cat, ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED if cat == "tagging" else T.INK)
        ax.text(xs["rank"], y, str(i + 1), ha="center", va="center", fontsize=T.FS["note"])
        ax.text(-0.6, y, g, ha="right", va="center", fontsize=T.FS["note"], fontstyle="italic")
        out.append((g, r["coloc_locus"] == "1", g in fm, g in DRUGMEND, cat))
    heads = [("lip", "lipid\nburden"), ("cv", "CV\nburden"), ("coloc", "coloc"), ("fm", "fine-\nmap"), ("drug", "drug/\nMend.")]
    for k, lab in heads:
        ax.text(xs[k], -0.75, lab, ha="center", va="bottom", fontsize=T.FS["note"], color=T.MUTED, linespacing=1.0)
    ax.text(xs["cat"], -0.75, "category", ha="left", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    ax.text(xs["rank"], -0.75, "rank", ha="center", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_xlim(-2.6, 8.1); ax.set_ylim(ng - 0.5, -0.5)
    ax.axis("off")
    return out, sorted(fm)


# ------------------------------------------------------------------ b, precision-recall
GRID = np.linspace(0.02, 1.0, 60)


def pr_on_grid(yv, sv):
    o = np.argsort(-sv, kind="mergesort"); yy = yv[o]
    npos = yy.sum()
    if npos < 2:
        return np.full(len(GRID), np.nan)
    tp = np.cumsum(yy); rec = tp / npos; prec = tp / np.arange(1, len(yy) + 1)
    m = yy == 1
    return np.interp(GRID, rec[m], prec[m], left=prec[m][0], right=prec[m][-1])


def draw_b(ax, y, si, sb, rng):
    n = len(y)

    def band(sv, B=200):
        mat = np.empty((B, len(GRID)))
        for b in range(B):
            idx = rng.integers(0, n, n)
            mat[b] = pr_on_grid(y[idx], sv[idx])
        return np.nanpercentile(mat, 2.5, axis=0), np.nanpercentile(mat, 97.5, axis=0)

    ap = {}
    for key, sv, col in (("burden", sb, BURDEN), ("full", si, T.INK)):
        ap[key] = average_precision_score(y, sv)
        lo, hi = band(sv)
        if key == "full":
            ax.fill_between(GRID, lo, hi, color=T.tint(T.INK, 0.14), lw=0, zorder=1)
        else:
            for b in (lo, hi):
                ax.plot(GRID, b, color=BURDEN, lw=0.5, ls=(0, (1.2, 1.2)), zorder=2)
        ax.plot(GRID, pr_on_grid(y, sv), color=col, lw=1.2, zorder=3 if key == "full" else 2.5)
    p_at = lambda sv, k: y[np.argsort(-sv, kind="mergesort")[:k]].sum() / k
    pk = {k: (p_at(si, k), p_at(sb, k)) for k in (10, 20, 50)}
    dap = average_precision_score(y, si) - average_precision_score(y, sb)
    r42 = np.random.default_rng(42)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    d = np.array([(lambda ix: average_precision_score(y[ix], si[ix]) - average_precision_score(y[ix], sb[ix]))
                  (np.concatenate([r42.choice(pos, pos.size), r42.choice(neg, neg.size)])) for _ in range(1000)])
    lo, hi = np.percentile(d, [2.5, 97.5])
    p = min(2 * min((d <= 0).mean(), (d >= 0).mean()), 1.0)
    au = (roc_auc_score(y, si), roc_auc_score(y, sb))
    notes = [("precision@k (full / burden):  @10  %.2f / %.2f  (tie-determined)" % pk[10], T.INK, T.FS["note"]),
             ("@20  %.2f / %.2f;  @50  %.2f / %.2f" % (pk[20] + pk[50]), T.INK, T.FS["note"]),
             (("$\\Delta$AP = %+.2f, bootstrap 95%% CI [%+.2f, %+.2f] includes 0 ($P$ = %.2f):" % (dap, lo, hi, p)).replace("-", "−"), T.INK, T.FS["note"]),
             ("not a robust improvement", T.INK, T.FS["note"]),
             ("AUROC %.2f vs %.2f (secondary; over-reassuring under class imbalance)" % au, T.MUTED, T.FS["small"])]
    ax.set_xlabel("Recall of established genes"); ax.set_ylabel("Precision")
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.03)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    T.clean(ax, left=True)
    return ap, pk, (dap, lo, hi, p), au, notes


# ------------------------------------------------------------------ c, per-gene ranks
def draw_c(ax, gene, y, si, sb):
    n = len(y)
    order = np.argsort(-si, kind="mergesort"); ri = np.empty(n, int); ri[order] = np.arange(1, n + 1)
    order = np.argsort(-sb, kind="mergesort"); rb = np.empty(n, int); rb[order] = np.arange(1, n + 1)
    pos = np.where(y == 1)[0]
    pos = pos[np.argsort(ri[pos])]
    gw = {r["gene"]: r for r in T.rows(os.path.join(SUP, "gene_ranking_scores_corrected.csv"))}
    rescue = {g for g in ("SORT1", "PSRC1", "CELSR2", "APOE", "APOC2") if gw.get(g, {}).get("coloc_w500") == "1"}
    assert len(rescue) == 5, sorted(rescue)
    miss = {"HMGCR", "MTTP"}
    npos = len(pos)
    ax.axvline(50, color=T.MUTED, ls=(0, (1.5, 1.5)), lw=T.LW["ref"], zorder=1)
    ax.text(50 * 1.12, npos - 0.1, "top 50", ha="left", va="center", fontsize=T.FS["small"], color=T.MUTED)
    out = []
    for j, idx in enumerate(pos):
        yy = j                      # as before: best integrated rank at the bottom
        g = gene[idx]; a, b = rb[idx], ri[idx]
        if g in rescue:
            ax.plot([a, b], [yy, yy], color=T.INK, lw=1.3, solid_capstyle="butt", zorder=2)
        else:
            ax.plot([a, b], [yy, yy], color=T.NS, lw=0.6, solid_capstyle="butt", zorder=2)
        ax.plot([a], [yy], ls="none", marker="o", ms=2.6, mfc=BURDEN, mec="white", mew=0.3, zorder=3)
        if g in miss:
            ax.plot([b], [yy], ls="none", marker="o", ms=3.0, mfc="white", mec=T.INK, mew=0.8, zorder=4)
        else:
            ax.plot([b], [yy], ls="none", marker="o", ms=3.2, mfc=T.INK, mec="white", mew=0.4, zorder=4)
        ax.text(0.62, yy, g, ha="right", va="center", fontsize=T.FS["small"], fontstyle="italic")
        out.append((g, int(a), int(b)))
    ax.set_xscale("log"); ax.set_xlim(0.7, 25000); ax.set_ylim(-0.7, npos + 0.3)
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000]))
    ax.xaxis.set_major_formatter(FixedFormatter(["1", "10", "100", "1,000", "10,000"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Genome-wide rank (log; left = better)")
    ax.set_yticks([])
    T.clean(ax)
    return out, sorted(rescue), sorted(miss)


# ------------------------------------------------------------------ d, label-permutation null
def draw_d(ax, y, si, rng):
    B = 10000
    null = np.empty(B); obs = roc_auc_score(y, si)
    for b in range(B):
        null[b] = roc_auc_score(rng.permutation(y), si)
    ax.hist(null, bins=40, color=T.NS, edgecolor="white", lw=0.3, zorder=2)
    ax.axvline(0.5, color=T.NULL, lw=T.LW["ref"], zorder=3)      # the null centre, drawn over the bars
    ax.axvline(obs, color=T.INK, lw=1.4, zorder=4)
    ymax = ax.get_ylim()[1]
    ax.annotate("observed\nAUROC %.2f" % obs, xy=(obs, ymax * 0.5), xytext=(obs - 0.13, ymax * 0.68), fontsize=T.FS["note"],
                fontweight="bold", ha="center", va="center", linespacing=1.05,
                arrowprops=dict(arrowstyle="-|>", color=T.INK, lw=0.6, mutation_scale=6, shrinkA=6, shrinkB=1))
    ax.text(0.665, ymax * 0.97, "grey = 10,000 random\ngene-label shuffles\n(null centred at 0.5)", fontsize=T.FS["note"],
            color=T.MUTED, va="top", linespacing=1.1)
    ax.text(0.665, ymax * 0.40, "no shuffle reached the\nobserved value (empirical", fontsize=T.FS["note"], va="top", linespacing=1.1)
    # the third line carries a power of ten: the exponent is set at 5 pt by hand, since mathtext would shrink it to 4.2 pt
    rend = ax.figure.canvas.get_renderer()
    x_pt = 0.0
    for txt, fs, rise in (("$P$ = 1.0×10", T.FS["note"], 0.0), ("−4", 5.0, 2.6), ("; AP likewise)", T.FS["note"], 0.0)):
        t = ax.text(0.665, ymax * 0.40, txt, fontsize=fs, va="baseline",
                    transform=offset_copy(ax.transData, fig=ax.figure, x=x_pt, y=-2 * 6.6 - 5.0 + rise, units="points"))
        x_pt += t.get_window_extent(rend).width * 72 / ax.figure.dpi - (0.9 if txt.endswith("10") else 0.0)
    ax.set_xlabel("AUROC under permuted gene labels"); ax.set_ylabel("Permutations")
    ax.set_xlim(0.42, 1.0); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    T.clean(ax, left=True)
    return obs, (null.min(), null.max(), float((null >= obs).mean()))


def build():
    grows = T.rows(os.path.join(FD, "gene_ranking_scores.csv"))
    gene = np.array([r["gene"] for r in grows])
    y = np.array([int(float(r["canonical"])) for r in grows])
    si = np.array([float(r["integrated"]) for r in grows])
    sb = np.array([float(r["burden_only"]) for r in grows])
    assert int(y.sum()) == 27 and len(y) == 18002
    rng = np.random.default_rng(1)

    P = T.Page(470.0)
    P.panel(2, 10, "a", "Top-ranked gene evidence architecture")
    sub(P, 2, 10, "dot size = burden −log10 P; square = coloc/fine-map; star = drug target or Mendelian")
    axa = P.ax(10, 42, 236, 150); a_out, fm = draw_a(axa)

    P.panel(262, 10, "b", "Precision–recall (primary)")
    sub(P, 262, 10, "27 established genes vs 18,002; biological calibration, not prospective validation")
    ap0 = (average_precision_score(y, si), average_precision_score(y, sb))
    P.key(272, 31, [("line", None, T.INK, "full score (AP %.2f)" % ap0[0]), ("line", None, BURDEN, "burden only (AP %.2f)" % ap0[1])], gap=8)
    P.text(272, 40.5, "bootstrap 95% CI: shaded (full score), dotted (burden only)", color=T.MUTED)
    axb = P.ax(298, 50, 200, 100); b_out = draw_b(axb, y, si, sb, rng)
    yy = 177.0
    for k, (txt, col, fs) in enumerate(b_out[4]):
        P.text(272, yy, txt, color=col, fontsize=fs)
        yy += 7.6 if k in (0, 2) else 9.0

    P.panel(2, 222, "c", "Rank of each canonical gene: burden vs integrated")
    x = rich(P, 11, 231.5, [("integrated rescues ", False), ("APOC2/APOE/CELSR2/PSRC1/SORT1", True), (" (no burden signal);", False)], color=T.MUTED)
    rich(P, 11, 239.5, [("misses ", False), ("HMGCR", True), (", ", False), ("MTTP", True)], color=T.MUTED)
    P.key(12, 252, [("o", BURDEN, "white", "burden rank"), ("o", T.INK, "white", "integrated rank"), ("line", None, T.INK, "rescued (regulatory)"),
                    ("o", "white", T.INK, "missed")], gap=7)
    axc = P.ax(40, 260, 206, 182); c_out = draw_c(axc, gene, y, si, sb)

    P.panel(262, 222, "d", "Enrichment vs a label-permutation null")
    sub(P, 262, 222, "observed ranking (black line) far exceeds random-label rankings (grey)")
    axd = P.ax(298, 248, 200, 194); d_out = draw_d(axd, y, si, rng)

    P.save(os.path.join(SRC_DIR, "figures", "fig3_convergence"))
    print("a:", a_out, "finemap", fm)
    print("b: AP", b_out[0], "p@k", b_out[1], "dAP", b_out[2], "AUROC", b_out[3])
    print("c:", c_out)
    print("d:", d_out)


if __name__ == "__main__":
    build()
