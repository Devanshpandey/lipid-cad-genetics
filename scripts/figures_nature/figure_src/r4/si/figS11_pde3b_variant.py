#!/usr/bin/env python3
"""Supplementary Figure (PDE3B variant-level burden architecture), redrawn in the r4 theme.

a  leave-one-variant-out: the pLoF-burden triglyceride effect with each qualifying variant removed in turn
b  single-variant triglyceride effects (95% CI from P) for the qualifying variants with >= 10 carriers

Same content as figure_src/make_figS11.py. Reads figure_data/si/pde3b_lovo.csv (see PROVENANCE_C.md).
"""
import math
import os
import sys

from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FD = os.path.join(SRC_DIR, "figure_data", "si")
TG = T.TRAIT["TG"]
FULL = -0.378


def load():
    rows = T.rows(os.path.join(FD, "pde3b_lovo.csv"))
    for r in rows:
        r["n_carriers"] = int(r["n_carriers"])
        for k in ("drop_beta", "single_beta", "single_P"):
            r[k] = T.num(r[k])
    rows.sort(key=lambda r: -r["n_carriers"])
    return rows


def draw_a(ax, rows):
    db = [r["drop_beta"] for r in rows]
    ax.axhspan(min(db), max(db), color=T.tint(TG, 0.16), lw=0, zorder=0)
    ax.axhline(FULL, color=T.INK, lw=T.LW["ref"], zorder=2)
    ax.text(len(rows) + 3.5, FULL, "full burden\n(%s)" % ("%.2f" % FULL).replace("-", "−"), ha="left", va="center",
            fontsize=T.FS["note"], linespacing=1.05)
    ax.plot(range(len(rows)), db, ls="none", marker="o", ms=2.6, mfc=TG, mec="white", mew=0.3, alpha=0.9, zorder=3)
    ax.plot([0], [db[0]], ls="none", marker="o", ms=T.MS + 0.6, mfc=TG, mec=T.INK, mew=0.7, zorder=4)
    ax.annotate("drop most common variant\n(780 carriers): %s" % ("%.2f" % db[0]).replace("-", "−"), (0, db[0]),
                textcoords="offset points", xytext=(9, 0), fontsize=T.FS["note"], va="center", linespacing=1.05)
    ax.set_xlim(-3, len(rows) + 2); ax.set_ylim(-0.42, -0.30)
    ax.set_xticks([0, 20, 40, 60, 80]); ax.set_yticks([-0.42, -0.40, -0.38, -0.36, -0.34, -0.32, -0.30])
    ax.set_xlabel("Qualifying pLoF variant removed (ranked by carrier count)")
    ax.set_ylabel("Triglyceride burden effect (SD)")
    T.clean(ax, left=True)
    return min(db), max(db), db[0], len(rows)


def draw_b(ax, rows):
    big = [r for r in rows if r["n_carriers"] >= 10 and not math.isnan(r["single_beta"]) and not math.isnan(r["single_P"])]
    big.sort(key=lambda r: r["n_carriers"])
    T.null(ax, 0.0)
    ax.axvline(FULL, color=T.MUTED, ls=(0, (1.5, 1.5)), lw=T.LW["ref"], zorder=1)
    ax.text(FULL, len(big) - 0.25, "full burden", ha="center", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    out = []
    for i, r in enumerate(big):
        b, p, n = r["single_beta"], r["single_P"], r["n_carriers"]
        se = abs(b) / norm.ppf(1 - p / 2) if p < 1 else abs(b)
        lo, hi = b - 1.96 * se, b + 1.96 * se
        T.point(ax, i, b, lo, hi, TG, "o", filled=p < 0.05)
        ax.text(1.02, i, "n = %d" % n, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
        out.append((r["variant"], n, round(b, 3), round(lo, 3), round(hi, 3), p < 0.05))
    ax.set_yticks(range(len(big))); ax.set_yticklabels([r["variant"].split(":")[1] for r in big])
    ax.set_ylim(-0.7, len(big) - 0.3)
    ax.set_xlim(-1.4, 0.65); ax.set_xticks([-1.0, -0.5, 0, 0.5])
    ax.set_xlabel("Single-variant triglyceride effect (SD, 95% CI)")
    ax.set_ylabel("Variant (chr11 position)")
    T.clean(ax)
    return out


def build():
    rows = load()
    P = T.Page(190.0)
    P.panel(2, 10, "a", "Leave-one-variant-out is stable")
    P.text(11, 19.5, "no single variant drives the burden signal", color=T.MUTED)
    axa = P.ax(40, 32, 174, 130); a_out = draw_a(axa, rows)
    P.panel(262, 10, "b", "Individual qualifying variants")
    P.text(271, 19.5, "variants with ≥10 carriers; filled = $P$ < 0.05", color=T.MUTED)
    axb = P.ax(326, 36, 150, 126); b_out = draw_b(axb, rows)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS11_pde3b_variant"))
    print("a: range %.3f to %.3f, drop-most-common %.3f, %d variants" % a_out)
    for r in b_out:
        print("b:", r)


if __name__ == "__main__":
    build()
