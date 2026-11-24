#!/usr/bin/env python3
"""Supplementary Figure (PDE3B carrier survival diagnostics), redrawn in the r4 theme.

a  adjusted cumulative incidence of CAD in pLoF carriers and non-carriers, with the carrier number at risk
b  power to detect a carrier CAD effect given 92 events

Same content as figure_src/make_figS12.py. Reads figure_data/si/pde3b_km.csv and pde3b_power.csv (see PROVENANCE_C.md).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FD = os.path.join(SRC_DIR, "figure_data", "si")
NONCARRIER = T.NULL       # mid grey; carriers are ink


def draw_a(P, ax, x0, w, ybase):
    km = T.rows(os.path.join(FD, "pde3b_km.csv"))
    for g, col in (("0", NONCARRIER), ("1", T.INK)):
        t = [float(x["time"]) for x in km if x["carrier"] == g]
        c = [float(x["cuminc"]) * 100 for x in km if x["carrier"] == g]
        ax.step(t, c, where="post", color=col, lw=1.1, zorder=3 if g == "1" else 2)
    ax.set_xlim(0, 14); ax.set_ylim(0, 10)
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12, 14]); ax.set_yticks([0, 2, 4, 6, 8, 10])
    ax.set_xlabel("Years from baseline"); ax.set_ylabel("Cumulative incidence of CAD (%)")
    T.clean(ax, left=True)
    nr = {float(x["time"]): int(x["n_risk"]) for x in km if x["carrier"] == "1"}
    ticks = [0, 4, 8, 12]
    P.text(x0 - 34, ybase, "Carriers at risk", color=T.MUTED)
    for t in ticks:
        P.text(x0 + w * t / 14.0, ybase + 8, format(nr.get(float(t)), ","), ha="center", color=T.MUTED)
    return {t: nr.get(float(t)) for t in ticks}


def draw_b(ax):
    pw = T.rows(os.path.join(FD, "pde3b_power.csv"))
    hr = [float(x["HR"]) for x in pw]; power = [float(x["power"]) for x in pw]
    ax.axvspan(0.745, 1.342, color=T.BAND, lw=0, zorder=0)
    ax.axhline(0.8, color=T.MUTED, ls=(0, (1.5, 1.5)), lw=T.LW["ref"], zorder=1)
    ax.text(0.515, 0.82, "80% power", ha="left", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    ax.text(1.0, 0.99, "under-powered zone\n(HR 0.75–1.34)", ha="center", va="top", fontsize=T.FS["note"], color=T.MUTED,
            linespacing=1.05)
    ax.axvline(1.028, ymax=0.76, color=T.INK, ls=(0, (3, 2)), lw=T.LW["ref"], zorder=2)
    ax.text(1.045, 0.55, "observed\nHR 1.03", ha="left", va="center", fontsize=T.FS["note"], linespacing=1.05)
    ax.plot(hr, power, color=T.INK, lw=1.2, zorder=3)
    ax.text(0.515, 0.04, "proportional hazards held\n(Schoenfeld carrier $P$ = 0.10)", ha="left", va="bottom", fontsize=T.FS["note"],
            linespacing=1.1)
    T.logx(ax, 0.5, 1.6, [0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.6])
    ax.set_ylim(0, 1.02); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("True carrier hazard ratio for CAD"); ax.set_ylabel("Power (92 events)")
    T.clean(ax, left=True)
    lo = max(h for h, p in zip(hr, power) if h < 1 and p >= 0.8)
    hi = min(h for h, p in zip(hr, power) if h > 1 and p >= 0.8)
    return len(pw), lo, hi


def build():
    P = T.Page(196.0)
    P.panel(2, 10, "a", "Carrier vs non-carrier CAD incidence")
    P.text(11, 19.5, "curves nearly identical (Cox HR 1.03, 92 carrier events)", color=T.MUTED)
    P.key(12, 32, [("line", None, NONCARRIER, "non-carrier"), ("line", None, T.INK, "pLoF carrier")], gap=8)
    axa = P.ax(46, 42, 190, 112); a_out = draw_a(P, axa, 46, 190, 182)
    P.panel(262, 10, "b", "Underpowered to exclude a modest effect")
    P.text(271, 19.5, "80% power only for HR ≤0.75 or ≥1.34", color=T.MUTED)
    axb = P.ax(306, 42, 194, 112); b_out = draw_b(axb)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS12_pde3b_survival"))
    print("a: carriers at risk", a_out)
    print("b: %d grid points; power >= 0.8 for HR <= %.3f and >= %.3f" % b_out)


if __name__ == "__main__":
    build()
