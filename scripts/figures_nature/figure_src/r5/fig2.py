#!/usr/bin/env python3
"""Fig. 2 (r5): translating the association into a trial effect, and conditional forecasts for the trials in progress.

a  the positive control: LDL-C score per mmol/L under four scalings of treated LDL-C, against CTT   (r4 Fig. 2c, top)
b  Lp(a) predicted for Lp(a)HORIZON with one or two slopes, against the trial                         (r4 Fig. 2c, rest)
c  HR across bands of predicted Lp(a), with the number of established-disease participants per band  (from the r4
   dose Supplementary Figure, with counts added)
d  forecasts as conditional scenarios: the primary ratio first, then each forecast beside its assumption (r4 Fig. 2d
   with the primary ratio and the assumptions added)

Reads only evidence/lpa_emulation; fits nothing.
"""
import os

from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

import common as C
from common import T
import fig2 as F2   # r4

LP, LD = T.TRAIT["Lp(a)"], T.TRAIT["LDL-C"]
N = lambda x: format(int(x), ",")


# ------------------------------------------------------------------ c: bands of predicted Lp(a)
def draw_c(ax, D):
    ax.axvspan(150, 320, color=T.BAND, lw=0, zorder=0)
    T.null(ax, 1.0, axis="y")
    for cohort, model, marker, key, dx in (("C standardized", "pre-registered", "o", "C_std", -3), ("A", "pre-registered (pooled)", "s", "A", 3)):
        bd = [r for r in D["bands"] if r["cohort"] == cohort and r["model"] == model and r["endpoint"] == "mace4"]
        assert [r["band"] for r in bd] == ["<50", "50-99", "100-149", ">=150"], [r["band"] for r in bd]
        xs = [T.num(r["predicted_mean"]) + dx for r in bd]
        for r, x in zip(bd, xs):
            hr, lo, hi = T.num(r["hr"]), T.num(r["hr_lo"]), T.num(r["hr_hi"])
            if r["band"] != "<50":
                ax.plot([x, x], [lo, hi], color=LP, lw=T.LW["ci"], solid_capstyle="butt", zorder=3)
            ax.plot([x], [hr], marker=marker, ms=T.MS, mfc=LP, mec="white", mew=0.5, ls="none", zorder=4)
            if key == "A":
                first = r["band"] == "<50"
                ax.text(4 if first else x, 0.925, ("n = " if first else "") + N(r["n"]), ha="left" if first else "center", va="center",
                        fontsize=T.FS["small"], color=T.MUTED)
            T.src("2", "c", "band %s (mean %.1f nmol/L)" % (r["band"], x - dx), key, hr, lo, hi, n=r["n"], events=r["events"])
    ax.text(316, 1.03, "above 150 nmol/L", fontsize=T.FS["note"], color=T.MUTED, ha="right", va="bottom")
    yb = 1.84
    ax.plot([F2.TRIAL_BASELINE[1], F2.TRIAL_BASELINE[2]], [yb, yb], color=T.TRIAL, lw=T.LW["ci"], solid_capstyle="butt", zorder=3)
    ax.plot([F2.TRIAL_BASELINE[0]], [yb], marker="D", ms=T.MS - 0.4, mfc=T.TRIAL, mec="white", mew=0.5, ls="none", zorder=4)
    ax.text(316, yb * 0.975, "Lp(a)HORIZON baseline, median and IQR", fontsize=T.FS["small"], color=T.MUTED, ha="right", va="top")
    ax.set_yscale("log"); ax.set_ylim(0.9, 1.95)
    ax.yaxis.set_major_locator(FixedLocator([1, 1.25, 1.5, 1.75])); ax.yaxis.set_major_formatter(FixedFormatter(["1", "1.25", "1.5", "1.75"]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0, 320); ax.set_xticks([0, 50, 100, 150, 200, 250, 300])
    ax.set_xlabel("Genetically predicted Lp(a), nmol/L (band mean)")
    ax.set_ylabel("HR against the band below 50 nmol/L")
    T.clean(ax, left=True)


# ------------------------------------------------------------------ d: forecasts as scenarios
FORECASTS = (("ACCLAIM-Lp(a), at risk", "acc_risk", "mace4", "o", "assumes 90% lowering from 250 nmol/L"),
             ("ACCLAIM-Lp(a), established", "acc_est", "mace4", "s", "assumes 90% lowering from 250 nmol/L"),
             ("OCEAN(a)-Outcomes", "ocean", "cor3", "s", "assumes 90% lowering from 275 nmol/L;\nonly 309 of 12,128 above its 200 nmol/L entry"))


def draw_d(P, ax, D, x_text, y_text):
    r = C.one(C.rows(os.path.join(C.EV, "forecast", "acclaim_stage_ratio.csv")), at_risk_cohort="acc_risk", endpoint="mace4", model="linear", scale="lifelong")
    ratio = (T.num(r["ratio"]), T.num(r["lo"]), T.num(r["hi"]))
    P.text(x_text, y_text, "Primary forecast, ACCLAIM-Lp(a): established-disease log HR", fontsize=T.FS["note"], fontweight="bold")
    P.text(x_text, y_text + 8.5, "as a share of the at-risk one, %s" % T.ci_text(*ratio), fontsize=T.FS["note"], fontweight="bold")
    P.text(x_text, y_text + 17, "1 = the same benefit after onset as before it", fontsize=T.FS["small"], color=T.MUTED)
    T.src("2", "d", "primary forecast: stage ratio of log HRs", "ACCLAIM-Lp(a)", *ratio)
    T.null(ax)
    y, yt, yl = 0.0, [], []
    for k, (head, cohort, endpoint, marker, note) in enumerate(FORECASTS):
        nl = note.count("\n") + 1
        span = 0.45 + nl * 0.72 + 0.55 + 1.0
        if k % 2 == 0:
            T.band(ax, y - 0.55, y + span + 0.45)
        ax.text(-0.02, y, head, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=T.FS["note"], fontweight="bold")
        ax.text(-0.02, y + 0.45, note, transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=T.FS["small"],
                color=T.MUTED, linespacing=1.05)
        y += 0.45 + nl * 0.72 + 0.55
        for model, lab in (("linear", "one slope"), ("spline", "second slope")):
            f = C.one(D["fc"], cohort=cohort, endpoint=endpoint, model=model)
            e, lo, hi = T.num(f["hr"]), T.num(f["lo"]), T.num(f["hi"])
            F2.mark(ax, y, e, lo, hi, LP, marker)
            F2.num_col(ax, 1.03, y, T.ci_text(e, lo, hi))
            yt.append(y); yl.append(lab)
            T.src("2", "d", "%s: %s" % (head, lab), marker, e, lo, hi)
            y += 1.0
    F2.num_col(ax, 1.03, -1.0, "HR (95% CI)", color=T.MUTED)
    ax.set_yticks(yt); ax.set_yticklabels(yl)
    ax.set_ylim(y - 0.4, -1.5)
    T.logx(ax, 0.5, 1.1, [0.5, 0.6, 0.7, 0.8, 1])
    ax.set_xlabel("HR for the assumed lowering"); T.clean(ax)


def build():
    C.relabel("2")
    D = F2.load()
    groups = F2.groups_c(D)
    P = T.Page(412.0)
    F2.stage_key(P, 10, 9, followup=False)
    P.key(10, 20, [("line", None, LD, "LDL-C score (positive control)"), ("line", None, LP, "Lp(a)")], gap=9)
    P.panel(2, 38, "a", "The positive control depends on the scaling")
    axa = P.ax(124, 52, 66, 80)
    F2.grouped(axa, groups[:1], "a", (0.46, 1.1), [0.5, 0.6, 0.7, 0.8, 1], "HR per mmol/L lower LDL-C")
    P.panel(262, 38, "b", "The predicted trial effect depends on the dose model")
    axb = P.ax(378, 52, 66, 118)
    F2.grouped(axb, groups[1:], "b", (0.46, 1.42), [0.5, 0.7, 1, 1.4], "HR for Lp(a)HORIZON")
    P.panel(2, 206, "c", "Few participants lie in the trial's range")
    axc = P.ax(44, 222, 196, 140); draw_c(axc, D)
    P.panel(262, 206, "d", "Conditional forecasts for the trials in progress")
    axd = P.ax(378, 256, 66, 124); draw_d(P, axd, D, 272, 222)
    P.save(os.path.join(C.SRC_DIR, "figures", "fig2_lpa_translation"))
    C.write("Fig2_source_data.csv")


if __name__ == "__main__":
    build()
    print("wrote Fig. 2 (r5)")
