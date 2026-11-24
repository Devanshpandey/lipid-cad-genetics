#!/usr/bin/env python3
"""Fig. 2: genetically predicted Lp(a) before and after disease onset, beside randomized trials;
and the Supplementary Figure on dose-response and per-144-nmol/L estimates.

a  stage contrast for Lp(a), the LDL-C score (positive control) and the non-lipid score     (fixed in advance)
b  Lp(a) by endpoint at two stages, beside the Lp(a)HORIZON components                      (fixed in advance)
c  translating the associations into trial effects                                          (post hoc)
d  conditional forecasts for the two Lp(a) trials in progress                               (committed before fitting)

Reads only evidence/lpa_emulation (copies of the LS6 result tables); fits nothing. Shape is the stage of the
population an estimate describes (circle, before a first event; square, established disease; triangle, after a first
event during follow-up); a diamond is a randomized trial; open markers have a 95% CI that includes 1.
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theme as T  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
EV = os.path.abspath(os.path.join(SRC_DIR, "..", "evidence", "lpa_emulation"))

# published trial results: Cho et al., N Engl J Med 2026, Table 2
TRIAL = {"mace4": (1.02, 0.90, 1.15), "mace3": (1.01, 0.88, 1.16), "revasc": (1.04, 0.80, 1.36),
         "mi": (0.99, 0.82, 1.18), "cvdeath": (0.93, 0.73, 1.19), "stroke": (1.06, 0.82, 1.37)}
TRIAL_BASELINE = (236.0, 204.0, 308.0)
TRIAL_ENTRY = 150.0
CTT = (0.78, 0.76, 0.80)
ENDPOINTS = [("mace4", "Four-point composite"), ("mace3", "Three-point composite"), ("revasc", "Revascularization"),
             ("mi", "Myocardial infarction"), ("cvdeath", "Cardiovascular death"), ("stroke", "Ischaemic stroke")]
STAGE = {"C_std": "o", "A": "s", "D": "^"}
HUE = {"x_lpa": T.TRAIT["Lp(a)"], "x_ldl": T.TRAIT["LDL-C"], "x_prs": T.TRAIT["non-lipid"]}


def one(table, **kw):
    r = [x for x in table if all(x[k] == v for k, v in kw.items())]
    assert len(r) == 1, (kw, len(r))
    return r[0]


def load():
    est = T.rows(os.path.join(EV, "emulation", "estimates.csv"))

    def pick(cohort, endpoint, exposure, horizon="full", model="joint"):
        r = one(est, cohort=cohort, endpoint=endpoint, exposure=exposure, horizon=horizon, model=model)
        return T.num(r["b"]), T.num(r["se"]), int(r["n"]), int(r["events"])
    D = dict(pick=pick)
    D["contrast"] = {r["exposure"]: r for r in T.rows(os.path.join(EV, "emulation", "stage_contrast.csv"))}
    D["cal2"] = T.rows(os.path.join(EV, "posthoc2", "calibration_statin.csv"))
    D["bands"] = T.rows(os.path.join(EV, "posthoc2", "linearity_bands.csv"))
    D["shape"] = T.rows(os.path.join(EV, "posthoc2", "linearity_shape.csv"))
    D["fc"] = [r for r in T.rows(os.path.join(EV, "forecast", "forecasts.csv")) if r["central"] == "TRUE" and r["scale"] == "lifelong"]
    D["mp"] = T.rows(os.path.join(EV, "matched_endpoint", "matched_predictions.csv"))
    return D


def mark(ax, y, est, lo, hi, colour, marker):
    T.point(ax, y, est, lo, hi, colour, marker, filled=T.excludes(lo, hi), ms=T.MS - 0.4 if marker == "D" else T.MS)


def num_col(ax, x, y, s, **kw):
    ax.text(x, y, s, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], **kw)


# ------------------------------------------------------------------ a
def draw_a(ax, D):
    groups = [("x_lpa", "Lp(a)", "per 144 nmol/L lower"), ("x_ldl", "LDL-C score", "per SD lower"), ("x_prs", "Non-lipid score", "per SD lower")]
    T.null(ax)
    y, ticks, labels = 0.0, [], []
    for k, (x, name, unit) in enumerate(groups):
        if k % 2 == 0:
            T.band(ax, y - 0.5, y + 2 * 0.8 + 0.5)
        for j, cohort in enumerate(("C_std", "A", "D")):
            b, s, n, ev = D["pick"](cohort, "mace4", x)
            hr, lo, hi = T.lower(b, s)
            yy = y + j * 0.8
            mark(ax, yy, hr, lo, hi, HUE[x], STAGE[cohort])
            num_col(ax, 1.03, yy, T.ci_text(hr, lo, hi))
            T.src("2", "a", name, cohort, hr, lo, hi, n=n, events=ev)
        c = D["contrast"][x]
        num_col(ax, 1.72, y + 0.8, "%.0f%% (%.0f–%.0f)" % (100 * T.num(c["ratio"]), 100 * T.num(c["ratio_lo"]), 100 * T.num(c["ratio_hi"])))
        ticks.append(y + 0.8); labels.append("%s\n%s" % (name, unit))
        y += 3 * 0.8 + 0.9
    num_col(ax, 1.03, -1.05, "HR (95% CI)", color=T.MUTED)
    num_col(ax, 1.72, -1.05, "Retained", color=T.MUTED)
    ax.set_yticks(ticks); ax.set_yticklabels(labels, linespacing=1.15)
    ax.set_ylim(y - 0.9 - 0.8 + 0.6, -1.4)
    T.logx(ax, 0.66, 1.04, [0.7, 0.8, 0.9, 1])
    ax.set_xlabel("HR, four-point composite"); T.clean(ax)


# ------------------------------------------------------------------ b
def draw_b(ax, D):
    T.null(ax)
    y, ticks = 0.0, []
    for k, (e, lab) in enumerate(ENDPOINTS):
        if k % 2 == 0:
            T.band(ax, y - 0.45, y + 1.3 + 0.45)
        for j, cohort in enumerate(("C_std", "A")):
            b, s, n, ev = D["pick"](cohort, e, "x_lpa")
            hr, lo, hi = T.lower(b, s)
            mark(ax, y + j * 0.65, hr, lo, hi, HUE["x_lpa"], STAGE[cohort])
            T.src("2", "b", lab, cohort, hr, lo, hi, events=ev)
        tr = TRIAL[e]
        mark(ax, y + 1.3, *tr, T.TRIAL, "D")
        T.src("2", "b", lab, "Lp(a)HORIZON", *tr)
        ticks.append(y + 0.65)
        y += 1.3 + 1.0
    ax.set_yticks(ticks); ax.set_yticklabels([l for _, l in ENDPOINTS])
    ax.set_ylim(y - 1.0 + 0.5, -0.5)
    T.logx(ax, 0.55, 1.45, [0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4])
    ax.set_xlabel("HR per 144 nmol/L lower Lp(a)"); T.clean(ax)


# ------------------------------------------------------------------ c and d: grouped forests with a number column
def grouped(ax, groups, panel, xlim, ticks, xlabel, col_x=1.03):
    T.null(ax)
    y, yt, yl = 0.0, [], []
    for k, (head, items) in enumerate(groups):
        if k % 2 == 0:
            T.band(ax, y - 0.55, y + len(items) + 0.45)
        ax.text(-0.02, y, head, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=T.FS["note"], fontweight="bold")
        y += 1.0
        for lab, est, lo, hi, colour, marker in items:
            mark(ax, y, est, lo, hi, colour, marker)
            num_col(ax, col_x, y, T.ci_text(est, lo, hi))
            yt.append(y); yl.append(lab)
            T.src("2", panel, "%s: %s" % (head, lab), marker, est, lo, hi)
            y += 1.0
    num_col(ax, col_x, -1.0, "HR (95% CI)", color=T.MUTED)
    ax.set_yticks(yt); ax.set_yticklabels(yl, color=T.INK)
    ax.set_ylim(y - 0.4, -1.5)
    T.logx(ax, xlim[0], xlim[1], ticks)
    ax.set_xlabel(xlabel); T.clean(ax)


def groups_c(D):
    def ldl(cohort, model, scale):
        r = one(D["cal2"], cohort=cohort, horizon="full", model=model, ldl_scale=scale)
        return 1 / T.num(r["genetic_hr_per_mmol_ldl"]), 1 / T.num(r["ldl_hi"]), 1 / T.num(r["ldl_lo"])
    A = lambda s: ldl("A", "pre-registered (pooled)", s)
    g = lambda endpoint, model, cmp_: tuple(T.num(one(D["mp"], endpoint=endpoint, model=model, scale="lifelong", compared_with=cmp_)[k]) for k in ("hr", "lo", "hi"))
    PRI, UREV = "trial primary, four-point", "trial urgent coronary revascularization"
    L, P = T.TRAIT["LDL-C"], T.TRAIT["Lp(a)"]
    return [("LDL-C score, per mmol/L lower", [
                ("statin trials (CTT)", *CTT, T.TRIAL, "D"),
                ("before a first event, LDL-C ÷ 0.6", *ldl("C standardized", "pre-registered", "LDL-C / 0.6 in statin users"), L, "o"),
                ("established, LDL-C as measured", *A("measured LDL-C"), L, "s"),
                ("established, LDL-C ÷ 0.7", *A("LDL-C / 0.7 in statin users"), L, "s"),
                ("established, LDL-C ÷ 0.6", *A("LDL-C / 0.6 in statin users"), L, "s"),
                ("established, LDL-C ÷ 0.5", *A("LDL-C / 0.5 in statin users"), L, "s")]),
            ("Lp(a), four-point composite", [
                ("Lp(a)HORIZON, observed", *TRIAL["mace4"], T.TRIAL, "D"),
                ("predicted, one slope", *g("mace4", "linear", PRI), P, "s"),
                ("predicted, second slope", *g("mace4", "spline", PRI), P, "s"),
                ("trial-proxy endpoint, one slope", *g("m4u7", "linear", PRI), P, "s"),
                ("trial-proxy endpoint, second slope", *g("m4u7", "spline", PRI), P, "s")]),
            ("Lp(a), urgent revascularization", [
                ("Lp(a)HORIZON, observed", *TRIAL["revasc"], T.TRIAL, "D"),
                ("acute-linked, one slope", *g("urev7", "linear", UREV), P, "s"),
                ("acute-linked, second slope", *g("urev7", "spline", UREV), P, "s")])]


def groups_d(D):
    out = []
    for head, cohort, endpoint, marker in (("ACCLAIM-Lp(a), at risk", "acc_risk", "mace4", "o"),
                                           ("ACCLAIM-Lp(a), established", "acc_est", "mace4", "s"),
                                           ("OCEAN(a)-Outcomes", "ocean", "cor3", "s")):
        items = []
        for model, lab in (("linear", "one slope"), ("spline", "second slope")):
            r = one(D["fc"], cohort=cohort, endpoint=endpoint, model=model)
            items.append((lab, T.num(r["hr"]), T.num(r["lo"]), T.num(r["hi"]), T.TRAIT["Lp(a)"], marker))
        out.append((head, items))
    return out


def stage_key(P, x, y, trial=True, followup=True):
    items = [("o", T.INK, "white", "before a first event"), ("s", T.INK, "white", "established disease at baseline")]
    if followup:
        items.append(("^", T.INK, "white", "after a first event during follow-up"))
    if trial:
        items.append(("D", T.TRIAL, "white", "randomized trial"))
    items.append(("s", "white", T.INK, "open: 95% CI includes 1"))
    return P.key(x, y, items, gap=9)


def main_figure(D, out):
    P = T.Page(455.0)
    stage_key(P, 10, 9)
    P.key(10, 20, [("line", None, T.TRAIT["Lp(a)"], "Lp(a)"), ("line", None, T.TRAIT["LDL-C"], "LDL-C score (positive control)"),
                   ("line", None, T.TRAIT["non-lipid"], "non-lipid coronary score")], gap=9)
    P.panel(2, 36, "a", "The associations persist after disease onset")
    axa = P.ax(80, 50, 72, 150); draw_a(axa, D)
    P.panel(262, 36, "b", "The persisting Lp(a) association is coronary")
    axb = P.ax(348, 50, 156, 150); draw_b(axb, D)
    P.panel(2, 232, "c", "Its size does not reliably predict the trial")
    axc = P.ax(124, 246, 66, 182)
    grouped(axc, groups_c(D), "c", (0.46, 1.42), [0.5, 0.7, 1, 1.4], "HR")
    P.panel(262, 232, "d", "Conditional forecasts for the trials in progress")
    axd = P.ax(366, 246, 82, 116)
    grouped(axd, groups_d(D), "d", (0.5, 1.1), [0.5, 0.6, 0.7, 0.8, 0.9, 1], "HR for the assumed lowering")
    P.save(os.path.join(SRC_DIR, "figures", "fig_lpa_trial_comparison"))


# ------------------------------------------------------------------ Supplementary Figure: dose-response and per-144 estimates
def si_dose(ax, D):
    ax.axvspan(TRIAL_ENTRY, 320, color=T.BAND, lw=0, zorder=0)
    T.null(ax, 1.0, axis="y")
    for cohort, model, marker, key in (("C standardized", "pre-registered", "o", "C_std"), ("A", "pre-registered (pooled)", "s", "A")):
        bd = [r for r in D["bands"] if r["cohort"] == cohort and r["model"] == model and r["endpoint"] == "mace4"]
        assert len(bd) == 4
        xs = [T.num(r["predicted_mean"]) for r in bd]
        b_lin = D["pick"](key, "mace4", "x_lpa")[0]
        ax.plot([xs[0], 320], [1, math.exp(b_lin * (320 - xs[0]) / 144)], color=T.TRAIT["Lp(a)"], lw=0.8,
                ls="-" if key == "A" else (0, (3, 2)), zorder=2)
        for r, x in zip(bd, xs):
            hr, lo, hi = T.num(r["hr"]), T.num(r["hr_lo"]), T.num(r["hr_hi"])
            if r["band"] != bd[0]["band"]:
                ax.plot([x, x], [lo, hi], color=T.TRAIT["Lp(a)"], lw=T.LW["ci"], solid_capstyle="butt", zorder=3)
            ax.plot([x], [hr], marker=marker, ms=T.MS, mfc=T.TRAIT["Lp(a)"], mec="white", mew=0.5, ls="none", zorder=4)
            T.src("S", "a", "band %s (mean %.1f nmol/L)" % (r["band"], x), key, hr, lo, hi, n=r["n"], events=r["events"])
    yb = 0.94
    ax.plot([TRIAL_BASELINE[1], TRIAL_BASELINE[2]], [yb, yb], color=T.TRIAL, lw=T.LW["ci"], solid_capstyle="butt", zorder=3)
    ax.plot([TRIAL_BASELINE[0]], [yb], marker="D", ms=T.MS - 0.4, mfc=T.TRIAL, mec="white", mew=0.5, ls="none", zorder=4)
    ax.text(TRIAL_BASELINE[1] - 5, yb, "trial baseline, median and IQR", fontsize=T.FS["note"], color=T.MUTED, ha="right", va="center")
    ax.text(TRIAL_ENTRY + 4, 1.9, "trial entry range", fontsize=T.FS["note"], color=T.MUTED, ha="left", va="top")
    ax.set_yscale("log"); ax.set_ylim(0.9, 1.95)
    from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
    ax.yaxis.set_major_locator(FixedLocator([1, 1.25, 1.5, 1.75])); ax.yaxis.set_major_formatter(FixedFormatter(["1", "1.25", "1.5", "1.75"]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0, 320); ax.set_xticks([0, 50, 100, 150, 200, 250, 300])
    ax.set_xlabel("Genetically predicted Lp(a), nmol/L (band mean)")
    ax.set_ylabel("HR against the lowest band")
    T.clean(ax, left=True)


def si_per144(ax, D):
    tr = TRIAL["mace4"]
    ax.axvspan(tr[1], tr[2], color=T.BAND, lw=0, zorder=0); T.null(ax)
    b, s, _, _ = D["pick"]("A", "mace4", "x_lpa")
    sp = one(D["shape"], cohort="A", model="pre-registered (pooled)", endpoint="mace4", test="spline, knot 150 nmol/L")
    rows_ = [("Lp(a)HORIZON, observed", *tr, T.TRIAL, "D"),
             ("lifelong difference, one slope", *T.lower(b, s), T.TRAIT["Lp(a)"], "s"),
             ("slope above 150 nmol/L", T.num(sp["hr_above_per_144_lower"]), T.num(sp["above_lo"]), T.num(sp["above_hi"]), T.TRAIT["Lp(a)"], "s")]
    for key, lab in (("measured LDL-C", "scaled to 4.7 years, LDL-C as measured"), ("LDL-C / 0.7 in statin users", "scaled, LDL-C ÷ 0.7"),
                     ("LDL-C / 0.6 in statin users", "scaled, LDL-C ÷ 0.6"), ("LDL-C / 0.5 in statin users", "scaled, LDL-C ÷ 0.5")):
        r = one(D["cal2"], cohort="A", horizon="full", model="pre-registered (pooled)", ldl_scale=key)
        rows_.append((lab, T.num(r["predicted_trial_hr"]), T.num(r["pred_lo"]), T.num(r["pred_hi"]), T.TRAIT["Lp(a)"], "s"))
    n = len(rows_)
    for i, (lab, e, lo, hi, c, m) in enumerate(rows_):
        mark(ax, i, e, lo, hi, c, m)
        T.src("S", "b", lab, m, e, lo, hi)
    ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows_])
    ax.set_ylim(n - 0.4, -0.6)
    T.logx(ax, 0.55, 1.3, [0.6, 0.7, 0.8, 0.9, 1, 1.15, 1.3])
    ax.set_xlabel("HR per 144 nmol/L lower Lp(a)"); T.clean(ax)


def si_figure(D):
    P = T.Page(200.0)
    stage_key(P, 10, 9, followup=False)
    P.panel(2, 30, "a", "Dose-response across genetically predicted Lp(a)")
    axa = P.ax(42, 44, 190, 124); si_dose(axa, D)
    P.panel(262, 30, "b", "The trial beside estimates per 144 nmol/L lower Lp(a)")
    axb = P.ax(400, 44, 104, 124); si_per144(axb, D)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS_lpa_dose_prediction"))


if __name__ == "__main__":
    D = load()
    main_figure(D, None)
    si_figure(D)
    T.write_src(os.path.join(EV, "Fig2_source_data.csv"))
    print("wrote Fig. 2, its Supplementary Figure and Fig2_source_data.csv")
