#!/usr/bin/env python3
"""Fig. 1: study design, and the two causal lipid axes.

a  the setting of the genetic evidence beside the setting of Lp(a)HORIZON, and the causal-analysis samples
b  genetic correlation between lipid pairs (LD-score regression)
c  LDL-C and ApoB modelled jointly under three perturbations: the coefficients invert
d  multivariable MR of LDL-C, triglycerides, HDL-C and Lp(a) against eight endpoints

Reads the deposited tables (source/supplementary) and the stage estimates (evidence/); fits nothing.
"""
import os
import sys

from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theme as T  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
EV = os.path.abspath(os.path.join(SRC_DIR, "..", "evidence", "lpa_emulation"))
N = lambda x: format(int(x), ",")


def stage_counts():
    est = T.rows(os.path.join(EV, "emulation", "estimates.csv"))
    out = {}
    for c in ("C_std", "A", "D"):
        r = [x for x in est if (x["cohort"], x["endpoint"], x["exposure"], x["horizon"], x["model"]) == (c, "mace4", "x_lpa", "full", "joint")]
        assert len(r) == 1, c
        out[c] = (int(r[0]["n"]), int(r[0]["events"]))
    return out


def cohort_table():
    t = {(r["category"], r["trait"]): r["N_or_cases"] for r in T.rows(os.path.join(SUP, "TableS1_cohort_phenotypes.csv"))}
    return int(t[("cohort", "genetically_inferred_European_unrelated")]), int(t[("binary_outcome", "CAD")])


# ------------------------------------------------------------------ a, design
def draw_a(P, x0, y0):
    ax = P.ax(x0, y0, 348, 112)
    ax.set_xlim(0, 348); ax.set_ylim(112, 0); ax.axis("off")
    box = lambda x, y, w, h, fc: ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor="none", zorder=1))
    txt = lambda x, y, s, **kw: ax.text(x, y, s, fontsize=kw.pop("fs", T.FS["note"]), color=kw.pop("c", T.INK), va=kw.pop("va", "center"), **kw)
    c1, c2, onset = (58, 196), (204, 348), 200
    txt((c1[0] + c1[1]) / 2, 4, "Before disease onset", ha="center", c=T.MUTED)
    txt((c2[0] + c2[1]) / 2, 4, "After disease onset", ha="center", c=T.MUTED)
    ax.plot([onset, onset], [9, 108], color=T.NULL, lw=T.LW["ref"], zorder=0)
    txt(onset, 111, "first myocardial infarction, ischaemic stroke or PAD", ha="center", c=T.MUTED, fs=T.FS["small"])
    for y, lab in ((19, "Genotype"), (48, "UK Biobank\ncohorts"), (91, "Trial")):
        txt(0, y, lab, ha="left", c=T.MUTED, linespacing=1.1)
    # lifelong genetic difference: one bar across both stages
    box(c1[0], 12, c2[1] - c1[0], 14, T.tint(T.TRAIT["Lp(a)"], 0.18))
    txt(c1[0] + 4, 19, "A lifelong genetic difference in Lp(a), LDL-C and non-lipid coronary risk", ha="left")
    ax.annotate("", xy=(c2[1] + 0.5, 19), xytext=(c2[1] - 6, 19), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=T.MUTED,
                mutation_scale=5), zorder=2)
    n = stage_counts()
    box(c1[0], 33, c1[1] - c1[0], 30, T.TINT)
    txt(c1[0] + 4, 42, "Free of vascular disease at baseline", ha="left")
    txt(c1[0] + 4, 53, "%s participants; %s events" % (N(n["C_std"][0]), N(n["C_std"][1])), ha="left", c=T.MUTED)
    box(c2[0], 33, c2[1] - c2[0], 19, T.TINT)
    txt(c2[0] + 4, 38.5, "Prior MI, ischaemic stroke or PAD at baseline", ha="left")
    txt(c2[0] + 4, 47.0, "%s participants; %s events" % (N(n["A"][0]), N(n["A"][1])), ha="left", c=T.MUTED)
    box(c2[0], 55, c2[1] - c2[0], 19, T.TINT)
    txt(c2[0] + 4, 60.5, "From 30 days after a first MI or ischaemic stroke", ha="left")
    txt(c2[0] + 4, 69.0, "%s participants; %s events" % (N(n["D"][0]), N(n["D"][1])), ha="left", c=T.MUTED)
    # the trial: a few years of treatment after onset
    box(c2[0], 83, 44, 16, T.tint(T.INK, 0.22))
    txt(c2[0] + 22, 91, "Lp(a)HORIZON", ha="center")
    txt(c2[0] + 48, 87.5, "8,322 patients; median 4.7 years", ha="left", c=T.MUTED)
    txt(c2[0] + 48, 95.5, "Lp(a) 144 nmol/L lower; HR 1.02", ha="left", c=T.MUTED)
    for k, (c, (nn, ev)) in enumerate(n.items()):
        T.src("1", "a", c, "participants; events", nn, ev, ev)


def draw_a_right(P, x0, y0):
    ax = P.ax(x0, y0, 146, 112)
    ax.set_xlim(0, 146); ax.set_ylim(112, 0); ax.axis("off")
    ax.add_patch(Rectangle((0, 12), 146, 87, facecolor=T.TINT, edgecolor="none"))
    n_all, n_cad = cohort_table()
    lines = [("Causal analyses", True, T.INK),
             ("UK Biobank: %s participants," % N(n_all), False, T.INK),
             ("%s CAD cases" % N(n_cad), False, T.INK),
             ("Multivariable MR, LD-score regression,", False, T.MUTED),
             ("colocalization, rare-variant burden", False, T.MUTED),
             ("External outcomes", True, T.INK),
             ("FinnGen and CARDIoGRAMplusC4D (CAD)", False, T.MUTED),
             ("MEGASTROKE (stroke subtypes)", False, T.MUTED),
             ("Million Veteran Program (African and", False, T.MUTED),
             ("Hispanic∕admixed-American ancestry)", False, T.MUTED)]
    y = 19
    for s, bold, c in lines:
        if bold and y > 20:
            y += 2.5
        ax.text(5, y, s, fontsize=T.FS["note"], fontweight="bold" if bold else "normal", color=c, va="center")
        y += 8.0


# ------------------------------------------------------------------ b, genetic correlation
def draw_b(ax):
    rg = {}
    for r in T.rows(os.path.join(SUP, "TableS2b_ldsc_rg_lipid_pairs.csv")):
        if r["rg"] and r["se"]:
            rg[frozenset((r["exposure"].replace("lipids_", ""), r["outcome"].replace("lipids_", "")))] = (float(r["rg"]), float(r["se"]))
    pairs = [(("LDL_C", "APOB"), "LDL-C and ApoB"), (("LDL_C", "TOT_CHOL"), "LDL-C and total cholesterol"),
             (("APOB", "TOT_CHOL"), "ApoB and total cholesterol"), (("HDL_C", "APOA1"), "HDL-C and ApoA1")]
    T.null(ax, 1.0)
    for i, (p, lab) in enumerate(pairs):
        y = len(pairs) - 1 - i
        v, se = rg[frozenset(p)]
        T.point(ax, y, v, v - 1.96 * se, v + 1.96 * se, T.INK, "s")
        ax.text(1.02, y, T.ci_text(v, v - 1.96 * se, v + 1.96 * se), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"])
        T.src("1", "b", lab, "r_g", v, v - 1.96 * se, v + 1.96 * se)
    ax.text(1.02, len(pairs) - 0.35, r"$r_g$ (95% CI)", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_yticks(range(len(pairs))); ax.set_yticklabels([l for _, l in pairs][::-1])
    ax.set_ylim(-0.6, len(pairs) - 0.4); ax.set_xlim(0.6, 1.45); ax.set_xticks([0.6, 0.8, 1.0, 1.2, 1.4])
    ax.set_xlabel(r"Genetic correlation, $r_g$"); T.clean(ax)


# ------------------------------------------------------------------ c, the LDL-C/ApoB inversion
def draw_c(ax):
    R = T.rows(os.path.join(SUP, "TableS3c_identifiability.csv"))
    pick = lambda arm, ex, model: next((r for r in R if r["arm"] == arm and r["exposure"] == ex and r["model"] == model), None)
    arms = [("biochemistry panel", "Biochemistry panel", (("LDL_C", "LDL-C", "LDL-C"), ("APOB", "ApoB", "ApoB"))),
            ("statin-naive GWAS", "Statin-naive participants", (("LDL_C", "LDL-C", "LDL-C"), ("APOB", "ApoB", "ApoB"))),
            ("NMR particle count", "NMR", (("LDL_C_NMR", "LDL-C", "LDL-C"), ("LDL_P", "LDL-P", "ApoB")))]
    T.null(ax, 1.0)
    y, ticks, labels, heads = 0.0, [], [], []
    for k, (arm, title, exps) in enumerate(arms):
        j0 = pick(arm, exps[0][0], "joint (MVMR)"); j1 = pick(arm, exps[1][0], "joint (MVMR)")
        heads.append((y - 0.62, "%s: %s instruments, conditional $F$ %.0f and %.0f" % (title, j0["n_instruments"], float(j0["condF_phenocov"]), float(j1["condF_phenocov"]))))
        if k % 2 == 0:
            T.band(ax, y - 1.0, y + 1.45)
        for i, (ex, lab, hue) in enumerate(exps):
            yy = y + i * 0.85
            u = pick(arm, ex, "univariable")
            if u is not None:
                T.point(ax, yy - 0.28, float(u["OR"]), float(u["lo"]), float(u["hi"]), T.NS, "o", ms=T.MS - 0.6, lw=0.8)
                T.src("1", "c", "%s, %s" % (title, lab), "univariable", float(u["OR"]), float(u["lo"]), float(u["hi"]))
            j = pick(arm, ex, "joint (MVMR)")
            o, lo, hi = float(j["OR"]), float(j["lo"]), float(j["hi"])
            T.point(ax, yy, o, lo, hi, T.TRAIT[hue], "s", filled=T.excludes(lo, hi))
            ax.text(1.02, yy, T.ci_text(o, lo, hi), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"])
            ticks.append(yy); labels.append(lab)
            T.src("1", "c", "%s, %s" % (title, lab), "joint", o, lo, hi)
        y += 2.45
    for yh, s in heads:
        ax.text(0.0, yh, s, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.text(1.02, -1.0, "Joint OR (95% CI)", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_yticks(ticks); ax.set_yticklabels(labels)
    T.logx(ax, 0.2, 5.0, [0.25, 0.5, 1, 2, 4])
    ax.set_ylim(y - 0.9, -1.25)
    ax.set_xlabel("Odds ratio for CAD per SD"); T.clean(ax)


# ------------------------------------------------------------------ d, multivariable MR
OUTCOMES = [("CAD", "CAD"), ("MI", "Myocardial infarction"), ("REVASC", "Revascularization"), ("MACE", "MACE"),
            ("PAD", "Peripheral artery disease"), ("CV_DEATH", "Cardiovascular death"), ("HF", "Heart failure"), ("STROKE", "Ischaemic stroke")]
EXPS = [("LDL_C", "LDL-C", "LDL-C"), ("TRIGLY", "Triglycerides", "TG"), ("HDL_C", "HDL-C", "HDL-C"), ("LPA", "Lp(a)", "Lp(a)")]


def draw_d(ax):
    R = T.rows(os.path.join(SUP, "TableS3_mvmr.csv"))
    assert len(R) == 32 and len({r["n_instruments"] for r in R}) == 1
    tab = {(r["exposure"], r["outcome"]): tuple(float(r[c]) for c in ("mvivw_OR", "mvivw_OR_lo95", "mvivw_OR_hi95", "mvivw_pval")) for r in R}
    bonf = 0.05 / len(R)
    T.null(ax, 1.0)
    y, ticks, state = 0.0, [], [0, 0, 0]
    step, gap = 0.78, 0.55
    for k, (oc, lab) in enumerate(OUTCOMES):
        top = y
        for ex, name, hue in EXPS:
            o, lo, hi, p = tab[(ex, oc)]
            s = 0 if p < bonf else (1 if p < 0.05 else 2)
            state[s] += 1
            c = T.TRAIT[hue]
            if s == 1:
                T.point(ax, y, o, lo, hi, T.tint(c, 0.45), "s", ci_colour=c)
            else:
                T.point(ax, y, o, lo, hi, c, "s", filled=(s == 0))
            T.src("1", "d", lab, name, o, lo, hi, p="%.3g" % p)
            y += step
        if k % 2 == 0:
            T.band(ax, top - 0.5, y - step + 0.5)
        ticks.append((top + y - step) / 2)
        y += gap
    assert state == [14, 6, 12], state
    ax.set_yticks(ticks); ax.set_yticklabels([l for _, l in OUTCOMES])
    ax.set_ylim(y - gap - step + 0.6, -0.6)
    T.logx(ax, 0.72, 1.6, [0.8, 0.9, 1, 1.2, 1.4, 1.6])
    ax.set_xlabel("Odds ratio per SD (161 instruments)"); T.clean(ax)
    return R[0]["n_instruments"]


def build(out):
    P = T.Page(470.0)
    P.panel(2, 10, "a", "Study design")
    draw_a(P, 6, 16)
    draw_a_right(P, 362, 16)
    P.panel(2, 146, "b", "LDL-C and ApoB are genetically indistinguishable")
    axb = P.ax(104, 156, 70, 62); draw_b(axb)
    P.panel(2, 246, "c", "Modelled jointly, they cannot be ranked")
    P.key(12, 258, [("s", T.INK, "white", "joint model"), ("s", "white", T.INK, "joint, CI includes 1"), ("o", T.NS, "white", "each fitted alone")], gap=8)
    axc = P.ax(58, 268, 116, 170); draw_c(axc)
    P.panel(262, 146, "d", "LDL/ApoB and Lp(a) carry independent effects")
    P.key(272, 160, [("s", T.TRAIT["LDL-C"], "white", "LDL-C"), ("s", T.TRAIT["TG"], "white", "Triglycerides"),
                     ("s", T.TRAIT["HDL-C"], "white", "HDL-C"), ("s", T.TRAIT["Lp(a)"], "white", "Lp(a)")], gap=8)
    P.key(272, 170, [("s", T.INK, "white", r"$P<1.6\times10^{-3}$"), ("s", T.tint(T.INK, 0.45), "white", r"$P<0.05$"),
                     ("s", "white", T.INK, r"$P\geq0.05$")], gap=8)
    axd = P.ax(352, 178, 152, 260); draw_d(axd)
    P.save(os.path.join(out, "fig1_architecture"))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SRC_DIR, "figures")
    build(out)
    print("wrote", os.path.join(out, "fig1_architecture.pdf"))
