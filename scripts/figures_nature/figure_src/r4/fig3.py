#!/usr/bin/env python3
"""Fig. 3: both axes track atherosclerotic disease, and the Lp(a) association is unchanged in statin users.

a  two-sample MR against MEGASTROKE subtypes (Lp(a) cis instrument; LDL-C genome-wide)
b  large-artery minus each other stroke endpoint, three Lp(a) instrument definitions
c  coronary minus stroke endpoints in UK Biobank (instrument bootstrap)
d  score x statin interaction on incident CAD under five adjustment strategies
e  the same interaction on an incident and on a prevalent-and-incident outcome
f  attenuation of each score's association by statin use

Reads the deposited tables (source/supplementary), figure_data/TableS35_statin_interaction.csv and the prevalent-
and-incident rows printed in _src_circgpm_supplement.tex; fits nothing.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theme as T  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
S35 = os.path.join(SRC_DIR, "figure_data", "TableS35_statin_interaction.csv")
SUPTEX = os.path.join(SRC_DIR, "_src_circgpm_supplement.tex")
LDL, LPA = T.TRAIT["LDL-C"], T.TRAIT["Lp(a)"]
R = lambda f: T.rows(os.path.join(SUP, f))

SUBTYPE = [("MEGASTROKE_LAS", "Large-artery atherosclerotic", True), ("MEGASTROKE_CES", "Cardioembolic", False),
           ("MEGASTROKE_SVS", "Small-vessel", False), ("MEGASTROKE_AIS", "Any ischaemic stroke", False),
           ("MEGASTROKE_AS", "Any stroke", False)]


def get(d, exp, out, method):
    r = [x for x in d if x["exposure"] == exp and x["outcome"] == out and x["method"] == method]
    assert len(r) == 1, (exp, out, method, len(r))
    return r[0]


def draw_a(ax):
    d = R("TableS36_stroke_subtype_mr.csv")
    T.null(ax)
    for i, (out, lab, athero) in enumerate(SUBTYPE):
        if athero:
            T.band(ax, i - 0.48, i + 0.48)
        for exp, dy, c in (("LPA", -0.17, LPA), ("LDL_C", 0.17, LDL)):
            r = get(d, exp, out, "IVW")
            o, lo, hi, p = float(r["OR_perSD"]), float(r["CI_lo95"]), float(r["CI_hi95"]), float(r["pval"])
            T.point(ax, i + dy, o, lo, hi, c, "s", filled=p < 0.01)
            T.src("3", "a", lab, exp + " IVW", o, lo, hi, p="%.3g" % p)
            for m, mk in (("weighted median", "o"), ("MR-Egger", "^")):
                rr = get(d, exp, out, m)
                ax.plot([float(rr["OR_perSD"])], [i + dy], marker=mk, ms=2.8, mfc="white", mec=c, mew=0.7, ls="none", zorder=5)
                T.src("3", "a", lab, "%s %s" % (exp, m), float(rr["OR_perSD"]), float(rr["CI_lo95"]), float(rr["CI_hi95"]))
    ax.set_yticks(range(len(SUBTYPE))); ax.set_yticklabels([s[1] for s in SUBTYPE])
    ax.set_ylim(len(SUBTYPE) - 0.5, -0.5)
    T.logx(ax, 0.86, 1.72, [0.9, 1, 1.2, 1.4, 1.6])
    ax.set_xlabel("Odds ratio per SD, MEGASTROKE"); T.clean(ax)


def draw_b(ax):
    con, het = R("TableS39_subtype_contrasts.csv"), R("TableS39b_subtype_heterogeneity.csv")
    m = re.search(r"cis \$r\^2<0\.01\$ \((\d+) instruments\)", open(SUPTEX).read())
    assert m, "instrument count not found"
    nv = int(m.group(1))
    order = [("LAS vs CES", "vs cardioembolic"), ("LAS vs SVS", "vs small-vessel"), ("LAS vs AIS", "vs any ischaemic"), ("LAS vs AS", "vs any stroke")]
    inst = [("cis r2<0.01", "s"), ("cis r2<0.001", "^"), ("rs10455872 only", "o")]
    T.null(ax, 0.0)
    for i, (c, lab) in enumerate(order):
        for k, (ins, mk) in enumerate(inst):
            r = next(x for x in con if x["contrast"] == c and x["instrument"] == ins)
            dd = float(r["delta_logOR"])
            seb = float(r["se_bootstrap"]) if r["se_bootstrap"] not in ("", "0") else 0.0
            se = seb if seb > 0 else float(r["se_independent"])
            T.point(ax, i + (k - 1) * 0.25, dd, dd - 1.96 * se, dd + 1.96 * se, LPA, mk, filled=T.excludes(dd - 1.96 * se, dd + 1.96 * se, 0.0))
            T.src("3", "b", lab, ins, dd, dd - 1.96 * se, dd + 1.96 * se)
    q = next(h for h in het if h["instrument"] == "cis r2<0.01")
    ax.text(1.0, 1.02, "Cochran $Q$ = %.1f, df %d, $P$ = %s" % (float(q["Q"]), int(float(q["df"])), T.p_text(float(q["P"]), 2)),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([o[1] for o in order])
    ax.set_ylim(len(order) - 0.55, -0.55)
    ax.set_xlim(-0.03, 0.42); ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlabel("Difference in log OR, large-artery minus other"); T.clean(ax)
    return nv


def draw_c(ax):
    d = R("TableS37_endpoint_contrasts.csv")
    lab = {"CAD vs STROKE": "CAD vs stroke", "MI vs STROKE": "MI vs stroke", "REVASC vs STROKE": "Revascularization vs stroke"}
    T.null(ax, 0.0)
    for i, r in enumerate(d):
        df, se, p = float(r["diff_logOR"]), float(r["se_bootstrap"]), float(r["p_bootstrap"])
        T.point(ax, i, df, df - 1.96 * se, df + 1.96 * se, LPA, "s")
        ax.text(df + 1.96 * se + 0.008, i, "$P$ = " + T.p_text(p), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
        T.src("3", "c", lab[r["contrast"]], "Lp(a)", df, df - 1.96 * se, df + 1.96 * se, p="%.3g" % p)
    ax.set_yticks(range(len(d))); ax.set_yticklabels([lab[r["contrast"]] for r in d])
    ax.set_ylim(len(d) - 0.5, -0.5)
    ax.set_xlim(-0.02, 0.3); ax.set_xticks([0, 0.1, 0.2, 0.3])
    ax.set_xlabel("Difference in log OR for Lp(a)"); T.clean(ax)


MODELS = [("primary (as reported)", "Primary model, this subset"), ("+ indication covariates", "+ indication covariates"),
          ("+ covariate interactions", "+ covariate interactions"), ("stabilized IPTW", "Stabilized IPTW"),
          ("both scores, mutually adjusted", "Both scores, mutually adjusted*")]


def draw_d(ax):
    rob = R("TableS34b_statin_robustness.csv")
    T.null(ax)
    ax.axvspan(0.855, 1.0, color=T.BAND, lw=0, zorder=0)
    for i, (m, lab) in enumerate(MODELS):
        for k, (term, c) in enumerate((("score_LDL:statin", LDL), ("score_LPA:statin", LPA))):
            r = next(x for x in rob if x["model"] == m and x["term"] == term)
            e, lo, hi = float(r["estimate"]), float(r["lo95"]), float(r["hi95"])
            T.point(ax, i + (k - 0.5) * 0.32, e, lo, hi, c, "s", filled=T.excludes(lo, hi))
            T.src("3", "d", lab, term, e, lo, hi)
    ax.set_yticks(range(len(MODELS))); ax.set_yticklabels([m[1] for m in MODELS])
    ax.set_ylim(len(MODELS) - 0.4, -0.6)
    ax.set_xlim(0.855, 1.065); ax.set_xticks([0.88, 0.92, 0.96, 1.0, 1.04])
    ax.text(0.02, 1.01, "shaded: weaker in statin users", transform=ax.transAxes, ha="left", va="bottom", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_xlabel("Ratio of HRs, score × statin"); T.clean(ax)


def strata():
    inc = {r["term"]: r for r in T.rows(S35)}
    tex = open(SUPTEX).read()

    def prevalent(score):
        m = re.search(r"CAD \(prevalent\+incident\) & " + re.escape(score) + r" score \$\\times\$ statin\s*& ratio of ORs ([\d.]+) \(([\d.]+)--([\d.]+)\)", tex)
        assert m, score
        return tuple(float(x) for x in m.groups())

    def incident(term):
        r = inc[term]
        return float(r["HR"]), float(r["lo95"]), float(r["hi95"])
    return [("Lp(a), prevalent and incident", prevalent("Lp(a)"), LPA, "^"),
            ("Lp(a), incident only", incident("Lp(a) score x baseline statin (ratio of HRs)"), LPA, "o"),
            ("LDL-C, prevalent and incident", prevalent("LDL-C"), LDL, "^"),
            ("LDL-C, incident only", incident("LDL-C score x baseline statin (ratio of HRs)"), LDL, "o")]


def draw_e(ax):
    st = strata()
    T.null(ax)
    for i, (lab, (e, lo, hi), c, mk) in enumerate(st):
        T.point(ax, i, e, lo, hi, c, mk, filled=T.excludes(lo, hi))
        T.src("3", "e", lab, mk, e, lo, hi)
    ax.annotate("reverse causation\nor selection", xy=(st[0][1][2], 0), xytext=(0.925, 0.5), fontsize=T.FS["note"], color=T.MUTED,
                ha="center", va="center", linespacing=1.0, arrowprops=dict(arrowstyle="-", color=T.MUTED, lw=0.4, shrinkA=0.5, shrinkB=1.5))
    ax.set_yticks(range(len(st))); ax.set_yticklabels([s[0] for s in st])
    ax.set_ylim(len(st) - 0.5, -0.5)
    ax.set_xlim(0.885, 1.075); ax.set_xticks([0.9, 0.95, 1.0, 1.05])
    ax.set_xlabel("Interaction with statin use (ratio)"); T.clean(ax)


def draw_f(ax):
    rob = R("TableS34b_statin_robustness.csv")
    g = lambda t: next(r for r in rob if r["model"] == "primary (as reported)" and r["term"] == t)
    att = lambda v: (1.0 - float(v)) * 100.0
    br = []
    for lab, term, main, c in (("LDL-C score", "score_LDL:statin", "score_LDL", LDL), ("Lp(a) score", "score_LPA:statin", "score_LPA", LPA)):
        r = g(term)
        br.append((lab, att(r["estimate"]), att(r["hi95"]), att(r["lo95"]), c))
    bound = br[1][3]
    ax.axvspan(bound, 14.0, color=T.BAND, lw=0, zorder=0)
    T.null(ax, 0.0)
    for i, (lab, e, lo, hi, c) in enumerate(br):
        T.point(ax, i, e, lo, hi, c, "s", filled=T.excludes(lo, hi, 0.0))
        ax.text(e, i - 0.3, "%.1f%%" % e, ha="center", va="bottom", fontsize=T.FS["note"])
        T.src("3", "f", lab, "attenuation, %", e, lo, hi)
    ax.text(13.6, 1.55, "excluded by the Lp(a) interval (> %.1f%%)" % bound, ha="right", va="center", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_yticks([0, 1]); ax.set_yticklabels([b[0] for b in br])
    ax.set_ylim(1.85, -0.6)
    ax.set_xlim(-3.0, 14.0); ax.set_xticks([0, 4, 8, 12])
    ax.set_xlabel("Attenuation by statin use (%)"); T.clean(ax)
    return bound


def build():
    P = T.Page(472.0)
    P.panel(2, 10, "a", "Stroke subtype MR")
    P.key(12, 22, [("line", None, LPA, r"Lp(a), $cis$ instrument"), ("line", None, LDL, "LDL-C, genome-wide")], gap=9)
    P.key(12, 31, [("s", T.INK, "white", "IVW"), ("o", "white", T.INK, "weighted median"), ("^", "white", T.INK, "MR-Egger"),
                   ("s", "white", T.INK, r"open IVW: $P\geq0.01$")], gap=8)
    axa = P.ax(92, 40, 158, 180); draw_a(axa)
    P.panel(262, 10, "b", "Between-subtype contrast for Lp(a)")
    nv = int(re.search(r"cis \$r\^2<0\.01\$ \((\d+) instruments\)", open(SUPTEX).read()).group(1))
    P.key(272, 22, [("s", LPA, "white", r"$cis$ $r^2<0.01$ (%d variants)" % nv), ("^", LPA, "white", r"$cis$ $r^2<0.001$"),
                    ("o", LPA, "white", "rs10455872 alone")], gap=8)
    axb = P.ax(340, 42, 164, 76); draw_b(axb)
    P.panel(262, 146, "c", "Coronary minus stroke endpoints")
    axc = P.ax(372, 156, 124, 50); draw_c(axc)
    P.panel(2, 252, "d", "Score × statin interaction by adjustment")
    P.key(12, 264, [("s", LDL, "white", "LDL-C score"), ("s", LPA, "white", "Lp(a) score"), ("s", "white", T.INK, "open: CI includes 1")], gap=8)
    axd = P.ax(112, 274, 138, 170); draw_d(axd)
    P.panel(262, 252, "e", "Incident versus prevalent outcome")
    P.key(272, 264, [("o", T.INK, "white", "incident, ratio of HRs"), ("^", T.INK, "white", "prevalent and incident, ratio of ORs")], gap=8)
    axe = P.ax(386, 272, 118, 62); draw_e(axe)
    P.panel(262, 368, "f", "Attenuation by statin use")
    axf = P.ax(330, 380, 174, 64); draw_f(axf)
    P.save(os.path.join(SRC_DIR, "figures", "fig_axes_specificity"))


if __name__ == "__main__":
    build()
    print("wrote Fig. 3")
