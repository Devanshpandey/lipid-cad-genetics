#!/usr/bin/env python3
"""
Figure 4 — Integrated common+rare evidence improves recovery of established lipid genes.
  (a) top-gene evidence architecture matrix (why genes rank; category + rank)
  (b) discrimination: ROC burden-only vs full score + layer-ablation AUROC inset
  (c) recall at practical screening depths (top 500), integrated vs burden
  (d) robustness: PDE3B rank across 108 specs + recovery AUROC by positive class

Outputs fig3_convergence.{pdf,png} (the manuscript's Figure 4 include).
Data: TableS8 (score components), roc_curve_points/recovery_curve (strengthen),
TableS20 (ablation), TableS22 (weight sensitivity), TableS23 (stratified AUROC).
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
SUP = os.path.join(HERE, "..", "supplementary")
TEAL, BURG, AMBER = LIPID["rare"], LIPID["outcome"], "#E08A2B"
CATCOL = {"established": TEAL, "candidate": AMBER, "tagging": GREY}


def hd(ax, letter, title, sub=None, lx=-0.02):
    ax.text(lx, 1.15, letter, transform=ax.transAxes, fontsize=13.5,
            fontweight="bold", va="top", ha="right", color=INK)
    ax.text(0.5, 1.145, title, transform=ax.transAxes, fontsize=9.6,
            fontweight="bold", va="top", ha="center", color=INK)
    if sub:
        ax.text(0.5, 1.065, sub, transform=ax.transAxes, fontsize=7.3,
                va="top", ha="center", color=SUBTLE)


fig = plt.figure(figsize=(13.2, 10.6))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.28, width_ratios=[1.12, 1.0])

# ======================================================= (a) evidence matrix
axA = fig.add_subplot(gs[0, 0])
rows = list(csv.DictReader(open(os.path.join(SUP, "TableS8_convergence_scores.csv"))))[:16]
CATEG = {  # classification per the manuscript narrative
 "LDLR":"established","PCSK9":"established","ANGPTL3":"established","APOB":"established",
 "LPA":"established","ABCA1":"established","APOA5":"established","APOC3":"established",
 "CETP":"established","LCAT":"established","APOA1":"established","LPL":"established",
 "PLG":"tagging","SLC22A2":"tagging","ZNF229":"tagging","ABO":"candidate","PDE3B":"candidate"}
DRUGMEND = {"LDLR","PCSK9","ANGPTL3","APOB","LPA","ABCA1","APOA5","APOC3","CETP","LCAT","APOA1","LPL"}
FINEMAP = {"PCSK9","LDLR"}
CAP = 60.0
cols = ["lipid\nburden", "CV\nburden", "coloc", "fine-\nmap", "drug/\nMend."]
ng = len(rows)
for i, r in enumerate(rows):
    g = r["gene"]; y = ng - 1 - i
    lb = min(float(r["lipid_burden_log10P"]), CAP); ob = float(r["outcome_burden_log10P"])
    axA.scatter(0, y, s=12 + 120 * lb / CAP, c=[TEAL], edgecolor="black", lw=0.4, zorder=3)
    axA.scatter(1, y, s=12 + 120 * min(ob, 13) / 13, c=[BURG], edgecolor="black", lw=0.4, zorder=3)
    if r["coloc_locus"] == "1":
        axA.scatter(2, y, s=95, marker="D", c=[TEAL], edgecolor="white", lw=0.8, zorder=3)
    if g in FINEMAP:
        axA.scatter(3, y, s=95, marker="D", c=[TEAL], edgecolor="white", lw=0.8, zorder=3)
    if g in DRUGMEND:
        axA.scatter(4, y, s=70, marker="*", c=[INK], edgecolor="none", zorder=3)
    cat = CATEG.get(g, "candidate")
    axA.text(5.4, y, cat, ha="left", va="center", fontsize=6.6, color=CATCOL[cat], fontweight="bold")
    axA.text(7.9, y, str(i + 1), ha="center", va="center", fontsize=7.0, color=INK)
    axA.text(-0.7, y, g, ha="right", va="center", fontsize=7.4, fontstyle="italic", color=INK)
for j, c in enumerate(cols):
    axA.text(j, ng - 0.3, c, ha="center", va="bottom", fontsize=6.6, fontweight="bold",
             color=INK, linespacing=0.9)
axA.text(5.4, ng - 0.3, "category", ha="left", va="bottom", fontsize=6.6, fontweight="bold", color=INK)
axA.text(7.9, ng - 0.3, "rank", ha="center", va="bottom", fontsize=6.6, fontweight="bold", color=INK)
axA.set_xlim(-3.0, 8.6); axA.set_ylim(-0.6, ng + 0.2)
axA.set_xticks([]); axA.set_yticks([])
for sp in axA.spines.values():
    sp.set_visible(False)
hd(axA, "a", "Top-ranked gene evidence architecture",
   "dot size = burden $-\\log_{10}P$; diamond = coloc/fine-map; star = established gene")

# ======================================================= (b) precision-recall
axB = fig.add_subplot(gs[0, 1])
# build PR curves + metrics from the full gene ranking (imbalance-appropriate)
grows = list(csv.DictReader(open(os.path.join(DATA, "gene_ranking_scores.csv"))))
for r in grows:
    r["canonical"] = int(float(r["canonical"]))
NPOS = sum(r["canonical"] for r in grows)


def pr_curve(scorekey):
    srt = sorted(grows, key=lambda r: -float(r[scorekey]))
    rec_x, prec_y, tp = [], [], 0
    for i, r in enumerate(srt, 1):
        if r["canonical"] == 1:
            tp += 1
            rec_x.append(tp / NPOS); prec_y.append(tp / i)
    ap = sum(prec_y) / NPOS
    p_at = {k: sum(x["canonical"] for x in srt[:k]) / k for k in (10, 25, 50)}
    return rec_x, prec_y, ap, p_at


for name, col, lab in (("burden_only", GREY, "burden only"), ("integrated", TEAL, "full score")):
    rx, py, ap, pat = pr_curve(name)
    axB.step(rx, py, where="post", color=col, lw=2.2, zorder=3,
             label=f"{lab} (PR-AUC {ap:.2f})")
    if name == "integrated":
        pat_i = pat
    else:
        pat_b = pat
axB.set_xlabel("recall of established genes"); axB.set_ylabel("precision")
axB.set_xlim(-0.02, 1.02); axB.set_ylim(0, 1.03)
axB.legend(loc="upper right", fontsize=7.0); despine(axB)
axB.text(0.03, 0.30,
         "precision@k (full / burden):\n"
         f"@10  {pat_i[10]:.2f} / {pat_b[10]:.2f}\n"
         f"@25  {pat_i[25]:.2f} / {pat_b[25]:.2f}\n"
         f"@50  {pat_i[50]:.2f} / {pat_b[50]:.2f}",
         transform=axB.transAxes, fontsize=6.4, va="top", color=INK)
axB.text(0.03, 0.03, "AUROC 0.97 vs 0.90 (secondary; inflated by class imbalance)",
         transform=axB.transAxes, fontsize=5.6, color=SUBTLE)
# layer-ablation inset (secondary AUROC diagnostic: which layer earns the gain)
# placed in the empty region below the high-precision part of the curves
# (axes-fraction from data box x[0.05,0.35] y[0.40,0.75] under xlim -0.02..1.02, ylim 0..1.03)
ins = axB.inset_axes([0.067, 0.388, 0.289, 0.340])
ab = list(csv.DictReader(open(os.path.join(SUP, "TableS20_layer_ablation.csv"))))
labs = ["burden", "b+coloc", "b+outc.", "full"]
au = [float(r["AUROC_canonical_recovery"]) for r in ab]
ins.bar(range(4), au, color=[GREY, TEAL, GREY, TEAL], width=0.7)
for xi, v in enumerate(au):
    ins.text(xi, v + 0.004, f"{v:.2f}", ha="center", fontsize=5.2, color=INK)
ins.set_xticks(range(4)); ins.set_xticklabels(labs, fontsize=5.2, rotation=25, ha="right")
ins.set_ylim(0.85, 1.0); ins.set_yticks([0.9, 1.0]); ins.tick_params(labelsize=5.0)
ins.set_title("ablation (AUROC): coloc earns the gain", fontsize=5.6, pad=2); despine(ins)
hd(axB, "b", "Precision-recall (biological calibration)",
   "27 established genes vs 18,002; not prospective validation")

# ======================================================= (c) recall at depth
axC = fig.add_subplot(gs[1, 0])
rec = list(csv.DictReader(open(os.path.join(DATA, "recovery_curve.csv"))))
for name, col, lab in (("burden", GREY, "burden only"), ("integrated", TEAL, "full score")):
    pts = [(int(r["k"]), float(r["recall"])) for r in rec if r["curve"] == name and int(r["k"]) <= 500]
    pts.sort()
    axC.plot([p[0] for p in pts], [p[1] for p in pts], color=col, lw=2.2, label=lab, zorder=3)
# mark practical depths
def recall_at(name, k):
    for r in rec:
        if r["curve"] == name and int(r["k"]) == k:
            return float(r["recall"])
    return None
for k in (10, 25, 50, 100, 250):
    v = recall_at("integrated", k)
    if v is not None:
        axC.scatter(k, v, s=20, color=TEAL, zorder=4)
# label only well-separated depths, above-right of each marker, to avoid crowding
for k in (50, 100, 250):
    v = recall_at("integrated", k)
    if v is not None:
        axC.annotate(f"{round(v*27)}/27 @{k}", (k, v), textcoords="offset points",
                     xytext=(7, 5), ha="left", fontsize=6.0, color=INK)
axC.set_xlabel("genes screened (top-$k$)"); axC.set_ylabel("recall of established genes")
axC.set_xlim(0, 500); axC.set_ylim(0, 1.02)
axC.legend(loc="lower right", fontsize=7.2); despine(axC)
axC.text(0.02, 0.97, "x-axis truncated at top 500\n(full curve reaches 1.0 genome-wide)",
         transform=axC.transAxes, fontsize=6.2, va="top", color=SUBTLE)
hd(axC, "c", "Recovery at practical screening depths",
   "integrated score reaches higher recall at every depth")

# ======================================================= (d) robustness
axD = fig.add_subplot(gs[1, 1])
gsd = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 1], wspace=0.55, width_ratios=[1.0, 1.0])
axD.axis("off")
axD1 = fig.add_subplot(gsd[0]); axD2 = fig.add_subplot(gsd[1])
ws = list(csv.DictReader(open(os.path.join(SUP, "TableS22_score_weight_sensitivity.csv"))))
ranks = np.array([int(r["pde3b_rank"]) for r in ws])
rng = np.random.default_rng(42)
axD1.scatter(rng.uniform(-0.20, 0.20, len(ranks)), ranks, s=15, color=CB["purple"],
             alpha=0.55, edgecolor="none", zorder=3)
axD1.boxplot([ranks], positions=[0], widths=0.52, showfliers=False,
             medianprops=dict(color=INK, lw=1.4), boxprops=dict(color=INK),
             whiskerprops=dict(color=INK), capprops=dict(color=INK))
axD1.axhline(18, color=GREY, ls=(0, (3, 3)), lw=0.8)
axD1.text(0.66, 18, "top 0.1%", fontsize=6.0, color=GREY, va="center", ha="right")
axD1.set_xticks([]); axD1.set_ylabel("PDE3B genome-wide rank", fontsize=7.4)
top = max(24, int(ranks.max()) + 3)
axD1.set_ylim(0.3, top); axD1.set_xlim(-0.72, 0.78)
axD1.invert_yaxis(); despine(axD1, bottom=False)
axD1.set_title(f"PDE3B rank ({len(ranks)} specs)", fontsize=7.2, pad=4)

st = list(csv.DictReader(open(os.path.join(SUP, "TableS23_stratified_auroc.csv"))))
labs = [s["stratum"].replace("approved/trial drug target", "approved/\ntrial")
        .replace("Mendelian dyslipidemia", "Mendelian").replace("all canonical", "all\ncanonical") for s in st]
vals = [float(s["AUROC"]) for s in st]; ns = [s["n"] for s in st]
axD2.bar(range(len(st)), vals, color=[AMBER, CB["blue"], TEAL], width=0.62)
for xi, (v, n) in enumerate(zip(vals, ns)):
    axD2.text(xi, v + 0.004, f"{v:.2f}\n(n={n})", ha="center", fontsize=6.0, color=INK)
axD2.set_xticks(range(len(st))); axD2.set_xticklabels(labs, fontsize=6.4)
axD2.set_ylim(0.85, 1.02); axD2.set_ylabel("AUROC", fontsize=7.4); despine(axD2)
axD2.set_title("recovery by target class", fontsize=7.2, pad=3)
axD.text(-0.02, 1.15, "d", transform=axD.transAxes, fontsize=13.5, fontweight="bold",
         va="top", ha="right", color=INK)
axD.text(0.5, 1.145, "Robustness of the score", transform=axD.transAxes, fontsize=9.6,
         fontweight="bold", va="top", ha="center", color=INK)

save(fig, "fig3_convergence")
