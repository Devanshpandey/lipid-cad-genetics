#!/usr/bin/env python3
"""
make_fig6_allelic.py — Rare-variant burden effects recapitulate lipid pharmacology.

Real rare-variant burden effect sizes (beta in SD units of the inverse-normal
trait, +/- 1.96*SE) for lipid drug-target genes, annotated with the matching
therapeutic. Shows genetic direction matches drug direction across the allelic
spectrum. Values from results/agent1_genetics/ burden summaries.

Output: fig6_allelic_series.{png,pdf}
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from figstyle import CB, INK, GREY, save  # house style applied on import

# gene, target lipid, burden beta (SD), SE, drug/therapeutic (matching direction)
G = [
 ("LDLR","LDL-C", 1.447,0.1050,"Statins (upregulate LDLR) lower LDL"),
 ("CETP","HDL-C", 0.914,0.0462,"CETP inhibitors raise HDL"),
 ("LPL","TG",     0.363,0.0316,"LPL pathway (APOC3/ANGPTL3i) lowers TG"),
 ("ABCA1","HDL-C",-0.343,0.0102,"Tangier gene (no approved drug)"),
 ("LPA","Lp(a)", -0.356,0.0077,"Pelacarsen (LPA ASO) lowers Lp(a)"),
 ("ANGPTL3","TG",-0.545,0.0225,"Evinacumab (ANGPTL3 mAb) lowers TG/LDL"),
 ("PCSK9","LDL-C",-0.608,0.0174,"Evolocumab (PCSK9i) lowers LDL"),
 ("LCAT","HDL-C",-1.028,0.0535,"LCAT deficiency (no approved drug)"),
 ("APOB","LDL-C",-1.118,0.0305,"Mipomersen (APOB ASO) lowers LDL"),
 ("APOC3","TG",  -1.178,0.0227,"Volanesorsen (APOC3 ASO) lowers TG"),
]
# has a directionally-matched approved/trial drug?
drug_match = {"LDLR":1,"CETP":1,"LPL":0,"ABCA1":0,"LPA":1,"ANGPTL3":1,"PCSK9":1,
              "LCAT":0,"APOB":1,"APOC3":1}

G = sorted(G, key=lambda r:r[2])
fig, ax = plt.subplots(figsize=(9.4, 4.8))
for i,(g,lip,b,se,drug) in enumerate(G):
    y=i; lo,hi=b-1.96*se,b+1.96*se
    c = CB["blue"] if b<0 else CB["red"]
    ax.plot([lo,hi],[y,y],color=c,lw=2.4,solid_capstyle="round",zorder=2)
    ax.scatter([b],[y],s=55,color=c,edgecolor="black",lw=0.5,zorder=3)
    ax.text(-2.55, y, f"{g}", ha="left", va="center", fontsize=9, fontweight="bold")
    ax.text(-2.05, y, f"({lip})", ha="left", va="center", fontsize=6.8, color=CB["grey"])
    # drug annotation on the right
    tick = r"  $\checkmark$" if drug_match[g] else ""
    ax.text(1.62, y, drug+tick, ha="left", va="center", fontsize=7.2,
            color=CB["green"] if drug_match[g] else CB["grey"])
ax.axvline(0,color=CB["ink"],ls=(0,(4,3)),lw=1)
ax.set_yticks([]); ax.set_ylim(-0.7,len(G)-0.3)
ax.set_xlim(-2.6,1.55)
ax.set_xlabel("Rare loss-of-function burden effect on target lipid (SD, 95% CI)")
ax.text(0.72, len(G)-0.45, r"raises lipid $\rightarrow$", fontsize=7.5, color=CB["red"], ha="center")
ax.text(-0.92, len(G)-0.45, r"$\leftarrow$ lowers lipid", fontsize=7.5, color=CB["blue"], ha="center")
ax.legend(handles=[Patch(color=CB["blue"],label="burden lowers lipid"),
                   Patch(color=CB["red"],label="burden raises lipid"),
                   Patch(color=CB["green"],label=r"$\checkmark$ matches drug direction")],
          loc="lower right", frameon=False, fontsize=7.5, bbox_to_anchor=(1.0,-0.02))
ax.text(0.0, 1.13, "Rare-variant allelic series recapitulates lipid pharmacology",
        transform=ax.transAxes, fontsize=11, fontweight="bold")
ax.text(0.0, 1.05, "loss-of-function direction matches the approved/trial therapeutic at each drug-target gene",
        transform=ax.transAxes, fontsize=8.2, color=GREY)
for ext in ("png","pdf"): fig.savefig(f"fig6_allelic_series.{ext}")
print("wrote fig6_allelic_series.png / .pdf")
