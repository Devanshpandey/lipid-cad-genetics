#!/usr/bin/env python3
"""
make_fig5_convergence.py — Integrated common+rare evidence matrix.

Real values from results/agent1_genetics/ (burden summaries, coloc, fine-mapping).
Genes grouped: recovered drug targets / Mendelian genes; nominated non-canonical
candidates; locus-tagging signals flagged (dagger). Columns are evidence layers.

Output: fig5_convergence.{png,pdf}
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from figstyle import CB, INK, GREY, save  # house style applied on import

# (gene, group, lipid_logP, lipid_trait, outcome_logP, coloc, finemap, druggable, flag)
# group: 0=recovered target/Mendelian, 1=nominated candidate, 2=locus-tagging(flag)
G = [
 ("LDLR",0,42.4,"LDL-C",13.1,1,1,1,""),
 ("PCSK9",0,267,"LDL-C",2.7,1,1,1,""),
 ("APOB",0,293,"LDL-C",1.8,0,0,1,""),
 ("APOC3",0,587,"TG",0.6,0,0,1,""),
 ("ANGPTL3",0,129,"TG",4.0,0,0,1,""),
 ("LPA",0,463,"Lp(a)",2.9,0,0,1,""),
 ("ABCA1",0,279,"ApoA1",2.1,0,0,0,""),
 ("CETP",0,86,"HDL-C",1.6,0,0,1,""),
 ("LCAT",0,81,"HDL-C",1.3,0,0,0,""),
 ("APOA5",0,67,"TG",1.3,0,0,0,""),
 # nominated candidates
 ("PDE3B",1,51,"TG",1.2,0,0,1,""),        # cilostazol target
 ("ABO",1,56,"LDL-C",4.9,0,0,0,""),        # thrombosis/CVD pleiotropy
 ("GSDMB",1,29,"HDL-C",4.4,0,0,0,""),      # 17q inflammation/CVD
 ("ANGPTL8",1,23,"TG",0.8,0,0,0,""),       # emerging TG target
 ("PLIN1",1,23,"HDL-C",2.8,0,0,0,""),      # lipodystrophy
 ("EFCAB13",1,32,"ApoB",3.0,0,0,0,""),
 # locus-tagging (flagged)
 ("SLC22A2",2,186,"Lp(a)",1.4,0,0,0,"†"),  # within LPA locus (chr6)
 ("ZNF229",2,52,"ApoB",1.4,0,0,0,"†"),     # within APOE locus (chr19)
]
cols = ["Rare lipid\nburden", "Rare outcome\nburden", "Coloc-\nalization",
        "Fine-\nmapped", "Drug\ntarget"]
group_names = {0:"Recovered drug targets / Mendelian genes",
               1:"Nominated non-canonical candidates",
               2:"Locus-tagging (†; not independent)"}

# build display rows with group separators
rows=[]; labels=[]; groups=[]
cur=None
for g in G:
    if g[1]!=cur: rows.append(None); labels.append(group_names[g[1]]); groups.append("hdr"); cur=g[1]
    rows.append(g); labels.append(g[0]+g[8]); groups.append("gene")

n=len(rows)
fig,ax=plt.subplots(figsize=(8.6, 0.34*n+1.2))
seqcmap=LinearSegmentedColormap.from_list("bu",["#f0f4fa","#3b6fb0","#0b3d91"])

def burden_color(lp):
    v=min(np.log10(max(lp,1))/np.log10(600),1.0)
    return seqcmap(v)

for i,(row,lab,grp) in enumerate(zip(rows,labels,groups)):
    y=n-1-i
    if grp=="hdr":
        ax.text(-2.15, y, lab, ha="left", va="center", fontsize=8.2,
                fontweight="bold", color=CB["ink"])
        continue
    gene,g,lipP,ltr,outP,col,fm,drug,flag=row
    novel = (g==1); art=(g==2)
    lc = CB["red"] if novel else (CB["grey"] if art else CB["ink"])
    ax.text(-2.15, y, lab, ha="left", va="center", fontsize=8.2,
            color=lc, fontweight="bold" if novel else "normal")
    # col0 lipid burden
    ax.add_patch(plt.Rectangle((0,y-0.4),1,0.8,color=burden_color(lipP),ec="white",lw=1))
    ax.text(0.5,y,f"{lipP:.0f}",ha="center",va="center",fontsize=6.5,
            color="white" if lipP>30 else CB["ink"])
    ax.text(1.02,y,ltr,ha="left",va="center",fontsize=5.6,color=CB["grey"])
    # col1 outcome burden
    ob=burden_color(outP) if outP>=5.6 else "#eef0f2"
    ax.add_patch(plt.Rectangle((1.9,y-0.4),1,0.8,color=ob,ec="white",lw=1))
    if outP>=1: ax.text(2.4,y,f"{outP:.1f}",ha="center",va="center",fontsize=6.2,
            color="white" if outP>=5.6 else CB["grey"])
    # binary cols: coloc, finemap, drug
    for j,val in enumerate([col,fm,drug]):
        cx=3.8+j*1.0
        if val:
            ax.scatter(cx+0.5,y,s=70,color=CB["green"] if j<2 else CB["orange"],
                       edgecolor="black",lw=0.5,zorder=3)
        else:
            ax.scatter(cx+0.5,y,s=22,facecolor="none",edgecolor="#cccccc",lw=0.8)

# column headers
xpos=[0.5,2.4,4.3,5.3,6.3]
for x,c in zip(xpos,cols):
    ax.text(x,n-0.2,c,ha="center",va="bottom",fontsize=7.3,fontweight="bold")
ax.set_xlim(-2.2,7.0); ax.set_ylim(-0.6,n+0.3)
ax.axis("off")
# legend
leg=[Patch(fc=seqcmap(0.8),label="burden $-\\log_{10}P$ (darker = stronger)"),
     Patch(fc=CB["green"],label="common-variant / fine-map evidence"),
     Patch(fc=CB["orange"],label="established drug target")]
ax.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,-0.11/(0.34*n/6)),
          ncol=3,frameon=False,fontsize=6.8)
ax.set_title("Discovery-only convergence score (rare burden + colocalization; "
             "drug-target status NOT used in the score)\n"
             "enriches known lipid targets (AUROC = 0.84; top-5 = 5/5; permutation "
             "P = 5×10$^{-5}$) and nominates new candidates",
             fontsize=8.8,pad=10)
for ext in ("png","pdf"): fig.savefig(f"fig5_convergence.{ext}")
print("wrote fig5_convergence.png / .pdf")
