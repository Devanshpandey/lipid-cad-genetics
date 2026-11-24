import csv, random, statistics, collections
random.seed(42)
R=list(csv.DictReader(open("gene_ranking_scores_corrected.csv")))
y=[1 if str(r["canonical"]).strip() in ("1","1.0") else 0 for r in R]
names=[r["gene"] for r in R]

def ap_ties(s, y, reps=200):
    """Average precision with ties broken at random, averaged over reps.
    The score has large exact ties (many genes at 6.00), and a fixed tie order
    can shift AP by more than the effect being measured."""
    idx=list(range(len(s))); out=[]
    for _ in range(reps):
        random.shuffle(idx)
        o=sorted(idx, key=lambda i:-s[i])
        tp=0; tot=sum(y); acc=0.0
        for k,i in enumerate(o,1):
            if y[i]: tp+=1; acc+=tp/k
        out.append(acc/tot)
    return statistics.mean(out), statistics.stdev(out)

def recall_at(s, y, k=20, reps=200):
    idx=list(range(len(s))); out=[]
    for _ in range(reps):
        random.shuffle(idx)
        o=sorted(idx, key=lambda i:-s[i])[:k]
        out.append(sum(y[i] for i in o))
    return statistics.mean(out), min(out), max(out)

cols=[("burden_only","burden only"),("integrated","published (buggy)"),
      ("integrated_w500","corrected +/-500kb"),("integrated_w100","corrected +/-100kb")]
print(f"{'score':22s} {'AP (tie-averaged)':>20s}   {'recall@20 of 27':>22s}")
vals={}
for c,lab in cols:
    s=[float(r[c]) for r in R]; vals[c]=s
    m,sd=ap_ties(s,y); rm,rlo,rhi=recall_at(s,y)
    print(f"  {lab:20s} {m:.3f} +/- {sd:.3f}      {rm:5.1f}  (range {rlo}-{rhi})")

print("\n=== top non-canonical genes (what the paper lists as PDE3B, ANGPTL8, PLIN1, ABO, GSDMB) ===")
for c,lab in cols[1:]:
    s=vals[c]
    o=sorted(range(len(R)), key=lambda i:(-s[i], names[i]))
    top=[names[i] for i in o if not y[i]][:8]
    print(f"  {lab:22s} {', '.join(top)}")
