import gzip, csv, collections, statistics, random
random.seed(42)
genes={}
with gzip.open("gencode.v19.gtf.gz","rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p=line.split("\t")
        if p[2]!="gene": continue
        c=p[0][3:] if p[0].startswith("chr") else p[0]
        if c in ("M","MT"): continue
        a=p[8]; i=a.find('gene_name "')
        if i<0: continue
        n=a[i+11:a.find('"',i+11)]; s,e=int(p[3]),int(p[4]); k=(n,c)
        genes[k]=(min(genes[k][0],s),max(genes[k][1],e)) if k in genes else (s,e)
by_name=collections.defaultdict(list)
for (n,c),(s,e) in genes.items():
    try: by_name[n].append((int(c),s,e))
    except ValueError: pass

R=list(csv.DictReader(open("gene_ranking_scores_corrected.csv")))
y=[1 if str(r["canonical"]).strip() in ("1","1.0") else 0 for r in R]
tot=sum(y)
def ap(sc,reps=120):
    idx=list(range(len(R))); out=[]
    for _ in range(reps):
        random.shuffle(idx); o=sorted(idx,key=lambda i:-sc[i]); tp=0; acc=0.0
        for k,i in enumerate(o,1):
            if y[i]: tp+=1; acc+=tp/k
        out.append(acc/tot)
    return statistics.mean(out)
def r20(sc,reps=120):
    idx=list(range(len(R))); out=[]
    for _ in range(reps):
        random.shuffle(idx); out.append(sum(y[i] for i in sorted(idx,key=lambda i:-sc[i])[:20]))
    return statistics.mean(out)

sent=sorted({(int(r["chr"]),int(float(r["bp"]))) for r in csv.DictReader(open("/path/to/scratch/TableS4_coloc_loci_NEW.csv"))})
bychr=collections.defaultdict(list)
for c,b in sent: bychr[c].append(b)
regions=[]
for c,v in bychr.items():
    v=sorted(v); cur=[]
    for b in v:
        if cur and b-cur[-1][1]<=500_000: cur[-1]=[cur[-1][0],b]
        else: cur.append([b,b])
    regions += [(c,lo,hi) for lo,hi in cur]
W=500_000
def flag(drop=None):
    out=[0.0]*len(R)
    for gi,r in enumerate(R):
        hit=0
        for c,s,e in by_name.get(r["gene"],[]):
            for b in bychr.get(c,[]):
                if drop is not None and drop[0]==c and drop[1]-W<=b<=drop[2]+W: continue
                if (b-W)<=e and (b+W)>=s: hit=1; break
            if hit: break
        out[gi]=hit
    return out
def score(fl):
    return [6.0*min(float(r["lip_lp"]),60.0)/60.0 + 2*fl[i] + 2*int(r["_I_out"]) if "_I_out" in r
            else 6.0*min(float(r["lip_lp"]),60.0)/60.0 + 2*fl[i] +
                 2*int(round((float(r["integrated"])-6.0*min(float(r["lip_lp"]),60.0)/60.0-2*float(r["coloc"]))/2))
            for i,r in enumerate(R)]
base=score(flag())
bap, br = ap(base), r20(base)
print(f"full model over {len(regions)} regions: AP={bap:.3f}  recall@20={br:.1f}\n")
aps=[];rs=[]
for reg in regions:
    s=score(flag(reg)); a=ap(s,40); rr=r20(s,40); aps.append(a); rs.append(rr)
print(f"leave-one-locus-out ({len(regions)} drops):")
print(f"  AP        min={min(aps):.3f}  median={statistics.median(aps):.3f}  max={max(aps):.3f}  (full {bap:.3f})")
print(f"  recall@20 min={min(rs):.1f}  median={statistics.median(rs):.1f}  max={max(rs):.1f}  (full {br:.1f})")
worst=sorted(zip(aps,rs,regions))[:3]
for a,rr,reg in worst:
    print(f"  largest AP drop: chr{reg[0]}:{reg[1]:,}-{reg[2]:,}  AP={a:.3f} (delta {a-bap:+.3f}) recall@20={rr:.1f}")
