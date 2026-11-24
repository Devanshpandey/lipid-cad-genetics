# Config: deCODE cis-pQTL x UKB GWAS colocalization for lipid effector genes
import os, io, gzip, math, time, urllib.request, subprocess
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor

CISDIR = os.environ.get("CISDIR", "/path/to/scratch/cispqtl")
os.makedirs(CISDIR, exist_ok=True)
GWASDIR = "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
# Token is a per-user deCODE credential tied to their data-use agreement and it expires.
# Never hard-code it (the previous value sat in this file in the corral archive).
TOKEN = os.environ.get("DECODE_TOKEN", "")
if not TOKEN: raise SystemExit("set DECODE_TOKEN in the environment")
BASE  = "https://download.decode.is/s3/download?token=%s&file=" % TOKEN
os.chdir(CISDIR)

# gene: (deCODE file, hg38 chrom/start/end, hg19 chrom/start/end, cognate lipid GWAS)
GENES = {
 "PCSK9":  dict(f="5231_79_PCSK9_PCSK9.txt.gz",  h38=("chr1",54700000,55400000),  h19=(1,55205000,55830000),  lipid="LDL_C"),
 "APOB":   dict(f="2797_56_APOB_Apo_B.txt.gz",   h38=("chr2",20700000,21350000),  h19=(2,20900000,21570000),  lipid="APOB"),
 "ANGPTL3":dict(f="10391_1_ANGPTL3_ANGL3.txt.gz",h38=("chr1",62300000,62900000),  h19=(1,62760000,63370000),  lipid="TRIGLY"),
 "APOC3":  dict(f="6461_54_APOC3_Apo_C_III.txt.gz",h38=("chr11",116500000,117150000),h19=(11,116400000,117010000),lipid="TRIGLY"),
}
GWAS = {"LDL_C":"lipids_LDL_C.txt.gz","APOB":"lipids_APOB.txt.gz","TRIGLY":"lipids_TRIGLY.txt.gz",
        "CAD":"outcomes_CAD.txt.gz","MI":"outcomes_MI.txt.gz"}
# coloc priors + prior SDs (quantitative traits 0.15; case-control log-OR 0.20)
P1,P2,P12 = 1e-4,1e-4,1e-5
W_QUANT, W_CC = 0.15, 0.20
CC_TRAITS = {"CAD","MI"}
print("config loaded; genes:", list(GENES))

# ==== CELL ====
# Download deCODE BGZF files to disk (wget -c: resumable, survives drops) IN PARALLEL,
# then decompress from the COMPLETE file with Python gzip (BGZF-safe) and keep the cis-window.
def download(gene, meta):
    dest = "%s.full.txt.gz" % gene
    url  = BASE + meta["f"]
    subprocess.run(["wget","-c","-q","--tries=30","--timeout=90","--waitretry=5","-O",dest,url])
    return gene, dest, (os.path.getsize(dest) if os.path.exists(dest) else 0)

def parse_window(gene, meta, path):
    chrom,start,end = meta["h38"]; rows=[]; seen=False
    with gzip.open(path,"rb") as gz:
        gz.readline()  # header
        for bline in gz:
            p = bline.split(b"\t")
            if len(p) < 12: continue
            c = p[0]
            if c == chrom.encode():
                seen=True; pos=int(p[1])
                if pos < start: continue
                if pos > end: break                     # past window on target chrom -> stop
                rows.append((p[3].decode(), pos, p[4].decode(), p[5].decode(),
                             float(p[6]), float(p[9]), float(p[8]), float(p[11])))
            elif seen:
                break                                   # left target chrom block -> stop
    df = pd.DataFrame(rows, columns=["rsid","pos38","ea","oa","beta","se","mlog10p","impmaf"])
    df = df[(df.rsid.str.startswith("rs")) & (df.impmaf>=0.01) & (df.se>0)]
    df = df[(df.ea.str.len()==1) & (df.oa.str.len()==1)]
    return df.sort_values("mlog10p", ascending=False).drop_duplicates("rsid").reset_index(drop=True)

t0=time.time(); print("downloading (parallel, resumable wget -c)...")
with ThreadPoolExecutor(max_workers=4) as ex:
    dl = list(ex.map(lambda kv: download(*kv), GENES.items()))
for g,path,sz in dl: print("  %-8s %6.0f MB  (expect ~953M)" % (g, sz/1e6))
pqtl={}
for g,path,sz in dl:
    if sz < 900_000_000:
        print("  %-8s INCOMPLETE (%.0f MB) - re-run this cell; wget -c will resume" % (g, sz/1e6)); continue
    pqtl[g] = parse_window(g, GENES[g], path)
    d=pqtl[g]
    if len(d):
        print("  %-8s %6d SNPs  lead=%s (min_log10p=%.1f)" % (g,len(d),d.iloc[0].rsid,d.iloc[0].mlog10p))
        d.to_csv("pqtl_%s.tsv"%g, sep="\t", index=False)
    else:
        print("  %-8s NO WINDOW ROWS" % g)
print("pQTL done in %.0fs. Full files kept as <gene>.full.txt.gz (rm to reclaim ~3.8GB)." % (time.time()-t0))

# ==== CELL ====
# Extract UKB GWAS cis-windows (hg19) with one parallel awk pass per file.
def awk_window(trait):
    f = os.path.join(GWASDIR, GWAS[trait])
    conds=[]
    for g,m in GENES.items():
        c,s,e = m["h19"]; conds.append("($1==%d && $2>=%d && $2<=%d)" % (c,s,e))
    cond = " || ".join(conds)
    prog = "zcat %s | awk 'NR==1{next} (%s){print $3,$1,$2,$4,$5,$9,$10,$12}'" % (f,cond)
    out = subprocess.run(["bash","-c",prog], capture_output=True, text=True).stdout
    df = pd.read_csv(io.StringIO(out), sep=" ", header=None,
                     names=["rsid","chrom","pos19","a0","a1","beta","se","log10p"])
    df = df[(df.se>0)].sort_values("log10p", ascending=False).drop_duplicates("rsid")
    return trait, df.reset_index(drop=True)

t0=time.time()
with ThreadPoolExecutor(max_workers=5) as ex:
    gwas = dict(ex.map(awk_window, GWAS.keys()))
for t,df in gwas.items(): print("%-7s %6d SNPs in loci" % (t,len(df)))
print("GWAS extract done in %.0fs" % (time.time()-t0))

# ==== CELL ====
# coloc.abf (Giambartolomei 2014), pure Python. Direction-agnostic (uses z^2).
def _lse(a):
    a=np.asarray(a,float); m=a.max(); return m+math.log(np.exp(a-m).sum())
def _ldiff(a,b):  # log(exp(a)-exp(b)), a>b
    return a + math.log1p(-math.exp(b-a))
def labf(beta, se, W):
    V = se*se; r = (W*W)/(W*W+V); z2 = (beta*beta)/V
    return 0.5*(np.log(1.0-r) + r*z2)
def coloc_abf(d1, d2, W1, W2):
    # d1,d2: DataFrames with rsid,beta,se ; merge by rsid
    m = d1.merge(d2, on="rsid", suffixes=("_1","_2"))
    if len(m)==0: return dict(nsnps=0,H0=np.nan,H1=np.nan,H2=np.nan,H3=np.nan,H4=np.nan,lead=None)
    l1 = labf(m.beta_1.values, m.se_1.values, W1)
    l2 = labf(m.beta_2.values, m.se_2.values, W2)
    L1,L2,L4 = _lse(l1), _lse(l2), _lse(l1+l2)
    s12 = L1+L2                                   # guard H3 when coloc is strong (s12 ~ L4)
    lH3 = -math.inf if L4 >= s12 else s12 + math.log1p(-math.exp(L4 - s12))
    lH = [0.0, math.log(P1)+L1, math.log(P2)+L2,
          math.log(P1)+math.log(P2)+lH3, math.log(P12)+L4]
    den = _lse(lH); pp = [math.exp(x-den) for x in lH]
    lead = m.loc[(l1+l2).argmax(),"rsid"]
    return dict(nsnps=len(m),H0=pp[0],H1=pp[1],H2=pp[2],H3=pp[3],H4=pp[4],lead=lead)
print("coloc.abf ready")

# ==== CELL ====
# Run protein<->lipid, protein<->CAD, protein<->MI coloc for each gene.
try: pqtl
except NameError: pqtl = {}
for g in GENES:
    if g not in pqtl or len(pqtl.get(g,[]))==0:
        try: pqtl[g] = pd.read_csv("pqtl_%s.tsv"%g, sep="\t")
        except Exception: pqtl[g] = pd.DataFrame(columns=["rsid","beta","se","mlog10p"])
rows=[]
for g,m in GENES.items():
    dp = pqtl[g][["rsid","beta","se"]]
    if len(dp)==0: continue
    for trait in [m["lipid"],"CAD","MI"]:
        dg = gwas[trait][["rsid","beta","se"]]
        W2 = W_CC if trait in CC_TRAITS else W_QUANT
        r = coloc_abf(dp, dg, W_QUANT, W2)
        rows.append(dict(gene=g, protein_lead=pqtl[g].iloc[0].rsid,
                         protein_mlog10p=round(pqtl[g].iloc[0].mlog10p,1),
                         comparison="pQTL_vs_%s"%trait, nsnps=r["nsnps"],
                         PP_H3=round(r["H3"],3), PP_H4=round(r["H4"],3), coloc_lead=r["lead"]))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(CISDIR,"cispqtl_coloc_results.csv"), index=False)
pd.set_option("display.width",160)
print(res.to_string(index=False))
print("\nsaved -> cispqtl_coloc_results.csv")
