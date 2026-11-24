#!/usr/bin/env python3
"""Recompute the convergence score's colocalization layer against the corrected
LD-clumped sentinel set, using a genuine GRCh37 gene annotation.

The published coloc flag was computed from the buggy May-24 export, whose
significant sentinels lay on chromosomes 1, 8 and 19 only; it fires for 8 genes
at 5 loci. The score formula itself is unchanged and is the one documented in
score_sensitivity.R:  integrated = 6*min(lip_lp,60)/60 + 2*coloc + 2*I_out
"""
import gzip, csv, collections, sys

# ---- 1. GRCh37 gene bodies from GENCODE v19 --------------------------------
genes = {}
with gzip.open("gencode.v19.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split("\t")
        if p[2] != "gene": continue
        chrom = p[0][3:] if p[0].startswith("chr") else p[0]
        if chrom in ("M","MT"): continue
        s, e = int(p[3]), int(p[4])
        attr = p[8]
        i = attr.find('gene_name "')
        if i < 0: continue
        name = attr[i+11:attr.find('"', i+11)]
        # merge multiple entries for one symbol on the same chromosome
        k = (name, chrom)
        if k in genes:
            a, b = genes[k]; genes[k] = (min(a, s), max(b, e))
        else:
            genes[k] = (s, e)
by_name = collections.defaultdict(list)
for (name, chrom), (s, e) in genes.items():
    by_name[name].append((chrom, s, e))
print(f"[annot] GENCODE v19: {len(genes)} gene records, {len(by_name)} distinct symbols", flush=True)

# ---- 2. corrected colocalized sentinels -------------------------------------
sent = sorted({(int(r["chr"]), int(float(r["bp"])))
               for r in csv.DictReader(open("/path/to/scratch/TableS4_coloc_loci_NEW.csv"))})
by_chr = collections.defaultdict(list)
for c, b in sent: by_chr[c].append(b)
print(f"[sent ] corrected LD-clumped colocalized sentinels: {len(sent)} on {len(by_chr)} chromosomes", flush=True)

def in_locus(name, pad):
    """gene body overlaps [sentinel-pad, sentinel+pad] for any colocalized sentinel"""
    for chrom, s, e in by_name.get(name, []):
        try: c = int(chrom)
        except ValueError: continue
        for b in by_chr.get(c, []):
            if (b - pad) <= e and (b + pad) >= s:
                return 1
    return 0

# ---- 3. rebuild the score ---------------------------------------------------
rows = list(csv.DictReader(open("gene_ranking_scores.csv")))
unmapped = [r["gene"] for r in rows if r["gene"] not in by_name]
print(f"[map  ] {len(rows)-len(unmapped)}/{len(rows)} scored genes mapped to GENCODE v19 "
      f"({len(unmapped)} unmapped)", flush=True)

def bt(lip): return 6.0 * min(float(lip), 60.0) / 60.0

# I_out is not stored; recover it exactly as score_sensitivity.R does
bad = 0
for r in rows:
    b = bt(r["lip_lp"]); c = float(r["coloc"])
    io = (float(r["integrated"]) - b - 2*c) / 2
    r["_I_out"] = int(round(io))
    if abs(io - round(io)) > 1e-6 or r["_I_out"] not in (0,1): bad += 1
print(f"[check] I_out recovered for all rows; non-integer/out-of-range: {bad}", flush=True)
recon_err = max(abs(bt(r["lip_lp"]) + 2*float(r["coloc"]) + 2*r["_I_out"] - float(r["integrated"])) for r in rows)
print(f"[check] formula reproduces published integrated score: max |err| = {recon_err:.3g}", flush=True)

for pad, tag in ((500_000, "w500"), (100_000, "w100"), (0, "body")):
    n = 0
    for r in rows:
        v = in_locus(r["gene"], pad); r["coloc_" + tag] = v; n += v
    print(f"[coloc] pad={pad//1000:>3d}kb -> coloc=1 for {n} of {len(rows)} genes "
          f"({sum(1 for r in rows if r['coloc_'+tag] and str(r['canonical']).strip() in ('1','1.0'))} of 27 canonical)", flush=True)
    for r in rows:
        r["integrated_" + tag] = round(bt(r["lip_lp"]) + 2*r["coloc_"+tag] + 2*r["_I_out"], 6)

cols = list(rows[0].keys())
with open("gene_ranking_scores_corrected.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
print("[out  ] wrote gene_ranking_scores_corrected.csv", flush=True)
