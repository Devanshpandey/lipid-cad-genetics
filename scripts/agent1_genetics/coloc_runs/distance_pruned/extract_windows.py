#!/usr/bin/env python3
"""Cut every summary-statistic file down to the union of coloc sentinel windows.

run_coloc.R reads whole genome-wide files; on a login node that is 15 x 467 MB of
decompression per pair. Extracting the union of +/-500 kb windows ONCE lets the
committed run_coloc.R be reused unchanged (same filenames, same REGENIE header),
so the analysis code path is identical to the cluster version.
"""
import gzip, os, glob
G = os.path.expandvars("/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged")
OUT = os.path.expandvars("$SCRATCH/coloc_dist/sumstat_windows")
os.makedirs(OUT, exist_ok=True)
W = 500_000

# union of all exposures' sentinel windows, merged per chromosome
iv = {}
for f in glob.glob("sentinels_*.tsv"):
    with open(f) as fh:
        fh.readline()
        for line in fh:
            rs, c, b, lp = line.split()
            c, b = int(c), int(b)
            iv.setdefault(c, []).append((b - W, b + W))
merged = {}
for c, v in iv.items():
    v.sort(); cur = []
    for lo, hi in v:
        if cur and lo <= cur[-1][1]: cur[-1] = (cur[-1][0], max(cur[-1][1], hi))
        else: cur.append((lo, hi))
    merged[c] = cur
tot = sum(hi-lo for v in merged.values() for lo,hi in v)
print(f"union: {sum(len(v) for v in merged.values())} merged windows, {tot/1e6:.0f} Mb", flush=True)

def inwin(c, b):
    v = merged.get(c)
    if not v: return False
    lo, hi = 0, len(v)-1
    while lo <= hi:
        m = (lo+hi)//2
        if b < v[m][0]: hi = m-1
        elif b > v[m][1]: lo = m+1
        else: return True
    return False

traits = [f"lipids_{t}" for t in ["LDL_C","HDL_C","TRIGLY","TOT_CHOL","APOA1","APOB"]] + \
         [f"outcomes_{o}" for o in ["CAD","MI","REVASC","STROKE","HF","CV_DEATH","MACE","AF","PAD"]]
for t in traits:
    dst = f"{OUT}/{t}.txt.gz"
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        print(f"  {t}: cached", flush=True); continue
    n = 0
    with gzip.open(f"{G}/{t}.txt.gz", "rt") as f, gzip.open(dst, "wt") as o:
        o.write(f.readline())
        for line in f:
            p = line.split(None, 3)
            if inwin(int(p[0]), int(p[1])):
                o.write(line); n += 1
    print(f"  {t}: {n} variants kept", flush=True)
print("EXTRACT_WINDOWS_DONE", flush=True)
