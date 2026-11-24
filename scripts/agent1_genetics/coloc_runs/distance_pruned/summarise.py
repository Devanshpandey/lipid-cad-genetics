#!/usr/bin/env python3
"""Corrected genome-wide coloc vs the shipped (broken) analysis."""
import csv, glob, os, collections
new = []
for f in glob.glob("out/*_coloc.csv"):
    new += list(csv.DictReader(open(f)))
old = list(csv.DictReader(open(os.path.expandvars(
    "/path/to/cad-genetics/results/agent1_genetics/coloc_results.csv"))))

def loci(rows):
    s = set()
    for r in rows:
        s.add((r["exposure"], int(r["chr"]), round(int(r["bp"])/1e6)))
    return s
def sig(rows):
    return [r for r in rows if float(r["PP.H4"]) > 0.8]

print(f"{'':22} {'SHIPPED (broken)':>18} {'CORRECTED':>12}")
print(f"{'tests':22} {len(old):>18} {len(new):>12}")
print(f"{'distinct exposure-loci':22} {len(loci(old)):>18} {len(loci(new)):>12}")
print(f"{'PP.H4 > 0.8':22} {len(sig(old)):>18} {len(sig(new)):>12}")
oc = {(r['exposure'],int(r['chr'])) for r in old}
nc = {(r['exposure'],int(r['chr'])) for r in new}
print(f"{'chromosomes covered':22} {len({c for _,c in oc}):>18} {len({c for _,c in nc}):>12}")

print("\n=== colocalized loci per exposure x outcome (corrected) ===")
cnt = collections.Counter((r["exposure"], r["outcome"]) for r in sig(new))
tested = collections.Counter((r["exposure"], r["outcome"]) for r in new)
for k in sorted(tested):
    print(f"  {k[0]:9} vs {k[1]:9} {tested[k]:4d} tested, {cnt.get(k,0):3d} colocalized")

print("\n=== newly colocalized loci NOT in the shipped analysis (top 25 by PP.H4) ===")
oldsig = loci(sig(old))
fresh = [r for r in sig(new) if (r["exposure"], int(r["chr"]), round(int(r["bp"])/1e6)) not in oldsig]
seen = set()
for r in sorted(fresh, key=lambda r: -float(r["PP.H4"])):
    k = (r["exposure"], int(r["chr"]), round(int(r["bp"])/1e6))
    if k in seen: continue
    seen.add(k)
    print(f"  {r['exposure']:9} {r['outcome']:9} chr{r['chr']}:{int(r['bp'])/1e6:8.2f}Mb "
          f"{r['sentinel_snp']:14} PP.H4={float(r['PP.H4']):.3f}")
    if len(seen) >= 25: break
print(f"\n  ({len(seen)} shown of {len({(r['exposure'],int(r['chr']),round(int(r['bp'])/1e6)) for r in fresh})} newly colocalized exposure-loci)")
