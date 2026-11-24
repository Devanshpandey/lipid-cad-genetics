#!/usr/bin/env python3
"""Aggregate a coloc out/ directory into one deposited CSV, WITH a header row.

The shipped TableS4_coloc_loci.csv has no header -- its first data record is
consumed as the column names -- so it is uninterpretable as deposited. Anything
written here carries an explicit header and an explicit sentinel_definition
column, so the two runs can never be confused for one another.
"""
import csv, glob, os, sys, collections

def load(d):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "out", "*_coloc.csv"))):
        with open(f) as fh:
            rows += list(csv.DictReader(fh))
    return rows

def summarise(rows, label):
    sig = [r for r in rows if str(r.get("coloc_sig", "")).strip().upper() in ("TRUE", "1")]
    loci = {r["sentinel_snp"] for r in sig}
    chrs = {r["chr"] for r in sig}
    pairs = {(r["exposure"], r["outcome"]) for r in rows}
    print(f"  {label:18s} tests={len(rows):6d}  PP.H4>0.8={len(sig):5d}  "
          f"loci={len(loci):4d}  chrom={len(chrs):3d}  pairs={len(pairs):3d}")
    return dict(tests=len(rows), sig=len(sig), loci=len(loci),
                chrom=len(chrs), pairs=len(pairs), loci_set=loci)

def write(rows, label, dst):
    if not rows:
        sys.exit(f"refusing to write an empty table for {label}")
    cols = list(rows[0].keys())
    if "sentinel_definition" not in cols:
        cols = cols + ["sentinel_definition"]
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            r = dict(r); r["sentinel_definition"] = label
            w.writerow(r)
    # deposited files must be re-readable with their own header
    with open(dst) as fh:
        back = list(csv.DictReader(fh))
    assert len(back) == len(rows) and "sentinel_snp" in back[0], "header round-trip failed"
    print(f"  wrote {dst}  ({len(rows)} rows, {len(cols)} cols, header verified)")

if __name__ == "__main__":
    S = os.path.expandvars("$SCRATCH")
    stats = {}
    for d, label in (("coloc_dist", "distance_pruned"), ("coloc_clump", "ld_clumped")):
        path = os.path.join(S, d)
        rows = load(path)
        if not rows:
            print(f"  {label:18s} (no results yet)"); continue
        stats[label] = summarise(rows, label)
        write(rows, label, os.path.join(S, f"TableS4_coloc_{label}.csv"))
    if len(stats) == 2:
        a, b = stats["distance_pruned"], stats["ld_clumped"]
        shared = a["loci_set"] & b["loci_set"]
        print(f"\n  concordance: {len(shared)} colocalized loci in BOTH definitions")
        print(f"    distance-pruned only: {len(a['loci_set']-b['loci_set'])}")
        print(f"    LD-clumped only     : {len(b['loci_set']-a['loci_set'])}")
