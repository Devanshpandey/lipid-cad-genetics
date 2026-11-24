#!/usr/bin/env python3
"""Assemble the LDL-C vs ApoB non-identifiability table from its actual sources.

The gencov-corrected conditional F values existed only as prose in
supplementary_information.tex - no machine-readable file carried them, and
TableS3b still ships the uncorrected gencov=0 column (7.48/7.87,
weak_instrument=TRUE). After this project's cis-MR episode the standing rule is
never to hand-type a number into a figure, so these are parsed back out of the
run logs that produced them and written to a table with provenance.
"""
import csv, os, re, glob
H = os.path.expanduser("~/cad-genetics")
E = f"{H}/exports/agent1"
out = []

def mv(path, arm, ninst, log, pair):
    rows = {r["exposure"]: r for r in csv.DictReader(open(path))}
    # phenocov F from the log: "F-statistic  48.23133  49.03864  171.6264"
    fs = None
    for line in open(log, errors="ignore"):
        m = re.search(r"F-statistic\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
        if m:
            cand = [float(x) for x in m.groups()]
            # the corrected run is the LARGER of the two F triples in the log
            fs = cand if (fs is None or cand[0] > fs[0]) else fs
    for i, ex in enumerate(pair):
        r = rows[ex]
        out.append(dict(arm=arm, n_instruments=ninst, exposure=ex, model="joint (MVMR)",
                        OR=float(r["mvivw_OR"]), lo=float(r["mvivw_OR_lo95"]),
                        hi=float(r["mvivw_OR_hi95"]),
                        condF_gencov0=float(r["cond_F"]),
                        condF_phenocov=(fs[i] if fs else ""),
                        source=f"{os.path.basename(path)} + {os.path.basename(log)}"))

mv(f"{E}/mvmr_apob_gencov/mvmr_apob_CAD.csv", "biochemistry panel", 151,
   f"{H}/apobgc_3448312.out", ["LDL_C", "APOB"])
mv(f"{E}/mvmr_apob_gencov/mvmr_apob_statinnaive_CAD.csv", "statin-naive GWAS", 148,
   f"{H}/apobgcn_3448345.out", ["LDL_C", "APOB"])

# NMR: ORs from the results CSV, both F values from the nmr log line
nmr = list(csv.DictReader(open(f"{E}/nmr_gencov/nmr_mvmr_results.csv")))
g0 = ph = None
for lg in glob.glob(f"{H}/nmrgc_*.out"):
    for line in open(lg, errors="ignore"):
        m = re.search(r"gencov=0 \(LDL_P / LDL_C_NMR\):\s*([\d.]+)\s*/\s*([\d.]+)", line)
        if m: g0 = [float(x) for x in m.groups()]
        m = re.search(r"phenocov \(LDL_P / LDL_C_NMR\):\s*([\d.]+)\s*/\s*([\d.]+)", line)
        if m: ph = [float(x) for x in m.groups()]
for i, ex in enumerate(["LDL_P", "LDL_C_NMR"]):
    for model in ("MVMR", "univariable"):
        r = next(x for x in nmr if x["exposure"] == ex and x["model"] == model)
        out.append(dict(arm="NMR particle count", n_instruments=int(r["nsnp"]), exposure=ex,
                        model=("joint (MVMR)" if model == "MVMR" else "univariable"),
                        OR=float(r["OR"]), lo=float(r["lo"]), hi=float(r["hi"]),
                        condF_gencov0=(g0[i] if (g0 and model == "MVMR") else ""),
                        condF_phenocov=(ph[i] if (ph and model == "MVMR") else ""),
                        source="nmr_mvmr_results.csv + nmrgc log"))

dst = f"{E}/mvmr_apob_gencov/TableS3c_identifiability.csv"
with open(dst, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)
print(f"wrote {dst}  ({len(out)} rows)")
for r in out:
    print(f"  {r['arm']:20} {r['exposure']:10} {r['model']:14} OR={r['OR']:.3f} "
          f"({r['lo']:.3f}-{r['hi']:.3f})  F: {r['condF_gencov0'] or '-'} -> {r['condF_phenocov'] or '-'}")
