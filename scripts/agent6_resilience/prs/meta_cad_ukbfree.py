#!/usr/bin/env python3
"""Fixed-effect inverse-variance meta-analysis of three UKB-free CAD GWAS,
matched on rsID with allele harmonisation, written directly in PRS-CS format.

  MVP EUR  (Verma 2024, phs002453, phecode 411.4)  OR + CI  -> beta, se
  CARDIoGRAMplusC4D 2015 (Nikpay)                  beta, se
  FinnGen R11 I9_CHD                               beta, se (rsids may be comma-joined)

Rules: rsIDs only; ACGT alleles only; strand-ambiguous (A/T, C/G) SNPs dropped in
every study; alleles aligned to the MVP effect allele (or to the first study
carrying the SNP); a SNP is kept if present in >= 2 studies. Output columns:
SNP A1 A2 BETA SE plus a per-SNP N_EFF and a summary JSON for the manifest.

Usage: meta_cad_ukbfree.py --mvp F --nikpay F --finngen F --out F.txt --n_finngen_cases N --n_finngen_controls N
"""
import argparse, json, math, sys
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--mvp", required=True); ap.add_argument("--nikpay", required=True); ap.add_argument("--finngen", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--n_finngen_cases", type=int, required=True); ap.add_argument("--n_finngen_controls", type=int, required=True)
ap.add_argument("--n_nikpay_cases", type=int, default=60801); ap.add_argument("--n_nikpay_controls", type=int, default=123504)
a = ap.parse_args()
ACGT = set("ACGT")
def clean(df):
    df = df[df.SNP.str.startswith("rs", na=False)].copy()
    df["A1"] = df.A1.str.upper(); df["A2"] = df.A2.str.upper()
    df = df[df.A1.isin(ACGT) & df.A2.isin(ACGT) & (df.A1 != df.A2)]
    amb = ((df.A1 == "A") & (df.A2 == "T")) | ((df.A1 == "T") & (df.A2 == "A")) | ((df.A1 == "C") & (df.A2 == "G")) | ((df.A1 == "G") & (df.A2 == "C"))
    df = df[~amb]
    df = df[np.isfinite(df.BETA) & np.isfinite(df.SE) & (df.SE > 0)]
    return df.drop_duplicates("SNP")
def neff(nc, nk): return 4.0 / (1.0 / nc + 1.0 / nk)

# --- MVP ---
m = pd.read_csv(a.mvp, sep="\t", usecols=["SNP_ID", "ref", "alt", "ea", "or", "ci", "num_cases", "num_controls"], dtype={"ci": str})
m = m[m["or"] > 0]
ci = m.ci.str.split(",", expand=True).astype(float)
m["BETA"] = np.log(m["or"]); m["SE"] = (np.log(ci[1]) - np.log(ci[0])) / (2 * 1.959964)
m["A1"] = m.ea; m["A2"] = np.where(m.ea == m.alt, m.ref, m.alt)
m = m.rename(columns={"SNP_ID": "SNP"})
n_mvp = neff(int(m.num_cases.iloc[0]), int(m.num_controls.iloc[0]))
m = clean(m[["SNP", "A1", "A2", "BETA", "SE"]]); m["NEFF"] = n_mvp
print(f"MVP: {len(m)} SNPs, Neff={n_mvp:.0f}", file=sys.stderr)

# --- Nikpay ---
k = pd.read_csv(a.nikpay, sep="\t", usecols=["markername", "effect_allele", "noneffect_allele", "beta", "se_dgc"], na_values=["."])
k = k.rename(columns={"markername": "SNP", "effect_allele": "A1", "noneffect_allele": "A2", "beta": "BETA", "se_dgc": "SE"})
k = clean(k); k["NEFF"] = neff(a.n_nikpay_cases, a.n_nikpay_controls)
print(f"Nikpay: {len(k)} SNPs", file=sys.stderr)

# --- FinnGen ---
f = pd.read_csv(a.finngen, sep="\t", usecols=["ref", "alt", "rsids", "beta", "sebeta"])
f = f.dropna(subset=["rsids"])
f = f.assign(SNP=f.rsids.str.split(",")).explode("SNP")
f = f.rename(columns={"alt": "A1", "ref": "A2", "beta": "BETA", "sebeta": "SE"})[["SNP", "A1", "A2", "BETA", "SE"]]
f = clean(f); f["NEFF"] = neff(a.n_finngen_cases, a.n_finngen_controls)
print(f"FinnGen: {len(f)} SNPs", file=sys.stderr)

# --- align to reference alleles (first study that has the SNP, MVP priority) ---
ref = pd.concat([m[["SNP", "A1", "A2"]], k[["SNP", "A1", "A2"]], f[["SNP", "A1", "A2"]]]).drop_duplicates("SNP").set_index("SNP")
def align(df, name):
    df = df.join(ref, on="SNP", rsuffix="_ref")
    same = (df.A1 == df.A1_ref) & (df.A2 == df.A2_ref)
    flip = (df.A1 == df.A2_ref) & (df.A2 == df.A1_ref)
    n_drop = int((~same & ~flip).sum())
    df = df[same | flip].copy()
    df.loc[flip[same | flip], "BETA"] *= -1
    df["w"] = 1.0 / df.SE ** 2
    print(f"{name}: {n_drop} allele-mismatch SNPs dropped", file=sys.stderr)
    return df[["SNP", "BETA", "w", "NEFF"]]
parts = [align(m, "MVP"), align(k, "Nikpay"), align(f, "FinnGen")]
allp = pd.concat(parts)
g = allp.assign(bw=allp.BETA * allp.w).groupby("SNP").agg(bw=("bw", "sum"), w=("w", "sum"), n=("w", "size"), NEFF=("NEFF", "sum"))
g = g[g.n >= 2]
g["BETA"] = g.bw / g.w; g["SE"] = 1.0 / np.sqrt(g.w)
out = g.join(ref).reset_index()[["SNP", "A1", "A2", "BETA", "SE", "NEFF", "n"]]
out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
summary = {"n_snps": int(len(out)), "n_in_3_studies": int((out.n == 3).sum()), "n_in_2_studies": int((out.n == 2).sum()),
           "neff_mvp": n_mvp, "neff_nikpay": neff(a.n_nikpay_cases, a.n_nikpay_controls), "neff_finngen": neff(a.n_finngen_cases, a.n_finngen_controls),
           "neff_total_max": float(out.NEFF.max()), "neff_median": float(out.NEFF.median()),
           "n_gwas_for_prscs": int(round(out.NEFF.median()))}
json.dump(summary, open(a.out + ".summary.json", "w"), indent=2)
print(json.dumps(summary), file=sys.stderr)
