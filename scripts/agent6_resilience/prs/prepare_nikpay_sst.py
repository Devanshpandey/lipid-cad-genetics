#!/usr/bin/env python3
"""Convert CARDIoGRAMplusC4D 2015 (Nikpay) additive sumstats to PRS-CS input.
Output columns: SNP A1 A2 BETA SE. Keeps rsIDs only, drops palindromic-ambiguous
strand SNPs with MAF > 0.4, drops SE <= 0. Prints a summary line for the manifest.
Usage: prepare_nikpay_sst.py IN OUT
"""
import sys, pandas as pd
src, dst = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, sep="\t", usecols=["markername", "effect_allele", "noneffect_allele", "effect_allele_freq", "beta", "se_dgc"], na_values=["."])
n0 = len(df)
df = df.rename(columns={"markername": "SNP", "effect_allele": "A1", "noneffect_allele": "A2", "beta": "BETA", "se_dgc": "SE"})
df = df[df["SNP"].str.startswith("rs", na=False)]
df = df.dropna(subset=["BETA", "SE"])
df = df[df["SE"] > 0]
df["A1"] = df["A1"].str.upper(); df["A2"] = df["A2"].str.upper()
df = df[df["A1"].isin(list("ACGT")) & df["A2"].isin(list("ACGT"))]
pal = ((df.A1 == "A") & (df.A2 == "T")) | ((df.A1 == "T") & (df.A2 == "A")) | ((df.A1 == "C") & (df.A2 == "G")) | ((df.A1 == "G") & (df.A2 == "C"))
amb = pal & (df["effect_allele_freq"].sub(0.5).abs() < 0.1)
df = df[~amb]
df = df.drop_duplicates("SNP")
df[["SNP", "A1", "A2", "BETA", "SE"]].to_csv(dst, sep="\t", index=False)
print(f"nikpay: {n0} rows in, {len(df)} SNPs out ({int(amb.sum())} ambiguous palindromes dropped)")
