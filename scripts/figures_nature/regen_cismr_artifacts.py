#!/usr/bin/env python3
"""Regenerate every cis-MR artifact from ONE source of truth.

The cis-MR numbers appear in six places (Methods counts, Results ORs, Table S32,
Table S40, the vendored figure CSV, Figure 3c). Transcribing them by hand is how
Lp(a)'s 1.15497 was once typed as 1.155 and printed as 1.16, contradicting the
Abstract. This reads the run output and emits the derived artifacts, so the only
human step is pasting blocks that were computed, not retyped.

Usage:
  regen_cismr_artifacts.py --run cis_mr_fullpanel_corrected_union.csv [--write-csv]
"""
import argparse, csv, math, os, sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
VENDORED = os.path.join(HERE, "Table_S33_drug_target_cis_MR.csv")

# annotation columns are editorial, not computed: carry them over from the
# existing vendored file keyed on (gene, exposure) so they cannot silently drift
def load_annotations(path):
    ann = {}
    if not os.path.exists(path):
        return ann
    with open(path) as f:
        for r in csv.DictReader(f):
            ann[(r["gene"], r["exposure"])] = (
                r.get("cognate_therapeutic", ""),
                r.get("direction_concordant", ""),
            )
    return ann

EXPO_LABEL = {"LDL_C": "LDL-C", "TRIGLY": "Triglycerides",
              "LPA": "Lp(a)", "HDL_C": "HDL-C", "APOB": "ApoB"}

def fmt_p(p):
    p = float(p)
    if p >= 1e-3:
        return f"{p:.3g}"
    e = int(math.floor(math.log10(p)))
    m = p / 10**e
    return f"${m:.0f}\\times10^{{{e}}}$"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="cis_mr_*.csv from run_cis_mr.R")
    ap.add_argument("--label", default="corrected, primary")
    ap.add_argument("--write-csv", action="store_true",
                    help="overwrite the vendored figure CSV (otherwise dry-run)")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.run)))
    if not rows:
        sys.exit(f"no rows in {a.run}")

    wide = OrderedDict()
    for r in rows:
        k = (r["gene"], r["exposure"])
        w = wide.setdefault(k, {"n": r["n_instruments"]})
        if r["n_instruments"] != w["n"]:
            sys.exit(f"FATAL: {k} has n={r['n_instruments']} for {r['outcome']} "
                     f"but n={w['n']} for another outcome -- harmonization broke")
        w[r["outcome"]] = r

    ann = load_annotations(VENDORED)
    missing = [(g, e) for (g, e) in wide
               if (g, EXPO_LABEL.get(e, e)) not in ann]
    if missing:
        print(f"!! no annotation carried over for: {missing}", file=sys.stderr)

    # ---- instrument-count sentence for Methods ----
    counts = sorted({(g, int(w["n"])) for (g, e), w in wide.items()},
                    key=lambda t: -t[1])
    print("=== Methods: instrument counts ===")
    print(", ".join(f"{n} at \\textit{{{g}}}" for g, n in counts))
    singles = [g for g, n in counts if n == 1]
    print(f"single-Wald targets: {singles}")

    # ---- Table S32 LaTeX rows ----
    print("\n=== Table S32 rows ===")
    for (g, e), w in wide.items():
        if "CAD" not in w or "MI" not in w:
            print(f"%% SKIPPED {g}/{e}: missing outcome arm", file=sys.stderr)
            continue
        c, m = w["CAD"], w["MI"]
        ther, dr = ann.get((g, EXPO_LABEL.get(e, e)), ("", ""))
        print(f"\\textit{{{g}}} & {EXPO_LABEL.get(e,e)} & {w['n']} & "
              f"{float(c['OR']):.2f} ({float(c['CI_lo']):.2f}--{float(c['CI_hi']):.2f}), {fmt_p(c['P'])} & "
              f"{float(m['OR']):.2f} ({float(m['CI_lo']):.2f}--{float(m['CI_hi']):.2f}), {fmt_p(m['P'])} & "
              f"{ther} & {dr}\\\\")

    # ---- vendored CSV for Figure 3c ----
    out = []
    for (g, e), w in wide.items():
        if "CAD" not in w or "MI" not in w:
            continue
        c, m = w["CAD"], w["MI"]
        el = EXPO_LABEL.get(e, e)
        ther, dr = ann.get((g, el), ("", ""))
        out.append({
            "gene": g, "exposure": el, "n_instruments": w["n"],
            "CAD_OR": f"{float(c['OR']):.4g}", "CAD_CI_low": f"{float(c['CI_lo']):.4g}",
            "CAD_CI_high": f"{float(c['CI_hi']):.4g}", "CAD_P": f"{float(c['P']):.2g}",
            "MI_OR": f"{float(m['OR']):.4g}", "MI_CI_low": f"{float(m['CI_lo']):.4g}",
            "MI_CI_high": f"{float(m['CI_hi']):.4g}", "MI_P": f"{float(m['P']):.2g}",
            "cognate_therapeutic": ther, "direction_concordant": dr,
            "ldl_ascertainment": f"6153 U 6177 ({a.label})",
        })
    print(f"\n=== vendored CSV: {len(out)} rows "
          f"({'WRITING' if a.write_csv else 'dry-run, pass --write-csv'}) ===")
    if a.write_csv:
        with open(VENDORED, "w", newline="") as f:
            wcsv = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            wcsv.writeheader(); wcsv.writerows(out)
        print(f"wrote {VENDORED}")

if __name__ == "__main__":
    main()
