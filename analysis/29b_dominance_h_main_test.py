"""Dominance coefficient h: main test on height, sensitivity analyses, literature heterozygotes.

Extends script 29 in three ways:
1. Main test restricted to height (reason given in 32b: only there are the patient values clinical SDS on percentile
   charts and thus in the same units as beta_het). This is a deviation from PRAEREG_stufe9.md part C and is reported
   as such.
2. Sensitivity analysis without clipped rows (|SDS| = 8): for patients beyond the scale the denominator of h is not
   reliable.
3. Literature h computed correctly: script 29 computes h_lit = (het mean - reference mean) / reference SD, which fails
   when the heterozygote value was already extracted in SDS units (marked in the note column, see HET_SDS_MARKER).
   Those rows are used directly here.
Also translates h for height into centimetres.

Inputs: 29_dominance_h.tsv, literature_patient_values.tsv. Output: 29b_dominance_h_main_test.tsv.
Call: python 29b_dominance_h_main_test.py
"""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr

import config

# Marker written by the literature extraction into the note column when the heterozygote value is already an SDS
# (raw data; the marker text is kept as extracted).
HET_SDS_MARKER = "HET_EINHEIT=SDS"

X = pd.read_csv(config.results("29_dominance_h.tsv"), sep="\t")
L = pd.read_csv(config.results("literature_patient_values.tsv"), sep="\t")
print(f"Pairs with h: {len(X)}; genes: {X.symbol.nunique()}")


def h1(d, name):
    if len(d) < 3:
        print(f"  {name:38s} n {len(d)} -- too few")
        return
    w5 = wilcoxon(d.h - 0.5, alternative="less").pvalue
    w0 = wilcoxon(d.h, alternative="greater").pvalue
    print(f"  {name:38s} n {len(d):3d}  median h {d.h.median():+.4f}  IQR {d.h.quantile(.25):+.4f}..{d.h.quantile(.75):+.4f}"
          f"  fraction>0 {(d.h > 0).mean():.2f}  p(h<0.5) {w5:.1e}  p(h>0) {w0:.4f}")


print("\n=== H1: median h ===")
h1(X, "all traits (as script 29)")
h1(X[X.trait == "height"], "MAIN TEST: height only")
h1(X[X.trait != "height"], "all except height")
print("\n  Sensitivity analyses:")
clipped = X.sds_patients.abs() >= 7.999
h1(X[~clipped], "without clipped rows (|SDS| = 8)")
h1(X[(X.trait == "height") & ~clipped], "height, without clipped rows")
h1(X[X.het_source == "AZ"], "AZ carrier effects only (not calibrated)")
h1(X[X.patients_n >= 8], "studies with n >= 8 patients only")
h1(X[(X.trait == "height") & (X.patients_n >= 8)], "height, n >= 8 patients")

print("\n=== H2: main product genes (W1 or W2b) ===")
for sub, label in [(X, "all"), (X[X.trait == "height"], "height only")]:
    a, b = sub[sub.main == True].h, sub[sub.main != True].h  # noqa: E712
    if len(a) >= 3 and len(b) >= 3:
        print(f"  {label:14s} main product {a.median():+.4f} (n {len(a)}) vs other {b.median():+.4f} (n {len(b)}); "
              f"one-sided MW p {mannwhitneyu(a, b, alternative='greater').pvalue:.3f}")

print("\n=== H3: h ~ s_het ===")
for sub, label in [(X, "all"), (X[X.trait == "height"], "height only")]:
    y = sub.dropna(subset=["post_mean"])
    if len(y) >= 5:
        r, p = spearmanr(y.h, y.post_mean)
        print(f"  {label:14s} rho {r:+.3f}, one-sided p {p / 2 if r > 0 else 1 - p / 2:.3f} (n {len(y)})")

print("\n=== Validation with heterozygotes from the literature ===")
M = L[L.het_literature_mean.notna()].copy()
M["sds_unit"] = M.note.astype(str).str.contains(HET_SDS_MARKER, case=False)
rows = []
for _, r in M.iterrows():
    if r.sds_unit:
        sds_het = r.het_literature_mean                      # already in SDS
        kind = "literature SDS"
    elif pd.notna(r.control_or_reference_mean) and r.reference_sd and r.reference_sd > 0:
        sds_het = (r.het_literature_mean - r.control_or_reference_mean) / r.reference_sd
        kind = "from raw value"
    else:
        continue
    hh = X[(X.symbol == r.symbol) & (X.trait == r.trait)].h
    rows.append(dict(symbol=r.symbol, trait=r.trait, het_n=r.het_n, sds_het=sds_het,
                     sds_patients=r.sds_patients, h_lit=sds_het / r.sds_patients,
                     h_ours=hh.iloc[0] if len(hh) else np.nan, kind=kind))
V = pd.DataFrame(rows)
print(V.round(3).to_string(index=False) if len(V) else "  none")
if len(V):
    print(f"\n  median h_lit {V.h_lit.median():+.3f} (n {len(V)}); median of our h for the same pairs "
          f"{V.h_ours.median():+.3f}")

print("\n=== What h means in numbers (height) ===")
Hh = X[(X.trait == "height") & X.b_het.notna()].copy()
Hh["cm_carrier"] = Hh.b_het * 9.3            # UK Biobank SD of height ~9.3 cm
Hh["cm_patient"] = Hh.sds_patients * 9.3
print(f"  median carrier effect {Hh.cm_carrier.median():+.2f} cm  versus  median patient {Hh.cm_patient.median():+.1f} cm")
print(f"  range carriers: {Hh.cm_carrier.min():+.2f} to {Hh.cm_carrier.max():+.2f} cm")
print(Hh.sort_values("h", ascending=False)[["symbol", "patients_n", "sds_patients", "b_het", "cm_carrier", "h"]]
      .head(8).round(3).to_string(index=False))
X.to_csv(config.results("29b_dominance_h_main_test.tsv"), sep="\t", index=False, float_format="%.4f")
