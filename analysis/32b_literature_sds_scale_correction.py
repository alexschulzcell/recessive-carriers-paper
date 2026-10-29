"""Scale correction of the patient SDS before computing the dominance coefficient h.

Deviation from PRAEREG_stufe9.md part C, explicitly marked and justified.

Problem: script 29 computes h = beta_het / sds_patients. beta_het comes from Genebass/AZ and is in SD units of a
rank-normalised (irnt) trait. In script 32, sds_patients was partly computed as (value - reference mean) / reference SD
on the raw scale, with an SD estimated from the normal range. For right-skewed laboratory values (GGT, bilirubin,
transaminases, creatinine, urate, triglycerides, cell counts) this is the wrong scale: it yields values such as
+53.7 SD (SLC25A13, GGT). In rank-normalised units such values are unreachable (in a cohort of 400k the maximum is
about |z| = 4.8). An inflated denominator pushes h towards 0 and thus in favour of hypothesis H1 ("median h clearly
below 0.5"); uncorrected, the test would be biased in our favour.

Correction:
1. Right-skewed traits are treated as log-normal: from reference mean m and reference SD s,
   sigma = sqrt(ln(1 + (s/m)^2)), mu = ln(m) - sigma^2/2, z = (ln(x) - mu) / sigma.
2. Approximately normal traits (height, haemoglobin, MCV) stay linear.
3. SDS reported by the study itself (height SDS on percentile charts) are always kept; this is the clinical standard
   and already in population SD units.
4. |z_hom| is capped at 8, and this is logged.
The main test for h is therefore restricted to height (clinical SDS, approximately normal trait, by far the most
genes); everything else is secondary. This restriction was not pre-specified and is reported as a deviation.

Input: 32_literature_values_merged.tsv. Output: literature_patient_values.tsv (read by 29, 29b, 45, 48, 49).
Call: python 32b_literature_sds_scale_correction.py
"""
import numpy as np
import pandas as pd

import config

LOG_TRAITS = {"bili", "ggt", "alt", "ast", "creat", "urate", "tg", "wbc", "lymph", "plt"}  # all other traits: linear
CAP = 8.0

L = pd.read_csv(config.results("32_literature_values_merged.tsv"), sep="\t")
L["sds_before_correction"] = L.sds_patients
L["sds_source"] = ""

# 1) keep directly reported SDS (there patient_mean is empty because nothing was recalculated)
direct = L.patient_mean.isna() & L.sds_patients.notna()
L.loc[direct, "sds_source"] = "SDS reported by the study"

# 2) recalculated values: check the scale
recalculated = ~direct & L.patient_mean.notna() & L.control_or_reference_mean.notna() & L.reference_sd.gt(0)
for i in L.index[recalculated]:
    t = L.at[i, "trait"]
    x = L.at[i, "patient_mean"]
    m = L.at[i, "control_or_reference_mean"]
    s = L.at[i, "reference_sd"]
    if t in LOG_TRAITS and x > 0 and m > 0:
        sig = np.sqrt(np.log1p((s / m) ** 2))
        mu = np.log(m) - sig ** 2 / 2
        L.at[i, "sds_patients"] = (np.log(x) - mu) / sig
        L.at[i, "sds_source"] = "recalculated, log-normal"
    else:
        L.at[i, "sds_patients"] = (x - m) / s
        L.at[i, "sds_source"] = "recalculated, linear"

print("Sources of the patient SDS:", L.sds_source.value_counts().to_dict())
unassigned = L[L.sds_source == ""]
if len(unassigned):
    print("without an assignable scale (dropped):", ", ".join(unassigned.symbol + ":" + unassigned.trait))
    L = L[L.sds_source != ""]

print("\nEffect of the correction (only rows with |sds_before_correction| > 8):")
g = L[L.sds_before_correction.abs() > CAP][["symbol", "trait", "patients_n", "sds_before_correction", "sds_patients",
                                            "sds_source"]]
print(g.sort_values("sds_before_correction", key=abs, ascending=False).round(2).to_string(index=False) if len(g) else "  none")

capped = L.sds_patients.abs() > CAP
print(f"\nafter correction still |SDS| > {CAP}: {int(capped.sum())}"
      + (f" ({', '.join(L.loc[capped, 'symbol'] + ':' + L.loc[capped, 'trait'])}) -> capped" if capped.sum() else ""))
L["sds_patients"] = L.sds_patients.clip(-CAP, CAP)

print("\nDistribution of sds_patients per trait (after correction):")
print(L.groupby("trait").sds_patients.agg(["count", "min", "median", "max"]).round(2).to_string())

out = config.results("literature_patient_values.tsv")
L.to_csv(out, sep="\t", index=False, float_format="%.4f")
print(f"\nwritten: {out.name} ({len(L)} rows, {L.symbol.nunique()} genes)")
print(f"Main-test subset height: {int((L.trait == 'height').sum())} genes; "
      f"|SDS| >= 0.5 among them: {int(((L.trait == 'height') & L.sds_patients.abs().ge(0.5)).sum())}")
