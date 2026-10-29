"""Merge and check the literature extraction of patient values (input to the dominance coefficient h, script 29).

The literature extraction was done in five batches (skeletal A/B, blood, immune, kidney/liver/lipids), one TSV each
(data folder literature_extraction/batch*.tsv). They are merged here, brought to the columns expected by scripts 29/32b
and filtered. Every number must come with a PMID and a verbatim quotation; rows without evidence are dropped and
counted. Missing patient SDS are recalculated from patient mean, reference mean and reference SD where possible;
for duplicated gene x trait pairs the study with the largest number of patients is kept.

Input: config.LITERATURE_DIR / batch*.tsv (original column headers, translated on reading).
Output: 32_literature_values_merged.tsv (read by 32b).
Call: python 32_literature_merge_extraction.py
"""
import numpy as np
import pandas as pd

import config

# Original column headers of the extraction batches -> English column names (raw input, kept as extracted)
BATCH_COLUMNS = {"merkmal": "trait", "patient_n": "patients_n", "sds_patient": "sds_patients",
                 "patient_mittel": "patient_mean", "kontroll_oder_referenz_mittel": "control_or_reference_mean",
                 "referenz_sd": "reference_sd", "het_literatur_mittel": "het_literature_mean", "zitat": "quotation",
                 "anmerkung": "note"}
# Original organ labels of the extraction batches -> English organ labels
BATCH_ORGANS = {"Blut": "blood", "Immun": "immune", "Niere": "kidney", "Leber": "liver", "Lipide": "lipids",
                "Glukose": "glucose", "Lunge": "lung", "Muskel": "muscle", "Nerven": "nervous", "Skelett": "skeletal",
                "Kreislauf": "circulatory"}

COLUMNS = ["symbol", "organ", "trait", "patients_n", "sds_patients", "patient_mean",
           "control_or_reference_mean", "reference_sd", "het_literature_mean", "het_n",
           "pmid", "quotation", "note"]
NUMERIC = ["patients_n", "sds_patients", "patient_mean", "control_or_reference_mean",
           "reference_sd", "het_literature_mean", "het_n"]

parts = []
for f in sorted(config.LITERATURE_DIR.glob("batch*.tsv"), key=lambda p: p.name):
    t = pd.read_csv(f, sep="\t", dtype=str)
    t.columns = [c.strip() for c in t.columns]
    t = t.rename(columns=BATCH_COLUMNS)
    missing = [c for c in COLUMNS if c not in t.columns]
    t["extraction_batch"] = f.stem
    print(f"{f.name:34s} {len(t):4d} rows" + (f"  MISSING COLUMNS: {missing}" if missing else ""))
    for c in missing:
        t[c] = np.nan
    parts.append(t[COLUMNS + ["extraction_batch"]])

if not parts:
    raise SystemExit("no extraction batches found")
L = pd.concat(parts, ignore_index=True)
n0 = len(L)
for c in NUMERIC:
    L[c] = pd.to_numeric(L[c].astype(str).str.replace(",", ".").str.strip().replace({"": None, "nan": None}),
                         errors="coerce")
L["trait"] = L.trait.astype(str).str.strip().str.lower()
L["symbol"] = L.symbol.astype(str).str.strip()
L["organ"] = L.organ.astype(str).str.strip().replace(BATCH_ORGANS)

# --- checks ---
dropped = {}
no_evidence = L.pmid.isna() | L.quotation.isna() | (L.quotation.astype(str).str.len() < 10)
dropped["no PMID or quotation"] = int(no_evidence.sum())
L = L[~no_evidence]

too_small = L.patients_n.fillna(0) < 3
dropped["fewer than 3 biallelic patients"] = int(too_small.sum())
L = L[~too_small]

# recalculate SDS where not reported directly
m = L.sds_patients.isna() & L.patient_mean.notna() & L.control_or_reference_mean.notna() & L.reference_sd.gt(0)
L.loc[m, "sds_patients"] = (L.loc[m, "patient_mean"] - L.loc[m, "control_or_reference_mean"]) / L.loc[m, "reference_sd"]
print(f"\nSDS recalculated from mean/reference: {int(m.sum())}")

no_sds = L.sds_patients.isna()
dropped["no SDS and no reference to recalculate"] = int(no_sds.sum())
addable = L[no_sds][["symbol", "trait", "patient_mean", "pmid"]]
L = L[~no_sds]

# duplicates: the largest n wins
L = L.sort_values("patients_n", ascending=False).drop_duplicates(["symbol", "trait"])

print("\ndropped:", dropped)
print(f"remaining: {len(L)} gene x trait pairs from {n0} extracted rows; "
      f"{L.symbol.nunique()} genes; with heterozygote values: {int(L.het_literature_mean.notna().sum())}")
print("per organ:", L.organ.value_counts().to_dict())
print("per trait:", L.trait.value_counts().to_dict())
print(f"\n|SDS| >= 0.5 (inclusion for h): {int(L.sds_patients.abs().ge(0.5).sum())}")
print("Distribution sds_patients:", L.sds_patients.describe()[["min", "25%", "50%", "75%", "max"]].round(2).to_dict())

if len(addable):
    print(f"\nRows with a patient value but without a reference (could be completed, {len(addable)}):")
    print(addable.to_string(index=False))

out = config.results("32_literature_values_merged.tsv")
L.to_csv(out, sep="\t", index=False, float_format="%.4f")
print(f"\nwritten: {out.name}")
print("\nTop and bottom by SDS:")
print(L.sort_values("sds_patients")[["symbol", "organ", "trait", "patients_n", "sds_patients", "pmid"]]
      .head(12).round(2).to_string(index=False))
print(L.sort_values("sds_patients")[["symbol", "organ", "trait", "patients_n", "sds_patients", "pmid"]]
      .tail(8).round(2).to_string(index=False))
