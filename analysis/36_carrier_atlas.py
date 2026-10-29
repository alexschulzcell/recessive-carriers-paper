"""Carrier phenotype atlas with a prospective prediction (Tables S7 and S5).

For every gene x trait cell of the gene-specific HPO direction map (direction_map.py, unchanged from script 31):
expected direction, observed signed z (UK Biobank pLoF / missense / synonymous), number of carriers, minimal detectable
effect (MDE, 80 % power, alpha 0.05), organ specificity (tissue, as in script 22), and a status: confirmed / opposite /
silent but measurable / not yet measurable (TET2 cells are excluded because of clonal haematopoiesis).
Prospective part: for the open cells, the predicted direction and the cohort size from which a uniform effect of
0.043 SD (median carrier effect on height, 0.40 cm / 9.3 cm) becomes detectable with 80 % power. The prediction table
is frozen with its date and SHA-256 (36_atlas_prediction_<date>.tsv; an existing file is never overwritten).

Inputs: direction map (HPO, Genebass, gene annotation), HPA tissue consensus (reference folder, hpa_tissue.zip).
Outputs: 36_carrier_atlas.tsv, 36_atlas_prediction_<today>.tsv.
Call: python 36_carrier_atlas.py
"""
import datetime
import hashlib
import zipfile

import numpy as np
import pandas as pd

import common
import config
import direction_map

dm = direction_map.build_direction_map()
GENES, Smat, TRAITS, RAW, annotation = dm.genes, dm.S, dm.traits, dm.raw, dm.annotation

TRAIT_TO_ORGAN = {t: o for o, v in common.ORGAN_TRAITS.items() for t in v}
TISSUES = {"blood": ["bone marrow"], "immune": ["bone marrow", "spleen", "lymph node"], "liver": ["liver"],
           "lipids": ["liver", "small intestine"], "kidney": ["kidney"], "muscle": ["skeletal muscle"], "lung": ["lung"],
           "circulatory": ["heart muscle", "blood vessel"], "nervous": ["cerebral cortex"], "glucose": ["pancreas"]}
# print order of the status categories (order of the original analysis)
STATUS_ORDER = ["excluded_clonal_haematopoiesis", "confirmed", "opposite", "silent_measurable", "not_yet_measurable"]

z = zipfile.ZipFile(config.REF_DIR / "hpa_tissue.zip")
hpa = pd.read_csv(z.open(z.namelist()[0]), sep="\t").pivot_table(index="Gene", columns="Tissue", values="nTPM", aggfunc="max")
hpa_median = hpa.median(axis=1)

ann = annotation.set_index("gene_id")
rows = []
for i, g in enumerate(GENES):
    for j, t in enumerate(TRAITS):
        s = Smat[i, j]
        if s == 0:
            continue
        o = TRAIT_TO_ORGAN[t]
        sp = (np.log2((hpa.loc[g, TISSUES[o]].max() + 1) / (hpa_median[g] + 1))
              if (o in TISSUES and g in hpa.index) else np.nan)
        rows.append(dict(gene_id=g, symbol=ann.symbol.get(g), trait=t, organ=o, expected_direction=int(s),
                         z_plof=RAW["plof"][i, j], z_missense=RAW["mislc"][i, j], z_synonymous=RAW["syn"][i, j],
                         carriers_plof=ann.carriers_plof.get(g), carrier_frequency_plof=ann.caf_plof.get(g),
                         tissue_specificity=sp))
X = pd.DataFrame(rows)
# z_a/2 1.96 + z_b 0.84; SE ~ 1/sqrt(carriers) for a rank-normalised trait
X["minimal_detectable_effect_sd"] = 2.80 / np.sqrt(X.carriers_plof.where(X.carriers_plof > 0))
best = np.where(np.abs(np.nan_to_num(X.z_plof)) >= np.abs(np.nan_to_num(X.z_missense)), X.z_plof, X.z_missense)
X["z_best"] = best
d = X.expected_direction * X.z_best
X["status"] = np.select([(np.abs(X.z_best) > 3) & (d > 0), (np.abs(X.z_best) > 3) & (d < 0),
                         X.minimal_detectable_effect_sd <= 0.2],
                        ["confirmed", "opposite", "silent_measurable"], "not_yet_measurable")
X.loc[X.symbol == "TET2", "status"] = "excluded_clonal_haematopoiesis"
X["specificity_tertile"] = X.groupby("organ").tissue_specificity.transform(
    lambda v: pd.qcut(v.rank(method="first"), 3, labels=[1, 2, 3]) if v.notna().sum() >= 6 else np.nan)

print("Cells per status:", X.status.value_counts().to_dict())
genes_per_status = X.groupby("status").gene_id.nunique()
print("Genes per status (at least one cell):",
      genes_per_status.reindex([s for s in STATUS_ORDER if s in genes_per_status.index]).to_dict())
k = X[np.isfinite(X.z_best) & (np.abs(X.z_best) > 2)]
print("\nCalibration (|z|>2): fraction in disease direction per specificity tertile")
for tt, x in k.groupby("specificity_tertile", observed=True):
    print(f"  tertile {tt}: {((x.expected_direction * x.z_best) > 0).mean():.3f} (n {len(x)})")

# Prospective: required carriers and cohort size. Expected effect uniformly 0.043 SD = median carrier effect on height
# (dominance analysis, script 29: -0.40 cm / 9.3 cm; not selected by significance). Effects of significant cells would be
# inflated by the winner's curse.
X["expected_effect_sd"] = 0.40 / 9.3
X["carriers_needed"] = np.ceil((2.80 / X.expected_effect_sd) ** 2)
X["cohort_size_needed"] = np.ceil(X.carriers_needed / (2 * X.carrier_frequency_plof.where(X.carrier_frequency_plof > 0)))
X["measurable_at_1_percent_today"] = X.carriers_plof >= X.carriers_needed
print(f"Cells in which a 1 % effect (0.043 SD) would already be measurable today: {int(X.measurable_at_1_percent_today.sum())}")
OPEN = X.status.isin(["not_yet_measurable", "silent_measurable"])
for label, N in [("0.5 million exomes", 5e5), ("1.4 million exomes", 1.4e6), ("5 million exomes", 5e6)]:
    print(f"  open cells (without |z|>3), testable with 80 % power from {label}: "
          f"{int((OPEN & (X.cohort_size_needed <= N)).sum())}")

cols = ["gene_id", "symbol", "organ", "trait", "expected_direction", "tissue_specificity", "specificity_tertile",
        "carriers_plof", "carrier_frequency_plof", "minimal_detectable_effect_sd", "z_plof", "z_missense", "z_synonymous",
        "status", "measurable_at_1_percent_today", "expected_effect_sd", "carriers_needed", "cohort_size_needed"]
(X[cols].sort_values(["organ", "trait", "symbol"], key=common.organ_column_key)
 .to_csv(config.results("36_carrier_atlas.tsv"), sep="\t", index=False, float_format="%.4g"))
today = datetime.date.today().isoformat()
V = X[OPEN][["gene_id", "symbol", "organ", "trait", "expected_direction", "specificity_tertile",
             "expected_effect_sd", "carriers_needed", "cohort_size_needed"]].sort_values("cohort_size_needed")
fv = config.results(f"36_atlas_prediction_{today}.tsv")
if fv.exists():
    print(f"Prediction {fv.name} already exists (frozen) -> NOT overwritten")
else:
    V.to_csv(fv, sep="\t", index=False, float_format="%.4g")
print(f"\nPrediction frozen: {fv.name}, {len(V)} cells, {V.gene_id.nunique()} genes, "
      f"SHA-256 {hashlib.sha256(fv.read_bytes()).hexdigest()}")
