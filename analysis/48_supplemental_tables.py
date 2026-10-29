"""Supplemental Tables S3-S7 in their published form, assembled from the pipeline outputs.

Only column/category renaming, an added plain-text trait name next to each trait code, and the English free-text notes
of Table S4 (from data/literature_notes_english.tsv); all numbers are copied unchanged as text.
The frozen prediction (Table S5) is kept byte-identical with its original column names
(config.FROZEN_PREDICTION); the English version written here is a copy with translated headers and organ labels.

Inputs (results folder): 36_carrier_atlas.tsv, 08_gene_x_phenotype_table.tsv, literature_patient_values.tsv,
47_st7a_gene_classes.tsv, 47_st7b_organ_to_mp_terms.tsv, 47_st7c_organ_to_icd10_blocks.tsv; config.FROZEN_PREDICTION;
data/literature_notes_english.tsv.
Outputs: results/supplemental_tables/Table_S{3A,3B,3C,4,5,6,7}_*.tsv
Call: python 48_supplemental_tables.py
"""
import pandas as pd

import common
import config

OUT = config.results("supplemental_tables")
OUT.mkdir(parents=True, exist_ok=True)
ENGLISH_TABLES = config.REPO_DIR / "data"

TRAIT_NAME = {"hb": "haemoglobin", "mcv": "mean corpuscular volume", "rbc": "red blood cell count", "plt": "platelet count",
              "wbc": "white blood cell count", "lymph": "lymphocyte count", "creat": "creatinine", "cysc": "cystatin C",
              "urate": "urate", "urea": "urea", "phosphate": "phosphate", "alt": "alanine aminotransferase",
              "ast": "aspartate aminotransferase", "ggt": "gamma-glutamyltransferase", "bili": "bilirubin",
              "albumin": "albumin", "ldl": "LDL cholesterol", "hdl": "HDL cholesterol", "tg": "triglycerides",
              "chol": "total cholesterol", "glucose": "glucose", "hba1c": "HbA1c", "fev1": "FEV1", "fvc": "FVC",
              "grip": "grip strength", "fluidint": "fluid intelligence", "reaction": "reaction time", "height": "height",
              "sitheight": "sitting height", "ebmd": "heel bone mineral density", "sbp": "systolic blood pressure"}
STATUS = {"confirmed", "opposite", "silent_measurable", "not_yet_measurable", "excluded_clonal_haematopoiesis"}
SDS_SOURCES = {"SDS reported by the study", "recalculated, log-normal", "recalculated, linear"}


def read_text(path):
    """Read a TSV with all values as text (numbers are copied exactly as written by the pipeline)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write(df, name):
    df.to_csv(OUT / name, sep="\t", index=False, lineterminator="\n")
    print("written:", name, df.shape)


def add_trait_name(df):
    """Keep the trait code and insert its plain-text name right after it."""
    df = df.copy()
    df.insert(df.columns.get_loc("trait") + 1, "trait_name", df["trait"].map(TRAIT_NAME))
    return df


# ST2 carrier atlas
t = add_trait_name(read_text(config.results("36_carrier_atlas.tsv")))
assert t["status"].isin(STATUS).all() and t["trait_name"].notna().all()
write(t, "Table_S7_carrier_atlas.tsv")

# ST3 frozen prediction: headers and organ labels translated on reading (same translation as
# common.load_frozen_prediction, but all values kept as text); the original file stays unchanged
t = read_text(config.FROZEN_PREDICTION).rename(columns=common._FROZEN_COLUMNS)
t["organ"] = t["organ"].map(lambda v: common._FROZEN_ORGANS.get(v, v))
t = add_trait_name(t)
write(t, "Table_S5_frozen_atlas_prediction_2026-10-20.tsv")

# ST4 skeletal gene x phenotype table
t = read_text(config.results("08_gene_x_phenotype_table.tsv"))
t["set"] = t["set"].map({"AR_nos": "recessive_dysplasia", "AD_nos": "dominant_dysplasia", "nos_other": "other_dysplasia"})
t = t.rename(columns={"set": "gene_set", "nos_inh": "nosology_inheritance", "nos_groups": "nosology_groups",
                      "hpo_ad": "dominant_disease_in_hpo", "oe_lof_upper": "loeuf"})
t.columns = [c.replace("qAR_", "q_within_recessive_").replace("_mislc", "_missense") for c in t.columns]
assert t["gene_set"].notna().all()
write(t, "Table_S6_skeletal_gene_x_phenotype.tsv")

# ST5 literature patient values
t = add_trait_name(read_text(config.results("literature_patient_values.tsv")))
notes = read_text(ENGLISH_TABLES / "literature_notes_english.tsv")
assert len(notes) == len(t) and list(notes["row"].astype(int)) == list(range(len(t)))
t["note"] = notes["note_en"].values
# Data correction: in one quotation the extraction had inserted a German explanatory prefix that is not part of the
# quoted source text; it is replaced by its English translation.
t["quotation"] = t["quotation"].str.replace("zwei norwegische Kinder:", "two Norwegian children:", regex=False)
assert t["sds_source"].isin(SDS_SOURCES).all() and (t["extraction_batch"] != "").all()
write(t, "Table_S4_literature_patient_values.tsv")

# ST7a-c definitions (written in English by script 47)
write(read_text(config.results("47_st7a_gene_classes.tsv")), "Table_S3A_gene_classes.tsv")
write(read_text(config.results("47_st7b_organ_to_mp_terms.tsv")), "Table_S3B_organ_to_MP_terms.tsv")
write(read_text(config.results("47_st7c_organ_to_icd10_blocks.tsv")), "Table_S3C_organ_to_ICD10_blocks.tsv")
