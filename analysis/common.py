"""Definitions shared by several analysis scripts: organ labels, their fixed order, and the frozen atlas prediction.

Organ order
-----------
Several permutation analyses loop over organs and draw random permutations inside that loop. The p-values are only
reproducible if the organs are visited in the same order as in the original analysis. ORGAN_ORDER fixes that order;
it is used wherever a script sorts organ labels (``sort_organs`` and ``organ_sort_key``). It is deliberately not the
alphabetical order of the English labels.
"""
import hashlib

import pandas as pd

import config

# Organ systems in the order in which they are defined in script 14 (HPO subtree per organ).
ORGAN_HPO_ROOTS = {
    "blood": ["HP:0001871"], "immune": ["HP:0002715"], "kidney": ["HP:0000077"], "liver": ["HP:0001392"],
    "lipids": ["HP:0003119"], "glucose": ["HP:0011014"], "lung": ["HP:0002086"], "muscle": ["HP:0003011"],
    "nervous": ["HP:0000707"], "skeletal": ["HP:0000924"], "circulatory": ["HP:0001626"],
}
# Traits of each organ (a priori pairs, script 14).
ORGAN_TRAITS = {
    "blood": ["hb", "mcv", "rbc", "plt"], "immune": ["wbc", "lymph"], "kidney": ["creat", "cysc", "urate", "urea", "phosphate"],
    "liver": ["alt", "ast", "ggt", "bili", "albumin"], "lipids": ["ldl", "hdl", "tg", "chol"], "glucose": ["glucose", "hba1c"],
    "lung": ["fev1", "fvc"], "muscle": ["grip"], "nervous": ["fluidint", "reaction"], "skeletal": ["height", "sitheight", "ebmd"],
    "circulatory": ["sbp"],
}

# Fixed sort order of organ labels (see module docstring).
ORGAN_ORDER = ("blood", "glucose", "immune", "circulatory", "liver", "lipids", "lung", "muscle", "nervous", "kidney", "skeletal")


def organ_sort_key(organ):
    return ORGAN_ORDER.index(organ)


def sort_organs(organs):
    """Sort organ labels in the fixed order ORGAN_ORDER."""
    return sorted(organs, key=organ_sort_key)


def organ_column_key(col):
    """``key`` argument for DataFrame.sort_values: organ columns in ORGAN_ORDER, all other columns unchanged."""
    return col.map(organ_sort_key) if col.name == "organ" else col


# The frozen prediction was written on 2026-10-20 and is kept byte-identical (SHA-256 checked below), with its
# original column names and organ labels. They are translated here when the file is read.
_FROZEN_COLUMNS = {"merkmal": "trait", "richtung_erwartet": "expected_direction", "spez_tertil": "specificity_tertile",
                   "effekt_erwartet_sd": "expected_effect_sd", "traeger_noetig": "carriers_needed",
                   "kohorte_noetig": "cohort_size_needed"}
_FROZEN_ORGANS = {"Blut": "blood", "Immun": "immune", "Niere": "kidney", "Leber": "liver", "Lipide": "lipids",
                  "Glukose": "glucose", "Lunge": "lung", "Muskel": "muscle", "Nerven": "nervous", "Skelett": "skeletal",
                  "Kreislauf": "circulatory"}


def load_frozen_prediction():
    """Read the frozen atlas prediction (Table S5) and print its SHA-256 prefix."""
    sha = hashlib.sha256(config.FROZEN_PREDICTION.read_bytes()).hexdigest()[:16]
    print(f"  SHA-256 of the prediction: {sha} (frozen: {config.FROZEN_PREDICTION_SHA256_PREFIX})")
    v = pd.read_csv(config.FROZEN_PREDICTION, sep="\t").rename(columns=_FROZEN_COLUMNS)
    v["organ"] = v["organ"].map(_FROZEN_ORGANS)
    return v
