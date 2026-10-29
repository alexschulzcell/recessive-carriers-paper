"""Directional replication in Genes & Health 44k (pre-specified, PRAEREG_vorzeichen.md, addendum A; test T10-8a).

The replication gene statistics of scripts 25/26 store only sum z^2, i.e. no signs. Here the REGENIE gene-burden files
of Genes & Health are read again and sign(BETA) * |BETA/SE| is kept.

The gene-specific direction map is not chosen again: it is built by direction_map.build_direction_map() exactly as in
script 31; only the z matrix is exchanged (UK Biobank -> Genes & Health). Statistic and permutation as in T10-7
(mean of expected sign * centred z over the cells of the map; permutation of the rows of the direction map, B = 1000).

Masks: A = LOFTEE high-confidence pLoF, C = pLoF + missense, D = synonymous (negative control).
Filters as in script 26: test ADD, AAF < 1 %, allele count >= 5, SE > 0, z^2 capped at 50.

Inputs: Genes & Health REGENIE gene-burden files (config.GH_BURDEN_DIR), direction map (HPO, script 03 annotation).
Output: 33_genes_and_health_directional.tsv (per mask x trait x expected direction).
Call: python 33_genes_and_health_directional.py
"""
import re

import numpy as np
import pandas as pd
from scipy.stats import binomtest

import config
import direction_map

# --- identical direction map as in script 31 ---
dm = direction_map.build_direction_map()
DIRECTIONS, GENES, Smat, TRAITS = dm.directions, dm.genes, dm.S, dm.traits
rng = np.random.default_rng(5)
B = 1000
print(f"\n--- Direction map taken over: {len(DIRECTIONS)} genes, {int((Smat != 0).sum())} cells ---", flush=True)

GH_TRAITS = {"Height": "height", "Hb": "hb", "MCV": "mcv", "RBC": "rbc", "Platelets": "plt", "WBC": "wbc",
             "Lymphocytes": "lymph", "creatinine": "creat", "Urate": "urate", "Urea": "urea",
             "Phosphate": "phosphate", "ALT": "alt", "AST": "ast", "GGT": "ggt", "Bilirubin": "bili",
             "Albumin": "albumin", "LDL-C": "ldl", "HDL-C": "hdl", "Triglycerides": "tg",
             "Total_cholesterol": "chol", "Glucose_non_fasting": "glucose", "HbA1c": "hba1c",
             "FEV1": "fev1", "Systolic_BP": "sbp"}
MASK = {"A_plof": "MASK_A.0.01", "C_plof_mis": "MASK_C.0.01", "D_syn": "MASK_D.0.01"}

zdata = {m: {} for m in MASK}
for f in sorted(config.GH_BURDEN_DIR.glob("*_genetests_*.regenie.gz")):
    mt = re.search(r"2024_05_08_(.+?)\.residual_GNH", f.name)
    if not mt or mt.group(1) not in GH_TRAITS:
        continue
    trait = GH_TRAITS[mt.group(1)]
    d = pd.read_csv(f, sep=" ", comment="#", usecols=["ID", "A1FREQ", "N", "TEST", "BETA", "SE"])
    d = d[d.TEST == "ADD"]
    d["gene_id"] = d.ID.str.extract(r"\((ENSG\d+)\)")[0]
    d["mask"] = d.ID.str.extract(r"\)\.(MASK_.+)$")[0]
    d["ac"] = d.A1FREQ * 2 * d.N
    for m, mk in MASK.items():
        x = d[(d["mask"] == mk) & (d.ac >= 5) & d.SE.gt(0)].dropna(subset=["gene_id"]).drop_duplicates("gene_id")
        z = np.sign(x.BETA.values) * np.clip(np.abs(x.BETA.values / x.SE.values), None, np.sqrt(50))
        zdata[m][trait] = pd.Series(z, index=x.gene_id.values)
    print(f"  {trait:10s} " + " ".join(f"{m}:{len(zdata[m][trait])}" for m in MASK), flush=True)

ALL_RECESSIVE = sorted(dm.recessive_genes)
RAW, CEN = {}, {}
for m in MASK:
    Z = pd.DataFrame({t: s.reindex(ALL_RECESSIVE) for t, s in zdata[m].items()}, index=ALL_RECESSIVE).reindex(columns=TRAITS)
    Zc = Z - Z.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    RAW[m] = Z.loc[GENES].values
    CEN[m] = Zc.loc[GENES].values


def stat(Z, Sm, mask):
    m = mask & (Sm != 0) & np.isfinite(Z)
    return np.nansum(np.where(m, Sm * Z, 0)) / max(m.sum(), 1), int(m.sum())


def test(m, drop=None, name=""):
    Z = CEN[m]
    mask = np.isfinite(Z)
    if drop is not None:
        mask = mask & ~drop[:, None]
    obs, n = stat(Z, Smat, mask)
    pm = np.array([stat(Z, Smat[rng.permutation(len(GENES))], mask)[0] for _ in range(B)])
    dz = (obs - pm.mean()) / pm.std()
    p = (1 + (pm >= obs).sum()) / (B + 1)
    print(f"  {name or m:32s} d {obs:+.3f} (perm {pm.mean():+.3f} +- {pm.std():.3f}) dz {dz:+5.1f}  p {p:.4f}  n {n}", flush=True)
    return obs, dz, p, n


print("\n=== T10-8a: Genes & Health 44k, directional, gene-specific HPO direction ===", flush=True)
R = {}
for m, label in [("A_plof", "mask A (pLoF)"), ("C_plof_mis", "mask C (pLoF+missense)"),
                 ("D_syn", "mask D (synonymous, control)")]:
    R[m] = test(m, name=label)

print("\n=== Robustness (without top genes by |z| in mask C) ===", flush=True)
max_abs = np.nanmax(np.where(Smat != 0, np.abs(np.nan_to_num(CEN["C_plof_mis"], nan=0.0)), np.nan), axis=1)
order = np.argsort(-np.nan_to_num(max_abs, nan=-1))
for m in ["A_plof", "C_plof_mis"]:
    for k in (10, 25):
        d = np.zeros(len(GENES), bool)
        d[order[:k]] = True
        test(m, drop=d, name=f"{m} without top {k}")

print("\n=== Fraction with expected sign (raw z, permutation of the direction maps) ===", flush=True)
for m in MASK:
    Z = RAW[m]
    for lo, hi, label in [(0, np.inf, "all"), (2, np.inf, "|z| > 2"), (3, np.inf, "|z| > 3")]:
        base = np.isfinite(Z) & (np.abs(Z) >= lo) & (np.abs(Z) < hi)

        def frac(Sm):
            k = base & (Sm != 0)
            return (((Sm * Z)[k] > 0).mean() if k.sum() else np.nan), int(k.sum())
        obs, n = frac(Smat)
        if n < 20:
            continue
        pm = np.array([frac(Smat[rng.permutation(len(GENES))])[0] for _ in range(B)])
        p = (1 + (pm >= obs).sum()) / (B + 1)
        bt = binomtest(int(round(obs * n)), n, 0.5, alternative="greater").pvalue
        print(f"  {m:12s} {label:9s} {obs:.3f} (n {n}; perm {pm.mean():.3f} +- {pm.std():.3f}, "
              f"dz {(obs - pm.mean()) / pm.std():+4.1f})  p_perm {p:.4f}  p_binom {bt:.1e}", flush=True)

# per trait x expected direction, all masks
rows = []
for m in MASK:
    Z = CEN[m]
    for j, t in enumerate(TRAITS):
        for sign in (+1, -1):
            k = (Smat[:, j] == sign) & np.isfinite(Z[:, j])
            if k.sum() >= 3:
                v = sign * Z[k, j]
                rows.append(dict(mask=m, trait=t, expected_sign=sign, n=int(k.sum()),
                                 d_mean=v.mean(), fraction_positive=(v > 0).mean()))
T = pd.DataFrame(rows)
T.to_csv(config.results("33_genes_and_health_directional.tsv"), sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print("\n=== per trait x expected direction (mask C) ===", flush=True)
print(T[T["mask"] == "C_plof_mis"].sort_values("d_mean", ascending=False).round(3).to_string(index=False), flush=True)
