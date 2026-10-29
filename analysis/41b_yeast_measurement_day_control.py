"""Measurement-day control for the Hillenmeyer replication of script 41.

Pre-specified (PRAEREG_hefe.md, addendum). Heterozygous experiments from the same day as a homozygous hit experiment
(hom z >= 3) but with a different compound are compared with heterozygous experiments with the same compound on other
days, each against a permutation of the hit sets between genes. Expectation: no signal on the same day with a different
compound (no day artefact).

Input: Hillenmeyer het/hom z-score tables in config.YEAST_HILLENMEYER_DIR (read via yeast_data.py). Output: log only.
Call: python 41b_yeast_measurement_day_control.py
"""
import numpy as np

import yeast_data

rng = np.random.default_rng(411)
B = 300
het, hom = yeast_data.read_hillenmeyer_z("het"), yeast_data.read_hillenmeyer_z("hom")
genes = sorted(set(het.dropna(how="all").index) & set(hom.dropna(how="all").index))
het, hom = het.loc[genes], hom.loc[genes]
# column names: "<yy_mm_dd_nn>:<compound>:..." -> compound and measurement day per experiment
het_compound = np.array([c.split(":")[1] for c in het.columns])
het_day = np.array([c.split(":")[0][:8] for c in het.columns])
hom_compound = np.array([c.split(":")[1] for c in hom.columns])
hom_day = np.array([c.split(":")[0][:8] for c in hom.columns])
V = het - het.mean(axis=0)
V = V.sub(V.mean(axis=1), axis=0)
V = V.values
fin = np.isfinite(V)
HIT = (hom.values >= 3)


def masks(hit):
    same_day = np.zeros(V.shape, bool)
    same_compound = np.zeros(V.shape, bool)
    for i in range(len(hit)):
        days, compounds = set(hom_day[hit[i]]), set(hom_compound[hit[i]])
        if not compounds:
            continue
        in_day, in_compound = np.isin(het_day, list(days)), np.isin(het_compound, list(compounds))
        same_day[i] = in_day & ~in_compound
        same_compound[i] = in_compound & ~in_day
    return same_day & fin, same_compound & fin


T, S = masks(HIT)
obs = (V[T].mean(), V[S].mean())
P = []
for b in range(B):
    t, s = masks(HIT[rng.permutation(len(HIT))])
    P.append((V[t].mean(), V[s].mean()))
P = np.array(P)
for k, label in enumerate(["same day, other compound", "same compound, other day"]):
    print(f"  {label:32s} n {[T, S][k].sum():7d}: {obs[k]:+.4f} (perm {P[:, k].mean():+.4f} +- {P[:, k].std():.4f}; "
          f"dz {(obs[k] - P[:, k].mean()) / P[:, k].std():+.1f})")
