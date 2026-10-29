"""Artefact control for the yeast condition diagonal of script 40 (exploratory, not pre-specified).

The HIP and HOP strains of a gene carry the same barcode, and HIP and HOP experiments were run in paired batches
(HIP 00n / HOP 00n+1). This control separates biology from a batch or barcode artefact: the heterozygous sensitivity is
compared in the same compound but a different batch (biology) and in the same batch but a different compound (artefact),
in addition to the diagonal itself, each against a permutation of the disease-condition sets between genes.

Input: Hoepfner HIP/HOP score tables in config.YEAST_HOEPFNER_DIR (read via yeast_data.py). Output: log only.
Call: python 40b_yeast_batch_control.py
"""
import re

import numpy as np

import yeast_data

rng = np.random.default_rng(41)
B = 200


def load(screen):
    """Scores averaged per experiment (compound number_concentration); {experiment: (compound number, batch)}."""
    d = yeast_data.read_hoepfner_scores(screen)
    meta = [re.search(r"Exp\. (\d+)_(.+?)_(HIP|HOP)_(\w+)", c).groups() for c in d.columns]
    d.columns = [f"{m[0]}_{m[1]}" for m in meta]
    return d.T.groupby(level=0).mean().T, {f"{m[0]}_{m[1]}": (m[0], m[3]) for m in meta}


hip, meta_hip = load("HIP")
hop, _ = load("HOP")
keys = sorted(set(hip.columns) & set(hop.columns))
genes = sorted(set(hip.dropna(how="all").index) & set(hop.dropna(how="all").index))
Sh, So = -hip.loc[genes, keys], -hop.loc[genes, keys]
compound = np.array([meta_hip[k][0] for k in keys])
batch = np.array([meta_hip[k][1] for k in keys])
M = (So >= 3).values & np.isfinite(Sh.values)
ok = M.sum(axis=1) >= 3
Sh, So, M = Sh[ok], So[ok], M[ok]
Vc = Sh - Sh.mean(axis=0)
Vc = Vc.sub(Vc.mean(axis=1), axis=0)
V = Vc.values
fin = np.isfinite(V)
print(f"Genes {len(M)}, experiments {len(keys)}, compounds {len(set(compound))}, batches {len(set(batch))}")


def masks(Mx):
    MS = np.zeros_like(Mx)
    MB = np.zeros_like(Mx)
    MSB = np.zeros_like(Mx)
    for i in range(len(Mx)):
        hit_compounds, hit_batches = set(compound[Mx[i]]), set(batch[Mx[i]])
        in_compound, in_batch = np.isin(compound, list(hit_compounds)), np.isin(batch, list(hit_batches))
        MS[i] = in_compound & ~in_batch & ~Mx[i]      # same compound, other batch
        MB[i] = in_batch & ~in_compound               # same batch, other compound
        MSB[i] = in_compound & in_batch & ~Mx[i]      # same compound, same batch, other dose (not a hit itself)
    return MS & fin, MB & fin, MSB & fin


MS, MB, MSB = masks(M)
obs = dict(diagonal=V[M & fin].mean(), compound_other_batch=V[MS].mean(), batch_other_compound=V[MB].mean(),
           compound_same_batch=V[MSB].mean())
n = dict(diagonal=(M & fin).sum(), compound_other_batch=MS.sum(), batch_other_compound=MB.sum(),
         compound_same_batch=MSB.sum())
perm = {k: [] for k in obs}
for b in range(B):
    P = M[rng.permutation(len(M))]
    a, bb, c = masks(P)
    for k, m in zip(obs, [P & fin, a, bb, c]):
        perm[k].append(V[m].mean())
for k in obs:
    pk = np.array(perm[k])
    print(f"  {k:26s} n {n[k]:7d}: {obs[k]:+.4f} (perm {pk.mean():+.4f} +- {pk.std():.4f}; dz {(obs[k] - pk.mean()) / pk.std():+.1f})")
