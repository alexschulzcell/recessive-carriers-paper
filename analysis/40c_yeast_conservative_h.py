"""Conservative yeast dominance coefficient h without shared experiments (exploratory, not pre-specified).

The disease compounds of a gene are defined from the HOP screen in a random half A of the batches (hit: HOP score
<= -3). The HIP and HOP effects are then measured in the experiments with the same compounds in the OTHER batches, so
that definition and measurement share neither a batch nor an experiment. h = delta HIP / delta HOP; five random splits.

Input: Hoepfner HIP/HOP score tables in config.YEAST_HOEPFNER_DIR (read via yeast_data.py). Output: log only.
Call: python 40c_yeast_conservative_h.py
"""
import re

import numpy as np
import pandas as pd

import yeast_data

rng = np.random.default_rng(42)


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
Sh = Sh - Sh.mean(axis=0)   # remove the experiment mean
So = So - So.mean(axis=0)
compound = np.array([meta_hip[k][0] for k in keys])
batch = np.array([meta_hip[k][1] for k in keys])
unique_batches = np.array(sorted(set(batch)))
res = []
for rep in range(5):
    half_a = set(rng.choice(unique_batches, len(unique_batches) // 2, replace=False))
    in_a = np.isin(batch, list(half_a))
    hit_a = (So.values >= 3) & in_a
    for i, g in enumerate(Sh.index):
        hit_compounds = set(compound[hit_a[i]])
        if len(hit_compounds) < 2:
            continue
        sel = np.isin(compound, list(hit_compounds)) & ~in_a
        oth = ~np.isin(compound, list(hit_compounds)) & ~in_a
        a, b = Sh.values[i], So.values[i]
        if sel.sum() < 2:
            continue
        dh = np.nanmean(a[sel]) - np.nanmean(a[oth])
        do = np.nanmean(b[sel]) - np.nanmean(b[oth])
        res.append((rep, g, dh, do))
R = pd.DataFrame(res, columns=["rep", "gene", "d_het", "d_hom"])
for rep, x in R.groupby("rep"):
    pos = x[x.d_hom > 0]
    print(f"Split {rep}: genes {len(x)}; delta het {x.d_het.mean():+.3f} (fraction > 0 {(x.d_het > 0).mean():.3f}), "
          f"delta hom {x.d_hom.mean():+.3f}; h = mean delta het / mean delta hom {x.d_het.mean() / x.d_hom.mean():.3f}; "
          f"median h per gene (delta hom > 0) {np.median(pos.d_het / pos.d_hom):.3f}")
