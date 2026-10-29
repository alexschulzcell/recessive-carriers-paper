"""Independent replication of the yeast condition diagonal with Hillenmeyer et al. 2008.

Pre-specified (PRAEREG_hefe.md, addendum; tests T16-R1 to T16-R3). Same design as script 40 with the heterozygous and
homozygous deletion screens of Hillenmeyer et al.: per gene, the disease compounds are those in which the homozygous strain
is sensitive (z >= 3); the heterozygous sensitivity in these compounds is compared with a permutation of the compound sets
between genes. T16-R3 estimates h as mean delta het / mean delta hom, and the per-gene diagonal is correlated with the
Hoepfner per-gene diagonal from script 40.

Inputs: Hillenmeyer het/hom z-score tables in config.YEAST_HILLENMEYER_DIR (read via yeast_data.py),
40_yeast_hoepfner_genes.tsv (script 40).
Output: 41_yeast_hillenmeyer_replication.tsv.
Call: python 41_yeast_hillenmeyer_replication.py
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import config
import yeast_data

rng = np.random.default_rng(41)
B = 1000
het, hom = yeast_data.read_hillenmeyer_z("het"), yeast_data.read_hillenmeyer_z("hom")


def experiment_meta(column):
    """(experiment ID with date, compound) from a column name such as '01_04_24_02:benomyl:6.9:um:...'."""
    return column.split(":")[0], column.split(":")[1]


def per_compound(d):
    compound = [experiment_meta(c)[1] for c in d.columns]
    return d.T.groupby(compound).mean().T


Hs, Os = per_compound(het), per_compound(hom)
compounds = sorted(set(Hs.columns) & set(Os.columns))
genes = sorted(set(Hs.dropna(how="all").index) & set(Os.dropna(how="all").index))
Hs, Os = Hs.loc[genes, compounds], Os.loc[genes, compounds]
print(f"shared compounds {len(compounds)}, genes with het+hom {len(genes)}")


def centre(X):
    Xc = X - X.mean(axis=0)
    return Xc.sub(Xc.mean(axis=1), axis=0)


def run(thr=3.0, name="primary", detail=False):
    M = (Os >= thr).values & np.isfinite(Hs.values)
    ok = M.sum(axis=1) >= 2
    h, o, M = Hs[ok], Os[ok], M[ok]
    V = centre(h).values
    R = h.values
    fin = np.isfinite(V)

    def stat(Mx):
        m = Mx & fin
        return V[m].mean(), (R[m] > 0).mean()
    o1, o2 = stat(M)
    P = np.array([stat(M[rng.permutation(len(M))]) for _ in range(B)])
    p1 = (1 + (P[:, 0] >= o1).sum()) / (B + 1)
    p2 = (1 + (P[:, 1] >= o2).sum()) / (B + 1)
    print(f"[{name}] genes {len(M)}, compounds {M.shape[1]}, disease cells {M.sum()}")
    print(f"   T16-R1 het sensitivity in hom disease compounds {o1:+.4f} (perm {P[:, 0].mean():+.4f} +- {P[:, 0].std():.4f}; "
          f"dz {(o1 - P[:, 0].mean()) / P[:, 0].std():+.1f}) p {p1:.4f}")
    print(f"   T16-R2 fraction het also sensitive {o2:.3f} (perm {P[:, 1].mean():.3f}) p {p2:.4f}", flush=True)
    out = dict(test=name, gene=len(M), cells=int(M.sum()), r1=o1, p1=p1, r2=o2, p2=p2)
    if detail:
        dh = np.array([np.nanmean(R[i][M[i] & fin[i]]) - np.nanmean(R[i][~M[i] & fin[i]]) for i in range(len(M))])
        do = np.array([np.nanmean(o.values[i][M[i] & fin[i]]) - np.nanmean(o.values[i][~M[i] & fin[i]]) for i in range(len(M))])
        print(f"   T16-R3 h = mean dhet / mean dhom = {np.nanmean(dh):.3f} / {np.nanmean(do):.3f} = {np.nanmean(dh) / np.nanmean(do):.3f}")
        gd = pd.Series([np.nanmean(V[i][M[i] & fin[i]]) for i in range(len(M))], index=h.index)
        hoepfner = pd.read_csv(config.results("40_yeast_hoepfner_genes.tsv"), sep="\t").set_index("gene").diagonal
        j = pd.concat([gd, hoepfner], axis=1, keys=["hill", "hoep"]).dropna()
        r = spearmanr(j.hill, j.hoep)
        print(f"   gene concordance with Hoepfner: rho {r.correlation:+.3f} (p {r.pvalue:.2g}, n {len(j)})")
        out.update(h=np.nanmean(dh) / np.nanmean(do), rho_hoepfner=r.correlation, p_rho=r.pvalue, n_rho=len(j))
    return out


rows = [run(detail=True), run(2.0, "threshold 2"), run(4.0, "threshold 4")]
# measurement days (yy_mm_dd) of the het and hom experiments; the measurement-day control itself is script 41b
het_days = pd.Series({c: experiment_meta(c)[0][:8] for c in het.columns})
hom_days = pd.Series({c: experiment_meta(c)[0][:8] for c in hom.columns})
print(f"\nMeasurement days: het {het_days.nunique()}, hom {hom_days.nunique()}, shared days {len(set(het_days) & set(hom_days))}")
pd.DataFrame(rows).to_csv(config.results("41_yeast_hillenmeyer_replication.tsv"), sep="\t", index=False, float_format="%.4f")
