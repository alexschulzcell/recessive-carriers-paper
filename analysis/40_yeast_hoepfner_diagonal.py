"""Yeast condition diagonal: heterozygous (HIP) versus homozygous (HOP) deletion strains, Hoepfner et al. 2014.

Pre-specified test (PRAEREG_hefe.md, tests T16-1 to T16-4). For every gene, the "disease conditions" are the experiments
in which its homozygous deletion strain is sensitive (HOP score <= -3). The test asks whether the heterozygous strain of
the same gene is also more sensitive in exactly these conditions (doubly centred HIP sensitivity, permutation of the
condition sets between genes). T16-3 estimates a yeast dominance coefficient h per gene, T16-4 asks whether genes with
few disease conditions (specific genes) show a stronger heterozygous diagonal. Sensitivity analyses vary the threshold,
use only the Ad. scores, and drop experiments in which more than 5 % of HOP strains are hits.

Inputs: Hoepfner HIP/HOP score tables in config.YEAST_HOEPFNER_DIR (read via yeast_data.py).
Outputs: 40_yeast_hoepfner_diagonal.tsv (one row per analysis) and 40_yeast_hoepfner_genes.tsv (per-gene values of the
primary analysis: gene, n_disease_conditions, diagonal, h_yeast; read by scripts 41 and 45).
Call: python 40_yeast_hoepfner_diagonal.py
"""
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

import config
import yeast_data

rng = np.random.default_rng(40)
B = 1000


def load(screen):
    """Score table with columns (experiment key = compound_concentration, score type)."""
    d = yeast_data.read_hoepfner_scores(screen)
    key = [re.search(r"Exp\. (.+?)_(HIP|HOP)_", c).group(1) for c in d.columns]
    score_type = [c.split(" for")[0] for c in d.columns]
    d.columns = pd.MultiIndex.from_arrays([key, score_type])
    return d


HIP, HOP = load("HIP"), load("HOP")


def prepare(ad_scores_only=False):
    hp, ho = HIP, HOP
    if ad_scores_only:
        hp = hp.loc[:, hp.columns.get_level_values(1) == "Ad. scores"]
        ho = ho.loc[:, ho.columns.get_level_values(1) == "Ad. scores"]
    hp = hp.T.groupby(level=0).mean().T   # mean per experiment key (compound_concentration)
    ho = ho.T.groupby(level=0).mean().T
    keys = sorted(set(hp.columns) & set(ho.columns))
    genes = sorted(set(hp.dropna(how="all").index) & set(ho.dropna(how="all").index))
    return -hp.loc[genes, keys], -ho.loc[genes, keys]   # sign flipped: sensitivity positive


def centre(X):
    Xc = X - X.mean(axis=0)
    return Xc.sub(Xc.mean(axis=1), axis=0)


def run(thr=3.0, ad_scores_only=False, max_frac=None, name="primary", detail=False):
    Shet, Shom = prepare(ad_scores_only)
    if max_frac is not None:
        frac = (Shom >= thr).mean(axis=0)
        keep = frac[frac <= max_frac].index
        Shet, Shom = Shet[keep], Shom[keep]
    M = (Shom >= thr).values & np.isfinite(Shet.values)          # disease conditions per gene
    ok = M.sum(axis=1) >= 3
    Shet, Shom, M = Shet[ok], Shom[ok], M[ok]
    Vc = centre(Shet).values
    Hraw = Shet.values
    fin = np.isfinite(Vc)

    def stat(Mx):
        m = Mx & fin
        return Vc[m].mean(), (Hraw[m] > 0).mean()
    o1, o2 = stat(M)
    perm = np.array([stat(M[rng.permutation(len(M))]) for _ in range(B)])
    p1 = (1 + (perm[:, 0] >= o1).sum()) / (B + 1)
    p2 = (1 + (perm[:, 1] >= o2).sum()) / (B + 1)
    print(f"[{name}] genes {len(M)}, experiments {M.shape[1]}, disease cells {M.sum()}:")
    print(f"   T16-1 HIP sensitivity in HOP disease conditions {o1:+.4f} (perm {perm[:, 0].mean():+.4f} +- {perm[:, 0].std():.4f}; "
          f"dz {(o1 - perm[:, 0].mean()) / perm[:, 0].std():+.1f}) p {p1:.4f}")
    print(f"   T16-2 fraction heterozygous also sensitive {o2:.3f} (perm {perm[:, 1].mean():.3f}) p {p2:.4f}", flush=True)
    res = dict(test=name, gene=len(M), cells=int(M.sum()), t1=o1, p1=p1, t2=o2, p2=p2)
    if detail:
        # T16-3 yeast h per gene
        gh, gd, nh = [], [], M.sum(axis=1)
        for i in range(len(M)):
            r = fin[i]
            a = np.nanmean(Hraw[i][M[i] & r]) - np.nanmean(Hraw[i][~M[i] & r])
            b = np.nanmean(Shom.values[i][M[i] & r]) - np.nanmean(Shom.values[i][~M[i] & r])
            gh.append(a / b if b > 0 else np.nan)
            gd.append(np.nanmean(Vc[i][M[i] & r]))
        gh, gd = np.array(gh), np.array(gd)
        print(f"   T16-3 yeast h: median {np.nanmedian(gh):.4f} (IQR {np.nanpercentile(gh, 25):.4f} to {np.nanpercentile(gh, 75):.4f}); "
              f"fraction h > 0: {(gh > 0).mean():.3f}")
        # T16-4 rule: specificity = few disease conditions
        r = spearmanr(gd, nh)
        tertile_labels = ["specific", "intermediate", "pleiotropic"]
        t = pd.qcut(pd.Series(nh).rank(method="first"), 3, labels=tertile_labels)
        per = {k: (np.nanmean(gd[t.values == k]), int((t.values == k).sum()), int(np.median(nh[t.values == k])))
               for k in tertile_labels}
        mw = mannwhitneyu(gd[t.values == "specific"], gd[t.values == "pleiotropic"], alternative="greater")
        print(f"   T16-4 rule: rho(gene diagonal, number of disease conditions) {r.correlation:+.3f} (p {r.pvalue:.2g}); "
              + "; ".join(f"{k} {v[0]:+.3f} (n {v[1]}, median {v[2]} cond.)" for k, v in per.items())
              + f" | MW specific > pleiotropic p {mw.pvalue:.2g}")
        res.update(h_median=np.nanmedian(gh), rho_specificity=r.correlation, p_rho=r.pvalue, p_mw=mw.pvalue)
        pd.DataFrame({"gene": Shet.index, "n_disease_conditions": nh, "diagonal": gd, "h_yeast": gh}).to_csv(
            config.results("40_yeast_hoepfner_genes.tsv"), sep="\t", index=False, float_format="%.4f")
    return res


rows = [run(detail=True)]
print("\n--- Sensitivity analyses ---")
rows.append(run(thr=2.0, name="threshold -2"))
rows.append(run(thr=4.0, name="threshold -4"))
rows.append(run(ad_scores_only=True, name="Ad. scores only"))
rows.append(run(max_frac=0.05, name="without experiments >5 % HOP hits"))
pd.DataFrame(rows).to_csv(config.results("40_yeast_hoepfner_diagonal.tsv"), sep="\t", index=False, float_format="%.4f")
