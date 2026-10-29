"""Shared code of the AlphaMissense dose analyses (scripts 42 and 43; PRAEREG_alphamissense.md, tests T17-1 to T17-3).

Rare variants (MAF < 1 %) of independent cohorts are assigned to recessive genes with a gene-specific expected direction
(direction_map.py) for the trait in question. Missense variants get their AlphaMissense pathogenicity score, pLoF
variants (gnomAD LOFTEE HC) are kept as reference. For every variant, d = expected sign x z (z oriented to the
AlphaMissense alternative allele, or to the LoF allele for pLoF). The dose test asks whether, within one cohort x gene x
trait cell, more pathogenic missense variants have a larger d (pooled within-cell regression slope, permutation of the
scores within cells).
"""
import numpy as np
import pandas as pd

CAP = np.sqrt(50)
CLASS_LABELS = ["benign", "ambiguous", "pathogen"]   # AlphaMissense classes (cut-offs 0.34 and 0.564)


def build_variant_table(sources, builds, annotation, directions, traits, alphamissense_dir):
    """Variants in direction genes with z, expected sign s, pLoF flag and AlphaMissense score (am, amflip)."""
    rows = []
    for coh, items in sources.items():
        for k, build, fn in items:
            if k not in traits:
                continue
            genes_ok = {g for g, r in directions.items() if k in r}
            df = fn()
            df = df[np.isfinite(df.z) & (df.maf < 0.01)]
            lof = annotation.lof19 if build == "hg19" else annotation.lof38
            n0 = len(rows)
            for ch, pos, ref, alt, z in zip(df.ch.astype(str), df.pos.astype(int), df.ref.astype(str), df.alt.astype(str), df.z):
                gs = [annotation.sym2id.get(s) for s in annotation.coding_genes(build, ch, pos)]
                gs = [g for g in gs if g in genes_ok]
                if not gs:
                    continue
                lof_gene = lof.get((ch, pos, ref, alt))
                lof_flip = False
                if lof_gene is None and lof.get((ch, pos, alt, ref)):
                    lof_gene = lof[(ch, pos, alt, ref)]
                    lof_flip = True
                for g in gs:
                    rows.append((coh, build, ch, pos, ref, alt, float(np.clip(z, -CAP, CAP)), g, k, directions[g][k],
                                 (lof_gene == g), lof_flip))
            print(f"{coh} {k}: {len(rows) - n0} variants in direction genes", flush=True)
    X = pd.DataFrame(rows, columns=["coh", "build", "ch", "pos", "ref", "alt", "z", "gene_id", "trait", "s", "plof", "lflip"])
    # stream AlphaMissense per build and join
    am_parts = []
    for build in builds:
        sub = X[X.build == build]
        keys = set(zip(sub.ch, sub.pos))
        f = alphamissense_dir / f"AlphaMissense_{build}.tsv.gz"
        for ch in pd.read_csv(f, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3, 8],
                              names=["ch", "pos", "am_ref", "am_alt", "am"], dtype={"ch": str}, chunksize=5_000_000):
            ch["ch"] = ch.ch.str.replace("chr", "", regex=False)
            m = [(a, b) in keys for a, b in zip(ch.ch, ch.pos)]
            am_parts.append(ch[m].assign(build=build))
        print(f"AlphaMissense {build}: matching rows so far {sum(len(p) for p in am_parts)}", flush=True)
    # overlapping genes: highest score
    AMd = pd.concat(am_parts).groupby(["build", "ch", "pos", "am_ref", "am_alt"], as_index=False).am.max()
    a = X.merge(AMd, left_on=["build", "ch", "pos", "ref", "alt"], right_on=["build", "ch", "pos", "am_ref", "am_alt"], how="left")
    b = X.merge(AMd, left_on=["build", "ch", "pos", "alt", "ref"], right_on=["build", "ch", "pos", "am_ref", "am_alt"], how="left")
    X["am"] = a.am.values
    flip = X.am.isna() & b.am.notna().values
    X.loc[flip, "am"] = b.am.values[flip.values]
    X["amflip"] = flip
    return X


def load_variant_table(cache_file, **kwargs):
    """Read the variant table from cache_file, or build it (build_variant_table) and write it there."""
    if cache_file.exists():
        return pd.read_csv(cache_file, sep="\t", dtype={"ch": str})
    X = build_variant_table(**kwargs)
    X.to_csv(cache_file, sep="\t", index=False)
    return X


def orient(X):
    """Orientation: missense to the AlphaMissense alternative allele, pLoF to the LoF allele; d = s * oriented z."""
    X["zo"] = np.where(X.amflip, -X.z, X.z)
    X.loc[X.plof, "zo"] = np.where(X.loc[X.plof, "lflip"], -X.loc[X.plof, "z"], X.loc[X.plof, "z"])
    X["d"] = X.s * X.zo
    X["cell"] = X.coh + "|" + X.gene_id + "|" + X.trait
    MS = X[X.am.notna() & ~X.plof].copy()
    print(f"\nMissense with AlphaMissense: {len(MS)} variants in {MS.cell.nunique()} cells; pLoF: {int(X.plof.sum())}")
    print("per source:", MS.groupby("coh").size().to_dict())
    return X, MS


def within_slope(df, col="d"):
    """Slope of col on am with a fixed effect per cell (pooled within regression)."""
    am_c = df.am - df.groupby("cell").am.transform("mean")
    y_c = df[col] - df.groupby("cell")[col].transform("mean")
    return (am_c * y_c).sum() / (am_c ** 2).sum()


def perm_within(df, rng, col="d", B=1000):
    """Observed within-cell slope, null distribution (scores permuted within cells) and one-sided p."""
    obs = within_slope(df, col)
    g = df.groupby("cell").indices
    am = df.am.values.copy()
    dist = np.empty(B)
    for b in range(B):
        a2 = am.copy()
        for ix in g.values():
            a2[ix] = am[rng.permutation(ix)]
        dist[b] = within_slope(df.assign(am=a2), col)
    return obs, dist, (1 + (dist >= obs).sum()) / (B + 1)


def dose_tests(X, MS, rng, directions):
    """Tests T17-1 to T17-3. Returns the result rows and the missense table restricted to cells with >= 2 variants."""
    MS = MS[MS.groupby("cell").am.transform("size") >= 2]
    print("\n=== T17-1 slope d ~ AlphaMissense within cell (pooled) ===")
    obs, dist, p = perm_within(MS, rng)
    print(f"  beta {obs:+.4f} (perm {dist.mean():+.4f} +- {dist.std():.4f}; dz {(obs - dist.mean()) / dist.std():+.1f}) p {p:.4f}  "
          f"(n {len(MS)}, cells {MS.cell.nunique()})")
    rows = [dict(test="T17-1 pooled", n=len(MS), beta=obs, p=p)]
    for coh, x in MS.groupby("coh"):
        o, d_, pp = perm_within(x, rng, B=500)
        print(f"  {coh:5s} beta {o:+.4f} (dz {(o - d_.mean()) / d_.std():+.1f}) p {pp:.4f} (n {len(x)})")
        rows.append(dict(test=f"T17-1 {coh}", n=len(x), beta=o, p=pp))
    x5 = MS[MS.groupby("cell").am.transform("size") >= 5]
    o, d_, pp = perm_within(x5, rng, B=500)
    print(f"  only cells >= 5 variants: beta {o:+.4f} p {pp:.4f} (n {len(x5)})")

    print("\n=== T17-2 classes (mean d, cell-centred and raw) ===")
    MS["am_class"] = pd.cut(MS.am, [-0.01, 0.34, 0.564, 1.01], labels=CLASS_LABELS)
    dc = MS.d - MS.groupby("cell").d.transform("mean")
    for k in CLASS_LABELS:
        m = MS.am_class == k
        print(f"  {k:10s} n {m.sum():6d}: d raw {MS.d[m].mean():+.4f} | cell-centred {dc[m].mean():+.4f} | "
              f"fraction d>0 {(MS.d[m] > 0).mean():.3f}")
    PL = X[X.plof & X.trait.notna()]
    print(f"  pLoF       n {len(PL):6d}: d raw {PL.d.mean():+.4f} | fraction d>0 {(PL.d > 0).mean():.3f}")
    b2 = MS[MS.am_class.isin(["benign", "pathogen"])].copy()
    b2["am"] = (b2.am_class == "pathogen").astype(float)
    b2 = b2[b2.groupby("cell").am.transform("nunique") == 2]
    o, d_, pp = perm_within(b2, rng, B=500)
    print(f"  pathogen - benign within cell: {o:+.4f} p {pp:.4f} (n {len(b2)}, cells {b2.cell.nunique()})")
    rows.append(dict(test="T17-2 pathogen-benign", n=len(b2), beta=o, p=pp))

    print("\n=== T17-3 negative control: direction maps permuted between genes ===")
    cells = MS[["gene_id", "trait"]].drop_duplicates()
    genes = sorted(set(cells.gene_id))
    null = []
    for b in range(200):
        perm = dict(zip(genes, rng.permutation(genes)))
        s2 = np.array([directions.get(perm[g], {}).get(t, 0) for g, t in zip(MS.gene_id, MS.trait)])
        y = MS.assign(dn=s2 * MS.zo)
        y = y[s2 != 0]
        null.append(within_slope(y, "dn") if len(y) > 100 else np.nan)
    null = np.array(null)
    print(f"  slope with foreign directions: mean {np.nanmean(null):+.4f} +- {np.nanstd(null):.4f}; observed {rows[0]['beta']:+.4f}")
    print(f"  undirected |z| ~ am (within cell): {within_slope(MS.assign(az=MS.zo.abs()), 'az'):+.4f}")
    return rows, MS
