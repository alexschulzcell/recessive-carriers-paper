"""Replication tests in cohorts without UK Biobank participants (shared by scripts 25, 35, 42, 43, 45 and 49).

Variant annotation
    pLoF   = gnomAD v2 LoF variant with LOFTEE "HC" (GeneBayes Zenodo file, GRCh37 positions plus pos_hg38), matched by
             position and alleles (both orientations).
    coding = variant in a coding exon +-2 bp (UCSC refGene NM_ transcripts, hg19 or hg38); mostly missense, synonymous
             variants are included.
    Only rare variants (MAF < 1 %) are used.

Gene statistics
    gene_stats_z2      undirected: sum of z^2 over the gene's variants ~ chi2(k), converted to a z^2 (capped at 50).
    signed_gene_stats  directed (Stouffer): sum z / sqrt(k), |z| <= sqrt(50), z oriented to the LoF allele for pLoF.

Tests (paper test IDs in brackets)
    t1_height_enrichment      [T1] height z2 of recessive skeletal-dysplasia genes (AR_nos) against matched controls.
    t2_organ_diagonal         [T2] organ diagonal: doubly centred z2 in the traits of each gene's disease organ,
                                   permutation of organ labels between genes.
    t3_producer_rule          [T3] pre-specified producer rule W1/W2b (classes from script 22) within organs.
    directional_tests         [T10-8b] gene-specific expected direction (direction_map.py), Stouffer gene z.
    frozen_prediction_test    [P1 / F5 / B5] frozen atlas prediction against an independent cohort.
    gene_scores, cohort_consistency  [P2 / F6 / B6] gene-level consistency of the directed carrier effect.
    skeletal_mechanism_test   [P4 / F7 / B7] matrix genes vs cilium + lysosome genes for height.
    dominance_h_height        [P5 / F8 / B8] dominance coefficient h for height.

Every test takes the random generator as an argument; the draws are made in the same order as in the original
analysis, so the permutation p-values are reproducible.
"""
import ast
import gzip
import re
import tarfile
from bisect import bisect_right
from collections import defaultdict
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2, combine_pvalues, fisher_exact, norm, spearmanr, wilcoxon

import common
import config

CAP = np.sqrt(50)


# ---------------------------------------------------------------- variant annotation
def cds_index(fn):
    """Coding intervals (+-2 bp) of NM_ transcripts per chromosome: {chrom: (sorted intervals, start positions)}."""
    iv = defaultdict(list)
    with gzip.open(fn, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if not c[1].startswith("NM_") or "_" in c[2]:
                continue
            ch, cs, ce, sym = c[2].replace("chr", ""), int(c[6]), int(c[7]), c[12]
            for s, e in zip(map(int, c[9].strip(",").split(",")), map(int, c[10].strip(",").split(","))):
                s2, e2 = max(s, cs), min(e, ce)
                if s2 < e2:
                    iv[ch].append((s2 - 2, e2 + 2, sym))
    return {ch: (sorted(v), [x[0] for x in sorted(v)]) for ch, v in iv.items()}


def _cds_genes(idx, ch, pos):
    if ch not in idx:
        return ()
    v, st = idx[ch]
    i = bisect_right(st, pos)
    return {g for s, e, g in v[max(0, i - 60):i] if s < pos <= e}


class VariantAnnotation:
    """pLoF variants (gnomAD v2 LOFTEE HC, GRCh37 and GRCh38) and coding-exon lookup (refGene hg19/hg38)."""

    def __init__(self, gene_symbols=True):
        # symbol -> gene_id from the gene annotation of script 03 (not needed by the download scripts 49a/49b)
        if gene_symbols:
            annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
            self.sym2id = annotation.drop_duplicates("symbol").set_index("symbol").gene_id.to_dict()
        self.lof19, self.lof38 = {}, {}
        lofs = pd.read_csv(config.REPLICATION_DIR / "annot" / "gnomad_v2_lofs.tsv.gz", sep="\t",
                           usecols=["ensg", "locus", "ref", "alt", "lof", "pos_hg38"])
        for r in lofs.itertuples():
            try:
                alts, lof_flags = ast.literal_eval(r.alt), ast.literal_eval(r.lof)
            except Exception:  # noqa: BLE001  (malformed list fields are skipped)
                continue
            ch, pos = r.locus.split(":")
            for a, flag in zip(alts, lof_flags):
                if flag == "HC":
                    self.lof19[(ch, int(pos), r.ref, a)] = r.ensg
                    if pd.notna(r.pos_hg38):
                        self.lof38[(ch, int(r.pos_hg38), r.ref, a)] = r.ensg
        print("gnomAD HC LoF variants:", len(self.lof19), "hg19,", len(self.lof38), "hg38", flush=True)
        self.cds = {"hg19": cds_index(config.REF_DIR / "refGene_hg19.txt.gz"),
                    "hg38": cds_index(config.REF_DIR / "refGene_hg38.txt.gz")}
        self.coding_genes = lru_cache(maxsize=None)(self._coding_genes)

    def _coding_genes(self, build, ch, pos):
        """Gene symbols whose coding exon (+-2 bp) contains the position."""
        return frozenset(_cds_genes(self.cds[build], ch, pos))

    def lof_dict(self, build):
        return self.lof19 if build == "hg19" else self.lof38


def gene_stats_z2(df, build, ann):
    """df: ch, pos, ref, alt, z, maf. Returns {'plof': Series z2, 'coding': Series z2} indexed by gene_id."""
    df = df[np.isfinite(df.z) & (df.maf < 0.01)]
    lofd = ann.lof_dict(build)
    pl, kd = defaultdict(list), defaultdict(list)
    for ch, pos, ref, alt, z in zip(df.ch.astype(str), df.pos.astype(int), df.ref.astype(str), df.alt.astype(str), df.z):
        g = lofd.get((ch, pos, ref, alt)) or lofd.get((ch, pos, alt, ref))
        if g:
            pl[g].append(z)
        for s in ann.coding_genes(build, ch, pos):
            gid = ann.sym2id.get(s)
            if gid:
                kd[gid].append(z)

    def agg(d):
        return pd.Series({g: min(norm.isf(max(chi2.sf(np.sum(np.square(v)), len(v)), 1e-300) / 2) ** 2, 50) for g, v in d.items()})
    return {"plof": agg(pl), "coding": agg(kd)}


def signed_gene_stats(df, build, ann):
    """Stouffer gene z (sum z / sqrt(k)); pLoF z oriented to the LoF allele. Returns {'plof': Series, 'coding': Series}."""
    df = df[np.isfinite(df.z) & (df.maf < 0.01)]
    lofd = ann.lof_dict(build)
    pl, kd = defaultdict(list), defaultdict(list)
    for ch, pos, ref, alt, z in zip(df.ch.astype(str), df.pos.astype(int), df.ref.astype(str), df.alt.astype(str), df.z):
        z = float(np.clip(z, -CAP, CAP))
        g = lofd.get((ch, pos, ref, alt))
        if g is None and lofd.get((ch, pos, alt, ref)):
            g = lofd[(ch, pos, alt, ref)]
            z = -z  # the file's alternative allele is the reference allele -> orient to the LoF allele
        if g:
            pl[g].append(z)
        for s in ann.coding_genes(build, ch, pos):
            gid = ann.sym2id.get(s)
            if gid:
                kd[gid].append(z)
    agg = lambda d: pd.Series({g: np.clip(np.sum(v) / np.sqrt(len(v)), -CAP, CAP) for g, v in d.items()})  # noqa: E731
    return {"plof": agg(pl), "coding": agg(kd)}


def load_gene_stats(sources, ann, cache_dir, signed=False):
    """Gene statistics per (cohort, trait), computed once and cached as TSV in cache_dir.

    Returns {(cohort, trait): {'plof': Series, 'coding': Series}}.
    """
    out = {}
    for coh, items in sources.items():
        for k, build, loader in items:
            f = cache_dir / f"{coh}_{k}.tsv.gz"
            if f.exists():
                g = pd.read_csv(f, sep="\t", index_col=0)
                out[(coh, k)] = {c: g[c].dropna() for c in g.columns}
                continue
            gs = (signed_gene_stats if signed else gene_stats_z2)(loader(), build, ann)
            pd.DataFrame(gs).to_csv(f, sep="\t")
            out[(coh, k)] = gs
            print(f"{coh} {k}: pLoF genes {len(gs['plof'])}, coding genes {len(gs['coding'])}", flush=True)
    return out


# ---------------------------------------------------------------- readers (one DataFrame: ch, pos, ref, alt, z, maf)
def read_genes_and_health(trait):
    d = pd.read_csv(config.REPLICATION_DIR / "gh" / f"{trait}.h.tsv.gz", sep="\t",
                    usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                             "effect_allele_frequency"], dtype={"chromosome": str})
    f = d.effect_allele_frequency
    return pd.DataFrame({"ch": d.chromosome, "pos": d.base_pair_location, "ref": d.other_allele, "alt": d.effect_allele,
                         "z": d.beta / d.standard_error, "maf": np.minimum(f, 1 - f)})


def read_finmetseq(name):
    d = pd.read_csv(config.REPLICATION_DIR / "finmetseq" / f"{name}.gz", sep="\t",
                    usecols=["CHROM", "BEG", "MARKER_ID", "MAF", "BETA", "SEBETA"], dtype={"CHROM": str})
    ra = d.MARKER_ID.str.split("_").str[1].str.split("/")
    return pd.DataFrame({"ch": d.CHROM, "pos": d.BEG, "ref": ra.str[0], "alt": ra.str[1], "z": d.BETA / d.SEBETA, "maf": d.MAF})


def read_exome_chip(path_or_obj):
    d = pd.read_csv(path_or_obj, sep="\t", na_values=["-"], dtype={"CHR": str, "CHROM": str}, low_memory=False).rename(columns={"CHROM": "CHR"})
    maf = pd.to_numeric(d.ExAC_MAF.astype(str).str.split(":").str[-1], errors="coerce")
    return pd.DataFrame({"ch": d.CHR, "pos": d.POS, "ref": d.REF, "alt": d.ALT, "z": d.beta / d.se, "maf": maf})


def read_blood_cell_consortium(tarname, member):
    with tarfile.open(config.REPLICATION_DIR / "bcx" / tarname) as t:
        return read_exome_chip(t.extractfile(member))


def read_ckdgen(name):
    d = pd.read_csv(config.REPLICATION_DIR / "ckdgen" / name, dtype={"chr": str})
    return pd.DataFrame({"ch": d.chr, "pos": d.bp_hg19, "ref": d.allele2, "alt": d.allele1,
                         "z": pd.to_numeric(d.beta, errors="coerce") / pd.to_numeric(d.se, errors="coerce"),
                         "maf": pd.to_numeric(d.ExAC_EUR_MAF, errors="coerce")})


def read_glgc(trait):
    d = pd.read_csv(config.EXOME_CHIP_DIR / f"glgc_exome_{trait}.txt", sep="\t")
    d[["ch", "pos"]] = d.SNP_hg19.str.split(":", expand=True)
    d = d[d.ch.str.match(r"^[0-9XY]+$")]
    f = d["Freq.A1.ESP.EUR"]
    return pd.DataFrame({"ch": d.ch, "pos": d.pos.astype(int), "ref": d.A2, "alt": d.A1, "z": d.Beta / d.SE, "maf": np.minimum(f, 1 - f)})


def independent_cohort_sources():
    """{cohort: [(trait, genome build, loader)]} for Genes & Health (GH), FinMetSeq (FMS) and the exome-chip consortia (CHIP)."""
    return {
        "GH": [(k, "hg38", lambda k=k: read_genes_and_health(k))
               for k in ["height", "hb", "mcv", "wbc", "plt", "creat", "alt", "bili", "ldl", "hdl", "tg", "hba1c"]],
        "FMS": [(k, "hg19", lambda f=f: read_finmetseq(f))
                for k, f in [("height", "Height"), ("creat", "Creatinine"), ("egfr", "eGFR"), ("ldl", "LDL-C"), ("hdl", "HDL-C"), ("tg", "Total_TG")]],
        "CHIP": [("height", "hg19", lambda: read_exome_chip(config.EXOME_CHIP_DIR / "giant_height_exomchip.txt.gz"))]
                + [(k, "hg19", lambda m=m: read_blood_cell_consortium("BCX_RBC_ALL.tar.gz", m))
                   for k, m in [("hb", "BCX_exomechip_HGB_ALL_20160420.txt"), ("mcv", "BCX_exomechip_MCV_ALL_20160420.txt"),
                                ("rbc", "BCX_exomechip_RBC_ALL_20160420.txt")]]
                + [("plt", "hg19", lambda: read_blood_cell_consortium("BCX_PLT_ALL.tar.gz", "BCX_exomechip_PLT_ALL_20160420.txt")),
                   ("wbc", "hg19", lambda: read_blood_cell_consortium("BCX_WBC_ALL.tar.gz", "BCX_exomechip_WBC_ALL_20160420.txt")),
                   ("lymph", "hg19", lambda: read_blood_cell_consortium("BCX_WBC_ALL.tar.gz", "BCX_exomechip_LYM_ALL_20160420.txt")),
                   ("egfr", "hg19", lambda: read_ckdgen("Published_eGFRcrea_overall_EA.csv.gz")),
                   ("egfrcys", "hg19", lambda: read_ckdgen("Published_eGFRcys_overall_EA.csv.gz"))]
                + [(t.lower(), "hg19", lambda t=t: read_glgc(t)) for t in ["LDL", "HDL", "TG"]],
    }


def read_finngen_bbj(cohort_dir, trait):
    """Pre-filtered FinnGen / BioBank Japan file (scripts 49a/49b): adds z, maf, and keeps beta and se."""
    d = pd.read_csv(cohort_dir / f"{trait}.tsv.gz", sep="\t", dtype={"ch": str})
    d = d[d.sebeta > 0]
    return pd.DataFrame({"ch": d.ch, "pos": d.pos, "ref": d.ref, "alt": d.alt, "z": d.beta / d.sebeta,
                         "maf": np.minimum(d.af_alt, 1 - d.af_alt), "beta": d.beta, "se": d.sebeta})


# ---------------------------------------------------------------- T1-T3 (undirected z2)
# Trait -> organ of the replication traits (a priori, as in script 14). eGFR is a kidney trait.
TRAIT_TO_ORGAN = {"height": "skeletal", "hb": "blood", "mcv": "blood", "rbc": "blood", "plt": "blood", "wbc": "immune",
                  "lymph": "immune", "creat": "kidney", "egfr": "kidney", "egfrcys": "kidney", "alt": "liver", "bili": "liver",
                  "ldl": "lipids", "hdl": "lipids", "tg": "lipids", "hba1c": "glucose"}


def matched_strata_annotation(annotation):
    """Protein-coding genes with LOEUF and CDS length; stratum = LOEUF decile x CDS-length quintile."""
    a = annotation[annotation.protein_coding & annotation.oe_lof_upper.notna() & annotation.cds_length.notna()].copy()
    a["stratum"] = pd.qcut(a.oe_lof_upper, 10, labels=False) * 10 + pd.qcut(a.cds_length, 5, labels=False)
    return a


def matched_test(T, C, col, rng, n_draws=5000):
    """Mean of col in T against stratum-matched random draws (with replacement) from C. Returns mean, dz, one-sided p."""
    by_stratum = {k: v[col].to_numpy() for k, v in C.groupby("stratum")}
    tot = np.zeros(n_draws)
    for k, n in T.stratum.value_counts().items():
        arr = by_stratum[k] if k in by_stratum else by_stratum[min(by_stratum, key=lambda x: abs(x - k))]
        tot += arr[rng.integers(0, len(arr), size=(n_draws, n))].sum(axis=1)
    tot /= len(T)
    o = T[col].mean()
    return o, (o - tot.mean()) / tot.std(), (1 + (tot >= o).sum()) / (n_draws + 1)


def t1_height_enrichment(stats, cohorts, strata, rng, rows):
    print("\nT1 height: AR_nos against matched controls")
    for coh in cohorts:
        for vs in ["plof", "coding"]:
            s = stats[(coh, "height")][vs]
            H = strata.merge(s.rename("z2"), left_on="gene_id", right_index=True)
            T = H[H.set == "AR_nos"]
            if len(T) < 20:
                print(f"  {coh} {vs}: too few AR genes ({len(T)})")
                continue
            for c in ["pool", "AR_nonskel"]:
                o, dz, p = matched_test(T, H[H.set == c], "z2", rng)
                print(f"  {coh:4s} {vs:6s} AR_nos n {len(T):3d} vs {c:10s}: z2 {o:.2f}, dz {dz:+.1f}, p {p:.4f}")
                rows.append(dict(test="T1", cohort=coh, set=vs, comparison=c, n=len(T), dz=dz, p=p))
            T2 = H[H.set == "AR_nonskel"]
            o, dz, p = matched_test(T2, H[H.set == "pool"], "z2", rng)
            print(f"  {coh:4s} {vs:6s} control AR_nonskel vs pool: dz {dz:+.1f}, p {p:.3f}")


def centred_z_matrix(stats, coh, vs, genes):
    """Genes x traits z2 of one cohort, traits with >= 50 genes, centred per trait and then per gene."""
    Z = pd.DataFrame({k: stats[(c, k)][vs] for (c, k) in stats if c == coh}).reindex(genes)
    Z = Z.loc[:, Z.notna().sum() >= 50]
    Zc = Z - Z.mean(axis=0)
    return Zc.sub(Zc.mean(axis=1), axis=0)


def t2_organ_diagonal(stats, cohorts, organ_members, trait_to_organ, rng, rows):
    print("\nT2 organ diagonal (permutation of organ labels, B=2000)")
    for coh in cohorts:
        for vs in ["plof", "coding"]:
            genes = sorted(set(organ_members.gene_id))
            Zc = centred_z_matrix(stats, coh, vs, genes)
            Zc = Zc[Zc.notna().sum(axis=1) >= 2]
            organs = common.sort_organs({trait_to_organ[t] for t in Zc.columns})
            g = list(Zc.index)
            Lm = np.array([[1.0 if x in set(organ_members[organ_members.organ == o].gene_id) else 0.0 for o in organs] for x in g])
            Cm = np.array([[1.0 if trait_to_organ[t] == o else 0.0 for t in Zc.columns] for o in organs])
            V = Zc.values
            valid = ~np.isnan(V)
            V0 = np.nan_to_num(V)
            stat = lambda Lx: (V0 * (((Lx @ Cm) > 0) & valid)).sum() / (((Lx @ Cm) > 0) & valid).sum()  # noqa: E731
            obs = stat(Lm)
            perm = np.array([stat(Lm[rng.permutation(len(Lm))]) for _ in range(2000)])
            dz, p = (obs - perm.mean()) / perm.std(), (1 + (perm >= obs).sum()) / 2001
            per_organ = []
            for o in organs:
                Lo = Lm[:, [organs.index(o)]]
                Co = Cm[[organs.index(o)], :]
                st = lambda Lx: (V0 * (((Lx @ Co) > 0) & valid)).sum() / max((((Lx @ Co) > 0) & valid).sum(), 1)  # noqa: E731
                ob = st(Lo)
                pe = np.array([st(Lo[rng.permutation(len(Lo))]) for _ in range(1000)])
                per_organ.append(f"{o} {(ob - pe.mean()) / pe.std():+.1f}")
            print(f"  {coh:4s} {vs:6s}: {len(g)} genes, traits {list(Zc.columns)} -> dz {dz:+.1f}, p {p:.4f} | per organ: {', '.join(per_organ)}")
            rows.append(dict(test="T2", cohort=coh, set=vs, comparison="organ_permutation", n=len(g), dz=dz, p=p))


def load_producer_units():
    """Gene x organ units of the producer rule (script 22, pLoF) without the skeletal organ."""
    U = pd.read_csv(config.results("22_main_product_units_plof.tsv"), sep="\t")[["gene_id", "organ", "W1", "W2b"]]
    return U[U.organ != "skeletal"]


def t3_producer_rule(stats, cohorts, U, trait_to_organ, rng, rows):
    print("\nT3 rule (W1, W2b; classes from script 22; without skeletal; permutation within organ, B=5000)")
    for coh in cohorts:
        for vs in ["plof", "coding"]:
            Zc = centred_z_matrix(stats, coh, vs, sorted(set(U.gene_id)))
            organ_traits = defaultdict(list)
            for t in Zc.columns:
                organ_traits[trait_to_organ[t]].append(t)
            vals = [Zc.loc[g, organ_traits[o]].mean() if o in organ_traits and g in Zc.index else np.nan for g, o in zip(U.gene_id, U.organ)]
            X = U.assign(value=vals).dropna(subset=["value"])
            for route in ["W1", "W2b"]:
                lab = X[route].values.astype(bool)
                w = X.value.values
                if lab.sum() < 5:
                    continue
                idx = [np.where(X.organ.values == o)[0] for o in X.organ.unique()]
                obs = w[lab].mean() - w[~lab].mean()
                d = np.empty(5000)
                for i in range(5000):
                    Lx = lab.copy()
                    for ix in idx:
                        Lx[ix] = rng.permutation(Lx[ix])
                    d[i] = w[Lx].mean() - w[~Lx].mean()
                p = (1 + (d >= obs).sum()) / 5001
                print(f"  {coh:4s} {vs:6s} {route}: diff {obs:+.3f} (n yes {int(lab.sum())}/{len(X)}; organs "
                      f"{common.sort_organs(X.organ.unique())}), p {p:.4f}")
                rows.append(dict(test=f"T3_{route}", cohort=coh, set=vs, comparison="within_organ", n=len(X), dz=obs, p=p))


# ---------------------------------------------------------------- T10-8b (directed, gene-specific direction)
FLIP = {"egfr": ("creat", -1), "egfrcys": ("cysc", -1)}   # sensitivity: eGFR as inverse creatinine / cystatin C


def directional_matrices(signed, dm, flip=False):
    """Raw and centred genes x traits matrices of signed gene z for one cohort x set ({trait: Series})."""
    all_recessive = sorted(dm.recessive_genes)
    d = {}
    for k, v in signed.items():
        if k in dm.traits:
            d[k] = v
        elif flip and k in FLIP and FLIP[k][0] not in signed:
            d[FLIP[k][0]] = FLIP[k][1] * v
    Z = pd.DataFrame({t: v.reindex(all_recessive) for t, v in d.items()}, index=all_recessive).reindex(columns=dm.traits)
    Zc = Z - Z.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    return Z.loc[dm.genes].values, Zc.loc[dm.genes].values


def _directional_stat(Z, Sm):
    m = (Sm != 0) & np.isfinite(Z)
    return np.nansum(np.where(m, Sm * Z, 0)) / max(m.sum(), 1), int(m.sum())


def directional_test(Zc, dm, name, rng, n_perm=1000):
    """Mean of expected sign x centred z over direction cells; permutation of the direction rows between genes."""
    obs, n = _directional_stat(Zc, dm.S)
    if n < 20:
        print(f"  {name:28s} too few cells (n {n})")
        return None
    pm = np.array([_directional_stat(Zc, dm.S[rng.permutation(len(dm.genes))])[0] for _ in range(n_perm)])
    dz = (obs - pm.mean()) / pm.std()
    p = (1 + (pm >= obs).sum()) / (n_perm + 1)
    print(f"  {name:28s} d {obs:+.3f} (perm {pm.mean():+.3f} +- {pm.std():.3f}) dz {dz:+5.1f}  p {p:.4f}  n {n}", flush=True)
    return dict(test=name, d=obs, dz=dz, p=p, n=n)


def directional_fraction(Z, dm, name, rng, lo=2, n_perm=1000):
    """Fraction of cells with |z| > lo whose sign matches the expectation; permutation as directional_test."""
    base = np.isfinite(Z) & (np.abs(Z) >= lo)

    def frac(Sm):
        k = base & (Sm != 0)
        return (((Sm * Z)[k] > 0).mean() if k.sum() else np.nan), int(k.sum())
    obs, n = frac(dm.S)
    if n < 20:
        return
    pm = np.array([frac(dm.S[rng.permutation(len(dm.genes))])[0] for _ in range(n_perm)])
    p = (1 + (pm >= obs).sum()) / (n_perm + 1)
    print(f"  {name:28s} fraction |z|>{lo}: {obs:.3f} (n {n}; perm {pm.mean():.3f}) p {p:.4f}", flush=True)


def directional_tests(signed_by_cohort, cohorts, dm, rng):
    """T10-8b for each cohort and set, then the eGFR sensitivity. signed_by_cohort: {(cohort, set): {trait: Series}}."""
    rows = []
    print("\n=== T10-8b: directed replication, variant based (Stouffer), gene-specific HPO direction ===")
    for coh in cohorts:
        for s in ("plof", "coding"):
            Z, Zc = directional_matrices(signed_by_cohort[(coh, s)], dm)
            r = directional_test(Zc, dm, f"{coh} {s}", rng)
            if r:
                rows.append(dict(r, cohort=coh, set=s, analysis="primary"))
            directional_fraction(Z, dm, f"{coh} {s}", rng)
    print("\n--- Sensitivity: eGFR -> creat/cysc (sign reversed) where creatinine is missing ---")
    for coh in cohorts:
        for s in ("plof", "coding"):
            Z, Zc = directional_matrices(signed_by_cohort[(coh, s)], dm, flip=True)
            r = directional_test(Zc, dm, f"{coh} {s} +eGFR", rng)
            if r:
                rows.append(dict(r, cohort=coh, set=s, analysis="with_eGFR"))
    R = pd.DataFrame(rows)
    for s in ("plof", "coding"):
        x = R[(R.analysis == "primary") & (R.set == s)]
        if len(x) > 1:
            print(f"\nFisher combination {s}: p {combine_pvalues(x.p, method='fisher')[1]:.2g} ({len(x)} sources)")
    return R


def signed_by_cohort_and_set(signed_stats):
    """{(cohort, trait): {set: Series}} -> {(cohort, set): {trait: Series}}."""
    out = {}
    for (coh, k), gs in signed_stats.items():
        for s in ("plof", "coding"):
            out.setdefault((coh, s), {})[k] = gs[s]
    return out


# ---------------------------------------------------------------- P1-P5 (cross-cohort tests on directed gene z)
def frozen_prediction_test(V, z_by_set, rng, rows):
    """P1: frozen atlas prediction against an independent cohort.

    V: frozen prediction (common.load_frozen_prediction). z_by_set: {'plof': {trait: Series}, 'plofmis': {...}}.
    Statistic: mean of expected direction x z over predicted cells; permutation of the expected directions within trait
    (B = 2000); secondary: fraction of correct signs among |z| > 2 (B = 500).
    """
    for m in ["plof", "plofmis"]:
        Z = z_by_set[m]
        x = V.copy()
        x["z"] = [Z[t].get(g, np.nan) if t in Z else np.nan for g, t in zip(x.gene_id, x.trait)]
        x = x.dropna(subset=["z"])
        x["d"] = x.expected_direction * x.z
        obs = x.d.mean()
        grp = x.groupby("trait").indices
        s = x.expected_direction.values
        z = x.z.values
        null = np.empty(2000)
        for b in range(2000):
            s2 = s.copy()
            for ix in grp.values():
                s2[ix] = s[rng.permutation(ix)]
            null[b] = (s2 * z).mean()
        p = (1 + (null >= obs).sum()) / 2001
        k = np.abs(z) > 2
        fr = (x.d[k] > 0).mean()
        nf = []
        for b in range(500):
            s2 = s.copy()
            for ix in grp.values():
                s2[ix] = s[rng.permutation(ix)]
            nf.append(((s2 * z)[k] > 0).mean())
        pf = (1 + (np.array(nf) >= fr).sum()) / 501
        print(f"  {m:8s} cells {len(x)} (genes {x.gene_id.nunique()}): mean d {obs:+.4f} (perm {null.mean():+.4f} +- {null.std():.4f}; "
              f"dz {(obs - null.mean()) / null.std():+.1f}) p {p:.4f}"
              f" | |z|>2: fraction correct {fr:.3f} (n {k.sum()}, perm {np.mean(nf):.3f}) p {pf:.4f}")
        for tt, y in x.groupby("specificity_tertile"):
            print(f"      specificity tertile {tt}: d {y.d.mean():+.4f} (n {len(y)})")
        rows.append(dict(test=f"P1 {m}", n=len(x), stat=obs, p=p))


def gene_scores(directions, get):
    """P2 gene score: mean of expected sign x z over the gene's direction cells (cells without z are skipped)."""
    out = {}
    for g, r in directions.items():
        v = [s * get(g, t) for t, s in r.items()]
        v = [a for a in v if np.isfinite(a)]
        if v:
            out[g] = np.mean(v)
    return pd.Series(out)


def ukb_gene_scores(dm):
    """P2 gene scores in UK Biobank (pLoF, raw signed Genebass z)."""
    gi = {g: i for i, g in enumerate(dm.genes)}
    ti = {t: j for j, t in enumerate(dm.traits)}
    return gene_scores(dm.directions, lambda g, t: dm.raw["plof"][gi[g], ti[t]] if t in ti else np.nan)


GH_BURDEN_TRAITS = {"Height": "height", "Hb": "hb", "MCV": "mcv", "RBC": "rbc", "Platelets": "plt", "WBC": "wbc",
                    "Lymphocytes": "lymph", "creatinine": "creat", "Urate": "urate", "Urea": "urea", "Phosphate": "phosphate",
                    "ALT": "alt", "AST": "ast", "GGT": "ggt", "Bilirubin": "bili", "Albumin": "albumin", "LDL-C": "ldl",
                    "HDL-C": "hdl", "Triglycerides": "tg", "Total_cholesterol": "chol", "Glucose_non_fasting": "glucose",
                    "HbA1c": "hba1c", "FEV1": "fev1", "Systolic_BP": "sbp"}


def genes_and_health_burden_z():
    """Signed gene-burden z of Genes & Health (REGENIE, mask A = LOFTEE HC pLoF, AAF < 1 %, allele count >= 5), capped at sqrt(50)."""
    Zg = {}
    for f in sorted(config.GH_BURDEN_DIR.glob("*_genetests_*.regenie.gz")):
        mt = re.search(r"2024_05_08_(.+?)\.residual_GNH", f.name)
        if not mt or mt.group(1) not in GH_BURDEN_TRAITS:
            continue
        d = pd.read_csv(f, sep=" ", comment="#", usecols=["ID", "A1FREQ", "N", "TEST", "BETA", "SE"])
        d = d[d.TEST == "ADD"]
        d["gene_id"] = d.ID.str.extract(r"\((ENSG\d+)\)")[0]
        d["mask"] = d.ID.str.extract(r"\)\.(MASK_.+)$")[0]
        d = d[(d["mask"] == "MASK_A.0.01") & (d.A1FREQ * 2 * d.N >= 5) & d.SE.gt(0)].drop_duplicates("gene_id")
        Zg[GH_BURDEN_TRAITS[mt.group(1)]] = pd.Series(np.clip(d.BETA / d.SE, -np.sqrt(50), np.sqrt(50)).values, index=d.gene_id.values)
    return Zg


def cohort_consistency(pairs, rng, rows):
    """P2: Spearman correlation of gene scores between two cohorts; one-sided permutation p (B = 2000)."""
    for a, b, na, nb in pairs:
        j = pd.concat([a, b], axis=1, keys=["a", "b"]).dropna()
        r = spearmanr(j.a, j.b).correlation
        null = np.array([spearmanr(j.a.values, rng.permutation(j.b.values)).correlation for _ in range(2000)])
        p = (1 + (null >= r).sum()) / 2001
        print(f"  {na}-{nb}: rho {r:+.3f} (n {len(j)}), p {p:.4f}")
        rows.append(dict(test=f"P2 {na}-{nb}", n=len(j), stat=r, p=p))


def skeletal_mechanism_test(height_z_by_set, rows):
    """P4: share of |z| > 1.96 for height among MATRIX genes vs cilium + lysosome genes (recessive dysplasia genes, Fisher)."""
    K = pd.read_csv(config.results("go_classes.tsv"), sep="\t")
    annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
    sym = annotation.drop_duplicates("gene_id").set_index("gene_id").symbol
    K = K[K.gene_id.isin(set(annotation[annotation.set == "AR_nos"].gene_id))]
    for m, Zm in height_z_by_set.items():
        x = K.assign(z=K.gene_id.map(Zm)).dropna(subset=["z"])
        x["hit"] = np.abs(x.z) > 1.96
        a = x[x.MATRIX == True]  # noqa: E712
        b = x[(x.CILIUM == True) | (x.LYSO == True)]  # noqa: E712
        o_, p_ = fisher_exact([[a.hit.sum(), (~a.hit).sum()], [b.hit.sum(), (~b.hit).sum()]], alternative="greater")
        print(f"  {m}: MATRIX {a.hit.mean():.3f} (n {len(a)}, z2 {np.mean(a.z ** 2):.2f}) vs cilium+lysosome {b.hit.mean():.3f} "
              f"(n {len(b)}, z2 {np.mean(b.z ** 2):.2f}) | OR {o_:.2f}, p {p_:.4f}")
        print("     MATRIX hits:", ", ".join(f"{sym.get(g, g)} {z:+.1f}" for g, z in x[(x.MATRIX == True) & x.hit][["gene_id", "z"]].values))  # noqa: E712
        rows.append(dict(test=f"P4 {m}", n=len(a) + len(b), stat=o_, p=p_))


def dominance_h_height(beta_het, rows):
    """P5: h = carrier effect (SD) / patient SDS for height, genes from the literature table (scripts 32/32b)."""
    annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
    L = pd.read_csv(config.results("literature_patient_values.tsv"), sep="\t")
    L = L[(L.trait == "height") & L.sds_patients.notna()].copy()
    L["gene_id"] = L.symbol.map(dict(zip(annotation.symbol, annotation.gene_id)))
    L["beta_het"] = [beta_het.get(g, np.nan) for g in L.gene_id]
    L = L.dropna(subset=["beta_het"]).drop_duplicates("symbol")
    L["h"] = L.beta_het / L.sds_patients
    h = L.h.values
    print(f"  genes {len(L)}: median h {np.median(h):+.4f} (IQR {np.percentile(h, 25):+.4f} to {np.percentile(h, 75):+.4f}); "
          f"fraction h>0 {(h > 0).mean():.2f}; sign test h>0 p {binomtest(int((h > 0).sum()), len(h), 0.5, alternative='greater').pvalue:.4f}; "
          f"Wilcoxon h<0.5 p {wilcoxon(h - 0.5, alternative='less').pvalue:.2g}")
    print(f"  median carrier effect {np.median(L.beta_het) * 9.3:+.2f} cm (UK Biobank SD 9.3 cm) against patients "
          f"{np.median(L.sds_patients) * 9.3:+.1f} cm")
    uk = pd.read_csv(config.results("29b_dominance_h_main_test.tsv"), sep="\t")
    if "h" in uk.columns and "symbol" in uk.columns:
        j = L.merge(uk[uk.trait == "height"][["symbol", "h"]], on="symbol", suffixes=("_cohort", "_ukb"))
        print(f"  gene concordance h cohort ~ UK Biobank: rho {spearmanr(j.h_cohort, j.h_ukb).correlation:+.3f} (n {len(j)})")
    rows.append(dict(test="P5 h height", n=len(L), stat=float(np.median(h)), p=np.nan))
