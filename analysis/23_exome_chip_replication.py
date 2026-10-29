"""Replication with freely available exome-chip consortia that predate UK Biobank exome sequencing (exploratory).

Sources: GIANT height (Marouli 2017, Height_All_add_SV, ~459k; exome-chip discovery across 147 studies before UK Biobank
exome sequencing) and GLGC lipids 2017 (LDL/HDL/TG, ~300k, before the UK Biobank biomarkers).
Variants are assigned to genes through coding exons (UCSC refGene hg19, NM_ transcripts, +-2 bp); only rare variants
(MAF < 1 %) are used. Exome-chip variants are mostly missense (no consequence annotation available), so they are compared
with the Genebass missense|LC burden.
Gene statistic: SKAT-like without LD (rare variants approximately independent): sum z^2 ~ chi2(k) -> p_gene ->
z2_gene = (Phi^-1(p/2))^2, capped at 50; directional (burden-like): Stouffer sum of the z_i.
Tests: (R1) height, AR_nos versus matched controls (LOEUF x CDS length) and versus AR_nonskel;
(R2) gene-level concordance with the Genebass missense burden for height;
(R3) organ diagonal for lipids: recessive genes of the lipid organ versus other recessive genes for LDL/HDL/TG.

Inputs: config.EXOME_CHIP_DIR (giant_height_exomchip.txt.gz, glgc_exome_{LDL,HDL,TG}.txt),
config.REF_DIR/refGene_hg19.txt.gz, gene_annotation.tsv.gz, go_classes.tsv, organ_fractions_recessive.tsv.gz,
config.GENEBASS_DIR/genebass_height_mislc.tsv.gz.
Outputs: 23_exome_chip_gene_stats_{height,ldl,hdl,tg}.tsv.
Call: python 23_exome_chip_replication.py
"""
import gzip
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2, spearmanr, fisher_exact

import config

rng = np.random.default_rng(17)

# --- coding-sequence intervals per chromosome (refGene hg19, NM_ transcripts, exons clipped to the CDS, +-2 bp)
intervals = defaultdict(list)
with gzip.open(config.REF_DIR / "refGene_hg19.txt.gz", "rt") as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if not c[1].startswith("NM_"):
            continue
        ch, cds_start, cds_end, sym = c[2].replace("chr", ""), int(c[6]), int(c[7]), c[12]
        for s, e in zip(map(int, c[9].strip(",").split(",")), map(int, c[10].strip(",").split(","))):
            s2, e2 = max(s, cds_start), min(e, cds_end)
            if s2 < e2:
                intervals[ch].append((s2 - 2, e2 + 2, sym))
index = {ch: (sorted(v), [x[0] for x in sorted(v)]) for ch, v in intervals.items()}


def genes_at(ch, pos):
    if ch not in index:
        return set()
    v, starts = index[ch]
    i = bisect_right(starts, pos)
    out = set()
    for j in range(max(0, i - 200), i):
        s, e, g = v[j]
        if s < pos <= e:
            out.add(g)
    return out


def gene_stats(df):
    """Per gene: number of variants k, z2 (chi2(k) combination, capped at 50) and directional Stouffer z."""
    df = df[np.isfinite(df.z)]
    per_gene = defaultdict(list)
    for ch, pos, z in zip(df.ch.astype(str), df.pos.astype(int), df.z):
        for g in genes_at(ch, pos):
            per_gene[g].append(z)
    out = []
    for g, zs in per_gene.items():
        zs = np.array(zs)
        k = len(zs)
        p = chi2.sf((zs ** 2).sum(), k)
        out.append((g, k, min(norm.isf(max(p, 1e-300) / 2) ** 2, 50), zs.sum() / np.sqrt(k)))
    return pd.DataFrame(out, columns=["symbol", "k", "z2", "zdir"]).set_index("symbol")


trait_stats = {}
g = pd.read_csv(config.EXOME_CHIP_DIR / "giant_height_exomchip.txt.gz", sep="\t", na_values=["-"])
g["maf"] = g.ExAC_MAF.astype(str).str.split(":").str[-1].astype(float, errors="ignore")
g["maf"] = pd.to_numeric(g.maf, errors="coerce")
g = g[(g.maf < 0.01) & g.se.gt(0)]
trait_stats["height"] = gene_stats(pd.DataFrame({"ch": g.CHR, "pos": g.POS, "z": g.beta / g.se}))
for t in ["LDL", "HDL", "TG"]:
    d = pd.read_csv(config.EXOME_CHIP_DIR / f"glgc_exome_{t}.txt", sep="\t")
    d[["ch", "pos"]] = d.SNP_hg19.str.split(":", expand=True)
    d = d[d.ch.str.match(r"^[0-9XY]+$")]
    maf = np.minimum(d["Freq.A1.ESP.EUR"], 1 - d["Freq.A1.ESP.EUR"])
    d = d[(maf < 0.01) & d.SE.gt(0)]
    trait_stats[t.lower()] = gene_stats(pd.DataFrame({"ch": d.ch, "pos": d.pos, "z": d.Beta / d.SE}))
for k, v in trait_stats.items():
    print(f"{k}: {len(v)} genes with rare chip variants, median {v.k.median():.0f} variants/gene")
    v.to_csv(config.results(f"23_exome_chip_gene_stats_{k}.tsv"), sep="\t")

A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding & A.oe_lof_upper.notna() & A.cds_length.notna()].copy()
A["stratum"] = pd.qcut(A.oe_lof_upper, 10, labels=False) * 10 + pd.qcut(A.cds_length, 5, labels=False)


def matched(T, C, col, n_boot=5000):
    """Mean of col in T versus control genes C drawn with replacement within the same LOEUF x length stratum."""
    by_stratum = {k: v[col].to_numpy() for k, v in C.groupby("stratum")}
    tot = np.zeros(n_boot)
    for k, n in T.stratum.value_counts().items():
        arr = by_stratum[k] if k in by_stratum else by_stratum[min(by_stratum, key=lambda x: abs(x - k))]
        tot += arr[rng.integers(0, len(arr), size=(n_boot, n))].sum(axis=1)
    tot /= len(T)
    o = T[col].mean()
    return o, tot.mean(), (o - tot.mean()) / tot.std(), (1 + (tot >= o).sum()) / (n_boot + 1)


print("\n(R1) Height, independent exome-chip cohort (GIANT):")
H = A.merge(trait_stats["height"], left_on="symbol", right_index=True)
H["sig"] = (H.z2 > norm.isf(0.025) ** 2).astype(float)
for t in ["AR_nos", "AD_nos", "nos_other", "AR_nonskel", "AR_skelsyndr"]:
    T = H[H.set == t]
    for c in ["pool", "AR_nonskel"]:
        if c == t:
            continue
        o, m, dz, p = matched(T, H[H.set == c], "z2")
        o2, m2, dz2, p2 = matched(T, H[H.set == c], "sig")
        print(f"  {t:12s} n {len(T):3d} vs {c:10s}: z2 {o:.2f} vs {m:.2f} (dz {dz:+.1f}, p {p:.4f}); "
              f"fraction p<0.05 {o2:.3f} vs {m2:.3f} (p {p2:.4f})")
# matrix genes
G = pd.read_csv(config.results("go_classes.tsv"), sep="\t")[["gene_id", "category"]]
H = H.merge(G, on="gene_id", how="left")
x = H[H.set == "AR_nos"]
print("  AR_nos by GO category (fraction p<0.05, n):",
      {k: (round(v.sig.mean(), 2), len(v)) for k, v in x.groupby("category")})
a, b = x[x.category == "MATRIX"], x[x.category.isin(["CILIUM", "LYSO"])]
p_fisher = fisher_exact([[a.sig.sum(), len(a) - a.sig.sum()], [b.sig.sum(), len(b) - b.sig.sum()]], alternative="greater")[1]
print(f"  MATRIX vs cilium+lysosome: {a.sig.mean():.2f} vs {b.sig.mean():.2f}, Fisher p {p_fisher:.3f}")

print("\n(R2) Concordance with Genebass missense (height):")
gb = pd.read_csv(config.GENEBASS_DIR / "genebass_height_mislc.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
gb["zgb"] = np.sign(gb.BETA_Burden) * norm.isf(gb.Pvalue_Burden.clip(1e-300) / 2)
y = H.merge(gb[["gene_id", "zgb"]], on="gene_id")
for t in ["AR_nos", "pool"]:
    yy = y[y.set == t]
    r, p = spearmanr(yy.zgb ** 2, yy.z2)
    r2, p2 = spearmanr(yy.zgb, yy.zdir)
    print(f"  {t}: z2 rho {r:+.3f} (p {p:.1e}), directional rho {r2:+.3f} (p {p2:.1e}), n {len(yy)}")
for s in ["SERPINH1", "COL27A1", "DDR2", "ADAMTS17", "NPR3", "PRKG2", "XYLT1", "SEC24D", "CHSY1", "SLC35D1", "CANT1",
          "ACAN", "NPR2", "IHH"]:
    r = y[y.symbol == s]
    if len(r):
        r = r.iloc[0]
        print(f"    {s:9s} chip k {int(r.k):3d} z2 {r.z2:6.1f} zdir {r.zdir:+.2f} | Genebass missense z {r.zgb:+.2f}")

print("\n(R3) Organ diagonal for lipids (GLGC):")
OF = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
members = OF[OF.frac >= 0.2]
lipid_genes = set(members[members.organ == "lipids"].gene_id)
all_recessive = set(members.gene_id)
for t in ["ldl", "hdl", "tg"]:
    L = A.merge(trait_stats[t], left_on="symbol", right_index=True)
    L = L[L.gene_id.isin(all_recessive)]
    a, b = L[L.gene_id.isin(lipid_genes)], L[~L.gene_id.isin(lipid_genes)]
    sig_a, sig_b = (a.z2 > 3.84).mean(), (b.z2 > 3.84).mean()
    perm = np.array([L.z2.sample(len(a), random_state=i).mean() for i in range(5000)])
    print(f"  {t}: lipid recessive genes n {len(a)} z2 {a.z2.mean():.2f} vs other recessive n {len(b)} {b.z2.mean():.2f}; "
          f"permutation p {(1 + (perm >= a.z2.mean()).sum()) / 5001:.4f}; "
          f"p<0.05 {sig_a:.2f} vs {sig_b:.2f}; top: {', '.join(a.sort_values('z2', ascending=False).symbol.head(5))}")
