"""Organ diagonal in heterozygous IMPC knockout mice.

Pre-specified test (PRAEREG_maus.md). For every gene with a 1:1 mouse ortholog, the IMPC release gives per zygosity and
MP top-level system whether any parameter of the system was significant. The organs of each recessive gene (from its
OMIM/HPO disease annotation, as in the human organ diagonal) are mapped to MP top-level systems (ORGAN_TO_MP). The test
asks whether heterozygous mice of a gene show phenotypes more often in the systems that match the organs of the human
disease (doubly centred hit matrix, permutation of the organ labels between genes). Homozygous mice serve as positive
control; a secondary analysis restricts to genes with both heterozygous and homozygous data and gives a rough mouse h.

Inputs: config.IMPC_RELEASE_DIR (statistical-results-ALL.csv.gz, HMD_HumanPhenotype.rpt), gene_annotation.tsv.gz,
organ_fractions_recessive.tsv.gz.
Output: 38_mouse_organ_diagonal.tsv.
Call: python 38_mouse_organ_diagonal.py
"""
import numpy as np
import pandas as pd

import config
import impc_mouse

# Pre-specified mapping of the organ systems to MP top-level systems (also read by scripts 38b and 47).
ORGAN_TO_MP = {"blood": ["hematopoietic system phenotype"], "immune": ["immune system phenotype"],
               "kidney": ["renal/urinary system phenotype"], "liver": ["liver/biliary system phenotype"],
               "lipids": ["homeostasis/metabolism phenotype"], "glucose": ["homeostasis/metabolism phenotype"],
               "lung": ["respiratory system phenotype"], "muscle": ["muscle phenotype"],
               "nervous": ["nervous system phenotype", "behavior/neurological phenotype"],
               "skeletal": ["skeleton phenotype", "limbs/digits/tail phenotype", "growth/size/body region phenotype"],
               "circulatory": ["cardiovascular system phenotype"]}
B = 10000


def load_significance():
    """IMPC: per gene x zygosity x top-level system, whether any parameter was significant."""
    cols = ["marker_symbol", "zygosity", "top_level_mp_term_name", "significant", "status"]
    parts = []
    for ch in pd.read_csv(impc_mouse.STATISTICAL_RESULTS, usecols=cols, chunksize=500_000, dtype=str):
        ch = ch[ch.zygosity.isin(["heterozygote", "homozygote"]) & ch.top_level_mp_term_name.notna() & (ch.status == "Successful")]
        ch = ch.assign(sys=ch.top_level_mp_term_name.str.split(",")).explode("sys")
        ch["sig"] = ch.significant.str.lower().eq("true")
        parts.append(ch.groupby(["marker_symbol", "zygosity", "sys"]).sig.max().reset_index())
    return pd.concat(parts).groupby(["marker_symbol", "zygosity", "sys"]).sig.max().reset_index()


def analyse(R, H, members, rng, zygosity, genes_only=None, label=""):
    x = R[(R.zygosity == zygosity)].merge(H[["mouse", "gene_id"]], left_on="marker_symbol", right_on="mouse")
    x = x[x.gene_id.isin(members.index)]
    if genes_only is not None:
        x = x[x.gene_id.isin(genes_only)]
    Z = x.pivot_table(index="gene_id", columns="sys", values="sig", aggfunc="max").astype(float)
    Z = Z[[c for c in Z.columns if c != "mammalian phenotype"]]
    Zc = Z - Z.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    genes = list(Zc.index)
    systems = list(Zc.columns)
    V = Zc.values
    labels = [members[g] for g in genes]

    def mask(labs):
        return impc_mouse.organ_cell_mask(labs, systems, V, ORGAN_TO_MP)
    m = mask(labels)
    obs = V[m].mean()
    raw_match = Z.values[m].mean()
    raw_other = Z.values[~m & np.isfinite(Z.values)].mean()
    perm = np.empty(B)
    for b in range(B):
        perm[b] = V[mask([labels[i] for i in rng.permutation(len(labels))])].mean()
    p = (1 + (perm >= obs).sum()) / (B + 1)
    print(f"  {label or zygosity:28s} genes {len(genes):4d}, matching cells {m.sum():4d}: hit rate matching {raw_match:.3f} "
          f"vs other {raw_other:.3f} | centred {obs:+.4f} (perm {perm.mean():+.4f} +- {perm.std():.4f}; "
          f"dz {(obs - perm.mean()) / perm.std():+.1f}) p {p:.4f}", flush=True)
    # per organ
    per = []
    for o in ORGAN_TO_MP:
        mo = mask([{o} & lab for lab in labels])
        if mo.sum() >= 10:
            per.append(f"{o} {Z.values[mo].mean():.2f}/{mo.sum()}")
    print("     per organ (hit rate matching / n):", ", ".join(per))
    return dict(test=label or zygosity, gene=len(genes), n=int(m.sum()), rate_matching=raw_match, rate_other=raw_other,
                stat=obs, p=p)


def main():
    rng = np.random.default_rng(38)
    R = load_significance()
    print("IMPC rows gene x zygosity x system:", len(R), "| genes het", R[R.zygosity == "heterozygote"].marker_symbol.nunique(),
          "hom", R[R.zygosity == "homozygote"].marker_symbol.nunique())
    H = impc_mouse.one_to_one_orthologs()
    members = impc_mouse.organ_membership()

    rows = []
    print("\n=== Primary: heterozygous ===")
    rows.append(analyse(R, H, members, rng, "heterozygote"))
    print("\n=== Positive control: homozygous ===")
    rows.append(analyse(R, H, members, rng, "homozygote"))
    both = set(R[R.zygosity == "heterozygote"].marker_symbol) & set(R[R.zygosity == "homozygote"].marker_symbol)
    genes_both = set(H[H.mouse.isin(both)].gene_id)
    print(f"\n=== Secondary: only genes with het and hom data (n {len(genes_both & set(members.index))}) ===")
    r1 = analyse(R, H, members, rng, "heterozygote", genes_both, "het (both available)")
    r2 = analyse(R, H, members, rng, "homozygote", genes_both, "hom (both available)")
    rows += [r1, r2]
    if r2["rate_matching"] - r2["rate_other"] != 0:
        print(f"  rough mouse h (excess matching het / hom): "
              f"{(r1['rate_matching'] - r1['rate_other']) / (r2['rate_matching'] - r2['rate_other']):.2f}")
    pd.DataFrame(rows).to_csv(config.results("38_mouse_organ_diagonal.tsv"), sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
