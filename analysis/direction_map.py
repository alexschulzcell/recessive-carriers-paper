"""Gene-specific expected direction of the carrier effect, derived from HPO (pre-specified, test family T10-7).

Every recessive gene gets its own expected direction per trait, taken from the direction-bearing HPO terms of its
OMIM diseases (clinical descriptions, independent of UK Biobank). Within one trait, different genes can therefore
expect opposite directions. Cells with contradicting directions are discarded.

The map is used by scripts 31 (UK Biobank test), 33 and 35 (Genes & Health, FinMetSeq, exome chip), 36 (carrier
atlas), 42/43 (AlphaMissense dose), 45 (cross-cohort consistency), 46 and 49/49c (FinnGen, BioBank Japan).

``build_direction_map()`` returns a namespace with
    directions   {gene_id: {trait: +1/-1}}
    genes        sorted gene_ids with at least one direction (rows of the matrices)
    traits       sorted traits of the HPO direction map (columns)
    S            genes x traits matrix of expected signs (+1/-1/0)
    raw          {mask: genes x traits matrix of signed UK Biobank burden z (Genebass)}
    centred      {mask: same, centred over all recessive genes, first per trait, then per gene}
    recessive_genes  set of recessive gene_ids (same gene set as scripts 14/20/30)
    annotation   gene annotation table (script 03)
    bidirectional_traits  traits with genes expected in both directions
Masks: plof (pLoF), mislc (missense|LC), syn (synonymous, negative control).
"""
from collections import defaultdict
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import norm

import config

# Pre-specified (PRAEREG_vorzeichen.md, part C): HPO term -> {trait: expected direction}
HPO_DIRECTION_MAP = {
    "HP:0001903": {"hb": -1},                                    # Anemia
    "HP:0001935": {"mcv": -1, "hb": -1},                          # Microcytic anemia
    "HP:0001889": {"mcv": +1, "hb": -1},                          # Macrocytic anemia
    "HP:0001873": {"plt": -1}, "HP:0001894": {"plt": +1},
    "HP:0001882": {"wbc": -1}, "HP:0001875": {"wbc": -1}, "HP:0001974": {"wbc": +1},
    "HP:0001888": {"lymph": -1},
    "HP:0003259": {"creat": +1, "cysc": +1},
    "HP:0012622": {"creat": +1, "cysc": +1, "urea": +1},
    "HP:0002905": {"phosphate": +1}, "HP:0002148": {"phosphate": -1},
    "HP:0002149": {"urate": +1}, "HP:0003537": {"urate": -1},
    "HP:0002904": {"bili": +1}, "HP:0002908": {"bili": +1},
    "HP:0008282": {"bili": +1}, "HP:0003573": {"bili": +1},
    "HP:0002910": {"alt": +1, "ast": +1}, "HP:0003234": {"ggt": +1},
    "HP:0003073": {"albumin": -1},
    "HP:0003233": {"hdl": -1}, "HP:0002155": {"tg": +1},
    "HP:0003124": {"chol": +1, "ldl": +1},
    "HP:0003563": {"ldl": -1, "chol": -1}, "HP:0008181": {"ldl": -1, "chol": -1},
    "HP:0003074": {"glucose": +1, "hba1c": +1}, "HP:0000857": {"glucose": +1, "hba1c": +1},
    "HP:0001943": {"glucose": -1},
    "HP:0002097": {"fev1": -1, "fvc": -1}, "HP:0006536": {"fev1": -1, "fvc": -1},
    "HP:0002113": {"fev1": -1, "fvc": -1}, "HP:0002206": {"fev1": -1, "fvc": -1},
    "HP:0002205": {"fev1": -1, "fvc": -1},
    "HP:0003701": {"grip": -1}, "HP:0003324": {"grip": -1}, "HP:0003202": {"grip": -1},
    "HP:0001249": {"fluidint": -1, "reaction": +1},
    "HP:0001263": {"fluidint": -1, "reaction": +1},
    "HP:0002376": {"fluidint": -1, "reaction": +1},
    "HP:0004322": {"height": -1, "sitheight": -1}, "HP:0003510": {"height": -1, "sitheight": -1},
    "HP:0008873": {"height": -1, "sitheight": -1}, "HP:0003498": {"height": -1, "sitheight": -1},
    "HP:0000098": {"height": +1, "sitheight": +1},
    "HP:0004349": {"ebmd": -1}, "HP:0000939": {"ebmd": -1},
    "HP:0011001": {"ebmd": +1}, "HP:0011002": {"ebmd": +1},
}
TRAITS = sorted({t for v in HPO_DIRECTION_MAP.values() for t in v})


def hpo_children(ref_dir):
    """{HPO term: set of direct children} from hp.obo (is_a relations)."""
    children = defaultdict(set)
    current = None
    for line in open(ref_dir / "hp.obo", encoding="utf-8"):
        if line.startswith("id: HP:"):
            current = line.split()[1]
        elif line.startswith("is_a:") and current:
            children[line.split()[1]].add(current)
    return children


def descendants(term, children):
    out, stack = {term}, [term]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def signed_ukb_z(annotation, trait, mask, index):
    """Signed Genebass burden z per gene (sign of BETA, |z| from the burden p, capped at sqrt(50)), carriers >= 10."""
    f = config.GENEBASS_DIR / f"genebass_{trait}_{mask}.tsv.gz"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
    g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
    g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)]  # noqa: E712
    g = g.drop_duplicates("gene_id").set_index("gene_id")
    z = np.clip(norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2), None, np.sqrt(50))
    return pd.Series(np.sign(g.BETA_Burden.values) * z, index=g.index).reindex(index)


def build_direction_map(ref_dir=config.REF_DIR):
    """Build the gene-specific direction map and the UK Biobank z matrices (prints a short summary)."""
    children = hpo_children(ref_dir)
    expected = defaultdict(lambda: defaultdict(list))
    for term, trait_signs in HPO_DIRECTION_MAP.items():
        for d in descendants(term, children):
            for trait, sign in trait_signs.items():
                expected[d][trait].append(sign)
    print(f"Direction-bearing HPO terms including descendants: {len(expected)}", flush=True)

    hpoa = pd.read_csv(ref_dir / "phenotype.hpoa", sep="\t", comment="#", dtype=str)
    hpoa = hpoa[hpoa.database_id.str.startswith("OMIM") & (hpoa.aspect == "P") & (hpoa.qualifier != "NOT")]
    disease_terms = hpoa.groupby("database_id").hpo_id.apply(set)
    g2d = pd.read_csv(ref_dir / "genes_to_disease.txt", sep="\t", dtype=str)
    g2d = g2d[g2d.association_type == "MENDELIAN"]

    annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
    organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
    recessive_genes = set(organ_fractions.gene_id)            # same recessive gene set as scripts 14/20/30
    sym2id = dict(zip(annotation.symbol, annotation.gene_id))
    gene_terms = g2d.groupby("gene_symbol").disease_id.apply(lambda ds: set().union(*[disease_terms.get(d, set()) for d in ds]))

    directions = {}
    n_conflicts = 0
    for symbol, terms in gene_terms.items():
        gid = sym2id.get(symbol)
        if gid is None or gid not in recessive_genes:
            continue
        acc = defaultdict(list)
        for t in terms:
            for trait, signs in expected.get(t, {}).items():
                acc[trait].extend(signs)
        r = {}
        for trait, signs in acc.items():
            u = set(signs)
            if len(u) == 1:
                r[trait] = u.pop()
            else:
                n_conflicts += 1
        if r:
            directions[gid] = r
    print(f"Recessive genes with a gene-specific direction: {len(directions)} (of {len(recessive_genes)}); "
          f"discarded contradicting cells: {n_conflicts}", flush=True)

    genes = sorted(directions)                                # matrix rows
    S = pd.DataFrame(0.0, index=genes, columns=TRAITS)
    for g, r in directions.items():
        for t, s in r.items():
            S.loc[g, t] = s

    all_recessive = sorted(recessive_genes)                   # centring over all recessive genes (as 14/20/30)
    raw, centred = {}, {}
    for mask in ["plof", "mislc", "syn"]:
        Z = pd.DataFrame({t: z for t in TRAITS if (z := signed_ukb_z(annotation, t, mask, all_recessive)) is not None},
                         index=all_recessive)
        Zc = Z - Z.mean(axis=0)
        Zc = Zc.sub(Zc.mean(axis=1), axis=0)
        raw[mask] = Z.reindex(columns=TRAITS).loc[genes].values
        centred[mask] = Zc.reindex(columns=TRAITS).loc[genes].values

    cells = defaultdict(int)
    for g, r in directions.items():
        for t, s in r.items():
            cells[(t, s)] += 1
    print("Cells per trait x direction:", {f"{t}{'+' if s > 0 else '-'}": n for (t, s), n in sorted(cells.items())}, flush=True)
    bidirectional = sorted({t for (t, s) in cells if (t, -s) in cells})
    print(f"Traits with genes in both directions (the generic loss-of-function direction cannot explain these): "
          f"{bidirectional}", flush=True)

    return SimpleNamespace(directions=directions, genes=genes, traits=TRAITS, S=S.values, raw=raw, centred=centred,
                           recessive_genes=recessive_genes, annotation=annotation, sym2id=sym2id,
                           bidirectional_traits=bidirectional)
