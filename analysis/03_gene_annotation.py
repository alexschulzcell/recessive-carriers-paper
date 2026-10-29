"""Gene annotation for all Genebass genes and definition of the gene sets.

For every gene in the Genebass gene QC this script collects: gene_id, symbol, nosology class (nos_class: AR_only, AD_only,
AD_AR, XL, other; script 02) and nosology groups, Mendelian status from HPO/OMIM (inheritance of the OMIM diseases;
skeletal involvement = any phenotype term below HP:0000924 abnormality of the skeletal system, HP:0004322 short stature
or HP:0000098 tall stature), gnomAD v2.1.1 LOEUF and CDS length, Genebass CAF per burden set (carriers ~ 2 * CAF * 394,841)
and the Genebass keep flags. From these it assigns each gene to one gene set: AR_nos (recessive skeletal dysplasia genes),
AD_nos, nos_other, AR_skelsyndr (recessive Mendelian genes outside the nosology with skeletal HPO terms), AR_nonskel
(recessive Mendelian genes without skeletal HPO terms), AD_mendel_other, mendel_mixed, and pool (no Mendelian gene).

Inputs: config.GENEBASS_DIR / genebass_qc_<set>.tsv.gz; nosology2023_genes.tsv (script 02); config.REF_DIR / hgnc.txt,
        gnomad_v2_lof.txt.bgz, phenotype.hpoa, genes_to_disease.txt, hp.obo
Output: gene_annotation.tsv.gz in RC_RESULTS_DIR (used by almost all later scripts)
Usage:  python 03_gene_annotation.py
"""
from collections import defaultdict

import pandas as pd

import config

REF = config.REF_DIR
N = 394841  # Genebass exomes

# --- Genebass gene universe and CAF
qc = {}
for s in ["plof", "mislc", "syn"]:
    q = pd.read_csv(config.GENEBASS_DIR / f"genebass_qc_{s}.tsv.gz", sep="\t")
    qc[s] = q.set_index("gene_id")
U = pd.DataFrame({"gene_id": sorted(set().union(*[set(q.index) for q in qc.values()]))})
U["symbol"] = U.gene_id.map(pd.concat([q.gene_symbol for q in qc.values()]).groupby(level=0).first())
for s, q in qc.items():
    U[f"caf_{s}"] = U.gene_id.map(q.CAF)
    U[f"carriers_{s}"] = (2 * U[f"caf_{s}"] * N).round()
    U[f"keep_{s}"] = U.gene_id.map(q.keep_gene_burden & q.keep_gene_coverage & q.keep_gene_n_var)

# --- HGNC: symbol / alias / previous symbol -> Ensembl
h = pd.read_csv(REF / "hgnc.txt", sep="\t", dtype=str,
                usecols=["symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id", "locus_group", "entrez_id"])
h = h[h.ensembl_gene_id.notna()]
sym2ens, alt2ens = {}, defaultdict(set)
for r in h.itertuples():
    sym2ens[r.symbol] = r.ensembl_gene_id
    for col in (r.alias_symbol, r.prev_symbol):
        if isinstance(col, str):
            for a in col.split("|"):
                alt2ens[a].add(r.ensembl_gene_id)
entrez2ens = dict(zip(h.entrez_id, h.ensembl_gene_id))
U["protein_coding"] = U.gene_id.isin(set(h[h.locus_group == "protein-coding gene"].ensembl_gene_id))
# manual corrections of nosology symbols (typos / outdated names); None = no protein-coding gene (loci, enhancers, miRNA)
manual = {"CCD134": "CCDC134", "DHPAT": "GNPAT", "FUCA": "FUCA1", "HSGNAT": "HGSNAT", "IFT40": None, "CGAF": None,
          "MAENLI": None, "HOXD": None, "ZRS": None, "MIR140": None, "OPG": "TNFRSF11B"}


def to_ensembl(g):
    """Ensembl id of a nosology gene symbol (approved symbol, else unique alias/previous symbol)."""
    g = manual.get(g, g)
    if g is None:
        return None
    if g in sym2ens:
        return sym2ens[g]
    a = alt2ens.get(g, set())
    return next(iter(a)) if len(a) == 1 else None


nos = pd.read_csv(config.results("nosology2023_genes.tsv"), sep="\t")
nos["gene_id"] = nos.gene.map(to_ensembl)
print("Nosology genes without Ensembl id:", sorted(nos[nos.gene_id.isna()].gene))
nosm = nos.dropna(subset=["gene_id"]).groupby("gene_id").agg(nos_symbol=("gene", "first"), nos_class=("nos_class", "first"),
                                                             nos_groups=("groups", "first"), nos_inh=("inh", "first"))
U = U.merge(nosm, left_on="gene_id", right_index=True, how="left")
print("Nosology genes in Genebass:", U.nos_class.value_counts().to_dict())

# --- HPO: skeletal subtree
parents = defaultdict(set)
cur = None
for line in open(REF / "hp.obo", encoding="utf-8"):
    if line.startswith("id: HP:"):
        cur = line.split()[1]
    elif line.startswith("is_a:") and cur:
        parents[cur].add(line.split()[1])
children = defaultdict(set)
for c, ps in parents.items():
    for p in ps:
        children[p].add(c)


def descendants(root):
    """HPO term and all its descendants."""
    out, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


SKELETAL = descendants("HP:0000924") | descendants("HP:0004322") | descendants("HP:0000098")
ha = pd.read_csv(REF / "phenotype.hpoa", sep="\t", comment="#", dtype=str)
ha = ha[ha.database_id.str.startswith("OMIM")]
inh = ha[ha.aspect == "I"].groupby("database_id").hpo_id.apply(set)
skel = ha[(ha.aspect == "P") & (ha.qualifier != "NOT")].groupby("database_id").hpo_id.apply(lambda s: len(set(s) & SKELETAL))
g2d = pd.read_csv(REF / "genes_to_disease.txt", sep="\t", dtype=str)
g2d = g2d[(g2d.association_type == "MENDELIAN") & g2d.disease_id.str.startswith("OMIM")]
g2d["gene_id"] = g2d.ncbi_gene_id.str.replace("NCBIGene:", "").map(entrez2ens)
rec = []
for gid, sub in g2d.dropna(subset=["gene_id"]).groupby("gene_id"):
    ar = ad = xl = 0
    sk = 0
    for d in sub.disease_id:
        i = inh.get(d, set())
        # HP:0000007 autosomal recessive, HP:0000006 autosomal dominant, X-linked terms
        ar |= "HP:0000007" in i
        ad |= "HP:0000006" in i
        xl |= bool(i & {"HP:0001417", "HP:0001419", "HP:0001423"})
        sk = max(sk, int(skel.get(d, 0)))
    rec.append((gid, int(ar), int(ad), int(xl), sk, len(sub)))
H = pd.DataFrame(rec, columns=["gene_id", "hpo_ar", "hpo_ad", "hpo_xl", "hpo_skel_terms", "hpo_n_dis"])
U = U.merge(H, on="gene_id", how="left")

# --- gnomAD v2.1.1 constraint (one row per gene: longest CDS)
c = pd.read_csv(REF / "gnomad_v2_lof.txt.bgz", sep="\t", compression="gzip",
                usecols=["gene_id", "oe_lof_upper", "cds_length", "gene_length", "pLI", "mis_z"])
c = c.sort_values("cds_length", ascending=False).drop_duplicates("gene_id")
U = U.merge(c, on="gene_id", how="left")

# --- gene sets
U["set"] = "pool"
U.loc[U.nos_class == "AR_only", "set"] = "AR_nos"
U.loc[U.nos_class == "AD_only", "set"] = "AD_nos"
U.loc[U.nos_class.isin(["AD_AR", "XL", "other"]), "set"] = "nos_other"
m = U.nos_class.isna() & (U.hpo_ar == 1) & (U.hpo_ad == 0) & (U.hpo_xl == 0)
U.loc[m & (U.hpo_skel_terms == 0), "set"] = "AR_nonskel"
U.loc[m & (U.hpo_skel_terms > 0), "set"] = "AR_skelsyndr"
U.loc[U.nos_class.isna() & (U.hpo_ad == 1) & (U.hpo_ar == 0), "set"] = "AD_mendel_other"
U.loc[(U.set == "pool") & U.hpo_n_dis.notna(), "set"] = "mendel_mixed"
U.to_csv(config.results("gene_annotation.tsv.gz"), sep="\t", index=False)
print(U.groupby("set").agg(n=("gene_id", "size"), loeuf=("oe_lof_upper", "median"), cds=("cds_length", "median"),
                           carriers_plof=("carriers_plof", "median"),
                           carriers_mis=("carriers_mislc", "median")).round(2).to_string())
