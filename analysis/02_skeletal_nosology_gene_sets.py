"""Parse the 2023 Nosology of genetic skeletal disorders into gene sets.

Source: Unger et al., Am J Med Genet A 2023 (PMID 36779427), author manuscript PMC10081954 as XML (downloaded by
00_download_reference_files.sh). Every table row of a disorder (NOS id) is split into one row per gene. Each gene is then
classified by the modes of inheritance of all its nosology entries: AR_only, AD_only, AD_AR, XL or other. The AR_only
genes are the recessive skeletal dysplasia genes (gene set AR_nos) tested in the skeletal analyses (scripts 03-12).

Input:  config.REF_DIR / nosology_pmc10081954.xml
Output: nosology2023_entries.tsv (one row per disorder x gene) and nosology2023_genes.tsv (per gene: inheritance modes,
        number of disorders, nosology groups, class nos_class), both in RC_RESULTS_DIR.
Usage:  python 02_skeletal_nosology_gene_sets.py
"""
import re
import xml.etree.ElementTree as ET

import pandas as pd

import config

root = ET.parse(config.REF_DIR / "nosology_pmc10081954.xml").getroot()
rows = []
group = None
for tw in root.iter("table-wrap"):
    for tr in tw.iter("tr"):
        c = ["".join(x.itertext()).strip() for x in tr]
        if len(c) < 5:
            continue
        if re.match(r"Group\s*\d+", c[0]):
            group = f"{c[0]}: {c[1]}"
            continue
        if not c[0].startswith("NOS"):
            continue
        inheritance, genes = c[2], c[3]
        for g in re.split(r"[,;/\s]+", genes):
            g = g.strip().strip("()")
            # gene symbols only (no chromosome bands such as 4p16.3)
            if re.fullmatch(r"[A-Z0-9][A-Z0-9\-]{1,14}", g) and not re.fullmatch(r"\d+[pq]?[\d.]*", g):
                # the NOS ids contain a non-breaking space, normalised to a plain space
                rows.append(dict(nos=c[0].replace("\u00a0", " "), group=group, disorder=c[1], inheritance=inheritance,
                                 gene=g, mim=c[4]))
E = pd.DataFrame(rows)
E.to_csv(config.results("nosology2023_entries.tsv"), sep="\t", index=False)


def classify(modes):
    """Gene class from the set of inheritance strings of all its nosology entries."""
    modes = set(modes)
    ad = any(re.search(r"\bAD\b|SMo|\bAD,", x) for x in modes)
    ar = any(re.search(r"\bAR\b", x) for x in modes)
    xl = any(re.search(r"XL", x) for x in modes)
    if ar and not ad and not xl:
        return "AR_only"
    if ad and not ar and not xl:
        return "AD_only"
    if ad and ar:
        return "AD_AR"
    if xl:
        return "XL"
    return "other"


G = E.groupby("gene").agg(inh=("inheritance", lambda s: "|".join(sorted(set(s)))), n_dis=("nos", "nunique"),
                          groups=("group", lambda s: "|".join(sorted(set(s)))))
G["nos_class"] = E.groupby("gene").inheritance.apply(classify)
G.reset_index().to_csv(config.results("nosology2023_genes.tsv"), sep="\t", index=False)
print(len(E), "entries,", len(G), "genes")
print(G.nos_class.value_counts().to_string())
print(E.inheritance.value_counts().head(20).to_string())
