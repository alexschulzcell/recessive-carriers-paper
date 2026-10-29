"""Table S3 (a priori definitions), written directly from the definitions in the analysis scripts.

47_st7a_gene_classes.tsv         gene classes MATRIX/LYSOSOME/CILIUM (script 10), SECRETED/CATABOLIC/CILIUM_ORGANS
                                 (script 21), W1 GO terms per organ (22), W2 tissues per organ (22), W4 producer cell types
                                 per organ (27); GO names from go-basic.obo, gene counts (script 10) from go_classes.tsv;
                                 the column 'obsolete' marks terms that are obsolete in the GO release used (then without
                                 descendants or annotations).
47_st7b_organ_to_mp_terms.tsv    organ -> MP top-level terms (ORGAN_TO_MP of script 38).
47_st7c_organ_to_icd10_blocks.tsv  organ block -> the 122 pre-selected ICD-10 codes (config.ICD10_SELECTION) with the
                                 pre-specified chapter ranges (PRAEREG_diagnosen.md); checks that every code lies within the
                                 range of its organ.
The definitions are read from the scripts with ast (no transcription, nothing is executed).
Script 48 copies these tables into the supplemental table folder.
Call: python 47_definition_tables.py
"""
import ast
import re
from pathlib import Path

import pandas as pd

import common
import config

HERE = Path(__file__).resolve().parent


def literal(script, name):
    """Literal assignment `name = ...` from a script (ast, without executing it)."""
    for n in ast.parse((HERE / script).read_text(encoding="utf-8")).body:
        if isinstance(n, ast.Assign) and any(getattr(t, "id", None) == name for t in n.targets):
            return ast.literal_eval(n.value)
    raise KeyError(f"{name} not in {script}")


def go_lists(script, name):
    """GO roots per class of a class assignment such as GO_CLASSES / MECHANISM_CLASSES in scripts 10/21."""
    source = (HERE / script).read_text(encoding="utf-8")
    for n in ast.parse(source).body:
        if isinstance(n, ast.Assign) and any(getattr(t, "id", None) == name for t in n.targets):
            return {ast.literal_eval(k): re.findall(r"GO:\d{7}", ast.get_source_segment(source, v))
                    for k, v in zip(n.value.keys, n.value.values)}
    raise KeyError(f"{name} not in {script}")


def write(df, name):
    df.to_csv(config.results(name), sep="\t", index=False, lineterminator="\n")


go_names, cur = {}, None
for line in open(config.REF_DIR / "go-basic.obo", encoding="utf-8"):
    if line.startswith("["):
        cur = None
    elif line.startswith("id: GO:"):
        cur = line.split()[1]
    elif line.startswith("name: ") and cur and cur not in go_names:
        go_names[cur] = line[6:].strip()
    elif line.startswith("namespace: ") and cur:
        go_names[cur + "_ns"] = line.split()[1]

# --- ST7a
GROUP_LABEL = {"LYSO": "LYSOSOME"}   # class label of script 10 as shown in the table
rows = []
n_class = pd.read_csv(config.results("go_classes.tsv"), sep="\t").category.value_counts()
for k, terms in go_lists("10_matrix_flux_go_classes.py", "GO_CLASSES").items():
    for t in terms:
        rows.append(dict(part="skeletal mechanism (script 10)", group=GROUP_LABEL.get(k, k), identifier=t,
                         name=go_names.get(t), ontology=go_names.get(t + "_ns"),
                         rule="GO term including descendants (is_a, part_of); precedence MATRIX > LYSOSOME > CILIUM",
                         n_genes_annotated=int(n_class.get(k, 0))))
for k, terms in go_lists("21_mechanism_across_organs.py", "MECHANISM_CLASSES").items():
    for t in terms:
        rows.append(dict(part="mechanism across organs (script 21, pre-specified)", group=k, identifier=t,
                         name=go_names.get(t), ontology=go_names.get(t + "_ns"),
                         rule="GO term including descendants; precedence SECRETED > CATABOLIC > CILIUM"))
for o, terms in literal("22_main_product_rule.py", "MAIN_PRODUCT_GO").items():
    for t in terms:
        rows.append(dict(part="W1 Gene Ontology main product (script 22, pre-specified)", group=o, identifier=t,
                         name=go_names.get(t), ontology=go_names.get(t + "_ns"),
                         rule="GO term including descendants"
                              + ("; skeleton = reference (MATRIX from script 10), not part of the test"
                                 if o == "skeletal" else "")))
for o, terms in literal("22_main_product_rule.py", "MAIN_PRODUCT_TISSUES").items():
    for t in terms:
        rows.append(dict(part="W2a/W2b tissue (script 22, pre-specified)", group=o, identifier=t, name=t,
                         ontology="HPA rna_tissue_consensus",
                         rule="maximum nTPM over the organ's tissues; specificity log2((nTPM+1)/(median over tissues+1)); "
                              "top tertile within organ"))
for o, terms in literal("27_producer_cell_type_w4.py", "PRODUCER_CELLS").items():
    for t in terms:
        rows.append(dict(part="W4 producer cell (script 27, pre-specified)", group=o, identifier=t, name=t,
                         ontology="HPA rna_single_cell_type",
                         rule="maximum nCPM over the organ's cell types; specificity log2((nCPM+1)/(median over cell types+1)); "
                              "top tertile within organ"))
T = pd.DataFrame(rows).astype({"n_genes_annotated": "Int64"})
T["obsolete"] = T.name.fillna("").str.startswith("obsolete ")
write(T, "47_st7a_gene_classes.tsv")
print("ST7a:", T.groupby("part", sort=False).size().to_dict(),
      "| GO names missing:", int(T[T.ontology != T.ontology].shape[0]),
      "| obsolete in go-basic.obo:", T[T.obsolete].identifier.tolist())

# --- ST7b
ORGAN_TO_MP = literal("38_mouse_organ_diagonal.py", "ORGAN_TO_MP")
T = pd.DataFrame([dict(organ=o, mp_top_level_term=m) for o, ms in ORGAN_TO_MP.items() for m in ms])
write(T, "47_st7b_organ_to_mp_terms.tsv")
print("ST7b:", len(T), "rows,", T.organ.nunique(), "organs,", T.mp_top_level_term.nunique(), "MP terms")

# --- ST7c
# pre-specified ICD-10 chapter ranges per organ (PRAEREG_diagnosen.md, fixed before loading the data)
BLOCK = {"blood": "D50-D77", "immune": "D80-D89", "kidney": "N00-N29", "liver": "K70-K83", "lipids": "E78",
         "glucose": "E10-E16", "lung": "J40-J99", "muscle": "G70-G73, M60-M63",
         "nervous": "G00-G69, G74-G99, F00-F03, F70-F79", "skeletal": "M15-M19, M40-M43, M80-M94, Q65-Q79",
         "circulatory": "I00-I99"}


def in_block(code, blocks):
    for r in blocks.split(", "):
        a, z = (r.split("-") + [r])[:2]
        if a[0] == code[0] and int(a[1:]) <= int(code[1:3]) <= int(z[1:]):
            return True
    return False


S = pd.read_csv(config.ICD10_SELECTION, sep="\t")
S["pre-specified_blocks"] = S.organ.map(BLOCK)
S["in_block"] = [in_block(c, b) for c, b in zip(S.icd, S["pre-specified_blocks"])]
S["description_short"] = S.description.str.extract(r"\((.*)\)$")[0]
T = (S[["organ", "pre-specified_blocks", "icd", "description_short", "n_cases", "analysis_id", "in_block"]]
     .sort_values(["organ", "icd"], key=common.organ_column_key)
     .rename(columns={"icd": "icd10", "description_short": "description"}))
write(T, "47_st7c_organ_to_icd10_blocks.tsv")
codes_per_organ = T.groupby("organ").size()
print("ST7c:", len(T), "codes;", codes_per_organ.reindex(common.sort_organs(codes_per_organ.index)).to_dict(),
      "| outside the block:", T[~T.in_block].icd10.tolist())
