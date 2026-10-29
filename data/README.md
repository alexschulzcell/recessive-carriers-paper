# Input data

All analyses use **public summary statistics**; no individual-level data were used. The raw inputs are not redistributed
here (they are large and subject to the terms of their providers). This file describes where each input comes from and
where the scripts expect it.

## Folder layout

The scripts read the raw inputs from one folder, set with the environment variable `RC_DATA_DIR` (default:
`data/raw` inside this repository; see `analysis/config.py`). Expected layout:

| Subfolder of `RC_DATA_DIR` | Content | How to obtain |
|---|---|---|
| `reference/` | HGNC, gnomAD v2.1.1 constraint, HPO ontology and annotations, Gene Ontology, Human Protein Atlas (tissue and single-cell type), UCSC refGene (hg19, hg38), 2023 skeletal nosology (PMC10081954, XML), Alliance orthology, SGD features | `bash analysis/00_download_reference_files.sh` |
| `genebass/` | UK Biobank gene-burden results (Genebass, 394,841 exomes), `genebass_<trait>_<set>.tsv.gz` with set = `plof`, `mislc`, `syn`, and gene QC `genebass_qc_<set>.tsv.gz` | `python analysis/01_download_genebass_skeletal.py`, `python analysis/13_download_genebass_organ_traits.py` (public API, several hours) |
| `genebass_icd10/` | Genebass results for the 122 pre-selected ICD-10 first-occurrence phenotypes (`data/icd10_selection.tsv`) | `python analysis/39a_download_genebass_icd10.py` |
| `azphewas/` | AstraZeneca PheWAS portal, gene-level results per trait and collapsing model (`az_<trait>_<model>.tsv.gz`; the portal returns only associations with P < 0.1) | `python analysis/06_download_azphewas.py` |
| `backman/` | UK Biobank single-variant exome statistics of Backman et al. 2021, one GWAS Catalog file per trait, saved as `<trait>.tsv.gz` | GWAS Catalog; accessions and URLs in `data/backman_accessions.tsv` |
| `replication/gh/`, `bcx/`, `ckdgen/`, `finmetseq/`, `annot/` | Genes & Health variant-level results (GWAS Catalog GCST90727297–GCST90727371), Blood Cell Consortium, CKDGen and FinMetSeq summary statistics, gnomAD v2 pLoF annotation, hg19-to-hg38 liftOver chain | `bash analysis/24_download_replication_cohorts.sh` |
| `replication/finngen/` | FinnGen R13 (height) and R12 (laboratory values) summary statistics, reduced to the variants needed | `python analysis/49a_download_finngen.py` |
| `replication/bbj/` | BioBank Japan summary statistics (NBDC Human Database, accession hum0197), reduced to the variants needed | `python analysis/49b_download_biobank_japan.py` |
| `genes_and_health_burden/` | Genes & Health 44k exome gene-burden results (REGENIE), one file per trait | `gs://genesandhealth_publicdatasets/results_44k_ExWAS/44kExWAS_quanttraits/regenie` (public bucket; CC BY-SA 4.0) |
| `exome_chip/` | GIANT exome-chip height (Marouli et al. 2017, discovery stage) and GLGC exome lipids (Liu et al. 2017) | GIANT consortium data files; https://csg.sph.umich.edu/willer/public/lipids2017/ |
| `genebayes/` | GeneBayes *s*<sub>het</sub> estimates (`shet_genebayes.tsv`) | Zeng et al. 2024, Zenodo record 10403680 |
| `impc_api_heterozygous/` | IMPC statistical results for heterozygous lines, body and bone procedures (`impc_het_Body.tsv.gz`, `impc_het_Bone.tsv.gz`) | IMPC statistical-result API (https://www.mousephenotype.org), zygosity heterozygote |
| `impc_release/` | IMPC latest release reports (`statistical-results-ALL.csv.gz`, `viability.csv.gz`, `genotype-phenotype-assertions-ALL.csv.gz`) and MGI orthology `HMD_HumanPhenotype.rpt` | https://www.mousephenotype.org/data/release (FTP reports); https://www.informatics.jax.org/downloads/reports/ |
| `alphamissense/` | AlphaMissense predictions, `AlphaMissense_hg19.tsv.gz` and `AlphaMissense_hg38.tsv.gz` | Zenodo record 8208688 (CC BY-NC-SA 4.0; not redistributed) |
| `yeast_hoepfner/` | HIP and HOP scores of Hoepfner et al. 2014 (`HIP_scores.txt`, `HOP_scores.txt`; Parquet copies are read if present) | Dryad, doi:10.5061/dryad.v5m8v (CC0) |
| `yeast_hillenmeyer/` | Web supplement of Hillenmeyer et al. 2008 (`het.z_result_nm.pub`, `hom.z_result_nm.pub`; Parquet copies are read if present) | Internet Archive snapshot of December 7, 2015 of the original supplement site |

Files shipped in this folder:

| File | Content |
|---|---|
| `icd10_selection.tsv` | The 122 ICD-10 codes selected before analysis, with their organ block (pre-specification `PRAEREG_diagnosen.md`) |
| `literature_extraction/batch*.tsv` | Patient values extracted from the literature by large-language-model agents, five batches, as extracted (column headers and free-text notes in German; read by script 32) |
| `literature_notes_english.tsv` | English versions of the free-text notes of the literature extraction (used for Table S4) |
| `backman_accessions.tsv` | GWAS Catalog accession and URL of each Backman et al. trait |
| `SHA256SUMS.txt` | SHA-256 checksums of every raw input in the layout above, as used for the paper |

## Checking the inputs

From inside `RC_DATA_DIR`:

```bash
sha256sum -c /path/to/repository/data/SHA256SUMS.txt
```

Providers update some resources (HPO, Gene Ontology, Alliance orthology, IMPC). A new download can therefore differ
from the files used for the paper; the checksums show which files changed. The versions used are: HPO release
2026-09-01, Gene Ontology release 2026-07-26, Human Protein Atlas v25.1, Alliance orthology release 9.0.0, Genebass and
AstraZeneca PheWAS portal downloaded October 19–20, 2026, FinnGen and BioBank Japan downloaded October 28, 2026.

## Terms of use

Each resource is subject to its own terms of use. Please cite the original studies (see the reference list of the paper)
and follow the acknowledgment requirements of UK Biobank, Genes & Health, FinnGen and BioBank Japan.
