# Carriers of recessive disease are a mild form of the same disease

Code and data companion to the paper by **Alexander Zentgraf**, the **Genes & Health Research Team** and
**Christian T. Thiel** (Institute of Human Genetics, Universitätsklinikum Erlangen, Friedrich-Alexander-Universität
Erlangen-Nürnberg).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23759254.svg)](https://doi.org/10.5281/zenodo.23759254)

Heterozygous carriers of recessive disease are considered silent. Using exome statistics from 394,841 UK Biobank
participants, we show that carriers of damaging variants in 2,470 genes that cause exclusively recessive disease deviate
in traits of the gene's disease organ, toward the gene's own disease. With analyses specified before the data were seen,
the direction replicated in Genes & Health, FinnGen and BioBank Japan. Carrier effects grow with predicted variant
damage, are strongest in genes specific to the organ's producer cell and, for height, reach about 1% of the patient
effect. In two yeast screens, heterozygous deletions are sensitive under the conditions that impair the homozygous
deletion.

All analyses use **public summary statistics**; no individual-level data were used.

## Contents

| Folder | Content |
|---|---|
| `analysis/` | All analysis scripts, numbered in the order in which they run, and `run_all.sh`, which runs them. Each script starts with a short description of its question, inputs and outputs. |
| `data/` | Description of the raw inputs and where to get them (`data/README.md`), checksums of the inputs used, and the small inputs that are shipped (ICD-10 selection, literature extraction). |
| `prespecification/` | The pre-specification documents (`PRAEREG_*.md`, German originals) and English translations (`english/`). |
| `tables/` | Supplemental Tables S1–S7 of the paper as tab-separated files, and the frozen atlas prediction in its original form. |
| `figures/` | Figures 1–6 (`main/`) and Figures S1–S12 (`supplemental/`) of the paper, as vector PDF and 600-dpi PNG. |

## Reproducing the results

1. Install Python 3.12 and the packages in `requirements.txt` (`pip install -r requirements.txt`). The shell scripts
   need Bash (on Windows, Git Bash works) and `curl`.
2. Download the raw inputs into one folder as described in `data/README.md`, and check them against
   `data/SHA256SUMS.txt`.
3. Run the pipeline:
   ```bash
   export RC_DATA_DIR=/path/to/raw-data        # folder with the raw inputs
   export RC_RESULTS_DIR=/path/to/results      # all results and logs are written here
   bash analysis/run_all.sh
   ```
   The run takes about two hours on a laptop. Every script writes its printed output to
   `$RC_RESULTS_DIR/<script name>.log`, and its tables to `$RC_RESULTS_DIR`.

Permutation tests use a fixed random seed per script, so the results are deterministic. A complete run of this
pipeline from the raw inputs reproduced every number reported in the paper.

## Which script gives which result

Step numbers in the Methods of the paper are the script numbers below.

| Result in the paper | Scripts |
|---|---|
| Gene sets and organ assignment of genes and traits (Methods; Figure S1) | `02`, `03`, `14` |
| Organ-matched carrier signal and its robustness (Figure 1B–D; Figure S2) | `14`, `20`, `20b`, `46` |
| Skeletal dysplasia genes and height (Figure 1E; Figures S1E and S11; Table S6) | `04`, `05`, `07`, `08`, `09`, `17` |
| Direction of carrier effects (Figure 2; Figure S3) | `30`, `30b`, `31` |
| Replication in Genes & Health, FinMetSeq and exome-chip data (Figure 3; Figure S4) | `23`, `25`, `26`, `33`, `35` |
| Replication in FinnGen and BioBank Japan (Figure 3; Figure S4) | `49a`, `49b`, `49`, `49c` |
| Cross-cohort consistency and yeast-human link (Figures 3E and 6F) | `45` |
| Carrier atlas, power and frozen prediction (Figure 3F; Figure S8; Tables S5 and S7) | `05`, `36` |
| Variant damage (AlphaMissense) and carrier effects (Figure 4A–E; Figure S5) | `42`, `43` |
| Dominance coefficient *h* (Figure 4F; Figure S6; Table S4) | `32`, `32b`, `29`, `29b` |
| Mouse, clinical diagnoses and constraint (Figure S7) | `15`, `38`, `38b`, `39`, `34`, `18` |
| Producer-cell rule (Figure 5A–C; Figure S9) | `21`, `21b`, `22`, `22b`, `27`, `37` |
| Growth plate (Figure 5D–F; Figure S12B,C) | `10`, `11`, `12`, `16` |
| Population-level summaries (Figure S12A) | `28` |
| Yeast (Figure 6; Figure S10) | `40`, `40b`, `40c`, `41`, `41b` |
| Supplemental tables (Tables S3–S7) | `47`, `48` |

Shared code: `config.py` (paths), `common.py` (organ definitions), `direction_map.py` (expected direction of each gene
from its disease), `replication_tests.py` (tests used for every replication cohort), `alphamissense_dose.py`,
`impc_mouse.py` and `yeast_data.py`.

Abbreviations used in the scripts: pLoF, predicted loss of function; `mislc`, the Genebass missense and low-confidence
pLoF set; `AR_nos`, genes of the 2023 skeletal nosology listed only with autosomal recessive inheritance; *d*<sub>z</sub>,
(observed − permutation mean)/permutation SD. Test labels (for example T10-7, F4, B4) are those of the
pre-specification documents and of Table S2.

## Pre-specification

For each confirmatory analysis, hypothesis, data, statistic and success criterion were written down before the numbers
were computed. The documents are in `prespecification/` and are listed with their checksums in Table S1. The documents
up to October 20, 2026 are dated internally; the documents for FinnGen (`PRAEREG_finngen.md`) and BioBank Japan
(`PRAEREG_bbj.md`) were committed to this repository on October 28, 2026, before any data of these cohorts were
downloaded. `PRAEREG_stufe19.md` is an excerpt limited to the parts used in the paper. All tests, including failed ones,
are listed in Table S2.

The frozen atlas prediction `tables/ATLAS_VORHERSAGE_2026-10-20.tsv` was written and hashed on October 20, 2026
(SHA-256 `7694042a7c18266b5367732cc2cf988fc46ae4992207682f3e9d54912eb26f67`) and is kept byte-identical, with its
original column names; Table S5 is the same table with English column names.

## Data sources

UK Biobank gene-burden statistics (Genebass), AstraZeneca PheWAS portal, UK Biobank single-variant exome statistics
(Backman et al. 2021), Genes & Health, FinnGen, BioBank Japan, FinMetSeq, exome-chip consortia (GIANT height discovery
stage, GLGC, Blood Cell Consortium, CKDGen), gnomAD, OMIM via the Human Phenotype Ontology, the 2023 skeletal nosology,
Human Protein Atlas, Gene Ontology, HGNC, RefSeq, AlphaMissense, GeneBayes, IMPC and MGI, Alliance of Genome Resources,
SGD, and the yeast chemogenomic screens of Hoepfner et al. 2014 and Hillenmeyer et al. 2008. We thank all participants
and the teams behind these resources. Each resource is subject to its own terms of use; the tables here contain only
aggregate statistics derived from UK Biobank summary statistics and values from the literature.

## Licence

- Code (`analysis/`): MIT licence, see `LICENSE`.
- Figures, tables, documents and data documentation: CC BY 4.0, see `LICENSE-DATA.md`.

## Citation

Please cite the paper (see `CITATION.cff`) and, for the code and data, the Zenodo record of this version:
https://doi.org/10.5281/zenodo.23759254 (all versions: https://doi.org/10.5281/zenodo.23527672).

## Versions

Version 2.0.0 (October 2026) accompanies the submitted manuscript and replaces earlier, withdrawn versions. The
repository history was rewritten for this version: the earlier commits contain only the pre-specification documents and
the frozen prediction, unchanged and with their original commit dates.
