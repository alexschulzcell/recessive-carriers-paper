"""Central path configuration for all analysis scripts.

Two environment variables control where the scripts look for inputs and where they write results:

    RC_DATA_DIR     folder with the raw inputs (default: <repository>/data/raw)
    RC_RESULTS_DIR  folder for all results, logs and intermediate files (default: <repository>/results)

The expected layout of RC_DATA_DIR is described in data/README.md. No script writes into RC_DATA_DIR
except the download scripts (01, 06, 13, 39a, 49a, 49b and the two shell download scripts).
"""
import os
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("RC_DATA_DIR", REPO_DIR / "data" / "raw"))
RESULTS_DIR = Path(os.environ.get("RC_RESULTS_DIR", REPO_DIR / "results"))

# raw inputs (see data/README.md for source, version and checksum of every file)
REF_DIR = DATA_DIR / "reference"                    # HPO, OMIM via HPO, HGNC, GO, HPA, refGene, nosology XML, Alliance, SGD
GENEBASS_DIR = DATA_DIR / "genebass"                # genebass_<trait>_<mask>.tsv.gz, genebass_qc_<mask>.tsv.gz
GENEBASS_ICD_DIR = DATA_DIR / "genebass_icd10"      # <ICD-10 code>_<mask>.tsv.gz
AZPHEWAS_DIR = DATA_DIR / "azphewas"                # az_<trait>_<model>.tsv.gz
GENEBAYES_FILE = DATA_DIR / "genebayes" / "shet_genebayes.tsv"
IMPC_API_DIR = DATA_DIR / "impc_api_heterozygous"   # impc_het_Body.tsv.gz, impc_het_Bone.tsv.gz
IMPC_RELEASE_DIR = DATA_DIR / "impc_release"        # statistical-results-ALL.csv.gz, viability.csv.gz
EXOME_CHIP_DIR = DATA_DIR / "exome_chip"            # GIANT height and GLGC lipid exome-chip files
REPLICATION_DIR = DATA_DIR / "replication"          # annot/, gh/, finmetseq/, bcx/, ckdgen/, finngen/, bbj/
GH_BURDEN_DIR = DATA_DIR / "genes_and_health_burden"  # Genes & Health REGENIE gene-burden files
YEAST_HOEPFNER_DIR = DATA_DIR / "yeast_hoepfner"
YEAST_HILLENMEYER_DIR = DATA_DIR / "yeast_hillenmeyer"
ALPHAMISSENSE_DIR = DATA_DIR / "alphamissense"
BACKMAN_DIR = DATA_DIR / "backman"
LITERATURE_DIR = REPO_DIR / "data" / "literature_extraction"  # batch*.tsv: literature extraction of patient values (shipped)

# frozen atlas prediction (Table S5), kept byte-identical in tables/ with its original column names
FROZEN_PREDICTION = REPO_DIR / "tables" / "ATLAS_VORHERSAGE_2026-10-20.tsv"
FROZEN_PREDICTION_SHA256_PREFIX = "7694042a7c18266b"

# pre-specified selection of 122 ICD-10 first-occurrence phenotypes (PRAEREG_diagnosen.md), shipped with the repository
ICD10_SELECTION = REPO_DIR / "data" / "icd10_selection.tsv"

# results
CACHE_DIR = RESULTS_DIR / "cache"


def results(name):
    """Path of a result file in RC_RESULTS_DIR (the folder is created if necessary)."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    return RESULTS_DIR / name


def cache(name):
    """Path of a cache folder in RC_RESULTS_DIR/cache (created if necessary)."""
    p = CACHE_DIR / name
    p.mkdir(parents=True, exist_ok=True)
    return p
