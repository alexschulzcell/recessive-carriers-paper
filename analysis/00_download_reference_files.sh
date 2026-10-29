#!/usr/bin/env bash
# Download the public reference files used by the whole pipeline (free, no login required).
# Target folder: $RC_DATA_DIR/reference (RC_DATA_DIR defaults to <repository>/data/raw, as in config.py).
# Files that already exist and are not empty are skipped. Download date and the versions of HPO, GO, the HPO
# annotations and the Alliance orthology file are written to VERSIONS.txt.
# Usage: bash 00_download_reference_files.sh
set -u
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="${RC_DATA_DIR:-$REPO_DIR/data/raw}"
Z="$DATA_DIR/reference"; mkdir -p "$Z"; cd "$Z"
get() { [ -s "$2" ] || { curl -sL --retry 3 -m 1800 -o "$2" "$1" && echo "ok $2 $(stat -c %s "$2")" || echo "FAILED $2"; sleep 2; }; }
get "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt" hgnc.txt
get "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz" gnomad_v2_lof.txt.bgz
get "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa" phenotype.hpoa
get "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt" genes_to_disease.txt
get "https://purl.obolibrary.org/obo/hp.obo" hp.obo
get "https://current.geneontology.org/annotations/goa_human.gaf.gz" goa_human.gaf.gz
get "https://current.geneontology.org/ontology/go-basic.obo" go-basic.obo
get "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip" hpa_tissue.zip
get "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip" hpa_singlecell.zip
get "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" refGene_hg19.txt.gz
get "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz" refGene_hg38.txt.gz
# Skeletal nosology 2023 (Unger et al., Am J Med Genet A 2023), author manuscript PMC10081954 as XML
get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=10081954" nosology_pmc10081954.xml
# Yeast-human orthologs (script 45, test P3); official current URLs.
# The analysis used Alliance release 9.0.0 (generated 2026-04-05) and the SGD_features.tab file served on 2026-10-16;
# a new download may return a newer version.
get "https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz" alliance_orthology.tsv.gz
get "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab" SGD_features.tab
{ echo "downloaded $(date -Iseconds)"; grep -m1 "^data-version" hp.obo; grep -m1 "^data-version" go-basic.obo; grep -m1 "#date\|#version" phenotype.hpoa; zcat alliance_orthology.tsv.gz | grep -m1 "Alliance Database Version"; } > VERSIONS.txt
ls -la
