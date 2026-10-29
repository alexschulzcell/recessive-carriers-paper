#!/usr/bin/env bash
# Download the raw summary statistics of the independent replication cohorts (free, no login required):
#   bcx/        Blood Cell Consortium exome-chip results (Chami, Eicher, Tajuddin 2016; archived GRASP copies)
#   ckdgen/     CKDGen exome-chip eGFRcrea/eGFRcys results (Li et al. 2017, European ancestry)
#   finmetseq/  FinMetSeq exome sequencing of about 19,000 Finns (GWAS Catalog GCST008673, hg19)
#   gh/         Genes & Health exome variant-level results (British South Asians; GWAS Catalog, harmonised GRCh38 files)
#   annot/      gnomAD v2 pLoF annotation with misannotation probabilities (Zenodo 10403680) and the hg19->hg38 liftOver chain
# Target folder: $RC_DATA_DIR/replication/<subfolder> (RC_DATA_DIR defaults to <repository>/data/raw, as in config.py).
# Files that already exist and are not empty are skipped.
# Usage: bash 24_download_replication_cohorts.sh
set -u
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="${RC_DATA_DIR:-$REPO_DIR/data/raw}"
Z="$DATA_DIR/replication"; mkdir -p "$Z"/bcx "$Z"/ckdgen "$Z"/finmetseq "$Z"/gh "$Z"/annot; cd "$Z"
get() { [ -s "$2" ] || { curl -sL --retry 3 -m 1800 -o "$2" "$1" && echo "ok $2 $(stat -c %s "$2")" || echo "FAILED $2"; sleep 3; }; }
WA="https://web.archive.org/web"
get "$WA/20161216051549id_/https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Chami/BCX_RBC_ALL.tar.gz" bcx/BCX_RBC_ALL.tar.gz
get "$WA/20161203144046id_/https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Eicher/BCX_PLT_ALL.tar.gz" bcx/BCX_PLT_ALL.tar.gz
get "$WA/20161216050429id_/https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Tajuddin/BCX_WBC_ALL.tar.gz" bcx/BCX_WBC_ALL.tar.gz
for f in Published_eGFRcrea_overall_EA Published_eGFRcys_overall_EA; do get "https://ckdgen.imbi.uni-freiburg.de/files/Li2017/$f.csv.gz" ckdgen/$f.csv.gz; done
FM="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008673"
for t in Height Creatinine eGFR LDL-C HDL-C Total_TG; do get "$FM/LockeAE_prePMID_${t}_sex-combined.gz" finmetseq/$t.gz; done
# Genes & Health: trait code -> last three digits of the GWAS Catalog accession GCST90727xxx
GH="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90727001-GCST90728000"
for p in height:327 hb:323 mcv:336 wbc:370 plt:343 creat:308 alt:300 bili:304 ldl:331 hdl:325 tg:360 hba1c:322; do
  k=${p%%:*}; n=${p##*:}; get "$GH/GCST90727$n/harmonised/GCST90727$n.h.tsv.gz" gh/$k.h.tsv.gz; done
get "https://zenodo.org/api/records/10403680/files/gnomad_lofs_with_misannotation_probabilities.tsv.gz/content" annot/gnomad_v2_lofs.tsv.gz
get "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" annot/hg19ToHg38.over.chain.gz
echo done
