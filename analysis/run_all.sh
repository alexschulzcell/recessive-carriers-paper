#!/usr/bin/env bash
# Runs all analyses of the paper in order, from the raw inputs to the results.
#
# Usage:
#   export RC_DATA_DIR=/path/to/raw-data        # layout: see data/README.md
#   export RC_RESULTS_DIR=/path/to/results      # all results, logs and intermediate files go here
#   bash analysis/run_all.sh [--download]
#
#   --download  also downloads Genebass (01, 13, 39a) and the AstraZeneca PheWAS portal (06) through their public APIs
#               (takes hours). Reference files and replication cohorts are downloaded separately with
#               00_download_reference_files.sh, 24_download_replication_cohorts.sh, 49a_download_finngen.py and
#               49b_download_biobank_japan.py; files that have to be obtained manually are listed in data/README.md.
#
# Each script writes its printed output to $RC_RESULTS_DIR/<script name>.log. The whole run takes about 2 hours.
set -eu
cd "$(dirname "$0")"
: "${RC_DATA_DIR:=$(cd .. && pwd)/data/raw}"
: "${RC_RESULTS_DIR:=$(cd .. && pwd)/results}"
export RC_DATA_DIR RC_RESULTS_DIR
mkdir -p "$RC_RESULTS_DIR"
PY=${PYTHON:-python}

# run <script> [args]: log name = script name (+ optional LOG_SUFFIX)
run() {
  local s=$1; shift
  local log="$RC_RESULTS_DIR/${s%.py}${LOG_SUFFIX:-}.log"
  echo "== $(date +%T) $s $*"
  "$PY" "$s" "$@" > "$log" 2>&1 || { echo "ERROR in $s (see $log)"; tail -5 "$log"; exit 1; }
}

if [ "${1:-}" = "--download" ]; then
  run 01_download_genebass_skeletal.py
  run 06_download_azphewas.py
  run 13_download_genebass_organ_traits.py
  run 39a_download_genebass_icd10.py
fi

# UK Biobank discovery: skeletal dysplasia genes
run 02_skeletal_nosology_gene_sets.py
run 03_gene_annotation.py
run 04_skeletal_enrichment.py 10
run 05_skeletal_robustness.py
run 07_skeletal_sensitivity.py
run 08_gene_x_phenotype_table_and_dose.py
run 09_kacser_burns_enzymes.py
run 10_matrix_flux_go_classes.py
run 11_selection_against_heterozygotes.py
run 12_prediction_new_dysplasia_genes.py
LOG_SUFFIX=_omim run 12_prediction_new_dysplasia_genes.py omim
# all organs
run 14_organ_trait_matrix.py
run 15_mouse_heterozygous_knockouts.py
run 16_body_proportion_and_igf1.py
run 17_non_european_replication.py
run 18_organ_matched_selection.py
run 20_organ_diagonal_robustness.py
run 20b_organ_diagonal_without_chip_genes.py
run 21_mechanism_across_organs.py
run 21b_skeletal_mechanism_corrected.py
run 22_main_product_rule.py
run 22b_main_product_confounders.py
# replication without UK Biobank
run 23_exome_chip_replication.py
run 25_independent_cohort_replication.py
run 26_genes_and_health_burden_replication.py
run 27_producer_cell_type_w4.py
run 28_population_impact.py
# dominance coefficient h
run 32_literature_merge_extraction.py
run 32b_literature_sds_scale_correction.py
run 29_dominance_coefficient_h.py
run 29b_dominance_h_main_test.py
# direction of the carrier effect
run 30_directional_organ_diagonal.py
run 30b_directional_post_hoc.py
run 31_gene_specific_direction.py
run 33_genes_and_health_directional.py
run 34_gnomad_constraint_prediction.py
run 35_directional_variant_replication.py
run 36_carrier_atlas.py
run 37_direction_specificity_genes_and_health.py
# mouse, diagnoses, yeast
run 38_mouse_organ_diagonal.py
run 38b_mouse_organ_diagonal_continuous.py
run 39_icd10_diagnoses_diagonal.py
run 40_yeast_hoepfner_diagonal.py
run 40b_yeast_batch_control.py
run 40c_yeast_conservative_h.py
run 41_yeast_hillenmeyer_replication.py
run 41b_yeast_measurement_day_control.py
# AlphaMissense dose, consistency checks, sensitivity
run 42_alphamissense_dose.py
run 43_alphamissense_backman.py
run 45_cross_cohort_and_yeast_human_consistency.py
run 46_clonal_hematopoiesis_sensitivity.py
# FinnGen and BioBank Japan
LOG_SUFFIX=_finngen run 49_finngen_bbj_replication.py FG
LOG_SUFFIX=_biobank_japan run 49_finngen_bbj_replication.py BBJ
run 49c_finngen_vs_biobank_japan.py
# supplemental tables (Tables S3-S7)
run 47_definition_tables.py
run 48_supplemental_tables.py
echo "== $(date +%T) done"
