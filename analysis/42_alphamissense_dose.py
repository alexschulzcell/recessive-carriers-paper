"""Dose curve within carriers: AlphaMissense pathogenicity versus carrier effect in independent cohorts.

Pre-specified test (PRAEREG_alphamissense.md, tests T17-1 to T17-3). Uses the variant-level summary statistics of the
independent cohorts of script 25 (Genes & Health, FinMetSeq, exome-chip consortia) and the gene-specific expected
directions of script 31 (direction_map.py). T17-1: within cohort x gene x trait cells, does the directional effect d
rise with the AlphaMissense score? T17-2: AlphaMissense classes and pathogen versus benign within cell. T17-3: negative
control with direction maps permuted between genes. See alphamissense_dose.py for the method.

Inputs: replication cohorts and annotation as in script 25 (replication_tests.py), HPO direction map inputs
(direction_map.py), config.ALPHAMISSENSE_DIR/AlphaMissense_{hg19,hg38}.tsv.gz.
Outputs: cache("alphamissense")/alphamissense_variants.tsv.gz (variant table; reused if present),
42_alphamissense_dose.tsv.
Call: python 42_alphamissense_dose.py
"""
import numpy as np
import pandas as pd

import alphamissense_dose as amd
import config
import direction_map
import replication_tests

annotation = replication_tests.VariantAnnotation()
sources = replication_tests.independent_cohort_sources()
dm = direction_map.build_direction_map()
rng = np.random.default_rng(42)

X = amd.load_variant_table(config.cache("alphamissense") / "alphamissense_variants.tsv.gz",
                           sources=sources, builds=["hg19", "hg38"], annotation=annotation, directions=dm.directions,
                           traits=dm.traits, alphamissense_dir=config.ALPHAMISSENSE_DIR)
X, MS = amd.orient(X)
rows, MS = amd.dose_tests(X, MS, rng, dm.directions)
pd.DataFrame(rows).to_csv(config.results("42_alphamissense_dose.tsv"), sep="\t", index=False, float_format="%.4f")
