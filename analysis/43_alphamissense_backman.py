"""AlphaMissense dose tests with UK Biobank single-variant results (Backman et al. 2021).

Pre-specified (PRAEREG_alphamissense.md, addendum). Same method as script 42 (tests T17-1 to T17-3, see
alphamissense_dose.py), applied to the Backman et al. single-variant exome association results (source BKM, GRCh38).
In addition, T17-4 (pre-specified primary test for gene specificity): in traits with genes expected in both directions,
the expected directions are permuted between genes within each trait; a generic loss-of-function direction cannot
produce the observed slope.

Inputs: config.BACKMAN_DIR/<trait>.tsv.gz, variant annotation as in script 25 (replication_tests.py), HPO direction map
inputs (direction_map.py), config.ALPHAMISSENSE_DIR/AlphaMissense_hg38.tsv.gz.
Outputs: cache("alphamissense")/alphamissense_variants_backman.tsv.gz (variant table; reused if present),
43_alphamissense_backman.tsv.
Call: python 43_alphamissense_backman.py
"""
import numpy as np
import pandas as pd

import alphamissense_dose as amd
import config
import direction_map
import replication_tests

annotation = replication_tests.VariantAnnotation()
dm = direction_map.build_direction_map()
rng = np.random.default_rng(43)


def read_backman(k):
    d = pd.read_csv(config.BACKMAN_DIR / f"{k}.tsv.gz", sep="\t",
                    usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                             "effect_allele_frequency"], dtype={"chromosome": str})
    f = d.effect_allele_frequency
    return pd.DataFrame({"ch": d.chromosome, "pos": d.base_pair_location, "ref": d.other_allele, "alt": d.effect_allele,
                         "z": d.beta / d.standard_error, "maf": np.minimum(f, 1 - f)})


sources = {"BKM": [(k, "hg38", (lambda k=k: read_backman(k)))
                   for k in sorted(p.name.split(".")[0] for p in config.BACKMAN_DIR.glob("*.tsv.gz"))]}

X = amd.load_variant_table(config.cache("alphamissense") / "alphamissense_variants_backman.tsv.gz",
                           sources=sources, builds=["hg38"], annotation=annotation, directions=dm.directions,
                           traits=dm.traits, alphamissense_dir=config.ALPHAMISSENSE_DIR)
X, MS = amd.orient(X)
rows, MS = amd.dose_tests(X, MS, rng, dm.directions)

print()
print("=== T17-4 (pre-specified, primary for gene specificity): bidirectional traits, directions permuted within trait ===")
bidirectional = [t for t, x in MS.groupby("trait") if set(x.s) == {-1, 1}]
Bd = MS[MS.trait.isin(bidirectional)].copy()
obs4 = amd.within_slope(Bd)
cl = Bd[["gene_id", "trait", "s"]].drop_duplicates(["gene_id", "trait"])
null = []
for b in range(500):
    c = cl.copy()
    c["s2"] = c.groupby("trait").s.transform(lambda v: rng.permutation(v.values))
    m = Bd.merge(c[["gene_id", "trait", "s2"]], on=["gene_id", "trait"])
    m["dn"] = m.s2 * m.zo
    null.append(amd.within_slope(m, "dn"))
null = np.array(null)
p4 = (1 + (null >= obs4).sum()) / 501
print(f"  traits {bidirectional}: n {len(Bd)}; beta {obs4:+.4f} vs perm {null.mean():+.4f} +- {null.std():.4f}; "
      f"dz {(obs4 - null.mean()) / null.std():+.1f}; p {p4:.4f}")
for (t, sg), x in Bd.groupby(["trait", "s"]):
    print(f"    {t:8s} expected {sg:+d}: slope zo~am {amd.within_slope(x.assign(zz=x.zo), 'zz'):+.4f} (n {len(x)})")
rows.append(dict(test="T17-4 bidirectional", n=len(Bd), beta=obs4, p=p4))
pd.DataFrame(rows).to_csv(config.results("43_alphamissense_backman.tsv"), sep="\t", index=False, float_format="%.4f")
