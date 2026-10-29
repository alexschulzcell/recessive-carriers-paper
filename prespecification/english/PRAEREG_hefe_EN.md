> English translation (21 October 2026) of `prespecification/PRAEREG_hefe.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specified: condition diagonal in yeast – does the heterozygous strain show its defect where the homozygous strain suffers? (stage 16, `work/40_*`)

Specified on 20 October 2026, **before** any computation. Known: format of the files, sign convention (negative = sensitive; TUB2 under benomyl rank 5/6,681), matrix sizes. Not known: any number of this test.

## Analogy to humans
| Human | Yeast |
|---|---|
| biallelic patient | homozygous deletion (HOP) |
| carrier | heterozygous deletion (HIP) |
| disease organ | conditions (compound × concentration) in which the homozygous deletion is strongly sensitive |
| organ diagonal | heterozygous sensitivity in exactly these conditions |

Advantage: an experiment with defined null alleles, isogenic, without problems of linkage disequilibrium, stratification or ascertainment.

## Data
Hoepfner et al. 2014 (Microbiol Res; Dryad doi:10.5061/dryad.v5m8v), `HIP_scores.txt`, `HOP_scores.txt` (genes × experiments, adjusted scores, negative = growth defect). Experiments paired via the key "Exp. <compound>_<concentration>" (HIP and HOP column of the same compound and dose). Genes: non-essential, with HIP and HOP data.

## Definitions
- Sensitivity s = −score (positive = sensitive).
- **"Disease conditions" of a gene:** paired experiments with HOP score ≤ −3. Genes with at least 3 such experiments.
- Heterozygous value: HIP sensitivity, double-centred (experiment mean, then gene mean), as in humans.

## Statistics
- **Primary (T16-1):** mean centred HIP sensitivity in the gene's disease conditions. Null: permute condition sets between genes (B = 1,000), one-sided. Criterion P < 0.05.
- **Direction (T16-2):** fraction of disease conditions with HIP score < 0 (heterozygote also sensitive), against the same permutation.
- **Dose (T16-3, descriptive):** "yeast h" = median over genes of (mean HIP sensitivity in disease conditions − otherwise) / (mean HOP sensitivity there − otherwise). Comparison with humans (h ≈ 0.01 for height, lower bound).
- **Rule (T16-4):** the human rule "half dose shows where the gene works specifically", transferred: genes with **few** disease conditions (condition-specific) against genes with many (pleiotropic). Tertiles of the number of HOP hits; prediction: diagonal (gene value from T16-1) larger in the specific tertile than in the pleiotropic one. Test: Spearman between gene diagonal and number of hits (expectation negative) plus tertile comparison (Mann–Whitney).
- Sensitivity: thresholds −2 and −4; "Ad." scores only; exclusion of experiments with > 5% HOP hits (unspecifically toxic).

## Interpretation, specified in advance
- T16-1 positive: the diagonal holds across species, in the causal experiment. Yeast–human convergence.
- T16-4 positive: the rule transfers as well.
- Null: the yeast condition diagonal is not an analogy to the organ (conditions are not tissues); to be reported as such.

---

# Addendum 20 October 2026: independent replication with Hillenmeyer et al. 2008 (Stanford, Giaever/Nislow; T16-R, `work/41`)

Specified **before** computing. Known: format (z-scores, **positive = sensitive**; ERG11 het under fluconazole rank 2/5,915, TUB3 hom under benomyl rank 8/4,742), het 5,915 genes × 726 experiments (332 conditions), hom 4,742 × 418 (161 conditions). Source: web supplement, preserved in the Internet Archive (snapshot 7 December 2015, `het.z_result_nm.pub`, `hom.z_result_nm.pub`), stored under `hillenmeyer_roh/`.

Unlike Hoepfner, HIP and HOP were measured **separately** (own experiments, own days). Pairing therefore by **compound** (condition name), not by experiment. For each gene × compound: mean z over the experiments of the compound.

- "Disease compounds" of a gene: hom mean z ≥ 3. Genes with ≥ 2 such compounds, with het and hom data.
- T16-R1: mean double-centred het sensitivity in the disease compounds; null: permute compound sets between genes (B = 1,000), one-sided. Criterion P < 0.05.
- T16-R2: fraction het z > 0 in disease compounds against the same permutation.
- T16-R3 (descriptive): h = mean Δhet / mean Δhom (z scale; as in script 40c not interpretable as fitness dominance).
- Sensitivity: thresholds 2 and 4; measurement-day control: het experiments from the same day as a hom hit experiment, but another compound (expectation null).
- Gene concordance with Hoepfner: Spearman of gene diagonals (genes in both data sets), expectation positive.
