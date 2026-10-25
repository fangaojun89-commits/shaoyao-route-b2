# GSE298414: within-source extraction-method sensitivity

25 October 2026. New-source gene comparisons, the source-label check and all six QC-triggered library omissions are complete. Saved-result checks passed without refitting the models. Existing canonical JGW03 results and candidate gates are unchanged.

## Design and interpretation

The analysis uses 72 PC12 libraries from GSE298414: two extraction solvents, two PR:GR ratios, two extraction methods and three reported ordinal levels, with three biological replicates in each of 24 conditions. All have recorded batch A. Twelve contrasts compare co-extraction (CEM) with post-extraction mixing (IEM) within a ratio, solvent and ordinal level. No untreated or cross-accession control is used. The source does not contain the canonical 1:1 CEM preparation JGW03.

The conflicting absolute dose labels remain unresolved; comparisons use low, middle and high ordering only. These are contrasts among deposited preparation conditions, not dose-equivalent exposures or estimates of clinical antispastic efficacy. The source shares the original research program and is not an independent phenotype-validation cohort.

The model retains 13,761 genes after expression and effective-length eligibility filtering. Gene TPMs were multiplied by the across-sample mean effective gene length and rescaled to the original library expected-count totals. The deposited estimates were not rounded. Primary robust edgeR quasi-likelihood tests and secondary robust limma-voom tests use the same 24-condition design and retained genes. BH adjustment is reported both within each contrast and across all 165,132 gene-by-contrast tests for each method. Pooled BH is a different family definition and need not be numerically more conservative within every individual contrast.

## Complete primary comparison table

Counts below are discoveries at BH q <= 0.05. The shared column uses the across-12-contrast family for both methods; it is descriptive cross-method overlap, not an independently calibrated test.

| Solvent | PR:GR | Ordinal level | edgeR within-contrast | edgeR across 12 | voom within-contrast | voom across 12 | Shared across 12 |
|---|---|---|---:|---:|---:|---:|---:|
| Water | 2:1 | Low | 1 | 8 | 1 | 3 | 3 |
| Water | 2:1 | Middle | 0 | 12 | 0 | 3 | 3 |
| Water | 2:1 | High | 6 | 22 | 6 | 6 | 6 |
| Water | 1:2 | Low | 0 | 1 | 0 | 0 | 0 |
| Water | 1:2 | Middle | 0 | 1 | 0 | 0 | 0 |
| Water | 1:2 | High | 0 | 2 | 0 | 0 | 0 |
| Ethanol | 2:1 | Low | 0 | 7 | 0 | 3 | 3 |
| Ethanol | 2:1 | Middle | 1 | 25 | 0 | 2 | 2 |
| Ethanol | 2:1 | High | 2,389 | 589 | 292 | 17 | 16 |
| Ethanol | 1:2 | Low | 1 | 7 | 1 | 1 | 1 |
| Ethanol | 1:2 | Middle | 0 | 4 | 0 | 0 | 0 |
| Ethanol | 1:2 | High | 189 | 124 | 207 | 64 | 64 |

Gene-effect ranks agree closely between methods across the 12 comparisons (Spearman rho 0.9637–0.9860), while significance membership is method-dependent. The largest discovery counts occur in the high-level ethanol conditions. The high-level ethanol 2:1 comparison requires particular caution because one CEM library is a pronounced expression outlier; this finding triggered the influence analysis below. Few discoveries at other conditions do not establish equivalence between preparations.

## Source-label sensitivity

GSM9013914 has an isolated CEM label in Table S3, while its JGE06 code and GEO operation order identify IEM. The full analysis uses the corroborated code/GEO interpretation and preserves the conflicting original label. Removing that library, retaining the full gene universe and imported length scaling, changes the affected low-level ethanol 2:1 edgeR discoveries across all 12 contrasts from 7 to 5. All five shared discoveries retain direction. Whole-profile effect-rank correlation is 0.8616; thus local quantitative effects are sensitive even though the shared discovery directions agree. Other comparison effect-rank correlations exceed 0.99999 except rounding at the fifth decimal; changes in shared dispersion and pooled BH can still alter their threshold membership.

## Expression QC and the high-level ethanol 2:1 comparison

The median within-condition mean log-expression correlation is 0.9899. The minimum pairwise correlation is 0.8885 in JGE02 high. GSM9013891 (JGE02 high, replicate 3) is extreme on PC1. Its filtered adjusted count total is 5.39 million, compared with 25.33 and 27.59 million for the other two CEM replicates. This establishes expression/depth heterogeneity, not its technical or biological cause. The sample remains in the primary result.

The post-fit influence analysis omitted each of all six libraries from the affected CEM/IEM comparison once, preserving the full retained gene set, 24-condition specification and complete 12-contrast BH family. This avoids selecting only the visually unusual sample for a favorable result. No replacement primary analysis was selected.

| Omitted library | Preparation | edgeR discoveries across 12 | voom discoveries across 12 | Shared, same direction | edgeR effect rho versus full |
|---|---|---:|---:|---:|---:|
| GSM9013879 | CEM replicate 1 | 2,362 | 1,004 | 927 | 0.9181 |
| GSM9013885 | CEM replicate 2 | 3,018 | 1,812 | 1,715 | 0.9433 |
| GSM9013891 | CEM replicate 3; PCA outlier | 27 | 24 | 23 | 0.5036 |
| GSM9013906 | IEM replicate 1 | 57 | 7 | 6 | 0.9671 |
| GSM9013912 | IEM replicate 2 | 51 | 2 | 2 | 0.9619 |
| GSM9013918 | IEM replicate 3 | 61 | 1 | 1 | 0.9619 |

Removing the PCA outlier reduces the primary edgeR discovery count from 589 to 27, retaining only nine of the original discoveries; those nine retain direction. The median absolute gene log2 fold change falls from 0.1832 to 0.0831. Omitting either of the other CEM replicates instead increases the contribution of the outlying library and produces thousands of discoveries. Omitting an IEM library also changes threshold membership substantially. Thus the apparent large high-level ethanol 2:1 response is strongly dependent on the small-sample composition and is not accepted as a stable preparation-wide effect. These checks do not identify the cause of the unusual library or justify deleting it.

## Manuscript-ready results and interpretation

An additional within-accession analysis compared co-extracted and post-extraction-mixed preparations across 72 PC12 libraries spanning two solvents, two PR:GR ratios and three reported ordinal levels. Twelve prespecified contrasts were evaluated over 13,761 genes. Gene-effect rankings were concordant between edgeR and limma-voom (Spearman rho, 0.964–0.986), but discovery counts depended strongly on the statistical method and sample composition. At the high ethanol level, the 2:1 contrast yielded 589 edgeR and 17 voom discoveries after BH adjustment across 165,132 gene-by-contrast tests, with 16 shared; the corresponding 1:2 contrast yielded 124 and 64, with 64 shared. One high-level 2:1 co-extract library was a pronounced expression outlier. Systematic omission of each of the six libraries in that comparison produced 27–3,018 edgeR discoveries. Omitting the outlying library reduced effect-rank agreement with the full analysis to 0.504 and retained nine of its 589 discoveries. These findings support explicit preparation- and library-dependent interpretation of the deposited perturbations. They do not establish formulation equivalence, a reproducible large extraction effect across preparations, or an antispastic mechanism. Absolute-dose discrepancies remain unresolved, and the canonical 1:1 co-extracted preparation is absent from this accession.

## Reproducible outputs and remaining work

- [Fixed analysis specification](../../analysis/preparation_extension/GSE298414_WITHIN_SOURCE_PREPARATION_ANALYSIS_20261025.md).
- [Source and label reconciliation](Route_B2_GSE298414_preparation_extension_20261025.md).
- Full matrices, models, all gene-level tests and QC: `data/preparation_extension/GSE298414/`.
- Saved-result checks: `review/review_receipt.json`; the checks verify complete-family BH values, saved voom coefficients, design rank, contrast direction and the exact omitted sample without fitting another model.
- Six-library influence results: `high_ethanol_six_library_influence/summary.tsv`, complete-family statistics and `review_receipt.json`. All six omissions and both methods' complete-family BH adjustments were checked.
- Plot: `review/preparation_QC_and_label_influence.pdf` and `.png`.

The earlier import attempts are retained and explicitly marked not accepted. The successful import both verifies the length-scaling formula and confirms that it changes raw counts. The deposited zero-effective-length assignments are retained in a separate table and excluded from inference through a structural availability filter.

Program-level interpretation, solvent/ratio/single-herb extensions, an exposure bridge and independent PSS validation remain open. The main article remains v38 pending scientific integration. This completed finite extraction-method analysis advances the preparation-sensitivity component; it does not complete the full preparation robustness requirement.
