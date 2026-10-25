# GSE298414 stratified preparation and program results

25 October 2026. Completed new contrasts and program comparisons, with independent readback verification. All earlier accepted gene results are retained.

## Scope and data reuse

The existing 72-library, 24-condition model now supplies 36 stratified contrasts: 12 extraction-method comparisons, 12 ethanol-minus-water comparisons, and 12 PR:GR 2:1-minus-1:2 comparisons. Each contrast compares three libraries per condition at the same reported ordinal level. The gene filter contains 13,761 genes. The new comparisons reuse the existing count import, normalization, dispersion, QL model and voom weights. Reconstructing the base weighted linear coefficients for the additional contrasts reproduced the old 12 voom coefficient and P-value matrices. No canonical JGW03 model was rerun.

These contrasts share samples. A contrast count is not an independent replication count. Between-stratum differences in statistical significance are not interaction tests. The source still lacks a within-accession 1:1 CEM preparation; absolute dose labels remain unresolved.

## Additional gene comparisons

BH values are reported within each contrast and across the gene-by-12-contrast family for each contrast class. The existing extraction family was not changed by the new solvent and ratio comparisons.

| Contrast | edgeR within q <= 0.05 | edgeR class q <= 0.05 | voom within q <= 0.05 | voom class q <= 0.05 | Shared class discoveries | Effect rho |
|---|---:|---:|---:|---:|---:|---:|
| JGE02_L1_minus_JGW02_L1 | 10 | 41 | 16 | 33 | 32 | 0.9813 |
| JGE02_L2_minus_JGW02_L2 | 4 | 50 | 10 | 46 | 46 | 0.9838 |
| JGE02_L3_minus_JGW02_L3 | 3012 | 1416 | 1541 | 454 | 441 | 0.9615 |
| JGE04_L1_minus_JGW04_L1 | 23 | 73 | 28 | 59 | 57 | 0.9865 |
| JGE04_L2_minus_JGW04_L2 | 41 | 104 | 42 | 77 | 77 | 0.9839 |
| JGE04_L3_minus_JGW04_L3 | 484 | 379 | 536 | 322 | 319 | 0.9845 |
| JGE06_L1_minus_JGW06_L1 | 14 | 59 | 14 | 44 | 44 | 0.9845 |
| JGE06_L2_minus_JGW06_L2 | 50 | 150 | 47 | 105 | 105 | 0.9865 |
| JGE06_L3_minus_JGW06_L3 | 630 | 445 | 711 | 384 | 376 | 0.9886 |
| JGE08_L1_minus_JGW08_L1 | 14 | 76 | 13 | 61 | 61 | 0.9853 |
| JGE08_L2_minus_JGW08_L2 | 5 | 66 | 4 | 47 | 47 | 0.9843 |
| JGE08_L3_minus_JGW08_L3 | 11 | 77 | 8 | 63 | 62 | 0.9849 |
| JGW02_L1_minus_JGW04_L1 | 0 | 3 | 0 | 3 | 2 | 0.9805 |
| JGW02_L2_minus_JGW04_L2 | 0 | 4 | 0 | 1 | 1 | 0.9831 |
| JGW02_L3_minus_JGW04_L3 | 0 | 5 | 0 | 2 | 2 | 0.9811 |
| JGW06_L1_minus_JGW08_L1 | 0 | 1 | 0 | 0 | 0 | 0.9794 |
| JGW06_L2_minus_JGW08_L2 | 0 | 2 | 0 | 0 | 0 | 0.9799 |
| JGW06_L3_minus_JGW08_L3 | 0 | 17 | 0 | 6 | 6 | 0.9844 |
| JGE02_L1_minus_JGE04_L1 | 0 | 4 | 0 | 2 | 2 | 0.9776 |
| JGE02_L2_minus_JGE04_L2 | 0 | 9 | 0 | 3 | 3 | 0.9821 |
| JGE02_L3_minus_JGE04_L3 | 2885 | 978 | 1335 | 219 | 217 | 0.9525 |
| JGE06_L1_minus_JGE08_L1 | 0 | 16 | 0 | 4 | 4 | 0.9871 |
| JGE06_L2_minus_JGE08_L2 | 2 | 27 | 0 | 11 | 11 | 0.9857 |
| JGE06_L3_minus_JGE08_L3 | 190 | 148 | 202 | 105 | 105 | 0.9861 |

The largest solvent and ratio comparisons both include JGE02 at the high level, which contains the previously identified expression outlier GSM9013891. Their discovery counts must not be used as stable whole-preparation effects without corresponding influence qualification. Other preparation contrasts also show differences, including high-level JGE06 versus JGW06, but these remain source-specific observations rather than clinical efficacy evidence.

## Resource coverage and program testing

Human ortholog mapping and label-independent representative-gene selection retained 11,356 background symbols; native rat mapping retained 11,870. Selection used highest mean TMM logCPM across all 72 libraries, with gene-ID tie breaking, matching the canonical mapping rule. Mapping was fixed for score influence analyses.

| Collection | Reference programs | Eligible programs | Test or score |
|---|---:|---:|---|
| Reactome97 | 2,836 | 1,469 | CAMERA; 10–500 expressed members |
| SynGO sensitivity | 249 | 69 | CAMERA; 10–500 expressed members |
| PROGENy | 14 | 14 | Signed ULM; at least 50 expressed targets |
| DoRothEA rat A/B | 117 | 107 | Signed ULM; at least five expressed targets |

CAMERA uses the saved human-mapped voom object and inter-gene correlation 0.01. It tests relative membership shifts; its direction is not a signed activation estimate. ULM uses frozen signed weights and the full mapped logCPM background. Library-level score contrasts use the same condition design and robust limma empirical Bayes with 95% intervals. ULM gene-regression P values are not used as treatment-effect P values. Collection-specific within-contrast and collection-by-contrast-class BH families are retained separately.

The following are **program–contrast findings**, not numbers of distinct programs or independent observations. Each row aggregates 12 related contrasts.

| Contrast class | Reactome findings | SynGO findings | PROGENy findings | DoRothEA findings |
|---|---:|---:|---:|---:|
| extraction | 514 | 25 | 1 | 10 |
| solvent | 445 | 11 | 7 | 44 |
| ratio | 1048 | 108 | 1 | 3 |

## Influential library and signed programs

In the high-level ethanol 2:1 CEM-minus-IEM comparison, three DoRothEA scores pass the collection-by-extraction-class BH threshold in the full analysis. All three lose that threshold after GSM9013891 is omitted. The table reports differences in inferred ULM score units, not expression fold changes or measured TF activity. Intervals are pointwise moderated 95% intervals and are not multiplicity-adjusted.

| TF | Full effect and 95% CI | Full class q | Effect after omission and 95% CI | Class q after omission |
|---|---|---:|---|---:|
| Foxo1 | 0.377 (0.203, 0.552) | 0.0134 | 0.093 (-0.002, 0.188) | 0.4001 |
| Hif1a | -0.381 (-0.551, -0.210) | 0.0121 | -0.191 (-0.343, -0.039) | 0.2565 |
| Smad3 | 0.409 (0.232, 0.586) | 0.0121 | 0.186 (0.041, 0.330) | 0.2565 |

All six libraries in that high-level comparison were omitted once; significant DoRothEA scores ranged from zero to 18. The separate disputed low-level label omission was also retained. These are fixed-score influence checks: membership, mapping, normalization and full-sample scores were held fixed, and only the library-level score model was refitted. They do not establish complete CAMERA or raw-to-result robustness. Hif1a and Smad3 retain nominal intervals excluding zero after the outlier omission, but neither retains the defined multiplicity-adjusted threshold. The threshold change should not be described as a reversal of direction or proof of no effect.

For high-level ethanol 1:2 CEM minus IEM, the inferred MAPK score difference is 0.570 (95% CI 0.377–0.763; class q approximately 3.1e-5). This is a preparation-context observation; it has not undergone a full influence analysis of its own comparison libraries and is not admitted as a disease-reversal candidate.

## Scientific interpretation

Preparation identity changes the deposited transcriptional and inferred-program responses, but the magnitude and statistical support are not uniformly stable across methods and libraries. The most conspicuous high-level 2:1 response and its three TF findings are substantially influenced by one library. These results strengthen the need to retain the exact recipe, solvent, extraction sequence and ordinal exposure labels. They do not connect an inferred program to antispastic efficacy or establish equivalence among preparations. S1 still admits zero programs; no downstream candidate gate is reopened.

## Reproducibility and next work

- [Analysis specification](../../analysis/preparation_extension/GSE298414_PROGRAM_AND_STRATIFIED_EXTENSION_20261025.md).
- [Initial extraction-method analysis and six-library gene influence](Route_B2_GSE298414_preparation_results_20261025.md).
- All outputs, coverage, complete test tables, fixed scores and readback receipt: `data/preparation_extension/GSE298414/program_and_stratified_extension_01`.
- Independent verification recomputed BH values and checked all 34,848 full/omission score effects directly from library means. The 144 program-family summaries were checked against complete tables; no model was fitted for these checks.

Single-herb comparisons, missing 1:1 CEM cells, the exposure bridge, independent PSS validation and the separate spinal pharmacology acquisition remain open. Whole-study completion is not claimed.
