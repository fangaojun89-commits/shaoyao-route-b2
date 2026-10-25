# Mixture and single herb transcriptional comparisons

25 October 2026. GSE298414 supplies a finite preparation comparison: 126 PC12 libraries, 42 conditions and 60 new mixture-minus-single-herb contrasts. Differences are substantial in several high-level ethanol conditions, including the separately extracted 1:1 mixture. The high-level ethanol 2:1 co-extract is also strongly influenced by GSM9013891. This evidence supports preparation-specific transcriptional differences. Because the comparisons do not equalize component exposure or identify an additive reference response, they do not estimate synergy, antagonism or clinical antispastic efficacy.

## Scope and methods

The fixed selection includes 36 single-herb libraries, 36 CEM libraries (2:1 and 1:2) and 54 IEM libraries (2:1, 1:1 and 1:2). All selected records carry batch A and three source-reported biological replicates per condition. Seventy-two compressed inputs were reused and 54 were unpacked from the existing verified official archive. No raw reads were downloaded for this extension. The previous 72-library/36-contrast analysis and all canonical models remain unchanged. The 72- and 126-library sets overlap and are not independent replications.

The new comparisons match solvent and source-reported ordinal level. The unresolved 20/100/500 versus 4/20/100 µg/mL discrepancy is preserved. The accession still lacks a 1:1 CEM condition. Source-label discrepancies for GSM9013914 and GSM9013917 are documented in the retained ledger; preparation codes and GEO protocols determine IEM and single licorice classification, respectively.

Gene-level RSEM TPM values were multiplied by mean effective length across the 126 libraries and rescaled to each original expected-count total (explicit tximport lengthScaledTPM). Positive effective lengths in every selected library and design-aware filterByExpr retained 13,898 genes. The categorical 42-condition model used TMM normalization, robust edgeR quasi-likelihood and robust limma-voom, with nominal voom 95% intervals. Residual degrees of freedom are 84. BH adjustments are computed within each contrast and separately across the 30 mixture-minus-peony and 30 mixture-minus-licorice contrasts, for each method. Class-wide BH can yield more discoveries than an individual contrast's BH; the procedures address different test families.

The planned GSM9013891 omission uses 125 libraries, the same filtered genes and the original 126-library length-scaled values, with recalculated TMM and gene models (residual df 83). It does not replace the primary results. Program tests reuse the saved primary voom fit and label-independent representative-symbol mapping. CAMERA uses an inter-gene correlation of 0.01 and 10–500 mapped expressed members. PROGENy and rat DoRothEA scores use ULM with center=FALSE and minimum 50 and five nonzero expressed targets, respectively; robust limma tests library-level score contrasts. Program multiplicity is separate by resource and comparator class. Three score omissions hold the primary mapping and scores fixed and omit GSM9013891, GSM9013914 or GSM9013917 individually; these are not CAMERA or raw-to-result omission analyses.

## Gene findings

The edgeR/voom effect-rank correlations across the 60 comparisons range from 0.952482 to 0.990224. At the high ordinal level in 70% ethanol, the 1:1 IEM-versus-peony comparison yields 1,885 edgeR and 1,964 voom class-adjusted findings, with 1,831 common genes; IEM-versus-licorice yields 1,406 and 1,410, with 1,330 common genes. Median absolute edgeR log2 fold changes are 0.177 and 0.172. These observations describe this source-specific comparison, not a superiority or synergy test. Library influence within those 1:1 comparison groups has not been exhaustively tested.

For high-level ethanol 2:1 CEM versus peony, omitting GSM9013891 changes edgeR/voom class-adjusted discovery counts from 2,215/1,373 to 604/630. Against licorice, counts change from 1,547/739 to 394/377. Full-versus-omission gene-effect rank correlations are 0.624441 and 0.563982. Of the original edgeR discoveries, 362 and 210 remain adjusted findings, including 359 and 205 with the same direction. Median absolute log2 fold changes decrease from 0.199 to 0.107 and from 0.182 to 0.107. The other 58 effect-rank correlations exceed 0.9997, although their discovery counts can change because dispersion, precision and multiplicity are estimated jointly.

The following complete table reports both primary methods for all 60 comparisons. Counts are gene–contrast findings, not independent genes or independent replications. [Full and omission summaries](Route_B2_GSE298414_single_herb_gene_summary_20261025.tsv) and [complete influence summary](Route_B2_GSE298414_single_herb_influence_20261025.tsv) retain both outcomes.

| contrast_id | edgeR_class05 | voom_class05 | shared_class05 | effect_spearman | median_abs_edgeR_logFC |
| --- | --- | --- | --- | --- | --- |
| JGW02_L1_minus_JGW01_L1 | 6 | 4 | 4 | 0.984566 | 0.0696267 |
| JGW02_L1_minus_JGW05_L1 | 28 | 24 | 23 | 0.983589 | 0.0784301 |
| JGW02_L2_minus_JGW01_L2 | 3 | 4 | 3 | 0.981348 | 0.0573166 |
| JGW02_L2_minus_JGW05_L2 | 11 | 8 | 7 | 0.984123 | 0.0761217 |
| JGW02_L3_minus_JGW01_L3 | 6 | 7 | 5 | 0.978424 | 0.0587035 |
| JGW02_L3_minus_JGW05_L3 | 64 | 59 | 57 | 0.986378 | 0.0930441 |
| JGW04_L1_minus_JGW01_L1 | 11 | 11 | 9 | 0.98092 | 0.0721454 |
| JGW04_L1_minus_JGW05_L1 | 10 | 5 | 3 | 0.980122 | 0.0665121 |
| JGW04_L2_minus_JGW01_L2 | 23 | 21 | 19 | 0.984926 | 0.0868184 |
| JGW04_L2_minus_JGW05_L2 | 6 | 4 | 4 | 0.980846 | 0.0752596 |
| JGW04_L3_minus_JGW01_L3 | 22 | 20 | 18 | 0.980282 | 0.0824377 |
| JGW04_L3_minus_JGW05_L3 | 8 | 7 | 5 | 0.977921 | 0.0747148 |
| JGW06_L1_minus_JGW01_L1 | 17 | 20 | 17 | 0.979888 | 0.0714956 |
| JGW06_L1_minus_JGW05_L1 | 1 | 4 | 0 | 0.976157 | 0.0521627 |
| JGW06_L2_minus_JGW01_L2 | 25 | 22 | 19 | 0.979831 | 0.0794673 |
| JGW06_L2_minus_JGW05_L2 | 6 | 2 | 2 | 0.980768 | 0.0686571 |
| JGW06_L3_minus_JGW01_L3 | 33 | 33 | 31 | 0.98376 | 0.0787034 |
| JGW06_L3_minus_JGW05_L3 | 23 | 23 | 21 | 0.985297 | 0.0868164 |
| JGW07_L1_minus_JGW01_L1 | 33 | 26 | 26 | 0.979548 | 0.0805186 |
| JGW07_L1_minus_JGW05_L1 | 5 | 1 | 0 | 0.975863 | 0.057127 |
| JGW07_L2_minus_JGW01_L2 | 43 | 39 | 36 | 0.980492 | 0.083399 |
| JGW07_L2_minus_JGW05_L2 | 15 | 5 | 4 | 0.978737 | 0.0614416 |
| JGW07_L3_minus_JGW01_L3 | 83 | 84 | 77 | 0.98504 | 0.0915172 |
| JGW07_L3_minus_JGW05_L3 | 5 | 1 | 1 | 0.983902 | 0.0602475 |
| JGW08_L1_minus_JGW01_L1 | 34 | 29 | 29 | 0.979699 | 0.0845321 |
| JGW08_L1_minus_JGW05_L1 | 1 | 2 | 1 | 0.980124 | 0.0648394 |
| JGW08_L2_minus_JGW01_L2 | 50 | 48 | 45 | 0.982414 | 0.0945893 |
| JGW08_L2_minus_JGW05_L2 | 10 | 3 | 3 | 0.981226 | 0.0735735 |
| JGW08_L3_minus_JGW01_L3 | 113 | 114 | 109 | 0.983795 | 0.098141 |
| JGW08_L3_minus_JGW05_L3 | 11 | 4 | 4 | 0.982624 | 0.0723385 |
| JGE02_L1_minus_JGE01_L1 | 6 | 4 | 1 | 0.977277 | 0.0543266 |
| JGE02_L1_minus_JGE05_L1 | 7 | 6 | 5 | 0.98222 | 0.0720038 |
| JGE02_L2_minus_JGE01_L2 | 4 | 3 | 1 | 0.978618 | 0.0735937 |
| JGE02_L2_minus_JGE05_L2 | 35 | 26 | 25 | 0.980432 | 0.0860355 |
| JGE02_L3_minus_JGE01_L3 | 2215 | 1373 | 1270 | 0.965295 | 0.198867 |
| JGE02_L3_minus_JGE05_L3 | 1547 | 739 | 685 | 0.952482 | 0.182331 |
| JGE04_L1_minus_JGE01_L1 | 26 | 26 | 22 | 0.979234 | 0.0756998 |
| JGE04_L1_minus_JGE05_L1 | 1 | 1 | 1 | 0.977958 | 0.0532766 |
| JGE04_L2_minus_JGE01_L2 | 66 | 66 | 61 | 0.981403 | 0.0797527 |
| JGE04_L2_minus_JGE05_L2 | 9 | 7 | 7 | 0.979873 | 0.0660896 |
| JGE04_L3_minus_JGE01_L3 | 1389 | 1395 | 1346 | 0.984769 | 0.144883 |
| JGE04_L3_minus_JGE05_L3 | 1 | 1 | 1 | 0.969385 | 0.0569583 |
| JGE06_L1_minus_JGE01_L1 | 33 | 34 | 28 | 0.97823 | 0.0715691 |
| JGE06_L1_minus_JGE05_L1 | 3 | 2 | 2 | 0.979408 | 0.0599494 |
| JGE06_L2_minus_JGE01_L2 | 21 | 17 | 16 | 0.97983 | 0.0673017 |
| JGE06_L2_minus_JGE05_L2 | 12 | 6 | 6 | 0.980845 | 0.0680222 |
| JGE06_L3_minus_JGE01_L3 | 219 | 210 | 201 | 0.981757 | 0.0861586 |
| JGE06_L3_minus_JGE05_L3 | 452 | 425 | 418 | 0.983522 | 0.108812 |
| JGE07_L1_minus_JGE01_L1 | 42 | 44 | 39 | 0.977925 | 0.0778284 |
| JGE07_L1_minus_JGE05_L1 | 7 | 4 | 3 | 0.982404 | 0.0760432 |
| JGE07_L2_minus_JGE01_L2 | 96 | 90 | 84 | 0.979727 | 0.0876263 |
| JGE07_L2_minus_JGE05_L2 | 26 | 18 | 17 | 0.984337 | 0.079252 |
| JGE07_L3_minus_JGE01_L3 | 1885 | 1964 | 1831 | 0.990224 | 0.176908 |
| JGE07_L3_minus_JGE05_L3 | 1406 | 1410 | 1330 | 0.988847 | 0.171788 |
| JGE08_L1_minus_JGE01_L1 | 39 | 39 | 33 | 0.980932 | 0.0774307 |
| JGE08_L1_minus_JGE05_L1 | 4 | 3 | 3 | 0.982961 | 0.0721416 |
| JGE08_L2_minus_JGE01_L2 | 42 | 36 | 34 | 0.980994 | 0.0777574 |
| JGE08_L2_minus_JGE05_L2 | 11 | 8 | 8 | 0.983688 | 0.0696273 |
| JGE08_L3_minus_JGE01_L3 | 1113 | 1114 | 1076 | 0.986553 | 0.129799 |
| JGE08_L3_minus_JGE05_L3 | 340 | 316 | 310 | 0.983457 | 0.107081 |

## Program coverage and all family findings

| collection | reference_programs | eligible_programs |
| --- | --- | --- |
| DoRothEA | 117 | 107 |
| PROGENy | 14 | 14 |
| Reactome97 | 2836 | 1475 |
| SynGO_sensitivity | 249 | 69 |

There are 240 program families (60 contrasts × four collections), 92,640 CAMERA tests and 7,260 primary signed-score tests. Counts below are program–contrast findings, not unique pathways or independent replications. [The complete family summary](Route_B2_GSE298414_single_herb_program_summary_20261025.tsv) includes null families.

| collection | contrast_class | tested_program_contrasts | class_FDR05 |
| --- | --- | --- | --- |
| DoRothEA | mixture_minus_GR | 3210 | 34 |
| DoRothEA | mixture_minus_PR | 3210 | 81 |
| PROGENy | mixture_minus_GR | 420 | 6 |
| PROGENy | mixture_minus_PR | 420 | 18 |
| Reactome97 | mixture_minus_GR | 44250 | 1741 |
| Reactome97 | mixture_minus_PR | 44250 | 2372 |
| SynGO_sensitivity | mixture_minus_GR | 2070 | 106 |
| SynGO_sensitivity | mixture_minus_PR | 2070 | 280 |

The complete 14-pathway PROGENy results for both high-level ethanol 1:1 comparisons are displayed as descriptive examples, chosen after the gene results. All other contrasts remain in the retained output. These are inferred transcriptional footprint scores rather than direct measurements of pathway activity.

| contrast_id | term_id | estimate | CI_low | CI_high | FDR_class |
| --- | --- | --- | --- | --- | --- |
| JGE07_L3_minus_JGE01_L3 | Androgen | -0.248414 | -0.499495 | 0.00266755 | 0.439742 |
| JGE07_L3_minus_JGE01_L3 | EGFR | -0.130654 | -0.488237 | 0.22693 | 0.916218 |
| JGE07_L3_minus_JGE01_L3 | Estrogen | -0.405614 | -0.642369 | -0.168858 | 0.0317733 |
| JGE07_L3_minus_JGE01_L3 | Hypoxia | -0.486709 | -0.761822 | -0.211595 | 0.0238879 |
| JGE07_L3_minus_JGE01_L3 | JAK-STAT | 0.067145 | -0.180669 | 0.314959 | 0.948007 |
| JGE07_L3_minus_JGE01_L3 | MAPK | -0.350878 | -0.512479 | -0.189277 | 0.00418084 |
| JGE07_L3_minus_JGE01_L3 | NFkB | -0.061368 | -0.244594 | 0.121858 | 0.932072 |
| JGE07_L3_minus_JGE01_L3 | PI3K | -0.47139 | -0.715367 | -0.227413 | 0.0121499 |
| JGE07_L3_minus_JGE01_L3 | TGFb | 0.0746401 | -0.167335 | 0.316615 | 0.945532 |
| JGE07_L3_minus_JGE01_L3 | TNFa | 0.0866562 | -0.110291 | 0.283603 | 0.916218 |
| JGE07_L3_minus_JGE01_L3 | Trail | 0.151049 | -0.0247561 | 0.326854 | 0.618582 |
| JGE07_L3_minus_JGE01_L3 | VEGF | 0.0339418 | -0.155697 | 0.223581 | 0.978196 |
| JGE07_L3_minus_JGE01_L3 | WNT | 0.10294 | -0.0620397 | 0.26792 | 0.901293 |
| JGE07_L3_minus_JGE01_L3 | p53 | 0.475785 | 0.248567 | 0.703004 | 0.004947 |
| JGE07_L3_minus_JGE05_L3 | Androgen | -0.275538 | -0.526619 | -0.0244566 | 0.371402 |
| JGE07_L3_minus_JGE05_L3 | EGFR | -0.222697 | -0.58028 | 0.134887 | 0.829572 |
| JGE07_L3_minus_JGE05_L3 | Estrogen | -0.376068 | -0.612824 | -0.139313 | 0.123304 |
| JGE07_L3_minus_JGE05_L3 | Hypoxia | 0.245109 | -0.0300048 | 0.520223 | 0.580373 |
| JGE07_L3_minus_JGE05_L3 | JAK-STAT | 0.199896 | -0.0479184 | 0.44771 | 0.665652 |
| JGE07_L3_minus_JGE05_L3 | MAPK | -0.459148 | -0.620749 | -0.297547 | 3.66874e-05 |
| JGE07_L3_minus_JGE05_L3 | NFkB | 0.0738861 | -0.10934 | 0.257113 | 0.901185 |
| JGE07_L3_minus_JGE05_L3 | PI3K | -0.506297 | -0.750274 | -0.26232 | 0.0113518 |
| JGE07_L3_minus_JGE05_L3 | TGFb | -0.266996 | -0.50897 | -0.0250208 | 0.371402 |
| JGE07_L3_minus_JGE05_L3 | TNFa | 0.146048 | -0.0508992 | 0.342994 | 0.728147 |
| JGE07_L3_minus_JGE05_L3 | Trail | 0.0212888 | -0.154516 | 0.197094 | 0.972589 |
| JGE07_L3_minus_JGE05_L3 | VEGF | 0.00749619 | -0.182143 | 0.197135 | 0.989018 |
| JGE07_L3_minus_JGE05_L3 | WNT | 0.232314 | 0.0673345 | 0.397294 | 0.219505 |
| JGE07_L3_minus_JGE05_L3 | p53 | 0.174742 | -0.052476 | 0.401961 | 0.719197 |

## Fixed score influence

| omission | collection | full_program_contrast_findings | retained_adjusted_findings | retained_same_direction |
| --- | --- | --- | --- | --- |
| omit_GSM9013891 | PROGENy | 24 | 17 | 17 |
| omit_GSM9013891 | DoRothEA | 115 | 88 | 88 |
| omit_GSM9013914 | PROGENy | 24 | 24 | 24 |
| omit_GSM9013914 | DoRothEA | 115 | 115 | 115 |
| omit_GSM9013917 | PROGENy | 24 | 23 | 23 |
| omit_GSM9013917 | DoRothEA | 115 | 115 | 115 |

Counts aggregate both comparator classes; their underlying BH families remain separate. Individual estimates, nominal intervals and adjusted values for every omission are retained. A stable sign alone is not persistence of adjusted evidence.

## Consequences for the study

Single-herb context is now evaluated for the available preparations. The comparison does not establish a classical decoction exposure bridge: the missing 1:1 CEM cell, conflicting absolute-dose descriptions and unmatched component concentrations remain unresolved. Shared cultures and comparators do not supply independent phenotype validation. S1 admission, canonical JGW03 and the candidate gate are unchanged.

## Reuse and verification

Independent retained-result checks reproduced within-contrast and class-wide BH adjustments, all gene summary counts, all 60 full-versus-omission effect comparisons, program coverage and all 240 program-family summaries. All 29,040 primary/omission signed effects were independently reconstructed from saved library-score means. These checks did not rerun scientific models and do not constitute clean-environment raw-data reconstruction.

Version 0.26.0 supplies all retained statistical tables, score matrices, coverage, sample and source identities, and the independent readback command in this repository. See [scope and instructions](RUN.md). Original model objects remain retained locally and are not part of the public payload. The earlier 72-library extension remains available separately.
