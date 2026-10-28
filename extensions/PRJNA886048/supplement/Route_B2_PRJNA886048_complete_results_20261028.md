# Complete injured-spinal pharmacological context: PRJNA886048

28 October 2026. Results below were computed from all 12 accepted libraries using the fixed unpaired baclofen-minus-saline contrast (six animals per group). This is retrospective pharmacological context in injured spinal tissue, not outcome-protected post-stroke spasticity validation. Candidate admission and S1 remain unchanged.

## Complete data and quality checks

All 24 original input files retain accepted source MD5 records. Each complete quantification was returned and accepted with archive and output SHA256, paired-read accounting, and 146,190 finite nonnegative transcript targets. The last two exact archives were delivered by the user; no additional assistant download was started. All 12 libraries were retained. Retained pairs ranged from 31,136,937 to 50,337,119 per library; pseudoalignment ranged from 84.22% to 91.52%, and post-filtering Q30 from 94.56% to 95.28%. The descriptive MDS plot shows within-group dispersion and overlapping treatment distributions; it was not used to remove samples.

## Gene and program results

Of 17,013 genes passing the fixed expression filter, none met BH FDR < 0.05 in either the primary robust edgeR quasi-likelihood analysis or the limma-voom sensitivity. Gene effect estimates were strongly rank-correlated between methods (Spearman rho = 0.9698). These findings do not establish absence of a drug response: they delimit gene-level evidence under this cohort, analysis and multiplicity family.

Competitive CAMERA tests identified 65 of 1,588 tested Reactome terms and four of 91 DoRothEA membership sets at collection-wise BH FDR < 0.05; none of the 14 PROGENy membership sets passed. The strongest Reactome membership signals included overlapping translation-related sets. These terms are correlated and must not be counted as independent mechanisms. Membership enrichment summarizes relative gene-level shifts and is distinct from signed footprint activity. None of 14 PROGENy or 91 DoRothEA signed score contrasts passed collection-wise BH FDR < 0.05. Complete effects, intervals, coverage and native q values remain available, including nonsignificant results.

## Prespecified cross-source comparison

All 14 PROGENy programs were compared with all 19 retained context/specification columns (266 cells), including nine core contexts and ten companions. Every comparison had 14 finite pairs, with no exact-zero or missing pairs. The table reports all columns without significance filtering. Rho describes ordering across the complete program panel; no agreement P value, pooled effect or clinical efficacy grade was computed. Sources retain their own effect scales and uncertainty definitions. Companion columns reuse data or represent different stages/doses and are not independent replications.

| Comparator | Role | Same sign | Opposite sign | Descriptive rho |
|---|---|---:|---:|---:|
| G16710_primary26_interaction | core_context | 7/14 | 7/14 | 0.147 |
| G128_raw_primary_all30_chronic | core_context | 7/14 | 7/14 | -0.064 |
| G128_raw_complete29_sensitivity_chronic | same_data_sensitivity | 7/14 | 7/14 | -0.213 |
| G128_adapter_primary_all30_chronic | same_data_sensitivity | 7/14 | 7/14 | -0.064 |
| G128_adapter_complete29_sensitivity_chronic | same_data_sensitivity | 7/14 | 7/14 | -0.174 |
| G277_v2a_chronic_combined | core_context | 8/14 | 6/14 | 0.222 |
| G277_v2a_early_combined | stage_context | 11/14 | 3/14 | 0.336 |
| G277_astrocyte_chronic_combined | core_context | 5/14 | 9/14 | -0.349 |
| G277_astrocyte_early_combined | stage_context | 8/14 | 6/14 | -0.169 |
| G277_microglia_chronic_combined | core_context | 11/14 | 3/14 | 0.244 |
| G277_microglia_early_combined | stage_context | 10/14 | 4/14 | 0.367 |
| JGW_PC12_high | core_context | 7/14 | 7/14 | -0.090 |
| JGW_PC12_trend | dose_context | 10/14 | 4/14 | -0.160 |
| JGW_C2C12_high | core_context | 6/14 | 8/14 | -0.178 |
| JGW_C2C12_trend | dose_context | 6/14 | 8/14 | -0.178 |
| JGW_HepG2_high | core_context | 9/14 | 5/14 | 0.358 |
| JGW_HepG2_trend | dose_context | 10/14 | 4/14 | 0.297 |
| G205_canonical_local_disease | injury_context | 7/14 | 7/14 | -0.143 |
| G205_canonical_intervention | core_context | 8/14 | 6/14 | -0.007 |

The primary phenotype-anchored GSE16710 and chronic bulk GSE128623 comparisons each had seven same-sign and seven opposite-sign programs. Cell-preparation and intervention comparisons varied in direction and ordering. This panel therefore does not establish a consistent, general reversal of the chronic post-stroke program by baclofen. It provides a complete, bounded pharmacological comparison to accompany the full-formula and injury results; it does not establish formulation equivalence, antispastic efficacy or a cell-intrinsic mechanism.

## Execution and reproducibility

The frozen driver completed and saved the gene models before stopping at a plotting interface incompatibility: edgeR 4.10.1 does not export plotMDS. A separately retained executed derivative changes only that call to the limma 3.68.4 generic, retaining the DGEList method. Input identity, gene import and fitted gene-model checksums were identical before and after continuation; no gene model was refitted. The original driver, failed-run log, executed derivative and continuation receipt are retained. A separate display-only MDS render expands margins to show all labels and adds group labels; the original plot is also retained. Statistical parameters and sample inclusion were unchanged.

Results: `R:/Shaoyao_storage/results/PRJNA886048_complete_cohort_20261028_01/`.

Complete descriptive panel: `R:/Shaoyao_storage/results/PRJNA886048_cross_source_comparison_20261028_01/`.

Execution and saved-content verification: `route_b2/results/diagnostics/PRJNA886048_COMPLETE_COHORT_EXECUTION_20261028.json`.

This supplement records verified new results. The verified results are integrated into the reviewed v46 main manuscript. Incorporation into the next public reproducibility release remains pending.
