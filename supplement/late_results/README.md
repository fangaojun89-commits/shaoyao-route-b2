# Completed conditional references and matched quantification: supplementary results

25 October 2026. Accepted late-stage reporting results, accompanying the current English manuscript.

This public reading package brings together final late-stage results that previously required several dated addenda. It complements Supplementary Tables S1–S18 and Figures S1–S7; A1–A7 are package-local identifiers, not renumbered journal tables. It does not replace the original files or add a new analysis. All seven TSV files preserve their source bytes, precision and missing-value labels.

The study's fixed spasticity anchor admitted no programs. The analyses below describe source-wide conditional references and same-library technical sensitivity; downstream candidate-restricted links remained unestimated. GSE128623 empirical BH thresholds, JGW03 unadjusted conditional P thresholds, and GSE205906 native within-family BH thresholds have different meanings and must not be pooled.

## Tables and reading order

| ID | File | Data rows |
|---|---|---:|
| A1 | [A1_G128_all_families.tsv](../../data/derived/G128_FAMILY_SUMMARY.tsv) | 448 |
| A2 | [A2_G128_family_overview.tsv](tables/A2_G128_family_overview.tsv) | 4 |
| A3 | [A3_G128_gene_counts.tsv](tables/A3_G128_gene_counts.tsv) | 12 |
| A4 | [A4_G128_chronic_programs.tsv](tables/A4_G128_chronic_programs.tsv) | 24 |
| A5 | [A5_G128_day28_reference_spaces.tsv](tables/A5_G128_day28_reference_spaces.tsv) | 36 |
| A6 | [A6_JGW_Linux_all_families.tsv](../../data/derived/JGW_FAMILY_SUMMARY.tsv) | 123 |
| A7 | [A7_G205_matched_alignment.tsv](tables/A7_G205_matched_alignment.tsv) | 22 |

A1 is the complete GSE128623 family-level reference table. A2–A5 are compact source-preserved views of its main comparisons, checked against A1. A6 contains the complete JGW03 separate-Linux family summary. A7 contains the complete matched HISAT2/STAR comparison summary. Family-feature events can count the same feature repeatedly across methods and contrasts; they are not sample sizes or unique molecular entities.

## GSE128623: completed reference spaces

Retrospective GSE128623 day/sham conditional references enumerated all 252 assignments within each specification and day, except predefined29 D4/sham, which had 126; other groups remained fixed. Inclusive tails counted the observed assignment once, with separate empirical BH families for each method, collection and contrast. Activity strengths used absolute moderated t from limma score models. A subsequent descriptive comparison evaluated signed effects among common programs with finite, nonzero effects in both specifications. Separately standardized effect magnitudes were not pooled. Separately completed six-group references comprised the observed assignment and 10,000 prespecified sampled assignments per processing/cohort specification, with 67 method/collection/contrast families per specification (268 in total). Saved-tail arithmetic, diagnostics and metadata checks authenticated their reporting outputs without rerunning completed models. The six-group and day/sham references condition on different label spaces; neither authenticates biological exchangeability or true FDR control.

The four specifications retain their original processing and cohort definitions. In A1–A5, the threshold is empirical BH q≤0.05 within each stored method/collection/contrast family. Counts are point-estimate threshold crossings, not interval-robust findings or effect magnitudes. `native_fraction_allocations_any_q_le_05` describes the same-library allocation space and is not an estimate of true FDR. The day-specific and six-group reference spaces answer different conditional questions.

### Family overview (A2)

| cell_id | sampled_family_n | nonzero_empirical_q_family_n | zero_empirical_q_family_n | feature_family_discovery_events_NOT_UNIQUE_FEATURES | unavailable_feature_family_events |
|---|---|---|---|---|---|
| raw_all30 | 67 | 19 | 48 | 12099 | 0 |
| raw_predefined29 | 67 | 27 | 40 | 21240 | 0 |
| adapter_all30 | 67 | 18 | 49 | 12136 | 0 |
| adapter_predefined29 | 67 | 27 | 40 | 21550 | 0 |


### Gene counts (A3)

| cell_id | method | tested_genes | chronic_combined | D28_vs_control | D42_vs_control | reference_allocation_n |
|---|---|---|---|---|---|---|
| raw_all30 | edgeR | 17070 | 1230 | 2773 | 0 | 10001 |
| raw_all30 | voom | 17070 | 1305 | 2923 | 0 | 10001 |
| raw_all30 | DESeq2 | 17070 | 957 | 2397 | 0 | 10001 |
| raw_predefined29 | edgeR | 17268 | 2359 | 3955 | 0 | 10001 |
| raw_predefined29 | voom | 17268 | 2046 | 3697 | 0 | 10001 |
| raw_predefined29 | DESeq2 | 17268 | 2225 | 4014 | 0 | 10001 |
| adapter_all30 | edgeR | 17032 | 1209 | 2803 | 0 | 10001 |
| adapter_all30 | voom | 17032 | 1310 | 2934 | 0 | 10001 |
| adapter_all30 | DESeq2 | 17032 | 955 | 2400 | 0 | 10001 |
| adapter_predefined29 | edgeR | 17218 | 2418 | 3973 | 0 | 10001 |
| adapter_predefined29 | voom | 17218 | 2144 | 3731 | 0 | 10001 |
| adapter_predefined29 | DESeq2 | 17218 | 2280 | 4043 | 0 | 10001 |


No D42 gene crossed this threshold in these specifications; this does not establish equivalence or absence of biological change. A4 preserves all chronic program methods, including zero-count methods. A5 compares D28 counts under the two reference spaces. Neither aggregate table establishes feature-wise overlaps across processing branches. Original native-model results and their denominators remain separate from these empirical-reference results.

## JGW03: complete separate-Linux reference

The summary below is generated directly from A6; it retains unavailable values in the reported accounting.

| Cell | Reference | Families | Method-feature records | Unavailable point P | Unavailable feature-assignment comparisons |
|---|---|---|---|---|---|
| PC12 | four_group_MC | 31 | 108723 | 138 | 77853 |
| PC12 | high_companion_exact | 12 | 46117 | 3 | 25 |
| C2C12 | four_group_MC | 29 | 113785 | 109 | 104022 |
| C2C12 | high_companion_exact | 11 | 48049 | 9 | 48 |
| HepG2 | four_group_MC | 29 | 118137 | 216 | 62699 |
| HepG2 | high_companion_exact | 11 | 49912 | 5 | 28 |

In the separate Linux four-dose reference, the high-dose/vehicle edgeR families had 2,728/13,325, 6,296/14,387, 2,572/15,004 genes with unadjusted conditional P<0.05 in PC12, C2C12 and HepG2, respectively. These descriptive counts were not multiplicity-adjusted discoveries and cannot rank cellular sensitivity because species, dose grids and tested universes differ. The complete 20-assignment companions had no strict P<0.05 values, as their minimum attainable value was 0.05. Missing required strengths remained explicit rather than reducing denominators.

## Statistical methods


### 1. Assignment populations and analysis units

Each cell line contributed 12 deposited culture libraries, with three libraries in each of four dose roles. The complete balanced four-role assignment population has 12!/(3!)^4=369,600 members. The frozen sampled space contained 10,000 distinct nonidentity assignments and the observed assignment once. The complete high-dose/vehicle companion permuted three high-dose and three vehicle roles among their six libraries, leaving the six low/medium libraries fixed; it contained choose(6,3)=20 assignments, including the observed assignment once. These are conditional references on the same libraries, not independent biological replications or evidence that source treatment labels were randomized.

### 2. Family scope and platform separation

The sampled space comprised 31 PC12, 29 C2C12 and 29 HepG2 method/collection/contrast families; companion counts were 12, 11 and 11. Native feature-key universes and row order were fixed from the observed Linux objects and retained under every assignment; statistic availability could vary across assignments. The exact field map controlled which high-dose, ordinal-dose and log-dose families existed; no extra sensitivity family was invented. PC12 Reactome mroast used 1,999 rotations in an equivalent reconstruction, not authentication of unrecorded historical draws. Linux sources used their own observed results and thresholds and were not numerically pooled with Windows companions.

### 3. Method-specific strengths and availability

Outer strengths were the unchanged native quasi-likelihood F for edgeR; absolute moderated t for limma-voom, mean-z program scores and ULM activities; absolute Wald statistics for DESeq2; minus the natural logarithm of native CAMERA P and of directional PC12 Reactome mroast P; and absolute normalized enrichment score for fgsea. The mroast Mixed result was ancillary rather than a newly added outer family. Strength availability required a finite native P in [0,1] and a strength other than NA/NaN; an unavailable adjusted P alone did not remove an otherwise available strength. Native P=0 yielded positive infinite minus-log-P strength without flooring. Tails compared stored strengths using inclusive exact greater-than-or-equal comparisons, without tolerance or rectifying native F roundoff.

### 4. Inclusive tails and missingness

For an available observed strength, b counted available nonidentity strengths meeting or exceeding it. With no missing required strength, the saved Monte Carlo point label-reference estimate was (1+b)/10,001, whereas the complete companion point P was (1+b)/20. The Monte Carlo estimate was distinct from the full-population tail probability targeted by the finite-population interval. If u required nonidentity strengths were unavailable, point P remained unavailable and the fixed-denominator bounds were (1+b)/N and (1+b+u)/N, with N=10,001 or 20. Unavailable observed strengths retained unavailable points and bounds. Counts of unavailable features and feature-by-assignment comparisons were retained without shrinking denominators.

### 5. Finite-population Monte Carlo uncertainty

For complete sampled feature references, under uniform sampling without replacement and a fixed strength for each assignment, B conditional on M followed a hypergeometric distribution with nonidentity population T=369,599, draws n=10,000 and M nonidentity-population tail hits. Equal-tail inversion retained integer M satisfying both Pr_M(B>=b)>=0.025 and Pr_M(B<=b)>=0.025 and transformed endpoints to (1+M)/369,600. Integer binary searches used log-tail evaluation, with an eight-machine-epsilon scale-dependent guard that conservatively enlarged the numerically accepted set near equality. These pointwise 95% intervals concern outer label-sampling error under those sampling and fixed-statistic conditions. They do not cover biological, inner rotation/multilevel, simultaneous or FDR uncertainty. In complete 20-assignment companions, both finite-population endpoints equal the census point P only for features with an available observed strength and all 19 nonidentity strengths available.

### 6. Multiplicity and interpretation

Native method-specific P/q values, candidate rules and S1 were retained. No new outer empirical BH family was computed in this JGW03 reference; conditional-P threshold counts are unadjusted descriptive quantities, not discoveries or estimates of true FDR. Complete outer companion enumeration does not eliminate inner stochastic uncertainty.




Read the family-level numbers and finite-population fields in A6 with these method definitions. Its pointwise reference intervals are not effect-size confidence intervals or simultaneous coverage. The earlier Windows reference is a separate historical implementation.

## GSE205906: matched alignment comparison

# GSE205906 matched alignment sensitivity

This addition reports the same nine libraries processed through HISAT2 (left) and STAR (right), with the same reference and fixed gene-fragment counting specification. Primary HISAT2 results and conditional-label references remain unchanged. The comparison is processing sensitivity, not independent biological replication.

## Complete method and collection comparison

Findings use strict q<0.05 within each method, collection and contrast. These are saved native/project within-family BH fields, not empirical conditional-label q values. In particular, mroast's project raw-P BH field is distinct from native mid-P FDR. SNI−sham is local disease; PLM−SNI is intervention. Eligible denominators remain branch-specific.

| Method | Collection | Contrast | Eligible HISAT2 / STAR | q<0.05 HISAT2 / STAR | Shared hits | Hit Jaccard | Paired continuous effects | Effect or NES Spearman rho |
|---|---|---|---:|---:|---:|---:|---:|---:|
| gene_edger | genes | SNI−sham | 14217 / 14238 | 0 / 0 | 0 | NA | 14022 | 0.9544 |
| gene_edger | genes | PLM−SNI | 14217 / 14238 | 62 / 79 | 60 | 0.7407 | 14022 | 0.9651 |
| gene_voom | genes | SNI−sham | 14217 / 14238 | 0 / 0 | 0 | NA | 14022 | 0.9527 |
| gene_voom | genes | PLM−SNI | 14217 / 14238 | 38 / 33 | 30 | 0.7317 | 14022 | 0.9639 |
| gene_deseq | genes | SNI−sham | 14217 / 14238 | 0 / 10 | 0 | 0.0000 | 14022 | 0.9544 |
| gene_deseq | genes | PLM−SNI | 14217 / 14238 | 453 / 469 | 411 | 0.8043 | 14022 | 0.9653 |
| unsigned_camera | Reactome | SNI−sham | 1546 / 1544 | 69 / 65 | 63 | 0.8873 | 0 | NA |
| unsigned_camera | Reactome | PLM−SNI | 1546 / 1544 | 267 / 282 | 257 | 0.8801 | 0 | NA |
| unsigned_mroast | Reactome | SNI−sham | 1546 / 1544 | 0 / 0 | 0 | NA | 0 | NA |
| unsigned_mroast | Reactome | PLM−SNI | 1546 / 1544 | 0 / 0 | 0 | NA | 0 | NA |
| unsigned_fgsea | Reactome | SNI−sham | 1546 / 1544 | 109 / 120 | 104 | 0.8320 | 1542 | 0.9634 |
| unsigned_fgsea | Reactome | PLM−SNI | 1546 / 1544 | 334 / 362 | 317 | 0.8364 | 1542 | 0.9851 |
| unsigned_camera | SynGO_sensitivity | SNI−sham | 78 / 78 | 31 / 31 | 31 | 1.0000 | 0 | NA |
| unsigned_camera | SynGO_sensitivity | PLM−SNI | 78 / 78 | 15 / 14 | 14 | 0.9333 | 0 | NA |
| unsigned_mroast | SynGO_sensitivity | SNI−sham | 78 / 78 | 0 / 0 | 0 | NA | 0 | NA |
| unsigned_mroast | SynGO_sensitivity | PLM−SNI | 78 / 78 | 0 / 0 | 0 | NA | 0 | NA |
| unsigned_fgsea | SynGO_sensitivity | SNI−sham | 78 / 78 | 39 / 34 | 33 | 0.8250 | 78 | 0.9881 |
| unsigned_fgsea | SynGO_sensitivity | PLM−SNI | 78 / 78 | 20 / 21 | 20 | 0.9524 | 78 | 0.9866 |
| activity_results | PROGENy | SNI−sham | 14 / 14 | 1 / 0 | 0 | 0.0000 | 14 | 0.9736 |
| activity_results | PROGENy | PLM−SNI | 14 / 14 | 0 / 3 | 0 | 0.0000 | 14 | 1.0000 |
| activity_results | DoRothEA_AB | SNI−sham | 81 / 80 | 0 / 0 | 0 | NA | 80 | 0.9659 |
| activity_results | DoRothEA_AB | PLM−SNI | 81 / 80 | 0 / 0 | 0 | NA | 80 | 0.9733 |

NA Jaccard denotes an empty union, not complete agreement. CAMERA and mroast expose categorical direction rather than a validated continuous effect/NES in this comparison; their continuous effect correlation is therefore NA. A high correlation does not establish equivalence, absence of confounding or invariant significance. Features are not biological replicates.

## Interpretation

Intervention gene findings overlapped substantially but did not coincide. DESeq2 local-disease findings changed from zero to 10; PROGENy q<0.05 counts changed from one to zero for disease and zero to three for intervention despite all 14 paired effect signs agreeing in both contrasts. These differences remain visible rather than being reduced to a claim of complete robustness. The constituent-combination context does not establish whole-formula efficacy, product/exposure equivalence or a PSS mechanism. Animal linkage, pooling, batch and dose histories remain unresolved.



In A7, `left` means HISAT2 and `right` means STAR. `local_disease` means SNI−sham and `intervention` means PLM−SNI. `NA` correlations/Jaccard values remain unavailable rather than becoming zero. Full saved-export aggregation was previously reconciled across 216,070 family-feature records with 528 assertions; the present package preserves the accepted 22-row summary and does not repeat alignment, counting or fitting. Missing historical stochastic internals are not retrospectively reconstructed.

## Provenance and scope of verification

[PROVENANCE.json](PROVENANCE.json) binds the seven source-preserved tables to their public paths, row counts and SHA-256 identities. A1 and A6 reuse the previously published data/derived files; the other five tables are supplied here. The repository-level [manifest](../../MANIFEST.json) covers this addition.

Before public packaging, the local reading assembly passed 653 scoped checks: 440 GSE128623 field reconciliations, seven source-file identities and 206 displayed-cell checks (24 JGW overview, six JGW high-dose and 176 matched-quantification numeric cells). Those checks concern the stated aggregates and displayed summaries; they do not independently reproduce the upstream experiments or every historical result. The original local receipt is identified by hash in the provenance file, rather than presented as a public download.

The accepted tables and existing reports, rather than newly calculated P values, are the basis of this package. No model, scoring, label allocation, alignment or counting was repeated for publication. This addition does not establish whole-study raw-to-result restoration or new biological validation. Original explanatory text and project aggregate tables follow the repository CC BY 4.0 notice; source-derived inputs and upstream resources retain their own terms. No raw reads, article PDFs, source gene-set memberships or raw failure messages are included.
