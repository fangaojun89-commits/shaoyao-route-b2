# Completed conditional references and matched quantification: supplementary results

25 October 2026. Accepted late-stage reporting results, accompanying the current English manuscript.

This public reading package brings together final late-stage results that previously required several dated addenda. It complements Supplementary Tables S1–S18 and Figures S1–S7; A1–A7 are package-local identifiers, not renumbered journal tables. It does not replace the original files or add a new analysis. All seven TSV files preserve their source bytes, precision and missing-value labels.

The study's fixed spasticity anchor admitted no programs. The analyses below describe source-wide conditional references and same-library technical sensitivity; downstream candidate-restricted links remained unestimated. GSE128623 empirical BH thresholds, JGW03 unadjusted conditional P thresholds, and GSE205906 native within-family BH thresholds have different meanings and must not be pooled.

## Tables and reading order

@@INDEX@@

A1 is the complete GSE128623 family-level reference table. A2–A5 are compact source-preserved views of its main comparisons, checked against A1. A6 contains the complete JGW03 separate-Linux family summary. A7 contains the complete matched HISAT2/STAR comparison summary. Family-feature events can count the same feature repeatedly across methods and contrasts; they are not sample sizes or unique molecular entities.

## GSE128623: completed reference spaces

Retrospective GSE128623 day/sham conditional references enumerated all 252 assignments within each specification and day, except predefined29 D4/sham, which had 126; other groups remained fixed. Inclusive tails counted the observed assignment once, with separate empirical BH families for each method, collection and contrast. Activity strengths used absolute moderated t from limma score models. A subsequent descriptive comparison evaluated signed effects among common programs with finite, nonzero effects in both specifications. Separately standardized effect magnitudes were not pooled. Separately completed six-group references comprised the observed assignment and 10,000 prespecified sampled assignments per processing/cohort specification, with 67 method/collection/contrast families per specification (268 in total). Saved-tail arithmetic, diagnostics and metadata checks authenticated their reporting outputs without rerunning completed models. The six-group and day/sham references condition on different label spaces; neither authenticates biological exchangeability or true FDR control.

The four specifications retain their original processing and cohort definitions. In A1–A5, the threshold is empirical BH q≤0.05 within each stored method/collection/contrast family. Counts are point-estimate threshold crossings, not interval-robust findings or effect magnitudes. `native_fraction_allocations_any_q_le_05` describes the same-library allocation space and is not an estimate of true FDR. The day-specific and six-group reference spaces answer different conditional questions.

### Family overview (A2)

@@A2@@


### Gene counts (A3)

@@A3@@


No D42 gene crossed this threshold in these specifications; this does not establish equivalence or absence of biological change. A4 preserves all chronic program methods, including zero-count methods. A5 compares D28 counts under the two reference spaces. Neither aggregate table establishes feature-wise overlaps across processing branches. Original native-model results and their denominators remain separate from these empirical-reference results.

## JGW03: complete separate-Linux reference

The summary below is generated directly from A6; it retains unavailable values in the reported accounting.

@@JGW@@

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

@@A7@@

NA Jaccard denotes an empty union, not complete agreement. CAMERA and mroast expose categorical direction rather than a validated continuous effect/NES in this comparison; their continuous effect correlation is therefore NA. A high correlation does not establish equivalence, absence of confounding or invariant significance. Features are not biological replicates.

## Interpretation

Intervention gene findings overlapped substantially but did not coincide. DESeq2 local-disease findings changed from zero to 10; PROGENy q<0.05 counts changed from one to zero for disease and zero to three for intervention despite all 14 paired effect signs agreeing in both contrasts. These differences remain visible rather than being reduced to a claim of complete robustness. The constituent-combination context does not establish whole-formula efficacy, product/exposure equivalence or a PSS mechanism. Animal linkage, pooling, batch and dose histories remain unresolved.



In A7, `left` means HISAT2 and `right` means STAR. `local_disease` means SNI−sham and `intervention` means PLM−SNI. `NA` correlations/Jaccard values remain unavailable rather than becoming zero. Full saved-export aggregation was previously reconciled across 216,070 family-feature records with 528 assertions; the present package preserves the accepted 22-row summary and does not repeat alignment, counting or fitting. Missing historical stochastic internals are not retrospectively reconstructed.

## Provenance and reporting scope

This report is generated solely from the seven public retained tables and a versioned interpretation template. [SOURCE_TABLES.json](SOURCE_TABLES.json) records the public input paths and exact copied identities; [SUMMARY_RECONCILIATION.tsv](SUMMARY_RECONCILIATION.tsv) records 440 field comparisons of A2–A5 against A1. [RECEIPT.json](RECEIPT.json) binds the generated files, code/template identities and actual Python version.

The A2/A3 tables, six-row JGW overview and complete 22-row matched-alignment display are generated from these inputs; A4/A5 remain available in full as linked TSV files. All seven TSVs retain their source bytes, precision, missing values and row order. The reading definitions are the retained interpretation of these fixed accepted results, not an automatic interpretation of arbitrary replacement data.

This command does not fit models, recalculate P/q values, draw labels, reconstruct upstream archive arithmetic, install software, contact external services, or establish whole-study cold restoration. Original text and project aggregate tables follow CC BY 4.0; original code follows MIT. Upstream source-specific terms remain applicable.
