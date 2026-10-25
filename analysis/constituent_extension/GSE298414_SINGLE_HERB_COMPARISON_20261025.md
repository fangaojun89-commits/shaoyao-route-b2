# Mixture versus single herb preparation comparisons

25 October 2026. This finite extension addresses the remaining single-herb preparation question using existing GSE298414 gene estimates. Earlier 72-library results and all canonical analyses remain unchanged. The new analysis estimates only contrasts not previously tested.

## Question and observations

Compare each available two-herb preparation with peony alone and licorice alone at the same source-reported ordinal level and solvent. The 126 extract libraries comprise 42 conditions with three source-reported culture replicates each: water/70% ethanol, three levels, and seven codes (01 peony, 05 licorice, 02/04 co-extracted 2:1/1:2, and 06/07/08 separately extracted then mixed 2:1/1:1/1:2). All have recorded batch A. Do not pool other accessions or use the untreated/drug libraries as matched vehicle or negative controls.

Five mixture codes × two solvents × three levels × two single-herb comparators gives 60 new contrasts, all positive mixture minus single herb. None is an earlier extraction, solvent or ratio contrast. The two classes, mixture-minus-peony and mixture-minus-licorice, each contain 30 contrasts. Source dose conflicts persist; the analysis uses ordinal levels without claiming resolved absolute concentration. Equal total extract labels do not match individual ingredient exposures. These contrasts do not identify synergy, antagonism or clinical efficacy. The absence of within-source 1:1 CEM remains unresolved.

## Analysis

Reuse 72 already extracted RSEM files and unpack only the additional 54 from the retained official archive. Use all 126 for a new 42-condition model required by these additional observations and estimands. Explicitly convert gene TPM and across-library mean effective lengths with tximport lengthScaledTPM, rescaled to original expected-count totals. Require positive effective length in all selected libraries and design-aware filterByExpr; retain fractional values. Fit robust edgeR QL with TMM and limma-voom with robust empirical Bayes, nominal voom 95% intervals, within-contrast BH and separate 30-contrast class-wide BH by method. Do not replace or recompute the existing 72-library analysis. Design rank is 42, with nominal residual df 84.

Preserve all 126 in the primary analysis. Because GSM9013891 was already influential in the previous preparation block, perform a planned omission for the new gene contrasts, retaining the primary gene filter and original length-scaled input values but recalculating TMM and model estimates for the 125 remaining libraries. Compare complete gene domains, effects, discovery overlap and directions; do not select the preferred result. This is a new-comparison influence analysis, not a rerun of the accepted 36 contrasts.

Use the same frozen resources and representative-symbol mapping rule as the preceding preparation extension. CAMERA on saved voom expression/weights uses inter-gene correlation 0.01 and 10–500 expressed members. ULM scores use logCPM with center=FALSE, minimum 50 PROGENy targets and five rat DoRothEA A/B targets. Fit library score contrasts under the 42-condition design and retain nominal intervals. Collection-by-comparator-class BH spans 30 contrasts separately. Retain full mappings and coverage. Fixed-score omissions cover GSM9013891 and both source-label discrepancies GSM9013914/GSM9013917; they do not imply complete raw-to-result or CAMERA robustness.

## Interpretation and verification

All 60 contrasts and all eligible programs are reported, including null results. Shared comparators and reused libraries are not independent replications. Differences between significance calls do not constitute interaction tests. Results remain preparation context; no candidate gate is changed. Check contrast direction from generated group means, complete input identities, finite estimates and multiplicity arithmetic, and independently reconstruct library-score effects and table summaries from retained outputs. No completed scientific model is rerun for reporting checks. Store bulk outputs on R; no cloud operation is required for this finite local comparison.
