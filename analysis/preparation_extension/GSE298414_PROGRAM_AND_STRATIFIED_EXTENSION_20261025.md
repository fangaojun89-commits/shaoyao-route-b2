# Reusing GSE298414 fits for stratified preparation and program comparisons

25 October 2026. Specified before the following new comparisons. This extension reuses the accepted new-source 72-library import, gene filter, TMM normalization, robust edgeR dispersion/QL fit and saved voom expression/precision weights. It does not reopen S1 or refit the canonical JGW03 source.

## Additional contrasts

Retain the original 12 CEM-minus-IEM contrasts. Add 12 ethanol-minus-water contrasts within extraction method, PR:GR ratio and ordinal level, and 12 PR:GR 2:1-minus-1:2 contrasts within extraction method, solvent and ordinal level. All coefficients are contrasts of the existing saturated 24-condition design. These 36 contrasts share libraries and are not independent replications. No interaction is claimed from a difference in significance between strata. Absolute-dose claims remain excluded.

For the 24 new gene contrasts, apply the saved edgeR QL fit. Reconstruct the base weighted linear coefficients from the saved voom expression/weights once for the new limma contrasts; do not recalculate the expression matrix, normalization, dispersion or voom weights. Preserve all existing tests. Report within-contrast BH and class-wide BH (each of the three 12-contrast classes separately), never reassign the existing extraction family by pooling all 36 tests.

## Programs and frozen resources

Use the existing frozen Reactome97, QuickGO SynGO sensitivity, PROGENy human top500 and DoRothEA rat A/B resources. Reuse the canonical gene annotation and its strict human ortholog mapping. As in the existing canonical source, where multiple retained gene IDs share a symbol, select the one with the highest mean TMM logCPM across all 72 libraries, resolving ties by raw gene ID. Record the selected mapping and freeze it for influence checks.

CAMERA tests use the saved human-mapped voom expression and weights, the same condition design, and fixed inter-gene correlation 0.01. Unsigned collections include programs with 10–500 expressed mapped members; retain coverage records for every reference program. Report membership direction and descriptive member-gene effects, without interpreting them as signed pathway activation. This saved-voom implementation is an explicitly documented new-source method; it is not an exact reproduction of the canonical DGEList CAMERA object.

Signed footprints use decoupleR ULM with center=FALSE, the full mapped expression background, original signed weights, and the existing canonical eligibility settings: at least 50 expressed targets for PROGENy and five for DoRothEA. Fit library-level score differences with the 24-condition design and robust limma empirical Bayes; report effects, 95% intervals and collection-specific BH. The ULM gene-regression P value is not the treatment-effect P value. Raw score magnitudes across sources are not treated as a shared measurement scale.

For each program collection, retain both within-contrast BH and collection-by-contrast-class BH across all 12 contrasts. No outcome-selected pathway subset is tested. Report all findings, including null and inconsistent results.

## Targeted influence

Retain the full analysis. For signed activities, use the same saved full-sample scores and remove each of the six high-level ethanol 2:1 comparison libraries once; separately remove the disputed JGE06 low-level library. Refit only the new library-level score contrasts, with identical program membership and background. These are fixed-score influence checks, not complete reruns of gene normalization, CAMERA or the whole pipeline. Their scope must be explicit in reporting. The unusual library is not deleted and no preferred omission replaces the primary result.

Single-herb comparisons, missing within-accession 1:1 CEM conditions, quantitative exposure bridging, independent PSS validation and other original gaps remain separate work. The current finite comparisons cannot resolve those requirements by themselves.
