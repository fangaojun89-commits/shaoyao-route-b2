# GSE298414 within-source preparation comparisons

25 October 2026, specified before differential-expression fitting. This implements the preparation-sensitivity component of SAP section 6.5 using newly acquired processed data. Existing accepted analyses are not refitted.

## Question and fixed scope

Compare co-extraction (CEM) with post-extraction mixing (IEM) in PC12 cells within GSE298414, at matching solvent, PR:GR ratio and reported ordinal level. Include the 72 libraries for codes JGW/JGE02, 04, 06 and 08. There are 24 conditions with three independent biological replicates per condition, as described in the source paper. Replicate numbers alone do not establish pairing across conditions. All selected libraries have recorded batch A. No untreated or cross-accession reference is borrowed.

Twelve contrasts are fixed in `R:/Shaoyao_storage/data/GSE298414_preparation/inputs_20261025/cem_iem_contrasts_supplement_verified.tsv`: CEM minus IEM for two ratios, two solvents and three levels. Model: independent libraries, `~0 + condition` (24 coefficients; 48 nominal residual degrees of freedom). No continuous-dose trend or absolute-dose equivalence is assumed. Source labels 20/100/500 and paper methods 4/20/100 remain unresolved.

## Import, filtering and inference

Import the original deposited RSEM gene tables with tximport (`type="rsem", txIn=FALSE, txOut=FALSE`). The installed gene-level RSEM importer ignores a requested abundance conversion and returns raw counts with a warning. Therefore explicitly use the exported `makeCountsFromAbundance` function: multiply each deposited gene TPM by its mean effective gene length across the 72 samples, then rescale each column to its original expected-count total. This is a length-scaled, count-scale estimate derived from the deposited rounded gene TPM and effective lengths, not a fresh transcript-level reconstruction. Verify the formula independently and verify that conversion actually changes the raw counts. Do not apply a second length offset or round estimates. Check sample ordering, gene identity and equality of original/adjusted column totals. Preserve the earlier raw expected-count matrix.

Primary model: edgeR TMM normalization, `filterByExpr` with the categorical design, robust dispersion estimation and robust quasi-likelihood fitting. Secondary method: limma-voom on the same filtered input and design, with robust empirical Bayes and 95% confidence intervals. Report all tested genes and all 12 contrasts; report BH within each contrast and BH across the complete gene-by-12-contrast family separately. Neither significance threshold nor discovery count constitutes evidence of equivalence when effects are not detected.

Retain sample-level library sizes, normalization factors, log-expression correlations and PCA coordinates. No outcome-driven exclusions. A fixed additional influence analysis removes only GSM9013914, the library with the conflicting published extraction label, and repeats the new-source comparison with the same retained gene universe. It has 71 libraries, 24 coefficients and 47 nominal residual degrees of freedom. This is not a correction of an existing accepted fit. Summarize effect correlation, sign agreement, shared discoveries and changes for the affected low-level ethanol 2:1 contrast; retain all sensitivity contrasts.

The matching annotation from the existing PC12 resource is reused only to label genes. Enrichment/activity analyses and solvent, ratio and single-herb extensions remain subsequent explicit work. S1, all candidate gates and canonical JGW03 results remain unchanged. These data are preparation context from the same research program, not independent disease-phenotype validation or therapeutic efficacy evidence.

## Sources

- [Source study](https://doi.org/10.1038/s41597-026-06759-6), Tables 2/3 and recovered Table S3.
- [Input/source reconciliation](../../supplement/preparation_extension/Route_B2_GSE298414_preparation_extension_20261025.md).
- [Official tximport vignette](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html), RSEM gene import and length-scaled abundance strategy.

Pre-fit input clarification: the deposited RSEM tables include 735 positive-count entries (178 genes; maximum expected count 5.16) with zero effective length and zero TPM. Record these entries, require zero TPM wherever effective length is zero, and exclude genes with a nonpositive effective length in any selected library before inference. This structural availability filter is combined with expression filtering and frozen for the omission analysis. The raw assignments remain preserved.
