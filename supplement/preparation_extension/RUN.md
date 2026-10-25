# GSE298414 preparation comparisons

This extension supplies all 36 retained gene contrasts for the 72-library preparation analysis, all four program collections, full fixed-score omission tables and the six-library gene-influence summaries. It accompanies article v39. The 36 contrasts reuse samples and do not constitute 36 independent replications.

The high-level ethanol 2:1 extraction comparison yielded 589 edgeR versus 17 voom findings across the extraction family. Omitting GSM9013891 reduced the edgeR count to 27 and retained nine original discoveries. All six affected libraries were examined; the primary analysis retains all 72. Foxo1, Hif1a and Smad3 score directions persisted after that omission but lost collection-by-extraction-class adjusted significance. Complete results, including other solvent/ratio comparisons and nonsignificant rows, are retained.

## Reuse the retained results

From the repository root, with Python, numpy and pandas installed:

```sh
python pipeline/preparation_readback.py --output /path/to/new/preparation-readback.json
```

The command verifies published payload hashes, recomputes BH arithmetic for the original 12 and additional 24 gene contrasts and all program tests, and checks 34,848 signed full/omission effects from library score means. It compares all 144 program-family summary rows. It does not refit gene models, infer scores from networks or claim raw-to-result restoration. The output file must not already exist. This command requires no credentials, cloud service or upstream network files.

## Original scientific analysis

Unchanged original sources and their hashes are in `analysis/preparation_extension/original/`. They retain historical absolute paths and must not be run as portable entry points. In dependency order: prepare GEO counts; reconcile article supplements; import 72 selected RSEM gene tables and fit the extraction model; run systematic high-level influence; extend the saved models to solvent/ratio contrasts and program analyses. The two dated specifications describe details. Original source acquisition includes all 141 deposited libraries; the accepted model uses the 72 listed in `full_72/sample_qc.tsv`.

The original archive identity and the four frozen-resource hashes are in `data/preparation_extension/GSE298414/PROVENANCE.json`. Source article DOI: 10.1038/s41597-026-06759-6. Retrieve the original source supplements and resource versions separately under their applicable terms. This release does not redistribute article files, gene-set memberships, networks or R model binaries. Consequently it is a complete retained-result/readback package for the stated contrasts, not a demonstrated clean-machine reconstruction of the scientific fits. Environment details for the new program extension are in its `sessionInfo.txt`.

## Interpretation and scope

GEO/supplement and paper absolute-dose values disagree, so comparisons retain ordinal levels. The source lacks 1:1 CEM. Two extraction-label conflicts are preserved in the 141-library source ledger; one lies in the 72-library block and has an omission sensitivity. The comparison uses positive effective length in every selected library and explicit gene-level lengthScaledTPM before TMM/edgeR and voom. It does not reconstruct transcript quantification. CAMERA is a membership test; ULM score differences are inferred activity, not direct TF measurements.

Class-wide BH families each contain 12 contrasts and remain separate for extraction, solvent and ratio, and for each program collection. Fixed-score omissions hold mapping and scores fixed; they are not full pipeline or CAMERA sensitivity. The six gene-influence affected-contrast tables include global q values from the original full 12-contrast fits; the rest of those omission fit matrices are not supplied here. Their global q values cannot be independently reconstructed from the affected-contrast subset alone.

Canonical JGW03/GSE128623/GSE205906/GSE277064 outputs and the zero-admission S1 result are unchanged. Single-herb comparisons, exposure bridging, independent paired PSS validation and whole-study cold restoration remain open. No model is run by the public readback command.

Original code is MIT; original reports, curation and statistical outputs are CC BY 4.0. Public source annotations retain their upstream provenance. No licence is granted for external resources by this repository; see THIRD_PARTY_NOTICES.md.
