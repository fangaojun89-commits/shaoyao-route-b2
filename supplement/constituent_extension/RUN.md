# Mixture and single-herb comparisons in GSE298414

Version 0.26.0 adds all 60 new mixture-minus-single-herb contrasts from 126 PC12 culture libraries (42 conditions, three libraries each). The earlier 72-library/36-contrast block is preserved. The two blocks share 72 libraries and are not independent replications. This finite preparation extension accompanies article v41.

## Contents and interpretation

The complete gene tables contain 13,898 genes in every contrast for both robust edgeR and limma-voom, for the primary 126-library model and the planned 125-library omission of GSM9013891. Gene filter, input identities, sample-level QC, PCA, all summaries and complete effect-influence comparisons are supplied. Local `rsem_path` columns are removed from three public metadata files; all other fields are unchanged and both original and published hashes are recorded.

The four program collections comprise 240 collection–contrast families: 92,640 CAMERA membership tests and 7,260 signed primary tests. Three fixed-score omissions add 21,780 signed tests. Public library-score matrices allow all 29,040 signed effect estimates to be reconstructed as group-mean differences. Full program coverage, representative-symbol choices, including nonsignificant results, and class/within-contrast multiplicity are retained.

At high ordinal level in ethanol, the 1:1 separately extracted mixture differs from both single-herb profiles. The 2:1 co-extract is substantially influenced by GSM9013891; both primary and omission outcomes are reported. Equal total-extract labels do not equalize component exposures, and these contrasts have no additive reference response. They do not estimate synergy, antagonism or clinical antispastic efficacy. Absolute-dose conflicts and the missing 1:1 co-extract remain unresolved. The complete [results report](Route_B2_GSE298414_single_herb_results_20261025.md) describes all comparisons, including the absence of exhaustive within-group influence testing for the 1:1 example.

## Recheck retained results

With Python, numpy and pandas installed, run from the repository root:

```sh
python pipeline/constituent_readback.py --output /path/to/new/constituent-readback.json
```

The output must not already exist. The command reads only published payloads and writes the specified new receipt. It checks every payload identity; 60 contrasts and both comparator classes; full and omission gene domains and BH values; both methods' complete summaries; all 60 gene-influence comparisons; program BH, coverage and all 240 family summaries; and all 29,040 signed effects from retained library scores. No credentials, cloud account or external networks are required. The command does not fit scientific models or establish clean-machine raw-to-result restoration. Confidence intervals are checked for ordered bounds containing their estimates, not independently refitted.

## Scientific sources

Unchanged historical scientific sources are in `analysis/constituent_extension/original/`; they contain local analysis paths and are not portable entry points. They document input selection, explicit gene-level lengthScaledTPM, TMM, the 42-condition model, robust edgeR/voom, program mapping and inference, and the independent readback. The finite [analysis specification](../../analysis/constituent_extension/GSE298414_SINGLE_HERB_COMPARISON_20261025.md) records decisions before fitting. The printed-matrix display calls in the gene source were wrapped in `invisible()` after execution; this only suppresses console output and no model was rerun.

`data/constituent_extension/GSE298414/PROVENANCE.json` records the official source archive, frozen resource identities, retained primary model identity and publication transformations. `input_receipt.json` records each selected archive member's identity. `gene_receipt.json` retains its historical pre-review status; the current review status is in PROVENANCE and the public verifier produces an independent current receipt. Session details are supplied for genes and programs. Readback does not regenerate score matrices from external resources or model binaries.

The fixed-score omissions hold primary mapping and scores fixed; they are not CAMERA or complete-pipeline sensitivities. Multiplicity is separate for each method/resource and the 30 contrasts in each comparator class. Shared comparators and libraries must not be counted as independent replications. Earlier canonical results and S1 candidate admission are unchanged. Independent paired PSS validation, exposure bridging and whole-study cold restoration remain open.

Original code is MIT; original statistical outputs, reports and curation are CC BY 4.0. Public GEO annotations keep their upstream provenance. Article PDFs, source supplemental workbooks, gene-set memberships, regulatory network edges and R model binaries are not redistributed. External resources retain their own terms; see THIRD_PARTY_NOTICES.md.
