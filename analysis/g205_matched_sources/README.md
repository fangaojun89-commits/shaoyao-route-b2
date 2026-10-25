# GSE205906 matched STAR/HISAT2 analysis: source guide

This guide connects the retained code for the completed nine-library processing sensitivity to the accepted 22-family comparison. Nineteen of the 20 selected sources were already distributed elsewhere in this repository; their existing copies are reused. The metadata/reservoir preparation source is newly supplied. All selected source bytes are unchanged. No alignment, counting, model, random-label reference or accepted-result check was replayed to create this guide.

## Scientific interpretation

The same nine libraries were processed with STAR and HISAT2 against Rnor6/Ensembl101 and compared under the fixed gene-fragment specification. This is a within-library processing sensitivity, not independent biological replication. The original HISAT2 analysis and candidate decisions remain unchanged. The complete 22-family summary is available in [the late-results supplement](../../supplement/late_results/README.md), including changes in thresholded findings. Native within-family BH values must not be substituted for conditional-label reference P or q values.

## Source order

1. `prepare_gse205906_matched_quantifier_inputs_20261015.py` binds retained sample metadata, prior QC and the existing read reservoir. Its CLI performs file operations; it is not a general-purpose dry-run entry.
2. `run_gse205906_matched_STAR_20261015.py` contains the fixed reference admission, index, library and matrix stages. It retains STAR, samtools and featureCounts identities and counting options. The one-library frontend and serial continuation preserve the original library function while changing scheduling boundaries. Their historical host/CPU/window checks describe the recorded run, not a portable configuration for a new machine.
3. `run_gse205906_matched_downstream_custody_v1_20261022.R` reads the accepted STAR matrix and calls `gse205906_matched_downstream_core_custody_v1_20261022.R`. The core imports the pinned original gene and program methods. Eleven caller/core/configuration/dependency identities were matched directly to the retained scientific completion record when this source selection was made.
4. `run_gse205906_matched_saved_comparison_custody_v2_20261022.R` compares the saved STAR and HISAT2 results through the original comparison module. The original saved-comparison entry is included as a pinned dependency. It is not an instruction to repeat the earlier HISAT2 preparation. The v2 entry retains the completed preparation and corrects a byte-size representation check without changing the comparison mathematics.
5. `verify_gse205906_saved_comparison_counts_20261022.py` independently reconciles exported eligibility, returned-result, q-threshold, overlap and sign counts. In the original completed check, it covered 216,070 unique family-feature records and 528 assertions across 22 families. These historical check totals are not a claim that this source-release step reran that check.

## Locate and verify the sources

[`PROVENANCE.json`](PROVENANCE.json) maps each original project-relative source role to its actual public location, with its size, SHA-256 and identity-evidence category. Read the `public_path` relative to the repository root. Existing shared source copies are deliberately not duplicated. The provenance record distinguishes identities matched to the accepted completion record from additional retained source bytes pinned by this selection; it does not treat every local file as independently authenticated historical execution.

From the repository root, the existing standard-library verifier checks the release payload without executing scientific sources:

```sh
python -I -B pipeline/sources.py verify
```

For an investigator who already holds the two original complete exports, the original saved-export checker accepts a directory containing `full_source_union.tsv` and `native_table_summary.tsv`:

```sh
python -I -B PATH_TO_THE_MAPPED_EXPORT_CHECKER.py /path/to/saved_exports
```

Replace the source path using its `public_path` in the provenance map. The seven-table public reading supplement supplies a reporting projection, not those two complete native exports; do not pass that projection as the checker input.

## Inputs and reconstruction scope

The complete nine-library HISAT2 strand-0 count matrix and all three complete STAR strand matrices are already available in the previously published [count-input package](../../data/count_inputs/README.md). Each has 32,883 gene rows and the same nine-library order; STAR strand 0 remains the primary counting specification. No new alignment or counting is needed to recover these accepted inputs. The three strand matrices are technical alternatives from the same libraries, not three independent datasets.

From the repository root, restore the exact original bytes into a new directory:

```sh
python -I -B pipeline/count_inputs.py materialize --output recovered-inputs
```

The command restores all eleven published count-package files, including both datasets; it neither fits models nor changes accepted results. `data/count_inputs/PROVENANCE.json` records each original size and SHA-256, compressed identity and relative destination. The four GSE205906 matrix paths beneath `recovered-inputs` are:

- `GSE205906_full_read_counts_fasp_v18e_20261005/formal/output/GSE205906_fixed9_gene_fragments_s0.tsv`
- `GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s0.tsv`
- `GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s1.tsv`
- `GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s2.tsv`

Raw FASTQs, the genome index, native fitted objects and all original process receipts are not bundled by this source guide. [Acquisition records](../../data/acquisition/README.md), [resource build](../resource_build/README.md) and [runtime bindings](../runtime_bindings/README.md) provide further reconstruction components. Count restoration alone does not create a portable execution adapter or restore the complete native scientific environment.

Historical absolute paths and receipt checks remain intact. Copying these sources into an arbitrary directory is insufficient to run the historical workflow. A future portable execution adapter would need explicit input, output and environment mappings and separate validation; this release supplies no such adapter and does not establish whole-study cold restoration. It makes the scientific implementation and its dependency relationships easier to inspect while preserving the completed calculations.

The retained checks support source and process identity, saved numerical components and descriptive export reconciliation. They do not reconstruct missing historical rotation draws, multilevel trajectories or every stochastic internal. The original S1 result and claim boundaries are unchanged.

Original code is covered by the repository MIT license. Source-derived data and gene-set resources retain their own terms; this guide does not extend their redistribution permissions.

## Source location table

| Original source role | Public file | Reused copy |
|---|---|---|
| `code/R/run_gse205906_matched_downstream_custody_v1_20261022.R` | [source](../../analysis/reference_sources/code/R/run_gse205906_matched_downstream_custody_v1_20261022.R) | yes |
| `code/R/gse205906_matched_downstream_core_custody_v1_20261022.R` | [source](../../analysis/reference_sources/code/R/gse205906_matched_downstream_core_custody_v1_20261022.R) | yes |
| `code/R/analyze_gse205906_bulk_model_v2_r4.R` | [source](../../analysis/study_sources/code/R/analyze_gse205906_bulk_model_v2_r4.R) | yes |
| `code/R/gse205906_bulk_model_v2_r4_lib.R` | [source](../../analysis/study_sources/code/R/gse205906_bulk_model_v2_r4_lib.R) | yes |
| `config/GSE205906_Ensembl101_full_model_v1_20261006.json` | [source](../../analysis/study_sources/config/GSE205906_Ensembl101_full_model_v1_20261006.json) | yes |
| `code/R/gse205906_ensembl101_input_loader_v2_20261006.R` | [source](../../analysis/study_sources/code/R/gse205906_ensembl101_input_loader_v2_20261006.R) | yes |
| `code/R/gse205906_ensembl101_effect_independent_prefit_v1_20261006.R` | [source](../../analysis/study_sources/code/R/gse205906_ensembl101_effect_independent_prefit_v1_20261006.R) | yes |
| `code/R/gse205906_ensembl101_annotation_adapter_v1_20261006.R` | [source](../../analysis/study_sources/code/R/gse205906_ensembl101_annotation_adapter_v1_20261006.R) | yes |
| `code/R/gse205906_ensembl101_nofit_RAM_worker_v2_20261006.R` | [source](../../analysis/study_sources/code/R/gse205906_ensembl101_nofit_RAM_worker_v2_20261006.R) | yes |
| `code/R/gse205906_resource_object_validation_v1.R` | [source](../../analysis/study_sources/code/R/gse205906_resource_object_validation_v1.R) | yes |
| `code/R/gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R` | [source](../../analysis/study_sources/code/R/gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R) | yes |
| `code/python/run_gse205906_matched_STAR_20261015.py` | [source](../../analysis/reference_sources/code/python/run_gse205906_matched_STAR_20261015.py) | yes |
| `code/python/run_gse205906_STAR_one_remaining_library_20261022.py` | [source](../../analysis/reference_sources/code/python/run_gse205906_STAR_one_remaining_library_20261022.py) | yes |
| `code/python/verify_gse205906_saved_comparison_counts_20261022.py` | [source](../../analysis/reference_sources/code/python/verify_gse205906_saved_comparison_counts_20261022.py) | yes |
| `code/R/gse205906_matched_saved_comparison_20261015.R` | [source](../../analysis/reference_sources/code/R/gse205906_matched_saved_comparison_20261015.R) | yes |
| `code/R/run_gse205906_matched_saved_comparison_custody_v2_20261022.R` | [source](../../analysis/reference_sources/code/R/run_gse205906_matched_saved_comparison_custody_v2_20261022.R) | yes |
| `code/R/run_gse205906_matched_saved_comparison_20261015.R` | [source](../../analysis/reference_sources/code/R/run_gse205906_matched_saved_comparison_20261015.R) | yes |
| `code/python/prepare_gse205906_matched_quantifier_inputs_20261015.py` | [source](../../analysis/g205_matched_sources/code/python/prepare_gse205906_matched_quantifier_inputs_20261015.py) | new |
| `code/python/continue_GSE205906_STAR_remaining_serial_20261022.py` | [source](../../analysis/reference_sources/code/python/continue_GSE205906_STAR_remaining_serial_20261022.py) | yes |
| `code/python/start_GSE205906_STAR_post_G128_prefix_v2_20261022.py` | [source](../../analysis/reference_sources/code/python/start_GSE205906_STAR_post_G128_prefix_v2_20261022.py) | yes |
