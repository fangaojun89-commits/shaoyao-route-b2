# Complete saved GSE205906 matched-comparison exports

These two unchanged tables make the previously accepted HISAT2-versus-STAR comparison directly inspectable. The archive contains `full_source_union.tsv` (216,070 unique family-feature rows across 22 families) and `native_table_summary.tsv`. Both SHA-256 identities match the accepted 22 October 2026 export-count receipt. They were retrieved and checked again on 26 October without fitting models, changing rows or recomputing scientific estimates.

From the repository root, with Python 3.12 or later:

```sh
python -I -B pipeline/sources.py verify
python -I -B -m tarfile -e data/g205_matched_exports/GSE205906_saved_exports_20261026.tar.gz recovered-g205-exports
python -I -B analysis/reference_sources/code/python/verify_gse205906_saved_comparison_counts_20261022.py recovered-g205-exports
```

Use a new, empty output directory. `PROVENANCE.json` records the archive and each uncompressed file identity. The last command reads existing tables and reconciles all 528 count, hit and sign assertions; it does not fit models. The returned receipt is in `VALIDATION.json`.

The tables contain project-generated gene/program statistics, identifiers, eligibility flags and native-result status fields. They contain no individual patient data, raw reads, gene-set memberships or regulatory edges. Original project-generated tables are covered by the repository CC BY 4.0 grant; upstream resource terms remain unchanged.

## Interpretation

Read the explicit left/right native presence, eligibility, uncertainty-status and continuous-quantity fields together. Historical columns remain byte-for-byte intact for provenance; their presence does not validate every interval. In particular, rows marked `NO_VALIDATED_NATIVE_QL_INTERVAL_OLD_SOURCE_COLUMNS_RETAINED_NOT_ADOPTED` do not supply an adopted edgeR interval. Missing or non-estimable values remain missing. Gene effects, activity coefficients and enrichment scores are distinct quantities.

This is a technical quantifier sensitivity analysis of the same nine libraries, not an independent biological replication. The existing scientific conclusions, original S1 gate and manuscript uncertainty rules are unchanged. The seven-table reading supplement remains the approved reporting projection. The exported-count check verifies descriptive reconciliation, not whole-study cold restoration or every historical stochastic internal.

See the [source and count-input guide](../../analysis/g205_matched_sources/README.md) for scientific code, dependencies and exact count restoration.
