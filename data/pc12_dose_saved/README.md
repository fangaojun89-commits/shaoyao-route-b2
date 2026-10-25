# Retained PC12 results for dose-encoding comparison

These 16 tables contain the exact field strings needed to compare the previously fitted ordinal and log-dose models for GSE289929. They cover the dose key, three gene-result tables and paired program/activity results for six additional families. They contain 46,307 table rows in total; this is not a sample count. The study design contains 12 libraries, three at each of four dose levels.

[PROVENANCE.json](PROVENANCE.json) identifies each original export by project-relative path and SHA-256, lists the selected columns, and records the projected file size/hash. Selection preserves row order and every retained field string. Repeated source gene-ID columns were checked for equality before selecting one. No numerical rounding, gene filtering, fitting, P-value calculation or multiple-testing adjustment was performed during projection. Gzip files use deterministic compression metadata.

The gene projections include only the gene identifier, coefficient, P value and FDR. Unsupported historical `ql_approx_*` fields are omitted. The program-score/activity projections retain their recorded nominal intervals. CAMERA tables retain direction and significance fields; no artificial CAMERA effect coefficient is supplied. No gene-set memberships, network edges, software archives or article PDFs are included.

These are project-generated statistical results derived from the public [GSE289929](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289929) submission. Original result tables/documentation are offered under the repository's CC BY 4.0 scope; upstream source data and annotations retain their provenance and applicable terms. Resource identifiers do not grant redistribution rights to the resource memberships. See [third-party notices](../../THIRD_PARTY_NOTICES.md).

[Run the comparison](../../supplement/pc12_dose/RUN.md). The original producing R analysis is retained [here](../../analysis/study_sources/code/R/analyze_jgw03_three_cell_lines.R); running it is not necessary for this saved-result comparison.
