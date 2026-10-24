# Frozen JGW03 analysis inputs

This directory contains the actual selected input tables and allocation registries used by the accepted separate Linux reference workflow, recovered from the retained bridge archive. All 12 files match both the archive inventory and the final Linux input-pin manifest. They are supplied unchanged, with SHA-256 identities in `PROVENANCE.json`.

| Accession | Context | Libraries | Gene rows before analysis filtering |
|---|---|---:|---:|
| [GSE227494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227494) | C2C12 myotubes | 12 | 42,807 |
| [GSE289929](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289929) | PC12 neural-like cells | 12 | 32,883 |
| [GSE299063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299063) | HepG2 off-tissue context | 12 | 60,675 |

Each accession has a count matrix, locked sample table, aligned gene annotation and `FULL_LABEL_REGISTRY.rds`. The matrices contain **RSEM expected counts**, including fractional values. They are not TPM/CPM and must not be silently rounded for every method. Original sources define method-specific handling. Four dose groups contain three submitted libraries each. Row/column order and mapping statuses are preserved, including unmapped genes.

The sample tables preserve historical curation fields verbatim. Their `statistical_unit` or `biological_unit_id` labels do not authenticate independent culture preparation, randomized assignment, paired culture blocks or exchangeability. Concentrations describe these cell-culture exposures; they do not establish clinically attainable exposure or human efficacy. This is the selected study snapshot, not a claim to reproduce the latest GEO record.

## Validate and place

Run the byte/content check before reading the frozen R objects:

```sh
python -I -B pipeline/jgw_inputs.py
Rscript --vanilla data/jgw03_inputs/verify_labels.R data/jgw03_inputs
```

The Python check verifies all 12 file identities, count validity, annotation alignment, sample order/groups and the published saved-summary identities. The R check reads the existing registries and verifies the observed assignment, all 10,000 distinct nonidentity assignments, the complete 20-assignment high-dose companion and fixed intermediate-dose labels. It generates no new allocations and executes no model. Each registry preserves its original seed and RNG states; the two reference spaces retain their distinct tail conventions.

For the historical Linux project layout, place each accession's count/sample/annotation files under `project/data/processed/jgw03_three_cell_lines/<accession>/`, and its registry under `intermediate/<accession>/`. The exact accepted relative namespaces are recorded per file. These are inputs to the published `JGW03_input_admission_v2_20261015.R`, preparation/core and Linux-reference sources, not a one-command full reconstruction. Frozen gene-set resources, fixed/prepared objects, software environments and other runtime records remain separate requirements. Do not rerun completed fits to validate this delivery.

## Sources and rights

The selected matrices derive from the public GEO submissions linked above. NCBI itself places no restrictions on GEO data use/distribution, while submitter rights may still apply; see the [official GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). The repository's code licence does not relicense source data or annotations. Original project curation retains its stated CC BY 4.0 scope. No article PDF or third-party gene-set membership collection is included here.

`data/derived/SAVED_OUTPUT_LINEAGE.json` separately connects the already public G128/JGW aggregate tables to retained accepted summaries. All 73 G128 small-evidence identities were rechecked without reopening fitted scientific objects. This saved-output check is not a new full native replay, proof of randomization or biological holdout validation.
