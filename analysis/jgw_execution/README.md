# JGW03 original primary execution entry

`pipeline/jgw_model.py` connects the unchanged original primary analysis and dose-summary scripts to their recorded GEO deposits, metadata, fixed configuration, resources and zero-admission S1 result. It provides a separate `check` command that runs only the original metadata and deposited-count readers. This closes the original-source input path in addition to the existing [prepared-table reader](../jgw_preparation/README.md).

## Check without repeating a scientific analysis

Use Python 3.12+ (standard library only), R 4.6.1 and the explicitly supplied existing R library. Direct package versions are checked against `CONTRACT.json`. The original primary Windows session and package table are retained under `sources/`; these are runtime evidence, not a complete portable environment lock. The helper requires jsonlite 2.0.0. The command installs nothing.

Place the following three original GEO deposits under an input root, preserving their bytes. Exact identities and official source URLs are in `CONTRACT.json` and the [acquisition catalogue](../../data/acquisition/INPUTS.json).

```text
input-root/data/raw/GSE289929/GSE289929_RAW.tar
input-root/data/raw/GSE227494/GSE227494_JGT_RCM_C2C12.txt.gz
input-root/data/raw/GSE299063/GSE299063_RAW.tar
```

These total 255,320,349 bytes and contain source-deposited processed RSEM results, not sequencing reads. They are not redistributed by this extension. Supply the unchanged eight resources from a completed [resource build](../resource_build/README.md).

```sh
python -I -B pipeline/jgw_model.py check \
  --input-root /path/to/input-root \
  --resource-root /path/to/completed-resource-build \
  --rscript /path/to/Rscript \
  --r-library /path/to/R-library \
  --output /path/to/new-jgw-input-check
```

Use a new directory on a data volume. The check stages approximately 300 MB of source files and resources, plus selected extracted archive members and receipts; allow at least 500 MB. It does not overwrite the source root or reuse an unverified output directory. The count reader selects the exact original archive members/columns and compares all counts, row/column order and locked sample fields with the accepted public input tables. Original GEO metadata checks retain their limitations. Historical metadata labels do not independently prove biological replication or source curation chronology.

The check parses the original program but binds only eight input/helper functions and the dataset specification. It does not bind or call `analyze_one`, construct a fixed layer, normalize, score, fit, draw labels or run a conditional reference. It writes a receipt, package information and selected input copies, not scientific results. `receipts/INPUT_CHECK.json` describes this finite check; `receipts/COMPLETE.json` includes source/input preservation and command records.

## Explicit future primary reconstruction

Replacing `check` with `run` performs the same input check and then executes these two unchanged original scripts, in order, in the fresh staged project:

1. `analyze_jgw03_three_cell_lines.R`: three separate culture analyses and their exploratory context summaries, with the original seed, configuration and S1 stopping rule.
2. `summarize_jgw03_dose_response.R`: descriptive four-level summaries from the newly generated model objects.

**The `run` branch performs scientific computation. It was not executed for this release.** Its availability is not proof that a new operating system, fresh library or complete study reconstruction has passed. Production time and peak resources have not been measured for this portable entry. Original accepted results were not rerun to test it.

The source script writes its historical primary output schema. This is not the whole current manuscript bundle: separate rotation reconstruction, 89/34-family Linux conditional references, later reporting projections and S3/S4 workbook changes remain separately versioned. Those stages must not be inferred from successful primary-script exit. The original broader output validator additionally expects a historical narrative report and later dose/report manifest products; it is not silently represented as a validated two-script pipeline. Do not substitute new output files for accepted results without scientific comparison.

## Preservation and source provenance

Raw source deposits, public source scripts, metadata and frozen resources are checked before staging and after execution. The original scripts receive a new project root; their resume environment is explicitly disabled to prevent accidental use of a partial checkpoint. Accepted results are outside that root.

The source metadata projections retain all fields, values and order for the selected accessions and 36 sample identifiers: 1,320 GEO metadata rows and three dataset-ledger rows. `PROVENANCE.json` identifies the complete original files and the exact projection rule. The earlier 99-row submitted-unit ledger, fixed configuration and S1 table are reused unchanged. Source metadata and count deposits retain their upstream terms; this extension introduces no new licensing permission over them.

See `VALIDATION.json` for actual executed checks. No cloud instance, complete model fit or new allocation was needed to prepare this entry.
