# JGW03 original input reader and portable unfitted bundles

This command connects all three frozen JGW03 count/annotation/sample tables and
their existing allocation registries to the unchanged original input reader and
resource checks. It stages the recorded relative layout and saves one portable
unfitted input bundle per accession. It reads the original kernel definitions,
but does not call fixed-layer construction, ULM scoring, model fitting, a label
allocation generator, or a reference dispatcher. Completed analyses are reused.

First rebuild the frozen resources using [build_resources.py](../resource_build/README.md).
Then use an existing R 4.6.1 installation with digest and jsonlite:

```sh
python -I -B pipeline/prepare_jgw_inputs.py \
  --resource-root /path/to/completed-resource-build \
  --rscript /path/to/Rscript \
  --r-library /path/to/R-library \
  --output /path/to/new-jgw-input-preparation
```

The explicit library is supplied only to the child process. No packages are
installed and no existing environment is modified. Python 3.12 or later is
required. Use a new output directory on a volume with at least 100 MB free.
The staged resources retain upstream terms and are not added to this repository.

```text
new-jgw-input-preparation/
  project/data/processed/  exact original table/resource layout
  project/intermediate/   copies of saved registries
  prepared/               three *_unfitted_inputs.rds bundles
  receipts/               invocation, runtime, source/input pins and readback QA
```

Each saved bundle includes fractional RSEM expected counts, aligned annotation,
the original locked sample metadata, original resource objects, portable source
records and the unchanged label registry. It is a new serialization of those
inputs, not an original historical fixed-layer or model checkpoint. Individual
sample scores, dispersion estimates, effects, P values and q values are absent.
The original bound fixed-layer function remains available for a separately
requested reconstruction; it would compute TMM and label-independent ULM scores.
This release deliberately does not execute that function.

The reader authenticates eight files per accession, including the resources.
Saved 10,001-row and 20-row label tables must align with the input column order.
The command verifies original sources and inputs before/after use, saved-object
content equality, and a separate direct TSV count/annotation/metadata readback.
These are input/preparation checks, not evidence of independent cultures,
randomized treatment allocation, biological validation, or a complete Linux
restoration. Original source curation and [upstream rights](../../THIRD_PARTY_NOTICES.md)
remain unchanged.

The supplied `negative_checks.R` exercises the original adapter's rejection of
a substituted formula ratio, duplicate sample ID and wrong resource taxon in
temporary memory. Use the same R library and locale as the main command:

```sh
Rscript --vanilla analysis/jgw_preparation/negative_checks.R /path/to/repository /path/to/completed-preparation
```
