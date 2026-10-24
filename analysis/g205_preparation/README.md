# GSE205906 original preparation, without model fitting

This entry turns the accepted HISAT2 strand-0 nine-library input into the original
ten design/filter/TMM objects (full data plus nine single-library omissions),
the corresponding four-collection membership objects, and an inspected original
consumer handle. It calls the unchanged 6 October Ensembl101 functions. It does
not call the handle's `execute()`, any reference `refit()`, or model `compute()`.
Accepted gene/program results are reused and remain untouched.

## Inputs and order

1. Run [`prepare_count_inputs.py`](../input_bridge/README.md) to obtain
   `prepared/G205_HISAT2_s0_inputs.rds`. The entry currently accepts the exact
   authenticated prepared bytes tested under the recorded Windows environment.
2. Run [`build_resources.py`](../resource_build/README.md) to recreate the original
   frozen resource bundle. Supply that output root; seven resource files are
   authenticated, including the original gene-set manifest supplied here.
3. Use the existing R 4.6.1 library with the versions in the unchanged full-model
   configuration. All ten R package versions are checked. One thread and
   BiocParallel SerialParam are used. This command installs nothing.

```sh
python -I -B pipeline/prepare_g205_model.py \
  --prepared /path/to/input-preparation/prepared/G205_HISAT2_s0_inputs.rds \
  --resource-root /path/to/resource-build \
  --rscript /path/to/Rscript \
  --r-library /path/to/R-library \
  --output /path/to/new-g205-preparation
```

The Python entry uses only the Python 3.12 standard library. The separate
`environment/G205_core` acquisition inventory documents the historical Linux
core; this preparation does not restore it or convert the anchor's `renv.lock`
into a universal environment. The original configuration is preserved. Its
PowerShell version is historical launcher metadata, not a dependency of this
Python-launched command.

```text
authenticated prepared inputs + exact sources/config + frozen resources
  -> original design, expression filters and TMM (10 phases)
  -> original program membership and eligibility (10 phases, 4 collections)
  -> original core/gene-reference/program-reference prepare + inspect
  -> saved preparation objects and QA; STOP before any fitted analysis

new-g205-preparation/
  derived/   new component serializations, prefit.rds, mapping.rds
  receipts/  RUN.json, runtime.json, sessionInfo.txt, stdout/stderr,
             PREPARATION_QA.json, PARENT.json
```

`PARENT.json` records the command, input/source/output hashes, timing, exit code
and before/after preservation. `PREPARATION_QA.json` reports observed retained
genes and eligible programs by omission. Those are preparation diagnostics,
not new effects, P values or candidate results. No random allocations are drawn;
the original future-model seed remains 20260817. Existing output directories
are refused. Failed directories remain available for diagnosis. Objects are
written only to the new directory and read back for exact serialized-content
equality. Process-local data.table self-reference pointers are not scientific data.

## Provenance and scope

The four component files and runtime JSON are explicitly **new** bindings. In
particular, `fixed_gene_ids.rds` contains IDs derived from the authenticated
count input; it is not a reconstruction of the historical `run_manifest.json`.
The original API field called `statistical_runtime_lock` binds the observed,
version-checked preparation runtime JSON, not historical Linux binary custody.
No historical success marker or model result is fabricated. Original source
assertions and their limits remain intact.

`gene_set_manifest.tsv` is the unchanged original 32,072-byte project manifest
(SHA-256 `b3a9089a6e136c5271cf6ac3999d10304dc2316fafd10f5c170d9f43084bbce1`).
It lists historical relative file identities. It neither redistributes the
listed third-party resources nor asserts that every historical file is part of
this repository. Resource reconstruction and upstream terms still apply.

The observations are nine deposited libraries. Animal pooling, batch and
administered dose remain unresolved. This preparation does not establish
independent biological replication, a spasticity holdout or therapeutic effects.

## Verification and remaining work

`VALIDATION.json` records actual tests of this entry. Additional read-only input
rejection checks can be run against a successful preparation directory:

```sh
Rscript --vanilla analysis/g205_preparation/negative_checks.R /path/to/repository /path/to/completed-preparation
```

Use the same explicit library and UTF-8 locale as the main command. A complete
Linux model-execution adapter, numerical model/refit comparison through that
adapter, and cold environment restoration are still untested. Reusing these
preparation objects alone is not proof of a complete historical model rerun.
