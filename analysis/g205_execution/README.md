# GSE205906 primary scientific execution entry

`pipeline/g205_model.py` connects the existing authenticated HISAT2 strand-0
input preparation to the unchanged original checked consumer. Paths are supplied
explicitly; the nine-library input, configurations, seeds, filtering, tests and
omission definitions are not changed. This is the primary Ensembl101 analysis,
not the matched STAR sensitivity or the separate label-allocation reference.

The entry has been checked with a real-input **plan only**, syntax checks and
small synthetic output/control-flow tests. No real-data fitting has been run
through this new adapter. End-to-end numerical equivalence to the accepted
historical results, including on a freshly restored Linux environment, remains
to be tested. `VALIDATION.json` distinguishes these scopes.

## Inputs and command

Use the exact `G205_HISAT2_s0_inputs.rds` from
[the count input bridge](../input_bridge/README.md), the unchanged output of
[the resource builder](../resource_build/README.md), and R 4.6.1 with the original
ten configured package versions. The prepared RDS is currently bound to the
tested Windows serialization bytes. Recreating a semantically equivalent RDS
with different bytes is not silently accepted. The gene-set manifest comes
from the release; all seven resource identities and original source/configuration
identities are checked. No packages are installed by this entry.

```sh
python -I -B pipeline/g205_model.py \
  --prepared /path/to/input-preparation/prepared/G205_HISAT2_s0_inputs.rds \
  --resource-root /path/to/resource-build \
  --rscript /path/to/Rscript \
  --r-library /path/to/R-library \
  --output /path/to/new-g205-plan
```

This default command checks input/source file identities and creates
`receipts/RUN.json` and `receipts/EXECUTION_PARENT.json`. It starts **no R process**
and does not establish runtime package compatibility. It returns
`PLAN_ONLY_NO_R_PROCESS`. To perform a new scientific execution, add `--execute`
and supply a different, new output directory. That flag performs real model
fitting and the original independent numerical reference refits. It must not be
used merely to inspect already accepted results.

## Execution and outputs

On explicit execution, `execute.R` evaluates the unchanged, hash-bound
[`prepare.R`](../g205_preparation/prepare.R). It creates the original ten
design/filter/TMM phases (full data plus nine single-library omissions), four
resource collections and the original consumer handle. It then calls
`handle$execute()`: the original independent gene reference, independent program
reference, core calculation, comparisons and finalization. The entry uses one
thread and the original serial BiocParallel backend. It does not implement a
second version of the statistical functions.

The original `gse205906_ens101_worker_write_plan_v1` serializes all 15 tables,
the complete checked R object and its file manifest. `write_result.R` writes
those exact planned bytes to an ordinary new directory and verifies each hash
by readback. It does not emulate historical Linux file descriptors or create
historical custody/completion receipts.

```text
new-g205-execution/
  derived/                 new component, prefit and mapping serializations
  results/                 15 TSVs, checked_payload.rds, file_manifest.json
  receipts/
    RUN.json               input/source bindings and chosen paths
    PREPARATION_QA.json     diagnostics of the preceding no-fit preparation stage
    EXECUTION_QA.json       original consumer receipt and output readback
    EXECUTION_PARENT.json   exit status, preservation checks and output hashes
    runtime.json, sessionInfo.txt, execution.stdout.txt, execution.stderr.txt
```

The final parent `PASS` requires successful original consumer/output readback,
zero R exit status and unchanged authenticated sources/inputs. It is acceptance
of this **new execution**, not proof of historical numerical equivalence. The
original file manifest deliberately retains its pending-parent-acceptance flags;
new adapter acceptance is reported separately. The preparation QA's zero-fit
counts describe that earlier stage, not the subsequent explicit execution.
Existing output directories are refused; failed execution directories remain
available for diagnosis.

The output tables are `gene_edger`, `gene_voom`, `gene_deseq`, `unsigned_camera`,
`unsigned_mroast`, `unsigned_fgsea`, `activity_results`, `mapping_coverage`,
`loao_long`, `loao_summary`, `design_out`, `contrast_out`, `filter_out`,
`annotation_metrics` and `gate`. TSV list columns retain quoted JSON vectors;
the checked RDS preserves their original classes and values.

The observation unit remains nine deposited libraries. Pooling, batch and
administered dose remain unresolved. Executing the frozen analysis does not
resolve those design limitations or establish an independent spasticity holdout.

## Lightweight adapter tests

These tests create only small synthetic objects; they perform no scientific fit.
Use a new output directory for each test and an R library containing jsonlite
and digest. Set a UTF-8 locale and `R_LIBS_USER` explicitly as appropriate.

```sh
python -I -B analysis/g205_execution/test_launcher.py /path/to/repository /path/to/new-launcher-test
Rscript --vanilla analysis/g205_execution/test_writer.R /path/to/repository /path/to/new-writer-test
```

The launcher test uses a mocked R subprocess and synthetic source/input files.
Its simulated `PASS` is a control-flow test, not a model result. The writer test
uses the unchanged original serializer on tiny synthetic tables and checks all
17 file hashes, RDS/TSV roundtrips, refusal to overwrite, wrong-hash rejection and
path traversal rejection. The separate real-input plan checks the actual frozen
inputs without rerunning the already completed preparation.
