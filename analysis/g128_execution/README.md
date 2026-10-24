# GSE128623 original analysis entry

This entry connects the published accepted count inputs and reconstructed frozen
resource bundle to the original 30-library and prespecified 29-library analysis.
It does not select new samples, generate new labels, or alter the accepted study
results. An explicit `run` command performs real fitting and should only be used
when a researcher intends to reproduce those analyses. Routine inspection uses
`check`; the project release test uses generated `synthetic` data.

## Inputs and execution order

```text
accepted gzip inputs --count_inputs.py materialize--> exact count/annotation/ledger files
eight lawful frozen source files --build_resources.py--> exact frozen resource RDS
count files + resource RDS + retained S1 gate + pinned R library
  --g128_model.py check--> original-reader, runtime and member checks (no fits)
  --g128_model.py run--> original analyzer --> independent validator --> comparison check

generated counts + generated memberships
  --g128_model.py synthetic--> same statistical kernels and independent validator
```

Use a new output directory for every invocation. Input files are read-only and
their identities are checked again after execution. The only copied input is a
separate byte-identical frozen-design ledger, required by the original driver's
distinct-path interface. Outputs remain under the new directory:

```text
run/
  frozen_design.tsv             # check/run only
  synthetic_inputs/             # synthetic mode only; visibly artificial genes/programs
  model_output/                 # absent in check mode
    primary_all30/
    complete29_sensitivity/
    comparison/
    ANALYZER_COMPLETE.tsv
    VALIDATION_REPORT.tsv
  receipts/
    RUN.json, COMMAND.json, stdout.txt, stderr.txt
    sessionInfo.txt, loaded_packages.tsv, R_VALIDATION.json, COMPLETE.json
```

## Environment

Python 3.12+ uses only its standard library; there is no Python dependency install.
Use an explicit R library. R 4.6.1 and the exact original versions are required:
data.table 1.18.4, edgeR 4.10.1, limma 3.68.4, DESeq2 1.52.0,
fgsea 1.38.0, jsonlite 2.0.0, yaml 2.3.12 and openssl 2.4.2, plus
their transitive dependencies (including BiocParallel). The anchor
[`renv.lock`](../GSE16710/renv.lock) records the existing Windows library;
the [Linux package guide](../../environment/G205_core/README.md) documents the
available core package archives and installation sequence. Neither proves that
an arbitrary new Linux installation has been restored. `sessionInfo.txt` records
the actual loaded environment and linked numerical libraries.

Formal execution retains the original Linux requirement. Windows supports the
original runtime's explicit test-only check for **check/synthetic modes only**.
Single-worker scheduling bounds desktop resource use; no statistical function,
contrast, filtering rule or significance threshold is replaced. The original
9,999 rotations, minimum 10 members, seed 20260819 and RNG kinds are fixed.
The separate count-fixture seed is 20261024. UTF-8 locale is set per subprocess;
global user settings are not changed.

## Commands

From the repository, after restoring counts and building the resources using the
two preceding guides:

```sh
python -I -B pipeline/g128_model.py check \
  --counts-root /work/count_inputs \
  --resources /work/resources/data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds \
  --cohort adapter --rscript /opt/R/bin/Rscript --r-library /work/R-library \
  --output /work/g128_adapter_check

python -I -B pipeline/g128_model.py synthetic \
  --rscript /opt/R/bin/Rscript --r-library /work/R-library \
  --output /work/g128_synthetic
```

The `synthetic` command rejects real-data input arguments. Its 360 generated gene
rows include ten all-zero rows, and its five artificial programs cover all four
resource collections. Accession-shaped sample labels are used solely because the
unchanged original kernel checks the fixed omitted sample/run identifier. No
generated value is a measurement from those public samples.

To intentionally reproduce the real model sequence **on Linux**, replace `check`
with `run`, retaining the input arguments and specifying a new output directory.
Choose `--cohort raw` or `--cohort adapter` explicitly; each invocation covers both
30 and 29 libraries. Do not execute these commands merely to inspect accepted
results. No real-data `run` is performed to validate this release.

The formal entry calls the authenticated original copy-isolated driver factory.
That historical compatibility layer changes exactly one validator binding to
avoid mutating expected-table order during comparison. Original scientific bodies
and comparison tolerances remain unchanged. The old source files' standalone CLI
still refuses its superseded incomplete contract; this entry uses their documented
function loader and the later driver, not a modification of that old CLI.

## Coverage and remaining work

The synthetic sequence exercises edgeR, limma-voom, DESeq2, activity scores,
CAMERA, mroast, fgsea, and gene/program leave-one-animal-out refits for both cohort
sizes, followed by independent reconstruction and the branch-comparison guard.
Its numerical validation uses the original tolerances. A successful synthetic
run is not a reproduction of actual scientific effect estimates or a Linux
platform check. Consult `VALIDATION.json` for what was executed for this release.

After a successful synthetic run, the separate negative checks deliberately
alter copies or in-memory values and verify rejection of an altered source,
wrong sample order, conflicting membership weight and altered comparison effect:

```sh
Rscript --vanilla analysis/g128_execution/negative_checks.R /path/to/repository /work/g128_synthetic
```

Use the same explicit R library via `R_LIBS_USER`. This writes only a new
`negative_checks/` directory under the synthetic run and performs no model fits.

The portable parent records commands and exact input/output hashes, but it does
not recreate the historical native held-file/custody system. Driver result flags
are preserved rather than turned into a claim of whole-study acceptance.
The 30/29 observed analyses are distinct from the much larger conditional-label
reference calculation; this entry does not launch or rerun that reference.
There is no mid-model resume mechanism in this original observed-model driver.
Interrupted synthetic tests may be repeated in a fresh directory; completed real
outputs should be reused instead of automatically rerunning them.

Whole-study cold-environment restoration, other branches' portable execution and
an independent phenotype-linked post-stroke-spasticity holdout remain unresolved.
The S1 gate still admits zero programs; this entry does not change candidate rules
or justify clinical efficacy or an antispasticity mechanism.
