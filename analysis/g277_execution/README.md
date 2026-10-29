# GSE277064 scientific execution and accepted-result reuse

This entry supplies the five-stage main G277 sequence: cell/time analysis,
deterministic fgsea refresh, exploratory summary, program figures and output
validation. It stages the three count matrices, both explicit metadata selections,
three configurations, frozen resource RDS and all source dependencies into a new
project. The earlier missing summary, validator and narrative-template dependency
are now supplied. `check` never evaluates a scientific source. `report-saved`
reuses the 169 pinned accepted result files and runs only summary/figures/validation
and the current gene reporting projection. It does not repeat normalization,
scores, fits, rotations or fgsea. `run` is for an intentional future reproduction
from the processed counts, with all five stages, and was subsequently executed once on Linux on 29 October 2026. The [complete comparison](../../verification/G277_20261029/README.md) records all 169 retained files and 24 gene-reporting families, including strict differences. Completed historical project analyses are retained.

## Commands

Use Python 3.12+ (standard library only), R 4.6.1 and the explicit R library whose
direct package versions are listed in `CONTRACT.json`. The existing Windows
environment is recorded in the executed preflight. The anchor's `renv.lock`
provides package provenance; it is not proof of a newly installed environment.
No package is installed automatically. On Linux, verify that en_US.UTF-8 exists before using the direct entry. The actual 29 October run used an explicit Rscript wrapper exporting LANG, LC_ALL and LC_CTYPE as C.UTF-8 before exec of the pinned Rscript, because en_US.UTF-8 was unavailable. The wrapper changed runtime encoding only; the scientific code and parameters were preserved. The code explicitly sets subprocess UTF-8,
the recorded Asia/Shanghai timezone and one computational thread. The original
seed 20260817, 9,999 rotations and primary fgsea flag TRUE are retained; the
subsequent original fgsea refresh has its own original tie-handling/seed logic.

First build the exact resource RDS using [resource reconstruction](../resource_build/README.md).
Then, from any directory:

```sh
python -I -B /repo/pipeline/g277_model.py check \
  --resources /work/resources/data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds \
  --rscript /opt/R/bin/Rscript --r-library /work/R-library --output /new/g277_check
```

This checks source/input identities, all direct package versions, syntax of all
five scripts, graphics availability, configuration and frozen resource versions.
Original count contents/designs and saved preprocessing were already tested by
the [input](../g277_preparation/README.md) and
[saved-input](../g277_saved_inputs/README.md) commands. The new check uses their
same pinned files; it creates no scientific `results` directory.

For intentional future production reproduction, replace `check` with `run`, use
another new directory and keep sufficient memory/disk available. This full
primary sequence can be lengthy. Native failures stop the sequence and retain
logs; it never restarts or replaces an existing directory automatically. A
zero exit from the historical validator establishes its enumerated checks,
not numerical agreement with every accepted result. A new environment or real
run still needs a direct comparison before adopting new scientific results.

For readback/reporting of the supplied project's accepted files:

```sh
python -I -B /repo/pipeline/g277_model.py report-saved \
  --saved-results /retained/results/GSE277064 \
  --rscript /opt/R/bin/Rscript --r-library /work/R-library --output /new/g277_report
```

`ACCEPTED_OUTPUTS.json` specifies all 169 required files and identities. These
large fitted/result files are not redistributed here; this mode requires the
retained file set. No model is called. Inputs are copied, the original source
directory remains read-only, and only the new copy is processed. The old
manuscript-ready file is left at its source and recreated at the new destination
from its exact frozen template, as required by the successor summary. That
template is a **historical narrative**, not a current inferential authority.

## Outputs and current reporting rule

```text
new-output/
  project/                 exact inputs/code, original historical result namespace
  receipts/                staging, commands, actual exits, package/session records
  reporting/
    view/                  current 15-column gene reporting tables (24 families)
    quarantine/            the three excluded legacy fields, traceability only
    schema/                exact field/source/row reconstruction bindings
    PROJECTION.json        comparison and lexical-preservation evidence
```

The original analysis emits `ql_se`, `ci_low` and `ci_high` at gene level. Those
legacy quantities do not establish standard edgeR QLF coefficient uncertainty.
The approved original pure reporting codec omits precisely these three columns
from the reader-facing gene tables, preserves all other strings and row order,
and retains excluded values in separate traceability files. `signed_sqrt_F`
remains a ranking statistic. This code does not infer a coefficient SE from it.
Method-specific activity-score intervals are untouched. The reporting projection
is applied after a future `run` as well as after `report-saved`; it never repairs
the historical files in place. The old template's gene-CI sentence remains
historical and is explicitly superseded by this rule.

The current summary and validator sources are retained local successors.
Their identities and the historical producer identities are both recorded in
the contract. They must not be described as byte-identical historical producers.
The three original statistical/plot sources retain their previously published
bytes. The reporting codec is also unchanged. The runtime receipt identifies
the actual scripts executed, without inventing historical native custody.

`report-saved` must reproduce every saved view/quarantine/schema hash and checks
every retained lexical cell independently. This is the approved reporting rule
behind S14/S16, not a reconstruction of Excel formatting or the whole paper.
The separate limma-voom/DESeq2 sensitivity, later presentation figures and other
branches are outside this five-stage primary sequence. Their accepted results
remain in place; this release does not claim complete cold reconstruction.

The retained S1 gate is checked explicitly (1,771 rows, zero admissions), because
the historical G277 sequence encodes an empty upstream candidate set. Source-wide
exploratory results cannot create an S3 admission or establish a PSS phenotype.
