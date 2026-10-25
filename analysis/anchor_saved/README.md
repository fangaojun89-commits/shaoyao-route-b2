# Reuse and check the completed experimental-spasticity anchor

This entry stages 30 retained input roles for the unchanged original
`validate_GSE16710_outputs.R`, runs its 51 output checks, then separately checks
all 1,771 recorded admission decisions against their component fields and
46,046 saved omission estimates. It reads stored permutation P/q joins but
does not draw labels, recompute activity scores or fit any model.

The gate admitted zero programs. This result stops the candidate-restricted
chain; it is not a test showing that the formula lacks clinical efficacy.

## Command and environment

Use Python 3.12 or later and an existing R 4.6.1 installation with data.table
1.18.4, jsonlite 2.0.0 and R.utils 2.13.0. Their dependencies must also be
available. The original [anchor lock](../GSE16710/renv.lock) records the full
project R inventory; this three-package check is not a complete environment
restoration. No Python third-party package is required. No package is installed
and no existing library is modified by this command.

```sh
python -I -B pipeline/anchor_saved.py \
  --rscript /path/to/Rscript \
  --r-library /path/to/R-library \
  --output /path/to/new-anchor-saved-check
```

The output must be new; allow 20 MB on the selected volume. Inputs and original
code remain immutable. The child uses an explicit library, UTF-8 locale, one
thread and the canonical retained QC record. Actual command, R session,
input/source pins, exit status and output hashes are recorded in `receipts/`.
An isolated output tree avoids writing into an accepted analysis directory.

```text
repository/data/anchor_saved_outputs/  retained derived tables and provenance
                 |
                 v
new-anchor-saved-check/project/       exact staged reader layout
                 |
                 +--> original 51-check R validator
                 +--> independent saved-gate/omission readback
                                      |
new-anchor-saved-check/receipts/       results, runtime, commands and identities
```

These inputs are previously generated results, not raw arrays. Raw-to-processed
and raw-to-model reconstruction are separate stages described in the
[anchor source guide](../GSE16710/README.md). No checkpoint or model is resumed
by this command. A second new directory can be used to check portability; the
original check table is deterministic, while timestamps and paths in receipts
are expected to differ. Complete clean-machine/Linux restoration remains
unverified.

## Content and interpretation

All 30 source roles match the retained internal output manifest. Twenty-eight
roles retain original bytes. Two fgsea result tables omit only `leadingEdge`;
every other field string and row order is preserved. This keeps derived
gene-membership lists outside the public selection while supplying every field
read by the original validator. The provenance gives original/public identities
and the explicit transformation. Some roles reuse already-public QC/gate files;
they are not duplicate scientific observations. No third-party source network,
gene-set membership or article PDF is added.

The separate reader verifies complete family keys, CI/sign/threshold gate
arithmetic, omission extrema, direction fractions, sign-flip counts, relative
changes and gate-summary/permutation joins. It does not regenerate original
effects, intervals, P values, q values or omitted-sample fits. The original R
validator also checks BH arithmetic from retained P values.

Several historical check names refer to unique biological units or outcome-blind
curation. Passing those checks verifies the stored identifiers/fields. It does
not authenticate animal independence, actual blinding or pre-outcome chronology.
The saved spasticity anchor is a spinal-injury model, not an independent stroke
cohort or a clinical validation set.

Original code is MIT; original curation/results are CC BY 4.0. Source-derived
facts and upstream resource identifiers retain their provenance and terms.
See [resource notices](../../THIRD_PARTY_NOTICES.md).
