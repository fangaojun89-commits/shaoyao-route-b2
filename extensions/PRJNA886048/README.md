# PRJNA886048: complete injured-spinal pharmacological context

This candidate contains the complete fixed six-baclofen versus six-saline analysis: 12 female mouse spinal-cord libraries, 17,013 tested genes, every program test and coverage row, 105 signed-score contrasts, per-library scores, and all 266 comparisons across 19 retained contexts. It is retrospective acute-SCI pharmacological context, not a protected post-stroke spasticity validation cohort. Source article: https://doi.org/10.1016/j.dib.2022.108809.

## Verified entry points

Read the retained results with Python standard library and base R:

```sh
python -I -B pipeline/verify.py --rscript Rscript
```

The fixed-input analysis entry, path-only configuration, restored dependencies and checkpoint recovery are documented in [REPRODUCE.md](REPRODUCE.md). It authenticates 63 accepted input/resource files. The three original scientific R modules and descriptive comparison script remain byte-identical to their executed sources. Both the installed-environment and newly restored-environment checkpoint exports matched all 13 complete TSV tables: exact text and row order, numeric tolerance 1e-9. All four existing model/import/program objects were reused, with zero import, gene-fit, CAMERA or ULM executions.

[VALIDATION.json](VALIDATION.json) records a real completed-run reuse and five synthetic rejection tests: existing output, forbidden statistical override, wrong resource, corrupted completed output and incomplete import checkpoint. Original accepted files and model objects were unchanged.

[FULL_ANALYSIS_VALIDATION.json](FULL_ANALYSIS_VALIDATION.json) records a subsequent real `--native-record-checks --run` from a clean committed export with the restored environments. Import, gene fitting, CAMERA and signed activity each executed once, with zero adopted or reused checkpoints. All 13 frozen tables matched in text/row order and within numeric tolerance 1e-9. A completed `--resume` then verified 44 data/checkpoint files without refitting; all 59 files in that validation run stayed byte-identical. This used the 12 retained quantification bundles, not newly processed raw reads.

## Reproduction plan and layout

```text
data/raw_ledger/       ENA accessions, URLs, sizes and MD5; raw reads remain external/read-only
data/derived/          frozen sample identity and comparison basis
data/output/           complete retained results and QC
analysis/original/     exact historical scientific modules and drivers
pipeline/             portable input checks, checkpoint execution and result validation
config/               example input paths; no statistical overrides
environment/          renv lock, dependency contracts and labelled session evidence
provenance/           input/checkpoint identities and source bindings
supplement/           complete methods and result context
pyproject.toml/uv.lock Python comparison dependency closure
MANIFEST.json         immutable payload hashes and sizes
```

```text
accepted quantification + frozen annotation/resources + sample identity
 -> pinned input/runtime checks -> import -> gene model
 -> CAMERA + signed ULM -> treatment contrasts -> 266 descriptive comparisons
 -> full-content readback + comparison to frozen tables -> output manifest
```

The execution output additionally separates `data/intermediate`, `data/derived`, `data/output`, `checkpoints`, `logs` and `provenance`. Raw inputs and accepted science directories are never overwritten. The fixed seed, explicit versions, actual commands, stage counters and output hashes are recorded. Completed checkpoints are verified before reuse; incomplete checkpoints stop for inspection.

## Environment and limits

`uv sync --locked --no-install-project` successfully created a new Python 3.12 environment. R 4.6.1 restoration to a new target library succeeded, and all 52 dependencies resolved only to that target or the R distribution's base/recommended packages. The old project library was not used during validation. This was a fresh project environment on the same Windows host, not a new operating-system installation. The original partial gene-export session and the later full namespace snapshot are explicitly labelled.

Fresh scientific-model execution from retained quantification inputs and completed-run recovery have now passed as recorded above. Full reconstruction from raw reads has not been demonstrated. Historical raw preparation sources retain workstation-specific paths. A separate [canonical native record/content check](pipeline/native/README.md) validates fixed inputs, arguments, orientation and actual outputs while allowing internally consistent path/time changes. Its bounded tests cover nine canonical historical libraries and invented relocations/errors; it is not remote-host attestation or a fresh native run. The analysis entry keeps historical identity checking by default; --native-record-checks binds actual validated files and forbids historical checkpoint adoption. Frozen-file annotation, sequence and membership reconstruction is integrated in pipeline/reference/run.py. A [portable native entry](pipeline/native/RUN.md) accepts explicit input/tool/index/output paths and retains the executed scientific body; its 17 bounded planning checks passed, but Linux preflight and actual native execution remain unverified. Freshly generated native inputs, raw-to-result execution and final release review remain outstanding. The candidate has not been tagged as a release, and the present validation does not establish complete end-to-end reproduction.

## Interpretation and licensing

No gene crossed BH q<0.05 in either gene method. CAMERA yielded 65 Reactome and four DoRothEA membership findings; signed footprint contrasts yielded none. The 19-context descriptive panel does not show consistent general reversal. These findings do not establish inactivity, equivalence, clinical efficacy or a mechanism. Companion contrasts reuse samples.

Original project code follows the parent repository MIT license; original reports and aggregate results follow its CC-BY-4.0 terms. External source data and resources retain their applicable terms. No gene-set memberships, article files, model binaries or cloud credentials are distributed. Public availability is not treated as permission to relicense third-party material.
