# GSE277064 executed reconstruction and complete comparison

The fixed selected-input five-stage sequence completed once on Linux on 29 October 2026. The native preflight and all five steps returned zero, with unchanged pinned inputs. Complete comparison covered all 169 retained result files,125 tables,5 RDS objects,33 graphics and 24 gene-reporting families.

Reported gene coefficients, native P/q values, program effects and supported score intervals passed the fixed comparisons. The original per-column R all.equal(tolerance=1e-9) criterion was retained. All strict failures remain in the files supplied here:38 checks in quarantined historical gene-interval fields,3 run-manifest/validation checks,10 metadata storage-type changes,11 internal-fit numerical differences and 3 refresh timestamps. Selected metadata values and all 5 design/contrast matrices agreed exactly. There is no claim that every component or historical file byte is identical.

All 24 new gene views preserve their source's retained lexical cells and row order. There are2,171,342 historical/new retained lexical differences; numerical agreement is recorded separately. The 33 graphical pairs were rendered and visually checked against the original plotting code and result sources. Four legacy program files retain long-subtitle cropping; the manuscript uses separately reviewed publication graphics.

## Inspect the evidence

From the repository root, use Python's standard library:

```sh
python -I -B verification/G277_20261029/verify.py
```

This checks the distributed evidence identities, coverage and strict-difference counts. It does not rerun a scientific model. The complete table and RDS component records are supplied without fitted model binaries. Distributed text uses LF line endings; cell values and order are preserved. The original data terms and third-party rights remain in force.

## Compare an independently generated result

Use compare_tables_RDS.R with a TSV containing name,kind,baseline,rebuilt columns and a new output directory. kind is table or RDS. Preserve original row order and the fixed tolerance. To inspect edgeR/limma S4 objects, use the package versions in the branch execution contract; missing classes or opaque objects are recorded rather than silently skipped.

```sh
Rscript --vanilla verification/G277_20261029/compare_tables_RDS.R comparison-spec.tsv new-comparison
```

The five-stage execution command and the UTF-8 runtime requirement are described in analysis/g277_execution/README.md. On the Linux environment used here, set ROUTE_B2_RSCRIPT_REAL to the pinned Rscript, make rscript_utf8.sh executable, and pass its absolute path as --rscript. The wrapper enforces C.UTF-8 after the Python entry constructs its child environment. Current summary/validation sources are retained adaptations with separate identities. Whole-ledger/config identities differ from older manifests; selected model values/designs were checked directly. This is a selected-input branch reconstruction with full comparison and retained differences, not whole-study cold restoration or independent phenotype validation.
