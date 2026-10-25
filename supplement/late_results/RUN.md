# Rebuild the late-results reading package

From the extracted repository directory, with Python 3.12 or later:

```sh
python -I -B pipeline/late_supplement.py --output output/late-results
```

Use a new output directory. The command needs only Python's standard library and the included public files; it performs no network access or scientific fits. This is a saved-result reporting stage: R, renv, scientific packages, raw reads, private caches and AutoDL are not used. There is no random-number generator or seed. The repository manifest and seven-table provenance pin the inputs; the receipt records code, template, Python version and output hashes. No dependency installation or additional Python lock file is needed for this standard-library-only command.

Read `output/late-results/SUPPLEMENTARY_RESULTS.md`. The directory also contains all seven unchanged TSV files, their source identities, 440 G128 summary-field comparisons and an execution receipt. `NA`, original precision and row order are retained in the copied TSVs; only the human-readable G205 correlations/Jaccards are rounded to four decimals, matching the accepted guide.

The dependency order is: release manifest → seven-table provenance and versioned reading template → seven accepted tables → summary reconciliation → generated reading package. Missing or altered inputs stop the command before publication of its output. Existing output directories are preserved. This short reporting stage needs no computation checkpoints; prior output directories remain separate.

A1/A6 are complete family summaries; A2–A5 are source-preserved compact views checked against A1; A7 is the matched-alignment summary. These totals overlap and are not unique genes or independent biological samples. No replacement analysis, new P value, new candidate, or external validation is implied.

The original code was adapted from the retained local reading-assembly logic to remove private path dependencies. The statistical definitions and seven accepted input tables are unchanged. Validation covers a fresh archive extraction, generated-table agreement with the previously published guide, exact table copying, and execution from an unrelated working directory. Tested on Windows with CPython 3.12.14; other supported-looking runtimes are not claimed as tested. Full scientific cold restoration remains separate and unverified.

## Observed verification

The command completed from a fresh archive extraction and an unrelated working directory, with isolated Python (`-I -B`). All seven table copies were byte-identical to the public inputs. The 44 displayed data rows (4 G128 overview, 12 G128 gene, 6 JGW overview and 22 G205 comparison rows) matched the previously published reading guide. All generated reader links resolved, and 440 G128 compact-view field comparisons matched the full family table. Existing output was preserved on a second invocation; a modified copied input was rejected before creating output. Those checks test reporting and input identity, not original scientific inference.
