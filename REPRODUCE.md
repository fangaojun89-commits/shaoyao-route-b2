# Reproduce and inspect Route B2

Start with the operation you need. All paths below are relative to the repository root. Use a fresh output directory and retain the version tag and `MANIFEST.json` with your results.

## 1. Read the completed scientific results

| Question | Complete public result or guide |
|---|---|
| Which programs passed the experimental-spasticity gate? | [Saved anchor results and original validator](analysis/anchor_saved/README.md): 1,771 gate decisions and 46,046 saved omission estimates |
| How do processing, cohort and conditional-label specifications affect post-stroke spinal findings? | [GSE128623 complete family tables and reading notes](supplement/late_results/README.md): 448 families |
| How stable are the whole-formula culture results across conditional label assignments? | [JGW03 Linux complete reference summary](supplement/late_results/README.md): 89 four-dose and 34 high-dose/vehicle families |
| How does the same-library alignment choice affect GSE205906? | [Complete STAR/HISAT2 exports](data/g205_matched_exports/README.md): 216,070 family-feature records across 22 families |
| What depends on extraction, solvent or herb ratio? | [GSE298414 preparation results](supplement/preparation_extension/RUN.md): 72 libraries and 36 stratified contrasts |
| How do mixtures compare with the single herbs? | [GSE298414 constituent comparisons](supplement/constituent_extension/RUN.md): 126 libraries and 60 contrasts |
| What was actually reconstructed and compared? | [Current branch-scope ledger](verification/RECONSTRUCTION_LEDGER_20261029.md), including G277 complete comparison and PRJ retained-input execution |
| Which deposited units were used? | [Dated submitted-unit ledger](data/ledger/README.md), supplemented by each later branch's sample and provenance tables |

The anchor admitted 0/1,771 programs. Source-wide exploratory results are reported in full; candidate-restricted downstream tests remain unevaluated. GSE298414's two analysis blocks share 72 libraries. STAR and HISAT2 describe the same nine GSE205906 libraries. These tables therefore represent distinct analytical views, rather than independent biological replications.

## 2. Recheck saved results without refitting

For the small GSE128623/JGW03 descriptive summaries, Python 3.12 or later and its standard library suffice:

```sh
python -I -B pipeline/sources.py verify
python -I -B pipeline/report.py --g128 data/derived/G128_FAMILY_SUMMARY.tsv --jgw data/derived/JGW_FAMILY_SUMMARY.tsv --output new-summary
```

Compare the three descriptive outputs with `expected/`; `RECEIPT.json` records their identities. The [anchor guide](analysis/anchor_saved/README.md) adds an R-based saved-result validator. The [matched-export guide](data/g205_matched_exports/README.md) gives the exact extraction and 528-assertion reconciliation commands. The preparation and constituent guides above supply their own numpy/pandas readbacks, including retained library-score effects. These operations use existing estimates; they do not repeat the original scientific fits.

## 3. Reconstruct inputs and choose a model branch

Obtain separately distributed resources under their upstream terms. Follow the [acquisition catalog](data/acquisition/README.md), [count restoration](data/count_inputs/README.md), [original count-input readers](analysis/input_bridge/README.md) and [frozen resource build](analysis/resource_build/README.md) in that order where applicable. Source archives and dependencies are not automatically downloaded by these commands.

| Branch | Entry and non-fitting operation | Explicit scientific execution and current validation |
|---|---|---|
| GSE16710 anchor | [Original sources and dependency order](analysis/GSE16710/README.md); `pipeline/anchor_saved.py` checks accepted results | Raw-to-model reconstruction remains a separate task; the saved check validates the completed gate |
| GSE128623 | [`pipeline/g128_model.py check`](analysis/g128_execution/README.md), selecting `raw` or `adapter` counts | `run` invokes the original all 30/predefined29 observed-model sequence; generated-data kernels and real-input checks passed. Portable real-model execution and later conditional-reference orchestration are separate scopes |
| JGW03 | [`pipeline/jgw_model.py check`](analysis/jgw_execution/README.md) reads deposited inputs and checks accepted membership | `run` exposes the original primary and dose-summary sequence; the new entry has not been validated by a repeat scientific fit. Later Linux conditional references have their own [retained sources](analysis/reference_sources/README.md) |
| GSE277064 | [`pipeline/g277_model.py check`](analysis/g277_execution/README.md) stages primary inputs; `report-saved` uses retained results | `run` executed the fixed five-stage primary sequence once on Linux, followed by the [complete retained-output comparison](verification/G277_20261029/README.md), with strict differences preserved. Later alternative-method orchestration is a separate scope |
| GSE205906 primary | [`pipeline/g205_model.py`](analysis/g205_execution/README.md) defaults to a plan using the prepared HISAT2 strand-0 object | `--execute` invokes the original checked consumer and serializer. Frozen-input planning and synthetic adapter checks passed; historical numerical equivalence through this new entry remains untested |
| GSE205906 matched STAR/HISAT2 | [`pipeline/g205_matched_model.py`](analysis/g205_matched_execution/README.md) defaults to an authenticated STAR s0 plan; [complete saved exports](data/g205_matched_exports/README.md) support result comparisons | `--prepare-only` performs deterministic preparation; `--execute` invokes the unchanged matched STAR kernel once. Synthetic/AST and six Python tests passed. New real-model execution and historical numerical equivalence remain unverified |
| GSE298414 extensions | Preparation/constituent guides in section 1 supply full results and no-fit readback | Original scientific sources retain historical paths; they document the calculations but are not portable raw-to-model launchers |

The observed-model commands are distinct from the later conditional-label reference runs. Reproducing one does not reproduce the other. The branch guides contain required arguments, exact input roles, dependencies and output contracts.

## 4. Restore the stated environment

The [GSE205906 Linux R-core guide](environment/G205_core/README.md) provides exact source identities, a serial installation order and `pipeline/restore_environment.py`. Its default is a read-only plan; `--execute` installs into a new isolated library on an existing compatible Linux/R installation. The new installer has synthetic control-flow tests; a fresh native installation was not performed for this release. Other branches have their own runtime requirements in the linked guides. Whole-study clean-machine restoration remains unverified.

## 5. Interpret the reproduction correctly

Keep native BH values, GSE128623 empirical BH descriptors and JGW03 conditional P values in their stated families and reference spaces. Report effects and missingness alongside thresholds. Preparation identity, component exposure, independent PSS validation and historical chronology retain the qualifications in the result guides. No independent outcome-protected PSS transcriptomic holdout has been established.

This public release covers the components listed above. The [PRJNA886048 extension](extensions/PRJNA886048/README.md) now contains the complete 12-library result and validated retained-input analysis on Windows and Linux. Its candidate status, per-column versus absolute-cell numerical comparison, raw-processing limits and unchanged inferential scope are documented in the extension. Historical version 0.30.0 remains unchanged.

Original code is MIT-licensed; project documentation and aggregate tables are CC BY 4.0. [Third-party notices](THIRD_PARTY_NOTICES.md) specify retained upstream terms. Cite the repository version and commit; a Zenodo DOI awaits finalized creator metadata.
