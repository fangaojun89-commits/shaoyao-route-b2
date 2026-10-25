# Shaoyao Route B2: reproducibility components

Version 0.27.1 clarifies that the complete GSE205906 HISAT2 and STAR count matrices are already public, and gives their exact restoration command and destinations in the [matched-source guide](analysis/g205_matched_sources/README.md). Scientific inputs, source code and results are unchanged.

Version 0.27.0 connects [the completed GSE205906 matched-alignment sources](analysis/g205_matched_sources/README.md): 20 unchanged source roles, 19 reused public copies and one newly supplied metadata-preparation source. Eleven identities are matched to the retained scientific completion record. This source-navigation release changes no scientific result and does not provide a portable cold-restore adapter.

Version 0.26.0 adds [complete constituent-comparison results and readback](supplement/constituent_extension/RUN.md): 126 libraries, 60 new gene contrasts, 240 program families and planned influence results. The portable check reconstructs retained statistics and library-score effects without refitting. Earlier preparation and canonical results are preserved.

Version 0.25.0 adds [GSE298414 preparation results and public readback](supplement/preparation_extension/RUN.md): 72 libraries, 36 gene contrasts, four program collections, source metadata and library-influence results. Unchanged original analysis sources are supplied; the portable command verifies retained statistics without refitting. This release does not claim cold restoration or new PSS validation.

Version 0.24.0 adds the [current six-preparation identity and exposure interpretation](data/preparation_identity/README.md). It distinguishes nominal GSE205906 group validity from unresolved actual exposure and qualifies JGW03 chemistry, while preserving historical source records and all accepted analyses.

Version 0.23.0 adds a [reproducible seven-family PC12 dose-encoding comparison](supplement/pc12_dose/RUN.md), using 16 projected retained-result tables. The summary and full paired tables reproduce the accepted article v36 comparison without refitting any model.

Version 0.22.0 adds a [portable seven-table reading-package command](supplement/late_results/RUN.md). It reconstructs the accepted displays and standalone table bundle from public retained results, using only the Python standard library. No scientific model is rerun.

The [completed late-results supplement](supplement/late_results/README.md) now exposes all seven accepted reporting tables (669 rows across overlapping families/views), including the complete 22-family GSE205906 matched-alignment comparison. Five tables are newly supplied and two reuse existing public files. No scientific analysis was rerun.

Version 0.20.0 connects the unchanged JGW03 primary and dose-summary scripts to their original deposited-count inputs, selected source metadata, fixed resources and S1 stopping result. The no-fit check reproduces all 1,636,380 count values and their order across 36 libraries and passes the original 51 metadata checks, retaining three source limitations. [Commands and executed scope](analysis/jgw_execution/README.md). The explicit production branch is supplied for future reconstruction and was not executed; completed scientific results are unchanged.

Version 0.19.0 publishes the retained protocol/SAP and a [22-item retrospective original-plan/current-method crosswalk](protocol/README.md). It discloses major amendments, unfulfilled broad-plan items, candidate-gate stopping and chronology limits. This is a documentation release: no completed analysis or scientific payload is changed.

Version 0.18.0 adds an explicit five-stage GSE277064 primary execution entry and a tested saved-result reporting mode. The latter reuses accepted outputs without fitting and reproduces all 24 approved gene-table projections, omitting only the three superseded uncertainty columns. Previously missing summary, validator and historical template dependencies are now supplied with distinct provenance. [Commands, current reporting rule and tested scope](analysis/g277_execution/README.md). Full production execution is available for an intentional future reproduction; it was not run for this release.

Version 0.17.0 connects GSE277064 retained preprocessing to original annotation/resource/input consumers. Five saved input views preserve counts, masks, normalization metadata and existing voom inputs, without fitted results, scores or network memberships. The command reuses v0.16.0 preparation and the verified resource build. No completed fit, score, filtering or normalization is repeated. [Command and scope](analysis/g277_saved_inputs/README.md).

Version 0.16.0 adds the previously missing GSE277064 source-input and original prefit-design entry: three unchanged GEO count files, an exact 63-title metadata selection, and the original 45-natural/18-trained library definitions with six reused D14 comparators. The command stops before normalization, scoring and fitting. See the [G277 preparation guide](analysis/g277_preparation/README.md). Earlier release descriptions retain their versioned scope; accepted analyses are not rerun.

Version 0.15.0 makes the completed experimental-spasticity anchor results reusable through the unchanged original 51-check validator and a separate readback of all 1,771 gate decisions and 46,046 saved omission estimates. Its 30 input roles match the retained internal output manifest; two explicit fgsea projections omit only derived leading-edge membership lists. No model, activity score or random allocation is recomputed. See the [saved-anchor command and scope](analysis/anchor_saved/README.md). The current gate remains zero admissions; whole-study raw-to-model and cold-environment restoration remain unverified.

Version 0.14.0 connects all three frozen JGW03 input tables, rebuilt resources and saved allocation registries to the unchanged original input reader. It supplies three portable unfitted bundles and checks their direct count/annotation/metadata readback without recomputing scores, fitting models or drawing labels. See the [JGW input preparation command](analysis/jgw_preparation/README.md). Earlier version descriptions below retain their original scope; whole-study cold restoration remains unverified.

Version 0.13.0 adds original GSE205906 deterministic preparation and consumer inspection without fitting: ten design/filter/TMM phases and their program memberships are connected to the original core and both reference preparers. Version 0.12.0 added an explicit GSE128623 model entry with separate no-fit input checks, full synthetic-kernel tests and a Linux-only real-model command. Accepted real-data analyses are not rerun for this release. Version 0.11.0 added a tested reconstruction of the frozen resource bundle from eight local source files, with all thirteen outputs matching accepted bytes and original consumer checks passing. The original UTF-8 locale is explicit. Version 0.10.0 supplies a portable input-preparation entry that calls the unchanged original GSE128623/GSE205906 readers, verifies the complete annotation conversion, and saves six unfitted input objects. Version 0.9.0 supplied their accepted complete count matrices and associated metadata. The repository also provides frozen JGW03 inputs for 36 libraries, selected scientific analysis sources, conditional-reference sources, the experimental-spasticity anchor, input and environment inventories, reporting code, and a dated data ledger for this secondary transcriptomic study. These components support inspection and reuse; the complete raw-to-model pipeline has not yet been restored. `.gitattributes` preserves pinned bytes on Windows checkouts. Two explicitly identified historical-source path derivatives have separate original/public hashes. Earlier tags remain in history.

## What is available

- `pipeline/g277_model.py`: separate no-fit staging, explicit five-stage primary execution, and accepted-result reporting modes. [Guide](analysis/g277_execution/README.md).

- `pipeline/connect_g277_saved_inputs.py`: original annotation/resource consumers for saved G277 preprocessing; [input views and scope](analysis/g277_saved_inputs/README.md).

- `pipeline/prepare_g277_inputs.py`: original G277 input admission and prefit designs, with unchanged [source-derived counts](data/g277_inputs/README.md) and explicit upstream terms; [command and output](analysis/g277_preparation/README.md).

- `pipeline/anchor_saved.py`: execute the original retained-output validator and separately check the complete saved gate/omission arithmetic; [inputs, environment, command and interpretation](analysis/anchor_saved/README.md).
- `data/anchor_saved_outputs/`: source-pinned retained anchor tables for this readback, with two explicitly documented column projections; no source networks or gene-set memberships.

- `pipeline/prepare_jgw_inputs.py`: original three-context JGW input admission with the preserved kernel definitions bound; no fixed-layer construction, scoring or fitting. [Requirements and outputs](analysis/jgw_preparation/README.md).
- `pipeline/jgw_model.py`: original deposited-count/metadata check and explicit primary/dose-summary sequence; [commands, runtime evidence and distinction from later conditional references](analysis/jgw_execution/README.md).

- `pipeline/prepare_g205_model.py`: original G205 preparation from accepted HISAT2 strand-0 inputs, with no model execution; [command, provenance and limits](analysis/g205_preparation/README.md).

- `pipeline/g128_model.py`: original G128 input/runtime check, full synthetic analysis and explicit Linux execution entry; [commands, tested scope and remaining requirements](analysis/g128_execution/README.md).

- `pipeline/build_resources.py`: rebuild the thirteen frozen resource outputs with unchanged original code; [source inputs, package/locale requirements and command](analysis/resource_build/README.md).

- `pipeline/prepare_count_inputs.py`: prepare and verify six unfitted G128/G205 input objects using the original readers; [requirements, command and next analysis stage](analysis/input_bridge/README.md).

- `data/count_inputs/`: eleven accepted G128/G205 count and metadata files in reversible gzip form; [content and placement guide](data/count_inputs/README.md). `pipeline/count_inputs.py` verifies and restores their exact original bytes without model execution.
- `pipeline/report.py`: original Python reporting code, standard library only.
- `data/derived/`: accepted GSE128623 (448 families) and JGW03 (123 families) aggregate reference tables.
- `expected/`: three fixed descriptive outputs for comparison.
- `MANIFEST.json`: SHA-256 identities of the release files.
- `data/ledger/`: the 15 October submitted-unit membership and reporting snapshots, their provenance, and the original metadata builder source.
- `pipeline/ledger.py`: export and verification of that dated metadata selection. See the [ledger guide](data/ledger/README.md) for later analyses and interpretation limits.
- `analysis/GSE16710/`: seven unchanged original R sources, frozen sample/QC metadata, resource/input inventories, dependency lock and the retained S1 gate. See the [anchor reconstruction guide](analysis/GSE16710/README.md) for dependency order, external inputs and the later companion sources still outside this release.
- `pipeline/anchor.py`: read-only source/metadata and enumerated external-input verification; never executes R or reruns a model.
- `analysis/study_sources/`: 38 unchanged R sources, five original configurations, 99 submitted-unit metadata rows and explicit GSE205906 source-role identities. See the [branch guide](analysis/study_sources/README.md) for dependency order, historical entry points and remaining reconstruction requirements.
- `pipeline/sources.py`: read-only release/source identity and embedded source-role verification.
- `analysis/reference_sources/`: a 151-source inventory spanning conditional references, anchor companions, matched quantification and saved-result checks; 124 source files are newly provided and 27 reuse prior public locations. [Guide and limits](analysis/reference_sources/README.md).
- `pipeline/references.py`: read-only verification of that inventory and the literal caller bindings of three source-archive omissions now supplied.
- `data/acquisition/`: 19 historical GEO files, 30 GSE128623 FASTQ identities, 18 GSE205906 FASTQ identities and four matching Rnor6/Ensembl101 reference assets. [Input placement and remaining access gaps](data/acquisition/README.md).
- `environment/G205_core/`: exact package acquisition identities, accepted installed versions/runtime inventory, serial build order, and a path-parameterized library checker. [Environment restoration sequence and validation scope](environment/G205_core/README.md).
- `pipeline/acquisition.py`: inspect these catalogs, check user-supplied input/package bytes, or emit a source-verified installation plan. No automatic download, installation or scientific computation.
- `analysis/runtime_bindings/`: 114 file pins from 12 retained GSE205906 records, including 18 occurrences resolved to public source/configuration files; four additional original R sources and separately labelled source assertions. [Runtime names and interpretation](analysis/runtime_bindings/README.md).
- `pipeline/bindings.py`: verify these source bindings or copy exact public source/configuration files into a fresh historical namespace layout without executing them.
- `data/jgw03_inputs/`: 12 unchanged accepted inputs across C2C12, PC12 and HepG2, including fractional RSEM expected counts and the original fixed allocation registries. [Content, placement and source terms](data/jgw03_inputs/README.md).
- `pipeline/jgw_inputs.py` and `data/jgw03_inputs/verify_labels.R`: file/content/order checks and read-only checks of the saved label assignments; no fitting or new random draws.
- `data/derived/SAVED_OUTPUT_LINEAGE.json`: provenance for the already public G128/JGW tables, including the rechecked 73-file G128 small-evidence chain and its retained interpretation limits.

The derived tables are family-level aggregate outputs. The ledger adds public submitted-array/library identifiers and qualified analysis membership. Neither contains individual-level clinical records, raw sequencing data or third-party gene-set memberships.

## Run

From the repository directory, using Python 3.12 or later:

```sh
python -I -B pipeline/report.py --g128 data/derived/G128_FAMILY_SUMMARY.tsv --jgw data/derived/JGW_FAMILY_SUMMARY.tsv --output output
```

Use a new output directory. Compare the three named outputs with `expected/`. `RECEIPT.json` records the runtime and file identities. The reporting code and data are unchanged from the privately verified companion, which ran under Python 3.12.14 and reproduced the expected outputs byte-for-byte. This release does not claim testing on all later Python versions.

Verify the dated ledger independently of the reporting command:

```sh
python -I -B pipeline/ledger.py verify data/ledger
```

The verifier checks 165 submitted identifiers, 18 overlapping analysis specifications and nine reporting contracts across eight accessions. The identifiers do not establish 165 independent animals or cultures. GSE85871 selected membership remains unknown, not zero.

## Interpretation

The study's specified experimental-spasticity gate admitted 0 of 1,771 programs. Downstream candidate tests remained unevaluated. The published tables describe separate source-wide conditional-reference analyses; they are not independent biological or clinical validation. G128 empirical-q descriptors and JGW conditional-P descriptors retain different definitions and must not be pooled. Feature-family counts are not independent samples or necessarily unique genes.

The command computes descriptive groups from already accepted tables. It does not refit models, regenerate random allocations, recalculate P values or demonstrate an antispasticity mechanism. The study was not publicly preregistered. No independent outcome-protected PSS transcriptomic holdout has been established.

## Scope of this release

The complete study source archive, most raw/intermediate data, manuscript, third-party gene-set resources and software dependency archives are not included in this release. Selected JGW processed inputs/allocation registries, accepted G128/G205 count matrices, unchanged G277 source counts, saved G277 preprocessing views and the saved anchor validation tables are included. G128/G205 input preparation and frozen-resource construction are executable and content-verified. The G128 observed-model sequence now has an explicit entry: its statistical kernels are tested on generated data, while real raw/adapter inputs pass no-fit checks. The portable Linux real-model driver has not been executed, and other branches still need equivalent execution entries. G277 now has an explicit primary sequence and tested accepted-result postprocessing; its full new model run and later alternative-method orchestration are not validated by this release. JGW now has an explicit original primary/dose-summary entry, with source-deposited counts checked against accepted inputs; its scientific run and subsequent conditional-reference orchestration are not verified by that check. G205 now has a tested original preparation/consumer-inspection command; it stops before statistical fits and does not recreate historical native custody. The saved-anchor command checks completed results, not raw-to-model reconstruction. Complete acquisition, dynamic execution manifests, environment provisioning and portable orchestration remain unfinished. The four GSE205906 reference assets have recovered official acquisition or gzip-derivation paths; the exact GTF bytes were revalidated in memory, while the FASTA is supported by the retained full-file receipt. The 70-package core installation passed previously on an existing Linux runtime; no new Linux restoration was executed for this release. Whole-study cold restoration is not established. The ledger is explicitly a 15 October metadata snapshot, not a complete log of later executions. Historical metadata wording does not independently establish animal/culture lineage or pre-outcome chronology.

## Licences and attribution

Original code: [MIT](LICENSE). Original explanatory documentation and the project's aggregate output tables: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), attributed to the Shaoyao Route B2 project contributors. This collective project credit is not a final scholarly author list. See [scope and resource notices](THIRD_PARTY_NOTICES.md).

Source accessions: [GSE128623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128623); JGW03 culture branches [GSE289929](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289929), [GSE227494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227494), and [GSE299063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299063). Their source data and the resources used upstream retain their own terms.

## Citation and archive

Cite this repository, version and commit for the reporting component. A Zenodo DOI is planned for a deposited version after creator metadata are finalized; no DOI is assigned or implied here. Final manuscript authors, affiliations and contact details are supplied separately by the project owner.
