# Shaoyao Route B2: reproducibility components

Version 0.10.0 adds a tested portable input-preparation entry that calls the unchanged original GSE128623/GSE205906 readers, verifies the complete annotation conversion, and saves six unfitted input objects. Version 0.9.0 supplied their accepted complete count matrices and associated metadata. The repository also provides frozen JGW03 inputs for 36 libraries, selected scientific analysis sources, conditional-reference sources, the experimental-spasticity anchor, input and environment inventories, reporting code, and a dated data ledger for this secondary transcriptomic study. These components support inspection and reuse; the complete raw-to-model pipeline has not yet been restored. `.gitattributes` preserves pinned bytes on Windows checkouts. Two explicitly identified historical-source path derivatives have separate original/public hashes. Earlier tags remain in history.

## What is available

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

The complete study source archive, most raw/intermediate data, manuscript, third-party gene-set resources and software dependency archives are not included in v0.10.0. Selected JGW processed inputs/allocation registries and accepted G128/G205 count matrices are included. G128/G205 input preparation is now executable and content-verified; scientific model execution remains a separate reconstruction requirement. Complete acquisition, dynamic execution manifests, environment provisioning and portable orchestration remain unfinished. The four GSE205906 reference assets have recovered official acquisition or gzip-derivation paths; the exact GTF bytes were revalidated in memory, while the FASTA is supported by the retained full-file receipt. The 70-package core installation passed previously on an existing Linux runtime; no new Linux restoration was executed for this release. Whole-study cold restoration is not established. The ledger is explicitly a 15 October metadata snapshot, not a complete log of later executions. Historical metadata wording does not independently establish animal/culture lineage or pre-outcome chronology.

## Licences and attribution

Original code: [MIT](LICENSE). Original explanatory documentation and the project's aggregate output tables: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), attributed to the Shaoyao Route B2 project contributors. This collective project credit is not a final scholarly author list. See [scope and resource notices](THIRD_PARTY_NOTICES.md).

Source accessions: [GSE128623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128623); JGW03 culture branches [GSE289929](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289929), [GSE227494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227494), and [GSE299063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299063). Their source data and the resources used upstream retain their own terms.

## Citation and archive

Cite this repository, version and commit for the reporting component. A Zenodo DOI is planned for a deposited version after creator metadata are finalized; no DOI is assigned or implied here. Final manuscript authors, affiliations and contact details are supplied separately by the project owner.
