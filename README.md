# Shaoyao Route B2: reproducibility components

Version 0.6.0 adds **71 technical input records and the exact 70-package GSE205906 Linux core environment** to the original scientific branch kernels, conditional-reference sources, experimental-spasticity anchor, accepted-table reporting companion and dated data ledger for a secondary transcriptomic study of Shaoyao–Gancao-related preparations and post-stroke spasticity. These are components of the study, not the complete raw-to-model pipeline. `.gitattributes` preserves pinned bytes on Windows checkouts; two explicitly identified historical-source path derivatives have separate original/public hashes. Earlier tags remain in history.

## What is available

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

The complete study source archive, raw/intermediate data, manuscript, third-party resources and software dependency archives are not included in v0.6.0. The source and acquisition inventories supply many previously absent identities, producers and readers, but complete acquisition, dynamic execution manifests, environment provisioning and portable orchestration remain unfinished. Four accepted reference assets still lack recovered original URLs in this selection. The 70-package core installation passed previously on an existing Linux runtime; this release checked the retained archives but did not execute a new Linux restoration. Whole-study cold restoration is not established. The ledger is explicitly a 15 October metadata snapshot, not a complete log of later executions. The anchor metadata preserves original curation language without independently establishing animal lineage or pre-outcome chronology.

## Licences and attribution

Original code: [MIT](LICENSE). Original explanatory documentation and the project's aggregate output tables: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), attributed to the Shaoyao Route B2 project contributors. This collective project credit is not a final scholarly author list. See [scope and resource notices](THIRD_PARTY_NOTICES.md).

Source accessions: [GSE128623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128623); JGW03 culture branches [GSE289929](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289929), [GSE227494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227494), and [GSE299063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299063). Their source data and the resources used upstream retain their own terms.

## Citation and archive

Cite this repository, version and commit for the reporting component. A Zenodo DOI is planned for a deposited version after creator metadata are finalized; no DOI is assigned or implied here. Final manuscript authors, affiliations and contact details are supplied separately by the project owner.
