# Submitted-unit data ledger

This is the **15 October 2026 metadata snapshot**, released on 24 October. It records the analysis membership established at that date and preserves uncertainty about biological units. It is not a complete inventory of later execution events.

## Contents and exact transformations

`ANALYSIS_MEMBERSHIP_20261015.tsv` is byte-identical to the retained 165-row membership file. `DATASET_REPORTING_20261015_public.tsv` retains all nine reporting rows and every field except `current_large_data_location`, which contained private storage locations. No scientific value, sample ID, group label, exclusion or uncertainty qualifier is changed. `PROVENANCE.json` records the original archive/file hashes and the public file hashes.

`historical_builder.py` is the original metadata-only builder, preserved byte-for-byte. Its 15 pinned source metadata/report inputs are not distributed in this release. The builder documents the historical transformation; do not expect that original builder to run using this directory alone. The public `pipeline/ledger.py verify` command is self-contained and checks the released ledger's identity and membership invariants. Its `export` command requires the exact retained source archive and does not rerun the historical models.

## Interpret the counts

| Source | Submitted identifiers | Qualification |
|---|---:|---|
| GSE16710 | 27 | 26 primary; GSM418682 omitted only from primary and retained in all27 sensitivity |
| GSE128623 | 30 | all30 primary; GSM3681230 omitted only from predefined29 sensitivity; raw/adapter arms reuse libraries |
| GSE277064 | 63 | 45 natural-course and 18 trained; six untrained day-14 comparators are reused, giving 69 branch uses |
| GSE289929 | 12 | PC12, own-accession vehicle controls |
| GSE227494 | 12 | C2C12, own-accession vehicle controls |
| GSE299063 | 12 | HepG2, own-accession vehicle controls |
| GSE205906 | 9 | Same nine libraries across the two recorded historical processing branches |
| GSE85871 | Unknown selected membership | 212 source profiles is an inventory count; candidate-dependent analysis was not run |

The first seven rows total 165 distinct submitted identifiers. They are not 165 independently authenticated biological replicates. The ledger retains the source labels and specifically qualifies unresolved animal linkage, pooling, batch, culture replication and the GSE205906 time field. Legacy claims explicitly labelled `NOT_ADOPTED` are preserved for transparency, not endorsed.

## Later work and evidence references

Subsequent GSE128623 conditional-reference reporting, JGW03 Linux references, GSE205906 STAR/HISAT2 matched comparison and GSE277064 method sensitivities are not additional independent cohorts. The dated membership table does not enumerate every later computational branch. This release's `data/derived` files contain later accepted GSE128623/JGW03 aggregate reporting results and are separately identified; do not date them to 15 October.

`dated_evidence_reference`, `reporting_evidence` and related fields retain original project-relative evidence identifiers. They are provenance labels, not working hyperlinks into this public repository. The referenced private reports have not been silently replaced. The original source bundle hash in PROVENANCE can bind this snapshot to the retained evidence, but a hash alone does not authenticate pre-outcome timing.

The snapshot records 27 resolved treated-library RIN values and nine unresolved vehicle RIN mappings in the original sidecar's checks; numeric RIN values are not republished here. No low-RIN refit is implied. No downstream candidate gate is promoted from unevaluated to passed.
