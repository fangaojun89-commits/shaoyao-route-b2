# Technical inputs and access records

`INPUTS.json` publishes 71 selected file identities projected from retained project records:

| Records | Role |
|---|---|
| 19 | Historical GEO source files, platform annotation and series metadata |
| 30 | Fixed GSE128623 single-end FASTQ inputs |
| 18 | Two FASTQ mates for each of nine GSE205906 libraries |
| 4 | Compressed/uncompressed GSE205906 Rnor6/Ensembl101 reference assets |

Download URLs and recorded expected sizes/checksums identify files; they are not newly tested availability claims or permission grants. The catalog contains no raw sequencing reads. The original technical-manifest hashes are recorded in `source_records`. GSE205906 animal-independence, dose and grouping assertions were deliberately excluded from this technical projection; use the qualified study design and manuscript for interpretation.

Place legally obtained files under an input root according to `path`, then run:

```sh
python -I -B pipeline/acquisition.py inspect
python -I -B pipeline/acquisition.py check-inputs --root /absolute/path/to/inputs
```

Checking streams file bytes without extracting, downloading or modifying them. Missing files and mismatches produce a nonzero exit code. FASTQ and reference placement paths are portable catalog aliases, not replacements for historical absolute paths in original analysis scripts. These need explicit binding when reconstructing a run. MD5 identifies the historical FASTQ delivery; it is not a security signature. GEO supplementary bundles may contain more libraries than the selected analysis; the sample selection remains in the existing branch metadata.

The four GSE205906 reference assets are tied to the accepted October asset manifest. Version 0.7.0 recovers the two official Ensembl release-101 compressed-file URLs and the two gzip-decompression relationships. They must be **Rnor6 / Ensembl101**, not the earlier mRatBN7.2 inventory. For each uncompressed record, `url` remains `null` because it is produced by decompressing the `derived_from.path` asset, rather than downloading a separate URL. Check the compressed hash first, decompress to a new path without modifying the source, then check the uncompressed size/hash.

The GTF URL was rechecked on 24 October using a TLS-verified 11,791,561-byte retrieval in memory; both compressed and decompressed SHA-256 identities exactly matched the retained assets. The FASTA URL and both identities come from the accepted 3 October full-file/gzip-CRC receipt; the 831,597,371-byte FASTA was not downloaded again. These are deliberately different evidence types. Do not substitute a newer assembly, annotation or freshly changed resource solely because its filename looks similar.

The GSE128623 published fold-change table is explicitly `published_context_only`; it does not replace the complete FASTQ-derived count matrix. The older GSE205906 public supplementary table likewise is a recorded source, not proof of the nine-library counting result. Existing source code and accepted producer manifests govern derived objects.

This selection does not close every scientific input. Frozen gene sets, orthology and platform transformations are enumerated in `analysis/GSE16710/`; their original versions remain required. Derived matrices, prepared objects, worker aliases and execution manifests still need complete producer-to-consumer binding. No model was rerun for this catalog.
