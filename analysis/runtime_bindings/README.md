# GSE205906 runtime names and immutable source identities

`G205_BINDINGS.json` projects file identities from 12 retained records across five GSE205906 reference/baseline runs. It contains 114 file pins, including 18 occurrences that resolve to public source/configuration files. Repeated pins describe shared dependencies, not independent runs or observations. Four previously undisclosed source files are supplied byte-for-byte from the sealed source archive. No models were executed to produce this map.

## Why filenames alone are insufficient

| Run or stage | Runtime name | Source identified by size and SHA-256 |
|---|---|---|
| Original label reference | `worker.R` | `gse205906_label_reference_20261014.R` |
| Canonical repair | `worker.R` | `gse205906_canonical_label_reference_20261015.R` |
| Gene-method reference | `worker.R` | `gse205906_gene_method_label_reference_20261015.R` |
| Unsigned reference | `worker.R` | `gse205906_unsigned_method_label_reference_v2_20261015.R` |
| Unsigned observed metadata completion | `finish_observed_bridge.R` | `finish_unsigned_observed_bridge_20261015.R` |
| Canonical independent lineage check | `lineage_checker_v2.R` | `independent_gse205906_canonical_lineage_v2_20261015.R` |
| Original program-baseline input inventory | `worker.R` | `gse205906_program_baselines_20261014.R` |
| Completed program-baseline recovery | `worker_v2.R` | `gse205906_program_baselines_20261014_v2.R` |

The program-baseline distinction is material: its original input manifest names the earlier worker, while `RECOVERY_COMPLETE.json` identifies the actual v2 recovery worker. Both are retained; the old inventory must not be treated as proof that the earlier worker produced the completed recovery. Similarly, the unsigned preparation contains a separately pinned metadata finisher; its presence is not evidence that an observed fit was repeated.

Two additional canonical `summary.R`/`checker.R` bindings come from literal source assertions and are labelled separately from runtime receipt observations. A source assertion defines a requirement; it is not itself independent evidence that a run executed.

## Use

```sh
python -I -B pipeline/bindings.py verify
python -I -B pipeline/bindings.py stage-sources --output /absolute/path/to/new-source-layout
```

The second command copies only the hash-matched public source/configuration files to their historical relative namespaces under a new directory. It preserves exact bytes, verifies staged identities and never executes R. Existing destinations are refused. Stage receipts identify whether each mapping was observed in a retained runtime record or required by a source assertion.

This is a source layout, not a runnable reconstructed analysis. External assets such as count matrices, prepared RDS objects, gene resources, observed checkpoints and allocation manifests are represented by identities but are not copied or generated. Source files retain their historical path assumptions and identity guards. A portable execution adapter must explicitly bind these remaining assets and cannot treat a different-byte relocated manifest as the historical original. In particular, the canonical repair depends on the accepted original label reference; do not start a new original enumeration merely to exercise this staging command.

The projection retains file roles and relative cloud namespaces, not host credentials, personal paths, raw error text or scientific arrays. Original record byte counts and hashes are included under `receipt_provenance`; the originals remain privately retained and were checked when making this projection. Neither this metadata map nor source staging authenticates animal exchangeability, proves whole-study cold restoration, or supplies an independent biological holdout.
