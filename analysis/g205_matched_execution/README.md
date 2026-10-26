# Portable matched STAR downstream entry

This entry connects the accepted public STAR strand-0 count preparation to the unchanged `g205m_fit()` scientific kernel. It fits only the STAR branch of the same nine-library quantifier sensitivity. It does not rerun alignment, counting, HISAT2, label references, or the saved comparison, and does not replace historical results.

## Inputs and command

Reuse `prepared/G205_STAR_s0_inputs.rds` from the [input bridge](../input_bridge/README.md). Its exact accepted identity is in `CONTRACT.json`; HISAT2, STAR s1/s2 and edited preparations are rejected. Reuse the frozen resource RDS and external orthology TSV under a [resource workspace](../resource_build/README.md). Resource memberships and their source snapshots are not redistributed here. All original source and resource identities are checked before creating output.

Use the existing pinned R 4.6.1 library, including the configuration's versions and BiocParallel 1.46.0. No installation occurs. Supply a new output directory, preferably on a data disk:

```sh
python -I -B pipeline/g205_matched_model.py --prepared /inputs/prepared/G205_STAR_s0_inputs.rds --resource-root /resources --rscript /R/bin/Rscript --r-library /R/library --output /new-plan
```

Default mode validates file identities and records a plan; it starts no R process. Add `--prepare-only` with a different new output directory to inspect runtime and run the original deterministic annotation/design/input preparation without fitting or drawing random values. Add `--execute` only for an intentional new STAR fit. It invokes the original matched kernel once, retaining its fixed RNG, single-thread backend, own expression filter, TMM, method order, program eligibility, 22 families, native missingness and P/q rules. Failed output directories are retained; there is no automatic retry or resume that could duplicate a fit.

Successful execution writes the seven complete native tables, full source domains, slots, coverage and native-not-returned ledger, plus a native RDS. A separate candidate RDS is saved during the original kernel's return callback before its final integrity check; it is never labelled accepted. A later export failure leaves completed native computation available for investigation without automatic refitting. Receipts describe the new process, exact inputs, current runtime, source identities and output hashes. Raw condition text is not persisted.

## Scientific mapping and provenance

The input bridge already ran the original input loader and checked all counts, sample order and GTF fields. `g205p_prepare()` starts from this authenticated loader result. Its contiguous preparation block—from annotation projection through orthology selection, resource assignment, original design and count/annotation validation—is AST-identical to `g205m_prepare()`. Only provenance/IO binding differs. It creates a distinctly named public-input packet rather than invoking or weakening the historical `g205m_matrix_binding()` contract. No historical COMPLETE, producer QC, process receipt or absolute-path namespace is fabricated.

The original engine, fit, runtime, source-domain constructor, factor-label adapter and TSV transport/writer functions are reused unchanged. The historical source files remain byte-for-byte intact. EdgeR uncertainty remains explicitly unreported where the matched kernel has no validated native QL interval. Nine libraries are the observation units; this is not an independent animal-level replication.

## Validation and limits

`test_structure.R` checks original source hashes, the full preparation AST block, synthetic 32,883-by-9 preparation and axes, contrast directions, no RNG/input mutation, negative inputs, strict one-to-one orthology, whole-graph copying, gzip RDS decoding, factor labels and typed-empty export. `test_entry.py` tests default no-R planning and refusal of wrong inputs/resources/sources, existing outputs and conflicting modes. Use explicit scratch on a data disk:

```sh
Rscript --vanilla analysis/g205_matched_execution/test_structure.R /path/to/repository /existing/test-scratch
# Set G205_TEST_SCRATCH to an existing data-disk scratch directory first.
python -I -B analysis/g205_matched_execution/test_entry.py
```

These tests do not fit real models. No new STAR fit or cross-platform numerical equivalence was established when adding this entry. `--execute` produces `NEW_NATIVE_RESULT_SAVED_NOT_NUMERICALLY_ACCEPTED`, not historical acceptance. A new result requires separate keyed comparisons against the right-hand STAR fields of the [complete saved exports](../../data/g205_matched_exports/README.md), preserving native unavailable values and family membership. The old saved-comparison driver also needs its original native tables/checkpoints; the union export is not a substitute for those file formats. The matched fit itself does not need raw FASTQ, genome indexes or historical process receipts.
