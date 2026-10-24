# Accepted GSE128623 and GSE205906 count inputs

These are complete retained gene-fragment count matrices used by the accepted analyses. The files were recovered from the existing computation, verified against its earlier receipts and checked for integer, finite, nonnegative values. No alignment, counting, model fit or label generation was repeated.

- GSE128623: original and adapter-processed 30-library matrices, plus each prespecified 29-library sensitivity matrix. The latter are exact column subsets omitting `GSM3681230`. The 30-library sample ledger and existing feature annotation are included. Original/adapter gene and sample axes agree.
- GSE205906: the original nine-library HISAT2 strand-0 matrix and all three STAR strand matrices. Strand 0 is the accepted primary specification. All matrices have 32,883 gene rows and the same nine-library order. The historical source/preparation overlay links labels to public accessions.

Files are gzip compressed with a zero timestamp. Decompression recovers the exact accepted original bytes; both compressed and original SHA-256 values appear in `PROVENANCE.json`. The retained runtime-relative paths identify their original input roles. The annotation can be much larger when decompressed because it preserves all original feature-coordinate fields.

```sh
python -I -B pipeline/count_inputs.py verify
python -I -B pipeline/count_inputs.py materialize --output recovered-inputs
```

Materialization requires a new directory, restores the eleven files under their retained relative paths and does not execute scientific code. It does not create a complete runnable study environment. Gene sets, additional prepared objects and full orchestration are separate requirements.

The original metadata is preserved, including dated preparation/QC flags and nominal biological-unit fields. Those fields do not establish current QC status, animal independence or verified administered dose. The 29-library subset and STAR strandedness alternatives retain their original roles; they are not additional independent cohorts. The existing paper and accepted results govern interpretation.

These source-derived public-data products retain applicable third-party rights; project curation is CC BY 4.0. See the [NCBI GEO data disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). The original project checker is MIT licensed.
