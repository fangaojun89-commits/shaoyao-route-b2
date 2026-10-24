# Prepare accepted counts with the original input readers

This entry converts the published GSE128623 and GSE205906 count inputs into six portable, unfiltered R input objects. It calls the unchanged historical input readers, preserving gene identifiers, annotation, count values, sample order and the predefined GSE128623 29-library subset. It supplies an executable next step after downloading the repository. Scientific fitting, normalization and new random allocation are not part of this command.

## Requirements and command

Use Python 3.12 or later, an existing R installation and an existing library containing `data.table` 1.18.4, `openssl` 2.4.2, `jsonlite` 2.0.0 and `digest` 0.6.39. The tested runtime and content checks are recorded in `VALIDATION.json`. These four packages serve input conversion; they are not a complete model-analysis environment.

Obtain the frozen [Ensembl 101 rat GTF](https://ftp.ensembl.org/pub/release-101/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.101.gtf.gz), or reuse an existing copy. `CONTRACT.json` specifies both compressed and uncompressed file identities and all six original R source identities. The command validates the full files before accepting them.

From the repository root:

```sh
python -I -B pipeline/prepare_count_inputs.py \
  --gtf-gz /path/to/Rattus_norvegicus.Rnor_6.0.101.gtf.gz \
  --rscript /path/to/Rscript \
  --r-library /path/to/existing/R-library \
  --output /path/to/new-input-workspace
```

Use ordinary command continuation syntax for your shell, or put the command on one line. Omit `--r-library` if the required packages are already on R's default library path. No installation or download occurs. The output must be a new directory. Allow approximately 100 MB of output space, preferably on a data disk. A failed workspace is retained for diagnosis and must not be treated as complete.

## Outputs and checks

| Output | Content and use |
|---|---|
| `original_inputs/` | Eleven restored files with the exact accepted bytes and historical relative names. |
| `derived/G205_gene_attributes.tsv` | All 32,883 gene records and all 14 fields extracted from the frozen GTF, including missing names, raw attributes and source ordinals. |
| `prepared/G128_raw_inputs.rds`, `G128_adapter_inputs.rds` | Original-reader counts, annotation, symbols and ledger; the 30-library design and contrast matrix; predefined 29-library column indices. |
| `prepared/G205_HISAT2_s0_inputs.rds`, `G205_STAR_s0_inputs.rds` | Original-loader inputs for the original and matched-primary quantification specifications, with the original transformation record. |
| `prepared/G205_STAR_s1_inputs.rds`, `G205_STAR_s2_inputs.rds` | Strand-diagnostic inputs, kept distinct from the matched-primary specification. |
| `receipts/` | The command, R logs/session information, content QA and final `PREPARE_COMPLETE.json`. |

The G128 analyzer and its independent validator must agree on the loaded counts, annotation, symbols and ledger. Their design values and axes must also agree; their local factor-variable names are not scientific differences. Both supplied 29-library matrices must be exact column subsets excluding GSM3681230. The same authenticated ledger is used for the reader's frozen-design comparison; this does not invoke or satisfy the historical execution driver's separate custody checks.

For G205, a streaming extraction limits memory use. Every extracted gene field is then compared with the original R GTF parser, applied in bounded chunks with source ordinals offset to their full-file positions. The original input loader checks the complete annotation join and renames the nine original matrix labels to their ledger GSM identifiers without changing order or counts. Fixed gene identifiers come from the authenticated complete HISAT2 matrix, not a fabricated historical execution manifest. All six saved objects are read back and compared with their in-memory inputs. The entry checks that R's random state is unchanged.

## Next analysis stage

These objects are **new input preparations**, not historical fitted checkpoints. They do not replace accepted model results. For inspection, use `x <- readRDS(...)`; G128 input components are in `x$input`, with the predefined subset in `x$prespecified_29_columns`; G205 components are in `x$inputs` and transformation details in `x$transformation`.

Future reconstruction must additionally supply the frozen program resources, the recorded branch-specific parameters, the appropriate full environment and a separately validated analysis adapter. The unchanged scientific kernels are described in [study_sources](../study_sources/README.md); historical source aliases and remaining pinned objects are described in [runtime_bindings](../runtime_bindings/README.md). The old G128 CLI is not a portable launcher, and the September G205 CLI/configuration must not be substituted for the accepted October analysis.

This preparation closes the count-to-original-reader step. It does not establish whole-study model reproduction, independent animal/culture lineage, or an independent PSS phenotype holdout. Technical quantification specifications retain the same source libraries and must not be counted as independent biological replications.
