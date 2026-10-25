# GSE277064 original input and prefit design preparation

This command connects the three unchanged count tables and fixed sample/title records to the original source input block. It also reuses only the original domain-design prefix preceding `DGEList`. The original file is not changed or sourced as a complete analysis. No filtering, normalization, annotation mapping, activity scoring, inferential model, new allocation or accepted-result write is performed.

Use R 4.6.1 and an explicit existing library with data.table, R.utils, limma, jsonlite and digest. The tested Windows package versions and session information are recorded in the execution receipts; this is not a new Linux environment restore. Python requires only its standard library (3.12 or later; the execution-tested version is recorded separately). Use a fresh output directory:

```sh
python -I -B pipeline/prepare_g277_inputs.py --rscript /path/to/Rscript --r-library /path/to/library --output /new/g277-prepared
```

The script checks pinned source/input identities, stages copies below `project/`, and writes one `prepared/GSE277064_unfitted_inputs.rds` plus five sample-membership TSVs. The RDS contains unchanged 63-column counts, the title-checked ledger and five prefit designs/contrast definitions. Natural-course designs have five groups and six contrasts per cell type. Trained designs have four groups and three contrasts for V2a and astrocytes, reusing the original D14 comparators. Counts are not duplicated into each domain object. The original count assertions are additionally exercised on small invented duplicate/negative/nonfinite inputs. Base R independently reads every supplied numeric count for exact comparison.

The saved preparation contains no fitted objects, gene symbols, resource memberships, TMM factors or inferred effects. The original annotation/resource construction and full observed/alternative-model execution remain separate stages; this command does not claim whole-study raw-to-result restoration. Input group checks do not authenticate animal independence or establish a spasticity phenotype.
