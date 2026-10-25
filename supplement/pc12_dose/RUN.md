# Reproduce the retained PC12 dose-encoding comparison

From the repository root, using Python with NumPy and pandas installed:

```sh
python -I -B pipeline/pc12_dose_comparison.py --output-dir pc12-dose-comparison-output
```

Select a new output directory. The command checks the 16 public projected-input hashes and preserves existing output by refusing to overwrite it. It needs neither a network connection, a private project directory, R, nor the upstream gene-set/network files. Tested with Python 3.12.14, NumPy 2.3.5 and pandas 3.0.1. These are reporting-runtime versions, not the original R fitting environment.

The command computes descriptive comparisons from saved estimates: directions, effect ranks, discovery overlap, fitted endpoint differences and supported interval zero-inclusion changes. It reuses original P values and FDR values. No model, dispersion, score, label reference or candidate gate is rerun.

Outputs are the seven-family TSV summary, seven compressed paired-feature tables under `paired/`, and a source-hash/scope receipt. The accepted [summary](Route_B2_PC12_dose_encoding_comparison_20261025.tsv) and [interpretation](README.md) are included for comparison. Local verification reproduced the summary byte-for-byte and all seven paired tables byte-for-byte after decompression; [VALIDATION.json](VALIDATION.json) records the checks. Gzip container bytes may vary between executions because of compression timestamps, without any change in table content.

The gene comparison retains 13,325 pairs, including 13,227 matching directions; ordinal/log-dose discoveries are 40/42, with all 40 shared. High-dose/trend-coherent counts are also 40/42, with 40 shared. These findings describe sensitivity within the same 12 libraries. They do not establish equivalence, exposure feasibility, independent replication or a molecular–spasticity association.

The code is a portable derivative of the project's dated saved-fit summary. Numerical comparison logic is retained; paths, public-input verification and output handling were adapted. Input provenance is recorded [here](../../data/pc12_dose_saved/PROVENANCE.json). The report excludes unsupported historical approximate gene intervals and does not calculate a confidence interval for the between-model difference.
