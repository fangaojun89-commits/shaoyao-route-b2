# Preparation responses relative to each single herb

This is the manuscript preparation display, with all five descriptive intersections among the 420 primary preparation-program pairs. Each pair has concordant effect directions and original class-adjusted q<=0.05 against both peony and licorice. All 420 pairs, original estimates, nominal 95% intervals and comparator-specific adjustment values are retained in the supplied source table.

The source is a deterministic pairing of the previously public GSE298414 constituent outputs. It reuses the same 126 libraries; this is not an independent cohort. The five illustrated pairs occur at the high ordinal level in 70% ethanol. JGE02 is a 2:1 co-extract, JGE07 a separately extracted 1:1 mixture and JGE08 a separately extracted 1:2 mixture. The co-extract pairs depend on the influential library; the display and manuscript retain the scope of the original omission analysis. The intersection is descriptive, selected after results, and is not an interaction or synergy test.

## Recreate the display

Use Python with the pandas/numpy/matplotlib versions recorded in ENVIRONMENT.json:

```sh
python -I -B supplement/preparation_display/plot.py --output new-preparation-display
```

The command refuses an existing output directory. It reads all 420 pairs, checks the complete selected set, verifies finite ordered intervals, and regenerates the five-pair plotting table, PDF and PNG. No count input, model fitting, pathway scoring or multiple-testing adjustment is performed. The plotting and selection body is unchanged from the executed local producer; only source/output paths are parameterized.

Compare plotted_data_all_five_pairs.tsv with expected_plotted_data.tsv by rows and values (or after normalizing CRLF/LF). The distributed source and expected TSV use LF line endings; all original cell values and row order are preserved. PDF metadata may differ by render time; figure content and the original numeric inputs, not PDF-byte equality, are the scientific reuse target. The plot's source identity and output counts are recorded in SOURCE_AND_DISPLAY.json after execution. Original project code uses MIT; original aggregate tables/documentation use CC BY 4.0 with external rights unchanged.
