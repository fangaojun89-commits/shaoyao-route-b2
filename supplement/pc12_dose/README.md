# PC12 dose encoding robustness comparison

25 October 2026. **Direction and effect ranks were closely concordant across the two retained dose encodings, while some discovery sets changed.** All 40 ordinal high-dose/trend-coherent genes were retained under log-dose encoding, which yielded 42 coherent genes. The ordinal model remains primary. This is a descriptive comparison of existing fits, not a new model run, independent validation or candidate-gate test.

## Design and scope

SAP section 10 specifies ordinal versus log-dose encoding for PC12. The canonical GSE289929 dataset contains 12 libraries, three at each of 0, 4, 20 and 100 µg/mL. The primary predictor is 0, 1, 2, 3; the alternative is log10(dose + 1). Both gene models retain the same 13,325 filtered genes, TMM normalization and robust edgeR quasi-likelihood framework. Their dispersions were estimated separately in the original execution. The producing source is [`analyze_jgw03_three_cell_lines.R`](../../analysis/study_sources/code/R/analyze_jgw03_three_cell_lines.R).

This post hoc descriptive summary covers every matched ordinal/log-dose export in seven retained families. Feature sets matched exactly; all compared numeric entries were finite. Methods without a retained matching log-dose fit are not represented as completed comparisons. The summary metrics were selected for this reporting step; no prospective numerical equivalence margin was defined.

## Paired findings

| Family | Paired features | Same direction | Effect Spearman rho | FDR <0.05 ordinal / log-dose / shared | Discovery Jaccard |
|---|---:|---:|---:|---:|---:|
| edgeR_genes | 13,325 | 13,227/13,325 | 0.999561 | 40 / 42 / 40 | 0.952381 |
| Reactome_score | 1,455 | 1,447/1,455 | 0.999736 | 0 / 0 / 0 | NA (empty union) |
| QuickGO_SynGO_score | 68 | 68/68 | 0.998702 | 0 / 0 / 0 | NA (empty union) |
| Reactome_CAMERA | 1,455 | 1,444/1,455 | NA | 37 / 40 / 37 | 0.925000 |
| QuickGO_SynGO_CAMERA | 68 | 68/68 | NA | 23 / 23 / 23 | 1.000000 |
| PROGENy | 14 | 14/14 | 0.982418 | 1 / 1 / 1 | 1.000000 |
| DoRothEA | 104 | 104/104 | 0.999157 | 7 / 9 / 7 | 0.777778 |

No feature significant under either encoding had an opposite effect direction between encodings. Across all genes, 98 changed sign; the full-set agreement was 99.265%. For genes the Kendall tau-b was 0.984574, and the top 20 absolute-effect genes overlapped 20/20. Corresponding top-20 overlaps were 18/20 for Reactome scores, 19/20 for QuickGO/SynGO scores and 20/20 for DoRothEA. PROGENy has only 14 features, so the 14/14 full-universe overlap is not evidence of top-ranked stability; its effect-rank correlation is informative instead.

Reactome CAMERA retained all 37 ordinal discoveries and added three; DoRothEA retained all seven and added two. CAMERA is a competitive gene-set test and does not provide the same estimand as the mean-z-score program model. The absence of FDR discoveries in the two score families therefore does not contradict the CAMERA counts. Counts from these dependent families are not pooled as independent evidence. An empty discovery union has undefined Jaccard, not perfect agreement.

## Effect scale and uncertainty

Per-step ordinal coefficients and per-unit log-dose coefficients have different units. The paired files multiply each coefficient by the observed predictor range: ordinal ×3 and log-dose ×2.0043213738. Supported program-score/activity intervals are rescaled in the same way. These quantities describe fitted endpoint differences over the observed dose grid; they are not the categorical high-dose-minus-vehicle contrast. Positive rescaling leaves direction and effect ranks unchanged.

For genes, the median absolute difference between the two fitted endpoint log2 effects was 0.002131, with 95th percentile 0.010655. Historical approximate edgeR QL intervals reconstructed from logFC and F remain excluded, consistent with the manuscript's existing interval review; no substitute gene interval is invented. The full TSV supplies family-specific scale differences, supported score/activity interval zero-inclusion changes and rank metrics. These nominal interval summaries are not multiplicity-adjusted and should not be equated with FDR membership. It is inappropriate to pool effect magnitudes across gene, program-score and activity scales. No test or confidence interval for the between-model difference is constructed because the two estimates reuse the same samples and their covariance is not available here.

## Interpretation and reporting

| Assumption | Base case | Alternative | Result change | Conclusion stability | Risk level | Recommended language |
|---|---|---|---|---|---|---|
| Dose-response functional form on the observed grid | Ordinal | log10(dose + 1) | Closely concordant directions/ranks; genes 40 to 42, Reactome CAMERA 37 to 40, DoRothEA 7 to 9; all ordinal discoveries retained | Conditionally stable | Moderate for thresholded discovery claims | Effect directions and ranks were closely concordant across the two encodings, with modest changes in thresholded discovery sets. |

The classification is descriptive rather than a predefined equivalence verdict. It applies only to the retained PC12 encodings and source-wide outputs. It does not establish a continuous tipping point, a concentration-response mechanism, clinical exposure feasibility, or an independent molecular–spasticity association. It does not promote the two additional genes into the frozen primary set. J2 and other candidate gates retain their existing zero-eligible-input status. No downgrade of the manuscript's bounded source-wide conclusion is required; any claim of identical discovery sets would be incorrect.

No additional fits were undertaken. A formal model-difference uncertainty analysis, nonlinear dose-response model or encoding-by-method factorial would require additional estimation and is outside this saved-result comparison. These optional analyses are not prerequisites to reporting the present findings. The planned dose-encoding comparison has been quantified; missing preparation/exposure and independent PSS validation evidence remain separate limitations.

## Reproduction and files

- [Rebuild instructions](RUN.md) use retained inputs and do not fit a scientific model.
- [Accepted summary](Route_B2_PC12_dose_encoding_comparison_20261025.tsv).
- [Input provenance and column projections](../../data/pc12_dose_saved/PROVENANCE.json).
- [Portable comparison script](../../pipeline/pc12_dose_comparison.py).
- Complete paired tables and a run receipt are created under the chosen output directory.
