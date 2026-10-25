# GSE277064 manuscript-ready text

## Methods

### Cell-type-resolved post-stroke cervical-spinal translatome analysis

GSE277064 was used to characterize the cell-type and temporal distribution of post-stroke cervical-spinal molecular programs, not to validate a post-stroke-spasticity phenotype. The accession comprises 63 RiboTag RNA-seq libraries generated from V2a interneurons, astrocytes, and microglia. The primary natural-course analysis was restricted to 45 libraries: control and post-stroke days 3, 7, 14, and 28, with three independent libraries per cell-type/time stratum. The 18 libraries collected 2, 4, or 8 h after rotarod training on post-stroke day 14 were excluded from all natural-course models. They were analyzed separately, within V2a interneurons and astrocytes, against the corresponding untrained day-14 group.

The three official GEO matrices contain RSEM expected counts generated against GRCm39/m32. Fractional expected counts arising from probabilistic allocation of multi-mapping reads were retained without rounding. Count-matrix columns were reconciled against the frozen sample ledger and official GEO sample titles before modeling. Ensembl gene identifiers were mapped with org.Mm.eg.db (version 3.23.0) under a fixed first-mapping rule, and the complete mapping table was retained.

Each cell type was analyzed independently. Genes were retained when they passed edgeR `filterByExpr` under the complete design and had at least 1 count per million in at least three libraries. Libraries were normalized with the trimmed mean of M values. Robust edgeR quasi-likelihood models were fitted as `~0 + timepoint`. The prespecified primary contrast was the mean of days 14 and 28 relative to control; day-14 and day-28 component contrasts were retained, and the mean of days 3 and 7 served as an early-stage control. Gene-level results are reported as log2 fold changes, quasi-likelihood statistics, 95% confidence intervals, and Benjamini-Hochberg false-discovery rates.

Reactome release 97 programs were transferred from human to mouse through the frozen strict one-to-one NCBI Gene Group map and restricted to 10–500 expressed features. `camera` was applied to all prespecified contrasts. The chronic-combined contrast was additionally evaluated with 9,999-rotation `mroast` and `fgseaMultilevel`; exact preranked-statistic ties were resolved deterministically by gene symbol with machine-epsilon perturbations. Mouse-native PROGENy 1.34.0 top-500 footprints and DoRothEA 1.23.0 confidence A/B regulons were scored with `decoupleR` univariate linear models and tested with the same design. Multiplicity was controlled within prespecified collections and contrasts. For exploratory all-program summaries, a deliberately conservative global adjustment covered program-by-cell-type chronic tests. Leave-one-library analyses quantified directional stability of chronic program effects. No sample was excluded on the basis of PCA position or desired effect direction.

### Sequential-evidence rule

GSE277064 could contribute S3 localization only for a program that had already passed S2 and whose disease direction had been frozen upstream. The frozen S1 analysis retained no programs; consequently, the eligible S2 input to GSE277064 was empty. S3 was therefore recorded as `NOT_IDENTIFIABLE_NO_S2_INPUT`, with no candidate-level test or pass count. All-program results below are exploratory descriptions of post-stroke cell-type/time biology and cannot create or rescue an S3 candidate.

Here, “prespecified” and “locked” denote an internally specified, locally frozen analysis workflow. The study was not publicly preregistered. This terminology clarification changes no sample-inclusion rule, contrast, method, threshold, numerical result, or S3 gate decision.

## Results

### Input integrity, quality control, and model estimability

All 63 count-matrix columns matched the frozen GSM ledger and official GEO titles exactly. The natural-course population contained 45 libraries, with the expected three biological replicates in each of 15 cell-type/time strata; no rotarod library entered a natural-course model. All three design matrices were full rank (rank 5; residual degrees of freedom 10). No library met the prespecified objective technical-warning criteria, and no sample was excluded.

After filtering, 15,511, 13,842, and 14,547 features were retained in V2a interneurons, astrocytes, and microglia, respectively. The corresponding biological coefficients of variation were 0.078, 0.084, and 0.069. In the chronic-combined contrast, gene-level BH FDR <=0.05 was reached by 6,299 V2a, 2,395 astrocyte, and 2,862 microglial features. These large within-cell-type responses indicate substantial post-stroke translatome remodeling but do not confer phenotype specificity.

### Exploratory cell-type and stage patterns

At the PROGENy footprint level, chronic V2a remodeling included increased Estrogen (standardized effect 0.534; global exploratory q=1.67e-4), MAPK (0.556; q=7.13e-4), and EGFR (0.442; q=9.83e-3) activity, together with reduced Hypoxia (-0.508; q=2.37e-4) and Androgen (-0.374; q=9.83e-3) activity. Estrogen, Hypoxia, and MAPK were not detected in the early-combined stage control (within-contrast q=0.63 for each), whereas the Androgen effect was already present early. This temporal contrast is compatible with, but does not establish, chronic-stage emergence of selected V2a programs.

Astrocytes showed lower chronic Androgen (-0.643; q=3.69e-3), TGF-beta (-0.557; q=8.14e-3), Estrogen (-0.384; q=1.19e-2), and MAPK (-0.236; q=1.30e-2) footprint activity and higher p53 (0.316; q=2.88e-2) and Hypoxia (0.226; q=4.82e-2) activity. Early-stage effects were heterogeneous, including an opposite-signed MAPK estimate; the chronic contrast therefore should not be summarized as a uniform astrocytic activation or suppression state.

Microglia showed marked chronic NF-kappaB (1.923; q=9.73e-10), TNF-alpha (1.460; q=3.72e-8), JAK-STAT (1.747; q=4.10e-6), and WNT (0.731; q=2.70e-5) activity. DoRothEA results were concordant for Nfkb1, Rela, Stat1, Stat2, and related inflammatory regulators. However, the corresponding PROGENy effects were already strong, and generally larger, in the day-3/day-7 stage control. These findings therefore mark persistent post-stroke immune remodeling rather than a chronic- or spasticity-specific program.

Across all locked collections, stringent exploratory consistency criteria were met by 296, 41, and 45 Reactome programs in V2a interneurons, astrocytes, and microglia, respectively; by 5, 6, and 4 PROGENy footprints; and by 18, 25, and 45 DoRothEA regulons. The Reactome counts include nested and overlapping pathways and must not be interpreted as independent discoveries. Every retained exploratory effect preserved its direction in at least 80% of leave-one-library analyses; the stricter summarized set contained no leave-one-library sign reversal. These quantities describe internal recurrence within GSE277064 only. They are neither candidate counts nor evidence grades.

### Rotarod boundary domain

The rotarod branch comprised, per cell type, three untrained day-14 libraries and nine trained libraries sampled at 2, 4, or 8 h. Relative to untrained day 14, within-contrast BH q<=0.05 was reached by 844, 417, and 327 Reactome programs in V2a interneurons at 2, 4, and 8 h, respectively, and by 1,090, 1,054, and 1,050 astrocyte programs. These extensive responses demonstrate that recent motor training materially alters the RiboTag translatome. Because the trained samples represent a recovery/activity perturbation without a spasticity endpoint, they were retained solely as a boundary-domain analysis and did not enter S1–S3 or E1–E4.

### Evidence conclusion

The analysis establishes a reproducible atlas of cell-type- and time-resolved cervical-spinal translatome changes after stroke. It does not validate post-stroke spasticity, Shaoyao–Gancao responsiveness, treatment efficacy, or causal mechanism. Because the internally specified and locally locked upstream candidate set is empty, the correct S3 result is not zero positive localizations but non-identifiability owing to the absence of eligible S2 input.
