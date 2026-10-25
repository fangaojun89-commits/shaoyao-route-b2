# Route B2 Statistical Analysis Plan

**Project:** Shaoyao–Gancao-related perturbations and molecular programs relevant to post-stroke spasticity  
**SAP version:** 1.0.0  
**Protocol date:** 2026-08-17  
**Status:** prospective protocol lock; no results are entered in this document  
**Machine-readable contrasts:** `../config/contrasts.yaml`

## 1. Study objective and inferential ceiling

This is a phenotype-anchored, perturbation-grounded and exposure-constrained computational prioritization study. The primary objective is to identify molecular programs that:

1. change with progression of experimentally observed spasticity in GSE16710;
2. recur during chronic post-stroke cervical-spinal remodeling in GSE128623, with cell-type/time localization assessed in GSE277064;
3. show a prespecified, directionally coherent response to canonical 1:1 hot-water co-extracted Jakyakgamcho-tang (JGT) in PC12 cells (GSE289929); and
4. receive independent directional support from paeoniflorin–liquiritin treatment in injured rat spinal cord (GSE205906), followed by empirical specificity calibration in GSE85871.

The estimands are within-study expression/program contrasts. They do **not** estimate clinical efficacy, a patient-level treatment effect, a post-stroke-spasticity-specific molecular signature, target engagement, or a causal mechanism. “Reversal” means transcriptomic directional opposition only. Zero retained programs is an admissible outcome.

## 2. Analysis populations, roles and statistical units

| Resource | Locked role | Statistical unit | Primary inclusion | Explicit non-use |
|---|---|---|---|---|
| GSE16710 | Spasticity-enriched cross-disease phenotype anchor | independent rat/array | all 27 arrays with recoverable group labels | not called PSS-specific |
| GSE128623 | Post-stroke cervical-spinal temporal recurrence | independent mouse/library | control, D4, D7, D14, D28 and D42 | no cross-study batch merging |
| GSE277064 | Held-out cell-type/time localization | independent mouse RiboTag library | natural-course control/D3/D7/D14/D28, stratified by V2a, astrocyte and microglia | rotarod samples are not pooled with natural course; cell types are not pooled |
| GSE289929 | Primary whole-formula perturbation | independent PC12 culture/library | same-accession vehicle and JGW03 1:1 hot-water CEM at 4, 20 and 100 µg/mL | ethanol extract is sensitivity only; no controls borrowed from another accession |
| GSE227494 | Secondary cellular-context analysis | independent C2C12 culture/library | same-accession vehicle and 1:1 hot-water CEM JGT | not external validation of GSE289929 |
| GSE299063 | Off-tissue/context negative control | independent HepG2 culture/library | same-batch matched vehicle and JGW03 only | not neural efficacy evidence |
| GSE205906 | Held-out diseased-spinal pharmacologic bridge | independent rat/library | sham, SNI and SNI+PF–LQ (PLM), n=3/group | repository processed genome annotation is not accepted for final inference |
| GSE85871 | Same-platform empirical compound null | independent culture replicate/profile | prespecified PF, liquiritin and glycyrrhizic acid plus eligible background compounds | not tissue relevance, PSS validation or direct-target evidence |

Cells, reads, genes, probes, image spots, pathway members and molecular-dynamics frames are never treated as biological replicates. If technical replicates exist, they are collapsed to the biological unit before inference. Repeated tissues/regions from the same animal require a paired term or animal random intercept; none may be relabeled as independent.

## 3. Data freeze, provenance and blinding

Before outcome modeling, `metadata/samples.tsv` must contain accession, GSM, BioSample/SRA where applicable, biological-unit ID, species, tissue/cell line, condition, time, dose, formula code, ratio, solvent, extraction mode, exposure duration, vehicle, batch and inclusion status. Duplicate BioSamples/runs across accessions are counted once. The sample ledger and checksums are frozen before program testing.

Sample labels are used to construct the design but candidate program names and drug-response results remain masked while technical QC exclusions are adjudicated. An exclusion requires a documented technical failure established without reference to the desired effect direction. PCA position alone is insufficient. All exclusions and protocol deviations are appended to `protocol/deviations.md`; the original ledger is never overwritten.

## 4. Common preprocessing and quality control

### 4.1 Microarray (GSE16710)

- Prefer raw CEL files; apply RMA background correction, quantile normalization and log2 summarization using a pinned Bioconductor environment.
- Fit arrayQualityMetrics-style diagnostics, RLE/NUSE where supported, sample correlations and PCA. A sample is removed only for a primary technical defect corroborated by at least one independent diagnostic; otherwise it remains and is examined by leave-one-animal influence analysis.
- Primary probe-to-gene mapping uses a current, version-locked annotation. Original platform annotation is a mandatory sensitivity analysis. Multiple probes per gene are not selected by outcome; program tests operate on probe-level statistics with duplicate-gene handling fixed in code, while displayed gene effects use the probe with the largest average expression, selected blind to group.

### 4.2 RNA-seq

- Raw counts are used; TPM/FPKM are not input to count-model differential expression.
- Genes are filtered with `edgeR::filterByExpr` using the complete design for that dataset. TMM normalization and edgeR quasi-likelihood (QL) models are primary; DESeq2 and limma-voom are prespecified sensitivity methods.
- Library size, detected features, composition, sample correlations, MDS/PCA, mean–variance behavior and model leverage are inspected. Batch is included only if recorded before analysis, estimable, and not aliased with condition. Surrogate-variable adjustment is not a primary analysis and, if justified, is reported only as sensitivity.
- There is no cross-study ComBat/ComBat-seq correction. Each study is analyzed separately and integrated only through effect direction, uncertainty, ranks and ortholog/pathway mappings.

### 4.3 Mandatory raw re-alignment of GSE205906

- FASTQ reads are reprocessed against **Rattus norvegicus mRatBN7.2** with a versioned, checksum-locked GTF. Reference and annotation release must be entered in the provenance ledger before execution; an unpinned reference is a hard failure.
- Adapter/quality assessment and strandedness diagnostics precede inference. The frozen local primary workflow uses **Rsubread 2.26.0 Subjunc** with a split gapped index and **featureCounts** exon-to-`gene_id` aggregation against the checksum-locked NCBI RefSeq annotation. Unstranded, forward-stranded and reverse-stranded assignment summaries are retained. A stranded mode is selected for a dataset only if (i) the same mode is dominant in at least 80% of libraries, (ii) the median dominant-strand assigned count is at least 80% of the unstranded assigned count, and (iii) the median minor/dominant assigned-count ratio is at most 0.20; otherwise unstranded counts are primary. This rule is applied without outcome labels. STAR two-pass plus featureCounts and Salmon remain pre-specified sensitivity analyses if a version-locked Linux/cloud runtime becomes available; unavailable sensitivity methods are reported as not run and are never imputed as agreement.
- A matched dual-reference rat-versus-human species screen is required because repository metadata contain a GRCh38 annotation concern. For each GSE205906 library, the same deterministic first 500,000 read pairs are uniquely mapped with identical Subjunc settings to checksum-locked mRatBN7.2 and GRCh38.p14 indexes. This is a gross species-identity diagnostic, not a contamination-rate estimator. Every library must have rat unique mapping ≥0.60, rat-minus-human unique mapping ≥0.10, and a rat/human mapping-rate ratio ≥1.10; the cohort median rat-minus-human difference must be ≥0.20. A failed/ambiguous prefix screen triggers full-library investigation and suspends all GSE205906 inference. Analysis also stops if group labels cannot be reconciled with SRA/BioSample metadata.
- The submitted processed matrix may be used only for a diagnostic concordance plot, not the final effect estimates.

### 4.4 Missingness and outliers

No expression values, group labels, times or doses are imputed. A sample with missing information required for its contrast is excluded from that contrast only. Robust empirical-Bayes fitting and influence diagnostics are used; a statistically influential but technically valid animal is retained in the primary analysis and removed only in a labeled sensitivity analysis.

## 5. Feature spaces and multiplicity families

The primary inferential unit is a predefined molecular program, not a DEG Venn intersection or post hoc PPI hub. Before unmasking drug data, versions and membership are frozen for:

- Reactome and SynGO programs in the neurospinal domain;
- PROGENy pathways and high-confidence DoRothEA regulons as separate regulator families.

Eligible gene sets contain 10–500 mapped genes in a dataset; the threshold is applied independently of outcomes. Ensembl high-confidence one-to-one orthologs are the primary gene-level bridge. Pathway-level transfer is primary for cross-species comparison; one-to-many mappings are sensitivity only. Mapping coverage is reported for every program.

Gene-level differential-expression tables are exploratory and use BH FDR within each dataset/contrast. Program-level tests use `camera`/`roast` as primary and `fgseaMultilevel` as sensitivity. Sample-level activity estimates used for effect sizes are fit with the same design as the gene model. PROGENy/DoRothEA activities are computed with a pinned `decoupleR` implementation and are adjusted within their own prespecified family.

Unless a section specifies otherwise:

- two-sided tests and 95% confidence intervals are reported;
- BH FDR is controlled at 0.05 within each explicitly named family;
- an unadjusted P value cannot pass a gate;
- all effect estimates, exact n, confidence intervals and adjusted P values are retained, irrespective of significance;
- sequential held-out families include only candidates admitted by the preceding locked gate, preventing rediscovery in the validation resource.

## 6. Dataset-specific estimands and models

### 6.1 GSE16710: phenotype-anchor interaction

Primary model (limma, robust empirical Bayes if diagnostics support it):

```text
expression ~ injury + time + injury:time
```

Reference levels are Sham and day 21. The primary estimand is:

```text
(Spi60 − Spi21) − (Sham60 − Sham21)
```

Supportive estimands are Spi21−Sham21, Spi60−Sham60, Spi60−Spi21 and Sham60−Sham21. The primary program family is all eligible locked Reactome/SynGO programs for the interaction. PROGENy and DoRothEA are separate families. The disease direction for each admitted program is frozen from the interaction coefficient before any JGT analysis.

Mandatory diagnostics/sensitivities are current versus original probe annotation, leave-one-animal influence, and the two time-specific injury contrasts. A program cannot pass S1 if its interaction direction changes under outcome-independent annotation or is driven by a single animal.

### 6.2 GSE128623: chronic versus acute post-stroke cervical cord

Primary model:

```text
counts ~ 0 + timepoint
```

The co-primary chronic estimands are `0.5×D28 + 0.5×D42 − control` and the two component contrasts D28−control and D42−control. The acute phase-control estimands are `0.5×D4 + 0.5×D7 − control`, D4−control and D7−control; D14−control is transitional.

Only S1 programs enter the held-out chronic family. Chronic support requires the combined chronic contrast to have the locked disease direction, BH FDR ≤0.05, a 95% CI excluding zero, and the same direction in both D28 and D42 component estimates. Acute significance may coexist, but acute effects cannot substitute for chronic support. All animals remain independent; timepoint is not modeled as repeated measurement.

### 6.3 GSE277064: cell-type/time-stratified RiboTag validation

V2a interneurons, astrocytes and microglia are analyzed in three separate models:

```text
counts ~ 0 + timepoint
```

Natural-course samples are control, D3, D7, D14 and D28. Primary held-out contrasts for S2 programs are D14−control, D28−control and `0.5×D14 + 0.5×D28 − control`; D3/D7 are stage controls. Multiplicity is controlled across admitted program × cell type × chronic contrast in one family. A cell type provides S3 localization only when the combined chronic estimate matches the locked disease direction, BH FDR ≤0.05, its 95% CI excludes zero, and D14/D28 do not have opposite signs. Effects, CIs and rank stability are emphasized because n=3/group.

Rotarod samples at 2/4/8 h on post-stroke day 14 are analyzed separately within V2a and astrocytes, using untrained D14 as the comparator only if metadata establish compatibility. They are boundary-domain recovery analyses and never determine S1–S3 or E1–E4.

### 6.4 GSE289929: canonical neural-like whole-formula perturbation

The primary subset is PC12, 1:1, hot-water, co-extracted JGW03, 24 h, with same-accession vehicle. Official GEO doses are 4, 20 and 100 µg/mL. The exact GSM sets are locked in `contrasts.yaml`.

Primary estimands are:

1. 100 µg/mL JGW03 versus matched vehicle;
2. a monotonic ordered-dose trend coded vehicle/low/mid/high = 0/1/2/3.

`log10(dose+1)` is a sensitivity encoding. If a verified, non-confounded batch exists, it is included; otherwise no batch term is invented. Only S2/S3 programs are tested. By default the expected therapeutic polarity is opposite to the locked disease direction. A “reinforcement” polarity is permitted only for a program designated as adaptive, with external justification recorded before the drug data are unmasked.

J1 requires a directionally coherent high-dose estimate. J2 additionally requires a coherent ordered-dose estimate; both estimates must have BH FDR ≤0.05 and 95% CIs excluding zero in their locked candidate family, and no replicate may be excluded for directional convenience. Ethanol JGT is sensitivity only.

### 6.5 JGT matched cellular-context controls

C2C12 (GSE227494) and HepG2 (GSE299063) are modeled separately. For every comparison, controls must match accession, cell line, experiment/batch, vehicle, exposure duration and other recoverable design factors. Controls are never borrowed across accessions. Canonical sample IDs are specified in `contrasts.yaml`; any disagreement with the frozen sample ledger stops that contrast.

High-dose and ordered-dose estimates are calculated as for PC12. C2C12 asks whether the response is compatible with a muscle context; HepG2 is an off-tissue/general-stress control. These resources can upgrade J3/J4 but are not independent replication of the whole-formula source. Ratios, water versus 70% ethanol, co-extracted versus separately extracted/mixed, and single herbs are factorial/stratified sensitivity analyses only. A co-extract deviation from a weighted single-herb expectation is called a non-additive transcriptomic response, never synergy.

### 6.6 GSE205906: independent injured-spinal PF–LQ perturbation

Primary edgeR QL model after raw rat re-alignment:

```text
counts ~ 0 + group
```

The primary article specifies 3.98 mg/mL liquiritin, whereas GEO labels the three PLM replicates as 3.99, 3.98 and 3.100 mg/mL. These values are treated as an unresolved repository/source conflict, not as three intended dose levels. The samples may support the prespecified directional PLM-group contrast after identity and raw-read QC, but no absolute-dose, dose-response, component-ratio or synergy inference is permitted.

The prespecified contrasts are SNI−sham (local disease effect) and PLM−SNI (intervention effect). Only E1 candidates are tested. Independent support requires: (a) SNI−sham does not contradict the locked disease direction; (b) PLM−SNI opposes both the locked disease direction and the local SNI effect; (c) the intervention estimate has BH FDR ≤0.05 and a 95% CI excluding zero. Given n=3/group, nominal P values alone are never described as confirmation; leave-one-animal and alignment-method stability are mandatory.

### 6.7 GSE85871: empirical specificity calibration

Paeoniflorin, liquiritin and glycyrrhizic acid are fixed before analysis. Each eligible compound signature must be constructed against its own metadata-compatible vehicle/reference. If a valid reference cannot be recovered, that compound is not tested.

For every E2 program, calculate a polarity-aligned reversal score for all eligible compounds. Signature magnitude is the root-mean-square moderated t statistic over the common tested genes. Each target compound is compared with its 20 nearest non-target compounds by magnitude; nearest 10 and 30 are sensitivity specifications. Empirical one-sided P values use the plus-one correction. BH adjustment covers E2 program × three prespecified compounds. Animal/culture and gene bootstraps (2,000 draws, fixed seed) yield 95% percentile intervals, rank distributions and top-5%/top-10% stability. This analysis addresses relative specificity only.

## 7. Sequential disease, formula and final evidence gates

The complete gate logic is evaluated exactly as configured and is never relaxed after results are seen.

| Grade | Prospective rule |
|---|---|
| S1 | GSE16710 interaction passes the locked program family and its sign is stable to annotation and influence analyses. |
| S2 | S1 plus GSE128623 chronic support; D28 and D42 have the same sign and the program is not supported solely by D4/D7. |
| S3 | S2 plus held-out GSE277064 chronic localization in at least one prespecified cell type under the multiplicity rule above. |
| J1 | A locked S2/S3 program has the prespecified high-dose PC12 JGW03 direction. |
| J2 | J1 plus a coherent PC12 ordered-dose trend, with both tests satisfying the locked FDR/CI rule. |
| J3 | J2 plus interpretable C2C12/HepG2 context discrimination rather than a stronger concordant off-tissue stress response. |
| J4 | J2/J3 remains directionally stable across admissible preparation/ratio sensitivity specifications. |
| E1 | S2 or S3 **and** J2. This is exploratory formula-responsive evidence. |
| E2 | E1 plus independent GSE205906 injured-spinal directional support. This is the minimum grade for a main-text candidate program. |
| E3 | E2 plus GSE85871 empirical specificity: empirical BH q≤0.05 and ≥80% bootstrap placement in the top 10% for at least one prespecified compound. |
| E4 | E3 plus separately curated, identity-matched in vivo exposure and direct biochemical/functional assay evidence. E4 is not inferred from transcriptomics. |

No weighted “total evidence score” is used. Evidence derived from overlapping samples, a shared experiment or a downstream transformation of the same model is counted once in the dependency matrix. Main-text biological candidates are limited to E2–E4; E1 is explicitly exploratory.

## 8. Confidence intervals, effect reporting and calibration

- Gene effects are log2 fold changes with model-based 95% CIs. Program/regulator effects are reported on their native standardized scale with model-based or animal-level bootstrap 95% CIs; scales are not pooled across incompatible assays.
- For very small groups, the moderated model estimate is primary and a leave-one-animal range is shown. A narrow gene-level P value is not allowed to overrule unstable sample-level behavior.
- Directional concordance is always accompanied by both component effects, not only a sign indicator.
- Empirical P values use `(1 + number at least as extreme)/(1 + number of valid null draws)`.
- Null calibration includes observed type-I tail rates and QQ plots; systematic anti-conservative behavior terminates the affected family.

## 9. Negative controls and falsification tests

1. **Stage control:** GSE128623 D4/D7 and GSE277064 D3/D7 test whether a candidate is merely an acute injury/inflammation signal.
2. **Off-tissue control:** matched HepG2 JGW03 tests whether reversal is a ubiquitous high-magnitude stress response.
3. **Compound null:** non-target GSE85871 compounds, matched on signature strength.
4. **Gene-set null:** 10,000 random sets matched on size and mean expression; network analyses, if any, additionally match node degree.
5. **Label null:** 10,000 design-respecting permutations, or exhaustive enumeration when fewer valid permutations exist. GSE16710 uses Freedman–Lane residual permutation for the interaction; RNA-seq permutations preserve batch/strata.
6. **Direction null:** random sign vectors and polarity permutations of equal dimension.
7. **Method baselines:** pathway findings must be compared with no-change and simple linear/rank baselines; a more complex model does not gain evidential status unless it improves out-of-sample calibration.

Null analyses estimate calibration and specificity; they do not replace the original biological contrast.

## 10. LOSO and multiverse robustness

The primary specification is fixed above. Sensitivity results cannot replace it. The multiverse contains only the prespecified choices in `contrasts.yaml`, including:

- current versus original GSE16710 probe annotation;
- edgeR QL versus DESeq2 versus limma-voom;
- camera/roast versus fgseaMultilevel;
- one-to-one ortholog versus pathway-only transfer;
- combined chronic contrast versus agreement of component timepoints;
- GSE277064 combined D14/D28 versus individual timepoints;
- PC12 ordinal dose versus `log10(dose+1)`;
- full versus objectively defined low-RIN exclusion;
- STAR/counting versus Salmon quantification for GSE205906;
- GSE85871 nearest-10/20/30 magnitude-matched nulls.

For each candidate, report selection frequency, sign frequency, median effect/rank, 95% range, Kendall’s tau between specifications, and top-k Jaccard overlap. A retained E2/E3 candidate must pass at least 80% of valid specifications and preserve the required direction in at least 90%. Specifications invalidated by rank deficiency or missing metadata are reported and excluded from the denominator; they are not silently replaced.

Leave-one-study-out (LOSO), leave-one-species-out and leave-one-timepoint-out analyses quantify dependence; they do not manufacture an alternate evidence grade when removal makes a gate structurally impossible. In that case the output is explicitly “not identifiable without source X,” and the remaining rank/direction stability is reported. A candidate dependent on a single animal, a single chronic timepoint or a single admissible preprocessing method is removed.

## 11. Stop rules

Analysis stops or is downgraded as follows:

1. Unrecoverable biological units, controls, group labels, duplicate identities or a rank-deficient design: remove the affected contrast; no guessed metadata.
2. Failure to resolve GSE205906 species/reference identity: no GSE205906 inference and therefore no E2 claim.
3. No S2 program: stop the drug-mechanism narrative; report the evidence map/negative result only.
4. No J2 program: do not call any program whole-formula responsive and do not proceed to E2 testing.
5. No E2 program: select no targets; redirect to a negative/resource report or stop result-oriented submission.
6. Empirical null not exceeded: do not claim compound specificity or E3.
7. Formula identity, direct assay or plausible measured exposure absent: E4 is unavailable and docking/MD is deleted.
8. GSE298414 remains excluded from absolute-dose inference unless its dose conflict is prospectively resolved; ordinal-only sensitivity cannot upgrade a gate.
9. Any virtual-perturbation model without context-matched real-perturbation calibration or superiority to simple baselines is deleted completely.
10. Spatial analysis is deleted unless an E2 candidate creates a prespecified cortical lesion-border/white-matter hypothesis with biological replicate slides.
11. Zero programs passing the complete rules is reported without changing thresholds, gene sets, directions or candidate definitions.

## 12. Reproducibility and reporting

All software versions, references, annotations, seeds, commands, QC exclusions and deviations are recorded. Primary scripts must read `config/contrasts.yaml`, validate every expected level and sample ID against the frozen ledger, and fail closed on a mismatch. Random procedures use the seed in that file. Raw data are read-only; all derived objects are reconstructible by workflow.

Results will report exact biological n, all prespecified contrasts, effect sizes, 95% CIs, BH-adjusted P values, null calibration, exclusions, failures and the complete candidate flow. Language is limited to association, recurrence, localization, directional modulation and computational prioritization. No analysis in this SAP can establish that SGD treats post-stroke spasticity.
