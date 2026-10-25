# Computational protocol: Route B2

Version: 0.1.0  
Date: 2026-08-17  
Status: pre-analysis draft; freeze after metadata audit and before outcome-bearing analysis

## 1. Study type and estimand class

This is a multi-source, retrospective computational evidence-triangulation study. It is exploratory and prioritization-oriented rather than confirmatory or causal. The inferential unit is the independent biological replicate defined by each source study.

The principal estimands are pathway- or regulon-level standardized effects within each source design. Cross-study evidence will be represented as a structured evidence vector and direction/stability profile, not as a pooled expression matrix or a single synthetic p value.

## 2. Primary research question

Which molecular programs associated with progression of experimentally measured spasticity recur in chronic post-stroke cervical-spinal remodeling and are directionally opposed by an exact whole-decoction perturbation in a prespecified neural-like cell context, with independent convergence under paeoniflorin–liquiritin treatment in injured spinal tissue?

## 3. Secondary research questions

1. Which candidate programs are localized to V2a interneuron, astrocyte, or microglial translatomes after stroke?
2. Are formula-responsive signatures consistent across dose while remaining distinguishable from generic transcriptional stress or broad compound responses?
3. Does a separate peripheral muscle/ECM branch recur across contracture, stroke-muscle, and in-vivo formula datasets without being conflated with reflex hyperexcitability?
4. For a surviving program, does spatial data add localization that cannot be answered by bulk or cell-type-resolved data?

## 4. Phenotype ontology and scope

Post-stroke spasticity is treated as three separable domains:

- involuntary motor activity/reflex hyperexcitability;
- impaired voluntary motor control;
- passive tissue and mechanical remodeling.

The primary paper addresses the first domain through neurospinal programs. Passive muscle/ECM remodeling is a separately labelled secondary branch. Generic stroke severity, infarct response, motor recovery, muscle atrophy, and contracture are not interchangeable with post-stroke spasticity.

## 5. Source roles

### 5.1 Core sources

- **GSE16710:** phenotype anchor; progressive spinal-cord-injury spasticity versus time-matched sham.
- **GSE128623:** independent post-stroke cervical-spinal bulk time course; chronic recurrence and acute-phase boundary control.
- **GSE277064:** post-stroke cervical-spinal RiboTag time course in V2a neurons, astrocytes, and microglia; cell-type and temporal localization.
- **GSE289929, GSE227494, GSE299063:** exact JGW03 whole-decoction perturbations in PC12, C2C12, and HepG2 cells, always contrasted with the labelled vehicle from the same accession and batch.
- **GSE205906:** injured-spinal paeoniflorin–liquiritin intervention; independent cross-model pharmacological convergence.

### 5.2 Secondary and calibration sources

- **GSE85871:** broad-compound empirical specificity/null analysis only.
- **GSE227495:** in-vivo whole-formula peripheral-muscle bridge only.
- **GSE31243 and PRJNA630495:** muscle/ECM cross-disease and acute post-stroke context only.
- **GSE270766 and GSE168934:** stroke recovery calibration; not evidence of spasticity.
- Spatial or single-cell resources enter only after the spatial decision gate in Section 12.

## 6. Eligibility and exclusion rules

A source is eligible for inferential analysis only if:

1. biological replicate identity can be recovered;
2. intervention, comparator, tissue/cell type, time, and species are unambiguous;
3. at least two independent biological units exist per modeled condition, with three preferred;
4. no required comparison is completely confounded with study/batch;
5. raw or appropriately processed data and a traceable annotation are available.

Sources with one library per condition, unresolved pooling, pseudo-controls assembled across studies, or incompatible intervention identity may be used descriptively but not for differential or validation claims.

## 7. Preprocessing and quality control

- Freeze accession-, sample-, and file-level metadata before reading outcome-bearing matrices.
- Verify sample counts against official repository records and associated primary papers.
- Compute file checksums and preserve raw files unchanged.
- Inspect library size, detected features, duplication where available, mapping/assignment summaries, sample–sample correlation, ordination, sex/time/batch labels, and outliers.
- Exclusions require a pre-specified technical criterion or documented metadata failure; they cannot be based on desired biological direction.
- For raw RNA-seq, align or quantify against a species-correct current reference and record release identifiers. Repository species/genome labels are verified rather than trusted automatically.
- Microarray probes are mapped with a frozen platform annotation; ambiguous many-to-many probes are removed or summarized using a pre-specified rule.

## 8. Within-source analyses

### 8.1 GSE16710 phenotype anchor

Use limma with factors for spasticity status, time, and their interaction. The primary contrast is:

`(Spasticity_day60 - Spasticity_day21) - (Sham_day60 - Sham_day21)`

The interaction isolates progression beyond time-related change in sham animals. Report standardized gene and gene-set effects, moderated statistics, false-discovery rates, confidence intervals where estimable, and leave-one-animal-out stability.

### 8.2 GSE128623 post-stroke cervical spinal cord

Use edgeR quasi-likelihood as the primary RNA-seq model, with DESeq2 or limma-voom as sensitivity analysis. The primary chronic estimand contrasts the average of days 28 and 42 with control, while preserving day as a modeled factor. Days 4 and 7 form an acute-phase boundary analysis; they are not combined with chronic samples.

### 8.3 GSE277064 cell-type-resolved post-stroke translatomes

Analyze each cell type with animal/library as the unit. Estimate time-specific effects against its cell-type-matched control, emphasizing days 14 and 28. Do not pool cell types or treat translated genes as equivalent to whole-tissue abundance. Exercise-response samples are secondary and analyzed within their original design.

### 8.4 Exact whole-formula perturbations

Each accession and cell line is modeled independently against its same-accession, same-batch vehicle. Ordered dose response is preferred; high-versus-vehicle and monotonic trend are co-primary descriptive estimands. When absolute dose metadata conflict, use only verified ordinal labels (low/mid/high) and flag the source until author clarification. HepG2 serves as tissue-context calibration and cannot establish neurospinal relevance.

### 8.5 GSE205906 injured-spinal intervention

Estimate injury versus sham and treatment versus untreated injury separately. A candidate is considered directionally opposed only when treatment moves the injury-associated program toward sham and the effect is not created by a single sample. Because this is peripheral nerve injury, it upgrades cross-model convergence but not stroke or spasticity specificity.

## 9. Gene-set and regulatory analyses

Primary inference is gene-set-first to reduce cross-species instability and target-selection bias. Pre-specified families are Reactome pathways, SynGO synaptic ontology, PROGENy pathway activities, and DoRothEA regulons. MatrisomeDB is reserved for the secondary passive-tissue branch.

Competitive gene-set testing with camera is primary; rotation testing with roast and rank-based fgsea are sensitivity analyses. Gene-set size is calculated after species and platform mapping. Multiple testing is controlled within each pre-specified family and contrast. Data-driven clusters may be reported as exploratory and must be re-tested in held-out sources.

## 10. Cross-species and cross-platform transfer

- Gene-level transfer uses traceable one-to-one Ensembl orthologs.
- Pathway/regulon effects are the primary transferable object.
- No expression matrices from different species, tissues, platforms, or studies are batch-corrected into a single pseudo-cohort.
- Directional recurrence, effect-size concordance, and rank-based overlap are reported separately. Rank overlap is not described as organ coupling, synchrony, resonance, or causal communication.

## 11. Evidence gates and candidate reporting

A program is classified through sequential gates rather than by a forced top-N list:

1. **S1 phenotype association:** stable signal in the GSE16710 progression interaction.
2. **S2 post-stroke recurrence:** directionally compatible chronic signal in GSE128623.
3. **S3 cell/time localization:** compatible signal in at least one pre-specified cell type and chronic/recovery time in GSE277064.
4. **J formula modulation:** coherent response to exact JGW03 in PC12, the prespecified neural-like formula-perturbation context. C2C12 and HepG2 are context analyses and cannot substitute for the PC12 gate.
5. **E2 intervention convergence:** injury-associated signal opposed by paeoniflorin–liquiritin in GSE205906.
6. **Specificity/exposure upgrades:** empirical specificity, measured exposure, or direct target evidence can upgrade but are not inferred from docking.

The main candidate set requires S1, S2, the locked PC12 formula-response gate (J2), and the independent injured-spinal intervention gate (E2). S3 upgrades a candidate with held-out cell/time localization but cannot substitute for S2. Sources that share a parent experiment cannot be counted as independent confirmations. The number of candidates may be zero.

## 12. Conditional advanced-method gates

### Spatial transcriptomics

Spatial analysis is performed only if a surviving program poses a localization question not resolved by GSE277064. Eligible spatial studies must have recoverable animal and slide identifiers. Spots are nested observations, not biological replicates. Multimodal measurements from the same animals count as one biological experiment.

### Virtual perturbation

Virtual knockout or perturbation is supplementary only. It requires a context-matched single-cell trajectory or multiomic reference, animal-level train/test separation, at least one context-matched real perturbation for calibration, and superiority over no-change, mean, linear, and additive baselines. Model output is a prediction, never a validation or causal intervention.

### Structural modeling

Docking and molecular dynamics are excluded from the default workflow. They may be added only after a candidate has relevant measured exposure, direct biochemical/functional support, a defensible binding site, and a validated experimental structure or explicitly qualified structure model. Even then, results support structural plausibility only.

## 13. Robustness and falsification

- leave-one-animal/sample-out stability;
- alternative count models and gene-set methods;
- ortholog-mapping sensitivity;
- acute-versus-chronic negative control;
- context negative controls and broad-compound empirical nulls;
- label permutations that preserve the original design;
- matched random gene-set and degree-preserving network nulls;
- leave-one-source-out evidence profiles;
- a multiverse/specification curve for all defensible preprocessing and modeling decisions;
- explicit dependency matrix to prevent double-counting related datasets.

## 14. Reporting and governance

The manuscript will report accession-level provenance, sample exclusions, all primary and sensitivity specifications, effect sizes, uncertainty, multiplicity control, and null findings. Code, environment locks, metadata corrections, and a machine-readable results ledger will accompany the paper, subject to repository licenses and data-use conditions.

No human-identifiable data are anticipated. Any future controlled-access or clinical data require separate ethics, consent/data-use, privacy, preregistration, and governance review.

## 15. Stopping and revision rules

- If biological replication cannot be established, the source is downgraded to descriptive use.
- If exact formula identity cannot be verified, the source cannot satisfy the formula gate.
- If no program passes the full gate, the paper reports the negative triangulation and identifies which bridge failed.
- A protocol amendment is permitted only with date, rationale, affected analyses, and a designation of pre-specified versus post-hoc work in `deviations.md`.
