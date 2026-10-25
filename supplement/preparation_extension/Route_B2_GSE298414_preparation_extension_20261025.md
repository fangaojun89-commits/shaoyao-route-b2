# GSE298414: preparation sensitivity inputs and source reconciliation

25 October 2026. **The deposited processed counts are available and structurally verified. The subsequent 12 extraction-method gene comparisons and influence checks are reported in the [results supplement](Route_B2_GSE298414_preparation_results_20261025.md).** The input-stage account below is retained; program-level comparisons remain pending. This extension addresses the preparation-sensitivity gap while retaining the accepted canonical JGW03 analyses and their candidate rules.

## Material obtained and reused

The official GEO archive contains 141 RSEM gene-result tables, all with the same ordered 32,883 gene identifiers. Each gzip member was decoded and its expected-count vector checked for finite, nonnegative values and a positive total. There are no exact duplicate expected-count vectors within these 141 libraries or between them and the existing 12-library canonical PC12 subset. GSM and BioSample identifiers also do not overlap with that subset. These checks exclude those forms of exact duplication; they do not prove broader biological independence.

The archive is 82,524,160 bytes, with local SHA-256 `eb08eea1f78621476f946fec97ff677c937dad48d25ca7f5d01c8182e23e41cc`. The assembled count matrix occupies approximately 4.85 MB compressed on R. The complete gene ordering matches the existing GSE289929 annotation exactly, so that annotation is reused. RSEM expected counts remain estimates and are not rounded to integer observed fragments. No FASTQ download, realignment or recounting is needed for this source.

| Source component | Libraries | Interpretation |
|---|---:|---|
| Fourteen extract codes, three reported levels, three replicates | 126 | Preparation comparisons |
| Untreated `Non-treat(N.T)` | 3 | Explicitly untreated; not automatically a matched DMSO vehicle |
| Four named comparator drugs | 12 | Source context; not assumed negative controls |
| Total | 141 | 47 source groups with three libraries each |

The archive-level expected-count totals range from 5,399,411.88 to 31,918,019.88, and 14,129–16,532 genes have positive estimated counts per library. TPM sums range from 999,999.05 to 1,000,000.86, consistent with rounded deposited TPM values. These are input summaries, not evidence that every sample has passed a complete biological quality assessment.

## Newly recovered sample table resolves the recorded batch

The original article's Table S3 was obtained through Europe PMC's supplementary-file service after direct publisher/PMC downloads failed. Its file size (27,080 bytes) agrees with the article XML and its contents identify the correct PC12 study. SHA-256: `9e819e6e0ccf3bee9f334ed452ae07e7c2d63b967c978e146e3d36b4ff13fab3`.

The table contains 162 PC12 rows: 141 GSE298414 and 21 GSE289929. **All 141 GSE298414 libraries are labelled batch A.** Every GSE298414 title and treatment string agrees with GEO. All 126 extract rows agree on ratio and solvent. This removes the missing recorded-batch field for within-GSE298414 contrasts; a constant batch indicator cannot be estimated separately and is unnecessary in a within-batch design. No controls are borrowed across accessions.

Two extraction labels in Table S3 conflict with the GEO operation order and article code definitions:

| Sample | Code | Table S3 extraction label | Code/GEO-supported interpretation |
|---|---|---|---|
| GSM9013914 | JGE06, low reported level, replicate 2 | CEM | IEM: Table 3 defines code 06 as post-extraction mixing and GEO states extraction followed by mixing |
| GSM9013917 | JGE05, low reported level, replicate 3 | IEM | The code belongs to the single-GR extraction series in Table 2; GEO states mixing followed by extraction |

Both original labels are retained alongside the interpretation and its evidence. This is a documented source adjudication, not an author-confirmed correction. Classification follows the corroborated recipe code and GEO treatment description. The first discrepancy affects one low-level CEM–IEM comparison; any conclusion depending on that contrast requires a targeted check without that library. The second concerns a single-herb record and does not enter the first CEM–IEM comparison block.

Table S5 was also recovered. Its cross-batch correlation summary does not establish that preparation effects are batch-free and is not used to alter source labels or substitute controls.

## Comparisons that can now be prepared

GSE298414 contains water/ethanol codes 01, 02, 04, 05, 06, 07 and 08. It **does not contain JGW03 or JGE03**. For ratios 2:1 and 1:2, both co-extraction (02/04) and post-extraction mixing (06/08) are present at all three reported levels and in both solvents. This forms a 72-library partial factorial: two ratios × two solvents × two extraction methods × three ordinal levels × three replicates.

The first finite comparison block is therefore 12 direct CEM-minus-IEM contrasts, one per ratio × solvent × ordinal level. It does not require the untreated samples to serve as vehicle controls. A categorical group model on that 24-cell subset has an estimable comparison structure; sample influence, expression filtering and uncertainty still require actual fitting. All contrasts are source-wide preparation context and cannot reopen S1 or upgrade a candidate gate.

The remaining single-herb, solvent and ratio comparisons must remain explicit subsequent parts of the preparation work. The absent within-accession 1:1 CEM cell prevents this source alone from supplying a complete three-ratio extraction factorial. No full preparation-robustness result is claimed at this stage.

## Absolute dose remains unresolved

GEO and Table S3 name the PC12 extract levels 20, 100 and 500 µg/mL. The paper's methods state 4, 20 and 100 µg/mL for PC12. The new table confirms the disagreement rather than resolving it. Preserve both sources and use only low/middle/high ordering for the proposed within-source comparison. Do not infer an absolute exposure, equate those levels with the canonical source, or convert them to a clinical dose. Recipe differences at matching reported levels remain observational contrasts of those deposited conditions, not proof of dose-equivalent pharmacology or synergy.

## Files and next analysis

- Original archive and prepared inputs: `R:/Shaoyao_storage/data/GSE298414_preparation/`.
- Current sample map: `inputs_20261025/sample_ledger_supplement_verified.tsv`.
- Complete retained label disagreements: `inputs_20261025/extraction_label_conflicts.tsv`.
- First contrast inventory: `inputs_20261025/cem_iem_contrasts_supplement_verified.tsv`.
- Source article, original GEO metadata, Tables S3/S5: `R:/Shaoyao_storage/literature/GSE298414_preparation_20261025/`.
- [Count preparation code](../../analysis/preparation_extension/original/prepare_GSE298414_existing_counts_20261025.py); [source reconciliation code](../../analysis/preparation_extension/original/reconcile_GSE298414_supplement_20261025.py).

Completed since this input review: the finite 72-library gene comparison block, the disputed-IEM-library sensitivity, and six systematic omissions prompted by an influential high-level ethanol CEM library. See the results supplement for the strong method/sample dependence. Next: program-level interpretation and the remaining preparation contrasts. Existing GSE289929/JGW03 results are retained. The main manuscript remains v38 pending integration of the new evidence.

Sources: [GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298414), [official GEO archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298414/suppl/GSE298414_RAW.tar), and Baek et al., *Multidimensional transcriptome dataset for systematic evaluation of Jakyakgamcho-tang-induced cell signatures*, DOI [10.1038/s41597-026-06759-6](https://doi.org/10.1038/s41597-026-06759-6), especially Tables 2/3/9 and Table S3. Original third-party materials remain internal; no new public redistribution occurred.
