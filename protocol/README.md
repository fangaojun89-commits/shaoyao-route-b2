# Protocol history and current-method concordance

This is a **retrospective reporting record**, prepared on 24 October 2026 (America/Tijuana; publication may be dated 25 October UTC). It connects the original broad plan, retained internal protocol/SAP and major amendments to the methods and interpretation of article v27. It does not change scientific inputs, estimates, thresholds or candidate gates, and it does not report a new analysis.

## Read the current disposition first

The [22-item crosswalk](CROSSWALK.tsv) distinguishes completed analyses, versioned method changes, unfulfilled items, conditional modules and unavailable validation. Each item names the historical location, actual/current disposition, outcome-visibility qualification, article section and public supporting file. [Provenance](PROVENANCE.json) identifies every copied original and the exact local article/review/addendum used. Local supporting files listed there are **not all redistributed**; the published support column identifies the material that is accessible in this repository. This is a selected record of consequential research decisions, not an exhaustive engineering incident log or a newly authenticated chronology of all executions.

| Topic | Original specification or wording | Current reporting |
|---|---|---|
| Registration | Local prospective-lock language and planned public deposit | Internally specified, not publicly preregistered; current archive is post-analysis |
| Statistical unit | Recovered independent animal/culture identity required | Submitted arrays/libraries modeled; unresolved animal/pool/batch linkage remains an explicit assumption, especially G205 |
| Candidate sequence | S1 through recurrence/formula/convergence gates | S1: 0/1,771; downstream candidate tests unevaluated; source-wide results exploratory |
| Orthology | Ensembl one-to-one | Frozen strict NCBI gene-group-derived one-to-one mapping |
| SynGO | Primary resource in initial plan | QuickGO assignedBy=SynGO sensitivity; no candidate admission |
| G205 processing | mRatBN7.2/RefSeq and Subjunc in historical SAP | Separate October Rnor6/Ensembl101 HISAT2 branch and matched STAR sensitivity; older branch retained |
| Gene uncertainty | Broad historical CI reporting | Unsupported edgeR gene intervals omitted; estimator-specific voom/DESeq2/activity intervals retained |
| JGW platform | Strict Windows/Linux numerical bridge | That bridge failed; separately approved Linux reference uses its own observed/null outputs and thresholds |
| Later label references | Broadly proposed sensitivity | Retrospective method-specific references, not biological replication or proof of exchangeability |
| Independent PSS validation | Outcome-protected phenotype-linked expression holdout | Not achieved within the documented search |

The original plan also named low-RIN exclusion, expanded preparation/G298414 comparisons and cross-study/species candidate-selection sensitivities. The current paper does **not** claim that all these broad-plan items were completed. RNA-integrity metadata review is not an exclusion sensitivity; completed source-wide method comparisons are not nonempty candidate-validation tests. Removal of a defining timepoint can make the original interaction unidentifiable, so a different contrast is not silently substituted. See D15-D18 for the exact scope.

## Historical copies and their limitations

`historical/` contains unchanged bytes copied from the retained project files at the time of this release. Original dates in titles, filesystem observations and embedded hashes are historical statements, not newly verified trusted timestamps. Some files were already consolidated or amended after outcomes. In particular:

- The chronology audit and event ledger describe the 18 August snapshot; the early deviations log ends on 29 August. They do not cover the complete later history.
- The audit's statement that an earlier G277 contrasts file was not found is a statement about that audit. Later public source selections have their own provenance. No retrospective rewrite of the audit is made.
- Historical labels such as `prospective protocol lock`, `held-out`, `independent rat/library` and the `prospective_rule` column must be read with D01-D05. They are not assertions of public preregistration, authenticated animal independence or successful validation.
- The old SAP and October proposal contain uncertainty specifications superseded by later reporting corrections. The [G277 current reporting rule](../analysis/g277_execution/README.md) excludes unsupported gene uncertainty fields while preserving the other values. Source-specific numerical intervals are not interchangeable.
- The historical audit recommends clean reruns and a signed/DOI archive. This release does not execute that recommendation, claim a signed release/DOI, or establish cold restoration. Existing results are retained under the user's no-rerun instruction.
- Historical Markdown links retain original project-relative addresses and may not resolve inside the public snapshot folder. Use this guide, the crosswalk and provenance for current navigation. Historical text was not rewritten merely to repair its links.

These are research records, not executable instructions. Archiving a plan does not establish that each proposed operation occurred. The current article is identified by hash, not included as a full manuscript in this release. No author identity, author correspondence, credentials, raw failure text, raw reads or third-party article PDF is added here.

## Two different uses of E labels

The historical claims document used **E1-E6 evidence classes** (association, recurrence, perturbation, systemic exposure, direct assay and clinical evidence). The executable candidate hierarchy uses **E1-E4 sequential gates** (integrated disease/formula support, injured-spinal convergence, compound specificity, exposure/direct-assay support). An E4 exposure *class* in literature is not a passed E4 candidate *gate*. The retained gate table carries the `route_b2_sequential_gate` namespace. Neither ontology changes in this release.

## What this resolves and what remains

The public reader can now inspect the actual original planning documents and the rationale/interpretation of major changes, rather than reconstructing these from current code defaults. Scientific output files, the v27 Markdown/Word manuscript, and all pre-existing scientific repository payloads remain unchanged.

The current work supports a bounded exploratory evidence-chain paper with a completed zero-admission result. A qualified independent PSS omics holdout remains unavailable, whole-study fresh-environment reconstruction remains unverified, and author information/declarations remain intentionally deferred. These limits do not require repeating accepted analyses or starting cloud computation for this documentation release.
