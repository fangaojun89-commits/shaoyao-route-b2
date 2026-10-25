# Route B2 analysis chronology and registration-status audit

**Audit snapshot:** 2026-08-18 08:16:07 CST (UTC+08:00)  
**Scope:** the supplied `route_b2` protocol/configuration files, deviation log, source code, run logs, manifests, summaries and machine-readable gate outputs.  
**Non-intervention statement:** this audit created only the two chronology documents requested. It did not edit, recompute or reinterpret any numerical result or gate output.

## 1. Core determination

Route B2 was **not shown to have been publicly preregistered**. No OSF/PROSPERO/AsPredicted/registry identifier, DOI, signed release, trusted timestamp, version-control commit, or other immutable public record was found in the supplied project. The correct status is therefore:

> **Internally specified and locally time-stamped, with branch-specific amendments; not publicly preregistered.**

This conclusion is also consistent with the current upstream redesign report and research plan, both of which explicitly state that outcome review had already begun and that retrospective preregistration must not be claimed. Those current reports were themselves updated after outcome availability, so they corroborate the terminology ceiling but do not create an earlier registration record.

The local record gives useful evidence that several core decisions were documented before the recorded inferential runs: the phenotype-anchor interaction, sequential S/J/E gate order, BH threshold, stability requirements, permissibility of zero candidates, stopping rules and claim ceiling. That evidence is scientifically useful but is not equivalent to public preregistration. Local `CreationTime`/`LastWriteTime` values are mutable; a SHA-256 hash identifies content but does not prove when that content existed.

The current SAP and configuration also cannot be described as untouched pre-analysis files:

- `statistical_analysis_plan.md` was created on 2026-08-17 but last modified at 06:40:57 on 2026-08-18, after the GSE16710 phenotype-anchor run had started at 06:39:08. The deviation log attributes this change to the still-unanalysed GSE128623/GSE205906 raw-processing branch. That is a reasonable **branch-specific pre-outcome implementation amendment**, but no earlier SAP snapshot survives to prove the exact diff.
- `contrasts.yaml` was last modified at 07:47:36, after the S1 null result, J-gate stop and S3 non-identifiability were known. It is a current consolidated configuration, not an unchanged preregistration artifact.
- `claims_and_evidence.md` was last modified at 07:47:39, also after those results. Its claim restrictions are appropriate, but its current form is a post-outcome governance revision.
- `deviations.md` was last modified at 08:01:43, after multiple branches had completed. It is a retrospective audit ledger whose entries state chronology relative to each affected branch; its statements are not independently timestamped proof.

The machine-readable event ledger is `../provenance/analysis_chronology.tsv`.

## 2. Evidence classes used in this audit

| Audit class | Meaning | Permissible interpretation |
|---|---|---|
| Internal design record | A local plan/configuration existed before the recorded run according to mutable filesystem metadata | “documented in local project files before the recorded run” |
| Branch-specific pre-outcome amendment | A change was documented after another branch had begun or completed but before outcome-bearing analysis of the affected branch, according to the deviation record | “amended before analysis of the affected branch” |
| Post-outcome implementation/reporting revision | A file was created or changed after relevant results were available | “added/updated after outcome inspection for reproducibility or interpretation” |
| Run/gate evidence | A log, result, validation or gate file records what was executed or adjudicated | Supports result chronology; it does not establish preregistration |
| Public preregistration | A public, immutable, independently time-stamped record predating outcome inspection | **Not demonstrated for Route B2** |

## 3. Reconstructed chronology

### 3.1 Internal design phase

Between 22:40 and 22:57 on 2026-08-17, the local project acquired the dataset inventory, protocol, thresholds, evidence-dependency ledger, formula/exposure records, SAP, contrasts configuration and 165-row sample ledger. The protocol and threshold files were last saved before any recorded inferential run. The current sample-ledger SHA-256 (`367894...bdf66`) later matched the GSE277064 run manifest.

On 2026-08-18, the gene-set configuration and frozen RDS object were finalized by 05:58:58, before the recorded GSE16710, JGW03 and GSE277064 modeling. The current frozen object SHA-256 (`0a8526...ba51`) matches the value in the GSE277064 manifest. This supports internal resource consistency, not public registration.

### 3.2 GSE16710 phenotype anchor

The QC log records preprocessing from 06:29:48 to 06:38:11. The primary model log records execution from 06:39:08 to 06:55:28. `S1_GATE.tsv` was written at 06:55:16 and contains the decisive result: 1,771 eligible programs were evaluated and zero passed.

The strongest pre-run chronology evidence for this branch is the earlier protocol and threshold file, not the current SAP as a whole. Two provenance qualifications are necessary:

1. the SAP was saved at 06:40:57, after this run had started; the deviation log says the change concerned unanalysed raw-RNA-seq branches, but no pre-change SAP snapshot survives; and
2. the current `analyze_GSE16710_spasticity_anchor.R` was saved at 07:03:41, after the core run finished. The GSE16710 file manifest was generated at 07:10–07:14 and hashes that later script. It is therefore a post-run integrity snapshot, not proof of the exact execution-time source.

These limitations do not by themselves change the S1 numerical result; they limit how strongly prospective specification and exact code provenance can be claimed.

### 3.3 JGW03 perturbation branch

File times show PC12, C2C12 and HepG2 outputs being generated from approximately 06:30 to 07:07. `J_GATE.tsv` was written at 07:07:37, after the zero-S1 result existed, and correctly reports J1–J4 as `NOT_IDENTIFIABLE_NO_S2_INPUT`, with zero eligible/evaluated inputs and `pass_n=NA`.

There is no dedicated start/end run log. The current main script was saved at 07:06:28, when some outputs already existed, and the output manifest was refreshed at 08:05:49. That manifest includes the 07:47 post-outcome `contrasts.yaml`. Accordingly:

- the three-cell-line source-wide results are exploratory/contextual;
- J1–J4 were not tested on an eligible disease-program family; and
- the current manifest is an integrity inventory, not a preregistration or exact execution-time manifest.

### 3.4 GSE277064 cell/time branch

The run log records execution from 06:56:23 to 07:18:48, beginning after S1 had already returned zero candidates. Its source-wide atlas is therefore a post-S1 exploratory analysis, not a held-out confirmation of an admitted program. `S3_GATE.tsv` correctly records `NOT_IDENTIFIABLE_NO_S2_INPUT`.

The run manifest is valuable but reveals version evolution. It records an execution-era contrasts SHA-256 of `110eb288...33d6`, whereas the current `contrasts.yaml` SHA-256 is `98ffd869...3492`; the corresponding earlier file content is absent. In addition, the current main script was saved at 07:23:04, after the 07:18:48 completion time, and the run manifest was later refreshed to hash that post-run script. Thus the manifest verifies the current analysis package, not an independently preserved execution-time source/configuration bundle.

### 3.5 Master gate ledger

The master-gate builder and ledger were created after S1, J and S3 outcomes were known. This is acceptable as transparent machine-readable adjudication because the ledger propagates the zero-S1 stopping rule and does not recode unevaluated gates as negative. It must be described as a **post-outcome adjudication/provenance ledger**, not as preregistration.

### 3.6 Raw GSE128623/GSE205906 branch

The switch from STAR/Salmon primary processing to Rsubread Subjunc/featureCounts, and the later paired-fragment/completion-contract correction, occurred after GSE16710/JGW03/GSE277064 work had begun or completed. The deviation log explicitly states that no outcome-bearing GSE128623/GSE205906 expression result had been inspected. Therefore these changes may be described as **branch-specific pre-outcome implementation amendments**, with two qualifications:

- the project as a whole was no longer outcome-blind; and
- the claim of affected-branch outcome blindness is documented internally, not independently time-stamped.

The amendments did not alter the biological contrasts, FDR threshold, statistical unit, sequential eligibility rule or claim ceiling. Any later GSE128623/GSE205906 result must cite the exact post-amendment source/configuration hash used in its final clean rerun.

## 4. Registration and reproducibility findings

### Finding 1 — no public preregistration

This is the principal terminology issue. “Preregistered,” “prospectively registered,” “registered analysis/subset,” and “as registered” are not supportable from the supplied record. The absence of a public immutable record does not invalidate the analysis, but it lowers protection against undisclosed analytic flexibility and must be stated as a limitation.

### Finding 2 — current rule files are version-consolidated, not pure pre-analysis snapshots

The protocol and threshold file are the strongest local evidence for core pre-run decisions. The SAP, contrasts and claim-governance files were subsequently edited. Because there is no Git history or preserved versioned copies, a reader cannot reconstruct every earlier content change from hashes alone.

### Finding 3 — exact execution-time source is unresolved for several completed runs

For GSE16710 and GSE277064, the current main-script save time is later than the recorded main-run completion. For JGW03, the script was saved while/after some outputs existed and there is no run log. Later manifests hash the current source. These are provenance gaps, not evidence that results were changed.

### Finding 4 — filesystem metadata and hashes answer different questions

- A local file time provides a useful chronology clue but can be changed by copying, restoring or editing. `S1_GATE.tsv` itself has a CreationTime later than its LastWriteTime, directly demonstrating that these fields are not immutable event timestamps.
- A SHA-256 hash establishes which bytes are being referenced. Without a trusted earlier publication or signature, it does not establish when those bytes existed.
- The MD5 values in `master_gate_source_provenance.tsv` currently match their files and are adequate for accidental-change detection, but the manuscript/repository should prefer SHA-256 for final provenance.

### Finding 5 — gate semantics remain correctly conservative

The chronology limitations do not justify relaxing or rerouting the hierarchy. S1 is a completed zero-pass result. S2/S3/J/E downstream gates lack eligible input and are not positive or negative tests. Source-wide post-S1 analyses remain exploratory and cannot originate candidates.

## 5. Terminology corrections required before manuscript freeze

The following replacements should be applied in the manuscript and supporting text without changing result tables:

| Avoid | Use |
|---|---|
| preregistered / prospectively registered | internally specified; documented in local project files |
| registered analysis/subset | internally specified analysis/subset |
| exactly as registered | exactly as specified in the internal SAP |
| registered upstream candidate set | internally specified upstream candidate set |
| registered sequential family | internally specified sequential family |
| deviation-registered | documented in the deviation log |

Current locations needing wording review include `results/master_gates/README.md`, `results/master_gates/master_gate_summary.tsv`, `results/GSE277064/programs/S3_GATE.tsv`, `results/GSE277064/manuscript_ready_methods_results.md`, `protocol/bulk_realign_downstream_models.md`, `results/GSE16710/GSE16710_analysis_report.md`, and the dataset-audit manuscript tables. These are terminology corrections only; the numerical outputs and gate decisions need not change.

## 6. Manuscript-ready wording

### Methods: analysis chronology and registration status

> The core phenotype-anchor contrast, sequential eligibility structure, false-discovery threshold, zero-candidate stopping rule, and claim ceiling were documented in local project files before the recorded phenotype-anchor modeling run. The study was not publicly preregistered. Resource and implementation changes were logged relative to inspection of outcome-bearing data in the affected branch; some branch-specific amendments occurred after results from other branches were available. SHA-256 hashes and a reconstructed chronology are provided as an internal audit trail, not as independent proof of public or immutable preregistration.

### Methods: version qualification

> The analysis package contains current version-consolidated protocol and configuration files. Because immutable pre-run snapshots were not retained for every revision, we distinguish core internally specified decisions from later implementation amendments and post-outcome reporting clarifications. Final numerical analyses are linked to run logs, gate outputs, validation files, and file hashes; where an exact execution-time source snapshot is unavailable, this is stated explicitly.

### Results: downstream gate status

> The phenotype-anchor gate evaluated a non-empty prespecified family and retained no programs. Consequently, downstream S2/S3/J/E gates had no eligible input and were not evaluated; source-wide analyses completed after this stopping point are reported only as exploratory or contextual results.

### Limitations

> The study lacked a public, immutable preregistration, and exact execution-time source snapshots were not retained for every completed branch. Local timestamps and hashes support internal reconstruction and content integrity but cannot independently establish prospective registration. We therefore avoid preregistration terminology, report branch-specific amendments transparently, and restrict post-gate source-wide analyses to exploratory interpretation.

### Short abstract-compatible form

> We used an internally specified, but not publicly preregistered, sequential computational framework with explicit stopping rules and a conservative claim ceiling.

## 7. Required remediation before final repository/publication

1. Freeze the final protocol, SAP, configurations and code in version control; create a signed release or public archival deposit with an immutable timestamp and DOI. This future deposit documents the final package but must not be backdated or described as pre-analysis preregistration.
2. Re-run every retained final analysis from the frozen release in a clean environment. Write the source/configuration SHA-256 values into the run manifest **at process start**, and write completion/output hashes at process end.
3. Preserve exact copies of the executed scripts and configuration beside each run manifest. Do not let a later manifest silently replace an execution-era hash with a current-source hash.
4. Give every future amendment an identifier, full timestamp, affected branch, triggering evidence, old/new rule, outcome-visibility declaration, author and cryptographic hash. Commit it before rerunning the affected branch.
5. Regenerate this chronology after all final reruns. The present TSV is a snapshot at 2026-08-18 08:16:07 CST and must not be assumed to cover later files.

## 8. Bottom line

The defensible statement is not “Route B2 was preregistered.” It is: **the core framework was documented internally before the recorded phenotype-anchor analysis, subsequent changes were partly branch-specific and logged, and the project provides a reconstructed local audit trail; however, no public immutable preregistration exists and exact execution-time source snapshots are incomplete.** This wording preserves the scientific value of the conservative zero-result design without overstating its registration status.
