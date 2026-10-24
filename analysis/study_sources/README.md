# Original sources for four study branches

This selection preserves 38 original R files and five configurations byte-for-byte from the sealed study source archive. `PROVENANCE.json` records their identities. It also includes a 99-row projection of the original submitted-unit metadata with all fields, values and ordering retained. These sources extend the public reporting and experimental-spasticity anchor components; they do not constitute a standalone launcher for the whole study.

## Branches and dependency order

| Branch | Included source sequence | Interpretation and remaining dependencies |
|---|---|---|
| GSE128623 | `analyze_gse128623_exact30_linux_model_v2.R` and validator → checked kernel loader → member/comparison guards → execution driver → copy-isolated driver and nonmutating comparison | The old model CLI has a deliberate NO_GO gate. The later driver imports its scientific kernels; running the old CLI is not the accepted execution route. Accepted counts and a tested original-reader preparation entry are now supplied in `data/count_inputs` and `analysis/input_bridge`; native execution context and full model/environment orchestration remain separate requirements. |
| GSE277064 | `analyze_GSE277064.R` → fgsea refresh and program plots; preset-method core → runner → comparison and independent verification | 63 submitted libraries comprise 45 natural-course and 18 trained libraries. Reuse of six natural-course comparators creates overlapping branch membership, not additional independent libraries. The alternative-method runner requires its original input/admission records and versioned environment. |
| JGW03 | `analyze_jgw03_three_cell_lines.R`; separate `reconstruct_pc12_rotations_equivalent_20261013.R` | Three distinct 12-library culture series retain their own controls and dose schemes. The PC12 reconstruction is an equivalent-producer comparison, not authentication of the historical original rotation producer. Full Linux conditional-reference controllers are outside this selection. |
| GSE205906 | 20 source roles in `G205_SOURCE_ROLES.json` → `gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R`, using the full October configuration | The October Rnor6/Ensembl101 pipeline binds older scientific kernels through checked adapters. The September CLI/configuration is a dependency, not a command to run against October counts. Native parent/descriptor context, input objects and environment are external. Later STAR matched-quantification extension sources are not in this selection. |

GSE128623 primary 30-library and predefined 29-library sensitivity analyses share libraries. Raw and adapter-treated analyses are technical specifications of those cohorts. The culture and GSE205906 branches are separate calibration sources, not independent PSS phenotype replications. This release changes no candidate gate or accepted result.

## Files and inputs

```text
analysis/study_sources/
  code/R/                       unchanged scientific sources
  config/                       unchanged historical configurations
  metadata/samples.tsv          four-accession submitted-unit projection
  G205_SOURCE_ROLES.json         exact bytes/SHA-256 for 20 source roles
  SHARED_EXTERNAL_INPUTS.json    identities, not redistributed resources
  PROVENANCE.json                selection and transformation record
```

The three entries in `SHARED_EXTERNAL_INPUTS.json` describe shared GEO metadata, dataset metadata and a frozen resource object. This is not a complete branch input manifest. Source-level input pins and assertions remain intact. Accepted G128/G205 counts and JGW03 input tables are supplied elsewhere in this repository, and `analysis/input_bridge` prepares new unfitted G128/G205 objects. Raw downloads, third-party gene memberships, native custody objects, fitted checkpoints and full environment archives are not supplied here. Input and result identities must be matched before any future reuse; an output directory's existence alone is not a valid checkpoint.

The global `contrasts.yaml` retains historical August planning text, including superseded GSE205906 reference/engine choices. It must not override the October Ensembl101 configuration or serve as a current cross-branch execution specification. Per-branch source and run provenance take precedence. In particular, the nine-library HISAT2 analysis and later matched STAR sensitivity use the same libraries and are not independent replication.

## Read-only verification

From the repository root, with Python 3.12 or later (standard library only):

```sh
python -I -B pipeline/sources.py verify
python -I -B pipeline/anchor.py verify
python -I -B pipeline/ledger.py verify data/ledger
```

The new verifier checks release hashes, all 43 selected source/configuration identities, the 20 GSE205906 roles against their literal pins in the original consumer/worker, and the 99 submitted-unit rows. It never calls R, downloads data, restores an environment, runs a model or edits results. Record the Git tag/commit and command output when using it. A successful check proves these identities and relationships, not scientific reproduction.

The scientific R files were parsed without evaluation during release QA. Preserve the original seeds, RNG, package checks and actual recorded run parameters when reconstructing a branch. A source default (for example a rotation count) must not be substituted for the parameter used by a recorded run. The anchor's `renv.lock` is not a universal lock for all Windows and Linux branches. No branch environment was installed or model rerun for this release.

Two unchanged sources contain historical `R:/Shaoyao_storage/` output-location checks: the GSE277064 preset-method runner and the PC12 equivalent reconstruction. They are recorded in `PROVENANCE.json`. They constrain portability; this release neither removes those checks nor presents those scripts as portable entry points. Any later portability adapter must be separately versioned and verified against the original source.

## Reconstruction status

Public source availability now permits inspection of these branch kernels and their configuration relationships. Full cold reconstruction still requires the remaining acquisition/reference/extension producers, complete branch input manifests, environment closure, and tested orchestration. A future reconstruction should use a fresh workspace with immutable raw inputs, separate intermediates and derived outputs, captured session information and commands, and input/code/environment-bound checkpoints. Do not evaluate these writers in an accepted-result directory. The original study's completed calculations were reused, not rerun, to prepare this release.
