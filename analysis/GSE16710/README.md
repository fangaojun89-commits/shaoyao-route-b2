# Experimental-spasticity anchor: original analysis source

This directory publishes seven original R scripts, unchanged from the sealed
23 October source archive. It covers the original GSE16710 preprocessing,
resource construction, primary analysis and its original sensitivities,
matched-set reference, precision reporting and output checks. The later
original-annotation/all-27 crossed cell and OLS baseline companions are separate
sources and are not included here. This is one branch of the larger study.

## What is included

| Path | Purpose |
|---|---|
| `code/R/preprocess_qc_GSE16710.R` | Raw-array QC, RMA and current/original probe annotation |
| `code/R/freeze_gene_sets.R` | Construct the frozen resource object from separately obtained resources |
| `code/R/analyze_GSE16710_spasticity_anchor.R` | Primary 26-array injury-by-time model; original annotation, all-27, leave-one-out and conditional-reference comparisons; S1 |
| `code/R/analyze_GSE16710_QuickGO_SynGO_sensitivity.R` | Separately scoped QuickGO assignedBy-SynGO sensitivity |
| `code/R/analyze_GSE16710_matched_random_gene_set_null.R` | Matched random gene-set reference |
| `code/R/summarize_GSE16710_precision.R` | Original effect-precision summary, sourced by the main analyzer |
| `code/R/validate_GSE16710_outputs.R` | Original checks on generated outputs |
| `metadata/samples.tsv` | Exactly the 27 GSE16710 rows selected from the original project ledger; all retained cells unchanged |
| `data/processed/GSE16710/qc_inclusion_decision.tsv` | Original retained technical adjudication, including the primary exclusion of GSM418682 |
| `expected/S1_GATE.tsv` | Unchanged retained result: 1,771 tested programs, zero admitted |
| `EXTERNAL_INPUTS.json` | Exact identities of 11 separately stored raw/processed inputs |
| `metadata/gene_set_manifest.tsv` | Original resource-file inventory; hashes and relative locations, not resource contents |
| `PROVENANCE.json` | Public/source identities and the one sample-row selection |

The archived QC text includes its historical outcome-blind designation. This
publication preserves that record; it does not establish independently verified
pre-outcome chronology. Similarly, a submitted identifier or a historical
biological-unit field does not independently authenticate animal lineage.

## Verify without fitting models

From the repository root:

```sh
python -I -B pipeline/anchor.py verify
python -I -B pipeline/anchor.py inspect-inputs --input-root /path/to/independent-anchor-workspace
```

The first command checks release hashes, source provenance, metadata and retained
gate counts. The second only reads the 11 named external files and checks their
identities; exit 3 means at least one is missing or different. Neither downloads,
installs, evaluates an R source, fits a model, modifies an input or creates output.
An input-identity PASS is not a whole-environment or whole-analysis PASS.

## Data and environment needed for reconstruction

Obtain the GSE16710 raw archive and series matrix, and GPL1355 annotation from
their [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16710) and
[platform](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1355) records.
Match the exact frozen files against `source_file_manifest.tsv` and
`EXTERNAL_INPUTS.json`; current upstream versions may differ. They are not
embedded in this repository.

The resource builder needs the original NCBI gene-group snapshot, Reactome 97
GMT, QuickGO assignedBy-SynGO snapshot, Matrisome masterlists and the specified
DoRothEA/PROGENy package data. Paths, identities and intended versions are in
`config/gene_sets.yaml`, `metadata/gene_set_manifest.tsv` and the builder itself.
Keep resource terms separate from the project's code licence. Replacing a missing
snapshot with a current network changes the scientific input and cannot be called
an exact reconstruction. The matching processed frozen RDS is an alternative
model input, not proof that the raw resource build was repeated.

The retained `renv.lock` specifies R 4.6.1 and 219 packages. It is an original
project lock, not a complete automatic installation recipe for this branch.
The preprocessing record additionally specifies rat2302cdf 2.18.0; the resource
configuration specifies DoRothEA 1.23.0 and PROGENy 1.34.0. Preserve exact package
versions and their source licences. `preprocessing_package_versions.tsv` records
the observed preprocessing packages. A fresh operating-system/toolchain restore
has not been demonstrated by this source release.

## Dependency order and future execution

Use a **fresh separate workspace** with this directory's structure; never point
the original writers at an existing study result tree. The R entries infer their
workspace root as two directories above `code/R`. They write their original
output names; the preprocessing script also extracts CEL files in its raw staging
directory. Keep original downloads immutable and use copies in that fresh tree.

1. Provision the exact R/native environment and lawfully obtained frozen inputs.
   Record platform, session information, actual commands and input identities.
2. Place the three raw GSE16710 files under `data/raw/GSE16710`, preserve the
   27-row metadata and the original QC decision. The QC decision is a retained
   curation input; it must not be silently recalculated or substituted.
3. In that fresh workspace, preprocessing is `Rscript --vanilla
   code/R/preprocess_qc_GSE16710.R`. The separate resource build is `Rscript
   --vanilla code/R/freeze_gene_sets.R` after all its frozen inputs are present.
4. Inspect raw/processed identities with the read-only Python command above.
   Run the main analyzer, then the QuickGO-SynGO sensitivity and matched-set
   reference using their same-named Rscript entries. The precision script is
   already sourced by the main analyzer. The original validator runs after the
   main and QuickGO-SynGO outputs exist.
5. Save actual runtime/session records, output identities and comparison results
   in the fresh workspace. Do not describe a future reconstruction as completed
   merely because these commands or a retained expected table are available.

The main analyzer sets seed 20260817; other random operations retain their
source-specific seeds and RNG behavior. Do not replace them with an inferred
global seed. Source code records the exact models, contrasts and statistical
families. This publication performed no model fitting or resampling.

## Interpretation and remaining scope

GSE16710 is an experimental spinal-injury spasticity anchor, not a stroke cohort.
The primary S1 result was empty. Downstream candidate tests remain unevaluated;
this does not demonstrate absence of treatment efficacy. The original validator
checks its declared outputs, not independent phenotype validity or every later
sensitivity. The published expected gate is an aggregate research result, not
third-party gene-set membership.

Other study-branch sources, all later anchor companions and complete cold
restoration remain to be delivered/established. Original code is MIT; project
curation and aggregate results are CC BY 4.0. Dependency metadata and all
externally acquired resources retain their respective provenance and terms.
