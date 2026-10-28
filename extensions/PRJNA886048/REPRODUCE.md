# Fixed-cohort execution and evidence boundaries

The entry binds the exact 12 independently accepted quantification bundles, all reference inputs, the sample map and unchanged scientific modules. It accepts only four path settings. Statistical settings stay in the original frozen modules. The raw FASTQ directory is not accessed. Inputs and the released package are read-only; all work goes to a new, non-overlapping directory.

## Environment

Python 3.12.14, with the exact comparison dependency closure in `uv.lock`:

```sh
uv sync --locked --no-install-project
```

R 4.6.1 uses `environment/renv.lock` (44 non-base packages; the contract also records eight base packages). The lock was captured from the actual installed package metadata, including the recorded Bioconductor/R-universe sources; it is not a reconstructed historical session. From this directory, in R:

```r
renv::restore(project=getwd(), lockfile="environment/renv.lock",
              library="renv/library", prompt=FALSE)
```

Bootstrap renv 1.2.4 separately if absent. Point `--r-library` to the restored library. The runner verifies R and every package version before any analysis. The current full namespace session and the historical partial gene-export session are both retained and labelled. R source compilation may require Rtools on Windows or a matching compiler/BLAS toolchain on Linux. The numeric comparison tolerance is 1e-9; graphics bytes are not treated as platform-independent.

## Inputs and commands

Copy `config/example.json` outside the release and supply the four actual input paths. Relative paths are resolved against that configuration file. No credentials or system connection addresses belong in it. `provenance/EXTERNAL_INPUTS.json` supplies expected sizes and hashes. The native root contains `linux_full_SRR21775846_01` through `linux_full_SRR21775857_01`, each with its receipt and four accepted outputs. No raw file download or hash pass is started by the entry.

```sh
uv run --locked python -I -B pipeline/run.py --config /path/inputs.json --check-only --rscript Rscript --r-library /path/R-library
uv run --locked python -I -B pipeline/run.py --config /path/inputs.json --run --output /path/new-run --rscript Rscript --r-library /path/R-library
uv run --locked python -I -B pipeline/run.py --config /path/inputs.json --resume --output /path/new-run --rscript Rscript --r-library /path/R-library
```

For the original data custodian only, existing completed checkpoints can be reused with `--adopt-completed-checkpoints /path/accepted-cohort-output` on `--check-only` or `--run`. All four objects must match `ACCEPTED_CHECKPOINTS.json`, which is tied to the same exact input manifest. These private model files are not distributed. Adoption exports existing results; it does not fit genes, run CAMERA or calculate new ULM scores. This option cannot accept arbitrary models or results from a different input set.

The untouched modules implement the original count import, model, membership and signed-score analyses. The public wrapper implements path resolution, identity, dependency checks, checkpoint state and export ordering. The MDS coordinate calculation follows the already reviewed limma namespace correction. The original comparison script is also unchanged.

## Directory and dependency contract

```text
release/analysis/original/       unchanged scientific modules
release/config/                 path-only example
release/environment/            renv lock and captured versions
release/provenance/             fixed input/checkpoint hashes
external/raw/                   original reads, never accessed by this entry
external/accepted-native/       read-only accepted quantification inputs
new-run/data/intermediate/      verified abundance copies
new-run/data/derived/           fixed sample map and comparison basis
new-run/data/output/            generated tables and comparisons
new-run/checkpoints/            bound RDS objects and stage receipts
new-run/logs/                   actual argv, stdout, stderr, status
new-run/provenance/             environment, binding and frozen output manifest
```

```text
manifest + input hashes + runtime closure -> preflight
accepted abundance + tx2gene -> import checkpoint -> gene model checkpoint
gene model + frozen memberships -> CAMERA checkpoint + signed-score checkpoint
existing/new checkpoints -> exported gene/program/score/QC/MDS tables
scores + complete fixed 19-context basis -> descriptive 266-cell comparison
all tables -> independent semantic readback + frozen numeric/text comparison
verified outputs -> COMPLETED manifest -> hash-verified reuse
```

The output identity binds input, code, lock and sample hashes but not workstation path strings. The seed is 20261025; the frozen statistical methods do not perform stochastic sampling. BLAS/OpenMP thread environment values are set to one for execution. Each R stage writes an in-progress marker before execution and commits a hash-bound checkpoint atomically. A verified checkpoint is reusable; an incomplete checkpoint stops for inspection. A partial descriptive-comparison directory also stops. Fresh runs refuse existing directories. A completed resume verifies all frozen data/checkpoint files and returns without invoking analysis. Original raw inputs, accepted output directories and source model objects are never overwritten.

## Raw reconstruction and resource access

`data/raw_ledger/ENA_quantification_manifest.tsv` records the 24 original read URLs, MD5s and sizes. Historical scientific preparation sources are under `analysis/original/`; their workstation paths are archival, not a portable raw-processing command. The native procedure used fastp 0.20.0, kallisto 0.46.1, Ensembl107 GRCm39 cDNA plus ncRNA, k=31, the original 250000-pair orientation audit, and the recorded index hash. Do not infer that this new accepted-quantification entry performs those upstream steps.

Obtain upstream Ensembl annotation, frozen strict one-to-one orthology, Reactome97 and the mouse PROGENy/DoRothEA resources through their original sources and applicable terms. Membership files and third-party article files are not redistributed. `EXTERNAL_INPUTS.json` prevents using a different resource snapshot silently. An independently verified portable resource acquisition recipe and raw-quantification entry remain release blockers; frozen-file mapping is now integrated below; users without the exact frozen files must not substitute a new version and claim the same analysis.

## Validation limits

R restoration to a new target library and Python restoration to a new virtual environment have now passed on the same Windows host; all 52 R dependencies resolved outside the old project library. Whole-machine restoration and full execution from raw reads have not been demonstrated. Completed-checkpoint export validation is distinct from a fresh model fit. Test receipts state which branch actually ran. No new PSS validation or change in the study's statistical/claim boundaries follows from this packaging work. Keep the candidate separate from the public release until remaining access/reconstruction and environment validation are addressed.

The exact retained native JSON pins include timestamps and command paths. They are suitable for authenticating the accepted bundles, but fresh native outputs need a separate semantic acceptance interface before the entry can claim end-to-end raw-read reproduction. The present validation is bounded to accepted-input and existing-checkpoint reuse.

## Integrated reference reconstruction (r03)

With the exact frozen compressed Ensembl107 inputs and four frozen resource tables available locally:

```sh
python -I -B pipeline/reference/run.py --source /path/frozen_ensembl107 --resources /path/frozen_gene_sets --out /path/new_reference
```

This command parses GTF, prepares the published cDNA/ncRNA universe and maps frozen memberships in sequence. It logs actual commands and stage output. The two mapping-reference files are copied from newly generated outputs because the mapper requires a shared directory. Existing output directories are refused, including partial attempts; inspect any failed attempt rather than overwrite it. This command neither downloads inputs nor builds the kallisto index. Raw reads and expression models are not accessed.

The source directory requires Mus_musculus.GRCm39.107.gtf.gz, Mus_musculus.GRCm39.cdna.all.fa.gz, Mus_musculus.GRCm39.ncrna.fa.gz and gtf_CHECKSUMS.txt, cdna_CHECKSUMS.txt, ncrna_CHECKSUMS.txt. Resources require ncbi_human_mouse_rat_one_to_one.tsv, reactome97_human_long.tsv, progeny_mouse_top500.tsv and dorothea_mouse_AB.tsv. All are hash-pinned by the component entries. Historical acquisition URLs for available records appear in provenance/REFERENCE_SOURCE_LOCATIONS.json; this is not a fresh availability check or a complete resource-acquisition workflow. Third-party resources are not bundled.

Use annotation/genes.tsv, sequences/kallisto_tx2gene.tsv and memberships outputs as reconstructed reference artifacts. The fixed-cohort configuration can point tx2gene to sequences/kallisto_tx2gene.tsv, mapped_edges to memberships/mapped_edges.tsv, and reference_terms to memberships/reference_coverage.tsv. Their exact required pins remain unchanged. The accepted-native root is still required separately. Fresh native quantification acceptance and automatic resource acquisition remain open.
