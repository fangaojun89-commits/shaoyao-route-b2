# Rebuild the frozen gene resources

This command runs the unchanged original `freeze_gene_sets.R` in a new workspace. It uses eight authenticated local source files and the recorded package versions, then requires all thirteen generated resource files—including the complete RDS bundle—to match the accepted bytes. The original G128 member constructors and G205 resource validator additionally check the rebuilt objects. No counts are fitted and no source file is modified.

## Inputs and environment

Place the eight files under a source root, preserving their paths in `CONTRACT.json`:

- NCBI human, mouse and rat protein-coding gene-group exports dated 18 August 2026, plus their download summary.
- The frozen Reactome release 97 GMT.
- The 17 August 2026 QuickGO annotations assigned by SynGO.
- The human and mouse Matrisome rev. 2014 master-list CSVs.

Their combined size is 11,463,526 bytes. The contract takes their identities from the original `gene_set_manifest.tsv`. Current live exports cannot substitute for these dated snapshots. External resource memberships are not redistributed by this repository; users need lawful access to these exact files and the package data. Source access/version/redistribution closure remains a separate requirement. In particular, this command does not silently fetch a newer upstream release.

Use R 4.6.1 and the ten dependency versions in `CONTRACT.json`. Eight versions are drawn from the retained anchor `renv.lock`; PROGENy 1.34.0 and DoRothEA 1.23.0 are additional requirements confirmed by the original builder and frozen bundle. The anchor lock alone does not include those two data packages. Python uses only its standard library; no Python third-party dependency lock is needed for this entry. No software installation occurs.

The original `environment/run_r.ps1` explicitly selected `English_United States.utf8`; its unchanged source is provided here as `original_run_r.ps1` for provenance. The new wrapper reproduces that locale only in its child processes. This matters: inheriting `C` produced ten matching TSVs but different GMT ordering and RDS encoding. Restoring the original locale produced all thirteen exact files. The tested platform is Windows 11 with the existing pinned R library. The `en_US.UTF-8` mapping on other systems is untested and still must pass the exact-output checks; do not relax those checks to accept a different resource.

## Run

From the repository root, on one line or using your shell's continuation syntax:

```sh
python -I -B pipeline/build_resources.py --source-root /path/to/original-project-root --rscript /path/to/Rscript --r-library /path/to/existing/R-library --output /path/to/new-resource-workspace
```

The output must not exist. A failure retains its workspace and logs. A successful `receipts/RESOURCE_BUILD_COMPLETE.json` authenticates all outputs and the original source inputs again. Existing completed outputs can be reused by matching every identity in that receipt and `CONTRACT.json`; this entry never refits a model or regenerates label assignments.

```text
new-resource-workspace/
  data/raw/gene_sets/          exact staged copies of the eight source inputs
  code/R/freeze_gene_sets.R    unchanged original builder
  data/processed/gene_sets/    thirteen verified derived resource files
  metadata/                   original builder's resource inventory
  results/diagnostics/        original builder's fourteen content checks
  receipts/                   commands, environment, logs and final identities
```

The dependency order is: source/package/locale verification → original resource construction → thirteen exact-output comparisons → G128/G205 original consumer checks → completion receipt. Commands, timestamps, runtime versions, session information, input/output hashes and the source contract are retained. Generated scientific resources stay in the chosen data workspace and are not automatically published.

## Connect to the prepared counts

First prepare G128/G205 counts with the separate [input bridge](../input_bridge/README.md). Then use the rebuilt `data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds` as the resource dependency alongside those six unfitted input objects. Its SHA-256 is `0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51`, the same frozen bundle required by the original G205 worker and retained G128 input record.

G128's analyzer and independent validator must generate the same full resource membership table. G205's original validator checks all ten RDS tables against their matching exports and checks the factor-label adapter. Export rounding is assessed with the original validator's relative tolerance of 1e-10 and zero absolute tolerance; this does not weaken the preceding byte-exact file requirements. Resource membership counts describe the full domain before expression filtering and do not establish tested or eligible programs.

The supplied entry completes resource reconstruction in an existing supported environment. A clean environment installation, complete cross-platform scientific model orchestration and independent phenotype holdout remain unverified. Preserve frozen results, package versions, random seeds and recorded execution parameters when implementing the next analysis stage. The historical model launchers should not be invoked merely to test this resource command.
