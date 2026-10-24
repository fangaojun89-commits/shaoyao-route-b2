# GSE205906 Linux R core: exact sources and accepted installation

`SOURCES.json` records the 70 exact R source archives, original/final acquisition URLs, byte counts, SHA-256 values, DESCRIPTION hashes, declared licences and accepted serial build order. `ACCEPTED_ENVIRONMENT.json` projects the completed 23 October installation check: all 70 installed versions matched, and eight requested root packages loaded. This is an environment check, not a new scientific analysis. Absolute workstation/container paths were omitted; package versions and runtime Conda build identities are retained.

The check reused an existing Linux R 4.6.1 runtime on Ubuntu 22.04.1. It does not establish clean operating-system/compiler provisioning, numerical equivalence of a new fit, restoration of the entire 219-package Windows lock, or restoration of every study branch. Runtime Conda identities are inventory evidence, not an explicit solvable environment file: original channels and complete system provisioning remain to be recovered.

## Restore sequence

1. Provision Linux R 4.6.1 and its native/compiler dependencies in a separate environment. Compare with the retained runtime package inventory. Do not alter an ongoing scientific environment.
2. Obtain each exact source archive using the recorded URL and store it as `filename` in one directory. Some Bioconductor sources were obtained at fixed Git commit URLs. Resolve an unavailable source without changing its expected hash or version. Package archives are not redistributed here; each retains its own licence.
3. Verify all source bytes before installation:

   ```sh
   python -I -B pipeline/acquisition.py check-packages --root /absolute/path/to/sources
   ```

4. Generate the source-verified installation plan with explicit paths:

   ```sh
   python -I -B pipeline/acquisition.py install-plan --root /absolute/path/to/sources --r /absolute/path/to/R --library /absolute/path/to/new-library > install-plan.json
   ```

   This command emits ordered argument arrays and isolation settings; it does not run them or validate the supplied R executable. Before execution, check `R --vanilla --slave -e 'cat(as.character(getRversion()))'` in the isolated environment and confirm 4.6.1. Create the new library, apply the recorded environment overrides, unset an inherited `R_HOME`, and execute each argument array directly, stopping on a nonzero exit. Keep per-package commands, exit codes and logs. The compiler and R must remain on a controlled PATH; use temporary storage on the same adequately provisioned volume.

5. After installation, check all installed package versions and the origins of loaded namespaces with the supplied acceptance script:

   ```sh
   /absolute/path/to/Rscript --vanilla environment/G205_core/verify_library.R /absolute/path/to/new-library /absolute/path/to/repository/environment/G205_core/SOURCES.json /absolute/path/to/new-acceptance-output
   ```

   The eight roots are limma, edgeR, DESeq2, decoupleR, fgsea, data.table, jsonlite and digest. The script saves its result and session information. Separately record native dependency resolution before any scientific reconstruction. A successful install command alone is not the final load check.

The catalog was checked against all 70 retained archives and their DESCRIPTION hashes, the accepted installed versions, and the original lock. The public source checker has been exercised on those archives; the emitted plan has not been executed as a new Linux restoration. The path-parameterized library checker is a new derivative of the accepted private checker and has been syntax-checked, not executed on a newly built Linux library. Completed scientific calculations were not repeated. These instructions expose a previously private environment dependency; they do not claim a fully portable whole-study pipeline.
