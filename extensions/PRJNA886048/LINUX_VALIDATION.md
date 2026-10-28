# Linux analysis validation from retained quantification inputs

The unmodified analysis entry at commit `2caf377c05e1cc2bf81c76c58ef12bde87fd095a` completed one full run in an isolated Linux environment. It used 63 accepted quantification/reference files for 12 libraries. Import, gene fitting, CAMERA and signed activity each executed once; no historical checkpoint was adopted. A completed-run resume verified 44 data/checkpoint files, fitted no models and left all 59 run files unchanged. [Machine-readable evidence](LINUX_ANALYSIS_VALIDATION.json) records the original comparison, the additional failed criterion and fixed-decision audit separately.

## Numerical agreement and its limits

All 13 tables passed the unchanged `base::all.equal(tolerance=1e-9, check.attributes=FALSE)` comparison per numeric column, with exact character content, dimensions and row correspondence. R uses a relative comparison where appropriate; this is not a maximum absolute error of 1e-9 for every cell. No tolerance or statistical parameter was changed.

An additional maximum-absolute-cell check passed 12 tables. The edgeR primary table exceeded that additional bound:

| Column | Maximum absolute difference | Cells above 1e-9 |
|---|---:|---:|
| logFC | 2.305888036069348e-9 | 19 |
| F | 7.154300085687737e-8 | 3783 |
| PValue | 4.008948006894997e-9 | 374 |
| FDR | 5.76918801620252e-9 | 1724 |

Both gene methods retained zero BH-significant genes; CAMERA retained the same 69 significant membership tests (65 Reactome and four DoRothEA); signed activity retained zero. The significant ID sets and effect directions were unchanged. The gene-method logFC Spearman correlation was 0.9697887148919758 on both platforms. These observations support stability of these fixed decisions in this run, not bitwise numerical identity or biological validation. Six tables were exact after line-ending normalization; none was byte-identical before normalization. The original failed extra-check receipt remains preserved and hash-bound in the evidence record.

## Reconstructing the validated environment

The [base explicit specification](environment/linux_base_explicit.txt) records 104 exact conda-forge artifacts for Linux x86_64, including Python 3.12.14, R 4.6.1 and uv. Its SHA-256 is `d79311167b2e29ee8c185983d97a5f2b0267676eb06e721c7cee8e3e170b8fca`. It contains package download identities, not third-party package binaries or a complete OS image. The [session capture](environment/linux_validation_sessionInfo.txt) records the actual R namespace and BLAS; only the private absolute prefix is replaced by `<isolated-prefix>`.

Use micromamba 2.3.3 with a new, empty work location and a separate cache. The tested bootstrap artifact was `https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.3.3-0.tar.bz2`, 6,515,385 bytes, SHA-256 `e7274528ceb9c20d048a428d6c22d7e02e268f8ffb762c4c365422347c8b8ba2`. Verify that digest before unpacking `bin/micromamba`. Create the new base prefix with `micromamba --no-rc create --prefix /new/work/runtime --file /path/environment/linux_base_explicit.txt --yes`. Use `micromamba --no-rc run --prefix /new/work/runtime` for subsequent commands so the prefix toolchain is activated.

From this extension, use the prefix's uv with `--no-config sync --locked --no-install-project --python /new/work/runtime/bin/python --no-python-downloads` and `UV_PROJECT_ENVIRONMENT=/new/work/python_env`. The existing `uv.lock` is unchanged. Restore `environment/renv.lock` into `/new/work/R_library` with a separate bootstrap library, project and renv cache; do not reuse an unrelated project library. The actual bootstrap was the official archive `https://cloud.r-project.org/src/contrib/Archive/renv/renv_1.2.4.tar.gz`, 1,438,389 bytes, SHA-256 `e63c637dc785d55848d9dbc6c9599378103803efd47c1f3f1f82057c00575e8c`. Its DESCRIPTION Package/Version was checked before installation. The old current-version URL returned 404 and was not silently replaced by a newer renv. Source compilation used the isolated prefix toolchain; no system package installation was performed.

Set `.libPaths` to the bootstrap library plus the new R distribution, load renv 1.2.4, and call `renv::restore(project="/new/work/R_restore_project", lockfile="/path/environment/renv.lock", library="/new/work/R_library", prompt=FALSE)`. All dependency versions are checked by the public preflight. Run the `--native-record-checks` check, fresh-analysis and completed-resume commands in [REPRODUCE.md](REPRODUCE.md) with these new interpreter/library locations. This is the environment configuration exercised here; availability of external package repositories can change.

## Preserved operator error and scope

The analysis completed before the original Python 3.10 supervisor stopped on an unsupported `hashlib.file_digest` call in its postprocessing inventory. That error record is retained. A verified Python 3.12 postcheck then inspected the completed outputs and invoked only completed-run resume; the scientific fit was not repeated to hide the supervisor error.

This validation restored isolated base and project dependencies on an existing Linux host. It did not install a fresh operating system, acquire or process raw FASTQs, repeat native quantification, test new samples, establish remote-host attestation or provide the protected PSS holdout. Frozen result tables and scientific modules remain unchanged. The candidate is not a tagged release.
