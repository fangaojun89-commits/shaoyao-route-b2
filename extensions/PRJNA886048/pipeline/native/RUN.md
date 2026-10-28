# Portable native entry: implemented, execution unverified

`run.py` plans or processes one new fixed-cohort library on Linux. It accepts only SRR21775846–SRR21775857 and the frozen ENA ledger. It contains no downloader, index builder, deletion, restart, or statistical model fitting. Do not use it to reprocess a completed library or bypass a partial historical stage with a different output directory. Existing accepted bundles should be consumed directly.

The scientific body in `scientific_core.py` is copied from the actually executed R3 native entry. `SCIENTIFIC_CORE_BINDING.json` records the source hash and normalized scientific-section hash. Only the function boundary and indentation change. The wrapper checks every generated command against the frozen five-stage plan before launching it.

## Paths and commands

Use Python 3.11 or later, without `-O`. Specify resolved absolute Linux paths. `--input-dir` contains both complete original files named `SRR21775856_1.fastq.gz` and `SRR21775856_2.fastq.gz` for this example. Prefix/suffix fragments are not supported. Input files and the package are read-only. The fastp 0.20.0 executable, kallisto 0.46.1 executable, Ensembl107/GRCm39/k31 index and transcript reference must match the exact identities enforced by the canonical record checker. Tools and index must already exist; downloading, installing and constructing them are separate inspected steps.

`--work-dir` and `--archive-dir` must already exist. Choose an inspected persistent filesystem for the archive directory. Work may use RAM only if its eventual retention and host shutdown are managed separately. This entry never shuts down a host or deletes temporary/intermediate files. The output is `WORK/linux_full_SRR21775856_01`; the archive is `ARCHIVES/SRR21775856_results_portable.tar.gz`. Keep inputs, package, output and archive directories separate.

From the extension directory:

```sh
python -I -B pipeline/native/run.py --run SRR21775856 \
  --input-dir /data/raw --fastp /tools/fastp --kallisto /tools/kallisto \
  --index /reference/Ensembl107_GRCm39_k31.idx \
  --tx2gene /reference/kallisto_tx2gene.tsv \
  --work-dir /work --archive-dir /persistent/results --plan
```

The plan reads the package ENA ledger only; it does not read raw files, hash tools or index, create directories, or launch native processes. It prints exact input sizes/MD5s, five tokenized commands and the historical resource budget. Replace `--plan` with `--preflight` for a read-only Linux check of paths, raw-file sizes, frozen tool/index/reference hashes, available disk and finite cgroup-v1/v2 memory headroom. Preflight does not hash raw reads or launch fastp/kallisto.

Only after inspecting those results, replace the mode with `--execute` to actually process that new library. The wrapper authenticates the package, refuses existing or partial output/archive paths, acquires an exclusive work-directory ownership lock, and hashes both complete raw files against the ledger before scientific processing. The commands retain 4 orientation threads and 8 fastp/full-quantification threads; cgroup quota/CPU affinity must admit at least eight cores. The historical memory budget requires 1.2 times compressed input bytes plus 18 GiB headroom; output space requires the same 1.2 envelope plus 1 GiB. A finite container memory limit is required; no host `MemAvailable` value is used. Resource budgets are conservative inherited admission rules, not measured portable-host peak usage.

## Results and failure handling

The unchanged science performs the first-250,000-pair three-direction audit, requires both stranded/unstranded pseudoalignment ratios strictly inside 0.45–0.55, runs the frozen fastp filters, then unstranded full kallisto quantification. Recorded source MD5s, stage PIDs/times/arguments, four output hashes and read accounting go to `receipt.json`. Output is checked by the canonical record/content gate, including the exact 146,190-target universe and finite nonnegative estimates. Successful local checking is not trusted-host attestation or acceptance by a separate data custodian.

A successful run creates a five-member archive (`receipt.json`, `fastp.json`, `fastp.html`, `quant/run_info.json`, `quant/abundance.tsv`), rereads every archived member and reports its size/SHA256. Uploading/downloading that archive and accepting it on R are separate operations. Native outputs can then be supplied to the analysis entry's explicit `--native-record-checks` mode; historical checkpoint adoption remains forbidden there.

Any existing output or failed stage stops for inspection. There is no automatic retry/resume and no cleanup. If interrupted while the exact child is still live, ownership stays in place; inspect its actual PID/command and files before considering recovery. An observation timeout does not prove completion or failure. Locks do not coordinate unrelated work directories or attest that a library was never processed elsewhere; the operator must first inspect and reuse any existing work.

## Verification boundary

Seventeen bounded local planning checks passed: generated command identity for nine canonical retained libraries after path relocation, normalized scientific-body identity, unknown-library/traversal/path-collision rejection, existing-output preservation, wrong-tool rejection and a real isolated `--plan` invocation without raw files or tools. Two additional synthetic archive-helper tests passed: five-member readback and preservation of a destination created just before archive publication. These 19 checks launched no native processes and read no original FASTQ. Linux cgroup preflight, the actual portable native execution and its real five-member archive readback have **not** been demonstrated. Synthetic test success must not be described as raw-read reproduction. The current 12 accepted libraries remain the source of the reported study results.
