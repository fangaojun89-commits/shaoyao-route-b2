# Canonical native record and numerical-content check

From the extension directory, with Python 3.11 or later:

```sh
python -I -B pipeline/native/verify.py --run SRR21775856 --native /path/to/native_bundle --tx2gene /path/to/kallisto_tx2gene.tsv
```

The bundle contains `receipt.json`, `fastp.json`, `fastp.html`, `quant/run_info.json` and `quant/abundance.tsv`. The transcript reference is the exact Ensembl107/GRCm39 reference produced by the sequence-preparation entry, SHA256 `b4bd74bd2e0056d950940686e31bcec66eacdcea2873d4a70cb1df431ee0fda5`. The source ENA ledger is authenticated from the package. Source reads remain outside the package and are not read by this command.

This read-only entry authenticates its three code files against the candidate manifest, compares the recorded original input MD5/size identities to the fixed ENA ledger, checks the fixed tool/index identities, five successful sequential stage records, exact fastp/kallisto scientific arguments, paired paths, output hashes and consistency with actual QC. Orientation requires the original 250,000-pair three-direction comparison with both stranded pseudoalignment fractions strictly between 0.45 and 0.55 of the unstranded result. Numerical checks require 146,190 unique targets matching the frozen reference, finite nonnegative abundance, paired-read accounting, estimated-count sum within max(1, pseudoaligned pairs × 1e-6), and TPM sum within 5 of one million. No clinical threshold or candidate-selection rule is introduced.

Paths, PIDs and valid sequential timezone-aware timestamps may vary. They are checked for internal agreement rather than compared to historical JSON bytes. A retained record is not independently signed execution evidence: `remote_execution_attested` remains false, and the command cannot prove that reported raw-file digests were measured on a trusted host. Use a separately inspected execution environment/runner to produce the records. No native processing, transfer, new acceptance receipt, checkpoint adoption, or model fitting occurs here.

Nine existing libraries (48–53 and 55–57) passed this canonical check. Libraries 46/47 used independent orientation evidence, and 54 has a documented continuation; the canonical entry stops on those historical variants and does not invalidate their existing historical acceptances. Sixteen invented invalid record/content combinations were rejected, including a coordinated change to the filter parameter and corresponding QC command/hash. A synthetic relocation with changed JSON hashes and times passed while preserving numerical content. This is bounded regression evidence, not a fresh raw-read execution on a second host.

The analysis entry retains historical exact-file checking by default. Its explicit --native-record-checks option binds all 63 actual input/resource files and checker code to subsequent checkpoints. Historical checkpoint adoption is rejected. Noncanonical 46/47/54 records must match all five historical file identities exactly. Inputs are checked again before final completion publication.

Add --native-record-checks to the documented pipeline/run.py --config ... --check-only command with the locked R/Python environments. This real read-only check passed on all 12 retained libraries, and historical checkpoint adoption was rejected. Full run/resume in this mode has not been executed. The existing frozen-table comparison remains required. Portable raw processing and independent fresh execution remain outstanding.
