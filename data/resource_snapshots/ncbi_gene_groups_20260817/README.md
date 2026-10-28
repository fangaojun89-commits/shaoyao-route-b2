# Frozen NCBI gene-group snapshot

This local release candidate retains the exact historical human, mouse and rat gene-group TSV exports and all 72 underlying response pages. It contains 70,154 valid gene records and accounts separately for three empty report slots. Retained acquisition dates are 17 August 2026 UTC; the filename suffix is a historical label. This is not a current NCBI release.

Run `python -I -B verify.py` from any working directory. Python 3.9 or later, standard library only; no network, package installation or random computation. The verifier checks every payload hash and reconstructs each combined TSV from the response pages. Manifest hashes establish consistency, not independent authenticity; retain a trusted manifest digest separately.

The portable summary changes only output_file paths, which now refer to sibling TSV files. The original summary hash is retained in PROVENANCE.json. The existing public resource builder requires the original summary hash and cannot yet accept this derivative without an explicitly bound adapter; do not overwrite its frozen contract or claim that this supplement alone closes end-to-end reproduction.

Raw HTTP headers are not distributed. Selected response metadata and original header hashes are retained; pagination token values, request identifiers and machine-specific paths are excluded. These metadata are project-held provenance, not signed NCBI attestations. Original source files remain unchanged. No model results are included.

This candidate has not been published. See SOURCE_NOTICE.md for attribution and source terms.
