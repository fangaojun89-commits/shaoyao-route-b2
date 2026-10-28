# PRJNA886048 and frozen-resource integration candidate

This local branch extends v0.30.0. It has not been published or tagged as a release. The PRJ extension reports the fixed six-versus-six acute injured-spinal baclofen context; it does not establish protected post-stroke spasticity validation, target engagement or efficacy of Shaoyao-Gancao Decoction.

## Saved results

From the repository root run:

    python -I -B extensions/PRJNA886048/pipeline/verify.py --rscript Rscript

This authenticates complete retained tables and checks their content without fitting models. Read extensions/PRJNA886048/REPRODUCE.md for the analysis entry, environment and checkpoint contracts. A separate canonical native record/content entry is documented in extensions/PRJNA886048/pipeline/native/README.md; it checks fixed parameters and source identities while allowing consistent runtime metadata changes. Its bounded tests are not fresh execution evidence. The historical exact-hash analysis entry must not be represented as accepting arbitrary reruns or granting historical checkpoint reuse for changed inputs.

## Frozen resources

The NCBI and QuickGO historical snapshots are in data/resource_snapshots, each with source attribution and provenance limits. NCBI has an offline page-reconstruction verifier. Raw HTTP headers, pagination tokens, machine logs and the machine-specific historical summary are excluded. Historical analysis code retains original path literals for source identity; use the portable entries for execution.

Acquire Reactome 97 from https://zenodo.org/records/21383214 and the rev.2014 human/mouse workbooks from https://sites.google.com/uic.edu/matrisome/matrisome-annotations/ . The accepted Matrisome CSV representation is the first worksheet parsed with pandas default header and index exported to CSV, LF newlines plus one additional terminal LF. Exact hashes in the contract must match. No Matrisome workbook or CSV is redistributed here.

Use extensions/frozen_resource_build/pipeline/prepare_inputs.py with --contract extensions/frozen_resource_build/analysis/resource_build/CONTRACT.json, --ncbi-snapshot data/resource_snapshots/ncbi_gene_groups_20260817, --quickgo pointing to the retained TSV, --reactome-gmt, --matrisome-human, --matrisome-mouse and a new --output. It validates eight inputs before copying; no latest-version fallback.

Then use extensions/frozen_resource_build/pipeline/build_resources.py with --source-root pointing to that output, --rscript, --r-library and a new --output. Its scientific builder and thirteen expected output identities are unchanged; the separately bound contract only accepts the path-portable NCBI summary. All thirteen files and original G128/G205 consumers passed a real existing-environment build. No whole-study clean-machine restoration is claimed.

Use the resulting four resource TSVs for the PRJ reference/mapping entry. PROGENy1.34.0 and DoRothEA1.23.0 source-package data were independently reacquired and reproduced the mouse inputs exactly; package tarballs are not included. Source-specific rights remain in effect.
