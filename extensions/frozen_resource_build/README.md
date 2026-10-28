# Portable summary resource-build contract variant

This is a separately identified local candidate, not public v0.30.0. PORTABLE_CONTRACT_BINDING.json binds the historical and portable summary hashes. Only output_file values in that metadata summary differ. The scientific R builder uses counts/API metadata, not this path column. All other inputs and all thirteen expected resource outputs retain their historical identities.

Run pipeline/build_resources.py with --source-root (root containing the eight input paths), --rscript, --r-library and a new --output directory. It authenticates every input and dependency and refuses an existing output. See the original resource construction documentation for exact environment requirements. No model fitting is performed. This candidate still requires all eight inputs; the NCBI supplement alone does not supply Reactome, QuickGO or Matrisome.

Do not alter the historical release contract to accept arbitrary summaries. This variant accepts only the explicitly bound portable bytes. Other metadata changes must fail. Current API exports are not substitutes for the frozen NCBI snapshot.
