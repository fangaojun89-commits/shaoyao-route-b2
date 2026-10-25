# Read-only original input admission; no build_fixed_layer(), fitting or scoring.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
repo <- spec$repo; out <- spec$output
stopifnot(as.character(getRversion()) == "4.6.1")
for (name in names(spec$source_pins)) {
  path <- file.path(repo, name); expected <- spec$source_pins[[name]]
  stopifnot(file.info(path)$size == expected$bytes,
    identical(digest::digest(file = path, algo = "sha256", serialize = FALSE), expected$sha256))
}
api <- new.env(parent = asNamespace("stats"))
sys.source(file.path(repo, spec$admission), envir = api)
# This binds the preserved original kernel definitions; it invokes no kernel.
kernel <- api$jgw03_v2_bind_kernel(file.path(repo, spec$core))
stopifnot(is.function(kernel$jgw03_build_fixed_layer), is.function(kernel$jgw03_fit_one_allocation))
runtime <- list(R = as.character(getRversion()), platform = R.version$platform,
  digest = as.character(utils::packageVersion("digest")),
  jsonlite = as.character(utils::packageVersion("jsonlite")), locale = Sys.getlocale())
jsonlite::write_json(runtime, file.path(out, "receipts/runtime.json"), auto_unbox = TRUE, pretty = TRUE)
summaries <- list()
for (accession in unlist(spec$datasets, use.names = FALSE)) {
  bound <- api$jgw03_v2_read_bound_inputs(spec$project, accession)
  labels <- readRDS(file.path(repo, "data/jgw03_inputs", accession, "FULL_LABEL_REGISTRY.rds"))
  stopifnot(identical(dim(labels$MC), c(10001L, 12L)),
    identical(dim(labels$high_companion), c(20L, 12L)),
    identical(colnames(labels$MC), as.character(bound$samples$sample_id)),
    identical(colnames(labels$high_companion), as.character(bound$samples$sample_id)),
    identical(as.integer(labels$MC[1L, ]), as.integer(bound$samples$ordered_dose)),
    identical(as.integer(labels$high_companion[1L, ]), as.integer(bound$samples$ordered_dose)))
  # Store scientific inputs with portable relative source records, never closures
  # or path-dependent replacements for historical fixed/fitted objects.
  inputs <- bound
  inputs$source_records$path <- NULL
  inputs$label_registry <- labels
  path <- file.path(out, "prepared", paste0(accession, "_unfitted_inputs.rds"))
  saveRDS(inputs, path, version = 3L, compress = "gzip")
  restored <- readRDS(path)
  stopifnot(identical(serialize(inputs, NULL, version = 3L), serialize(restored, NULL, version = 3L)))
  # Independent TSV-to-matrix and annotation readback, without calling the adapter.
  source_root <- file.path(spec$project, "data/processed/jgw03_three_cell_lines", accession)
  tab <- utils::read.delim(file.path(source_root, "locked_expected_counts.tsv.gz"),
    check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(tab[, -1L, drop = FALSE]); rownames(counts) <- tab$raw_gene_id
  stopifnot(identical(restored$counts, counts),
    identical(restored$annotation, utils::read.delim(file.path(source_root, "gene_annotation.tsv"),
      check.names = FALSE, stringsAsFactors = FALSE)),
    identical(restored$samples, utils::read.delim(file.path(source_root, "locked_samples.tsv"),
      check.names = FALSE, stringsAsFactors = FALSE)), identical(restored$label_registry, labels))
  summaries[[accession]] <- list(count_rows = nrow(inputs$counts), libraries = ncol(inputs$counts),
    fractional_count_cells = sum(inputs$counts != floor(inputs$counts)),
    source_files = nrow(inputs$source_records),
    reactome_terms = length(inputs$resources$collections$Reactome97$sets),
    syngo_proxy_terms = length(inputs$resources$collections$QuickGO_SynGO_sensitivity$sets),
    MC_assignments = nrow(labels$MC), high_companion_assignments = nrow(labels$high_companion),
    exact_saved_object_content_readback = TRUE, independent_count_metadata_readback = TRUE)
}
writeLines(capture.output(sessionInfo()), file.path(out, "receipts/sessionInfo.txt"), useBytes = TRUE)
jsonlite::write_json(list(status = "PASS_ORIGINAL_JGW_INPUT_READER_WITHOUT_SCORING_OR_FITS",
  datasets = summaries, original_kernel_definitions_bound = TRUE,
  scientific_model_fits = 0L, scores_recomputed = 0L, random_allocations_generated = 0L,
  fixed_layers_recreated = FALSE, historical_execution_recreated = FALSE,
  scope = "New portable unfitted input bundles; original counts, annotation, metadata, resources and saved labels"),
  file.path(out, "receipts/PREPARATION_QA.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
