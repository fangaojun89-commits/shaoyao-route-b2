# Explicit metadata-reader adapter for the immutable canonical entry.
# Same six CLI arguments. Only the job admission's four string cell IDs are
# normalized from jsonlite's simplifyVector=FALSE list to a character vector.
# Counts, scores, observed objects, labels, kernels, files and RNG are untouched.
compat_args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(compat_args) == 6L)
compat_source <- normalizePath(compat_args[[1L]], winslash = "/", mustWork = TRUE)
compat_job <- normalizePath(compat_args[[2L]], winslash = "/", mustWork = TRUE)
compat_entry <- file.path(compat_source, "run_GSE128623_reference_shard_20261015.R")
stopifnot(identical(digest::digest(file = compat_entry, algo = "sha256"),
  "21a850894037b5e7bcb487909f15200f007ebecbdf149200b317b5670166c88f"))
compat_environment <- new.env(parent = globalenv())
compat_environment$readRDS <- local({
  expected_job <- compat_job
  expected_cells <- c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
  function(file, ...) {
    value <- base::readRDS(file, ...)
    if (is.character(file) && length(file) == 1L &&
      identical(normalizePath(file, winslash = "/", mustWork = TRUE), expected_job)) {
      cell_ids <- value$context$admission$cell_ids
      stopifnot((is.list(cell_ids) || is.character(cell_ids)), length(cell_ids) == 4L,
        all(vapply(cell_ids, function(id) is.character(id) && length(id) == 1L && !is.na(id), logical(1))),
        identical(unname(unlist(cell_ids, use.names = FALSE)), expected_cells))
      value$context$admission$cell_ids <- expected_cells
    }
    value
  }
})
sys.source(compat_entry, envir = compat_environment)
