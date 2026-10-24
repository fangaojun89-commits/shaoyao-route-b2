#!/usr/bin/env Rscript
# Read-only cached-score bridge. No model, expression preprocessing, or worker source.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
library(digest)
library(jsonlite)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, quote = "", na.strings = "NA")
}
identity <- function(path) {
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
       bytes = unname(file.info(path)$size),
       sha256 = digest(file = path, algo = "sha256", serialize = FALSE))
}

arguments <- commandArgs(trailingOnly = TRUE)
assert(length(arguments) == 3L, "Usage: bridge PROJECT_ROOT OUTPUT_ROOT NEW_RECEIPT")
project_root <- normalizePath(arguments[1], winslash = "/", mustWork = TRUE)
output_root <- normalizePath(arguments[2], winslash = "/", mustWork = TRUE)
receipt_path <- arguments[3]
assert(startsWith(receipt_path, "R:/Shaoyao_storage/cache/"), "Receipt must be on R")
assert(!file.exists(receipt_path), "Refusing existing receipt")
manifest_path <- file.path(output_root, "INPUT_MANIFEST.tsv")
export_path <- file.path(output_root, "SCORE_EXPORT_FOR_VERIFICATION.tsv")
design_path <- file.path(output_root, "SAMPLE_DESIGN.tsv")
manifest <- read_tsv(manifest_path)
families <- c("Reactome", "PROGENy", "DoRothEA")
counts <- c(Reactome = 1656L, PROGENy = 14L, DoRothEA = 101L)
score_rows <- manifest[match(paste0("score_", families), manifest$role), , drop = FALSE]
assert(!anyNA(score_rows$role), "Missing score roles")
score_paths <- file.path(project_root, score_rows$relative_path)
assert(all(startsWith(normalizePath(score_paths, winslash = "/", mustWork = TRUE),
                      paste0(project_root, "/data/processed/GSE16710/"))), "Unexpected score location")
paths <- c(score_paths, manifest_path, export_path, design_path)
before <- lapply(paths, identity)
assert(identical(vapply(before[1:3], `[[`, "", "sha256"), score_rows$sha256),
       "Score manifest hash mismatch")
assert(all(vapply(before[1:3], `[[`, 0, "bytes") == score_rows$bytes),
       "Score manifest byte mismatch")
export <- read_tsv(export_path)
design <- read_tsv(design_path)
assert(nrow(design) == 26L && !anyDuplicated(design$sample_id), "Invalid sample IDs")
assert(identical(names(export), c("collection", "feature_id", design$sample_id)),
       "Export sample columns differ")
assert(nrow(export) == 1771L && !anyDuplicated(paste(export$collection, export$feature_id)),
       "Export feature IDs differ")
comparisons <- lapply(seq_along(families), function(index) {
  family <- families[index]
  cached <- readRDS(score_paths[index])
  observed <- export[export$collection == family, , drop = FALSE]
  assert(is.matrix(cached) && is.numeric(cached) &&
           identical(dim(cached), c(counts[[family]], 26L)), "Wrong cached score dimensions")
  assert(identical(rownames(cached), observed$feature_id) &&
           identical(colnames(cached), design$sample_id), "Cached score order differs")
  exported <- as.matrix(observed[, design$sample_id, drop = FALSE])
  assert(is.numeric(exported) && all(is.finite(cached)) && all(is.finite(exported)),
         "Nonfinite or nonnumeric score bridge")
  differences <- abs(cached - exported)
  assert(all(differences <= 1e-12 * (1 + abs(cached))), "Score values exceed export precision")
  list(collection = family, features = nrow(cached), samples = ncol(cached),
       compared_elements = length(cached), maximum_absolute_difference = max(differences),
       tolerance = "1e-12 * (1 + abs(cached))", value_identity = "MATCH_WITHIN_EXPORT_PRECISION")
})
after <- lapply(paths, identity)
assert(identical(before, after), "Bridge input changed during read-only comparison")
receipt <- list(status = "PASS_READ_ONLY_RDS_EXPORT_BINDING", inputs = before,
                inputs_unchanged = identical(before, after), comparisons = comparisons,
                total_elements = sum(vapply(comparisons, `[[`, 0L, "compared_elements")),
                R_version = R.version.string,
                packages = list(digest = as.character(packageVersion("digest")),
                                jsonlite = as.character(packageVersion("jsonlite"))),
                model_fitted = FALSE, original_worker_sourced = FALSE,
                original_CEL_or_expression_matrix_read = FALSE)
write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE,
           digits = 17, null = "null", na = "null")
cat(receipt$status, receipt$total_elements, "elements; no fitting\n")
