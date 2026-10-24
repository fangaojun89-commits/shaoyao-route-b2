# Recover only the four independently pinned, already fitted ID-2 sentinels.
# No model, download, count, source-file replacement or old-partial deletion.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
stopifnot(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01"))
pin_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stopifnot(startsWith(pin_path, paste0(root, "/")),
  identical(digest::digest(file = pin_path, algo = "sha256"),
    "667e784de22005581683b8a6174c895be940f002ef5ac41acb2e8a8e40d3fdce"))
preparer <- new.env(parent = globalenv())
sys.source(file.path(root, "source/prepare_GSE128623_reference_observed_20261015.R"), envir = preparer)
runtime <- preparer$g128p_runtime()
prepared <- preparer$g128p_load(root)
bound <- preparer$g128p_bind(root)
transport_path <- file.path(root, "source/GSE128623_archive_transport_v2_20261015.R")
transport_sha <- "622f09f17fe23b6cc2c0d933672b431376a09998a47de8cdb6d2e83550acd753"
stopifnot(identical(digest::digest(file = transport_path, algo = "sha256"), transport_sha))
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
stopifnot(length(self) == 1L)
execution_pins <- list(transport_source = list(path = transport_path, sha256 = transport_sha),
  entry_source = list(path = self, sha256 = digest::digest(file = self, algo = "sha256")),
  original_source_manifest_sha256 = prepared$source_manifest_sha256)
transport <- new.env(parent = baseenv())
sys.source(transport_path, envir = transport)
diagnostic <- jsonlite::read_json(pin_path, simplifyVector = FALSE)
stopifnot(diagnostic$run_id == basename(root), diagnostic$diagnostic_actual_exit_code == 0L,
  length(diagnostic$partials) == 4L)
jobs <- jsonlite::read_json(file.path(root, "JOBS_PREPARED.json"), simplifyVector = FALSE)$jobs
output <- file.path(root, "saved_sentinel_recovery_02")
stopifnot(!file.exists(output), dir.create(output))
records <- list()
for (cell in preparer$g128p_cells()) {
  job_record <- jobs[[paste0(cell, "/six_group")]]
  preparer$g128p_check_record(job_record$path, job_record)
  job <- readRDS(job_record$path)
  observed <- bound$binding$pack(job$observed_bundle)
  labels <- stats::setNames(unname(job$allocations[2L, ]), colnames(job$allocations))
  expected <- bound$worker_environment$g128w_archive_identity(job$context, 2L,
    bound$worker_environment$g128w_row_hash(names(labels), labels))
  relative <- paste0("allocations/", cell, "/six_group/allocation_00002.rds.partial")
  candidates <- Filter(function(record) identical(record$relative_path, relative), diagnostic$partials)
  stopifnot(length(candidates) == 1L)
  source_pin <- candidates[[1L]]
  partial <- file.path(root, relative)
  stopifnot(identical(digest::digest(readRDS(partial)$packed, algo = "sha256"),
    source_pin$decoded_stable_packed_sha256))
  final <- sub("[.]partial$", "", partial)
  result <- transport$g128t_recover_partial(partial, final, root, expected,
    observed$keys, observed$families, bound$binding, execution_pins,
    list(bytes = source_pin$bytes, sha256 = source_pin$sha256), job$context$observed_evidence)
  verified <- bound$worker_environment$g128w_read_completed(final, expected,
    observed$keys, observed$families, bound$binding, job$context$observed_evidence)
  stopifnot(identical(result$sha256, verified$sha256), file.exists(partial),
    identical(digest::digest(file = partial, algo = "sha256"), source_pin$sha256))
  records[[cell]] <- list(final = verified, original_partial = source_pin,
    decoded_packed_sha256 = digest::digest(readRDS(final)$packed, algo = "sha256"),
    original_partial_preserved = TRUE, model_calls = 0L)
  jsonlite::write_json(records[[cell]], file.path(output, paste0(cell, ".json")),
    auto_unbox = TRUE, pretty = TRUE, digits = 17)
  rm(job, observed)
  invisible(gc())
}
preparer$g128p_verify_manifest(prepared$manifest)
transport$g128t_check_execution(execution_pins)
jsonlite::write_json(list(status = "FOUR_SAVED_ID2_SENTINELS_RECOVERED_ZERO_MODEL_CALLS",
  model_calls = 0L, records = records, additional_execution_sources = execution_pins,
  independent_diagnostic = preparer$g128p_pin(pin_path), runtime = runtime,
  full_reference_complete = FALSE), file.path(output, "RECOVERY_RECEIPT.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 17)
cat("Recovered four saved native sentinels without fitting; original partials preserved.\n")
