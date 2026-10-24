# Explicit metadata/transport-v2 entry; statistical kernels stay canonical.
# Its additional execution pins supplement, not replace, the original source manifest.
# SOURCE_DIR JOB_RDS JOB_SHA256 FIRST_ID LAST_ID WORKER_ID
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
source_directory <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
job_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stopifnot(grepl("^/root/autodl-tmp/gse128623_label_reference_[A-Za-z0-9_-]+/source$", source_directory),
  grepl("^/root/autodl-tmp/gse128623_label_reference_[A-Za-z0-9_-]+/jobs/[A-Za-z0-9_-]+[.]rds$", job_path),
  grepl("^[0-9a-f]{64}$", args[[3L]]),
  all(grepl("^[1-9][0-9]*$", args[4:5])), grepl("^[A-Za-z0-9_-]+$", args[[6L]]))
first_id <- as.integer(args[[4L]])
last_id <- as.integer(args[[5L]])
stopifnot(!anyNA(c(first_id, last_id)), first_id <= last_id,
  identical(digest::digest(file = job_path, algo = "sha256"), args[[3L]]))
job <- readRDS(job_path)
cell_ids <- job$context$admission$cell_ids
stopifnot((is.list(cell_ids) || is.character(cell_ids)), length(cell_ids) == 4L,
  all(vapply(cell_ids, function(id) is.character(id) && length(id) == 1L && !is.na(id), logical(1))),
  identical(unname(unlist(cell_ids, use.names = FALSE)),
    c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")))
job$context$admission$cell_ids <- unname(unlist(cell_ids, use.names = FALSE))
stopifnot(identical(digest::digest(file = job_path, algo = "sha256"), args[[3L]]),
  is.list(job), all(c("context", "inputs", "observed_bundle", "allocations") %in% names(job)),
  identical(job$context$admission$status, "ALL_FOUR_OBSERVED_AND_RESOURCE_PASS"),
  identical(job$context$admission$cell_ids,
    c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")))
run_root <- normalizePath(job$context$run_root, winslash = "/", mustWork = TRUE)
stopifnot(identical(source_directory, file.path(run_root, "source")),
  startsWith(job_path, paste0(run_root, "/jobs/")),
  identical(digest::digest(job$inputs, algo = "sha256"), job$context$inputs_object_sha256),
  identical(digest::digest(job$observed_bundle, algo = "sha256"), job$context$observed_bundle_object_sha256))
observed_file <- normalizePath(job$context$observed_evidence$file, winslash = "/", mustWork = TRUE)
stopifnot(startsWith(observed_file, paste0(run_root, "/observed/")),
  identical(digest::digest(file = observed_file, algo = "sha256"), job$context$observed_evidence$sha256))
preserved_observed <- readRDS(observed_file)
stopifnot(identical(digest::digest(preserved_observed, algo = "sha256"), job$context$observed_bundle_object_sha256),
  identical(digest::digest(file = observed_file, algo = "sha256"), job$context$observed_evidence$sha256))
rm(preserved_observed)
worker_path <- file.path(source_directory, "GSE128623_reference_shard_core_20261015.R")
worker_sha256 <- digest::digest(file = worker_path, algo = "sha256")
stopifnot(identical(worker_sha256, job$context$worker_core_sha256),
  identical(as.character(getRversion()), "4.6.1"))
source(worker_path, local = TRUE)
source_paths <- setNames(file.path(source_directory, c(
  "GSE128623_gene_camera_label_core_20261015.R", "GSE128623_activity_label_core_20261015.R",
  "GSE128623_unsigned_label_core_20261015.R", "GSE128623_label_bundle_core_v2_20261015.R",
  "GSE128623_reference_archive_core_20261015.R")), names(g128w_source_pins()))
bindings <- g128w_bind_sources(source_paths)
# Preserve the canonical source binding; additional I/O source pins are below.
source_identity <- c(bindings$pins, worker = worker_sha256,
                      entry = job$context$entry_source_sha256)
entry_path <- file.path(source_directory, "run_GSE128623_reference_shard_20261015.R")
stopifnot(identical(digest::digest(file = entry_path, algo = "sha256"), job$context$entry_source_sha256),
  identical(digest::digest(source_identity, algo = "sha256"), job$context$source_manifest_sha256),
  all(vapply(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"),
    function(key) identical(Sys.getenv(key), "1"), logical(1))))

transport_path <- file.path(source_directory, "GSE128623_archive_transport_v2_20261015.R")
transport_expected_sha <- "622f09f17fe23b6cc2c0d933672b431376a09998a47de8cdb6d2e83550acd753"
stopifnot(identical(digest::digest(file = transport_path, algo = "sha256"), transport_expected_sha))
invoked <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(invoked) == 1L)
actual_entry <- normalizePath(sub("^--file=", "", invoked), winslash = "/", mustWork = TRUE)
stopifnot(identical(actual_entry, file.path(source_directory, "run_GSE128623_reference_shard_transport_v2_20261015.R")))
execution_pins <- list(transport_source = list(path = transport_path, sha256 = transport_expected_sha),
  entry_source = list(path = actual_entry, sha256 = digest::digest(file = actual_entry, algo = "sha256")),
  original_source_manifest_sha256 = job$context$source_manifest_sha256)
transport_environment <- new.env(parent = baseenv())
sys.source(transport_path, envir = transport_environment)
bindings$write <- transport_environment$g128t_make_writer(bindings, execution_pins)
data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
stopifnot(is.null(BiocParallel::bpRNGseed(BiocParallel::bpparam())))
job$context$worker_id <- args[[6L]]
worker_output <- file.path(run_root, "workers", job$context$worker_id)
if (!dir.exists(dirname(worker_output))) stopifnot(dir.create(dirname(worker_output)))
stopifnot(!file.exists(worker_output), dir.create(worker_output))
writeLines(c(commandArgs(), paste("job_sha256", args[[3L]]), paste("pid", Sys.getpid())),
  file.path(worker_output, "COMMAND.txt"))
start <- Sys.time()
result <- g128w_run_shard(job$inputs, job$observed_bundle, job$allocations,
  seq.int(first_id, last_id), job$context, bindings, synthetic = FALSE)
stopifnot(identical(digest::digest(file = job_path, algo = "sha256"), args[[3L]]),
  identical(digest::digest(file = worker_path, algo = "sha256"), worker_sha256),
  identical(digest::digest(file = entry_path, algo = "sha256"), job$context$entry_source_sha256),
  identical(digest::digest(file = observed_file, algo = "sha256"), job$context$observed_evidence$sha256),
  identical(digest::digest(job$inputs, algo = "sha256"), job$context$inputs_object_sha256))
stopifnot(identical(digest::digest(file = transport_path, algo = "sha256"), transport_expected_sha),
  identical(digest::digest(file = actual_entry, algo = "sha256"), execution_pins$entry_source$sha256))
saveRDS(list(shard = result, job_sha256 = args[[3L]], sources = source_identity,
  additional_execution_sources = execution_pins,
  started_utc = format(start, tz = "UTC", usetz = TRUE),
  ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  actual_process_exit_confirmed_here = FALSE, full_run_complete = FALSE),
  file.path(worker_output, "SHARD_RECEIPT.rds"), version = 3L, compress = "gzip")
write.table(result$files, file.path(worker_output, "ARCHIVE_MANIFEST.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(worker_output, "SESSION_INFO.txt"))
cat("Requested shard returned successfully; OS exit and full coverage remain separate checks.\n")
