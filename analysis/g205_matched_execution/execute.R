# Explicit runner. No model call unless the Python parent supplied execute mode.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
stopifnot(identical(spec$schema, "G205_MATCHED_PORTABLE_RUN/v1"), spec$mode %in% c("prepare", "execute"))
out <- spec$output; repo <- spec$repo; phase <- "AUTHENTICATE"
read_raw <- function(path) readBin(path, "raw", n = file.info(path)$size)
pin <- function(raw) list(bytes = length(raw), sha256 = digest::digest(raw, "sha256", serialize = FALSE))
held <- function(path, expected) {
  raw <- read_raw(path); actual <- pin(raw)
  stopifnot(actual$bytes == expected$bytes, identical(actual$sha256, expected$sha256)); raw
}
source_raw <- lapply(names(spec$source_pins), function(p) held(file.path(repo, p), spec$source_pins[[p]]))
names(source_raw) <- names(spec$source_pins)
inputs_raw <- lapply(names(spec$input_pins), function(p) held(p, spec$input_pins[[p]]))
names(inputs_raw) <- names(spec$input_pins)
definitions <- function(raw) {
  env <- new.env(parent = baseenv())
  for (node in as.list(parse(text = rawToChar(raw), keep.source = FALSE))) {
    stopifnot(is.call(node), identical(node[[1L]], as.name("<-")), length(node) == 3L,
      is.symbol(node[[2L]]), is.call(node[[3L]]), identical(node[[3L]][[1L]], as.name("function")))
    eval(node, env)
  }
  env
}
json <- function(value, path) jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
write_once_rds <- function(value, path) {
  stopifnot(!file.exists(path), !file.exists(paste0(path, ".partial")))
  saveRDS(value, paste0(path, ".partial"), compress = "gzip", version = 3L)
  stopifnot(file.rename(paste0(path, ".partial"), path))
}
status <- tryCatch(withCallingHandlers({
  core <- definitions(source_raw[[spec$contract$auxiliary_roles$core]])
  adapter <- definitions(source_raw[["analysis/g205_matched_execution/adapter.R"]])
  sources <- lapply(spec$contract$engine_roles, function(p) source_raw[[p]])
  pure <- new.env(parent = baseenv())
  core$g205m_import(source_raw[[spec$contract$auxiliary_roles$resource_adapter]], "gse205906_resource_label_adapter_v1", pure)
  core$g205m_import(source_raw[[spec$contract$auxiliary_roles$writer]],
    c("gse205906_ens101_worker_need_v1", "gse205906_ens101_worker_tsv_v1"), pure)
  core$g205m_import(source_raw[[spec$contract$auxiliary_roles$caller]],
    c("g205run_need", "g205run_source_universe", "g205run_domains", "g205run_tsv_view"), pure)
  phase <- "ORIGINAL_PREPARATION"
  # Decode the already hash-authenticated bytes; never reread an unchecked RDS.
  converted <- adapter$g205p_decode_gzip_rds(inputs_raw[[spec$prepared]])
  resources <- pure$gse205906_resource_label_adapter_v1(adapter$g205p_decode_gzip_rds(inputs_raw[[spec$resource_paths$frozen_resource_RDS]]))
  con <- textConnection(rawToChar(inputs_raw[[spec$resource_paths$orthology]]), local = TRUE)
  orthology <- tryCatch(utils::read.delim(con, sep = "\t", quote = '"', colClasses = "character",
    na.strings = "NA", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE), finally = close(con))
  binding <- list(schema = "G205_MATCHED_PUBLIC_INPUT_BINDING/v1", prepared_input = spec$prepared_identity,
    input_pins = spec$input_pins, source_pins = spec$source_pins, historical_native_custody_recreated = FALSE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_before <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  rng_before <- RNGkind()
  before <- serialize(list(converted, resources, orthology), NULL, version = 3L)
  prepared <- adapter$g205p_prepare(core, sources, converted, resources, orthology, binding)
  runtime <- adapter$g205p_runtime_check(prepared$engine$cfg)
  data.table::setDTthreads(1L)
  stopifnot(identical(before, serialize(list(converted, resources, orthology), NULL, version = 3L)),
    identical(rng_before, RNGkind()), identical(had_seed, exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)),
    !had_seed || identical(seed_before, get(".Random.seed", envir = .GlobalEnv)))
  json(runtime, file.path(out, "receipts/runtime.json"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "receipts/sessionInfo.txt"))
  qa <- list(status = "PREPARED_ORIGINAL_FUNCTIONS_NO_FIT", scientific_model_runs = 0L,
    prepared_scientific_hash = prepared$packet$sha256, preparation_RNG_unchanged = TRUE,
    source_input_genes = nrow(prepared$packet$counts), libraries = ncol(prepared$packet$counts),
    historical_native_custody_recreated = FALSE, independent_numeric_acceptance = FALSE,
    numerical_equivalence_to_historical_results_established = FALSE)
  if (spec$mode == "execute") {
    phase <- "ORIGINAL_MATCHED_FIT_ONCE"
    callback_n <- 0L
    candidate <- function(value, packet, checks) {
      stopifnot(callback_n == 0L); callback_n <<- callback_n + 1L
      write_once_rds(list(schema = "G205_MATCHED_PORTABLE_CANDIDATE/v1", result = value,
        prepared_packet = packet, scientific_acceptance = FALSE), file.path(out, "results/candidate.rds"))
    }
    result <- core$g205m_fit(prepared, on_candidate = candidate)
    stopifnot(callback_n == 1L, all(unlist(result$integrity_checks)))
    # Save completed native computation before any derived export; never retry fits.
    write_once_rds(result, file.path(out, "results/native_result.rds"))
    phase <- "EXPORT_AND_READBACK"
    domains <- pure$g205run_domains(prepared$packet, result, resources)
    tables <- c(result$tables, list(source_domains = domains, slots = result$slots,
      mapping_coverage = result$mapping_coverage, native_not_returned_programs = result$native_not_returned_programs))
    for (role in names(tables)) {
      bytes <- pure$gse205906_ens101_worker_tsv_v1(pure$g205run_tsv_view(tables[[role]]))
      path <- file.path(out, "results", paste0(role, ".tsv"))
      stopifnot(!file.exists(path)); writeBin(bytes, path)
      stopifnot(identical(bytes, read_raw(path)))
    }
    restored <- readRDS(file.path(out, "results/native_result.rds"))
    stopifnot(identical(serialize(restored$tables, NULL, version = 3L), serialize(result$tables, NULL, version = 3L)),
      identical(serialize(restored$slots, NULL, version = 3L), serialize(result$slots, NULL, version = 3L)))
    qa$status <- "NEW_NATIVE_RESULT_SAVED_NOT_NUMERICALLY_ACCEPTED"
    qa$scientific_model_runs <- 1L; qa$slots <- nrow(result$slots)
    qa$native_fit_invocations <- 1L; qa$source_domain_rows <- nrow(domains)
    qa$table_readback_passed <- TRUE
  }
  phase <- "FINAL_INPUT_IDENTITY"
  for (p in names(spec$input_pins)) invisible(held(p, spec$input_pins[[p]]))
  for (p in names(spec$source_pins)) invisible(held(file.path(repo, p), spec$source_pins[[p]]))
  json(qa, file.path(out, "receipts/EXECUTION_QA.json"))
  0L
}, warning = function(w) invokeRestart("muffleWarning"), message = function(m) invokeRestart("muffleMessage")),
error = function(e) {
  json(list(status = "FAILED_NEW_EXECUTION", phase = phase, error_class = class(e),
    raw_error_text_saved = FALSE, automatic_fit_retry = FALSE), file.path(out, "receipts/FAILURE.json"))
  1L
})
quit(status = status)
