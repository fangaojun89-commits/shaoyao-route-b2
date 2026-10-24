#!/usr/bin/env Rscript
# Metadata/readback closure only. Reuses the already completed observed RDS;
# never imports a driver/fit function or invokes an unsigned method.
# Usage: Rscript --vanilla THIS.R EXISTING_RUN
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
library(data.table)
library(digest)
library(jsonlite)
require_true <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
identity <- function(path) list(path = normalizePath(path, mustWork = TRUE),
  bytes = unname(file.info(path)$size), sha256 = digest(file = path, algo = "sha256", serialize = FALSE))
verify_identity <- function(expected) {
  actual <- identity(expected$path)
  require_true(actual$sha256 == expected$sha256 && actual$bytes == expected$bytes,
    paste("Preserved input identity", basename(expected$path)))
  actual
}
arguments <- commandArgs(trailingOnly = TRUE)
require_true(length(arguments) == 1L, "Only the existing failed-preparation run is accepted")
run <- normalizePath(arguments[1L], mustWork = TRUE)
worker_path <- file.path(run, "worker.R")
expected <- c(worker.R = "1eb606c6fea739521e786863d5def7d184d3741bb647168a58671af5bfa4a540",
  definition.md = "5af50b1cf337744a6cbb08e331f1457c274f87c789ed75a3ae4a3119d5721215",
  INPUT_MANIFEST.tsv = "591439fffc147ef32154e2309b61e1401dc9e78cd8a8afc722dd4daf40236155",
  RUNTIME_BEFORE.json = "d2ae1f6ea9196c5bbbd0840b474b447aaa110601df44346036755bb85845c9ae",
  "inputs/state.rds" = "afe159fbcc698b1b1f9ac8962fa587d261d45155abafe42a93a74090ae1d12bd",
  "allocations/allocation_0001.rds" = "942a102b24f0602b4974d8004b6a4e76ee9c940696450be6a59c678b7910478d")
before <- lapply(names(expected), function(name) {
  actual <- identity(file.path(run, name))
  require_true(actual$sha256 == expected[name], paste("Exact failed-attempt input", name))
  actual
})
names(before) <- names(expected)
outputs <- file.path(run, c("OBSERVED_COMPATIBILITY.tsv", "OBSERVED.tsv",
  "OBSERVED_DETAILS_mroast.tsv", "OBSERVED_DETAILS_fgsea.tsv", "sessionInfo.txt",
  "FINISHER_RUNTIME.json", "OBSERVED_BRIDGE_TRANSPORT_RECOVERY.json",
  "OBSERVED_METADATA_CLOSURE_MANIFEST.tsv", "PREPARED.json"))
require_true(!any(file.exists(outputs)), "All closure products must be new; do not overwrite a prior attempt")

# Import only the pinned V2's explicitly named readback/validation declarations.
# fit_unsigned, prepare, run_worker, selftest, stochastic_record, and every model
# entrypoint are excluded. No source() or evaluation of top-level driver occurs.
functions <- c("need", "utc", "pin", "verify", "jwrite", "twrite", "backend_snapshot",
  "kernel_contract", "runtime_snapshot", "native_columns", "finite_diagnostic",
  "validate_result", "observed_bridge")
constants <- c("groups", "contrasts", "collections", "key_columns", "numeric_columns",
  "seed_base", "claim_role", "versions")
scope <- new.env(parent = globalenv())
found <- character()
for (expression in as.list(parse(worker_path, keep.source = FALSE))) {
  if (!is.call(expression) || !identical(expression[[1L]], as.name("<-")) || !is.symbol(expression[[2L]])) next
  name <- as.character(expression[[2L]])
  is_function <- is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function"))
  if ((name %in% functions && is_function) || (name %in% constants && !is_function)) {
    require_true(!name %in% found, "Unique declaration import")
    eval(expression, scope); found <- c(found, name)
  }
}
require_true(setequal(found, c(functions, constants)), "Complete allowlisted declarations only")
rownames(scope$contrasts) <- scope$groups
require_true(as.character(getRversion()) == "4.6.1", "Same existing R runtime")
for (package in names(scope$versions)) require_true(as.character(packageVersion(package)) ==
  scope$versions[package], "Same pinned package version")
require_true(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")), "Same RNG kind")
setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
require_true(is.null(BiocParallel::bpRNGseed(BiocParallel::bpparam())) &&
  BiocParallel::bpnworkers(BiocParallel::bpparam()) == 1L, "Same serial backend")
manifest <- fread(file.path(run, "INPUT_MANIFEST.tsv"), na.strings = "NA")
require_true(all(c("role", "path", "bytes", "sha256") %in% names(manifest)) &&
  !anyDuplicated(manifest$role), "Existing immutable manifest schema")
input_pins <- lapply(seq_len(nrow(manifest)), function(i) verify_identity(as.list(manifest[i])))
upstream_path <- manifest$path[manifest$role == "canonical_input_manifest"]
require_true(length(upstream_path) == 1L, "Unique bound upstream manifest")
upstream <- fread(upstream_path)
upstream_pins <- lapply(seq_len(nrow(upstream)), function(i) verify_identity(as.list(upstream[i])))
runtime_before <- fromJSON(file.path(run, "RUNTIME_BEFORE.json"), simplifyVector = FALSE)
require_true(scope$runtime_snapshot()$kernel_contract_sha256 == runtime_before$kernel_contract_sha256,
  "Kernel/default/backend identity matches the completed fit process")
state_path <- file.path(run, "inputs", "state.rds")
observed_path <- file.path(run, "allocations", "allocation_0001.rds")
state <- readRDS(state_path)
observed <- readRDS(observed_path)
require_true(state$source_binding == expected[["INPUT_MANIFEST.tsv"]] &&
  identical(dim(state$labels), c(1680L, 9L)) && state$spaces$observed == 1L &&
  identical(unname(state$labels[1L, ]), rep(scope$groups, each = 3L)) &&
  all(observed$canonical_specification$orientation == 1) &&
  is.finite(observed$elapsed_seconds) && observed$elapsed_seconds > 0,
  "Bound completed original observed assignment and saved elapsed time")
scope$validate_result(observed, state, 1L)
require_true(nrow(observed$result) == 6496L &&
  nrow(observed$native_payload$mroast) == 3248L && nrow(observed$native_payload$fgsea) == 3248L,
  "Eight complete already-fitted unsigned families")

targets <- setNames(lapply(c("native_mroast", "native_fgsea"), function(role) {
  selected <- manifest$path[manifest$role == role]
  require_true(length(selected) == 1L, "Unique original native bridge target")
  normalizePath(selected, mustWork = TRUE)
}), c("mroast", "fgsea"))

# Lexical override is limited to these two read-only native TSVs. Their fields
# use a normal CSV quote layer: outer quotes and doubled interior quote bytes.
# A standards-compliant CSV read removes that one transport layer. JSON then
# decodes the ordered leadingEdge array. No raw quote substitution is applied.
allowed_targets <- unname(unlist(targets))
scope$fread <- local({
  allowed <- allowed_targets
  function(file, ...) {
    path <- normalizePath(file, mustWork = TRUE)
    require_true(path %in% allowed, "CSV transport parser restricted to the native targets")
    as.data.table(utils::read.delim(path, header = TRUE, sep = "\t", quote = "\"",
      comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
      na.strings = "NA", fill = FALSE))
  }
})
bridge <- scope$observed_bridge(observed, targets)
require_true(nrow(bridge) == 35L && all(bridge$comparisons == 3248L) &&
  all(bridge$exact_NA_mask) && all(bridge$accepted), "Complete original native numeric/text/list compatibility")
leading_edge <- bridge[method == "fgsea" & field == "leadingEdge"]
require_true(nrow(leading_edge) == 1L && leading_edge$accepted &&
  leading_edge$comparisons == 3248L, "Every ordered native leadingEdge matches the saved RDS")

# Finish only the missing metadata/export steps. All existing source/state/RDS
# files are rechecked before and after; no fit function was imported or called.
this_script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
finisher_pin <- identity(this_script)
invisible(lapply(c(before, input_pins, upstream_pins), verify_identity))
scope$twrite(bridge, file.path(run, "OBSERVED_COMPATIBILITY.tsv"))
scope$twrite(observed$result, file.path(run, "OBSERVED.tsv"))
for (method in c("mroast", "fgsea")) scope$twrite(observed$native_payload[[method]],
  file.path(run, paste0("OBSERVED_DETAILS_", method, ".tsv")))
capture.output(sessionInfo(), file = file.path(run, "sessionInfo.txt"))
scope$jwrite(scope$runtime_snapshot(), file.path(run, "FINISHER_RUNTIME.json"))
scope$jwrite(list(status = "OBSERVED_NATIVE_BRIDGE_METADATA_ONLY_RECOVERY",
  time_utc = scope$utc(), original_prepare_terminal_exit = 1L,
  original_failed_attempt_worker = before[["worker.R"]],
  original_failed_attempt_definition = before[["definition.md"]],
  original_input_manifest_unchanged = before[["INPUT_MANIFEST.tsv"]],
  original_state = before[["inputs/state.rds"]], original_observed_checkpoint = before[["allocations/allocation_0001.rds"]],
  finisher = finisher_pin, imported_declarations = functions,
  native_targets = lapply(targets, identity),
  parser = list(transport = "TSV with standard CSV outer quoting and doubled interior quotes",
    transport_decoder = "base utils::read.delim sep=TAB quote=doublequote comment.char=empty fill=FALSE",
    list_decoder = "jsonlite::fromJSON after exactly one CSV quote layer; compare ordered character vectors",
    blanket_raw_quote_replacement = FALSE, frozen_worker_source_modified = FALSE,
    frozen_INPUT_MANIFEST_modified = FALSE),
  ordered_leadingEdge_rows_verified = 3248L, native_field_rows = 35L,
  observed_fit_elapsed_seconds_from_saved_RDS = observed$elapsed_seconds,
  original_fit_process_peak_memory_kb = NULL,
  original_fit_process_peak_status = "NOT_CAPTURED_ORIGINAL_FIT_PROCESS_ALREADY_ENDED",
  fit_or_stochastic_method_calls_by_finisher = 0L,
  scientific_parameters_or_tolerances_changed = FALSE,
  source_and_input_bytes_unchanged = TRUE,
  saved_full_precision_RDS_native_payload_compared_with_original_TSV_targets = TRUE,
  later_export_readback_is_a_separate_transport_measurement = TRUE,
  authenticated_exchangeability = FALSE, candidate_gate_eligible = FALSE,
  actual_finisher_terminal_exit_is_recorded_by_calling_process = TRUE),
  file.path(run, "OBSERVED_BRIDGE_TRANSPORT_RECOVERY.json"))
closure_inputs <- c(before, list(finisher = finisher_pin), lapply(targets, identity),
  list(transport_recovery = identity(file.path(run, "OBSERVED_BRIDGE_TRANSPORT_RECOVERY.json")),
    finisher_runtime = identity(file.path(run, "FINISHER_RUNTIME.json"))))
closure_manifest <- rbindlist(lapply(names(closure_inputs), function(role) {
  value <- closure_inputs[[role]]
  data.table(role = role, path = value$path, bytes = value$bytes, sha256 = value$sha256)
}))
scope$twrite(closure_manifest, file.path(run, "OBSERVED_METADATA_CLOSURE_MANIFEST.tsv"))
find_identity <- function(role) {
  i <- which(manifest$role == role)
  require_true(length(i) == 1L, "Unique existing manifest role")
  as.list(manifest[i, .(path, bytes, sha256)])
}
scope$jwrite(list(status = "PREPARED_FULL_UNSIGNED_NATIVE_PAYLOAD_BRIDGE", completed_utc = scope$utc(),
  derived = identity(state_path), observed = identity(file.path(run, "OBSERVED.tsv")),
  observed_checkpoint = identity(observed_path), worker = find_identity("worker"),
  definition = find_identity("definition"), input_manifest = identity(file.path(run, "INPUT_MANIFEST.tsv")),
  native_bridge = identity(file.path(run, "OBSERVED_COMPATIBILITY.tsv")),
  all_labels = 1680L, companion_labels = c(20L, 20L), family_n = 8L,
  result_rows = 6496L, retained_genes = 14217L,
  program_n = c(Reactome = 1546L, SynGO_sensitivity = 78L),
  observed_seconds = observed$elapsed_seconds, process_peak_memory_kb = NULL,
  process_peak_memory_status = "NOT_CAPTURED_ORIGINAL_FIT_PROCESS_ALREADY_ENDED",
  resource_acceptance_not_inferred_from_missing_peak = TRUE,
  packages = as.list(scope$versions), RNG_kind = RNGkind(), runtime = find_identity("runtime"),
  canonical_payload = TRUE, restored_native_Direction_or_leadingEdge = FALSE,
  inner_MC_uncertainty_eliminated_by_outer_enumeration = FALSE,
  authenticated_exchangeability = FALSE, candidate_gate_eligible = FALSE,
  historical_runtime_authenticated = FALSE, primary_results_modified = FALSE,
  metadata_finisher = finisher_pin,
  metadata_closure_manifest = identity(file.path(run, "OBSERVED_METADATA_CLOSURE_MANIFEST.tsv")),
  metadata_closure_receipt = identity(file.path(run, "OBSERVED_BRIDGE_TRANSPORT_RECOVERY.json")),
  original_prepare_terminal_exit = 1L, observed_fit_repeated = FALSE), file.path(run, "PREPARED.json"))
invisible(lapply(c(before, input_pins, upstream_pins, list(finisher_pin)), verify_identity))
cat("PREPARED_FULL_UNSIGNED_NATIVE_PAYLOAD_BRIDGE; metadata-only closure; no fit calls\n")
