# Original deterministic preparation only. No execute(), refit() or compute().
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
repo <- spec$repo; out <- spec$output
data.table::setDTthreads(1L)
pin <- function(raw) list(bytes = length(raw), sha256 = digest::digest(raw, "sha256", serialize = FALSE))
read_raw <- function(path) readBin(path, "raw", n = file.info(path)$size)
capture <- function(relative) {
  raw <- read_raw(file.path(repo, relative))
  expected <- spec$source_pins[[relative]]
  stopifnot(!is.null(expected), identical(pin(raw)$sha256, expected$sha256), length(raw) == expected$bytes)
  raw
}
import <- function(raw) {
  env <- new.env(parent = baseenv())
  env$head <- utils::head; env$tail <- utils::tail; env$setNames <- stats::setNames
  for (node in as.list(parse(text = rawToChar(raw), keep.source = FALSE))) {
    stopifnot(is.call(node), identical(node[[1L]], as.name("<-")), is.symbol(node[[2L]]),
      is.call(node[[3L]]), identical(node[[3L]][[1L]], as.name("function")))
    eval(node, env)
  }
  lockEnvironment(env, bindings = TRUE); env
}
roles <- jsonlite::fromJSON(file.path(repo, "analysis/study_sources/G205_SOURCE_ROLES.json"), simplifyVector = FALSE)$roles
bundle <- lapply(roles, function(r) capture(paste0("analysis/study_sources/", r$path)))
names(bundle) <- vapply(roles, `[[`, "", "role")
worker <- import(capture(spec$worker))
full_raw <- capture(spec$full_config); baseline_raw <- capture(spec$baseline_config)
api <- worker$gse205906_ens101_worker_sources_v1(bundle, full_raw, baseline_raw)
full <- api$config$gse205906_ens101_check_full_model_config_v1(full_raw, baseline_raw)
stopifnot(as.character(getRversion()) == full$software$R)
package_names <- setdiff(names(full$software), c("R", "PowerShell"))
observed <- setNames(lapply(package_names, function(p) as.character(utils::packageVersion(p))), package_names)
stopifnot(identical(observed, full$software[package_names]))
worker$gse205906_ens101_worker_serial_backend_v3(full)
writeLines(capture.output(sessionInfo()), file.path(out, "receipts/sessionInfo.txt"), useBytes = TRUE)
runtime <- list(scope = "CURRENT_PREPARATION_RUNTIME_VERSION_CHECK_NOT_HISTORICAL_LINUX_CUSTODY",
  R = as.character(getRversion()), platform = R.version$platform, packages = observed,
  threads = data.table::getDTthreads(), backend = class(BiocParallel::bpparam()), RNGkind = RNGkind())
jsonlite::write_json(runtime, file.path(out, "receipts/runtime.json"), auto_unbox = TRUE, pretty = TRUE)
input_raw <- read_raw(spec$prepared)
stopifnot(identical(pin(input_raw)$sha256, spec$prepared_identity$sha256), length(input_raw) == spec$prepared_identity$bytes)
prepared <- readRDS(spec$prepared); inputs <- prepared$inputs
stopifnot(identical(names(inputs), c("count_matrix", "gene_annotation", "sample_ledger", "fixed_gene_ids")))
input_before <- serialize(inputs, NULL, version = 3L)
# These are explicitly NEW serializations of authenticated prepared components.
# The fixed ID list is a new derived list, never the historical run_manifest.json.
component_names <- c("count_matrix", "gene_annotation", "sample_ledger", "fixed_gene_ids")
component_pins <- lapply(component_names, function(k) {
  raw <- serialize(inputs[[k]], NULL, version = 3L)
  writeBin(raw, file.path(out, "derived", paste0(k, ".rds"))); pin(raw)
})
names(component_pins) <- c("new_counts", "new_annotation", "new_metadata", "fixed_gene_manifest")
input_pins <- c(component_pins, list(new_config = pin(full_raw), new_prefit_SOURCE = pin(bundle$prefit),
  statistical_runtime_lock = pin(read_raw(file.path(out, "receipts/runtime.json")))))
view <- api$config$gse205906_ens101_prefit_view_v1(full, baseline_raw)
prefit <- api$prefit$gse205906_ens101_effect_independent_prefit_v1(inputs, view,
  bundle[c("model", "helpers", "annotation_adapter")], input_pins)
stopifnot(identical(prefit$status, "PASS_TEN_DESIGN_FILTER_TMM_ONLY"), prefit$passing_phase_n == 10L,
  !prefit$gene_dispersion_or_effect_fit_called, !prefit$program_scoring_or_testing_called)
resource_raw <- lapply(spec$resource_paths, function(p) read_raw(p))
historical_pins <- worker$gse205906_ens101_worker_inputs_v1()
for (k in names(resource_raw)) stopifnot(identical(pin(resource_raw[[k]])$sha256, historical_pins[[k]]$sha256),
  length(resource_raw[[k]]) == historical_pins[[k]]$bytes)
resource <- worker$gse205906_ens101_worker_resources_v1(resource_raw, api)
mapping <- api$mapping$gse205906_ens101_mapping_v2(prefit, resource$resources, resource$orthology,
  worker$gse205906_ens101_worker_mapping_cfg_v1(full), bundle[c("model", "prefit")], resource$pins)
stopifnot(identical(mapping$status, "PASS_TEN_MEMBERSHIP_ONLY"), mapping$passing_phase_n == 10L,
  !mapping$model_fit_called, !mapping$program_scoring_or_testing_called)
rawfive <- inputs[c("count_matrix", "gene_annotation", "sample_ledger")]
rawfive$gene_resources <- resource$resources; rawfive$orthology_table <- resource$orthology
consumer <- api$consumer$load_gse205906_ens101_checked_consumer_v1(full_raw, baseline_raw,
  bundle[seq_len(15L)], prefit, mapping)
core <- worker$gse205906_ens101_worker_import_v1(bundle$core)
handle <- consumer$prepare(rawfive, core$gse205906_ens101_raw_input_identities_v1(rawfive))
receipt <- handle$inspect()
stopifnot(identical(receipt$last_stage, "PREPARED"), !isTRUE(receipt$output_released))
stopifnot(identical(input_before, serialize(inputs, NULL, version = 3L)))
# Persist objects, never live execution closures. Completed scientific outputs stay untouched.
saveRDS(prefit, file.path(out, "derived/prefit.rds"), version = 3L, compress = "gzip")
saveRDS(mapping, file.path(out, "derived/mapping.rds"), version = 3L, compress = "gzip")
# data.table internal self-reference pointers are process-local; R serialization
# preserves the complete scientific content but deliberately not pointer identity.
stopifnot(identical(serialize(prefit, NULL, version = 3L),
    serialize(readRDS(file.path(out, "derived/prefit.rds")), NULL, version = 3L)),
  identical(serialize(mapping, NULL, version = 3L),
    serialize(readRDS(file.path(out, "derived/mapping.rds")), NULL, version = 3L)))
summary <- list(status = "PASS_ORIGINAL_G205_PREPARATION_WITHOUT_MODEL_FITS",
  prepared_input = spec$prepared_identity, phase_count = length(prefit$runs),
  retained_gene_n = lapply(prefit$runs, `[[`, "retained_gene_n"),
  eligible_programs = lapply(mapping$runs, `[[`, "eligible_program_n"),
  resource_validation = resource$summary, consumer_receipt = receipt,
  original_input_objects_unchanged = TRUE, saved_object_serialized_content_equal = TRUE,
  new_real_data_model_fits = 0L, reference_refits = 0L, score_or_test_calls = 0L,
  historical_run_manifest_recreated = FALSE,
  provenance_scope = "New preparation from authenticated public HISAT2 s0 inputs; not original native run custody",
  scientific_scope = "Nine deposited libraries; animal pooling, batch and administered dose remain unresolved",
  future_execution = "Consumer prepared and inspected only; execute() deliberately not called")
jsonlite::write_json(summary, file.path(out, "receipts/PREPARATION_QA.json"), auto_unbox = TRUE,
  pretty = TRUE, null = "null", digits = NA)
writeLines(capture.output(sessionInfo()), file.path(out, "receipts/sessionInfo.txt"), useBytes = TRUE)
