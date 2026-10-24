# Separate executor, NOT an edit/relabel of the old CLI or scientific identity.
# No driver/finish/parallel expansion. Definitions only when sourced.
# CLI: --execute CONFIG ACCESSION four_group_MC FIRST LAST PREVIOUS_PREFIX
#      PREVIOUS_BYTES PREVIOUS_SHA LEGACY_ACCEPTANCE LEGACY_SHA EXECUTOR_ACCEPTANCE EXECUTOR_SHA
# EXECUTOR_ACCEPTANCE is an independently produced ONE-JOB permission. Its
# evidence pins authenticate review/resources/the previous native boundary;
# this entry never creates its own execution permission or claims native exit.
jgwsp_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwsp_old_sha <- "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d"
jgwsp_helper_sha <- "74bb8dc0d3cdb16bb9373fdf5a9e6153af7298d5c104276e58c3247fc41b972e"
jgwsp_root <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03"

# Deliberate copy of old jgwl_open. AST must differ at ONLY the entry assignment:
# entry here is the LEGACY SCIENTIFIC CONTRACT, not the actual running CLI.
# Every other original guard is retained. The real executor is pinned separately.
jgwsp_legacy_scientific_context <- function(config_path, accession) {
  config_pin <- jgwl_pin(config_path)
  cfg <- jsonlite::read_json(config_pin$path, simplifyVector = FALSE)
  project <- normalizePath(cfg$project_root, winslash = "/", mustWork = TRUE)
  run <- normalizePath(cfg$run_root, winslash = "/", mustWork = TRUE)
  jgwl_need(identical(cfg$schema, "jgw03-separate-linux-reference-config-v2") &&
    identical(run, "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03") &&
    identical(project, paste0(run, "/project")) && identical(config_pin$path, paste0(run, "/CONFIG.json")) &&
    accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    setequal(names(cfg$cells), c("GSE289929", "GSE227494", "GSE299063")),
    "Exact separate Linux root/config/cell contract required")
  relative <- names(cfg$source_pins)
  jgwl_need(length(relative) > 0L && !anyDuplicated(relative) &&
    all(nzchar(relative)) && all(!grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))", relative)),
    "Project-relative source manifest required")
  sources <- lapply(seq_along(relative), function(i) {
    rec <- cfg$source_pins[[i]]; rec$path <- file.path(project, relative[i])
    pin <- jgwl_record(rec)
    jgwl_need(jgwl_inside(pin$path, project), "Source escapes staged project")
    pin
  }); names(sources) <- relative
  required <- c("code/R/prepare_JGW03_reference_first_nonidentity_20261015.R",
    "code/R/run_JGW03_full_reference_single_20261015.R",
    "code/R/prepare_JGW03_observed_label_pilot_v3_20261015.R",
    "code/R/JGW03_label_reference_core_20261015.R", "code/R/JGW03_input_admission_v2_20261015.R",
    "code/R/JGW03_reference_field_map_20261015.R", "code/R/JGW03_reference_compact_tail_20261015.R",
    "code/R/JGW03_finite_population_mc_v2_20261015.R", "code/R/analyze_jgw03_three_cell_lines.R",
    "code/R/run_JGW03_full_reference_linux_20261017.R",
    "code/R/run_JGW03_full_reference_linux_v3_20261017.R")
  jgwl_need(all(required %in% relative), "Required immutable source absent from manifest")
  jgwl_need(is.character(cfg$checker_source) && length(cfg$checker_source) == 1L &&
    cfg$checker_source %in% relative, "Independent checker must be source-pinned")
  jgwl_need(identical(sources[["code/R/run_JGW03_full_reference_linux_20261017.R"]]$sha256,
    "57e3643a06b5c17f66c46eb658d150665c9bc55a1b4cc025ceb3e21d6d2cdc8f"), "Preserved v1 entry differs")
  definition <- jgwl_record(cfg$compatibility_definition); approval <- jgwl_record(cfg$approval)
  jgwl_need(identical(definition$sha256,
    "5c06cb8e054c4b2e3a261c2250944c09939d6b6851a65043222ac2360de3aff8"), "Approved compatibility definition differs")
  origin <- jgwl2_origin_records(cfg$reused_first_PC12)
  entry <- jgwl_pin(file.path(project, "code/R/run_JGW03_full_reference_linux_v3_20261017.R"))
  jgwl_need(identical(entry, sources[[tail(required, 1L)]]), "Executing entry differs from manifest")
  # Both old entries are source-safe. Their main/bootstrap is NOT invoked.
  e <- new.env(parent = asNamespace("stats"))
  sys.source(sources[[required[1L]]]$path, envir = e)
  modules <- e$jgwx_modules(project)
  jgwl_need(identical(sources[[required[2L]]]$sha256,
    "f1de9c41003e0a96c6e08e3367a3e34f3a9ce898cce7a0be694bf9d4965ff4e5"),
    "Original full worker differs")
  sys.source(sources[[required[2L]]]$path, envir = e)
  runtime <- jgwl_runtime(modules$helper, cfg$expected_runtime)
  inputs <- jgwl_records(cfg$input_pins)
  jgwl_need(all(vapply(inputs, function(x) jgwl_inside(x$path, project), logical(1))),
    "Scientific inputs must belong to staged project")
  cell <- cfg$cells[[accession]]
  ref <- lapply(cell[c("registry", "windows_observed", "windows_nonidentity", "windows_fixed")], jgwl_record)
  jgwl_need(length(ref) == 4L && all(vapply(ref, function(x) jgwl_inside(x$path, run), logical(1))),
    "Explicit reference copies must belong to separate run")
  receipts <- jgwl_records(cell$windows_receipts)
  jgwl_need(all(c("prepared", "observed_complete", "nonidentity_complete", "nonidentity_saved") %in%
    names(receipts)) && all(vapply(receipts, function(x) jgwl_inside(x$path, run), logical(1))),
    "Windows provenance receipts required")
  metadata <- lapply(receipts, function(x) jsonlite::read_json(x$path, simplifyVector = FALSE))
  equal_pin <- function(left, right) left$bytes == right$bytes && identical(left$sha256, right$sha256)
  wp <- metadata$prepared; wo <- metadata$observed_complete
  wn <- metadata$nonidentity_complete; ws <- metadata$nonidentity_saved
  jgwl_need(identical(wo$accession, accession) && identical(wn$accession, accession) &&
    identical(wo$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_ONLY") && isTRUE(wo$bridge$all_fields_pass) &&
    identical(wn$status, "FIRST_NONIDENTITY_SAVED_AND_COMPACTED_ONLY") && wn$allocation_id == 2 &&
    equal_pin(ref$registry, wp$label_registry) && equal_pin(ref$windows_observed, wo$result) &&
    equal_pin(ref$windows_fixed, wp$observed_files$fixed) && equal_pin(ref$windows_nonidentity, wn$native) &&
    equal_pin(ref$windows_nonidentity, ws$native) &&
    identical(wn$prepared_sha256, wp$prepared$sha256) &&
    identical(ws$prepared_sha256, wp$prepared$sha256) && identical(ws$allocation_key, wn$allocation_key),
    "Windows reference-to-receipt association differs")
  # A serialized registry is imported byte-for-byte, never regenerated here.
  registry <- readRDS(ref$registry$path)
  keys <- e$jgwx_labels(registry)
  jgwl_need(identical(keys$four_group_MC[2L], wn$allocation_key), "Fixed Windows ID2 key differs")
  expected_inputs <- modules$admission$jgw03_v2_expected_files(accession)
  ledger_path <- file.path(project, "metadata/samples.tsv")
  mandatory <- c(file.path(project, expected_inputs$relative_path), ledger_path)
  jgwl_need(all(mandatory %in% vapply(inputs, `[[`, character(1L), "path")),
    "Eight frozen scientific inputs and independent sample ledger must be manifest-bound")
  before <- list(config = config_pin, sources = sources, inputs = inputs, references = ref,
    windows_receipts = receipts, entry = entry, runtime = runtime,
    compatibility_definition = definition, approval = approval, original_PC12 = origin)
  list(cfg = cfg, project = project, run = run, accession = accession, e = e, modules = modules,
    h = modules$helper, ref = ref, metadata = metadata, registry = registry, keys = keys,
    cell_dir = file.path(run, "intermediate", accession), before = before)
}

jgwsp_same_pin <- function(left, right) {
  is.list(left) && is.list(right) && identical(left$path, right$path) &&
    is.numeric(left$bytes) && length(left$bytes) == 1L && is.finite(left$bytes) &&
    is.numeric(right$bytes) && length(right$bytes) == 1L && is.finite(right$bytes) &&
    left$bytes == right$bytes && identical(left$sha256, right$sha256)
}

# Pure permission binding. Authentication/rehashing of its evidence is done by
# main before the first new allocation. Legacy acceptance alone is insufficient.
jgwsp_permission <- function(value, z, actual, helper, previous, legacy, first, last) {
  jgwsp_need(is.list(value) &&
    identical(value$schema, "jgw03-incremental-prefix-executor-acceptance-v1") &&
    identical(value$status, "ACCEPTED_FOR_ONE_BOUNDED_INCREMENTAL_RUN") &&
    isTRUE(value$resources_admitted) && isTRUE(value$independent_wiring_review_accepted) &&
    isTRUE(value$saved_prefix_boundary_comparison_passed) &&
    isTRUE(value$previous_boundary_authenticated_by_caller) &&
    identical(value$maximum_parallel_single_core_workers, 3L) &&
    is.numeric(value$per_worker_memory_bytes) && length(value$per_worker_memory_bytes) == 1L &&
    is.finite(value$per_worker_memory_bytes) && value$per_worker_memory_bytes == 3*1024^3,
    "SeparateReviewedThreeWorkerExecutorAcceptanceRequired")
  jgwsp_need(jgwsp_same_pin(value$actual_executor, actual) &&
    jgwsp_same_pin(value$legacy_scientific_entry, z$before$entry) &&
    jgwsp_same_pin(value$config, z$before$config) && jgwsp_same_pin(value$legacy_acceptance, legacy) &&
    jgwsp_same_pin(value$increment_helper, helper) && identical(helper$sha256, jgwsp_helper_sha) &&
    identical(value$job$accession, z$accession) && identical(value$job$space, "four_group_MC") &&
    is.numeric(value$job$first) && length(value$job$first) == 1L &&
    is.numeric(value$job$last) && length(value$job$last) == 1L &&
    value$job$first == first && value$job$last == last &&
    jgwsp_same_pin(value$job$previous_prefix, previous), "OneJobExactSourceAndPrefixPermission")
  jgwsp_need(is.list(value$evidence) && length(value$evidence) == 3L &&
    !anyDuplicated(names(value$evidence)) && setequal(names(value$evidence),
    c("review", "resources", "previous_boundary")), "ExplicitIndependentEvidencePinsRequired")
  invisible(TRUE)
}

jgwsp_range_guard <- function(previous, b, first, last) {
  whole <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x == floor(x)
  jgwsp_need(whole(first) && whole(last) && first >= 3 && first <= last &&
    last <= 10001 && last-first+1 <= 50 &&
    is.list(previous) && identical(previous$schema, "jgw03-fixed-registry-tail-prefix-v1") &&
    identical(previous$space, "four_group_MC") && whole(previous$through_allocation_id) &&
    first == previous$through_allocation_id + 1 &&
    identical(previous$scientific_identity_sha256, b$identity_sha256) &&
    identical(previous$full_registry, b$context$label_registry) &&
    identical(previous$observed, b$context$observed_files$observed) &&
    identical(previous$state$allocation_keys, b$keys$four_group_MC) &&
    identical(previous$state$processed_nonidentity, as.integer(first-2)) &&
    identical(previous$provisional_P_computed, FALSE), "StrictNextUnchangedMCSpanAtMost50")
  ids <- seq.int(2L, previous$through_allocation_id)
  pins <- previous$completed_compact_pins
  jgwsp_need(is.list(pins) && identical(names(pins), as.character(ids)) &&
    all(vapply(ids, function(i) is.list(pins[[as.character(i)]]) &&
      identical(pins[[as.character(i)]]$allocation_key, b$keys$four_group_MC[i]), logical(1L))),
    "CompletePreexistingPrefixKeyPinsBeforeNewNativeCall")
  invisible(TRUE)
}

# Wiring only: the old native callback remains the sole producer. Completed
# draws are reused by that callback; incomplete leaves stop without auto-refit.
jgwsp_range <- function(z, b, previous, first, last, actual, fixed, advance,
                        leaf_exists, get_draw, read_archive, pin, write_prefix) {
  jgwsp_range_guard(previous, b, first, last)
  old_identity <- b$identity
  old_identity_sha <- b$identity_sha256
  b$entry <- actual  # Actual STARTED/archive actor, NOT the scientific contract.
  # Inspect only existing leaves in this bounded requested range before any fit.
  # The original reader rejects partial leaves. A complete OLD actor is not a
  # successor-owned checkpoint, even though its scientific identity is the same.
  for (i in seq.int(first, last)) {
    key <- b$keys$four_group_MC[i]
    if (leaf_exists(file.path(b$store, key))) {
      prior <- read_archive(b, key)
      jgwsp_need(is.list(prior) && identical(prior$scientific_identity_sha256, old_identity_sha) &&
        identical(prior$allocation_key, key) && jgwsp_same_pin(prior$execution_entry, actual),
        "PreexistingSuccessorRangeMustHaveActualNewExecutor")
    }
  }
  records <- vector("list", last-first+1)
  for (i in seq.int(first, last)) {
    labels <- stats::setNames(as.numeric(b$registry$MC[i, ]), colnames(b$registry$MC))
    b$modules$kernel$jgw03_check_allocation(labels, b$domain$sample_ids)
    key <- b$keys$four_group_MC[i]
    existed <- leaf_exists(file.path(b$store, key))
    draw <- get_draw(b, key, labels, fixed, z$cell_dir, TRUE)
    jgwsp_need(identical(draw$scientific_identity_sha256, old_identity_sha) &&
      identical(draw$allocation_key, key) && jgwsp_same_pin(draw$execution_entry, actual) &&
      is.integer(draw$new_native_calls) && length(draw$new_native_calls) == 1L &&
      draw$new_native_calls %in% c(0L, 1L),
      "EverySuccessorRangeArchiveHasActualNewExecutor")
    records[[i-first+1]] <- list(allocation_id = i, allocation_key = key,
      execution_entry = draw$execution_entry, reused_completed_archive = existed,
      new_native_calls_this_invocation = if (existed) 0L else draw$new_native_calls)
    if (i %% 10L == 0L) invisible(gc())
  }
  prefix <- advance(previous, b, last, read_archive, b$h$jgw03c_tail_update, pin)
  jgwsp_need(identical(b$identity, old_identity) && identical(b$identity_sha256, old_identity_sha),
    "FrozenScientificIdentityNeverRelabelled")
  target <- file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", last))
  product <- write_prefix(prefix, target)
  list(prefix = prefix, product = product, allocation_records = records,
    legacy_scientific_identity_sha256 = old_identity_sha, actual_executor = actual)
}

jgwsp_main <- function(args) {
  jgwsp_need(length(args) == 13L && identical(args[1L], "--execute") &&
    args[3L] %in% c("GSE289929", "GSE227494", "GSE299063") &&
    identical(args[4L], "four_group_MC") && all(grepl("^[0-9]+$", args[c(5L, 6L, 8L)])),
    "ExplicitOneJobMCExecutionArgumentsRequired")
  old_path <- file.path(jgwsp_root, "project/code/R/run_JGW03_full_reference_linux_v3_20261017.R")
  jgwsp_need(identical(digest::digest(file = old_path, algo = "sha256"), jgwsp_old_sha),
    "ImmutableLegacyScientificEntryRequired")
  old <- new.env(parent = asNamespace("stats")); sys.source(old_path, old)
  actual <- old$jgwl_entry()  # Here --file is truthfully the NEW executor.
  jgwsp_need(old$jgwl_inside(actual$path, jgwsp_root) &&
    !old$jgwl_inside(actual$path, file.path(jgwsp_root, "project")), "NewExecutorOutsideFrozenProject")
  loader <- jgwsp_legacy_scientific_context
  environment(loader) <- old  # Bind original imports; never replace old$jgwl_open.
  z <- loader(args[2L], args[3L])
  jgwsp_need(identical(z$before$entry$sha256, jgwsp_old_sha), "LegacyContractNotActualCLI")
  legacy <- old$jgwl_accept(z, args[10L], args[11L])
  permission_pin <- old$jgwl_record(list(path = args[12L], bytes = unname(file.info(args[12L])$size), sha256 = args[13L]))
  permission <- jsonlite::read_json(permission_pin$path, simplifyVector = FALSE)
  previous_pin <- old$jgwl_record(list(path = args[7L], bytes = as.numeric(args[8L]), sha256 = args[9L]))
  helper <- old$jgwl_record(permission$increment_helper)
  first <- as.numeric(args[5L]); last <- as.numeric(args[6L])
  jgwsp_permission(permission, z, actual, helper, previous_pin, legacy, first, last)
  evidence <- lapply(permission$evidence, old$jgwl_record)
  jgwsp_need(old$jgwl_inside(helper$path, z$run) &&
    all(vapply(evidence, function(x) old$jgwl_inside(x$path, z$run), logical(1L))), "RunBoundNewExecutorEvidence")
  h <- new.env(parent = baseenv()); sys.source(helper$path, h)
  opened <- old$jgwl_boot(z); b <- opened$b
  previous <- readRDS(previous_pin$path)
  jgwsp_range_guard(previous, b, first, last)
  expected_previous <- file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", first-1))
  jgwsp_need(identical(previous_pin$path, expected_previous), "CanonicalPreviousPrefixForSameCell")
  target <- file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", last))
  jgwsp_need(!file.exists(target) && !file.exists(paste0(target, ".partial")), "NoExistingTargetOrAutomaticRetry")
  output <- permission$output_leaf
  parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  jgwsp_need(old$jgwl_inside(parent, file.path(z$run, "provenance")) &&
    identical(output, file.path(parent, basename(output))) && !basename(output) %in% c(".", "..") &&
    !file.exists(output) && dir.create(output), "ExclusiveNewDiagnosticReceiptLeaf")
  sources <- list(actual, z$before$entry, legacy, permission_pin, helper, previous_pin)
  old$jgwl_json(list(status = "NEW_EXECUTOR_ONE_RANGE_STARTED_NOT_NATIVE_EXIT",
    actual_executor = actual, legacy_scientific_entry = z$before$entry,
    scientific_identity_sha256 = b$identity_sha256, first = first, last = last,
    previous_prefix = previous_pin, executor_acceptance = permission_pin), file.path(output, "STARTED.json"))
  completed <- FALSE
  on.exit(if (!completed) old$jgwl_json(list(status = "STOP_RETAINED_NO_AUTOMATIC_REFIT",
    actual_executor = actual, actual_native_exit_observed_here = FALSE), file.path(output, "FAILED.json")), add = TRUE)
  fixed <- readRDS(b$context$fixed$path)
  result <- jgwsp_range(z, b, previous, first, last, actual, fixed,
    h$jgw03_advance_saved_tail_prefix, dir.exists, b$e$jgwf_get_draw,
    b$e$jgwf_read_archive, b$e$jgwx_pin, b$e$jgwx_atomic_plain)
  old$jgwl_after(z)
  invisible(lapply(c(sources, evidence), old$jgwl_record))
  receipt <- list(schema = "jgw03-incremental-prefix-executor-result-v1",
    status = "BOUNDED_INCREMENTAL_PREFIX_PRODUCER_COMPLETE_NOT_INDEPENDENT_ACCEPTANCE",
    accession = z$accession, space = "four_group_MC", first = first, last = last,
    actual_executor = actual, legacy_scientific_entry = z$before$entry,
    scientific_identity_sha256 = b$identity_sha256, previous_prefix = previous_pin,
    new_prefix = c(list(path = target), result$product), allocation_records = result$allocation_records,
    executor_acceptance = permission_pin, legacy_acceptance = legacy, evidence = evidence,
    old_CLI_success_asserted = FALSE, actual_native_exit_observed_here = FALSE,
    previous_boundary_native_exit_authenticated_by_external_driver_not_revalidated_here = TRUE,
    scientific_acceptance = FALSE, checker_requires_exact_entry_segment_mapping = TRUE,
    existing_scientific_source_config_or_acceptance_modified = FALSE,
    maximum_parallel_single_core_workers = 3L, R = R.version.string,
    session_info = utils::capture.output(utils::sessionInfo()))
  old$jgwl_json(receipt, file.path(output, "R_COMPLETE.json"))
  completed <- TRUE
  0L
}

if (sys.nframe() == 0L) {
  exit <- tryCatch(jgwsp_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_INCREMENTAL_SUCCESSOR_STOP_RETAINED_NO_AUTOMATIC_REFIT")
    1L
  })
  quit(status = exit)
}
