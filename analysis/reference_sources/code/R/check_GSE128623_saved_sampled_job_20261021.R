# Saved-output QA only. CLI: ROOT CELL NEW_OUTPUT. No archive/model/RNG replay.
# The fixed RUN_CONTEXT reuses the accepted 918-shard preflight; only relevant
# small metadata and the selected saved result tables are opened here.

g128sj_need <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
}

# Pure comparison: original 1e-12 relative tolerance, preserving NA versus NaN
# and signed infinities. This never replaces or repairs a saved value.
g128sj_equal <- function(left, right, label) {
  if (is.data.frame(left) || is.data.frame(right)) {
    g128sj_need(is.data.frame(left) && is.data.frame(right) &&
      identical(names(left), names(right)) && nrow(left) == nrow(right), label)
    for (name in names(left)) g128sj_equal(left[[name]], right[[name]], label)
  } else if (is.list(left) || is.list(right)) {
    g128sj_need(is.list(left) && is.list(right) && length(left) == length(right) &&
      identical(names(left), names(right)), label)
    for (i in seq_along(left)) g128sj_equal(left[[i]], right[[i]], label)
  } else if (is.numeric(left) || is.numeric(right)) {
    g128sj_need(is.numeric(left) && is.numeric(right) && length(left) == length(right) &&
      identical(names(left), names(right)) && identical(is.na(left), is.na(right)) &&
      identical(is.nan(left), is.nan(right)) &&
      identical(is.infinite(left), is.infinite(right)), label)
    infinite <- is.infinite(left)
    finite <- is.finite(left)
    g128sj_need(all(left[infinite] == right[infinite]) &&
      all(abs(left[finite] - right[finite]) <= 1e-12 * (1 + abs(right[finite]))), label)
  } else {
    g128sj_need(identical(left, right), label)
  }
  invisible(TRUE)
}

g128sj_ids <- function(ids, total, label) {
  g128sj_need(is.numeric(ids) && !anyNA(ids) &&
    identical(as.numeric(ids), as.numeric(seq_len(total))), label)
  invisible(TRUE)
}

# Pure manifest assembly. Inputs are already pinned metadata, not archives.
# Registry-row association is inherited from the unchanged full preflight pins.
g128sj_expected <- function(tasks, manifests, spec, root) {
  g128sj_need(length(tasks) > 0L && length(tasks) == length(manifests), "ExpectedManifestDomain")
  pieces <- vector("list", length(tasks))
  for (i in seq_along(tasks)) {
    task <- tasks[[i]]
    part <- manifests[[i]]
    ids <- seq.int(as.integer(task$first), as.integer(task$last))
    g128sj_need(identical(task$key, spec$key) && min(ids) >= 1L &&
      max(ids) <= spec$required_n && length(ids) <= 50L && is.data.frame(part) &&
      all(c("allocation_id", "file", "sha256", "bytes", "allocation_row_sha256",
        "reused_checkpoint") %in% names(part)), "SelectedManifestSchemaAndTask")
    g128sj_need(identical(as.numeric(part$allocation_id), as.numeric(ids)) &&
      identical(part$file, file.path(root, "allocations", spec$cell_id, "six_group",
        sprintf("allocation_%05d.rds", ids))) &&
      is.character(part$sha256) && !anyNA(part$sha256) &&
      all(grepl("^[a-f0-9]{64}$", part$sha256)) &&
      is.character(part$allocation_row_sha256) && !anyNA(part$allocation_row_sha256) &&
      all(grepl("^[a-f0-9]{64}$", part$allocation_row_sha256)) &&
      is.numeric(part$bytes) && !anyNA(part$bytes) &&
      all(is.finite(part$bytes) & part$bytes > 0 & part$bytes == floor(part$bytes)) &&
      is.logical(part$reused_checkpoint) && !anyNA(part$reused_checkpoint),
      "SelectedManifestOrderedPathsHashesSizes")
    pieces[[i]] <- part[, c("allocation_id", "file", "sha256", "bytes",
      "allocation_row_sha256", "reused_checkpoint"), drop = FALSE]
  }
  expected <- do.call(rbind, pieces)
  expected <- expected[order(expected$allocation_id), , drop = FALSE]
  rownames(expected) <- NULL
  g128sj_ids(expected$allocation_id, spec$required_n, "CompleteManifestIDsOnce")
  g128sj_need(!anyDuplicated(expected$file) && !anyDuplicated(expected$allocation_row_sha256),
    "UniqueExpectedArchiveAndRegistryRows")
  expected
}

# Pure saved-table validation, suitable for invented fixtures with smaller spec.
# mc_fun and summary_fun are supplied frozen deterministic functions, not fitters.
g128sj_check_saved <- function(producer, independent, lineage, compact, mc,
                              expected, spec, template, mc_fun, summary_fun) {
  total <- spec$required_n
  family_n <- spec$family_n
  g128sj_need(identical(spec$reference_space, "six_group") && identical(spec$kind, "sampled") &&
    is.list(independent) && identical(independent$status,
      "PASS_EXPECTED_ARCHIVE_BYTES_IDENTITIES_DIAGNOSTICS_AND_TAILS") &&
    independent$checked_allocation_n == total && independent$checked_family_n == family_n &&
    identical(independent$RNG_final_consumption_verified, FALSE) &&
    identical(independent$full_268_real_family_acceptance, FALSE) &&
    identical(independent$real_execution_authority_granted, FALSE), "SavedIndependentScope")
  receipts <- independent$file_receipts
  receipt_fields <- c("allocation_id", "file", "sha256", "bytes")
  g128sj_need(is.data.frame(receipts) && identical(names(receipts), receipt_fields),
    "SavedFileReceiptSchema")
  g128sj_ids(receipts$allocation_id, total, "SavedFileReceiptCompleteIDs")
  g128sj_need(identical(receipts$file, expected$file) &&
    identical(receipts$sha256, expected$sha256) &&
    is.numeric(receipts$bytes) && !anyNA(receipts$bytes) &&
    identical(as.numeric(receipts$bytes), as.numeric(expected$bytes)),
    "SavedFileReceiptsMatchCurrentPinnedManifests")
  tail <- independent$tail
  g128sj_need(is.data.frame(tail) && nrow(tail) > 0L &&
    is.data.frame(template) && nrow(template) == family_n &&
    !anyDuplicated(template$family_id) &&
    all(c(names(template), "feature_id", "required_total_n", "reference_kind",
      "empirical_p", "empirical_BH_q") %in% names(tail)) &&
    all(vapply(tail[c(names(template), "feature_id")], is.character, logical(1L))) &&
    !anyNA(tail[c(names(template), "feature_id")]) &&
    !anyDuplicated(tail[c(names(template), "feature_id")]) &&
    setequal(tail$family_id, template$family_id) &&
    all(tail$required_total_n == total) && all(tail$reference_kind == "sampled"),
    "SavedTailCompleteFamilyKindAndDenominator")
  at <- match(tail$family_id, template$family_id)
  for (name in names(template)) g128sj_need(
    identical(tail[[name]], template[[name]][at]), "SavedTailFamilyAssociation")
  g128sj_equal(producer, tail, "SavedProducerIndependentTailAgreement")
  family <- independent$empirical_family_diagnostics
  g128sj_need(is.data.frame(family) && all(c(names(template), "n", "actual_BH_adjustment_n") %in%
    names(family)) && nrow(family) == family_n &&
    all(vapply(names(template), function(name) identical(family[[name]], template[[name]]), logical(1L))),
    "SavedFamilyTemplate")  # Saved row names are not scientific family keys.
  feature_n <- vapply(template$family_id, function(key) sum(tail$family_id == key), integer(1L))
  adjustment_n <- vapply(template$family_id, function(key)
    sum(!is.na(tail$empirical_p[tail$family_id == key])), integer(1L))
  g128sj_need(identical(as.numeric(family$n), as.numeric(feature_n)) &&
    identical(as.numeric(family$actual_BH_adjustment_n), as.numeric(adjustment_n)),
    "SavedFamilyFeatureAndBHCounts")
  native <- independent$native_diagnostics
  g128sj_need(is.data.frame(native) && nrow(native) == total * family_n &&
    all(c("allocation_id", names(template), "n", "native_q_any_le_05") %in% names(native)) &&
    !anyNA(native[c("allocation_id", "family_id")]) &&
    setequal(native$family_id, template$family_id) &&
    is.logical(native$native_q_any_le_05) && !anyNA(native$native_q_any_le_05),
    "SavedNativeDiagnosticsSchema")
  for (key in template$family_id) g128sj_ids(native$allocation_id[native$family_id == key],
    total, "SavedNativeEveryFamilyEveryAllocationOnce")
  at <- match(native$family_id, template$family_id)
  for (name in names(template)) g128sj_need(
    identical(native[[name]], template[[name]][at]), "SavedNativeFamilyAssociation")
  g128sj_need(is.numeric(native$n) && !anyNA(native$n) &&
    all(native$n == feature_n[at]), "SavedNativeFeatureDomain")
  g128sj_need(is.data.frame(lineage) && identical(names(lineage),
    c("allocation_id", "mode", "origin", "lineage_status")), "SavedLineageSchema")
  g128sj_ids(lineage$allocation_id, total, "SavedLineageCompleteIDs")
  mode <- ifelse(expected$reused_checkpoint, "legacy", "transport_fresh")
  origin <- ifelse(expected$reused_checkpoint, "ORIGINAL_MAIN01_SUCCESSFUL_FIRST_CREATION",
    "CURRENT_MAIN02_FIRST_CREATION_ACTUAL_EXIT_ZERO")
  for (id in c(1L, 2L, 101L)) if (id <= total && expected$reused_checkpoint[id]) {
    mode[id] <- if (id == 2L) "transport_recovered" else if (id == 101L) "transport_fresh" else "legacy"
    origin[id] <- switch(as.character(id),
      "1" = "PRESERVED_OBSERVED_NO_FIT_VIEW_MAIN02_REUSE_AND_ACCEPTED_OBSERVED_BUNDLE",
      "2" = "PINNED_ZERO_FIT_RECOVERY_AND_AUTHENTICATED_OLD_PARTIAL",
      "101" = "PINNED_FRESH_ID101_MEASURED_EXIT_ZERO")
  }
  g128sj_need(identical(lineage$mode, mode) && identical(lineage$origin, origin) &&
    identical(lineage$lineage_status, ifelse(mode == "legacy", "PASS_LEGACY_IDENTITY_ONLY",
      "PASS_SUPPLEMENTARY_TRANSPORT_LINEAGE_ONLY")), "SavedLineageExpectedScope")
  g128sj_need(is.data.frame(compact) && identical(names(compact), c("allocation_id", "warning_n",
    "message_n", "native_evidence_kind", "full_native_fields_independently_accepted")),
    "SavedCompactScopeSchema")
  g128sj_ids(compact$allocation_id, total, "SavedCompactCompleteIDs")
  kind <- rep("compact_statistics_P_and_diagnostics_only", total)
  kind[intersect(c(2L, 101L, 1001L, 5001L, 10001L), seq_len(total))] <- "fixed_six_group_native_sentinel"
  kind[1L] <- "preserved_observed_bundle"
  g128sj_need(identical(compact$native_evidence_kind, kind) &&
    is.logical(compact$full_native_fields_independently_accepted) &&
    !anyNA(compact$full_native_fields_independently_accepted) &&
    !any(compact$full_native_fields_independently_accepted), "SavedCompactAcceptanceNotPromoted")
  for (name in c("warning_n", "message_n")) g128sj_need(is.numeric(compact[[name]]) &&
    !anyNA(compact[[name]]) && all(is.finite(compact[[name]]) & compact[[name]] >= 0 &
      compact[[name]] == floor(compact[[name]])), "SavedConditionCounts")
  cohort <- if (spec$cohort == "all30") "all30" else "pre29"
  expected_mc <- mc_fun(tail, cohort)
  g128sj_need(is.list(mc) && is.list(mc$context) && identical(mc$context$cohort, cohort) &&
    is.data.frame(mc$intervals) &&
    identical(mc$intervals$original_empirical_p, tail$empirical_p) &&
    identical(mc$intervals$original_empirical_BH_q, tail$empirical_BH_q),
    "SavedMCOriginalPAndQNotChanged")
  for (name in c("theta_hat", "theta_lower", "theta_upper")) g128sj_need(
    is.numeric(mc$intervals[[name]]) && length(mc$intervals[[name]]) == nrow(tail) &&
      all(is.na(mc$intervals[[name]]) | (is.finite(mc$intervals[[name]]) &
        mc$intervals[[name]] >= 0 & mc$intervals[[name]] <= 1)), "SavedMCProbabilityRange")
  g128sj_equal(mc, expected_mc, "SavedMCMatchesFrozenDeterministicSidecar")
  summary <- summary_fun(list(spec = spec), independent)
  g128sj_need(nrow(summary) == family_n && identical(summary$family_id, template$family_id),
    "SavedAuditSummaryDomain")
  summary$saved_output_QA <- "PASS_SAVED_TABLES_AND_PINNED_MANIFEST_BINDING_ONLY"
  summary
}

g128sj_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  g128sj_need(length(args) == 3L && args[2L] %in%
    c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29"), "OneFixedSampledCell")
  g128sj_need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE) &&
    as.character(getRversion()) == "4.6.1" &&
    as.character(utils::packageVersion("digest")) == "0.6.39" &&
    as.character(utils::packageVersion("jsonlite")) == "2.0.0", "PinnedSavedQARuntime")
  g128sj_need(identical(Sys.getenv("LANG"), "C.UTF-8") && identical(Sys.getenv("LC_ALL"), "C.UTF-8") &&
    all(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")) == "1"),
    "PinnedUTF8SingleThreadEnvironment")
  root <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
  g128sj_need(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01"), "FixedRunRoot")
  cell <- args[2L]
  output <- gsub("\\", "/", args[3L], fixed = TRUE)
  g128sj_need(identical(output, file.path(root, "saved_sampled_output_QA_20261021_r01", cell)) &&
    !file.exists(output) && dir.exists(dirname(output)), "NewSpecificSavedQAOutput")
  reader_path <- file.path(root, "parallel_finalization_source_20261021_r01",
    "finalize_GSE128623_adapter_cell_parallel_20261021.R")
  g128sj_need(identical(digest::digest(file = reader_path, algo = "sha256"),
    "18b4289b31d0ed85381359ecb5ca7fd20814446d42ccd4b72ebee074f0fe16df"), "FrozenReaderDefinitions")
  reader <- new.env(parent = asNamespace("stats"))
  sys.source(reader_path, envir = reader)
  original_path <- file.path(root, "finalization_preparation_20261015_r01",
    "finalize_GSE128623_label_reference_20261015.R")
  original <- reader$g128pc_load(original_path,
    "2d3298699771c82a983edf4f599c0cfaeb0e210f52f824f95934d9e0220c2240")
  self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  g128sj_need(length(self) == 1L, "PhysicalSavedQAEntry")
  context_path <- file.path(root, "parallel_finalization_20261021_r01/adapter_all30/RUN_CONTEXT.json")
  context_pin <- original$g128f_pin(context_path)
  g128sj_need(identical(context_pin$sha256,
    "e1be75a269fed58bdffe0efaf9e63c77431b7b5a4767de23741caa4775a08247"), "AcceptedFullPreflightContext")
  context <- original$g128f_json(context_path)$context
  g128sj_need(identical(context$root, root) && identical(context$cell_id, "adapter_all30") &&
    identical(context$schema, "gse128623-original-arithmetic-cell-context-v1") &&
    identical(context$scientific_model_calls, 0L), "AcceptedPreflightContextScope")
  # Both JSON arrays and mixed-named-list JSON objects are lists of pin records.
  pins <- reader$g128pc_metadata_union(c(unname(context$metadata_records), unname(context$source_records)))
  records <- list(context_pin, original$g128f_pin(self))
  bound <- function(path) {
    pin <- pins[[path]]
    g128sj_need(!is.null(pin), "RelevantMetadataMustBeInAcceptedPreflight")
    pin <- original$g128f_verify(pin)
    records[[path]] <<- pin
    pin
  }
  read_json <- function(relative) original$g128f_json(bound(file.path(root, relative))$path)
  contract <- read_json("full_finalization_execution_20261020_r01/CONTRACT.json")
  g128sj_need(identical(digest::digest(file = file.path(root,
    "full_finalization_execution_20261020_r01/CONTRACT.json"), algo = "sha256"),
    "35297c8beb13e50e571888f488f7e6a1419606e682cdbf3af00f86c5f325bbb5"), "OriginalContract")
  for (pin in unname(context$source_records)) {
    original$g128f_verify(pin)
    records[[pin$path]] <- pin
  }
  helpers <- original$g128f_load_helpers(contract$sources)
  config <- read_json("DISPATCH_CONFIG_v2.json")
  for (name in names(config$source_pins)) bound(file.path(root, "source", name))
  bound(file.path(root, "INPUT_MANIFEST.tsv"))
  read_json("JOBS_PREPARED.json")
  admission <- read_json("ADMISSION.json")
  read_json("saved_sentinel_recovery_02/RECOVERY_RECEIPT.json")
  plan <- read_json("dispatch/main_02/PLAN.json")
  terminal <- read_json("dispatch/main_02/DISPATCH_EXIT.json")
  original$g128f_require_terminal(terminal)
  g128sj_need(length(plan$tasks) == 918L && plan$required_archives == 44792L,
    "AcceptedWholeRunMetadataScope")
  specs <- original$g128f_specification()
  spec <- specs[specs$cell_id == cell & specs$reference_space == "six_group", , drop = FALSE]
  g128sj_need(nrow(spec) == 1L && spec$required_n == 10001L && spec$family_n == 67L,
    "Frozen10001By67Domain")
  registry <- read_json("registry/PREPARED.json")
  registry_pin <- bound(file.path(root, "registry", registry$cohorts[[spec$cohort]]$sampled_space$file))
  g128sj_need(identical(registry_pin$sha256,
    registry$cohorts[[spec$cohort]]$sampled_space$sha256), "FrozenSampledRegistry")
  observed <- admission$observed_receipts[[cell]]
  bound(observed$path)  # Do not open the scientific observed bundle or old partial.
  tasks <- Filter(function(task) identical(task$key, spec$key), plan$tasks)
  manifests <- lapply(tasks, function(task) {
    folder <- file.path(root, "dispatch/main_02", task$worker)
    bound(file.path(folder, "STARTED.json"))
    bound(file.path(folder, "PROCESS_EXIT.json"))
    original$g128f_table(bound(file.path(root, "workers", task$worker, "ARCHIVE_MANIFEST.tsv"))$path)
  })
  # Recheck selected legacy/fresh creation metadata already authenticated in the
  # accepted preflight, without opening any SHARD_RECEIPT/model/partial RDS.
  for (pin in pins) {
    suffix <- substring(pin$path, nchar(root) + 2L)
    selected_creation <- startsWith(suffix, paste0("workers/", cell, "_six_group_")) ||
      startsWith(suffix, paste0("dispatch/main_01/", cell, "_six_group_")) ||
      startsWith(suffix, paste0("workers/fresh_transport_", cell, "_00101/")) ||
      startsWith(suffix, paste0("logs/fresh_transport_", cell, "_101_"))
    if (selected_creation && !grepl("[.]rds$", pin$path)) bound(pin$path)
  }
  expected <- g128sj_expected(tasks, manifests, spec, root)
  leaf <- if (startsWith(cell, "raw_")) file.path(root,
    "full_reference_finalization_20261020_r01", paste0(cell, "__six_group")) else
    file.path(root, "parallel_finalization_20261021_r01", cell, paste0(cell, "__six_group"))
  g128sj_need(identical(normalizePath(leaf, winslash = "/", mustWork = TRUE), leaf), "FixedSavedLeaf")
  filenames <- c("PRODUCER_TAIL.rds", "INDEPENDENT_ARITHMETIC.rds", "TRANSPORT_LINEAGE.tsv",
    "COMPACT_DIAGNOSTICS_SCOPE.tsv", "SELECTED_OUTER_MC_SIDECAR.rds")
  output_pins <- lapply(file.path(leaf, filenames), original$g128f_pin)
  producer <- readRDS(output_pins[[1L]]$path)
  independent <- readRDS(output_pins[[2L]]$path)
  lineage <- original$g128f_table(output_pins[[3L]]$path)
  compact <- original$g128f_table(output_pins[[4L]]$path)
  mc <- readRDS(output_pins[[5L]]$path)
  summary <- g128sj_check_saved(producer, independent, lineage, compact, mc, expected, spec,
    helpers$independent$g128i_family_template(), helpers$mc$g128mc_summarize, reader$g128pc_summary)
  extra_names <- c("FAMILY_REFERENCE_SUMMARY.tsv", "JOB_COMPLETE.json")
  for (name in extra_names) if (file.exists(file.path(leaf, name))) {
    pin <- original$g128f_pin(file.path(leaf, name))
    if (name == "FAMILY_REFERENCE_SUMMARY.tsv") {
      saved <- original$g128f_table(pin$path)
      g128sj_equal(saved, summary[setdiff(names(summary), "saved_output_QA")], "AdapterSavedSummaryAgreement")
    } else {
      job <- original$g128f_json(pin$path)
      g128sj_need(identical(job$key, spec$key) && job$checked_allocation_n == 10001L &&
        job$family_n == 67L && identical(job$scientific_model_calls, 0L) &&
        identical(job$final_scientific_acceptance, FALSE), "AdapterSavedJobReceiptScope")
    }
    output_pins[[length(output_pins) + 1L]] <- pin
  }
  for (pin in c(records, output_pins)) original$g128f_verify(pin)
  g128sj_need(dir.create(output), "ExclusiveSavedQAOutput")
  summary_path <- file.path(output, "SAVED_FAMILY_AUDIT.tsv")
  original$g128f_write_new(summary, summary_path, FALSE)
  receipt <- list(schema = "gse128623-saved-sampled-output-qa-v1",
    status = "PASS_SELECTED_SAVED_SAMPLE_JOB_OUTPUTS_AND_PINNED_METADATA_ONLY",
    cell_id = cell, job_key = spec$key, saved_leaf = leaf, checked_allocation_n = 10001L,
    checked_family_n = 67L, saved_input_records = output_pins, metadata_records = records,
    accepted_preflight_context = context_pin, archive_payloads_opened = 0L,
    registry_row_binding = "REUSED_ACCEPTED_PREFLIGHT_WITH_UNCHANGED_MANIFEST_AND_REGISTRY_PINS",
    producer_independent_saved_tails_agree = TRUE, native_source_BH_diagnostic_ID_coverage = TRUE,
    worker_fit_diagnostics_checked_here = FALSE, deterministic_saved_MC_checked = TRUE,
    model_calls = 0L, RNG_replay_calls = 0L, input_or_saved_results_modified = FALSE,
    native_process_exit_verified_here = FALSE, original_finalizer_CLI_success_asserted = FALSE,
    full_native_acceptance_granted = FALSE,
    final_scientific_acceptance = FALSE, whole_24_job_acceptance_still_required = TRUE,
    audit_summary = original$g128f_pin(summary_path), R = R.version.string,
    completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  reader$g128pc_write_json(receipt, file.path(output, "SAVED_OUTPUT_QA.json"))
  invisible(receipt)
}

if (sys.nframe() == 0L) g128sj_main()
