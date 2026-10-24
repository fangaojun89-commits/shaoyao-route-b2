# Bounded one-cell/one-space worker. No I/O or model call occurs on import.
# The controller admits the full four-cell run separately; this is not authority.

g128w_need <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

g128w_source_pins <- function() {
  c(gene_camera = "dcb6b4ba4cdc33093df91bd3efcd55297a44339eb4be273c7aa943558ac6b613",
    activity = "12bb16ccfb85cc147180a721519fac47e414eeadfe1b315d997649599179c12b",
    unsigned = "f2efd9e8ba2160dcf01b66771c12345b8c71596a41937d554a745c2b06f0def6",
    bundle = "8f57128488ff832a12466b1359eac16ebb819f23e534c66c5b1d8fe804a1a151",
    archive = "ff4784ea1dc7da2418638eed49031788783c9ae01da712633d5bfbaf9bc4ad24")
}

g128w_bind_sources <- function(paths) {
  expected <- g128w_source_pins()
  g128w_need(is.character(paths) && identical(names(paths), names(expected)) &&
               all(file.exists(paths)), "Supply all five ordered, pinned source files.")
  actual <- vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))
  g128w_need(identical(unname(actual), unname(expected)), "Native/helper source pin changed.")
  environments <- lapply(paths, function(path) {
    environment <- new.env(parent = baseenv())
    sys.source(path, environment)
    environment
  })
  kernels <- list(gene_camera = environments$gene_camera$g128lc_fit_gene_camera,
                  activity = environments$activity$g128a_fit_activity,
                  unsigned = environments$unsigned$g128_unsigned_chronic)
  list(mode = "native_source_pinned", paths = paths, pins = actual,
       fit = environments$bundle$g128b_fit_bundle,
       pack = environments$archive$g128r_pack_bundle,
       validate = environments$archive$g128r_validate_packed,
       write = environments$archive$g128r_write_archive, kernels = kernels)
}

g128w_row_hash <- function(sample_ids, labels) {
  digest::digest(paste(c(paste(sample_ids, collapse = "\t"),
                         paste(unname(labels), collapse = "\t")), collapse = "\n"),
                 algo = "sha256", serialize = FALSE)
}

g128w_validate_allocations <- function(allocations, metadata, cohort, target, synthetic) {
  groups <- c("control", "D4", "D7", "D14", "D28", "D42")
  capacities <- if (cohort == "all30") rep(5L, 6L) else c(5L, 4L, rep(5L, 4L))
  g128w_need(is.matrix(allocations) && is.character(allocations) &&
               !anyNA(allocations) && ncol(allocations) == sum(capacities) &&
               identical(colnames(allocations), as.character(metadata$sample_id)) &&
               identical(rownames(allocations), as.character(seq_len(nrow(allocations)))) &&
               nrow(allocations) >= 2L, "Allocation sample/row identities or dimensions differ.")
  observed <- as.character(metadata$group)
  g128w_need(identical(unname(allocations[1L, ]), observed) &&
               all(allocations %in% groups), "Allocation 1 must be the unchanged observed labels.")
  codes <- matrix(match(allocations, groups), nrow = nrow(allocations))
  g128w_need(all(apply(codes, 1L, function(row) identical(tabulate(row, 6L), capacities))),
             "Allocation changed fixed group capacities.")
  keys <- apply(allocations, 1L, paste, collapse = "\t")
  g128w_need(!anyDuplicated(keys), "Repeated allocation, including an extra observed row.")
  expected_n <- 10001L
  if (!is.null(target)) {
    day <- sub("_vs_control$", "", target)
    changing <- observed %in% c("control", day)
    g128w_need(all(vapply(which(!changing), function(column) {
      all(allocations[, column] == observed[[column]])
    }, logical(1))), "Day companion changed nuisance-group labels.")
    expected_n <- if (cohort == "predefined29" && day == "D4") 126L else 252L
  }
  g128w_need(synthetic || nrow(allocations) == expected_n,
             "Real worker requires the complete frozen allocation space/list.")
  invisible(TRUE)
}

g128w_compact_diagnostics <- function(bundle) {
  diagnostics <- bundle$gene_camera$deseq_diagnostics
  if (is.null(diagnostics)) return(list(native_diagnostics_available = FALSE))
  metadata <- diagnostics$gene_metadata
  diagnostics$gene_metadata <- NULL
  beta <- if (is.data.frame(metadata) && "betaConv" %in% names(metadata)) metadata$betaConv else NULL
  if (!is.null(beta)) names(beta) <- bundle$gene_camera$retained_gene_ids
  weights <- bundle$gene_camera$native$voom$weights
  list(native_diagnostics_available = TRUE, DESeq2 = diagnostics, beta_convergence = beta,
       TMM_samples = bundle$gene_camera$native$dge$samples,
       nonfinite_voom_weight_n = sum(!is.finite(weights)),
       nonpositive_voom_weight_n = sum(is.finite(weights) & weights <= 0),
       integration = bundle$integration_diagnostics)
}

g128w_archive_identity <- function(context, allocation_id, allocation_row_sha256) {
  list(run_id = context$run_id, cell_id = context$cell_id,
       reference_space = context$reference_space, allocation_id = as.integer(allocation_id),
       input_manifest_sha256 = context$input_manifest_sha256,
       source_manifest_sha256 = context$source_manifest_sha256,
       allocation_row_sha256 = allocation_row_sha256)
}

g128w_read_completed <- function(path, expected_identity, observed_keys, observed_families, bindings,
                                  expected_observed_evidence = NULL) {
  payload <- readRDS(path)
  g128w_need(identical(payload$schema_version, "gse128623-statistic-native-p-archive-v1") &&
               identical(payload$identity[names(expected_identity)], expected_identity),
             "Existing checkpoint belongs to another run/input/source/allocation.")
  bindings$validate(payload$packed)
  g128w_need(identical(payload$packed$keys, observed_keys) &&
               identical(payload$packed$families, observed_families) &&
               identical(payload$identity$packed_sha256, digest::digest(payload$packed, algo = "sha256")),
             "Checkpoint domain/content fingerprint changed; do not refit or repair.")
  g128w_need(is.list(payload$packed$worker_diagnostics) &&
               is.list(payload$packed$native_evidence), "Checkpoint diagnostics/native-evidence record missing.")
  diagnostic <- payload$packed$worker_diagnostics
  g128w_need(identical(diagnostic$allocation_id, expected_identity$allocation_id) &&
               identical(paste(diagnostic$processing, diagnostic$cohort, sep = "_"), expected_identity$cell_id) &&
               is.character(diagnostic$labels) && !anyNA(diagnostic$labels) &&
               identical(g128w_row_hash(names(diagnostic$labels), diagnostic$labels), expected_identity$allocation_row_sha256),
             "Checkpoint diagnostic labels/allocation/cell differ from admitted identities.")
  evidence <- payload$packed$native_evidence
  is_sentinel <- identical(expected_identity$reference_space, "six_group") &&
    expected_identity$allocation_id %in% c(2L, 101L, 1001L, 5001L, 10001L)
  if (expected_identity$allocation_id == 1L) {
    g128w_need(identical(evidence$kind, "preserved_observed_bundle") &&
                 identical(evidence$reference, expected_observed_evidence), "Preserved observed evidence reference changed.")
  } else if (is_sentinel) {
    g128w_need(identical(evidence$kind, "fixed_six_group_native_sentinel") &&
                 all(vapply(evidence$bundle[c("gene_camera", "activity", "unsigned")], is.list, logical(1))),
               "Required fixed native sentinel evidence is absent.")
    repacked <- bindings$pack(evidence$bundle, target_contrast = NULL)
    for (field in c("keys", "families", "statistic", "native_p", "native_diagnostics", "mroast_auxiliary", "rng_calls")) {
      g128w_need(identical(repacked[[field]], payload$packed[[field]]), "Native sentinel and packed evidence disagree.")
    }
  } else g128w_need(identical(evidence$kind, "compact_statistics_P_and_diagnostics_only"), "Unexpected native-evidence kind.")
  list(file = normalizePath(path, winslash = "/", mustWork = TRUE),
       sha256 = digest::digest(file = path, algo = "sha256"), bytes = unname(file.info(path)$size))
}

g128w_run_shard <- function(inputs, observed_bundle, allocations, allocation_ids,
                             context, bindings, synthetic = FALSE) {
  g128w_need(is.logical(synthetic) && length(synthetic) == 1L && !is.na(synthetic), "Explicit worker mode required.")
  g128w_need(is.list(context) && context$cohort %in% c("all30", "predefined29") &&
               context$processing %in% c("raw", "adapter") &&
               identical(context$cell_id, paste(context$processing, context$cohort, sep = "_")),
             "Explicit processing/cohort cell identity required.")
  for (field in c("run_id", "cell_id", "reference_space", "worker_id")) {
    g128w_need(is.character(context[[field]]) && length(context[[field]]) == 1L &&
                 !is.na(context[[field]]) && grepl("^[A-Za-z0-9_-]+$", context[[field]]),
               paste("Invalid worker identity:", field))
  }
  for (field in c("input_manifest_sha256", "source_manifest_sha256", "allocation_matrix_sha256",
                  "inputs_object_sha256", "observed_bundle_object_sha256")) {
    g128w_need(is.character(context[[field]]) && length(context[[field]]) == 1L &&
                 !is.na(context[[field]]) && grepl("^[0-9a-f]{64}$", context[[field]]),
               paste("Invalid worker coupling hash:", field))
  }
  target <- if (identical(context$reference_space, "six_group")) NULL else {
    sub("_exact$", "", context$reference_space)
  }
  if (!is.null(target)) g128w_need(target %in% paste0(c("D4", "D7", "D14", "D28", "D42"), "_vs_control") &&
    identical(context$reference_space, paste0(target, "_exact")), "Unknown fixed reference space.")
  root <- normalizePath(context$run_root, winslash = "/", mustWork = TRUE)
  root_pattern <- if (synthetic) "^R:/Shaoyao_storage/cache/GSE128623_[A-Za-z0-9_-]+$" else
    "^/root/autodl-tmp/gse128623_label_reference_[A-Za-z0-9_-]+$"
  g128w_need(grepl(root_pattern, root), "Worker output must stay in its explicit owned run root.")
  g128w_need(is.list(bindings) && all(vapply(bindings[c("fit", "pack", "validate", "write")], is.function, logical(1))),
             "Missing explicit method/archive bindings.")
  if (!synthetic) {
    g128w_need(identical(bindings$mode, "native_source_pinned") &&
                 identical(context$admission$status, "ALL_FOUR_OBSERVED_AND_RESOURCE_PASS") &&
                 identical(context$admission$run_id, context$run_id) &&
                 identical(context$admission$input_manifest_sha256, context$input_manifest_sha256) &&
                 identical(context$admission$source_manifest_sha256, context$source_manifest_sha256),
               "Real execution needs the separately reviewed full-run observed/resource admission.")
  }
  g128w_need(identical(digest::digest(allocations, algo = "sha256"), context$allocation_matrix_sha256),
             "Allocation matrix differs from the explicitly admitted registry-derived object.")
  g128w_need(identical(digest::digest(inputs, algo = "sha256"), context$inputs_object_sha256) &&
               identical(digest::digest(observed_bundle, algo = "sha256"), context$observed_bundle_object_sha256),
             "Input/observed objects differ from their declared caller binding.")
  g128w_validate_allocations(allocations, inputs$metadata, context$cohort, target, synthetic)
  original_labels <- setNames(as.character(inputs$metadata$group), as.character(inputs$metadata$sample_id))
  observed_diagnostic <- observed_bundle$integration_diagnostics
  g128w_need(identical(observed_diagnostic$labels, original_labels) &&
               identical(observed_diagnostic$sample_ids, names(original_labels)) &&
               identical(observed_diagnostic$cohort, context$cohort),
             "Preserved observed bundle is not this cohort's unchanged named labels/sample order.")
  g128w_need(is.numeric(allocation_ids) && length(allocation_ids) > 0L && !anyNA(allocation_ids) &&
               all(is.finite(allocation_ids) & allocation_ids == as.integer(allocation_ids) &
                     allocation_ids >= 1L & allocation_ids <= nrow(allocations)) &&
               !anyDuplicated(allocation_ids), "Invalid/repeated shard allocation IDs.")
  observed_pack <- bindings$pack(observed_bundle, target_contrast = target)
  bindings$validate(observed_pack)
  expected_ids <- observed_bundle$gene_camera$retained_gene_ids
  baseline_E <- observed_bundle$gene_camera$native$voom$E
  observed_counts <- observed_bundle$gene_camera$native$dge$counts
  g128w_need(is.matrix(inputs$counts) &&
               identical(colnames(inputs$counts), names(original_labels)) &&
               all(expected_ids %in% rownames(inputs$counts)) &&
               identical(dimnames(observed_counts), list(expected_ids, names(original_labels))) &&
               identical(as.numeric(observed_counts),
                 as.numeric(round(inputs$counts[expected_ids, , drop = FALSE]))) &&
               identical(dimnames(baseline_E), dimnames(observed_counts)),
             "Preserved observed native count/E domain does not bind this full count input.")
  leaf <- file.path(root, "allocations", context$cell_id, context$reference_space)
  leases <- file.path(root, "leases", context$cell_id, context$reference_space)
  for (directory in c(leaf, leases)) {
    if (!dir.exists(directory)) g128w_need(dir.create(directory, recursive = TRUE), "Could not create owned worker leaf.")
    resolved <- normalizePath(directory, winslash = "/", mustWork = TRUE)
    g128w_need(startsWith(resolved, paste0(root, "/")), "Resolved worker leaf escapes its run root.")
  }
  rows <- vector("list", length(allocation_ids))
  for (index in seq_along(allocation_ids)) {
    allocation_id <- as.integer(allocation_ids[[index]])
    labels <- setNames(unname(allocations[allocation_id, ]), colnames(allocations))
    row_sha <- g128w_row_hash(names(labels), labels)
    expected_identity <- g128w_archive_identity(context, allocation_id, row_sha)
    path <- file.path(leaf, sprintf("allocation_%05d.rds", allocation_id))
    reused <- file.exists(path)
    if (reused) {
      archive <- g128w_read_completed(path, expected_identity, observed_pack$keys, observed_pack$families,
                                      bindings, context$observed_evidence)
    } else {
      lease <- file.path(leases, sprintf("allocation_%05d", allocation_id))
      g128w_need(dir.create(lease, showWarnings = FALSE),
                 "Allocation already has an owner/incomplete attempt; inspect it, never auto-refit.")
      saveRDS(list(run_id = context$run_id, worker_id = context$worker_id,
        allocation_id = allocation_id, pid = Sys.getpid(), allocation_row_sha256 = row_sha),
        file.path(lease, "OWNER.rds"), version = 3L)
      conditions <- list(warning = character(), message = character())
      start <- proc.time()[["elapsed"]]
      # Reuse the already preserved observed fit. Null rows perform their own fit.
      bundle <- if (allocation_id == 1L) observed_bundle else withCallingHandlers(
        do.call(bindings$fit, c(inputs, list(expected_retained_gene_ids = expected_ids,
          baseline_voom_E = baseline_E, labels = labels, cohort = context$cohort,
          kernels = bindings$kernels, target_contrast = target))),
        warning = function(condition) conditions$warning <<- c(conditions$warning, conditionMessage(condition)),
        message = function(condition) conditions$message <<- c(conditions$message, conditionMessage(condition)))
      elapsed <- proc.time()[["elapsed"]] - start
      packed <- bindings$pack(bundle, target_contrast = target)
      g128w_need(identical(packed$keys, observed_pack$keys) &&
                   identical(packed$families, observed_pack$families), "Allocated native result domain changed.")
      packed$worker_diagnostics <- list(labels = labels, allocation_id = allocation_id,
        cohort = context$cohort, processing = context$processing,
        observed_fit_reused = allocation_id == 1L, elapsed_seconds = elapsed,
        unsigned_trace_scope = if (!is.null(target) && length(packed$rng_calls))
          "RETAINED_OBSERVED_CHRONIC_TRACE_NOT_DAY_FIT" else if (!is.null(target))
          "NO_UNSIGNED_DAY_FIT_OR_TRACE" else "NATIVE_CHRONIC_REFERENCE_TRACE",
        peak_memory_measured = FALSE, conditions = conditions,
        native = g128w_compact_diagnostics(bundle))
      is_sentinel <- is.null(target) && allocation_id %in% c(2L, 101L, 1001L, 5001L, 10001L)
      packed$native_evidence <- if (allocation_id == 1L) {
        list(kind = "preserved_observed_bundle", reference = context$observed_evidence)
      } else if (is_sentinel) list(kind = "fixed_six_group_native_sentinel", bundle = bundle) else {
        list(kind = "compact_statistics_P_and_diagnostics_only")
      }
      expected_identity$packed_sha256 <- digest::digest(packed, algo = "sha256")
      archive <- bindings$write(packed, path, root, expected_identity)
      # One atomic archive contains all compact evidence; no metadata codec
      # retry can accidentally refit a completed allocation.
      archive <- g128w_read_completed(path, expected_identity[names(expected_identity) != "packed_sha256"],
                                      observed_pack$keys, observed_pack$families, bindings, context$observed_evidence)
    }
    rows[[index]] <- data.frame(file = archive$file, sha256 = archive$sha256,
      allocation_id = allocation_id, allocation_row_sha256 = row_sha,
      bytes = archive$bytes, reused_checkpoint = reused, stringsAsFactors = FALSE)
  }
  if (identical(bindings$mode, "native_source_pinned")) {
    after <- vapply(bindings$paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))
    g128w_need(identical(after, bindings$pins), "Source files changed during shard; do not accept it.")
  }
  g128w_need(identical(digest::digest(inputs, algo = "sha256"), context$inputs_object_sha256) &&
               identical(digest::digest(observed_bundle, algo = "sha256"), context$observed_bundle_object_sha256),
             "Input/observed objects changed during shard.")
  list(status = "REQUESTED_SHARD_ARCHIVES_WRITTEN_NOT_FULL_RUN_ACCEPTANCE", files = do.call(rbind, rows),
       run_id = context$run_id, cell_id = context$cell_id, reference_space = context$reference_space,
       requested_ids = as.integer(allocation_ids), total_required_ids = nrow(allocations),
       independent_archive_check_complete = FALSE, real_authority_granted_by_helper = FALSE)
}
