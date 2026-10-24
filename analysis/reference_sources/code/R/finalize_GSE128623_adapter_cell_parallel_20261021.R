# No fits: one not-yet-started adapter six_group, original arithmetic/inspections.
# CLI RUN_ROOT CELL OUTPUT. Original run/contract/cores/helpers remain immutable.
# All 20 accepted exact companions are reused later, never recomputed here.

g128pc_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
g128pc_sha <- function(path) digest::digest(file = path, algo = "sha256")
g128pc_load <- function(path, sha) {
  g128pc_need(identical(g128pc_sha(path), sha), "PinnedFunctionSource")
  expressions <- parse(path, keep.source = FALSE)
  functions <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    length(x) == 3L && is.symbol(x[[2L]]) && is.call(x[[3L]]) &&
    identical(x[[3L]][[1L]], as.name("function")), expressions)
  env <- new.env(parent = asNamespace("stats"))
  for (item in functions) eval(item, env)
  env
}

g128pc_inspector <- function(finalizer, bindings) {
  found <- list()
  walk <- function(node) {
    if (is.call(node) && identical(node[[1L]], as.name("<-")) &&
      length(node) == 3L && identical(node[[2L]], as.name("inspections"))) {
      found[[length(found) + 1L]] <<- node[[3L]]
    }
    # The pinned closure is an assignment in nested for/braced statement blocks.
    # Do not walk expressions/formals/subscripts: they can have missing AST args.
    if (is.call(node) && identical(node[[1L]], as.name("{")))
      for (i in seq.int(2L, length(node))) walk(node[[i]])
    if (is.call(node) && identical(node[[1L]], as.name("for"))) walk(node[[4L]])
  }
  walk(body(finalizer$g128f_run))
  g128pc_need(length(found) == 1L && is.call(found[[1L]]) &&
    identical(found[[1L]][[1L]], as.name("function")), "OneExactOriginalInspectionClosure")
  eval(found[[1L]], list2env(bindings, parent = finalizer))
}

g128pc_summary <- function(index, independent) {
  tail <- independent$tail
  native <- independent$native_diagnostics
  families <- independent$empirical_family_diagnostics
  families$cell_id <- index$spec$cell_id
  families$reference_space <- index$spec$reference_space
  families$required_allocation_n <- index$spec$required_n
  families$empirical_unavailable_n <- vapply(families$family_id, function(key)
    sum(is.na(tail$empirical_p[tail$family_id == key])), integer(1L))
  families$empirical_q_le_05_n <- vapply(families$family_id, function(key)
    sum(tail$empirical_BH_q[tail$family_id == key] <= .05, na.rm = TRUE), integer(1L))
  families$native_fraction_allocations_any_q_le_05 <- vapply(families$family_id, function(key)
    mean(native$native_q_any_le_05[native$family_id == key]), numeric(1L))
  families
}

g128pc_metadata_union <- function(records) {
  answer <- list()
  for (record in records) {
    g128pc_need(is.list(record) && all(c("path", "bytes", "sha256") %in% names(record)),
      "ExplicitMetadataPin")
    old <- answer[[record$path]]
    if (!is.null(old)) g128pc_need(identical(old, record), "ConflictingMetadataIdentity")
    answer[[record$path]] <- record
  }
  answer
}

g128pc_write_json <- function(value, path) {
  g128pc_need(!file.exists(path), "NoReceiptOverwrite")
  temporary <- paste0(path, ".writing_", Sys.getpid())
  g128pc_need(!file.exists(temporary), "NoPartialOverwrite")
  jsonlite::write_json(value, temporary, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA)
  g128pc_need(file.rename(temporary, path), "AtomicReceiptPublication")
}

g128pc_one <- function(admitted, index, bundle, bundle_pin, output, context_sha, finalizer) {
  key <- index$spec$key
  leaf <- file.path(output, paste(index$spec$cell_id, index$spec$reference_space, sep = "__"))
  complete_path <- file.path(leaf, "JOB_COMPLETE.json")
  g128pc_need(!file.exists(complete_path), "CompletedJobMustBeReusedNotReplayed")
  g128pc_need(!dir.exists(leaf), "UncommittedJobRetainedNoAutomaticReplay")
  target <- if (index$spec$kind == "sampled") NULL else sub("_exact$", "", index$spec$reference_space)
  cached <- admitted$helpers$archive$g128r_pack_bundle(bundle, target_contrast = target)
  first <- admitted$helpers$independent$g128i_authenticated_rds(
    index$expected$file[1L], index$expected$sha256[1L])$payload
  for (field in c("keys", "families", "statistic", "native_p", "native_diagnostics",
    "mroast_auxiliary", "rng_calls", "scope")) {
    g128pc_need(identical(first$packed[[field]], cached[[field]]), "OriginalObservedRepackBinding")
  }
  inspections <- g128pc_inspector(finalizer, list(index = index, target = target,
    bundle_pin = bundle_pin, admitted = admitted))
  result <- finalizer$g128f_one_job(index, first$packed, admitted$root, admitted$helpers, inspections)
  g128pc_need(dir.create(leaf), "ExclusiveCellJobOutput")
  finalizer$g128f_write_new(result$producer_tail, file.path(leaf, "PRODUCER_TAIL.rds"))
  finalizer$g128f_write_new(result$independent, file.path(leaf, "INDEPENDENT_ARITHMETIC.rds"))
  lineage <- do.call(rbind, lapply(result$inspections, "[[", "lineage"))
  compact <- do.call(rbind, lapply(result$inspections, "[[", "compact"))
  finalizer$g128f_write_new(lineage, file.path(leaf, "TRANSPORT_LINEAGE.tsv"), FALSE)
  finalizer$g128f_write_new(compact, file.path(leaf, "COMPACT_DIAGNOSTICS_SCOPE.tsv"), FALSE)
  if (index$spec$kind == "sampled") {
    mc <- admitted$helpers$mc$g128mc_summarize(result$independent$tail,
      if (index$spec$cohort == "all30") "all30" else "pre29")
    finalizer$g128f_write_new(mc, file.path(leaf, "SELECTED_OUTER_MC_SIDECAR.rds"))
  }
  summary <- g128pc_summary(index, result$independent)
  finalizer$g128f_write_new(summary, file.path(leaf, "FAMILY_REFERENCE_SUMMARY.tsv"), FALSE)
  extra <- g128pc_metadata_union(unlist(lapply(result$inspections, "[[", "records"), recursive = FALSE))
  for (record in extra) finalizer$g128f_verify(record)
  files <- list.files(leaf, full.names = TRUE)
  receipt <- list(schema = "gse128623-original-arithmetic-cell-job-v1",
    status = "ORIGINAL_ONE_JOB_ARITHMETIC_AND_TRANSPORT_WRITTEN",
    key = key, context_sha256 = context_sha, checked_allocation_n = index$spec$required_n,
    family_n = index$spec$family_n,
    output_records = lapply(files, finalizer$g128f_pin), extra_metadata_records = extra,
    original_inspections_extracted_unchanged = TRUE, original_one_job_function_used = TRUE,
    scientific_model_calls = 0L, original_results_modified = FALSE,
    RNG_final_consumption_verified = FALSE, final_scientific_acceptance = FALSE)
  g128pc_write_json(receipt, complete_path)
  receipt
}

g128pc_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  g128pc_need(length(args) == 3L && args[2L] %in% c("adapter_all30", "adapter_predefined29"),
    "OneExplicitUnstartedAdapterCell")
  g128pc_need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE) &&
    as.character(getRversion()) == "4.6.1" &&
    as.character(utils::packageVersion("digest")) == "0.6.39" &&
    as.character(utils::packageVersion("jsonlite")) == "2.0.0", "ExistingPinnedRuntime")
  g128pc_need(identical(Sys.getenv("LANG"), "C.UTF-8") && identical(Sys.getenv("LC_ALL"), "C.UTF-8") &&
    all(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")) == "1"),
    "OriginalLocaleAndSingleThreads")
  root <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
  g128pc_need(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01"), "ExactRun")
  original <- file.path(root, "finalization_preparation_20261015_r01/finalize_GSE128623_label_reference_20261015.R")
  finalizer <- g128pc_load(original, "2d3298699771c82a983edf4f599c0cfaeb0e210f52f824f95934d9e0220c2240")
  contract_path <- file.path(root, "full_finalization_execution_20261020_r01/CONTRACT.json")
  g128pc_need(g128pc_sha(contract_path) == "35297c8beb13e50e571888f488f7e6a1419606e682cdbf3af00f86c5f325bbb5",
    "OriginalTerminalBoundContract")
  admitted <- finalizer$g128f_preflight(root, finalizer$g128f_json(contract_path))
  admitted$helpers <- finalizer$g128f_load_helpers(admitted$contract$sources)
  cell <- args[2L]
  selected <- Filter(function(index) identical(index$spec$cell_id, cell) &&
    identical(index$spec$reference_space, "six_group"), admitted$jobs)
  g128pc_need(length(selected) == 1L && selected[[1L]]$spec$family_n == 67L &&
    selected[[1L]]$spec$required_n == 10001L && identical(selected[[1L]]$spec$kind, "sampled"),
    "OneFullOriginalSampledJobDomain")
  original_output <- file.path(root, "full_reference_finalization_20261020_r01")
  g128pc_need(!any(dir.exists(vapply(selected, function(x) file.path(original_output,
    paste(x$spec$cell_id, x$spec$reference_space, sep = "__")), character(1L)))),
    "OriginalHasNotPublishedAdapterCell")
  output <- gsub("\\", "/", args[3L], fixed = TRUE)
  parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  g128pc_need(identical(parent, file.path(root, "parallel_finalization_20261021_r01")) &&
    identical(basename(output), cell), "ExactIndependentCellOutput")
  g128pc_need(!file.exists(output), "NewSingleJobOutputOnlyNoReplay")
  self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  g128pc_need(length(self) == 1L, "ExplicitPhysicalCellEntry")
  source_records <- c(list(finalizer$g128f_pin(self), finalizer$g128f_pin(original),
    finalizer$g128f_pin(contract_path)), lapply(admitted$contract$sources, finalizer$g128f_verify))
  observed_pin <- finalizer$g128f_verify(selected[[1L]]$observed)
  observed_receipt <- finalizer$g128f_json(observed_pin$path)
  g128pc_need(identical(observed_receipt$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED"),
    "OriginalAcceptedObservedReceipt")
  bundle_pin <- finalizer$g128f_verify(observed_receipt$observed_bundle)
  metadata <- g128pc_metadata_union(c(admitted$records, list(observed_pin, bundle_pin)))
  runtime <- list(R = R.version.string, libraries = .libPaths(), locale = Sys.getlocale(),
    thread_environment = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")))
  context <- list(schema = "gse128623-original-arithmetic-cell-context-v1", root = root, cell_id = cell,
    source_records = source_records, metadata_records = metadata, runtime = runtime,
    job_keys = names(selected), original_sources_not_modified = TRUE, scientific_model_calls = 0L)
  context_sha <- digest::digest(context, algo = "sha256")
  context_path <- file.path(output, "RUN_CONTEXT.json")
  g128pc_need(dir.create(output), "CreateSingleJobOutput")
  g128pc_write_json(list(context_sha256 = context_sha, context = context), context_path)
  writeLines(capture.output(sessionInfo()), file.path(output, "SESSION_INFO.txt"))
  previous_options <- options(digits = 17)
  on.exit(options(previous_options), add = TRUE)
  bundle <- readRDS(bundle_pin$path)
  results <- lapply(selected, function(index) {
    result <- g128pc_one(admitted, index, bundle, bundle_pin, output, context_sha, finalizer)
    invisible(gc()); result
  })
  for (record in c(metadata, source_records)) finalizer$g128f_verify(record)
  for (result in results) for (record in result$extra_metadata_records) finalizer$g128f_verify(record)
  summary_path <- file.path(output, paste0(cell, "__six_group"), "FAMILY_REFERENCE_SUMMARY.tsv")
  summary <- finalizer$g128f_table(summary_path)
  g128pc_need(nrow(summary) == 67L && !anyDuplicated(summary$family_id), "OriginalJobFamilySummary")
  receipt <- list(schema = "gse128623-original-arithmetic-cell-complete-v1",
    status = "ONE_ORIGINAL_SAMPLED_JOB_WRITTEN_REQUIRES_NATIVE_EXIT_AND_JOIN_QA",
    cell_id = cell, context_sha256 = context_sha, source_records = source_records,
    job_records = lapply(selected, function(index) finalizer$g128f_pin(file.path(output,
      paste(index$spec$cell_id, index$spec$reference_space, sep = "__"), "JOB_COMPLETE.json"))),
    job_n = 1L, family_n = 67L,
    required_archive_n = sum(vapply(selected, function(x) x$spec$required_n, numeric(1L))),
    summary = finalizer$g128f_pin(summary_path), metadata_records = metadata,
    scientific_model_calls = 0L, original_results_modified = FALSE,
    RNG_final_consumption_verified = FALSE, final_scientific_acceptance = FALSE,
    all_20_already_accepted_exact_companions_to_be_reused = TRUE,
    original_parent_handoff_and_whole_24_job_acceptance_still_required = TRUE,
    completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  g128pc_write_json(receipt, file.path(output, "CELL_COMPLETE.json"))
  invisible(receipt)
}

if (sys.nframe() == 0L) g128pc_main()

