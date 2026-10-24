# Saved-only completion of all 20 already-approved exact day/control spaces.
# RUN_ROOT CONTRACT_JSON CONTRACT_SHA256 NEW_OUTPUT
# Reuse the unchanged full finalizer's declarations; never call its full-run
# preflight/run or manufacture a terminal for the still-running six-group jobs.

g128x_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  stopifnot(length(arguments) == 4L, identical(as.character(getRversion()), "4.6.1"),
    requireNamespace("digest", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
  sha <- function(path) digest::digest(file = path, algo = "sha256")
  root <- normalizePath(arguments[1L], winslash = "/", mustWork = TRUE)
  stopifnot(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01"),
    identical(sha(arguments[2L]), arguments[3L]))
  contract <- jsonlite::read_json(arguments[2L], simplifyVector = FALSE)
  self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  stopifnot(length(self) == 1L, identical(sha(self), contract$entry_sha256),
    identical(contract$schema, "gse128623-complete-exact-companion-contract-v1"),
    identical(contract$run_root, root), length(contract$metadata) > 7L,
    identical(contract$finalizer$sha256,
      "2d3298699771c82a983edf4f599c0cfaeb0e210f52f824f95934d9e0220c2240"),
    identical(sha(contract$finalizer$path), contract$finalizer$sha256))
  finalizer <- new.env(parent = asNamespace("stats"))
  sys.source(contract$finalizer$path, envir = finalizer)
  need <- finalizer$g128f_need
  records <- lapply(contract$metadata, finalizer$g128f_verify)
  paths <- vapply(records, `[[`, character(1L), "path")
  need(!anyDuplicated(paths) && all(startsWith(paths, paste0(root, "/"))),
    "Exact metadata pins must be unique and inside this run.")
  names(records) <- paths
  bound <- function(path) {
    need(path %in% names(records), "Required creation/registry metadata was not caller-pinned.")
    records[[path]]
  }
  read_json <- function(relative) finalizer$g128f_json(bound(file.path(root, relative))$path)
  fixed <- c(PREPARED.rds = "bb58eee81848b9bcac46ba36f57dffd4722e268ddfb52ef5b3a09ba1176b9a65",
    INPUT_MANIFEST.tsv = "1d87d05c46bf18a80cb04636e369b3888847b205a9f81c030af3957193a7cc56",
    JOBS_PREPARED.json = "e07be471f4cadd9d515f9f44ef71475e780e96b0eb79dcb75a95db810b7c8750",
    DISPATCH_CONFIG_v2.json = "847a049b410828b0b920a9f355ef93e34d9d51d30b68cb4559bf1059f20d44a8",
    ADMISSION.json = "2fe4115f6b496c5a00f47806bd3dc8270a05f993d8fdb1d60251fa746f429307",
    "registry/PREPARED.json" = "63be688e2d548b2917479a9ecc6e1dba4b3ab37542ac8c4aca9244e78c3ef5d0")
  for (name in names(fixed)) need(identical(bound(file.path(root, name))$sha256,
    unname(fixed[[name]])), "Original frozen run metadata differs.")
  prepared <- readRDS(bound(file.path(root, "PREPARED.rds"))$path)
  need(identical(prepared$status, "INPUTS_FROZEN_NO_MODELS") &&
    identical(prepared$run_root, root) &&
    identical(prepared$input_manifest_sha256, unname(fixed[["INPUT_MANIFEST.tsv"]])) &&
    identical(prepared$source_manifest_sha256,
      "77578bfd27acca61e79ea2290242ee6ac4f1215bc10fc3d54f6143d55752622e") &&
    identical(digest::digest(prepared$source_identity, algo = "sha256"), prepared$source_manifest_sha256),
    "Preserved preparation/input/source identity differs.")
  config <- read_json("DISPATCH_CONFIG_v2.json")
  jobs <- read_json("JOBS_PREPARED.json")
  plan <- read_json("dispatch/main_02/PLAN.json")
  registry <- read_json("registry/PREPARED.json")
  admission <- read_json("ADMISSION.json")
  need(identical(plan$config_sha256, unname(fixed[["DISPATCH_CONFIG_v2.json"]])) &&
    identical(plan$controller_sha256, "ba55d2ec1d7aa742ccf0c310802271d8a4e72d224db235c90a75de2db5434a40") &&
    length(plan$tasks) == 918L && plan$required_archives == 44792L &&
    identical(config$jobs_prepared_sha256, unname(fixed[["JOBS_PREPARED.json"]])) &&
    identical(jobs$status, "JOBS_PREPARED_NOT_LAUNCHED") &&
    identical(jobs$admission$sha256, unname(fixed[["ADMISSION.json"]])) &&
    identical(admission$status, "ALL_FOUR_OBSERVED_AND_RESOURCE_PASS") &&
    isTRUE(admission$resource_acceptance) && identical(admission$run_id, basename(root)),
    "Original plan/job/admission coupling differs.")
  for (name in names(config$source_pins)) need(identical(bound(file.path(root, "source", name))$sha256,
    config$source_pins[[name]]), "Original deployed execution source differs.")
  specification <- finalizer$g128f_specification()
  specification <- specification[specification$kind == "exhaustive", , drop = FALSE]
  need(nrow(specification) == 20L && sum(specification$required_n) == 4788L &&
    all(specification$family_n == 9L), "Require all 20 complete exact spaces.")
  indices <- list()
  for (row in seq_len(nrow(specification))) {
    spec <- specification[row, ]
    job <- jobs$jobs[[spec$key]]
    need(identical(job$path, file.path(root, "jobs", paste0(spec$cell_id, "_", spec$reference_space, ".rds"))) &&
      file.exists(job$path) && file.info(job$path)$size == job$bytes, "Opaque original job path/size differs.")
    cohort <- registry$cohorts[[spec$cohort]]
    item <- cohort$exact_day_control_spaces[[sub("_vs_control_exact$", "", spec$reference_space)]]
    registry_pin <- bound(file.path(root, "registry", item$file))
    need(identical(registry_pin$sha256, item$sha256), "Frozen exact label registry differs.")
    table <- utils::read.delim(registry_pin$path, stringsAsFactors = FALSE,
      colClasses = "character", check.names = FALSE)
    samples <- as.character(unlist(cohort$sample_ids, use.names = FALSE))
    need(nrow(table) == spec$required_n &&
      identical(table$allocation_id, as.character(seq_len(spec$required_n))) &&
      identical(names(table), c("allocation_id", "is_observed", samples)) &&
      length(samples) == if (spec$cohort == "all30") 30L else 29L,
      "Complete ordered registry/sample domain differs.")
    labels <- as.matrix(table[, samples, drop = FALSE])
    hashes <- apply(labels, 1L, function(values) digest::digest(paste(c(paste(samples, collapse = "\t"),
      paste(unname(values), collapse = "\t")), collapse = "\n"), algo = "sha256", serialize = FALSE))
    need(!anyNA(labels) && !anyDuplicated(hashes), "Missing/repeated exact labels.")
    expected <- data.frame(file = file.path(root, "allocations", spec$cell_id, spec$reference_space,
      sprintf("allocation_%05d.rds", seq_len(spec$required_n))), sha256 = NA_character_,
      allocation_id = seq_len(spec$required_n), allocation_row_sha256 = unname(hashes),
      bytes = NA_real_, reused_checkpoint = NA, creator_worker = NA_character_, stringsAsFactors = FALSE)
    indices[[spec$key]] <- list(spec = spec, expected = expected, labels = labels, sample_ids = samples,
      fixed = list(run_id = basename(root), cell_id = spec$cell_id, reference_space = spec$reference_space,
        input_manifest_sha256 = prepared$input_manifest_sha256, source_manifest_sha256 = prepared$source_manifest_sha256))
  }
  tasks <- Filter(function(task) task$key %in% names(indices), plan$tasks)
  need(length(tasks) == 114L && !anyDuplicated(vapply(tasks, `[[`, character(1L), "worker")),
    "Complete exact-space shard set required.")
  for (task in tasks) {
    index <- indices[[task$key]]
    ids <- seq.int(as.integer(task$first), as.integer(task$last))
    need(min(ids) >= 1L && max(ids) <= index$spec$required_n && length(ids) <= 50L &&
      all(is.na(index$expected$sha256[ids])) &&
      identical(task$worker, paste0("main_02_", index$spec$cell_id, "_", index$spec$reference_space,
        "_", sprintf("%05d", min(ids)), "_", sprintf("%05d", max(ids)))) &&
      identical(task$job, jobs$jobs[[task$key]]$path) &&
      identical(task$job_sha256, jobs$jobs[[task$key]]$sha256), "Task range/job identity differs.")
    dispatch <- file.path(root, "dispatch/main_02", task$worker)
    started <- finalizer$g128f_json(bound(file.path(dispatch, "STARTED.json"))$path)
    ended <- finalizer$g128f_json(bound(file.path(dispatch, "PROCESS_EXIT.json"))$path)
    command <- as.character(unlist(started$command, use.names = FALSE))
    need(identical(started$task, task) && ended$exit_code == 0L && identical(command,
      c("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript", "--vanilla",
        file.path(root, "source/run_GSE128623_reference_shard_transport_v2_20261015.R"),
        file.path(root, "source"), task$job, task$job_sha256,
        as.character(task$first), as.character(task$last), task$worker)), "Actual selected shard exit/entry differs.")
    worker <- file.path(root, "workers", task$worker)
    manifest <- finalizer$g128f_table(bound(file.path(worker, "ARCHIVE_MANIFEST.tsv"))$path)
    need(identical(as.integer(manifest$allocation_id), ids) &&
      identical(manifest$file, index$expected$file[ids]) &&
      identical(manifest$allocation_row_sha256, index$expected$allocation_row_sha256[ids]) &&
      is.logical(manifest$reused_checkpoint) && !anyNA(manifest$reused_checkpoint) &&
      all(grepl("^[a-f0-9]{64}$", manifest$sha256)) &&
      all(file.exists(manifest$file)) && all(file.info(manifest$file)$size == manifest$bytes),
      "Selected final archive/registry coverage differs.")
    finalizer$g128f_bind_shard_receipt(readRDS(bound(file.path(worker, "SHARD_RECEIPT.rds"))$path),
      task, index, manifest, prepared$source_identity)
    index$expected[ids, c("sha256", "bytes", "reused_checkpoint")] <- manifest[, c("sha256", "bytes", "reused_checkpoint")]
    index$expected$creator_worker[ids] <- task$worker
    indices[[task$key]] <- index
  }
  need(all(vapply(indices, function(index) !anyNA(index$expected), logical(1L))),
    "No missing exact-space allocation or shortened denominator is allowed.")
  # Independently bind every reused legacy range before any result is exported.
  for (index in indices) for (first in unique(((which(index$expected$reused_checkpoint) - 1L) %/% 50L) * 50L + 1L)) {
    last <- min(first + 49L, index$spec$required_n)
    worker <- paste(index$spec$cell_id, index$spec$reference_space, sprintf("%05d", first), sprintf("%05d", last), sep = "_")
    old <- finalizer$g128f_table(bound(file.path(root, "workers", worker, "ARCHIVE_MANIFEST.tsv"))$path)
    dispatch <- file.path(root, "dispatch/main_01", worker)
    start <- finalizer$g128f_json(bound(file.path(dispatch, "STARTED.json"))$path)
    end <- finalizer$g128f_json(bound(file.path(dispatch, "PROCESS_EXIT.json"))$path)
    ids <- which(index$expected$reused_checkpoint & index$expected$allocation_id >= first & index$expected$allocation_id <= last)
    at <- match(ids, old$allocation_id)
    command <- as.character(unlist(start$command, use.names = FALSE))
    need(!anyNA(at) && !any(old$reused_checkpoint[at]) &&
      identical(old$file[at], index$expected$file[ids]) && identical(old$sha256[at], index$expected$sha256[ids]) &&
      identical(old$allocation_row_sha256[at], index$expected$allocation_row_sha256[ids]) && end$exit_code == 0L &&
      identical(start$task$key, index$spec$key) && identical(start$task$worker, worker) &&
      identical(start$task$job_sha256, jobs$jobs[[index$spec$key]]$sha256) &&
      start$task$first == first && start$task$last == last && identical(command,
        c("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript", "--vanilla",
          file.path(root, "source/run_GSE128623_reference_shard_json_compat_20261015.R"),
          file.path(root, "source"), jobs$jobs[[index$spec$key]]$path, jobs$jobs[[index$spec$key]]$sha256,
          as.character(first), as.character(last), worker)), "Legacy first-creation range/exit differs.")
  }
  helpers <- finalizer$g128f_load_helpers(contract$sources)
  output <- file.path(normalizePath(dirname(arguments[4L]), winslash = "/", mustWork = TRUE), basename(arguments[4L]))
  need(startsWith(output, paste0(root, "/exact_companion_readback_")) &&
    grepl("^[A-Za-z0-9_-]+$", basename(output)) && !file.exists(output), "Use a fresh exact-companion output.")
  need(dir.create(output), "Could not create exact-companion output.")
  writeLines(capture.output(sessionInfo()), file.path(output, "SESSION_INFO.txt"))
  old_options <- options(digits = 17); on.exit(options(old_options), add = TRUE)
  lineage_context <- list(root = root, helpers = helpers)
  summaries <- list(); output_records <- list(); accepted <- list()
  for (cell in unique(specification$cell_id)) {
    observed_pin <- bound(admission$observed_receipts[[cell]]$path)
    need(identical(observed_pin$sha256, admission$observed_receipts[[cell]]$sha256), "Accepted observed receipt differs.")
    observed <- finalizer$g128f_json(observed_pin$path)
    need(identical(observed$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED"), "Observed bridge status differs.")
    bundle_pin <- bound(observed$observed_bundle$path)
    need(identical(bundle_pin$sha256, observed$observed_bundle$sha256) && bundle_pin$bytes == observed$observed_bundle$bytes,
      "Accepted observed bundle differs.")
    bundle <- readRDS(bundle_pin$path)
    for (index in Filter(function(index) identical(index$spec$cell_id, cell), indices)) {
      target <- sub("_exact$", "", index$spec$reference_space)
      cached <- helpers$archive$g128r_pack_bundle(bundle, target_contrast = target)
      first <- helpers$independent$g128i_authenticated_rds(index$expected$file[1L], index$expected$sha256[1L])$payload
      for (field in c("keys", "families", "statistic", "native_p", "native_diagnostics", "mroast_auxiliary", "rng_calls", "scope"))
        need(identical(first$packed[[field]], cached[[field]]), "Accepted observed no-fit repack differs.")
      inspect <- function(payload, row, id) {
        diagnostic <- payload$packed$worker_diagnostics
        native <- payload$packed$native_evidence
        retained <- id == 1L && length(payload$packed$rng_calls) > 0L
        need(is.list(diagnostic) && identical(as.character(diagnostic$labels), unname(index$labels[id, ])) &&
          identical(names(diagnostic$labels), index$sample_ids) && diagnostic$allocation_id == id &&
          identical(diagnostic$cohort, index$spec$cohort) && identical(diagnostic$processing, index$spec$processing) &&
          identical(diagnostic$observed_fit_reused, id == 1L) &&
          identical(diagnostic$unsigned_trace_scope, if (retained) "RETAINED_OBSERVED_CHRONIC_TRACE_NOT_DAY_FIT" else
            "NO_UNSIGNED_DAY_FIT_OR_TRACE") &&
          identical(native$kind, if (id == 1L) "preserved_observed_bundle" else "compact_statistics_P_and_diagnostics_only"),
          "Saved exact-space labels/evidence/unsigned scope differs.")
        if (id == 1L) need(identical(native$reference$file, bundle_pin$path) &&
          identical(native$reference$sha256, bundle_pin$sha256), "Observed native reference differs.")
        lineage <- finalizer$g128f_transport_expectation(index, row, payload, lineage_context)
        list(lineage = lineage$row, compact = data.frame(allocation_id = id,
          warning_n = length(diagnostic$conditions$warning), message_n = length(diagnostic$conditions$message),
          native_evidence_kind = native$kind, full_native_fields_independently_accepted = FALSE,
          stringsAsFactors = FALSE))
      }
      result <- finalizer$g128f_one_job(index, first$packed, root, helpers, inspect)
      leaf <- file.path(output, paste(index$spec$cell_id, index$spec$reference_space, sep = "__"))
      need(dir.create(leaf), "Exact job output exists.")
      finalizer$g128f_write_new(result$producer_tail, file.path(leaf, "PRODUCER_TAIL.rds"))
      finalizer$g128f_write_new(result$independent, file.path(leaf, "INDEPENDENT_ARITHMETIC.rds"))
      finalizer$g128f_write_new(do.call(rbind, lapply(result$inspections, `[[`, "lineage")), file.path(leaf, "TRANSPORT_LINEAGE.tsv"), FALSE)
      finalizer$g128f_write_new(do.call(rbind, lapply(result$inspections, `[[`, "compact")), file.path(leaf, "COMPACT_DIAGNOSTICS_SCOPE.tsv"), FALSE)
      tail <- result$independent$tail; native <- result$independent$native_diagnostics
      families <- result$independent$empirical_family_diagnostics
      families$cell_id <- index$spec$cell_id; families$reference_space <- index$spec$reference_space
      families$required_allocation_n <- index$spec$required_n
      families$empirical_unavailable_n <- vapply(families$family_id, function(key) sum(is.na(tail$empirical_p[tail$family_id == key])), integer(1L))
      families$empirical_q_le_05_n <- vapply(families$family_id, function(key) sum(tail$empirical_BH_q[tail$family_id == key] <= .05, na.rm = TRUE), integer(1L))
      families$native_fraction_allocations_any_q_le_05 <- vapply(families$family_id, function(key) mean(native$native_q_any_le_05[native$family_id == key]), numeric(1L))
      summaries[[index$spec$key]] <- families
      output_records <- c(output_records, lapply(list.files(leaf, full.names = TRUE), finalizer$g128f_pin))
      accepted[[index$spec$key]] <- list(allocations = index$spec$required_n, families = 9L,
        feature_rows = nrow(tail), producer_independent_tails_and_BH_agree = TRUE)
      cat(jsonlite::toJSON(list(job_complete = index$spec$key, allocations = index$spec$required_n), auto_unbox = TRUE), "\n")
      rm(result, first, cached, tail, native); invisible(gc())
    }
    rm(bundle); invisible(gc())
  }
  for (record in c(records, contract$sources, list(contract$finalizer))) finalizer$g128f_verify(record)
  need(identical(sha(arguments[2L]), arguments[3L]) && identical(sha(self), contract$entry_sha256), "Caller/entry changed.")
  summary <- do.call(rbind, summaries); rownames(summary) <- NULL
  need(nrow(summary) == 180L && length(accepted) == 20L &&
    sum(vapply(accepted, `[[`, integer(1L), "allocations")) == 4788L, "Complete exact family/allocation inventory differs.")
  summary_path <- file.path(output, "FAMILY_REFERENCE_SUMMARY.tsv")
  finalizer$g128f_write_new(summary, summary_path, FALSE)
  receipt <- list(schema = "gse128623-complete-exact-companion-readback-v1",
    status = "PASS_ALL_20_EXACT_COMPANION_ARITHMETIC_AND_LINEAGE_ONLY",
    completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), run_root = root,
    job_n = 20L, archive_n = 4788L, family_slots = 180L, jobs = accepted,
    selected_114_actual_shards_exit_zero = TRUE, all_selected_metadata_and_sources_unchanged = TRUE,
    caller_contract = finalizer$g128f_pin(arguments[2L]), entry = finalizer$g128f_pin(self),
    unchanged_finalizer = contract$finalizer, helper_sources = contract$sources,
    metadata_records = records, output_records = c(output_records, list(finalizer$g128f_pin(summary_path))),
    source_native_observed_repack_fields_match = TRUE, archive_bytes_authenticated_each_arithmetic_pass = TRUE,
    producer_independent_complete_tails_and_BH_agree = TRUE, scientific_model_calls = 0L,
    original_results_modified = FALSE, exact_companions_receive_outer_MC_intervals = FALSE,
    full_24_job_finalization_invoked = FALSE, six_group_reference_accepted = FALSE,
    full_native_fields_or_final_RNG_consumption_independently_reestimated = FALSE,
    biological_exchangeability_authenticated = FALSE, candidate_gate_changed = FALSE,
    whole_study_complete = FALSE, R = R.version.string, library_paths = .libPaths(), command = commandArgs())
  jsonlite::write_json(receipt, file.path(output, "COMPANION_READBACK_COMPLETE.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  invisible(receipt)
}

if (sys.nframe() == 0L) g128x_main()
