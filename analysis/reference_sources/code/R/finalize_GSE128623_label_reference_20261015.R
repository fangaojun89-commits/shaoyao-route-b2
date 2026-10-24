# No fits: connect the existing 24-job archives to producer/independent arithmetic.
# CLI: RUN_ROOT CONTRACT_JSON CONTRACT_SHA256 NEW_OUTPUT
# A caller-pinned contract supplies named metadata pins, not new scientific rules.

g128f_need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

g128f_sha <- function(path) digest::digest(file = path, algo = "sha256")

g128f_pin <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = g128f_sha(path))
}

g128f_verify <- function(record) {
  g128f_need(is.list(record) && all(c("path", "bytes", "sha256") %in% names(record)),
    "Explicit path/bytes/SHA record required.")
  actual <- g128f_pin(record$path)
  g128f_need(identical(actual$sha256, record$sha256) && actual$bytes == record$bytes,
    paste("Pinned file differs:", basename(record$path)))
  actual
}

g128f_json <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)

g128f_table <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
    quote = "\"", comment.char = "", na.strings = "NA")
}

g128f_specification <- function() {
  cells <- c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
  spaces <- c("six_group", paste0(c("D4", "D7", "D14", "D28", "D42"), "_vs_control_exact"))
  result <- do.call(rbind, lapply(cells, function(cell) {
    data.frame(cell_id = cell, reference_space = spaces, stringsAsFactors = FALSE)
  }))
  result$cohort <- ifelse(grepl("predefined29$", result$cell_id), "predefined29", "all30")
  result$processing <- ifelse(startsWith(result$cell_id, "adapter_"), "adapter", "raw")
  result$required_n <- ifelse(result$reference_space == "six_group", 10001L,
    ifelse(result$cohort == "predefined29" & result$reference_space == "D4_vs_control_exact", 126L, 252L))
  result$kind <- ifelse(result$reference_space == "six_group", "sampled", "exhaustive")
  result$family_n <- ifelse(result$kind == "sampled", 67L, 9L)
  result$key <- paste(result$cell_id, result$reference_space, sep = "/")
  result
}

g128f_helper_pins <- function() {
  c(archive = "ff4784ea1dc7da2418638eed49031788783c9ae01da712633d5bfbaf9bc4ad24",
    independent = "e83fafb8580b15f200746527da81d232c36709695dd62ee31a72fab8e1f42c09",
    lineage = "5e48ebca6bd63600ff7ee233fd7e63a320fc8a4d37ab06d143a82bc7b347b33c",
    mc = "c27d12f9a3fa26be4e3a7ff282ef085df05ff05cc191b67b921413b0188c3cc8")
}

g128f_load_helpers <- function(source_records) {
  required <- g128f_helper_pins()
  g128f_need(is.list(source_records) && setequal(names(source_records), names(required)),
    "Supply precisely archive, independent, lineage and selected mc source pins.")
  environments <- lapply(names(required), function(role) {
    record <- g128f_verify(source_records[[role]])
    g128f_need(identical(record$sha256, unname(required[[role]])), "Frozen helper source differs.")
    target <- new.env(parent = baseenv())
    sys.source(record$path, envir = target)
    target
  })
  names(environments) <- names(required)
  environments
}

g128f_bind_shard_receipt <- function(receipt, task, index, manifest, source_identity) {
  shard <- receipt$shard
  g128f_need(is.list(receipt) && is.list(shard) &&
    identical(receipt$job_sha256, task$job_sha256) && identical(receipt$sources, source_identity) &&
    identical(receipt$actual_process_exit_confirmed_here, FALSE) && identical(receipt$full_run_complete, FALSE) &&
    identical(shard$status, "REQUESTED_SHARD_ARCHIVES_WRITTEN_NOT_FULL_RUN_ACCEPTANCE") &&
    identical(shard$run_id, index$fixed$run_id) && identical(shard$cell_id, index$spec$cell_id) &&
    identical(shard$reference_space, index$spec$reference_space) &&
    identical(as.integer(shard$requested_ids), as.integer(manifest$allocation_id)) &&
    shard$total_required_ids == index$spec$required_n &&
    identical(shard$independent_archive_check_complete, FALSE), "Shard receipt job/domain/source differs.")
  original <- shard$files
  for (field in c("file", "sha256", "allocation_row_sha256", "reused_checkpoint")) {
    g128f_need(identical(original[[field]], manifest[[field]]), "Shard manifest/receipt differs.")
  }
  g128f_need(identical(as.integer(original$allocation_id), as.integer(manifest$allocation_id)) &&
    identical(as.numeric(original$bytes), as.numeric(manifest$bytes)) &&
    identical(receipt$additional_execution_sources$transport_source$sha256,
      "622f09f17fe23b6cc2c0d933672b431376a09998a47de8cdb6d2e83550acd753") &&
    identical(receipt$additional_execution_sources$entry_source$sha256,
      "aa5739825a03abf2c47b718e227ba954764df73c340158e88f8616c358dc4a21") &&
    identical(receipt$additional_execution_sources$original_source_manifest_sha256,
      index$fixed$source_manifest_sha256), "Shard additional transport/entry source differs.")
  invisible(TRUE)
}

g128f_require_terminal <- function(terminal) {
  g128f_need(identical(terminal$status, "ALL_REQUESTED_SHARDS_EXITED_ZERO_NEEDS_INDEPENDENT_CHECK") &&
    terminal$planned_shards == 918L && terminal$started_shards == 918L && terminal$finished_shards == 918L &&
    identical(terminal$whole_run_independently_accepted, FALSE),
    "All 918 current shards must actually terminate successfully before finalization.")
  invisible(TRUE)
}

g128f_preflight <- function(root, contract) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  g128f_need(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01") &&
    identical(contract$schema_version, "gse128623-finalization-contract-v1") &&
    identical(contract$run_root, root), "This CLI is bound to the one existing full run.")
  locations <- c(prepared = "PREPARED.rds", input_manifest = "INPUT_MANIFEST.tsv",
    jobs = "JOBS_PREPARED.json", config = "DISPATCH_CONFIG_v2.json",
    admission = "ADMISSION.json", plan = "dispatch/main_02/PLAN.json",
    terminal = "dispatch/main_02/DISPATCH_EXIT.json",
    recovery = "saved_sentinel_recovery_02/RECOVERY_RECEIPT.json")
  g128f_need(setequal(names(contract$metadata), names(locations)), "Eight explicit metadata pins required.")
  records <- lapply(names(locations), function(role) {
    record <- contract$metadata[[role]]
    g128f_need(identical(record$path, file.path(root, locations[[role]])), "Unexpected metadata path.")
    g128f_verify(record)
  })
  names(records) <- names(locations)
  fixed <- c(prepared = "bb58eee81848b9bcac46ba36f57dffd4722e268ddfb52ef5b3a09ba1176b9a65",
    input_manifest = "1d87d05c46bf18a80cb04636e369b3888847b205a9f81c030af3957193a7cc56",
    jobs = "e07be471f4cadd9d515f9f44ef71475e780e96b0eb79dcb75a95db810b7c8750",
    config = "847a049b410828b0b920a9f355ef93e34d9d51d30b68cb4559bf1059f20d44a8",
    admission = "2fe4115f6b496c5a00f47806bd3dc8270a05f993d8fdb1d60251fa746f429307",
    recovery = "f1a920deb0872388b2e336ffced942cf81988f7f8a7f1715772f52189b827303")
  for (role in names(fixed)) g128f_need(identical(records[[role]]$sha256, unname(fixed[[role]])),
    paste("Frozen run metadata differs:", role))
  config <- g128f_json(records$config$path)
  jobs <- g128f_json(records$jobs$path)
  plan <- g128f_json(records$plan$path)
  terminal <- g128f_json(records$terminal$path)
  g128f_require_terminal(terminal)
  g128f_need(identical(plan$config_sha256, records$config$sha256) &&
    identical(plan$controller_sha256, "ba55d2ec1d7aa742ccf0c310802271d8a4e72d224db235c90a75de2db5434a40") &&
    plan$required_archives == 44792L && length(plan$tasks) == 918L &&
    identical(config$jobs_prepared_sha256, records$jobs$sha256), "Controller plan/config binding differs.")
  admission <- g128f_json(records$admission$path)
  g128f_need(identical(admission$status, "ALL_FOUR_OBSERVED_AND_RESOURCE_PASS") &&
    isTRUE(admission$resource_acceptance) && identical(admission$run_id, basename(root)) &&
    identical(jobs$status, "JOBS_PREPARED_NOT_LAUNCHED") &&
    identical(jobs$admission$sha256, records$admission$sha256), "Original job admission differs.")
  # Only the preparation metadata RDS is decoded; job/count/model RDS are not read.
  prepared <- readRDS(records$prepared$path)
  g128f_need(identical(prepared$status, "INPUTS_FROZEN_NO_MODELS") &&
    identical(prepared$run_root, root) &&
    identical(prepared$input_manifest_sha256, records$input_manifest$sha256) &&
    identical(prepared$source_manifest_sha256, "77578bfd27acca61e79ea2290242ee6ac4f1215bc10fc3d54f6143d55752622e") &&
    identical(digest::digest(prepared$source_identity, algo = "sha256"), prepared$source_manifest_sha256),
    "Prepared input/source identity differs.")
  for (file in names(config$source_pins)) {
    record <- g128f_pin(file.path(root, "source", file))
    g128f_need(identical(record$sha256, config$source_pins[[file]]), "Deployed execution source drift.")
    records[[paste0("source/", file)]] <- record
  }
  registry_path <- file.path(root, "registry/PREPARED.json")
  g128f_need(identical(g128f_sha(registry_path),
    "63be688e2d548b2917479a9ecc6e1dba4b3ab37542ac8c4aca9244e78c3ef5d0"), "Frozen registry differs.")
  records$registry <- g128f_pin(registry_path)
  registry <- g128f_json(registry_path)
  specification <- g128f_specification()
  g128f_need(setequal(names(jobs$jobs), specification$key) && sum(specification$required_n) == 44792L,
    "Full fixed 24-job domain required.")
  indices <- vector("list", nrow(specification))
  names(indices) <- specification$key
  for (row in seq_len(nrow(specification))) {
    spec <- specification[row, ]
    job <- jobs$jobs[[spec$key]]
    g128f_need(identical(job$path, file.path(root, "jobs", paste0(spec$cell_id, "_", spec$reference_space, ".rds"))) &&
      file.exists(job$path) && file.info(job$path)$size == job$bytes, "Pinned job path/size differs.")
    cohort <- registry$cohorts[[spec$cohort]]
    item <- if (spec$kind == "sampled") cohort$sampled_space else
      cohort$exact_day_control_spaces[[sub("_vs_control_exact$", "", spec$reference_space)]]
    path <- file.path(root, "registry", item$file)
    g128f_need(identical(g128f_sha(path), item$sha256), "Frozen allocation registry byte identity differs.")
    records[[paste0("registry/", spec$cohort, "/", spec$reference_space)]] <- g128f_pin(path)
    table <- utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
    sample_ids <- as.character(unlist(cohort$sample_ids, use.names = FALSE))
    g128f_need(nrow(table) == spec$required_n && identical(table$allocation_id, as.character(seq_len(spec$required_n))) &&
      all(sample_ids %in% names(table)) && length(sample_ids) == if (spec$cohort == "all30") 30L else 29L,
      "Allocation IDs/sample domain differs.")
    labels <- as.matrix(table[, sample_ids, drop = FALSE])
    row_hashes <- apply(labels, 1L, function(values) digest::digest(
      paste(c(paste(sample_ids, collapse = "\t"), paste(unname(values), collapse = "\t")), collapse = "\n"),
      algo = "sha256", serialize = FALSE))
    g128f_need(!anyNA(labels) && !anyDuplicated(row_hashes), "Missing/repeated registry label row.")
    expected <- data.frame(file = file.path(root, "allocations", spec$cell_id, spec$reference_space,
      sprintf("allocation_%05d.rds", seq_len(spec$required_n))), sha256 = NA_character_,
      allocation_id = seq_len(spec$required_n), allocation_row_sha256 = unname(row_hashes),
      bytes = NA_real_, reused_checkpoint = NA, creator_worker = NA_character_, stringsAsFactors = FALSE)
    indices[[spec$key]] <- list(spec = spec, expected = expected, labels = labels,
      sample_ids = sample_ids, fixed = list(run_id = basename(root), cell_id = spec$cell_id,
        reference_space = spec$reference_space, input_manifest_sha256 = prepared$input_manifest_sha256,
        source_manifest_sha256 = prepared$source_manifest_sha256), observed = admission$observed_receipts[[spec$cell_id]])
  }
  g128f_need(!anyDuplicated(vapply(plan$tasks, function(task) task$worker, character(1L))), "Repeated shard worker ID.")
  for (task in plan$tasks) {
    index <- indices[[task$key]]
    g128f_need(!is.null(index), "Unexpected planned cell/space.")
    ids <- seq.int(as.integer(task$first), as.integer(task$last))
    g128f_need(identical(task$worker, paste0("main_02_", index$spec$cell_id, "_", index$spec$reference_space,
      "_", sprintf("%05d", min(ids)), "_", sprintf("%05d", max(ids)))) &&
      length(ids) <= 50L && min(ids) >= 1L && max(ids) <= index$spec$required_n &&
      all(is.na(index$expected$sha256[ids])) &&
      identical(task$job, jobs$jobs[[task$key]]$path) && identical(task$job_sha256, jobs$jobs[[task$key]]$sha256),
      "Shard task identity/range/job differs or overlaps.")
    start_path <- file.path(root, "dispatch/main_02", task$worker, "STARTED.json")
    exit_path <- file.path(dirname(start_path), "PROCESS_EXIT.json")
    start <- g128f_json(start_path)
    ended <- g128f_json(exit_path)
    command <- as.character(unlist(start$command, use.names = FALSE))
    g128f_need(identical(start$task, task) && ended$exit_code == 0L && length(command) == 9L &&
      identical(command[1:3], c("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript", "--vanilla",
        file.path(root, "source/run_GSE128623_reference_shard_transport_v2_20261015.R"))) &&
      identical(command[4:9], c(file.path(root, "source"), task$job, task$job_sha256,
        as.character(task$first), as.character(task$last), task$worker)), "Actual shard terminal/entry binding differs.")
    path <- file.path(root, "workers", task$worker, "ARCHIVE_MANIFEST.tsv")
    manifest <- g128f_table(path)
    g128f_need(identical(as.integer(manifest$allocation_id), ids) &&
      identical(manifest$file, index$expected$file[ids]) &&
      identical(manifest$allocation_row_sha256, index$expected$allocation_row_sha256[ids]) &&
      is.logical(manifest$reused_checkpoint) && !anyNA(manifest$reused_checkpoint) &&
      all(grepl("^[a-f0-9]{64}$", manifest$sha256)) &&
      all(file.exists(manifest$file)) && all(file.info(manifest$file)$size == manifest$bytes),
      "Required final archive path/coverage/registry/size differs.")
    shard_path <- file.path(dirname(path), "SHARD_RECEIPT.rds")
    g128f_bind_shard_receipt(readRDS(shard_path), task, index, manifest, prepared$source_identity)
    index$expected[ids, c("sha256", "bytes", "reused_checkpoint")] <-
      manifest[, c("sha256", "bytes", "reused_checkpoint")]
    index$expected$creator_worker[ids] <- task$worker
    indices[[task$key]] <- index
    for (meta in c(start_path, exit_path, path, shard_path)) records[[meta]] <- g128f_pin(meta)
  }
  g128f_need(all(vapply(indices, function(index) !anyNA(index$expected), logical(1L))),
    "The complete 44,792 archive index is not ready; no partial scientific output.")
  # Authenticate external creation metadata before any scientific output exists.
  # Old model payloads and saved partials are not decoded in this preflight.
  recovery <- g128f_json(records$recovery$path)
  g128f_need(recovery$model_calls == 0L &&
    identical(recovery$status, "FOUR_SAVED_ID2_SENTINELS_RECOVERED_ZERO_MODEL_CALLS"), "Zero-fit recovery differs.")
  for (cell in unique(specification$cell_id)) {
    observed <- g128f_verify(admission$observed_receipts[[cell]])
    receipt <- g128f_json(observed$path)
    g128f_need(identical(receipt$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED") &&
      file.exists(receipt$observed_bundle$path) &&
      file.info(receipt$observed_bundle$path)$size == receipt$observed_bundle$bytes,
      "Accepted observed cache is not ready.")
    records[[observed$path]] <- observed
    full <- indices[[paste(cell, "six_group", sep = "/")]]$expected
    saved <- recovery$records[[cell]]
    partial <- file.path(root, saved$original_partial$relative_path)
    g128f_need(full$reused_checkpoint[2L] && identical(saved$final$file, full$file[2L]) &&
      identical(saved$final$sha256, full$sha256[2L]) && saved$model_calls == 0L &&
      isTRUE(saved$original_partial_preserved) && file.exists(partial) &&
      file.info(partial)$size == saved$original_partial$bytes, "Recovered saved ID2 metadata not ready.")
    worker <- paste0("fresh_transport_", cell, "_00101")
    path <- file.path(root, "workers", worker, "ARCHIVE_MANIFEST.tsv")
    fresh <- g128f_table(path)
    exit_path <- file.path(root, "logs", paste0("fresh_transport_", cell, "_101_01.resources.json"))
    declared <- Filter(function(record) identical(record$path, exit_path), config$fresh_fixed_sentinel_101_processes)
    g128f_need(full$reused_checkpoint[101L] && nrow(fresh) == 1L && fresh$allocation_id == 101L &&
      !fresh$reused_checkpoint && identical(fresh$file, full$file[101L]) &&
      identical(fresh$sha256, full$sha256[101L]) && length(declared) == 1L &&
      identical(g128f_verify(declared[[1L]])$path, exit_path) && g128f_json(exit_path)$exit_code == 0L,
      "Measured fixed ID101 first-creation metadata not ready.")
    for (meta in c(path, exit_path)) records[[meta]] <- g128f_pin(meta)
  }
  for (index in indices) {
    rows <- which(index$expected$reused_checkpoint)
    if (index$spec$kind == "sampled") rows <- setdiff(rows, c(1L, 2L, 101L))
    for (first in unique(((rows - 1L) %/% 50L) * 50L + 1L)) {
      last <- min(first + 49L, index$spec$required_n)
      worker <- paste(index$spec$cell_id, index$spec$reference_space,
        sprintf("%05d", first), sprintf("%05d", last), sep = "_")
      path <- file.path(root, "workers", worker, "ARCHIVE_MANIFEST.tsv")
      prior <- g128f_table(path)
      selected <- rows[rows >= first & rows <= last]
      at <- match(selected, prior$allocation_id)
      exit_path <- file.path(root, "dispatch/main_01", worker, "PROCESS_EXIT.json")
      start_path <- file.path(dirname(exit_path), "STARTED.json")
      start <- g128f_json(start_path)
      command <- as.character(unlist(start$command, use.names = FALSE))
      g128f_need(!anyNA(at) && !any(prior$reused_checkpoint[at]) &&
        identical(prior$file[at], index$expected$file[selected]) &&
        identical(prior$sha256[at], index$expected$sha256[selected]) &&
        identical(prior$allocation_row_sha256[at], index$expected$allocation_row_sha256[selected]) &&
        g128f_json(exit_path)$exit_code == 0L && identical(start$task$key, index$spec$key) &&
        identical(start$task$job_sha256, jobs$jobs[[index$spec$key]]$sha256) &&
        identical(start$task$worker, worker) && start$task$first == first && start$task$last == last &&
        length(command) == 9L && identical(command[3L],
          file.path(root, "source/run_GSE128623_reference_shard_json_compat_20261015.R")) &&
        identical(command[5:9], c(jobs$jobs[[index$spec$key]]$path, jobs$jobs[[index$spec$key]]$sha256,
          as.character(first), as.character(last), worker)), "Reused legacy first-creation metadata not ready.")
      for (meta in c(path, exit_path, start_path)) records[[meta]] <- g128f_pin(meta)
    }
  }
  list(root = root, jobs = indices, records = records, contract = contract,
    prepared = prepared, config = config, recovery = recovery,
    synthetic = FALSE, all_required_terminal_and_file_metadata_ready = TRUE)
}

g128f_transport_expectation <- function(index, row, payload, prepared) {
  root <- prepared$root
  id <- as.integer(row$allocation_id)
  expected <- c(index$fixed, list(allocation_id = id, allocation_row_sha256 = row$allocation_row_sha256))
  expected <- expected[c("run_id", "cell_id", "reference_space", "allocation_id",
    "input_manifest_sha256", "source_manifest_sha256", "allocation_row_sha256")]
  pins <- list(transport_source_sha256 = "622f09f17fe23b6cc2c0d933672b431376a09998a47de8cdb6d2e83550acd753",
    fresh_entry_sha256 = "aa5739825a03abf2c47b718e227ba954764df73c340158e88f8616c358dc4a21",
    recovery_entry_sha256 = "772641571d8260a9d75765abd95c04a8ebee6c1d968e9e36a08157ce1c2c309b",
    original_source_manifest_sha256 = index$fixed$source_manifest_sha256)
  extra_pins <- list()
  expected_recovery <- NULL
  full <- index$spec$kind == "sampled"
  mode <- "transport_fresh"
  origin <- "CURRENT_MAIN02_FIRST_CREATION_ACTUAL_EXIT_ZERO"
  if (isTRUE(row$reused_checkpoint)) {
    if (full && id == 2L) {
      mode <- "transport_recovered"
      origin <- "PINNED_ZERO_FIT_RECOVERY_AND_AUTHENTICATED_OLD_PARTIAL"
      record <- prepared$recovery$records[[index$spec$cell_id]]
      g128f_need(record$model_calls == 0L && isTRUE(record$original_partial_preserved) &&
        identical(record$final$file, row$file) && identical(record$final$sha256, row$sha256), "Recovered sentinel final binding differs.")
      path <- file.path(root, record$original_partial$relative_path)
      old <- prepared$helpers$independent$g128i_authenticated_rds(path, record$original_partial$sha256)$payload
      g128f_need(identical(old$identity[names(expected)], expected) &&
        identical(digest::digest(old$packed, algo = "sha256"), record$decoded_packed_sha256), "Original saved sentinel differs.")
      expected_recovery <- list(path = path, bytes = record$original_partial$bytes,
        sha256 = record$original_partial$sha256, old_fresh_packed_sha256 = old$identity$packed_sha256,
        old_fresh_hash_not_reasserted_as_decoded_hash = TRUE, model_calls = 0L,
        original_partial_preserved = TRUE, decoded_packed_sha256 = record$decoded_packed_sha256)
      extra_pins[[path]] <- g128f_pin(path)
    } else if (full && id == 101L) {
      origin <- "PINNED_FRESH_ID101_MEASURED_EXIT_ZERO"
      worker <- paste0("fresh_transport_", index$spec$cell_id, "_00101")
      path <- file.path(root, "workers", worker, "ARCHIVE_MANIFEST.tsv")
      prior <- g128f_table(path)
      g128f_need(nrow(prior) == 1L && prior$allocation_id == id && !prior$reused_checkpoint &&
        identical(prior$sha256, row$sha256) && identical(prior$file, row$file), "Fresh fixed sentinel manifest differs.")
      exit_path <- file.path(root, "logs", paste0("fresh_transport_", index$spec$cell_id, "_101_01.resources.json"))
      declared <- Filter(function(record) identical(record$path, exit_path), prepared$config$fresh_fixed_sentinel_101_processes)
      g128f_need(length(declared) == 1L && identical(g128f_verify(declared[[1L]])$path, exit_path) &&
        g128f_json(exit_path)$exit_code == 0L, "Fresh fixed sentinel actual exit differs.")
      for (meta in c(path, exit_path)) extra_pins[[meta]] <- g128f_pin(meta)
    } else if (full && id == 1L) {
      mode <- "legacy"
      origin <- "PRESERVED_OBSERVED_NO_FIT_VIEW_MAIN02_REUSE_AND_ACCEPTED_OBSERVED_BUNDLE"
      # Its first enclosing main_01 shard later failed at ID2. Do not rewrite that
      # failure as exit0: complete current reuse + exact observed repack binds ID1.
    } else {
      mode <- "legacy"
      origin <- "ORIGINAL_MAIN01_SUCCESSFUL_FIRST_CREATION"
      first <- ((id - 1L) %/% 50L) * 50L + 1L
      last <- min(first + 49L, index$spec$required_n)
      worker <- paste(index$spec$cell_id, index$spec$reference_space,
        sprintf("%05d", first), sprintf("%05d", last), sep = "_")
      path <- file.path(root, "workers", worker, "ARCHIVE_MANIFEST.tsv")
      prior <- g128f_table(path)
      prior <- prior[prior$allocation_id == id, , drop = FALSE]
      exit_path <- file.path(root, "dispatch/main_01", worker, "PROCESS_EXIT.json")
      g128f_need(nrow(prior) == 1L && !prior$reused_checkpoint &&
        identical(prior$file, row$file) && identical(prior$sha256, row$sha256) &&
        identical(prior$allocation_row_sha256, row$allocation_row_sha256) &&
        g128f_json(exit_path)$exit_code == 0L, "Original successful legacy creation differs.")
      for (meta in c(path, exit_path)) extra_pins[[meta]] <- g128f_pin(meta)
    }
  }
  checked <- prepared$helpers$lineage$g128l_check(payload, expected, mode, pins, expected_recovery)
  list(row = data.frame(allocation_id = id, mode = mode, origin = origin,
    lineage_status = checked$status, stringsAsFactors = FALSE), records = extra_pins)
}

g128f_one_job <- function(index, observed_pack, root, helpers, inspect_payload) {
  expected <- index$expected
  total <- index$spec$required_n
  g128f_need(nrow(expected) == total && identical(as.integer(expected$allocation_id), seq_len(total)),
    "Full expected ordered IDs required.")
  state <- helpers$archive$g128r_new_reference(observed_pack, total - 1L, index$spec$kind)
  inspections <- vector("list", total)
  for (id in seq_len(total)) {
    payload <- helpers$independent$g128i_authenticated_rds(expected$file[id], expected$sha256[id])$payload
    inspections[[id]] <- inspect_payload(payload, expected[id, , drop = FALSE], id)
    if (id == 1L) g128f_need(identical(payload$packed, observed_pack), "Observed pack differs.") else {
      state <- helpers$archive$g128r_add_reference(state, payload$packed, id)
    }
  }
  producer <- helpers$archive$g128r_finish_reference(state)
  target <- if (index$spec$kind == "sampled") NULL else sub("_exact$", "", index$spec$reference_space)
  independent <- helpers$independent$g128i_check_archives(expected, index$fixed, observed_pack,
    seq_len(total), index$spec$kind, root, target_contrast = target)
  g128f_need(identical(names(producer), names(independent$tail)) && nrow(producer) == nrow(independent$tail),
    "Producer/independent tail table structure differs.")
  for (field in names(producer)) {
    left <- producer[[field]]
    right <- independent$tail[[field]]
    if (is.numeric(left)) {
      g128f_need(identical(is.na(left), is.na(right)) && identical(is.nan(left), is.nan(right)) &&
        all(abs(left[!is.na(left)] - right[!is.na(right)]) <= 1e-12 * (1 + abs(right[!is.na(right)]))),
        paste("Producer/independent arithmetic differs:", field))
    } else g128f_need(identical(left, right), paste("Producer/independent key/status differs:", field))
  }
  g128f_need(independent$checked_allocation_n == total && independent$checked_family_n == index$spec$family_n,
    "Independent complete family/ID domain differs.")
  list(producer_tail = producer, independent = independent, inspections = inspections)
}

g128f_write_new <- function(object, path, rds = TRUE) {
  g128f_need(!file.exists(path), "Output exists; no overwrite.")
  if (rds) saveRDS(object, path, version = 3L, compress = "gzip") else {
    utils::write.table(object, path, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA")
  }
}

g128f_run <- function(admitted, output) {
  g128f_need(isTRUE(admitted$all_required_terminal_and_file_metadata_ready) &&
    !isTRUE(admitted$synthetic), "Only the actual full preflight result can enter run finalization.")
  output <- gsub("\\", "/", output, fixed = TRUE)
  parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  g128f_need(!file.exists(output) &&
    (startsWith(parent, paste0(admitted$root, "/")) || identical(parent, admitted$root) ||
      identical(parent, "R:/Shaoyao_storage/cache")) &&
    grepl("^[A-Za-z0-9_-]+$", basename(output)), "Use an explicit new cloud/R output directory.")
  admitted$helpers <- g128f_load_helpers(admitted$contract$sources)
  sources <- lapply(admitted$contract$sources, g128f_verify)
  # No output is created until every actual terminal and every final filename,
  # range, registry row and byte-size precondition has passed.
  g128f_need(dir.create(output), "Could not create new output.")
  previous_options <- options(digits = 17)
  on.exit(options(previous_options), add = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(output, "SESSION_INFO.txt"))
  summaries <- list()
  output_records <- list()
  evidence <- list()
  for (cell in unique(vapply(admitted$jobs, function(index) index$spec$cell_id, character(1L)))) {
    selected <- Filter(function(index) identical(index$spec$cell_id, cell), admitted$jobs)
    observed_pin <- g128f_verify(selected[[1L]]$observed)
    observed_receipt <- g128f_json(observed_pin$path)
    g128f_need(identical(observed_receipt$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED"),
      "Accepted observed receipt differs.")
    bundle_pin <- g128f_verify(observed_receipt$observed_bundle)
    bundle <- readRDS(bundle_pin$path)
    admitted$records[[observed_pin$path]] <- observed_pin
    admitted$records[[bundle_pin$path]] <- bundle_pin
    for (index in selected) {
      target <- if (index$spec$kind == "sampled") NULL else sub("_exact$", "", index$spec$reference_space)
      cached <- admitted$helpers$archive$g128r_pack_bundle(bundle, target_contrast = target)
      first <- admitted$helpers$independent$g128i_authenticated_rds(index$expected$file[1L], index$expected$sha256[1L])$payload
      for (field in c("keys", "families", "statistic", "native_p", "native_diagnostics", "mroast_auxiliary", "rng_calls", "scope")) {
        g128f_need(identical(first$packed[[field]], cached[[field]]), "Preserved observed no-fit repack differs.")
      }
      inspections <- function(payload, row, id) {
        diagnostic <- payload$packed$worker_diagnostics
        g128f_need(is.list(diagnostic) && identical(as.character(diagnostic$labels), unname(index$labels[id, ])) &&
          identical(names(diagnostic$labels), index$sample_ids) && diagnostic$allocation_id == id &&
          identical(diagnostic$cohort, index$spec$cohort) && identical(diagnostic$processing, index$spec$processing) &&
          identical(diagnostic$observed_fit_reused, id == 1L), "Compact labels/cohort/observed identity differs.")
        retained <- !is.null(target) && id == 1L && length(payload$packed$rng_calls) > 0L
        expected_scope <- if (retained) "RETAINED_OBSERVED_CHRONIC_TRACE_NOT_DAY_FIT" else
          if (!is.null(target)) "NO_UNSIGNED_DAY_FIT_OR_TRACE" else "NATIVE_CHRONIC_REFERENCE_TRACE"
        g128f_need(identical(diagnostic$unsigned_trace_scope, expected_scope), "Unsigned trace scope differs.")
        native <- payload$packed$native_evidence
        sentinel <- is.null(target) && id %in% c(2L, 101L, 1001L, 5001L, 10001L)
        kind <- if (id == 1L) "preserved_observed_bundle" else if (sentinel)
          "fixed_six_group_native_sentinel" else "compact_statistics_P_and_diagnostics_only"
        g128f_need(is.list(native) && identical(native$kind, kind), "Fixed native evidence kind differs.")
        if (id == 1L) g128f_need(identical(native$reference$sha256, bundle_pin$sha256) &&
          identical(native$reference$file, bundle_pin$path), "Observed native source differs.")
        if (sentinel) g128f_need(is.list(native$bundle) &&
          all(c("gene_camera", "activity", "unsigned", "integration_diagnostics") %in% names(native$bundle)),
          "Prespecified full native sentinel bundle missing.")
        lineage <- g128f_transport_expectation(index, row, payload, admitted)
        warning_n <- length(diagnostic$conditions$warning)
        message_n <- length(diagnostic$conditions$message)
        list(lineage = lineage$row, records = lineage$records,
          compact = data.frame(allocation_id = id, warning_n = warning_n, message_n = message_n,
            native_evidence_kind = kind, full_native_fields_independently_accepted = FALSE,
            stringsAsFactors = FALSE))
      }
      result <- g128f_one_job(index, first$packed, admitted$root, admitted$helpers, inspections)
      leaf <- file.path(output, paste(index$spec$cell_id, index$spec$reference_space, sep = "__"))
      g128f_need(dir.create(leaf), "Job output exists.")
      g128f_write_new(result$producer_tail, file.path(leaf, "PRODUCER_TAIL.rds"))
      g128f_write_new(result$independent, file.path(leaf, "INDEPENDENT_ARITHMETIC.rds"))
      lineage <- do.call(rbind, lapply(result$inspections, `[[`, "lineage"))
      compact <- do.call(rbind, lapply(result$inspections, `[[`, "compact"))
      g128f_write_new(lineage, file.path(leaf, "TRANSPORT_LINEAGE.tsv"), FALSE)
      g128f_write_new(compact, file.path(leaf, "COMPACT_DIAGNOSTICS_SCOPE.tsv"), FALSE)
      for (inspection in result$inspections) admitted$records <- c(admitted$records, inspection$records)
      tail <- result$independent$tail
      native <- result$independent$native_diagnostics
      families <- result$independent$empirical_family_diagnostics
      families$cell_id <- index$spec$cell_id
      families$reference_space <- index$spec$reference_space
      families$required_allocation_n <- index$spec$required_n
      families$empirical_unavailable_n <- vapply(families$family_id, function(key)
        sum(is.na(tail$empirical_p[tail$family_id == key])), integer(1L))
      families$empirical_q_le_05_n <- vapply(families$family_id, function(key)
        sum(tail$empirical_BH_q[tail$family_id == key] <= .05, na.rm = TRUE), integer(1L))
      families$native_fraction_allocations_any_q_le_05 <- vapply(families$family_id, function(key)
        mean(native$native_q_any_le_05[native$family_id == key]), numeric(1L))
      summaries[[index$spec$key]] <- families
      if (index$spec$kind == "sampled") {
        mc <- admitted$helpers$mc$g128mc_summarize(tail,
          if (index$spec$cohort == "all30") "all30" else "pre29")
        g128f_write_new(mc, file.path(leaf, "SELECTED_OUTER_MC_SIDECAR.rds"))
      }
      paths <- list.files(leaf, full.names = TRUE, recursive = FALSE)
      output_records <- c(output_records, lapply(paths, g128f_pin))
      evidence[[index$spec$key]] <- list(checked_allocations = index$spec$required_n,
        family_n = index$spec$family_n, arithmetic_and_transport = TRUE,
        exact_companion_outer_MC = if (index$spec$kind == "exhaustive") "NOT_APPLICABLE" else "SELECTED_SIDECAR",
        full_native_bundle_fields_and_RNG_consumption_accepted = FALSE)
      rm(result, first, cached)
      invisible(gc())
    }
    rm(bundle)
    invisible(gc())
  }
  # Metadata/source drift is rechecked; archive bytes were authenticated during
  # both arithmetic passes, not falsely described as a third full-disk rehash.
  for (record in c(admitted$records, sources)) g128f_verify(record)
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  g128f_need(nrow(summary) == 448L &&
    sum(summary$reference_space == "six_group") == 268L, "Final family inventory differs.")
  g128f_write_new(summary, file.path(output, "FAMILY_REFERENCE_SUMMARY.tsv"), FALSE)
  receipt <- list(status = "FULL_24_JOB_ARITHMETIC_AND_TRANSPORT_COMPLETE_NATIVE_ACCEPTANCE_PENDING",
    completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), run_root = admitted$root,
    jobs = evidence, job_n = 24L, required_archive_n = 44792L, native_six_group_family_n = 268L,
    companion_family_slots = 180L, producer_independent_tails_and_BH_agree = TRUE,
    all_required_metadata_and_sources_unchanged = TRUE, archive_bytes_authenticated_each_arithmetic_pass = TRUE,
    archives_rehashed_again_after_all_jobs = FALSE, raw_count_expression_and_job_payloads_rehashed_here = FALSE,
    source_records = sources,
    metadata_records = admitted$records, output_records = c(output_records, list(g128f_pin(file.path(output, "FAMILY_REFERENCE_SUMMARY.tsv")))),
    scientific_model_calls = 0L, original_results_modified = FALSE, candidate_gate_eligible = FALSE,
    R = R.version.string, library_paths = .libPaths(), command = commandArgs(),
    full_native_evidence_independently_accepted = FALSE, final_scientific_acceptance = FALSE,
    unresolved = c("Full native sentinel/extra model fields and final RNG consumption need separate bounded acceptance",
      "Conditional same-data reference does not authenticate biological exchangeability or true FDR",
      "Pointwise outer MC intervals do not quantify inner mroast/fgsea simulation error"))
  path <- file.path(output, "FINALIZATION_RECEIPT.json")
  g128f_need(!file.exists(path), "Final receipt exists.")
  jsonlite::write_json(receipt, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  invisible(receipt)
}

g128f_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  g128f_need(length(arguments) == 4L, "RUN_ROOT CONTRACT_JSON CONTRACT_SHA256 NEW_OUTPUT")
  g128f_need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE),
    "Use the existing digest/jsonlite environment; do not install here.")
  g128f_need(identical(as.character(getRversion()), "4.6.1"), "Pinned R4.6.1 required.")
  contract_path <- normalizePath(arguments[2L], winslash = "/", mustWork = TRUE)
  g128f_need(identical(g128f_sha(contract_path), arguments[3L]), "Caller contract SHA differs.")
  contract <- g128f_json(contract_path)
  admitted <- g128f_preflight(arguments[1L], contract)
  # Fresh output and source/contract pinning are local to this no-fit step.
  admitted$records$caller_contract <- g128f_pin(contract_path)
  self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  g128f_need(length(self) == 1L, "Exact finalizer source path required.")
  admitted$records$finalizer_source <- g128f_pin(self)
  g128f_run(admitted, arguments[4L])
}

if (sys.nframe() == 0L) g128f_main()
