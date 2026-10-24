# Saved-only ordered merge; source-safe. No fit, labels, RNG, finish-P or BH.
# CLI (15 args): --execute CONFIG ACCESSION four_group_MC FIRST LAST PREFIX
# PREFIX_BYTES PREFIX_SHA LEGACY_ACCEPT LEGACY_SHA SOURCE452 SHA452 PERMIT PERMIT_SHA
jgw11sm_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgw11sm_root <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03"
jgw11sm_old_sha <- "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d"
jgw11sm_context_sha <- "452acfa10dbbaa29631ad41877feab120ba421dec6852bd42ab5180ef9b1b56c"
jgw11sm_helper_sha <- "74bb8dc0d3cdb16bb9373fdf5a9e6153af7298d5c104276e58c3247fc41b972e"
jgw11sm_fit_sha <- "38fd274dcf96af76f3a4056dbd3021d396490da67da7064da8c8fe4c8357cd02"
jgw11sm_whole <- function(x, minimum = 0) is.numeric(x) && length(x) == 1L &&
  !is.na(x) && is.finite(x) && x >= minimum && x == floor(x)
jgw11sm_same_pin <- function(a, b) is.list(a) && is.list(b) &&
  identical(a$path, b$path) && jgw11sm_whole(a$bytes, 1) && jgw11sm_whole(b$bytes, 1) &&
  a$bytes == b$bytes && identical(a$sha256, b$sha256)

jgw11sm_permission <- function(p, z, actual, context, helper, fit, previous, legacy, first, last) {
  jgw11sm_need(is.list(p) && identical(p$schema, "jgw03-eleven-worker-saved-fit-shard-merge-acceptance-v1") &&
    identical(p$status, "ACCEPTED_FOR_ONE_ORDERED_SAVED_SHARD_MERGE") &&
    isTRUE(p$resources_admitted) && isTRUE(p$independent_wiring_review_accepted) &&
    isTRUE(p$previous_boundary_authenticated_by_caller) &&
    isTRUE(p$single_owner_disjoint_partitions_authenticated_by_caller) &&
    identical(p$maximum_parallel_single_core_workers, 11L) &&
    jgw11sm_whole(p$per_worker_memory_bytes, 1) && p$per_worker_memory_bytes == 3*1024^3,
    "SeparateReviewedElevenWorkerSavedMergePermission")
  jgw11sm_need(jgw11sm_same_pin(p$actual_merge, actual) &&
    jgw11sm_same_pin(p$legacy_scientific_entry, z$before$entry) &&
    jgw11sm_same_pin(p$config, z$before$config) && jgw11sm_same_pin(p$legacy_acceptance, legacy) &&
    jgw11sm_same_pin(p$successor_context_source, context) && identical(context$sha256, jgw11sm_context_sha) &&
    jgw11sm_same_pin(p$increment_helper, helper) && identical(helper$sha256, jgw11sm_helper_sha) &&
    jgw11sm_same_pin(p$fit_executor, fit) && identical(fit$sha256, jgw11sm_fit_sha) &&
    identical(p$job$accession, z$accession) && identical(p$job$space, "four_group_MC") &&
    jgw11sm_whole(first, 3) && jgw11sm_whole(last, first) && last <= 10001 && last-first+1 <= 50 &&
    jgw11sm_whole(p$job$first, 3) && jgw11sm_whole(p$job$last, first) &&
    p$job$first == first && p$job$last == last &&
    jgw11sm_same_pin(p$job$previous_prefix, previous), "ExactMergeSourcesCellRangeAndPreviousPrefix")
  jgw11sm_need(is.list(p$evidence) && length(p$evidence) == 3L && !anyDuplicated(names(p$evidence)) &&
    setequal(names(p$evidence), c("review", "resources", "previous_boundary")) &&
    is.list(p$shards) && length(p$shards) %in% c(1L, 2L, 3L, 4L), "ExplicitMergeEvidenceAndShardPins")
  invisible(TRUE)
}

# Every supplied record is small authenticated JSON. Native exit is a direct
# parent observation, never inferred from the R_COMPLETE status string.
jgw11sm_shards <- function(shards, z, b, first, last, previous, fit, parts, record, read_json) {
  jgw11sm_need(length(shards) == length(parts), "ExactFixedShardPartitionCount")
  records <- list()
  for (index in seq_along(parts)) {
    pins <- shards[[index]]
    jgw11sm_need(is.list(pins) && setequal(names(pins), c("r_complete", "native_exit")) &&
      !anyDuplicated(names(pins)), "TwoExplicitShardReceiptPins")
    rp <- record(pins$r_complete); np <- record(pins$native_exit)
    r <- read_json(rp$path); n <- read_json(np$path)
    span <- parts[[index]]
    jgw11sm_need(identical(r$schema, "jgw03-eleven-worker-fit-only-shard-result-v1") &&
      identical(r$status, "FIT_ONLY_SHARD_COMPLETE_NOT_PREFIX_OR_INDEPENDENT_ACCEPTANCE") &&
      identical(r$accession, z$accession) && identical(r$space, "four_group_MC") &&
      jgw11sm_whole(r$first, first) && jgw11sm_whole(r$last, r$first) &&
      r$first == span[["first"]] && r$last == span[["last"]] &&
      jgw11sm_whole(r$round_first, first) && r$round_first == first &&
      jgw11sm_whole(r$round_last, first) && r$round_last == last &&
      jgw11sm_same_pin(r$actual_executor, fit) &&
      jgw11sm_same_pin(r$legacy_scientific_entry, z$before$entry) &&
      identical(r$scientific_identity_sha256, b$identity_sha256) &&
      jgw11sm_same_pin(r$previous_prefix, previous) && identical(r$prefix_written, FALSE) &&
      identical(r$actual_native_exit_observed_here, FALSE), "ExactCompletedFitShardNotPrefix")
    jgw11sm_need(identical(n$schema, "jgw03-eleven-worker-child-native-exit-v1") &&
      identical(n$actual_native_exit_code, 0L) && identical(n$actual_R_exit_code, 0L) &&
      is.null(n$observer_failure_class_only) && isTRUE(n$memory_budget_met) &&
      identical(n$raw_streams_persisted, FALSE) && identical(n$automatic_native_retry, FALSE) &&
      jgw11sm_whole(n$actual_pid, 1) && identical(n$actual_R_pid, n$actual_pid) &&
      jgw11sm_whole(n$start_ticks, 1) && is.list(n$R_argv) && length(n$R_argv) > 0L &&
      jgw11sm_same_pin(n$actual_executor, fit) && jgw11sm_same_pin(n$r_complete, rp),
      "DirectParentStrictZeroAndExactCompletedShardRequired")
    command <- read_json(record(n$command)$path); started <- read_json(record(n$started)$path)
    job <- command$job
    jgw11sm_need(identical(command$entry_sha256, fit$sha256) && identical(command$R_argv, n$R_argv) &&
      identical(n$job, job) &&
      identical(job$mode, "fit_only") && identical(job$space, "four_group_MC") &&
      identical(job$accession, z$accession) && job$first == r$first && job$last == r$last &&
      job$round_first == first && job$round_last == last && job$shard_index == index &&
      identical(started$job, job) && identical(started$entry_sha256, fit$sha256) &&
      identical(started$pid, n$actual_pid) && identical(started$start_ticks, n$start_ticks) &&
      jgw11sm_whole(started$parent, 1),
      "ActualShardParentCommandAndStartedIdentity")
    ids <- seq.int(span[["first"]], span[["last"]])
    jgw11sm_need(is.list(r$allocation_records) && length(r$allocation_records) == length(ids),
      "ExactShardAllocationRecordCount")
    for (position in seq_along(ids)) {
      i <- ids[position]; key <- b$keys$four_group_MC[i]; one <- r$allocation_records[[position]]
      jgw11sm_need(jgw11sm_whole(one$allocation_id, first) && one$allocation_id == i &&
        identical(one$allocation_key, key) && jgw11sm_same_pin(one$execution_entry, fit) &&
        identical(one$prior_leaf_existed_at_preflight, FALSE) &&
        identical(one$new_native_calls_recorded_in_archive, 1L), "OrderedUniqueShardKeyAndHonestFitActor")
      cp <- record(one$completion_metadata)
      jgw11sm_need(identical(cp$path, file.path(b$store, key, "COMPLETE.json")), "CanonicalSmallCompletionMetadata")
      c <- read_json(cp$path)
      jgw11sm_need(identical(c$status, "IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT") &&
        identical(c$allocation_key, key) && identical(c$scientific_identity_sha256, b$identity_sha256) &&
        identical(c$new_native_calls, 1L) && is.list(c$archive) &&
        identical(names(c$archive), c("bytes", "sha256")) && jgw11sm_whole(c$archive$bytes, 1) &&
        is.character(c$archive$sha256) && length(c$archive$sha256) == 1L &&
        !is.na(c$archive$sha256) && grepl("^[a-f0-9]{64}$", c$archive$sha256), "OriginalCompletedArchivePin")
      records[[as.character(i)]] <- list(allocation_id = i, allocation_key = key,
        execution_entry = fit, completion_metadata = cp, archive = c$archive)
    }
  }
  jgw11sm_need(identical(names(records), as.character(seq.int(first, last))), "CompleteOrderedUniqueRoundIDs")
  records
}

jgw11sm_merge <- function(previous, b, first, last, records, fit, guard, advance, reader, update, pin, writer) {
  guard(previous, b, first, last)
  jgw11sm_need(identical(names(records), as.character(seq.int(first, last))), "ExactlyOneRoundBeforeReadback")
  old_identity <- b$identity; old_sha <- b$identity_sha256
  read_once <- function(model, key) {
    id <- match(key, b$keys$four_group_MC)
    jgw11sm_need(!is.na(id) && id >= first && id <= last, "NoPriorOrFutureArchiveRead")
    value <- reader(model, key)
    jgw11sm_need(jgw11sm_same_pin(value$execution_entry, fit) &&
      identical(value$new_native_calls, 1L), "SavedArchiveHasActualFitOnlyExecutor")
    value
  }
  # Frozen helper streams each new key exactly once in frozen registry order.
  bound_pin <- function(path) {
    id <- match(basename(dirname(path)), b$keys$four_group_MC)
    jgw11sm_need(!is.na(id) && id >= first && id <= last && basename(path) == "COMPACT.rds",
      "PinOnlyThisRoundsArchives")
    got <- pin(path); expected <- records[[as.character(id)]]$archive
    jgw11sm_need(identical(names(got), c("bytes", "sha256")) &&
      got$bytes == expected$bytes && identical(got$sha256, expected$sha256), "ArchiveMatchesAuthenticatedCompletion")
    got
  }
  prefix <- advance(previous, b, last, read_once, update, bound_pin)
  jgw11sm_need(identical(b$identity, old_identity) && identical(b$identity_sha256, old_sha), "ScientificIdentityPreserved")
  target <- file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", last))
  list(product = writer(prefix, target), target = target, prefix = prefix)
}

jgw11sm_main <- function(args) {
  jgw11sm_need(length(args) == 15L && identical(args[1L], "--execute") &&
    args[3L] %in% c("GSE289929", "GSE227494", "GSE299063") && identical(args[4L], "four_group_MC") &&
    all(grepl("^[0-9]+$", args[c(5L, 6L, 8L)])) && identical(args[13L], jgw11sm_context_sha), "ExplicitOneSavedMergeJob")
  old_path <- file.path(jgw11sm_root, "project/code/R/run_JGW03_full_reference_linux_v3_20261017.R")
  jgw11sm_need(identical(digest::digest(file = old_path, algo = "sha256"), jgw11sm_old_sha), "FixedLegacyScience")
  old <- new.env(parent = asNamespace("stats")); sys.source(old_path, old)
  actual <- old$jgwl_entry()
  jgw11sm_need(old$jgwl_inside(actual$path, jgw11sm_root) &&
    !old$jgwl_inside(actual$path, file.path(jgw11sm_root, "project")), "HonestSeparateMergeEntry")
  context <- old$jgwl_pin(args[12L])
  jgw11sm_need(identical(context$sha256, jgw11sm_context_sha) && old$jgwl_inside(context$path, jgw11sm_root), "Fixed452Context")
  candidate <- new.env(parent = asNamespace("stats")); sys.source(context$path, candidate)
  loader <- candidate$jgwsp_legacy_scientific_context; environment(loader) <- old
  z <- loader(args[2L], args[3L]); legacy <- old$jgwl_accept(z, args[10L], args[11L])
  previous_pin <- old$jgwl_record(list(path = args[7L], bytes = as.numeric(args[8L]), sha256 = args[9L]))
  permit_pin <- old$jgwl_record(list(path = args[14L], bytes = unname(file.info(args[14L])$size), sha256 = args[15L]))
  p <- jsonlite::read_json(permit_pin$path, simplifyVector = FALSE)
  helper <- old$jgwl_record(p$increment_helper); fit <- old$jgwl_record(p$fit_executor)
  first <- as.numeric(args[5L]); last <- as.numeric(args[6L])
  jgw11sm_permission(p, z, actual, context, helper, fit, previous_pin, legacy, first, last)
  evidence <- lapply(p$evidence, old$jgwl_record)
  jgw11sm_need(all(vapply(c(list(helper, fit), evidence), function(x) old$jgwl_inside(x$path, z$run), logical(1L))), "RunBoundMergeSourcesEvidence")
  producer <- new.env(parent = asNamespace("stats")); sys.source(fit$path, producer)
  h <- new.env(parent = baseenv()); sys.source(helper$path, h)
  b <- old$jgwl_boot(z)$b; previous <- readRDS(previous_pin$path)
  candidate$jgwsp_range_guard(previous, b, first, last)
  jgw11sm_need(identical(previous_pin$path, file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", first-1))), "CanonicalPreviousSameCellPrefix")
  target <- file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", last))
  jgw11sm_need(!file.exists(target) && !file.exists(paste0(target, ".partial")), "FreshPrefixNoAutomaticRetry")
  record <- function(x) { a <- old$jgwl_record(x); jgw11sm_need(old$jgwl_inside(a$path, z$run), "RunBoundShardMetadata"); a }
  read_json <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)
  records <- jgw11sm_shards(p$shards, z, b, first, last, previous_pin, fit,
    producer$jgw11sh_partitions(z$accession, first, last), record, read_json)
  output <- p$output_leaf; parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  jgw11sm_need(old$jgwl_inside(parent, file.path(z$run, "provenance")) &&
    identical(output, file.path(parent, basename(output))) && !basename(output) %in% c(".", "..") &&
    !file.exists(output) && dir.create(output), "ExclusiveMergeReceiptLeaf")
  old$jgwl_json(list(status = "SAVED_SHARD_MERGE_STARTED_NOT_NATIVE_EXIT", actual_merge = actual,
    previous_prefix = previous_pin, executor_acceptance = permit_pin), file.path(output, "STARTED.json"))
  complete <- FALSE
  on.exit(if (!complete) old$jgwl_json(list(status = "STOP_RETAINED_NO_AUTOMATIC_REMERGE",
    actual_native_exit_observed_here = FALSE), file.path(output, "FAILED.json")), add = TRUE)
  result <- jgw11sm_merge(previous, b, first, last, records, fit, candidate$jgwsp_range_guard,
    h$jgw03_advance_saved_tail_prefix, b$e$jgwf_read_archive, b$h$jgw03c_tail_update, b$e$jgwx_pin,
    b$e$jgwx_atomic_plain)
  old$jgwl_after(z)
  invisible(lapply(c(list(actual, context, helper, fit, legacy, previous_pin, permit_pin), evidence,
    unlist(lapply(p$shards, function(x) unname(x)), recursive = FALSE),
    lapply(records, `[[`, "completion_metadata")), old$jgwl_record))
  old$jgwl_json(list(schema = "jgw03-eleven-worker-saved-fit-shard-merge-result-v1",
    status = "ORDERED_SAVED_SHARD_PREFIX_COMPLETE_NOT_INDEPENDENT_ACCEPTANCE", accession = z$accession,
    space = "four_group_MC", first = first, last = last, actual_merge = actual, fit_executor = fit,
    legacy_scientific_entry = z$before$entry, successor_context_source = context,
    scientific_identity_sha256 = b$identity_sha256, previous_prefix = previous_pin,
    new_prefix = c(list(path = result$target), result$product), allocation_records = unname(records),
    executor_acceptance = permit_pin, legacy_acceptance = legacy, evidence = evidence, shards = p$shards,
    new_fits = 0L, allocation_RNG_calls = 0L, new_archives_read = last-first+1,
    prior_archives_read = 0L, original_CLI_success_asserted = FALSE,
    actual_native_exit_observed_here = FALSE, scientific_acceptance = FALSE), file.path(output, "R_COMPLETE.json"))
  complete <- TRUE
  0L
}
if (sys.nframe() == 0L) {
  exit <- tryCatch(jgw11sm_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_SAVED_SHARD_MERGE_STOP_RETAINED_NO_AUTOMATIC_RETRY"); 1L
  })
  quit(status = exit)
}
