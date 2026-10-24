# Source-safe thin producer. No scheduler, fork, prefix write, finish-P or BH.
# CLI: --execute CONFIG ACCESSION four_group_MC FIRST LAST PREVIOUS_PREFIX
#      PREVIOUS_BYTES PREVIOUS_SHA LEGACY_ACCEPT LEGACY_SHA SOURCE452 SHA452
#      SHARD_PERMIT SHARD_PERMIT_SHA
# FIRST/LAST are this shard. The permit binds the whole <=50-ID next round.
# Actual entry is this new file; legacy scientific identity remains unchanged.
jgw11sh_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgw11sh_root <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03"
jgw11sh_old_sha <- "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d"
jgw11sh_context_sha <- "452acfa10dbbaa29631ad41877feab120ba421dec6852bd42ab5180ef9b1b56c"
jgw11sh_merge_sha <- "74bb8dc0d3cdb16bb9373fdf5a9e6153af7298d5c104276e58c3247fc41b972e"
jgw11sh_whole <- function(x, minimum = 0) is.numeric(x) && length(x) == 1L &&
  !is.na(x) && is.finite(x) && x >= minimum && x == floor(x)
jgw11sh_same_pin <- function(a, b) is.list(a) && is.list(b) &&
  identical(a$path, b$path) && jgw11sh_whole(a$bytes, 1) && jgw11sh_whole(b$bytes, 1) &&
  a$bytes == b$bytes && identical(a$sha256, b$sha256)

# Deterministic partition, not new labels/sampling. Last one-ID round has one
# shard, never an empty worker. The external owner must issue disjoint permits.
jgw11sh_partitions <- function(accession, first, last) {
  jgw11sh_need(accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    jgw11sh_whole(first, 3) && jgw11sh_whole(last, first) && last <= 10001 &&
    last - first + 1 <= 50, "OneBoundedFutureRound")
  n <- last - first + 1
  count <- as.integer(min(if (accession == "GSE227494") 3L else 4L, n))
  if (count == 1L) return(list(c(first = first, last = last)))
  sizes <- rep(floor(n / count), count)
  remainder <- n %% count
  if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1
  ends <- first - 1 + cumsum(sizes)
  starts <- c(first, head(ends, -1L) + 1)
  lapply(seq_len(count), function(i) c(first = starts[i], last = ends[i]))
}

jgw11sh_permission <- function(p, z, actual, context, previous, legacy, first, last) {
  jgw11sh_need(is.list(p) && identical(p$schema, "jgw03-eleven-worker-fit-only-shard-acceptance-v1") &&
    identical(p$status, "ACCEPTED_FOR_ONE_UNUSED_FIT_ONLY_SHARD") &&
    isTRUE(p$resources_admitted) && isTRUE(p$independent_wiring_review_accepted) &&
    isTRUE(p$previous_boundary_authenticated_by_caller) &&
    isTRUE(p$single_owner_disjoint_partitions_authenticated_by_caller) &&
    identical(p$maximum_parallel_single_core_workers, 11L) &&
    jgw11sh_whole(p$per_worker_memory_bytes, 1) && p$per_worker_memory_bytes == 3*1024^3,
    "NewElevenSingleCoreWorkerPermissionRequired")
  jgw11sh_need(jgw11sh_same_pin(p$actual_executor, actual) &&
    jgw11sh_same_pin(p$legacy_scientific_entry, z$before$entry) &&
    jgw11sh_same_pin(p$successor_context_source, context) && identical(context$sha256, jgw11sh_context_sha) &&
    jgw11sh_same_pin(p$config, z$before$config) && jgw11sh_same_pin(p$legacy_acceptance, legacy) &&
    jgw11sh_same_pin(p$job$previous_prefix, previous) &&
    identical(p$job$accession, z$accession) && identical(p$job$space, "four_group_MC") &&
    jgw11sh_whole(p$job$first, 3) && jgw11sh_whole(p$job$last, p$job$first) &&
    p$job$first == first && p$job$last == last &&
    jgw11sh_whole(p$job$round_first, 3) && jgw11sh_whole(p$job$round_last, p$job$round_first) &&
    jgw11sh_whole(p$job$shard_index, 1), "ExactShardSourceCellRangePreviousPrefix")
  parts <- jgw11sh_partitions(z$accession, p$job$round_first, p$job$round_last)
  index <- p$job$shard_index
  jgw11sh_need(index <= length(parts) && first == parts[[index]][["first"]] &&
    last == parts[[index]][["last"]], "OnlyFixedNonoverlappingPartition")
  jgw11sh_need(is.list(p$evidence) && length(p$evidence) == 3L && !anyDuplicated(names(p$evidence)) &&
    setequal(names(p$evidence), c("review", "resources", "previous_boundary")), "ThreeExplicitShardEvidencePins")
  invisible(TRUE)
}

# All own keys/labels are preflighted before the first callback. Never call the
# old tail-prefix builder or the 452 range runner. Existing leaves (complete or
# partial, any actor) are retained and rejected here, not refitted or relabelled.
# The original get_draw callback owns each new leaf via its original dir.create.
jgw11sh_fit <- function(z, b, previous, round_first, round_last, first, last,
                      actual, fixed, context_range_guard, leaf_exists, get_draw, pin_complete) {
  context_range_guard(previous, b, round_first, round_last)
  jgw11sh_need(jgw11sh_whole(first, round_first) && jgw11sh_whole(last, first) && last <= round_last &&
    is.function(leaf_exists) && is.function(get_draw) && is.function(pin_complete),
    "OnlyUnusedBoundedShardCallbacks")
  keys <- b$keys$four_group_MC
  jgw11sh_need(is.character(keys) && length(keys) == 10001L && !anyNA(keys) &&
    !anyDuplicated(keys) && all(nzchar(keys)), "CompleteUnchangedFrozenKeys")
  ids <- seq.int(first, last)
  labels <- lapply(ids, function(i) {
    allocation <- stats::setNames(as.numeric(b$registry$MC[i, ]), colnames(b$registry$MC))
    b$modules$kernel$jgw03_check_allocation(allocation, b$domain$sample_ids)
    jgw11sh_need(!leaf_exists(file.path(b$store, keys[i])), "EveryShardKeyMustBeUnusedBeforeAnyFit")
    allocation
  })
  identity <- b$identity
  identity_sha <- b$identity_sha256
  b$entry <- actual
  records <- vector("list", length(ids))
  for (position in seq_along(ids)) {
    i <- ids[position]
    key <- keys[i]
    jgw11sh_need(!leaf_exists(file.path(b$store, key)), "NoRacingOrPreexistingKeyBeforeCallback")
    # Exactly the frozen callback/signature. It resets the original seed and
    # calls the unchanged kernel once; no internal parallel backend is changed.
    draw <- get_draw(b, key, labels[[position]], fixed, z$cell_dir, TRUE)
    jgw11sh_need(is.list(draw) && identical(draw$scientific_identity_sha256, identity_sha) &&
      identical(draw$allocation_key, key) && jgw11sh_same_pin(draw$execution_entry, actual) &&
      identical(draw$new_native_calls, 1L), "OneNewOriginalCallbackArchiveWithHonestActor")
    # Pin only the small original COMPLETE.json, whose archive bytes/SHA were
    # just validated by get_draw. No additional scientific RDS read/hash here.
    completion_path <- file.path(b$store, key, "COMPLETE.json")
    completion <- pin_complete(completion_path)
    jgw11sh_need(is.list(completion) && identical(names(completion), c("bytes", "sha256")) &&
      jgw11sh_whole(completion$bytes, 1) && is.character(completion$sha256) &&
      length(completion$sha256) == 1L && !is.na(completion$sha256) &&
      grepl("^[0-9a-f]{64}$", completion$sha256), "OriginalSmallCompletePinSchema")
    records[[position]] <- list(allocation_id = i, allocation_key = key,
      execution_entry = draw$execution_entry, prior_leaf_existed_at_preflight = FALSE,
      new_native_calls_recorded_in_archive = draw$new_native_calls,
      completion_metadata = c(list(path = completion_path), completion))
    if (i %% 10L == 0L) invisible(gc())
  }
  jgw11sh_need(identical(b$identity, identity) && identical(b$identity_sha256, identity_sha),
    "ScientificIdentityNeverRelabelled")
  list(allocation_records = records, scientific_identity_sha256 = identity_sha,
    actual_executor = actual, prefix_written = FALSE)
}

jgw11sh_main <- function(args) {
  jgw11sh_need(length(args) == 15L && identical(args[1L], "--execute") &&
    args[3L] %in% c("GSE289929", "GSE227494", "GSE299063") &&
    identical(args[4L], "four_group_MC") && all(grepl("^[0-9]+$", args[c(5L, 6L, 8L)])) &&
    identical(args[13L], jgw11sh_context_sha), "ExplicitOneFitOnlyShardArguments")
  old_path <- file.path(jgw11sh_root, "project/code/R/run_JGW03_full_reference_linux_v3_20261017.R")
  jgw11sh_need(identical(digest::digest(file = old_path, algo = "sha256"), jgw11sh_old_sha), "FixedLegacyScientificSource")
  old <- new.env(parent = asNamespace("stats"))
  sys.source(old_path, old)
  actual <- old$jgwl_entry()
  jgw11sh_need(old$jgwl_inside(actual$path, jgw11sh_root) &&
    !old$jgwl_inside(actual$path, file.path(jgw11sh_root, "project")), "HonestNewShardExecutor")
  context <- old$jgwl_pin(args[12L])
  jgw11sh_need(identical(context$sha256, jgw11sh_context_sha) &&
    old$jgwl_inside(context$path, jgw11sh_root), "Unchanged452ContextSource")
  candidate <- new.env(parent = asNamespace("stats"))
  sys.source(context$path, candidate)
  loader <- candidate$jgwsp_legacy_scientific_context
  environment(loader) <- old
  z <- loader(args[2L], args[3L])
  legacy <- old$jgwl_accept(z, args[10L], args[11L])
  previous_pin <- old$jgwl_record(list(path = args[7L], bytes = as.numeric(args[8L]), sha256 = args[9L]))
  permit_pin <- old$jgwl_record(list(path = args[14L], bytes = unname(file.info(args[14L])$size), sha256 = args[15L]))
  permit <- jsonlite::read_json(permit_pin$path, simplifyVector = FALSE)
  first <- as.numeric(args[5L])
  last <- as.numeric(args[6L])
  jgw11sh_permission(permit, z, actual, context, previous_pin, legacy, first, last)
  evidence <- lapply(permit$evidence, old$jgwl_record)
  jgw11sh_need(all(vapply(evidence, function(x) old$jgwl_inside(x$path, z$run), logical(1L))),
    "RunBoundIndependentShardEvidence")
  b <- old$jgwl_boot(z)$b
  previous <- readRDS(previous_pin$path)
  round_first <- permit$job$round_first
  round_last <- permit$job$round_last
  candidate$jgwsp_range_guard(previous, b, round_first, round_last)
  jgw11sh_need(identical(previous_pin$path,
    file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", round_first-1))),
    "SameCellCanonicalPreRoundPrefix")
  round_target <- file.path(b$store, sprintf("four_group_MC_TAIL_PREFIX_%05d.rds", round_last))
  jgw11sh_need(!file.exists(round_target) && !file.exists(paste0(round_target, ".partial")),
    "RoundPrefixAlreadyCompleteNoNewFit")
  output <- permit$output_leaf
  parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  jgw11sh_need(old$jgwl_inside(parent, file.path(z$run, "provenance")) &&
    identical(output, file.path(parent, basename(output))) && !basename(output) %in% c(".", "..") &&
    !file.exists(output) && dir.create(output), "ExclusiveNewShardReceiptLeaf")
  old$jgwl_json(list(status = "FIT_ONLY_SHARD_STARTED_NOT_NATIVE_EXIT", actual_executor = actual,
    legacy_scientific_entry = z$before$entry, scientific_identity_sha256 = b$identity_sha256,
    first = first, last = last, round_first = round_first, round_last = round_last,
    previous_prefix = previous_pin, executor_acceptance = permit_pin), file.path(output, "STARTED.json"))
  complete <- FALSE
  on.exit(if (!complete) old$jgwl_json(list(status = "STOP_RETAINED_NO_AUTOMATIC_REFIT",
    actual_executor = actual, actual_native_exit_observed_here = FALSE), file.path(output, "FAILED.json")), add = TRUE)
  fixed <- readRDS(b$context$fixed$path)
  result <- jgw11sh_fit(z, b, previous, round_first, round_last, first, last, actual, fixed,
    candidate$jgwsp_range_guard, function(path) file.exists(path) || dir.exists(path),
    b$e$jgwf_get_draw, b$e$jgwx_pin)
  old$jgwl_after(z)
  invisible(lapply(c(list(actual, z$before$entry, context, legacy, previous_pin, permit_pin), evidence), old$jgwl_record))
  old$jgwl_json(list(schema = "jgw03-eleven-worker-fit-only-shard-result-v1",
    status = "FIT_ONLY_SHARD_COMPLETE_NOT_PREFIX_OR_INDEPENDENT_ACCEPTANCE",
    accession = z$accession, space = "four_group_MC", first = first, last = last,
    round_first = round_first, round_last = round_last, actual_executor = actual,
    legacy_scientific_entry = z$before$entry, successor_context_source = context,
    scientific_identity_sha256 = result$scientific_identity_sha256,
    previous_prefix = previous_pin, allocation_records = result$allocation_records,
    executor_acceptance = permit_pin, evidence = evidence, legacy_acceptance = legacy,
    maximum_parallel_single_core_workers = 11L, internal_backend = "SerialParam", internal_workers = 1L,
    prefix_written = FALSE, original_or_452_CLI_success_asserted = FALSE,
    actual_native_exit_observed_here = FALSE, scientific_acceptance = FALSE,
    frozen_merge_helper_sha256_for_external_ordered_readback = jgw11sh_merge_sha,
    R = R.version.string, session_info = utils::capture.output(utils::sessionInfo())),
    file.path(output, "R_COMPLETE.json"))
  complete <- TRUE
  0L
}

if (sys.nframe() == 0L) {
  exit <- tryCatch(jgw11sh_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_FIT_ONLY_SHARD_STOP_RETAINED_NO_AUTOMATIC_REFIT")
    1L
  })
  quit(status = exit)
}
