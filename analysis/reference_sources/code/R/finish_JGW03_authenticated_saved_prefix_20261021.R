# Saved-only final arithmetic. Existing checkpoints/archives are never changed.
# No old tail-prefix rebuild, archive stream, native model or RNG is invoked.
# CLI: --execute CONFIG ACCESSION SPACE PREFIX BYTES SHA LEGACY_ACCEPT LEGACY_SHA
#      SOURCE452 SOURCE452_SHA PERMIT PERMIT_SHA
jgwsf_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwsf_root <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03"
jgwsf_old_sha <- "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d"
jgwsf_context_sha <- "452acfa10dbbaa29631ad41877feab120ba421dec6852bd42ab5180ef9b1b56c"
jgwsf_helper_sha <- "74bb8dc0d3cdb16bb9373fdf5a9e6153af7298d5c104276e58c3247fc41b972e"

jgwsf_same_pin <- function(a, b) {
  is.list(a) && is.list(b) && identical(a$path, b$path) &&
    is.numeric(a$bytes) && length(a$bytes) == 1L && is.finite(a$bytes) &&
    is.numeric(b$bytes) && length(b$bytes) == 1L && is.finite(b$bytes) &&
    a$bytes == b$bytes && identical(a$sha256, b$sha256)
}

jgwsf_permission <- function(p, z, actual, context, helper, previous, legacy, space) {
  jgwsf_need(is.list(p) && identical(p$schema, "jgw03-authenticated-saved-prefix-finish-acceptance-v1") &&
    identical(p$status, "ACCEPTED_FOR_ONE_SAVED_PREFIX_FINISH") && isTRUE(p$resources_admitted) &&
    isTRUE(p$independent_wiring_review_accepted) && isTRUE(p$previous_boundary_authenticated_by_caller) &&
    identical(p$maximum_parallel_single_core_workers, 3L) &&
    is.numeric(p$per_worker_memory_bytes) && length(p$per_worker_memory_bytes) == 1L &&
    is.finite(p$per_worker_memory_bytes) && p$per_worker_memory_bytes == 3*1024^3,
    "SeparateReviewedSavedFinishPermission")
  jgwsf_need(jgwsf_same_pin(p$actual_finish, actual) && jgwsf_same_pin(p$legacy_scientific_entry, z$before$entry) &&
    jgwsf_same_pin(p$successor_context_source, context) && identical(context$sha256, jgwsf_context_sha) &&
    jgwsf_same_pin(p$increment_helper, helper) && identical(helper$sha256, jgwsf_helper_sha) &&
    jgwsf_same_pin(p$config, z$before$config) && jgwsf_same_pin(p$legacy_acceptance, legacy) &&
    identical(p$job$accession, z$accession) && identical(p$job$space, space) &&
    jgwsf_same_pin(p$job$previous_prefix, previous), "ExactFinishSourcesCellSpaceAndPrefix")
  jgwsf_need(is.list(p$evidence) && length(p$evidence) == 3L && !anyDuplicated(names(p$evidence)) &&
    setequal(names(p$evidence), c("review", "resources", "previous_boundary")), "ThreeFinishEvidencePins")
  invisible(TRUE)
}

# Pure structural/domain check only. The caller authenticates the canonical
# prefix bytes/SHA and actual producer exit. No historical archive is re-read.
jgwsf_prefix <- function(prefix, b, space, template) {
  whole <- function(x, minimum = 0) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x >= minimum && x == floor(x)
  sha <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{64}$", x)
  total <- if (identical(space, "four_group_MC")) 10001L else 20L
  keys <- b$keys[[space]]
  fields <- c("schema", "space", "through_allocation_id", "scientific_identity_sha256",
    "full_registry", "observed", "completed_compact_pins", "state", "provisional_P_computed")
  jgwsf_need(space %in% c("four_group_MC", "high_companion_exact") &&
    is.list(prefix) && identical(names(prefix), fields) &&
    identical(prefix$schema, "jgw03-fixed-registry-tail-prefix-v1") && identical(prefix$space, space) &&
    whole(prefix$through_allocation_id) && prefix$through_allocation_id == total &&
    identical(prefix$scientific_identity_sha256, b$identity_sha256) &&
    identical(prefix$full_registry, b$context$label_registry) &&
    identical(prefix$observed, b$context$observed_files$observed) && identical(prefix$provisional_P_computed, FALSE) &&
    is.character(keys) && length(keys) == total && !anyNA(keys) && !anyDuplicated(keys), "AuthenticatedCompleteFinalPrefix")
  # Numeric/integer through is semantic equality only; never coerce the prefix.
  pins <- prefix$completed_compact_pins
  ids <- seq.int(2L, total)
  jgwsf_need(is.list(pins) && identical(names(pins), as.character(ids)), "EveryOrderedFinalPrefixPinID")
  for (i in ids) {
    one <- pins[[as.character(i)]]
    jgwsf_need(is.list(one) && identical(names(one), c("allocation_key", "archive")) &&
      identical(one$allocation_key, keys[i]) && is.list(one$archive) &&
      identical(names(one$archive), c("bytes", "sha256")) && whole(one$archive$bytes, 1) && sha(one$archive$sha256),
      "FinalPrefixKeyAndOriginalArchivePinSchema")
  }
  state <- prefix$state
  changing <- c("processed_nonidentity", "families", "execution_metadata_missing_allocations")
  fixed <- setdiff(names(template), changing)
  jgwsf_need(is.list(state) && identical(names(state), names(template)) &&
    identical(state[fixed], template[fixed]) && identical(state$processed_nonidentity, total-1L) &&
    is.integer(state$execution_metadata_missing_allocations) && whole(state$execution_metadata_missing_allocations) &&
    state$execution_metadata_missing_allocations <= total &&
    identical(names(state$families), names(template$families)), "FinalStateFixedDomainAndCompleteCount")
  for (id in names(template$families)) {
    x <- state$families[[id]]; expected <- template$families[[id]]
    hits <- x$tail_hits_nonidentity; missing <- x$null_unavailable
    jgwsf_need(is.list(x) && identical(names(x), names(expected)) &&
      identical(x$observed, expected$observed) && identical(x$observed_available, expected$observed_available) &&
      is.integer(hits) && is.integer(missing) && length(hits) == length(expected$observed) &&
      length(missing) == length(hits) && !anyNA(hits) && !anyNA(missing) &&
      all(hits >= 0L & missing >= 0L & hits+missing <= total-1L) &&
      all(hits[!expected$observed_available] == 0L), "FinalOriginalIntegerCountersAndObservedValues")
  }
  invisible(TRUE)
}

jgwsf_finish <- function(prefix, b, space, actual, tail_start, tail_finish, finite, finite_sha, writer) {
  template <- tail_start(b$observed, b$domain, space, b$keys[[space]])
  jgwsf_prefix(prefix, b, space, template)
  answer <- tail_finish(prefix$state, finite, finite_sha)
  product <- list(summary = answer, tail_checkpoint = prefix$state,
    scientific_identity_sha256 = b$identity_sha256, full_registry = b$context$registry,
    execution_entry = actual)
  target <- file.path(b$store, paste0(space, "_TAIL_SUMMARY.rds"))
  list(product = writer(product, target), target = target)
}

jgwsf_main <- function(args) {
  jgwsf_need(length(args) == 13L && identical(args[1L], "--execute") &&
    args[3L] %in% c("GSE289929", "GSE227494", "GSE299063") &&
    args[4L] %in% c("four_group_MC", "high_companion_exact") && grepl("^[0-9]+$", args[6L]) &&
    identical(args[11L], jgwsf_context_sha), "OneExplicitSavedFinishJob")
  old_path <- file.path(jgwsf_root, "project/code/R/run_JGW03_full_reference_linux_v3_20261017.R")
  jgwsf_need(identical(digest::digest(file = old_path, algo = "sha256"), jgwsf_old_sha), "FixedLegacyScientificSource")
  old <- new.env(parent = asNamespace("stats")); sys.source(old_path, old)
  actual <- old$jgwl_entry()
  jgwsf_need(old$jgwl_inside(actual$path, jgwsf_root) &&
    !old$jgwl_inside(actual$path, file.path(jgwsf_root, "project")), "HonestSeparateSavedFinishExecutor")
  context <- old$jgwl_pin(args[10L])
  jgwsf_need(identical(context$sha256, jgwsf_context_sha) && old$jgwl_inside(context$path, jgwsf_root), "Fixed452ContextLoader")
  candidate <- new.env(parent = asNamespace("stats")); sys.source(context$path, candidate)
  loader <- candidate$jgwsp_legacy_scientific_context; environment(loader) <- old
  z <- loader(args[2L], args[3L]); legacy <- old$jgwl_accept(z, args[8L], args[9L])
  prefix_pin <- old$jgwl_record(list(path = args[5L], bytes = as.numeric(args[6L]), sha256 = args[7L]))
  permit_pin <- old$jgwl_record(list(path = args[12L], bytes = unname(file.info(args[12L])$size), sha256 = args[13L]))
  permit <- jsonlite::read_json(permit_pin$path, simplifyVector = FALSE)
  helper <- old$jgwl_record(permit$increment_helper)
  jgwsf_permission(permit, z, actual, context, helper, prefix_pin, legacy, args[4L])
  evidence <- lapply(permit$evidence, old$jgwl_record)
  jgwsf_need(old$jgwl_inside(helper$path, z$run) &&
    all(vapply(evidence, function(x) old$jgwl_inside(x$path, z$run), logical(1L))), "RunBoundFinishEvidence")
  b <- old$jgwl_boot(z)$b
  space <- args[4L]; total <- if (space == "four_group_MC") 10001L else 20L
  jgwsf_need(identical(prefix_pin$path, file.path(b$store, sprintf("%s_TAIL_PREFIX_%05d.rds", space, total))),
    "SameCellCanonicalFinalPrefix")
  target <- file.path(b$store, paste0(space, "_TAIL_SUMMARY.rds"))
  jgwsf_need(!file.exists(target) && !file.exists(paste0(target, ".partial")), "FreshSummaryNoAutomaticRetry")
  output <- permit$output_leaf
  parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  jgwsf_need(old$jgwl_inside(parent, file.path(z$run, "provenance")) &&
    identical(output, file.path(parent, basename(output))) && !basename(output) %in% c(".", "..") &&
    !file.exists(output) && dir.create(output), "ExclusiveSavedFinishReceiptLeaf")
  old$jgwl_json(list(status = "SAVED_FINISH_STARTED_NOT_NATIVE_EXIT", actual_finish = actual,
    legacy_scientific_entry = z$before$entry, final_prefix = prefix_pin), file.path(output, "STARTED.json"))
  completed <- FALSE
  on.exit(if (!completed) old$jgwl_json(list(status = "STOP_RETAINED_NO_AUTOMATIC_REFINISH",
    actual_native_exit_observed_here = FALSE), file.path(output, "FAILED.json")), add = TRUE)
  prefix <- readRDS(prefix_pin$path)
  result <- jgwsf_finish(prefix, b, space, actual, b$h$jgw03c_tail_start, b$h$jgw03c_tail_finish,
    list(MC = b$h$jgw03_finite_mc_interval, companion = b$h$jgw03_exact20_companion_interval),
    b$modules$pins$finite$sha256, b$e$jgwx_atomic_plain)
  old$jgwl_after(z)
  invisible(lapply(c(list(actual, z$before$entry, context, helper, legacy, prefix_pin, permit_pin), evidence), old$jgwl_record))
  old$jgwl_json(list(schema = "jgw03-authenticated-saved-prefix-finish-result-v1",
    status = "SAVED_PREFIX_ARITHMETIC_COMPLETE_NOT_FULL_ARCHIVE_ACCEPTANCE", accession = z$accession,
    space = space, actual_finish = actual, legacy_scientific_entry = z$before$entry,
    successor_context_source = context, scientific_identity_sha256 = b$identity_sha256,
    final_prefix = prefix_pin, summary = c(list(path = result$target), result$product),
    executor_acceptance = permit_pin, legacy_acceptance = legacy, evidence = evidence,
    original_CLI_success_asserted = FALSE, actual_native_exit_observed_here = FALSE,
    new_fits = 0L, allocation_RNG_calls = 0L, archives_reread = 0L, prefix_modified = FALSE,
    scientific_acceptance = FALSE, full_archive_independent_checker_still_required = TRUE,
    R = R.version.string, session_info = utils::capture.output(utils::sessionInfo())), file.path(output, "R_COMPLETE.json"))
  completed <- TRUE
  0L
}

if (sys.nframe() == 0L) {
  status <- tryCatch(jgwsf_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_SAVED_FINISH_STOP_RETAINED_NO_AUTOMATIC_RETRY"); 1L
  })
  quit(status = status)
}
