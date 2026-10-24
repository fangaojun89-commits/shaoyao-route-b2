# Definitions only. No files, packages, native fits, RNG, finish-P or BH calls.
# Caller MUST authenticate the previous checkpoint and its upstream scientific
# identity/pins. This helper checks structure/continuity, not historical truth,
# and does not re-read/re-pin old archives. Only new IDs reach the callbacks.
# Minimal b: identity_sha256, keys, context$label_registry,
# context$observed_files$observed, store. No CLI/bootstrap or scientific RDS.
# Callbacks: read_archive(b, key), tail_update(state, compact), pin(path).
# pin() must return the ORIGINAL jgwx_pin shape: list(bytes, sha256), no path.

jgw03_advance_saved_tail_prefix <- function(previous, b, last,
                                           read_archive, tail_update, pin) {
  need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
  scalar_count <- function(x, minimum = 0) is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x >= minimum && x == floor(x)
  sha <- function(x) is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^[a-f0-9]{64}$", x)
  valid_pin <- function(x) is.list(x) && identical(names(x), c("bytes", "sha256")) &&
    scalar_count(x$bytes, 1) && sha(x$sha256)
  need(is.function(read_archive) && is.function(tail_update) && is.function(pin),
    "ExplicitSavedOnlyCallbacksRequired")
  fields <- c("schema", "space", "through_allocation_id", "scientific_identity_sha256",
    "full_registry", "observed", "completed_compact_pins", "state", "provisional_P_computed")
  need(is.list(previous) && identical(names(previous), fields) &&
    identical(previous$schema, "jgw03-fixed-registry-tail-prefix-v1") &&
    previous$space %in% c("four_group_MC", "high_companion_exact") &&
    identical(previous$provisional_P_computed, FALSE), "OriginalPrefixSchemaAndScope")
  need(is.list(b) && sha(b$identity_sha256) &&
    identical(previous$scientific_identity_sha256, b$identity_sha256) &&
    identical(previous$full_registry, b$context$label_registry) &&
    identical(previous$observed, b$context$observed_files$observed) &&
    is.character(b$store) && length(b$store) == 1L && !is.na(b$store) && nzchar(b$store),
    "AuthenticatedParentScientificIdentityRegistryObservedAndStore")
  space <- previous$space
  keys <- b$keys[[space]]
  total <- if (space == "four_group_MC") 10001L else 20L
  need(is.character(keys) && length(keys) == total && !anyNA(keys) &&
    all(nzchar(keys)) && !anyDuplicated(keys), "CompleteDistinctFrozenRegistryKeys")
  old_last <- previous$through_allocation_id
  need(scalar_count(old_last, 1) && old_last <= total && scalar_count(last, 1) &&
    last > old_last && last <= total, "StrictNewBoundedPrefixNoDuplicateOrBackwardAdvance")
  state_fields <- c("schema", "domain_sha256", "space", "total_N", "allocation_keys",
    "processed_nonidentity", "families", "execution_metadata_missing_allocations", "definition")
  validate_state <- function(state, through) {
    need(is.list(state) && identical(names(state), state_fields) &&
      identical(state$schema, "jgw03-tail-accumulator-v1") && sha(state$domain_sha256) &&
      identical(state$space, space) && identical(state$total_N, total) &&
      identical(state$allocation_keys, keys) && is.integer(state$processed_nonidentity) &&
      identical(state$processed_nonidentity, as.integer(through - 1)) &&
      is.integer(state$execution_metadata_missing_allocations) &&
      scalar_count(state$execution_metadata_missing_allocations) &&
      state$execution_metadata_missing_allocations <= through &&
      is.list(state$families) && length(state$families) > 0L &&
      !is.null(names(state$families)) && !anyDuplicated(names(state$families)) &&
      all(nzchar(names(state$families))) && is.list(state$definition), "OriginalTailStateDomainAndCount")
    for (family in state$families) {
      need(is.list(family) && identical(names(family), c("observed", "observed_available",
        "tail_hits_nonidentity", "null_unavailable")) && is.numeric(family$observed) &&
        length(family$observed) > 0L && is.logical(family$observed_available) &&
        !anyNA(family$observed_available) && length(family$observed_available) == length(family$observed),
        "OriginalTailFamilyShape")
      hits <- family$tail_hits_nonidentity
      unavailable <- family$null_unavailable
      need(is.integer(hits) && is.integer(unavailable) &&
        length(hits) == length(family$observed) && length(unavailable) == length(hits) &&
        !anyNA(hits) && !anyNA(unavailable) && all(hits >= 0L & unavailable >= 0L) &&
        all(hits + unavailable <= through - 1) && all(hits[!family$observed_available] == 0L),
        "OriginalIntegerTailCounters")
    }
    invisible(TRUE)
  }
  state <- previous$state
  validate_state(state, old_last)
  pins <- previous$completed_compact_pins
  old_ids <- if (old_last >= 2) seq.int(2L, old_last) else integer()
  expected_names <- if (length(old_ids)) as.character(old_ids) else NULL
  need(is.list(pins) && length(pins) == length(old_ids) &&
    identical(names(pins), expected_names), "CompleteOrderedNonduplicateOldPrefixPinIDs")
  for (i in old_ids) {
    record <- pins[[as.character(i)]]
    need(is.list(record) && identical(names(record), c("allocation_key", "archive")) &&
      identical(record$allocation_key, keys[i]) && valid_pin(record$archive), "OldPrefixPinKeyAndSchema")
  }
  for (i in seq.int(old_last + 1, last)) {
    archive <- read_archive(b, keys[i])
    need(is.list(archive) && identical(archive$schema, "jgw03-full-reference-allocation-v1") &&
      identical(archive$scientific_identity_sha256, b$identity_sha256) &&
      identical(archive$allocation_key, keys[i]) && is.list(archive$compact) &&
      identical(archive$compact$allocation_key, keys[i]) &&
      identical(archive$compact$domain_sha256, state$domain_sha256) &&
      identical(archive$compact$definition, state$definition), "NextAuthenticatedSavedArchiveRequired")
    updated <- tail_update(state, archive$compact)
    need(is.list(updated) && identical(updated$processed_nonidentity, as.integer(i - 1)) &&
      identical(updated[setdiff(state_fields, c("families", "processed_nonidentity",
        "execution_metadata_missing_allocations"))], state[setdiff(state_fields,
        c("families", "processed_nonidentity", "execution_metadata_missing_allocations"))]) &&
      identical(names(updated$families), names(state$families)) &&
      all(vapply(names(state$families), function(id)
        identical(updated$families[[id]]$observed, state$families[[id]]$observed) &&
        identical(updated$families[[id]]$observed_available, state$families[[id]]$observed_available), logical(1L))),
      "TailCallbackPreservesFrozenDomainAndAdvancesExactlyOnce")
    state <- updated
    archive_pin <- pin(file.path(b$store, keys[i], "COMPACT.rds"))
    need(valid_pin(archive_pin), "NewSavedArchiveOriginalPinSchema")
    pins[[as.character(i)]] <- list(allocation_key = keys[i], archive = archive_pin)
  }
  validate_state(state, last)
  list(schema = "jgw03-fixed-registry-tail-prefix-v1", space = space,
    through_allocation_id = last, scientific_identity_sha256 = b$identity_sha256,
    full_registry = b$context$label_registry, observed = b$context$observed_files$observed,
    completed_compact_pins = pins, state = state, provisional_P_computed = FALSE)
}
