# Functions only. Supplement to the frozen numerical archive reader, not a fit,
# data repair, archive authorizer or replacement numerical/biological acceptance.
# Caller authenticates the archive bytes and derives expected mode from the
# historical creation manifests/recovery receipt, NEVER from this payload alone.
g128l_check <- function(payload, expected_identity, expected_mode, pins,
                        expected_recovery = NULL) {
  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  hash <- function(x) is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^[0-9a-f]{64}$", x)
  scalar_integer <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x >= 0 && x == floor(x)
  named_list <- function(x, fields) is.list(x) && !anyDuplicated(names(x)) &&
    setequal(names(x), fields)
  base_fields <- c("run_id", "cell_id", "reference_space", "allocation_id",
    "input_manifest_sha256", "source_manifest_sha256", "allocation_row_sha256")
  pin_fields <- c("transport_source_sha256", "fresh_entry_sha256",
    "recovery_entry_sha256", "original_source_manifest_sha256")
  transport_fields <- c("transport_source_sha256", "transport_entry_sha256",
    "transport_original_source_manifest_sha256", "transport_old_fresh_packed_sha256",
    "transport_initial_serialization_sha256", "transport_digest_definition",
    "transport_content_check")
  need(is.character(expected_mode) && length(expected_mode) == 1L && !is.na(expected_mode) &&
    expected_mode %in% c("legacy", "transport_fresh", "transport_recovered"),
    "Explicit externally authenticated creation mode required")
  need(named_list(expected_identity, base_fields) &&
    all(vapply(expected_identity[c("run_id", "cell_id", "reference_space")],
      function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x), logical(1))) &&
    scalar_integer(expected_identity$allocation_id) && expected_identity$allocation_id >= 1 &&
    all(vapply(expected_identity[grepl("sha256$", names(expected_identity))], hash, logical(1))),
    "Explicit complete original identity required")
  need(named_list(pins, pin_fields) && all(vapply(pins, hash, logical(1))) &&
    identical(expected_identity$source_manifest_sha256, pins$original_source_manifest_sha256),
    "Expected original/effective source pins required")
  need(named_list(payload, c("schema_version", "identity", "packed", "real_execution_authority_granted")) &&
    identical(payload$schema_version, "gse128623-statistic-native-p-archive-v1") &&
    identical(payload$real_execution_authority_granted, FALSE) && is.list(payload$packed),
    "Archive envelope/schema differs")
  identity <- payload$identity
  need(is.list(identity) && !anyDuplicated(names(identity)) &&
    identical(identity[names(expected_identity)], expected_identity) && hash(identity$packed_sha256) &&
    identical(identity$packed_sha256, digest::digest(payload$packed, algo = "sha256")),
    "Original identity or stable decoded packed digest differs")
  recovered <- identical(expected_mode, "transport_recovered")
  expected_fields <- c(base_fields, "packed_sha256",
    if (!identical(expected_mode, "legacy")) transport_fields,
    if (recovered) "transport_recovery")
  need(setequal(names(identity), expected_fields),
    "Missing, partial, unknown or incorrectly legacy-classified identity fields")
  if (identical(expected_mode, "legacy")) {
    need(is.null(expected_recovery), "Legacy archive cannot carry a recovery expectation")
    return(list(status = "PASS_LEGACY_IDENTITY_ONLY", mode = expected_mode,
      decoded_packed_digest_checked = TRUE, supplementary_fields_applicable = FALSE,
      numerical_or_native_acceptance_granted = FALSE))
  }
  for (field in transport_fields[grepl("sha256$", transport_fields)])
    need(hash(identity[[field]]), "Malformed supplementary hash")
  need(identical(identity$transport_source_sha256, pins$transport_source_sha256) &&
    identical(identity$transport_entry_sha256,
      pins[[if (recovered) "recovery_entry_sha256" else "fresh_entry_sha256"]]) &&
    identical(identity$transport_original_source_manifest_sha256, pins$original_source_manifest_sha256) &&
    identical(identity$transport_digest_definition,
      "digest(decoded packed after one v3 serialize/unserialize); not fresh-RAM raw-byte canonicality"),
    "Supplementary source/entry/digest-definition coupling differs")
  count_fields <- c("ordinary_leaf_n", "function_n", "environment_address_exempt_n",
                     "external_pointer_address_exempt_n")
  need(named_list(identity$transport_content_check, count_fields) &&
    all(vapply(identity$transport_content_check, scalar_integer, logical(1))) &&
    identity$transport_content_check$ordinary_leaf_n > 0,
    "Malformed transport comparison counters")
  if (recovered) {
    fields <- c("path", "bytes", "sha256", "old_fresh_packed_sha256",
      "old_fresh_hash_not_reasserted_as_decoded_hash", "model_calls", "original_partial_preserved")
    record <- identity$transport_recovery
    need(named_list(expected_recovery, c(fields, "decoded_packed_sha256")) &&
      named_list(record, fields) &&
      identical(expected_identity$reference_space, "six_group") &&
      expected_identity$allocation_id == 2L &&
      is.character(record$path) && length(record$path) == 1L && !is.na(record$path) &&
      identical(record$path, paste0("/root/autodl-tmp/", expected_identity$run_id,
        "/allocations/", expected_identity$cell_id, "/six_group/allocation_00002.rds.partial")) &&
      identical(record$path, expected_recovery$path) &&
      scalar_integer(record$bytes) && record$bytes > 0 &&
      scalar_integer(expected_recovery$bytes) && as.numeric(record$bytes) == as.numeric(expected_recovery$bytes),
      "Recovered archive must match an externally verified original partial")
    for (field in c("sha256", "old_fresh_packed_sha256")) need(hash(record[[field]]) &&
      identical(record[[field]], expected_recovery[[field]]), "Recovery hash differs")
    need(hash(expected_recovery$decoded_packed_sha256) &&
      identical(identity$packed_sha256, expected_recovery$decoded_packed_sha256) &&
      identical(identity$transport_old_fresh_packed_sha256, record$old_fresh_packed_sha256) &&
      identical(record$old_fresh_hash_not_reasserted_as_decoded_hash, TRUE) &&
      identical(record$original_partial_preserved, TRUE) &&
      scalar_integer(record$model_calls) && record$model_calls == 0 &&
      identical(expected_recovery$old_fresh_hash_not_reasserted_as_decoded_hash, TRUE) &&
      identical(expected_recovery$original_partial_preserved, TRUE) &&
      scalar_integer(expected_recovery$model_calls) && expected_recovery$model_calls == 0,
      "Recovered stable/old digest or zero-fit history differs")
  } else need(is.null(expected_recovery), "Fresh archive cannot carry a recovery expectation")
  list(status = "PASS_SUPPLEMENTARY_TRANSPORT_LINEAGE_ONLY", mode = expected_mode,
    decoded_packed_digest_checked = TRUE, expected_effective_source_pins_checked = TRUE,
    original_partial_binding_checked_against_caller = recovered,
    initial_serialization_hash = "FORMAT_AND_DECLARED_DEFINITION_ONLY_RAW_BUFFER_NOT_RETAINED",
    fresh_RAM_hash = if (recovered) "BOUND_TO_CALLER_AUTHENTICATED_OLD_PARTIAL" else
      "FORMAT_AND_PRODUCER_PROCEDURE_ONLY_ORIGINAL_FRESH_RAM_UNAVAILABLE",
    lexical_environment_binding_values_independently_compared = FALSE,
    numerical_or_native_acceptance_granted = FALSE)
}
