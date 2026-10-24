# Independent transport-only successor. Declarations only; no fit or driver.
# Factory returns the original four-argument writer interface.
# Frozen native kernels, packer, statistics, tails and source manifest stay intact.

g128t_need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
g128t_sha <- function(path) digest::digest(file = path, algo = "sha256")
g128t_hash_ok <- function(value) is.character(value) && length(value) == 1L &&
  !is.na(value) && grepl("^[0-9a-f]{64}$", value)

g128t_check_execution <- function(pins) {
  g128t_need(is.list(pins) && all(c("transport_source", "entry_source",
    "original_source_manifest_sha256") %in% names(pins)) &&
    g128t_hash_ok(pins$original_source_manifest_sha256), "Explicit additional execution pins required")
  for (name in c("transport_source", "entry_source")) {
    record <- pins[[name]]
    g128t_need(is.list(record) && is.character(record$path) && length(record$path) == 1L &&
      file.exists(record$path) && g128t_hash_ok(record$sha256) &&
      identical(g128t_sha(record$path), record$sha256), "Additional transport/entry source pin mismatch")
  }
  invisible(TRUE)
}

g128t_location <- function(path, root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  g128t_need(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01") ||
    grepl("^R:/Shaoyao_storage/cache/GSE128623_[A-Za-z0-9_-]+$", root), "Use an explicit owned run root")
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  g128t_need((identical(parent, root) || startsWith(parent, paste0(root, "/"))) &&
    grepl("^[A-Za-z0-9_-]+[.]rds$", basename(path)), "Archive destination escapes its run")
  list(root = root, final = file.path(parent, basename(path)))
}

g128t_check_identity <- function(identity, pins) {
  required <- c("run_id", "cell_id", "reference_space", "allocation_id", "input_manifest_sha256",
    "source_manifest_sha256", "allocation_row_sha256", "packed_sha256")
  g128t_need(is.list(identity) && all(required %in% names(identity)), "Incomplete original archive identity")
  for (field in required[grepl("sha256$", required)]) g128t_need(g128t_hash_ok(identity[[field]]),
    "Invalid original archive identity hash")
  for (field in c("run_id", "cell_id", "reference_space")) g128t_need(is.character(identity[[field]]) &&
    length(identity[[field]]) == 1L && !is.na(identity[[field]]) && nzchar(identity[[field]]), "Invalid run/cell/space")
  g128t_need(length(identity$allocation_id) == 1L && is.finite(identity$allocation_id) &&
    identity$allocation_id >= 1L && identity$allocation_id == as.integer(identity$allocation_id) &&
    identical(identity$source_manifest_sha256, pins$original_source_manifest_sha256),
    "Allocation or unchanged original source-manifest coupling differs")
  invisible(TRUE)
}

# Compare every ordinary field, attribute, matrix and S4 slot recursively.
# Closures retain exact formals/body/attributes but lexical environment addresses
# and external pointer addresses are not a meaningful serialization equality.
# Their complete graph is still transported by R's single v3 serialization.
g128t_compare_content <- function(before, after) {
  counts <- new.env(parent = emptyenv())
  counts$ordinary_leaf_n <- 0L; counts$function_n <- 0L
  counts$environment_address_exempt_n <- 0L; counts$external_pointer_address_exempt_n <- 0L
  compare <- function(a, b, where) {
    g128t_need(identical(typeof(a), typeof(b)) && identical(class(a), class(b)),
      paste("Transport type/class changed:", where))
    if (is.environment(a)) {
      g128t_need(identical(environmentIsLocked(a), environmentIsLocked(b)), "Environment lock state changed")
      counts$environment_address_exempt_n <- counts$environment_address_exempt_n + 1L
      return(invisible(TRUE))
    }
    if (typeof(a) == "externalptr") {
      counts$external_pointer_address_exempt_n <- counts$external_pointer_address_exempt_n + 1L
      return(invisible(TRUE))
    }
    if (is.function(a)) {
      g128t_need(identical(formals(a), formals(b)) && identical(body(a), body(b)),
        paste("Closure code/formals changed:", where))
      counts$function_n <- counts$function_n + 1L
    } else if (isS4(a)) {
      slots <- methods::slotNames(a)
      g128t_need(identical(slots, methods::slotNames(b)) &&
        isTRUE(methods::validObject(a, test = TRUE)) && isTRUE(methods::validObject(b, test = TRUE)),
        paste("S4 slot schema/validity changed:", where))
      for (slot in slots) compare(methods::slot(a, slot), methods::slot(b, slot), paste0(where, "@", slot))
    } else if (is.call(a) || is.expression(a)) {
      # Formula/language attributes may include a lexical .Environment. Walk
      # the exact AST first, then apply the same attribute/environment rules.
      g128t_need(length(a) == length(b) && identical(names(a), names(b)), paste("Language shape changed:", where))
      for (index in seq_along(a)) compare(base::.subset2(a, index), base::.subset2(b, index), paste0(where, "/AST/", index))
    } else if (is.list(a) || is.pairlist(a)) {
      g128t_need(length(a) == length(b) && identical(names(a), names(b)), paste("List shape changed:", where))
      # Inspect stored members without S3 [[ dispatch (numeric_version's [[
      # deliberately preserves its wrapper). Exact class/attributes stay checked.
      for (index in seq_along(a)) compare(base::.subset2(a, index), base::.subset2(b, index), paste0(where, "[[", index, "]]"))
    } else {
      # Atomic identities include numeric precision, NA/NaN/Inf and exact signs.
      g128t_need(identical(a, b), paste("Ordinary field changed:", where))
      counts$ordinary_leaf_n <- counts$ordinary_leaf_n + 1L
    }
    aa <- attributes(a); bb <- attributes(b)
    if (isS4(a)) { aa <- aa[setdiff(names(aa), methods::slotNames(a))]; bb <- bb[setdiff(names(bb), methods::slotNames(b))] }
    if (inherits(a, "data.table")) {
      aa <- aa[setdiff(names(aa), ".internal.selfref")]; bb <- bb[setdiff(names(bb), ".internal.selfref")]
    }
    g128t_need(identical(names(aa), names(bb)), paste("Attribute schema changed:", where))
    for (name in names(aa)) compare(aa[[name]], bb[[name]], paste0(where, "/attr/", name))
    invisible(TRUE)
  }
  compare(before, after, "packed")
  as.list(counts)
}

g128t_native_check <- function(packed, bindings) {
  bindings$validate(packed)
  evidence <- packed$native_evidence
  g128t_need(is.list(packed$worker_diagnostics) && is.list(evidence), "Native/worker evidence missing")
  if (!identical(evidence$kind, "fixed_six_group_native_sentinel")) return(invisible(TRUE))
  bundle <- evidence$bundle
  g128t_need(all(vapply(bundle[c("gene_camera", "activity", "unsigned")], is.list, logical(1))),
    "Full native sentinel bundle missing")
  repacked <- bindings$pack(bundle, target_contrast = NULL)
  for (field in c("keys", "families", "statistic", "native_p", "native_diagnostics", "mroast_auxiliary", "rng_calls", "scope"))
    g128t_need(identical(repacked[[field]], packed[[field]]), paste("Native repack differs:", field))
  labels <- packed$worker_diagnostics$labels
  g128t_need(identical(bundle$integration_diagnostics$labels, labels) &&
    identical(bundle$integration_diagnostics$sample_ids, names(labels)) &&
    identical(colnames(bundle$gene_camera$native$dge$counts), names(labels)) &&
    identical(rownames(bundle$gene_camera$native$dge$counts), bundle$gene_camera$retained_gene_ids) &&
    identical(dimnames(bundle$gene_camera$native$voom$E), dimnames(bundle$gene_camera$native$dge$counts)),
    "Native sentinel sample/gene/label coupling differs")
  invisible(TRUE)
}

g128t_publish <- function(packed, path, run_root, identity, bindings, pins, recovery = NULL) {
  g128t_check_execution(pins); g128t_check_identity(identity, pins)
  location <- g128t_location(path, run_root)
  final <- location$final; partial <- paste0(final, ".transport.partial")
  old_partial <- paste0(final, ".partial")
  g128t_need(!file.exists(final) && !file.exists(partial), "Final/new transport partial exists; never overwrite")
  if (file.exists(old_partial)) g128t_need(is.list(recovery) &&
    identical(normalizePath(old_partial, winslash = "/"), recovery$path) &&
    identical(g128t_sha(old_partial), recovery$sha256), "Old partial requires separately proven zero-fit recovery")
  g128t_native_check(packed, bindings)
  # Exactly one explicit serialization/decode normalizes the fresh closure graph.
  # Do not claim that its bytes equal a new serialization of the original RAM.
  frozen_raw <- serialize(packed, NULL, version = 3L)
  frozen <- unserialize(frozen_raw)
  content_check <- g128t_compare_content(packed, frozen)
  g128t_native_check(frozen, bindings)
  old_fresh_hash <- identity$packed_sha256
  identity$packed_sha256 <- digest::digest(frozen, algo = "sha256")
  identity$transport_source_sha256 <- pins$transport_source$sha256
  identity$transport_entry_sha256 <- pins$entry_source$sha256
  identity$transport_original_source_manifest_sha256 <- pins$original_source_manifest_sha256
  identity$transport_old_fresh_packed_sha256 <- old_fresh_hash
  identity$transport_initial_serialization_sha256 <- digest::digest(frozen_raw, algo = "sha256", serialize = FALSE)
  identity$transport_digest_definition <- "digest(decoded packed after one v3 serialize/unserialize); not fresh-RAM raw-byte canonicality"
  identity$transport_content_check <- content_check
  identity$transport_recovery <- recovery
  payload <- list(schema_version = "gse128623-statistic-native-p-archive-v1", identity = identity,
    packed = frozen, real_execution_authority_granted = FALSE)
  saveRDS(payload, partial, compress = "gzip", version = 3L)
  check <- function(file) {
    restored <- readRDS(file)
    g128t_need(identical(restored$schema_version, payload$schema_version) &&
      identical(restored$identity, identity) && identical(restored$real_execution_authority_granted, FALSE) &&
      identical(digest::digest(restored$packed, algo = "sha256"), identity$packed_sha256),
      "Decoded archive schema/identity/stable-packed digest changed; preserve partial")
    g128t_compare_content(frozen, restored$packed)
    g128t_native_check(restored$packed, bindings)
    invisible(TRUE)
  }
  check(partial)
  g128t_check_execution(pins)
  g128t_need(file.link(partial, final), "Exclusive archive publication failed; preserve transport partial")
  check(final)
  if (!is.null(recovery)) g128t_need(identical(g128t_sha(recovery$path), recovery$sha256) &&
    unname(file.info(recovery$path)$size) == recovery$bytes, "Original partial changed during recovery")
  cleanup <- unlink(partial, recursive = FALSE)
  list(file = final, bytes = unname(file.info(final)$size), sha256 = g128t_sha(final),
    own_partial_removed = identical(cleanup, 0L))
}

g128t_make_writer <- function(original_bindings, execution_pins) {
  g128t_need(is.list(original_bindings) && is.function(original_bindings$pack) &&
    is.function(original_bindings$validate), "Original pinned pack/validate bindings required")
  g128t_check_execution(execution_pins)
  force(original_bindings); force(execution_pins)
  function(packed, path, run_root, identity) g128t_publish(packed, path, run_root, identity,
    original_bindings, execution_pins)
}

# Caller supplies exact independently reviewed partial pin and expected coupling.
# This function never runs a model or deletes/replaces the original .partial.
g128t_recover_partial <- function(partial_path, final_path, run_root, expected_identity,
                                  observed_keys, observed_families, original_bindings,
                                  execution_pins, expected_partial,
                                  expected_observed_evidence = NULL) {
  g128t_check_execution(execution_pins)
  location <- g128t_location(final_path, run_root)
  partial_path <- normalizePath(partial_path, winslash = "/", mustWork = TRUE)
  g128t_need(identical(partial_path, paste0(location$final, ".partial")) &&
    is.list(expected_partial) && g128t_hash_ok(expected_partial$sha256) &&
    identical(g128t_sha(partial_path), expected_partial$sha256) &&
    unname(file.info(partial_path)$size) == expected_partial$bytes,
    "Exact preserved original partial identity required")
  payload <- readRDS(partial_path)
  required_coupling <- c("run_id", "cell_id", "reference_space", "allocation_id",
    "input_manifest_sha256", "source_manifest_sha256", "allocation_row_sha256")
  g128t_need(is.list(expected_identity) && all(required_coupling %in% names(expected_identity)) &&
    !anyDuplicated(names(expected_identity)), "All seven caller-reviewed coupling fields required")
  g128t_need(identical(names(payload), c("schema_version", "identity", "packed", "real_execution_authority_granted")) &&
    identical(payload$schema_version, "gse128623-statistic-native-p-archive-v1") &&
    identical(payload$real_execution_authority_granted, FALSE) &&
    identical(payload$identity[names(expected_identity)], expected_identity), "Original partial coupling/schema differs")
  g128t_check_identity(payload$identity, execution_pins)
  packed <- payload$packed
  original_bindings$validate(packed)
  g128t_need(identical(packed$keys, observed_keys) && identical(packed$families, observed_families),
    "Original partial full feature/family domain differs")
  diagnostic <- packed$worker_diagnostics
  labels <- diagnostic$labels
  row_hash <- digest::digest(paste(c(paste(names(labels), collapse = "\t"),
    paste(unname(labels), collapse = "\t")), collapse = "\n"), algo = "sha256", serialize = FALSE)
  g128t_need(identical(diagnostic$allocation_id, expected_identity$allocation_id) &&
    identical(paste(diagnostic$processing, diagnostic$cohort, sep = "_"), expected_identity$cell_id) &&
    is.character(labels) && !anyNA(labels) && !is.null(names(labels)) && !anyDuplicated(names(labels)) &&
    identical(row_hash, expected_identity$allocation_row_sha256) &&
    identical(expected_identity$reference_space, "six_group") &&
    expected_identity$allocation_id %in% c(2L, 101L, 1001L, 5001L, 10001L) &&
    identical(packed$native_evidence$kind, "fixed_six_group_native_sentinel"),
    "Preserved sentinel allocation/label/native-evidence coupling differs")
  g128t_native_check(packed, original_bindings)
  recovery <- list(path = partial_path, bytes = unname(file.info(partial_path)$size),
    sha256 = expected_partial$sha256, old_fresh_packed_sha256 = payload$identity$packed_sha256,
    old_fresh_hash_not_reasserted_as_decoded_hash = TRUE, model_calls = 0L,
    original_partial_preserved = TRUE)
  g128t_publish(packed, location$final, location$root, payload$identity,
    original_bindings, execution_pins, recovery = recovery)
}
