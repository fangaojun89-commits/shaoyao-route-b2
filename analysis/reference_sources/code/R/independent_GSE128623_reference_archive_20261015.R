# Independent schema-v1 archive reader: definitions only, no producer imports.
# Contract source: GSE128623_FULL_LABEL_REFERENCE_EXECUTION_PROPOSAL_20261015.md.
# Wire producer reviewed at SHA256 920928115af859f608b7d03349e1e7979cfbca353d87c3ba011cd0b0b3cefca5.
# Required allocation IDs and identity hashes must come from the frozen caller
# registry. This helper neither authorizes a run nor certifies 268 real families.

g128i_require <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

g128i_hash_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}

g128i_family_template <- function(target_contrast = NULL) {
  tests <- c("chronic_combined", "D28_vs_control", "D42_vs_control", "acute_combined",
             "D4_vs_control", "D7_vs_control", "D14_vs_control")
  collections <- c("DoRothEA_AB", "PROGENy", "Reactome", "SynGO_assignedBy_proxy_sensitivity")
  specification <- list(edgeR = "GENE", voom = "GENE", DESeq2 = "GENE",
    activity = collections, CAMERA = collections[3:4],
    mroast = collections[3:4], fgsea = collections[3:4])
  blocks <- lapply(names(specification), function(method) {
    selected_tests <- if (method %in% c("mroast", "fgsea")) tests[1L] else tests
    do.call(rbind, lapply(specification[[method]], function(collection) {
      data.frame(method = method, collection = collection, contrast = selected_tests,
        stringsAsFactors = FALSE)
    }))
  })
  families <- do.call(rbind, blocks)
  families$family_id <- with(families, paste(method, collection, contrast, sep = "::"))
  if (!is.null(target_contrast)) {
    g128i_require(is.character(target_contrast) && length(target_contrast) == 1L &&
      !is.na(target_contrast) && target_contrast %in% tests[c(2L, 3L, 5L, 6L, 7L)],
      "Only an explicitly named fixed day/control companion is supported.")
    families <- families[families$contrast == target_contrast, , drop = FALSE]
  }
  rownames(families) <- NULL
  families
}

# Independent BH algebra. NA/NaN positions are retained; the adjustment count
# is the actual nonmissing vector length, matching native p.adjust defaults.
g128i_bh <- function(probability) {
  result <- probability
  present <- which(!is.na(probability))
  n <- length(present)
  if (n == 0L) return(result)
  ordering <- order(probability[present], decreasing = TRUE)
  result[present[ordering]] <- pmin(1,
    cummin((n / seq.int(n, 1L)) * probability[present[ordering]]))
  result
}

g128i_source_q <- function(probability, method) {
  if (identical(method, "mroast")) {
    corrected <- probability - 0.5 / 9999
    corrected <- pmax(corrected, 0)
    return(pmax(probability, g128i_bh(corrected)))
  }
  g128i_bh(probability)
}

g128i_evidence_available <- function(statistic, probability, method) {
  available <- is.finite(statistic)
  zero_p_transform <- method %in% c("CAMERA", "mroast") &
    !is.na(probability) & probability == 0 & is.infinite(statistic) & statistic > 0
  available | zero_p_transform
}

g128i_native_initial_state <- function() {
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) assign(".Random.seed", previous_seed, .GlobalEnv) else {
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20260819L)
  get(".Random.seed", .GlobalEnv, inherits = FALSE)
}

g128i_check_rng <- function(records, keys, initial_state, observed_records = NULL,
                            retained_chronic = FALSE) {
  unsigned_collections <- c("Reactome", "SynGO_assignedBy_proxy_sensitivity")
  collection_order <- unique(keys$collection[keys$method == "mroast"])
  fgsea_collections <- unique(keys$collection[keys$method == "fgsea"])
  g128i_require(all(collection_order %in% unsigned_collections) &&
    all(fgsea_collections %in% collection_order), "Unsigned collection/input domains differ.")
  expected <- list()
  for (collection in collection_order) for (method in c("mroast", "fgsea")) {
    selected <- keys$method == method & keys$collection == collection
    if (any(selected)) expected[[length(expected) + 1L]] <- list(
      collection = collection, method = if (method == "fgsea") "fgseaMultilevel" else method,
      feature_ids = keys$feature_id[selected])
  }
  if (retained_chronic) {
    g128i_require(is.list(records), "Invalid retained chronic RNG trace.")
    roster <- vapply(records, function(call) {
      g128i_require(is.list(call) && is.character(call$collection) && length(call$collection) == 1L &&
        is.character(call$method) && length(call$method) == 1L, "Invalid retained RNG roster.")
      paste(call$collection, call$method, sep = "::")
    }, character(1L))
    retained_collections <- unique(vapply(records, `[[`, character(1L), "collection"))
    g128i_require(all(retained_collections %in% unsigned_collections),
      "Unknown collection in retained chronic trace.")
    allowed <- as.vector(t(outer(retained_collections, c("mroast", "fgseaMultilevel"),
      paste, sep = "::")))
    positions <- match(roster, allowed)
    g128i_require(!anyNA(positions) && !anyDuplicated(positions) &&
      all(diff(positions) > 0), "Retained chronic RNG roster/order changed.")
    expected <- lapply(records, function(call) list(collection = call$collection,
      method = call$method, feature_ids = call$program_ids))
  }
  g128i_require(is.list(records) && length(records) == length(expected),
    "Unsigned RNG call coverage differs from the archived method roster.")
  if (length(expected) == 0L) return(invisible(TRUE))
  rng_kind <- c("Mersenne-Twister", "Inversion", "Rejection")
  chronic <- c(control = -1, D4 = 0, D7 = 0, D14 = 0, D28 = .5, D42 = .5)
  ranking_anchor <- NULL
  for (index in seq_along(expected)) {
    call <- records[[index]]
    required <- c("collection", "method", "seed", "initial_state", "final_state",
      "rng_kind", "final_rng_kind", "contrast", "program_ids", "ranking", "backend")
    g128i_require(is.list(call) && !anyDuplicated(names(call)) &&
      all(required %in% names(call)), "Incomplete RNG call record.")
    g128i_require(identical(call$collection, expected[[index]]$collection) &&
      identical(call$method, expected[[index]]$method) &&
      is.numeric(call$seed) && length(call$seed) == 1L && !is.na(call$seed) &&
      call$seed == 20260819 && identical(call$rng_kind, rng_kind) &&
      identical(call$final_rng_kind, rng_kind), "Unsigned RNG method/seed/kind changed.")
    g128i_require(is.integer(call$initial_state) && identical(call$initial_state, initial_state) &&
      is.integer(call$final_state) && length(call$final_state) == length(initial_state) &&
      !anyNA(call$final_state) && call$final_state[1L] == initial_state[1L] &&
      call$final_state[2L] >= 0L && call$final_state[2L] <= 624L,
      "Unsigned RNG initial state or final-state structure differs.")
    g128i_require(is.numeric(call$contrast) &&
      identical(names(call$contrast), names(chronic)) &&
      identical(as.numeric(call$contrast), as.numeric(chronic)), "Unsigned contrast changed.")
    programs <- call$program_ids
    g128i_require(is.character(programs) && !anyNA(programs) &&
      !anyDuplicated(programs) && all(nzchar(trimws(programs))) &&
      length(programs) > 0L && length(programs) == length(expected[[index]]$feature_ids) &&
      setequal(programs, expected[[index]]$feature_ids), "RNG input program roster differs.")
    if (identical(call$method, "mroast")) g128i_require(
      identical(programs, expected[[index]]$feature_ids), "mroast input/program row order changed.")
    if (!is.null(observed_records)) g128i_require(length(records) == length(observed_records) &&
      identical(call$collection, observed_records[[index]]$collection) &&
      identical(call$method, observed_records[[index]]$method) &&
      identical(programs, observed_records[[index]]$program_ids),
      "Unsigned call roster/input program order changed across allocations.")
    g128i_require(is.list(call$ranking) &&
      identical(call$ranking$rule,
        "max_abs_signed_sqrt_QLF_F_per_symbol_then_score_desc_symbol_asc") &&
      is.numeric(call$ranking$count) && length(call$ranking$count) == 1L &&
      is.finite(call$ranking$count) && call$ranking$count >= 0 &&
      call$ranking$count == floor(call$ranking$count) && g128i_hash_string(call$ranking$sha256),
      "Invalid unsigned ranking metadata.")
    current_ranking <- list(count = as.numeric(call$ranking$count),
      sha256 = as.character(call$ranking$sha256))
    if (is.null(ranking_anchor)) ranking_anchor <- current_ranking else g128i_require(
      identical(current_ranking, ranking_anchor), "Within-allocation ranking metadata disagree.")
    backend <- call$backend
    g128i_require(is.list(backend) && is.character(backend$class) &&
      length(backend$class) == 1L && backend$class == "SerialParam" &&
      is.numeric(backend$workers) && length(backend$workers) == 1L &&
      !is.na(backend$workers) && backend$workers == 1 &&
      "RNGseed" %in% names(backend) && is.null(backend$RNGseed),
      "Unsigned backend differs from SerialParam/RNGseed NULL.")
  }
  for (collection in unsigned_collections) {
    roast <- which(vapply(records, function(call) identical(call$collection, collection) &&
      identical(call$method, "mroast"), logical(1L)))
    fg <- which(vapply(records, function(call) identical(call$collection, collection) &&
      identical(call$method, "fgseaMultilevel"), logical(1L)))
    if (length(fg)) g128i_require(length(roast) == 1L &&
      all(records[[fg]]$program_ids %in% records[[roast]]$program_ids),
      "fgsea RNG programs are not a subset of the original mroast input domain.")
  }
  invisible(TRUE)
}

g128i_check_packed <- function(packed, template, initial_state, observed_pack = NULL,
                               allow_retained_chronic_trace = FALSE) {
  key_names <- c("method", "collection", "contrast", "feature_id", "family_id")
  required <- c("keys", "statistic", "native_p", "families", "native_diagnostics",
    "mroast_auxiliary", "rng_calls")
  g128i_require(is.list(packed) && !anyDuplicated(names(packed)) &&
    all(required %in% names(packed)), "Required packed fields are missing.")
  keys <- packed$keys
  g128i_require(is.data.frame(keys) && identical(names(keys), key_names) &&
    all(vapply(keys, is.character, logical(1L))) && !anyNA(keys) && !anyDuplicated(keys) &&
    all(nzchar(trimws(unlist(keys, use.names = FALSE)))), "Invalid ordered feature keys.")
  g128i_require(is.data.frame(packed$families) && identical(packed$families, template) &&
    all(keys$family_id %in% template$family_id) &&
    identical(with(keys, paste(method, collection, contrast, sep = "::")), keys$family_id),
    "Family template/feature association differs from the declared scope.")
  n <- nrow(keys)
  for (field in c("statistic", "native_p")) g128i_require(
    is.double(packed[[field]]) && is.null(attributes(packed[[field]])) &&
      length(packed[[field]]) == n, "Statistic/P must be full plain double vectors.")
  p <- packed$native_p
  statistic <- packed$statistic
  g128i_require(all(is.na(p) | (is.finite(p) & p >= 0 & p <= 1)), "Native P out of range.")
  transformed <- which(keys$method %in% c("CAMERA", "mroast"))
  reconstructed <- -log10(p[transformed])
  actual <- statistic[transformed]
  g128i_require(identical(is.na(actual), is.na(reconstructed)) &&
    identical(is.nan(actual), is.nan(reconstructed)) &&
    identical(is.infinite(actual), is.infinite(reconstructed)) &&
    all(actual[is.infinite(actual)] == reconstructed[is.infinite(reconstructed)]),
    "P-transform statistic missing/NaN/infinite positions differ from native P.")
  finite <- is.finite(reconstructed)
  g128i_require(all(abs(actual[finite] - reconstructed[finite]) <=
    1e-12 * (1 + abs(reconstructed[finite]))), "Finite P-transform statistic differs from native P.")
  absolute <- which(keys$method %in% c("voom", "DESeq2", "activity", "fgsea") & is.finite(statistic))
  g128i_require(all(statistic[absolute] >= 0), "A declared absolute outer statistic is negative.")
  if (!is.null(observed_pack)) g128i_require(identical(keys, observed_pack$keys) &&
    identical(packed$families, observed_pack$families), "Draw key/order/family domain changed.")
  available <- g128i_evidence_available(statistic, p, keys$method)
  diagnostics <- lapply(seq_len(nrow(template)), function(index) {
    selected <- which(keys$family_id == template$family_id[index])
    probability <- p[selected]
    q <- g128i_source_q(probability, template$method[index])
    t <- statistic[selected]
    data.frame(template[index, , drop = FALSE], n = length(selected),
      available_statistic_n = sum(available[selected]), native_p_available_n = sum(!is.na(probability)),
      native_q_available_n = sum(!is.na(q)), native_p_lt_05_n = sum(probability < .05, na.rm = TRUE),
      native_q_le_05_n = sum(q <= .05, na.rm = TRUE), native_q_any_le_05 = any(q <= .05, na.rm = TRUE),
      native_p_zero_n = sum(probability == 0, na.rm = TRUE),
      positive_infinite_statistic_n = sum(is.infinite(t) & t > 0),
      negative_edge_F_n = if (template$method[index] == "edgeR") sum(t < 0, na.rm = TRUE) else 0L,
      stringsAsFactors = FALSE)
  })
  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  g128i_require(identical(packed$native_diagnostics, diagnostics),
    "Stored native diagnostics differ from independent source-BH reconstruction.")
  mroast_keys <- keys[keys$method == "mroast", , drop = FALSE]
  rownames(mroast_keys) <- NULL
  if (nrow(mroast_keys) == 0L) {
    g128i_require(is.null(packed$mroast_auxiliary), "Unexpected mroast auxiliary data.")
  } else {
    auxiliary <- packed$mroast_auxiliary
    extra <- c("native_fdr", "mixed_p", "mixed_fdr")
    g128i_require(is.data.frame(auxiliary) &&
      identical(names(auxiliary), c(key_names, extra)) &&
      identical(auxiliary[key_names], mroast_keys), "mroast auxiliary key/order differs.")
    for (field in extra) g128i_require(is.double(auxiliary[[field]]) &&
      is.null(attributes(auxiliary[[field]])) && all(is.na(auxiliary[[field]]) |
        (is.finite(auxiliary[[field]]) & auxiliary[[field]] >= 0 & auxiliary[[field]] <= 1)),
      "Invalid native mroast auxiliary probabilities.")
  }
  companion <- nrow(template) == 9L
  retained <- companion && length(packed$rng_calls) > 0L
  if (retained) g128i_require(allow_retained_chronic_trace &&
    is.list(packed$worker_diagnostics) &&
    identical(packed$worker_diagnostics$unsigned_trace_scope,
      "RETAINED_OBSERVED_CHRONIC_TRACE_NOT_DAY_FIT"),
    "Companion retained chronic trace requires the observed-only explicit scope marker.")
  g128i_check_rng(packed$rng_calls, keys, initial_state,
    if (is.null(observed_pack)) NULL else observed_pack$rng_calls, retained_chronic = retained)
  list(available = available, diagnostics = diagnostics,
    unchecked_fields = setdiff(names(packed), c(required, "scope")),
    rng_trace_scope = if (retained) "RETAINED_CHRONIC_SCHEMA_ONLY_NOT_COMPANION_EFFECTIVENESS" else
      if (companion) "NO_UNSIGNED_CALLS_FOR_COMPANION" else "DECLARED_UNSIGNED_SCHEMA_AND_INITIAL_STATE",
    limited_worker_fields_checked = if (retained) "worker_diagnostics.unsigned_trace_scope" else character())
}

g128i_authenticated_rds <- function(path, sha256) {
  size <- unname(file.info(path)$size)
  g128i_require(length(size) == 1L && is.finite(size) && size > 0,
    "Missing or empty expected archive.")
  bytes <- readBin(path, what = "raw", n = size)
  g128i_require(length(bytes) == size &&
    identical(digest::digest(bytes, algo = "sha256", serialize = FALSE), sha256),
    "Archive bytes do not match the explicit expected SHA-256.")
  # Decode the same authenticated buffer, not a second open of a mutable file.
  connection <- rawConnection(bytes, open = "rb")
  input <- if (length(bytes) >= 2L && identical(bytes[1:2], as.raw(c(31, 139))))
    gzcon(connection) else connection
  on.exit(close(input), add = TRUE)
  list(payload = readRDS(input), bytes = size)
}

g128i_check_archives <- function(expected, fixed_identity, observed_pack,
                                 required_ids, reference_kind, run_root,
                                 target_contrast = NULL) {
  required_columns <- c("file", "sha256", "allocation_id", "allocation_row_sha256")
  g128i_require(requireNamespace("digest", quietly = TRUE), "Existing digest package required.")
  g128i_require(is.numeric(required_ids) && length(required_ids) >= 2L &&
    !anyNA(required_ids) && identical(as.numeric(required_ids), as.numeric(seq_along(required_ids))),
    "Explicit required IDs must cover observed 1 and every consecutive required allocation.")
  total <- length(required_ids)
  g128i_require(identical(reference_kind, "sampled") || identical(reference_kind, "exhaustive"),
    "Reference kind must be sampled or exhaustive.")
  g128i_require(is.data.frame(expected) && !anyDuplicated(names(expected)) &&
    all(required_columns %in% names(expected)) && nrow(expected) == total &&
    is.numeric(expected$allocation_id) && !anyNA(expected$allocation_id) &&
    !anyDuplicated(expected$allocation_id) &&
    setequal(expected$allocation_id, required_ids), "Expected allocation coverage differs.")
  g128i_require(is.character(expected$file) && !anyNA(expected$file) &&
    all(nzchar(expected$file)) && all(grepl("^([A-Za-z]:[/\\\\]|/)", expected$file)) &&
    all(vapply(expected$sha256, g128i_hash_string, logical(1L))) &&
    all(vapply(expected$allocation_row_sha256, g128i_hash_string, logical(1L))) &&
    !anyDuplicated(expected$allocation_row_sha256), "Invalid explicit file/hash/label-row manifest.")
  identity_fields <- c("run_id", "cell_id", "reference_space",
    "input_manifest_sha256", "source_manifest_sha256")
  g128i_require(is.list(fixed_identity) && !anyDuplicated(names(fixed_identity)) &&
    identical(sort(names(fixed_identity)), sort(identity_fields)), "Explicit fixed identity fields required.")
  for (field in identity_fields[1:3]) g128i_require(is.character(fixed_identity[[field]]) &&
    length(fixed_identity[[field]]) == 1L && !is.na(fixed_identity[[field]]) &&
    nzchar(trimws(fixed_identity[[field]])), "Invalid fixed run/cell/reference identity.")
  for (field in identity_fields[4:5]) g128i_require(g128i_hash_string(fixed_identity[[field]]),
    "Invalid fixed input/source manifest hash.")
  root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
  g128i_require(grepl("^([Rr]:/Shaoyao_storage/cache/GSE128623_[A-Za-z0-9_-]+|/root/autodl-tmp/gse128623_label_reference_[A-Za-z0-9_-]+)$", root),
    "Use an explicit specific R-cache/cloud archive root, not a broad or C root.")
  paths <- vapply(expected$file, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
  g128i_require(!anyDuplicated(paths) && all(startsWith(paths, paste0(root, "/"))) &&
    all(grepl("[.]rds$", paths)), "Archive paths escape the run root or alias another expected file.")
  template <- g128i_family_template(target_contrast)
  initial_state <- g128i_native_initial_state()
  observed_check <- g128i_check_packed(observed_pack, template, initial_state,
    allow_retained_chronic_trace = !is.null(target_contrast))
  observed_t <- observed_pack$statistic
  observed_valid <- observed_check$available
  hits <- as.integer(observed_valid)
  missing <- integer(length(observed_t))
  all_diagnostics <- vector("list", total)
  receipts <- vector("list", total)
  unchecked <- vector("list", total)
  for (allocation in seq_len(total)) {
    row <- match(allocation, expected$allocation_id)
    readback <- g128i_authenticated_rds(paths[row], expected$sha256[row])
    payload <- readback$payload
    g128i_require(is.list(payload) && !anyDuplicated(names(payload)) &&
      identical(payload$schema_version, "gse128623-statistic-native-p-archive-v1") &&
      identical(payload$real_execution_authority_granted, FALSE), "Unexpected archive envelope/schema.")
    identity <- payload$identity
    g128i_require(is.list(identity) && !anyDuplicated(names(identity)) &&
      all(c(identity_fields, "allocation_id", "allocation_row_sha256") %in% names(identity)),
      "Incomplete archive identity.")
    for (field in identity_fields) g128i_require(identical(identity[[field]], fixed_identity[[field]]),
      paste("Archive identity differs:", field))
    g128i_require(is.numeric(identity$allocation_id) && length(identity$allocation_id) == 1L &&
      !is.na(identity$allocation_id) && identity$allocation_id == allocation &&
      identical(identity$allocation_row_sha256, expected$allocation_row_sha256[row]),
      "Archive allocation ID/label-row hash differs from the expected registry row.")
    if ("packed_sha256" %in% names(identity)) g128i_require(
      g128i_hash_string(identity$packed_sha256) &&
      identical(identity$packed_sha256, digest::digest(payload$packed, algo = "sha256")),
      "Archive packed-content self hash differs.")
    if ("packed_sha256" %in% names(expected)) g128i_require(
      g128i_hash_string(expected$packed_sha256[row]) &&
      identical(expected$packed_sha256[row], digest::digest(payload$packed, algo = "sha256")),
      "Expected packed-content hash differs.")
    checked <- g128i_check_packed(payload$packed, template, initial_state, observed_pack,
      allow_retained_chronic_trace = allocation == 1L && !is.null(target_contrast))
    if (allocation == 1L) g128i_require(identical(payload$packed, observed_pack),
      "Observed archive is not the explicitly supplied observed packed object.") else {
      valid <- checked$available
      extreme <- rep(FALSE, length(observed_t))
      finite <- which(observed_valid & is.finite(observed_t) & valid)
      extreme[finite] <- payload$packed$statistic[finite] >=
        observed_t[finite] - 1e-10 * (1 + abs(observed_t[finite]))
      infinite <- which(observed_valid & is.infinite(observed_t) & observed_t > 0 & valid)
      extreme[infinite] <- is.infinite(payload$packed$statistic[infinite]) &
        payload$packed$statistic[infinite] > 0
      hits <- hits + as.integer(extreme)
      missing <- missing + as.integer(!valid)
    }
    all_diagnostics[[allocation]] <- data.frame(allocation_id = allocation,
      checked$diagnostics, stringsAsFactors = FALSE)
    receipts[[allocation]] <- data.frame(allocation_id = allocation, file = paths[row],
      sha256 = expected$sha256[row], bytes = readback$bytes, stringsAsFactors = FALSE)
    unchecked[[allocation]] <- list(allocation_id = allocation,
      packed_fields_not_accepted = checked$unchecked_fields,
      envelope_fields_not_accepted = setdiff(names(payload), c("schema_version", "identity", "packed",
        "real_execution_authority_granted")),
      identity_fields_not_accepted = setdiff(names(identity), c(identity_fields, "allocation_id",
        "allocation_row_sha256", "packed_sha256")), rng_trace_scope = checked$rng_trace_scope,
      limited_worker_fields_checked = checked$limited_worker_fields_checked)
  }
  p <- rep(NA_real_, length(observed_t))
  completed <- observed_valid & missing == 0L
  p[completed] <- hits[completed] / total
  q <- rep(NA_real_, length(p))
  adjustment_n <- integer(nrow(template))
  for (index in seq_len(nrow(template))) {
    selected <- which(observed_pack$keys$family_id == template$family_id[index])
    q[selected] <- g128i_bh(p[selected])
    adjustment_n[index] <- sum(!is.na(p[selected]))
  }
  tail <- data.frame(observed_pack$keys, empirical_p = p, empirical_BH_q = q,
    hits_including_observed = hits, required_total_n = total, missing_null_n = missing,
    observed_unavailable = !observed_valid,
    available_total_n = total - missing - as.integer(!observed_valid),
    reference_kind = reference_kind, mc_uncertainty_computed = FALSE,
    stringsAsFactors = FALSE, check.names = FALSE)
  list(status = "PASS_EXPECTED_ARCHIVE_BYTES_IDENTITIES_DIAGNOSTICS_AND_TAILS",
    tail = tail, native_diagnostics = do.call(rbind, all_diagnostics),
    empirical_family_diagnostics = data.frame(template, n = vapply(template$family_id,
      function(family) sum(observed_pack$keys$family_id == family), integer(1L)),
      actual_BH_adjustment_n = adjustment_n, stringsAsFactors = FALSE),
    file_receipts = do.call(rbind, receipts), unchecked_fields = unchecked,
    checked_allocation_n = total, checked_family_n = nrow(template),
    mc_uncertainty = "NOT_IMPLEMENTED", RNG_final_consumption_verified = FALSE,
    mroast_auxiliary_acceptance = "ORDERED_KEYS_TYPES_RANGES_ONLY_NOT_NATIVE_REFIT",
    full_268_real_family_acceptance = FALSE, real_execution_authority_granted = FALSE)
}
