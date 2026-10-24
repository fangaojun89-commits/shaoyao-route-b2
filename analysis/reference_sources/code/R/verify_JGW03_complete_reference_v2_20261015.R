#!/usr/bin/env Rscript
# Saved-only full-space readback. No producer functions, native fits or RNG.
# The parent must additionally verify Rterm/Rscript/R are absent before launch.
jgwz_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwz_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
jgwz_digest <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 3L)
jgwz_counts <- function(x) c(n = length(x), NA_only = sum(is.na(x) & !is.nan(x)),
  "NaN" = sum(is.nan(x)), positive_Inf = sum(is.infinite(x) & x > 0),
  negative_Inf = sum(is.infinite(x) & x < 0), zero = sum(!is.na(x) & x == 0))
jgwz_invalid_p <- function(x) sum(!is.na(x) & (!is.finite(x) | x < 0 | x > 1))

jgwz_validate_compact <- function(draw, domain, key, check = jgwz_need) {
  check(identical(draw$schema, "jgw03-native-compact-v1") &&
    identical(draw$accession, domain$accession) && identical(draw$domain_sha256, domain$sha256) &&
    identical(draw$definition, domain$definition) && identical(draw$allocation_key, key) &&
    identical(names(draw$families), names(domain$families)) &&
    isTRUE(draw$native_table_numeric_and_dynamic_text_only) && identical(draw$DDS_or_fits_included, FALSE),
    "Full native compact identity/domain/order differs")
  if (identical(draw$execution_metadata$status, "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA")) {
    actual <- draw$execution_metadata$actual
    check(is.numeric(actual$seed) && identical(as.numeric(actual$seed), 20260817) &&
      identical(actual$RNGkind, c("Mersenne-Twister", "Inversion", "Rejection")) &&
      identical(actual$backend, "SerialParam") && identical(as.numeric(actual$workers), 1) &&
      identical(as.numeric(actual$mroast_nrot), 1999) && identical(as.numeric(actual$fgsea_eps), 0),
      "Saved caller-supplied fixed seed/backend/inner parameter fields differ")
  }
  totals <- c(rows = 0, numeric = 0, dynamic = 0)
  for (id in names(domain$families)) {
    spec <- domain$families[[id]]; value <- draw$families[[id]]; n <- length(spec$keys)
    check(is.character(spec$keys) && !anyNA(spec$keys) && all(nzchar(spec$keys)) && !anyDuplicated(spec$keys) &&
      identical(names(value$native_numeric), spec$numeric_names) &&
      identical(vapply(value$native_numeric, typeof, character(1)), spec$numeric_types) &&
      all(lengths(value$native_numeric) == n) && identical(names(value$native_text), spec$dynamic_names) &&
      all(lengths(value$native_text) == n) &&
      all(vapply(value$native_text, function(x) is.character(x) || is.factor(x) ||
        (is.logical(x) && all(is.na(x))), logical(1))), "Native key/vector/type schema differs")
    p <- value$native_numeric[[spec$map$native_P]]; q <- value$native_numeric[[spec$map$native_q]]
    check(is.numeric(p) && is.numeric(q) && jgwz_invalid_p(p) == 0L && jgwz_invalid_p(q) == 0L,
      "Native P/q invalid probability; unavailable values cannot be imputed")
    if (spec$map$method %in% c("CAMERA", "mroast1999_equivalent_not_historical_authentication")) {
      outer <- p; k <- which(!is.na(p)); outer[k] <- -log(p[k])
    } else {
      outer <- value$native_numeric[[spec$map$native_statistic]]
      check(is.numeric(outer), "Native statistic absent")
      if (spec$map$method != "edgeR_QL") outer <- abs(outer)
    }
    available <- is.finite(p) & p >= 0 & p <= 1 & !is.na(outer)
    check(identical(outer, value$outer_statistic) && identical(available, value$outer_available),
      "Independent outer-strength/availability arithmetic differs")
    for (column in names(value$native_numeric)) check(
      identical(jgwz_counts(value$native_numeric[[column]]), value$diagnostics$native_numeric_counts[[column]]),
      "Native NA/NaN/Inf/zero counts differ")
    for (column in names(value$native_text)) {
      x <- value$native_text[[column]]
      check(identical(c(n = length(x), unavailable = sum(is.na(x)),
        empty = sum(!is.na(x) & as.character(x) == "")), value$diagnostics$native_text_counts[[column]]),
        "Dynamic-text missingness differs")
    }
    check(identical(jgwz_counts(outer), value$diagnostics$outer_counts) &&
      value$diagnostics$outer_unavailable_n == sum(!available) &&
      value$diagnostics$raw_P_zero_n == sum(!is.na(p) & p == 0) &&
      value$diagnostics$raw_P_invalid_n == 0 && value$diagnostics$raw_q_invalid_n == 0 &&
      value$diagnostics$native_q_unavailable_n == sum(is.na(q)) &&
      value$diagnostics$native_q_le_005_known_n == sum(!is.na(q) & is.finite(q) & q >= 0 & q <= .05),
      "Native P/q and outer diagnostic arithmetic differs")
    totals <- totals + c(rows = n, numeric = sum(lengths(value$native_numeric)),
      dynamic = sum(lengths(value$native_text)))
  }
  totals
}

# A separate direct extended-real comparison, not the producer's tail helpers.
jgwz_increment <- function(observed, draw) {
  jgwz_need(is.numeric(observed$outer_statistic) && is.numeric(draw$outer_statistic) &&
    length(observed$outer_statistic) == length(draw$outer_statistic) &&
    is.logical(observed$outer_available) && is.logical(draw$outer_available) &&
    !anyNA(observed$outer_available) && !anyNA(draw$outer_available) &&
    length(observed$outer_available) == length(observed$outer_statistic) &&
    length(draw$outer_available) == length(draw$outer_statistic), "Tail comparison domain differs")
  comparison <- draw$outer_statistic >= observed$outer_statistic
  known <- observed$outer_available & draw$outer_available
  jgwz_need(!anyNA(comparison[known]), "Available extended-real comparison undefined")
  comparison[!known] <- FALSE
  list(hits = as.integer(comparison), unknown = as.integer(!draw$outer_available))
}

# Check the saved integer endpoints by acceptance and immediate-neighbor
# rejection. Monotonic hypergeometric tails identify the unique integer bounds;
# no producer binary search or interval function is used as an oracle.
jgwz_mc_endpoints <- function(b, lower, upper, n = 10000L, T = 369599L, confidence = .95) {
  plain_integer <- function(x) typeof(x) %in% c("integer", "double") && !is.object(x) &&
    is.null(dim(x)) && length(x) == 1L && is.finite(x) && x == floor(x)
  jgwz_need(all(vapply(list(b, n, T), plain_integer, logical(1))) && b >= 0 && b <= n && n >= 1 && n <= T &&
    T <= .Machine$integer.max - 1 && is.numeric(confidence) && !is.object(confidence) &&
    is.null(dim(confidence)) && length(confidence) == 1L && is.finite(confidence) && confidence > 0 && confidence < 1 &&
    is.double(lower) && is.double(upper) && !is.object(lower) && !is.object(upper) &&
    is.null(dim(lower)) && is.null(dim(upper)) && length(lower) == 1L && length(upper) == 1L &&
    is.finite(lower) && is.finite(upper), "Complete finite-population endpoint inputs required")
  lo <- round(lower * (T + 1) - 1); hi <- round(upper * (T + 1) - 1)
  jgwz_need(lo >= 0 && hi <= T && lo <= hi &&
    identical(lower, (1 + lo) / (T + 1)) && identical(upper, (1 + hi) / (T + 1)),
    "Saved P endpoints are not the integer-population transformation")
  if (n == T) {
    jgwz_need(lo == b && hi == b, "Census uncertainty must be an exact point")
  } else {
    cutoff <- log1p(-confidence) - log(2)
    cutoff <- cutoff - 8 * .Machine$double.eps * max(1, abs(cutoff))
    upper_tail <- function(M) stats::phyper(b - 1, M, T - M, n, lower.tail = FALSE, log.p = TRUE)
    lower_tail <- function(M) stats::phyper(b, M, T - M, n, lower.tail = TRUE, log.p = TRUE)
    jgwz_need(upper_tail(lo) >= cutoff && lower_tail(lo) >= cutoff &&
      upper_tail(hi) >= cutoff && lower_tail(hi) >= cutoff &&
      (lo == 0 || upper_tail(lo - 1) < cutoff) && (hi == T || lower_tail(hi + 1) < cutoff),
      "Hypergeometric endpoint or adjacent integer acceptance differs")
  }
  c(M_lower = lo, M_upper = hi)
}

jgwz_check_summary <- function(table, hits, unknown, observed_available, space,
    cache = new.env(parent = emptyenv()), check = jgwz_need) {
  check(space %in% c("four_group_MC", "high_companion_exact"), "Known complete space required")
  required <- if (space == "four_group_MC") 10000L else 19L
  check(is.integer(hits) && is.integer(unknown) && is.logical(observed_available) &&
    length(hits) > 0L && length(hits) == length(unknown) && length(hits) == length(observed_available) &&
    !anyNA(hits) && !anyNA(unknown) && !anyNA(observed_available) &&
    all(hits >= 0 & unknown >= 0 & hits + unknown <= required), "Complete fixed hit/missingness counts required")
  lower <- (1 + hits) / (required + 1)
  upper <- (1 + hits + unknown) / (required + 1)
  lower[!observed_available] <- upper[!observed_available] <- NA_real_
  point <- lower; point[unknown > 0L | !observed_available] <- NA_real_
  wanted <- data.frame(tail_hits_nonidentity = hits, null_unavailable_n = unknown,
    required_nonidentity_n = required, total_denominator = required + 1,
    observed_available = observed_available,
    observed_self_hit = ifelse(observed_available, 1L, NA_integer_), conditional_P = point,
    missingness_P_lower = lower, missingness_P_upper = upper,
    reference_complete_for_feature = observed_available & unknown == 0L, check.names = FALSE)
  # The frozen producer appends upper first through right-to-left chained
  # assignment, in both spaces. Match its schema without rewriting saved data.
  expected_names <- c(names(wanted), "finite_population_P_upper", "finite_population_P_lower")
  check(is.data.frame(table) && identical(names(table), expected_names) && nrow(table) == length(hits),
    "Complete saved summary column schema/order differs")
  for (column in names(wanted)) check(identical(table[[column]], wanted[[column]]),
    "Independent fixed-denominator summary field differs")
  eligible <- which(wanted$reference_complete_for_feature)
  missing <- which(!wanted$reference_complete_for_feature)
  check(is.double(table$finite_population_P_lower) && is.double(table$finite_population_P_upper) &&
    all(is.na(table$finite_population_P_lower[missing]) & !is.nan(table$finite_population_P_lower[missing])) &&
    all(is.na(table$finite_population_P_upper[missing]) & !is.nan(table$finite_population_P_upper[missing])),
    "Unavailable-feature uncertainty must remain NA, never imputed")
  if (space == "high_companion_exact") {
    check(identical(table$finite_population_P_lower[eligible], point[eligible]) &&
      identical(table$finite_population_P_upper[eligible], point[eligible]),
      "Complete companion interval must equal exact total hits over twenty")
  } else {
    # Cache only after checking each distinct b's actual saved endpoint pair.
    for (b in unique(hits[eligible])) {
      rows <- eligible[hits[eligible] == b]
      lo <- table$finite_population_P_lower[rows]; hi <- table$finite_population_P_upper[rows]
      check(length(unique(lo)) == 1L && length(unique(hi)) == 1L && !anyNA(lo) && !anyNA(hi),
        "Same count must retain the same finite-MC endpoint pair")
      key <- as.character(b); pair <- c(lo[1L], hi[1L])
      if (exists(key, cache, inherits = FALSE)) check(identical(get(key, cache), pair),
        "Cross-family finite-MC endpoint cache differs") else {
        jgwz_mc_endpoints(as.integer(b), pair[1L], pair[2L]); assign(key, pair, cache)
      }
    }
  }
  invisible(TRUE)
}

jgwz_main <- function(args) {
  jgwz_need(length(args) == 4L, "RUN_ROOT SPACE STAGE_EXIT_JSON NEW_R_OUTPUT required")
  root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE); space <- args[[2L]]
  stage_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
  out <- gsub("\\\\", "/", args[[4L]])
  jgwz_need(identical(root, "R:/Shaoyao_storage/cache/JGW03_label_reference_20261015r01") &&
    space %in% c("four_group_MC", "high_companion_exact") && startsWith(stage_path, paste0(root, "/")) &&
    startsWith(out, "R:/Shaoyao_storage/cache/") && !grepl("(^|/)\\.\\.(/|$)", out) &&
    !file.exists(out) && dir.exists(dirname(out)), "Explicit fixed root/space, stage and fresh R leaf required")
  checks <- 0; check <- function(ok, label) { jgwz_need(ok, label); checks <<- checks + 1 }
  check(identical(as.character(getRversion()), "4.6.1"), "Existing project R 4.6.1 required")
  pins <- new.env(hash = TRUE, parent = emptyenv())
  bind <- function(path, expected = NULL) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    actual <- list(path = path, bytes = unname(file.info(path)$size), sha256 = jgwz_sha(path))
    if (!is.null(expected)) check(is.numeric(expected$bytes) && length(expected$bytes) == 1L &&
      actual$bytes == expected$bytes && identical(actual$sha256, expected$sha256), "File bytes/SHA binding differs")
    if (exists(path, pins, inherits = FALSE)) check(identical(get(path, pins), actual), "Bound file changed")
    assign(path, actual, pins); actual
  }
  json <- function(path, expected = NULL) { bind(path, expected); jsonlite::read_json(path, simplifyVector = FALSE) }
  load <- function(path, expected = NULL) { bind(path, expected); readRDS(path) }
  self <- grep("^--file=", commandArgs(), value = TRUE)
  check(length(self) == 1L, "Explicit saved-only reader file required")
  self <- bind(sub("^--file=", "", self))
  stage <- json(stage_path); plan <- json(file.path(dirname(stage_path), "FIXED_STAGE_PLAN.json"))
  stage_name <- if (space == "high_companion_exact") "companions" else "remainder"
  worker_sha <- "f1de9c41003e0a96c6e08e3367a3e34f3a9ce898cce7a0be694bf9d4965ff4e5"
  finite_sha <- "349a5de42e47c41c7211c9f64c0f012d3658c72d840f5b75e4deb8912d52db6b"
  check(identical(stage$status, "FIXED_STAGE_ACTUAL_R_EXITS_ZERO_NOT_INDEPENDENT_ACCEPTANCE") &&
    identical(stage$stage, stage_name) && identical(plan$stage, stage_name) &&
    identical(plan$run_root, root) && plan$maximum_R_children == 1 &&
    stage$actual_stage_exit_code == 0 && stage$planned_R_processes == 6 &&
    length(stage$completed_processes) == 6L && length(plan$jobs) == 6L &&
    identical(stage$automatic_next_stage_launched, FALSE) && isTRUE(stage$single_R_worker) &&
    identical(stage$entry_sha256, worker_sha) && identical(plan$entry_sha256, worker_sha) &&
    identical(stage$launcher_sha256, plan$launcher_sha256), "Corresponding six actual stage exits required")
  bind(plan$entry, list(bytes = 16249, sha256 = worker_sha))
  launcher <- bind("<PROJECT_ROOT>/code/powershell/run_JGW03_full_reference_measured_20261015.ps1")
  check(identical(launcher$sha256, stage$launcher_sha256), "Actual stage launcher identity differs")
  cells <- c("GSE289929", "GSE227494", "GSE299063")
  for (j in seq_len(6L)) {
    cell <- cells[(j + 1L) %/% 2L]; mode <- if (j %% 2L == 1L) "run" else "finish"
    job <- plan$jobs[[j]]; row <- stage$completed_processes[[j]]; process <- json(row$process_receipt)
    check(row$job_number == j && identical(row$accession, cell) && row$actual_R_exit_code == 0 &&
      identical(job$cell, cell) && identical(job$mode, mode) && identical(job$space, space) &&
      process$actual_R_exit_code == 0 && process$job_number == j && identical(process$accession, cell) &&
      identical(process$mode, mode) && identical(process$space, space) &&
      identical(normalizePath(process$output_leaf, winslash = "/", mustWork = TRUE),
        file.path(root, paste0(cell, "_reference_01"))) &&
      identical(process$entry_sha256, worker_sha) && identical(process$launcher_sha256, stage$launcher_sha256) &&
      isTRUE(process$source_after_unchanged) && identical(process$startup_metadata_write_failed, FALSE) &&
      identical(process$raw_stdout_stderr_persisted, FALSE), "Actual cell/run/finish completion differs")
    if (mode == "run") {
      first <- if (space == "high_companion_exact") 1 else 6
      last <- if (space == "high_companion_exact") 20 else 10001
      check(job$first == first && job$last == last && process$first == first && process$last == last,
        "Complete fixed run range differs")
    }
  }
  # Only after the six process receipts pass may saved scientific compacts open.
  total_N <- if (space == "high_companion_exact") 20L else 10001L
  cache <- new.env(hash = TRUE, parent = emptyenv()); summaries <- list()
  for (accession in cells) {
    leaf <- file.path(root, paste0(accession, "_reference_01"))
    prepared <- json(file.path(leaf, "PREPARED.json")); context <- load(file.path(leaf, "PREPARED.rds"), prepared$prepared)
    check(identical(context$schema, "jgw03-first-nonidentity-preparation-v1") &&
      identical(context$accession, accession) && context$full_MC_nonidentity_draws == 10000L &&
      context$full_companion_including_observed == 20L && context$full_global_method_family_n == 89L &&
      context$full_companion_method_family_n == 34L, "Frozen full preparation scope differs")
    for (record in context$sources) bind(record$path, record)
    bind(context$entry$path, context$entry)
    for (record in context$observed_files) bind(record$path, record) # fixed file hashed only
    registry <- load(file.path(leaf, "FULL_LABEL_REGISTRY.rds"), context$label_registry)
    domain <- load(context$observed_files$domain$path, context$observed_files$domain)
    observed <- load(context$observed_files$observed$path, context$observed_files$observed)
    domain_body <- domain; domain_body$sha256 <- NULL
    full_families <- if (accession == "GSE289929") 31L else 29L
    check(identical(domain$schema, "jgw03-native-compact-domain-v1") && identical(domain$accession, accession) &&
      identical(jgwz_digest(domain_body), domain$sha256) && identical(domain$sha256, context$domain_sha256) &&
      length(domain$families) == full_families && identical(names(domain$observed_labels), domain$sample_ids) &&
      length(domain$sample_ids) == 12L && !anyDuplicated(domain$sample_ids), "Accepted fixed full domain differs")
    matrix <- if (space == "four_group_MC") registry$MC else registry$high_companion
    check(is.matrix(matrix) && identical(dim(matrix), c(total_N, 12L)) &&
      identical(colnames(matrix), domain$sample_ids), "Complete registry shape/sample order differs")
    keys <- unname(apply(matrix, 1L, paste0, collapse = ""))
    check(length(unique(keys)) == total_N && !anyNA(keys) && all(grepl("^[0-3]{12}$", keys)) &&
      identical(keys[1L], observed$allocation_key) &&
      identical(keys[1L], paste(as.integer(domain$observed_labels), collapse = "")) &&
      all(apply(matrix, 1L, function(x) identical(tabulate(x + 1L, 4L), rep(3L, 4L)))),
      "Distinct ordered four-by-three labels/identity differ")
    selected <- names(domain$families)
    if (space == "high_companion_exact") {
      nuisance <- which(!domain$observed_labels %in% c(0, 3))
      check(all(vapply(seq_len(total_N), function(i) identical(as.integer(matrix[i, nuisance]),
        as.integer(domain$observed_labels[nuisance])), logical(1))), "Companion nuisance roles changed")
      selected <- selected[vapply(domain$families, function(x) x$map$contrast == "high_vs_vehicle", logical(1))]
    }
    check(length(selected) == if (space == "four_group_MC") full_families else
      if (accession == "GSE289929") 12L else 11L, "Complete method-family selection differs")
    identity <- list(prepared_sha256 = prepared$prepared$sha256, registry = context$label_registry,
      domain_sha256 = domain$sha256, sources = context$sources, protected = context$protected,
      observed_files = context$observed_files, preparation_entry = context$entry, runtime = context$runtime,
      fixed_ram_sha256_is_checked_by_admission_before_every_new_fit = TRUE)
    identity_sha <- jgwz_digest(identity); store <- file.path(leaf, "full_reference_compacts_01")
    prefix <- load(file.path(store, sprintf("%s_TAIL_PREFIX_%05d.rds", space, total_N)))
    product <- load(file.path(store, paste0(space, "_TAIL_SUMMARY.rds")))
    check(identical(prefix$schema, "jgw03-fixed-registry-tail-prefix-v1") && identical(prefix$space, space) &&
      prefix$through_allocation_id == total_N && identical(prefix$scientific_identity_sha256, identity_sha) &&
      identical(prefix$full_registry, context$label_registry) && identical(prefix$observed, context$observed_files$observed) &&
      identical(prefix$provisional_P_computed, FALSE) &&
      length(prefix$completed_compact_pins) == total_N - 1L &&
      identical(names(prefix$completed_compact_pins), as.character(seq.int(2L, total_N))) &&
      identical(product$tail_checkpoint, prefix$state) && identical(product$scientific_identity_sha256, identity_sha) &&
      identical(product$full_registry, context$label_registry) && identical(product$execution_entry$sha256, worker_sha),
      "Full prefix/summary/archive lineage differs")
    state <- prefix$state; summary <- product$summary
    check(identical(state$schema, "jgw03-tail-accumulator-v1") && identical(state$space, space) &&
      identical(state$domain_sha256, domain$sha256) && state$total_N == total_N &&
      state$processed_nonidentity == total_N - 1L && identical(state$allocation_keys, keys) &&
      identical(state$definition, domain$definition) && identical(names(state$families), selected),
      "Full tail state domain/order/denominator differs")
    counters <- lapply(selected, function(id) {
      origin <- observed$families[[id]]
      check(all((origin$outer_statistic >= origin$outer_statistic)[origin$outer_available]),
        "Available observed value must self-hit including Inf ties")
      list(hits = integer(length(origin$outer_statistic)), unknown = integer(length(origin$outer_statistic)))
    }); names(counters) <- selected
    metadata_missing <- 0L; totals <- c(rows = 0, numeric = 0, dynamic = 0)
    for (i in seq_len(total_N)) {
      archive_leaf <- file.path(store, keys[i]); meta <- json(file.path(archive_leaf, "COMPLETE.json"))
      archive <- load(file.path(archive_leaf, "COMPACT.rds"), meta$archive)
      check(identical(meta$status, "IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT") &&
        identical(meta$allocation_key, keys[i]) && identical(meta$scientific_identity_sha256, identity_sha) &&
        identical(archive$schema, "jgw03-full-reference-allocation-v1") && identical(archive$accession, accession) &&
        identical(archive$scientific_identity_sha256, identity_sha) && identical(archive$allocation_key, keys[i]) &&
        identical(archive$compact_sha256, jgwz_digest(archive$compact)) &&
        identical(archive$execution_entry$sha256, worker_sha), "Saved allocation binding differs")
      draw <- archive$compact; totals <- totals + jgwz_validate_compact(draw, domain, keys[i], check)
      if (i == 1L) check(identical(draw, observed) && archive$new_native_calls == 0L,
        "Original label must exactly reuse accepted compact without fit") else {
        pointer <- prefix$completed_compact_pins[[as.character(i)]]
        check(identical(pointer$allocation_key, keys[i]) &&
          pointer$archive$bytes == meta$archive$bytes &&
          identical(pointer$archive$sha256, meta$archive$sha256), "Full ordered prefix archive pin differs")
        for (id in selected) {
          add <- jgwz_increment(observed$families[[id]], draw$families[[id]])
          counters[[id]]$hits <- counters[[id]]$hits + add$hits
          counters[[id]]$unknown <- counters[[id]]$unknown + add$unknown
        }
      }
      metadata_missing <- metadata_missing + as.integer(
        !identical(draw$execution_metadata$status, "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA"))
      rm(draw, archive); if (i %% 25L == 0L) invisible(gc())
    }
    check(identical(summary$schema, "jgw03-tail-summary-v1") && identical(summary$space, space) &&
      identical(summary$domain_sha256, domain$sha256) && summary$processed_nonidentity == total_N - 1L &&
      identical(names(summary$families), selected) && identical(summary$finite_mc_source_sha256, finite_sha) &&
      summary$execution_metadata_missing_allocations == metadata_missing &&
      state$execution_metadata_missing_allocations == metadata_missing &&
      identical(summary$execution_metadata_coverage_complete, metadata_missing == 0L) &&
      identical(summary$execution_metadata_is_caller_supplied_not_independently_authenticated, TRUE) &&
      identical(summary$native_P_q_changed, FALSE) && identical(summary$empirical_BH_computed, FALSE) &&
      identical(summary$uncertainty_scope, "Pointwise conditional finite-population MC only; no inner-MC/biological/FDR coverage") &&
      identical(summary$reference_acceptance, "Arithmetic only, not independent native or scientific acceptance"),
      "Complete summary metadata/coverage/scope differs")
    for (id in selected) {
      origin <- observed$families[[id]]; saved <- state$families[[id]]; count <- counters[[id]]
      check(identical(saved$observed, origin$outer_statistic) && identical(saved$observed_available, origin$outer_available) &&
        identical(saved$tail_hits_nonidentity, count$hits) && identical(saved$null_unavailable, count$unknown),
        "Independent full inclusive tail/unknown accumulator differs")
      jgwz_check_summary(summary$families[[id]], count$hits, count$unknown, origin$outer_available, space, cache, check)
    }
    summaries[[accession]] <- list(allocation_n = total_N, complete_space_family_n = length(selected),
      archived_full_native_family_n_per_draw = full_families, numeric_values_readback = totals[["numeric"]],
      dynamic_values_readback = totals[["dynamic"]], feature_rows_readback = totals[["rows"]],
      execution_metadata_missing_allocations = metadata_missing,
      tail_hits_nonidentity = sum(vapply(counters, function(x) sum(x$hits), numeric(1))),
      null_unavailable_feature_draws = sum(vapply(counters, function(x) sum(x$unknown), numeric(1))),
      original_label_exact_zero_fit_reuse = TRUE, scientific_identity_sha256 = identity_sha)
    rm(prefix, product, state, summary, counters, observed, registry, domain); invisible(gc())
  }
  paths <- sort(ls(pins, all.names = TRUE)); manifest <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    record <- get(paths[i], pins); check(unname(file.info(paths[i])$size) == record$bytes &&
      identical(jgwz_sha(paths[i]), record$sha256), "Bound input/source changed during full readback")
    manifest[[i]] <- data.frame(path = paths[i], bytes = record$bytes, sha256 = record$sha256)
  }
  dir.create(out)
  manifest_path <- file.path(out, "INPUT_PIN_MANIFEST.tsv.gz")
  connection <- gzfile(manifest_path, "wt")
  tryCatch(utils::write.table(do.call(rbind, manifest), connection, sep = "\t", quote = TRUE,
    row.names = FALSE), finally = close(connection))
  receipt <- list(status = "PASS_INDEPENDENT_SAVED_COMPLETE_SPACE_READBACK_ONLY", space = space,
    UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), checks = checks, source = self,
    runtime = list(R = as.character(getRversion()), stats = as.character(utils::packageVersion("stats")),
      digest = as.character(utils::packageVersion("digest")), jsonlite = as.character(utils::packageVersion("jsonlite"))),
    stage = get(normalizePath(stage_path, winslash = "/"), pins), corresponding_six_actual_exits_zero = TRUE,
    cells = summaries, bound_files = length(paths), all_bound_pins_unchanged = TRUE,
    pin_manifest = list(path = manifest_path, bytes = unname(file.info(manifest_path)$size), sha256 = jgwz_sha(manifest_path)),
    distinct_MC_hit_counts_with_independent_endpoint_checks = length(ls(cache, all.names = TRUE)),
    endpoint_oracle = "Integer endpoint/neighbor acceptance using log phyper and accepted conservative eight-epsilon guard; no producer interval/tail import",
    native_model_P_q_are_saved_readback_not_refitted = TRUE, new_native_fit_calls = 0L,
    native_P_q_or_BH_recomputed = FALSE, unavailable_values_imputed = FALSE, denominator_reduced = FALSE,
    finite_MC_scope = "Conditional fixed-per-assignment pointwise finite population only; no inner-MC/biological/simultaneous/FDR coverage",
    full_native_algorithm_or_biological_acceptance = FALSE, original_records_removed = 0L)
  jsonlite::write_json(receipt, file.path(out, "INDEPENDENT_COMPLETE_REFERENCE_READBACK.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
  cat(receipt$status, space, checks, "\n"); invisible(0L)
}
if (sys.nframe() == 0L) {
  result <- tryCatch(jgwz_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_COMPLETE_SPACE_READBACK_STOP_NO_RAW_CONDITION_PERSISTED"); 1L
  })
  quit(status = result)
}
