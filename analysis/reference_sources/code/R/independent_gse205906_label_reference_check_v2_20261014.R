#!/usr/bin/env Rscript
# Independent arithmetic check of completed reference objects; never fits models.
# Usage: Rscript --vanilla THIS.R --selftest
#        Rscript --vanilla THIS.R RUN_DIRECTORY NEW_DIAGNOSTIC_DIRECTORY
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

require_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

file_identity <- function(path) {
  require_true(file.exists(path), paste("Missing input", basename(path)))
  list(path = normalizePath(path, mustWork = TRUE),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

verify_identity <- function(identity) {
  actual <- file_identity(identity$path)
  require_true(identical(actual$sha256, identity$sha256),
    paste("Input identity", basename(identity$path)))
  if (!is.null(identity$bytes)) require_true(actual$bytes == identity$bytes,
    paste("Input bytes", basename(identity$path)))
  actual
}

near <- function(actual, expected, tolerance = 1e-12) {
  length(actual) == length(expected) && all(is.finite(actual)) &&
    all(is.finite(expected)) &&
    all(abs(actual - expected) <= tolerance * (1 + abs(expected)))
}

# Sorted monotone adjustment, independent of the producer's p.adjust call.
independent_bh <- function(probabilities) {
  count <- length(probabilities)
  ordering <- order(probabilities)
  adjusted_sorted <- pmin(1, rev(cummin(rev(
    probabilities[ordering] * count / seq_len(count)))))
  answer <- numeric(count)
  answer[ordering] <- adjusted_sorted
  answer
}

# Explicit type-7 interpolation rather than the producer's quantile invocation.
linear_quantiles <- function(values, probabilities) {
  ordered <- sort(values)
  positions <- 1 + (length(ordered) - 1) * probabilities
  lower <- floor(positions)
  fraction <- positions - lower
  ordered[lower] + fraction *
    (ordered[pmin(length(ordered), lower + 1L)] - ordered[lower])
}

label_signature <- function(labels) apply(labels, 1L, paste0, collapse = ":")

# Independently enumerate 3^9 words and select balanced words, not nested combn.
independent_label_space <- function() {
  words <- as.matrix(expand.grid(rep(list(c("sham", "SNI", "PLM")), 9L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  words[rowSums(words == "sham") == 3L & rowSums(words == "SNI") == 3L &
    rowSums(words == "PLM") == 3L, , drop = FALSE]
}

valid_statistic <- function(statistic, collection, method) {
  is.finite(statistic) & (statistic >= 0 |
    (collection == "genes" & method == "edgeR_robust_QL" & statistic >= -1e-10))
}

synthetic_selftest <- function() {
  label_space <- independent_label_space()
  require_true(nrow(label_space) == 1680L &&
    !anyDuplicated(label_signature(label_space)), "Synthetic balanced space")
  require_true(sum(rowSums(label_space[, 7:9, drop = FALSE] == "PLM") == 3L) == 20L &&
    sum(rowSums(label_space[, 1:3, drop = FALSE] == "sham") == 3L) == 20L,
    "Synthetic companion spaces")
  require_true(near(independent_bh(c(.01, .01, .04, 1)), c(.02, .02, 4 * .04 / 3, 1)),
    "Synthetic BH ties")
  require_true(near(linear_quantiles(c(8, 0, 4, 2), c(0, .25, .5, .95, 1)),
    c(0, 1.5, 3, 7.4, 8)), "Synthetic interpolation")
  statistics <- rbind(c(0, -3e-12, 0, 1), c(2, 2, 1, 3))
  observed <- statistics[, 1L]
  tolerant <- rowSums(statistics >= observed - 1e-10 * (1 + abs(observed)))
  require_true(identical(as.integer(tolerant), c(4L, 3L)), "Synthetic tolerant tails")
  require_true(near(tolerant / 4, (1 + rowSums(
    statistics[, -1L, drop = FALSE] >= observed - 1e-10 * (1 + abs(observed)))) / 4),
    "Synthetic exhaustive/plus-one identity")
  require_true(all(valid_statistic(c(-3e-12, 0), "genes", "edgeR_robust_QL")) &&
    !valid_statistic(-2e-10, "genes", "edgeR_robust_QL") &&
    !valid_statistic(-3e-12, "PROGENy", "ULM_robust_limma"),
    "Synthetic method-specific floating boundary")
  array_counts <- unlist(jsonlite::fromJSON("[20,20]", simplifyVector = FALSE), use.names = TRUE)
  object_counts <- unlist(jsonlite::fromJSON(
    '{"local_disease":20,"intervention":20}', simplifyVector = FALSE), use.names = TRUE)
  require_true(is.null(names(array_counts)) && identical(as.numeric(array_counts), c(20, 20)) &&
    identical(names(object_counts), c("local_disease", "intervention")) &&
    identical(as.numeric(object_counts), c(20, 20)), "Synthetic pairwise JSON schemas")
  cat("PASS independent arithmetic selftest; no scientific inputs read\n")
}

compare_tables <- function(actual, expected, key_columns, label) {
  require_true(setequal(names(actual), names(expected)), paste(label, "columns"))
  actual <- as.data.frame(actual[, names(expected), drop = FALSE])
  expected <- as.data.frame(expected)
  make_keys <- function(table) do.call(paste, c(table[key_columns], sep = "\034"))
  actual_keys <- make_keys(actual)
  expected_keys <- make_keys(expected)
  require_true(nrow(actual) == nrow(expected) && !anyDuplicated(actual_keys) &&
    !anyDuplicated(expected_keys) && setequal(actual_keys, expected_keys),
    paste(label, "complete key domain"))
  actual <- actual[match(expected_keys, actual_keys), , drop = FALSE]
  for (column in names(expected)) {
    if (is.numeric(expected[[column]])) {
      require_true(near(actual[[column]], expected[[column]]), paste(label, column))
    } else {
      require_true(identical(as.character(actual[[column]]), as.character(expected[[column]])),
        paste(label, column))
    }
  }
  list(rows = nrow(expected), columns = ncol(expected),
    compared_values = nrow(expected) * ncol(expected))
}

main <- function(run_path, output_path) {
  require_true(requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("jsonlite", quietly = TRUE) &&
    requireNamespace("digest", quietly = TRUE), "Existing check dependencies required")
  data.table::setDTthreads(1L)
  run_path <- normalizePath(run_path, mustWork = TRUE)
  require_true(!dir.exists(output_path) && !file.exists(output_path),
    "Diagnostic directory must be new")
  completion_path <- file.path(run_path, "ENUMERATION_COMPLETE.json")
  require_true(file.exists(completion_path), "Await completed summary; no partial result checks")
  read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
  read_tsv <- function(path) as.data.frame(data.table::fread(path, na.strings = "NA",
    check.names = FALSE))
  completion <- read_json(completion_path)
  prepared_path <- file.path(run_path, "PREPARED.json")
  prepared <- read_json(prepared_path)
  preparation_pins <- list(prepared$derived, prepared$observed, prepared$worker,
    prepared$definition, prepared$input_manifest, prepared$observed_checkpoint)
  before_pins <- c(list(file_identity(completion_path), file_identity(prepared_path)),
    lapply(preparation_pins, verify_identity), lapply(completion$outputs, verify_identity),
    list(verify_identity(completion$summary_source)))
  verify_identity(completion$preparation)
  input_manifest <- read_tsv(prepared$input_manifest$path)
  require_true(all(c("path", "sha256") %in% names(input_manifest)) &&
    !anyDuplicated(input_manifest$path), "Explicit source input inventory")
  source_pins <- lapply(seq_len(nrow(input_manifest)), function(index) {
    verify_identity(as.list(input_manifest[index, , drop = FALSE]))
  })
  before_pins <- c(before_pins, source_pins)
  continuation_paths <- list.files(run_path,
    pattern = "^CONTINUATION_[1-4]_STARTED\\.json$", full.names = TRUE)
  continuation_sources <- list()
  for (path in continuation_paths) {
    continuation <- read_json(path)
    require_true(continuation$status == "VALIDATION_ONLY_CONTINUATION" &&
      isTRUE(continuation$all_model_kernels_unchanged) &&
      continuation$data_table_threads == 1L &&
      continuation$frozen_inferential_worker$sha256 == prepared$worker$sha256 &&
      continuation$derived$sha256 == prepared$derived$sha256,
      "Validation-only continuation identity and declared boundary")
    continuation_sources[[length(continuation_sources) + 1L]] <-
      verify_identity(continuation$continuation_source)
    before_pins <- c(before_pins, list(file_identity(path)),
      list(verify_identity(continuation$frozen_inferential_worker)),
      list(verify_identity(continuation$derived)),
      list(verify_identity(continuation$continuation_source)))
  }
  # Only structural state fields are used; no counts, E, scores, or model calls.
  state <- readRDS(prepared$derived$path)
  labels <- state$labels
  spaces <- state$spaces
  expected_keys <- state$expected_keys
  source_binding <- state$source_binding
  rm(state)
  invisible(gc())
  require_true(identical(dim(labels), c(1680L, 9L)), "Stored label matrix")
  independent_labels <- independent_label_space()
  signatures <- label_signature(labels)
  require_true(!anyDuplicated(signatures) &&
    setequal(signatures, label_signature(independent_labels)), "All legal allocations exactly once")
  observed_labels <- rep(c("sham", "SNI", "PLM"), each = 3L)
  observed_id <- which(signatures == paste(observed_labels, collapse = ":"))
  companion <- list(
    local_disease = which(rowSums(labels[, 7:9, drop = FALSE] == "PLM") == 3L),
    intervention = which(rowSums(labels[, 1:3, drop = FALSE] == "sham") == 3L))
  require_true(length(observed_id) == 1L && observed_id == prepared$observed_id &&
    identical(as.integer(spaces$global), seq_len(1680L)) &&
    identical(as.integer(spaces$observed), as.integer(observed_id)) &&
    all(vapply(names(companion), function(name) {
      length(companion[[name]]) == 20L &&
        identical(as.integer(spaces[[name]]), as.integer(companion[[name]]))
    }, logical(1))) && identical(intersect(companion[[1L]], companion[[2L]]), observed_id),
    "Observed and companion identities")
  allocation_table <- read_tsv(file.path(run_path, "ALLOCATIONS.tsv"))
  require_true(nrow(allocation_table) == 1680L && ncol(allocation_table) == 10L &&
    identical(as.integer(allocation_table[[1L]]), seq_len(1680L)) &&
    identical(unname(as.matrix(allocation_table[, -1L])), unname(labels)),
    "Allocation TSV/state agreement")
  before_pins <- c(before_pins, list(file_identity(file.path(run_path, "ALLOCATIONS.tsv"))))
  filter_proof_path <- file.path(run_path, "FILTER_INVARIANCE.tsv")
  filter_proof <- read_tsv(filter_proof_path)
  require_true(nrow(filter_proof) == 1680L &&
    identical(as.integer(filter_proof$allocation_id), seq_len(1680L)) &&
    all(filter_proof$retained_gene_n == 14217L) && all(filter_proof$same_mask),
    "All-allocation retained-mask preparation receipt")
  before_pins <- c(before_pins, list(file_identity(filter_proof_path)))
  observed <- read_tsv(prepared$observed$path)
  key_columns <- c("collection", "method", "contrast", "feature")
  result_columns <- c(key_columns, "effect", "statistic", "p_value", "q_value")
  require_true(identical(names(observed), result_columns) &&
    identical(observed[key_columns], expected_keys) && !anyNA(observed) &&
    setequal(unique(observed$contrast), names(companion)), "Observed schema and full keys")
  family_labels <- do.call(paste, c(observed[c("collection", "method", "contrast")], sep = "::"))
  family_names <- unique(family_labels)
  family_indices <- lapply(family_names, function(name) which(family_labels == name))
  require_true(length(family_indices) == 10L && length(unique(observed$collection)) == 5L &&
    sum(observed$collection == "genes") == 2L * prepared$retained_gene_n &&
    prepared$library_n == 9L && prepared$retained_gene_n == 14217L &&
    prepared$signed_program_n == 95L && prepared$all_labels_filter_verified == 1680L,
    "Fixed source-wide collection denominators")
  manifest <- read_tsv(completion$outputs$allocation_manifest$path)
  require_true(nrow(manifest) == 1680L &&
    identical(as.integer(manifest$allocation_id), seq_len(1680L)), "All allocation manifests")
  feature_n <- nrow(observed)
  # Two dense matrices plus a largest-family working slice stay below 4 GiB.
  estimated_peak_bytes <- 8 * 1680 * (2 * feature_n + 2 * max(lengths(family_indices)))
  require_true(estimated_peak_bytes < 3.2 * 1024^3, "Independent matrix memory bound")
  statistic_matrix <- matrix(NA_real_, feature_n, 1680L)
  probability_matrix <- matrix(NA_real_, feature_n, 1680L)
  bh_rejections <- matrix(0L, length(family_indices), 1680L)
  allocation_pins <- vector("list", 1680L)
  negative_counts <- integer(1680L)
  minima <- numeric(1680L)
  for (allocation_id in seq_len(1680L)) {
    path <- file.path(run_path, "allocations", sprintf("allocation_%04d.rds", allocation_id))
    allocation_pins[[allocation_id]] <- file_identity(path)
    require_true(allocation_pins[[allocation_id]]$sha256 == manifest$sha256[allocation_id] &&
      allocation_pins[[allocation_id]]$bytes == manifest$bytes[allocation_id], "Manifest checkpoint pin")
    checkpoint <- readRDS(path)
    result <- checkpoint$result
    require_true(identical(checkpoint$id, allocation_id) &&
      identical(checkpoint$labels, labels[allocation_id, ]) &&
      identical(checkpoint$source_binding, source_binding) &&
      data.table::is.data.table(result) && identical(names(result), result_columns) &&
      identical(as.data.frame(result[, key_columns, with = FALSE]), expected_keys) &&
      !anyNA(result), "Allocation identity/source/domain")
    require_true(all(is.finite(result$effect)) &&
      all(valid_statistic(result$statistic, result$collection, result$method)) &&
      all(is.finite(result$p_value) & result$p_value > 0 & result$p_value <= 1) &&
      all(is.finite(result$q_value) & result$q_value > 0 & result$q_value <= 1),
      "Allocation finite and method-specific statistic/P/q boundaries")
    negative_counts[allocation_id] <- sum(result$statistic < 0)
    minima[allocation_id] <- min(result$statistic)
    require_true(all(abs(result$p_value[result$statistic < 0] - 1) <= 1e-12),
      "Near-zero negative gene QLF has native P=1")
    require_true(near(checkpoint$elapsed_seconds, manifest$elapsed_seconds[allocation_id]) &&
      negative_counts[allocation_id] == manifest$negative_native_statistic_n[allocation_id] &&
      near(minima[allocation_id], manifest$minimum_native_statistic[allocation_id]),
      "Allocation diagnostic manifest fields")
    for (family_id in seq_along(family_indices)) {
      index <- family_indices[[family_id]]
      independent_q <- independent_bh(result$p_value[index])
      require_true(near(result$q_value[index], independent_q), "Complete native BH family")
      bh_rejections[family_id, allocation_id] <- sum(independent_q < .05)
    }
    statistic_matrix[, allocation_id] <- result$statistic
    probability_matrix[, allocation_id] <- result$p_value
    if (allocation_id == observed_id) {
      for (field in c("effect", "statistic", "p_value", "q_value")) {
        require_true(near(result[[field]], observed[[field]]), paste("Observed checkpoint", field))
      }
    }
    if (allocation_id %% 300L == 0L) cat("independently loaded", allocation_id, "of 1680\n")
  }
  rm(checkpoint, result)
  invisible(gc())
  feature_expected <- observed
  feature_expected$tail_global_raw <- feature_expected$tail_global_tolerant <- integer(feature_n)
  feature_expected$tail_pairwise_raw <- feature_expected$tail_pairwise_tolerant <- integer(feature_n)
  chunk_groups <- split(seq_len(feature_n), ceiling(seq_len(feature_n) / 384L))
  for (index in chunk_groups) {
    values <- statistic_matrix[index, , drop = FALSE]
    threshold <- observed$statistic[index]
    tolerant_threshold <- threshold - 1e-10 * (1 + abs(threshold))
    feature_expected$tail_global_raw[index] <- rowSums(values >= threshold)
    feature_expected$tail_global_tolerant[index] <- rowSums(values >= tolerant_threshold)
    for (contrast in names(companion)) {
      position <- which(observed$contrast[index] == contrast)
      target <- index[position]
      feature_expected$tail_pairwise_raw[target] <- rowSums(
        values[position, companion[[contrast]], drop = FALSE] >= threshold[position])
      feature_expected$tail_pairwise_tolerant[target] <- rowSums(
        values[position, companion[[contrast]], drop = FALSE] >= tolerant_threshold[position])
    }
  }
  feature_expected$global_N <- 1680L
  feature_expected$pairwise_N <- 20L
  feature_expected$exact_P_global <- feature_expected$tail_global_tolerant / 1680
  feature_expected$exact_P_pairwise <- feature_expected$tail_pairwise_tolerant / 20
  feature_expected$candidate_gate_eligible <- FALSE
  feature_expected$empirical_BH_global <- feature_expected$empirical_BH_pairwise <- numeric(feature_n)
  for (index in family_indices) {
    feature_expected$empirical_BH_global[index] <- independent_bh(feature_expected$exact_P_global[index])
    feature_expected$empirical_BH_pairwise[index] <- independent_bh(feature_expected$exact_P_pairwise[index])
  }
  feature_check <- compare_tables(read_tsv(completion$outputs$features$path),
    feature_expected, key_columns, "All feature tails/empirical P/BH")
  require_true(all(feature_expected$tail_global_tolerant >= 1L) &&
    all(feature_expected$tail_pairwise_tolerant >= 1L), "Observed included once in every tail")
  allocation_expected_rows <- list()
  summary_expected_rows <- list()
  qq_probabilities <- c(.001, .01, .05, .10, .25, .50, .75, .90, .99)
  for (family_id in seq_along(family_indices)) {
    index <- family_indices[[family_id]]
    context <- observed[index[1L], c("collection", "method", "contrast"), drop = FALSE]
    values <- probability_matrix[index, , drop = FALSE]
    quantiles <- vapply(seq_len(1680L), function(column) {
      linear_quantiles(values[, column], qq_probabilities)
    }, numeric(length(qq_probabilities)))
    metrics <- data.frame(allocation_id = seq_len(1680L),
      collection = context$collection, method = context$method, contrast = context$contrast,
      features = length(index), tail_01 = colMeans(values < .01),
      tail_05 = colMeans(values < .05), tail_10 = colMeans(values < .10),
      BH_rejections = bh_rejections[family_id, ],
      any_BH_rejection = bh_rejections[family_id, ] > 0L,
      minimum_P = apply(values, 2L, min))
    for (probability_id in seq_along(qq_probabilities)) {
      metrics[[paste0("P_quantile_", qq_probabilities[probability_id])]] <- quantiles[probability_id, ]
    }
    for (space in c("global", "pairwise_third_group_fixed")) {
      allocations <- if (space == "global") seq_len(1680L) else companion[[context$contrast]]
      space_metrics <- metrics[allocations, , drop = FALSE]
      space_metrics$space <- space
      allocation_expected_rows[[length(allocation_expected_rows) + 1L]] <- space_metrics
      median_01 <- linear_quantiles(space_metrics$tail_01, .5)
      median_05 <- linear_quantiles(space_metrics$tail_05, .5)
      median_10 <- linear_quantiles(space_metrics$tail_10, .5)
      warning_levels <- sum(c(median_01 > .02, median_05 > .10, median_10 > .20))
      summary_expected_rows[[length(summary_expected_rows) + 1L]] <- data.frame(
        space = space, collection = context$collection, method = context$method,
        contrast = context$contrast, allocations = length(allocations), feature_n = length(index),
        median_tail_01 = median_01, median_tail_05 = median_05, median_tail_10 = median_10,
        p05_tail_05 = linear_quantiles(space_metrics$tail_05, .05),
        p95_tail_05 = linear_quantiles(space_metrics$tail_05, .95),
        median_BH_rejections = linear_quantiles(space_metrics$BH_rejections, .5),
        maximum_BH_rejections = max(space_metrics$BH_rejections),
        fraction_allocations_any_BH = mean(space_metrics$any_BH_rejection),
        warning_levels = warning_levels, conditional_tail_warning = warning_levels >= 2L,
        interpretation = "CONDITIONAL_RANDOM_LABEL_REFERENCE_NOT_AUTHENTICATED_TYPE_I_OR_FDR")
    }
    rm(values, quantiles)
    invisible(gc())
  }
  allocation_expected <- do.call(rbind, allocation_expected_rows)
  summary_expected <- do.call(rbind, summary_expected_rows)
  allocation_check <- compare_tables(read_tsv(completion$outputs$allocation_families$path),
    allocation_expected, c("allocation_id", "space", "collection", "method", "contrast"),
    "All allocation-family diagnostics/QQ")
  summary_check <- compare_tables(read_tsv(completion$outputs$family_summary$path),
    summary_expected, c("space", "collection", "method", "contrast"), "All family summaries")
  swap_rows <- list()
  family_asymmetry_rows <- list()
  for (contrast in names(companion)) {
    permutation <- if (contrast == "local_disease") c(sham = "SNI", SNI = "sham", PLM = "PLM") else
      c(sham = "sham", SNI = "PLM", PLM = "SNI")
    swapped_labels <- matrix(unname(permutation[labels]), nrow(labels), ncol(labels))
    partners <- match(label_signature(swapped_labels), signatures)
    require_true(!anyNA(partners) && identical(as.integer(partners[partners]), seq_len(1680L)) &&
      all(partners != seq_len(1680L)) && setequal(partners[companion[[contrast]]], companion[[contrast]]),
      "Independent two-sided swap involution and companion closure")
    first <- which(seq_len(1680L) < partners)
    index <- which(observed$contrast == contrast)
    maxima <- numeric(length(first))
    above <- integer(length(first))
    contrast_families <- unique(family_labels[index])
    family_diagnostics <- setNames(lapply(contrast_families, function(name) {
      list(above = 0L, pair_maxima = numeric(length(first)),
        pair_above = rep(FALSE, length(first)), absolute_maximum = 0,
        negative_n = 0L, minimum_native = Inf)
    }), contrast_families)
    for (block in split(index, ceiling(seq_along(index) / 384L))) {
      left <- statistic_matrix[block, first, drop = FALSE]
      right <- statistic_matrix[block, partners[first], drop = FALSE]
      scaled <- abs(left - right) / (1 + pmax(abs(left), abs(right)))
      maxima <- pmax(maxima, apply(scaled, 2L, max))
      above <- above + colSums(scaled > 1e-10)
      # Read only the matrices already loaded for the accepted arithmetic check.
      # These diagnostic strata do not alter the pre-fixed symmetry gate.
      for (name in unique(family_labels[block])) {
        position <- which(family_labels[block] == name)
        family_scaled <- scaled[position, , drop = FALSE]
        diagnostic <- family_diagnostics[[name]]
        diagnostic$above <- diagnostic$above + sum(family_scaled > 1e-10)
        diagnostic$pair_maxima <- pmax(diagnostic$pair_maxima, apply(family_scaled, 2L, max))
        diagnostic$pair_above <- diagnostic$pair_above | colSums(family_scaled > 1e-10) > 0L
        diagnostic$absolute_maximum <- max(diagnostic$absolute_maximum,
          abs(left[position, , drop = FALSE] - right[position, , drop = FALSE]))
        diagnostic$negative_n <- diagnostic$negative_n +
          sum(left[position, , drop = FALSE] < 0) + sum(right[position, , drop = FALSE] < 0)
        diagnostic$minimum_native <- min(diagnostic$minimum_native,
          left[position, , drop = FALSE], right[position, , drop = FALSE])
        family_diagnostics[[name]] <- diagnostic
      }
    }
    swap_rows[[contrast]] <- data.frame(allocation_id = first,
      partner_id = partners[first], contrast = contrast,
      maximum_scaled_difference = maxima, above_fixed_tie_tolerance_n = above)
    for (name in contrast_families) {
      family_index <- which(family_labels == name)
      context <- observed[family_index[1L], c("collection", "method", "contrast"), drop = FALSE]
      diagnostic <- family_diagnostics[[name]]
      pair_quantiles <- linear_quantiles(diagnostic$pair_maxima, c(.5, .95, .99))
      family_asymmetry_rows[[length(family_asymmetry_rows) + 1L]] <- data.frame(
        collection = context$collection, method = context$method, contrast = contrast,
        features = length(family_index), allocations = 1680L, swap_pairs = length(first),
        feature_pair_comparisons = length(family_index) * length(first),
        above_fixed_tolerance_n = diagnostic$above,
        fraction_above_fixed_tolerance = diagnostic$above / (length(family_index) * length(first)),
        swap_pairs_above_fixed_tolerance_n = sum(diagnostic$pair_above),
        maximum_scaled_difference = max(diagnostic$pair_maxima),
        median_pair_maximum_scaled_difference = pair_quantiles[1L],
        p95_pair_maximum_scaled_difference = pair_quantiles[2L],
        p99_pair_maximum_scaled_difference = pair_quantiles[3L],
        maximum_absolute_statistic_difference = diagnostic$absolute_maximum,
        negative_native_statistic_n = diagnostic$negative_n,
        minimum_native_statistic = diagnostic$minimum_native,
        global_tolerance_changed_feature_n = sum(
          feature_expected$tail_global_raw[family_index] != feature_expected$tail_global_tolerant[family_index]),
        pairwise_tolerance_changed_feature_n = sum(
          feature_expected$tail_pairwise_raw[family_index] != feature_expected$tail_pairwise_tolerant[family_index]),
        interpretation = "READ_ONLY_NUMERICAL_DIAGNOSTIC_NOT_A_GATE_RELAXATION_OR_BIOLOGICAL_VALIDATION")
    }
  }
  swap_expected <- do.call(rbind, swap_rows)
  family_asymmetry <- do.call(rbind, family_asymmetry_rows)
  require_true(nrow(family_asymmetry) == 10L &&
    sum(family_asymmetry$above_fixed_tolerance_n) == sum(swap_expected$above_fixed_tie_tolerance_n) &&
    near(max(family_asymmetry$maximum_scaled_difference), max(swap_expected$maximum_scaled_difference)) &&
    sum(family_asymmetry$negative_native_statistic_n) == sum(negative_counts),
    "Method-stratified diagnostic closure")
  swap_check <- compare_tables(read_tsv(completion$outputs$label_swaps$path), swap_expected,
    c("allocation_id", "partner_id", "contrast"), "All two-sided swaps")
  resolution_status <- if (any(feature_expected$exact_P_pairwise < .10) ||
    any(feature_expected$exact_P_global < 1 / 840)) "NEEDS_NUMERICAL_SYMMETRY_REVIEW" else
    "TWO_SIDED_RESOLUTION_RESPECTED"
  numerical_acceptance <- resolution_status == "TWO_SIDED_RESOLUTION_RESPECTED" &&
    sum(swap_expected$above_fixed_tie_tolerance_n) == 0L
  # jsonlite auto_unbox serializes the producer's named length-two R vector as
  # [20,20]. Accept that exact positional schema or the equivalent named object;
  # the independently reconstructed, contrast-labelled spaces remain decisive.
  completion_pair_counts <- unlist(completion$pairwise_allocations, use.names = TRUE)
  require_true(is.numeric(completion_pair_counts) && length(completion_pair_counts) == 2L &&
    all(is.finite(completion_pair_counts)), "Pairwise completion count schema")
  completion_pair_schema <- if (is.null(names(completion_pair_counts))) "array_local_disease_intervention" else
    "named_contrast_object"
  if (is.null(names(completion_pair_counts))) {
    names(completion_pair_counts) <- c("local_disease", "intervention")
  }
  require_true(setequal(names(completion_pair_counts), names(companion)),
    "Pairwise completion contrast labels")
  require_true(completion$allocations == 1680L &&
    completion_pair_counts[["local_disease"]] == 20L &&
    completion_pair_counts[["intervention"]] == 20L &&
    identical(completion$scoped_numerical_acceptance, numerical_acceptance) &&
    identical(completion$two_sided_resolution_status, resolution_status) &&
    near(completion$swap_maximum_scaled_difference, max(swap_expected$maximum_scaled_difference)) &&
    completion$swap_above_tolerance_comparisons == sum(swap_expected$above_fixed_tie_tolerance_n) &&
    completion$gene_QLF_roundoff_negative_count == sum(negative_counts) &&
    near(completion$minimum_native_statistic, min(minima)) &&
    completion$global_conditional_warning_family_n == sum(
      summary_expected$conditional_tail_warning[summary_expected$space == "global"]) &&
    isTRUE(completion$observed_counted_once) && isTRUE(completion$no_valid_allocation_omitted) &&
    identical(completion$candidate_gate_eligible, FALSE) &&
    identical(completion$authenticated_design_exchangeability, FALSE) &&
    identical(completion$other_sources_and_method_controls_complete, FALSE) &&
    identical(completion$primary_results_modified, FALSE), "Complete receipt arithmetic and scope")
  require_true(completion$status == if (numerical_acceptance)
    "ENUMERATION_COMPLETE_SCOPED_CONDITIONAL_REFERENCE" else
    "NUMERICAL_REVIEW_REQUIRED_ENUMERATION_PRESERVED", "Complete receipt status")
  # Recheck every original/source/product byte identity after the read-only check.
  invisible(lapply(c(before_pins, allocation_pins), verify_identity))
  require_true(dir.create(output_path, recursive = TRUE), "New diagnostic directory")
  output_path <- normalizePath(output_path, mustWork = TRUE)
  data.table::fwrite(summary_expected, file.path(output_path, "INDEPENDENT_FAMILY_SUMMARY.tsv"),
    sep = "\t", quote = TRUE, na = "NA")
  asymmetry_path <- file.path(output_path, "INDEPENDENT_METHOD_SWAP_DIAGNOSTIC.tsv")
  data.table::fwrite(family_asymmetry, asymmetry_path, sep = "\t", quote = TRUE, na = "NA")
  this_script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
  receipt <- list(status = "PASS_INDEPENDENT_COMPLETED_REFERENCE_ARITHMETIC",
    completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    checker = file_identity(this_script), run = run_path,
    completed_reference = file_identity(completion_path), preparation = file_identity(prepared_path),
    output_pins = completion$outputs, external_input_n = nrow(input_manifest),
    validation_only_continuation_receipt_n = length(continuation_paths),
    actual_continuation_source_pins = continuation_sources,
    allocation_n = 1680L, observed_id = observed_id,
    companion_n = c(local_disease = 20L, intervention = 20L),
    completion_pair_count_JSON_schema = completion_pair_schema,
    full_feature_n = feature_n, full_native_family_n = length(family_indices),
    expected_peak_dense_bytes = estimated_peak_bytes,
    feature_check = feature_check, allocation_family_check = allocation_check,
    family_summary_check = summary_check, swap_check = swap_check,
    all_native_BH_families_checked = 1680L * length(family_indices),
    gene_QLF_negative_n = sum(negative_counts), minimum_native_statistic = min(minima),
    method_stratified_swap_diagnostic = file_identity(asymmetry_path),
    diagnostic_does_not_relax_numerical_gate = TRUE,
    two_sided_resolution_status = resolution_status,
    producer_scoped_numerical_acceptance = numerical_acceptance,
    all_source_and_allocation_pins_unchanged = TRUE,
    native_values_clipped_or_changed = FALSE, models_or_reference_fits_executed = FALSE,
    authenticated_design_exchangeability = FALSE, candidate_gate_eligible = FALSE,
    overall_scientific_provenance = "NEEDS REVIEW",
    scope = "Completed fixed-nine-library conditional label reference only; no independent biological validation, authenticated type-I/FDR control, candidate admission, other-source control completion, or whole-project acceptance")
  jsonlite::write_json(receipt, file.path(output_path, "INDEPENDENT_LABEL_REFERENCE_CHECK.json"),
    pretty = TRUE, auto_unbox = TRUE, digits = 17, na = "null")
  capture.output(sessionInfo(), file = file.path(output_path, "sessionInfo.txt"))
  cat(receipt$status, ";", feature_n, "feature rows; all 1680 allocations checked\n")
}

arguments <- commandArgs(trailingOnly = TRUE)
if (identical(arguments, "--selftest")) {
  synthetic_selftest()
} else {
  require_true(length(arguments) == 2L, "Use --selftest OR explicit run and new diagnostic paths")
  main(arguments[1L], arguments[2L])
}
