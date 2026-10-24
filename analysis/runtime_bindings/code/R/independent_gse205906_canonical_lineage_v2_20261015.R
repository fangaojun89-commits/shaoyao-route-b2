#!/usr/bin/env Rscript
# Read-only lineage/arithmetic verification. No edgeR/limma fit is called.
# Usage: Rscript --vanilla THIS.R --selftest
#        Rscript --vanilla THIS.R OLD_RUN NEW_RUN NEW_DIAGNOSTIC_DIRECTORY
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

require_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

identity <- function(path) {
  require_true(file.exists(path), paste("Missing input", basename(path)))
  list(path = normalizePath(path, mustWork = TRUE), bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

verify_identity <- function(expected) {
  actual <- identity(expected$path)
  require_true(actual$sha256 == expected$sha256 && actual$bytes == expected$bytes,
    paste("Frozen identity", basename(expected$path)))
  actual
}

near <- function(actual, expected) {
  length(actual) == length(expected) && all(is.finite(actual)) && all(is.finite(expected)) &&
    all(abs(actual - expected) <= 1e-12 * (1 + abs(expected)))
}

independent_bh <- function(probabilities) {
  count <- length(probabilities)
  ordering <- order(probabilities)
  adjusted <- pmin(1, rev(cummin(rev(probabilities[ordering] * count / seq_len(count)))))
  answer <- numeric(count)
  answer[ordering] <- adjusted
  answer
}

independent_labels <- function() {
  words <- as.matrix(expand.grid(rep(list(c("sham", "SNI", "PLM")), 9L),
    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
  words[rowSums(words == "sham") == 3L & rowSums(words == "SNI") == 3L &
    rowSums(words == "PLM") == 3L, , drop = FALSE]
}

signature <- function(labels) apply(labels, 1L, paste0, collapse = ":")

# First-appearance cells reconstruct the same ordering without model.matrix or
# the producer's order-of-first-position implementation.
independent_specification <- function(labels) {
  groups <- c("sham", "SNI", "PLM")
  cells <- unique(unname(labels))
  require_true(length(cells) == 3L && all(tabulate(match(labels, groups), 3L) == 3L),
    "Balanced named cells")
  design <- outer(unname(labels), cells, "==") * 1
  rownames(design) <- names(labels)
  colnames(design) <- paste0("cell_", seq_len(3L))
  scientific_weights <- list(
    local_disease = c(sham = -1, SNI = 1, PLM = 0),
    intervention = c(sham = 0, SNI = -1, PLM = 1))
  vectors <- matrix(0, 3L, 2L,
    dimnames = list(colnames(design), names(scientific_weights)))
  directions <- setNames(numeric(2L), names(scientific_weights))
  for (contrast in names(scientific_weights)) {
    weights <- unname(scientific_weights[[contrast]][cells])
    last <- max(which(weights != 0))
    directions[contrast] <- if (weights[last] > 0) 1 else -1
    vectors[, contrast] <- weights * directions[contrast]
    require_true(identical(as.numeric(design %*% vectors[, contrast] * directions[contrast]),
      as.numeric(scientific_weights[[contrast]][labels])), "Restored sample-wise scientific contrast")
  }
  list(design = design, contrasts = vectors, orientation = directions,
    permutation = match(cells, groups))
}

check_full_geometry <- function(labels) {
  keys <- signature(labels)
  require_true(identical(dim(labels), c(1680L, 9L)) && !anyDuplicated(keys) &&
    setequal(keys, signature(independent_labels())), "Complete independent balanced space")
  observed <- which(keys == paste(rep(c("sham", "SNI", "PLM"), each = 3L), collapse = ":"))
  companion <- list(local_disease = which(rowSums(labels[, 7:9, drop = FALSE] == "PLM") == 3L),
    intervention = which(rowSums(labels[, 1:3, drop = FALSE] == "sham") == 3L))
  require_true(length(observed) == 1L && all(lengths(companion) == 20L) &&
    identical(intersect(companion[[1L]], companion[[2L]]), observed), "Observed/20/20 spaces")
  partitions <- character(1680L)
  for (allocation_id in seq_len(1680L)) {
    specification <- independent_specification(labels[allocation_id, ])
    partitions[allocation_id] <- paste(max.col(specification$design), collapse = ":")
    for (contrast in names(companion)) {
      pair <- if (contrast == "local_disease") c("sham", "SNI") else c("SNI", "PLM")
      translation <- setNames(c("sham", "SNI", "PLM"), c("sham", "SNI", "PLM"))
      translation[pair] <- rev(pair)
      partner <- match(paste(unname(translation[labels[allocation_id, ]]), collapse = ":"), keys)
      require_true(!is.na(partner) && partner != allocation_id, "Swap partner identity")
      other <- independent_specification(labels[partner, ])
      require_true(identical(specification$design, other$design, num.eq = FALSE) &&
        identical(specification$contrasts[, contrast], other$contrasts[, contrast]) &&
        specification$orientation[contrast] == -other$orientation[contrast],
        "All canonical swap inputs and direction")
    }
  }
  require_true(length(unique(partitions)) == 280L && all(table(partitions) == 6L),
    "Full-space sixfold labelled multiplicity; no space reduction")
  list(observed = observed, companion = companion, partitions = 280L, partner_checks = 3360L)
}

selftest <- function() {
  labels <- independent_labels()
  colnames(labels) <- paste0("GSM", 6235292:6235300)
  geometry <- check_full_geometry(labels)
  require_true(near(independent_bh(c(.01, .01, .04, 1)), c(.02, .02, 4 * .04 / 3, 1)),
    "Independent BH ties")
  positive_zero <- c(0, 1)
  negative_zero <- positive_zero * c(-1, 1)
  require_true(!identical(positive_zero, negative_zero, num.eq = FALSE), "Bit-level field comparison")
  cat("PASS independent full-space algebra/geometry selftest; no model/input files read\n")
  invisible(geometry)
}

read_manifest <- function(path) {
  value <- data.table::fread(path, check.names = FALSE)
  require_true(identical(as.integer(value$allocation_id), seq_len(1680L)) &&
    !anyNA(value$sha256) && !anyDuplicated(value$allocation_id), "Complete ordered SHA manifest")
  value
}

assert_result <- function(value, expected_keys) {
  keys <- c("collection", "method", "contrast", "feature")
  numbers <- c("effect", "statistic", "p_value", "q_value")
  require_true(data.table::is.data.table(value) && identical(names(value), c(keys, numbers)) &&
    nrow(value) == 31872L && identical(as.data.frame(value[, keys, with = FALSE]), expected_keys) &&
    !anyNA(value), "Full combined native schema/keys")
  for (field in numbers) require_true(all(is.finite(value[[field]])), "Finite native fields")
  require_true(all(value$p_value > 0 & value$p_value <= 1 & value$q_value > 0 & value$q_value <= 1),
    "Native P/q range")
  allowed_negative <- value$collection == "genes" & value$method == "edgeR_robust_QL" &
    value$statistic >= -1e-10
  require_true(all(value$statistic >= 0 | allowed_negative), "Unchanged method-specific F roundoff bound")
}

repair_impact <- function(old, fresh, allocation_id) {
  rows <- lapply(c("local_disease", "intervention"), function(contrast) {
    position <- which(old$collection == "genes" & old$contrast == contrast)
    require_true(length(position) == 14217L &&
      all(fresh$collection[position] == "genes" & fresh$contrast[position] == contrast),
      "Full source-aligned gene family")
    previous <- old[position]
    current <- fresh[position]
    require_true(near(previous$q_value, independent_bh(previous$p_value)) &&
      near(current$q_value, independent_bh(current$p_value)), "Independent full-gene BH for old/new")
    data.frame(allocation_id = allocation_id, contrast = contrast, genes = length(position),
      maximum_absolute_effect_difference = max(abs(current$effect - previous$effect)),
      maximum_scaled_effect_difference = max(abs(current$effect - previous$effect) /
        (1 + pmax(abs(current$effect), abs(previous$effect)))),
      maximum_absolute_F_difference = max(abs(current$statistic - previous$statistic)),
      maximum_scaled_F_difference = max(abs(current$statistic - previous$statistic) /
        (1 + pmax(abs(current$statistic), abs(previous$statistic)))),
      maximum_absolute_log10_P_difference = max(abs(log10(current$p_value) - log10(previous$p_value))),
      maximum_absolute_log10_q_difference = max(abs(log10(current$q_value) - log10(previous$q_value))),
      old_BH_rejections = sum(previous$q_value < .05), new_BH_rejections = sum(current$q_value < .05),
      BH_status_changed_n = sum((previous$q_value < .05) != (current$q_value < .05)),
      minimum_new_F = min(current$statistic), negative_new_F_n = sum(current$statistic < 0))
  })
  do.call(rbind, rows)
}

main <- function(old_run, new_run, output_directory) {
  require_true(requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("digest", quietly = TRUE),
    "Existing read-only check dependencies")
  data.table::setDTthreads(1L)
  old_run <- normalizePath(old_run, mustWork = TRUE)
  new_run <- normalizePath(new_run, mustWork = TRUE)
  require_true(old_run != new_run && !file.exists(output_directory) && !dir.exists(output_directory),
    "Distinct immutable runs and new diagnostic destination")
  read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
  new_complete_path <- file.path(new_run, "ENUMERATION_COMPLETE.json")
  require_true(file.exists(new_complete_path), "Wait for complete canonical enumeration/summary")
  old_complete <- read_json(file.path(old_run, "ENUMERATION_COMPLETE.json"))
  new_complete <- read_json(new_complete_path)
  old_prepared <- read_json(file.path(old_run, "PREPARED.json"))
  new_prepared <- read_json(file.path(new_run, "PREPARED.json"))
  before <- lapply(c(file.path(old_run, c("PREPARED.json", "ENUMERATION_COMPLETE.json")),
    file.path(new_run, c("PREPARED.json", "ENUMERATION_COMPLETE.json"))), identity)
  for (prepared in list(old_prepared, new_prepared)) {
    for (field in c("derived", "observed", "worker", "definition", "input_manifest", "observed_checkpoint")) {
      before[[length(before) + 1L]] <- verify_identity(prepared[[field]])
    }
    inventory <- data.table::fread(prepared$input_manifest$path)
    require_true(!anyDuplicated(inventory$path), "Explicit unique preparation input inventory")
    for (row in seq_len(nrow(inventory))) {
      before[[length(before) + 1L]] <- verify_identity(as.list(inventory[row]))
    }
  }
  for (completed in list(old_complete, new_complete)) {
    for (item in completed$outputs) before[[length(before) + 1L]] <- verify_identity(item)
  }
  old_manifest <- read_manifest(old_complete$outputs$allocation_manifest$path)
  new_manifest <- read_manifest(new_complete$outputs$allocation_manifest$path)
  old_state <- readRDS(old_prepared$derived$path)
  new_state <- readRDS(new_prepared$derived$path)
  for (field in setdiff(names(old_state), "source_binding")) {
    require_true(identical(old_state[[field]], new_state[[field]], num.eq = FALSE),
      paste("Exact preserved derived payload", field))
  }
  require_true(identical(new_state$original_source_binding, old_state$source_binding) &&
    !identical(new_state$source_binding, old_state$source_binding) &&
    normalizePath(new_state$original_run, mustWork = TRUE) == old_run &&
    new_state$original_manifest$sha256 == identity(old_complete$outputs$allocation_manifest$path)$sha256 &&
    new_prepared$original_preparation$sha256 == identity(file.path(old_run, "PREPARED.json"))$sha256 &&
    isTRUE(new_prepared$canonical_gene_independent_refit) &&
    isTRUE(new_prepared$same_label_non_gene_checkpoint_reuse), "Explicit old/new lineage bindings")
  verify_identity(new_state$original_manifest)
  new_inputs <- data.table::fread(new_prepared$input_manifest$path)
  require_true(any(new_inputs$sha256 == old_prepared$derived$sha256) &&
    any(new_inputs$sha256 == new_state$original_manifest$sha256) &&
    any(new_inputs$sha256 == new_prepared$worker$sha256) &&
    any(new_inputs$sha256 == new_prepared$definition$sha256), "Binding payload includes source and implementation identities")
  # V2 metadata-only repair: a TSV roundtrip cannot recreate the complete R
  # serialization graph merely by coercing bytes to double. Reconstruct the
  # original typed construction from pinned files; first require every value
  # to equal the new manifest, then retain the exact digest assertion below.
  binding_base <- data.table::fread(old_prepared$input_manifest$path)
  binding_extra_paths <- list(original_prepared = file.path(old_run, "PREPARED.json"),
    original_derived = old_prepared$derived$path,
    original_allocation_manifest = file.path(old_run, "ALLOCATION_MANIFEST.tsv"),
    original_completion = file.path(old_run, "ENUMERATION_COMPLETE.json"),
    canonical_worker = file.path(new_run, "worker.R"),
    canonical_definition = file.path(new_run, "definition.md"),
    unchanged_summary_source = file.path(new_run, "summary.R"),
    unchanged_checker_source = file.path(new_run, "checker.R"))
  binding_extra_rows <- data.table::rbindlist(lapply(names(binding_extra_paths), function(role) {
    value <- identity(binding_extra_paths[[role]])
    data.table::data.table(role = role, path = value$path, bytes = value$bytes, sha256 = value$sha256)
  }))
  binding_inputs <- data.table::rbindlist(list(binding_base, binding_extra_rows))
  require_true(identical(names(binding_inputs), names(new_inputs)) &&
    nrow(binding_inputs) == nrow(new_inputs) &&
    all(vapply(names(new_inputs), function(field) {
      all(binding_inputs[[field]] == new_inputs[[field]])
    }, logical(1))), "Original binding constructor equals every pinned manifest value")
  recomputed_binding <- digest::digest(list(inputs = binding_inputs,
    previous_source_binding = old_state$source_binding,
    numerical_method = "canonical_per_allocation_independent_refit_v1"), algo = "sha256")
  require_true(identical(recomputed_binding, new_state$source_binding), "Independently reconstructed source-binding payload")
  geometry <- check_full_geometry(new_state$labels)
  require_true(geometry$observed == new_state$spaces$observed &&
    identical(as.integer(new_state$spaces$global), seq_len(1680L)) &&
    all(vapply(names(geometry$companion), function(contrast) {
      identical(geometry$companion[[contrast]], new_state$spaces[[contrast]])
    }, logical(1))), "Saved spaces versus independent geometry")
  labels <- new_state$labels
  expected_keys <- new_state$expected_keys
  old_binding <- old_state$source_binding
  new_binding <- new_state$source_binding
  rm(old_state, new_state, binding_inputs)
  invisible(gc())
  identities <- vector("list", 3360L)
  impact_rows <- vector("list", 1680L)
  transition_rows <- vector("list", 1680L)
  observed_exact <- NULL
  for (allocation_id in seq_len(1680L)) {
    old_path <- file.path(old_run, "allocations", sprintf("allocation_%04d.rds", allocation_id))
    new_path <- file.path(new_run, "allocations", sprintf("allocation_%04d.rds", allocation_id))
    old_pin <- identity(old_path)
    new_pin <- identity(new_path)
    require_true(old_pin$sha256 == old_manifest$sha256[allocation_id] &&
      old_pin$bytes == old_manifest$bytes[allocation_id] &&
      new_pin$sha256 == new_manifest$sha256[allocation_id] &&
      new_pin$bytes == new_manifest$bytes[allocation_id], "All paired old/new checkpoint SHA and bytes")
    identities[[2L * allocation_id - 1L]] <- old_pin
    identities[[2L * allocation_id]] <- new_pin
    previous <- readRDS(old_path)
    current <- readRDS(new_path)
    require_true(identical(previous$id, allocation_id) && identical(current$id, allocation_id) &&
      identical(previous$labels, labels[allocation_id, ]) && identical(current$labels, labels[allocation_id, ]) &&
      identical(previous$source_binding, old_binding) && identical(current$source_binding, new_binding) &&
      isTRUE(current$new_independent_gene_refit) && isTRUE(current$same_label_non_gene_exact),
      "All paired IDs, original labels, bindings and producer declarations")
    require_true(identical(current$original_allocation$path, old_pin$path) &&
      current$original_allocation$bytes == old_pin$bytes && current$original_allocation$sha256 == old_pin$sha256,
      "Per-checkpoint original allocation identity")
    specification <- independent_specification(labels[allocation_id, ])
    require_true(identical(names(current$canonical_specification), names(specification)), "Canonical specification schema")
    for (field in names(specification)) {
      # Numeric +0/-0 in contrast metadata denote the same exact coefficient;
      # bit identity remains mandatory separately for every reused data field.
      require_true(identical(current$canonical_specification[[field]], specification[[field]]),
        paste("All saved canonical metadata versus independent reconstruction", field))
    }
    assert_result(previous$result, expected_keys)
    assert_result(current$result, expected_keys)
    non_gene <- which(previous$result$collection != "genes")
    require_true(length(non_gene) == 3438L, "Full signed/unsigned non-gene denominator")
    for (field in names(previous$result)) {
      require_true(identical(previous$result[[field]][non_gene], current$result[[field]][non_gene], num.eq = FALSE),
        paste("Bit-exact same-label non-gene field", field))
    }
    expected_impact <- repair_impact(previous$result, current$result, allocation_id)
    recorded_impact <- as.data.frame(current$repair_impact)
    require_true(identical(names(recorded_impact), names(expected_impact)) && nrow(recorded_impact) == 2L &&
      setequal(recorded_impact$contrast, expected_impact$contrast), "Complete recorded repair-impact schema")
    recorded_impact <- recorded_impact[match(expected_impact$contrast, recorded_impact$contrast), , drop = FALSE]
    for (field in names(expected_impact)) {
      if (is.numeric(expected_impact[[field]])) {
        require_true(near(recorded_impact[[field]], expected_impact[[field]]), paste("Repair impact", field))
      } else require_true(identical(recorded_impact[[field]], expected_impact[[field]]), paste("Repair impact", field))
    }
    impact_rows[[allocation_id]] <- expected_impact
    transitions <- lapply(expected_impact$contrast, function(contrast) {
      position <- which(previous$result$collection == "genes" & previous$result$contrast == contrast)
      old_selected <- previous$result$q_value[position] < .05
      new_selected <- current$result$q_value[position] < .05
      data.frame(allocation_id = allocation_id, contrast = contrast,
        newly_below_q005_n = sum(!old_selected & new_selected),
        no_longer_below_q005_n = sum(old_selected & !new_selected))
    })
    transition_rows[[allocation_id]] <- do.call(rbind, transitions)
    if (allocation_id == geometry$observed) {
      observed_exact <- setNames(vapply(c("effect", "statistic", "p_value", "q_value"), function(field) {
        identical(previous$result[[field]], current$result[[field]], num.eq = FALSE)
      }, logical(1)), c("effect", "statistic", "p_value", "q_value"))
      require_true(all(observed_exact), "Actual exact observed bridge; no accumulated compatibility tolerance")
    }
    if (allocation_id %% 300L == 0L) cat("independent lineage", allocation_id, "of 1680\n")
  }
  impact <- do.call(rbind, impact_rows)
  transitions <- do.call(rbind, transition_rows)
  require_true(nrow(impact) == 3360L && nrow(transitions) == 3360L &&
    all(impact$BH_status_changed_n == transitions$newly_below_q005_n + transitions$no_longer_below_q005_n) &&
    all(impact$new_BH_rejections - impact$old_BH_rejections ==
      transitions$newly_below_q005_n - transitions$no_longer_below_q005_n), "All directional BH transition closure")
  aggregate <- do.call(rbind, lapply(c("local_disease", "intervention"), function(contrast) {
    index <- which(impact$contrast == contrast)
    data.frame(contrast = contrast, allocations = length(index), genes_per_allocation = 14217L,
      feature_allocation_pairs = length(index) * 14217L,
      allocations_with_BH_status_change = sum(impact$BH_status_changed_n[index] > 0),
      changed_feature_allocation_pairs = sum(impact$BH_status_changed_n[index]),
      newly_below_q005_feature_allocation_pairs = sum(transitions$newly_below_q005_n[index]),
      no_longer_below_q005_feature_allocation_pairs = sum(transitions$no_longer_below_q005_n[index]),
      old_below_q005_feature_allocation_pairs = sum(impact$old_BH_rejections[index]),
      new_below_q005_feature_allocation_pairs = sum(impact$new_BH_rejections[index]),
      maximum_absolute_gene_effect_difference = max(impact$maximum_absolute_effect_difference[index]),
      maximum_scaled_gene_effect_difference = max(impact$maximum_scaled_effect_difference[index]),
      maximum_absolute_gene_F_difference = max(impact$maximum_absolute_F_difference[index]),
      maximum_scaled_gene_F_difference = max(impact$maximum_scaled_F_difference[index]),
      maximum_absolute_log10_P_difference = max(impact$maximum_absolute_log10_P_difference[index]),
      maximum_absolute_log10_q_difference = max(impact$maximum_absolute_log10_q_difference[index]))
  }))
  invisible(lapply(c(before, identities), verify_identity))
  require_true(dir.create(output_directory, recursive = TRUE), "New diagnostic directory")
  output_directory <- normalizePath(output_directory, mustWork = TRUE)
  impact_path <- file.path(output_directory, "INDEPENDENT_REPAIR_IMPACT.tsv")
  transition_path <- file.path(output_directory, "INDEPENDENT_BH_TRANSITIONS.tsv")
  aggregate_path <- file.path(output_directory, "INDEPENDENT_REPAIR_IMPACT_SUMMARY.tsv")
  data.table::fwrite(impact, impact_path, sep = "\t", quote = TRUE, na = "NA")
  data.table::fwrite(transitions, transition_path, sep = "\t", quote = TRUE, na = "NA")
  data.table::fwrite(aggregate, aggregate_path, sep = "\t", quote = TRUE, na = "NA")
  this_script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
  receipt <- list(status = "PASS_INDEPENDENT_CANONICAL_LINEAGE_AND_REPAIR_IMPACT",
    completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), checker = identity(this_script),
    old_completion = identity(file.path(old_run, "ENUMERATION_COMPLETE.json")),
    new_completion = identity(new_complete_path), old_preparation = identity(file.path(old_run, "PREPARED.json")),
    new_preparation = identity(file.path(new_run, "PREPARED.json")),
    old_run = old_run, new_run = new_run,
    old_worker = old_prepared$worker, new_worker = new_prepared$worker,
    old_definition = old_prepared$definition, new_definition = new_prepared$definition,
    old_source_binding = old_binding, new_source_binding = new_binding,
    independently_reconstructed_source_binding = recomputed_binding,
    source_binding_reconstruction = "Pinned old manifest plus eight role-ordered typed pin rows, rbindlist; all values equal new pinned manifest before exact digest test",
    allocation_pair_n = 1680L, allocation_SHA_checks = 3360L,
    same_label_non_gene_rows_per_allocation = 3438L, all_non_gene_fields_bit_exact = TRUE,
    full_combined_rows_per_allocation = 31872L, new_gene_family_n = 14217L,
    complete_old_new_gene_BH_families_checked = 6720L,
    canonical_partner_checks = geometry$partner_checks, full_labelled_space_preserved = TRUE,
    observed_fields_bit_exact = as.list(observed_exact), repair_impact_rows = nrow(impact),
    BH_transition_feature_allocation_total = sum(impact$BH_status_changed_n),
    maximum_scaled_gene_F_difference = max(impact$maximum_scaled_F_difference),
    maximum_absolute_log10_P_difference = max(impact$maximum_absolute_log10_P_difference),
    maximum_absolute_log10_q_difference = max(impact$maximum_absolute_log10_q_difference),
    outputs = list(repair_impact = identity(impact_path), BH_transitions = identity(transition_path),
      repair_impact_summary = identity(aggregate_path)),
    all_before_after_pins_unchanged = TRUE, models_or_new_reference_fits_executed = FALSE,
    source_binding_and_complete_metadata_checked_not_an_independent_refit = TRUE,
    producer_scoped_numerical_acceptance = new_complete$scoped_numerical_acceptance,
    full_space_arithmetic_checker_is_separate = TRUE, camera_warnings_resolved = FALSE,
    authenticated_design_exchangeability = FALSE, overall_scientific_provenance = "NEEDS REVIEW",
    scope = "Canonical lineage, exact same-label non-gene preservation and recorded repair impact only; no biological validation, authenticated type-I/FDR control, primary replacement or project completion")
  jsonlite::write_json(receipt, file.path(output_directory, "INDEPENDENT_CANONICAL_LINEAGE_CHECK.json"),
    pretty = TRUE, auto_unbox = TRUE, digits = 17, na = "null")
  capture.output(sessionInfo(), file = file.path(output_directory, "sessionInfo.txt"))
  cat(receipt$status, "; all 1680 old/new allocation pairs checked\n")
}

arguments <- commandArgs(trailingOnly = TRUE)
if (identical(arguments, "--selftest")) selftest() else {
  require_true(length(arguments) == 3L, "Use --selftest OR explicit old/new run and new diagnostic paths")
  main(arguments[1L], arguments[2L], arguments[3L])
}
