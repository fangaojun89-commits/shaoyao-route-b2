#!/usr/bin/env Rscript
# Complete-space, read-only verification of cloud allocation objects and summaries.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
pin <- function(path) list(path = normalizePath(path, mustWork = TRUE),
  bytes = unname(file.info(path)$size), sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
jwrite <- function(value, path) {
  need(!file.exists(path), paste("Existing summary", basename(path)))
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, digits = 17, na = "null")
}
twrite <- function(value, path) {
  need(!file.exists(path), paste("Existing summary", basename(path)))
  data.table::fwrite(value, path, sep = "\t", quote = TRUE, na = "NA")
}
read_table <- function(path) data.table::fread(path, na.strings = "NA", check.names = FALSE)
args <- commandArgs(trailingOnly = TRUE)
need(length(args) == 1L, "Explicit run directory required")
run <- normalizePath(args[1], mustWork = TRUE)
prepared <- jsonlite::fromJSON(file.path(run, "PREPARED.json"), simplifyVector = FALSE)
for (item in list(prepared$derived, prepared$observed, prepared$worker, prepared$definition,
  prepared$input_manifest, prepared$observed_checkpoint))
  need(identical(pin(item$path)$sha256, item$sha256), "Frozen preparation pin")
inputs <- read_table(prepared$input_manifest$path)
for (i in seq_len(nrow(inputs))) need(pin(inputs$path[i])$sha256 == inputs$sha256[i], "Frozen external source")
state <- readRDS(prepared$derived$path)
observed <- read_table(prepared$observed$path)
keys <- c("collection", "method", "contrast", "feature")
need(identical(as.data.frame(observed[, ..keys]), state$expected_keys), "Observed key domain")
observed_checkpoint <- readRDS(prepared$observed_checkpoint$path)
for (field in c("effect", "statistic", "p_value", "q_value")) need(
  all(abs(observed[[field]] - observed_checkpoint$result[[field]]) <= 1e-12 * (1 + abs(observed[[field]]))),
  paste("Pinned observed checkpoint value", field))
family <- paste(observed$collection, observed$method, observed$contrast, sep = "::")
families <- split(seq_len(nrow(observed)), family)
tail_raw_global <- tail_tol_global <- tail_raw_pair <- tail_tol_pair <- integer(nrow(observed))
allocation_rows <- list(); manifests <- list(); swap_rows <- list()
qq_probs <- c(.001, .01, .05, .10, .25, .50, .75, .90, .99)
tolerance <- 1e-10 * (1 + abs(observed$statistic))
label_keys <- apply(state$labels, 1, paste, collapse = ":")
read_allocation <- function(id) {
  path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", id))
  need(file.exists(path), paste("Missing legal allocation", id))
  value <- readRDS(path); result <- value$result
  need(identical(value$id, id) && identical(value$labels, state$labels[id, ]) &&
    identical(value$source_binding, state$source_binding), "Allocation source/labels")
  need(data.table::is.data.table(result) && identical(names(result), names(observed)) &&
    identical(as.data.frame(result[, ..keys]), state$expected_keys) && !anyNA(result), "Complete allocation schema/domain")
  for (field in c("effect", "statistic", "p_value", "q_value")) need(all(is.finite(result[[field]])), "Finite allocation")
  permitted_roundoff <- result$collection == "genes" & result$method == "edgeR_robust_QL" &
    result$statistic >= -1e-10
  need(all((result$statistic >= 0 | permitted_roundoff) & result$p_value > 0 & result$p_value <= 1 &
    result$q_value > 0 & result$q_value <= 1), "Allocation statistic/P/q range")
  for (index in families) need(all(abs(p.adjust(result$p_value[index], "BH") - result$q_value[index]) <= 1e-12), "Allocation BH family")
  value
}
for (i in state$spaces$global) {
  value <- read_allocation(i); result <- value$result
  path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", i))
  identity <- pin(path)
  manifests[[i]] <- data.table::data.table(allocation_id = i, bytes = identity$bytes,
    sha256 = identity$sha256, elapsed_seconds = value$elapsed_seconds,
    negative_native_statistic_n = sum(result$statistic < 0), minimum_native_statistic = min(result$statistic))
  raw <- result$statistic >= observed$statistic
  tolerant <- result$statistic >= observed$statistic - tolerance
  tail_raw_global <- tail_raw_global + raw
  tail_tol_global <- tail_tol_global + tolerant
  paired <- (observed$contrast == "local_disease" & i %in% state$spaces$local_disease) |
    (observed$contrast == "intervention" & i %in% state$spaces$intervention)
  tail_raw_pair <- tail_raw_pair + (raw & paired)
  tail_tol_pair <- tail_tol_pair + (tolerant & paired)
  for (index in families) {
    spaces <- "global"
    contrast <- observed$contrast[index[1]]
    if (i %in% state$spaces[[contrast]]) spaces <- c(spaces, "pairwise_third_group_fixed")
    p <- result$p_value[index]; q <- result$q_value[index]
    quantiles <- stats::quantile(p, qq_probs, names = FALSE, type = 7)
    for (space in spaces) {
      row <- data.table::data.table(allocation_id = i, space = space,
        collection = observed$collection[index[1]], method = observed$method[index[1]],
        contrast = contrast, features = length(index), tail_01 = mean(p < .01),
        tail_05 = mean(p < .05), tail_10 = mean(p < .10), BH_rejections = sum(q < .05),
        any_BH_rejection = any(q < .05), minimum_P = min(p))
      for (k in seq_along(qq_probs)) row[[paste0("P_quantile_", qq_probs[k])]] <- quantiles[k]
      allocation_rows[[length(allocation_rows) + 1L]] <- row
    }
  }
  for (contrast in c("local_disease", "intervention")) {
    swapped <- state$labels[i, ]
    pair <- if (contrast == "local_disease") c("sham", "SNI") else c("SNI", "PLM")
    swapped[state$labels[i, ] == pair[1]] <- pair[2]
    swapped[state$labels[i, ] == pair[2]] <- pair[1]
    partner_id <- match(paste(swapped, collapse = ":"), label_keys)
    need(!is.na(partner_id), "Swap closure")
    if (i < partner_id) {
      partner <- read_allocation(partner_id)$result
      index <- which(observed$contrast == contrast)
      difference <- abs(result$statistic[index] - partner$statistic[index]) /
        (1 + pmax(abs(result$statistic[index]), abs(partner$statistic[index])))
      swap_rows[[length(swap_rows) + 1L]] <- data.table::data.table(allocation_id = i,
        partner_id = partner_id, contrast = contrast, maximum_scaled_difference = max(difference),
        above_fixed_tie_tolerance_n = sum(difference > 1e-10))
    }
  }
  if (i %% 200 == 0L) cat("verified", i, "of 1680\n")
}
need(all(tail_tol_global >= 1L & tail_tol_pair >= 1L), "Observed counted in every reference")
feature <- data.table::copy(observed)
feature[, `:=`(tail_global_raw = tail_raw_global, tail_global_tolerant = tail_tol_global,
  tail_pairwise_raw = tail_raw_pair, tail_pairwise_tolerant = tail_tol_pair,
  global_N = 1680L, pairwise_N = 20L,
  exact_P_global = tail_tol_global / 1680, exact_P_pairwise = tail_tol_pair / 20,
  candidate_gate_eligible = FALSE)]
feature[, empirical_BH_global := p.adjust(exact_P_global, "BH"), by = .(collection, method, contrast)]
feature[, empirical_BH_pairwise := p.adjust(exact_P_pairwise, "BH"), by = .(collection, method, contrast)]
allocation <- data.table::rbindlist(allocation_rows)
summary <- allocation[, .(allocations = .N, feature_n = unique(features),
  median_tail_01 = median(tail_01), median_tail_05 = median(tail_05), median_tail_10 = median(tail_10),
  p05_tail_05 = as.numeric(quantile(tail_05, .05)), p95_tail_05 = as.numeric(quantile(tail_05, .95)),
  median_BH_rejections = median(BH_rejections), maximum_BH_rejections = max(BH_rejections),
  fraction_allocations_any_BH = mean(any_BH_rejection)), by = .(space, collection, method, contrast)]
summary[, warning_levels := (median_tail_01 > .02) + (median_tail_05 > .10) + (median_tail_10 > .20)]
summary[, conditional_tail_warning := warning_levels >= 2L]
summary[, interpretation := "CONDITIONAL_RANDOM_LABEL_REFERENCE_NOT_AUTHENTICATED_TYPE_I_OR_FDR"]
swap <- data.table::rbindlist(swap_rows)
need(nrow(swap) == 1680L, "Both contrast swap pair inventories")
resolution_status <- if (any(feature$exact_P_pairwise < .10) || any(feature$exact_P_global < 1/840))
  "NEEDS_NUMERICAL_SYMMETRY_REVIEW" else "TWO_SIDED_RESOLUTION_RESPECTED"
numerical_acceptance <- resolution_status == "TWO_SIDED_RESOLUTION_RESPECTED" &&
  sum(swap$above_fixed_tie_tolerance_n) == 0L
for (item in list(prepared$derived, prepared$observed, prepared$worker, prepared$definition,
  prepared$input_manifest, prepared$observed_checkpoint))
  need(identical(pin(item$path)$sha256, item$sha256), "Frozen preparation pin after enumeration")
for (i in seq_len(nrow(inputs))) need(pin(inputs$path[i])$sha256 == inputs$sha256[i], "External source unchanged after enumeration")
twrite(feature, file.path(run, "FEATURE_EMPIRICAL_REFERENCE.tsv"))
twrite(allocation, file.path(run, "ALLOCATION_FAMILY_DIAGNOSTICS.tsv"))
twrite(summary, file.path(run, "FAMILY_REFERENCE_SUMMARY.tsv"))
twrite(swap, file.path(run, "LABEL_SWAP_NUMERICAL_CHECK.tsv"))
twrite(data.table::rbindlist(manifests), file.path(run, "ALLOCATION_MANIFEST.tsv"))
outputs <- setNames(lapply(file.path(run, c("FEATURE_EMPIRICAL_REFERENCE.tsv", "ALLOCATION_FAMILY_DIAGNOSTICS.tsv",
  "FAMILY_REFERENCE_SUMMARY.tsv", "LABEL_SWAP_NUMERICAL_CHECK.tsv", "ALLOCATION_MANIFEST.tsv")), pin),
  c("features", "allocation_families", "family_summary", "label_swaps", "allocation_manifest"))
jwrite(list(status = if (numerical_acceptance) "ENUMERATION_COMPLETE_SCOPED_CONDITIONAL_REFERENCE" else
    "NUMERICAL_REVIEW_REQUIRED_ENUMERATION_PRESERVED", scoped_numerical_acceptance = numerical_acceptance,
  allocations = 1680L,
  pairwise_allocations = c(local_disease = 20L, intervention = 20L),
  completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  preparation = pin(file.path(run, "PREPARED.json")), outputs = outputs,
  summary_source = pin(sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1])),
  observed_counted_once = TRUE, no_valid_allocation_omitted = TRUE,
  two_sided_resolution_status = resolution_status,
  swap_maximum_scaled_difference = max(swap$maximum_scaled_difference),
  swap_above_tolerance_comparisons = sum(swap$above_fixed_tie_tolerance_n),
  gene_QLF_roundoff_negative_count = sum(vapply(manifests, function(x) x$negative_native_statistic_n, integer(1))),
  minimum_native_statistic = min(vapply(manifests, function(x) x$minimum_native_statistic, numeric(1))),
  global_conditional_warning_family_n = sum(summary[space == "global"]$conditional_tail_warning),
  candidate_gate_eligible = FALSE, authenticated_design_exchangeability = FALSE,
  other_sources_and_method_controls_complete = FALSE, primary_results_modified = FALSE),
  file.path(run, "ENUMERATION_COMPLETE.json"))
cat("ENUMERATION_COMPLETE scoped reference only; numerical acceptance:", numerical_acceptance,
  ";", resolution_status, "\n")
