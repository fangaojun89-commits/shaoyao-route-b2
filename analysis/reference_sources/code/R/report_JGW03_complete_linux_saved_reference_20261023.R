#!/usr/bin/env Rscript
# Saved-summary reporting only. No fits, RNG, tail completion, P/q adjustment,
# or archive scan. CLI: BINDING_JSON FRESH_OUTPUT_DIRECTORY.
jgwlrep_need <- function(ok, rule) {
  if (!isTRUE(ok)) stop(rule, call. = FALSE)
}
jgwlrep_cells <- function() c("GSE289929", "GSE227494", "GSE299063")
jgwlrep_columns <- function() c("tail_hits_nonidentity", "null_unavailable_n",
  "required_nonidentity_n", "total_denominator", "observed_available", "observed_self_hit",
  "conditional_P", "missingness_P_lower", "missingness_P_upper", "reference_complete_for_feature",
  "finite_population_P_upper", "finite_population_P_lower")
jgwlrep_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
jgwlrep_pin <- function(path) {
  jgwlrep_need(is.character(path) && length(path) == 1L && !is.na(path) &&
    file.exists(path) && !dir.exists(path), "ExistingExplicitFile")
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  jgwlrep_need(identical(gsub("\\\\", "/", path), resolved), "CanonicalExplicitFile")
  ancestors <- resolved
  while (!identical(dirname(tail(ancestors, 1L)), tail(ancestors, 1L))) {
    ancestors <- c(ancestors, dirname(tail(ancestors, 1L)))
  }
  links <- Sys.readlink(ancestors)
  jgwlrep_need(all(is.na(links) | !nzchar(links)), "NoSymlinkInput")
  list(path = resolved, bytes = unname(file.info(resolved)$size), sha256 = jgwlrep_sha(resolved))
}
jgwlrep_bind <- function(record, pins) {
  jgwlrep_need(is.list(record) && is.numeric(record$bytes) && length(record$bytes) == 1L &&
    is.finite(record$bytes) && record$bytes >= 0 && record$bytes == floor(record$bytes) &&
    is.character(record$sha256) && length(record$sha256) == 1L &&
    !is.na(record$sha256) && grepl("^[0-9a-f]{64}$", record$sha256), "ExplicitSizeHashRecord")
  actual <- jgwlrep_pin(record$path)
  jgwlrep_need(actual$bytes == record$bytes && identical(actual$sha256, record$sha256), "InputPinMatch")
  if (exists(actual$path, envir = pins, inherits = FALSE)) {
    jgwlrep_need(identical(get(actual$path, envir = pins), actual), "RepeatedPinMatch")
  }
  assign(actual$path, actual, envir = pins)
  actual
}
jgwlrep_manifest_match <- function(record, manifest) {
  jgwlrep_need(!grepl("^[Rr]:", record$path), "NoWindowsRScientificBinding")
  rows <- which(manifest$path == record$path)
  jgwlrep_need(length(rows) == 1L && manifest$bytes[rows] == record$bytes &&
    identical(manifest$sha256[rows], record$sha256), "ExactlyOneAcceptedManifestAssociation")
  invisible(TRUE)
}
jgwlrep_write_tsv <- function(table, connection) {
  # Seventeen significant digits preserve stored doubles through a text
  # round-trip; write.table's default numeric formatting is not that contract.
  numeric_columns <- vapply(table, is.numeric, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], function(x) sprintf("%.17g", x))
  utils::write.table(table, connection, sep = "\t", quote = TRUE, row.names = FALSE,
    na = "NA", qmethod = "double")
}

# Validate small saved objects, then return stored values verbatim. Arithmetic
# expressions below are assertions only; none replace a saved P or interval.
jgwlrep_cell <- function(domain, product, accession, space, evidence) {
  full_n <- if (accession == "GSE289929") 31L else 29L
  exact_n <- if (accession == "GSE289929") 12L else 11L
  total <- if (space == "four_group_MC") 10001L else 20L
  jgwlrep_need(accession %in% jgwlrep_cells() &&
    space %in% c("four_group_MC", "high_companion_exact"), "KnownCellAndSpace")
  body <- domain
  body$sha256 <- NULL
  jgwlrep_need(identical(domain$schema, "jgw03-native-compact-domain-v1") &&
    identical(domain$accession, accession) && length(domain$families) == full_n &&
    identical(digest::digest(body, "sha256", serializeVersion = 3L), domain$sha256) &&
    identical(domain$definition$native_P_q_changed, FALSE) &&
    identical(domain$definition$empirical_BH_computed, FALSE) &&
    identical(domain$definition$Mixed_outer_family_added, FALSE), "FullOwnLinuxDomain")
  ids <- domain$sample_ids
  labels <- domain$observed_labels
  jgwlrep_need(is.character(ids) && length(ids) == 12L && !anyNA(ids) && !anyDuplicated(ids) &&
    all(nzchar(ids)) && is.numeric(labels) && identical(names(labels), ids) && !anyNA(labels) &&
    all(labels %in% 0:3) && identical(tabulate(labels + 1L, 4L), rep(3L, 4L)), "TwelveLibraryGeometry")
  family_ids <- names(domain$families)
  jgwlrep_need(is.character(family_ids) && !anyNA(family_ids) && !anyDuplicated(family_ids), "UniqueFamilyIdentity")
  for (id in family_ids) {
    map <- domain$families[[id]]$map
    jgwlrep_need(identical(map$accession, accession) && identical(map$space, "four_group_MC") &&
      map$null_draw_n == 10000 && isTRUE(map$observed_once) &&
      identical(id, paste(map$method, map$collection, map$contrast, sep = "::")), "SourceMethodCollectionContrastIdentity")
  }
  selected <- family_ids
  if (space == "high_companion_exact") {
    selected <- selected[vapply(domain$families, function(x) identical(x$map$contrast, "high_vs_vehicle"), logical(1))]
  }
  state <- product$tail_checkpoint
  summary <- product$summary
  jgwlrep_need(length(selected) == if (space == "four_group_MC") full_n else exact_n,
    "CompleteSelectedFamilyCount")
  jgwlrep_need(identical(product$scientific_identity_sha256, evidence$scientific_identity_sha256) &&
    evidence$spaces[[space]]$allocation_n == total &&
    evidence$spaces[[space]]$selected_family_n == length(selected) &&
    identical(summary$schema, "jgw03-tail-summary-v1") && identical(summary$space, space) &&
    identical(summary$domain_sha256, domain$sha256) && summary$processed_nonidentity == total - 1L &&
    identical(names(summary$families), selected) && identical(names(state$families), selected) &&
    identical(state$schema, "jgw03-tail-accumulator-v1") && identical(state$space, space) &&
    identical(state$domain_sha256, domain$sha256) && state$total_N == total &&
    state$processed_nonidentity == total - 1L && identical(state$definition, domain$definition) &&
    identical(summary$native_P_q_changed, FALSE) && identical(summary$empirical_BH_computed, FALSE),
    "AcceptedSavedSummaryIdentityOrderAndDenominator")
  keys <- state$allocation_keys
  jgwlrep_need(is.character(keys) && length(keys) == total && !anyNA(keys) && !anyDuplicated(keys) &&
    all(grepl("^[0-3]{12}$", keys)) && identical(keys[1L], paste(as.integer(labels), collapse = "")), "FixedCompleteOrderedKeys")
  allocations <- lapply(strsplit(keys, "", fixed = TRUE), as.integer)
  nuisance <- which(!labels %in% c(0, 3))
  jgwlrep_need(all(vapply(allocations, function(x) identical(tabulate(x + 1L, 4L), rep(3L, 4L)) &&
    (space == "four_group_MC" || identical(x[nuisance], as.integer(labels[nuisance]))), logical(1))), "FrozenLabelGeometry")
  features <- families <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    id <- selected[i]
    spec <- domain$families[[id]]
    tab <- summary$families[[id]]
    saved <- state$families[[id]]
    n <- length(spec$keys)
    jgwlrep_need(is.character(spec$keys) && n > 0L && !anyNA(spec$keys) &&
      all(nzchar(spec$keys)) && !anyDuplicated(spec$keys) &&
      spec$map$native_key %in% c("raw_gene_id", "feature", "term_id") &&
      is.data.frame(tab) && nrow(tab) == n && length(names(tab)) == 12L &&
      identical(names(tab)[seq_len(10L)], jgwlrep_columns()[seq_len(10L)]) &&
      setequal(names(tab)[11:12], jgwlrep_columns()[11:12]), "FeatureKeysAndSavedTwelveColumns")
    # Older generators create the final two fields in either assignment order.
    # The export uses this declared canonical column order; values are untouched.
    tab <- tab[jgwlrep_columns()]
    hits <- saved$tail_hits_nonidentity
    missing <- saved$null_unavailable
    available <- saved$observed_available
    jgwlrep_need(is.integer(hits) && is.integer(missing) && is.logical(available) &&
      length(hits) == n && length(missing) == n && length(available) == n &&
      !anyNA(hits) && !anyNA(missing) && !anyNA(available) &&
      all(hits >= 0L & missing >= 0L & hits + missing <= total - 1L) &&
      all(hits[!available] == 0L) && is.numeric(saved$observed) && length(saved$observed) == n &&
      all(!is.na(saved$observed[available])), "SavedObservedCountersAndMissingness")
    complete <- available & missing == 0L
    jgwlrep_need(identical(tab$tail_hits_nonidentity, hits) && identical(tab$null_unavailable_n, missing) &&
      identical(tab$observed_available, available) && identical(tab$reference_complete_for_feature, complete) &&
      is.numeric(tab$required_nonidentity_n) && all(tab$required_nonidentity_n == total - 1L) &&
      is.numeric(tab$total_denominator) && all(tab$total_denominator == total) &&
      identical(tab$observed_self_hit, ifelse(available, 1L, NA_integer_)), "VerbatimSavedCounterFields")
    probability_columns <- c("conditional_P", "missingness_P_lower", "missingness_P_upper",
      "finite_population_P_upper", "finite_population_P_lower")
    jgwlrep_need(all(vapply(tab[probability_columns], function(x) is.numeric(x) &&
      !any(is.nan(x)) && all(is.na(x) | (is.finite(x) & x >= 0 & x <= 1)), logical(1))) &&
      identical(is.na(tab$conditional_P), !complete) &&
      identical(is.na(tab$missingness_P_lower), !available) &&
      identical(is.na(tab$missingness_P_upper), !available) &&
      identical(is.na(tab$finite_population_P_lower), !complete) &&
      identical(is.na(tab$finite_population_P_upper), !complete), "StoredProbabilityAndMissingnessStates")
    jgwlrep_need(all(tab$conditional_P[complete] == (1 + hits[complete]) / total) &&
      all(tab$missingness_P_lower[available] == (1 + hits[available]) / total) &&
      all(tab$missingness_P_upper[available] == (1 + hits[available] + missing[available]) / total) &&
      all(tab$finite_population_P_lower[complete] <= tab$finite_population_P_upper[complete]), "StoredArithmeticAssertionsOnly")
    if (space == "high_companion_exact") {
      jgwlrep_need(identical(tab$finite_population_P_lower, tab$conditional_P) &&
        identical(tab$finite_population_P_upper, tab$conditional_P), "StoredExactCensusEndpoints")
    }
    info <- data.frame(accession = accession, space = space, family_id = id,
      method = spec$map$method, collection = spec$map$collection, contrast = spec$map$contrast,
      stringsAsFactors = FALSE)
    features[[i]] <- cbind(info[rep(1L, n), , drop = FALSE], native_key_field = spec$map$native_key,
      key_order = seq_len(n), feature_key = spec$keys, observed_strength = saved$observed, tab, row.names = NULL)
    known <- tab$conditional_P[!is.na(tab$conditional_P)]
    range_or_na <- function(x, fun) if (all(is.na(x))) NA_real_ else fun(x, na.rm = TRUE)
    families[[i]] <- cbind(info, feature_n = n, observed_missing_n = sum(!available),
      null_missing_feature_n = sum(missing > 0L), conditional_P_available_n = length(known),
      conditional_P_unavailable_n = sum(is.na(tab$conditional_P)),
      null_unavailable_feature_draws = sum(missing), conditional_P_min = range_or_na(known, min),
      conditional_P_max = range_or_na(known, max), conditional_P_strict_lt_005_n = sum(known < .05),
      conditional_P_equal_005_n = sum(known == .05), conditional_P_le_005_n = sum(known <= .05),
      pointwise_finite_population_lower_min = range_or_na(tab$finite_population_P_lower, min),
      pointwise_finite_population_lower_max = range_or_na(tab$finite_population_P_lower, max),
      pointwise_finite_population_upper_min = range_or_na(tab$finite_population_P_upper, min),
      pointwise_finite_population_upper_max = range_or_na(tab$finite_population_P_upper, max),
      reference_N = total, nonidentity_N = total - 1L,
      uncertainty_scope = "Pointwise conditional finite-population reference; not effect-size, inner-MC, biological, simultaneous or FDR coverage")
  }
  list(families = do.call(rbind, families), features = do.call(rbind, features))
}

jgwlrep_main <- function(args) {
  jgwlrep_need(length(args) == 2L, "BindingAndFreshOutputRequired")
  pins <- new.env(parent = emptyenv())
  binding_pin <- jgwlrep_bind(jgwlrep_pin(args[1L]), pins)
  binding <- jsonlite::read_json(binding_pin$path, simplifyVector = FALSE)
  jgwlrep_need(identical(binding$schema, "jgw03-complete-linux-saved-report-binding-v1") &&
    identical(names(binding$cells), jgwlrep_cells()), "ExplicitOrderedThreeCellBinding")
  accepted_pin <- jgwlrep_bind(binding$accepted_readback, pins)
  manifest_pin <- jgwlrep_bind(binding$accepted_manifest, pins)
  exit_pin <- jgwlrep_bind(binding$actual_readback_parent_exit, pins)
  accepted <- jsonlite::read_json(accepted_pin$path, simplifyVector = FALSE)
  terminal <- jsonlite::read_json(exit_pin$path, simplifyVector = FALSE)
  jgwlrep_need(identical(accepted$status, "PASS_INDEPENDENT_COMPLETE_LINUX_SAVED_REFERENCE_ONLY") &&
    identical(accepted$source$sha256, "fd1bb54a1798cf51431dcc49044e3d6af039ed7b6b72549b6c7cbccf6c29e38f") &&
    identical(accepted$pin_manifest$path, manifest_pin$path) &&
    accepted$pin_manifest$bytes == manifest_pin$bytes &&
    identical(accepted$pin_manifest$sha256, manifest_pin$sha256) && accepted$total_MC_families == 89L &&
    accepted$total_companion_families == 34L && identical(names(accepted$cells), jgwlrep_cells()) &&
    accepted$new_native_fit_calls == 0L && isTRUE(accepted$all_bound_pins_unchanged) &&
    identical(accepted$unavailable_values_imputed, FALSE) && identical(accepted$denominator_reduced, FALSE) &&
    identical(accepted$native_BH_recomputed, FALSE) && identical(accepted$outer_empirical_BH_computed, FALSE),
    "ActualAcceptedV9ReadbackScope")
  jgwlrep_need(terminal$actual_native_exit_code == 0 && is.numeric(terminal$child_pid) &&
    length(terminal$child_pid) == 1L && terminal$child_pid > 0 && is.numeric(terminal$birth) &&
    length(terminal$birth) == 1L && terminal$birth > 0 &&
    identical(terminal$output, dirname(accepted_pin$path)) &&
    identical(terminal$raw_streams_persisted, FALSE) && identical(terminal$automatic_retry, FALSE),
    "GenuineSavedReaderParentExitNotComputeParent")
  connection <- gzfile(manifest_pin$path, "rt")
  manifest <- tryCatch(utils::read.delim(connection, quote = '"', check.names = FALSE,
    stringsAsFactors = FALSE), finally = close(connection))
  jgwlrep_need(identical(names(manifest), c("path", "bytes", "sha256")) &&
    nrow(manifest) == accepted$bound_files && !anyNA(manifest) && !anyDuplicated(manifest$path), "AcceptedManifestSchema")
  results <- list()
  for (accession in jgwlrep_cells()) {
    records <- binding$cells[[accession]]
    jgwlrep_need(setequal(names(records), c("domain", "four_group_MC", "high_companion_exact")), "ExactlyThreeCellFiles")
    for (record in records) jgwlrep_manifest_match(jgwlrep_bind(record, pins), manifest)
    domain <- readRDS(records$domain$path)
    for (space in c("four_group_MC", "high_companion_exact")) {
      results[[paste(accession, space)]] <- jgwlrep_cell(domain, readRDS(records[[space]]$path),
        accession, space, accepted$cells[[accession]])
    }
  }
  families <- do.call(rbind, lapply(results, `[[`, "families"))
  features <- do.call(rbind, lapply(results, `[[`, "features"))
  rownames(families) <- rownames(features) <- NULL
  jgwlrep_need(nrow(families) == 123L && sum(families$space == "four_group_MC") == 89L &&
    sum(families$space == "high_companion_exact") == 34L, "All123SeparateFamilies")
  self <- grep("^--file=", commandArgs(), value = TRUE)
  jgwlrep_need(length(self) == 1L, "ExplicitSourceEntry")
  self <- jgwlrep_bind(jgwlrep_pin(normalizePath(sub("^--file=", "", self), winslash = "/", mustWork = TRUE)), pins)
  before <- mget(ls(pins, all.names = TRUE), envir = pins)
  unchanged <- function() all(vapply(before, function(x) identical(jgwlrep_pin(x$path), x), logical(1)))
  out <- file.path(normalizePath(dirname(args[2L]), winslash = "/", mustWork = TRUE), basename(args[2L]))
  jgwlrep_need(!file.exists(out) && !basename(out) %in% c(".", "..") && unchanged(), "FreshOutputAndStableInputs")
  jgwlrep_need(dir.create(out, showWarnings = FALSE), "ExclusiveOutput")
  family_path <- file.path(out, "FAMILY_REFERENCE_SUMMARY.tsv")
  feature_path <- file.path(out, "FEATURE_SAVED_REFERENCE.tsv.gz")
  jgwlrep_write_tsv(families, family_path)
  connection <- gzfile(feature_path, "wt")
  tryCatch(jgwlrep_write_tsv(features, connection), finally = close(connection))
  report_path <- file.path(out, "READING_NOTES.md")
  writeLines(c("# Separate Linux JGW03 saved conditional references", "",
    "All 89 four-group MC and 34 complete high/vehicle companion families are exported without selecting features or families by P. Three original cellular contexts remain separate. Outputs describe the same 36 libraries, not new biological replicates.", "",
    "Stored conditional P, missingness bounds and pointwise finite-population endpoints are copied verbatim in a declared canonical column order (finite_population_P_upper then finite_population_P_lower). Four-group MC samples 10,000 distinct nonidentity assignments without replacement from the 369,599 nonidentity members of the balanced 369,600-assignment population, and adds observed once. Its conditional_P is a saved point label-reference estimate, not exhaustive four-group enumeration. The separate companion is the complete 20-assignment census (19 nonidentity plus observed once). Unavailable required strengths retain unavailable point P and fixed-denominator bounds. Companion P<0.05 is structurally unattainable; P=0.05 is a descriptive boundary, not a native discovery.", "",
    "No native P/q is present in these summary tables. P<, = and <=0.05 counts are unadjusted conditional-reference descriptors, not native or empirical-BH discoveries. Pointwise label-reference bounds are not effect-size confidence intervals or inner-MC, biological, simultaneous or FDR coverage. No models, RNG, P/q adjustment or endpoints were generated by this exporter. Arithmetic assertions do not replace saved values.", "",
    "Original Windows companions and native observed results remain separate and unchanged. S1 and candidate gates are unchanged. Full numerical readback is not independent biological validation or full native-algorithm certification."), report_path, useBytes = TRUE)
  jgwlrep_need(unchanged(), "InputsUnchangedAfterExport")
  receipt <- list(status = "PASS_SEPARATE_LINUX_SAVED_SUMMARY_DESCRIPTIVE_EXPORT_ONLY",
    UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), source = self, inputs = before,
    outputs = lapply(c(family_path, feature_path, report_path), jgwlrep_pin),
    family_n = 123L, MC_family_n = 89L, companion_family_n = 34L,
    feature_rows = nrow(features), all_input_pins_unchanged = TRUE, saved_twelve_column_values_copied = TRUE,
    canonical_saved_column_order = jgwlrep_columns(), full_balanced_four_group_population = 369600L,
    numeric_text_precision = "17 significant digits; no replacement of stored doubles",
    MC_sampling = "10000 distinct nonidentity assignments without replacement; observed added once",
    new_native_fit_calls = 0L, RNG_calls = 0L, P_q_BH_or_endpoints_recomputed = FALSE,
    arithmetic_assertions_only = TRUE, archives_scanned = 0L, unavailable_values_imputed = FALSE,
    denominator_reduced = FALSE, original_S1_or_candidate_gates_changed = FALSE,
    biological_or_whole_study_acceptance = FALSE,
    runtime = list(R = as.character(getRversion()), digest = as.character(utils::packageVersion("digest")),
      jsonlite = as.character(utils::packageVersion("jsonlite"))))
  jsonlite::write_json(receipt, file.path(out, "REPORT_COMPLETE.json"), pretty = TRUE,
    auto_unbox = TRUE, digits = NA, null = "null", na = "null")
  cat(receipt$status, "\n")
  invisible(0L)
}
if (sys.nframe() == 0L) {
  code <- tryCatch(jgwlrep_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_LINUX_REPORT_STOP_NO_RAW_CONDITION_PERSISTED")
    1L
  })
  quit(status = code)
}
