# Functions only: no files, model calls, RNG, package attachment or driver.
# Retrospective implementation, 15 October 2026; not a claim that these outer
# scores/availability rules were already specified in the original SAP.

jgw03c_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)

jgw03c_definition <- function() list(
  schema = "jgw03-native-compact-tail-v1",
  outer = c(edgeR_QL = "F unchanged", limma_voom_trend_TRUE = "abs(t)",
    DESeq2_native_IF_Cooks = "abs(stat)", CAMERA = "-log(PValue), derived P-strength",
    mean_z_limma = "abs(t)", ULM_activity_limma = "abs(t)",
    fgseaMultilevel = "abs(NES)",
    mroast1999_equivalent_not_historical_authentication = "-log(directional PValue), derived P-strength"),
  tail = "inclusive exact >=, no tolerance/rectification; signed infinities compared explicitly",
  availability = "native P finite in [0,1] and outer not NA/NaN; native q NA alone does not remove a statistic",
  MC = "10000 distinct nonidentity draws; observed separate; (1+b)/10001",
  companion = "19 nonidentity high/vehicle assignments plus observed once; total_hits/20",
  unavailable = "fixed denominator bounds; point P NA if observed or any required draw is unavailable",
  finite_mc_source_sha256 = "349a5de42e47c41c7211c9f64c0f012d3658c72d840f5b75e4deb8912d52db6b",
  native_P_q_changed = FALSE, Mixed_outer_family_added = FALSE,
  empirical_BH_computed = FALSE, real_reference_acceptance = FALSE)

jgw03c_table <- function(result, family) {
  label <- family$kernel_label
  method <- family$method
  collection <- family$collection
  value <- switch(method,
    edgeR_QL = result$gene_edgeR[[label]],
    limma_voom_trend_TRUE = result$gene_voom[[label]]$value,
    DESeq2_native_IF_Cooks = result$gene_DESeq2[[label]]$value,
    CAMERA = result$programs[[collection]][[label]],
    mean_z_limma = result$programs[[collection]]$score_fit[[label]],
    ULM_activity_limma = result$activities[[collection]][[label]],
    fgseaMultilevel = result$programs[[collection]][[paste0("fgsea_", label)]],
    mroast1999_equivalent_not_historical_authentication = result$programs$Reactome97$roast[[label]])
  jgw03c_need(is.data.frame(value) && nrow(value) > 0L, "Native result table absent or empty")
  # The saved QL test exposes actual df; never reconstruct an SE/CI from F.
  if (method == "edgeR_QL") {
    degrees <- result$native_edgeR_tests[[label]]$df.total
    if (length(degrees) == 1L) degrees <- rep(degrees, nrow(value))
    jgw03c_need(is.numeric(degrees) && length(degrees) == nrow(value), "Saved native QL df domain differs")
    if ("ql_df_total" %in% names(value)) jgw03c_need(
      identical(value$ql_df_total, degrees), "Saved QL df duplicated inconsistently") else
      value$ql_df_total <- degrees
  }
  value
}

jgw03c_labels <- function(result, sample_ids = NULL) {
  labels <- result$labels
  jgw03c_need(is.numeric(labels) && length(labels) == 12L &&
    !anyNA(labels) && all(is.finite(labels)) && all(labels %in% 0:3) &&
    !is.null(names(labels)) && !anyNA(names(labels)) &&
    all(nzchar(names(labels))) && !anyDuplicated(names(labels)) &&
    identical(tabulate(labels + 1L, 4L), rep(3L, 4L)), "Complete named four-by-three labels required")
  if (!is.null(sample_ids)) jgw03c_need(identical(names(labels), sample_ids), "Sample-label order differs")
  labels
}

# This one-time domain stores keys/order/schema, not matrices or fitted objects.
# The parent's admitted source/registry identity remains an external obligation.
jgw03c_make_domain <- function(observed_result, field_map) {
  accession <- observed_result$accession
  jgw03c_need(is.character(accession) && length(accession) == 1L &&
    accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    is.data.frame(field_map) && nrow(field_map) == 123L, "One admitted JGW03 result and full field map required")
  selected <- field_map[field_map$accession == accession & field_map$space == "four_group_MC", , drop = FALSE]
  jgw03c_need(nrow(selected) == if (accession == "GSE289929") 31L else 29L,
    "Complete accession-specific MC family domain required")
  ids <- paste(selected$method, selected$collection, selected$contrast, sep = "::")
  jgw03c_need(!anyDuplicated(ids), "Duplicate family domain")
  labels <- jgw03c_labels(observed_result)
  annotation_names <- c("raw_gene_id", "base_gene_id", "gene_symbol", "annotation_status", "human_symbol")
  excluded <- c(annotation_names, "ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high")
  families <- lapply(seq_len(nrow(selected)), function(i) {
    row <- as.list(selected[i, , drop = FALSE])
    table <- jgw03c_table(observed_result, row)
    key <- table[[row$native_key]]
    jgw03c_need(is.character(key) && !anyNA(key) && all(nzchar(key)) &&
      !anyDuplicated(key), "Complete unique native key vector required")
    indices <- which(!names(table) %in% c(excluded, row$native_key))
    numeric_indices <- indices[vapply(table[indices], is.numeric, logical(1))]
    dynamic_indices <- indices[names(table)[indices] %in% c("Direction", "leadingEdge")]
    static_indices <- setdiff(indices, c(numeric_indices, dynamic_indices))
    jgw03c_need(!anyDuplicated(names(table)[numeric_indices]) &&
      !anyDuplicated(names(table)[dynamic_indices]) &&
      all(c(row$native_P, row$native_q) %in% names(table)[numeric_indices]),
      "Native inference field schema differs")
    list(map = row, keys = key, names = names(table),
      numeric_indices = numeric_indices, numeric_names = names(table)[numeric_indices],
      numeric_types = vapply(table[numeric_indices], typeof, character(1)),
      dynamic_indices = dynamic_indices, dynamic_names = names(table)[dynamic_indices],
      static_indices = static_indices, static_values = table[static_indices])
  })
  names(families) <- ids
  domain <- list(schema = "jgw03-native-compact-domain-v1", accession = accession,
    sample_ids = names(labels), observed_labels = labels, families = families,
    definition = jgw03c_definition())
  domain$sha256 <- digest::digest(domain, algo = "sha256", serializeVersion = 3L)
  domain
}

jgw03c_numeric_counts <- function(value) c(n = length(value),
  NA_only = sum(is.na(value) & !is.nan(value)), "NaN" = sum(is.nan(value)),
  positive_Inf = sum(is.infinite(value) & value > 0),
  negative_Inf = sum(is.infinite(value) & value < 0), zero = sum(!is.na(value) & value == 0))

jgw03c_execution_metadata <- function(value) {
  if (is.null(value)) return(list(status = "NOT_PROVIDED",
    actual_seed_backend_authenticated_here = FALSE, inner_MC_uncertainty_covered = FALSE))
  required <- c("seed", "RNGkind", "backend", "workers", "mroast_nrot", "fgsea_eps", "provenance")
  jgw03c_need(is.list(value) && all(required %in% names(value)) &&
    is.numeric(value$seed) && identical(as.numeric(value$seed), 20260817) &&
    identical(value$RNGkind, c("Mersenne-Twister", "Inversion", "Rejection")) &&
    identical(value$backend, "SerialParam") && identical(as.numeric(value$workers), 1) &&
    identical(as.numeric(value$mroast_nrot), 1999) && identical(as.numeric(value$fgsea_eps), 0) &&
    is.character(value$provenance) && length(value$provenance) == 1L &&
    !is.na(value$provenance) && nzchar(value$provenance),
    "Caller-supplied actual execution metadata must retain the approved seed/backend/inner parameters")
  plain <- function(x) {
    if (isS4(x) || is.environment(x) || is.function(x) || is.language(x)) return(FALSE)
    if (is.list(x)) return(all(vapply(x, plain, logical(1))))
    is.null(x) || is.atomic(x)
  }
  jgw03c_need(plain(value), "Execution metadata must not contain fitted objects/environments")
  list(status = "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA", actual = value,
    actual_seed_backend_authenticated_here = FALSE,
    source_reset_each_inner_kernel = TRUE, inner_MC_uncertainty_covered = FALSE)
}

# Getter only, from the already fitted/saved DDS. Never call DESeq/results.
# The native function converted res to a data.frame, losing result metadata;
# acknowledge that rather than recreating independent-filtering choices.
jgw03c_saved_deseq_metadata <- function(result) {
  answer <- lapply(c("high", "trend"), function(label) {
    saved <- result$gene_DESeq2[[label]]$native_locals
    fit_type <- if (is.null(saved$dds)) list(status = "SAVED_DDS_NOT_RETAINED", value = NULL) else
      tryCatch({
        value <- attr(DESeq2::dispersionFunction(saved$dds), "fitType")
        list(status = if (is.null(value)) "FIT_TYPE_ATTRIBUTE_NOT_RETAINED" else "SAVED_FITTED_DDS_ATTRIBUTE",
          value = value)
      }, error = function(condition) list(status = "SAVED_ATTRIBUTE_GETTER_UNAVAILABLE",
        value = NULL, exception_class_only = class(condition)[1L]))
    list(fitType = fit_type, independentFiltering_metadata =
      if (is.data.frame(saved$res)) "NOT_RETAINED_NATIVE_RES_CONVERTED_TO_DATA_FRAME" else
        "NOT_CAPTURED_BY_THIS_COMPACT_HELPER_NO_RESULTS_RECOMPUTATION",
      new_fit_or_results_call = FALSE)
  })
  names(answer) <- c("high", "trend")
  answer
}

# Saved native tables only. All numeric columns and NA/NaN/Inf states retain
# their original vector values/types; dynamic text is retained, not fabricated.
# No keys, annotations, DDS, voom E or fitted object is repeated in this result.
jgw03c_extract <- function(result, domain, execution_metadata = NULL) {
  jgw03c_need(identical(domain$schema, "jgw03-native-compact-domain-v1") &&
    identical(result$accession, domain$accession) &&
    identical(domain$definition, jgw03c_definition()), "Accession/domain/definition mismatch")
  check_domain <- domain; check_domain$sha256 <- NULL
  jgw03c_need(identical(domain$sha256, digest::digest(check_domain, algo = "sha256", serializeVersion = 3L)),
    "Compact domain content differs")
  labels <- jgw03c_labels(result, domain$sample_ids)
  families <- lapply(domain$families, function(spec) {
    table <- jgw03c_table(result, spec$map)
    jgw03c_need(identical(names(table), spec$names) &&
      identical(table[[spec$map$native_key]], spec$keys) &&
      identical(vapply(table[spec$numeric_indices], typeof, character(1)), spec$numeric_types) &&
      identical(table[spec$static_indices], spec$static_values),
      "Native key/order/schema/static text differs; no reordering or row dropping allowed")
    native <- as.list(table[spec$numeric_indices])
    text <- as.list(table[spec$dynamic_indices])
    jgw03c_need(all(vapply(text, function(x) is.character(x) || is.factor(x) ||
      (is.logical(x) && all(is.na(x))), logical(1))),
      "Native dynamic text must retain source text/factor or all-missing logical encoding")
    raw_p <- native[[spec$map$native_P]]
    raw_q <- native[[spec$map$native_q]]
    method <- spec$map$method
    p_strength <- method %in% c("CAMERA", "mroast1999_equivalent_not_historical_authentication")
    if (p_strength) {
      outer <- raw_p
      valid <- is.finite(raw_p) & raw_p >= 0 & raw_p <= 1
      outer[valid] <- -log(raw_p[valid])  # Real P=0 stays Inf; never floored.
      outer[!is.na(raw_p) & !valid] <- NA_real_
    } else {
      outer <- native[[spec$map$native_statistic]]
      jgw03c_need(is.numeric(outer), "Required native statistic absent")
      if (method != "edgeR_QL") outer <- abs(outer)
    }
    available <- is.finite(raw_p) & raw_p >= 0 & raw_p <= 1 & !is.na(outer)
    list(native_numeric = native, native_text = text, outer_statistic = outer,
      outer_available = available,
      diagnostics = list(native_numeric_counts = lapply(native, jgw03c_numeric_counts),
        native_text_counts = lapply(text, function(x) c(n = length(x), unavailable = sum(is.na(x)),
          empty = sum(!is.na(x) & as.character(x) == ""))),
        outer_counts = jgw03c_numeric_counts(outer), outer_unavailable_n = sum(!available),
        raw_P_zero_n = sum(!is.na(raw_p) & raw_p == 0),
        raw_P_invalid_n = sum(!is.na(raw_p) & (!is.finite(raw_p) | raw_p < 0 | raw_p > 1)),
        raw_q_invalid_n = sum(!is.na(raw_q) & (!is.finite(raw_q) | raw_q < 0 | raw_q > 1)),
        native_q_unavailable_n = sum(is.na(raw_q)),
        native_q_le_005_known_n = sum(!is.na(raw_q) & is.finite(raw_q) & raw_q >= 0 & raw_q <= .05)))
  })
  list(schema = "jgw03-native-compact-v1", accession = domain$accession,
    domain_sha256 = domain$sha256, allocation_key = paste(as.integer(labels), collapse = ""),
    families = families, execution_metadata = jgw03c_execution_metadata(execution_metadata),
    saved_DESeq2_metadata = jgw03c_saved_deseq_metadata(result),
    native_table_numeric_and_dynamic_text_only = TRUE,
    DDS_or_fits_included = FALSE, definition = jgw03c_definition())
}

# Extended-real inclusive comparison, not abs(a-b) or a tolerance subtraction.
jgw03c_ge <- function(draw, observed) {
  jgw03c_need(is.numeric(draw) && is.numeric(observed) && length(draw) == length(observed),
    "Equal-length numeric strengths required")
  answer <- rep(NA, length(draw))
  known <- !is.na(draw) & !is.na(observed)
  finite <- known & is.finite(draw) & is.finite(observed)
  answer[finite] <- draw[finite] >= observed[finite]
  extended <- known & !finite
  answer[extended] <- draw[extended] >= observed[extended]  # +/-Inf same-sign ties are TRUE.
  answer
}

# allocation_keys are supplied from the separately pinned frozen registry,
# identity first. This helper authenticates geometry/order, not their RNG origin.
jgw03c_tail_start <- function(observed, domain, space, allocation_keys) {
  jgw03c_need(space %in% c("four_group_MC", "high_companion_exact") &&
    identical(observed$domain_sha256, domain$sha256) &&
    identical(observed$definition, domain$definition) &&
    identical(observed$allocation_key, paste(as.integer(domain$observed_labels), collapse = "")),
    "Tail observed identity/domain/space mismatch")
  total <- if (space == "four_group_MC") 10001L else 20L
  jgw03c_need(is.character(allocation_keys) && length(allocation_keys) == total &&
    !anyNA(allocation_keys) && all(grepl("^[0-3]{12}$", allocation_keys)) &&
    !anyDuplicated(allocation_keys) && allocation_keys[1L] == observed$allocation_key,
    "Complete distinct registry keys, observed first, required")
  labels <- lapply(strsplit(allocation_keys, "", fixed = TRUE), as.integer)
  jgw03c_need(all(vapply(labels, function(x) identical(tabulate(x + 1L, 4L), rep(3L, 4L)), logical(1))),
    "Every allocation must retain four groups of three")
  selected <- seq_along(domain$families)
  if (space == "high_companion_exact") {
    nuisance <- which(!domain$observed_labels %in% c(0, 3))
    original <- as.integer(domain$observed_labels)
    jgw03c_need(all(vapply(labels, function(x) identical(x[nuisance], original[nuisance]), logical(1))),
      "High companion may not change low/medium roles")
    selected <- which(vapply(domain$families, function(x) x$map$contrast == "high_vs_vehicle", logical(1)))
  }
  ids <- names(domain$families)[selected]
  jgw03c_need(identical(names(observed$families), names(domain$families)), "Observed family order differs")
  rows <- lapply(ids, function(id) {
    value <- observed$families[[id]]
    list(observed = value$outer_statistic, observed_available = value$outer_available,
      tail_hits_nonidentity = integer(length(value$outer_statistic)),
      null_unavailable = integer(length(value$outer_statistic)))
  })
  names(rows) <- ids
  list(schema = "jgw03-tail-accumulator-v1", domain_sha256 = domain$sha256,
    space = space, total_N = total, allocation_keys = allocation_keys,
    processed_nonidentity = 0L, families = rows,
    execution_metadata_missing_allocations = as.integer(
      !identical(observed$execution_metadata$status, "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA")),
    definition = jgw03c_definition())
}

jgw03c_tail_update <- function(state, draw) {
  next_id <- state$processed_nonidentity + 2L
  jgw03c_need(identical(state$schema, "jgw03-tail-accumulator-v1") &&
    identical(draw$domain_sha256, state$domain_sha256) && next_id <= state$total_N &&
    identical(draw$definition, state$definition) &&
    identical(draw$allocation_key, state$allocation_keys[next_id]) &&
    all(names(state$families) %in% names(draw$families)), "Next frozen nonidentity allocation required")
  for (id in names(state$families)) {
    saved <- state$families[[id]]
    value <- draw$families[[id]]
    jgw03c_need(is.logical(value$outer_available) && !anyNA(value$outer_available) &&
      length(value$outer_available) == length(saved$observed) &&
      is.logical(saved$observed_available) && !anyNA(saved$observed_available), "Tail availability domain differs")
    compare <- jgw03c_ge(value$outer_statistic, saved$observed)
    known <- value$outer_available & saved$observed_available
    jgw03c_need(!anyNA(compare[known]), "Available strengths must be comparable")
    # Unknowns are counted separately, never silently treated as a non-hit.
    saved$tail_hits_nonidentity <- saved$tail_hits_nonidentity + as.integer(known & !is.na(compare) & compare)
    saved$null_unavailable <- saved$null_unavailable + as.integer(!value$outer_available)
    state$families[[id]] <- saved
  }
  state$processed_nonidentity <- state$processed_nonidentity + 1L
  state$execution_metadata_missing_allocations <- state$execution_metadata_missing_allocations +
    as.integer(!identical(draw$execution_metadata$status, "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA"))
  state
}

# Vector arithmetic is exposed separately for small invented boundary checks.
# Completed draw count is mandatory; no denominator shrinks to complete cases.
jgw03c_tail_arithmetic <- function(hits_nonidentity, unavailable_nonidentity,
    observed_available, space, processed_nonidentity) {
  jgw03c_need(space %in% c("four_group_MC", "high_companion_exact"), "Declared tail space required")
  required <- if (space == "four_group_MC") 10000L else 19L
  integer_vector <- function(x) is.numeric(x) && length(x) > 0L && !anyNA(x) &&
    all(is.finite(x)) && all(x >= 0 & x == floor(x))
  jgw03c_need(integer_vector(hits_nonidentity) && integer_vector(unavailable_nonidentity) &&
    is.logical(observed_available) && !anyNA(observed_available) &&
    length(hits_nonidentity) == length(unavailable_nonidentity) &&
    length(hits_nonidentity) == length(observed_available) &&
    length(processed_nonidentity) == 1L && integer_vector(processed_nonidentity) &&
    processed_nonidentity == required &&
    all(hits_nonidentity + unavailable_nonidentity <= required),
    "Complete fixed draw denominator and valid hit/unavailable counts required")
  N <- required + 1
  lower <- (1 + hits_nonidentity) / N
  upper <- (1 + hits_nonidentity + unavailable_nonidentity) / N
  lower[!observed_available] <- upper[!observed_available] <- NA_real_
  point <- lower
  point[unavailable_nonidentity > 0 | !observed_available] <- NA_real_
  data.frame(tail_hits_nonidentity = hits_nonidentity,
    null_unavailable_n = unavailable_nonidentity, required_nonidentity_n = required,
    total_denominator = N, observed_available = observed_available,
    observed_self_hit = ifelse(observed_available, 1L, NA_integer_),
    conditional_P = point, missingness_P_lower = lower, missingness_P_upper = upper,
    reference_complete_for_feature = observed_available & unavailable_nonidentity == 0L,
    check.names = FALSE)
}

# finite_mc must be the parent's source-pinned accepted v2 function pair.
# The supplied SHA is recorded/asserted, not a substitute for parent source pins.
jgw03c_tail_finish <- function(state, finite_mc, finite_mc_source_sha256) {
  jgw03c_need(identical(finite_mc_source_sha256, jgw03c_definition()$finite_mc_source_sha256) &&
    is.list(finite_mc) && is.function(finite_mc$MC) && is.function(finite_mc$companion),
    "Accepted source-pinned finite-population v2 callback pair required")
  rows <- lapply(state$families, function(value) {
    answer <- jgw03c_tail_arithmetic(value$tail_hits_nonidentity, value$null_unavailable,
      value$observed_available, state$space, state$processed_nonidentity)
    answer$finite_population_P_lower <- answer$finite_population_P_upper <- NA_real_
    eligible <- which(answer$reference_complete_for_feature)
    if (length(eligible)) {
      if (state$space == "four_group_MC") {
        interval <- finite_mc$MC(value$tail_hits_nonidentity[eligible], draws = 10000L,
          nonidentity_population = 369599L, confidence = .95)$intervals
        jgw03c_need(identical(as.numeric(interval$SAP_MC_P), answer$conditional_P[eligible]),
          "Finite-population helper/SAP point convention differs")
      } else {
        interval <- finite_mc$companion(value$tail_hits_nonidentity[eligible] + 1L, confidence = .95)$intervals
        jgw03c_need(identical(as.numeric(interval$exact_companion_P), answer$conditional_P[eligible]),
          "Exact companion convention differs")
      }
      answer$finite_population_P_lower[eligible] <- interval$conditional_full_P_lower
      answer$finite_population_P_upper[eligible] <- interval$conditional_full_P_upper
    }
    answer
  })
  list(schema = "jgw03-tail-summary-v1", domain_sha256 = state$domain_sha256,
    space = state$space, processed_nonidentity = state$processed_nonidentity,
    families = rows, finite_mc_source_sha256 = finite_mc_source_sha256,
    execution_metadata_missing_allocations = state$execution_metadata_missing_allocations,
    execution_metadata_coverage_complete = state$execution_metadata_missing_allocations == 0L,
    execution_metadata_is_caller_supplied_not_independently_authenticated = TRUE,
    native_P_q_changed = FALSE, empirical_BH_computed = FALSE,
    uncertainty_scope = "Pointwise conditional finite-population MC only; no inner-MC/biological/FDR coverage",
    reference_acceptance = "Arithmetic only, not independent native or scientific acceptance")
}
