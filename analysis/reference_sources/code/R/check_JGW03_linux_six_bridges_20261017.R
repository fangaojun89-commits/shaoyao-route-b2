#!/usr/bin/env Rscript
# Saved-only independent readback. Sourcing defines functions, never fits.
# CLI: Rscript --vanilla THIS_FILE BINDING_JSON OUT
# Binding schema jgw03-linux-six-bridge-checker-binding-v1:
# config, handoff_source, controller_started, controller_exit: {path,bytes,sha256};
# controller_record_dir; controller_identity={pid,start_ticks}; jobs (nine, in
# cell/prepare/bridge1/bridge2 order): {accession,mode,allocation,attempt,
# controller_command,controller_exit,command,started,process_exit}, each of the
# last five a complete file record. Controller start_ticks is an externally
# observed binding, NOT a field claimed to exist in its STARTED.json.
# The producer comparison/extractor are never imported or invoked. Only the
# source-pinned table getter and pure metadata registry/map/definition are used.

jgw6_need <- function(ok, rule) {
  if (!isTRUE(ok)) stop(paste0("JGW6_RULE:", rule), call. = FALSE)
}
jgw6_hash <- function(value) digest::digest(value, "sha256", serializeVersion = 3L)
jgw6_pin <- function(path) {
  jgw6_need(is.character(path) && length(path) == 1L && !is.na(path) &&
    file.exists(path) && !dir.exists(path), "regular_file")
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
jgw6_record <- function(record) {
  jgw6_need(is.list(record) && is.numeric(record$bytes) && length(record$bytes) == 1L &&
    is.finite(record$bytes) && record$bytes >= 0 && record$bytes == floor(record$bytes) &&
    is.character(record$sha256) && length(record$sha256) == 1L &&
    !is.na(record$sha256) && grepl("^[0-9a-f]{64}$", record$sha256), "file_record")
  actual <- jgw6_pin(record$path)
  jgw6_need(actual$bytes == record$bytes && identical(actual$sha256, record$sha256), "file_identity")
  actual
}
jgw6_same_pin <- function(a, b) {
  is.list(a) && is.list(b) && isTRUE(a$bytes == b$bytes) && identical(a$sha256, b$sha256)
}
jgw6_json <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)
jgw6_json_same <- function(a, b) {
  # JSON numbers/arrays do not preserve R integer/double or named-list classes.
  identical(jsonlite::toJSON(a, auto_unbox = TRUE, null = "null", na = "null", digits = NA),
    jsonlite::toJSON(b, auto_unbox = TRUE, null = "null", na = "null", digits = NA))
}
jgw6_import <- function(path, functions, envir) {
  expressions <- parse(path, keep.source = FALSE)
  found <- character()
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))) {
      name <- as.character(expression[[2L]])
      if (name %in% functions) { eval(expression, envir); found <- c(found, name) }
    }
  }
  jgw6_need(setequal(found, functions) && !anyDuplicated(found), "pure_getter_definitions")
  invisible(TRUE)
}
jgw6_plain <- function(table) {
  jgw6_need(is.data.frame(table), "table_required")
  answer <- lapply(seq_along(table), function(i) table[[i]])
  names(answer) <- names(table)
  structure(answer, class = "data.frame", row.names = attr(table, "row.names"))
}
jgw6_probability <- function(name) {
  name %in% c("PValue", "FDR", "P.Value", "adj.P.Val", "pvalue", "padj", "pval",
    "PValue.Mixed", "FDR.Mixed", "FDR_high", "FDR_trend") ||
    grepl("(^|_)(PValue|FDR)$", name)
}
jgw6_direction <- function(name) {
  name %in% c("logFC", "effect", "t", "stat", "log2FoldChange", "ES", "NES", "F", "CI_low", "CI_high") ||
    grepl("_(effect|t|CI_low|CI_high)$", name)
}

# Reference is Windows; tolerance applies only to finite cross-platform values.
# Missingness, infinity sign, zero, direction and threshold masks are exact.
jgw6_vector <- function(actual, reference, probability = FALSE, direction = FALSE,
                        absolute = 1e-10, relative = 1e-8) {
  shape <- identical(typeof(actual), typeof(reference)) && identical(class(actual), class(reference)) &&
    length(actual) == length(reference) && identical(names(actual), names(reference))
  numeric <- is.numeric(actual) && is.numeric(reference)
  masks <- shape && identical(is.na(actual), is.na(reference)) &&
    identical(is.nan(actual), is.nan(reference))
  finite <- if (numeric && shape) is.finite(actual) & is.finite(reference) else logical()
  max_abs <- max_log <- NA_real_
  zeros <- directions <- thresholds <- log_ok <- range_ok <- TRUE
  values <- FALSE
  if (numeric && shape) {
    infinite <- !is.na(reference) & !is.finite(reference)
    inf_ok <- identical(is.infinite(actual), is.infinite(reference)) &&
      all(actual[infinite] == reference[infinite])
    difference <- abs(actual[finite] - reference[finite])
    if (length(difference)) max_abs <- max(difference)
    values <- masks && inf_ok && all(difference <= absolute + relative * abs(reference[finite]))
    zeros <- identical(!is.na(actual) & actual == 0, !is.na(reference) & reference == 0)
    if (direction) directions <- all(sign(actual[finite]) == sign(reference[finite]))
    if (probability) {
      range_ok <- all(is.na(actual) | (is.finite(actual) & actual >= 0 & actual <= 1)) &&
        all(is.na(reference) | (is.finite(reference) & reference >= 0 & reference <= 1))
      thresholds <- identical(actual <= .05, reference <= .05) &&
        identical(actual < .05, reference < .05)
      positive <- finite & actual > 0 & reference > 0
      if (any(positive)) max_log <- max(abs(log10(actual[positive]) - log10(reference[positive])))
      log_ok <- !any(positive) || max_log <= 1e-6
    }
  } else if (shape) values <- identical(actual, reference)
  list(pass = isTRUE(values && zeros && directions && thresholds && log_ok && range_ok),
    shape_types = shape, missing_masks = masks, zero_mask = zeros, direction_mask = directions,
    threshold_masks = thresholds, probability_range = range_ok, positive_log_bound = log_ok,
    n = length(reference), compared_finite = sum(finite), unavailable = sum(is.na(reference)),
    max_abs = max_abs, max_log10P = max_log)
}
jgw6_table <- function(actual, reference, endpoint) {
  a <- jgw6_plain(actual); b <- jgw6_plain(reference)
  excluded <- c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high")
  a <- structure(unclass(a)[!names(a) %in% excluded], class = "data.frame", row.names = attr(a, "row.names"))
  b <- structure(unclass(b)[!names(b) %in% excluded], class = "data.frame", row.names = attr(b, "row.names"))
  key <- if ("raw_gene_id" %in% names(b)) "raw_gene_id" else if ("term_id" %in% names(b)) "term_id" else "feature"
  schema <- identical(names(a), names(b)) && nrow(a) == nrow(b) && key %in% names(b) &&
    is.character(b[[key]]) && !anyNA(b[[key]]) && all(nzchar(b[[key]])) &&
    !anyDuplicated(b[[key]]) && identical(a[[key]], b[[key]])
  if (!schema) return(data.frame(endpoint = endpoint, column = "SCHEMA_KEYS_ORDER", pass = FALSE,
    shape_types = FALSE, missing_masks = NA, zero_mask = NA, direction_mask = NA,
    threshold_masks = NA, probability_range = NA, positive_log_bound = NA,
    n = nrow(b), compared_finite = 0L, unavailable = 0L, max_abs = NA_real_, max_log10P = NA_real_))
  rows <- lapply(seq_along(b), function(i) {
    name <- names(b)[i]
    result <- jgw6_vector(a[[i]], b[[i]], jgw6_probability(name), jgw6_direction(name))
    data.frame(endpoint = endpoint, column = paste0(name, "[", i, "]"),
      as.data.frame(result, stringsAsFactors = FALSE), check.names = FALSE)
  })
  do.call(rbind, rows)
}
jgw6_fixed <- function(actual, reference) {
  check_metadata <- function(a, b, dose = FALSE) {
    a <- jgw6_plain(a); b <- jgw6_plain(b)
    if (!identical(names(a), names(b)) || !identical(rownames(a), rownames(b))) return(FALSE)
    all(vapply(seq_along(b), function(i) {
      if (dose && names(b)[i] == "log10_dose_plus_1")
        jgw6_vector(a[[i]], b[[i]], absolute = 1e-12, relative = 1e-12)$pass else identical(a[[i]], b[[i]])
    }, logical(1)))
  }
  list(accession = identical(actual$accession, reference$accession),
    keep = identical(actual$keep, reference$keep),
    full_counts = identical(actual$full_counts, reference$full_counts),
    full_annotation = check_metadata(actual$full_annotation, reference$full_annotation),
    retained_counts = identical(actual$y$counts, reference$y$counts),
    retained_gene_sample_order = identical(dimnames(actual$y$counts), dimnames(reference$y$counts)),
    metadata = check_metadata(actual$original_samples, reference$original_samples, TRUE),
    annotation = check_metadata(actual$annotation, reference$annotation),
    library_sizes = isTRUE(all.equal(actual$y$samples$lib.size, reference$y$samples$lib.size,
      tolerance = 1e-12, check.attributes = FALSE)),
    TMM_factors = isTRUE(all.equal(actual$y$samples$norm.factors, reference$y$samples$norm.factors,
      tolerance = 1e-12, check.attributes = FALSE)))
}
jgw6_model <- function(actual, reference, labels, accession) {
  same_metadata <- function(a, b) {
    a <- jgw6_plain(a); b <- jgw6_plain(b)
    identical(names(a), names(b)) && identical(rownames(a), rownames(b)) &&
      all(vapply(seq_along(b), function(i) if (names(b)[i] == "log10_dose_plus_1")
        jgw6_vector(a[[i]], b[[i]], absolute = 1e-12, relative = 1e-12)$pass else
        identical(a[[i]], b[[i]]), logical(1)))
  }
  expected_dose <- (if (accession == "GSE289929") c(0, 4, 20, 100) else c(0, 20, 100, 500))[labels + 1L]
  samples <- actual$model_metadata
  expected_high <- outer(unname(labels), 0:3, "==") * 1
  groups <- c("vehicle", "low", "medium", "high")
  designs <- actual$designs
  row_order <- all(vapply(c("high", "trend", "logtrend"), function(name)
    identical(rownames(designs[[name]]), rownames(samples)), logical(1)))
  metadata_ok <- same_metadata(samples, reference$model_metadata) && nrow(samples) == 12L &&
    identical(as.character(samples$sample_id), names(labels)) &&
    identical(as.numeric(samples$ordered_dose), unname(labels)) &&
    identical(as.numeric(samples$dose_numeric), unname(expected_dose)) &&
    identical(as.character(samples$dose_group), groups[labels + 1L]) &&
    identical(levels(samples$dose_group), groups)
  high <- is.matrix(designs$high) && identical(dim(designs$high), c(12L, 4L)) &&
    identical(colnames(designs$high), groups) && all(designs$high == expected_high)
  trend <- is.matrix(designs$trend) && identical(dim(designs$trend), c(12L, 2L)) &&
    identical(colnames(designs$trend), c("(Intercept)", "ordered_dose")) &&
    all(designs$trend[, 1L] == 1) && all(designs$trend[, 2L] == labels)
  logtrend <- is.matrix(designs$logtrend) && identical(dim(designs$logtrend), c(12L, 2L)) &&
    identical(colnames(designs$logtrend), c("(Intercept)", "log10_dose_plus_1")) &&
    all(designs$logtrend[, 1L] == 1) &&
    all(abs(designs$logtrend[, 2L] - log10(expected_dose + 1)) <= 1e-12 * (1 + abs(log10(expected_dose + 1))))
  contrast <- is.matrix(designs$contrast_high) && identical(dim(designs$contrast_high), c(4L, 1L)) &&
    identical(rownames(designs$contrast_high), groups) && all(designs$contrast_high[, 1L] == c(-1, 0, 0, 1))
  list(metadata = metadata_ok, design_row_order = row_order, high = high, trend = trend,
    logtrend = logtrend, high_contrast = contrast)
}
jgw6_labels <- function(value, sample_ids) {
  is.numeric(value) && length(value) == 12L && !anyNA(value) && all(is.finite(value)) &&
    all(value %in% 0:3) && identical(names(value), sample_ids) && length(sample_ids) == 12L &&
    !anyNA(sample_ids) && all(nzchar(sample_ids)) && !anyDuplicated(sample_ids) &&
    identical(tabulate(value + 1L, 4L), rep(3L, 4L))
}
jgw6_registry <- function(registry, observed) {
  jgw6_need(jgw6_labels(observed, names(observed)), "observed_labels")
  mc <- registry$MC; day <- registry$high_companion
  check <- function(x, n) is.matrix(x) && is.numeric(x) && identical(dim(x), c(n, 12L)) &&
    identical(colnames(x), names(observed)) && !anyNA(x) && all(x %in% 0:3) &&
    all(rowSums(x == 0) == 3L & rowSums(x == 1) == 3L & rowSums(x == 2) == 3L & rowSums(x == 3) == 3L)
  jgw6_need(check(mc, 10001L) && check(day, 20L), "complete_registry_dimensions_groups")
  keys <- unname(apply(mc, 1L, paste0, collapse = ""))
  day_keys <- unname(apply(day, 1L, paste0, collapse = ""))
  key <- paste0(as.integer(observed), collapse = "")
  fixed <- which(observed %in% c(1, 2))
  jgw6_need(!anyDuplicated(keys) && !anyDuplicated(day_keys) && keys[1L] == key &&
    keys[2L] != key && day_keys[1L] == key && sum(day_keys == key) == 1L &&
    all(vapply(fixed, function(i) all(day[, i] == observed[i]), logical(1))) &&
    registry$identity_row == 1 && registry$MC_nonidentity_draw_n == 10000 && registry$full_space_n == 369600,
    "registry_unique_identity_companion_scope")
  list(keys = keys[1:2], MC_rows = nrow(mc), companion_rows = nrow(day),
    MC_unique = length(unique(keys)), companion_unique = length(unique(day_keys)))
}

# Independently select saved fields; never evaluate result_selector strings.
jgw6_native_table <- function(result, row) {
  label <- row$kernel_label; collection <- row$collection
  table <- switch(row$method,
    edgeR_QL = result$gene_edgeR[[label]],
    limma_voom_trend_TRUE = result$gene_voom[[label]]$value,
    DESeq2_native_IF_Cooks = result$gene_DESeq2[[label]]$value,
    CAMERA = result$programs[[collection]][[label]],
    mean_z_limma = result$programs[[collection]]$score_fit[[label]],
    ULM_activity_limma = result$activities[[collection]][[label]],
    fgseaMultilevel = result$programs[[collection]][[paste0("fgsea_", label)]],
    mroast1999_equivalent_not_historical_authentication = result$programs$Reactome97$roast[[label]])
  table <- jgw6_plain(table)
  jgw6_need(nrow(table) > 0L, "complete_native_family")
  if (identical(row$method, "edgeR_QL")) {
    degrees <- result$native_edgeR_tests[[label]]$df.total
    if (length(degrees) == 1L) degrees <- rep(degrees, nrow(table))
    jgw6_need(is.numeric(degrees) && length(degrees) == nrow(table), "saved_ql_df")
    if ("ql_df_total" %in% names(table)) jgw6_need(identical(table$ql_df_total, degrees), "ql_df_duplicate") else
      table$ql_df_total <- degrees
  }
  table
}
jgw6_counts <- function(value) c(n = length(value), NA_only = sum(is.na(value) & !is.nan(value)),
  "NaN" = sum(is.nan(value)), positive_Inf = sum(is.infinite(value) & value > 0),
  negative_Inf = sum(is.infinite(value) & value < 0), zero = sum(!is.na(value) & value == 0))
jgw6_outer <- function(table, row) {
  p <- table[[row$native_P]]; q <- table[[row$native_q]]
  jgw6_need(is.numeric(p) && is.numeric(q), "numeric_probability_vectors")
  valid <- is.finite(p) & p >= 0 & p <= 1
  if (row$method %in% c("CAMERA", "mroast1999_equivalent_not_historical_authentication")) {
    strength <- p
    strength[valid] <- -log(p[valid])
    strength[!is.na(p) & !valid] <- NA_real_
  } else {
    strength <- table[[row$native_statistic]]
    jgw6_need(is.numeric(strength), "native_statistic_vector")
    if (row$method != "edgeR_QL") strength <- abs(strength)
  }
  list(statistic = strength, available = valid & !is.na(strength))
}
jgw6_compact <- function(result, compact, domain, field_map, definition) {
  selected <- field_map[field_map$accession == result$accession & field_map$space == "four_group_MC", , drop = FALSE]
  n <- if (result$accession == "GSE289929") 31L else 29L
  ids <- paste(selected$method, selected$collection, selected$contrast, sep = "::")
  copy <- domain; copy$sha256 <- NULL
  jgw6_need(nrow(selected) == n && !anyDuplicated(ids) &&
    identical(domain$schema, "jgw03-native-compact-domain-v1") &&
    identical(domain$accession, result$accession) && identical(domain$definition, definition) &&
    identical(domain$sha256, jgw6_hash(copy)) && identical(names(domain$families), ids) &&
    identical(compact$schema, "jgw03-native-compact-v1") && identical(compact$accession, result$accession) &&
    identical(compact$domain_sha256, domain$sha256) && identical(compact$definition, definition) &&
    identical(names(compact$families), ids) && jgw6_labels(result$labels, domain$sample_ids) &&
    identical(compact$allocation_key, paste0(as.integer(result$labels), collapse = "")) &&
    isTRUE(compact$native_table_numeric_and_dynamic_text_only) && identical(compact$DDS_or_fits_included, FALSE),
    "linux_compact_identity_full_family_domain")
  excluded <- c("raw_gene_id", "base_gene_id", "gene_symbol", "annotation_status", "human_symbol",
    "ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high")
  answer <- lapply(seq_len(n), function(i) {
    row <- as.list(selected[i, , drop = FALSE]); table <- jgw6_native_table(result, row)
    spec <- domain$families[[i]]; family <- compact$families[[i]]
    keys <- table[[row$native_key]]
    jgw6_need(is.character(keys) && !anyNA(keys) && all(nzchar(keys)) && !anyDuplicated(keys), "native_family_keys")
    indices <- which(!names(table) %in% c(excluded, row$native_key))
    numeric_indices <- indices[vapply(indices, function(j) is.numeric(table[[j]]), logical(1))]
    text_indices <- indices[names(table)[indices] %in% c("Direction", "leadingEdge")]
    static_indices <- setdiff(indices, c(numeric_indices, text_indices))
    numeric <- as.list(table[numeric_indices]); text <- as.list(table[text_indices])
    outer <- jgw6_outer(table, row); p <- table[[row$native_P]]; q <- table[[row$native_q]]
    schema <- identical(spec$map, row) && identical(spec$keys, keys) && identical(spec$names, names(table)) &&
      identical(spec$numeric_indices, numeric_indices) && identical(spec$numeric_names, names(numeric)) &&
      identical(spec$numeric_types, vapply(numeric, typeof, character(1))) &&
      identical(spec$dynamic_indices, text_indices) && identical(spec$dynamic_names, names(text)) &&
      identical(spec$static_indices, static_indices) && identical(spec$static_values, table[static_indices])
    diagnostics <- list(native_numeric_counts = lapply(numeric, jgw6_counts),
      native_text_counts = lapply(text, function(x) c(n = length(x), unavailable = sum(is.na(x)),
        empty = sum(!is.na(x) & as.character(x) == ""))),
      outer_counts = jgw6_counts(outer$statistic), outer_unavailable_n = sum(!outer$available),
      raw_P_zero_n = sum(!is.na(p) & p == 0),
      raw_P_invalid_n = sum(!is.na(p) & (!is.finite(p) | p < 0 | p > 1)),
      raw_q_invalid_n = sum(!is.na(q) & (!is.finite(q) | q < 0 | q > 1)),
      native_q_unavailable_n = sum(is.na(q)),
      native_q_le_005_known_n = sum(!is.na(q) & is.finite(q) & q >= 0 & q <= .05))
    flags <- c(schema = schema, numeric_exact = identical(family$native_numeric, numeric),
      dynamic_text_exact = identical(family$native_text, text),
      outer_exact = identical(family$outer_statistic, outer$statistic),
      availability_exact = identical(family$outer_available, outer$available),
      diagnostics_exact = identical(family$diagnostics, diagnostics),
      native_probability_range = diagnostics$raw_P_invalid_n == 0L && diagnostics$raw_q_invalid_n == 0L)
    data.frame(family = ids[i], n = length(keys), as.list(flags), pass = all(flags), check.names = FALSE)
  })
  do.call(rbind, answer)
}
jgw6_deseq_get <- function(dds) {
  list(fitType = attr(DESeq2::dispersionFunction(dds), "fitType"),
    betaConv = S4Vectors::mcols(dds)$betaConv)
}
jgw6_deseq <- function(result, getter = jgw6_deseq_get) {
  answer <- lapply(c("high", "trend"), function(label) {
    saved <- result$gene_DESeq2[[label]]$native_locals
    jgw6_need(!is.null(saved$dds), "saved_DDS_required_no_reconstruction")
    fields <- getter(saved$dds); flags <- fields$betaConv
    jgw6_need(is.null(flags) || is.logical(flags), "native_betaConv_type")
    list(fitType = fields$fitType,
      betaConv = if (is.null(flags)) list(status = "NOT_RETAINED", false_n = NULL, unavailable_n = NULL) else
        list(status = "SAVED_NATIVE_FLAG_COUNTS", n = length(flags),
          false_n = sum(!is.na(flags) & !flags), unavailable_n = sum(is.na(flags))),
      native_results_metadata = "Original res converted to data.frame; no new results() reconstruction")
  })
  names(answer) <- c("high", "trend"); answer
}
jgw6_saved_deseq <- function(result, diagnostics) {
  answer <- lapply(c("high", "trend"), function(label) {
    value <- diagnostics[[label]]$fitType
    list(fitType = list(status = if (is.null(value)) "FIT_TYPE_ATTRIBUTE_NOT_RETAINED" else
      "SAVED_FITTED_DDS_ATTRIBUTE", value = value), independentFiltering_metadata =
      if (is.data.frame(result$gene_DESeq2[[label]]$native_locals$res))
        "NOT_RETAINED_NATIVE_RES_CONVERTED_TO_DATA_FRAME" else
        "NOT_CAPTURED_BY_THIS_COMPACT_HELPER_NO_RESULTS_RECOMPUTATION", new_fit_or_results_call = FALSE)
  })
  names(answer) <- c("high", "trend"); answer
}
jgw6_roast <- function(result) {
  if (result$accession != "GSE289929") return(list(applicable = FALSE, pass = TRUE))
  tables <- result$programs$Reactome97$roast
  jgw6_need(identical(names(tables), c("high", "trend")), "two_PC12_roast_tables")
  flags <- lapply(tables, function(table) {
    check <- function(p, q) is.numeric(p) && is.numeric(q) && length(p) == length(q) &&
      all(is.finite(p) & is.finite(q) & p >= 0 & p <= 1 & q >= 0 & q <= 1) &&
      all(abs(q - pmax(stats::p.adjust(p - .5 / 2000, method = "BH"), p)) <= 1e-12)
    c(directional = check(table$PValue, table$FDR), mixed = check(table$PValue.Mixed, table$FDR.Mixed))
  })
  list(applicable = TRUE, tables = flags, pass = all(unlist(flags)), nrot = 1999L,
    rule = "Saved native mid-P FDR arithmetic only; no new P, q, fit or reference")
}
jgw6_execution <- function(compact, saved, accession, prepared_sha) {
  x <- compact$execution_metadata; actual <- x$actual
  expected <- saved$execution; expected$RNGkind <- unlist(expected$RNGkind, use.names = FALSE)
  jgw6_need(identical(x$status, "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA") &&
    identical(x$actual_seed_backend_authenticated_here, FALSE) &&
    identical(x$inner_MC_uncertainty_covered, FALSE) && isTRUE(x$source_reset_each_inner_kernel) &&
    jgw6_json_same(actual, expected) && actual$seed == 20260817 && actual$workers == 1 &&
    identical(actual$backend, "SerialParam") && actual$mroast_nrot == 1999 && actual$fgsea_eps == 0 &&
    identical(actual$RNGkind, c("Mersenne-Twister", "Inversion", "Rejection")) &&
    identical(actual$mroast_executed, accession == "GSE289929") &&
    identical(actual$mroast_applicable, accession == "GSE289929") &&
    (if (accession == "GSE289929") isTRUE(actual$actual_mroast_nrot == 1999) else is.null(actual$actual_mroast_nrot)) &&
    identical(actual$provenance, paste("Separate Linux native bridge; prepared", prepared_sha)) &&
    is.character(actual$RNG_after_sha256) && grepl("^[0-9a-f]{64}$", actual$RNG_after_sha256) &&
    identical(actual$historical_backend_authenticated, FALSE), "execution_metadata_scope")
  invisible(TRUE)
}

# The pre-fit handoff may have a fresh, explicitly pinned successor. Do not
# mistake an obsolete controller PID/directory for a scientific input identity.
jgw6_controller_paths <- function(binding) {
  records <- binding$controller_record_dir
  jgw6_need(is.character(records) && length(records) == 1L && !is.na(records) &&
    identical(dirname(records), "/root/autodl-tmp") &&
    grepl("^JGW03_linux_handoff_20261017r[0-9]+$", basename(records)) &&
    identical(binding$controller_started$path, paste0(records, "/STARTED.json")) &&
    identical(binding$controller_exit$path, paste0(records, "/CONTROLLER_EXIT.json")), "explicit_controller_record_scope")
  records
}
jgw6_controller <- function(binding, began, ended, handoff) {
  records <- jgw6_controller_paths(binding)
  positive_int <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x > 0 && x == floor(x)
  identity <- binding$controller_identity
  jgw6_need(is.list(identity) && positive_int(identity$pid) && positive_int(identity$start_ticks) &&
    positive_int(began$actual_controller_pid) && began$actual_controller_pid == identity$pid &&
    identical(ended$status, "SIX_BRIDGES_EXITED_ZERO_PENDING_INDEPENDENT_ACCEPTANCE") &&
    ended$actual_controller_exit == 0 && is.null(ended$exception_class_only) &&
    identical(ended$raw_conditions_persisted, FALSE) && identical(ended$full_reference_dispatch, FALSE) &&
    identical(ended$independent_six_bridge_or_full_resource_acceptance_claimed, FALSE) &&
    length(ended$commands) == 9L && jgw6_same_pin(began$source, handoff) &&
    identical(began$full_reference_dispatch, FALSE), "real_controller_terminal")
  # The birth tick is retained as the main operator's external observation;
  # original controller STARTED.json does not claim to have recorded it.
  records
}

# This validates actual terminal records, not their acceptance-like status alone.
jgw6_process <- function(command, started, terminal, controller_command, controller_exit,
                         job, config, entry, supervisor, project) {
  argc <- c("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript", "--vanilla",
    file.path(project, "code/R/run_JGW03_full_reference_linux_20261017.R"), config$path,
    job$mode, job$accession, if (!is.null(job$allocation)) as.character(job$allocation))
  control <- unlist(controller_command$argv, use.names = FALSE)
  expected_tail <- c("-B", supervisor$path, "--config-sha256", config$sha256,
    "--mode", job$mode, "--accession", job$accession, "--attempt", job$attempt,
    if (!is.null(job$allocation)) c("--allocation", as.character(job$allocation)))
  jgw6_need(length(control) == length(expected_tail) + 1L && identical(control[-1L], expected_tail) &&
    identical(unlist(command$argv, use.names = FALSE), argc) &&
    jgw6_same_pin(command$config, config) && jgw6_same_pin(command$entry, entry) &&
    jgw6_same_pin(command$supervisor, supervisor) && jgw6_same_pin(controller_command$supervisor, supervisor) &&
    identical(controller_command$config_sha256, config$sha256) &&
    identical(command$mode, job$mode) && identical(command$accession, job$accession) &&
    identical(command$allocation, job$allocation) && identical(started$config_sha256, config$sha256), "process_command_scope")
  positive_int <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x > 0 && x == floor(x)
  jgw6_need(positive_int(started$actual_R_pid) && positive_int(started$start_ticks) &&
    terminal$actual_R_pid == started$actual_R_pid && terminal$actual_R_exit_code == 0 &&
    positive_int(terminal$peak_child_RSS_bytes) && terminal$peak_child_RSS_bytes <= 3 * 1024^3 &&
    isTRUE(terminal$bridge_memory_budget_met) && is.null(terminal$observer_failure_class_only) &&
    is.numeric(terminal$elapsed_seconds) && is.finite(terminal$elapsed_seconds) && terminal$elapsed_seconds >= 0 &&
    is.numeric(terminal$minimum_direct_memory_headroom_bytes) && is.finite(terminal$minimum_direct_memory_headroom_bytes) &&
    identical(terminal$raw_streams_persisted, FALSE) && isTRUE(terminal$no_automatic_retry_or_full_dispatch) &&
    identical(controller_exit$attempt, job$attempt) && positive_int(controller_exit$actual_supervisor_pid) &&
    controller_exit$actual_supervisor_exit == 0 && identical(controller_exit$raw_streams_persisted, FALSE) &&
    identical(controller_exit$automatic_retry, FALSE), "actual_process_exit_resource")
  for (item in list(terminal, controller_exit)) for (stream in c("stdout", "stderr")) {
    n <- item[[paste0(stream, "_bytes")]]; sha <- item[[paste0(stream, "_sha256")]]
    jgw6_need(is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0 && n == floor(n) &&
      is.character(sha) && length(sha) == 1L && grepl("^[0-9a-f]{64}$", sha), "safe_stream_metadata")
  }
  data.frame(accession = job$accession, mode = job$mode,
    allocation = if (is.null(job$allocation)) NA_integer_ else job$allocation,
    actual_R_pid = started$actual_R_pid, start_ticks = started$start_ticks,
    actual_R_exit_code = terminal$actual_R_exit_code, peak_RSS_bytes = terminal$peak_child_RSS_bytes,
    elapsed_seconds = terminal$elapsed_seconds, pass = TRUE)
}

jgw6_run <- function(binding_path, output) {
  jgw6_need(!file.exists(output) && dir.exists(dirname(output)), "fresh_output_existing_parent")
  tracked <- list()
  track <- function(rec) {
    pin <- jgw6_record(rec); tracked[[pin$path]] <<- pin; pin
  }
  read_meta <- function(rec) { pin <- track(rec); jgw6_json(pin$path) }
  local_meta <- function(path) read_meta(jgw6_pin(path))
  binding_pin <- track(jgw6_pin(binding_path)); binding <- jgw6_json(binding_pin$path)
  jgw6_need(identical(binding$schema, "jgw03-linux-six-bridge-checker-binding-v1"), "explicit_checker_binding")
  config_pin <- track(binding$config); cfg <- jgw6_json(config_pin$path)
  run <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r01"
  project <- paste0(run, "/project")
  records <- jgw6_controller_paths(binding)
  jgw6_need(identical(cfg$run_root, run) && identical(cfg$project_root, project) &&
    identical(config_pin$path, paste0(run, "/CONFIG.json")), "fixed_separate_linux_root")
  handoff <- track(binding$handoff_source)
  jgw6_need(identical(handoff$path, "/root/autodl-tmp/handoff_JGW03_linux_six_bridges_20261017_r01.py") &&
    identical(handoff$sha256, "245d9aa32ad8611c588caa116418c234b362fa866fc31a111490385ead02025d"), "handoff_source")
  began <- read_meta(binding$controller_started); ended <- read_meta(binding$controller_exit)
  jgw6_controller(binding, began, ended, handoff)
  source_names <- names(cfg$source_pins)
  jgw6_need(length(source_names) > 0L && !anyDuplicated(source_names) && all(nzchar(source_names)) &&
    all(!grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))", source_names)), "source_relative_paths")
  sources <- lapply(source_names, function(name) {
    rec <- cfg$source_pins[[name]]; rec$path <- file.path(project, name); track(rec)
  }); names(sources) <- source_names
  required <- c(
    "code/R/run_JGW03_full_reference_linux_20261017.R" = "57e3643a06b5c17f66c46eb658d150665c9bc55a1b4cc025ceb3e21d6d2cdc8f",
    "code/python/run_JGW03_linux_bridge_20261017.py" = "f388ecf2574d935865fb46d6cd0fa926284ed0fb8690ec2a043d5c0318d85e16",
    "code/R/prepare_JGW03_observed_label_pilot_v3_20261015.R" = "fe01d32bf79239cc51a5b4b75140c3ad209b31501adc94f141623c9ecf991f38",
    "code/R/run_JGW03_full_reference_single_20261015.R" = "f1de9c41003e0a96c6e08e3367a3e34f3a9ce898cce7a0be694bf9d4965ff4e5",
    "code/R/JGW03_label_reference_core_20261015.R" = "e3939f371f02b24cfdfefee0a9cd9e39bdf47c8920caa44f58434a93f7a2974e",
    "code/R/JGW03_reference_field_map_20261015.R" = "c349eab9c5f3718d95a06c0b51e44584f4adb60ecbfb2f2d75235af7bfd7da42",
    "code/R/JGW03_reference_compact_tail_20261015.R" = "1074cc5001899cf4f2500c72a52f6f218b1cbb763982c6a0a9b43bfeab82ce8d",
    "code/R/JGW03_input_admission_v2_20261015.R" = "63a24a37bac60e1e1c99cd11e9a69c6a11ab8dbb70c44766d06f2e33b87e9c0c",
    "code/R/JGW03_finite_population_mc_v2_20261015.R" = "349a5de42e47c41c7211c9f64c0f012d3658c72d840f5b75e4deb8912d52db6b",
    "code/R/analyze_jgw03_three_cell_lines.R" = "5df6378ed4dadfb4050d91e252a63af4f44d2be80d3070d889d1e0e007796775",
    "code/R/prepare_JGW03_reference_first_nonidentity_20261015.R" = "99e4c614c0ddcbe6fd1970f6fc1fb4029731b812e14f92a6aafce43d63e9b82c")
  jgw6_need(all(names(required) %in% names(sources)) && all(vapply(names(required), function(name)
    identical(sources[[name]]$sha256, required[[name]]), logical(1))), "reviewed_getter_map_sources")
  entry <- sources[[names(required)[1L]]]; supervisor <- sources[[names(required)[2L]]]
  inputs <- lapply(cfg$input_pins, track)
  jgw6_need(length(inputs) > 0L && all(vapply(inputs, function(x) startsWith(x$path, paste0(project, "/")), logical(1))), "staged_inputs")
  getter <- new.env(parent = baseenv())
  jgw6_import(sources[[names(required)[3L]]]$path, c("jgwp_need", "jgwp_fresh_tables", "jgwp_expected_tables"), getter)
  jgw6_import(sources[[names(required)[5L]]]$path, c("jgw03_require", "jgw03_family_registry"), getter)
  jgw6_import(sources[[names(required)[6L]]]$path, "jgw03_reference_field_map", getter)
  jgw6_import(sources[[names(required)[7L]]]$path, "jgw03c_definition", getter)
  registry_map <- getter$jgw03_family_registry(); field_map <- getter$jgw03_reference_field_map(registry_map)
  jgw6_need(nrow(field_map) == 123L && sum(field_map$space == "four_group_MC") == 89L &&
    sum(field_map$space == "high_companion_exact") == 34L, "full_89_34_metadata_scope")
  cells <- c(GSE289929 = "PC12", GSE227494 = "C2C12", GSE299063 = "HepG2")
  jgw6_need(setequal(names(cfg$cells), names(cells)) && length(binding$jobs) == 9L, "nine_jobs_three_cells")
  process_rows <- list(); job_number <- 0L
  for (accession in names(cells)) for (phase in c("prepare", "bridge1", "bridge2")) {
    job_number <- job_number + 1L; job <- binding$jobs[[job_number]]
    mode <- if (phase == "prepare") "prepare" else "bridge"
    allocation <- if (phase == "prepare") NULL else as.numeric(sub("bridge", "", phase))
    attempt <- paste0(cells[[accession]], "_", phase, "_01")
    folder <- paste0(run, "/provenance/", attempt)
    jgw6_need(identical(job$accession, accession) && identical(job$mode, mode) &&
      identical(job$allocation, allocation) && identical(job$attempt, attempt), "fixed_job_order_scope")
    paths <- c(controller_command = paste0(records, "/", attempt, "_COMMAND.json"),
      controller_exit = paste0(records, "/", attempt, "_EXIT.json"), command = paste0(folder, "/COMMAND.json"),
      started = paste0(folder, "/STARTED.json"), process_exit = paste0(folder, "/PROCESS_EXIT.json"))
    jgw6_need(all(vapply(names(paths), function(name) identical(job[[name]]$path, paths[[name]]), logical(1))), "job_record_paths")
    item <- lapply(job[names(paths)], read_meta)
    jgw6_need(jgw6_json_same(ended$commands[[job_number]], item$controller_exit), "controller_job_terminal_binding")
    process_rows[[job_number]] <- jgw6_process(item$command, item$started, item$process_exit,
      item$controller_command, item$controller_exit, job, config_pin, entry, supervisor, project)
    binary <- track(jgw6_pin(unlist(item$command$argv, use.names = FALSE)[1L]))
    jgw6_need(jgw6_same_pin(binary, item$command$Rscript), "actual_R_binary_pin")
    # Safe scalar resource provenance; acceptance here is finite RSS, not future reservation.
    local_meta(paste0(folder, "/ADMISSION.json"))
  }
  processes <- do.call(rbind, process_rows)
  jgw6_need(length(unique(processes$actual_R_pid)) == 9L, "nine_distinct_actual_R_pids")
  all_columns <- list(); all_families <- list(); summaries <- list(); sample_ids <- character()
  for (accession in names(cells)) {
    cell <- cfg$cells[[accession]]; leaf <- paste0(run, "/intermediate/", accession)
    references <- lapply(cell[c("registry", "windows_observed", "windows_nonidentity", "windows_fixed")], track)
    windows <- lapply(cell$windows_receipts, read_meta)
    jgw6_need(all(c("prepared", "observed_complete", "nonidentity_complete", "nonidentity_saved",
      "observed_process", "nonidentity_process") %in% names(windows)), "complete_windows_receipts")
    wp <- windows$prepared; wo <- windows$observed_complete; wn <- windows$nonidentity_complete; ws <- windows$nonidentity_saved
    jgw6_need(identical(wo$accession, accession) && identical(wn$accession, accession) &&
      identical(wo$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_ONLY") && isTRUE(wo$bridge$all_fields_pass) &&
      identical(wn$status, "FIRST_NONIDENTITY_SAVED_AND_COMPACTED_ONLY") && wn$allocation_id == 2 &&
      jgw6_same_pin(references$registry, wp$label_registry) && jgw6_same_pin(references$windows_observed, wo$result) &&
      jgw6_same_pin(references$windows_fixed, wp$observed_files$fixed) &&
      jgw6_same_pin(references$windows_nonidentity, wn$native) && jgw6_same_pin(references$windows_nonidentity, ws$native) &&
      identical(wn$prepared_sha256, wp$prepared$sha256) && identical(ws$prepared_sha256, wp$prepared$sha256) &&
      identical(wn$allocation_key, ws$allocation_key) && windows$observed_process$actual_R_exit_code == 0 &&
      windows$nonidentity_process$actual_R_exit_code == 0, "windows_reference_receipt_lineage")
    preparation_started <- local_meta(paste0(leaf, "/PREPARATION_STARTED.json"))
    preparation <- local_meta(paste0(leaf, "/PREPARED.json")); prepared_pin <- track(preparation$prepared)
    context <- readRDS(prepared_pin$path)
    jgw6_need(identical(context$schema, "jgw03-separate-linux-preparation-v1") &&
      identical(context$accession, accession) && identical(context$root, project) && identical(context$run_root, run) &&
      identical(preparation$config_sha256, config_pin$sha256) && identical(preparation$entry_sha256, entry$sha256) &&
      jgw6_same_pin(context$binding$config, config_pin) && jgw6_same_pin(context$entry, entry) &&
      jgw6_json_same(preparation_started$binding, context$binding) &&
      jgw6_json_same(context$binding$sources, sources) && jgw6_json_same(context$binding$inputs, inputs) &&
      jgw6_json_same(context$binding$references, references) &&
      jgw6_json_same(context$binding$windows_receipts, lapply(cell$windows_receipts, jgw6_record)) &&
      jgw6_same_pin(context$binding$entry, entry) &&
      context$full_MC_nonidentity_draws == 10000 && context$full_companion_including_observed == 20 &&
      context$full_MC_family_n == 89 && context$full_companion_family_n == 34 && context$native_allocation_calls == 0 &&
      isTRUE(context$independent_bridge_and_resource_acceptance_required), "linux_preparation_identity")
    runtime_fields <- c("R", "packages", "platform", "BLAS", "backend", "workers", "seed", "RNGkind")
    jgw6_need(jgw6_json_same(context$runtime[runtime_fields], cfg$expected_runtime[runtime_fields]) &&
      identical(context$runtime$R, "4.6.1") && grepl("linux", context$runtime$platform) &&
      all(unlist(context$runtime$thread_environment, use.names = FALSE) == "1"), "recorded_linux_runtime")
    for (name in c("fixed", "registry", "session_info")) {
      track(context[[name]]); jgw6_need(jgw6_same_pin(context[[name]], preparation[[name]]), "prepared_product_lineage")
    }
    track(context$fixed_bridge)
    for (rec in context$sources) track(rec)
    registry <- readRDS(context$registry$path)
    jgw6_need(jgw6_same_pin(context$registry, references$registry), "original_registry_bytes_unchanged")
    registry_check <- jgw6_registry(registry, context$observed_labels)
    jgw6_need(identical(registry_check$keys[2L], wn$allocation_key), "fixed_id2_key")
    fixed <- readRDS(context$fixed$path); old_fixed <- readRDS(references$windows_fixed$path)
    fixed_checks <- jgw6_fixed(fixed, old_fixed)
    jgw6_need(jgw6_labels(fixed$observed_labels, names(context$observed_labels)) &&
      identical(fixed$observed_labels, context$observed_labels) &&
      identical(colnames(fixed$y$counts), names(context$observed_labels)) &&
      identical(as.character(fixed$original_samples$sample_id), names(context$observed_labels)), "fixed_sample_alignment")
    sample_ids <- c(sample_ids, names(context$observed_labels))
    domain <- NULL; cell_summary <- list(fixed = fixed_checks, registry = registry_check, bridges = list())
    for (id in 1:2) {
      bridge_leaf <- paste0(leaf, "/bridge_", id)
      start <- local_meta(paste0(bridge_leaf, "/FIT_STARTED.json"))
      saved <- local_meta(paste0(bridge_leaf, "/FIT_SAVED.json"))
      bridge <- local_meta(paste0(bridge_leaf, "/BRIDGE_COMPLETE.json"))
      jgw6_need(identical(bridge$status, "LINUX_NATIVE_BRIDGE_PASS_ONLY") && bridge$allocation_id == id &&
        identical(bridge$accession, accession) && identical(bridge$allocation_key, registry_check$keys[id]) &&
        identical(bridge$prepared_sha256, prepared_pin$sha256) && identical(bridge$config_sha256, config_pin$sha256) &&
        identical(bridge$entry_sha256, entry$sha256) && identical(start$prepared_sha256, prepared_pin$sha256) &&
        identical(start$allocation_key, registry_check$keys[id]) && jgw6_same_pin(start$entry, entry) &&
        identical(saved$allocation_key, registry_check$keys[id]) && identical(saved$prepared_sha256, prepared_pin$sha256) &&
        identical(saved$config_sha256, config_pin$sha256) && identical(saved$entry_sha256, entry$sha256) &&
        jgw6_same_pin(saved$native, bridge$native) && jgw6_same_pin(saved$session_info, bridge$session_info) &&
        isTRUE(bridge$all_fields_pass) && identical(bridge$full_execution_admitted_here, FALSE) &&
        identical(bridge$original_Windows_result_or_gate_changed, FALSE) && isTRUE(bridge$tail_uses_no_tolerance),
        "bridge_native_prepared_config_lineage")
      for (name in c("native", "compact", "checks", "checks_saved", "session_info")) track(bridge[[name]])
      check_saved <- jgw6_json(bridge$checks_saved$path)
      jgw6_need(jgw6_same_pin(check_saved$file, bridge$checks) && identical(check_saved$native_sha256, bridge$native$sha256) &&
        identical(check_saved$prepared_sha256, prepared_pin$sha256), "producer_diagnostic_byte_lineage_not_numeric_oracle")
      if (id == 1L) { track(bridge$domain); domain <- readRDS(bridge$domain$path) }
      reference <- if (id == 1L) references$windows_observed else references$windows_nonidentity
      jgw6_need(jgw6_same_pin(bridge$Windows_reference_only_never_used_as_Linux_result, reference) &&
        !identical(bridge$native$path, reference$path), "linux_native_not_windows_observed")
      native <- readRDS(bridge$native$path); old <- readRDS(reference$path); compact <- readRDS(bridge$compact$path)
      labels <- setNames(as.numeric(registry$MC[id, ]), colnames(registry$MC))
      jgw6_need(identical(native$accession, accession) && identical(old$accession, accession) &&
        identical(native$labels, labels) && identical(old$labels, labels) &&
        identical(domain$observed_labels, context$observed_labels), "native_actual_allocation_identity")
      fresh_tables <- getter$jgwp_fresh_tables(fixed, native); old_tables <- getter$jgwp_fresh_tables(old_fixed, old)
      jgw6_need(length(old_tables) > 0L && setequal(names(fresh_tables), names(old_tables)) &&
        !anyDuplicated(names(fresh_tables)) && !anyDuplicated(names(old_tables)) &&
        setequal(names(old_tables), basename(getter$jgwp_expected_tables(accession))), "all_native_and_score_endpoint_coverage")
      column_rows <- lapply(names(old_tables), function(name) jgw6_table(fresh_tables[[name]], old_tables[[name]], name))
      # Schema-failure rows have fewer diagnostic fields; explicit union fill.
      cols <- unique(unlist(lapply(column_rows, names), use.names = FALSE))
      column_rows <- lapply(column_rows, function(x) { for (name in setdiff(cols, names(x))) x[[name]] <- NA; x[, cols, drop = FALSE] })
      columns <- do.call(rbind, column_rows); columns$accession <- accession; columns$allocation_id <- id
      families <- jgw6_compact(native, compact, domain, field_map, getter$jgw03c_definition())
      families$accession <- accession; families$allocation_id <- id
      native_deseq <- jgw6_deseq(native); windows_deseq <- jgw6_deseq(old)
      model_checks <- jgw6_model(native, old, labels, accession)
      roast_checks <- jgw6_roast(native)
      desequ <- jgw6_json_same(native_deseq, windows_deseq) &&
        jgw6_json_same(native_deseq, bridge$DESeq2_native_diagnostics) &&
        jgw6_json_same(windows_deseq, bridge$Windows_DESeq2_diagnostics) &&
        identical(compact$saved_DESeq2_metadata, jgw6_saved_deseq(native, native_deseq))
      jgw6_execution(compact, saved, accession, prepared_pin$sha256)
      all_columns[[paste(accession, id)]] <- columns; all_families[[paste(accession, id)]] <- families
      cell_summary$bridges[[as.character(id)]] <- list(all_tables_pass = all(columns$pass),
        endpoint_n = length(old_tables), columns_checked = nrow(columns), native_families = nrow(families),
        all_native_compact_vectors_exact = all(families$pass), DESeq2_diagnostics_exact = desequ,
        DESeq2_native_diagnostics = native_deseq, saved_model_metadata_design = model_checks,
        saved_native_mroast_FDR_arithmetic = roast_checks,
        original_approximate_QL_interval_fields_excluded_only =
          c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high"))
      rm(native, old, compact, fresh_tables, old_tables); invisible(gc())
    }
    summaries[[accession]] <- cell_summary
    rm(fixed, old_fixed, registry, context, domain); invisible(gc())
  }
  columns <- do.call(rbind, all_columns); families <- do.call(rbind, all_families)
  unchanged <- all(vapply(tracked, function(rec) identical(jgw6_record(rec), rec), logical(1)))
  pass <- unchanged && length(sample_ids) == 36L && !anyDuplicated(sample_ids) && all(columns$pass) &&
    all(families$pass) && nrow(families) == 178L && all(vapply(summaries, function(x)
      all(unlist(x$fixed)) && all(vapply(x$bridges, function(b) isTRUE(b$DESeq2_diagnostics_exact) &&
        all(unlist(b$saved_model_metadata_design)) && isTRUE(b$saved_native_mroast_FDR_arithmetic$pass), logical(1))), logical(1)))
  jgw6_need(dir.create(output, recursive = FALSE), "exclusive_new_output")
  utils::write.table(columns, file.path(output, "COLUMN_CHECKS.tsv"), sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  utils::write.table(families, file.path(output, "FAMILY_CHECKS.tsv"), sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  utils::write.table(processes, file.path(output, "PROCESS_CHECKS.tsv"), sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  result <- list(schema = "jgw03-independent-linux-six-bridge-readback-v1",
    status = if (pass) "SIX_NATIVE_VECTOR_BRIDGES_ONLY_NOT_FULL_REFERENCE_OR_RESOURCE_ACCEPTANCE" else "FAIL_SAVED_BRIDGE_READBACK",
    validation_status = if (pass) "PASS" else "FAIL", all_checks_pass = pass,
    binding = binding_pin, checker_source = {
      arg <- grep("^--file=", commandArgs(), value = TRUE)
      if (length(arg) == 1L) jgw6_pin(sub("^--file=", "", arg)) else NULL
    }, controller_identity_external_start_ticks_not_historical_STARTED_field = binding$controller_identity,
    processes = processes, cells = summaries, MC_family_n = 89L, exact_family_n = 34L,
    six_saved_native_objects = 6L, two_allocation_family_rows = nrow(families),
    sources_inputs_native_and_metadata_unchanged = unchanged,
    before_after_file_pins = unname(tracked),
    before_validation = list(libraries = 36L, unique_library_ids = length(unique(sample_ids)), unique_subjects = "NOT_AUTHENTICATED"),
    after_validation = list(original_libraries_retained = 36L, removed_libraries = 0L, input_modified = FALSE),
    clinical_QC = "NOT_APPLICABLE_TO_SAVED_LIBRARY_BRIDGE_SCOPE",
    warnings = c("Library observations are not authenticated independent subjects or randomized culture units",
      "Recorded batch does not authenticate exchangeability or paired cultures",
      "Same seed/backend and saved final RNG hash do not prove all internal RNG consumption",
      "Inner mroast/fgsea Monte Carlo uncertainty remains outside this saved-vector bridge",
      "Windows references are comparison objects, never Linux observed/null results"),
    full_reference_accepted = FALSE, full_dispatch_or_resource_acceptance = FALSE,
    fits_or_results_recomputed = FALSE, raw_conditions_persisted = FALSE,
    output_files = lapply(c("COLUMN_CHECKS.tsv", "FAMILY_CHECKS.tsv", "PROCESS_CHECKS.tsv"),
      function(name) jgw6_pin(file.path(output, name))))
  jsonlite::write_json(result, file.path(output, "VALIDATION.json"), pretty = TRUE,
    auto_unbox = TRUE, digits = NA, null = "null", na = "null")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  code <- tryCatch({
    jgw6_need(length(args) == 2L, "CLI_BINDING_JSON_OUT")
    answer <- jgw6_run(args[1L], args[2L])
    cat(jsonlite::toJSON(list(status = answer$status, all_checks_pass = answer$all_checks_pass,
      raw_conditions_persisted = FALSE), auto_unbox = TRUE), "\n")
    if (isTRUE(answer$all_checks_pass)) 0L else 2L
  }, error = function(condition) {
    # Never print/persist raw native condition text or a call trace.
    cat(jsonlite::toJSON(list(status = "CHECKER_STOP_RETAIN_INPUTS", exception_class_only = class(condition)[1L],
      raw_conditions_persisted = FALSE), auto_unbox = TRUE), "\n"); 2L
  })
  quit(save = "no", status = code)
}
