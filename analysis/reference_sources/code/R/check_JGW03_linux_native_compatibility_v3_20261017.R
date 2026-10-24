#!/usr/bin/env Rscript
# Saved-only, one native object. CLI: CONFIG ACCESSION ID FRESH_OUT.
# This is a new, explicitly non-equivalence contract. Original abs/rel/log-P
# bounds are diagnostics, never enlarged or re-labelled as their old PASS.
# No fitting, results(), stochastic call, installation or dispatch occurs.

jgw2_need <- function(ok, rule) {
  if (!isTRUE(ok)) stop(paste0("JGW2_RULE:", rule), call. = FALSE)
}
jgw2_pin <- function(path) {
  jgw2_need(is.character(path) && length(path) == 1L && !is.na(path) &&
    file.exists(path) && !dir.exists(path), "regular_file")
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE), bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
jgw2_record <- function(record) {
  jgw2_need(is.list(record) && is.numeric(record$bytes) && length(record$bytes) == 1L &&
    is.finite(record$bytes) && record$bytes >= 0 && record$bytes == floor(record$bytes) &&
    is.character(record$sha256) && length(record$sha256) == 1L &&
    !is.na(record$sha256) && grepl("^[0-9a-f]{64}$", record$sha256), "file_record")
  actual <- jgw2_pin(record$path)
  jgw2_need(actual$bytes == record$bytes && identical(actual$sha256, record$sha256), "file_identity")
  actual
}
jgw2_same_pin <- function(a, b) is.list(a) && is.list(b) &&
  isTRUE(a$bytes == b$bytes) && identical(a$sha256, b$sha256)
jgw2_json <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)
jgw2_plain <- function(table) {
  jgw2_need(is.data.frame(table), "table")
  value <- lapply(seq_along(table), function(i) table[[i]])
  names(value) <- names(table)
  structure(value, class = "data.frame", row.names = attr(table, "row.names"))
}
jgw2_import <- function(path, wanted, envir) {
  found <- character()
  for (expression in parse(path, keep.source = FALSE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))) {
      name <- as.character(expression[[2L]])
      if (name %in% wanted) { eval(expression, envir); found <- c(found, name) }
    }
  }
  jgw2_need(setequal(found, wanted) && !anyDuplicated(found), "pure_definitions")
  invisible(TRUE)
}
jgw2_probability <- function(name) name %in%
  c("pvalue", "padj", "pval", "P.Value", "adj.P.Val", "PValue.Mixed", "FDR.Mixed") ||
  grepl("(^|_)(PValue|FDR)($|_)", name)
jgw2_direction <- function(name) name %in%
  c("logFC", "effect", "t", "stat", "log2FoldChange", "ES", "NES", "F", "CI_low", "CI_high") ||
  grepl("_(effect|t|CI_low|CI_high)$", name)
jgw2_utf8 <- function(value) enc2utf8(as.character(value))

# Predeclared arithmetic tolerances are scaled machine precision, not fitted to
# any observed Windows/Linux discrepancy. Cross-platform old bounds are separate.
jgw2_arithmetic_policy <- function() list(absolute = 64 * .Machine$double.eps,
  relative = 512 * .Machine$double.eps, positive_probability_log = 4096 * .Machine$double.eps,
  rationale = "Forward arithmetic from already saved native sufficient fields; no numerical refitting",
  old_crossplatform_absolute = 1e-10, old_crossplatform_relative = 1e-8,
  old_crossplatform_log10P = 1e-6)
jgw2_math <- function(actual, expected, probability = FALSE) {
  policy <- jgw2_arithmetic_policy()
  shape <- is.numeric(actual) && is.numeric(expected) && length(actual) == length(expected)
  masks <- shape && identical(is.na(actual), is.na(expected)) && identical(is.nan(actual), is.nan(expected)) &&
    identical(is.infinite(actual), is.infinite(expected))
  finite <- if (shape) is.finite(actual) & is.finite(expected) else logical()
  unavailable <- if (shape) !is.na(expected) & !is.finite(expected) else logical()
  delta <- if (any(finite)) abs(actual[finite] - expected[finite]) else numeric()
  values <- masks && all(actual[unavailable] == expected[unavailable]) &&
    all(delta <= policy$absolute + policy$relative * abs(expected[finite]))
  log_delta <- numeric()
  if (probability && shape) {
    known <- !is.na(actual)
    values <- values && all(is.finite(actual[known]) & actual[known] >= 0 & actual[known] <= 1) &&
      identical(!is.na(actual) & actual == 0, !is.na(expected) & expected == 0)
    positive <- finite & actual > 0 & expected > 0
    log_delta <- abs(log(actual[positive]) - log(expected[positive]))
    values <- values && all(log_delta <= policy$positive_probability_log * (1 + abs(log(expected[positive]))))
  }
  list(pass = isTRUE(values), n = length(expected), missing_masks_exact = masks,
    max_abs = if (length(delta)) max(delta) else NA_real_,
    max_natural_logP = if (length(log_delta)) max(log_delta) else NA_real_)
}

# Base-R BH formula, preserving NA/NaN positions. The default native call's
# adjustment universe is its nonmissing P vector, not another method/family.
jgw2_bh <- function(p) {
  jgw2_need(is.numeric(p) && all(is.na(p) | (is.finite(p) & p >= 0 & p <= 1)), "BH_probability_domain")
  answer <- as.numeric(p); names(answer) <- names(p)
  available <- which(!is.na(p)); n <- length(available)
  if (!n) return(answer)
  order <- order(p[available], seq_along(available))
  sorted <- p[available][order]
  adjusted <- pmin(1, rev(cummin(rev((n / seq_len(n)) * sorted))))
  answer[available[order]] <- adjusted
  answer
}
jgw2_midp <- function(p) {
  jgw2_need(is.numeric(p) && all(is.finite(p) & p >= .5 / 2000 & p <= 1), "native_mroast_2000_midP_domain")
  pmax(jgw2_bh(p - .5 / 2000), p)
}
jgw2_compare_column <- function(actual, reference, name) {
  types <- identical(typeof(actual), typeof(reference)) && identical(class(actual), class(reference)) &&
    length(actual) == length(reference)
  masks <- types && identical(is.na(actual), is.na(reference)) && identical(is.nan(actual), is.nan(reference))
  numeric <- is.numeric(actual) && is.numeric(reference)
  infinite <- zero <- direction <- threshold <- range <- TRUE
  text_equal <- TRUE; above <- above_log <- 0L; max_abs <- max_log <- NA_real_
  if (types && numeric) {
    finite <- is.finite(actual) & is.finite(reference)
    other <- !is.na(reference) & !is.finite(reference)
    infinite <- identical(is.infinite(actual), is.infinite(reference)) && all(actual[other] == reference[other])
    zero <- identical(!is.na(actual) & actual == 0, !is.na(reference) & reference == 0)
    if (jgw2_direction(name)) direction <- all(sign(actual[finite]) == sign(reference[finite]))
    delta <- abs(actual[finite] - reference[finite])
    above <- sum(delta > 1e-10 + 1e-8 * abs(reference[finite]))
    if (length(delta)) max_abs <- max(delta)
    if (jgw2_probability(name)) {
      range <- all(is.na(actual) | (is.finite(actual) & actual >= 0 & actual <= 1)) &&
        all(is.na(reference) | (is.finite(reference) & reference >= 0 & reference <= 1))
      threshold <- identical(actual <= .05, reference <= .05) && identical(actual < .05, reference < .05)
      positive <- finite & actual > 0 & reference > 0
      differences <- abs(log10(actual[positive]) - log10(reference[positive]))
      above_log <- sum(differences > 1e-6)
      if (length(differences)) max_log <- max(differences)
    }
  } else if (types) {
    if (is.character(actual) || is.factor(actual)) {
      text_equal <- identical(jgw2_utf8(actual), jgw2_utf8(reference)) &&
        (!is.factor(actual) || identical(jgw2_utf8(levels(actual)), jgw2_utf8(levels(reference))))
    } else text_equal <- identical(actual, reference)
  }
  hard <- isTRUE(types && masks && infinite && zero && direction && threshold && range && text_equal)
  list(hard_pass = hard, types_exact = types, missing_NaN_exact = masks, infinity_exact = infinite,
    zero_mask_exact = zero, scientific_direction_exact = direction, both_005_masks_exact = threshold,
    probability_range = range, UTF8_text_exact = text_equal, n = length(reference),
    original_absrel_exceeded_n = above, original_log10P_exceeded_n = above_log,
    max_abs = max_abs, max_log10P = max_log,
    original_aligned_continuous_bound_pass = hard && above == 0L && above_log == 0L)
}
jgw2_compare_table <- function(actual, reference, endpoint) {
  a <- jgw2_plain(actual); b <- jgw2_plain(reference)
  exclude <- c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high")
  a <- structure(unclass(a)[!names(a) %in% exclude], class = "data.frame", row.names = attr(a, "row.names"))
  b <- structure(unclass(b)[!names(b) %in% exclude], class = "data.frame", row.names = attr(b, "row.names"))
  key <- if ("raw_gene_id" %in% names(b)) "raw_gene_id" else if ("term_id" %in% names(b)) "term_id" else "feature"
  jgw2_need(identical(names(a), names(b)) && nrow(a) == nrow(b) && key %in% names(b) &&
    is.character(a[[key]]) && is.character(b[[key]]), "full_endpoint_schema")
  ak <- jgw2_utf8(a[[key]]); bk <- jgw2_utf8(b[[key]])
  jgw2_need(!anyNA(ak) && !anyNA(bk) && all(nzchar(ak)) && all(nzchar(bk)) &&
    !anyDuplicated(ak) && !anyDuplicated(bk) && setequal(ak, bk), "complete_key_bijection")
  alignment <- match(bk, ak)
  rows <- lapply(seq_along(b), function(i) {
    result <- jgw2_compare_column(a[[i]][alignment], b[[i]], names(b)[i])
    data.frame(endpoint = endpoint, column = paste0(names(b)[i], "[", i, "]"),
      original_key_order_exact = identical(ak, bk), reordered_read_view_n = sum(alignment != seq_along(bk)),
      as.data.frame(result, stringsAsFactors = FALSE), check.names = FALSE)
  })
  do.call(rbind, rows)
}

jgw2_arithmetic_row <- function(method, contrast, field, actual, expected, probability = FALSE) {
  data.frame(method = method, contrast = contrast, field = field,
    as.data.frame(jgw2_math(actual, expected, probability)), stringsAsFactors = FALSE)
}
jgw2_limma <- function(table, fit, label) {
  table <- jgw2_plain(table); index <- if (label == "high") 1L else 2L
  jgw2_need(is.matrix(fit$coefficients) && is.matrix(fit$stdev.unscaled) &&
    identical(rownames(fit$coefficients), table$raw_gene_id) &&
    identical(dimnames(fit$coefficients), dimnames(fit$stdev.unscaled)) && index <= ncol(fit$coefficients),
    "saved_voom_fit_gene_coefficient_domain")
  effect <- fit$coefficients[, index]; unscaled <- fit$stdev.unscaled[, index]
  variance <- fit$s2.post; degrees <- fit$df.total
  if (length(variance) == 1L) variance <- rep(variance, nrow(table))
  if (length(degrees) == 1L) degrees <- rep(degrees, nrow(table))
  jgw2_need(length(variance) == nrow(table) && length(degrees) == nrow(table) &&
    all(is.finite(variance) & variance >= 0) && all(!is.na(degrees) & degrees > 0), "saved_voom_variance_df")
  se <- unscaled * sqrt(variance)
  statistic <- effect / unscaled / sqrt(variance)
  p <- 2 * stats::pt(abs(statistic), df = degrees, lower.tail = FALSE)
  critical <- stats::qt(.975, df = degrees)
  expected <- list(logFC = effect, t = statistic, P.Value = p,
    adj.P.Val = jgw2_bh(p), CI_low = effect - critical * se, CI_high = effect + critical * se)
  do.call(rbind, lapply(names(expected), function(name) jgw2_arithmetic_row("limma_voom_saved_fit", label,
    name, table[[name]], unname(expected[[name]]), name %in% c("P.Value", "adj.P.Val"))))
}
jgw2_ql_saved <- function(result, fixed, label) {
  table <- result$gene_edgeR[[label]]; test <- result$native_edgeR_tests[[label]]
  jgw2_need(is.data.frame(test$table) && identical(rownames(test$table), table$raw_gene_id) &&
    identical(table$raw_gene_id, fixed$annotation$raw_gene_id) &&
    identical(result$native_edgeR_fits[[label]]$design, result$designs[[label]]), "saved_QL_test_design_domain")
  fields <- intersect(c("logFC", "logCPM", "F", "PValue"), names(test$table))
  jgw2_need(all(c("logFC", "F", "PValue") %in% fields), "native_QL_test_columns")
  rows <- lapply(fields, function(name) jgw2_arithmetic_row("edgeR_saved_test_extraction", label,
    name, table[[name]], test$table[[name]], name == "PValue"))
  rows[[length(rows) + 1L]] <- jgw2_arithmetic_row("edgeR_BH", label, "FDR", table$FDR, jgw2_bh(table$PValue), TRUE)
  do.call(rbind, rows)
}

# Independent transcription of the installed 1.52.0 selection rule (text pins
# are authenticated by the CLI). No DESeq2 adjustment/results function is called.
jgw2_independent_filter <- function(base_mean, p, alpha = .05) {
  jgw2_need(is.numeric(base_mean) && is.numeric(p) && length(base_mean) == length(p) &&
    length(p) > 0L && all(is.finite(base_mean) & base_mean >= 0) && identical(alpha, .05), "native_IF_inputs_alpha")
  lower <- mean(base_mean == 0)
  theta <- seq(lower, if (lower < .95) .95 else 1, length.out = 50L)
  cutoffs <- stats::quantile(base_mean, probs = theta)
  adjusted <- matrix(NA_real_, length(p), 50L)
  for (j in seq_len(50L)) {
    use <- base_mean >= cutoffs[j]
    if (any(use)) adjusted[use, j] <- jgw2_bh(p[use])
  }
  rejections <- colSums(adjusted < alpha, na.rm = TRUE)
  # This is the original deterministic P-adjustment selection smoother, not
  # another expression model or native fit. Its input is saved P/baseMean only.
  smooth <- stats::lowess(theta, rejections, f = 1 / 5)
  if (max(rejections) <= 10) chosen <- 1L else {
    residual <- if (all(rejections == 0)) 0 else rejections[rejections > 0] - smooth$y[rejections > 0]
    maximum <- max(smooth$y)
    threshold <- maximum - sqrt(mean(residual^2))
    chosen <- if (any(rejections > threshold)) which(rejections > threshold)[1L] else
      if (any(rejections > .9 * maximum)) which(rejections > .9 * maximum)[1L] else
        if (any(rejections > .8 * maximum)) which(rejections > .8 * maximum)[1L] else 1L
  }
  list(padj = adjusted[, chosen], selected_index = chosen,
    filter_threshold = unname(cutoffs[chosen]), filter_theta = theta[chosen],
    rejection_counts = unname(rejections), alpha = alpha, adjusted_family_n = sum(base_mean >= cutoffs[chosen]),
    scope = "Independent saved-P native IF arithmetic; no results() or refit")
}
jgw2_deseq_fields <- function(dds) {
  list(counts = DESeq2::counts(dds), size_factors = DESeq2::sizeFactors(dds),
    normalization_factors = DESeq2::normalizationFactors(dds),
    columns = as.data.frame(SummarizedExperiment::colData(dds)),
    row_metadata = as.data.frame(S4Vectors::mcols(dds)),
    coefficient_names = DESeq2::resultsNames(dds),
    dispersion_design = attr(dds, "dispModelMatrix"),
    test = attr(dds, "test"), beta_prior = attr(dds, "betaPrior"),
    model_matrix_type = attr(dds, "modelMatrixType"),
    fit_type = attr(DESeq2::dispersionFunction(dds), "fitType"))
}
jgw2_deseq_expected <- function(fields, labels, label) {
  counts <- fields$counts; meta <- fields$row_metadata; samples <- fields$columns
  jgw2_need(is.matrix(counts) && is.numeric(counts) && ncol(counts) == 12L &&
    identical(colnames(counts), names(labels)) && identical(rownames(meta), rownames(counts)) &&
    identical(rownames(samples), colnames(counts)) && all(is.finite(counts) & counts >= 0 & counts == round(counts)) &&
    identical(fields$test, "Wald") && identical(fields$beta_prior, FALSE) &&
    identical(fields$model_matrix_type, "standard") && is.null(fields$normalization_factors) &&
    all(!is.na(labels) & labels %in% 0:3), "saved_DDS_native_default_domain")
  sf <- fields$size_factors
  jgw2_need(is.numeric(sf) && length(sf) == 12L && all(is.finite(sf) & sf > 0), "saved_sizeFactors")
  normal_mean <- rowMeans(sweep(counts, 2L, sf, "/"))
  expected_design <- if (label == "high") cbind(1, as.numeric(labels == 1), as.numeric(labels == 2), as.numeric(labels == 3)) else
    cbind(1, unname(labels))
  dm <- fields$dispersion_design
  jgw2_need(is.matrix(dm) && identical(dim(dm), dim(expected_design)) &&
    identical(rownames(dm), names(labels)) && all(dm == expected_design) &&
    !"tDegreesFreedom" %in% names(meta) &&
    all(c("baseMean", "allZero", "maxCooks", "betaConv") %in% names(meta)) &&
    is.logical(meta$allZero) && !anyNA(meta$allZero) &&
    identical(meta$allZero, unname(rowSums(counts) == 0)) &&
    (is.null(meta$replace) || !any(meta$replace, na.rm = TRUE)), "saved_DDS_design_and_no_replacement")
  if (label == "high") {
    groups <- c("vehicle", "low", "medium", "high")
    jgw2_need(is.factor(samples$dose_group) && identical(levels(samples$dose_group), groups) &&
      identical(as.character(samples$dose_group), groups[labels + 1L]), "DDS_four_group_factor")
    coefficient <- "dose_group_high_vs_vehicle"
  } else {
    jgw2_need(identical(as.numeric(samples$ordered_dose), unname(labels)), "DDS_ordered_dose")
    coefficient <- "ordered_dose"
  }
  names <- c(coefficient, paste0("SE_", coefficient), paste0("WaldStatistic_", coefficient), paste0("WaldPvalue_", coefficient))
  jgw2_need(coefficient %in% fields$coefficient_names && all(names %in% names(meta)), "saved_DDS_target_coefficient")
  effect <- meta[[names[1L]]]; se <- meta[[names[2L]]]; statistic <- meta[[names[3L]]]; p <- meta[[names[4L]]]
  zero <- rep(FALSE, nrow(counts))
  if (label == "high") zero <- rowSums(counts[, labels %in% c(0, 3), drop = FALSE]) == 0 & !meta$allZero
  # cleanContrast's standard base-reference branch copies saved coefficient/SE,
  # then resets only contrast-zero LFC/stat/P. No contrast refit is reproduced.
  export_effect <- effect; export_statistic <- statistic; export_p <- p
  export_effect[zero] <- 0; export_statistic[zero] <- 0; export_p[zero] <- 1
  cooks_cutoff <- stats::qf(.99, ncol(dm), nrow(dm) - ncol(dm))
  cooks <- meta$maxCooks > cooks_cutoff
  # The source's don't-filter exception is for a TWO-level factor only. These
  # frozen designs have four factor levels or one numeric dose, so it is absent.
  export_p[which(cooks)] <- NA_real_
  filtering <- jgw2_independent_filter(meta$baseMean, export_p, .05)
  list(baseMean = meta$baseMean, log2FoldChange = export_effect, lfcSE = se,
    stat = export_statistic, pvalue = export_p, padj = filtering$padj,
    raw_effect = effect, raw_se = se, raw_statistic = statistic, raw_p = p,
    normalized_mean = normal_mean, contrast_zero = zero, Cook_outlier = cooks,
    Cook_cutoff = cooks_cutoff, filtering = filtering,
    fit_type = fields$fit_type, betaConv = meta$betaConv,
    gene_ids = rownames(counts), coefficient = coefficient)
}
jgw2_deseq <- function(export, saved_res, fields, fixed, labels, label) {
  expected <- jgw2_deseq_expected(fields, labels, label)
  jgw2_need(identical(export$raw_gene_id, expected$gene_ids) &&
    identical(rownames(saved_res), expected$gene_ids) &&
    identical(dimnames(fields$counts), dimnames(fixed$full_counts[fixed$keep, , drop = FALSE])) &&
    all(fields$counts == round(fixed$full_counts[fixed$keep, , drop = FALSE])), "DDS_export_count_gene_sample_alignment")
  rows <- lapply(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), function(name) {
    jgw2_need(identical(export[[name]], saved_res[[name]]), "unmodified_native_export_equals_saved_res")
    jgw2_arithmetic_row("DESeq2_saved_native_rule", label, name, export[[name]], expected[[name]],
      name %in% c("pvalue", "padj"))
  })
  rows[[length(rows) + 1L]] <- jgw2_arithmetic_row("DESeq2_normalization", label, "baseMean",
    fields$row_metadata$baseMean, unname(expected$normalized_mean))
  rows[[length(rows) + 1L]] <- jgw2_arithmetic_row("DESeq2_raw_Wald", label, "statistic_ratio",
    expected$raw_statistic, expected$raw_effect / expected$raw_se)
  normal_p <- 2 * stats::pnorm(abs(expected$raw_statistic), lower.tail = FALSE)
  rows[[length(rows) + 1L]] <- jgw2_arithmetic_row("DESeq2_raw_Wald", label, "normal_two_sided_P",
    expected$raw_p, normal_p, TRUE)
  list(checks = do.call(rbind, rows), diagnostics = list(fitType = expected$fit_type,
    betaConv_n = length(expected$betaConv), betaConv_missing = sum(is.na(expected$betaConv)),
    betaConv_false = sum(!is.na(expected$betaConv) & !expected$betaConv),
    contrast_zero_n = sum(expected$contrast_zero), Cook_outlier_n = sum(expected$Cook_outlier, na.rm = TRUE),
    Cook_unavailable_n = sum(is.na(expected$Cook_outlier)), Cook_cutoff = expected$Cook_cutoff,
    independent_filter = expected$filtering[setdiff(names(expected$filtering), "padj")]))
}

jgw2_family_adjustments <- function(native, field_map, old) {
  rows <- field_map[field_map$accession == native$accession & field_map$space == "four_group_MC", , drop = FALSE]
  answer <- lapply(seq_len(nrow(rows)), function(i) {
    row <- as.list(rows[i, , drop = FALSE])
    if (row$method == "DESeq2_native_IF_Cooks") return(NULL) # independently handled with saved DDS
    table <- old$jgw6_native_table(native, row)
    expected <- if (row$method == "mroast1999_equivalent_not_historical_authentication")
      jgw2_midp(table[[row$native_P]]) else jgw2_bh(table[[row$native_P]])
    value <- jgw2_arithmetic_row(row$method, row$contrast, paste(row$collection, row$native_q, sep = ":"),
      table[[row$native_q]], expected, TRUE)
    if (row$method == "mroast1999_equivalent_not_historical_authentication") value <- rbind(value,
      jgw2_arithmetic_row(row$method, row$contrast, paste(row$collection, "FDR.Mixed", sep = ":"),
        table$FDR.Mixed, jgw2_midp(table$PValue.Mixed), TRUE))
    value
  })
  do.call(rbind, answer)
}
jgw2_installed_source <- function(path, name) {
  syntax <- parse(path, keep.source = FALSE)
  jgw2_need(length(syntax) == 1L && is.call(syntax[[1L]]) &&
    identical(syntax[[1L]][[1L]], as.name("function")), "one_installed_function_text")
  actual <- get0(name, envir = asNamespace("DESeq2"), inherits = FALSE)
  jgw2_need(is.function(actual), "installed_function_present")
  clean <- function(x) parse(text = deparse(x, control = "all"), keep.source = FALSE)[[1L]]
  identical(clean(formals(actual)), clean(syntax[[1L]][[2L]])) &&
    identical(clean(body(actual)), clean(syntax[[1L]][[3L]]))
}

jgw2_run <- function(config_path, accession, id, output) {
  jgw2_need(accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    is.numeric(id) && length(id) == 1L && !is.na(id) && id %in% 1:2 &&
    !file.exists(output) && dir.exists(dirname(output)), "one_object_new_output")
  run <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03"
  jgw2_need(identical(config_path, paste0(run, "/CONFIG.json")), "explicit_v2_config")
  tracked <- list(); warning_classes <- character()
  track <- function(rec) {
    pin <- jgw2_record(rec); tracked[[pin$path]] <<- pin; pin
  }
  read_meta <- function(path) {
    pin <- track(jgw2_pin(path)); jgw2_json(pin$path)
  }
  read_saved <- function(rec) {
    pin <- track(rec)
    withCallingHandlers(readRDS(pin$path), warning = function(condition) {
      warning_classes <<- unique(c(warning_classes, class(condition)[1L]))
      invokeRestart("muffleWarning") # classes are retained, never raw text
    })
  }
  cfg_pin <- track(jgw2_pin(config_path)); cfg <- jgw2_json(config_path)
  project <- paste0(run, "/project")
  jgw2_need(identical(cfg$schema, "jgw03-separate-linux-reference-config-v2") &&
    identical(cfg$run_root, run) && identical(cfg$project_root, project) &&
    setequal(names(cfg$cells), c("GSE289929", "GSE227494", "GSE299063")), "v2_scientific_scope")
  definition <- track(cfg$compatibility_definition); approval <- track(cfg$approval)
  jgw2_need(identical(definition$sha256, "5c06cb8e054c4b2e3a261c2250944c09939d6b6851a65043222ac2360de3aff8"),
    "approved_non_equivalence_contract")
  source_names <- names(cfg$source_pins)
  jgw2_need(length(source_names) > 0L && !anyDuplicated(source_names) && all(nzchar(source_names)) &&
    all(!grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))", source_names)), "source_relative_paths")
  sources <- lapply(source_names, function(name) {
    rec <- cfg$source_pins[[name]]; rec$path <- file.path(project, name); track(rec)
  }); names(sources) <- source_names
  known <- c(
    "code/R/check_JGW03_linux_six_bridges_20261017.R" = "ab2c820a86e0acd22e2c0a90bfd8884e8251e919b2565c418b8de6da973f007d",
    "code/R/prepare_JGW03_observed_label_pilot_v3_20261015.R" = "fe01d32bf79239cc51a5b4b75140c3ad209b31501adc94f141623c9ecf991f38",
    "code/R/JGW03_label_reference_core_20261015.R" = "e3939f371f02b24cfdfefee0a9cd9e39bdf47c8920caa44f58434a93f7a2974e",
    "code/R/JGW03_reference_field_map_20261015.R" = "c349eab9c5f3718d95a06c0b51e44584f4adb60ecbfb2f2d75235af7bfd7da42",
    "code/R/JGW03_reference_compact_tail_20261015.R" = "1074cc5001899cf4f2500c72a52f6f218b1cbb763982c6a0a9b43bfeab82ce8d",
    "code/R/analyze_jgw03_three_cell_lines.R" = "5df6378ed4dadfb4050d91e252a63af4f44d2be80d3070d889d1e0e007796775",
    "code/R/JGW03_input_admission_v2_20261015.R" = "63a24a37bac60e1e1c99cd11e9a69c6a11ab8dbb70c44766d06f2e33b87e9c0c",
    "code/R/run_JGW03_full_reference_linux_20261017.R" = "57e3643a06b5c17f66c46eb658d150665c9bc55a1b4cc025ceb3e21d6d2cdc8f")
  producer_name <- "code/R/run_JGW03_full_reference_linux_v3_20261017.R"
  checker_name <- "code/R/check_JGW03_linux_native_compatibility_v3_20261017.R"
  jgw2_need(all(c(names(known), producer_name, checker_name) %in% source_names) &&
    all(vapply(names(known), function(name) identical(sources[[name]]$sha256, known[[name]]), logical(1))) &&
    identical(cfg$checker_source, checker_name), "original_science_and_versioned_entry_pins")
  script <- grep("^--file=", commandArgs(), value = TRUE)
  jgw2_need(length(script) == 1L && identical(jgw2_pin(sub("^--file=", "", script)), sources[[checker_name]]),
    "actual_checker_source")
  adjustment_sha <- c(pvalueAdjustment = "29c82f8d183cdedbbf2d50a1a8e1219721c155c587476eb15dd50ce377712249",
    results = "7d0738387acf17c13a3173d91300803fbd99b8d39878b8f9cc99b889ccfc1ef0",
    cleanContrast = "6a88ce5c9dc780ea14e24155bccdd87e9208748cf8e0c15cacaba112a611246b",
    filtered_p = "d48b88d38d2d4380bf39cd388b81cf5a605270c01ea72ad11cfa0a8fe2f4ef0d")
  jgw2_need(setequal(names(cfg$native_adjustment_sources), names(adjustment_sha)), "four_installed_adjustment_sources")
  adjustment_sources <- lapply(names(adjustment_sha), function(name) {
    rec <- track(cfg$native_adjustment_sources[[name]])
    jgw2_need(identical(rec$sha256, adjustment_sha[[name]]) && jgw2_installed_source(rec$path, name),
      "native_adjustment_source_version_and_live_definition")
    rec
  }); names(adjustment_sources) <- names(adjustment_sha)
  jgw2_need(identical(as.character(utils::packageVersion("DESeq2")), "1.52.0"), "actual_native_DESeq2_version")
  # Reuse only already independent structural helpers, never its old table
  # comparator/run or the producer comparison/extraction functions.
  old <- new.env(parent = baseenv())
  jgw2_import(sources[[names(known)[1L]]]$path, c("jgw6_need", "jgw6_hash", "jgw6_plain",
    "jgw6_probability", "jgw6_direction", "jgw6_vector", "jgw6_fixed", "jgw6_model", "jgw6_labels",
    "jgw6_registry", "jgw6_native_table", "jgw6_counts", "jgw6_outer", "jgw6_compact",
    "jgw6_json_same", "jgw6_saved_deseq"), old)
  getter <- new.env(parent = baseenv())
  jgw2_import(sources[[names(known)[2L]]]$path, c("jgwp_need", "jgwp_fresh_tables", "jgwp_expected_tables"), getter)
  jgw2_import(sources[[names(known)[3L]]]$path, c("jgw03_require", "jgw03_family_registry"), getter)
  jgw2_import(sources[[names(known)[4L]]]$path, "jgw03_reference_field_map", getter)
  jgw2_import(sources[[names(known)[5L]]]$path, "jgw03c_definition", getter)
  family_map <- getter$jgw03_reference_field_map(getter$jgw03_family_registry())
  jgw2_need(nrow(family_map) == 123L && sum(family_map$space == "four_group_MC") == 89L &&
    sum(family_map$space == "high_companion_exact") == 34L, "unchanged_global_family_registry")
  inputs <- lapply(cfg$input_pins, track)
  cell <- cfg$cells[[accession]]
  references <- lapply(cell[c("registry", "windows_observed", "windows_nonidentity", "windows_fixed")], track)
  windows_receipts <- lapply(cell$windows_receipts, function(rec) { pin <- track(rec); jgw2_json(pin$path) })
  wp <- windows_receipts$prepared; wo <- windows_receipts$observed_complete
  wn <- windows_receipts$nonidentity_complete; ws <- windows_receipts$nonidentity_saved
  jgw2_need(jgw2_same_pin(wp$label_registry, references$registry) &&
    jgw2_same_pin(wo$result, references$windows_observed) &&
    jgw2_same_pin(wp$observed_files$fixed, references$windows_fixed) &&
    jgw2_same_pin(wn$native, references$windows_nonidentity) && jgw2_same_pin(ws$native, references$windows_nonidentity) &&
    identical(wo$accession, accession) && identical(wn$accession, accession) && wn$allocation_id == 2 &&
    identical(wn$prepared_sha256, wp$prepared$sha256) && identical(ws$prepared_sha256, wp$prepared$sha256) &&
    identical(wn$allocation_key, ws$allocation_key), "Windows_comparison_lineage")
  leaf <- paste0(run, "/intermediate/", accession)
  preparation_started <- read_meta(paste0(leaf, "/PREPARATION_STARTED.json"))
  prepared <- read_meta(paste0(leaf, "/PREPARED.json")); prepared_pin <- track(prepared$prepared)
  context <- read_saved(prepared_pin)
  jgw2_need(identical(context$schema, "jgw03-separate-linux-preparation-v2") &&
    identical(context$accession, accession) && identical(context$root, project) && identical(context$run_root, run) &&
    identical(prepared$config_sha256, cfg_pin$sha256) && identical(prepared$entry_sha256, sources[[producer_name]]$sha256) &&
    jgw2_same_pin(context$binding$config, cfg_pin) && jgw2_same_pin(context$entry, sources[[producer_name]]) &&
    old$jgw6_json_same(context$binding$sources, sources) && old$jgw6_json_same(context$binding$inputs, inputs) &&
    old$jgw6_json_same(context$binding$references, references) &&
    old$jgw6_json_same(preparation_started$binding, context$binding) &&
    context$full_MC_nonidentity_draws == 10000 && context$full_companion_including_observed == 20 &&
    context$full_MC_family_n == 89 && context$full_companion_family_n == 34 && context$native_allocation_calls == 0,
    "Linux_v2_preparation_input_identity")
  for (name in c("fixed", "registry", "session_info")) {
    track(context[[name]]); jgw2_need(jgw2_same_pin(context[[name]], prepared[[name]]), "prepared_product_pin")
  }
  track(context$fixed_bridge)
  for (rec in context$sources) track(rec)
  runtime_fields <- c("R", "packages", "platform", "BLAS", "backend", "workers", "seed", "RNGkind")
  jgw2_need(old$jgw6_json_same(context$runtime[runtime_fields], cfg$expected_runtime[runtime_fields]) &&
    identical(context$runtime$R, "4.6.1") && grepl("linux", context$runtime$platform) &&
    all(unlist(context$runtime$thread_environment, use.names = FALSE) == "1"), "saved_native_runtime")
  registry <- read_saved(context$registry)
  jgw2_need(jgw2_same_pin(context$registry, references$registry), "original_registry_byte_identity")
  registry_scope <- old$jgw6_registry(registry, context$observed_labels)
  labels <- setNames(as.numeric(registry$MC[id, ]), colnames(registry$MC))
  fixed <- read_saved(context$fixed); windows_fixed <- read_saved(references$windows_fixed)
  fixed_checks <- old$jgw6_fixed(fixed, windows_fixed)
  bridge_leaf <- paste0(leaf, "/bridge_", id)
  saved <- read_meta(paste0(bridge_leaf, "/FIT_SAVED.json"))
  complete <- read_meta(paste0(bridge_leaf, "/BRIDGE_COMPLETE.json"))
  jgw2_need(identical(saved$allocation_key, registry_scope$keys[id]) &&
    identical(saved$prepared_sha256, prepared_pin$sha256) && identical(saved$config_sha256, cfg_pin$sha256) &&
    identical(saved$entry_sha256, sources[[producer_name]]$sha256) &&
    identical(complete$accession, accession) && complete$allocation_id == id &&
    identical(complete$allocation_key, saved$allocation_key) &&
    identical(complete$prepared_sha256, prepared_pin$sha256) && identical(complete$config_sha256, cfg_pin$sha256) &&
    identical(complete$entry_sha256, sources[[producer_name]]$sha256) &&
    jgw2_same_pin(saved$native, complete$native), "native_v2_saved_complete_lineage")
  # Both producer statuses may be inspected: independent hard flags decide this
  # result. An independent PASS never changes the producer/old strict status.
  jgw2_need(complete$status %in% c("LINUX_NATIVE_COMPATIBILITY_PASS_V2", "LINUX_NATIVE_COMPATIBILITY_REVIEW_REQUIRED_V2"),
    "v2_producer_receipt_schema")
  for (name in c("native", "compact", "checks", "semantic_checks", "aligned_strict_checks", "checks_saved", "session_info"))
    track(complete[[name]])
  track(saved$session_info)
  reused <- accession == "GSE289929" && id == 1L
  origin_summary <- NULL
  if (reused) {
    origin <- cfg$reused_first_PC12
    jgw2_need(identical(origin$origin_run_root, "/root/autodl-tmp/JGW03_full_reference_linux_20261017r01"), "explicit_r01_origin")
    origin_records <- lapply(origin[setdiff(names(origin), "origin_run_root")], track)
    om <- lapply(origin_records[c("config", "prepared_receipt", "fit_saved", "bridge_complete", "process_exit",
      "process_command", "process_started")],
      function(rec) jgw2_json(rec$path))
    op <- read_saved(origin_records$prepared)
    jgw2_need(identical(saved$kind, "REUSED_ORIGINAL_LINUX_NATIVE") && saved$new_native_calls == 0 &&
      jgw2_same_pin(saved$native, origin_records$native) &&
      jgw2_same_pin(saved$origin$records$native, origin_records$native) &&
      jgw2_same_pin(om$prepared_receipt$prepared, origin_records$prepared) &&
      jgw2_same_pin(om$fit_saved$native, origin_records$native) && jgw2_same_pin(om$bridge_complete$native, origin_records$native) &&
      identical(om$fit_saved$config_sha256, origin_records$config$sha256) &&
      identical(om$fit_saved$prepared_sha256, origin_records$prepared$sha256) &&
      identical(om$fit_saved$entry_sha256, known[["code/R/run_JGW03_full_reference_linux_20261017.R"]]) &&
      identical(om$bridge_complete$status, "LINUX_NATIVE_BRIDGE_REVIEW_REQUIRED") &&
      om$process_exit$actual_R_exit_code == 2 && saved$origin$actual_original_R_exit_code == 2 &&
      jgw2_same_pin(om$process_command$config, origin_records$config) &&
      identical(om$process_command$entry$sha256, om$fit_saved$entry_sha256) &&
      identical(om$process_started$config_sha256, origin_records$config$sha256) &&
      is.numeric(om$process_started$start_ticks) && om$process_started$start_ticks > 0 &&
      om$process_started$actual_R_pid == om$process_exit$actual_R_pid &&
      saved$origin$actual_original_R_pid == om$process_exit$actual_R_pid &&
      saved$origin$actual_original_R_start_ticks == om$process_started$start_ticks &&
      identical(saved$origin$source_generation_entry_sha256, om$fit_saved$entry_sha256) &&
      jgw2_same_pin(saved$origin$new_linkage_entry, sources[[producer_name]]) &&
      old$jgw6_json_same(saved$execution, om$fit_saved$execution) &&
      old$jgw6_json_same(op$runtime, context$runtime) &&
      identical(op$registry$sha256, context$registry$sha256), "preserved_exit2_no_refit_lineage")
    origin_fixed <- read_saved(origin_records$fixed)
    science <- setdiff(names(fixed), "input_admission_v2")
    jgw2_need(identical(science, setdiff(names(origin_fixed), "input_admission_v2")) &&
      identical(fixed[science], origin_fixed[science]), "reused_native_exact_original_fixed_values")
    for (name in names(om$config$source_pins)) jgw2_need(jgw2_same_pin(om$config$source_pins[[name]], cfg$source_pins[[name]]),
      "reused_original_scientific_source_identity")
    origin_summary <- list(kind = saved$kind, actual_original_R_exit_code = 2L, original_strict_status = om$bridge_complete$status,
      origin_records = origin_records, new_native_calls = 0L)
    rm(op, origin_fixed)
  } else {
    jgw2_need(identical(saved$kind, "NEW_LINUX_NATIVE_V2") && saved$new_native_calls == 1 &&
      identical(saved$native$path, paste0(bridge_leaf, "/NATIVE_CAPSULE.rds")) &&
      identical(saved$source_generation_entry_sha256, sources[[producer_name]]$sha256), "new_native_own_v2_identity")
    started <- read_meta(paste0(bridge_leaf, "/FIT_STARTED.json"))
    jgw2_need(identical(started$prepared_sha256, prepared_pin$sha256) &&
      identical(started$allocation_key, saved$allocation_key) && jgw2_same_pin(started$entry, sources[[producer_name]]), "new_fit_started_identity")
  }
  native <- read_saved(saved$native)
  reference <- if (id == 1L) references$windows_observed else references$windows_nonidentity
  windows <- read_saved(reference)
  jgw2_need(identical(native$accession, accession) && identical(windows$accession, accession) &&
    identical(native$labels, labels) && identical(windows$labels, labels) &&
    identical(native$native_source_sha256, known[["code/R/analyze_jgw03_three_cell_lines.R"]]), "saved_native_label_source_identity")
  measured_library_counts <- list(before = length(windows$labels), after = length(native$labels),
    fixed_count_columns = ncol(fixed$y$counts), unique_native_sample_ids = length(unique(names(native$labels))),
    same_library_id_order = identical(names(native$labels), names(windows$labels)), deleted = 0L)
  model_checks <- old$jgw6_model(native, windows, labels, accession)
  if (id == 1L) domain <- read_saved(complete$domain) else {
    first <- read_meta(paste0(leaf, "/bridge_1/BRIDGE_COMPLETE.json"))
    jgw2_need(identical(first$config_sha256, cfg_pin$sha256) && identical(first$prepared_sha256, prepared_pin$sha256), "same_Linux_observed_domain")
    domain <- read_saved(first$domain)
  }
  compact <- read_saved(complete$compact)
  family_checks <- old$jgw6_compact(native, compact, domain, family_map, getter$jgw03c_definition())
  jgw2_need(identical(domain$observed_labels, context$observed_labels), "Linux_observed_not_Windows_threshold")
  execution <- saved$execution; execution$RNGkind <- unlist(execution$RNGkind, use.names = FALSE)
  jgw2_need(old$jgw6_json_same(compact$execution_metadata$actual, execution) && execution$seed == 20260817 &&
    execution$workers == 1 && identical(execution$backend, "SerialParam") && execution$mroast_nrot == 1999 && execution$fgsea_eps == 0 &&
    identical(execution$RNGkind, c("Mersenne-Twister", "Inversion", "Rejection")) &&
    identical(execution$mroast_executed, accession == "GSE289929") &&
    identical(execution$historical_backend_authenticated, FALSE) &&
    identical(compact$execution_metadata$inner_MC_uncertainty_covered, FALSE), "saved_execution_not_RNG_consumption_proof")
  previous_locale <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", previous_locale), add = TRUE)
  locale <- Sys.setlocale("LC_CTYPE", "C.UTF-8")
  jgw2_need(nzchar(locale), "UTF8_comparison_read_view_available")
  a <- getter$jgwp_fresh_tables(fixed, native); b <- getter$jgwp_fresh_tables(windows_fixed, windows)
  expected_endpoints <- basename(getter$jgwp_expected_tables(accession))
  jgw2_need(setequal(names(a), expected_endpoints) && setequal(names(b), expected_endpoints) &&
    !anyDuplicated(names(a)) && !anyDuplicated(names(b)), "complete_native_auxiliary_endpoint_roster")
  columns <- do.call(rbind, lapply(names(b), function(name) jgw2_compare_table(a[[name]], b[[name]], name)))
  arithmetic <- list(jgw2_family_adjustments(native, family_map, old))
  for (label in c("high", "trend", "logtrend")) arithmetic[[length(arithmetic) + 1L]] <- jgw2_ql_saved(native, fixed, label)
  diagnostics <- list(); windows_diagnostics <- list()
  for (label in c("high", "trend")) {
    arithmetic[[length(arithmetic) + 1L]] <- jgw2_limma(native$gene_voom[[label]]$value,
      native$gene_voom[[label]]$native_locals$fit, label)
    fields <- jgw2_deseq_fields(native$gene_DESeq2[[label]]$native_locals$dds)
    checked <- jgw2_deseq(native$gene_DESeq2[[label]]$value, native$gene_DESeq2[[label]]$native_locals$res,
      fields, fixed, labels, label)
    arithmetic[[length(arithmetic) + 1L]] <- checked$checks
    diagnostics[[label]] <- checked$diagnostics
    prior <- jgw2_deseq_fields(windows$gene_DESeq2[[label]]$native_locals$dds)
    windows_diagnostics[[label]] <- list(fitType = prior$fit_type, betaConv_n = length(prior$row_metadata$betaConv),
      betaConv_missing = sum(is.na(prior$row_metadata$betaConv)),
      betaConv_false = sum(!is.na(prior$row_metadata$betaConv) & !prior$row_metadata$betaConv))
    diagnostics[[label]]$Windows_fitType_same <- identical(fields$fit_type, prior$fit_type)
    diagnostics[[label]]$Windows_betaConv_flags_same <- identical(fields$row_metadata$betaConv, prior$row_metadata$betaConv)
    jgw2_need(is.character(fields$fit_type) && length(fields$fit_type) == 1L && nzchar(fields$fit_type) &&
      is.logical(fields$row_metadata$betaConv) && length(fields$row_metadata$betaConv) == nrow(fields$counts), "native_convergence_diagnostic_domain")
  }
  arithmetic <- do.call(rbind, arithmetic)
  # The original compact's saved DESeq metadata is checked against actual DDS
  # getter values, not the producer's copied flags or a fresh results() object.
  native_diagnostics <- lapply(diagnostics, function(x) list(fitType = x$fitType,
    betaConv = list(status = "SAVED_NATIVE_FLAG_COUNTS", n = x$betaConv_n,
      false_n = x$betaConv_false, unavailable_n = x$betaConv_missing),
    native_results_metadata = "Original res converted to data.frame; no new results() reconstruction"))
  compact_DDS_exact <- identical(compact$saved_DESeq2_metadata, old$jgw6_saved_deseq(native, native_diagnostics))
  unchanged <- all(vapply(tracked, function(rec) identical(jgw2_record(rec), rec), logical(1)))
  passed <- all(unlist(fixed_checks)) && all(unlist(model_checks)) && all(columns$hard_pass) &&
    all(family_checks$pass) && all(arithmetic$pass) && compact_DDS_exact && unchanged
  jgw2_need(dir.create(output, recursive = FALSE), "exclusive_checker_output")
  tables <- list(COLUMN_COMPATIBILITY = columns, NATIVE_ARITHMETIC = arithmetic, FAMILY_COMPACT_CHECKS = family_checks)
  for (name in names(tables)) utils::write.table(tables[[name]], file.path(output, paste0(name, ".tsv")),
    sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  receipt <- list(schema = "jgw03-saved-native-compatibility-v2-independent-v1",
    status = if (passed) "PASS_SAVED_NATIVE_COMPATIBILITY_V2_NOT_EQUIVALENCE" else "FAIL_SAVED_NATIVE_COMPATIBILITY_V2",
    all_hard_invariants_pass = passed, accession = accession, allocation_id = as.integer(id), allocation_key = saved$allocation_key,
    config_sha256 = cfg_pin$sha256, producer_entry_sha256 = sources[[producer_name]]$sha256,
    prepared_sha256 = prepared_pin$sha256, native_sha256 = saved$native$sha256, checker_sha256 = sources[[checker_name]]$sha256,
    config = cfg_pin, checker = sources[[checker_name]], producer = sources[[producer_name]], prepared = prepared_pin,
    native = saved$native, comparison_only_Windows = reference, compatibility_definition = definition, approval = approval,
    native_generation = if (reused) origin_summary else list(kind = saved$kind, new_native_calls = 1L),
    original_producer_strict_status_or_flags_not_changed = list(producer_status = complete$status,
      strict_original_order_all_fields_pass = complete$strict_original_order_all_fields_pass,
      strict_aligned_UTF8_all_fields_pass = complete$strict_aligned_UTF8_all_fields_pass),
    original_absrel_exceeded_field_rows = sum(columns$original_absrel_exceeded_n > 0),
    original_logP_exceeded_field_rows = sum(columns$original_log10P_exceeded_n > 0),
    own_aligned_original_continuous_bounds_all_pass = all(columns$original_aligned_continuous_bound_pass),
    continuous_bounds_are_diagnostic_not_new_admission = TRUE,
    arithmetic_policy = jgw2_arithmetic_policy(), installed_native_adjustment_sources = adjustment_sources,
    fixed_checks = fixed_checks, model_checks = model_checks, native_DESeq2_diagnostics = diagnostics,
    Windows_DESeq2_diagnostics = windows_diagnostics, compact_saved_DDS_metadata_exact = compact_DDS_exact,
    endpoint_n = length(expected_endpoints), column_n = nrow(columns),
    probability_column_n = sum(vapply(sub("\\[[0-9]+\\]$", "", columns$column), jgw2_probability, logical(1))),
    global_source_MC_family_n = 89L, global_source_exact_family_n = 34L,
    this_object_full_MC_families_checked = nrow(family_checks),
    this_cell_source_exact_families = sum(family_map$accession == accession & family_map$space == "high_companion_exact"),
    registry = registry_scope, comparison_LC_CTYPE = locale,
    before_after_pins = unname(tracked), sources_inputs_native_unchanged = unchanged,
    before_after_libraries = measured_library_counts,
    warning_classes_only = warning_classes, raw_conditions_persisted = FALSE,
    no_native_fit_or_results_invoked = TRUE, original_values_or_tolerances_changed = FALSE,
    native_process_terminal_and_resource_acceptance = "SEPARATE_MAIN_AGGREGATION_REQUIRED",
    internal_validation_limits = c("Gene voom saved coefficients/SE/t/P/CI independently recomputed; no fitting",
      "DESeq2 saved coefficient/Wald/normal P/Cook/zero/independent filtering independently recomputed",
      "edgeR saved test extraction and original BH checked; Poisson-bound QLF internals not independently refitted",
      "All methods' existing adjustment domains checked; program/activity moderation internals were not retained",
      "CAMERA/fgsea stochastic internal P algorithm and complete RNG consumption are not re-executed"),
    warnings = c("Library observations do not authenticate independent subjects, batches or randomized exchangeability",
      "Two fixed labels do not establish full-null-tail equivalence; Linux observed/domain/null must remain together",
      "Native fit diagnostics are retained; no field or feature is excluded to obtain compatibility"),
    full_reference_accepted = FALSE, full_dispatch_or_resource_accepted = FALSE, clinical_QC = "NOT_APPLICABLE",
    outputs = lapply(names(tables), function(name) jgw2_pin(file.path(output, paste0(name, ".tsv")))))
  jsonlite::write_json(receipt, file.path(output, "COMPATIBILITY.json"), pretty = TRUE,
    auto_unbox = TRUE, digits = NA, null = "null", na = "null")
  invisible(receipt)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  code <- tryCatch({
    jgw2_need(length(args) == 4L && grepl("^[12]$", args[3L]), "CLI_CONFIG_ACCESSION_ID_OUT")
    receipt <- jgw2_run(args[1L], args[2L], as.numeric(args[3L]), args[4L])
    cat(jsonlite::toJSON(receipt[c("status", "accession", "allocation_id", "all_hard_invariants_pass")], auto_unbox = TRUE), "\n")
    if (receipt$all_hard_invariants_pass) 0L else 2L
  }, error = function(condition) {
    cat(jsonlite::toJSON(list(status = "STOP_SAVED_COMPATIBILITY_INPUTS_RETAINED",
      exception_class_only = class(condition)[1L], raw_conditions_persisted = FALSE), auto_unbox = TRUE), "\n"); 2L
  })
  quit(save = "no", status = code)
}

