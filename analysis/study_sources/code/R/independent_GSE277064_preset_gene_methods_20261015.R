#!/usr/bin/env Rscript
# Separate-pass native-object readback only. No producer helpers or model fits.
# Usage: PROJECT_ROOT COMPLETE_RUN_ROOT FRESH_RUN_SUBDIRECTORY
# Or --self-test, which uses invented vectors only and reads no scientific file.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

require_true <- function(value, label) {
  if (!isTRUE(value)) stop(label, call. = FALSE)
}

manual_bh <- function(probabilities) {
  require_true(is.numeric(probabilities) &&
    all(is.na(probabilities) | (is.finite(probabilities) &
      probabilities >= 0 & probabilities <= 1)), "Invalid BH probabilities")
  answer <- rep(NA_real_, length(probabilities))
  tested <- which(!is.na(probabilities))
  if (!length(tested)) return(answer)
  ordering <- order(probabilities[tested], method = "radix")
  sorted <- probabilities[tested[ordering]]
  adjusted <- pmin(1, rev(cummin(rev(sorted * length(tested) / seq_along(sorted)))))
  answer[tested[ordering]] <- adjusted
  answer
}

numeric_comparison <- function(actual, expected, probability = FALSE,
                               signed = FALSE) {
  require_true(length(actual) > 0L && length(actual) == length(expected),
    "Numeric comparison length differs or is zero")
  require_true((is.numeric(actual) || (is.logical(actual) && all(is.na(actual)))) &&
    (is.numeric(expected) || (is.logical(expected) && all(is.na(expected)))),
    "Expected numeric or wholly missing transport column")
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  missing_equal <- identical(is.na(actual), is.na(expected))
  infinite_equal <- identical(is.infinite(actual), is.infinite(expected))
  infinite_positions <- which(is.infinite(expected) & !is.na(expected))
  infinite_sign_equal <- !length(infinite_positions) ||
    identical(sign(actual[infinite_positions]), sign(expected[infinite_positions]))
  finite <- which(is.finite(actual) & is.finite(expected))
  differences <- abs(actual[finite] - expected[finite])
  bound <- 1e-10 + 1e-8 * abs(expected[finite])
  zero_equal <- TRUE
  signs_equal <- TRUE
  log_delta <- NA_real_
  thresholds_equal <- TRUE
  if (probability) {
    require_true(all(is.na(actual) | (is.finite(actual) & actual >= 0 & actual <= 1)) &&
      all(is.na(expected) | (is.finite(expected) & expected >= 0 & expected <= 1)),
      "Probability outside native domain")
    zero_equal <- identical(actual[finite] == 0, expected[finite] == 0)
    positive <- finite[actual[finite] > 0 & expected[finite] > 0]
    log_delta <- if (length(positive))
      max(abs(log10(actual[positive]) - log10(expected[positive]))) else 0
    thresholds_equal <- identical(actual[finite] <= .05, expected[finite] <= .05)
  }
  if (signed) {
    zero_equal <- identical(actual[finite] == 0, expected[finite] == 0)
    signs_equal <- identical(sign(actual[finite]), sign(expected[finite]))
  }
  list(pass = missing_equal && infinite_equal && infinite_sign_equal &&
      all(differences <= bound) && zero_equal && signs_equal && thresholds_equal &&
      (!probability || log_delta <= 1e-6),
    compared_n = length(actual), finite_pair_n = length(finite),
    expected_missing_n = sum(is.na(expected)), missing_equal = missing_equal,
    infinite_equal = infinite_equal && infinite_sign_equal,
    zero_equal = zero_equal, signs_equal = signs_equal,
    threshold_005_equal = thresholds_equal,
    maximum_absolute_difference = if (length(finite)) max(differences) else NA_real_,
    maximum_log10_probability_difference = log_delta)
}

file_identity <- function(path) {
  require_true(file.exists(path) && !dir.exists(path), "Bound file absent")
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

authenticate_file <- function(binding, label) {
  actual <- file_identity(binding$path)
  require_true(identical(actual$sha256, binding$sha256) &&
    actual$bytes == binding$bytes, paste("File binding mismatch:", label))
  actual
}

read_transport <- function(path) {
  connection <- if (endsWith(path, ".gz")) gzfile(path, open = "rt") else
    file(path, open = "rt")
  on.exit(close(connection))
  read.delim(connection, sep = "\t", quote = '"', header = TRUE,
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA",
    comment.char = "", encoding = "UTF-8")
}

locked_contrasts <- function() {
  answer <- cbind(chronic_combined = c(-1, 0, 0, .5, .5),
    D14_vs_control = c(-1, 0, 0, 1, 0),
    D28_vs_control = c(-1, 0, 0, 0, 1),
    early_combined = c(-1, .5, .5, 0, 0),
    D3_vs_control = c(-1, 1, 0, 0, 0),
    D7_vs_control = c(-1, 0, 1, 0, 0))
  rownames(answer) <- c("control", "D3", "D7", "D14", "D28")
  answer
}

family_summary <- function(result, contrast_names) {
  do.call(rbind, lapply(contrast_names, function(contrast_name) {
    rows <- result[result$contrast == contrast_name, , drop = FALSE]
    data.frame(method = unique(rows$method), contrast = contrast_name,
      input_gene_n = nrow(rows), finite_effect_n = sum(is.finite(rows$log2_effect)),
      finite_se_n = sum(is.finite(rows$standard_error)),
      finite_ci_n = sum(is.finite(rows$ci_lower) & is.finite(rows$ci_upper)),
      finite_p_n = sum(is.finite(rows$p_value)), finite_q_n = sum(is.finite(rows$bh_q)),
      native_missing_p_n = sum(!is.finite(rows$p_value)),
      eligible_inference_n = sum(rows$inference_eligible),
      q_le_005_eligible_n = sum(rows$inference_eligible & rows$bh_q <= .05, na.rm = TRUE),
      cell_type = unique(rows$cell_type))
  }))
}

run_readback <- function(project_root, run_root, output) {
  require_true(as.character(getRversion()) == "4.6.1", "R version differs")
  required_versions <- c(limma = "3.68.4", DESeq2 = "1.52.0", edgeR = "4.10.1",
    digest = "0.6.39", jsonlite = "2.0.0", data.table = "1.18.4")
  for (package in names(required_versions)) {
    require_true(as.character(utils::packageVersion(package)) == required_versions[package],
      paste("Package version differs:", package))
  }
  data.table::setDTthreads(1L)
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  run_root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
  output <- gsub("\\\\", "/", output)
  require_true(identical(normalizePath(dirname(output), winslash = "/", mustWork = TRUE), run_root) &&
    grepl("^independent_readback_[0-9]+$", basename(output)) && !file.exists(output),
    "Use a fresh independently named run subdirectory")
  completion_path <- file.path(run_root, "COMPLETE.json")
  require_true(file.exists(completion_path) && !file.exists(file.path(run_root, "FAILED.json")),
    "Readback requires completed, nonfailed run")
  completion_binding <- file_identity(completion_path)
  completion <- jsonlite::read_json(completion_path, simplifyVector = FALSE)
  require_true(identical(completion$status,
    "ADDITIONAL_GENE_METHOD_FITS_COMPLETE_INDEPENDENT_CHECK_PENDING") &&
    completion$libraries == 45L && completion$cells == 3L &&
    completion$fitted_methods_per_cell == 2L && completion$contrasts_per_method_cell == 6L &&
    completion$result_families == 36L && isTRUE(completion$all_source_pins_unchanged) &&
    !isTRUE(completion$original_edgeR_refitted) && !isTRUE(completion$preprocessing_recomputed) &&
    !isTRUE(completion$gene_universes_changed) && !isTRUE(completion$candidate_gate_eligible) &&
    !isTRUE(completion$original_primary_or_gates_changed), "Completion scope differs")
  require_true(completion$core$sha256 ==
    "e90ec660f5fa0b4f93b3d36879c9e94e5b4afec167cf7b0f04c8566b236e713b" &&
    completion$driver$sha256 ==
    "95ba4e8389dc86ad61ec508ab31d16e666319b3d65887210d87c54a5193c011c" &&
    completion$preparation$sha256 ==
    "19958d0cd1ef0dbf421a0ac90fb3456be05f912d8a5bbd66cb13bda850335daa",
    "Reviewed source identity differs")
  authenticate_file(completion$core, "core")
  authenticate_file(completion$driver, "driver")
  authenticate_file(completion$preparation, "preparation")
  preparation <- jsonlite::read_json(completion$preparation$path, simplifyVector = FALSE)
  authenticate_file(preparation$readiness, "readiness")
  readiness <- jsonlite::read_json(preparation$readiness$path, simplifyVector = FALSE)
  require_true(length(readiness$inputs) == 26L && readiness$natural_unique_library_n == 45L,
    "Fixed scientific input inventory differs")
  source_before <- lapply(names(readiness$inputs), function(key) {
    actual <- authenticate_file(readiness$inputs[[key]], key)
    data.frame(kind = "scientific_input", key = key, path = actual$path,
      bytes = actual$bytes, sha256 = actual$sha256)
  })
  expected_products <- c(as.vector(outer(c("astrocyte", "microglia", "v2a"),
    c("saved_voom", "integer_DESeq2"), paste, sep = "_")))
  product_keys <- c(paste0(expected_products, "_checkpoint"),
    paste0(expected_products, "_table"), "family_summary", "session_info")
  require_true(length(completion$products) == 14L &&
    setequal(names(completion$products), product_keys), "Complete six-fit product inventory required")
  product_before <- lapply(names(completion$products), function(key) {
    actual <- authenticate_file(completion$products[[key]], key)
    require_true(identical(dirname(actual$path), run_root), "Product outside this run")
    data.frame(kind = "run_product", key = key, path = actual$path,
      bytes = actual$bytes, sha256 = actual$sha256)
  })
  started <- jsonlite::read_json(file.path(run_root, "STARTED.json"), simplifyVector = FALSE)
  require_true(started$core$sha256 == completion$core$sha256 &&
    started$driver$sha256 == completion$driver$sha256 && started$seed == 20260817L &&
    identical(unlist(started$actual_RNG_kind, use.names = FALSE),
      c("Mersenne-Twister", "Inversion", "Rejection")), "Run start provenance differs")
  require_true(dir.create(output, recursive = FALSE), "Cannot create fresh readback directory")
  comparisons <- list()
  family_rows <- list()
  fit_rows <- list()
  scalar_n <- 0
  compare_table <- function(actual, expected, scope) {
    require_true(is.data.frame(actual) && is.data.frame(expected) &&
      identical(names(actual), names(expected)) && nrow(actual) == nrow(expected),
      paste("Table schema differs:", scope))
    for (field in names(expected)) {
      if (is.numeric(expected[[field]])) {
        result <- numeric_comparison(actual[[field]], expected[[field]],
          probability = field %in% c("p_value", "bh_q"),
          signed = field %in% c("log2_effect", "statistic", "ci_lower", "ci_upper"))
        record <- c(list(scope = scope, field = field, kind = "numeric"), result)
      } else {
        equal <- identical(actual[[field]], expected[[field]])
        record <- list(scope = scope, field = field, kind = "exact_key_or_flag",
          pass = equal, compared_n = length(expected[[field]]),
          finite_pair_n = NA_integer_, expected_missing_n = sum(is.na(expected[[field]])),
          missing_equal = identical(is.na(actual[[field]]), is.na(expected[[field]])),
          infinite_equal = NA, zero_equal = NA, signs_equal = NA,
          threshold_005_equal = NA, maximum_absolute_difference = NA_real_,
          maximum_log10_probability_difference = NA_real_)
      }
      comparisons[[length(comparisons) + 1L]] <<- as.data.frame(record)
      scalar_n <<- scalar_n + length(expected[[field]])
      require_true(isTRUE(record$pass), paste("Native readback mismatch:", scope, field))
    }
    invisible(TRUE)
  }
  compare_vector <- function(actual, expected, scope, probability = FALSE) {
    result <- numeric_comparison(actual, expected, probability = probability)
    comparisons[[length(comparisons) + 1L]] <<-
      as.data.frame(c(list(scope = scope, field = "native_relation", kind = "numeric"), result))
    scalar_n <<- scalar_n + length(expected)
    require_true(isTRUE(result$pass), paste("Native internal relation differs:", scope))
  }
  contrast_matrix <- locked_contrasts()
  levels <- rownames(contrast_matrix)
  for (cell in c("astrocyte", "microglia", "v2a")) {
    source_key <- paste0("results/GSE277064/objects/", cell, "_analysis.rds")
    original <- readRDS(readiness$inputs[[source_key]]$path)
    counts <- original$dge$counts
    ids <- rownames(counts)
    groups <- as.character(original$sample_metadata$group)
    samples <- colnames(counts)
    expected_design <- t(vapply(groups, function(group) as.numeric(levels == group), numeric(5)))
    dimnames(expected_design) <- list(samples, levels)
    require_true(is.matrix(counts) && is.numeric(counts) && ncol(counts) == 15L &&
      nrow(counts) == readiness$cell_types[[cell]]$filtered_gene_rows &&
      !anyNA(ids) && !anyDuplicated(ids) && !anyNA(samples) && !anyDuplicated(samples) &&
      all(is.finite(counts)) && all(counts >= 0) && identical(original$cell_type, cell) &&
      identical(original$domain, "natural_course"), "Source count schema differs")
    require_true(identical(samples, as.character(original$sample_metadata$source_column_or_member)) &&
      identical(samples, unlist(readiness$cell_types[[cell]]$source_columns, use.names = FALSE)) &&
      identical(as.character(original$sample_metadata$sample_id),
        unlist(readiness$cell_types[[cell]]$source_GSM_ids, use.names = FALSE)) &&
      !anyNA(groups) && identical(as.integer(table(factor(groups, levels = levels))), rep(3L, 5L)),
      "Source sample identity, order or group capacities differ")
    require_true(identical(dim(original$design), c(15L, 5L)) &&
      identical(rownames(original$design), samples) && identical(colnames(original$design), levels) &&
      identical(as.numeric(original$design), as.numeric(expected_design)) &&
      qr(original$design)$rank == 5L &&
      identical(rownames(original$contrasts), levels) &&
      identical(colnames(original$contrasts), colnames(contrast_matrix)) &&
      identical(as.numeric(original$contrasts), as.numeric(contrast_matrix)),
      "Saved design or fixed six contrast definitions differ")
    require_true(identical(dimnames(original$voom$E), dimnames(counts)) &&
      identical(dim(original$voom$weights), dim(counts)) &&
      all(is.finite(original$voom$E)) && all(is.finite(original$voom$weights)) &&
      all(original$voom$weights > 0), "Saved voom layer differs")
    for (method in c("saved_voom", "integer_DESeq2")) {
      key <- paste(cell, method, sep = "_")
      checkpoint <- readRDS(completion$products[[paste0(key, "_checkpoint")]]$path)
      require_true(identical(checkpoint$cell_type, cell) && identical(checkpoint$method, method) &&
        checkpoint$source$sha256 == readiness$inputs[[source_key]]$sha256 &&
        checkpoint$core$sha256 == completion$core$sha256 &&
        checkpoint$driver$sha256 == completion$driver$sha256 &&
        isTRUE(checkpoint$original_counts_unchanged) && !isTRUE(checkpoint$candidate_gate_eligible) &&
        is.finite(checkpoint$elapsed_seconds) && checkpoint$elapsed_seconds >= 0,
        "Saved native checkpoint provenance differs")
      stored <- checkpoint$scientific_result
      native_fit <- stored$fit
      expected <- list()
      actual_fit_type <- NA_character_
      beta_nonconverged <- cooks_above <- missing_cooks <- dispersion_outliers <- NA_integer_
      if (method == "saved_voom") {
        require_true(inherits(native_fit, "MArrayLM") &&
          identical(dim(native_fit$coefficients), c(length(ids), 6L)) &&
          identical(rownames(native_fit$coefficients), ids) &&
          identical(colnames(native_fit$coefficients), colnames(contrast_matrix)) &&
          identical(as.numeric(native_fit$design), as.numeric(expected_design)) &&
          identical(rownames(native_fit$design), samples) &&
          identical(colnames(native_fit$design), levels) &&
          identical(as.numeric(native_fit$contrasts), as.numeric(contrast_matrix)) &&
          identical(rownames(native_fit$contrasts), levels) &&
          identical(colnames(native_fit$contrasts), colnames(contrast_matrix)) &&
          all(native_fit$df.residual == 10), "Saved limma native keys/design differ")
        for (field in c("stdev.unscaled", "t", "p.value")) {
          require_true(is.matrix(native_fit[[field]]) &&
            identical(dim(native_fit[[field]]), c(length(ids), 6L)) &&
            identical(rownames(native_fit[[field]]), ids) &&
            identical(colnames(native_fit[[field]]), colnames(contrast_matrix)),
            paste("Native limma matrix keys differ:", field))
        }
        require_true(isTRUE(stored$diagnostics$checkpoint_reused) &&
          !isTRUE(stored$diagnostics$voom_recomputed) &&
          isTRUE(stored$diagnostics$robust_ebayes) && !isTRUE(stored$diagnostics$trend),
          "Saved voom invocation diagnostics differ")
        se <- sweep(native_fit$stdev.unscaled, 1L,
          sqrt(rep(native_fit$s2.post, length.out = length(ids))), `*`)
        degrees <- rep(native_fit$df.total, length.out = length(ids))
        compare_vector(as.vector(native_fit$t), as.vector(native_fit$coefficients / se),
          paste(key, "t_equals_effect_over_SE"))
        derived_p <- 2 * stats::pt(-abs(native_fit$t), df = degrees)
        compare_vector(as.vector(native_fit$p.value), as.vector(derived_p),
          paste(key, "two_sided_native_t_P"), probability = TRUE)
        for (index in seq_len(6L)) {
          effect <- unname(native_fit$coefficients[, index])
          standard_error <- unname(se[, index])
          p <- unname(native_fit$p.value[, index])
          half <- stats::qt(.975, df = degrees) * standard_error
          expected[[index]] <- data.frame(gene_id = ids, method = "saved_voom_limma",
            contrast = colnames(contrast_matrix)[index], log2_effect = effect,
            standard_error = standard_error, statistic = unname(native_fit$t[, index]),
            degrees_freedom = unname(degrees), ci_lower = effect - half,
            ci_upper = effect + half, ci_kind = "nominal_95_percent_moderated_t",
            p_value = p, bh_q = manual_bh(p),
            inference_eligible = is.finite(effect) & is.finite(standard_error) &
              standard_error > 0 & is.finite(degrees) & degrees > 0 & is.finite(p),
            cell_type = cell)
        }
      } else {
        require_true(inherits(native_fit, "DESeqDataSet"), "Expected saved DESeqDataSet")
        rounded <- round(counts)
        storage.mode(rounded) <- "integer"
        native_counts <- SummarizedExperiment::assay(native_fit, "counts")
        native_model <- attr(native_fit, "modelMatrix")
        native_metadata <- as.data.frame(SummarizedExperiment::colData(native_fit))
        coefficient_names <- DESeq2::resultsNames(native_fit)
        require_true(identical(native_counts, rounded) &&
          identical(rownames(native_fit), ids) && identical(colnames(native_fit), samples) &&
          identical(as.character(native_metadata$group), groups) &&
          identical(rownames(native_metadata), samples) &&
          identical(dim(native_model), c(15L, 5L)) &&
          identical(as.numeric(native_model), as.numeric(expected_design)) &&
          identical(sub("^group_?", "", coefficient_names), levels) &&
          identical(colnames(stats::coef(native_fit)), coefficient_names) &&
          is.null(DESeq2::normalizationFactors(native_fit)) &&
          !"replaceCounts" %in% SummarizedExperiment::assayNames(native_fit),
          "Native DESeq counts, design, coefficient mapping or replacement differs")
        compare_table(stored$coefficient_mapping,
          data.frame(saved_group = levels, native_coefficient = coefficient_names),
          paste(key, "coefficient_mapping"))
        gene_metadata <- as.data.frame(S4Vectors::mcols(native_fit))
        require_true(all(c("betaConv", "allZero", "maxCooks", "dispOutlier") %in%
          names(gene_metadata)), "Missing native DESeq diagnostics")
        factors <- DESeq2::sizeFactors(native_fit)
        require_true(all(is.finite(factors)) && all(factors > 0), "Invalid native size factors")
        compare_vector(unlist(stored$diagnostics$size_factors[samples], use.names = FALSE),
          unname(factors), paste(key, "native_size_factors"))
        cutoff <- stats::qf(.99, 5L, 10L)
        actual_fit_type <- attr(DESeq2::dispersionFunction(native_fit), "fitType")
        beta_nonconverged <- sum(!gene_metadata$betaConv, na.rm = TRUE)
        cooks_above <- sum(!is.na(gene_metadata$maxCooks) & gene_metadata$maxCooks > cutoff)
        missing_cooks <- sum(is.na(gene_metadata$maxCooks))
        dispersion_outliers <- sum(gene_metadata$dispOutlier, na.rm = TRUE)
        require_true(actual_fit_type %in% c("parametric", "local") &&
          identical(stored$diagnostics$requested_fit_type, "parametric") &&
          identical(stored$diagnostics$actual_fit_type, actual_fit_type) &&
          isTRUE(stored$diagnostics$native_fitted_model_matrix_verified) &&
          identical(stored$diagnostics$rounding_rule, "R round ties to even in RAM; source remains fractional") &&
          stored$diagnostics$fractional_entry_n == sum(counts != rounded) &&
          stored$diagnostics$maximum_rounding_change == max(abs(counts - rounded)) &&
          stored$diagnostics$rounded_all_zero_n == sum(rowSums(rounded) == 0) &&
          stored$diagnostics$beta_nonconvergence_n == beta_nonconverged &&
          !isTRUE(stored$diagnostics$replacement_performed) &&
          !isTRUE(stored$diagnostics$independent_filtering) &&
          !isTRUE(stored$diagnostics$gene_universe_changed) &&
          identical(unname(gene_metadata$allZero), unname(rowSums(rounded) == 0)),
          "Native DESeq fitting/rounding diagnostics differ")
        compare_vector(stored$diagnostics$cook_cutoff, cutoff, paste(key, "standard_Cooks_cutoff"))
        for (index in seq_len(6L)) {
          # This accessor recomputes a contrast from saved fit parameters, not a fit.
          native <- as.data.frame(DESeq2::results(native_fit,
            contrast = unname(contrast_matrix[, index]), alpha = .05,
            pAdjustMethod = "BH", independentFiltering = FALSE,
            cooksCutoff = TRUE, parallel = FALSE))
          require_true(identical(rownames(native), ids), "Native DESeq result order differs")
          p <- native$pvalue
          compare_vector(native$padj, manual_bh(p), paste(key, index, "native_BH"), probability = TRUE)
          finite_p <- which(is.finite(p))
          if (length(finite_p)) {
            compare_vector(p[finite_p], 2 * stats::pnorm(-abs(native$stat[finite_p])),
              paste(key, index, "finite_native_Wald_P"), probability = TRUE)
          }
          above <- !is.na(gene_metadata$maxCooks) & gene_metadata$maxCooks > cutoff
          missing_p <- !is.finite(p)
          eligible <- !is.na(gene_metadata$betaConv) & gene_metadata$betaConv &
            is.finite(p) & is.finite(native$padj) & is.finite(native$log2FoldChange) &
            is.finite(native$lfcSE) & native$lfcSE > 0
          half <- stats::qnorm(.975) * native$lfcSE
          expected[[index]] <- data.frame(gene_id = ids, method = "DESeq2_integer_approximation",
            contrast = colnames(contrast_matrix)[index], base_mean = native$baseMean,
            log2_effect = native$log2FoldChange, standard_error = native$lfcSE,
            statistic = native$stat, degrees_freedom = NA_real_,
            ci_lower = native$log2FoldChange - half, ci_upper = native$log2FoldChange + half,
            ci_kind = "nominal_95_percent_unshrunk_Wald", p_value = p, bh_q = native$padj,
            beta_converged = gene_metadata$betaConv, all_zero = gene_metadata$allZero,
            max_cooks = gene_metadata$maxCooks, standard_cook_cutoff = cutoff,
            above_standard_cook_threshold = above, cook_value_unknown = is.na(gene_metadata$maxCooks),
            dispersion_outlier = gene_metadata$dispOutlier,
            native_missing_p = missing_p, native_missing_q = !is.finite(native$padj),
            missing_p_all_zero_compatible = missing_p & gene_metadata$allZero,
            missing_p_cooks_compatible = missing_p & above & !gene_metadata$allZero,
            missing_p_convergence_compatible = missing_p & !is.na(gene_metadata$betaConv) &
              !gene_metadata$betaConv,
            missing_p_other_or_unresolved = missing_p & !gene_metadata$allZero & !above &
              (is.na(gene_metadata$betaConv) | gene_metadata$betaConv),
            inference_eligible = eligible, cell_type = cell)
        }
      }
      expected_table <- do.call(rbind, expected)
      require_true(nrow(expected_table) == length(ids) * 6L &&
        !anyDuplicated(expected_table[c("gene_id", "method", "contrast")]),
        "Incomplete or duplicate expected native families")
      compare_table(stored$results, expected_table, paste(key, "RDS_native_readback"))
      transport <- read_transport(completion$products[[paste0(key, "_table")]]$path)
      compare_table(transport, expected_table, paste(key, "TSV_transport_readback"))
      expected_summary <- family_summary(expected_table, colnames(contrast_matrix))
      compare_table(stored$summary, expected_summary, paste(key, "saved_family_summary"))
      family_rows[[key]] <- expected_summary
      fit_rows[[key]] <- data.frame(cell_type = cell, method = method,
        retained_gene_n = length(ids), library_observation_n = 15L, group_n = 5L,
        libraries_per_group = 3L, design_rank = 5L, residual_df = 10L,
        contrast_n = 6L, gene_contrast_rows = nrow(expected_table),
        actual_dispersion_fit = actual_fit_type, beta_nonconverged_n = beta_nonconverged,
        above_standard_Cooks_n = cooks_above, unknown_maxCooks_n = missing_cooks,
        dispersion_outlier_n = dispersion_outliers,
        finite_p_n = sum(is.finite(expected_table$p_value)),
        missing_p_n = sum(is.na(expected_table$p_value)),
        eligible_n = sum(expected_table$inference_eligible),
        captured_warning_n = length(checkpoint$warning_messages),
        captured_message_n = length(checkpoint$messages), elapsed_seconds = checkpoint$elapsed_seconds,
        source_counts_not_changed = TRUE, no_model_refit_by_checker = TRUE)
      rm(checkpoint, stored, native_fit, expected, expected_table, transport)
      invisible(gc(FALSE))
      cat("READBACK_COMPLETE", key, "six complete families\n")
    }
    authenticate_file(readiness$inputs[[source_key]], source_key)
    rm(original, counts)
    invisible(gc(FALSE))
  }
  all_families <- do.call(rbind, family_rows)
  rownames(all_families) <- NULL
  require_true(nrow(all_families) == 36L, "Exactly 36 complete families required")
  compare_table(read_transport(completion$products$family_summary$path), all_families,
    "complete_36_row_family_summary")
  for (key in names(readiness$inputs)) authenticate_file(readiness$inputs[[key]], key)
  for (key in names(completion$products)) authenticate_file(completion$products[[key]], key)
  authenticate_file(completion$core, "core_closure")
  authenticate_file(completion$driver, "driver_closure")
  authenticate_file(completion$preparation, "preparation_closure")
  authenticate_file(preparation$readiness, "readiness_closure")
  authenticate_file(completion_binding, "completion_closure")
  comparisons <- do.call(rbind, comparisons)
  fits <- do.call(rbind, fit_rows)
  file_bindings <- do.call(rbind, c(source_before, product_before))
  tables <- list(FIELD_COMPARISONS = comparisons, FAMILY_COUNTS = all_families,
    FIT_DIAGNOSTICS = fits, FILE_BINDINGS = file_bindings)
  output_bindings <- lapply(names(tables), function(name) {
    path <- file.path(output, paste0(name, ".tsv"))
    write.table(tables[[name]], path, sep = "\t", quote = TRUE, row.names = FALSE,
      na = "NA", fileEncoding = "UTF-8")
    file_identity(path)
  })
  names(output_bindings) <- names(tables)
  checker_path <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
  receipt <- list(status = "PASS_SCOPED_NATIVE_RESULT_READBACK_NOT_WHOLE_GOAL_ACCEPTANCE",
    validation_status = "PASS", completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    checker = file_identity(checker_path), completion = completion_binding,
    core = completion$core, driver = completion$driver,
    preparation = completion$preparation, readiness = preparation$readiness,
    model_refit_n = 0L, native_model_objects_read_n = 6L, original_cell_objects_read_n = 3L,
    method_cell_fits = 6L, contrast_families = 36L,
    gene_contrast_rows = sum(fits$gene_contrast_rows),
    numeric_and_exact_field_comparisons = nrow(comparisons), compared_scalar_n = scalar_n,
    all_field_comparisons_pass = all(comparisons$pass),
    scientific_input_pins_unchanged_n = 26L, protected_master_gate_pin_unchanged = TRUE,
    before_library_observations = 45L, after_library_observations = 45L,
    original_records_removed = 0L, original_files_modified = 0L,
    RDS_native_and_TSV_transport_separately_compared = TRUE,
    arithmetic_tolerance = list(absolute = 1e-10, relative_to_expected = 1e-8,
      positive_probability_log10 = 1e-6, missing_pattern = "exact",
      probability_zero_and_005_threshold = "exact", signed_zero_and_direction = "exact",
      all_missing_TSV_numeric_type = "logical-all-NA accepted as missing numeric transport"),
    BH = "Manual BH over native finite-P rows within each cell/method/contrast; missing rows retained, no denominator or eligibility substitution",
    scientific_status = "EXPLORATORY_SAME_DATA_SENSITIVITY_WITH_UNRESOLVED_PROVENANCE_LIMITATIONS",
    limitations = c("No independent biological cohort, animal/pooling or batch authentication",
      "No recomputation of regression or moderation; native saved-fit readback only",
      "Saved voom E/weights use is source-bound invocation provenance, not independent regression reconstruction",
      "RSEM rounding is an integer approximation, not uncertainty propagation or length-offset restoration",
      "Native missingness/nonconvergence/Cooks diagnostics preserved, not repaired",
      "No validation of old edgeR uncertainty, S14/S16, S1/S3 candidate gates or mechanistic/efficacy claims",
      "No cold environment restoration, author review, independent human certification or whole-goal completion"),
    candidate_gate_updated = FALSE, original_edgeR_or_gates_changed = FALSE,
    whole_goal_complete = FALSE, outputs = output_bindings,
    actual_R = R.version.string,
    actual_packages = as.list(vapply(names(required_versions), function(package)
      as.character(utils::packageVersion(package)), character(1))))
  jsonlite::write_json(receipt, file.path(output, "INDEPENDENT_READBACK_RECEIPT.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
  cat("PASS_SCOPED_NATIVE_RESULT_READBACK; six native objects, 36 families; zero model fits\n")
  invisible(receipt)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (identical(arguments, "--self-test")) {
  # Invented arithmetic/missingness tests; not native model validation.
  checks <- c(
    manual_BH = identical(manual_bh(c(.01, NA_real_, .04, .03)), c(.03, NA_real_, .04, .04)),
    all_missing_BH = identical(manual_bh(c(NA_real_, NA_real_)), c(NA_real_, NA_real_)),
    missing_mask_rejected = !numeric_comparison(c(1, NA_real_), c(NA_real_, 1))$pass,
    signed_cross_zero_rejected = !numeric_comparison(1e-13, -1e-13, signed = TRUE)$pass,
    probability_tiny_ratio_rejected = !numeric_comparison(1e-100, 1e-90, probability = TRUE)$pass,
    probability_zero_rejected = !numeric_comparison(0, 1e-100, probability = TRUE)$pass,
    probability_threshold_rejected = !numeric_comparison(.05000000001, .05, probability = TRUE)$pass,
    transport_roundoff_allowed = numeric_comparison(1.23456789012346, 1.234567890123456)$pass,
    all_missing_transport = numeric_comparison(c(NA, NA), c(NA_real_, NA_real_))$pass,
    identical_inf = numeric_comparison(c(Inf, -Inf), c(Inf, -Inf))$pass,
    opposite_inf_rejected = !numeric_comparison(Inf, -Inf)$pass)
  require_true(all(checks), "Independent invented arithmetic test failed")
  print(checks)
  cat("PASS_INVENTED_READBACK_ARITHMETIC; no scientific inputs read and no fits\n")
} else {
  require_true(length(arguments) == 3L, "PROJECT_ROOT COMPLETE_RUN_ROOT FRESH_RUN_SUBDIRECTORY required")
  run_readback(arguments[1L], arguments[2L], arguments[3L])
}
