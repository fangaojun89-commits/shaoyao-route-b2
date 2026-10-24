# Functions only: this file does not read scientific inputs or start a fit.
# Limited same-data exploratory sensitivities; primary edgeR is untouched.
gene_method_contract <- function() {
  group_levels <- c("control", "D3", "D7", "D14", "D28")
  contrast_matrix <- matrix(c(
    -1, 0, 0, .5, .5,
    -1, 0, 0, 1, 0,
    -1, 0, 0, 0, 1,
    -1, .5, .5, 0, 0,
    -1, 1, 0, 0, 0,
    -1, 0, 1, 0, 0
  ), nrow = 5L, dimnames = list(group_levels, c(
    "chronic_combined", "D14_vs_control", "D28_vs_control",
    "early_combined", "D3_vs_control", "D7_vs_control"
  )))
  list(group_levels = group_levels, contrasts = contrast_matrix, alpha = .05)
}

assert_gene_method <- function(condition, label) {
  if (!isTRUE(condition)) stop(label, call. = FALSE)
}

validate_gene_method_inputs <- function(counts, groups, design, contrasts,
                                        voom_checkpoint = NULL) {
  contract <- gene_method_contract()
  assert_gene_method(is.matrix(counts) && is.numeric(counts), "Counts must be a numeric matrix")
  assert_gene_method(nrow(counts) > 0L && ncol(counts) == 15L, "Expected 15 selected libraries")
  assert_gene_method(all(is.finite(counts)) && all(counts >= 0), "Invalid unnormalized counts")
  assert_gene_method(!is.null(rownames(counts)) && !is.null(colnames(counts)) &&
                       !anyNA(rownames(counts)) && !anyNA(colnames(counts)) &&
                       all(nzchar(rownames(counts))) && all(nzchar(colnames(counts))) &&
                       !anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)),
                     "Invalid count keys")
  group_factor <- factor(as.character(groups), levels = contract$group_levels)
  assert_gene_method(!anyNA(group_factor) && length(group_factor) == ncol(counts) &&
                       identical(as.integer(table(group_factor)), rep(3L, 5L)),
                     "Locked group population differs")
  rebuilt <- model.matrix(~ 0 + group_factor)
  assert_gene_method(is.matrix(design) && identical(dim(design), dim(rebuilt)) &&
                       identical(rownames(design), colnames(counts)) &&
                       identical(colnames(design), contract$group_levels) &&
                       identical(as.numeric(design), as.numeric(rebuilt)) && qr(design)$rank == 5L,
                     "Locked design or sample order differs")
  assert_gene_method(is.matrix(contrasts) && identical(dim(contrasts), dim(contract$contrasts)) &&
                       identical(rownames(contrasts), rownames(contract$contrasts)) &&
                       identical(colnames(contrasts), colnames(contract$contrasts)) &&
                       identical(as.numeric(contrasts), as.numeric(contract$contrasts)),
                     "Locked contrasts differ")
  assert_gene_method(all(colSums(counts) > 0), "Zero-count library")
  if (!is.null(voom_checkpoint)) {
    assert_gene_method(identical(dimnames(voom_checkpoint$E), dimnames(counts)) &&
                         identical(dim(voom_checkpoint$weights), dim(counts)) &&
                         all(is.finite(voom_checkpoint$E)) &&
                         all(is.finite(voom_checkpoint$weights)) &&
                         all(voom_checkpoint$weights > 0), "Invalid saved voom checkpoint")
    if (!is.null(dimnames(voom_checkpoint$weights))) {
      assert_gene_method(identical(dimnames(voom_checkpoint$weights), dimnames(counts)),
                         "Named saved precision-weight keys differ")
    }
  }
  invisible(group_factor)
}

gene_family_summary <- function(result_table) {
  keys <- unique(result_table[, c("method", "contrast"), drop = FALSE])
  summaries <- lapply(seq_len(nrow(keys)), function(index) {
    selected <- result_table[result_table$method == keys$method[index] &
                               result_table$contrast == keys$contrast[index], , drop = FALSE]
    data.frame(
      method = keys$method[index], contrast = keys$contrast[index],
      input_gene_n = nrow(selected), finite_effect_n = sum(is.finite(selected$log2_effect)),
      finite_se_n = sum(is.finite(selected$standard_error)),
      finite_ci_n = sum(is.finite(selected$ci_lower) & is.finite(selected$ci_upper)),
      finite_p_n = sum(is.finite(selected$p_value)), finite_q_n = sum(is.finite(selected$bh_q)),
      native_missing_p_n = sum(!is.finite(selected$p_value)),
      eligible_inference_n = sum(selected$inference_eligible),
      q_le_005_eligible_n = sum(selected$inference_eligible & selected$bh_q <= .05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summaries)
}

fit_saved_voom_genes <- function(counts, groups, design, contrasts, voom_checkpoint) {
  validate_gene_method_inputs(counts, groups, design, contrasts, voom_checkpoint)
  fit <- limma::lmFit(voom_checkpoint$E, design = design, weights = voom_checkpoint$weights)
  fit <- limma::contrasts.fit(fit, contrasts = contrasts)
  fit <- limma::eBayes(fit, robust = TRUE, trend = FALSE)
  assert_gene_method(identical(dim(fit$coefficients), c(nrow(counts), ncol(contrasts))) &&
                       identical(rownames(fit$coefficients), rownames(counts)) &&
                       identical(colnames(fit$coefficients), colnames(contrasts)),
                     "Voom fit result keys differ")
  degrees_freedom <- rep(fit$df.total, length.out = nrow(counts))
  posterior_variance <- rep(fit$s2.post, length.out = nrow(counts))
  records <- lapply(seq_len(ncol(contrasts)), function(index) {
    effect <- fit$coefficients[, index]
    standard_error <- fit$stdev.unscaled[, index] * sqrt(posterior_variance)
    half_width <- qt(.975, df = degrees_freedom) * standard_error
    p_value <- fit$p.value[, index]
    valid <- is.finite(effect) & is.finite(standard_error) & standard_error > 0 &
      is.finite(degrees_freedom) & degrees_freedom > 0 & is.finite(p_value)
    data.frame(
      gene_id = rownames(counts), method = "saved_voom_limma", contrast = colnames(contrasts)[index],
      log2_effect = unname(effect), standard_error = unname(standard_error),
      statistic = unname(fit$t[, index]), degrees_freedom = unname(degrees_freedom),
      ci_lower = unname(effect - half_width), ci_upper = unname(effect + half_width),
      ci_kind = "nominal_95_percent_moderated_t", p_value = unname(p_value),
      bh_q = p.adjust(p_value, method = "BH"), inference_eligible = unname(valid),
      stringsAsFactors = FALSE
    )
  })
  result_table <- do.call(rbind, records)
  list(results = result_table, summary = gene_family_summary(result_table), fit = fit,
       diagnostics = list(checkpoint_reused = TRUE, voom_recomputed = FALSE,
                          robust_ebayes = TRUE, trend = FALSE, bh_family = "cell_by_contrast_by_method"))
}

assemble_deseq_gene_rows <- function(native, gene_ids, contrast_name, beta_converged,
                                      all_zero, max_cooks, cook_cutoff, dispersion_outlier) {
  assert_gene_method(identical(rownames(native), gene_ids), "DESeq2 result row order differs")
  assert_gene_method(all(vapply(list(beta_converged, all_zero, max_cooks, dispersion_outlier),
                               length, integer(1)) == length(gene_ids)), "DESeq2 diagnostic length differs")
  native_p <- native$pvalue
  native_q <- native$padj
  above_cook <- !is.na(max_cooks) & max_cooks > cook_cutoff
  missing_p <- !is.finite(native_p)
  eligible <- !is.na(beta_converged) & beta_converged & is.finite(native_p) &
    is.finite(native_q) & is.finite(native$log2FoldChange) & is.finite(native$lfcSE) & native$lfcSE > 0
  half_width <- qnorm(.975) * native$lfcSE
  # Preserve native numbers, including NA; reasons are compatible diagnostics,
  # not assertions that a high Cook's value always caused the missing result.
  data.frame(
    gene_id = gene_ids, method = "DESeq2_integer_approximation", contrast = contrast_name,
    base_mean = native$baseMean, log2_effect = native$log2FoldChange,
    standard_error = native$lfcSE, statistic = native$stat, degrees_freedom = NA_real_,
    ci_lower = native$log2FoldChange - half_width, ci_upper = native$log2FoldChange + half_width,
    ci_kind = "nominal_95_percent_unshrunk_Wald", p_value = native_p, bh_q = native_q,
    beta_converged = beta_converged, all_zero = all_zero, max_cooks = max_cooks,
    standard_cook_cutoff = cook_cutoff, above_standard_cook_threshold = above_cook,
    cook_value_unknown = is.na(max_cooks),
    dispersion_outlier = dispersion_outlier, native_missing_p = missing_p,
    native_missing_q = !is.finite(native_q),
    missing_p_all_zero_compatible = missing_p & all_zero,
    missing_p_cooks_compatible = missing_p & above_cook & !all_zero,
    missing_p_convergence_compatible = missing_p & !is.na(beta_converged) & !beta_converged,
    missing_p_other_or_unresolved = missing_p & !all_zero & !above_cook &
      (is.na(beta_converged) | beta_converged),
    inference_eligible = eligible, stringsAsFactors = FALSE
  )
}

fit_integer_deseq_genes <- function(counts, groups, design, contrasts) {
  group_factor <- validate_gene_method_inputs(counts, groups, design, contrasts)
  rounded <- round(counts)
  assert_gene_method(all(rounded <= .Machine$integer.max), "DESeq2 integer overflow")
  assert_gene_method(!any(rowSums(rounded) == 0), "Rounding created an all-zero selected gene")
  storage.mode(rounded) <- "integer"
  sample_metadata <- data.frame(group = group_factor, row.names = colnames(counts))
  reconstructed <- model.matrix(~ 0 + group, data = sample_metadata)
  assert_gene_method(identical(dim(reconstructed), dim(design)) &&
                       identical(as.numeric(reconstructed), as.numeric(design)) &&
                       identical(rownames(reconstructed), rownames(design)), "DESeq2 design differs")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = rounded, colData = sample_metadata,
                                       design = ~ 0 + group)
  dds <- DESeq2::DESeq(dds, test = "Wald", betaPrior = FALSE, fitType = "parametric",
                       sfType = "ratio", minReplicatesForReplace = Inf, parallel = FALSE,
                       quiet = FALSE)
  coefficient_names <- DESeq2::resultsNames(dds)
  coefficient_semantics <- sub("^group_?", "", coefficient_names)
  assert_gene_method(identical(coefficient_semantics, gene_method_contract()$group_levels),
                     "DESeq2 coefficient semantics/order differ")
  native_model_matrix <- attr(dds, "modelMatrix")
  assert_gene_method(is.matrix(native_model_matrix) &&
                       identical(dim(native_model_matrix), dim(design)) &&
                       identical(as.numeric(native_model_matrix), as.numeric(design)),
                     "Native DESeq2 fitted matrix differs")
  assert_gene_method(identical(colnames(stats::coef(dds)), coefficient_names),
                     "DESeq2 fitted coefficient order differs")
  gene_metadata <- as.data.frame(S4Vectors::mcols(dds))
  required <- c("betaConv", "allZero", "maxCooks", "dispOutlier")
  assert_gene_method(all(required %in% names(gene_metadata)), "DESeq2 diagnostic fields missing")
  size_factors <- DESeq2::sizeFactors(dds)
  assert_gene_method(all(is.finite(size_factors)) && all(size_factors > 0), "Invalid DESeq2 size factors")
  actual_fit_type <- attr(DESeq2::dispersionFunction(dds), "fitType")
  assert_gene_method(actual_fit_type %in% c("parametric", "local"), "Unexpected dispersion-fit substitution")
  cook_cutoff <- qf(.99, ncol(design), nrow(design) - ncol(design))
  records <- lapply(seq_len(ncol(contrasts)), function(index) {
    native <- as.data.frame(DESeq2::results(
      dds, contrast = unname(contrasts[, index]), alpha = .05, pAdjustMethod = "BH",
      independentFiltering = FALSE, cooksCutoff = TRUE, parallel = FALSE
    ))
    assemble_deseq_gene_rows(native, rownames(counts), colnames(contrasts)[index],
                             gene_metadata$betaConv, gene_metadata$allZero,
                             gene_metadata$maxCooks, cook_cutoff, gene_metadata$dispOutlier)
  })
  result_table <- do.call(rbind, records)
  assert_gene_method(!"replaceCounts" %in% SummarizedExperiment::assayNames(dds),
                     "Unexpected DESeq2 outlier replacement")
  list(
    results = result_table, summary = gene_family_summary(result_table), fit = dds,
    coefficient_mapping = data.frame(saved_group = colnames(design),
                                     native_coefficient = coefficient_names),
    diagnostics = list(
      requested_fit_type = "parametric", actual_fit_type = actual_fit_type,
      native_fitted_model_matrix_verified = TRUE,
      rounding_rule = "R round ties to even in RAM; source remains fractional",
      fractional_entry_n = sum(counts != rounded),
      maximum_rounding_change = max(abs(counts - rounded)),
      rounded_all_zero_n = sum(rowSums(rounded) == 0),
      size_factors = setNames(as.list(size_factors), colnames(counts)),
      replacement_performed = FALSE, independent_filtering = FALSE,
      cook_cutoff = cook_cutoff, beta_nonconvergence_n = sum(!gene_metadata$betaConv, na.rm = TRUE),
      gene_universe_changed = FALSE, bh_family = "cell_by_contrast_by_method"
    )
  )
}
