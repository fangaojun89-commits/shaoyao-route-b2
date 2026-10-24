# Functions only. Caller-held NEW matched checkpoint, never the old checked_payload.
# No file IO, fit dispatcher, source import, RNG call or acceptance-file creation.
# This bounded component does not authenticate filesystem/producer custody and
# does not replay CAMERA, mroast, fgsea, mapping or activity-score construction.

g205n_need <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
}

g205n_equal <- function(observed, expected, label) {
  g205n_need(is.numeric(observed) && is.numeric(expected) &&
    length(observed) == length(expected), paste0(label, ": numeric domain"))
  g205n_need(identical(is.na(observed), is.na(expected)) &&
    identical(is.nan(observed), is.nan(expected)), paste0(label, ": missingness"))
  finite <- is.finite(expected)
  g205n_need(identical(is.finite(observed), finite), paste0(label, ": finiteness"))
  infinite <- !finite & !is.na(expected)
  g205n_need(all(observed[infinite] == expected[infinite]), paste0(label, ": infinity sign"))
  # Fixed tolerance, never widened from observed discrepancies.
  g205n_need(all(abs(observed[finite] - expected[finite]) <=
    1e-12 + 1e-10 * abs(expected[finite])), paste0(label, ": values"))
  invisible(TRUE)
}

g205n_expand <- function(value, n, label) {
  g205n_need(is.numeric(value) && length(value) %in% c(1L, n), label)
  if (length(value) == 1L) rep(value, n) else value
}

g205n_probability <- function(value, label) {
  g205n_need(is.numeric(value) &&
    all(is.na(value) | (is.finite(value) & value >= 0 & value <= 1)), label)
  invisible(TRUE)
}

g205n_rows <- function(table, contrast, ids, id_column, collection = NULL) {
  table <- as.data.frame(table)
  g205n_need(!anyDuplicated(names(table)) &&
    all(c(id_column, "contrast") %in% names(table)), "Table keys missing")
  g205n_need(!anyNA(table$contrast) &&
    all(table$contrast %in% c("local_disease", "intervention")), "Unexpected contrast")
  selected <- table$contrast == contrast
  if (!is.null(collection)) {
    g205n_need("collection" %in% names(table) && !anyNA(table$collection), "Collection missing")
    selected <- selected & table$collection == collection
  }
  answer <- table[selected, , drop = FALSE]
  g205n_need(identical(as.character(answer[[id_column]]), ids) &&
    !anyNA(answer[[id_column]]) && !anyDuplicated(answer[[id_column]]), "Native ordered keys differ")
  answer
}

g205n_limma <- function(table, fit, ids, effect_column, id_column,
                        collection = NULL) {
  contrasts <- c("local_disease", "intervention")
  if (!is.null(collection)) g205n_need(!anyNA(table$collection) &&
    all(table$collection %in% c("PROGENy", "DoRothEA_AB")), "Unexpected activity collection")
  for (field in c("coefficients", "stdev.unscaled", "t", "p.value")) {
    value <- fit[[field]]
    g205n_need(is.matrix(value) && identical(rownames(value), ids) &&
      identical(colnames(value), contrasts), paste0("Saved limma ", field, " keys"))
  }
  posterior <- g205n_expand(fit$s2.post, length(ids), "Saved posterior variance domain")
  degrees <- g205n_expand(fit$df.total, length(ids), "Saved moderated df domain")
  for (j in seq_along(contrasts)) {
    row <- g205n_rows(table, contrasts[[j]], ids, id_column, collection)
    method <- if (id_column == "gene_id") "limma_voom_robust_eBayes" else "decoupleR_ULM_then_limma"
    g205n_need("method" %in% names(row) && !anyNA(row$method) && all(row$method == method), "Limma method label differs")
    effect <- fit$coefficients[, j]
    error <- fit$stdev.unscaled[, j] * sqrt(posterior)
    critical <- stats::qt(.975, degrees)
    expected <- list(se = error, ci_low = effect - critical * error,
      ci_high = effect + critical * error, t = fit$t[, j], df_total = degrees,
      p_value = fit$p.value[, j])
    g205n_probability(fit$p.value[, j], "Saved limma P domain")
    expected[[effect_column]] <- effect
    q_column <- if (id_column == "gene_id") "q_value" else "q_within_collection_contrast"
    expected[[q_column]] <- stats::p.adjust(fit$p.value[, j], "BH")
    for (field in names(expected)) g205n_equal(row[[field]], unname(expected[[field]]), field)
    g205n_equal(unname(fit$t[, j]), unname(effect / error), "Saved limma effect/SE and t")
    g205n_equal(unname(fit$p.value[, j]), unname(2 * stats::pt(-abs(fit$t[, j]), degrees)),
      "Saved limma moderated two-sided P")
  }
  invisible(TRUE)
}

g205n_edger <- function(table, tests, ids) {
  g205n_need(identical(names(tests), c("local_disease", "intervention")), "Saved QL contrast keys")
  for (contrast in names(tests)) {
    test <- tests[[contrast]]
    native <- as.data.frame(test$table)
    g205n_need(identical(rownames(native), ids), "Saved QL gene keys")
    row <- g205n_rows(table, contrast, ids, "gene_id")
    g205n_need("method" %in% names(row) && !anyNA(row$method) &&
      all(row$method == "edgeR_QL_robust"), "QL method label differs")
    for (field in c("logFC", "logCPM", "F", "PValue")) {
      g205n_equal(row[[field]], native[[field]], paste0("QL ", field))
    }
    q <- stats::p.adjust(native$PValue, "BH")
    g205n_probability(native$PValue, "Saved QL P domain")
    g205n_equal(row$FDR, q, "QL FDR")
    g205n_equal(row$q_value, q, "QL q")
    degrees <- if (is.null(test$df.total)) rep(NA_real_, length(ids)) else
      g205n_expand(test$df.total, length(ids), "Saved QL df domain")
    g205n_equal(row$df_total, unname(degrees), "QL df")
    for (field in c("se", "ci_low", "ci_high", "se_df")) {
      g205n_equal(row[[field]], rep(NA_real_, nrow(row)), "Unvalidated QL interval introduced")
    }
    g205n_need(is.character(row$se_status) && length(row$se_status) == nrow(row) &&
      all(row$se_status == "NOT_REPORTED_NO_VALIDATED_NATIVE_QL_INTERVAL") &&
      !anyNA(row$se_status), "QL interval status differs")
  }
  invisible(TRUE)
}

g205n_deseq_result <- function(table, native, contrast) {
  native <- as.data.frame(native)
  g205n_probability(native$pvalue, "Saved DESeq P domain")
  g205n_probability(native$padj, "Saved DESeq q domain")
  row <- g205n_rows(table, contrast, rownames(native), "gene_id")
  g205n_need("method" %in% names(row) && !anyNA(row$method) &&
    all(row$method == "DESeq2_Wald_sensitivity"), "DESeq method label differs")
  for (field in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
    g205n_equal(row[[field]], native[[field]], paste0("DESeq ", field))
  }
  g205n_equal(row$p_value, native$pvalue, "DESeq P alias")
  g205n_equal(row$q_value, native$padj, "DESeq q alias")
  g205n_equal(row$ci_low, native$log2FoldChange - stats::qnorm(.975) * native$lfcSE, "DESeq lower")
  g205n_equal(row$ci_high, native$log2FoldChange + stats::qnorm(.975) * native$lfcSE, "DESeq upper")
  invisible(TRUE)
}

g205n_unsigned_q <- function(table, kind) {
  table <- as.data.frame(table)
  g205n_need(all(c("collection", "contrast", "program_id", "p_value",
    "q_within_collection_contrast") %in% names(table)), "Unsigned fields missing")
  g205n_need(!anyNA(table[c("collection", "contrast", "program_id")]) &&
    all(table$collection %in% c("Reactome", "SynGO_sensitivity")) &&
    all(table$contrast %in% c("local_disease", "intervention")) &&
    !anyDuplicated(table[c("collection", "contrast", "program_id")]), "Unsigned family keys differ")
  p_column <- if (kind == "fgsea") "pval" else "PValue"
  for (collection in c("Reactome", "SynGO_sensitivity")) {
    for (contrast in c("local_disease", "intervention")) {
      row <- table[table$collection == collection & table$contrast == contrast, , drop = FALSE]
      p <- row[[p_column]]
      g205n_need(is.numeric(p) && all(is.na(p) | (is.finite(p) & p >= 0 & p <= 1)), "Unsigned P domain")
      g205n_equal(row$p_value, p, "Unsigned native P alias")
      g205n_equal(row$q_within_collection_contrast, stats::p.adjust(p, "BH"), "Unsigned family BH")
      if (kind == "fgsea") g205n_equal(row$q_within_collection_contrast, row$padj, "fgsea native q alias")
    }
  }
  # Native mroast FDR/FDR.Mixed are NOT reconstructed as ordinary raw-P BH.
  invisible(TRUE)
}

g205n_readback <- function(checkpoint) {
  g205n_need(identical(checkpoint$schema, "GSE205906_MATCHED_SAVED_NATIVE_CHECKPOINT/v1") &&
    isTRUE(checkpoint$statistical_fit_executed) && !is.null(checkpoint$result), "New fitted checkpoint required")
  packet <- checkpoint$prepared_packet
  result <- checkpoint$result
  g205n_need(identical(result$schema, "gse205906-matched-downstream-result-v1") &&
    identical(result$binding, packet$binding) &&
    identical(result$prepared_sha256, packet$sha256) &&
    identical(result$candidate_admission, FALSE), "Saved prepared/result association differs")
  g205n_need(identical(as.character(getRversion()), result$runtime$R), "Saved R version required")
  for (package in c("edgeR", "limma", "DESeq2")) {
    g205n_need(requireNamespace(package, quietly = TRUE) &&
      identical(as.character(utils::packageVersion(package)), result$runtime$packages[[package]]), "Saved package version required")
  }
  counts <- packet$counts
  base <- result$native$base_fit
  contrasts <- packet$design$contrast
  g205n_need(is.matrix(contrasts) && identical(rownames(contrasts), c("sham", "SNI", "PLM")) &&
    identical(colnames(contrasts), c("local_disease", "intervention")) &&
    identical(as.numeric(contrasts), c(-1, 1, 0, 0, -1, 1)), "Frozen contrast directions differ")
  g205n_need(is.matrix(counts) && nrow(counts) == 32883L && ncol(counts) == 9L &&
    identical(colnames(counts), paste0("GSM", 6235292:6235300)) &&
    identical(colnames(counts), rownames(packet$design$design)) &&
    all(is.finite(counts)) && all(counts >= 0 & counts == floor(counts)), "Full new input domain differs")
  expected_design <- diag(3L)[rep(seq_len(3L), each = 3L), , drop = FALSE]
  g205n_need(identical(colnames(packet$design$design), rownames(contrasts)) &&
    identical(dim(packet$design$design), c(9L, 3L)) &&
    identical(as.numeric(packet$design$design), as.numeric(expected_design)), "Frozen no-intercept group design differs")
  y0 <- edgeR::DGEList(counts = counts)
  expected_filter <- edgeR::filterByExpr(y0, design = packet$design$design)
  expected_cpm <- rowSums(edgeR::cpm(y0, log = FALSE) >= 1) >= 3L
  g205n_need(identical(base$keep_filter_by_expr, expected_filter) &&
    identical(base$keep_cpm, expected_cpm) &&
    identical(base$keep, expected_filter & expected_cpm), "Own full-design filtering differs")
  retained <- counts[base$keep, , drop = FALSE]
  g205n_need(identical(base$y$counts, retained), "Own retained counts differ")
  g205n_equal(base$y$samples$lib.size, unname(colSums(retained)), "Reset retained library sizes")
  ids <- rownames(retained)
  g205n_edger(result$tables$gene_edger, result$native$QL_tests, ids)
  g205n_need(identical(result$native$voom_fit$contrasts, contrasts), "Saved voom contrast matrix differs")
  g205n_limma(result$tables$gene_voom, result$native$voom_fit, ids, "logFC", "gene_id")
  dds <- result$native$DESeq_fit
  g205n_need(identical(rownames(dds), ids) && identical(colnames(dds), colnames(retained)) &&
    identical(DESeq2::counts(dds), retained), "Saved DDS count/gene/sample domain differs")
  names_dds <- DESeq2::resultsNames(dds)
  fit_type <- attr(DESeq2::dispersionFunction(dds), "fitType")
  g205n_need(is.character(fit_type) && length(fit_type) == 1L &&
    fit_type %in% c("parametric", "local") &&
    identical(fit_type, result$diagnostics$actual_DESeq_fit_type) &&
    is.character(result$tables$gene_deseq$actual_fit_type_after_DESeq) &&
    length(result$tables$gene_deseq$actual_fit_type_after_DESeq) == nrow(result$tables$gene_deseq) &&
    !anyNA(result$tables$gene_deseq$actual_fit_type_after_DESeq) &&
    all(result$tables$gene_deseq$actual_fit_type_after_DESeq == fit_type), "Saved DESeq fit type differs")
  expected_names <- paste0("group", rownames(contrasts))
  g205n_need(all(expected_names %in% names_dds), "Saved DESeq coefficient domain differs")
  for (j in seq_len(ncol(contrasts))) {
    coefficients <- stats::setNames(numeric(length(names_dds)), names_dds)
    coefficients[expected_names] <- contrasts[, j]
    # Saved-DDS extraction; no DESeq(), iterative fitting or contrast-SE shortcut.
    native <- DESeq2::results(result$native$DESeq_fit, contrast = coefficients,
      alpha = .05, independentFiltering = FALSE, cooksCutoff = FALSE)
    g205n_deseq_result(result$tables$gene_deseq, native, colnames(contrasts)[[j]])
  }
  g205n_need(is.character(result$tables$activity_results$collection) &&
    !anyNA(result$tables$activity_results$collection) &&
    all(result$tables$activity_results$collection %in% c("PROGENy", "DoRothEA_AB")), "Activity collection domain differs")
  for (collection in c("PROGENy", "DoRothEA_AB")) {
    scores <- result$own_activity_scores[[collection]]
    fit <- result$native$activity_fits[[collection]]$fit
    if (!nrow(scores)) {
      g205n_need(is.null(fit) && !any(result$tables$activity_results$collection == collection), "Empty activity family fitted")
    } else {
      g205n_need(identical(fit$contrasts, contrasts), "Saved activity contrast matrix differs")
      g205n_limma(result$tables$activity_results, fit, rownames(scores), "effect", "program_id", collection)
    }
  }
  for (kind in c("camera", "mroast", "fgsea")) {
    g205n_unsigned_q(result$tables[[paste0("unsigned_", kind)]], kind)
  }
  list(status = "HELD_SAVED_NATIVE_COMPONENT_NUMERICS_VERIFIED",
    relative_tolerance = 1e-10, absolute_tolerance = 1e-12,
    filtering_and_reset_library_sizes_verified = TRUE,
    gene_native_extraction_and_applicable_intervals_verified = TRUE,
    activity_saved_fit_extraction_verified = TRUE,
    unsigned_native_P_alias_and_family_BH_only_verified = TRUE,
    unsigned_P_replay_mapping_activity_score_and_TMM_reproduction = "NOT_PERFORMED",
    external_custody_TSV_transport_and_actual_parent_exit_authenticated = FALSE,
    iterative_model_fits = 0L, candidate_admission = FALSE,
    full_independent_numeric_acceptance = FALSE)
}
