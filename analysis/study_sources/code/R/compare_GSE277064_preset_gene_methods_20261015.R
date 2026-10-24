#!/usr/bin/env Rscript
# Descriptive comparison of existing edgeR gene tables with accepted native
# readbacks. No model, moderation, original QLF-derived SE or CI is evaluated.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
assert_ok <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
file_id <- function(path) list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
  bytes = unname(file.info(path)$size),
  sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
verify_id <- function(binding) {
  actual <- file_id(binding$path)
  assert_ok(actual$sha256 == binding$sha256 && actual$bytes == binding$bytes, "Input identity changed")
}
read_selected <- function(path, selected = NULL) {
  connection <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  header <- scan(text = readLines(connection, n = 1L, warn = FALSE), what = character(),
    sep = "\t", quote = '"', quiet = TRUE)
  close(connection)
  connection <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  column_classes <- rep(NA_character_, length(header))
  if (!is.null(selected)) {
    assert_ok(all(selected %in% header), "Selected table fields absent")
    column_classes[!header %in% selected] <- "NULL"
  }
  read.delim(connection, sep = "\t", quote = '"', na.strings = "NA",
    comment.char = "", check.names = FALSE, colClasses = column_classes,
    stringsAsFactors = FALSE, encoding = "UTF-8")
}
manual_bh <- function(p) {
  answer <- rep(NA_real_, length(p))
  valid <- which(is.finite(p))
  if (!length(valid)) return(answer)
  permutation <- order(p[valid], method = "radix")
  answer[valid[permutation]] <- pmin(1, rev(cummin(rev(
    p[valid[permutation]] * length(valid) / seq_along(valid)))))
  answer
}

arguments <- commandArgs(trailingOnly = TRUE)
assert_ok(length(arguments) == 3L, "PROJECT_ROOT COMPLETED_RUN FRESH_COMPARISON_SUBDIRECTORY required")
root <- normalizePath(arguments[1L], winslash = "/", mustWork = TRUE)
run_root <- normalizePath(arguments[2L], winslash = "/", mustWork = TRUE)
output <- gsub("\\\\", "/", arguments[3L])
assert_ok(!file.exists(output) &&
  identical(normalizePath(dirname(output), winslash = "/", mustWork = TRUE), run_root) &&
  grepl("^method_comparison_[0-9]+$", basename(output)), "Use a fresh comparison subdirectory")
receipt_path <- file.path(run_root, "independent_readback_01", "INDEPENDENT_READBACK_RECEIPT.json")
independent <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
assert_ok(independent$status == "PASS_SCOPED_NATIVE_RESULT_READBACK_NOT_WHOLE_GOAL_ACCEPTANCE" &&
  isTRUE(independent$all_field_comparisons_pass) && independent$model_refit_n == 0L &&
  independent$contrast_families == 36L, "Accepted scoped native readback required")
complete_path <- file.path(run_root, "COMPLETE.json")
verify_id(independent$completion)
completion <- jsonlite::read_json(complete_path, simplifyVector = FALSE)
readiness <- jsonlite::read_json(independent$readiness$path, simplifyVector = FALSE)
for (binding in readiness$inputs) verify_id(binding)
inputs <- list(file_id(receipt_path), file_id(complete_path), file_id(independent$readiness$path))
families <- pairs <- list()
contrast_names <- c("chronic_combined", "D14_vs_control", "D28_vs_control", "early_combined",
  "D3_vs_control", "D7_vs_control")
for (cell in c("astrocyte", "microglia", "v2a")) {
  new_tables <- lapply(c("saved_voom", "integer_DESeq2"), function(method) {
    binding <- completion$products[[paste0(cell, "_", method, "_table")]]
    verify_id(binding)
    inputs[[length(inputs) + 1L]] <<- file_id(binding$path)
    read_selected(binding$path, c("gene_id", "method", "contrast", "log2_effect",
      "p_value", "bh_q", "inference_eligible", "cell_type"))
  })
  names(new_tables) <- c("saved_voom_limma", "DESeq2_integer_approximation")
  summary_path <- file.path(root, "results/GSE277064/gene_level", paste0(cell, "_contrast_summary.tsv"))
  original_summary <- read_selected(summary_path)
  inputs[[length(inputs) + 1L]] <- file_id(summary_path)
  for (contrast_name in contrast_names) {
    edge_path <- file.path(root, "results/GSE277064/gene_level", cell,
      paste0(contrast_name, ".tsv.gz"))
    inputs[[length(inputs) + 1L]] <- file_id(edge_path)
    # Explicit projection excludes all historical ql_se/ci_low/ci_high columns.
    edge <- read_selected(edge_path, c("gene_id_versioned", "logFC", "PValue", "FDR",
      "contrast", "cell_type", "domain"))
    assert_ok(!anyDuplicated(edge$gene_id_versioned) && !anyNA(edge$gene_id_versioned) &&
      all(edge$contrast == contrast_name) && all(edge$cell_type == cell) &&
      all(edge$domain == "natural_course"), "Existing edgeR family identity differs")
    assert_ok(all(is.na(edge$PValue) | (is.finite(edge$PValue) & edge$PValue >= 0 & edge$PValue <= 1)) &&
      all(is.na(edge$FDR) | (is.finite(edge$FDR) & edge$FDR >= 0 & edge$FDR <= 1)),
      "Existing edgeR probability domain differs")
    adjusted <- manual_bh(edge$PValue)
    assert_ok(identical(is.na(adjusted), is.na(edge$FDR)) &&
      all(abs(edge$FDR - adjusted) <= 1e-10 + 1e-8 * abs(adjusted), na.rm = TRUE),
      "Existing edgeR FDR is inconsistent with its native finite-P BH family")
    old_summary <- original_summary[original_summary$contrast == contrast_name, , drop = FALSE]
    assert_ok(nrow(old_summary) == 1L && old_summary$tested_features == nrow(edge) &&
      old_summary$FDR_0_05 == sum(edge$FDR <= .05, na.rm = TRUE) &&
      old_summary$FDR_0_05_up == sum(edge$FDR <= .05 & edge$logFC > 0, na.rm = TRUE) &&
      old_summary$FDR_0_05_down == sum(edge$FDR <= .05 & edge$logFC < 0, na.rm = TRUE),
      "Existing edgeR summary disagrees with current saved gene table")
    tables <- list(edgeR_QL = data.frame(gene_id = edge$gene_id_versioned,
      log2_effect = edge$logFC, p_value = edge$PValue, bh_q = edge$FDR,
      inference_eligible = is.finite(edge$logFC) & is.finite(edge$PValue) & is.finite(edge$FDR)))
    for (method in names(new_tables)) {
      selected <- new_tables[[method]][new_tables[[method]]$contrast == contrast_name, , drop = FALSE]
      assert_ok(identical(selected$gene_id, edge$gene_id_versioned) &&
        all(selected$cell_type == cell) && all(selected$method == method),
        "Gene universe/order or new family labels differ")
      tables[[method]] <- selected[c("gene_id", "log2_effect", "p_value", "bh_q", "inference_eligible")]
    }
    for (method in names(tables)) {
      selected <- tables[[method]]
      native_sig <- is.finite(selected$bh_q) & selected$bh_q <= .05
      eligible_sig <- selected$inference_eligible & native_sig
      families[[length(families) + 1L]] <- data.frame(cell_type = cell, contrast = contrast_name,
        method = method, retained_gene_n = nrow(selected), finite_effect_n = sum(is.finite(selected$log2_effect)),
        finite_p_n = sum(is.finite(selected$p_value)), finite_q_n = sum(is.finite(selected$bh_q)),
        missing_p_n = sum(is.na(selected$p_value)), missing_q_n = sum(is.na(selected$bh_q)),
        inference_eligible_n = sum(selected$inference_eligible),
        native_q_le_005_n = sum(native_sig), eligible_q_le_005_n = sum(eligible_sig),
        eligible_q_le_005_up_n = sum(eligible_sig & selected$log2_effect > 0),
        eligible_q_le_005_down_n = sum(eligible_sig & selected$log2_effect < 0),
        median_absolute_log2_effect = stats::median(abs(selected$log2_effect), na.rm = TRUE))
    }
    method_pairs <- combn(names(tables), 2L, simplify = FALSE)
    for (method_pair in method_pairs) {
      left <- tables[[method_pair[1L]]]
      right <- tables[[method_pair[2L]]]
      paired <- is.finite(left$log2_effect) & is.finite(right$log2_effect)
      nonzero <- paired & left$log2_effect != 0 & right$log2_effect != 0
      direction_n <- sum(sign(left$log2_effect[nonzero]) == sign(right$log2_effect[nonzero]))
      significant <- function(value, eligible) value$gene_id[is.finite(value$bh_q) &
        value$bh_q <= .05 & (!eligible | value$inference_eligible)]
      native_left <- significant(left, FALSE)
      native_right <- significant(right, FALSE)
      eligible_left <- significant(left, TRUE)
      eligible_right <- significant(right, TRUE)
      native_union <- union(native_left, native_right)
      eligible_union <- union(eligible_left, eligible_right)
      native_intersection <- intersect(native_left, native_right)
      eligible_intersection <- intersect(eligible_left, eligible_right)
      rho <- if (sum(paired) >= 2L && stats::sd(left$log2_effect[paired]) > 0 &&
        stats::sd(right$log2_effect[paired]) > 0) stats::cor(left$log2_effect[paired],
          right$log2_effect[paired], method = "spearman") else NA_real_
      pairs[[length(pairs) + 1L]] <- data.frame(cell_type = cell, contrast = contrast_name,
        left_method = method_pair[1L], right_method = method_pair[2L],
        common_gene_n = nrow(left), paired_finite_effect_n = sum(paired),
        log2_effect_spearman = rho, paired_nonzero_effect_n = sum(nonzero),
        direction_agreement_n = direction_n,
        direction_agreement_fraction = if (sum(nonzero)) direction_n / sum(nonzero) else NA_real_,
        left_zero_effect_n = sum(paired & left$log2_effect == 0),
        right_zero_effect_n = sum(paired & right$log2_effect == 0),
        both_zero_effect_n = sum(paired & left$log2_effect == 0 & right$log2_effect == 0),
        native_qset_left_n = length(native_left), native_qset_right_n = length(native_right),
        native_qset_intersection_n = length(native_intersection), native_qset_union_n = length(native_union),
        native_qset_jaccard = if (length(native_union)) length(native_intersection) / length(native_union) else NA_real_,
        eligible_qset_left_n = length(eligible_left), eligible_qset_right_n = length(eligible_right),
        eligible_qset_intersection_n = length(eligible_intersection), eligible_qset_union_n = length(eligible_union),
        eligible_qset_jaccard = if (length(eligible_union)) length(eligible_intersection) / length(eligible_union) else NA_real_)
    }
  }
}
family_table <- do.call(rbind, families)
pair_table <- do.call(rbind, pairs)
assert_ok(nrow(family_table) == 54L && nrow(pair_table) == 54L &&
  length(unique(paste(pair_table$cell_type, pair_table$contrast))) == 18L,
  "Complete three-method by 18-contrast coverage required")
for (binding in inputs) verify_id(binding)
for (binding in readiness$inputs) verify_id(binding)
assert_ok(dir.create(output, recursive = FALSE), "Cannot create fresh comparison output")
outputs <- list()
for (name in c("METHOD_FAMILY_COUNTS", "PAIRWISE_METHOD_COMPARISONS")) {
  value <- if (name == "METHOD_FAMILY_COUNTS") family_table else pair_table
  path <- file.path(output, paste0(name, ".tsv"))
  write.table(value, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  outputs[[name]] <- file_id(path)
}
script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
jsonlite::write_json(list(status = "DESCRIPTIVE_SAME_DATA_THREE_METHOD_COMPARISON_COMPLETE",
  script = file_id(script), independent_native_readback = file_id(receipt_path), inputs = inputs,
  outputs = outputs, method_family_rows = 54L, added_method_family_rows = 36L,
  distinct_cell_contrasts = 18L, pairwise_comparison_rows = 54L,
  method_refits = 0L, old_edgeR_SE_or_CI_used = FALSE, all_inputs_unchanged = TRUE,
  methods = "Original robust edgeR QL; saved voom robust limma; rounded-RSEM DESeq2",
  log2_effect_spearman = "All paired finite gene effects; average ranks for ties; not significance-selected",
  direction_fraction = "Agreement among paired finite nonzero effects; zero counts reported separately",
  qsets = "Native within-cell/method/contrast BH q<=0.05, native and eligibility-restricted sets reported separately; empty union Jaccard NA",
  interpretation = "Descriptive source-wide same-data method sensitivity, not independent replication, true-FDR calibration, equivalence, mechanism or efficacy",
  library_observations = 45L, authenticated_independent_animal_n = NA_integer_,
  original_primary_or_gates_changed = FALSE, whole_goal_complete = FALSE),
  file.path(output, "COMPARISON_RECEIPT.json"), auto_unbox = TRUE, pretty = TRUE,
  digits = NA, na = "null", null = "null")
cat("DESCRIPTIVE_THREE_METHOD_COMPARISON_COMPLETE; 54 method families, 54 pairs, 18 cell contrasts; zero fits\n")
