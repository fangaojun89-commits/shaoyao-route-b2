# Functions only. No file reads/writes, source(), options(), RNG reset or driver.
# Source-faithful extraction from analyze_gse128623_exact30_linux_model_v2.R
# SHA256 e251c6ae9bb70cd944ebb1d781895878ee80963946f6e8ad94a52a83713add50.
# Call only with explicitly admitted in-memory inputs. This preparation is not
# permission to run real label references. No canonical reparameterization or
# allocation/partition cache is introduced here.

g128lc_contract <- function() {
  groups <- c("control", "D4", "D7", "D14", "D28", "D42")
  contrasts <- cbind(
    chronic_combined = c(-1, 0, 0, 0, .5, .5),
    D28_vs_control = c(-1, 0, 0, 0, 1, 0),
    D42_vs_control = c(-1, 0, 0, 0, 0, 1),
    acute_combined = c(-1, .5, .5, 0, 0, 0),
    D4_vs_control = c(-1, 1, 0, 0, 0, 0),
    D7_vs_control = c(-1, 0, 1, 0, 0, 0),
    D14_vs_control = c(-1, 0, 0, 1, 0, 0)
  )
  rownames(contrasts) <- groups
  list(
    groups = groups, contrasts = contrasts,
    source_sha256 = "e251c6ae9bb70cd944ebb1d781895878ee80963946f6e8ad94a52a83713add50",
    R_version = "4.6.1",
    packages = c(edgeR = "4.10.1", limma = "3.68.4", DESeq2 = "1.52.0",
                 data.table = "1.18.4"),
    camera_collections = c("Reactome", "SynGO_assignedBy_proxy_sensitivity"),
    minimum_members = 10L, maximum_members = 500L,
    edge_ci_status = "NOT_AVAILABLE_FROM_STANDARD_EDGER_QLF_OUTPUT"
  )
}

g128lc_assert <- function(condition, label) {
  if (!isTRUE(condition)) stop(label, call. = FALSE)
}

g128lc_require_runtime <- function() {
  contract <- g128lc_contract()
  g128lc_assert(as.character(getRversion()) == contract$R_version, "R version differs")
  for (package in names(contract$packages)) {
    g128lc_assert(requireNamespace(package, quietly = TRUE), paste("Missing package", package))
    g128lc_assert(as.character(utils::packageVersion(package)) == contract$packages[[package]],
                 paste("Package version differs", package))
  }
  invisible(TRUE)
}

g128lc_validate_inputs <- function(counts, metadata, annotation, labels,
                                   annotation_feature_id = "gene_id",
                                   gene_symbol_column = "gene_symbol") {
  contract <- g128lc_contract()
  g128lc_assert(is.matrix(counts) && is.numeric(counts) && nrow(counts) > 1L &&
                 ncol(counts) %in% c(29L, 30L), "Numeric full 29/30-library count matrix required")
  for (keys in dimnames(counts)) {
    g128lc_assert(!is.null(keys) && !anyNA(keys) && all(nzchar(keys)) && !anyDuplicated(keys),
                 "Missing, empty or duplicate matrix keys")
  }
  g128lc_assert(length(dimnames(counts)) == 2L && !is.null(rownames(counts)) &&
                 !is.null(colnames(counts)), "Matrix row/column keys required")
  g128lc_assert(all(is.finite(counts)) && all(counts >= 0) &&
                 all(abs(counts - round(counts)) <= .Machine$double.eps^.5),
               "Counts must obey source finite nonnegative near-integer rule")
  # Same read-input rounding as the source, not fractional-RSEM approximation.
  rounded <- round(counts)
  g128lc_assert(all(colSums(rounded) > 0), "Zero library")
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  g128lc_assert(all(c("sample_id", "group") %in% names(metadata)) &&
                 nrow(metadata) == ncol(counts) &&
                 identical(as.character(metadata$sample_id), colnames(counts)),
               "Metadata must match fixed sample order exactly")
  g128lc_assert(is.character(labels) && !anyNA(labels) &&
                 identical(names(labels), as.character(metadata$sample_id)),
               "Character allocation labels must carry exact fixed sample IDs in order")
  baseline <- factor(as.character(metadata$group), levels = contract$groups)
  allocated <- factor(as.character(labels), levels = contract$groups)
  expected <- if (ncol(counts) == 30L) rep(5L, 6L) else c(5L, 4L, 5L, 5L, 5L, 5L)
  g128lc_assert(!anyNA(baseline) && identical(as.integer(table(baseline)), expected),
               "Frozen cohort group populations differ")
  g128lc_assert(length(allocated) == ncol(counts) && !anyNA(allocated) &&
                 identical(as.integer(table(allocated)), expected),
               "Allocation must preserve all fixed group populations")
  g128lc_assert(all(c(annotation_feature_id, gene_symbol_column) %in% names(annotation)),
               "Annotation columns missing")
  ids <- trimws(as.character(annotation[[annotation_feature_id]]))
  g128lc_assert(nrow(annotation) == nrow(counts) && identical(ids, rownames(counts)),
               "Annotation gene order differs")
  symbols <- trimws(as.character(annotation[[gene_symbol_column]]))
  symbols[is.na(symbols) | !nzchar(symbols)] <- NA_character_
  list(counts = rounded, metadata = metadata, symbols = symbols, labels = allocated)
}

g128lc_design <- function(labels) {
  group <- factor(as.character(labels), levels = g128lc_contract()$groups)
  g128lc_assert(!anyNA(group), "Unknown allocation label")
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- g128lc_contract()$groups
  g128lc_assert(qr(design)$rank == 6L, "Allocation design is not full rank")
  design
}

g128lc_filter <- function(counts, design, groups) {
  y0 <- edgeR::DGEList(counts = counts)
  by_expr <- edgeR::filterByExpr(y0, design = design)
  min_n <- min(table(groups))
  by_cpm <- rowSums(edgeR::cpm(y0, log = FALSE) >= 1) >= min_n
  keep <- by_expr & by_cpm
  g128lc_assert(sum(keep) >= 2L, "Too few retained features")
  list(y0 = y0, keep = keep, min_group_n = as.integer(min_n),
       by_expr_n = sum(by_expr), cpm_n = sum(by_cpm))
}

g128lc_edge_tables <- function(fit, contrast_matrix, gene_ids) {
  bundles <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
    name <- colnames(contrast_matrix)[[index]]
    test <- edgeR::glmQLFTest(fit, contrast = contrast_matrix[, index])
    tab <- data.table::as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table)
    g128lc_assert(nrow(tab) == length(gene_ids), paste("edgeR row count", name))
    table <- data.table::data.table(
      gene_id = gene_ids, contrast = name, logFC = tab$logFC,
      logCPM = tab$logCPM, QLF_F = tab$F, df1 = 1,
      df2 = as.numeric(test$df.total), PValue = tab$PValue,
      BH_q = stats::p.adjust(tab$PValue, method = "BH"),
      tested_feature_count = length(gene_ids), ci_low = NA_real_, ci_high = NA_real_,
      ci_status = g128lc_contract()$edge_ci_status, method = "edgeR_TMM_robust_QL_primary",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
    list(table = table, test = test)
  })
  tests <- lapply(bundles, `[[`, "test")
  names(tests) <- colnames(contrast_matrix)
  list(table = data.table::rbindlist(lapply(bundles, `[[`, "table")), tests = tests)
}

g128lc_limma_tables <- function(voom, design, contrast_matrix, gene_ids) {
  fit <- limma::lmFit(voom, design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), robust = TRUE)
  rows <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
    se <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
    df <- fit$df.total
    critical <- stats::qt(.975, df = df)
    effect <- fit$coefficients[, index]
    p <- fit$p.value[, index]
    data.table::data.table(
      gene_id = gene_ids, contrast = colnames(contrast_matrix)[[index]],
      logFC = effect, SE = se, df = df, t = fit$t[, index],
      ci_low = effect - critical * se, ci_high = effect + critical * se,
      PValue = p, BH_q = stats::p.adjust(p, method = "BH"),
      tested_feature_count = length(gene_ids), method = "limma_voom_robust_sensitivity",
      ci_status = "METHOD_SPECIFIC_LIMMA_VOOM_95_PERCENT",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  })
  list(table = data.table::rbindlist(rows), fit = fit)
}

g128lc_deseq_tables <- function(counts, sample_ids, labels, keep, contrast_matrix) {
  group <- factor(as.character(labels), levels = g128lc_contract()$groups)
  col_data <- data.frame(group = group, row.names = as.character(sample_ids))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts[keep, , drop = FALSE]), colData = col_data, design = ~ 0 + group
  )
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", betaPrior = FALSE,
                       minReplicatesForReplace = Inf, quiet = TRUE)
  g128lc_assert(identical(DESeq2::resultsNames(dds), paste0("group", g128lc_contract()$groups)),
               "DESeq2 coefficient order differs")
  bundles <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
    result <- DESeq2::results(dds, contrast = contrast_matrix[, index], alpha = .05,
                              independentFiltering = FALSE, cooksCutoff = FALSE)
    effect <- as.numeric(result$log2FoldChange)
    se <- as.numeric(result$lfcSE)
    p <- as.numeric(result$pvalue)
    z <- stats::qnorm(.975)
    table <- data.table::data.table(
      gene_id = rownames(result), contrast = colnames(contrast_matrix)[[index]],
      logFC = effect, lfcSE = se, WaldStatistic = as.numeric(result$stat),
      ci_low = effect - z * se, ci_high = effect + z * se,
      PValue = p, BH_q = stats::p.adjust(p, method = "BH"),
      tested_feature_count = nrow(result), method = "DESeq2_Wald_sensitivity",
      ci_status = "METHOD_SPECIFIC_DESEQ2_WALD_95_PERCENT",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
    list(table = table, native = result)
  })
  native_results <- lapply(bundles, `[[`, "native")
  names(native_results) <- colnames(contrast_matrix)
  # Diagnostics do not alter native P/q, impute missing values or drop rows.
  gene_metadata <- as.data.frame(S4Vectors::mcols(dds))
  beta <- if ("betaConv" %in% names(gene_metadata)) gene_metadata$betaConv else
    rep(NA, nrow(dds))
  size_factors <- DESeq2::sizeFactors(dds)
  dispersion <- DESeq2::dispersions(dds)
  list(
    table = data.table::rbindlist(lapply(bundles, `[[`, "table")), fit = dds, results = native_results,
    diagnostics = list(
      requested_fit_type = "parametric",
      actual_fit_type = attr(DESeq2::dispersionFunction(dds), "fitType"),
      beta_convergence_available = "betaConv" %in% names(gene_metadata),
      beta_nonconverged_n = sum(!beta, na.rm = TRUE), beta_missing_n = sum(is.na(beta)),
      size_factors = size_factors,
      nonfinite_size_factor_n = sum(!is.finite(size_factors)),
      nonpositive_size_factor_n = sum(is.finite(size_factors) & size_factors <= 0),
      nonfinite_dispersion_n = sum(!is.finite(dispersion)),
      gene_metadata = gene_metadata,
      cooks_cutoff_applied = FALSE, independent_filtering = FALSE,
      min_replicates_for_replace = Inf,
      replacement_assay_present = "replaceCounts" %in% SummarizedExperiment::assayNames(dds)
    )
  )
}

g128lc_camera_indices <- function(symbols, members, minimum_members = 10L) {
  g128lc_assert(is.numeric(minimum_members) && length(minimum_members) == 1L &&
                 is.finite(minimum_members) && minimum_members == 10,
               "Fixed minimum CAMERA members is exactly ten")
  members <- data.table::copy(data.table::as.data.table(members))
  required <- c("collection", "program_id", "program_name", "gene_symbol", "test_class")
  g128lc_assert(all(required %in% names(members)), "CAMERA member columns missing")
  g128lc_assert(all(vapply(members[, ..required], is.character, logical(1))),
               "CAMERA member identifiers must be character")
  unsigned <- members[grepl("unsigned", test_class)]
  g128lc_assert(all(unique(unsigned$collection) %in% g128lc_contract()$camera_collections),
               "Unsupported unsigned collection")
  symbol_to_row <- split(seq_along(symbols), symbols)
  result <- lapply(unique(unsigned$collection), function(collection_name) {
    block <- unsigned[collection == collection_name]
    registry <- unique(block[, .(program_id, program_name)])
    g128lc_assert(!anyNA(registry$program_id) && all(nzchar(registry$program_id)) &&
                   !anyDuplicated(registry$program_id), "Ambiguous CAMERA registry keys")
    index_list <- lapply(registry$program_id, function(id) {
      genes <- unique(block[program_id == id, gene_symbol])
      unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    })
    names(index_list) <- registry$program_id
    eligible <- lengths(index_list) >= minimum_members & lengths(index_list) <= 500L
    list(collection = collection_name, registry = registry, eligible = eligible,
         indices = index_list[eligible], all_index_lengths = lengths(index_list))
  })
  names(result) <- unique(unsigned$collection)
  result
}

g128lc_camera_tables <- function(voom, symbols, design, contrast_matrix, members,
                                  minimum_members = 10L) {
  indices <- g128lc_camera_indices(symbols, members, minimum_members)
  rows <- list()
  native <- list()
  cursor <- 1L
  for (block in indices) {
    if (length(block$indices) == 0L) next
    for (j in seq_len(ncol(contrast_matrix))) {
      payload <- limma::camera(
        voom, index = block$indices, design = design, contrast = contrast_matrix[, j],
        inter.gene.cor = .01, sort = FALSE
      )
      result <- data.table::as.data.table(payload, keep.rownames = "program_id")
      result[, `:=`(
        collection = block$collection, contrast = colnames(contrast_matrix)[[j]],
        BH_q = stats::p.adjust(PValue, "BH"), method = "camera_inter_gene_cor_0.01",
        package_version = as.character(utils::packageVersion("limma")),
        inter_gene_correlation = .01, sort_rule = "input_program_id_order",
        formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
        s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
      )]
      rows[[cursor]] <- result
      native[[paste(block$collection, colnames(contrast_matrix)[[j]], sep = "::")]] <- payload
      cursor <- cursor + 1L
    }
  }
  list(table = if (length(rows)) data.table::rbindlist(rows, fill = TRUE) else
         data.table::data.table(), indices = indices, native = native)
}

g128lc_fit_gene_camera <- function(counts, metadata, annotation, members, labels,
                                   annotation_feature_id = "gene_id",
                                   gene_symbol_column = "gene_symbol") {
  g128lc_require_runtime()
  input <- g128lc_validate_inputs(counts, metadata, annotation, labels,
                                   annotation_feature_id, gene_symbol_column)
  design <- g128lc_design(input$labels)
  contrast_matrix <- g128lc_contract()$contrasts
  filtering <- g128lc_filter(input$counts, design, input$labels)
  y <- edgeR::normLibSizes(filtering$y0[filtering$keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  ql_fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  voom <- limma::voom(y, design, plot = FALSE)
  gene_ids <- rownames(y)
  edge <- g128lc_edge_tables(ql_fit, contrast_matrix, gene_ids)
  limma <- g128lc_limma_tables(voom, design, contrast_matrix, gene_ids)
  deseq <- g128lc_deseq_tables(input$counts, input$metadata$sample_id, input$labels,
                               filtering$keep, contrast_matrix)
  camera <- g128lc_camera_tables(voom, input$symbols[filtering$keep], design,
                                  contrast_matrix, members)
  list(
    tables = list(edgeR = edge$table, voom = limma$table, DESeq2 = deseq$table,
                  CAMERA = camera$table),
    native = list(dge = y, ql_fit = ql_fit, ql_tests = edge$tests, voom = voom,
                  limma_fit = limma$fit, deseq_fit = deseq$fit,
                  deseq_results = deseq$results, camera = camera$native),
    filtering = filtering, design = design, contrasts = contrast_matrix,
    retained_gene_ids = gene_ids, retained_symbols = input$symbols[filtering$keep],
    camera_indices = camera$indices, deseq_diagnostics = deseq$diagnostics,
    labels = as.character(input$labels),
    scope = list(candidate_generation_allowed = FALSE, no_historical_fit_identity_claim = TRUE,
                 allocation_cache_used = FALSE, real_execution_authorization_claimed = FALSE)
  )
}
