#!/usr/bin/env Rscript

# Fresh exact30 model-v2 implementation. It does not source, transform, or
# execute any v4/v7 analyzer. The frozen v2 contract is deliberately NO-GO;
# command-line execution therefore refuses before reading outcome counts until
# an independently reviewed successor supplies every future pin.

V2A_CANDIDATE <- "GSE128623_exact30_linux_model_v2"
V2A_GROUPS <- c("control", "D4", "D7", "D14", "D28", "D42")
V2A_CONTRASTS <- list(
  chronic_combined = c(-1, 0, 0, 0, 0.5, 0.5),
  D28_vs_control = c(-1, 0, 0, 0, 1, 0),
  D42_vs_control = c(-1, 0, 0, 0, 0, 1),
  acute_combined = c(-1, 0.5, 0.5, 0, 0, 0),
  D4_vs_control = c(-1, 1, 0, 0, 0, 0),
  D7_vs_control = c(-1, 0, 1, 0, 0, 0),
  D14_vs_control = c(-1, 0, 0, 1, 0, 0)
)
V2A_REQUIRED_PACKAGES <- c(
  data.table = "1.18.4", edgeR = "4.10.1", limma = "3.68.4",
  DESeq2 = "1.52.0", fgsea = "1.38.0", jsonlite = "2.0.0",
  yaml = "2.3.12", openssl = "2.4.2"
)
V2A_EDGE_CI_STATUS <- "NOT_AVAILABLE_FROM_STANDARD_EDGER_QLF_OUTPUT"
V2A_ENRICHMENT_SEED <- 20260819L
V2A_RNG_KIND <- c("Mersenne-Twister", "Inversion", "Rejection")
V2A_RESULT_FLAGS <- list(
  formal_sequential_inference = FALSE,
  candidate_generation_allowed = FALSE,
  s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
)

v2a_stop <- function(...) stop(sprintf(...), call. = FALSE)
v2a_assert <- function(value, ...) if (!isTRUE(value)) v2a_stop(...)

v2a_reset_enrichment_rng <- function() {
  RNGkind(kind = V2A_RNG_KIND[[1L]], normal.kind = V2A_RNG_KIND[[2L]],
          sample.kind = V2A_RNG_KIND[[3L]])
  set.seed(V2A_ENRICHMENT_SEED)
  invisible(TRUE)
}

v2a_require_runtime <- function(test_only = FALSE) {
  required_r <- "4.6.1"
  v2a_assert(as.character(getRversion()) == required_r,
             "R_VERSION_MISMATCH:%s", as.character(getRversion()))
  for (package in names(V2A_REQUIRED_PACKAGES)) {
    v2a_assert(requireNamespace(package, quietly = TRUE),
               "MISSING_PACKAGE:%s", package)
    observed <- as.character(utils::packageVersion(package))
    v2a_assert(observed == V2A_REQUIRED_PACKAGES[[package]],
               "PACKAGE_VERSION_MISMATCH:%s:%s", package, observed)
  }
  if (!test_only) {
    v2a_assert(identical(Sys.info()[["sysname"]], "Linux"), "LINUX_ONLY")
  }
  invisible(TRUE)
}

v2a_file_identity <- function(path) {
  v2a_assert(length(path) == 1L && nzchar(path) && file.exists(path) && !dir.exists(path),
             "REGULAR_FILE_REQUIRED:%s", path)
  info <- file.info(path)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  sha <- as.character(openssl::sha256(connection))
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
       bytes = unname(as.numeric(info$size)), sha256 = tolower(sha))
}

v2a_assert_pin <- function(path, pin, label) {
  v2a_assert(is.list(pin) && identical(sort(names(pin)), c("bytes", "path", "sha256")),
             "PIN_SCHEMA:%s", label)
  v2a_assert(is.character(pin$path) && length(pin$path) == 1L && nzchar(pin$path),
             "PIN_PATH:%s", label)
  v2a_assert(is.numeric(pin$bytes) && length(pin$bytes) == 1L && is.finite(pin$bytes) &&
               pin$bytes >= 0 && pin$bytes == floor(pin$bytes), "PIN_BYTES:%s", label)
  v2a_assert(is.character(pin$sha256) && length(pin$sha256) == 1L &&
               grepl("^[0-9a-f]{64}$", pin$sha256), "PIN_SHA256:%s", label)
  observed <- v2a_file_identity(path)
  v2a_assert(observed$bytes == pin$bytes && observed$sha256 == pin$sha256,
             "PIN_IDENTITY_MISMATCH:%s", label)
  observed
}

v2a_read_json <- function(path, label) {
  v2a_assert(file.exists(path) && !dir.exists(path), "JSON_MISSING:%s", label)
  value <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  v2a_assert(is.list(value) && !is.null(names(value)), "JSON_OBJECT_REQUIRED:%s", label)
  value
}

v2a_exact_names <- function(value, expected, label) {
  v2a_assert(is.list(value) && identical(sort(names(value)), sort(expected)),
             "JSON_KEYS_MISMATCH:%s", label)
}

v2a_parse_args <- function(args) {
  result <- list(contract = NULL, descriptor = NULL, staging = NULL, test_only = FALSE)
  for (arg in args) {
    if (startsWith(arg, "--contract=")) result$contract <- substring(arg, 12L)
    else if (startsWith(arg, "--descriptor=")) result$descriptor <- substring(arg, 14L)
    else if (startsWith(arg, "--staging=")) result$staging <- substring(arg, 11L)
    else if (identical(arg, "--test-only")) result$test_only <- TRUE
    else v2a_stop("UNKNOWN_ARGUMENT:%s", arg)
  }
  v2a_assert(!is.null(result$contract) && !is.null(result$descriptor) && !is.null(result$staging),
             "REQUIRED_ARGUMENT_MISSING")
  result
}

v2a_validate_contract_gate <- function(contract, test_only = FALSE) {
  v2a_assert(identical(contract$schema,
                       "route-b2-gse128623-exact30-linux-model-v2-contract-v1"),
             "CONTRACT_SCHEMA")
  v2a_assert(identical(contract$candidate_id, V2A_CANDIDATE), "CONTRACT_CANDIDATE")
  v2a_assert(identical(contract$claim_boundary$s2_status,
                       "NOT_IDENTIFIABLE_NO_S1_INPUT"), "S2_BOUNDARY")
  v2a_assert(identical(contract$claim_boundary$candidate_generation_allowed, FALSE),
             "CANDIDATE_BOUNDARY")
  v2a_assert(identical(contract$batch_policy$batch_adjustment_allowed, FALSE) &&
               identical(contract$batch_policy$combat_allowed, FALSE), "BATCH_POLICY")
  if (!test_only) {
    v2a_assert(identical(contract$production_go, TRUE),
               "PRODUCTION_NO_GO_CONTRACT_REFUSES_BEFORE_COUNTS")
    v2a_assert(identical(contract$methodological_amendment$effective, TRUE),
               "CI_AMENDMENT_NOT_EFFECTIVE")
  }
  invisible(TRUE)
}

v2a_validate_future_pin_set <- function(contract) {
  pins <- contract$admission$future_pins
  required <- c("qc_successor_release", "qc_final_closure", "count_matrix",
                "feature_annotation", "sample_ledger", "strandedness_decision",
                "linux_runtime_lock")
  v2a_assert(is.list(pins) && all(required %in% names(pins)), "FUTURE_PIN_SET_SCHEMA")
  for (name in required) {
    pin <- pins[[name]]
    v2a_assert(is.list(pin), "FUTURE_PIN_OBJECT:%s", name)
    for (field in c("path", "bytes", "sha256")) {
      v2a_assert(!is.null(pin[[field]]), "FUTURE_PIN_NULL:%s:%s", name, field)
    }
    v2a_assert(!identical(pin$path, "PENDING") && !identical(pin$sha256, "PENDING"),
               "FUTURE_PIN_PENDING:%s", name)
  }
  v2a_assert(length(pins$qc_final_closure$accepted_schema_versions) > 0L,
             "QC_ACCEPTED_SCHEMA_EMPTY")
  v2a_assert(length(pins$qc_successor_release$accepted_candidate_ids) > 0L &&
               length(pins$qc_successor_release$accepted_statuses) > 0L,
             "QC_SUCCESSOR_ACCEPTED_ID_OR_STATUS_EMPTY")
  invisible(TRUE)
}

v2a_read_s1_zero <- function(path) {
  required <- c("collection", "feature_id", "feature_name", "S1_pass", "direction",
                "primary_effect", "primary_ci_low", "primary_ci_high", "primary_q_value",
                "permutation_q_value_BH", "loo_direction_fraction",
                "original_annotation_direction_stable", "all27_direction_stable")
  tab <- data.table::fread(path, sep = "\t", header = TRUE, na.strings = c("NA", ""))
  v2a_assert(identical(names(tab), required), "S1_SCHEMA")
  v2a_assert(nrow(tab) == 1771L && !anyDuplicated(tab[, .(collection, feature_id)]),
             "S1_KEYSET")
  values <- toupper(as.character(tab$S1_pass))
  v2a_assert(all(values %in% c("TRUE", "FALSE")) && sum(values == "TRUE") == 0L,
             "S1_NONZERO_OR_INVALID")
  invisible(tab)
}

v2a_validate_sample_ledger <- function(ledger, frozen_design) {
  required <- c("sample_order", "sample_id", "biological_unit_id", "biosample",
                "sra_experiment", "run_accession", "sex", "group")
  v2a_assert(identical(names(ledger), required), "SAMPLE_LEDGER_SCHEMA")
  v2a_assert(nrow(ledger) == 30L && identical(as.integer(ledger$sample_order), 1:30),
             "SAMPLE_LEDGER_EXACT30_ORDER")
  for (name in c("sample_id", "biological_unit_id", "biosample", "sra_experiment",
                 "run_accession")) {
    values <- as.character(ledger[[name]])
    v2a_assert(all(nzchar(values)) && !anyDuplicated(values), "IDENTITY_NOT_UNIQUE:%s", name)
  }
  v2a_assert(all(ledger$sex == "female"), "COHORT_NOT_ALL_FEMALE")
  v2a_assert(identical(sort(unique(ledger$group)), sort(V2A_GROUPS)), "GROUP_SET")
  counts <- table(factor(ledger$group, levels = V2A_GROUPS))
  v2a_assert(all(counts == 5L), "GROUP_COUNTS_NOT_SIX_BY_FIVE")
  v2a_assert(identical(names(frozen_design), required) && nrow(frozen_design) == 30L,
             "FROZEN_DESIGN_SCHEMA")
  for (name in required) {
    v2a_assert(identical(as.character(ledger[[name]]), as.character(frozen_design[[name]])),
               "SAMPLE_LEDGER_DRIFT:%s", name)
  }
  invisible(TRUE)
}

v2a_read_inputs <- function(count_path, annotation_path, ledger_path, frozen_design_path,
                            count_feature_id, annotation_feature_id, gene_symbol_column) {
  ledger <- data.table::fread(ledger_path, sep = "\t", header = TRUE,
                              colClasses = "character")
  frozen <- data.table::fread(frozen_design_path, sep = "\t", header = TRUE,
                              colClasses = "character")
  v2a_validate_sample_ledger(ledger, frozen)
  counts_dt <- data.table::fread(count_path, sep = "\t", header = TRUE,
                                 na.strings = c("NA", "NaN", ""), check.names = FALSE)
  v2a_assert(identical(names(counts_dt)[1L], count_feature_id), "COUNT_FEATURE_COLUMN")
  v2a_assert(identical(names(counts_dt)[-1L], as.character(ledger$sample_id)),
             "COUNT_COLUMNS_NOT_EXACT_LEDGER_ORDER")
  gene_id <- trimws(as.character(counts_dt[[1L]]))
  v2a_assert(length(gene_id) > 0L && all(nzchar(gene_id)) && !anyDuplicated(gene_id),
             "COUNT_FEATURE_IDS")
  counts <- as.matrix(counts_dt[, -1L, with = FALSE])
  storage.mode(counts) <- "double"
  v2a_assert(all(is.finite(counts)) && all(counts >= 0) &&
               all(abs(counts - round(counts)) <= .Machine$double.eps^0.5),
             "COUNTS_NOT_FINITE_NONNEGATIVE_INTEGERS")
  counts <- round(counts)
  rownames(counts) <- gene_id
  colnames(counts) <- as.character(ledger$sample_id)
  v2a_assert(all(colSums(counts) > 0) && sum(counts) > 0, "ZERO_LIBRARY_OR_MATRIX")
  annotation <- data.table::fread(annotation_path, sep = "\t", header = TRUE,
                                  colClasses = "character", na.strings = c("NA", ""))
  v2a_assert(annotation_feature_id %in% names(annotation) && gene_symbol_column %in% names(annotation),
             "ANNOTATION_COLUMNS")
  ann_id <- trimws(as.character(annotation[[annotation_feature_id]]))
  v2a_assert(nrow(annotation) == nrow(counts) && identical(ann_id, rownames(counts)) &&
               all(nzchar(ann_id)) && !anyDuplicated(ann_id), "ANNOTATION_ORDER")
  symbols <- trimws(as.character(annotation[[gene_symbol_column]]))
  symbols[is.na(symbols) | !nzchar(symbols)] <- NA_character_
  list(counts = counts, annotation = annotation, symbols = symbols, ledger = ledger)
}

v2a_design <- function(groups) {
  group <- factor(as.character(groups), levels = V2A_GROUPS)
  v2a_assert(!anyNA(group), "UNKNOWN_GROUP")
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- V2A_GROUPS
  v2a_assert(identical(colnames(design), V2A_GROUPS) && qr(design)$rank == 6L,
             "DESIGN_RANK_OR_ORDER")
  design
}

v2a_contrast_matrix <- function() {
  matrix <- do.call(cbind, V2A_CONTRASTS)
  rownames(matrix) <- V2A_GROUPS
  v2a_assert(all(is.finite(matrix)) && all(abs(colSums(matrix)) < 1e-12),
             "CONTRAST_INVALID")
  matrix
}

v2a_filter <- function(counts, design, groups) {
  y0 <- edgeR::DGEList(counts = counts)
  by_expr <- edgeR::filterByExpr(y0, design = design)
  min_n <- min(table(groups))
  by_cpm <- rowSums(edgeR::cpm(y0, log = FALSE) >= 1) >= min_n
  keep <- by_expr & by_cpm
  v2a_assert(sum(keep) >= 2L, "TOO_FEW_RETAINED_FEATURES")
  list(y0 = y0, keep = keep, min_group_n = as.integer(min_n),
       by_expr_n = sum(by_expr), cpm_n = sum(by_cpm))
}

v2a_edge_tables <- function(fit, contrast_matrix, gene_ids) {
  data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(index) {
    name <- colnames(contrast_matrix)[[index]]
    test <- edgeR::glmQLFTest(fit, contrast = contrast_matrix[, index])
    tab <- data.table::as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table)
    v2a_assert(nrow(tab) == length(gene_ids), "EDGER_ROW_COUNT:%s", name)
    data.table::data.table(
      gene_id = gene_ids, contrast = name, logFC = tab$logFC,
      logCPM = tab$logCPM, QLF_F = tab$F, df1 = 1,
      df2 = as.numeric(test$df.total), PValue = tab$PValue,
      BH_q = stats::p.adjust(tab$PValue, method = "BH"),
      tested_feature_count = length(gene_ids), ci_low = NA_real_, ci_high = NA_real_,
      ci_status = V2A_EDGE_CI_STATUS, method = "edgeR_TMM_robust_QL_primary",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2a_limma_tables <- function(voom, design, contrast_matrix, gene_ids) {
  fit <- limma::lmFit(voom, design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), robust = TRUE)
  data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(index) {
    se <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
    df <- fit$df.total
    critical <- stats::qt(0.975, df = df)
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
  }))
}

v2a_deseq_tables <- function(counts, ledger, keep, contrast_matrix) {
  group <- factor(as.character(ledger$group), levels = V2A_GROUPS)
  col_data <- data.frame(group = group, row.names = as.character(ledger$sample_id))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts[keep, , drop = FALSE]), colData = col_data, design = ~ 0 + group
  )
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", betaPrior = FALSE,
                       minReplicatesForReplace = Inf, quiet = TRUE)
  result_names <- DESeq2::resultsNames(dds)
  expected <- paste0("group", V2A_GROUPS)
  v2a_assert(identical(result_names, expected), "DESEQ_COEFFICIENT_ORDER")
  data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(index) {
    result <- DESeq2::results(dds, contrast = contrast_matrix[, index], alpha = 0.05,
                              independentFiltering = FALSE, cooksCutoff = FALSE)
    effect <- as.numeric(result$log2FoldChange)
    se <- as.numeric(result$lfcSE)
    p <- as.numeric(result$pvalue)
    z <- stats::qnorm(0.975)
    data.table::data.table(
      gene_id = rownames(result), contrast = colnames(contrast_matrix)[[index]],
      logFC = effect, lfcSE = se, WaldStatistic = as.numeric(result$stat),
      ci_low = effect - z * se, ci_high = effect + z * se,
      PValue = p, BH_q = stats::p.adjust(p, method = "BH"),
      tested_feature_count = nrow(result), method = "DESeq2_Wald_sensitivity",
      ci_status = "METHOD_SPECIFIC_DESEQ2_WALD_95_PERCENT",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2a_collapse_symbol <- function(expression, symbols) {
  valid <- !is.na(symbols) & nzchar(symbols)
  v2a_assert(any(valid), "NO_ANNOTATED_SYMBOLS")
  sums <- rowsum(expression[valid, , drop = FALSE], group = symbols[valid], reorder = FALSE)
  denominator <- as.numeric(table(factor(symbols[valid], levels = rownames(sums))))
  sums / denominator
}

v2a_scale_rows <- function(matrix) {
  center <- rowMeans(matrix)
  spread <- apply(matrix, 1L, stats::sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  sweep(sweep(matrix, 1L, center, "-"), 1L, spread, "/")
}

v2a_program_members <- function(resources) {
  required <- c("reactome", "syngo_sensitivity", "ortholog_one_to_one",
                "progeny_mouse", "dorothea_mouse")
  v2a_assert(is.list(resources) && all(required %in% names(resources)),
             "GENE_SET_RESOURCE_SCHEMA")
  orthology <- data.table::as.data.table(resources$ortholog_one_to_one)
  mouse_map <- unique(orthology[other_tax_id == 10090L & !is.na(human_symbol) &
                                  !is.na(other_symbol),
                                .(human_symbol = as.character(human_symbol),
                                  gene_symbol = as.character(other_symbol))])
  map_human <- function(tab, id_col, name_col, collection) {
    tab <- data.table::as.data.table(tab)
    source <- data.table::data.table(
      program_id = as.character(tab[[id_col]]), program_name = as.character(tab[[name_col]]),
      human_symbol = as.character(tab$gene_symbol)
    )
    mapped <- merge(source, mouse_map, by = "human_symbol", allow.cartesian = TRUE)
    unique(mapped[, .(collection = collection, program_id, program_name, gene_symbol,
                       weight = 1, test_class = "unsigned_and_activity")])
  }
  reactome <- data.table::as.data.table(resources$reactome)
  v2a_assert(identical(unique(as.character(reactome$release)), "97"), "REACTOME_RELEASE")
  members <- data.table::rbindlist(list(
    map_human(reactome, "term_id", "term_name", "Reactome"),
    map_human(resources$syngo_sensitivity, "term_id", "term_name",
              "SynGO_assignedBy_proxy_sensitivity"),
    data.table::as.data.table(resources$progeny_mouse)[, .(
      collection = "PROGENy", program_id = as.character(pathway),
      program_name = as.character(pathway), gene_symbol = as.character(gene_symbol),
      weight = as.numeric(weight), test_class = "activity_only")],
    data.table::as.data.table(resources$dorothea_mouse)[, .(
      collection = "DoRothEA_AB", program_id = as.character(tf),
      program_name = as.character(tf), gene_symbol = as.character(target),
      weight = as.numeric(mor), test_class = "activity_only")]
  ), fill = TRUE)
  members <- unique(members[
    nzchar(program_id) & nzchar(program_name) & nzchar(collection) &
      nzchar(gene_symbol) & is.finite(weight)
  ], by = c("collection", "program_id", "gene_symbol"))
  conflict <- members[, .(n_weight = data.table::uniqueN(weight)),
                      by = .(collection, program_id, gene_symbol)]
  v2a_assert(all(conflict$n_weight == 1L), "PROGRAM_WEIGHT_CONFLICT")
  members
}

v2a_activity_tables <- function(voom, symbols, design, contrast_matrix, members,
                                minimum_members = 10L) {
  expression <- v2a_collapse_symbol(voom$E, symbols)
  z_expression <- v2a_scale_rows(expression)
  registry <- unique(members[, .(collection, program_id, program_name)])
  data.table::setorder(registry, collection, program_id)
  scores <- lapply(seq_len(nrow(registry)), function(index) {
    row <- registry[index]
    subset <- members[collection == row$collection & program_id == row$program_id &
                        gene_symbol %in% rownames(z_expression)]
    subset <- unique(subset, by = "gene_symbol")
    if (nrow(subset) < minimum_members || nrow(subset) > 500L) return(NULL)
    score <- colSums(z_expression[subset$gene_symbol, , drop = FALSE] * subset$weight) /
      sum(abs(subset$weight))
    if (!all(is.finite(score)) || !is.finite(stats::sd(score)) || stats::sd(score) == 0) return(NULL)
    list(key = row, score = as.numeric(scale(score)), members = nrow(subset))
  })
  scores <- Filter(Negate(is.null), scores)
  if (length(scores) == 0L) return(data.table::data.table())
  score_matrix <- do.call(rbind, lapply(scores, `[[`, "score"))
  rownames(score_matrix) <- vapply(scores, function(x) paste(x$key$collection, x$key$program_id,
                                                             sep = "::"), character(1))
  fit <- limma::lmFit(score_matrix, design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), robust = TRUE)
  data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(index) {
    se <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
    critical <- stats::qt(0.975, df = fit$df.total)
    effect <- fit$coefficients[, index]
    p <- fit$p.value[, index]
    keys <- data.table::tstrsplit(rownames(score_matrix), "::", fixed = TRUE)
    data.table::data.table(
      collection = keys[[1L]], program_id = keys[[2L]],
      contrast = colnames(contrast_matrix)[[index]], effect = effect, SE = se,
      ci_low = effect - critical * se, ci_high = effect + critical * se,
      PValue = p, BH_q = stats::ave(p, keys[[1L]], FUN = function(x) stats::p.adjust(x, "BH")),
      tested_program_count = length(effect), method = "weighted_voom_activity_descriptive",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2a_camera_tables <- function(voom, symbols, design, contrast_matrix, members,
                              minimum_members = 10L) {
  unsigned <- members[grepl("unsigned", test_class)]
  symbol_to_row <- split(seq_along(symbols), symbols)
  collections <- unique(unsigned$collection)
  rows <- list()
  cursor <- 1L
  for (collection_name in collections) {
    block <- unsigned[collection == collection_name]
    registry <- unique(block[, .(program_id, program_name)])
    index_list <- lapply(registry$program_id, function(id) {
      genes <- unique(block[program_id == id, gene_symbol])
      unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    })
    names(index_list) <- registry$program_id
    eligible <- lengths(index_list) >= minimum_members & lengths(index_list) <= 500L
    index_list <- index_list[eligible]
    if (length(index_list) == 0L) next
    for (j in seq_len(ncol(contrast_matrix))) {
      result <- data.table::as.data.table(limma::camera(
        voom, index = index_list, design = design, contrast = contrast_matrix[, j],
        inter.gene.cor = 0.01, sort = FALSE
      ), keep.rownames = "program_id")
      result[, `:=`(
        collection = collection_name, contrast = colnames(contrast_matrix)[[j]],
        BH_q = stats::p.adjust(PValue, "BH"), method = "camera_inter_gene_cor_0.01",
        package_version = as.character(utils::packageVersion("limma")),
        inter_gene_correlation = 0.01, sort_rule = "input_program_id_order",
        formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
        s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
      )]
      rows[[cursor]] <- result
      cursor <- cursor + 1L
    }
  }
  if (length(rows) == 0L) data.table::data.table() else data.table::rbindlist(rows, fill = TRUE)
}

v2a_chronic_program_sensitivities <- function(voom, symbols, design, contrast_matrix,
                                              members, edge_chronic, nrot = 9999L,
                                              minimum_members = 10L) {
  unsigned <- members[grepl("unsigned", test_class)]
  symbol_to_row <- split(seq_along(symbols), symbols)
  contrast <- contrast_matrix[, "chronic_combined"]
  roast_rows <- list()
  fgsea_rows <- list()
  cursor <- 1L
  rank_by_symbol <- data.table::data.table(
    symbol = symbols, gene_id = edge_chronic$gene_id,
    score = sign(edge_chronic$logFC) * sqrt(pmax(edge_chronic$QLF_F, 0))
  )
  rank_by_symbol <- rank_by_symbol[!is.na(symbol) & nzchar(symbol) & is.finite(score)]
  rank_by_symbol[, abs_score := abs(score)]
  data.table::setorder(rank_by_symbol, symbol, -abs_score, gene_id)
  rank_by_symbol <- unique(rank_by_symbol, by = "symbol")
  data.table::setorder(rank_by_symbol, -score, symbol)
  ranks <- rank_by_symbol$score
  names(ranks) <- rank_by_symbol$symbol
  ranking_rule <- "max_abs_signed_sqrt_QLF_F_per_symbol_then_score_desc_symbol_asc"
  ranking_payload <- paste0(paste(sprintf("%s\t%.17g", names(ranks), ranks),
                                  collapse = "\n"), "\n")
  ranking_sha256 <- tolower(as.character(openssl::sha256(charToRaw(ranking_payload))))
  for (collection_name in unique(unsigned$collection)) {
    block <- unsigned[collection == collection_name]
    registry <- unique(block[, .(program_id, program_name)])
    index_list <- lapply(registry$program_id, function(id) {
      genes <- unique(block[program_id == id, gene_symbol])
      unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    })
    names(index_list) <- registry$program_id
    eligible <- lengths(index_list) >= minimum_members & lengths(index_list) <= 500L
    index_list <- index_list[eligible]
    if (length(index_list) == 0L) next
    v2a_reset_enrichment_rng()
    roast <- data.table::as.data.table(limma::mroast(
      voom, index = index_list, design = design, contrast = contrast,
      nrot = nrot, set.statistic = "mean", midp = TRUE, sort = "none"
    ), keep.rownames = "program_id")
    roast[, `:=`(
      collection = collection_name, contrast = "chronic_combined",
      BH_q = pmax(PValue, stats::p.adjust(pmax(PValue - 0.5 / nrot, 0), "BH")),
      method = sprintf("mroast_mean_%d_rotations", nrot),
      package_version = as.character(utils::packageVersion("limma")),
      rng_seed = V2A_ENRICHMENT_SEED, rng_kind = V2A_RNG_KIND[[1L]],
      rng_normal_kind = V2A_RNG_KIND[[2L]], rng_sample_kind = V2A_RNG_KIND[[3L]],
      rotation_count = nrot, set_statistic = "mean", midp = TRUE,
      sort_rule = "input_program_id_order",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )]
    roast_rows[[cursor]] <- roast
    paths <- lapply(index_list, function(index) unique(symbols[index]))
    paths <- lapply(paths, function(x) intersect(x, names(ranks)))
    paths <- paths[lengths(paths) >= minimum_members & lengths(paths) <= 500L]
    if (length(paths) > 0L) {
      v2a_reset_enrichment_rng()
      fg <- data.table::as.data.table(fgsea::fgseaMultilevel(
        pathways = paths, stats = ranks, minSize = minimum_members, maxSize = 500L,
        eps = 0, scoreType = "std", nPermSimple = 1000L
      ))
      data.table::setnames(fg, "pathway", "program_id")
      fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ","),
                                 character(1L))]
      fg[, `:=`(
        collection = collection_name, contrast = "chronic_combined",
        BH_q = stats::p.adjust(pval, "BH"), method = "fgseaMultilevel_sensitivity",
        package_version = as.character(utils::packageVersion("fgsea")),
        rng_seed = V2A_ENRICHMENT_SEED, rng_kind = V2A_RNG_KIND[[1L]],
        rng_normal_kind = V2A_RNG_KIND[[2L]], rng_sample_kind = V2A_RNG_KIND[[3L]],
        ranking_rule = ranking_rule, ranking_count = length(ranks),
        ranking_sha256 = ranking_sha256, min_size = minimum_members,
        max_size = 500L, eps = 0, score_type = "std", n_perm_simple = 1000L,
        tie_rule = "score_desc_then_symbol_asc",
        formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
        s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
      )]
      fgsea_rows[[cursor]] <- fg
    }
    cursor <- cursor + 1L
  }
  list(
    mroast = if (length(roast_rows)) data.table::rbindlist(roast_rows, fill = TRUE) else data.table::data.table(),
    fgsea = if (length(fgsea_rows)) data.table::rbindlist(fgsea_rows, fill = TRUE) else data.table::data.table()
  )
}

v2a_gene_loao <- function(counts, ledger, full_gene_ids) {
  rows <- vector("list", ncol(counts))
  chronic <- V2A_CONTRASTS$chronic_combined
  for (omit in seq_len(ncol(counts))) {
    retained <- seq_len(ncol(counts)) != omit
    counts_i <- counts[, retained, drop = FALSE]
    ledger_i <- ledger[retained]
    groups_i <- factor(ledger_i$group, levels = V2A_GROUPS)
    design_i <- v2a_design(groups_i)
    filtering <- v2a_filter(counts_i, design_i, groups_i)
    y <- edgeR::normLibSizes(filtering$y0[filtering$keep, , keep.lib.sizes = FALSE],
                             method = "TMM")
    y <- edgeR::estimateDisp(y, design_i, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design_i, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, contrast = chronic)
    tab <- data.table::as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table)
    observed <- data.table::data.table(
      gene_id = rownames(y), eligible = TRUE, ineligibility_reason = "ELIGIBLE",
      logFC = tab$logFC, QLF_F = tab$F, df2 = as.numeric(test$df.total),
      PValue = tab$PValue, BH_q = stats::p.adjust(tab$PValue, "BH")
    )
    base <- data.table::data.table(gene_id = full_gene_ids)
    joined <- merge(base, observed, by = "gene_id", all.x = TRUE, sort = FALSE)
    joined[is.na(eligible), `:=`(eligible = FALSE,
                                 ineligibility_reason = "FILTERED_IN_OMISSION_REFIT")]
    joined[, `:=`(
      omitted_sample_id = as.character(ledger$sample_id[[omit]]),
      omitted_biological_unit_id = as.character(ledger$biological_unit_id[[omit]]),
      omitted_group = as.character(ledger$group[[omit]]), retained_sample_n = sum(retained),
      retained_gene_n = sum(filtering$keep), recomputed_filter = TRUE,
      recomputed_TMM = TRUE, recomputed_dispersion = TRUE, recomputed_QL_fit = TRUE,
      contrast = "chronic_combined", ci_status = V2A_EDGE_CI_STATUS,
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )]
    rows[[omit]] <- joined
  }
  result <- data.table::rbindlist(rows, fill = TRUE)
  v2a_assert(nrow(result) == length(full_gene_ids) * ncol(counts) &&
               !anyDuplicated(result[, .(gene_id, omitted_sample_id)]), "GENE_LOAO_COVERAGE")
  result
}

v2a_program_loao <- function(counts, ledger, symbols, members, minimum_members = 10L) {
  registry <- unique(members[, .(collection, program_id, program_name)])
  data.table::setorder(registry, collection, program_id)
  rows <- vector("list", ncol(counts))
  chronic <- matrix(V2A_CONTRASTS$chronic_combined, ncol = 1L,
                    dimnames = list(V2A_GROUPS, "chronic_combined"))
  for (omit in seq_len(ncol(counts))) {
    retained <- seq_len(ncol(counts)) != omit
    counts_i <- counts[, retained, drop = FALSE]
    ledger_i <- ledger[retained]
    groups_i <- factor(ledger_i$group, levels = V2A_GROUPS)
    design_i <- v2a_design(groups_i)
    filtering <- v2a_filter(counts_i, design_i, groups_i)
    y <- edgeR::normLibSizes(filtering$y0[filtering$keep, , keep.lib.sizes = FALSE],
                             method = "TMM")
    y <- edgeR::estimateDisp(y, design_i, robust = TRUE)
    voom_i <- limma::voom(y, design_i, plot = FALSE)
    expression <- v2a_collapse_symbol(voom_i$E, symbols[filtering$keep])
    z_expression <- v2a_scale_rows(expression)
    fold <- data.table::rbindlist(lapply(seq_len(nrow(registry)), function(index) {
      row <- registry[index]
      subset <- members[collection == row$collection & program_id == row$program_id &
                          gene_symbol %in% rownames(z_expression)]
      subset <- unique(subset, by = "gene_symbol")
      eligible <- nrow(subset) >= minimum_members && nrow(subset) <= 500L
      reason <- if (nrow(subset) < minimum_members) "TOO_FEW_RETAINED_MEMBER_GENES" else
        if (nrow(subset) > 500L) "TOO_MANY_RETAINED_MEMBER_GENES" else "ELIGIBLE"
      effect <- se <- ci_low <- ci_high <- p <- NA_real_
      if (eligible) {
        score <- colSums(z_expression[subset$gene_symbol, , drop = FALSE] * subset$weight) /
          sum(abs(subset$weight))
        if (!all(is.finite(score)) || !is.finite(stats::sd(score)) || stats::sd(score) == 0) {
          eligible <- FALSE
          reason <- "NONFINITE_OR_ZERO_SCORE_VARIANCE"
        } else {
          score <- matrix(as.numeric(scale(score)), nrow = 1L)
          fit <- limma::lmFit(score, design_i)
          fit <- limma::eBayes(limma::contrasts.fit(fit, chronic), robust = TRUE)
          effect <- fit$coefficients[[1L]]
          se <- fit$stdev.unscaled[[1L]] * sqrt(fit$s2.post[[1L]])
          critical <- stats::qt(0.975, df = fit$df.total[[1L]])
          ci_low <- effect - critical * se
          ci_high <- effect + critical * se
          p <- fit$p.value[[1L]]
        }
      }
      data.table::data.table(
        collection = row$collection, program_id = row$program_id,
        program_name = row$program_name, eligible = eligible,
        ineligibility_reason = reason, tested_feature_count = nrow(subset),
        effect = effect, SE = se, ci_low = ci_low, ci_high = ci_high, PValue = p,
        omitted_sample_id = as.character(ledger$sample_id[[omit]]),
        omitted_biological_unit_id = as.character(ledger$biological_unit_id[[omit]]),
        omitted_group = as.character(ledger$group[[omit]]), retained_sample_n = sum(retained),
        retained_gene_n = sum(filtering$keep), recomputed_filter = TRUE,
        recomputed_TMM = TRUE, recomputed_dispersion = TRUE, recomputed_voom = TRUE,
        recomputed_gene_z_scores = TRUE, recomputed_activity_scores = TRUE,
        recomputed_limma_fit = eligible, contrast = "chronic_combined",
        formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
        s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
      )
    }), fill = TRUE)
    fold[, BH_q := stats::ave(PValue, collection, FUN = function(x) stats::p.adjust(x, "BH"))]
    rows[[omit]] <- fold
  }
  result <- data.table::rbindlist(rows, fill = TRUE)
  v2a_assert(nrow(result) == nrow(registry) * ncol(counts) &&
               !anyDuplicated(result[, .(collection, program_id, omitted_sample_id)]),
             "PROGRAM_LOAO_COVERAGE")
  result
}

v2a_diagnostics <- function(counts, ledger, filtering, y, ql_fit, voom, design) {
  raw_sizes <- colSums(counts)
  retained_sizes <- colSums(counts[filtering$keep, , drop = FALSE])
  library <- data.table::data.table(
    sample_id = colnames(counts), biological_unit_id = ledger$biological_unit_id,
    group = ledger$group, raw_library_size = raw_sizes,
    retained_library_size = retained_sizes, TMM_norm_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors
  )
  mds <- limma::plotMDS(y, top = min(500L, nrow(y)), plot = FALSE)
  composition <- data.table::data.table(
    sample_id = colnames(counts), group = ledger$group,
    MDS1 = mds$x, MDS2 = mds$y, distance_method = "leading_logFC"
  )
  dispersion <- data.table::data.table(
    gene_id = rownames(y), AveLogCPM = ql_fit$AveLogCPM,
    tagwise_dispersion = y$tagwise.dispersion,
    BCV = sqrt(y$tagwise.dispersion), QL_s2_post = ql_fit$s2.post,
    QL_df_prior = ql_fit$df.prior, QL_df_residual = ql_fit$df.residual
  )
  mean_variance <- data.table::data.table(
    gene_id = rownames(y), voom_Amean = voom$Amean,
    mean_voom_weight = rowMeans(voom$weights), QL_s2_post = ql_fit$s2.post
  )
  fit0 <- limma::lmFit(voom, design)
  fitted <- fit0$coefficients %*% t(design)
  residual <- voom$E - fitted
  leverage <- diag(design %*% solve(crossprod(design), t(design)))
  influence <- data.table::data.table(
    sample_id = colnames(counts), group = ledger$group, design_leverage = leverage,
    residual_RMSE = sqrt(colMeans(residual^2)), max_abs_residual = apply(abs(residual), 2L, max)
  )
  prefilter <- data.table::data.table(
    input_feature_count = nrow(counts), all_zero_feature_count = sum(rowSums(counts) == 0),
    filterByExpr_pass_count = filtering$by_expr_n,
    CPM_rule_pass_count = filtering$cpm_n, retained_feature_count = sum(filtering$keep),
    smallest_group_n = filtering$min_group_n
  )
  list(library = library, composition = composition, dispersion = dispersion,
       mean_variance = mean_variance, influence = influence, prefilter = prefilter)
}

v2a_write_table <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(table, path, sep = "\t", quote = FALSE, na = "NA")
}

v2a_fit_branch <- function(input, selected, branch, members, root, nrot = 9999L,
                           minimum_members = 10L) {
  counts <- input$counts[, selected, drop = FALSE]
  ledger <- input$ledger[selected]
  symbols <- input$symbols
  groups <- factor(ledger$group, levels = V2A_GROUPS)
  expected <- if (branch == "primary_all30") rep(5L, 6L) else c(5L, 4L, 5L, 5L, 5L, 5L)
  v2a_assert(identical(as.integer(table(groups)), expected), "BRANCH_GROUP_COUNTS:%s", branch)
  design <- v2a_design(groups)
  contrast_matrix <- v2a_contrast_matrix()
  filtering <- v2a_filter(counts, design, groups)
  y <- edgeR::normLibSizes(filtering$y0[filtering$keep, , keep.lib.sizes = FALSE],
                           method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  ql_fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  voom <- limma::voom(y, design, plot = FALSE)
  gene_ids <- rownames(y)
  edge <- v2a_edge_tables(ql_fit, contrast_matrix, gene_ids)
  limma <- v2a_limma_tables(voom, design, contrast_matrix, gene_ids)
  deseq <- v2a_deseq_tables(counts, ledger, filtering$keep, contrast_matrix)
  diagnostics <- v2a_diagnostics(counts, ledger, filtering, y, ql_fit, voom, design)
  activity <- v2a_activity_tables(voom, symbols[filtering$keep], design, contrast_matrix,
                                  members, minimum_members)
  camera <- v2a_camera_tables(voom, symbols[filtering$keep], design, contrast_matrix,
                              members, minimum_members)
  chronic <- edge[contrast == "chronic_combined"]
  sensitivities <- v2a_chronic_program_sensitivities(
    voom, symbols[filtering$keep], design, contrast_matrix, members, chronic,
    nrot = nrot, minimum_members = minimum_members
  )
  gene_loao <- v2a_gene_loao(counts, ledger, rownames(counts))
  program_loao <- v2a_program_loao(counts, ledger, symbols, members, minimum_members)
  branch_root <- file.path(root, branch)
  outputs <- list(
    "gene_level/edgeR_primary.tsv.gz" = edge,
    "gene_level/limma_voom_sensitivity.tsv.gz" = limma,
    "gene_level/DESeq2_Wald_sensitivity.tsv.gz" = deseq,
    "diagnostics/library_sizes_TMM.tsv" = diagnostics$library,
    "diagnostics/MDS_coordinates.tsv" = diagnostics$composition,
    "diagnostics/dispersion_BCV_QL.tsv.gz" = diagnostics$dispersion,
    "diagnostics/mean_variance.tsv.gz" = diagnostics$mean_variance,
    "diagnostics/sample_influence.tsv" = diagnostics$influence,
    "diagnostics/prefilter_summary.tsv" = diagnostics$prefilter,
    "programs/camera.tsv.gz" = camera,
    "programs/mroast_chronic.tsv.gz" = sensitivities$mroast,
    "programs/fgsea_chronic.tsv.gz" = sensitivities$fgsea,
    "programs/activity_all_contrasts.tsv.gz" = activity,
    "loao/gene_chronic_true_refit.tsv.gz" = gene_loao,
    "loao/program_chronic_true_refit.tsv.gz" = program_loao
  )
  for (relative in names(outputs)) v2a_write_table(outputs[[relative]], file.path(branch_root, relative))
  list(edge = edge, retained = gene_ids, output_paths = file.path(branch, names(outputs)))
}

v2a_compare_branches <- function(primary, complete29, root) {
  merged <- merge(
    primary$edge[, .(gene_id, contrast, logFC_all30 = logFC, q_all30 = BH_q)],
    complete29$edge[, .(gene_id, contrast, logFC_complete29 = logFC, q_complete29 = BH_q)],
    by = c("gene_id", "contrast"), all = TRUE
  )
  summary <- merged[, .(
    all30_tested = sum(is.finite(logFC_all30)), complete29_tested = sum(is.finite(logFC_complete29)),
    common_tested = sum(is.finite(logFC_all30) & is.finite(logFC_complete29)),
    sign_concordance = mean(sign(logFC_all30) == sign(logFC_complete29), na.rm = TRUE),
    spearman_logFC = suppressWarnings(stats::cor(logFC_all30, logFC_complete29,
                                                  method = "spearman", use = "complete.obs"))
  ), by = contrast]
  v2a_write_table(merged, file.path(root, "comparison/all30_vs_complete29_gene.tsv.gz"))
  v2a_write_table(summary, file.path(root, "comparison/all30_vs_complete29_summary.tsv"))
  c("comparison/all30_vs_complete29_gene.tsv.gz", "comparison/all30_vs_complete29_summary.tsv")
}

v2a_execute_model <- function(input, gene_set_path, staging, nrot = 9999L,
                              minimum_members = 10L) {
  v2a_assert(!file.exists(staging), "STAGING_ALREADY_EXISTS")
  v2a_assert(dir.create(staging, recursive = FALSE, showWarnings = FALSE),
             "STAGING_CREATE_FAILED")
  resources <- readRDS(gene_set_path)
  gene_set_identity <- v2a_file_identity(gene_set_path)
  members <- v2a_program_members(resources)
  primary_selected <- rep(TRUE, nrow(input$ledger))
  sensitivity_selected <- input$ledger$sample_id != "GSM3681230"
  v2a_assert(sum(sensitivity_selected) == 29L &&
               identical(input$ledger$run_accession[!sensitivity_selected], "SRR8757507"),
             "COMPLETE29_EXCLUSION")
  primary <- v2a_fit_branch(input, primary_selected, "primary_all30", members, staging,
                            nrot, minimum_members)
  complete29 <- v2a_fit_branch(input, sensitivity_selected, "complete29_sensitivity",
                               members, staging, nrot, minimum_members)
  comparison <- v2a_compare_branches(primary, complete29, staging)
  marker <- data.table::data.table(
    candidate_id = V2A_CANDIDATE, status = "ANALYZER_COMPLETE_PENDING_INDEPENDENT_VALIDATION",
    primary_sample_n = 30L, complete29_sample_n = 29L,
    excluded_sample_id = "GSM3681230", excluded_run_accession = "SRR8757507",
    gene_set_bytes = gene_set_identity$bytes,
    gene_set_sha256 = gene_set_identity$sha256,
    enrichment_seed = V2A_ENRICHMENT_SEED,
    rng_kind = V2A_RNG_KIND[[1L]], rng_normal_kind = V2A_RNG_KIND[[2L]],
    rng_sample_kind = V2A_RNG_KIND[[3L]], mroast_rotation_count = nrot,
    edgeR_ci_status = V2A_EDGE_CI_STATUS, formal_sequential_inference = FALSE,
    candidate_generation_allowed = FALSE, s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT",
    batch_adjustment_performed = FALSE, combat_performed = FALSE
  )
  v2a_write_table(marker, file.path(staging, "ANALYZER_COMPLETE.tsv"))
  invisible(c(primary$output_paths, complete29$output_paths, comparison, "ANALYZER_COMPLETE.tsv"))
}

v2a_main <- function() {
  args <- v2a_parse_args(commandArgs(trailingOnly = TRUE))
  v2a_require_runtime(test_only = args$test_only)
  contract <- v2a_read_json(args$contract, "contract")
  v2a_validate_contract_gate(contract, test_only = args$test_only)
  # Even test execution needs a separate synthetic harness; the command-line
  # production path may never use --test-only as authorization.
  v2a_assert(!args$test_only, "TEST_ONLY_NOT_A_PRODUCTION_ENTRY")
  v2a_validate_future_pin_set(contract)
  v2a_stop("V2_DESCRIPTOR_EXECUTION_REQUIRES_A_FUTURE_REVIEWED_CONTRACT_SUCCESSOR")
}

if (sys.nframe() == 0L) v2a_main()
