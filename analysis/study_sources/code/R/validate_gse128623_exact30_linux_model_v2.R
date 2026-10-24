#!/usr/bin/env Rscript

# Independent reconstruction validator for exact30 model v2. This file does not
# source or import the analyzer or any rejected v4/v7 statistical core.

V2V_CANDIDATE <- "GSE128623_exact30_linux_model_v2"
V2V_GROUPS <- c("control", "D4", "D7", "D14", "D28", "D42")
V2V_CONTRASTS <- cbind(
  chronic_combined = c(-1, 0, 0, 0, 0.5, 0.5),
  D28_vs_control = c(-1, 0, 0, 0, 1, 0),
  D42_vs_control = c(-1, 0, 0, 0, 0, 1),
  acute_combined = c(-1, 0.5, 0.5, 0, 0, 0),
  D4_vs_control = c(-1, 1, 0, 0, 0, 0),
  D7_vs_control = c(-1, 0, 1, 0, 0, 0),
  D14_vs_control = c(-1, 0, 0, 1, 0, 0)
)
rownames(V2V_CONTRASTS) <- V2V_GROUPS
V2V_VERSIONS <- c(
  data.table = "1.18.4", edgeR = "4.10.1", limma = "3.68.4",
  DESeq2 = "1.52.0", fgsea = "1.38.0", jsonlite = "2.0.0",
  yaml = "2.3.12", openssl = "2.4.2"
)
V2V_EDGE_CI <- "NOT_AVAILABLE_FROM_STANDARD_EDGER_QLF_OUTPUT"
V2V_ENRICHMENT_SEED <- 20260819L
V2V_RNG_KIND <- c("Mersenne-Twister", "Inversion", "Rejection")

v2v_stop <- function(...) stop(sprintf(...), call. = FALSE)
v2v_require <- function(x, ...) if (!isTRUE(x)) v2v_stop(...)

v2v_reset_enrichment_rng <- function() {
  RNGkind(kind = V2V_RNG_KIND[[1L]], normal.kind = V2V_RNG_KIND[[2L]],
          sample.kind = V2V_RNG_KIND[[3L]])
  set.seed(V2V_ENRICHMENT_SEED)
  invisible(TRUE)
}

v2v_file_identity <- function(path) {
  v2v_require(file.exists(path) && !dir.exists(path), "IDENTITY_FILE:%s", path)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  list(bytes = unname(as.numeric(file.info(path)$size)),
       sha256 = tolower(as.character(openssl::sha256(connection))))
}

v2v_runtime <- function(test_only = FALSE) {
  v2v_require(as.character(getRversion()) == "4.6.1", "R_VERSION")
  for (p in names(V2V_VERSIONS)) {
    v2v_require(requireNamespace(p, quietly = TRUE), "PACKAGE_MISSING:%s", p)
    v2v_require(as.character(utils::packageVersion(p)) == V2V_VERSIONS[[p]],
                "PACKAGE_VERSION:%s", p)
  }
  if (!test_only) v2v_require(identical(Sys.info()[["sysname"]], "Linux"), "LINUX_ONLY")
}

v2v_args <- function(x) {
  out <- list(contract = NULL, descriptor = NULL, staging = NULL, test_only = FALSE)
  for (value in x) {
    if (startsWith(value, "--contract=")) out$contract <- substring(value, 12L)
    else if (startsWith(value, "--descriptor=")) out$descriptor <- substring(value, 14L)
    else if (startsWith(value, "--staging=")) out$staging <- substring(value, 11L)
    else if (identical(value, "--test-only")) out$test_only <- TRUE
    else v2v_stop("UNKNOWN_ARGUMENT:%s", value)
  }
  v2v_require(all(vapply(out[1:3], function(z) !is.null(z) && nzchar(z), logical(1))),
                "ARGUMENT_MISSING")
  out
}

v2v_json <- function(path) {
  v2v_require(file.exists(path) && !dir.exists(path), "JSON_MISSING:%s", path)
  value <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  v2v_require(is.list(value) && !is.null(names(value)), "JSON_OBJECT:%s", path)
  value
}

v2v_contract_gate <- function(contract, test_only = FALSE) {
  v2v_require(identical(contract$schema,
                        "route-b2-gse128623-exact30-linux-model-v2-contract-v1"),
                "CONTRACT_SCHEMA")
  v2v_require(identical(contract$candidate_id, V2V_CANDIDATE), "CANDIDATE")
  v2v_require(identical(contract$claim_boundary$candidate_generation_allowed, FALSE) &&
                  identical(contract$claim_boundary$formal_sequential_inference, FALSE) &&
                  identical(contract$claim_boundary$s2_status,
                            "NOT_IDENTIFIABLE_NO_S1_INPUT"), "CLAIM_BOUNDARY")
  v2v_require(identical(contract$batch_policy$batch_adjustment_allowed, FALSE) &&
                  identical(contract$batch_policy$combat_allowed, FALSE), "BATCH_BOUNDARY")
  if (!test_only) {
    v2v_require(identical(contract$production_go, TRUE), "PRODUCTION_NO_GO")
    v2v_require(identical(contract$methodological_amendment$effective, TRUE),
                  "CI_AMENDMENT_NOT_EFFECTIVE")
  }
}

v2v_read_input <- function(count_path, annotation_path, ledger_path,
                           count_id_col, annotation_id_col, symbol_col) {
  ledger <- data.table::fread(ledger_path, sep = "\t", header = TRUE, colClasses = "character")
  required <- c("sample_order", "sample_id", "biological_unit_id", "biosample",
                "sra_experiment", "run_accession", "sex", "group")
  v2v_require(identical(names(ledger), required) && nrow(ledger) == 30L,
                "LEDGER_SCHEMA_OR_N")
  v2v_require(identical(as.integer(ledger$sample_order), 1:30), "LEDGER_ORDER")
  for (field in required[2:6]) {
    v2v_require(all(nzchar(ledger[[field]])) && !anyDuplicated(ledger[[field]]),
                  "LEDGER_IDENTITY:%s", field)
  }
  v2v_require(all(ledger$sex == "female"), "SEX")
  v2v_require(identical(as.integer(table(factor(ledger$group, levels = V2V_GROUPS))),
                          rep(5L, 6L)), "GROUP_COUNTS")
  count_table <- data.table::fread(count_path, sep = "\t", header = TRUE,
                                   na.strings = c("NA", "NaN", ""), check.names = FALSE)
  v2v_require(identical(names(count_table)[1L], count_id_col), "COUNT_ID_COLUMN")
  v2v_require(identical(names(count_table)[-1L], ledger$sample_id), "COUNT_SAMPLE_ORDER")
  ids <- trimws(as.character(count_table[[1L]]))
  v2v_require(all(nzchar(ids)) && !anyDuplicated(ids), "COUNT_IDS")
  counts <- as.matrix(count_table[, -1L, with = FALSE])
  storage.mode(counts) <- "double"
  v2v_require(all(is.finite(counts)) && all(counts >= 0) &&
                  all(abs(counts - round(counts)) <= .Machine$double.eps^0.5), "COUNT_VALUES")
  counts <- round(counts)
  rownames(counts) <- ids
  colnames(counts) <- ledger$sample_id
  v2v_require(all(colSums(counts) > 0) && sum(counts) > 0, "ZERO_COUNTS")
  annotation <- data.table::fread(annotation_path, sep = "\t", header = TRUE,
                                  colClasses = "character", na.strings = c("NA", ""))
  v2v_require(all(c(annotation_id_col, symbol_col) %in% names(annotation)),
                "ANNOTATION_COLUMNS")
  ann_ids <- trimws(annotation[[annotation_id_col]])
  v2v_require(identical(ann_ids, ids) && !anyDuplicated(ann_ids), "ANNOTATION_ORDER")
  symbols <- trimws(annotation[[symbol_col]])
  symbols[is.na(symbols) | !nzchar(symbols)] <- NA_character_
  list(counts = counts, ledger = ledger, annotation = annotation, symbols = symbols)
}

v2v_s1_zero <- function(path) {
  expected <- c("collection", "feature_id", "feature_name", "S1_pass", "direction",
                "primary_effect", "primary_ci_low", "primary_ci_high", "primary_q_value",
                "permutation_q_value_BH", "loo_direction_fraction",
                "original_annotation_direction_stable", "all27_direction_stable")
  value <- data.table::fread(path, sep = "\t", header = TRUE)
  v2v_require(identical(names(value), expected) && nrow(value) == 1771L &&
                  !anyDuplicated(value[, .(collection, feature_id)]), "S1_SCHEMA_KEYS")
  pass <- toupper(as.character(value$S1_pass))
  v2v_require(all(pass %in% c("TRUE", "FALSE")) && sum(pass == "TRUE") == 0L,
                "S1_NOT_ZERO")
}

v2v_design <- function(group) {
  factor_group <- factor(group, levels = V2V_GROUPS)
  v2v_require(!anyNA(factor_group), "GROUP_LEVEL")
  design <- stats::model.matrix(~ 0 + factor_group)
  colnames(design) <- V2V_GROUPS
  v2v_require(qr(design)$rank == ncol(design) && identical(colnames(design), V2V_GROUPS),
                "DESIGN")
  design
}

v2v_branch_fit <- function(counts, ledger) {
  group <- factor(ledger$group, levels = V2V_GROUPS)
  design <- v2v_design(group)
  y0 <- edgeR::DGEList(counts = counts)
  by_expr <- edgeR::filterByExpr(y0, design = design)
  minimum_group_n <- min(table(group))
  cpm_rule <- rowSums(edgeR::cpm(y0, log = FALSE) >= 1) >= minimum_group_n
  keep <- by_expr & cpm_rule
  v2v_require(sum(keep) >= 2L, "FILTER_RETAINED")
  y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  ql <- edgeR::glmQLFit(y, design, robust = TRUE)
  voom <- limma::voom(y, design, plot = FALSE)
  list(group = group, design = design, y0 = y0, by_expr = by_expr,
       cpm_rule = cpm_rule, minimum_group_n = minimum_group_n, keep = keep,
       y = y, ql = ql, voom = voom)
}

v2v_edge <- function(fit) {
  data.table::rbindlist(lapply(seq_len(ncol(V2V_CONTRASTS)), function(j) {
    test <- edgeR::glmQLFTest(fit$ql, contrast = V2V_CONTRASTS[, j])
    table <- data.table::as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table)
    data.table::data.table(
      gene_id = rownames(fit$y), contrast = colnames(V2V_CONTRASTS)[[j]],
      logFC = table$logFC, logCPM = table$logCPM, QLF_F = table$F,
      df1 = 1, df2 = as.numeric(test$df.total), PValue = table$PValue,
      BH_q = stats::p.adjust(table$PValue, "BH"),
      tested_feature_count = nrow(fit$y), ci_low = NA_real_, ci_high = NA_real_,
      ci_status = V2V_EDGE_CI, method = "edgeR_TMM_robust_QL_primary",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2v_limma <- function(fit) {
  model <- limma::lmFit(fit$voom, fit$design)
  model <- limma::eBayes(limma::contrasts.fit(model, V2V_CONTRASTS), robust = TRUE)
  data.table::rbindlist(lapply(seq_len(ncol(V2V_CONTRASTS)), function(j) {
    se <- model$stdev.unscaled[, j] * sqrt(model$s2.post)
    effect <- model$coefficients[, j]
    critical <- stats::qt(0.975, model$df.total)
    p <- model$p.value[, j]
    data.table::data.table(
      gene_id = rownames(fit$y), contrast = colnames(V2V_CONTRASTS)[[j]],
      logFC = effect, SE = se, df = model$df.total, t = model$t[, j],
      ci_low = effect - critical * se, ci_high = effect + critical * se,
      PValue = p, BH_q = stats::p.adjust(p, "BH"),
      tested_feature_count = nrow(fit$y), method = "limma_voom_robust_sensitivity",
      ci_status = "METHOD_SPECIFIC_LIMMA_VOOM_95_PERCENT",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2v_deseq <- function(counts, ledger, keep) {
  col_data <- data.frame(group = factor(ledger$group, levels = V2V_GROUPS),
                         row.names = ledger$sample_id)
  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts[keep, , drop = FALSE]), col_data,
                                        design = ~ 0 + group)
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", betaPrior = FALSE,
                       minReplicatesForReplace = Inf, quiet = TRUE)
  v2v_require(identical(DESeq2::resultsNames(dds), paste0("group", V2V_GROUPS)),
                "DESEQ_NAMES")
  data.table::rbindlist(lapply(seq_len(ncol(V2V_CONTRASTS)), function(j) {
    result <- DESeq2::results(dds, contrast = V2V_CONTRASTS[, j],
                              independentFiltering = FALSE, cooksCutoff = FALSE)
    effect <- as.numeric(result$log2FoldChange)
    se <- as.numeric(result$lfcSE)
    p <- as.numeric(result$pvalue)
    data.table::data.table(
      gene_id = rownames(result), contrast = colnames(V2V_CONTRASTS)[[j]],
      logFC = effect, lfcSE = se, WaldStatistic = as.numeric(result$stat),
      ci_low = effect - stats::qnorm(0.975) * se,
      ci_high = effect + stats::qnorm(0.975) * se,
      PValue = p, BH_q = stats::p.adjust(p, "BH"),
      tested_feature_count = nrow(result), method = "DESeq2_Wald_sensitivity",
      ci_status = "METHOD_SPECIFIC_DESEQ2_WALD_95_PERCENT",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2v_numeric_equal <- function(observed, expected, tolerance = 2e-7) {
  if (length(observed) != length(expected)) return(FALSE)
  both_na <- is.na(observed) & is.na(expected)
  finite <- is.finite(observed) & is.finite(expected)
  all(both_na | (finite & abs(observed - expected) <= tolerance *
                   pmax(1, abs(observed), abs(expected))))
}

v2v_compare_table <- function(path, expected, keys, numeric, exact, label) {
  observed <- data.table::fread(path, sep = "\t", header = TRUE,
                                na.strings = c("NA", "NaN", ""))
  v2v_require(identical(names(observed), names(expected)), "TABLE_SCHEMA:%s", label)
  v2v_require(nrow(observed) == nrow(expected) && !anyDuplicated(observed[, ..keys]),
                "TABLE_KEYS:%s", label)
  data.table::setorderv(observed, keys)
  data.table::setorderv(expected, keys)
  for (field in exact) {
    v2v_require(identical(as.character(observed[[field]]), as.character(expected[[field]])),
                  "TABLE_EXACT:%s:%s", label, field)
  }
  for (field in numeric) {
    v2v_require(v2v_numeric_equal(as.numeric(observed[[field]]), as.numeric(expected[[field]])),
                  "TABLE_NUMERIC:%s:%s", label, field)
  }
  invisible(TRUE)
}

v2v_compare_reconstruction <- function(path, expected, keys, label) {
  numeric <- names(expected)[vapply(expected, is.numeric, logical(1L))]
  exact <- setdiff(names(expected), numeric)
  v2v_compare_table(path, expected, keys, numeric, exact, label)
}

v2v_members <- function(resources) {
  need <- c("reactome", "syngo_sensitivity", "ortholog_one_to_one",
            "progeny_mouse", "dorothea_mouse")
  v2v_require(all(need %in% names(resources)), "GENE_SET_SCHEMA")
  ortho <- data.table::as.data.table(resources$ortholog_one_to_one)
  mapping <- unique(ortho[other_tax_id == 10090L & !is.na(human_symbol) &
                            !is.na(other_symbol),
                          .(human_symbol = as.character(human_symbol),
                            gene_symbol = as.character(other_symbol))])
  human <- function(tab, id, name, collection_value) {
    tab <- data.table::as.data.table(tab)
    left <- data.table::data.table(program_id = as.character(tab[[id]]),
                                   program_name = as.character(tab[[name]]),
                                   human_symbol = as.character(tab$gene_symbol))
    mapped <- merge(left, mapping, by = "human_symbol", allow.cartesian = TRUE)
    unique(mapped[, .(collection = collection_value, program_id, program_name,
                       gene_symbol, weight = 1, test_class = "unsigned_and_activity")])
  }
  all <- data.table::rbindlist(list(
    human(resources$reactome, "term_id", "term_name", "Reactome"),
    human(resources$syngo_sensitivity, "term_id", "term_name",
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
  unique(all[nzchar(collection) & nzchar(program_id) & nzchar(program_name) &
               nzchar(gene_symbol) & is.finite(weight)],
         by = c("collection", "program_id", "gene_symbol"))
}

v2v_symbol_expression <- function(voom, symbols) {
  valid <- !is.na(symbols) & nzchar(symbols)
  sums <- rowsum(voom$E[valid, , drop = FALSE], symbols[valid], reorder = FALSE)
  counts <- as.numeric(table(factor(symbols[valid], levels = rownames(sums))))
  expression <- sums / counts
  center <- rowMeans(expression)
  spread <- apply(expression, 1L, stats::sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  sweep(sweep(expression, 1L, center, "-"), 1L, spread, "/")
}

v2v_activity <- function(fit, symbols, members, minimum_members) {
  z <- v2v_symbol_expression(fit$voom, symbols[fit$keep])
  registry <- unique(members[, .(collection, program_id)])
  data.table::setorder(registry, collection, program_id)
  records <- lapply(seq_len(nrow(registry)), function(i) {
    key <- registry[i]
    block <- unique(members[collection == key$collection & program_id == key$program_id &
                              gene_symbol %in% rownames(z)], by = "gene_symbol")
    if (nrow(block) < minimum_members || nrow(block) > 500L) return(NULL)
    score <- colSums(z[block$gene_symbol, , drop = FALSE] * block$weight) / sum(abs(block$weight))
    if (!all(is.finite(score)) || stats::sd(score) == 0) return(NULL)
    list(key = key, score = as.numeric(scale(score)))
  })
  records <- Filter(Negate(is.null), records)
  if (!length(records)) return(data.table::data.table())
  matrix <- do.call(rbind, lapply(records, `[[`, "score"))
  rownames(matrix) <- vapply(records, function(x) paste(x$key$collection, x$key$program_id,
                                                        sep = "::"), character(1))
  model <- limma::lmFit(matrix, fit$design)
  model <- limma::eBayes(limma::contrasts.fit(model, V2V_CONTRASTS), robust = TRUE)
  data.table::rbindlist(lapply(seq_len(ncol(V2V_CONTRASTS)), function(j) {
    se <- model$stdev.unscaled[, j] * sqrt(model$s2.post)
    effect <- model$coefficients[, j]
    critical <- stats::qt(0.975, model$df.total)
    p <- model$p.value[, j]
    keys <- data.table::tstrsplit(rownames(matrix), "::", fixed = TRUE)
    data.table::data.table(
      collection = keys[[1]], program_id = keys[[2]],
      contrast = colnames(V2V_CONTRASTS)[[j]], effect = effect, SE = se,
      ci_low = effect - critical * se, ci_high = effect + critical * se,
      PValue = p, BH_q = stats::ave(p, keys[[1]], FUN = function(x) stats::p.adjust(x, "BH")),
      tested_program_count = length(effect), method = "weighted_voom_activity_descriptive",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )
  }))
}

v2v_diagnostics <- function(counts, ledger, fit) {
  raw_sizes <- colSums(counts)
  retained_sizes <- colSums(counts[fit$keep, , drop = FALSE])
  library <- data.table::data.table(
    sample_id = colnames(counts), biological_unit_id = ledger$biological_unit_id,
    group = ledger$group, raw_library_size = raw_sizes,
    retained_library_size = retained_sizes, TMM_norm_factor = fit$y$samples$norm.factors,
    effective_library_size = fit$y$samples$lib.size * fit$y$samples$norm.factors
  )
  mds <- limma::plotMDS(fit$y, top = min(500L, nrow(fit$y)), plot = FALSE)
  composition <- data.table::data.table(
    sample_id = colnames(counts), group = ledger$group,
    MDS1 = mds$x, MDS2 = mds$y, distance_method = "leading_logFC"
  )
  dispersion <- data.table::data.table(
    gene_id = rownames(fit$y), AveLogCPM = fit$ql$AveLogCPM,
    tagwise_dispersion = fit$y$tagwise.dispersion,
    BCV = sqrt(fit$y$tagwise.dispersion), QL_s2_post = fit$ql$s2.post,
    QL_df_prior = fit$ql$df.prior, QL_df_residual = fit$ql$df.residual
  )
  mean_variance <- data.table::data.table(
    gene_id = rownames(fit$y), voom_Amean = fit$voom$Amean,
    mean_voom_weight = rowMeans(fit$voom$weights), QL_s2_post = fit$ql$s2.post
  )
  model <- limma::lmFit(fit$voom, fit$design)
  fitted <- model$coefficients %*% t(fit$design)
  residual <- fit$voom$E - fitted
  leverage <- diag(fit$design %*% solve(crossprod(fit$design), t(fit$design)))
  influence <- data.table::data.table(
    sample_id = colnames(counts), group = ledger$group, design_leverage = leverage,
    residual_RMSE = sqrt(colMeans(residual^2)),
    max_abs_residual = apply(abs(residual), 2L, max)
  )
  prefilter <- data.table::data.table(
    input_feature_count = nrow(counts), all_zero_feature_count = sum(rowSums(counts) == 0),
    filterByExpr_pass_count = sum(fit$by_expr), CPM_rule_pass_count = sum(fit$cpm_rule),
    retained_feature_count = sum(fit$keep), smallest_group_n = fit$minimum_group_n
  )
  list(library = library, composition = composition, dispersion = dispersion,
       mean_variance = mean_variance, influence = influence, prefilter = prefilter)
}

v2v_camera <- function(fit, symbols, members, minimum_members) {
  unsigned <- members[grepl("unsigned", test_class)]
  symbol_to_row <- split(seq_along(symbols), symbols)
  rows <- list()
  cursor <- 1L
  for (collection_name in unique(unsigned$collection)) {
    block <- unsigned[collection == collection_name]
    registry <- unique(block[, .(program_id, program_name)])
    index_list <- lapply(registry$program_id, function(id) {
      genes <- unique(block[program_id == id, gene_symbol])
      unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    })
    names(index_list) <- registry$program_id
    index_list <- index_list[lengths(index_list) >= minimum_members & lengths(index_list) <= 500L]
    if (!length(index_list)) next
    for (j in seq_len(ncol(V2V_CONTRASTS))) {
      result <- data.table::as.data.table(limma::camera(
        fit$voom, index = index_list, design = fit$design,
        contrast = V2V_CONTRASTS[, j], inter.gene.cor = 0.01, sort = FALSE
      ), keep.rownames = "program_id")
      result[, `:=`(
        collection = collection_name, contrast = colnames(V2V_CONTRASTS)[[j]],
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
  if (!length(rows)) data.table::data.table() else data.table::rbindlist(rows, fill = TRUE)
}

v2v_ranking <- function(symbols, edge_chronic) {
  table <- data.table::data.table(
    symbol = symbols, gene_id = edge_chronic$gene_id,
    score = sign(edge_chronic$logFC) * sqrt(pmax(edge_chronic$QLF_F, 0))
  )
  table <- table[!is.na(symbol) & nzchar(symbol) & is.finite(score)]
  table[, abs_score := abs(score)]
  data.table::setorder(table, symbol, -abs_score, gene_id)
  table <- unique(table, by = "symbol")
  data.table::setorder(table, -score, symbol)
  ranks <- table$score
  names(ranks) <- table$symbol
  payload <- paste0(paste(sprintf("%s\t%.17g", names(ranks), ranks), collapse = "\n"), "\n")
  list(
    values = ranks,
    rule = "max_abs_signed_sqrt_QLF_F_per_symbol_then_score_desc_symbol_asc",
    sha256 = tolower(as.character(openssl::sha256(charToRaw(payload))))
  )
}

v2v_program_sensitivities <- function(fit, symbols, members, edge_chronic,
                                      nrot, minimum_members) {
  unsigned <- members[grepl("unsigned", test_class)]
  symbol_to_row <- split(seq_along(symbols), symbols)
  rank_profile <- v2v_ranking(symbols, edge_chronic)
  ranks <- rank_profile$values
  roast_rows <- list()
  fgsea_rows <- list()
  cursor <- 1L
  for (collection_name in unique(unsigned$collection)) {
    block <- unsigned[collection == collection_name]
    registry <- unique(block[, .(program_id, program_name)])
    index_list <- lapply(registry$program_id, function(id) {
      genes <- unique(block[program_id == id, gene_symbol])
      unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    })
    names(index_list) <- registry$program_id
    index_list <- index_list[lengths(index_list) >= minimum_members & lengths(index_list) <= 500L]
    if (!length(index_list)) next
    v2v_reset_enrichment_rng()
    roast <- data.table::as.data.table(limma::mroast(
      fit$voom, index = index_list, design = fit$design,
      contrast = V2V_CONTRASTS[, "chronic_combined"], nrot = nrot,
      set.statistic = "mean", midp = TRUE, sort = "none"
    ), keep.rownames = "program_id")
    roast[, `:=`(
      collection = collection_name, contrast = "chronic_combined",
      BH_q = pmax(PValue, stats::p.adjust(pmax(PValue - 0.5 / nrot, 0), "BH")),
      method = sprintf("mroast_mean_%d_rotations", nrot),
      package_version = as.character(utils::packageVersion("limma")),
      rng_seed = V2V_ENRICHMENT_SEED, rng_kind = V2V_RNG_KIND[[1L]],
      rng_normal_kind = V2V_RNG_KIND[[2L]], rng_sample_kind = V2V_RNG_KIND[[3L]],
      rotation_count = nrot, set_statistic = "mean", midp = TRUE,
      sort_rule = "input_program_id_order",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )]
    roast_rows[[cursor]] <- roast
    paths <- lapply(index_list, function(index) unique(symbols[index]))
    paths <- lapply(paths, function(x) intersect(x, names(ranks)))
    paths <- paths[lengths(paths) >= minimum_members & lengths(paths) <= 500L]
    if (length(paths)) {
      v2v_reset_enrichment_rng()
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
        rng_seed = V2V_ENRICHMENT_SEED, rng_kind = V2V_RNG_KIND[[1L]],
        rng_normal_kind = V2V_RNG_KIND[[2L]], rng_sample_kind = V2V_RNG_KIND[[3L]],
        ranking_rule = rank_profile$rule, ranking_count = length(ranks),
        ranking_sha256 = rank_profile$sha256, min_size = minimum_members,
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
    mroast = if (length(roast_rows)) data.table::rbindlist(roast_rows, fill = TRUE) else
      data.table::data.table(),
    fgsea = if (length(fgsea_rows)) data.table::rbindlist(fgsea_rows, fill = TRUE) else
      data.table::data.table()
  )
}

v2v_gene_loao_keys_and_stats <- function(counts, ledger) {
  output <- vector("list", ncol(counts))
  for (omit in seq_len(ncol(counts))) {
    keep_sample <- seq_len(ncol(counts)) != omit
    fit <- v2v_branch_fit(counts[, keep_sample, drop = FALSE], ledger[keep_sample])
    test <- edgeR::glmQLFTest(fit$ql, contrast = V2V_CONTRASTS[, "chronic_combined"])
    table <- data.table::as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table)
    hit <- data.table::data.table(
      gene_id = rownames(fit$y), eligible = TRUE,
      ineligibility_reason = "ELIGIBLE", logFC = table$logFC, QLF_F = table$F,
      df2 = as.numeric(test$df.total), PValue = table$PValue,
      BH_q = stats::p.adjust(table$PValue, "BH")
    )
    all <- merge(data.table::data.table(gene_id = rownames(counts)), hit,
                 by = "gene_id", all.x = TRUE, sort = FALSE)
    all[is.na(eligible), `:=`(eligible = FALSE,
                              ineligibility_reason = "FILTERED_IN_OMISSION_REFIT")]
    all[, `:=`(
      omitted_sample_id = ledger$sample_id[[omit]],
      omitted_biological_unit_id = ledger$biological_unit_id[[omit]],
      omitted_group = ledger$group[[omit]], retained_sample_n = sum(keep_sample),
      retained_gene_n = sum(fit$keep), recomputed_filter = TRUE,
      recomputed_TMM = TRUE, recomputed_dispersion = TRUE, recomputed_QL_fit = TRUE,
      contrast = "chronic_combined", ci_status = V2V_EDGE_CI,
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
    )]
    output[[omit]] <- all
  }
  data.table::rbindlist(output, fill = TRUE)
}

v2v_program_loao_stats <- function(counts, ledger, symbols, members, minimum_members) {
  registry <- unique(members[, .(collection, program_id, program_name)])
  data.table::setorder(registry, collection, program_id)
  output <- vector("list", ncol(counts))
  chronic <- matrix(V2V_CONTRASTS[, "chronic_combined"], ncol = 1L,
                    dimnames = list(V2V_GROUPS, "chronic_combined"))
  for (omit in seq_len(ncol(counts))) {
    retained <- seq_len(ncol(counts)) != omit
    fit <- v2v_branch_fit(counts[, retained, drop = FALSE], ledger[retained])
    z <- v2v_symbol_expression(fit$voom, symbols[fit$keep])
    fold <- data.table::rbindlist(lapply(seq_len(nrow(registry)), function(index) {
      row <- registry[index]
      block <- members[collection == row$collection & program_id == row$program_id &
                         gene_symbol %in% rownames(z)]
      block <- unique(block, by = "gene_symbol")
      eligible <- nrow(block) >= minimum_members && nrow(block) <= 500L
      reason <- if (nrow(block) < minimum_members) "TOO_FEW_RETAINED_MEMBER_GENES" else
        if (nrow(block) > 500L) "TOO_MANY_RETAINED_MEMBER_GENES" else "ELIGIBLE"
      effect <- se <- ci_low <- ci_high <- p <- NA_real_
      if (eligible) {
        score <- colSums(z[block$gene_symbol, , drop = FALSE] * block$weight) /
          sum(abs(block$weight))
        if (!all(is.finite(score)) || !is.finite(stats::sd(score)) || stats::sd(score) == 0) {
          eligible <- FALSE
          reason <- "NONFINITE_OR_ZERO_SCORE_VARIANCE"
        } else {
          score <- matrix(as.numeric(scale(score)), nrow = 1L)
          model <- limma::lmFit(score, fit$design)
          model <- limma::eBayes(limma::contrasts.fit(model, chronic), robust = TRUE)
          effect <- model$coefficients[[1L]]
          se <- model$stdev.unscaled[[1L]] * sqrt(model$s2.post[[1L]])
          critical <- stats::qt(0.975, model$df.total[[1L]])
          ci_low <- effect - critical * se
          ci_high <- effect + critical * se
          p <- model$p.value[[1L]]
        }
      }
      data.table::data.table(
        collection = row$collection, program_id = row$program_id,
        program_name = row$program_name, eligible = eligible,
        ineligibility_reason = reason, tested_feature_count = nrow(block),
        effect = effect, SE = se, ci_low = ci_low, ci_high = ci_high, PValue = p,
        omitted_sample_id = ledger$sample_id[[omit]],
        omitted_biological_unit_id = ledger$biological_unit_id[[omit]],
        omitted_group = ledger$group[[omit]], retained_sample_n = sum(retained),
        retained_gene_n = sum(fit$keep), recomputed_filter = TRUE,
        recomputed_TMM = TRUE, recomputed_dispersion = TRUE, recomputed_voom = TRUE,
        recomputed_gene_z_scores = TRUE, recomputed_activity_scores = TRUE,
        recomputed_limma_fit = eligible, contrast = "chronic_combined",
        formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
        s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
      )
    }), fill = TRUE)
    fold[, BH_q := stats::ave(PValue, collection, FUN = function(x) stats::p.adjust(x, "BH"))]
    output[[omit]] <- fold
  }
  result <- data.table::rbindlist(output, fill = TRUE)
  v2v_require(nrow(result) == nrow(registry) * ncol(counts) &&
                  !anyDuplicated(result[, .(collection, program_id, omitted_sample_id)]),
                "PROGRAM_LOAO_RECONSTRUCTION_COVERAGE")
  result
}

v2v_validate_branch <- function(input, selected, branch, members, root,
                                nrot = 9999L, minimum_members = 10L) {
  counts <- input$counts[, selected, drop = FALSE]
  ledger <- input$ledger[selected]
  fit <- v2v_branch_fit(counts, ledger)
  edge <- v2v_edge(fit)
  limma <- v2v_limma(fit)
  deseq <- v2v_deseq(counts, ledger, fit$keep)
  branch_root <- file.path(root, branch)
  v2v_compare_table(file.path(branch_root, "gene_level/edgeR_primary.tsv.gz"), edge,
                    c("gene_id", "contrast"), c("logFC", "logCPM", "QLF_F", "df1",
                    "df2", "PValue", "BH_q", "tested_feature_count"),
                    c("gene_id", "contrast", "ci_status", "method",
                      "formal_sequential_inference", "candidate_generation_allowed", "s2_status"),
                    paste0(branch, ":edgeR"))
  observed_edge <- data.table::fread(file.path(branch_root, "gene_level/edgeR_primary.tsv.gz"))
  v2v_require(all(is.na(observed_edge$ci_low)) && all(is.na(observed_edge$ci_high)) &&
                  all(observed_edge$ci_status == V2V_EDGE_CI), "EDGER_CI_MUST_BE_ABSENT")
  v2v_compare_table(file.path(branch_root, "gene_level/limma_voom_sensitivity.tsv.gz"), limma,
                    c("gene_id", "contrast"), c("logFC", "SE", "df", "t", "ci_low",
                    "ci_high", "PValue", "BH_q", "tested_feature_count"),
                    c("gene_id", "contrast", "method", "ci_status",
                      "formal_sequential_inference", "candidate_generation_allowed", "s2_status"),
                    paste0(branch, ":limma"))
  v2v_compare_table(file.path(branch_root, "gene_level/DESeq2_Wald_sensitivity.tsv.gz"), deseq,
                    c("gene_id", "contrast"), c("logFC", "lfcSE", "WaldStatistic", "ci_low",
                    "ci_high", "PValue", "BH_q", "tested_feature_count"),
                    c("gene_id", "contrast", "method", "ci_status",
                      "formal_sequential_inference", "candidate_generation_allowed", "s2_status"),
                    paste0(branch, ":DESeq2"))
  activity <- v2v_activity(fit, input$symbols, members, minimum_members)
  v2v_compare_table(file.path(branch_root, "programs/activity_all_contrasts.tsv.gz"), activity,
                    c("collection", "program_id", "contrast"),
                    c("effect", "SE", "ci_low", "ci_high", "PValue", "BH_q",
                      "tested_program_count"),
                    c("collection", "program_id", "contrast", "method",
                      "formal_sequential_inference", "candidate_generation_allowed", "s2_status"),
                    paste0(branch, ":activity"))
  diagnostics <- v2v_diagnostics(counts, ledger, fit)
  diagnostic_files <- c(
    library = "diagnostics/library_sizes_TMM.tsv",
    composition = "diagnostics/MDS_coordinates.tsv",
    dispersion = "diagnostics/dispersion_BCV_QL.tsv.gz",
    mean_variance = "diagnostics/mean_variance.tsv.gz",
    influence = "diagnostics/sample_influence.tsv",
    prefilter = "diagnostics/prefilter_summary.tsv"
  )
  diagnostic_keys <- list(library = "sample_id", composition = "sample_id",
                          dispersion = "gene_id", mean_variance = "gene_id",
                          influence = "sample_id", prefilter = "input_feature_count")
  for (name in names(diagnostic_files)) {
    v2v_compare_reconstruction(file.path(branch_root, diagnostic_files[[name]]),
                               diagnostics[[name]], diagnostic_keys[[name]],
                               paste0(branch, ":diagnostic:", name))
  }
  camera <- v2v_camera(fit, input$symbols[fit$keep], members, minimum_members)
  v2v_require(nrow(camera) > 0L, "CAMERA_EMPTY:%s", branch)
  v2v_compare_reconstruction(file.path(branch_root, "programs/camera.tsv.gz"), camera,
                             c("collection", "program_id", "contrast"),
                             paste0(branch, ":camera"))
  sensitivities <- v2v_program_sensitivities(
    fit, input$symbols[fit$keep], members, edge[contrast == "chronic_combined"],
    nrot, minimum_members
  )
  v2v_require(nrow(sensitivities$mroast) > 0L && nrow(sensitivities$fgsea) > 0L,
                "PROGRAM_SENSITIVITY_EMPTY:%s", branch)
  v2v_compare_reconstruction(file.path(branch_root, "programs/mroast_chronic.tsv.gz"),
                             sensitivities$mroast,
                             c("collection", "program_id", "contrast"),
                             paste0(branch, ":mroast"))
  v2v_compare_reconstruction(file.path(branch_root, "programs/fgsea_chronic.tsv.gz"),
                             sensitivities$fgsea,
                             c("collection", "program_id", "contrast"),
                             paste0(branch, ":fgsea"))
  gene_loao <- v2v_gene_loao_keys_and_stats(counts, ledger)
  v2v_compare_table(file.path(branch_root, "loao/gene_chronic_true_refit.tsv.gz"), gene_loao,
                    c("gene_id", "omitted_sample_id"),
                    c("logFC", "QLF_F", "df2", "PValue", "BH_q", "retained_sample_n",
                      "retained_gene_n"),
                    c("gene_id", "omitted_sample_id", "omitted_biological_unit_id",
                      "omitted_group", "eligible", "ineligibility_reason", "recomputed_filter",
                      "recomputed_TMM", "recomputed_dispersion", "recomputed_QL_fit",
                      "contrast", "ci_status", "formal_sequential_inference",
                      "candidate_generation_allowed", "s2_status"), paste0(branch, ":gene_LOAO"))
  program_loao <- v2v_program_loao_stats(counts, ledger, input$symbols, members,
                                         minimum_members)
  v2v_compare_reconstruction(file.path(branch_root,
                                       "loao/program_chronic_true_refit.tsv.gz"),
                             program_loao,
                             c("collection", "program_id", "omitted_sample_id"),
                             paste0(branch, ":program_LOAO"))
  list(edge = edge, retained = rownames(fit$y))
}

v2v_validate_output <- function(input, gene_set_path, staging, nrot = 9999L,
                                minimum_members = 10L) {
  required_marker <- file.path(staging, "ANALYZER_COMPLETE.tsv")
  v2v_require(file.exists(required_marker), "ANALYZER_MARKER_MISSING")
  marker <- data.table::fread(required_marker)
  gene_set_identity <- v2v_file_identity(gene_set_path)
  v2v_require(nrow(marker) == 1L && marker$candidate_id[[1L]] == V2V_CANDIDATE &&
                  marker$status[[1L]] == "ANALYZER_COMPLETE_PENDING_INDEPENDENT_VALIDATION" &&
                  as.numeric(marker$gene_set_bytes[[1L]]) == gene_set_identity$bytes &&
                  marker$gene_set_sha256[[1L]] == gene_set_identity$sha256 &&
                  as.integer(marker$enrichment_seed[[1L]]) == V2V_ENRICHMENT_SEED &&
                  marker$rng_kind[[1L]] == V2V_RNG_KIND[[1L]] &&
                  marker$rng_normal_kind[[1L]] == V2V_RNG_KIND[[2L]] &&
                  marker$rng_sample_kind[[1L]] == V2V_RNG_KIND[[3L]] &&
                  as.integer(marker$mroast_rotation_count[[1L]]) == nrot &&
                  marker$edgeR_ci_status[[1L]] == V2V_EDGE_CI &&
                  !as.logical(marker$candidate_generation_allowed[[1L]]) &&
                  !as.logical(marker$formal_sequential_inference[[1L]]) &&
                  !as.logical(marker$batch_adjustment_performed[[1L]]) &&
                  !as.logical(marker$combat_performed[[1L]]), "ANALYZER_MARKER")
  resources <- readRDS(gene_set_path)
  members <- v2v_members(resources)
  primary <- v2v_validate_branch(input, rep(TRUE, 30L), "primary_all30", members,
                                 staging, nrot, minimum_members)
  selected <- input$ledger$sample_id != "GSM3681230"
  v2v_require(sum(selected) == 29L && input$ledger$run_accession[!selected] == "SRR8757507",
                "COMPLETE29_EXCLUSION")
  sensitivity <- v2v_validate_branch(input, selected, "complete29_sensitivity", members,
                                     staging, nrot, minimum_members)
  observed_comparison <- data.table::fread(file.path(staging,
                                                     "comparison/all30_vs_complete29_gene.tsv.gz"))
  expected_keys <- merge(primary$edge[, .(gene_id, contrast)],
                         sensitivity$edge[, .(gene_id, contrast)],
                         by = c("gene_id", "contrast"), all = TRUE)
  v2v_require(nrow(observed_comparison) == nrow(expected_keys) &&
                  !anyDuplicated(observed_comparison[, .(gene_id, contrast)]) &&
                  setequal(paste(observed_comparison$gene_id, observed_comparison$contrast),
                           paste(expected_keys$gene_id, expected_keys$contrast)),
                "BRANCH_COMPARISON_KEYS")
  report <- data.table::data.table(
    candidate_id = V2V_CANDIDATE, status = "PASS_INDEPENDENT_RECONSTRUCTION",
    primary_edge_key_count = nrow(primary$edge),
    complete29_edge_key_count = nrow(sensitivity$edge),
    primary_retained_feature_count = length(primary$retained),
    complete29_retained_feature_count = length(sensitivity$retained),
    edgeR_ci_status = V2V_EDGE_CI, batch_adjustment_performed = FALSE,
    combat_performed = FALSE, formal_sequential_inference = FALSE,
    candidate_generation_allowed = FALSE, s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
  )
  data.table::fwrite(report, file.path(staging, "VALIDATION_REPORT.tsv"), sep = "\t",
                     quote = FALSE, na = "NA")
  invisible(report)
}

v2v_main <- function() {
  args <- v2v_args(commandArgs(trailingOnly = TRUE))
  v2v_runtime(args$test_only)
  contract <- v2v_json(args$contract)
  v2v_contract_gate(contract, args$test_only)
  v2v_require(!args$test_only, "TEST_ONLY_NOT_A_PRODUCTION_ENTRY")
  v2v_stop("V2_VALIDATOR_REQUIRES_FUTURE_REVIEWED_CONTRACT_SUCCESSOR")
}

if (sys.nframe() == 0L) v2v_main()
