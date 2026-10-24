#!/usr/bin/env Rscript

# Canonical JGW03 perturbation analysis for Route B2
#
# Datasets are deliberately analyzed separately:
#   * GSE289929, rat PC12 (primary neural-like perturbation)
#   * GSE227494, mouse C2C12 (muscle-context calibration)
#   * GSE299063, human HepG2 batch A (off-tissue calibration)
#
# The script fails closed when the frozen GSM ledger or accession-specific
# vehicle sets disagree with config/contrasts.yaml.  It never borrows controls
# across accessions.  Its outputs support transcriptional-response and context
# calibration statements only; they cannot establish efficacy, target
# engagement, direct binding, or a causal mechanism.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(decoupleR)
  library(fgsea)
  library(yaml)
  library(ggplot2)
  library(pheatmap)
})

options(stringsAsFactors = FALSE, width = 140)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  normalizePath(file.path(getwd(), "route_b2"), mustWork = TRUE)
candidate_file <- if (length(args) >= 2L && nzchar(args[[2]])) normalizePath(args[[2]], mustWork = TRUE) else NA_character_

processed_root <- file.path(project_root, "data", "processed", "jgw03_three_cell_lines")
results_root <- file.path(project_root, "results", "jgw03_three_cell_lines")
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "cross_lineage"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "provenance"), recursive = TRUE, showWarnings = FALSE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
assert_true <- function(x, fmt, ...) if (!isTRUE(x)) stopf(fmt, ...)

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, file = con, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

config <- yaml::read_yaml(file.path(project_root, "config", "contrasts.yaml"))
seed <- as.integer(config$global$random_seed)
set.seed(seed)

ledger <- utils::read.delim(
  file.path(project_root, "metadata", "samples.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
geo_long <- utils::read.delim(
  file.path(project_root, "metadata", "geo_series_metadata_long.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
datasets_ledger <- utils::read.delim(
  file.path(project_root, "metadata", "datasets.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
orthology <- utils::read.delim(
  file.path(project_root, "data", "processed", "gene_sets", "ncbi_human_mouse_rat_one_to_one.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)
reactome_long <- utils::read.delim(
  file.path(project_root, "data", "processed", "gene_sets", "reactome97_human_long.tsv"),
  check.names = FALSE, na.strings = c("", "NA")
)

dataset_specs <- list(
  GSE289929 = list(
    cell_line = "PC12", organism = "Rattus norvegicus", tax_id = 10116L,
    raw_type = "rsem_tar", max_dose = 100, reference = "RN6",
    tar = file.path(project_root, "data", "raw", "GSE289929", "GSE289929_RAW.tar"),
    config_path = c("canonical_subset", "gsm_sets"),
    config_vehicle = "vehicle", config_treatment = c("dose_100_ug_ml", "dose_20_ug_ml", "dose_4_ug_ml"),
    interpretation = "primary_neural_like_whole_formula_perturbation"
  ),
  GSE227494 = list(
    cell_line = "C2C12", organism = "Mus musculus", tax_id = 10090L,
    raw_type = "matrix", max_dose = 500, reference = "GRCm38.86",
    matrix = file.path(project_root, "data", "raw", "GSE227494", "GSE227494_JGT_RCM_C2C12.txt.gz"),
    config_path = c("canonical_subset", "gsm_sets"),
    config_vehicle = "vehicle", config_treatment = "jgt_all_doses",
    interpretation = "muscle_context_calibration_not_external_replication"
  ),
  GSE299063 = list(
    cell_line = "HepG2", organism = "Homo sapiens", tax_id = 9606L,
    raw_type = "rsem_tar", max_dose = 500, reference = "hg38",
    tar = file.path(project_root, "data", "raw", "GSE299063", "GSE299063_RAW.tar"),
    config_path = c("canonical_batch_A_subset", "gsm_sets"),
    config_vehicle = "matched_vehicle", config_treatment = "jgt_all_doses",
    interpretation = "off_tissue_context_calibration_not_neural_efficacy"
  )
)

get_nested <- function(x, keys) {
  for (key in keys) x <- x[[key]]
  x
}

normalize_sample_table <- function(gse, spec) {
  x <- ledger[ledger$dataset_id == gse & ledger$include_primary == "TRUE", , drop = FALSE]
  x <- x[order(as.integer(x$sample_order)), , drop = FALSE]
  assert_true(nrow(x) == 12L, "%s must contain exactly 12 primary libraries; observed %d", gse, nrow(x))
  assert_true(!anyDuplicated(x$sample_id), "%s contains duplicate GSM identifiers", gse)
  assert_true(all(x$tissue_or_cell_line == spec$cell_line), "%s cell-line mismatch", gse)
  assert_true(length(unique(x$comparator_set)) == 1L, "%s has more than one comparator set", gse)
  assert_true(length(unique(x$batch)) == 1L, "%s canonical subset spans unresolved batches", gse)
  assert_true(all(as.numeric(x$exposure_duration_hour) == 24), "%s exposure-duration mismatch", gse)
  assert_true(sum(x$condition == "vehicle") == 3L, "%s must have three vehicle libraries", gse)
  assert_true(sum(x$condition == "treatment") == 9L, "%s must have nine JGW03 libraries", gse)
  assert_true(all(x$intervention_code[x$condition == "treatment"] == "JGW03"), "%s treatment code mismatch", gse)
  assert_true(all(x$formula_ratio[x$condition == "treatment"] == "1:1"), "%s ratio mismatch", gse)
  assert_true(all(x$solvent[x$condition == "treatment"] == "water"), "%s solvent mismatch", gse)
  assert_true(all(x$extraction_method[x$condition == "treatment"] == "CEM"), "%s extraction-mode mismatch", gse)

  ds_cfg <- config$datasets[[gse]]
  gsm_cfg <- get_nested(ds_cfg, spec$config_path)
  cfg_vehicle <- unlist(gsm_cfg[[spec$config_vehicle]], use.names = FALSE)
  cfg_treatment <- unlist(gsm_cfg[spec$config_treatment], use.names = FALSE)
  assert_true(setequal(cfg_vehicle, x$sample_id[x$condition == "vehicle"]),
              "%s vehicle GSMs disagree between samples.tsv and contrasts.yaml", gse)
  assert_true(setequal(cfg_treatment, x$sample_id[x$condition == "treatment"]),
              "%s treatment GSMs disagree between samples.tsv and contrasts.yaml", gse)

  x$dose_numeric <- ifelse(x$condition == "vehicle", 0, safe_numeric(x$dose_value))
  assert_true(!anyNA(x$dose_numeric), "%s contains unresolved numeric doses", gse)
  expected_doses <- if (gse == "GSE289929") c(0, 4, 20, 100) else c(0, 20, 100, 500)
  assert_true(identical(sort(unique(x$dose_numeric)), expected_doses), "%s dose sequence mismatch", gse)
  assert_true(all(table(x$dose_numeric) == 3L), "%s does not have n=3 at every dose", gse)

  x$dose_group <- factor(
    ifelse(x$condition == "vehicle", "vehicle",
           ifelse(x$dose_level == "low", "low",
                  ifelse(x$dose_level == "medium", "medium", "high"))),
    levels = c("vehicle", "low", "medium", "high")
  )
  x$ordered_dose <- as.numeric(x$dose_group) - 1
  x$log10_dose_plus_1 <- log10(x$dose_numeric + 1)
  assert_true(!anyNA(x$dose_group), "%s dose-level mapping is incomplete", gse)
  assert_true(max(x$dose_numeric) == spec$max_dose, "%s high dose differs from specification", gse)
  x
}

validate_geo_metadata <- function(gse, samples, spec) {
  g <- geo_long[geo_long$dataset_id == gse & geo_long$sample_id %in% samples$sample_id, , drop = FALSE]
  assert_true(length(unique(g$sample_id)) == nrow(samples), "%s official series metadata misses locked GSMs", gse)
  rows <- list()
  add <- function(check_id, status, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(dataset_id = gse, check_id = check_id,
                                             status = status, detail = detail)
  }
  add("locked_gsm_count", if (length(unique(g$sample_id)) == 12L) "PASS" else "FAIL",
      sprintf("12 expected; %d found", length(unique(g$sample_id))))

  source_rows <- g[g$field == "Sample_source_name_ch1", , drop = FALSE]
  source_ok <- all(source_rows$value == spec$cell_line) && length(unique(source_rows$sample_id)) == 12L
  add("cell_line", if (source_ok) "PASS" else "FAIL",
      sprintf("Official GEO source_name expected %s", spec$cell_line))

  org_rows <- g[g$field == "Sample_organism_ch1", , drop = FALSE]
  org_ok <- all(org_rows$value == spec$organism) && length(unique(org_rows$sample_id)) == 12L
  add("organism", if (org_ok) "PASS" else "FAIL",
      sprintf("Official GEO organism expected %s", spec$organism))

  for (i in seq_len(nrow(samples))) {
    gsm <- samples$sample_id[[i]]
    vals <- paste(g$value[g$sample_id == gsm & g$field %in% c("Sample_title", "Sample_characteristics_ch1", "Sample_description")],
                  collapse = " | ")
    if (samples$condition[[i]] == "vehicle") {
      ok <- grepl("vehicle", vals, ignore.case = TRUE)
      detail <- sprintf("%s must be an accession-matched vehicle", gsm)
    } else {
      dose <- format(samples$dose_numeric[[i]], scientific = FALSE, trim = TRUE)
      ok <- grepl("JGT|JGW03", vals, ignore.case = TRUE) && grepl(dose, vals, fixed = TRUE)
      detail <- sprintf("%s must identify JGT/JGW03 at %s ug/mL", gsm, dose)
    }
    add(paste0("sample_identity_", gsm), if (ok) "PASS" else "FAIL", detail)
  }

  geo_time <- g[g$field == "Sample_characteristics_ch1" & grepl("time:", g$value, ignore.case = TRUE), , drop = FALSE]
  if (nrow(geo_time) && all(grepl("24", geo_time$value))) {
    add("exposure_duration", "PASS", "24 h recorded in official sample characteristics")
  } else {
    ds <- datasets_ledger[datasets_ledger$dataset_id == gse, , drop = FALSE]
    ds_ok <- nrow(ds) == 1L && grepl("24 h", ds$timepoints, fixed = TRUE)
    add("exposure_duration", if (ds_ok) "PASS_WITH_LIMITATION" else "FAIL",
        "24 h is locked from the primary study/data ledger but is not present as a GEO sample characteristic")
  }

  processing <- unique(g$value[g$field == "Sample_data_processing"])
  if (gse == "GSE299063" && any(grepl("Rn6", processing, ignore.case = TRUE))) {
    add("reference_annotation", "PASS_WITH_LIMITATION",
        "GEO processing text contains an Rn6 copy error, whereas organism=Homo sapiens, Assembly=hg38 and the deposited files contain ENSG identifiers")
  } else {
    add("reference_annotation", "PASS", sprintf("Official processing/assembly metadata support %s", spec$reference))
  }
  out <- do.call(rbind, rows)
  assert_true(!any(out$status == "FAIL"), "%s failed official metadata validation", gse)
  out
}

extract_rsem_counts <- function(gse, samples, spec, input_dir) {
  assert_true(file.exists(spec$tar), "%s raw archive is missing", gse)
  members <- utils::untar(spec$tar, list = TRUE)
  selected <- character(nrow(samples))
  for (i in seq_len(nrow(samples))) {
    gsm <- samples$sample_id[[i]]
    hit <- grep(paste0("^", gsm, "_.*\\.STAR\\.RSEM\\.genes\\.results\\.gz$"), members, value = TRUE)
    assert_true(length(hit) == 1L, "%s must have one RSEM member for %s; found %d", gse, gsm, length(hit))
    assert_true(basename(hit) == samples$source_column_or_member[[i]],
                "%s archive member for %s disagrees with samples.tsv", gse, gsm)
    selected[[i]] <- hit
  }
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(spec$tar, files = selected, exdir = input_dir)

  tables <- lapply(selected, function(member) {
    p <- file.path(input_dir, member)
    z <- utils::read.delim(gzfile(p), check.names = FALSE)
    required <- c("gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM")
    assert_true(all(required %in% names(z)), "%s has a malformed RSEM member: %s", gse, member)
    z
  })
  ids <- lapply(tables, `[[`, "gene_id")
  assert_true(all(vapply(ids[-1L], identical, logical(1), ids[[1L]])), "%s RSEM gene rows/order differ across samples", gse)
  counts <- do.call(cbind, lapply(tables, `[[`, "expected_count"))
  rownames(counts) <- ids[[1L]]
  colnames(counts) <- samples$sample_id
  feature <- tables[[1L]][, c("gene_id", "length", "effective_length"), drop = FALSE]
  list(counts = counts, feature = feature, members = selected)
}

read_c2c12_counts <- function(gse, samples, spec) {
  assert_true(file.exists(spec$matrix), "%s processed expected-count matrix is missing", gse)
  # GEO deposited this combined matrix as whitespace-delimited text despite
  # the .txt.gz suffix; an explicit tab separator yields zero data columns.
  x <- utils::read.table(gzfile(spec$matrix), header = TRUE, sep = "", row.names = 1L,
                         check.names = FALSE, quote = "", comment.char = "")
  cols <- samples$source_column_or_member
  assert_true(all(cols %in% colnames(x)), "%s matrix does not contain every locked source column", gse)
  counts <- as.matrix(x[, cols, drop = FALSE])
  storage.mode(counts) <- "numeric"
  colnames(counts) <- samples$sample_id
  list(counts = counts,
       feature = data.frame(gene_id = rownames(counts), length = NA_real_, effective_length = NA_real_),
       members = cols)
}

annotate_features <- function(raw_ids, tax_id) {
  base <- sub("\\.[0-9]+$", "", raw_ids)
  if (tax_id == 10116L) {
    symbol <- AnnotationDbi::mapIds(org.Rn.eg.db, keys = unique(base), column = "SYMBOL",
                                    keytype = "ENSEMBL", multiVals = "first")
    symbol <- unname(symbol[base])
  } else if (tax_id == 9606L) {
    symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(base), column = "SYMBOL",
                                    keytype = "ENSEMBL", multiVals = "first")
    symbol <- unname(symbol[base])
  } else if (tax_id == 10090L) {
    # The GSE227494 combined matrix is deposited with mouse gene symbols.
    valid <- base %in% AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
    symbol <- ifelse(valid, base, NA_character_)
  } else stopf("Unsupported taxon: %s", tax_id)
  data.frame(raw_gene_id = raw_ids, base_gene_id = base, gene_symbol = symbol,
             annotation_status = ifelse(is.na(symbol), "unmapped", "mapped"))
}

map_to_human <- function(annotation, tax_id) {
  x <- annotation
  x$human_symbol <- NA_character_
  if (tax_id == 9606L) {
    x$human_symbol <- toupper(x$gene_symbol)
  } else {
    m <- orthology[orthology$other_tax_id == tax_id, c("human_symbol", "other_symbol"), drop = FALSE]
    m <- m[!is.na(m$human_symbol) & !is.na(m$other_symbol), , drop = FALSE]
    m$other_key <- toupper(m$other_symbol)
    m <- m[!duplicated(m$other_key), , drop = FALSE]
    lookup <- setNames(toupper(m$human_symbol), m$other_key)
    x$human_symbol <- unname(lookup[toupper(x$gene_symbol)])
  }
  x
}

robust_z <- function(x) {
  s <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - stats::median(x, na.rm = TRUE)) / s
}

adjudicate_qc_metrics <- function(qc) {
  # Correlation and PCA are treatment-sensitive expression-pattern screens,
  # not stand-alone technical failures. A removal-eligible screen therefore
  # requires an extreme library-size diagnostic plus an independent pattern
  # diagnostic. Even then, exclusion requires source-level evidence.
  qc$library_size_screen_flag <- abs(qc$library_robust_z) > 3.5
  qc$expression_pattern_screen_flag <- qc$correlation_robust_z < -3.5 |
    qc$pca_distance_robust_z > 3.5
  qc$corroborated_technical_failure <- qc$library_size_screen_flag &
    qc$expression_pattern_screen_flag
  if ("objective_qc_flag" %in% names(qc)) qc$objective_qc_flag <- NULL
  qc$excluded_from_primary <- FALSE
  qc$adjudication <- ifelse(
    qc$corroborated_technical_failure,
    "requires_source_level_review_but_retained_in_primary",
    ifelse(qc$expression_pattern_screen_flag,
           "expression_pattern_screen_only_treatment_sensitive_retained",
           "no_corroborated_technical_failure_retained")
  )
  qc
}

plot_qc <- function(gse, y, logcpm, samples, fit_high, result_dir) {
  dir.create(file.path(result_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  group_cols <- c(vehicle = "#4D4D4D", low = "#56B4E9", medium = "#E69F00", high = "#D55E00")
  sample_df <- data.frame(sample_id = samples$sample_id, dose_group = samples$dose_group,
                          dose_numeric = samples$dose_numeric)

  p_lib <- ggplot(sample_df, aes(x = sample_id, y = colSums(y$counts), fill = dose_group)) +
    geom_col() + scale_fill_manual(values = group_cols) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",", accuracy = 1)) +
    labs(x = NULL, y = "Deposited expected-count library size", title = paste(gse, "library sizes")) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(result_dir, "figures", paste0(gse, "_library_sizes.png")), p_lib, width = 8, height = 5, dpi = 180)

  vars <- apply(logcpm, 1L, stats::var)
  top <- order(vars, decreasing = TRUE)[seq_len(min(5000L, length(vars)))]
  pca <- stats::prcomp(t(logcpm[top, , drop = FALSE]), scale. = FALSE)
  pv <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  pca_df <- cbind(sample_df, PC1 = pca$x[, 1L], PC2 = pca$x[, 2L])
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = dose_group, label = sample_id)) +
    geom_point(size = 3) + geom_text(vjust = -0.8, size = 2.5, check_overlap = TRUE) +
    scale_color_manual(values = group_cols) + theme_bw(base_size = 10) +
    labs(title = paste(gse, "PCA (top variable filtered genes)"),
         x = sprintf("PC1 (%.1f%%)", pv[[1L]]), y = sprintf("PC2 (%.1f%%)", pv[[2L]]))
  ggsave(file.path(result_dir, "figures", paste0(gse, "_PCA.png")), p_pca, width = 7, height = 5.5, dpi = 180)

  cor_mat <- stats::cor(logcpm, method = "pearson")
  ann <- data.frame(dose_group = samples$dose_group, row.names = samples$sample_id)
  grDevices::png(file.path(result_dir, "figures", paste0(gse, "_sample_correlation.png")),
                 width = 1800, height = 1600, res = 180)
  pheatmap::pheatmap(cor_mat, annotation_col = ann, annotation_row = ann,
                     main = paste(gse, "sample correlation"), border_color = NA)
  grDevices::dev.off()

  grDevices::png(file.path(result_dir, "figures", paste0(gse, "_BCV.png")),
                 width = 1500, height = 1200, res = 180)
  edgeR::plotBCV(y, main = paste(gse, "biological coefficient of variation"))
  grDevices::dev.off()
  grDevices::png(file.path(result_dir, "figures", paste0(gse, "_QL_dispersion.png")),
                 width = 1500, height = 1200, res = 180)
  edgeR::plotQLDisp(fit_high, main = paste(gse, "quasi-likelihood dispersion"))
  grDevices::dev.off()

  med_cor <- apply(cor_mat, 2L, function(z) stats::median(z[names(z) != names(which.max(z))], na.rm = TRUE))
  # Recompute median off-diagonal explicitly because self-correlation names can be duplicated in generic code.
  med_cor <- vapply(seq_len(ncol(cor_mat)), function(j) stats::median(cor_mat[-j, j]), numeric(1))
  pca_dist <- sqrt(rowSums(scale(pca$x[, seq_len(min(5L, ncol(pca$x))), drop = FALSE])^2))
  qc <- data.frame(
    dataset_id = gse,
    sample_id = samples$sample_id,
    dose_group = as.character(samples$dose_group),
    dose_ug_ml = samples$dose_numeric,
    raw_library_size = colSums(y$counts),
    tmm_norm_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    detected_raw_ge_1 = colSums(y$counts >= 1),
    detected_filtered_ge_1 = colSums(y$counts >= 1),
    median_pairwise_correlation = med_cor,
    PC1 = pca$x[, 1L], PC2 = pca$x[, 2L], pca_robust_distance_5pc = pca_dist,
    library_robust_z = robust_z(log10(colSums(y$counts) + 1)),
    correlation_robust_z = robust_z(med_cor),
    pca_distance_robust_z = robust_z(pca_dist)
  )
  adjudicate_qc_metrics(qc)
}

ql_ci <- function(tab, test) {
  f <- pmax(tab$F, 0)
  se <- ifelse(f > 0, abs(tab$logFC) / sqrt(f), NA_real_)
  df <- test$df.total
  if (length(df) == 1L) df <- rep(df, nrow(tab))
  crit <- stats::qt(0.975, df = df)
  data.frame(ql_approx_se = se,
             ql_approx_CI_low = tab$logFC - crit * se,
             ql_approx_CI_high = tab$logFC + crit * se,
             ql_df_total = df)
}

edge_table <- function(test, annotation, contrast_name) {
  tab <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  assert_true(identical(rownames(tab), annotation$raw_gene_id), "Gene table order mismatch in %s", contrast_name)
  ci <- ql_ci(tab, test)
  out <- cbind(annotation, contrast = contrast_name, tab, ci)
  out$FDR <- stats::p.adjust(out$PValue, method = "BH")
  out
}

limma_voom_table <- function(y, design, contrast = NULL, coef = NULL, annotation, contrast_name) {
  v <- limma::voom(y, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  if (!is.null(contrast)) {
    fit <- limma::contrasts.fit(fit, contrast)
    coef <- 1L
  }
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  tab <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = 0.95)
  assert_true(identical(rownames(tab), annotation$raw_gene_id), "limma-voom row mismatch in %s", contrast_name)
  data.frame(annotation, contrast = contrast_name,
             logFC = tab$logFC, AveExpr = tab$AveExpr, t = tab$t,
             P.Value = tab$P.Value, adj.P.Val = tab$adj.P.Val,
             CI_low = tab$CI.L, CI_high = tab$CI.R, check.names = FALSE)
}

deseq_table <- function(counts, keep, samples, design_type, annotation, contrast_name) {
  rounded <- round(counts[keep, , drop = FALSE])
  coldata <- data.frame(row.names = samples$sample_id,
                        dose_group = samples$dose_group,
                        ordered_dose = samples$ordered_dose)
  if (design_type == "categorical") {
    dds <- DESeq2::DESeqDataSetFromMatrix(rounded, coldata, design = ~ dose_group)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
    res <- DESeq2::results(dds, contrast = c("dose_group", "high", "vehicle"), alpha = 0.05)
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(rounded, coldata, design = ~ ordered_dose)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
    res <- DESeq2::results(dds, name = "ordered_dose", alpha = 0.05)
  }
  res <- as.data.frame(res)
  ann <- annotation[keep, , drop = FALSE]
  assert_true(identical(rownames(res), ann$raw_gene_id), "DESeq2 row mismatch in %s", contrast_name)
  data.frame(ann, contrast = contrast_name, rounded_expected_counts_sensitivity = TRUE,
             baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
             lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
}

collapse_expression_by_symbol <- function(logcpm, annotation, symbol_col) {
  symbols <- annotation[[symbol_col]]
  valid <- !is.na(symbols) & nzchar(symbols)
  idx <- which(valid)
  means <- rowMeans(logcpm[idx, , drop = FALSE])
  ord <- order(symbols[idx], -means, annotation$raw_gene_id[idx])
  idx <- idx[ord]
  selected <- idx[!duplicated(symbols[idx])]
  out <- logcpm[selected, , drop = FALSE]
  rownames(out) <- symbols[selected]
  list(expression = out, selected_index = selected,
       mapping = data.frame(raw_gene_id = annotation$raw_gene_id[selected],
                            gene_symbol = annotation$gene_symbol[selected],
                            mapped_symbol = symbols[selected],
                            mean_logCPM = rowMeans(out)))
}

fit_activity_matrix <- function(activity, samples, family_name) {
  stopifnot(identical(colnames(activity), samples$sample_id))
  design_high <- model.matrix(~ 0 + dose_group, data = samples)
  colnames(design_high) <- levels(samples$dose_group)
  con_high <- limma::makeContrasts(high - vehicle, levels = design_high)
  design_trend <- model.matrix(~ ordered_dose, data = samples)
  design_log <- model.matrix(~ log10_dose_plus_1, data = samples)

  fit_one <- function(design, contrast = NULL, coef = NULL, label) {
    fit <- limma::lmFit(activity, design)
    if (!is.null(contrast)) {
      fit <- limma::contrasts.fit(fit, contrast)
      coef <- 1L
    }
    fit <- limma::eBayes(fit, robust = TRUE, trend = FALSE)
    tab <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = 0.95)
    data.frame(feature = rownames(tab), family = family_name, contrast = label,
               effect = tab$logFC, CI_low = tab$CI.L, CI_high = tab$CI.R,
               t = tab$t, PValue = tab$P.Value,
               FDR = stats::p.adjust(tab$P.Value, method = "BH"))
  }
  list(
    high = fit_one(design_high, contrast = con_high, label = "high_vs_vehicle"),
    trend = fit_one(design_trend, coef = 2L, label = "ordered_dose_per_step"),
    logtrend = fit_one(design_log, coef = 2L, label = "log10_dose_plus_1")
  )
}

run_decoupler_activity <- function(expr, network, source_col, target_col, mor_col, samples, family_name, minsize) {
  net <- data.frame(source = network[[source_col]], target = network[[target_col]], mor = network[[mor_col]])
  net <- net[complete.cases(net) & net$target %in% rownames(expr), , drop = FALSE]
  assert_true(nrow(net) > 0L, "No mapped edges for %s", family_name)
  act_long <- decoupleR::run_ulm(mat = expr, network = net, .source = source,
                                 .target = target, .mor = mor, minsize = minsize)
  nms <- names(act_long)
  source_name <- intersect(c("source", ".source"), nms)[1L]
  condition_name <- intersect(c("condition", ".condition"), nms)[1L]
  score_name <- intersect(c("score", "estimate", "statistic"), nms)[1L]
  assert_true(!anyNA(c(source_name, condition_name, score_name)), "Unexpected decoupleR output for %s", family_name)
  features <- sort(unique(act_long[[source_name]]))
  activity <- matrix(NA_real_, nrow = length(features), ncol = ncol(expr),
                     dimnames = list(features, colnames(expr)))
  for (i in seq_len(nrow(act_long))) {
    activity[as.character(act_long[[source_name]][[i]]), as.character(act_long[[condition_name]][[i]])] <-
      as.numeric(act_long[[score_name]][[i]])
  }
  activity <- activity[, samples$sample_id, drop = FALSE]
  assert_true(!anyNA(activity), "Incomplete activity matrix for %s", family_name)
  fits <- fit_activity_matrix(activity, samples, family_name)
  list(activity = activity, long = act_long, fits = fits,
       mapped_edges = nrow(net), mapped_sources = length(unique(net$source)))
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(z) unique(toupper(z[-c(1L, 2L)])))
  names(sets) <- vapply(parts, `[[`, character(1), 1L)
  sets
}

fit_program_scores <- function(scores, samples, collection) {
  res <- fit_activity_matrix(scores, samples, collection)
  res
}

make_program_scores <- function(expr_human, sets) {
  gene_sd <- apply(expr_human, 1L, stats::sd)
  expr_human <- expr_human[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
  z <- t(scale(t(expr_human)))
  scores <- t(vapply(sets, function(gs) {
    ix <- intersect(gs, rownames(z))
    colMeans(z[ix, , drop = FALSE])
  }, numeric(ncol(z))))
  rownames(scores) <- names(sets)
  colnames(scores) <- colnames(z)
  scores <- t(scale(t(scores)))
  scores
}

run_camera_collection <- function(gse, collection, y_human, expr_human, sets, term_meta,
                                  samples, designs, tests, result_dir, run_roast = FALSE) {
  mapped_sets <- lapply(sets, intersect, y = rownames(y_human))
  sizes <- lengths(mapped_sets)
  eligible <- sizes >= 10L & sizes <= 500L
  mapped_sets <- mapped_sets[eligible]
  assert_true(length(mapped_sets) > 0L, "%s has no eligible %s programs", gse, collection)
  index <- limma::ids2indices(mapped_sets, identifiers = rownames(y_human), remove.empty = TRUE)

  cam_high <- limma::camera(y_human, index = index, design = designs$high,
                            contrast = designs$contrast_high, inter.gene.cor = 0.01,
                            use.ranks = FALSE, sort = FALSE)
  cam_trend <- limma::camera(y_human, index = index, design = designs$trend,
                             contrast = 2L, inter.gene.cor = 0.01,
                             use.ranks = FALSE, sort = FALSE)
  cam_log <- limma::camera(y_human, index = index, design = designs$logtrend,
                           contrast = 2L, inter.gene.cor = 0.01,
                           use.ranks = FALSE, sort = FALSE)
  standardize_camera <- function(x, label) {
    x <- data.frame(term_id = rownames(x), x, row.names = NULL, check.names = FALSE)
    x$collection <- collection
    x$contrast <- label
    x$FDR <- stats::p.adjust(x$PValue, method = "BH")
    merge(x, term_meta, by = "term_id", all.x = TRUE, sort = FALSE)
  }
  cam_high <- standardize_camera(cam_high, "high_vs_vehicle")
  cam_trend <- standardize_camera(cam_trend, "ordered_dose_per_step")
  cam_log <- standardize_camera(cam_log, "log10_dose_plus_1")
  write_tsv(cam_high, file.path(result_dir, "programs", paste0(collection, "_camera_high_vs_vehicle.tsv")))
  write_tsv(cam_trend, file.path(result_dir, "programs", paste0(collection, "_camera_ordered_trend.tsv")))
  write_tsv(cam_log, file.path(result_dir, "programs", paste0(collection, "_camera_logdose_trend.tsv")))

  scores <- make_program_scores(expr_human, mapped_sets)
  score_fit <- fit_program_scores(scores, samples, paste0(collection, "_mean_z_score"))
  write_tsv(data.frame(term_id = rownames(scores), scores, check.names = FALSE),
            file.path(result_dir, "programs", paste0(collection, "_sample_scores.tsv")))
  for (nm in names(score_fit)) {
    write_tsv(score_fit[[nm]], file.path(result_dir, "programs", paste0(collection, "_score_", nm, ".tsv")))
  }

  score_high <- score_fit$high
  score_trend <- score_fit$trend
  names(score_high)[names(score_high) != "feature"] <- paste0("high_score_", names(score_high)[names(score_high) != "feature"])
  names(score_trend)[names(score_trend) != "feature"] <- paste0("trend_score_", names(score_trend)[names(score_trend) != "feature"])
  dual <- merge(cam_high, cam_trend, by = c("term_id", "collection", "term_name"),
                suffixes = c("_high", "_trend"), all = TRUE)
  dual <- merge(dual, score_high, by.x = "term_id", by.y = "feature", all.x = TRUE)
  dual <- merge(dual, score_trend, by.x = "term_id", by.y = "feature", all.x = TRUE)
  dual$high_CI_excludes_zero <- with(dual, high_score_CI_low * high_score_CI_high > 0)
  dual$trend_CI_excludes_zero <- with(dual, trend_score_CI_low * trend_score_CI_high > 0)
  dual$high_trend_direction_coherent <- with(dual, sign(high_score_effect) == sign(trend_score_effect))
  dual$internal_dual_response_rule <- with(
    dual,
    FDR_high <= 0.05 & FDR_trend <= 0.05 & high_CI_excludes_zero &
      trend_CI_excludes_zero & high_trend_direction_coherent
  )
  dual$gate_note <- "Internal dual-response rule only; J1/J2 additionally require an upstream S2/S3 candidate and prespecified polarity opposite the locked disease direction"
  write_tsv(dual, file.path(result_dir, "programs", paste0(collection, "_high_trend_dual.tsv")))

  # Rank-based enrichment is a sensitivity analysis and does not replace camera.
  fgsea_one <- function(test, label) {
    stat <- sign(test$table$logFC) * sqrt(pmax(test$table$F, 0))
    names(stat) <- rownames(test$table)
    stat <- sort(stat[is.finite(stat)], decreasing = TRUE)
    set.seed(seed)
    x <- suppressWarnings(fgsea::fgseaMultilevel(pathways = mapped_sets, stats = stat,
                                                  minSize = 10, maxSize = 500, eps = 0))
    x <- as.data.frame(x)
    if ("leadingEdge" %in% names(x)) x$leadingEdge <- vapply(x$leadingEdge, paste, collapse = ";", character(1))
    names(x)[names(x) == "pathway"] <- "term_id"
    x$collection <- collection
    x$contrast <- label
    merge(x, term_meta, by = "term_id", all.x = TRUE)
  }
  fg_high <- fgsea_one(tests$high, "high_vs_vehicle")
  fg_trend <- fgsea_one(tests$trend, "ordered_dose_per_step")
  write_tsv(fg_high, file.path(result_dir, "programs", paste0(collection, "_fgsea_high_vs_vehicle.tsv")))
  write_tsv(fg_trend, file.path(result_dir, "programs", paste0(collection, "_fgsea_ordered_trend.tsv")))

  roast_out <- NULL
  if (run_roast) {
    # 1,999 rotations provide a prespecified sensitivity with adequate tail
    # resolution while keeping the all-program calculation reproducible.
    set.seed(seed)
    ro_high <- limma::mroast(y_human, index = index, design = designs$high,
                             contrast = designs$contrast_high, nrot = 1999,
                             set.statistic = "mean", gene.weights = NULL, sort = "none")
    set.seed(seed)
    ro_trend <- limma::mroast(y_human, index = index, design = designs$trend,
                              contrast = 2L, nrot = 1999,
                              set.statistic = "mean", gene.weights = NULL, sort = "none")
    ro_high <- data.frame(term_id = rownames(ro_high), ro_high, row.names = NULL, check.names = FALSE)
    ro_trend <- data.frame(term_id = rownames(ro_trend), ro_trend, row.names = NULL, check.names = FALSE)
    write_tsv(ro_high, file.path(result_dir, "programs", paste0(collection, "_mroast1999_high_vs_vehicle.tsv")))
    write_tsv(ro_trend, file.path(result_dir, "programs", paste0(collection, "_mroast1999_ordered_trend.tsv")))
    roast_out <- list(high = ro_high, trend = ro_trend)
  }

  list(high = cam_high, trend = cam_trend, logtrend = cam_log, dual = dual,
       scores = scores, score_fit = score_fit, fgsea_high = fg_high,
       fgsea_trend = fg_trend, roast = roast_out,
       mapping = data.frame(term_id = names(mapped_sets), mapped_genes = lengths(mapped_sets)))
}

analysis_summary_columns <- c(
  "dataset_id", "cell_line", "tax_id", "n_vehicle", "n_low", "n_medium", "n_high",
  "doses_ug_ml", "genes_input", "genes_filterByExpr", "native_symbol_mapping_n",
  "strict_human_mapping_n", "metadata_limitations_n", "edgeR_high_FDR05_n",
  "edgeR_trend_FDR05_n", "edgeR_high_and_trend_FDR05_coherent_n",
  "reactome_eligible_n", "reactome_internal_dual_n", "syngo_sensitivity_eligible_n",
  "syngo_internal_dual_n", "progeny_high_FDR05_n", "progeny_trend_FDR05_n",
  "dorothea_high_FDR05_n", "dorothea_trend_FDR05_n",
  "expression_pattern_screen_flags_n", "corroborated_technical_failure_n",
  "samples_excluded_n"
)

order_analysis_summary <- function(frame) {
  missing_columns <- setdiff(analysis_summary_columns, names(frame))
  extra_columns <- setdiff(names(frame), analysis_summary_columns)
  if (length(missing_columns) || length(extra_columns)) {
    stop(sprintf(
      "JGW03 analysis-summary schema drift: missing=%s; extra=%s",
      paste(missing_columns, collapse = ","), paste(extra_columns, collapse = ",")
    ))
  }
  frame[, analysis_summary_columns, drop = FALSE]
}

analyze_one <- function(gse, spec) {
  message("Analyzing ", gse, " (", spec$cell_line, ")")
  result_dir <- file.path(results_root, gse)
  input_dir <- file.path(processed_root, gse, "selected_rsem")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  for (d in c("genes", "programs", "activities", "figures", "qc", "provenance"))
    dir.create(file.path(result_dir, d), recursive = TRUE, showWarnings = FALSE)

  samples <- normalize_sample_table(gse, spec)
  metadata_validation <- validate_geo_metadata(gse, samples, spec)
  write_tsv(metadata_validation, file.path(result_dir, "qc", "official_metadata_validation.tsv"))

  inp <- if (spec$raw_type == "rsem_tar")
    extract_rsem_counts(gse, samples, spec, input_dir) else
      read_c2c12_counts(gse, samples, spec)
  counts <- inp$counts
  assert_true(identical(colnames(counts), samples$sample_id), "%s count/sample order mismatch", gse)
  assert_true(all(is.finite(counts)) && all(counts >= 0), "%s has invalid expected counts", gse)
  assert_true(!anyDuplicated(rownames(counts)), "%s has duplicate raw gene identifiers", gse)

  annotation <- annotate_features(rownames(counts), spec$tax_id)
  annotation <- map_to_human(annotation, spec$tax_id)
  rownames(annotation) <- annotation$raw_gene_id
  write_tsv(annotation, file.path(processed_root, gse, "gene_annotation.tsv"))
  write_tsv_gz(data.frame(raw_gene_id = rownames(counts), counts, check.names = FALSE),
               file.path(processed_root, gse, "locked_expected_counts.tsv.gz"))
  write_tsv(samples, file.path(processed_root, gse, "locked_samples.tsv"))

  design_high <- model.matrix(~ 0 + dose_group, data = samples)
  colnames(design_high) <- levels(samples$dose_group)
  design_trend <- model.matrix(~ ordered_dose, data = samples)
  design_logtrend <- model.matrix(~ log10_dose_plus_1, data = samples)
  contrast_high <- limma::makeContrasts(high - vehicle, levels = design_high)
  assert_true(qr(design_high)$rank == ncol(design_high), "%s categorical design is rank deficient", gse)
  assert_true(qr(design_trend)$rank == ncol(design_trend), "%s trend design is rank deficient", gse)

  y0 <- edgeR::DGEList(counts = counts, genes = annotation)
  keep <- edgeR::filterByExpr(y0, design = design_high)
  assert_true(sum(keep) >= 1000L, "%s retained fewer than 1,000 expressed genes", gse)
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::normLibSizes(y, method = "TMM")

  y_high <- edgeR::estimateDisp(y, design_high, robust = TRUE)
  fit_high <- edgeR::glmQLFit(y_high, design_high, robust = TRUE)
  test_high <- edgeR::glmQLFTest(fit_high, contrast = contrast_high)
  y_trend <- edgeR::estimateDisp(y, design_trend, robust = TRUE)
  fit_trend <- edgeR::glmQLFit(y_trend, design_trend, robust = TRUE)
  test_trend <- edgeR::glmQLFTest(fit_trend, coef = 2L)
  y_log <- edgeR::estimateDisp(y, design_logtrend, robust = TRUE)
  fit_log <- edgeR::glmQLFit(y_log, design_logtrend, robust = TRUE)
  test_log <- edgeR::glmQLFTest(fit_log, coef = 2L)

  ann_keep <- annotation[keep, , drop = FALSE]
  gene_high <- edge_table(test_high, ann_keep, "high_vs_vehicle")
  gene_trend <- edge_table(test_trend, ann_keep, "ordered_dose_per_step")
  gene_log <- edge_table(test_log, ann_keep, "log10_dose_plus_1")
  write_tsv_gz(gene_high, file.path(result_dir, "genes", "edgeR_QL_high_vs_vehicle_all_genes.tsv.gz"))
  write_tsv_gz(gene_trend, file.path(result_dir, "genes", "edgeR_QL_ordered_trend_all_genes.tsv.gz"))
  write_tsv_gz(gene_log, file.path(result_dir, "genes", "edgeR_QL_logdose_trend_all_genes.tsv.gz"))

  voom_high <- limma_voom_table(y, design_high, contrast = contrast_high,
                                 annotation = ann_keep, contrast_name = "high_vs_vehicle")
  voom_trend <- limma_voom_table(y, design_trend, coef = 2L,
                                  annotation = ann_keep, contrast_name = "ordered_dose_per_step")
  write_tsv_gz(voom_high, file.path(result_dir, "genes", "limma_voom_high_vs_vehicle_all_genes.tsv.gz"))
  write_tsv_gz(voom_trend, file.path(result_dir, "genes", "limma_voom_ordered_trend_all_genes.tsv.gz"))

  deseq_high <- deseq_table(counts, keep, samples, "categorical", annotation, "high_vs_vehicle")
  deseq_trend <- deseq_table(counts, keep, samples, "trend", annotation, "ordered_dose_per_step")
  write_tsv_gz(deseq_high, file.path(result_dir, "genes", "DESeq2_rounded_high_vs_vehicle_all_genes.tsv.gz"))
  write_tsv_gz(deseq_trend, file.path(result_dir, "genes", "DESeq2_rounded_ordered_trend_all_genes.tsv.gz"))

  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  write_tsv_gz(data.frame(raw_gene_id = rownames(logcpm), logcpm, check.names = FALSE),
               file.path(processed_root, gse, "TMM_logCPM_filtered.tsv.gz"))
  qc <- plot_qc(gse, y_high, logcpm, samples, fit_high, result_dir)
  qc$genes_input <- nrow(counts)
  qc$genes_retained_filterByExpr <- sum(keep)
  qc$mapped_native_symbol <- sum(!is.na(annotation$gene_symbol))
  qc$mapped_strict_human_symbol <- sum(!is.na(annotation$human_symbol))
  write_tsv(qc, file.path(result_dir, "qc", "sample_qc_metrics.tsv"))

  method_summary <- data.frame(
    dataset_id = gse,
    cell_line = spec$cell_line,
    contrast = c("high_vs_vehicle", "ordered_dose_per_step"),
    edgeR_FDR05_up = c(sum(gene_high$FDR <= 0.05 & gene_high$logFC > 0),
                       sum(gene_trend$FDR <= 0.05 & gene_trend$logFC > 0)),
    edgeR_FDR05_down = c(sum(gene_high$FDR <= 0.05 & gene_high$logFC < 0),
                         sum(gene_trend$FDR <= 0.05 & gene_trend$logFC < 0)),
    edgeR_limma_logFC_spearman = c(
      stats::cor(gene_high$logFC, voom_high$logFC, method = "spearman", use = "complete.obs"),
      stats::cor(gene_trend$logFC, voom_trend$logFC, method = "spearman", use = "complete.obs")
    ),
    edgeR_DESeq2_logFC_spearman = c(
      stats::cor(gene_high$logFC, deseq_high$log2FoldChange, method = "spearman", use = "complete.obs"),
      stats::cor(gene_trend$logFC, deseq_trend$log2FoldChange, method = "spearman", use = "complete.obs")
    )
  )
  write_tsv(method_summary, file.path(result_dir, "qc", "method_concordance.tsv"))

  # Human-mapped feature space for cross-species Reactome/SynGO comparisons.
  human_collapsed <- collapse_expression_by_symbol(logcpm, ann_keep, "human_symbol")
  native_collapsed <- collapse_expression_by_symbol(logcpm, ann_keep, "gene_symbol")
  write_tsv(human_collapsed$mapping, file.path(result_dir, "provenance", "human_symbol_feature_selection.tsv"))
  write_tsv(native_collapsed$mapping, file.path(result_dir, "provenance", "native_symbol_feature_selection.tsv"))
  y_human <- y[human_collapsed$selected_index, , keep.lib.sizes = FALSE]
  rownames(y_human$counts) <- rownames(human_collapsed$expression)
  y_human$genes <- data.frame(human_symbol = rownames(y_human))
  y_human <- edgeR::normLibSizes(y_human, method = "TMM")
  # camera.DGEList and mroast.DGEList require an attached dispersion estimate.
  # The complete categorical design is used for this common program-testing
  # object; trend-specific edgeR tests below retain their own dispersion fit.
  y_human <- edgeR::estimateDisp(y_human, design_high, robust = TRUE)
  expr_human <- human_collapsed$expression
  expr_native <- native_collapsed$expression

  reactome_meta <- unique(reactome_long[, c("term_id", "term_name")])
  reactome_sets <- split(toupper(reactome_long$gene_symbol), reactome_long$term_id)
  reactome_sets <- lapply(reactome_sets, unique)
  designs <- list(high = design_high, trend = design_trend, logtrend = design_logtrend,
                  contrast_high = contrast_high)
  tests_human <- list(
    high = edgeR::glmQLFTest(edgeR::glmQLFit(edgeR::estimateDisp(y_human, design_high, robust = TRUE),
                                              design_high, robust = TRUE), contrast = contrast_high),
    trend = edgeR::glmQLFTest(edgeR::glmQLFit(edgeR::estimateDisp(y_human, design_trend, robust = TRUE),
                                               design_trend, robust = TRUE), coef = 2L)
  )
  reactome <- run_camera_collection(gse, "Reactome97", y_human, expr_human,
                                    reactome_sets, reactome_meta, samples, designs,
                                    tests_human, result_dir, run_roast = FALSE)

  syngo_path <- file.path(project_root, "data", "processed", "gene_sets",
                          "quickgo_syngo_mammalian_direct_human.gmt")
  syngo_sets <- read_gmt(syngo_path)
  syngo_meta <- data.frame(term_id = names(syngo_sets), term_name = names(syngo_sets))
  syngo <- run_camera_collection(gse, "QuickGO_SynGO_sensitivity", y_human, expr_human,
                                 syngo_sets, syngo_meta, samples, designs,
                                 tests_human, result_dir, run_roast = FALSE)

  if (spec$tax_id == 10090L) {
    progeny_path <- file.path(project_root, "data", "processed", "gene_sets", "progeny_mouse_top500.tsv")
    dorothea_path <- file.path(project_root, "data", "processed", "gene_sets", "dorothea_mouse_AB.tsv")
    progeny_expr <- expr_native
    dorothea_expr <- expr_native
  } else if (spec$tax_id == 10116L) {
    progeny_path <- file.path(project_root, "data", "processed", "gene_sets", "progeny_human_top500.tsv")
    dorothea_path <- file.path(project_root, "data", "processed", "gene_sets", "dorothea_rat_AB_one_to_one.tsv")
    progeny_expr <- expr_human
    dorothea_expr <- expr_native
  } else {
    progeny_path <- file.path(project_root, "data", "processed", "gene_sets", "progeny_human_top500.tsv")
    dorothea_path <- file.path(project_root, "data", "processed", "gene_sets", "dorothea_human_AB.tsv")
    progeny_expr <- expr_native
    dorothea_expr <- expr_native
  }
  progeny_net <- utils::read.delim(progeny_path, check.names = FALSE)
  dorothea_net <- utils::read.delim(dorothea_path, check.names = FALSE)
  progeny <- run_decoupler_activity(progeny_expr, progeny_net, "pathway", "gene_symbol", "weight",
                                    samples, "PROGENy_top500_ULM", minsize = 50L)
  dorothea <- run_decoupler_activity(dorothea_expr, dorothea_net, "tf", "target", "mor",
                                     samples, "DoRothEA_AB_ULM", minsize = 5L)
  write_tsv(data.frame(feature = rownames(progeny$activity), progeny$activity, check.names = FALSE),
            file.path(result_dir, "activities", "PROGENy_sample_activity.tsv"))
  write_tsv(data.frame(feature = rownames(dorothea$activity), dorothea$activity, check.names = FALSE),
            file.path(result_dir, "activities", "DoRothEA_sample_activity.tsv"))
  for (nm in names(progeny$fits))
    write_tsv(progeny$fits[[nm]], file.path(result_dir, "activities", paste0("PROGENy_", nm, ".tsv")))
  for (nm in names(dorothea$fits))
    write_tsv(dorothea$fits[[nm]], file.path(result_dir, "activities", paste0("DoRothEA_", nm, ".tsv")))

  dose_means <- do.call(rbind, lapply(levels(samples$dose_group), function(d) {
    ids <- samples$sample_id[samples$dose_group == d]
    data.frame(dataset_id = gse, dose_group = d,
               dose_ug_ml = unique(samples$dose_numeric[samples$dose_group == d]),
               mean_library_size = mean(colSums(counts[, ids, drop = FALSE])),
               median_library_size = stats::median(colSums(counts[, ids, drop = FALSE])))
  }))
  write_tsv(dose_means, file.path(result_dir, "qc", "dose_level_summary.tsv"))

  provenance <- data.frame(
    dataset_id = gse,
    source = if (spec$raw_type == "rsem_tar") spec$tar else spec$matrix,
    source_size_bytes = file.info(if (spec$raw_type == "rsem_tar") spec$tar else spec$matrix)$size,
    source_md5 = unname(tools::md5sum(if (spec$raw_type == "rsem_tar") spec$tar else spec$matrix)),
    selected_members_or_columns = paste(inp$members, collapse = ";"),
    reference = spec$reference,
    count_measure = "RSEM expected_count (fractional values retained for edgeR/limma; rounded only for labeled DESeq2 sensitivity)",
    interpretation = spec$interpretation
  )
  write_tsv(provenance, file.path(result_dir, "provenance", "input_provenance.tsv"))

  summary <- data.frame(
    dataset_id = gse, cell_line = spec$cell_line, tax_id = spec$tax_id,
    n_vehicle = sum(samples$dose_group == "vehicle"), n_low = sum(samples$dose_group == "low"),
    n_medium = sum(samples$dose_group == "medium"), n_high = sum(samples$dose_group == "high"),
    doses_ug_ml = paste(sort(unique(samples$dose_numeric)), collapse = ","),
    genes_input = nrow(counts), genes_filterByExpr = sum(keep),
    native_symbol_mapping_n = sum(!is.na(annotation$gene_symbol)),
    strict_human_mapping_n = sum(!is.na(annotation$human_symbol)),
    expression_pattern_screen_flags_n = sum(qc$expression_pattern_screen_flag),
    corroborated_technical_failure_n = sum(qc$corroborated_technical_failure),
    samples_excluded_n = sum(qc$excluded_from_primary),
    metadata_limitations_n = sum(metadata_validation$status == "PASS_WITH_LIMITATION"),
    edgeR_high_FDR05_n = sum(gene_high$FDR <= 0.05, na.rm = TRUE),
    edgeR_trend_FDR05_n = sum(gene_trend$FDR <= 0.05, na.rm = TRUE),
    edgeR_high_and_trend_FDR05_coherent_n = sum(
      gene_high$FDR <= 0.05 & gene_trend$FDR <= 0.05 &
        sign(gene_high$logFC) == sign(gene_trend$logFC), na.rm = TRUE),
    reactome_eligible_n = nrow(reactome$dual),
    reactome_internal_dual_n = sum(reactome$dual$internal_dual_response_rule, na.rm = TRUE),
    syngo_sensitivity_eligible_n = nrow(syngo$dual),
    syngo_internal_dual_n = sum(syngo$dual$internal_dual_response_rule, na.rm = TRUE),
    progeny_high_FDR05_n = sum(progeny$fits$high$FDR <= 0.05, na.rm = TRUE),
    progeny_trend_FDR05_n = sum(progeny$fits$trend$FDR <= 0.05, na.rm = TRUE),
    dorothea_high_FDR05_n = sum(dorothea$fits$high$FDR <= 0.05, na.rm = TRUE),
    dorothea_trend_FDR05_n = sum(dorothea$fits$trend$FDR <= 0.05, na.rm = TRUE)
  )
  summary <- order_analysis_summary(summary)
  write_tsv(summary, file.path(result_dir, "analysis_summary.tsv"))

  saveRDS(list(samples = samples, annotation = annotation, keep = keep, y = y,
               designs = designs, tests = list(high = test_high, trend = test_trend, logtrend = test_log),
               reactome = reactome, syngo = syngo, progeny = progeny,
               dorothea = dorothea, summary = summary),
          file.path(processed_root, gse, "analysis_objects.rds"), compress = "xz")
  list(gse = gse, spec = spec, samples = samples, summary = summary, qc = qc,
       metadata_validation = metadata_validation, gene_high = gene_high,
       gene_trend = gene_trend, reactome = reactome, syngo = syngo,
       progeny = progeny, dorothea = dorothea)
}

resume_enabled <- tolower(Sys.getenv("ROUTE_B2_JGW03_RESUME", "false")) %in% c("1", "true", "yes")
load_completed_one <- function(gse, spec) {
  object_path <- file.path(processed_root, gse, "analysis_objects.rds")
  if (!resume_enabled || !file.exists(object_path)) return(NULL)
  message("Reusing completed checkpoint for ", gse)
  obj <- readRDS(object_path)
  result_dir <- file.path(results_root, gse)
  qc <- utils::read.delim(file.path(result_dir, "qc", "sample_qc_metrics.tsv"), check.names = FALSE)
  qc <- adjudicate_qc_metrics(qc)
  write_tsv(qc, file.path(result_dir, "qc", "sample_qc_metrics.tsv"))
  summary <- obj$summary
  if ("objective_sample_flags_n" %in% names(summary)) summary$objective_sample_flags_n <- NULL
  summary$expression_pattern_screen_flags_n <- sum(qc$expression_pattern_screen_flag)
  summary$corroborated_technical_failure_n <- sum(qc$corroborated_technical_failure)
  summary$samples_excluded_n <- sum(qc$excluded_from_primary)
  summary <- order_analysis_summary(summary)
  write_tsv(summary, file.path(result_dir, "analysis_summary.tsv"))
  list(
    gse = gse, spec = spec, samples = obj$samples, summary = summary,
    qc = qc,
    metadata_validation = utils::read.delim(file.path(result_dir, "qc", "official_metadata_validation.tsv"), check.names = FALSE),
    gene_high = NULL, gene_trend = NULL, reactome = obj$reactome, syngo = obj$syngo,
    progeny = obj$progeny, dorothea = obj$dorothea
  )
}

analyses <- lapply(names(dataset_specs), function(gse) {
  cached <- load_completed_one(gse, dataset_specs[[gse]])
  if (is.null(cached)) analyze_one(gse, dataset_specs[[gse]]) else cached
})
names(analyses) <- names(dataset_specs)

merge_context_family <- function(extractor, family_label) {
  get_dual <- function(a) extractor(a)
  pc <- get_dual(analyses$GSE289929)
  c2 <- get_dual(analyses$GSE227494)
  hp <- get_dual(analyses$GSE299063)
  select_cols <- function(x, prefix) {
    keep <- intersect(c("term_id", "feature", "term_name", "high_score_effect", "high_score_CI_low",
                        "high_score_CI_high", "trend_score_effect", "trend_score_CI_low",
                        "trend_score_CI_high", "FDR_high", "FDR_trend",
                        "internal_dual_response_rule", "effect_high", "effect_trend",
                        "FDR_high_activity", "FDR_trend_activity"), names(x))
    z <- x[, keep, drop = FALSE]
    id <- intersect(c("term_id", "feature"), names(z))[1L]
    names(z)[names(z) != id] <- paste0(prefix, "_", names(z)[names(z) != id])
    names(z)[names(z) == id] <- "feature_id"
    z
  }
  out <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE),
                list(select_cols(pc, "PC12"), select_cols(c2, "C2C12"), select_cols(hp, "HepG2")))
  out$family <- family_label
  if (all(c("PC12_high_score_effect", "C2C12_high_score_effect", "HepG2_high_score_effect") %in% names(out))) {
    out$C2C12_same_direction_high <- sign(out$C2C12_high_score_effect) == sign(out$PC12_high_score_effect)
    out$HepG2_same_direction_high <- sign(out$HepG2_high_score_effect) == sign(out$PC12_high_score_effect)
    out$C2C12_same_direction_trend <- sign(out$C2C12_trend_score_effect) == sign(out$PC12_trend_score_effect)
    out$HepG2_same_direction_trend <- sign(out$HepG2_trend_score_effect) == sign(out$PC12_trend_score_effect)
    out$C2C12_to_PC12_abs_high_ratio <- abs(out$C2C12_high_score_effect) / abs(out$PC12_high_score_effect)
    out$HepG2_to_PC12_abs_high_ratio <- abs(out$HepG2_high_score_effect) / abs(out$PC12_high_score_effect)
    out$C2C12_to_PC12_abs_trend_ratio <- abs(out$C2C12_trend_score_effect) / abs(out$PC12_trend_score_effect)
    out$HepG2_to_PC12_abs_trend_ratio <- abs(out$HepG2_trend_score_effect) / abs(out$PC12_trend_score_effect)
    out$descriptive_hepg2_stronger_concordant_both <- with(
      out, HepG2_same_direction_high & HepG2_same_direction_trend &
        HepG2_to_PC12_abs_high_ratio >= 1 & HepG2_to_PC12_abs_trend_ratio >= 1
    )
    out$context_rule_status <- "DESCRIPTIVE_ONLY: dose maxima differ (PC12 100; C2C12/HepG2 500 ug/mL), so amplitude ratios do not by themselves pass/fail J3"
  }
  out
}

reactome_context <- merge_context_family(function(a) a$reactome$dual, "Reactome97")
syngo_context <- merge_context_family(function(a) a$syngo$dual, "QuickGO_SynGO_sensitivity")
write_tsv(reactome_context, file.path(results_root, "cross_lineage", "Reactome_context_comparison.tsv"))
write_tsv(syngo_context, file.path(results_root, "cross_lineage", "QuickGO_SynGO_context_comparison.tsv"))

activity_dual <- function(a, which) {
  z <- a[[which]]$fits
  h <- z$high[, c("feature", "effect", "CI_low", "CI_high", "FDR")]
  t <- z$trend[, c("feature", "effect", "CI_low", "CI_high", "FDR")]
  names(h)[-1L] <- paste0(names(h)[-1L], "_high")
  names(t)[-1L] <- paste0(names(t)[-1L], "_trend")
  merge(h, t, by = "feature", all = TRUE)
}
progeny_context <- merge_context_family(function(a) activity_dual(a, "progeny"), "PROGENy_top500_ULM")
dorothea_context <- merge_context_family(function(a) activity_dual(a, "dorothea"), "DoRothEA_AB_ULM")
write_tsv(progeny_context, file.path(results_root, "cross_lineage", "PROGENy_context_comparison.tsv"))
write_tsv(dorothea_context, file.path(results_root, "cross_lineage", "DoRothEA_context_comparison.tsv"))

cross_cor <- function(context, effect_stem, family) {
  cols <- paste0(c("PC12", "C2C12", "HepG2"), "_", effect_stem)
  if (!all(cols %in% names(context))) return(NULL)
  pairs <- list(c("PC12", "C2C12"), c("PC12", "HepG2"), c("C2C12", "HepG2"))
  do.call(rbind, lapply(pairs, function(p) {
    x <- context[[paste0(p[[1]], "_", effect_stem)]]
    y <- context[[paste0(p[[2]], "_", effect_stem)]]
    ok <- complete.cases(x, y)
    n_ok <- sum(ok)
    data.frame(family = family, effect = effect_stem, lineage_1 = p[[1]], lineage_2 = p[[2]],
               n_common = n_ok,
               spearman_rho = if (n_ok >= 2L) stats::cor(x[ok], y[ok], method = "spearman") else NA_real_,
               kendall_tau = if (n_ok >= 2L) stats::cor(x[ok], y[ok], method = "kendall") else NA_real_)
  }))
}
cross_correlations <- do.call(rbind, Filter(Negate(is.null), list(
  cross_cor(reactome_context, "high_score_effect", "Reactome97"),
  cross_cor(reactome_context, "trend_score_effect", "Reactome97"),
  cross_cor(syngo_context, "high_score_effect", "QuickGO_SynGO_sensitivity"),
  cross_cor(syngo_context, "trend_score_effect", "QuickGO_SynGO_sensitivity"),
  cross_cor(progeny_context, "effect_high", "PROGENy_top500_ULM"),
  cross_cor(progeny_context, "effect_trend", "PROGENy_top500_ULM"),
  cross_cor(dorothea_context, "effect_high", "DoRothEA_AB_ULM"),
  cross_cor(dorothea_context, "effect_trend", "DoRothEA_AB_ULM")
)))
write_tsv(cross_correlations, file.path(results_root, "cross_lineage", "effect_correlations.tsv"))

# J gates cannot be called from drug data alone.  If a compatible upstream
# candidate file is supplied, it is joined without changing its disease
# direction; otherwise the required status is explicitly pending.
j_gate <- data.frame(
  gate = c("J1", "J2", "J3", "J4"),
  status = c(
    "NOT_EVALUABLE_PENDING_S2_S3_AND_LOCKED_DISEASE_DIRECTION",
    "NOT_EVALUABLE_PENDING_S2_S3_AND_LOCKED_DISEASE_DIRECTION",
    "NOT_EVALUABLE_PENDING_J2_AND_A_PROSPECTIVELY_QUANTIFIED_CONTEXT_RULE",
    "NOT_ASSESSED_IN_THREE_CANONICAL_CELL_LINE_SUBTASK"
  ),
  reason = c(
    "High-dose PC12 response is available, but J1 also requires membership in S2/S3 and the expected polarity relative to the locked disease direction.",
    "High-dose and ordered-dose PC12 results are available, but J2 also requires membership in S2/S3 and the expected polarity relative to the locked disease direction.",
    "C2C12/HepG2 context estimates are available; high doses differ (100 versus 500 ug/mL), and contrasts.yaml does not define a numeric J3 amplitude threshold.",
    "J4 requires preparation/ratio sensitivity analyses outside this canonical three-cell-line task."
  )
)

s1_gate_path <- file.path(project_root, "results", "GSE16710", "program_level", "S1_GATE.tsv")
s1_zero <- FALSE
s1_provenance <- data.frame(
  source_file = s1_gate_path, source_rows = NA_integer_, S1_pass_n = NA_integer_,
  S1_missing_n = NA_integer_, source_md5 = NA_character_, stop_rule = NA_character_
)
if (file.exists(s1_gate_path)) {
  s1 <- utils::read.delim(s1_gate_path, check.names = FALSE)
  assert_true("S1_pass" %in% names(s1), "S1 gate file lacks S1_pass: %s", s1_gate_path)
  s1_bool <- if (is.logical(s1$S1_pass)) s1$S1_pass else
    ifelse(toupper(as.character(s1$S1_pass)) == "TRUE", TRUE,
           ifelse(toupper(as.character(s1$S1_pass)) == "FALSE", FALSE, NA))
  s1_zero <- !any(s1_bool, na.rm = TRUE) && !anyNA(s1_bool)
  s1_provenance <- data.frame(
    source_file = s1_gate_path, source_rows = nrow(s1),
    S1_pass_n = sum(s1_bool, na.rm = TRUE), S1_missing_n = sum(is.na(s1_bool)),
    source_md5 = unname(tools::md5sum(s1_gate_path)),
    stop_rule = if (s1_zero) "STOP_NO_S2" else "not_triggered"
  )
}
write_tsv(s1_provenance, file.path(results_root, "provenance", "upstream_S1_gate_provenance.tsv"))

if (s1_zero) {
  j_gate$status <- "NOT_IDENTIFIABLE_NO_S2_INPUT"
  j_gate$reason <- c(
    "STOP_NO_S2: zero S1 programs can enter S2; high-dose PC12 data are retained as exploratory background but J1 is not evaluated.",
    "STOP_NO_S2: zero S1 programs can enter S2; ordered-dose PC12 data are retained as exploratory background but J2 is neither positive nor negative.",
    "J3 is downstream of J2 and is not identifiable when the sequential S2 input is empty; C2C12/HepG2 remain descriptive context calibrators.",
    "J4 is downstream of J2/J3 and is not identifiable when the sequential S2 input is empty."
  )
} else if (!is.na(candidate_file)) {
  candidates <- utils::read.delim(candidate_file, check.names = FALSE)
  required <- c("term_id", "disease_direction", "S2_or_S3")
  assert_true(all(required %in% names(candidates)),
              "Candidate file must contain: %s", paste(required, collapse = ", "))
  pc <- analyses$GSE289929$reactome$dual
  joined <- merge(candidates, pc, by = "term_id", all.x = TRUE)
  dir_sign <- ifelse(tolower(joined$disease_direction) %in% c("up", "positive", "1"), 1,
                     ifelse(tolower(joined$disease_direction) %in% c("down", "negative", "-1"), -1, NA))
  joined$expected_drug_sign <- -dir_sign
  joined$J1_pass <- joined$S2_or_S3 & joined$FDR_high <= 0.05 &
    joined$high_CI_excludes_zero & sign(joined$high_score_effect) == joined$expected_drug_sign
  joined$J2_pass <- joined$J1_pass & joined$FDR_trend <= 0.05 &
    joined$trend_CI_excludes_zero & joined$high_trend_direction_coherent &
    sign(joined$trend_score_effect) == joined$expected_drug_sign
  write_tsv(joined, file.path(results_root, "cross_lineage", "candidate_conditioned_J1_J2.tsv"))
  j_gate$status[j_gate$gate == "J1"] <- if (any(joined$J1_pass, na.rm = TRUE)) "EVALUATED_WITH_PASSING_CANDIDATES" else "EVALUATED_ZERO_PASS"
  j_gate$status[j_gate$gate == "J2"] <- if (any(joined$J2_pass, na.rm = TRUE)) "EVALUATED_WITH_PASSING_CANDIDATES" else "EVALUATED_ZERO_PASS"
  j_gate$reason[j_gate$gate == "J1"] <- sprintf("%d of %d supplied candidates pass J1.", sum(joined$J1_pass, na.rm = TRUE), nrow(joined))
  j_gate$reason[j_gate$gate == "J2"] <- sprintf("%d of %d supplied candidates pass J2.", sum(joined$J2_pass, na.rm = TRUE), nrow(joined))
}
write_tsv(j_gate, file.path(results_root, "cross_lineage", "J_gate_status.tsv"))

j_gate_machine <- data.frame(
  gate = j_gate$gate,
  eligible_input_n = if (s1_zero) 0L else NA_integer_,
  evaluated_n = if (s1_zero) 0L else NA_integer_,
  pass_n = NA_integer_,
  status = j_gate$status,
  stop_rule = if (s1_zero) "STOP_NO_S2" else "not_triggered_or_pending",
  upstream_source = if (file.exists(s1_gate_path)) s1_gate_path else NA_character_,
  interpretation = j_gate$reason
)
write_tsv(j_gate_machine, file.path(results_root, "J_GATE.tsv"))

all_summaries <- do.call(rbind, lapply(analyses, `[[`, "summary"))
write_tsv(all_summaries, file.path(results_root, "analysis_summary_all_cell_lines.tsv"))

session <- capture.output(sessionInfo())
writeLines(session, file.path(results_root, "provenance", "sessionInfo.txt"))
pkg_versions <- data.frame(
  package = c("R", "edgeR", "limma", "DESeq2", "decoupleR", "fgsea", "AnnotationDbi",
              "org.Rn.eg.db", "org.Mm.eg.db", "org.Hs.eg.db", "yaml", "ggplot2"),
  version = c(as.character(getRversion()), vapply(c("edgeR", "limma", "DESeq2", "decoupleR", "fgsea",
                                                   "AnnotationDbi", "org.Rn.eg.db", "org.Mm.eg.db",
                                                   "org.Hs.eg.db", "yaml", "ggplot2"),
                                                 function(p) as.character(utils::packageVersion(p)), character(1)))
)
write_tsv(pkg_versions, file.path(results_root, "provenance", "package_versions.tsv"))

report <- c(
  "# Canonical JGW03 three-cell-line perturbation analysis",
  "",
  sprintf("Generated with seed %d. Each accession was modeled separately; no control was borrowed across accessions.", seed),
  "",
  "## Inferential boundary",
  "",
  "These analyses quantify transcriptional responses in PC12, C2C12 and HepG2 cultures. They do not establish clinical efficacy, target engagement, direct binding, or a causal mechanism. C2C12 and HepG2 are context calibrators, not independent replication of PC12.",
  "",
  "## Locked sample design",
  "",
  paste(capture.output(print(all_summaries[, c("dataset_id", "cell_line", "doses_ug_ml", "genes_filterByExpr", "corroborated_technical_failure_n", "samples_excluded_n")], row.names = FALSE)), collapse = "\n"),
  "",
  "## J-gate status",
  "",
  paste(capture.output(print(j_gate, row.names = FALSE)), collapse = "\n"),
  "",
  "The sequential STOP_NO_S2 rule is enforced when the frozen S1 input contains zero passing programs. In that state, J2 is not evaluated and must not be described as either positive or negative; all three-cell-line results are exploratory/contextual.",
  "",
  "## Important source limitation",
  "",
  "The GSE299063 GEO processing text contains an apparent Rn6 copy error. The official organism is Homo sapiens, the assembly field is hg38, and the deposited RSEM gene identifiers are ENSG; this discrepancy is reported rather than silently corrected.",
  "",
  "## Context-comparison limitation",
  "",
  "PC12 high dose is 100 ug/mL, whereas C2C12 and HepG2 high doses are 500 ug/mL. Standardized effects and correlations are therefore descriptive context diagnostics; they are not direct potency comparisons and do not alone pass or fail J3."
)
writeLines(report, file.path(results_root, "README_results.md"))

message("Completed canonical JGW03 three-cell-line analysis. Results: ", results_root)
