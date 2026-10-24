#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 120)

required <- c(
  "affy", "affyPLM", "arrayQualityMetrics", "Biobase", "AnnotationDbi",
  "rat2302.db", "org.Rn.eg.db", "data.table", "matrixStats", "ggplot2", "jsonlite"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(affy)
  library(affyPLM)
  library(Biobase)
  library(AnnotationDbi)
  library(data.table)
  library(matrixStats)
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

raw_dir <- file.path(project_root, "data", "raw", "GSE16710")
processed_dir <- file.path(project_root, "data", "processed", "GSE16710")
result_dir <- file.path(project_root, "results", "GSE16710")
qc_dir <- file.path(result_dir, "qc")
figure_dir <- file.path(result_dir, "figures")
log_dir <- file.path(result_dir, "logs")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "preprocess_qc_GSE16710.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("GSE16710 preprocessing/QC started:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
set.seed(20260817)

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), file = path, sep = "\t", quote = FALSE, na = "NA")
}

robust_z <- function(x) {
  m <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, center = m, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

read_gpl_annotation <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  begin <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  if (is.na(begin) || is.na(end) || end <= begin + 1L) stop("Malformed GPL1355 annotation file")
  data.table::fread(text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"),
                    sep = "\t", header = TRUE, quote = "", fill = TRUE, data.table = TRUE)
}

split_multi <- function(dt, value_col, separator = "///") {
  out <- dt[, .(value = trimws(unlist(strsplit(as.character(get(value_col)), separator, fixed = TRUE)))), by = ID]
  out <- out[!is.na(value) & nzchar(value) & value != "---"]
  unique(out)
}

ledger_path <- file.path(project_root, "metadata", "samples.tsv")
ledger_all <- fread(ledger_path, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
ledger <- ledger_all[dataset_id == "GSE16710" & include_primary == TRUE]
ledger[, sample_order_int__ := as.integer(sample_order)]
setorder(ledger, sample_order_int__)
ledger[, sample_order_int__ := NULL]

expected_counts <- c(Sham_21 = 6L, Sham_60 = 5L, Spi_21 = 8L, Spi_60 = 8L)
observed_counts <- table(paste(ledger$experimental_arm, ledger$time_day, sep = "_"))
if (nrow(ledger) != 27L || !identical(as.integer(observed_counts[names(expected_counts)]), unname(expected_counts))) {
  stop("Frozen sample ledger does not match the prespecified 27-array group structure")
}
if (anyDuplicated(ledger$sample_id) || anyDuplicated(ledger$biological_unit_id)) {
  stop("Duplicate sample or biological-unit identifiers in GSE16710 ledger")
}

tar_file <- file.path(raw_dir, "GSE16710_RAW.tar")
series_file <- file.path(raw_dir, "GSE16710_series_matrix.txt.gz")
gpl_file <- file.path(raw_dir, "GPL1355.annot.gz")
stopifnot(file.exists(tar_file), file.exists(series_file), file.exists(gpl_file))

tar_members <- utils::untar(tar_file, list = TRUE)
expected_members <- ledger$source_column_or_member
if (!setequal(tar_members, expected_members)) {
  stop("CEL members in GSE16710_RAW.tar do not exactly match the frozen sample ledger")
}

temp_cel <- tempfile("GSE16710_CEL_")
dir.create(temp_cel, recursive = TRUE)
on.exit(unlink(temp_cel, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(tar_file, files = expected_members, exdir = temp_cel)
cel_files <- file.path(temp_cel, expected_members)
if (!all(file.exists(cel_files))) stop("Not all CEL files were extracted")

cat("Reading", length(cel_files), "official CEL files...\n")
abatch <- affy::ReadAffy(filenames = cel_files, sampleNames = ledger$sample_id)
if (!identical(sampleNames(abatch), ledger$sample_id)) stop("AffyBatch sample order differs from ledger")
if (annotation(abatch) != "rat2302") stop("Unexpected Affymetrix annotation/CDF: ", annotation(abatch))

cat("Computing raw-intensity metrics and RMA...\n")
pm_log2 <- log2(affy::pm(abatch))
raw_metrics <- data.table(
  sample_id = colnames(pm_log2),
  raw_log2_median = matrixStats::colMedians(pm_log2),
  raw_log2_iqr = matrixStats::colIQRs(pm_log2),
  raw_log2_mad = matrixStats::colMads(pm_log2)
)
rm(pm_log2)
gc(verbose = FALSE)

rma_eset <- affy::rma(abatch, background = TRUE, normalize = TRUE, normalize.method = "quantiles")
rma_expr <- Biobase::exprs(rma_eset)
colnames(rma_expr) <- ledger$sample_id
if (ncol(rma_expr) != 27L || any(!is.finite(rma_expr))) stop("Invalid RMA expression matrix")

cat("Computing RNA-degradation diagnostic...\n")
rna_deg <- tryCatch(affy::AffyRNAdeg(abatch), error = function(e) e)
if (inherits(rna_deg, "error")) {
  warning("AffyRNAdeg failed: ", conditionMessage(rna_deg))
  rna_deg_slope <- rep(NA_real_, ncol(rma_expr))
} else {
  deg_summary <- affy::summaryAffyRNAdeg(rna_deg)
  capture.output(str(rna_deg), file = file.path(log_dir, "AffyRNAdeg_structure.txt"))
  if (is.list(rna_deg) && "slope" %in% names(rna_deg)) {
    rna_deg_slope <- as.numeric(rna_deg$slope)
  } else if (is.list(deg_summary) && "slope" %in% names(deg_summary)) {
    rna_deg_slope <- as.numeric(deg_summary$slope)
  } else if (is.matrix(deg_summary) || is.data.frame(deg_summary)) {
    slope_col <- grep("slope", colnames(deg_summary), ignore.case = TRUE, value = TRUE)[1]
    rna_deg_slope <- if (!is.na(slope_col)) as.numeric(deg_summary[, slope_col]) else rep(NA_real_, ncol(rma_expr))
  } else {
    rna_deg_slope <- rep(NA_real_, ncol(rma_expr))
  }
  grDevices::pdf(file.path(figure_dir, "GSE16710_RNA_degradation.pdf"), width = 9, height = 7)
  try(affy::plotAffyRNAdeg(rna_deg, cols = seq_len(ncol(rma_expr))), silent = TRUE)
  grDevices::dev.off()
}

cat("Fitting probe-level models for RLE/NUSE...\n")
plm <- tryCatch(affyPLM::fitPLM(abatch, background = TRUE, normalize = TRUE, verbosity.level = 0), error = function(e) e)
if (inherits(plm, "error")) {
  warning("fitPLM failed: ", conditionMessage(plm))
  rle_mat <- nuse_mat <- matrix(NA_real_, nrow = 1L, ncol = ncol(rma_expr), dimnames = list(NULL, colnames(rma_expr)))
} else {
  rle_mat <- affyPLM::RLE(plm, type = "values")
  nuse_mat <- affyPLM::NUSE(plm, type = "values")
  grDevices::pdf(file.path(figure_dir, "GSE16710_RLE_NUSE.pdf"), width = 12, height = 8)
  par(mfrow = c(2, 1), mar = c(7, 4, 2, 1))
  boxplot(rle_mat, las = 2, outline = FALSE, main = "Relative log expression (RLE)", ylab = "RLE")
  abline(h = 0, col = "red", lty = 2)
  boxplot(nuse_mat, las = 2, outline = FALSE, main = "Normalized unscaled standard error (NUSE)", ylab = "NUSE")
  abline(h = 1, col = "red", lty = 2)
  grDevices::dev.off()
}

sample_cor <- stats::cor(rma_expr, method = "pearson")
diag(sample_cor) <- NA_real_
median_cor <- apply(sample_cor, 2L, stats::median, na.rm = TRUE)

pca <- stats::prcomp(t(rma_expr), center = TRUE, scale. = FALSE)
n_pc <- min(4L, ncol(pca$x))
pca_distance <- sqrt(rowSums(scale(pca$x[, seq_len(n_pc), drop = FALSE])^2))

qc <- merge(ledger[, .(sample_id, biological_unit_id, experimental_arm, time_day, replicate)], raw_metrics,
            by = "sample_id", sort = FALSE)
qc[, rna_degradation_slope := rna_deg_slope[match(sample_id, colnames(rma_expr))]]
qc[, rle_median := matrixStats::colMedians(rle_mat)[match(sample_id, colnames(rma_expr))]]
qc[, rle_iqr := matrixStats::colIQRs(rle_mat)[match(sample_id, colnames(rma_expr))]]
qc[, nuse_median := matrixStats::colMedians(nuse_mat)[match(sample_id, colnames(rma_expr))]]
qc[, nuse_iqr := matrixStats::colIQRs(nuse_mat)[match(sample_id, colnames(rma_expr))]]
qc[, median_interarray_correlation := median_cor[match(sample_id, names(median_cor))]]
qc[, pca_robust_distance := pca_distance[match(sample_id, names(pca_distance))]]

for (metric in c("raw_log2_median", "raw_log2_iqr", "rna_degradation_slope", "rle_median", "rle_iqr",
                 "nuse_median", "nuse_iqr", "median_interarray_correlation", "pca_robust_distance")) {
  qc[, paste0(metric, "_robust_z") := robust_z(get(metric))]
}

qc[, flag_raw_intensity := abs(raw_log2_median_robust_z) > 3.5 | abs(raw_log2_iqr_robust_z) > 3.5]
qc[, flag_degradation := abs(rna_degradation_slope_robust_z) > 3.5]
qc[, flag_probe_model := abs(rle_median_robust_z) > 3.5 | rle_iqr_robust_z > 3.5 |
                              abs(nuse_median_robust_z) > 3.5 | nuse_iqr_robust_z > 3.5]
qc[, flag_correlation := median_interarray_correlation_robust_z < -3.5]
qc[, flag_pca_only := pca_robust_distance_robust_z > 3.5]
qc[, independent_technical_domains_flagged := rowSums(.SD, na.rm = TRUE),
   .SDcols = c("flag_raw_intensity", "flag_degradation", "flag_probe_model", "flag_correlation")]
qc[, flag_for_outcome_blind_review := independent_technical_domains_flagged >= 2L]
qc[, automatic_exclusion := FALSE]
qc[, qc_note := fifelse(flag_for_outcome_blind_review,
                        "manual_review_required_before_modeling; no automatic exclusion",
                        "retain_unless_external_technical_failure_is_documented")]
qc[, sample_order_int__ := match(sample_id, ledger$sample_id)]
setorder(qc, sample_order_int__)
qc[, sample_order_int__ := NULL]
write_tsv(qc, file.path(qc_dir, "array_qc_metrics_and_flags.tsv"))
write_tsv(as.data.table(sample_cor, keep.rownames = "sample_id"), file.path(qc_dir, "sample_correlation_matrix.tsv"))

pca_dt <- data.table(sample_id = rownames(pca$x), pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE])
pca_dt <- merge(pca_dt, ledger[, .(sample_id, experimental_arm, time_day)], by = "sample_id", sort = FALSE)
write_tsv(pca_dt, file.path(qc_dir, "pca_scores.tsv"))
var_exp <- pca$sdev^2 / sum(pca$sdev^2)
write_tsv(data.table(PC = paste0("PC", seq_along(var_exp)), variance_explained = var_exp,
                     cumulative_variance = cumsum(var_exp)), file.path(qc_dir, "pca_variance.tsv"))

p_pca <- ggplot(pca_dt, aes(x = PC1, y = PC2, colour = experimental_arm, shape = factor(time_day), label = sample_id)) +
  geom_point(size = 3) + geom_text(size = 2.4, vjust = -0.8, check_overlap = TRUE) +
  theme_bw(base_size = 11) + labs(shape = "Day", colour = "Injury", title = "GSE16710 RMA PCA")
ggsave(file.path(figure_dir, "GSE16710_RMA_PCA.pdf"), p_pca, width = 8, height = 6)
ggsave(file.path(figure_dir, "GSE16710_RMA_PCA.png"), p_pca, width = 8, height = 6, dpi = 300)

cor_long <- as.data.table(as.table(sample_cor))
setnames(cor_long, c("sample_1", "sample_2", "correlation"))
p_cor <- ggplot(cor_long, aes(sample_1, sample_2, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c(limits = c(min(sample_cor, na.rm = TRUE), 1), na.value = "white") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank()) +
  labs(title = "GSE16710 inter-array correlation")
ggsave(file.path(figure_dir, "GSE16710_sample_correlation.pdf"), p_cor, width = 9, height = 8)

cat("Building current Bioconductor probe annotation...\n")
probe_ids <- rownames(rma_expr)
current_map <- AnnotationDbi::select(rat2302.db::rat2302.db, keys = probe_ids,
                                     keytype = "PROBEID", columns = c("ENTREZID"))
current_map <- as.data.table(current_map)
setnames(current_map, c("PROBEID", "ENTREZID"), c("probe_id", "entrez_id"))
current_map <- current_map[!is.na(entrez_id) & nzchar(entrez_id)]
current_map[, gene_symbol := AnnotationDbi::mapIds(org.Rn.eg.db::org.Rn.eg.db,
                                                   keys = entrez_id, keytype = "ENTREZID",
                                                   column = "SYMBOL", multiVals = "first")]
current_map <- unique(current_map[!is.na(gene_symbol) & nzchar(gene_symbol)])
current_map[, annotation_source := paste0("rat2302.db_", as.character(packageVersion("rat2302.db")),
                                          "+org.Rn.eg.db_", as.character(packageVersion("org.Rn.eg.db")))]
write_tsv(current_map, file.path(processed_dir, "probe_annotation_current_bioconductor.tsv"))

cat("Parsing official GEO GPL1355 annotation (2016 snapshot)...\n")
gpl <- read_gpl_annotation(gpl_file)
if (!all(c("ID", "Gene symbol", "Gene ID") %in% names(gpl))) stop("Required GPL1355 annotation columns absent")
original_symbol_map <- split_multi(gpl[, .(ID, `Gene symbol`)], "Gene symbol")
setnames(original_symbol_map, c("ID", "value"), c("probe_id", "gene_symbol"))
original_symbol_map[, annotation_source := "GEO_GPL1355_annotation_2016-08-09"]
original_symbol_map <- unique(original_symbol_map[probe_id %in% probe_ids])
write_tsv(original_symbol_map, file.path(processed_dir, "probe_annotation_original_GPL1355_2016.tsv"))

probe_mean <- rowMeans(rma_expr)
make_gene_matrix <- function(map_dt, label) {
  x <- copy(map_dt)[probe_id %in% rownames(rma_expr)]
  x[, probe_mean_expression := probe_mean[match(probe_id, names(probe_mean))]]
  setorder(x, gene_symbol, -probe_mean_expression, probe_id)
  rep_map <- x[, .SD[1L], by = gene_symbol]
  mat <- rma_expr[rep_map$probe_id, , drop = FALSE]
  rownames(mat) <- rep_map$gene_symbol
  write_tsv(rep_map, file.path(processed_dir, paste0("gene_representative_probe_", label, ".tsv")))
  saveRDS(mat, file.path(processed_dir, paste0("rma_gene_expression_", label, ".rds")), compress = "xz")
  write_tsv(as.data.table(mat, keep.rownames = "gene_symbol"),
            file.path(processed_dir, paste0("rma_gene_expression_", label, ".tsv.gz")))
  list(map = x, representative = rep_map, matrix = mat)
}

current_gene <- make_gene_matrix(current_map, "current_bioconductor")
original_gene <- make_gene_matrix(original_symbol_map, "original_GPL1355_2016")

saveRDS(rma_expr, file.path(processed_dir, "rma_probe_expression.rds"), compress = "xz")
write_tsv(as.data.table(rma_expr, keep.rownames = "probe_id"),
          file.path(processed_dir, "rma_probe_expression.tsv.gz"))
write_tsv(ledger, file.path(processed_dir, "sample_ledger_verified.tsv"))

cat("Comparing raw-CEL RMA with author-submitted processed matrix...\n")
series_lines <- readLines(gzfile(series_file), warn = FALSE)
tab_begin <- match("!series_matrix_table_begin", series_lines)
tab_end <- match("!series_matrix_table_end", series_lines)
if (is.na(tab_begin) || is.na(tab_end)) stop("Malformed GSE16710 series matrix")
series_dt <- fread(text = paste(series_lines[(tab_begin + 1L):(tab_end - 1L)], collapse = "\n"),
                   sep = "\t", header = TRUE, quote = "\"", data.table = TRUE)
setnames(series_dt, 1L, "probe_id")
series_samples <- setdiff(names(series_dt), "probe_id")
if (!identical(series_samples, ledger$sample_id)) stop("Series matrix sample order differs from ledger")
common_probe <- intersect(series_dt$probe_id, rownames(rma_expr))
series_mat <- as.matrix(series_dt[match(common_probe, probe_id), ..series_samples])
storage.mode(series_mat) <- "double"
rownames(series_mat) <- common_probe
raw_rma_common <- rma_expr[common_probe, series_samples, drop = FALSE]
series_concordance <- data.table(
  sample_id = series_samples,
  pearson_r = vapply(series_samples, function(s) cor(series_mat[, s], raw_rma_common[, s]), numeric(1)),
  median_absolute_difference = vapply(series_samples, function(s) median(abs(series_mat[, s] - raw_rma_common[, s])), numeric(1))
)
write_tsv(series_concordance, file.path(qc_dir, "author_series_vs_raw_RMA_concordance.tsv"))

cat("Generating arrayQualityMetrics report (failure is logged but does not stop core QC)...\n")
pData(rma_eset) <- data.frame(
  sample_id = ledger$sample_id,
  group = paste(ledger$experimental_arm, ledger$time_day, sep = "_"),
  row.names = ledger$sample_id,
  stringsAsFactors = FALSE
)
aqm_dir <- file.path(qc_dir, "arrayQualityMetrics_report")
aqm_error <- NULL
tryCatch(
  arrayQualityMetrics::arrayQualityMetrics(
    expressionset = rma_eset,
    outdir = aqm_dir,
    force = TRUE,
    do.logtransform = FALSE,
    intgroup = "group"
  ),
  error = function(e) aqm_error <<- conditionMessage(e)
)
writeLines(if (is.null(aqm_error)) "PASS" else paste("FAILED:", aqm_error),
           file.path(qc_dir, "arrayQualityMetrics_status.txt"))

source_files <- data.table(
  file = c("GSE16710_RAW.tar", "GSE16710_series_matrix.txt.gz", "GPL1355.annot.gz"),
  bytes = file.info(c(tar_file, series_file, gpl_file))$size,
  md5 = unname(tools::md5sum(c(tar_file, series_file, gpl_file)))
)
write_tsv(source_files, file.path(processed_dir, "source_file_manifest.tsv"))

package_names <- unique(c(required, "rat2302cdf"))
package_versions <- data.table(
  package = package_names,
  version = vapply(package_names, function(x) as.character(packageVersion(x)), character(1))
)
write_tsv(package_versions, file.path(processed_dir, "preprocessing_package_versions.tsv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_preprocess_qc.txt"))

summary_obj <- list(
  dataset = "GSE16710",
  stage = "outcome_blind_preprocessing_and_qc",
  n_arrays = ncol(rma_expr),
  n_probe_sets = nrow(rma_expr),
  group_counts = as.list(as.integer(observed_counts[names(expected_counts)])),
  group_count_names = names(expected_counts),
  n_current_annotated_probe_gene_pairs = nrow(current_map),
  n_current_representative_genes = nrow(current_gene$matrix),
  n_original_annotated_probe_gene_pairs = nrow(original_symbol_map),
  n_original_representative_genes = nrow(original_gene$matrix),
  n_arrays_requiring_outcome_blind_review = sum(qc$flag_for_outcome_blind_review),
  automatic_exclusions = 0L,
  arrayQualityMetrics_status = if (is.null(aqm_error)) "PASS" else paste("FAILED", aqm_error),
  completed_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
)
jsonlite::write_json(summary_obj, file.path(qc_dir, "preprocessing_qc_summary.json"), pretty = TRUE, auto_unbox = TRUE)

cat("GSE16710 preprocessing/QC completed:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
