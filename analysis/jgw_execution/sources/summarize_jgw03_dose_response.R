#!/usr/bin/env Rscript

# Derive transparent four-level dose summaries from completed canonical JGW03
# models. These are descriptive diagnostics alongside, not replacements for,
# the prespecified edgeR ordered-dose test.

suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else
  normalizePath(file.path(getwd(), "route_b2"), mustWork = TRUE)
processed_root <- file.path(project_root, "data", "processed", "jgw03_three_cell_lines")
results_root <- file.path(project_root, "results", "jgw03_three_cell_lines")

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

summarize_matrix <- function(mat, samples, feature_label = "feature_id") {
  stopifnot(identical(colnames(mat), samples$sample_id))
  groups <- c("vehicle", "low", "medium", "high")
  means <- sapply(groups, function(g) rowMeans(mat[, samples$dose_group == g, drop = FALSE]))
  ses <- sapply(groups, function(g) {
    z <- mat[, samples$dose_group == g, drop = FALSE]
    apply(z, 1L, stats::sd) / sqrt(ncol(z))
  })
  colnames(means) <- paste0("mean_", groups)
  colnames(ses) <- paste0("SE_", groups)
  d <- t(apply(means, 1L, diff))
  overall <- means[, "mean_high"] - means[, "mean_vehicle"]
  expected_sign <- sign(overall)
  monotone_in_overall_direction <- vapply(seq_len(nrow(d)), function(i) {
    if (expected_sign[[i]] == 0) return(FALSE)
    all(d[i, ] * expected_sign[[i]] >= -1e-12)
  }, logical(1))
  strict_monotone_in_overall_direction <- vapply(seq_len(nrow(d)), function(i) {
    if (expected_sign[[i]] == 0) return(FALSE)
    all(d[i, ] * expected_sign[[i]] > 0)
  }, logical(1))
  rho <- apply(means, 1L, function(z) suppressWarnings(stats::cor(0:3, z, method = "spearman")))
  out <- data.frame(feature = rownames(mat), means, ses,
                    delta_low_minus_vehicle = d[, 1L],
                    delta_medium_minus_low = d[, 2L],
                    delta_high_minus_medium = d[, 3L],
                    high_minus_vehicle = overall,
                    spearman_group_mean_vs_ordinal = rho,
                    monotone_in_overall_direction = monotone_in_overall_direction,
                    strict_monotone_in_overall_direction = strict_monotone_in_overall_direction,
                    check.names = FALSE)
  names(out)[names(out) == "feature"] <- feature_label
  out
}

for (gse in c("GSE289929", "GSE227494", "GSE299063")) {
  obj_path <- file.path(processed_root, gse, "analysis_objects.rds")
  if (!file.exists(obj_path)) stop("Missing checkpoint: ", obj_path)
  obj <- readRDS(obj_path)
  samples <- obj$samples
  result_dir <- file.path(results_root, gse, "dose_response")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  dose_key <- unique(samples[, c("dose_group", "ordered_dose", "dose_numeric")])
  dose_key <- dose_key[order(dose_key$ordered_dose), ]
  names(dose_key)[names(dose_key) == "dose_numeric"] <- "dose_ug_ml"
  write_tsv(dose_key, file.path(result_dir, "dose_key.tsv"))

  logcpm <- edgeR::cpm(obj$y, log = TRUE, prior.count = 2)
  gene <- summarize_matrix(logcpm, samples, "raw_gene_id")
  write_tsv_gz(gene, file.path(result_dir, "gene_TMM_logCPM_four_level_diagnostics.tsv.gz"))

  reactome <- summarize_matrix(obj$reactome$scores, samples, "term_id")
  rstat <- obj$reactome$dual[, c("term_id", "term_name", "FDR_high", "FDR_trend",
                                 "high_score_effect", "high_score_CI_low", "high_score_CI_high",
                                 "trend_score_effect", "trend_score_CI_low", "trend_score_CI_high",
                                 "internal_dual_response_rule")]
  reactome <- merge(reactome, rstat, by = "term_id", all.x = TRUE, sort = FALSE)
  write_tsv(reactome, file.path(result_dir, "Reactome97_four_level_diagnostics.tsv"))

  syngo <- summarize_matrix(obj$syngo$scores, samples, "term_id")
  sstat <- obj$syngo$dual[, c("term_id", "term_name", "FDR_high", "FDR_trend",
                              "high_score_effect", "high_score_CI_low", "high_score_CI_high",
                              "trend_score_effect", "trend_score_CI_low", "trend_score_CI_high",
                              "internal_dual_response_rule")]
  syngo <- merge(syngo, sstat, by = "term_id", all.x = TRUE, sort = FALSE)
  write_tsv(syngo, file.path(result_dir, "QuickGO_SynGO_sensitivity_four_level_diagnostics.tsv"))

  add_activity_stats <- function(activity, fits) {
    z <- summarize_matrix(activity, samples, "feature")
    h <- fits$high[, c("feature", "effect", "CI_low", "CI_high", "FDR")]
    t <- fits$trend[, c("feature", "effect", "CI_low", "CI_high", "FDR")]
    names(h)[-1L] <- paste0("high_", names(h)[-1L])
    names(t)[-1L] <- paste0("trend_", names(t)[-1L])
    merge(merge(z, h, by = "feature", all.x = TRUE), t, by = "feature", all.x = TRUE)
  }
  progeny <- add_activity_stats(obj$progeny$activity, obj$progeny$fits)
  dorothea <- add_activity_stats(obj$dorothea$activity, obj$dorothea$fits)
  write_tsv(progeny, file.path(result_dir, "PROGENy_four_level_diagnostics.tsv"))
  write_tsv(dorothea, file.path(result_dir, "DoRothEA_four_level_diagnostics.tsv"))
}

message("Dose-response diagnostics completed: ", results_root)
