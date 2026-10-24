#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(file_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(digest)
})

result_root <- file.path(project_root, "results", "GSE277064")
random_seed <- 20260817L

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

break_rank_stat_ties <- function(stats) {
  stats <- stats[is.finite(stats)]
  original_ties <- sum(duplicated(stats))
  iteration <- 0L
  while (anyDuplicated(stats) && iteration < 12L) {
    iteration <- iteration + 1L
    tie_groups <- split(seq_along(stats), stats)
    tie_groups <- tie_groups[lengths(tie_groups) > 1L]
    for (indices in tie_groups) {
      ordered_indices <- indices[order(names(stats)[indices])]
      base_value <- stats[ordered_indices[1L]]
      step <- (2 ^ iteration) * .Machine$double.eps * max(1, abs(base_value))
      centered_rank <- seq_along(ordered_indices) - mean(seq_along(ordered_indices))
      stats[ordered_indices] <- base_value + centered_rank * step
    }
  }
  assert_true(!anyDuplicated(stats), "Unable to deterministically resolve preranked-statistic ties")
  attr(stats, "original_ties") <- original_ties
  attr(stats, "iterations") <- iteration
  stats
}

refresh_one <- function(cell_type, seed_offset) {
  object_path <- file.path(result_root, "objects", paste0(cell_type, "_analysis.rds"))
  object <- readRDS(object_path)
  assert_true(object$domain == "natural_course", "Unexpected analysis object domain")
  assert_true(identical(object$cell_type, cell_type), "Analysis object cell type mismatch")

  gene_stats <- fread(file.path(result_root, "gene_level", cell_type, "chronic_combined.tsv.gz"))
  stats_table <- gene_stats[!is.na(symbol) & nzchar(symbol),
                            .(rank_stat = mean(signed_sqrt_F, na.rm = TRUE)), by = symbol]
  stats_vector <- stats_table$rank_stat
  names(stats_vector) <- stats_table$symbol
  stats_vector <- break_rank_stat_ties(stats_vector)
  tie_count_before <- attr(stats_vector, "original_ties")
  tie_iterations <- attr(stats_vector, "iterations")
  stats_vector <- sort(stats_vector, decreasing = TRUE)

  pathway_symbols <- lapply(object$reactome_indices, function(index) {
    unique(na.omit(object$voom$genes$symbol[index]))
  })
  set.seed(random_seed + seed_offset)
  result <- as.data.table(fgseaMultilevel(
    pathways = pathway_symbols,
    stats = stats_vector,
    minSize = 10L,
    maxSize = 500L,
    eps = 0,
    nproc = 1
  ))
  result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
  setnames(result, c("pathway", "pval", "padj", "size"),
           c("program_id", "p_value", "q_within_celltype_contrast", "tested_feature_n"))
  result[, contrast := "chronic_combined"]
  result[, method := "fgseaMultilevel_deterministic_ties"]

  name_map <- unique(object$reactome_camera[, .(program_id, term_name)])
  result <- merge(result, name_map, by = "program_id", all.x = TRUE, sort = FALSE)
  result[, collection := "Reactome"]
  result[, cell_type := cell_type]
  result[, domain := "natural_course"]
  setcolorder(result, c(
    "program_id", "p_value", "q_within_celltype_contrast", "log2err", "ES", "NES",
    "tested_feature_n", "leadingEdge", "contrast", "method", "term_name", "collection",
    "cell_type", "domain"
  ))
  write_tsv(result, file.path(result_root, "programs", paste0(cell_type, "_Reactome_fgseaMultilevel.tsv.gz")))

  object$reactome_fgsea <- result
  object$fgsea_tie_break <- list(
    rule = "exact ties broken deterministically by gene symbol with machine-epsilon perturbations",
    ties_before = tie_count_before,
    ties_after = sum(duplicated(stats_vector)),
    iterations = tie_iterations,
    refreshed_on = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  saveRDS(object, object_path, compress = "xz")

  list(
    result = result,
    audit = data.table(
      cell_type = cell_type,
      ranked_symbols = length(stats_vector),
      exact_ties_before = tie_count_before,
      exact_ties_after = sum(duplicated(stats_vector)),
      tie_break_iterations = tie_iterations,
      pathways_tested = nrow(result),
      seed = random_seed + seed_offset
    )
  )
}

cell_types <- c("v2a", "astrocyte", "microglia")
refreshed <- Map(refresh_one, cell_types, rep(100L, length(cell_types)))
all_results <- rbindlist(lapply(refreshed, function(x) x$result), fill = TRUE)
all_results[, q_global_chronic_all_locked_precompute := p.adjust(p_value, method = "BH")]
write_tsv(all_results, file.path(result_root, "programs", "all_celltypes_Reactome_fgseaMultilevel.tsv.gz"))
write_tsv(rbindlist(lapply(refreshed, function(x) x$audit)),
          file.path(result_root, "diagnostics", "fgsea_deterministic_tie_audit.tsv"))

manifest_path <- file.path(result_root, "audit", "run_manifest.tsv")
manifest <- fread(manifest_path)
manifest[manifest_key == "script_sha256", value := digest(
  file.path(project_root, "code", "R", "analyze_GSE277064.R"),
  algo = "sha256", file = TRUE
)]
manifest <- manifest[!manifest_key %chin% c("fgsea_refresh_script_sha256", "fgsea_tie_break_rule")]
manifest <- rbind(
  manifest,
  data.table(
    manifest_key = c("fgsea_refresh_script_sha256", "fgsea_tie_break_rule"),
    value = c(
      digest(script_path, algo = "sha256", file = TRUE),
      "exact_ties_resolved_by_gene_symbol_machine_epsilon_before_fgseaMultilevel"
    )
  )
)
write_tsv(manifest, manifest_path)
writeLines(capture.output(sessionInfo()),
           file.path(result_root, "diagnostics", "fgsea_refresh_sessionInfo.txt"), useBytes = TRUE)

cat("Refreshed deterministic fgsea results for", paste(cell_types, collapse = ", "), "\n")

