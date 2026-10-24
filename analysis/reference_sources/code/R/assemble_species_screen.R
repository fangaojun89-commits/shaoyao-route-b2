#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L) stop("assemble_species_screen.R takes no arguments", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
manifest <- data.table::fread(file.path(project, "metadata", "GSE205906_sra_samples.tsv"))
data.table::setorder(manifest, sample_order)
if (nrow(manifest) != 9L || anyDuplicated(manifest$sample_id) ||
    anyDuplicated(manifest$biological_unit_id) || anyDuplicated(manifest$biosample) ||
    anyDuplicated(manifest$run_accession) ||
    !identical(as.integer(manifest$sample_order), seq_len(9L))) {
  stop("GSE205906 manifest is not the frozen nine-library one-to-one design", call. = FALSE)
}
root <- file.path(project, "data", "processed", "GSE205906", "species_screen")

records <- lapply(manifest$sample_id, function(sample_id) {
  path <- file.path(root, sample_id, "COMPLETE.tsv")
  metrics_path <- file.path(root, sample_id, "mapping_metrics_long.tsv")
  prefix_summary_path <- file.path(root, sample_id, "fastq_prefix_summary.json")
  if (!file.exists(path)) stop("Missing species screen: ", sample_id, call. = FALSE)
  if (!file.exists(metrics_path) || file.info(metrics_path)$size <= 0 ||
      !file.exists(prefix_summary_path) || file.info(prefix_summary_path)$size <= 0) {
    stop("Missing species-screen supporting output: ", sample_id, call. = FALSE)
  }
  value <- data.table::fread(path)
  if (nrow(value) != 1L) stop("Malformed species screen: ", sample_id, call. = FALSE)
  required <- c(
    "dataset_id", "sample_id", "mapping_unit", "total_reads_or_pairs",
    "rat_unique_mapping_rate", "human_unique_mapping_rate", "rat_minus_human_rate",
    "rat_to_human_rate_ratio", "threshold_rat_rate", "threshold_rate_difference",
    "threshold_rate_ratio", "status"
  )
  if (!all(required %in% names(value))) stop("Incomplete species screen: ", sample_id, call. = FALSE)
  value
})
screen <- data.table::rbindlist(records, fill = TRUE)
if (!identical(screen$sample_id, manifest$sample_id)) stop("Species-screen sample order mismatch")

config <- yaml::read_yaml(file.path(project, "config", "contrasts.yaml"))
rules <- config$datasets$GSE205906$raw_processing
threshold <- as.numeric(
  rules$species_screen_minimum_cohort_median_rate_difference
)
min_rat <- as.numeric(rules$species_screen_minimum_rat_unique_mapping_rate)
min_difference <- as.numeric(rules$species_screen_minimum_rat_minus_human_rate)
min_ratio <- as.numeric(rules$species_screen_minimum_rat_to_human_rate_ratio)
expected_prefix_n <- as.numeric(rules$species_screen_fragments_per_library)

numeric_columns <- c(
  "total_reads_or_pairs", "rat_unique_mapping_rate", "human_unique_mapping_rate",
  "rat_minus_human_rate", "rat_to_human_rate_ratio", "threshold_rat_rate",
  "threshold_rate_difference", "threshold_rate_ratio"
)
for (column in numeric_columns) data.table::set(screen, j = column, value = as.numeric(screen[[column]]))
finite_columns <- setdiff(numeric_columns, "rat_to_human_rate_ratio")
if (any(!is.finite(as.matrix(screen[, ..finite_columns]))) ||
    any(!(is.finite(screen$rat_to_human_rate_ratio) | is.infinite(screen$rat_to_human_rate_ratio))) ||
    any(screen$dataset_id != "GSE205906") || any(screen$mapping_unit != "fragments") ||
    any(!screen$status %in% c("PASS", "FAIL_SPECIES_IDENTITY")) ||
    any(screen$total_reads_or_pairs != expected_prefix_n) ||
    any(abs(screen$rat_minus_human_rate -
      (screen$rat_unique_mapping_rate - screen$human_unique_mapping_rate)) > 1e-12) ||
    any(screen$threshold_rat_rate != min_rat) ||
    any(screen$threshold_rate_difference != min_difference) ||
    any(screen$threshold_rate_ratio != min_ratio)) {
  stop("Species-screen records violate the frozen metric/provenance contract", call. = FALSE)
}
expected_ratio <- ifelse(
  screen$human_unique_mapping_rate == 0,
  Inf,
  screen$rat_unique_mapping_rate / screen$human_unique_mapping_rate
)
if (any((is.finite(expected_ratio) & abs(screen$rat_to_human_rate_ratio - expected_ratio) > 1e-12) |
        (is.infinite(expected_ratio) & !is.infinite(screen$rat_to_human_rate_ratio)))) {
  stop("Species-screen ratio is internally inconsistent", call. = FALSE)
}
recomputed_sample_pass <- screen$rat_unique_mapping_rate >= min_rat &
  screen$rat_minus_human_rate >= min_difference &
  screen$rat_to_human_rate_ratio >= min_ratio
if (any(screen$status == "PASS" & !recomputed_sample_pass)) {
  stop("A PASS species-screen marker fails the frozen sample-level rules", call. = FALSE)
}
median_difference <- stats::median(screen$rat_minus_human_rate)
all_samples_pass <- all(screen$status == "PASS")
cohort_pass <- all_samples_pass && median_difference >= threshold
cohort <- data.table(
  dataset_id = "GSE205906",
  n_libraries = nrow(screen),
  all_sample_level_rules_pass = all_samples_pass,
  median_rat_minus_human_rate = median_difference,
  threshold_median_rat_minus_human_rate = threshold,
  status = if (cohort_pass) "PASS" else "FAIL_SPECIES_IDENTITY"
)
data.table::fwrite(screen, file.path(root, "all_samples_species_screen.tsv"), sep = "\t")
data.table::fwrite(cohort, file.path(root, "COHORT_COMPLETE.tsv"), sep = "\t")
if (!cohort_pass) stop("Cohort species screen failed; GSE205906 inference suspended", call. = FALSE)
message("PASS GSE205906 cohort rat-versus-human species screen")
