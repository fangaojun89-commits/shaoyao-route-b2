#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rsubread)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    paste(
      "Usage: species_screen_one_sample.R <sample_id> <SINGLE|PAIRED>",
      "<rat_index_prefix> <human_index_prefix> <prefix_fastq1.gz>",
      "<prefix_fastq2.gz|-> <output_dir> <bam_work_dir>"
    ),
    call. = FALSE
  )
}
sample_id <- args[[1L]]
layout <- args[[2L]]
rat_index <- args[[3L]]
human_index <- args[[4L]]
fastq1 <- args[[5L]]
fastq2 <- args[[6L]]
output_dir <- args[[7L]]
bam_work_dir <- args[[8L]]

if (!grepl("^GSM[0-9]+$", sample_id)) stop("Invalid sample ID", call. = FALSE)
if (!layout %in% c("SINGLE", "PAIRED")) stop("Invalid layout", call. = FALSE)
if (!file.exists(fastq1)) stop("Missing diagnostic FASTQ1", call. = FALSE)
if (layout == "PAIRED" && !file.exists(fastq2)) stop("Missing diagnostic FASTQ2", call. = FALSE)
if (layout == "SINGLE" && fastq2 != "-") stop("Single-end screen must use '-' for FASTQ2", call. = FALSE)

validate_index <- function(prefix) {
  files <- list.files(dirname(prefix), pattern = paste0("^", basename(prefix), "\\."), full.names = TRUE)
  core <- files[grepl("(\\.files$|\\.reads$|\\.[0-9]+\\.b\\.(array|tab)$)", files)]
  if (length(core) < 4L || any(file.info(core)$size <= 0)) {
    stop("Missing/incomplete index: ", prefix, call. = FALSE)
  }
}
validate_index(rat_index)
validate_index(human_index)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(project, "config", "contrasts.yaml"))
rules <- config$datasets$GSE205906$raw_processing
min_rat <- as.numeric(rules$species_screen_minimum_rat_unique_mapping_rate)
min_difference <- as.numeric(rules$species_screen_minimum_rat_minus_human_rate)
min_ratio <- as.numeric(rules$species_screen_minimum_rat_to_human_rate_ratio)

marker <- file.path(output_dir, "COMPLETE.tsv")
if (file.exists(marker)) {
  old <- utils::read.delim(marker, check.names = FALSE)
  required_columns <- c(
    "dataset_id", "sample_id", "diagnostic_method", "mapping_unit", "total_reads_or_pairs",
    "rat_unique_mapping_rate", "human_unique_mapping_rate",
    "rat_minus_human_rate", "rat_to_human_rate_ratio",
    "threshold_rat_rate", "threshold_rate_difference", "threshold_rate_ratio", "status"
  )
  valid_old <- nrow(old) == 1L && all(required_columns %in% names(old))
  if (valid_old) {
    numeric_values <- unlist(old[1L, c(
      "total_reads_or_pairs", "rat_unique_mapping_rate", "human_unique_mapping_rate",
      "rat_minus_human_rate", "rat_to_human_rate_ratio", "threshold_rat_rate",
      "threshold_rate_difference", "threshold_rate_ratio"
    )], use.names = TRUE)
    numeric_values <- suppressWarnings(as.numeric(numeric_values))
    names(numeric_values) <- c(
      "total", "rat_rate", "human_rate", "difference", "ratio",
      "recorded_min_rat", "recorded_min_difference", "recorded_min_ratio"
    )
    finite_required <- all(is.finite(numeric_values[c(
      "total", "rat_rate", "human_rate", "difference",
      "recorded_min_rat", "recorded_min_difference", "recorded_min_ratio"
    )]))
    ratio_valid <- is.finite(numeric_values[["ratio"]]) || is.infinite(numeric_values[["ratio"]])
    thresholds_match <- isTRUE(all.equal(numeric_values[["recorded_min_rat"]], min_rat)) &&
      isTRUE(all.equal(numeric_values[["recorded_min_difference"]], min_difference)) &&
      isTRUE(all.equal(numeric_values[["recorded_min_ratio"]], min_ratio))
    recomputed_pass <- finite_required && ratio_valid && numeric_values[["total"]] > 0 &&
      numeric_values[["rat_rate"]] >= min_rat &&
      numeric_values[["difference"]] >= min_difference &&
      numeric_values[["ratio"]] >= min_ratio
    expected_ratio <- if (numeric_values[["human_rate"]] == 0) {
      Inf
    } else {
      numeric_values[["rat_rate"]] / numeric_values[["human_rate"]]
    }
    metrics_coherent <- abs(numeric_values[["difference"]] -
      (numeric_values[["rat_rate"]] - numeric_values[["human_rate"]])) <= 1e-12 &&
      ((is.finite(expected_ratio) &&
        abs(numeric_values[["ratio"]] - expected_ratio) <= 1e-12) ||
        (is.infinite(expected_ratio) && is.infinite(numeric_values[["ratio"]])))
    valid_old <- old$dataset_id[[1L]] == "GSE205906" &&
      old$sample_id[[1L]] == sample_id &&
      old$diagnostic_method[[1L]] == "same_deterministic_prefix_dual_reference_unique_Subjunc" &&
      old$mapping_unit[[1L]] == (if (layout == "PAIRED") "fragments" else "reads") &&
      old$status[[1L]] == "PASS" && thresholds_match && metrics_coherent && recomputed_pass &&
      file.exists(file.path(output_dir, "mapping_metrics_long.tsv")) &&
      file.info(file.path(output_dir, "mapping_metrics_long.tsv"))$size > 0
  }
  if (isTRUE(valid_old)) {
    message("PASS existing species screen: ", sample_id)
    quit(save = "no", status = 0L)
  }
  stop("Existing species screen marker/supporting output is invalid: ", marker, call. = FALSE)
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Partial species-screen output exists: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bam_work_dir, recursive = TRUE, showWarnings = FALSE)

map_one <- function(label, index) {
  bam <- file.path(bam_work_dir, paste0(sample_id, ".species_", label, ".bam"))
  result <- Rsubread::subjunc(
    index = index,
    readfile1 = fastq1,
    readfile2 = if (layout == "PAIRED") fastq2 else NULL,
    input_format = "gzFASTQ",
    output_format = "BAM",
    output_file = bam,
    phredOffset = 33,
    unique = TRUE,
    nBestLocations = 1,
    nthreads = 4,
    sortReadsByCoordinates = FALSE,
    useAnnotation = FALSE
  )
  if (!file.exists(bam) || file.info(bam)$size <= 0) stop("Missing ", label, " screen BAM")
  values <- as.numeric(result[, 1L])
  names(values) <- rownames(result)
  metric_unit <- if (layout == "PAIRED") "fragments" else "reads"
  required <- paste0(c("Total_", "Mapped_", "Uniquely_mapped_"), metric_unit)
  if (!all(required %in% names(values))) stop("Unexpected Subjunc mapping summary")
  selected <- values[required]
  names(selected) <- c("total", "mapped", "unique")
  if (any(!is.finite(selected)) || selected[["total"]] <= 0 ||
      any(selected < 0) || selected[["mapped"]] > selected[["total"]] ||
      selected[["unique"]] > selected[["mapped"]]) {
    stop("Invalid ", label, " Subjunc mapping totals", call. = FALSE)
  }
  list(label = label, bam = bam, values = values, selected = selected, metric_unit = metric_unit)
}

rat <- map_one("rat", rat_index)
human <- map_one("human", human_index)
if (rat$selected[["total"]] != human$selected[["total"]]) {
  stop("Rat and human screens used different read totals", call. = FALSE)
}
rat_rate <- rat$selected[["unique"]] / rat$selected[["total"]]
human_rate <- human$selected[["unique"]] / human$selected[["total"]]
difference <- rat_rate - human_rate
ratio <- if (human_rate == 0) Inf else rat_rate / human_rate
pass <- rat_rate >= min_rat && difference >= min_difference && ratio >= min_ratio

long <- do.call(
  rbind,
  lapply(list(rat, human), function(x) {
    data.frame(
      sample_id = sample_id,
      reference_species = x$label,
      metric = names(x$values),
      value = as.numeric(x$values),
      stringsAsFactors = FALSE
    )
  })
)
utils::write.table(long, file.path(output_dir, "mapping_metrics_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

record <- data.frame(
  dataset_id = "GSE205906",
  sample_id = sample_id,
  diagnostic_method = "same_deterministic_prefix_dual_reference_unique_Subjunc",
  mapping_unit = rat$metric_unit,
  total_reads_or_pairs = rat$selected[["total"]],
  rat_unique_mapping_rate = rat_rate,
  human_unique_mapping_rate = human_rate,
  rat_minus_human_rate = difference,
  rat_to_human_rate_ratio = ratio,
  threshold_rat_rate = min_rat,
  threshold_rate_difference = min_difference,
  threshold_rate_ratio = min_ratio,
  status = if (pass) "PASS" else "FAIL_SPECIES_IDENTITY",
  stringsAsFactors = FALSE
)
utils::write.table(record, marker, sep = "\t", quote = FALSE, row.names = FALSE)
if (!pass) {
  stop("Species identity screen failed for ", sample_id, "; GSE205906 inference suspended", call. = FALSE)
}
message("PASS rat-versus-human species screen: ", sample_id)
