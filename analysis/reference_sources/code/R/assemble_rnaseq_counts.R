#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("GSE128623", "GSE205906")) {
  stop("Usage: assemble_rnaseq_counts.R <GSE128623|GSE205906>", call. = FALSE)
}
dataset_id <- args[[1L]]

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)

manifest_path <- file.path(project, "metadata", paste0(dataset_id, "_sra_samples.tsv"))
input_root <- file.path(project, "data", "processed", dataset_id, "raw_realign")
output_dir <- file.path(project, "data", "processed", dataset_id, "assembled_counts")
if (!file.exists(manifest_path)) stop("Missing sample manifest", call. = FALSE)
manifest <- data.table::fread(manifest_path)
data.table::setorder(manifest, sample_order)
config <- yaml::read_yaml(file.path(project, "config", "contrasts.yaml"))
raw_rules <- config$datasets[[dataset_id]]$raw_processing
dominant_mode_fraction_min <- as.numeric(raw_rules$strandedness_same_dominant_mode_min_fraction)
dominant_over_unstranded_min <- as.numeric(raw_rules$strandedness_median_dominant_over_unstranded_min)
minor_over_dominant_max <- as.numeric(raw_rules$strandedness_median_minor_over_dominant_max)
if (any(!is.finite(c(
      dominant_mode_fraction_min, dominant_over_unstranded_min, minor_over_dominant_max
    ))) || dominant_mode_fraction_min < 0 || dominant_mode_fraction_min > 1 ||
    dominant_over_unstranded_min < 0 || dominant_over_unstranded_min > 1 ||
    minor_over_dominant_max < 0 || minor_over_dominant_max > 1) {
  stop("Missing/invalid frozen strandedness thresholds", call. = FALSE)
}
expected_n <- if (dataset_id == "GSE128623") 30L else 9L
expected_layout <- if (dataset_id == "GSE128623") "SINGLE" else "PAIRED"
expected_organism <- if (dataset_id == "GSE128623") "Mus musculus" else "Rattus norvegicus"
required_manifest_columns <- c(
  "dataset_id", "sample_order", "sample_id", "biological_unit_id", "statistical_unit",
  "organism", "experimental_arm", "condition", "time_day", "biosample", "sra_experiment",
  "run_accession", "library_layout", "library_strategy", "include_primary"
)
if (!all(required_manifest_columns %in% names(manifest)) || nrow(manifest) != expected_n ||
    anyDuplicated(manifest$sample_id) || anyDuplicated(manifest$biological_unit_id) ||
    anyDuplicated(manifest$biosample) || anyDuplicated(manifest$sra_experiment) ||
    anyDuplicated(manifest$run_accession) ||
    !identical(as.integer(manifest$sample_order), seq_len(expected_n)) ||
    any(manifest$dataset_id != dataset_id) || any(manifest$library_layout != expected_layout) ||
    any(manifest$organism != expected_organism) || any(manifest$library_strategy != "RNA-Seq")) {
  stop("Sample manifest violates the frozen one-row-per-biological-library design", call. = FALSE)
}

if (dataset_id == "GSE128623") {
  observed_groups <- manifest[, .N, by = .(experimental_arm, time_day)]
  expected_groups <- data.table(
    experimental_arm = c("control", rep("stroke", 5L)),
    time_day = c(0L, 4L, 7L, 14L, 28L, 42L),
    N = rep(5L, 6L)
  )
  data.table::setorder(observed_groups, time_day)
  data.table::setorder(expected_groups, time_day)
  if (!identical(observed_groups, expected_groups) || any(manifest$include_primary != "TRUE")) {
    stop("GSE128623 manifest groups differ from the frozen 6 x 5 design", call. = FALSE)
  }
} else {
  observed_groups <- manifest[, .N, by = experimental_arm]
  data.table::setorder(observed_groups, experimental_arm)
  expected_groups <- data.table(experimental_arm = c("JA", "PLM", "SNI"), N = rep(3L, 3L))
  if (!identical(observed_groups, expected_groups) || any(manifest$include_primary != "CONDITIONAL")) {
    stop("GSE205906 manifest groups differ from the frozen 3 x 3 design", call. = FALSE)
  }
  cohort_path <- file.path(project, "data", "processed", "GSE205906", "species_screen", "COHORT_COMPLETE.tsv")
  if (!file.exists(cohort_path)) {
    stop("GSE205906 count assembly is prohibited before the cohort species gate", call. = FALSE)
  }
  cohort <- data.table::fread(cohort_path)
  threshold <- as.numeric(
    config$datasets$GSE205906$raw_processing$species_screen_minimum_cohort_median_rate_difference
  )
  if (nrow(cohort) != 1L || cohort$dataset_id[[1L]] != "GSE205906" ||
      cohort$n_libraries[[1L]] != expected_n || cohort$status[[1L]] != "PASS" ||
      !isTRUE(as.logical(cohort$all_sample_level_rules_pass[[1L]])) ||
      !is.finite(cohort$median_rat_minus_human_rate[[1L]]) ||
      cohort$threshold_median_rat_minus_human_rate[[1L]] != threshold ||
      cohort$median_rat_minus_human_rate[[1L]] < threshold) {
    stop("GSE205906 cohort species gate is missing, malformed, or non-PASS", call. = FALSE)
  }
}

assembly_marker <- file.path(output_dir, "COMPLETE.tsv")
expected_assembly_outputs <- c(
  "assignment_qc_long.tsv", "assignment_qc_by_sample.tsv", "strand_decision.tsv",
  "gene_annotation.tsv.gz", "gene_counts_primary.tsv.gz", "library_qc.tsv",
  "gene_counts_primary.rds"
)
if (file.exists(assembly_marker)) {
  completed <- data.table::fread(assembly_marker)
  output_manifest_path <- file.path(output_dir, "OUTPUT_MANIFEST.tsv")
  if (nrow(completed) != 1L || completed$dataset_id[[1L]] != dataset_id ||
      completed$status[[1L]] != "PASS" || !file.exists(output_manifest_path)) {
    stop("Existing assembly completion state is malformed", call. = FALSE)
  }
  output_manifest <- data.table::fread(output_manifest_path)
  if (nrow(output_manifest) != length(expected_assembly_outputs) ||
      anyDuplicated(output_manifest$filename) ||
      !setequal(output_manifest$filename, expected_assembly_outputs)) {
    stop("Existing assembly output manifest has the wrong file set", call. = FALSE)
  }
  for (i in seq_len(nrow(output_manifest))) {
    filename <- output_manifest$filename[[i]]
    if (basename(filename) != filename || !grepl("^[0-9a-fA-F]{32}$", output_manifest$md5[[i]])) {
      stop("Unsafe/malformed existing assembly output manifest", call. = FALSE)
    }
    path <- file.path(output_dir, filename)
    if (!file.exists(path) || file.info(path)$size != output_manifest$bytes[[i]] ||
        unname(tools::md5sum(path)) != tolower(output_manifest$md5[[i]])) {
      stop("Existing assembly output fails byte/MD5 validation: ", path, call. = FALSE)
    }
  }
  message("PASS existing validated assembly: ", dataset_id)
  quit(save = "no", status = 0L)
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Partial unmarked assembly exists; preserve/archive before retrying: ", output_dir, call. = FALSE)
}

qc_rows <- list()
for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[[i]]
  sample_dir <- file.path(input_root, sample_id)
  marker <- file.path(sample_dir, "COMPLETE.tsv")
  if (!file.exists(marker)) stop("Incomplete sample: ", sample_id, call. = FALSE)
  completed <- data.table::fread(marker)
  if (nrow(completed) != 1L || completed$status[[1L]] != "PASS") {
    stop("Invalid completion marker: ", marker, call. = FALSE)
  }
  for (mode in 0:2) {
    stat_path <- file.path(sample_dir, paste0("assignment_strand", mode, ".tsv"))
    if (!file.exists(stat_path)) stop("Missing assignment summary: ", stat_path, call. = FALSE)
    stat <- data.table::fread(stat_path)
    if (!"Status" %in% names(stat) || ncol(stat) < 2L) stop("Malformed: ", stat_path, call. = FALSE)
    count_col <- names(stat)[[2L]]
    stat_values <- as.numeric(stat[[count_col]])
    if (any(!is.finite(stat_values)) || any(stat_values < 0) ||
        any(stat_values != floor(stat_values))) {
      stop("Invalid assignment counts: ", stat_path, call. = FALSE)
    }
    assigned <- stat[Status == "Assigned", get(count_col)]
    if (length(assigned) != 1L || !is.finite(assigned)) stop("Invalid Assigned row: ", stat_path)
    total <- sum(stat_values)
    if (total <= 0 || assigned > total) stop("Invalid assignment total: ", stat_path, call. = FALSE)
    qc_rows[[length(qc_rows) + 1L]] <- data.table(
      dataset_id = dataset_id,
      sample_order = manifest$sample_order[[i]],
      sample_id = sample_id,
      strand_mode = mode,
      assigned = as.numeric(assigned),
      total = as.numeric(total),
      assigned_fraction = as.numeric(assigned) / as.numeric(total)
    )
  }
}
strand_qc <- data.table::rbindlist(qc_rows)
wide <- data.table::dcast(strand_qc, sample_order + sample_id ~ strand_mode, value.var = "assigned")
data.table::setnames(wide, c("0", "1", "2"), c("assigned_0", "assigned_1", "assigned_2"))
wide[, dominant_mode := fifelse(assigned_1 > assigned_2, 1L, fifelse(assigned_2 > assigned_1, 2L, NA_integer_))]
wide[, dominant_assigned := pmax(assigned_1, assigned_2)]
wide[, minor_assigned := pmin(assigned_1, assigned_2)]
if (any(wide$assigned_0 <= 0) || any(wide$dominant_assigned <= 0)) {
  stop("Strandedness diagnostics contain zero denominators", call. = FALSE)
}
wide[, dominant_over_unstranded := dominant_assigned / assigned_0]
wide[, minor_over_dominant := minor_assigned / dominant_assigned]

mode_frequencies <- table(wide$dominant_mode, useNA = "no")
if (length(mode_frequencies) == 0L) stop("No non-tied strand diagnostics", call. = FALSE)
candidate_mode <- as.integer(names(which.max(mode_frequencies)))
dominance_fraction <- as.numeric(max(mode_frequencies)) / nrow(wide)
median_dominant_unstranded <- stats::median(wide$dominant_over_unstranded, na.rm = TRUE)
median_minor_dominant <- stats::median(wide$minor_over_dominant, na.rm = TRUE)
if (!is.finite(dominance_fraction) || !is.finite(median_dominant_unstranded) ||
    !is.finite(median_minor_dominant)) {
  stop("Strandedness decision metrics are not finite", call. = FALSE)
}
selected_mode <- if (
  dominance_fraction >= dominant_mode_fraction_min &&
  median_dominant_unstranded >= dominant_over_unstranded_min &&
  median_minor_dominant <= minor_over_dominant_max
) candidate_mode else 0L

decision <- data.table(
  dataset_id = dataset_id,
  n_libraries = nrow(wide),
  candidate_stranded_mode = candidate_mode,
  dominant_mode_fraction = dominance_fraction,
  median_dominant_over_unstranded = median_dominant_unstranded,
  median_minor_over_dominant = median_minor_dominant,
  threshold_dominant_mode_fraction = dominant_mode_fraction_min,
  threshold_dominant_over_unstranded = dominant_over_unstranded_min,
  threshold_minor_over_dominant = minor_over_dominant_max,
  selected_primary_strand_mode = selected_mode,
  decision = if (selected_mode == 0L) "UNSTRANDED" else paste0("STRANDED_MODE_", selected_mode),
  labels_used_for_decision = FALSE,
  status = "PASS"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(strand_qc, file.path(output_dir, "assignment_qc_long.tsv"), sep = "\t")
data.table::fwrite(wide, file.path(output_dir, "assignment_qc_by_sample.tsv"), sep = "\t")
data.table::fwrite(decision, file.path(output_dir, "strand_decision.tsv"), sep = "\t")

annotation <- NULL
gene_ids <- NULL
count_columns <- list()
library_qc <- list()
for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[[i]]
  count_path <- file.path(
    input_root,
    sample_id,
    paste0("counts_strand", selected_mode, ".tsv.gz")
  )
  counts <- data.table::fread(count_path)
  if (!all(c("GeneID", "count") %in% names(counts))) {
    stop("Count table lacks GeneID/count: ", count_path, call. = FALSE)
  }
  if (anyDuplicated(counts$GeneID)) stop("Duplicated GeneID in ", count_path, call. = FALSE)
  if (i == 1L) {
    gene_ids <- counts$GeneID
    annotation_columns <- setdiff(names(counts), "count")
    annotation <- counts[, ..annotation_columns]
  } else if (!identical(counts$GeneID, gene_ids)) {
    stop("Gene order/identity differs for sample ", sample_id, call. = FALSE)
  } else if (!identical(as.data.frame(counts[, ..annotation_columns]), as.data.frame(annotation))) {
    stop("Gene annotation differs for sample ", sample_id, call. = FALSE)
  }
  values <- as.numeric(counts$count)
  if (any(!is.finite(values)) || any(values < 0) || any(values != floor(values))) {
    stop("Non-finite, negative, or non-integer counts in ", count_path, call. = FALSE)
  }
  count_columns[[sample_id]] <- values
  expected_assigned <- wide[
    wide$sample_id == manifest$sample_id[[i]],
    get(paste0("assigned_", selected_mode))
  ]
  if (length(expected_assigned) != 1L || sum(values) != expected_assigned) {
    stop("Count-table sum differs from featureCounts Assigned total for ", sample_id, call. = FALSE)
  }
  library_qc[[i]] <- data.table(
    dataset_id = dataset_id,
    sample_order = manifest$sample_order[[i]],
    sample_id = sample_id,
    selected_strand_mode = selected_mode,
    assigned_gene_counts = sum(values),
    detected_genes_count_ge_1 = sum(values >= 1),
    detected_genes_count_ge_10 = sum(values >= 10)
  )
}

count_matrix <- do.call(cbind, count_columns)
rownames(count_matrix) <- gene_ids
storage.mode(count_matrix) <- "integer"
count_table <- data.table(GeneID = gene_ids)
for (sample_id in names(count_columns)) count_table[[sample_id]] <- count_columns[[sample_id]]

data.table::fwrite(annotation, file.path(output_dir, "gene_annotation.tsv.gz"), sep = "\t")
data.table::fwrite(count_table, file.path(output_dir, "gene_counts_primary.tsv.gz"), sep = "\t")
data.table::fwrite(data.table::rbindlist(library_qc), file.path(output_dir, "library_qc.tsv"), sep = "\t")
saveRDS(
  list(
    dataset_id = dataset_id,
    counts = count_matrix,
    gene_annotation = annotation,
    samples = manifest,
    strand_decision = decision,
    assignment_qc = wide
  ),
  file.path(output_dir, "gene_counts_primary.rds"),
  compress = "xz"
)

output_files <- list.files(output_dir, full.names = TRUE)
output_files <- output_files[!basename(output_files) %in% c("OUTPUT_MANIFEST.tsv", "COMPLETE.tsv")]
if (!setequal(basename(output_files), expected_assembly_outputs)) {
  stop("Assembly produced an unexpected output file set", call. = FALSE)
}
output_manifest <- data.table(
  filename = basename(output_files),
  bytes = file.info(output_files)$size,
  md5 = unname(tools::md5sum(output_files))
)
data.table::fwrite(output_manifest, file.path(output_dir, "OUTPUT_MANIFEST.tsv"), sep = "\t")
completion <- data.table(
  dataset_id = dataset_id,
  n_libraries = nrow(manifest),
  n_genes = nrow(count_matrix),
  selected_primary_strand_mode = selected_mode,
  status = "PASS"
)
data.table::fwrite(completion, assembly_marker, sep = "\t")
message("PASS assembled ", dataset_id, " using strand mode ", selected_mode)
