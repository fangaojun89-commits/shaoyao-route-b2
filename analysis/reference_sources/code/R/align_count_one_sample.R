#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Rsubread))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(
    paste(
      "Usage: align_count_one_sample.R",
      "<dataset_id> <sample_id> <SINGLE|PAIRED> <index_prefix>",
      "<annotation.gtf.gz> <fastq1.gz> <fastq2.gz|-> <sample_output_dir>",
      "<bam_work_dir>"
    ),
    call. = FALSE
  )
}

dataset_id <- args[[1L]]
sample_id <- args[[2L]]
layout <- args[[3L]]
index_prefix <- args[[4L]]
gtf <- args[[5L]]
fastq1 <- args[[6L]]
fastq2 <- args[[7L]]
output_dir <- args[[8L]]
bam_work_dir <- args[[9L]]

if (!grepl("^GSE[0-9]+$", dataset_id)) stop("Invalid dataset ID", call. = FALSE)
if (!grepl("^GSM[0-9]+$", sample_id)) stop("Invalid sample ID", call. = FALSE)
if (!layout %in% c("SINGLE", "PAIRED")) stop("Invalid library layout", call. = FALSE)
if (!file.exists(gtf)) stop("Missing GTF: ", gtf, call. = FALSE)
if (!file.exists(fastq1)) stop("Missing FASTQ1: ", fastq1, call. = FALSE)
if (layout == "PAIRED" && !file.exists(fastq2)) stop("Missing FASTQ2: ", fastq2, call. = FALSE)
if (layout == "SINGLE" && fastq2 != "-") stop("Single-end library must use '-' for FASTQ2", call. = FALSE)

index_dir <- dirname(index_prefix)
index_stem <- basename(index_prefix)
index_files <- list.files(index_dir, pattern = paste0("^", index_stem, "\\."), full.names = TRUE)
core_index_files <- index_files[grepl("(\\.files$|\\.reads$|\\.[0-9]+\\.b\\.(array|tab)$)", index_files)]
if (length(core_index_files) < 4L || any(file.info(core_index_files)$size <= 0)) {
  stop("Missing or incomplete Rsubread index: ", index_prefix, call. = FALSE)
}

marker <- file.path(output_dir, "COMPLETE.tsv")
expected_outputs <- c(
  "mapping_summary.tsv",
  paste0("counts_strand", 0:2, ".tsv.gz"),
  paste0("assignment_strand", 0:2, ".tsv"),
  "counts_all_strands.rds",
  "session_info.txt"
)
if (file.exists(marker)) {
  completed <- utils::read.delim(marker, check.names = FALSE)
  output_manifest_path <- file.path(output_dir, "OUTPUT_MANIFEST.tsv")
  required_completion_columns <- c(
    "dataset_id", "sample_id", "layout", "rsubread_version", "assigned_strand0",
    "assigned_strand1", "assigned_strand2", "assigned_fraction_strand0", "status"
  )
  valid_completion <- nrow(completed) == 1L &&
    all(required_completion_columns %in% names(completed)) &&
    completed$dataset_id[[1L]] == dataset_id && completed$sample_id[[1L]] == sample_id &&
    completed$layout[[1L]] == layout && completed$rsubread_version[[1L]] == "2.26.0" &&
    completed$status[[1L]] == "PASS" && file.exists(output_manifest_path)
  if (valid_completion) {
    assigned_values <- suppressWarnings(as.numeric(unlist(completed[1L, c(
      "assigned_strand0", "assigned_strand1", "assigned_strand2", "assigned_fraction_strand0"
    )], use.names = FALSE)))
    valid_completion <- all(is.finite(assigned_values)) && assigned_values[[1L]] > 0 &&
      all(assigned_values[2:3] >= 0) && assigned_values[[4L]] > 0 && assigned_values[[4L]] <= 1
  }
  if (valid_completion) {
    output_manifest <- utils::read.delim(output_manifest_path, check.names = FALSE)
    valid_completion <- nrow(output_manifest) == length(expected_outputs) &&
      all(c("filename", "bytes", "md5") %in% names(output_manifest)) &&
      !anyDuplicated(output_manifest$filename) && setequal(output_manifest$filename, expected_outputs) &&
      all(basename(output_manifest$filename) == output_manifest$filename) &&
      all(grepl("^[0-9a-fA-F]{32}$", output_manifest$md5))
  }
  if (valid_completion) {
    paths <- file.path(output_dir, output_manifest$filename)
    valid_completion <- all(file.exists(paths)) &&
      all(file.info(paths)$size == as.numeric(output_manifest$bytes)) &&
      all(unname(tools::md5sum(paths)) == tolower(output_manifest$md5))
  }
  if (!isTRUE(valid_completion)) {
    stop("Existing completion marker/output manifest fails validation: ", output_dir, call. = FALSE)
  }
  message("PASS existing completed sample: ", dataset_id, "/", sample_id)
  quit(save = "no", status = 0L)
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop(
    "Partial unmarked sample output exists; preserve/archive before retrying: ",
    output_dir,
    call. = FALSE
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bam_work_dir, recursive = TRUE, showWarnings = FALSE)

bam <- file.path(bam_work_dir, paste0(sample_id, ".subjunc.bam"))
started <- Sys.time()
message("Aligning ", dataset_id, "/", sample_id, " (", layout, ")")

mapping <- Rsubread::subjunc(
  index = index_prefix,
  readfile1 = fastq1,
  readfile2 = if (layout == "PAIRED") fastq2 else NULL,
  input_format = "gzFASTQ",
  output_format = "BAM",
  output_file = bam,
  phredOffset = 33,
  unique = FALSE,
  nBestLocations = 1,
  nthreads = 4,
  sortReadsByCoordinates = FALSE,
  useAnnotation = TRUE,
  annot.ext = gtf,
  isGTF = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id"
)
if (!file.exists(bam) || file.info(bam)$size <= 0) stop("Subjunc did not create a BAM", call. = FALSE)

mapping_df <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping_df$metric <- rownames(mapping_df)
mapping_df <- mapping_df[, c("metric", setdiff(names(mapping_df), "metric")), drop = FALSE]
rownames(mapping_df) <- NULL
mapping_df$dataset_id <- dataset_id
mapping_df$sample_id <- sample_id
utils::write.table(
  mapping_df,
  file.path(output_dir, "mapping_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

run_feature_counts <- function(strand_mode) {
  Rsubread::featureCounts(
    files = bam,
    annot.ext = gtf,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    GTF.attrType.extra = c("gene", "gene_biotype"),
    useMetaFeatures = TRUE,
    allowMultiOverlap = FALSE,
    minOverlap = 1,
    countMultiMappingReads = FALSE,
    fraction = FALSE,
    minMQS = 10,
    primaryOnly = TRUE,
    strandSpecific = strand_mode,
    isPairedEnd = layout == "PAIRED",
    countReadPairs = layout == "PAIRED",
    requireBothEndsMapped = layout == "PAIRED",
    countChimericFragments = FALSE,
    autosort = TRUE,
    nthreads = 4,
    tmpDir = bam_work_dir,
    verbose = FALSE
  )
}

fc_objects <- vector("list", 3L)
names(fc_objects) <- paste0("strand", 0:2)
for (strand_mode in 0:2) {
  message("Counting ", sample_id, " with strandSpecific=", strand_mode)
  fc <- run_feature_counts(strand_mode)
  if (ncol(fc$counts) != 1L || nrow(fc$counts) == 0L) {
    stop("Unexpected featureCounts dimensions for strand mode ", strand_mode, call. = FALSE)
  }
  count_table <- cbind(
    fc$annotation,
    count = as.numeric(fc$counts[, 1L])
  )
  count_path <- file.path(output_dir, paste0("counts_strand", strand_mode, ".tsv.gz"))
  con <- gzfile(count_path, open = "wt", encoding = "UTF-8")
  utils::write.table(count_table, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)

  stat <- as.data.frame(fc$stat, stringsAsFactors = FALSE)
  names(stat)[1L] <- "Status"
  stat$dataset_id <- dataset_id
  stat$sample_id <- sample_id
  stat$strand_mode <- strand_mode
  utils::write.table(
    stat,
    file.path(output_dir, paste0("assignment_strand", strand_mode, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  fc_objects[[strand_mode + 1L]] <- list(
    counts = setNames(as.numeric(fc$counts[, 1L]), rownames(fc$counts)),
    annotation = fc$annotation,
    stat = stat
  )
}

saveRDS(
  list(
    dataset_id = dataset_id,
    sample_id = sample_id,
    layout = layout,
    reference_index = index_prefix,
    annotation_gtf = gtf,
    mapping = mapping_df,
    feature_counts = fc_objects
  ),
  file.path(output_dir, "counts_all_strands.rds"),
  compress = "xz"
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

assigned <- vapply(
  fc_objects,
  function(x) {
    idx <- which(x$stat$Status == "Assigned")
    if (length(idx) != 1L) return(NA_real_)
    as.numeric(x$stat[idx, 2L])
  },
  numeric(1L)
)
totals <- vapply(
  fc_objects,
  function(x) sum(as.numeric(x$stat[, 2L]), na.rm = TRUE),
  numeric(1L)
)
if (any(!is.finite(assigned)) || any(!is.finite(totals)) || any(totals <= 0) ||
    assigned[[1L]] <= 0 || any(assigned < 0) || any(assigned > totals)) {
  stop(
    "featureCounts returned an invalid or zero unstranded assignment; retaining FASTQ/BAM for adjudication",
    call. = FALSE
  )
}

manifest_files <- list.files(output_dir, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest_files <- setdiff(manifest_files, marker)
output_manifest <- data.frame(
  filename = basename(manifest_files),
  bytes = file.info(manifest_files)$size,
  md5 = unname(tools::md5sum(manifest_files)),
  stringsAsFactors = FALSE
)
utils::write.table(
  output_manifest,
  file.path(output_dir, "OUTPUT_MANIFEST.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

finished <- Sys.time()
record <- data.frame(
  dataset_id = dataset_id,
  sample_id = sample_id,
  layout = layout,
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  rsubread_version = as.character(utils::packageVersion("Rsubread")),
  assigned_strand0 = assigned[[1L]],
  assigned_strand1 = assigned[[2L]],
  assigned_strand2 = assigned[[3L]],
  assigned_fraction_strand0 = assigned[[1L]] / totals[[1L]],
  bam_bytes_before_cache_cleanup = file.info(bam)$size,
  started_at = format(started, tz = "UTC", usetz = TRUE),
  finished_at = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_minutes = as.numeric(difftime(finished, started, units = "mins")),
  status = "PASS",
  stringsAsFactors = FALSE
)
utils::write.table(record, marker, sep = "\t", quote = FALSE, row.names = FALSE)
message("PASS aligned and counted ", dataset_id, "/", sample_id)
