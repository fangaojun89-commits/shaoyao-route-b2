# Executes only original metadata and deposited-count readers, never models.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
contract <- spec$contract
stopifnot(as.character(getRversion()) == contract$R)
versions <- vapply(names(contract$packages), function(p) as.character(utils::packageVersion(p)), "")
stopifnot(identical(unname(versions), unname(unlist(contract$packages))))
primary <- file.path(spec$project, contract$sequence[[1L]])
expressions <- parse(primary, keep.source = FALSE)
e <- new.env(parent = globalenv())
e$project_root <- spec$project
e$config <- yaml::read_yaml(file.path(spec$project, "config/contrasts.yaml"))
stopifnot(as.integer(e$config$global$random_seed) == contract$seed)
e$ledger <- utils::read.delim(file.path(spec$project, "metadata/samples.tsv"), check.names = FALSE, na.strings = c("", "NA"))
e$geo_long <- utils::read.delim(file.path(spec$project, "metadata/geo_series_metadata_long.tsv"), check.names = FALSE, na.strings = c("", "NA"))
e$datasets_ledger <- utils::read.delim(file.path(spec$project, "metadata/datasets.tsv"), check.names = FALSE, na.strings = c("", "NA"))
names_to_bind <- c("stopf", "assert_true", "safe_numeric", "get_nested", "normalize_sample_table",
                   "validate_geo_metadata", "extract_rsem_counts", "read_c2c12_counts", "dataset_specs")
bound <- character()
for (x in expressions) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]])) {
    name <- as.character(x[[2L]])
    if (name %in% names_to_bind) {
      stopifnot(!name %in% bound)
      if (name != "dataset_specs") stopifnot(is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")))
      eval(x, e); bound <- c(bound, name)
    }
  }
}
stopifnot(setequal(bound, names_to_bind), !exists("analyze_one", envir = e, inherits = FALSE))
s1 <- utils::read.delim(file.path(spec$project, "results/GSE16710/program_level/S1_GATE.tsv"), check.names = FALSE)
stopifnot(nrow(s1) == 1771L, !anyNA(s1$S1_pass), all(toupper(as.character(s1$S1_pass)) == "FALSE"))
rows <- list()
for (gse in names(e$dataset_specs)) {
  ds <- e$dataset_specs[[gse]]
  samples <- e$normalize_sample_table(gse, ds)
  metadata <- e$validate_geo_metadata(gse, samples, ds)
  extract_dir <- file.path(spec$output, "input_check", gse)
  original <- if (ds$raw_type == "rsem_tar") e$extract_rsem_counts(gse, samples, ds, extract_dir) else e$read_c2c12_counts(gse, samples, ds)
  reference <- file.path(spec$repo, "data/jgw03_inputs", gse)
  count_table <- utils::read.delim(file.path(reference, "locked_expected_counts.tsv.gz"), check.names = FALSE)
  counts <- as.matrix(count_table[, -1L, drop = FALSE]); rownames(counts) <- count_table$raw_gene_id
  stopifnot(identical(dimnames(original$counts), dimnames(counts)),
            identical(dim(original$counts), dim(counts)), all(original$counts == counts))
  saved <- utils::read.delim(file.path(reference, "locked_samples.tsv"), check.names = FALSE)
  stopifnot(identical(names(samples), names(saved)), nrow(saved) == 12L)
  for (field in names(saved)) {
    left <- samples[[field]]; right <- saved[[field]]
    stopifnot(identical(is.na(left), is.na(right)))
    ok <- !is.na(left)
    if (field == "log10_dose_plus_1") stopifnot(all(abs(left[ok]-right[ok]) <= 1e-12*(1+abs(left[ok])))) else
      stopifnot(identical(as.character(left[ok]), as.character(right[ok])))
  }
  rows[[gse]] <- list(libraries = ncol(counts), genes = nrow(counts),
    compared_count_values = length(counts), exact_count_values_and_order = TRUE,
    sample_fields_match = TRUE, metadata_checks = nrow(metadata),
    metadata_limitation_checks = sum(metadata$status == "PASS_WITH_LIMITATION"))
}
for (source in contract$sequence) parse(file.path(spec$project, source), keep.source = FALSE)
capture.output(sessionInfo(), file = file.path(spec$output, "receipts/check_sessionInfo.txt"))
jsonlite::write_json(list(status = "PASS_ORIGINAL_SOURCE_READERS_NO_FITS", datasets = rows,
  R = as.character(getRversion()), packages = as.list(versions), locale = Sys.getlocale(),
  seed_in_original_config = contract$seed, original_functions_bound = bound,
  S1_programs = 1771L, S1_admissions = 0L, scientific_fits = 0L,
  score_or_normalization_runs = 0L, new_label_allocations = 0L,
  primary_writer_executed = FALSE, scope = "Source-deposited processed counts to accepted input tables; not FASTQ quantification or model validation"),
  file.path(spec$output, "receipts/INPUT_CHECK.json"), auto_unbox = TRUE, pretty = TRUE)
