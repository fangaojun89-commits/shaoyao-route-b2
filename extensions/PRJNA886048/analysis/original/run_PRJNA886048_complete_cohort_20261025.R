# Assemble the fixed complete cohort and call the already checked analysis modules.
# --check-only reports readiness without staging files or fitting biological models.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(3L, 4L))
project <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(project, "environment/R-library/4.6"), "C:/R/Library", .libPaths()))
data_root <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
output <- args[3]
check_only <- length(args) == 4L && identical(args[4], "--check-only")
if (length(args) == 4L) stopifnot(check_only)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256")
sample_path <- file.path(project, "manuscript/supplement/Route_B2_baclofen_SCI_sample_map_20261025.tsv")
samples <- read.delim(sample_path, stringsAsFactors = FALSE, check.names = FALSE)
ena_path <- file.path(dirname(dirname(data_root)),
                      "literature/SCI_baclofen_PRJNA886048_20261025/ENA_quantification_manifest.tsv")
ena <- read.delim(ena_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(!anyDuplicated(ena$run_accession), all(samples$run %in% ena$run_accession))
source(file.path(project, "code/R/analyze_PRJNA886048_gene_contrast_20261025.R"))
invisible(prjna886048_design(samples))
stopifnot(setequal(samples$run, paste0("SRR217758", 46:57)))
native <- setNames(file.path(data_root, paste0("linux_full_", samples$run, "_01")), samples$run)
ready <- file.exists(file.path(native, "receipt.json"))
inventory <- data.frame(run = samples$run, treatment = samples$treatment,
                        status = ifelse(ready, "RECEIPT_PRESENT", "PENDING_QUANTIFICATION"))
qc <- list()
receipts <- list()
expected_index <- "8345052ff38d43971d50cbff8f5b767ad438cfb6e868c820d836bc056d2b1f73"
commands_complete <- function(commands) {
  succeeded <- vapply(commands, function(z) identical(z$exit_code, 0L), logical(1))
  if (all(succeeded)) return(TRUE)
  # Preserve the interrupted kallisto record when its identical command later succeeded.
  interrupted <- which(!succeeded)
  if (length(interrupted) != 1L || interrupted == length(commands)) return(FALSE)
  prior <- commands[[interrupted]]
  final <- commands[[length(commands)]]
  identical(prior$stage, "KALLISTO") && is.null(prior$exit_code) &&
    identical(final$stage, "KALLISTO_RESUME_R2") && identical(final$exit_code, 0L) &&
    identical(prior$command, final$command)
}
for (i in which(ready)) {
  run <- samples$run[i]
  d <- native[[run]]
  r <- read_json(file.path(d, "receipt.json"))
  if (!identical(r$status, "COMPLETE_FULL_LIBRARY_QUANTIFICATION")) {
    ready[i] <- FALSE
    inventory$status[i] <- r$status
    next
  }
  stopifnot(identical(r$library, run), identical(r$index_sha256, expected_index),
            identical(r$orientation, "unstranded"),
            commands_complete(r$commands))
  submitted <- ena[match(run, ena$run_accession), , drop = FALSE]
  for (mate in 1:2) {
    input <- r$input_checks[[paste0(run, "_", mate, ".fastq.gz")]]
    stopifnot(identical(input$md5, strsplit(submitted$fastq_md5, ";", fixed = TRUE)[[1]][mate]),
              input$bytes == as.numeric(strsplit(submitted$fastq_bytes, ";", fixed = TRUE)[[1]][mate]))
  }
  for (name in names(r$outputs)) {
    p <- file.path(d, name)
    stopifnot(file.exists(p), file.info(p)$size == r$outputs[[name]]$bytes,
              identical(sha(p), r$outputs[[name]]$sha256))
  }
  f <- read_json(file.path(d, "fastp.json"))$summary
  k <- read_json(file.path(d, "quant/run_info.json"))
  stopifnot(k$n_targets == 146190L, f$after_filtering$total_reads == 2 * k$n_processed,
            k$n_pseudoaligned <= k$n_processed, k$n_processed > 0)
  inventory$status[i] <- "COMPLETE_OUTPUT_HASHES_AND_READ_COUNTS_VERIFIED"
  qc[[run]] <- data.frame(run = run, input_pairs = f$before_filtering$total_reads / 2,
                          retained_pairs = k$n_processed, pseudoaligned_pairs = k$n_pseudoaligned,
                          pseudoaligned_fraction = k$n_pseudoaligned / k$n_processed,
                          post_q30_fraction = f$after_filtering$q30_rate,
                          post_gc_fraction = f$after_filtering$gc_content)
  receipts[[run]] <- r
}
if (check_only || !all(ready)) {
  cat(jsonlite::toJSON(list(status = if (all(ready)) "COMPLETE_COHORT_READY_FOR_QC_AND_FIT" else
                             "WAITING_FOR_REMAINING_LIBRARIES", completed = sum(ready), total = 12L,
                           inventory = inventory, biological_models_run = 0L),
                       auto_unbox = TRUE, pretty = TRUE), "\n")
  quit(status = if (check_only) 0L else 2L)
}
map_path <- file.path(data_root, "reference_preparation/kallisto_tx2gene.tsv")
edge_path <- file.path(data_root, "program_mapping_20261025/mapped_edges.tsv")
term_path <- file.path(data_root, "program_mapping_20261025/reference_coverage.tsv")
module_paths <- file.path(project, "code/R", c("analyze_PRJNA886048_gene_contrast_20261025.R",
                           "analyze_PRJNA886048_programs_20261025.R",
                           "analyze_PRJNA886048_signed_activity_20261025.R"))
bound_paths <- c(sample_path, ena_path, map_path, edge_path, term_path, module_paths,
                 file.path(project, "code/R/run_PRJNA886048_complete_cohort_20261025.R"),
                 file.path(native, "quant/abundance.tsv"))
identity <- setNames(vapply(bound_paths, sha, character(1)), bound_paths)
identity_path <- file.path(output, "input_and_module_identity.rds")
if (dir.exists(output)) {
  stopifnot(file.exists(identity_path), identical(readRDS(identity_path), identity))
  if (file.exists(file.path(output, "receipt.json"))) {
    cat("REUSED_COMPLETED_ANALYSIS; no models rerun\n")
    quit(status = 0L)
  }
} else {
  dir.create(output, recursive = TRUE)
  saveRDS(identity, identity_path)
}
quant_root <- file.path(output, "quantification_inputs")
for (run in samples$run) {
  to <- file.path(quant_root, run)
  if (!dir.exists(to)) dir.create(to, recursive = TRUE)
  destination <- file.path(to, "abundance.tsv")
  if (!file.exists(destination)) {
    stopifnot(file.copy(file.path(native[[run]], "quant/abundance.tsv"), destination))
  }
  stopifnot(identical(sha(destination), sha(file.path(native[[run]], "quant/abundance.tsv"))))
}
write_tsv <- function(x, name) write.table(x, file.path(output, name), sep = "\t",
                                          quote = FALSE, row.names = FALSE, na = "NA")
write_tsv(samples, "samples.tsv")
write_tsv(do.call(rbind, qc), "complete_library_qc.tsv")
import_path <- file.path(output, "gene_import_objects.rds")
if (file.exists(import_path)) imported <- readRDS(import_path) else {
  imported <- prjna886048_import(samples, quant_root, map_path)
  saveRDS(imported, import_path)
}
write_tsv(imported$inputs, "quantification_input_hashes.tsv")
fit_path <- file.path(output, "genes/gene_model_objects.rds")
if (file.exists(fit_path)) fit <- readRDS(fit_path) else {
  fit <- prjna886048_fit_genes(imported$length_adjusted$counts, samples)
}
prjna886048_save_genes(fit, file.path(output, "genes"))
# All libraries remain in the analysis. MDS is descriptive QC, not an exclusion rule.
grDevices::pdf(file.path(output, "complete_cohort_MDS.pdf"), width = 8, height = 6)
edgeR::plotMDS(fit$dge, labels = samples$run,
              col = ifelse(samples$treatment == "Baclofen", "#1F77B4", "#E1812C"))
grDevices::dev.off()
source(file.path(project, "code/R/analyze_PRJNA886048_programs_20261025.R"))
source(file.path(project, "code/R/analyze_PRJNA886048_signed_activity_20261025.R"))
edges <- read.delim(edge_path, stringsAsFactors = FALSE, check.names = FALSE)
terms <- read.delim(term_path, stringsAsFactors = FALSE, check.names = FALSE)
program_path <- file.path(output, "program_membership_objects.rds")
activity_path <- file.path(output, "signed_activity_objects.rds")
if (file.exists(program_path)) programs <- readRDS(program_path) else {
  programs <- prjna886048_camera(fit, edges, reference_terms = terms)
  saveRDS(programs, program_path)
}
if (file.exists(activity_path)) activity <- readRDS(activity_path) else {
  activity <- prjna886048_signed_activity(fit, edges, samples, terms)
  saveRDS(activity, activity_path)
}
write_tsv(programs$camera, "camera_complete_results.tsv")
write_tsv(programs$coverage, "camera_complete_coverage.tsv")
write_tsv(activity$contrasts, "signed_activity_contrasts.tsv")
write_tsv(activity$coverage, "signed_activity_coverage.tsv")
for (collection in names(activity$scores)) {
  mat <- activity$scores[[collection]]
  write_tsv(data.frame(term_id = rownames(mat), mat, check.names = FALSE, row.names = NULL),
            paste0(collection, "_per_library_scores.tsv"))
}
jsonlite::write_json(list(status = "COMPLETE_FIXED_COHORT_ANALYSIS", utc = format(Sys.time(), tz = "UTC"),
                         libraries = samples$run, group_sizes = c(vehicle = 6L, baclofen = 6L),
                         sample_map_sha256 = sha(sample_path), tx2gene_sha256 = sha(map_path),
                         mapped_edges_sha256 = sha(edge_path), reference_terms_sha256 = sha(term_path),
                         tested_genes = nrow(fit$primary_edgeR),
                         scientific_scope = "Retrospective pharmacological context; S1 unchanged; not protected PSS validation",
                         previously_completed_source_models_rerun = FALSE),
                     file.path(output, "receipt.json"), auto_unbox = TRUE, pretty = TRUE)
