# Path-independent orchestration; the three scientific modules remain byte-identical.
# Called only by the Python entry after payload/input identity verification.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
if (nzchar(args[4])) {
  libraries <- strsplit(args[4], .Platform$path.sep, fixed = TRUE)[[1]]
  .libPaths(c(normalizePath(libraries, winslash = "/", mustWork = TRUE), .libPaths()))
}
stopifnot(requireNamespace("jsonlite", quietly = TRUE), requireNamespace("digest", quietly = TRUE))
read_json <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
sha <- function(p) digest::digest(file = p, algo = "sha256")
env <- read_json(file.path(root, "environment/R_dependency_contract.json"))
stopifnot(identical(as.character(getRversion()), env$R))
for (p in names(env$packages)) {
  observed <- as.character(utils::packageVersion(p))
  if (!identical(observed, as.character(package_version(env$packages[[p]])))) {
    stop("R dependency mismatch: ", p, "; expected ", env$packages[[p]], "; found ", observed)
  }
}
for (p in unlist(env$direct)) stopifnot(requireNamespace(p, quietly = TRUE))
if (args[2] == "--environment-only") {
  cat(jsonlite::toJSON(list(status = "PASS_R_DEPENDENCY_CLOSURE", R = as.character(getRversion()),
                           package_count = length(env$packages), models_fitted = 0L), auto_unbox = TRUE), "\n")
  quit(status = 0L)
}
job <- read_json(args[2])
out <- normalizePath(args[3], winslash = "/", mustWork = TRUE)
stopifnot(identical(job$mode, "run") || identical(job$mode, "resume"))
Sys.setenv(TZ = "UTC")
Sys.setlocale("LC_NUMERIC", "C")
set.seed(20261025L)
for (name in c("gene_contrast", "programs", "signed_activity")) {
  source(file.path(root, "analysis/original", paste0("analyze_PRJNA886048_", name, "_20261025.R")))
}
samples <- read.delim(file.path(root, "data/derived/samples.tsv"), stringsAsFactors = FALSE)
prjna886048_design(samples)
write_tsv <- function(x, name) {
  p <- file.path(out, "data/output", name)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
counts <- c(import_executed = 0L, gene_fit_executed = 0L, camera_executed = 0L, activity_executed = 0L,
            adopted_checkpoints = 0L, reused_checkpoints = 0L)
checkpoint <- function(key, counter, expr) {
  p <- file.path(out, "checkpoints", paste0(key, ".rds"))
  receipt <- paste0(p, ".json")
  pending <- paste0(p, ".pending")
  if (file.exists(p) && file.exists(receipt)) {
    z <- read_json(receipt)
    stopifnot(identical(z$binding, job$binding), identical(z$sha256, sha(p)))
    # A lingering pending marker after an atomic verified commit is harmless.
    counts["reused_checkpoints"] <<- counts["reused_checkpoints"] + 1L
    return(readRDS(p))
  }
  if (file.exists(p) || file.exists(receipt) || file.exists(pending) || file.exists(paste0(p, ".partial"))) {
    stop("Incomplete checkpoint requires inspection; no automatic stage repetition: ", key)
  }
  writeLines(job$binding, pending)
  if (!is.null(job$adopt[[key]])) {
    a <- job$adopt[[key]]
    stopifnot(identical(sha(a$path), a$sha256))
    value <- readRDS(a$path)
    stopifnot(file.copy(a$path, paste0(p, ".partial")))
    counts["adopted_checkpoints"] <<- counts["adopted_checkpoints"] + 1L
  } else {
    value <- force(expr)
    counts[counter] <<- counts[counter] + 1L
    saveRDS(value, paste0(p, ".partial"))
  }
  stopifnot(file.rename(paste0(p, ".partial"), p))
  jsonlite::write_json(list(binding = job$binding, sha256 = sha(p), key = key),
                      paste0(receipt, ".partial"), auto_unbox = TRUE, pretty = TRUE)
  stopifnot(file.rename(paste0(receipt, ".partial"), receipt))
  value
}
quant <- file.path(out, "data/intermediate/quantification_inputs")
for (run in samples$run) {
  src <- file.path(job$inputs$native_root, paste0("linux_full_", run, "_01"), "quant/abundance.tsv")
  dst <- file.path(quant, run, "abundance.tsv")
  dir.create(dirname(dst), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(dst)) stopifnot(file.copy(src, dst))
  stopifnot(identical(sha(dst), sha(src)))
}
imported <- checkpoint("import", "import_executed", prjna886048_import(samples, quant, job$inputs$tx2gene))
fit <- checkpoint("genes", "gene_fit_executed", prjna886048_fit_genes(imported$length_adjusted$counts, samples))
stopifnot(identical(colnames(fit$voom$E), samples$run))
edges <- read.delim(job$inputs$mapped_edges, check.names = FALSE, stringsAsFactors = FALSE)
terms <- read.delim(job$inputs$reference_terms, check.names = FALSE, stringsAsFactors = FALSE)
programs <- checkpoint("programs", "camera_executed", prjna886048_camera(fit, edges, reference_terms = terms))
activity <- checkpoint("activity", "activity_executed", prjna886048_signed_activity(fit, edges, samples, terms))
for (key in c("primary_edgeR", "sensitivity_voom", "gene_filter")) write_tsv(fit[[key]], paste0("genes/", key, ".tsv"))
write_tsv(programs$camera, "camera_complete_results.tsv")
write_tsv(programs$coverage, "camera_complete_coverage.tsv")
write_tsv(activity$contrasts, "signed_activity_contrasts.tsv")
write_tsv(activity$coverage, "signed_activity_coverage.tsv")
for (collection in names(activity$scores)) {
  mat <- activity$scores[[collection]]
  write_tsv(data.frame(term_id = rownames(mat), mat, check.names = FALSE, row.names = NULL),
            paste0(collection, "_per_library_scores.tsv"))
}
write_tsv(samples, "samples.tsv")
# Exactly the presentation-only calculation used in the reviewed saved-MDS renderer.
mds <- limma::plotMDS(fit$dge, plot = FALSE)
write_tsv(data.frame(run = samples$run, treatment = samples$treatment, x = mds$x, y = mds$y),
          "complete_cohort_MDS_coordinates.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "provenance/sessionInfo.txt"))
jsonlite::write_json(list(status = "PORTABLE_R_STAGES_COMPLETE", binding = job$binding,
                         counts = as.list(counts), random_seed = 20261025L,
                         deterministic_methods = "No stochastic sampling in the frozen analysis modules",
                         R = as.character(getRversion()), scope = "Acute SCI context; not protected PSS validation"),
                     file.path(out, "provenance/R_execution.json"), auto_unbox = TRUE, pretty = TRUE)
cat(jsonlite::toJSON(as.list(counts), auto_unbox = TRUE), "\n")
