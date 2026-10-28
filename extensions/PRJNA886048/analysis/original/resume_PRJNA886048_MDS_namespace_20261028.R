# Resume the already fitted cohort with one recorded plotting-namespace correction.
# The frozen driver and input identity remain unchanged; the executed derivative
# is retained separately, together with unchanged pre/post model hashes.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
project <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(project, "environment/R-library/4.6"), "C:/R/Library", .libPaths()))
output <- normalizePath(args[3], winslash = "/", mustWork = TRUE)
driver <- file.path(project, "code/R/run_PRJNA886048_complete_cohort_20261025.R")
sha <- function(p) digest::digest(file = p, algo = "sha256")
stopifnot(identical(sha(driver), "f6ea807cdb7d83109e6637b7f73e8398e04df64525e8eb1b62d76161487dca30"))
checkpoints <- file.path(output, c("input_and_module_identity.rds", "gene_import_objects.rds",
                                  "genes/gene_model_objects.rds"))
stopifnot(all(file.exists(checkpoints)), !file.exists(file.path(output, "receipt.json")))
before <- setNames(vapply(checkpoints, sha, character(1)), checkpoints)
stopifnot("plotMDS" %in% getNamespaceExports("limma"),
          is.function(getS3method("plotMDS", "DGEList", envir = asNamespace("edgeR"))))
lines <- readLines(driver, warn = FALSE)
at <- grep("edgeR::plotMDS(", lines, fixed = TRUE)
stopifnot(length(at) == 1L)
lines[at] <- sub("edgeR::plotMDS(", "limma::plotMDS(", lines[at], fixed = TRUE)
executed <- file.path(output, "executed_driver_MDS_namespace_20261028.R")
stopifnot(!file.exists(executed))
writeLines(lines, executed, useBytes = TRUE)
record <- list(status = "STARTED_CHECKPOINT_RESUME_PLOT_NAMESPACE_ONLY",
               original_driver = driver, original_driver_sha256 = sha(driver),
               executed_driver = executed, executed_driver_sha256 = sha(executed),
               change = "edgeR::plotMDS -> limma::plotMDS; DGEList dispatch retained",
               edgeR_version = as.character(packageVersion("edgeR")),
               limma_version = as.character(packageVersion("limma")),
               checkpoint_hashes_before = as.list(before), gene_models_refitted = FALSE,
               statistical_parameters_changed = FALSE)
record_path <- file.path(output, "MDS_namespace_resume_receipt_20261028.json")
jsonlite::write_json(record, record_path, auto_unbox = TRUE, pretty = TRUE)
sys.source(executed, envir = new.env(parent = globalenv()))
after <- setNames(vapply(checkpoints, sha, character(1)), checkpoints)
stopifnot(identical(before, after), file.exists(file.path(output, "receipt.json")))
record$status <- "PASS_COHORT_RESUMED_WITH_UNCHANGED_GENE_CHECKPOINTS"
record$checkpoint_hashes_after <- as.list(after)
record$complete_receipt_sha256 <- sha(file.path(output, "receipt.json"))
record$completed_utc <- format(Sys.time(), tz = "UTC")
jsonlite::write_json(record, record_path, auto_unbox = TRUE, pretty = TRUE)
cat(record$status, "\n")
