# In-memory original-adapter checks. No input modification, scores or fits.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
prepared <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
api <- new.env(parent = asNamespace("stats"))
sys.source(file.path(repo, "analysis/reference_sources/code/R/JGW03_input_admission_v2_20261015.R"), envir = api)
input <- readRDS(file.path(prepared, "prepared/GSE289929_unfitted_inputs.rds"))
rejects <- function(expression, expected) {
  error <- tryCatch({ force(expression); NULL }, error = function(e) conditionMessage(e))
  stopifnot(is.character(error), length(error) == 1L, grepl(expected, error, fixed = TRUE))
  invisible(TRUE)
}
bad_samples <- input$samples
bad_samples$formula_ratio[bad_samples$condition == "treatment"] <- "3:1"
rejects(api$jgw03_v2_validate_metadata(bad_samples, input$accession), "Frozen treatment metadata differs")
bad_samples <- input$samples
bad_samples$sample_id[2L] <- bad_samples$sample_id[1L]
rejects(api$jgw03_v2_validate_metadata(bad_samples, input$accession), "Unique ordered sample IDs required")
bad_resources <- input$resources
bad_resources$progeny$tax_id[1L] <- 10090L
rejects(api$jgw03_v2_validate_resources(bad_resources, input$accession), "Resource species/release differs")
cat("PASS: preparation identity rejects a substituted formula ratio, duplicate sample ID and resource taxon; no scores or fits\n")
