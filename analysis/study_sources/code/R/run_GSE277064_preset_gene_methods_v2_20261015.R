#!/usr/bin/env Rscript
# Real-data entry prepared for the separately requested limited addition.
# Do not invoke before the owner's explicit new-fit decision and source review.
# No primary edgeR fit, preprocessing, program test or candidate gate is called.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
library(data.table)
library(jsonlite)
library(digest)

assert <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
identity_file <- function(path) {
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
       bytes = unname(file.info(path)$size),
       sha256 = digest(file = path, algo = "sha256", serialize = FALSE))
}
verify_file <- function(binding) {
  actual <- identity_file(binding$path)
  assert(actual$bytes == binding$bytes && actual$sha256 == binding$sha256,
         paste("Pinned source changed:", basename(binding$path)))
  invisible(actual)
}
write_json_new <- function(value, path) {
  assert(!file.exists(path), "Never overwrite a result/receipt")
  write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = NA,
             na = "null", null = "null")
}
save_checkpoint <- function(value, path) {
  assert(!file.exists(path) && !file.exists(paste0(path, ".partial")),
         "Preserve earlier checkpoints")
  temporary <- paste0(path, ".partial")
  saveRDS(value, temporary, compress = "gzip")
  assert(file.rename(temporary, path), "Checkpoint commit failed")
  identity_file(path)
}
capture_method <- function(expression) {
  conditions <- new.env(parent = emptyenv())
  conditions$warnings <- conditions$messages <- character()
  value <- withCallingHandlers(force(expression),
    warning = function(condition) {
      conditions$warnings <- c(conditions$warnings, conditionMessage(condition))
    },
    message = function(condition) {
      conditions$messages <- c(conditions$messages, conditionMessage(condition))
    })
  list(value = value, warnings = unique(conditions$warnings), messages = unique(conditions$messages))
}
assert_complete_rows <- function(result, gene_ids, contrasts) {
  assert(nrow(result) == length(gene_ids) * ncol(contrasts), "All six complete gene families required")
  assert(!anyDuplicated(result[, c("gene_id", "method", "contrast")]), "Duplicate output keys")
  for (contrast in colnames(contrasts)) {
    selected <- result[result$contrast == contrast, , drop = FALSE]
    assert(identical(selected$gene_id, gene_ids), "All source genes in original order")
    assert(all(is.na(selected$p_value) | (is.finite(selected$p_value) & selected$p_value >= 0 & selected$p_value <= 1)) &&
             all(is.na(selected$bh_q) | (is.finite(selected$bh_q) & selected$bh_q >= 0 & selected$bh_q <= 1)),
           "Native probabilities must remain legal or missing")
    assert(all(!selected$inference_eligible |
                 (is.finite(selected$log2_effect) & is.finite(selected$standard_error) &
                    selected$standard_error > 0 & is.finite(selected$p_value) & is.finite(selected$bh_q))),
           "No inference from non-testable/non-converged rows")
  }
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
assert(length(args) == 2L, "Explicit PROJECT_ROOT FRESH_R_DRIVE_OUTPUT required")
root <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
output <- gsub("\\\\", "/", args[2L])
assert(grepl("^R:/Shaoyao_storage/(data|cache)/[A-Za-z0-9_-]+$", output),
       "Added scientific outputs belong to a specifically named R-drive directory")
assert(!file.exists(output), "Use a fresh independent output, not an old run")
setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20260817L)
assert(as.character(getRversion()) == "4.6.1", "Use the reviewed R version")
versions <- c(edgeR = "4.10.1", limma = "3.68.4", DESeq2 = "1.52.0",
              data.table = "1.18.4", jsonlite = "2.0.0", digest = "0.6.39")
for (package in names(versions)) {
  assert(as.character(packageVersion(package)) == versions[package], "Installed package version differs")
}
preparation_path <- file.path(root, "results/diagnostics/GSE277064_preset_gene_method_preparation_20261013/PREPARATION_RECEIPT_20261014.json")
assert(identity_file(preparation_path)$sha256 == "19958d0cd1ef0dbf421a0ac90fb3456be05f912d8a5bbd66cb13bda850335daa",
       "Use the exact reviewed preparation receipt")
preparation <- fromJSON(preparation_path, simplifyVector = FALSE)
core_path <- file.path(root, "code/R/GSE277064_gene_method_core_v2_20261015.R")
assert(identity_file(core_path)$sha256 == "e90ec660f5fa0b4f93b3d36879c9e94e5b4afec167cf7b0f04c8566b236e713b",
       "Use the reviewed key-attribute-only corrected method core")
verify_file(preparation$core)
verify_file(preparation$test)
verify_file(preparation$readiness)
verify_file(preparation$independent_review)
verify_file(preparation$proposal)
for (binding in preparation$synthetic_files) verify_file(binding)
assert(preparation$synthetic_check_n == 84L && preparation$observed_final_synthetic_exit_code == 0L &&
         !isTRUE(preparation$real_fits_run), "Bound preparation and synthetic scope")
readiness <- fromJSON(preparation$readiness$path, simplifyVector = FALSE)
assert(length(readiness$inputs) == 26L && readiness$natural_unique_library_n == 45L,
       "Expected frozen natural-course input inventory")
for (binding in readiness$inputs) verify_file(binding)
source(core_path, local = TRUE)
core_before <- identity_file(core_path)
entry_path <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
entry_before <- identity_file(entry_path)

assert(dir.create(output, recursive = FALSE), "Cannot create explicit fresh output")
write_json_new(list(status = "ADDITIONAL_GENE_METHOD_RUN_STARTED_NOT_ACCEPTED",
                    started_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                    process_id = Sys.getpid(), preparation = identity_file(preparation_path),
                    driver = entry_before, core = core_before, seed = 20260817L,
                    actual_R = R.version.string, actual_packages = as.list(versions),
                    actual_library_paths = .libPaths(), actual_RNG_kind = RNGkind(),
                    source_raw_readonly = TRUE, original_primary_or_gates_changed = FALSE),
               file.path(output, "STARTED.json"))

tryCatch({
  products <- list()
  summaries <- list()
  for (cell_type in c("astrocyte", "microglia", "v2a")) {
    object_path <- file.path(root, "results/GSE277064/objects", paste0(cell_type, "_analysis.rds"))
    object_key <- paste0("results/GSE277064/objects/", cell_type, "_analysis.rds")
    verify_file(readiness$inputs[[object_key]])
    original <- readRDS(object_path)
    assert(identical(original$cell_type, cell_type) && identical(original$domain, "natural_course"),
           "Only saved natural-course cell-specific models are eligible")
    counts <- original$dge$counts
    groups <- as.character(original$sample_metadata$group)
    assert(identical(colnames(counts), as.character(original$sample_metadata$source_column_or_member)),
           "Saved library/metadata order must remain identical")
    expected_n <- readiness$cell_types[[cell_type]]$filtered_gene_rows
    assert(nrow(counts) == expected_n && ncol(counts) == 15L, "Fixed cell-specific gene universe")
    validate_gene_method_inputs(counts, groups, original$design, original$contrasts, original$voom)
    count_hash <- digest(counts, algo = "sha256")
    for (method in c("saved_voom", "integer_DESeq2")) {
      started <- proc.time()[["elapsed"]]
      captured <- if (method == "saved_voom") {
        capture_method(fit_saved_voom_genes(counts, groups, original$design, original$contrasts, original$voom))
      } else {
        capture_method(fit_integer_deseq_genes(counts, groups, original$design, original$contrasts))
      }
      elapsed <- proc.time()[["elapsed"]] - started
      assert_complete_rows(captured$value$results, rownames(counts), original$contrasts)
      assert(identical(digest(counts, algo = "sha256"), count_hash), "Unnormalized source counts changed in RAM")
      captured$value$results$cell_type <- cell_type
      captured$value$summary$cell_type <- cell_type
      checkpoint <- list(cell_type = cell_type, method = method,
                         source = readiness$inputs[[object_key]], driver = entry_before, core = core_before,
                         scientific_result = captured$value, warning_messages = captured$warnings,
                         messages = captured$messages, elapsed_seconds = elapsed,
                         original_counts_unchanged = TRUE, candidate_gate_eligible = FALSE)
      key <- paste(cell_type, method, sep = "_")
      products[[paste0(key, "_checkpoint")]] <- save_checkpoint(checkpoint, file.path(output, paste0(key, ".rds")))
      table_path <- file.path(output, paste0(key, "_genes.tsv.gz"))
      assert(!file.exists(table_path), "Preserve earlier result table")
      fwrite(captured$value$results, table_path, sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
      products[[paste0(key, "_table")]] <- identity_file(table_path)
      summaries[[key]] <- captured$value$summary
      cat("COMPLETED", cell_type, method, "six families; seconds", elapsed, "\n")
      rm(captured, checkpoint)
      gc(FALSE)
    }
    verify_file(readiness$inputs[[object_key]])
    rm(original, counts)
    gc(FALSE)
  }
  summary_path <- file.path(output, "GENE_FAMILY_SUMMARY.tsv")
  summary <- rbindlist(summaries)
  assert(nrow(summary) == 36L, "Three cells by two methods by six contrasts")
  fwrite(summary, summary_path, sep = "\t", quote = TRUE, na = "NA")
  products$family_summary <- identity_file(summary_path)
  for (binding in readiness$inputs) verify_file(binding)
  verify_file(core_before)
  verify_file(entry_before)
  session_path <- file.path(output, "SESSION_INFO.txt")
  writeLines(capture.output(sessionInfo()), session_path)
  products$session_info <- identity_file(session_path)
  write_json_new(list(status = "ADDITIONAL_GENE_METHOD_FITS_COMPLETE_INDEPENDENT_CHECK_PENDING",
                      completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                      libraries = 45L, cells = 3L, fitted_methods_per_cell = 2L,
                      contrasts_per_method_cell = 6L, result_families = 36L,
                      preparation = identity_file(preparation_path), core = core_before, driver = entry_before,
                      products = products, all_source_pins_unchanged = TRUE,
                      original_edgeR_refitted = FALSE, preprocessing_recomputed = FALSE,
                      gene_universes_changed = FALSE, candidate_gate_eligible = FALSE,
                      original_primary_or_gates_changed = FALSE,
                      historical_runtime_authenticated = FALSE,
                      independent_real_result_check_passed = FALSE,
                      whole_goal_complete = FALSE, overall_scientific_status = "NEEDS REVIEW"),
                 file.path(output, "COMPLETE.json"))
  cat("ADDITIONAL_GENE_METHOD_FITS_COMPLETE; independent results check pending\n")
}, error = function(condition) {
  write_json_new(list(status = "ADDITIONAL_GENE_METHOD_RUN_FAILED_PARTIAL_OUTPUTS_PRESERVED",
                      condition_class = class(condition), error = conditionMessage(condition),
                      original_primary_or_gates_changed = FALSE, independent_real_result_check_passed = FALSE,
                      whole_goal_complete = FALSE), file.path(output, "FAILED.json"))
  stop(condition)
})
