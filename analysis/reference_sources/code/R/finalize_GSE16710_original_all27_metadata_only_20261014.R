# Verify already computed tables and finish metadata. No matrix, activity or model is rebuilt.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1L]], mustWork = TRUE, winslash = "/")
out <- normalizePath(args[[2L]], mustWork = TRUE, winslash = "/")
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
sha <- function(p) tolower(unname(tools::sha256sum(p)))
stopifnot(!file.exists(file.path(out, "COMPLETE.json")))
registry_path <- file.path(out, "REGISTRY.json")
registry <- fromJSON(registry_path)
paths <- file.path(root, registry$inputs$relative_path)
observed_hashes <- vapply(paths, sha, character(1))
stopifnot(identical(registry$inputs$sha256, unname(observed_hashes)),
          identical(as.numeric(registry$inputs$bytes), as.numeric(file.info(paths)$size)),
          registry$helper_sha256 == sha(registry$helper_path))
compatibility <- fread(file.path(out, "ORIGINAL26_COMPATIBILITY.tsv"))
stopifnot(nrow(compatibility) == 41L, all(compatibility$pass),
          all(is.finite(compatibility$maximum_absolute_difference[!is.na(compatibility$maximum_absolute_difference)])))
files <- c("reactome_camera_original_all27_sensitivity.tsv", "reactome_activity_original_all27_sensitivity.tsv",
           "progeny_original_all27_sensitivity.tsv", "dorothea_original_all27_sensitivity.tsv")
layers <- c("CAMERA", "Reactome_mean_z", "PROGENy_ULM", "DoRothEA_ULM")
expected_contrasts <- c("spasticity_progression_interaction", "injury_at_day21", "injury_at_day60", "progression_within_spi", "time_within_sham")
family <- fread(file.path(out, "FAMILY_SUMMARY.tsv"))
validation <- list()
for (i in seq_along(files)) {
  x <- fread(file.path(out, files[[i]]))
  stopifnot(!anyDuplicated(x[, .(feature_id, contrast)]),
            setequal(unique(x$contrast), if (i == 1L) expected_contrasts[[1L]] else expected_contrasts),
            all(is.finite(x$p_value)), all(is.finite(x$q_value_BH)),
            all(x$p_value >= 0 & x$p_value <= 1), all(x$q_value_BH >= 0 & x$q_value_BH <= 1))
  if (i == 1L) {
    stopifnot(all(x$direction %in% c("Up", "Down")), all(x$n_genes > 0), all(x$n_genes == as.integer(x$n_genes)))
  } else {
    fields <- c("estimate", "standard_error", "ci_low", "ci_high", "moderated_t", "df_total", "B_statistic")
    stopifnot(all(vapply(x[, ..fields], function(z) all(is.finite(z)), logical(1))),
              all(x$standard_error > 0), all(x$df_total > 0), all(x$ci_low <= x$estimate), all(x$estimate <= x$ci_high),
              max(abs(x$moderated_t - x$estimate / x$standard_error)) <= 1e-10,
              max(abs(x$ci_low - (x$estimate - qt(.975, x$df_total) * x$standard_error))) <= 1e-10,
              max(abs(x$ci_high - (x$estimate + qt(.975, x$df_total) * x$standard_error))) <= 1e-10)
  }
  for (contrast_name in unique(x$contrast)) {
    z <- x[contrast == contrast_name]
    q_delta <- max(abs(z$q_value_BH - p.adjust(z$p_value, "BH")))
    stopifnot(q_delta <= 1e-12)
    row <- family[layer == layers[[i]] & contrast == contrast_name]
    stopifnot(nrow(row) == 1L, row$tested == nrow(z), row$q_le_005 == sum(z$q_value_BH <= .05),
              abs(row$minimum_q - min(z$q_value_BH)) <= 1e-12)
    validation[[length(validation) + 1L]] <- data.table(layer = layers[[i]], contrast = contrast_name, rows = nrow(z),
      finite_fields = TRUE, native_BH_recalculation_max_error = q_delta, interval_and_t_consistent = TRUE)
  }
}
checks <- rbindlist(validation)
stopifnot(nrow(checks) == 16L, nrow(family) == 16L)
design <- fread(file.path(out, "DESIGN_ALL27.tsv"))
stopifnot(identical(design$sample_id, registry$sample_order_all27), qr(as.matrix(design[, -"sample_id"]))$rank == 4L)
fwrite(checks, file.path(out, "FINITE_FIELDS_AND_BH_CHECKS.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "METADATA_FINALIZER_SESSION.txt"))
stopifnot(identical(registry$inputs$sha256, unname(vapply(paths, sha, character(1)))))
output_names <- c(files, "FAMILY_SUMMARY.tsv", "ORIGINAL26_COMPATIBILITY.tsv", "ORIGINAL_MAPPING_COVERAGE.tsv",
                  "DESIGN_ALL27.tsv", "SESSION_INFO.txt", "FINITE_FIELDS_AND_BH_CHECKS.tsv", "METADATA_FINALIZER_SESSION.txt")
output_paths <- file.path(out, output_names)
helper <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), mustWork = TRUE, winslash = "/")
receipt <- list(status = "COMPLETE_EXPLORATORY_ORIGINAL_ALL27_CELL", utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  registry_sha256 = sha(registry_path), helper_sha256 = registry$helper_sha256, inputs_unchanged = TRUE,
  original26_comparisons_passed = TRUE, comparison_field_rows = nrow(compatibility), native_family_summary = family,
  output_identities = data.table(path = output_paths, bytes = file.info(output_paths)$size, sha256 = vapply(output_paths, sha, character(1))),
  native_fit_driver_exit_code = 1L, native_driver_failure_phase = "final_input_hash_vector_attribute_comparison",
  input_hash_values_all_equal = TRUE, input_bytes_all_equal = TRUE, failure_is_named_vs_unnamed_vector = TRUE,
  native_failure_preserved = list(log = file.path(out, "RUN.log"), sha256 = sha(file.path(out, "RUN.log"))),
  metadata_only_finalizer = list(path = helper, sha256 = sha(helper), invocation = commandArgs()),
  finalizer_model_rerun = FALSE, finalizer_activity_or_raw_matrix_read = FALSE,
  source_derived_matrix_reproduction_only = TRUE, raw_RMA_recomputed = FALSE, new_cell_model_fitted = TRUE,
  primary_or_gate_changed = FALSE, independent_validation = FALSE, whole_multiverse_or_goal_complete = FALSE)
write_json(receipt, file.path(out, "COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", digits = NA)
cat("METADATA_FINALIZATION_PASS: all20 input hashes/bytes unchanged; 41 original26 fields passed; all16 native BH/finite/CI families verified; no model rerun\n")
