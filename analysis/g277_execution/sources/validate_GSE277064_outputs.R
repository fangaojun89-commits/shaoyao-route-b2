#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(file_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

result_root <- file.path(project_root, "results", "GSE277064")
checks <- list()
add_check <- function(id, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.table(
    check_id = id,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

audit <- fread(file.path(result_root, "audit", "analysis_audit_summary.tsv"))
add_check("input_audit_all_pass", all(audit$status == "PASS"), paste(unique(audit$status), collapse = ","), "PASS")

models <- fread(file.path(result_root, "diagnostics", "all_celltypes_model_summary.tsv"))
add_check("three_natural_models", nrow(models) == 3L, nrow(models), 3L)
add_check("natural_models_n15", all(models$n_samples == 15L), paste(models$n_samples, collapse = ","), "15,15,15")
add_check("natural_models_full_rank", all(models$design_rank == models$design_columns),
          paste(models$design_rank, models$design_columns, sep = "/", collapse = ","), "5/5 each")
add_check("natural_models_residual_df10", all(models$residual_df == 10L), paste(models$residual_df, collapse = ","), "10,10,10")

qc <- fread(file.path(result_root, "diagnostics", "all_celltypes_sample_QC.tsv"))
add_check("natural_QC_45_samples", nrow(qc) == 45L, nrow(qc), 45L)
add_check("no_technical_exclusions", !any(qc$exclusion), sum(qc$exclusion), 0L)
add_check("no_objective_metric_flags", !any(qc$objective_metric_flag), sum(qc$objective_metric_flag), 0L)

activity <- fread(file.path(result_root, "programs", "all_celltypes_activity_contrasts.tsv.gz"))
camera <- fread(file.path(result_root, "programs", "all_celltypes_Reactome_camera.tsv.gz"))
roast <- fread(file.path(result_root, "programs", "all_celltypes_Reactome_mroast.tsv.gz"))
fgsea <- fread(file.path(result_root, "programs", "all_celltypes_Reactome_fgseaMultilevel.tsv.gz"))
loso <- fread(file.path(result_root, "programs", "all_celltypes_chronic_activity_LOSO.tsv.gz"))
internal <- fread(file.path(result_root, "programs", "chronic_internal_consistency_not_S3.tsv.gz"))

add_check("activity_row_count", nrow(activity) == 29322L, nrow(activity), 29322L)
add_check("camera_row_count", nrow(camera) == 27654L, nrow(camera), 27654L)
add_check("mroast_chronic_row_count", nrow(roast) == 4609L, nrow(roast), 4609L)
add_check("fgsea_chronic_row_count", nrow(fgsea) == 4609L, nrow(fgsea), 4609L)
add_check("LOSO_program_row_count", nrow(loso) == 4887L, nrow(loso), 4887L)
add_check("internal_program_row_count", nrow(internal) == 4887L, nrow(internal), 4887L)

probability_columns <- list(
  activity = activity$p_value,
  camera = camera$p_value,
  roast = roast$p_value,
  fgsea = fgsea$p_value
)
for (table_name in names(probability_columns)) {
  values <- probability_columns[[table_name]]
  add_check(
    paste0(table_name, "_p_values_valid"),
    all(is.finite(values) & values >= 0 & values <= 1),
    sum(!is.finite(values) | values < 0 | values > 1),
    0L
  )
}

tie_audit <- fread(file.path(result_root, "diagnostics", "fgsea_deterministic_tie_audit.tsv"))
add_check("fgsea_ties_resolved", all(tie_audit$exact_ties_after == 0L),
          paste(tie_audit$exact_ties_after, collapse = ","), "0,0,0")
add_check("fgsea_method_labeled_deterministic",
          identical(unique(fgsea$method), "fgseaMultilevel_deterministic_ties"),
          paste(unique(fgsea$method), collapse = ","), "fgseaMultilevel_deterministic_ties")

gate <- fread(file.path(result_root, "programs", "S3_GATE.tsv"))
add_check("S3_gate_single_row", nrow(gate) == 1L, nrow(gate), 1L)
add_check("S3_eligible_input_zero", gate$eligible_input_n == 0L, gate$eligible_input_n, 0L)
add_check("S3_pass_is_NA", is.na(gate$passed_candidate_n), gate$passed_candidate_n, "NA")
gate_reason <- as.character(gate$reason)
gate_status_and_provenance_ok <- gate$status == "NOT_IDENTIFIABLE_NO_S2_INPUT" &&
  grepl("internally specified and locally locked", gate_reason, fixed = TRUE) &&
  grepl("not a public preregistration", gate_reason, fixed = TRUE) &&
  !grepl("registered sequential", gate_reason, fixed = TRUE)
add_check("S3_not_identifiable_status_internal_lock", gate_status_and_provenance_ok,
          paste(gate$status, gate_reason, sep = " | "),
          "NOT_IDENTIFIABLE_NO_S2_INPUT with internal-lock/no-public-preregistration provenance")
add_check("all_internal_rows_not_S3", all(internal$S3_status == "NOT_IDENTIFIABLE_NO_S2_INPUT"),
          paste(unique(internal$S3_status), collapse = ","), "NOT_IDENTIFIABLE_NO_S2_INPUT")
add_check("no_internal_gate_eligibility", !any(internal$upstream_disease_direction_checked),
          sum(internal$upstream_disease_direction_checked), 0L)

for (cell_type in c("v2a", "astrocyte", "microglia")) {
  object <- readRDS(file.path(result_root, "objects", paste0(cell_type, "_analysis.rds")))
  add_check(
    paste0(cell_type, "_object_15_samples"),
    nrow(object$sample_metadata) == 15L && ncol(object$dge$counts) == 15L,
    paste(nrow(object$sample_metadata), ncol(object$dge$counts), sep = "/"),
    "15/15"
  )
  add_check(
    paste0(cell_type, "_object_design_rank"),
    qr(object$design)$rank == 5L,
    qr(object$design)$rank,
    5L
  )
  add_check(
    paste0(cell_type, "_fgsea_object_refreshed"),
    !is.null(object$fgsea_tie_break) && object$fgsea_tie_break$ties_after == 0L,
    if (is.null(object$fgsea_tie_break)) "missing" else object$fgsea_tie_break$ties_after,
    0L
  )
}

for (cell_type in c("v2a", "astrocyte")) {
  object <- readRDS(file.path(result_root, "boundary_rotarod", "objects", paste0(cell_type, "_analysis.rds")))
  add_check(
    paste0(cell_type, "_boundary_object_12_samples"),
    nrow(object$sample_metadata) == 12L && ncol(object$dge$counts) == 12L,
    paste(nrow(object$sample_metadata), ncol(object$dge$counts), sep = "/"),
    "12/12"
  )
  add_check(
    paste0(cell_type, "_boundary_separate_domain"),
    identical(object$domain, "rotarod_boundary"), object$domain, "rotarod_boundary"
  )
}

validation_output_path <- file.path(result_root, "diagnostics", "output_validation.tsv")
all_files <- list.files(result_root, recursive = TRUE, full.names = TRUE)
file_sizes <- file.info(all_files)$size
prospective_file_count <- length(all_files) + as.integer(!file.exists(validation_output_path))
add_check("result_file_set_169_nonempty", prospective_file_count == 169L && all(file_sizes > 0),
          paste0("files=", prospective_file_count, ";zero_bytes=", sum(file_sizes == 0)),
          "files=169;zero_bytes=0")

manifest <- fread(file.path(result_root, "audit", "run_manifest.tsv"))
validator_sha <- digest(script_path, algo = "sha256", file = TRUE)
if (any(manifest$manifest_key == "output_validation_script_sha256")) {
  manifest[manifest_key == "output_validation_script_sha256", value := validator_sha]
} else {
  manifest <- rbind(manifest, data.table(
    manifest_key = "output_validation_script_sha256",
    value = validator_sha
  ))
}
fwrite(manifest, file.path(result_root, "audit", "run_manifest.tsv"),
       sep = "\t", quote = FALSE, na = "NA")
manifest_value <- function(key) manifest[manifest_key == key, value][1L]
main_sha <- digest(file.path(project_root, "code", "R", "analyze_GSE277064.R"), algo = "sha256", file = TRUE)
refresh_sha <- digest(file.path(project_root, "code", "R", "refresh_GSE277064_fgsea.R"), algo = "sha256", file = TRUE)
summary_sha <- digest(file.path(project_root, "code", "R", "summarize_GSE277064_exploratory.R"), algo = "sha256", file = TRUE)
add_check("main_script_hash_matches", identical(manifest_value("script_sha256"), main_sha),
          manifest_value("script_sha256"), main_sha)
add_check("fgsea_refresh_hash_matches", identical(manifest_value("fgsea_refresh_script_sha256"), refresh_sha),
          manifest_value("fgsea_refresh_script_sha256"), refresh_sha)
add_check("summary_script_hash_matches", identical(manifest_value("exploratory_summary_script_sha256"), summary_sha),
          manifest_value("exploratory_summary_script_sha256"), summary_sha)
add_check("roast_rotations_9999", identical(manifest_value("roast_nrot"), "9999"),
          manifest_value("roast_nrot"), "9999")

validation <- rbindlist(checks)
validation[, status := fifelse(passed, "PASS", "FAIL")]
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
write_tsv(validation, validation_output_path)

final_files <- list.files(result_root, recursive = TRUE, full.names = TRUE)
final_file_sizes <- file.info(final_files)$size
if (length(final_files) != 169L || any(final_file_sizes <= 0)) {
  stop("GSE277064 post-validation result inventory failed", call. = FALSE)
}

if (!all(validation$passed)) {
  print(validation[passed == FALSE])
  stop("GSE277064 output validation failed", call. = FALSE)
}

cat("GSE277064 output validation:", nrow(validation), "of", nrow(validation), "checks PASS\n")
