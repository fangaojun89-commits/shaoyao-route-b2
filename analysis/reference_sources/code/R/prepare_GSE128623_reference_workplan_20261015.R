# Deterministic work map from the already frozen label registry. No model/data fit.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
project <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
registry_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- gsub("\\", "/", args[[3L]], fixed = TRUE)
stopifnot(grepl("^R:/Shaoyao_storage/cache/GSE128623_[A-Za-z0-9_-]+$", output),
          !file.exists(output), dir.exists(dirname(output)))
prepared <- file.path(registry_root, "PREPARED.json")
prepared_sha <- "63be688e2d548b2917479a9ecc6e1dba4b3ab37542ac8c4aca9244e78c3ef5d0"
stopifnot(identical(digest::digest(file = prepared, algo = "sha256"), prepared_sha))
receipt <- jsonlite::read_json(prepared, simplifyVector = FALSE)
proposal <- file.path(project, "protocol/GSE128623_FULL_LABEL_REFERENCE_EXECUTION_PROPOSAL_20261015.md")
stopifnot(identical(digest::digest(file = proposal, algo = "sha256"),
  "15fec7caee23ac4eeb12af6a7363c5590aa14a001300a8686aa203be3b5096b4"))
rows <- list()
for (processing in c("raw", "adapter")) for (cohort in c("all30", "predefined29")) {
  record <- receipt$cohorts[[cohort]]
  spaces <- c(list(six_group = record$sampled_space),
    setNames(record$exact_day_control_spaces,
      paste0(names(record$exact_day_control_spaces), "_vs_control_exact")))
  for (space in names(spaces)) {
    item <- spaces[[space]]
    path <- file.path(registry_root, item$file)
    stopifnot(identical(digest::digest(file = path, algo = "sha256"), item$sha256))
    full <- identical(space, "six_group")
    n <- as.integer(item$rows)
    stopifnot(n == if (full) 10001L else if (cohort == "predefined29" &&
      identical(space, "D4_vs_control_exact")) 126L else 252L)
    rows[[length(rows) + 1L]] <- data.frame(cell_id = paste(processing, cohort, sep = "_"),
      processing = processing, cohort = cohort, reference_space = space,
      kind = if (full) "sampled" else "exhaustive",
      target_contrast = if (full) "NOT_APPLICABLE_FULL_SEVEN" else sub("_exact$", "", space),
      allocation_file = item$file, allocation_file_sha256 = item$sha256,
      assignment_order_sha256 = item$assignment_order_sha256,
      required_archive_n = n, observed_allocation_id = 1L,
      observed_reference_view_reuses_preserved_cell_fit = TRUE,
      nonidentity_joint_gene_camera_calls = n - 1L,
      nonidentity_joint_activity_calls = n - 1L,
      nonidentity_unsigned_calls = if (full) n - 1L else 0L,
      family_n = if (full) 67L else 9L,
      stringsAsFactors = FALSE)
  }
}
workplan <- do.call(rbind, rows)
stopifnot(nrow(workplan) == 24L, !anyDuplicated(workplan[c("cell_id", "reference_space")]),
  length(unique(workplan$cell_id)) == 4L, all(table(workplan$cell_id) == 6L),
  sum(workplan$required_archive_n) == 44792L,
  sum(workplan$nonidentity_joint_gene_camera_calls) == 44768L,
  sum(workplan$nonidentity_joint_activity_calls) == 44768L,
  sum(workplan$nonidentity_unsigned_calls) == 40000L,
  sum(workplan$family_n[workplan$kind == "sampled"]) == 268L,
  sum(workplan$family_n[workplan$kind == "exhaustive"]) == 180L)
for (cohort in c("all30", "predefined29")) for (space in unique(workplan$reference_space)) {
  paired <- workplan[workplan$cohort == cohort & workplan$reference_space == space, ]
  stopifnot(nrow(paired) == 2L, length(unique(paired$allocation_file_sha256)) == 1L,
    length(unique(paired$assignment_order_sha256)) == 1L)
}
stopifnot(dir.create(output))
write.table(workplan, file.path(output, "REFERENCE_WORKPLAN.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(list(status = "REGISTRY_DERIVED_WORKPLAN_ONLY_NOT_ADMITTED_JOBS",
  registry_prepared_sha256 = prepared_sha,
  plan_sha256 = digest::digest(file = file.path(output, "REFERENCE_WORKPLAN.tsv"), algo = "sha256"),
  processing_cohort_cells = 4L, reference_jobs = 24L, planned_atomic_archives = 44792L,
  nonidentity_joint_gene_camera_calls = 44768L, nonidentity_joint_activity_calls = 44768L,
  nonidentity_unsigned_calls = 40000L, preserved_observed_cell_fits_required = 4L,
  observed_reference_views = 24L, six_group_native_families = 268L,
  companion_reference_family_slots = 180L, actual_model_calls = 0L,
  job_input_preparation_complete = FALSE, observed_bridge_complete = FALSE,
  resource_budget_accepted = FALSE, scientific_execution_authority_granted = FALSE),
  file.path(output, "WORKPLAN_PREPARED.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(output, "SESSION_INFO.txt"))
cat("PREPARED 24 registry-derived work-map rows; 44792 planned archives; zero fits\n")
