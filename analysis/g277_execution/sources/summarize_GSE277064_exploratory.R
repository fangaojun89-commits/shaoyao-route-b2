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
program_root <- file.path(result_root, "programs")
manuscript_template_path <- file.path(
  project_root, "code", "templates", "GSE277064_manuscript_ready_methods_results.md"
)
manuscript_output_path <- file.path(result_root, "manuscript_ready_methods_results.md")

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

activity <- fread(file.path(program_root, "all_celltypes_activity_contrasts.tsv.gz"))
camera <- fread(file.path(program_root, "all_celltypes_Reactome_camera.tsv.gz"))
roast <- fread(file.path(program_root, "all_celltypes_Reactome_mroast.tsv.gz"))
fgsea <- fread(file.path(program_root, "all_celltypes_Reactome_fgseaMultilevel.tsv.gz"))
loso <- fread(file.path(program_root, "all_celltypes_chronic_activity_LOSO.tsv.gz"))
internal <- fread(file.path(program_root, "chronic_internal_consistency_not_S3.tsv.gz"))
s3_gate <- fread(file.path(program_root, "S3_GATE.tsv"))

stopifnot(
  nrow(s3_gate) == 1L,
  s3_gate$eligible_input_n == 0L,
  s3_gate$status == "NOT_IDENTIFIABLE_NO_S2_INPUT"
)

direction_code <- function(x) fifelse(x == "Up", 1L, fifelse(x == "Down", -1L, NA_integer_))

reactome_base <- internal[collection == "Reactome"]
reactome_base <- merge(
  reactome_base,
  camera[contrast == "chronic_combined", .(
    cell_type, program_id, term_name,
    camera_direction_full = direction,
    camera_p_full = p_value,
    camera_q_global = q_global_chronic_all_locked_precompute
  )],
  by = c("cell_type", "program_id"), all.x = TRUE
)
reactome_base <- merge(
  reactome_base,
  roast[contrast == "chronic_combined", .(
    cell_type, program_id,
    mroast_direction = direction,
    mroast_p = p_value,
    mroast_q_global = q_global_chronic_all_locked_precompute,
    mroast_p_mixed = p_value_mixed
  )],
  by = c("cell_type", "program_id"), all.x = TRUE
)
reactome_base <- merge(
  reactome_base,
  fgsea[contrast == "chronic_combined", .(
    cell_type, program_id,
    fgsea_ES = ES,
    fgsea_NES = NES,
    fgsea_p = p_value,
    fgsea_q_global = q_global_chronic_all_locked_precompute
  )],
  by = c("cell_type", "program_id"), all.x = TRUE
)
reactome_base <- merge(
  reactome_base,
  loso[collection == "Reactome", .(
    cell_type, program_id,
    loso_direction_fraction = direction_fraction,
    loso_sign_flip_any = sign_flip_any,
    loso_leave_one_min = leave_one_min,
    loso_leave_one_max = leave_one_max,
    loso_max_abs_delta = max_abs_delta
  )],
  by = c("cell_type", "program_id"), all.x = TRUE
)
reactome_base[, activity_sign := sign(effect_chronic_combined)]
reactome_base[, camera_sign := direction_code(camera_direction_full)]
reactome_base[, mroast_sign := direction_code(mroast_direction)]
reactome_base[, fgsea_sign := sign(fgsea_NES)]
reactome_base[, all_method_direction_match :=
                 activity_sign == camera_sign & activity_sign == mroast_sign & activity_sign == fgsea_sign]
reactome_base[, maximum_global_q := pmax(
  q_global_chronic_all_locked_precompute_chronic_combined,
  camera_q_global,
  mroast_q_global,
  fgsea_q_global,
  na.rm = FALSE
)]
reactome_base[, strict_all_locked_exploratory_support :=
                 q_global_chronic_all_locked_precompute_chronic_combined <= 0.05 &
                 combined_ci_excludes_zero &
                 !D14_D28_opposite_signs &
                 camera_q_global <= 0.05 &
                 mroast_q_global <= 0.05 &
                 fgsea_q_global <= 0.05 &
                 all_method_direction_match &
                 loso_direction_fraction >= 0.80]
reactome_base[, gate_eligible := FALSE]
reactome_base[, S3_status := "NOT_IDENTIFIABLE_NO_S2_INPUT"]
reactome_base[, absolute_chronic_effect := abs(effect_chronic_combined)]
setorder(reactome_base, cell_type, -strict_all_locked_exploratory_support, maximum_global_q, -absolute_chronic_effect)
write_tsv(reactome_base, file.path(program_root, "exploratory_Reactome_multimethod_consensus_not_S3.tsv.gz"))

activity_base <- internal[collection != "Reactome"]
activity_base <- merge(
  activity_base,
  loso[collection != "Reactome", .(
    collection, cell_type, program_id,
    loso_direction_fraction = direction_fraction,
    loso_sign_flip_any = sign_flip_any,
    loso_leave_one_min = leave_one_min,
    loso_leave_one_max = leave_one_max,
    loso_max_abs_delta = max_abs_delta
  )],
  by = c("collection", "cell_type", "program_id"), all.x = TRUE
)
activity_base[, strict_all_locked_exploratory_support :=
                q_global_chronic_all_locked_precompute_chronic_combined <= 0.05 &
                combined_ci_excludes_zero &
                !D14_D28_opposite_signs &
                loso_direction_fraction >= 0.80]
activity_base[, gate_eligible := FALSE]
activity_base[, S3_status := "NOT_IDENTIFIABLE_NO_S2_INPUT"]
activity_base[, absolute_chronic_effect := abs(effect_chronic_combined)]
setorder(activity_base, collection, cell_type, -strict_all_locked_exploratory_support,
         q_global_chronic_all_locked_precompute_chronic_combined, -absolute_chronic_effect)
write_tsv(activity_base, file.path(program_root, "exploratory_activity_consistency_not_S3.tsv.gz"))

stage_table <- merge(
  activity[contrast == "chronic_combined", .(
    collection, cell_type, program_id,
    chronic_effect = effect,
    chronic_ci_low = ci_low,
    chronic_ci_high = ci_high,
    chronic_p = p_value,
    chronic_q_global = q_global_chronic_all_locked_precompute
  )],
  activity[contrast == "early_combined", .(
    collection, cell_type, program_id,
    early_effect = effect,
    early_ci_low = ci_low,
    early_ci_high = ci_high,
    early_p = p_value,
    early_q_within = q_within_celltype_contrast
  )],
  by = c("collection", "cell_type", "program_id"), all = TRUE
)
stage_table[, chronic_early_same_direction := sign(chronic_effect) == sign(early_effect)]
stage_table[, absolute_chronic_to_early_ratio := abs(chronic_effect) / pmax(abs(early_effect), 1e-8)]
stage_table[, gate_eligible := FALSE]
stage_table[, S3_status := "NOT_IDENTIFIABLE_NO_S2_INPUT"]
write_tsv(stage_table, file.path(program_root, "exploratory_chronic_vs_early_activity_not_S3.tsv.gz"))

consensus_counts <- rbindlist(list(
  reactome_base[, .(
    programs_tested = .N,
    strict_all_locked_exploratory_support_n = sum(strict_all_locked_exploratory_support, na.rm = TRUE),
    strict_support_with_no_LOSO_sign_flip_n = sum(strict_all_locked_exploratory_support & !loso_sign_flip_any, na.rm = TRUE)
  ), by = .(collection, cell_type)],
  activity_base[, .(
    programs_tested = .N,
    strict_all_locked_exploratory_support_n = sum(strict_all_locked_exploratory_support, na.rm = TRUE),
    strict_support_with_no_LOSO_sign_flip_n = sum(strict_all_locked_exploratory_support & !loso_sign_flip_any, na.rm = TRUE)
  ), by = .(collection, cell_type)]
), fill = TRUE)
consensus_counts[, eligible_S2_input_n := 0L]
consensus_counts[, S3_pass_n := NA_integer_]
consensus_counts[, S3_status := "NOT_IDENTIFIABLE_NO_S2_INPUT"]
write_tsv(consensus_counts, file.path(program_root, "exploratory_consensus_counts_not_S3.tsv"))

top_reactome <- reactome_base[strict_all_locked_exploratory_support == TRUE,
                               head(.SD, 20L), by = cell_type]
top_activity <- activity_base[strict_all_locked_exploratory_support == TRUE,
                              head(.SD, 20L), by = .(collection, cell_type)]
write_tsv(top_reactome, file.path(program_root, "exploratory_top20_Reactome_per_celltype_not_S3.tsv"))
write_tsv(top_activity, file.path(program_root, "exploratory_top20_activity_per_collection_celltype_not_S3.tsv"))

boundary_activity <- fread(file.path(result_root, "boundary_rotarod", "programs", "all_celltypes_activity_contrasts.tsv.gz"))
boundary_counts <- boundary_activity[, .(
  programs_tested = .N,
  within_contrast_q_0_05 = sum(q_within_celltype_contrast <= 0.05, na.rm = TRUE),
  positive_q_0_05 = sum(q_within_celltype_contrast <= 0.05 & effect > 0, na.rm = TRUE),
  negative_q_0_05 = sum(q_within_celltype_contrast <= 0.05 & effect < 0, na.rm = TRUE)
), by = .(collection, cell_type, contrast)]
boundary_counts[, evidence_role := "rotarod_boundary_only_not_S_or_E_gate"]
write_tsv(boundary_counts, file.path(result_root, "boundary_rotarod", "programs", "boundary_activity_result_counts.tsv"))

warning_adjudication <- data.table(
  warning = c(
    "AnnotationDbi one-to-many mapping message",
    "fgsea exact preranked-statistic ties",
    "empty mroast/fgsea files in rotarod branch",
    "analysis-plan registration terminology"
  ),
  adjudication = c(
    "Expected for some Ensembl identifiers; a fixed multiVals=first rule was applied and the full mapping table is retained",
    "Resolved by deterministic gene-symbol machine-epsilon tie breaking; final fgsea tables were refreshed and the before/after counts are retained",
    "Resolved by explicit NOT_RUN status tables because these sensitivity methods are restricted to the natural-course chronic-combined contrast",
    "The workflow was internally specified and locally locked; there was no public preregistration. Correcting this terminology does not alter any analysis or gate decision"
  ),
  unresolved = FALSE
)
write_tsv(warning_adjudication, file.path(result_root, "diagnostics", "warning_adjudication.tsv"))

report_lines <- c(
  "# GSE277064 exploratory result synopsis (not S3)",
  "",
  "The frozen S1 set contains zero programs; therefore the eligible S2 input is empty and S3 is formally NOT_IDENTIFIABLE_NO_S2_INPUT. The tables below summarize all-program post-stroke cell-type/time behavior only and cannot create candidates.",
  "",
  "## Analysis-plan provenance",
  "",
  "Prespecified and locked denote the internally specified, locally frozen workflow; they do not imply public preregistration. The terminology correction changes no sample inclusion, contrast, method, threshold, numerical result or S3 gate decision.",
  "",
  "## Strict all-locked exploratory counts",
  "",
  paste(capture.output(print(consensus_counts)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "Strict exploratory support requires all-program global BH q<=0.05, a chronic activity CI excluding zero, non-opposite D14/D28 activity signs and >=80% leave-one-library direction stability; Reactome additionally requires camera, 9,999-rotation mroast and deterministic-tie fgsea direction/FDR agreement. This deliberately stringent descriptor is not an evidence gate and is not a post-stroke-spasticity validation.",
  "",
  "Rotarod results are maintained in a separate boundary branch and cannot enter S1-S3 or E1-E4."
)
writeLines(report_lines, file.path(result_root, "exploratory_result_synopsis_not_S3.md"), useBytes = TRUE)

if (!file.exists(manuscript_template_path)) {
  stop("The release-frozen GSE277064 manuscript-ready template is missing", call. = FALSE)
}
if (file.exists(manuscript_output_path)) {
  stop("Refusing to overwrite an existing GSE277064 manuscript-ready result", call. = FALSE)
}
if (!isTRUE(file.copy(manuscript_template_path, manuscript_output_path, overwrite = FALSE))) {
  stop("Could not materialize the release-frozen GSE277064 manuscript-ready template", call. = FALSE)
}

manifest_path <- file.path(result_root, "audit", "run_manifest.tsv")
manifest <- fread(manifest_path)
tracked_manifest_keys <- c(
  "exploratory_summary_script_sha256", "manuscript_ready_text_sha256",
  "manuscript_ready_text_source_sha256", "manuscript_ready_text_role",
  "analysis_summary_sha256", "S3_gate_sha256", "analysis_plan_provenance",
  "public_preregistration", "terminology_correction_scope", "result_file_count"
)
manifest <- manifest[!manifest_key %chin% tracked_manifest_keys]
manifest <- rbind(manifest, data.table(
  manifest_key = tracked_manifest_keys,
  value = c(
    digest(script_path, algo = "sha256", file = TRUE),
    digest(manuscript_output_path, algo = "sha256", file = TRUE),
    digest(manuscript_template_path, algo = "sha256", file = TRUE),
    "release-frozen narrative template; not a statistical input or independent analytical result",
    digest(file.path(result_root, "analysis_summary.md"), algo = "sha256", file = TRUE),
    digest(file.path(result_root, "programs", "S3_GATE.tsv"), algo = "sha256", file = TRUE),
    "internally specified and locally locked",
    "none",
    "wording only; no analysis or S3 gate decision changed",
    as.character(length(list.files(result_root, recursive = TRUE, full.names = TRUE)))
  )
))
write_tsv(manifest, manifest_path)

cat("Wrote exploratory summaries with S3 disabled by empty S2 input.\n")
