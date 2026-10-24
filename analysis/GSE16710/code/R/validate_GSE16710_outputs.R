#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
result_dir <- file.path(project_root, "results", "GSE16710")
program_dir <- file.path(result_dir, "program_level")
gene_dir <- file.path(result_dir, "gene_level")
diagnostic_dir <- file.path(result_dir, "diagnostics")
processed_dir <- file.path(project_root, "data", "processed", "GSE16710")
adjudication_source <- Sys.getenv("GSE16710_ADJUDICATION_SOURCE", unset = "CANONICAL_POST_QC_RECORD")
if (!adjudication_source %in% c("CANONICAL_POST_QC_RECORD", "RELEASE_FROZEN")) {
  stop("Unsupported GSE16710_ADJUDICATION_SOURCE")
}
adjudication_dir <- if (adjudication_source == "RELEASE_FROZEN") {
  file.path(project_root, "data", "frozen_adjudication", "GSE16710")
} else {
  processed_dir
}

checks <- list()
add_check <- function(name, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    status = if (isTRUE(pass)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

ledger <- fread(file.path(processed_dir, "sample_ledger_verified.tsv"))
decision <- fread(file.path(adjudication_dir, "qc_inclusion_decision.tsv"))
qc <- fread(file.path(result_dir, "qc", "array_qc_metrics_and_flags.tsv"))
model <- fread(file.path(diagnostic_dir, "model_estimability.tsv"))
gene_counts <- fread(file.path(gene_dir, "gene_interaction_summary_counts.tsv"))
camera <- fread(file.path(program_dir, "reactome_camera_primary_all_contrasts.tsv"))
mroast <- fread(file.path(program_dir, "reactome_mroast_interaction_sensitivity.tsv"))
fgsea <- fread(file.path(program_dir, "reactome_fgsea_interaction_sensitivity.tsv"))
react_s1 <- fread(file.path(program_dir, "S1_gate_Reactome.tsv"))
prog_s1 <- fread(file.path(program_dir, "S1_gate_PROGENy.tsv"))
doro_s1 <- fread(file.path(program_dir, "S1_gate_DoRothEA.tsv"))
s1 <- fread(file.path(program_dir, "S1_gate_all_collections.tsv"))
s1_machine <- fread(file.path(program_dir, "S1_GATE.tsv"))
perm <- fread(file.path(diagnostic_dir, "program_Freedman_Lane_10000.tsv"))
cal <- fread(file.path(diagnostic_dir, "program_Freedman_Lane_null_calibration.tsv"))
all27 <- fread(file.path(program_dir, "reactome_camera_all27_sensitivity.tsv"))
orig <- fread(file.path(program_dir, "reactome_camera_original_annotation_sensitivity.tsv"))
all27_progeny <- fread(file.path(program_dir, "progeny_all27_sensitivity.tsv"))
all27_dorothea <- fread(file.path(program_dir, "dorothea_all27_sensitivity.tsv"))
orig_progeny <- fread(file.path(program_dir, "progeny_original_annotation_sensitivity.tsv"))
orig_dorothea <- fread(file.path(program_dir, "dorothea_original_annotation_sensitivity.tsv"))
loo_reactome <- fread(file.path(diagnostic_dir, "reactome_leave_one_animal_out_long.tsv.gz"))
loo_progeny <- fread(file.path(diagnostic_dir, "progeny_leave_one_animal_out_long.tsv.gz"))
loo_dorothea <- fread(file.path(diagnostic_dir, "dorothea_leave_one_animal_out_long.tsv.gz"))
precision <- fread(file.path(diagnostic_dir, "interaction_precision_summary.tsv"))
quickgo_status <- fread(file.path(diagnostic_dir, "quickgo_assignedBy_SynGO_sensitivity_status.tsv"))
quickgo_camera <- fread(file.path(program_dir, "quickgo_assignedBy_SynGO_camera_sensitivity.tsv"))
quickgo_mroast <- fread(file.path(program_dir, "quickgo_assignedBy_SynGO_mroast_sensitivity.tsv"))
quickgo_fgsea <- fread(file.path(program_dir, "quickgo_assignedBy_SynGO_fgsea_sensitivity.tsv"))
quickgo_activity <- fread(file.path(program_dir, "quickgo_assignedBy_SynGO_activity_sensitivity.tsv"))

add_check("official_array_count", nrow(ledger) == 27L, nrow(ledger), 27L)
add_check("unique_biological_units", uniqueN(ledger$biological_unit_id) == 27L,
          uniqueN(ledger$biological_unit_id), 27L)
add_check("frozen_primary_sample_count", sum(decision$primary_analysis_decision == "retain") == 26L,
          sum(decision$primary_analysis_decision == "retain"), 26L)
primary_ids <- decision[primary_analysis_decision == "retain", sample_id]
primary_cells <- ledger[sample_id %in% primary_ids,
                        .N, by = .(experimental_arm, time_day)][order(experimental_arm, time_day)]
observed_cells <- paste0(primary_cells$experimental_arm, primary_cells$time_day, "=", primary_cells$N,
                         collapse = ";")
add_check("primary_cell_sizes", identical(observed_cells, "Sham21=6;Sham60=4;Spi21=8;Spi60=8"),
          observed_cells, "Sham21=6;Sham60=4;Spi21=8;Spi60=8")
add_check("outcome_blind_exclusion_identity",
          identical(decision[primary_analysis_decision != "retain", sample_id], "GSM418682"),
          paste(decision[primary_analysis_decision != "retain", sample_id], collapse = ";"), "GSM418682")
add_check("excluded_array_corroborated_QC_domains",
          qc[sample_id == "GSM418682", independent_technical_domains_flagged] >= 2L,
          qc[sample_id == "GSM418682", independent_technical_domains_flagged], ">=2")
add_check("GSM418682_exclusion_direction_free",
          nrow(decision) == 27L &&
            identical(decision[primary_analysis_decision != "retain", sample_id], "GSM418682") &&
            all(toupper(as.character(decision$outcome_information_used)) == "FALSE") &&
            all(decision$decision_stage == "outcome_blind_qc"),
          paste0("outcome_information_used=", paste(unique(decision$outcome_information_used), collapse = ";")),
          "FALSE for all frozen decision rows; internal process record, not external proof")
add_check("primary_model_full_rank", model[analysis == "primary_outcome_blind_QC_set", model_estimable] == TRUE,
          model[analysis == "primary_outcome_blind_QC_set", paste0(design_rank, "/", design_columns)], "4/4")
add_check("primary_residual_df", model[analysis == "primary_outcome_blind_QC_set", residual_df] == 22L,
          model[analysis == "primary_outcome_blind_QC_set", residual_df], 22L)
add_check("gene_test_count", gene_counts$n_tested == 15375L, gene_counts$n_tested, 15375L)
add_check("gene_BH_hits_zero", gene_counts$n_BH_q_le_0_05 == 0L, gene_counts$n_BH_q_le_0_05, 0L)

cam_int <- camera[contrast == "spasticity_progression_interaction"]
camera_q_recalc <- p.adjust(cam_int$p_value, method = "BH")
add_check("reactome_camera_program_count", nrow(cam_int) == 1656L, nrow(cam_int), 1656L)
add_check("reactome_camera_BH_reproducible", isTRUE(all.equal(camera_q_recalc, cam_int$q_value_BH, tolerance = 1e-12)),
          max(abs(camera_q_recalc - cam_int$q_value_BH)), "<=1e-12")
add_check("reactome_camera_primary_hits_zero", sum(cam_int$q_value_BH <= 0.05) == 0L,
          sum(cam_int$q_value_BH <= 0.05), 0L)
add_check("reactome_camera_min_q", abs(min(cam_int$q_value_BH) - 0.105978073883113) < 1e-12,
          min(cam_int$q_value_BH), 0.105978073883113)
add_check("mroast_directional_hits_zero", sum(mroast$FDR <= 0.05) == 0L,
          sum(mroast$FDR <= 0.05), 0L)
add_check("mroast_mixed_hits_zero", sum(mroast$FDR.Mixed <= 0.05) == 0L,
          sum(mroast$FDR.Mixed <= 0.05), 0L)
add_check("fgsea_sensitivity_hits_recorded", sum(fgsea$padj <= 0.05) == 108L,
          sum(fgsea$padj <= 0.05), 108L)

add_check("S1_reactome_count", nrow(react_s1) == 1656L, nrow(react_s1), 1656L)
add_check("S1_progeny_count", nrow(prog_s1) == 14L, nrow(prog_s1), 14L)
add_check("S1_dorothea_count", nrow(doro_s1) == 101L, nrow(doro_s1), 101L)
add_check("S1_total_count", nrow(s1) == 1771L, nrow(s1), 1771L)
add_check("S1_pass_zero", sum(s1$S1_pass, na.rm = TRUE) == 0L, sum(s1$S1_pass, na.rm = TRUE), 0L)
add_check("S1_GATE_machine_readable_matches", isTRUE(all.equal(s1, s1_machine, check.attributes = FALSE)),
          nrow(s1_machine), paste0(nrow(s1), " identical rows"))
add_check("fgsea_not_used_as_S1_gate", !"fgsea_padj" %in% names(s1),
          paste(names(s1), collapse = ";"), "no fgsea gate column")

add_check("Freedman_Lane_program_count", nrow(perm) == 1771L, nrow(perm), 1771L)
add_check("Freedman_Lane_10000_valid_each", all(perm$valid_permutations == 10000L),
          paste(range(perm$valid_permutations), collapse = "-"), "10000-10000")
add_check("Freedman_Lane_BH_hits_zero", sum(perm$empirical_q_value_BH <= 0.05) == 0L,
          sum(perm$empirical_q_value_BH <= 0.05), 0L)
add_check("Freedman_Lane_not_anti_conservative", identical(cal$anti_conservative_flag, FALSE),
          cal$anti_conservative_flag, FALSE)
add_check("all27_camera_influence_recorded", sum(all27$q_value_BH <= 0.05) == 15L,
          sum(all27$q_value_BH <= 0.05), 15L)
add_check("original_annotation_camera_hits_zero", sum(orig$q_value_BH <= 0.05) == 0L,
          sum(orig$q_value_BH <= 0.05), 0L)
all27_prog_int <- all27_progeny[contrast == "spasticity_progression_interaction"]
all27_doro_int <- all27_dorothea[contrast == "spasticity_progression_interaction"]
orig_prog_int <- orig_progeny[contrast == "spasticity_progression_interaction"]
orig_doro_int <- orig_dorothea[contrast == "spasticity_progression_interaction"]
add_check("all27_PROGENy_complete", nrow(all27_prog_int) == 14L,
          nrow(all27_prog_int), 14L)
add_check("all27_DoRothEA_complete", nrow(all27_doro_int) == 101L,
          nrow(all27_doro_int), 101L)
add_check("original_annotation_PROGENy_complete", nrow(orig_prog_int) == 14L,
          nrow(orig_prog_int), 14L)
add_check("original_annotation_DoRothEA_mapped_count", nrow(orig_doro_int) == 100L,
          nrow(orig_doro_int), 100L)
add_check("activity_all27_and_original_hits_zero",
          sum(all27_prog_int$q_value_BH <= 0.05, na.rm = TRUE) == 0L &&
            sum(all27_doro_int$q_value_BH <= 0.05, na.rm = TRUE) == 0L &&
            sum(orig_prog_int$q_value_BH <= 0.05, na.rm = TRUE) == 0L &&
            sum(orig_doro_int$q_value_BH <= 0.05, na.rm = TRUE) == 0L,
          paste(c(sum(all27_prog_int$q_value_BH <= 0.05, na.rm = TRUE),
                  sum(all27_doro_int$q_value_BH <= 0.05, na.rm = TRUE),
                  sum(orig_prog_int$q_value_BH <= 0.05, na.rm = TRUE),
                  sum(orig_doro_int$q_value_BH <= 0.05, na.rm = TRUE)), collapse = ";"),
          "0;0;0;0")
add_check("LOO_Reactome_complete",
          nrow(loo_reactome) == 1656L * 26L &&
            uniqueN(loo_reactome$feature_id) == 1656L &&
            uniqueN(loo_reactome$omitted_sample_id) == 26L,
          paste(nrow(loo_reactome), uniqueN(loo_reactome$feature_id),
                uniqueN(loo_reactome$omitted_sample_id), sep = ";"),
          "43056;1656;26")
add_check("LOO_PROGENy_complete",
          nrow(loo_progeny) == 14L * 26L &&
            uniqueN(loo_progeny$feature_id) == 14L &&
            uniqueN(loo_progeny$omitted_sample_id) == 26L,
          paste(nrow(loo_progeny), uniqueN(loo_progeny$feature_id),
                uniqueN(loo_progeny$omitted_sample_id), sep = ";"),
          "364;14;26")
add_check("LOO_DoRothEA_complete",
          nrow(loo_dorothea) == 101L * 26L &&
            uniqueN(loo_dorothea$feature_id) == 101L &&
            uniqueN(loo_dorothea$omitted_sample_id) == 26L,
          paste(nrow(loo_dorothea), uniqueN(loo_dorothea$feature_id),
                uniqueN(loo_dorothea$omitted_sample_id), sep = ";"),
          "2626;101;26")

add_check("precision_families_complete", nrow(precision) == 4L &&
            setequal(precision$family, c("representative_gene", "Reactome", "PROGENy", "DoRothEA")),
          paste(precision$family, collapse = ";"),
          "representative_gene;Reactome;PROGENy;DoRothEA")
add_check("precision_CI_widths_finite_positive",
          all(is.finite(precision$ci_width_median) & precision$ci_width_median > 0),
          paste(signif(precision$ci_width_median, 6), collapse = ";"), "all finite and >0")
add_check("precision_not_interpreted_as_equivalence",
          all(precision$precision_interpretation ==
                "descriptive_precision_only_not_equivalence_or_evidence_of_no_effect"),
          paste(unique(precision$precision_interpretation), collapse = ";"),
          "descriptive precision only; not equivalence")

quickgo_tables <- list(quickgo_camera, quickgo_mroast, quickgo_fgsea, quickgo_activity)
quickgo_gate_flags_ok <- all(vapply(quickgo_tables, function(x) {
  all(x$sensitivity_only) && all(!x$candidate_gate_eligible) && all(!x$S1_entry_permitted)
}, logical(1)))
add_check("QuickGO_SynGO_proxy_completed", nrow(quickgo_status) == 1L &&
            quickgo_status$analysis_status == "COMPLETED_SENSITIVITY_ONLY" &&
            quickgo_status$n_eligible_sets == 78L,
          paste(quickgo_status$analysis_status, quickgo_status$n_eligible_sets, sep = ";"),
          "COMPLETED_SENSITIVITY_ONLY;78")
add_check("QuickGO_SynGO_proxy_not_portal_release",
          identical(quickgo_status$equivalent_to_SynGO_portal_1_3_bulk, FALSE),
          quickgo_status$equivalent_to_SynGO_portal_1_3_bulk, FALSE)
add_check("QuickGO_SynGO_all_outputs_sensitivity_only", quickgo_gate_flags_ok,
          quickgo_gate_flags_ok, TRUE)
add_check("QuickGO_SynGO_camera_hits_zero", sum(quickgo_camera$q_value_BH <= 0.05) == 0L,
          sum(quickgo_camera$q_value_BH <= 0.05), 0L)
add_check("QuickGO_SynGO_mroast_directional_hits_zero", sum(quickgo_mroast$FDR <= 0.05) == 0L,
          sum(quickgo_mroast$FDR <= 0.05), 0L)
add_check("QuickGO_SynGO_mroast_mixed_hits_zero", sum(quickgo_mroast$FDR.Mixed <= 0.05) == 0L,
          sum(quickgo_mroast$FDR.Mixed <= 0.05), 0L)
add_check("QuickGO_SynGO_activity_hits_zero", sum(quickgo_activity$q_value_BH <= 0.05) == 0L,
          sum(quickgo_activity$q_value_BH <= 0.05), 0L)
add_check("QuickGO_SynGO_fgsea_sensitivity_hits_recorded", sum(quickgo_fgsea$padj <= 0.05) == 14L,
          sum(quickgo_fgsea$padj <= 0.05), 14L)
add_check("QuickGO_SynGO_absent_from_S1_GATE",
          !any(s1$collection %chin% c("QuickGO", "QuickGO_SynGO_sensitivity",
                                     "QuickGO_assignedBy_SynGO")),
          paste(sort(unique(s1$collection)), collapse = ";"),
          "Reactome;PROGENy;DoRothEA only")

checks_dt <- rbindlist(checks)
fwrite(checks_dt, file.path(diagnostic_dir, "GSE16710_output_validation.tsv"), sep = "\t", quote = FALSE)

summary <- list(
  checks_total = nrow(checks_dt),
  checks_pass = sum(checks_dt$status == "PASS"),
  checks_fail = sum(checks_dt$status == "FAIL"),
  overall_status = if (all(checks_dt$status == "PASS")) "PASS" else "FAIL",
  validated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
)
write_json(summary, file.path(diagnostic_dir, "GSE16710_output_validation_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(checks_dt)
if (any(checks_dt$status == "FAIL")) stop("One or more GSE16710 validation checks failed")
