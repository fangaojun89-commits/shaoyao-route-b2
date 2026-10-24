#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
result_dir <- file.path(project_root, "results", "GSE16710")

gene <- fread(file.path(result_dir, "gene_level", "limma_gene_interaction_primary.tsv.gz"))
react <- fread(file.path(result_dir, "program_level", "reactome_mean_z_activity_all_contrasts.tsv"))[
  contrast == "spasticity_progression_interaction"]
prog <- fread(file.path(result_dir, "program_level", "progeny_ulm_activity_all_contrasts.tsv"))[
  contrast == "spasticity_progression_interaction"]
doro <- fread(file.path(result_dir, "program_level", "dorothea_ulm_activity_all_contrasts.tsv"))[
  contrast == "spasticity_progression_interaction"]

summarise_ci <- function(x, family, scale) {
  width <- x$ci_high - x$ci_low
  data.table(
    family = family,
    effect_scale = scale,
    n_effects = sum(is.finite(width)),
    ci_width_min = min(width, na.rm = TRUE),
    ci_width_q25 = quantile(width, 0.25, na.rm = TRUE, names = FALSE),
    ci_width_median = median(width, na.rm = TRUE),
    ci_width_q75 = quantile(width, 0.75, na.rm = TRUE, names = FALSE),
    ci_width_max = max(width, na.rm = TRUE),
    median_half_width = median(width, na.rm = TRUE) / 2,
    precision_interpretation = "descriptive_precision_only_not_equivalence_or_evidence_of_no_effect"
  )
}

out <- rbindlist(list(
  summarise_ci(gene, "representative_gene", "interaction_log2_expression"),
  summarise_ci(react, "Reactome", "mean_gene_z_score_program_activity"),
  summarise_ci(prog, "PROGENy", "decoupleR_ULM_activity_score"),
  summarise_ci(doro, "DoRothEA", "decoupleR_ULM_activity_score")
))

fwrite(out, file.path(result_dir, "diagnostics", "interaction_precision_summary.tsv"),
       sep = "\t", quote = FALSE, na = "NA")
print(out)
