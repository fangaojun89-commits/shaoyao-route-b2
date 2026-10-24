#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(matrixStats)
  library(fgsea)
  library(BiocParallel)
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
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
gene_set_dir <- file.path(project_root, "data", "processed", "gene_sets")
program_dir <- file.path(project_root, "results", "GSE16710", "program_level")
diagnostic_dir <- file.path(project_root, "results", "GSE16710", "diagnostics")

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")

make_design <- function(meta) {
  meta <- copy(meta)
  meta[, injury := factor(experimental_arm, levels = c("Sham", "Spi"))]
  meta[, time_factor := factor(as.character(time_day), levels = c("21", "60"))]
  x <- model.matrix(~ injury * time_factor, meta)
  rownames(x) <- meta$sample_id
  if (qr(x)$rank != 4L) stop("QuickGO-SynGO sensitivity design is not full rank")
  x
}

zscore_rows <- function(x) {
  s <- matrixStats::rowSds(x)
  keep <- is.finite(s) & s > 0
  z <- sweep(x[keep, , drop = FALSE], 1L, rowMeans(x[keep, , drop = FALSE]), "-")
  sweep(z, 1L, s[keep], "/")
}

activity_stats <- function(activity, design) {
  fit <- lmFit(activity, design)
  fit <- contrasts.fit(fit, c(0, 0, 0, 1))
  fit <- eBayes(fit, robust = TRUE)
  se <- fit$stdev.unscaled[, 1L] * sqrt(fit$s2.post)
  crit <- qt(0.975, df = rep_len(fit$df.total, nrow(activity)))
  data.table(
    feature_id = rownames(activity),
    estimate = fit$coefficients[, 1L],
    standard_error = se,
    ci_low = fit$coefficients[, 1L] - crit * se,
    ci_high = fit$coefficients[, 1L] + crit * se,
    moderated_t = fit$t[, 1L],
    p_value = fit$p.value[, 1L],
    q_value_BH = p.adjust(fit$p.value[, 1L], "BH")
  )
}

probe_all <- readRDS(file.path(processed_dir, "rma_probe_expression.rds"))
gene_all <- readRDS(file.path(processed_dir, "rma_gene_expression_current_bioconductor.rds"))
ledger <- fread(file.path(processed_dir, "sample_ledger_verified.tsv"))
decision <- fread(file.path(adjudication_dir, "qc_inclusion_decision.tsv"))
probe_map <- fread(file.path(processed_dir, "probe_annotation_current_bioconductor.tsv"))
frozen <- readRDS(file.path(gene_set_dir, "gene_set_resources_frozen_2026-08-18.rds"))

ledger <- merge(ledger, decision, by = "sample_id", all.x = TRUE, sort = FALSE)
ledger[, order__ := as.integer(sample_order)]
setorder(ledger, order__)
primary_samples <- ledger[primary_analysis_decision == "retain", sample_id]
meta <- ledger[match(primary_samples, sample_id)]
design <- make_design(meta)
probe <- probe_all[, primary_samples, drop = FALSE]
gene <- gene_all[, primary_samples, drop = FALSE]

ortho <- frozen$ortholog_one_to_one[other_tax_id == 10116,
                                    .(human_symbol, rat_symbol = other_symbol)]
syngo <- merge(frozen$syngo_sensitivity[, .(term_id, term_name, human_symbol = gene_symbol)],
               ortho, by = "human_symbol", allow.cartesian = TRUE)
syngo <- unique(syngo)
coverage <- syngo[, .(
  source_human_genes = uniqueN(human_symbol),
  one_to_one_rat_genes = uniqueN(rat_symbol),
  mapped_expression_genes = uniqueN(rat_symbol[rat_symbol %in% rownames(gene)])
), by = .(term_id, term_name)]
coverage[, `:=`(
  eligible_10_to_500 = mapped_expression_genes >= 10L & mapped_expression_genes <= 500L,
  source_identity = "EMBL_EBI_QuickGO_assignedBy_SynGO_snapshot_2026-08-17",
  equivalent_to_SynGO_portal_1_3_bulk = FALSE,
  candidate_gate_eligible = FALSE,
  S1_entry_permitted = FALSE
)]
write_tsv(coverage, file.path(program_dir, "quickgo_assignedBy_SynGO_mapping_coverage.tsv"))
write_tsv(syngo, file.path(program_dir, "quickgo_assignedBy_SynGO_human_to_rat_mapping.tsv.gz"))

eligible <- coverage[eligible_10_to_500 == TRUE, term_id]
if (!length(eligible)) {
  status <- data.table(
    source_identity = "EMBL_EBI_QuickGO_assignedBy_SynGO_snapshot_2026-08-17",
    analysis_status = "NOT_RUN_NO_ELIGIBLE_SETS_AFTER_MAPPING",
    candidate_gate_eligible = FALSE,
    S1_entry_permitted = FALSE
  )
  write_tsv(status, file.path(diagnostic_dir, "quickgo_assignedBy_SynGO_sensitivity_status.tsv"))
  quit(save = "no", status = 0L)
}

pathways <- lapply(eligible, function(id) sort(unique(syngo[term_id == id & rat_symbol %in% rownames(gene), rat_symbol])))
names(pathways) <- eligible
probe_map[, probe_index := match(probe_id, rownames(probe))]
index <- lapply(pathways, function(g) sort(unique(probe_map[gene_symbol %in% g & !is.na(probe_index), probe_index])))

camera <- limma::camera(probe, index = index, design = design, contrast = c(0, 0, 0, 1),
                        inter.gene.cor = 0.01, use.ranks = FALSE, allow.neg.cor = FALSE,
                        sort = FALSE, directional = TRUE)
camera <- as.data.table(camera, keep.rownames = "feature_id")
setnames(camera, c("NGenes", "Direction", "PValue", "FDR"),
         c("n_probe_sets", "direction", "p_value", "camera_reported_FDR"))
camera[, `:=`(
  q_value_BH = p.adjust(p_value, "BH"),
  sensitivity_only = TRUE,
  candidate_gate_eligible = FALSE,
  S1_entry_permitted = FALSE
)]
camera <- merge(coverage[, .(feature_id = term_id, term_name)], camera, by = "feature_id", all.y = TRUE)
write_tsv(camera, file.path(program_dir, "quickgo_assignedBy_SynGO_camera_sensitivity.tsv"))

set.seed(20260817)
roast <- limma::mroast(probe, index = index, design = design, contrast = c(0, 0, 0, 1),
                       set.statistic = "mean", nrot = 9999, approx.zscore = TRUE,
                       adjust.method = "BH", midp = TRUE, sort = "none")
roast <- as.data.table(roast, keep.rownames = "feature_id")
roast[, `:=`(sensitivity_only = TRUE, candidate_gate_eligible = FALSE, S1_entry_permitted = FALSE)]
write_tsv(roast, file.path(program_dir, "quickgo_assignedBy_SynGO_mroast_sensitivity.tsv"))

z <- zscore_rows(gene)
pathways_z <- lapply(pathways, intersect, y = rownames(z))
activity <- t(vapply(pathways_z, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(z))))
colnames(activity) <- colnames(z)
act <- activity_stats(activity, design)
act[, `:=`(sensitivity_only = TRUE, candidate_gate_eligible = FALSE, S1_entry_permitted = FALSE)]
act <- merge(coverage[, .(feature_id = term_id, term_name)], act, by = "feature_id", all.y = TRUE)
write_tsv(act, file.path(program_dir, "quickgo_assignedBy_SynGO_activity_sensitivity.tsv"))

gene_stats <- fread(file.path(project_root, "results", "GSE16710", "gene_level", "limma_gene_interaction_primary.tsv.gz"))
ranks <- setNames(gene_stats$moderated_t, gene_stats$gene_symbol)
ranks <- sort(ranks[is.finite(ranks) & !is.na(names(ranks))], decreasing = TRUE)
set.seed(20260817)
fg <- fgsea::fgseaMultilevel(pathways, ranks, minSize = 10, maxSize = 500, eps = 0,
                             scoreType = "std", nPermSimple = 10000,
                             BPPARAM = BiocParallel::SerialParam())
fg <- as.data.table(fg)
if ("leadingEdge" %in% names(fg)) fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
setnames(fg, "pathway", "feature_id")
fg[, `:=`(sensitivity_only = TRUE, candidate_gate_eligible = FALSE, S1_entry_permitted = FALSE)]
write_tsv(fg, file.path(program_dir, "quickgo_assignedBy_SynGO_fgsea_sensitivity.tsv"))

status <- data.table(
  source_identity = "EMBL_EBI_QuickGO_assignedBy_SynGO_snapshot_2026-08-17",
  analysis_status = "COMPLETED_SENSITIVITY_ONLY",
  n_eligible_sets = length(eligible),
  camera_BH_hits = sum(camera$q_value_BH <= 0.05, na.rm = TRUE),
  mroast_directional_BH_hits = sum(roast$FDR <= 0.05, na.rm = TRUE),
  fgsea_BH_hits = sum(fg$padj <= 0.05, na.rm = TRUE),
  activity_BH_hits = sum(act$q_value_BH <= 0.05, na.rm = TRUE),
  equivalent_to_SynGO_portal_1_3_bulk = FALSE,
  candidate_gate_eligible = FALSE,
  S1_entry_permitted = FALSE,
  interpretation = "proxy_sensitivity_only_never_enters_S1"
)
write_tsv(status, file.path(diagnostic_dir, "quickgo_assignedBy_SynGO_sensitivity_status.tsv"))
print(status)
