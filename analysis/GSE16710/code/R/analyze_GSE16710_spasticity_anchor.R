#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 140)

required <- c(
  "limma", "data.table", "matrixStats", "fgsea", "decoupleR", "BiocParallel",
  "ggplot2", "yaml", "jsonlite"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(matrixStats)
  library(fgsea)
  library(decoupleR)
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
if (length(script_arg) != 1L) stop("Cannot determine script path")
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
result_dir <- file.path(project_root, "results", "GSE16710")
gene_result_dir <- file.path(result_dir, "gene_level")
program_result_dir <- file.path(result_dir, "program_level")
diagnostic_dir <- file.path(result_dir, "diagnostics")
figure_dir <- file.path(result_dir, "figures")
log_dir <- file.path(result_dir, "logs")
for (d in c(gene_result_dir, program_result_dir, diagnostic_dir, figure_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "analyze_GSE16710_spasticity_anchor.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

set.seed(20260817)
cat("GSE16710 phenotype-anchor modeling started:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), file = path, sep = "\t", quote = FALSE, na = "NA")
}

safe_sign <- function(x, tol = 1e-12) {
  out <- rep(NA_integer_, length(x))
  out[is.finite(x) & x > tol] <- 1L
  out[is.finite(x) & x < -tol] <- -1L
  out[is.finite(x) & abs(x) <= tol] <- 0L
  out
}

zscore_rows <- function(x) {
  center <- rowMeans(x)
  spread <- matrixStats::rowSds(x)
  keep <- is.finite(spread) & spread > 0 & is.finite(center)
  z <- sweep(x[keep, , drop = FALSE], 1L, center[keep], "-")
  z <- sweep(z, 1L, spread[keep], "/")
  z
}

make_design <- function(meta) {
  meta <- copy(meta)
  meta[, injury := factor(experimental_arm, levels = c("Sham", "Spi"))]
  meta[, time_factor := factor(as.character(time_day), levels = c("21", "60"))]
  design <- stats::model.matrix(~ injury * time_factor, data = meta)
  expected <- c("(Intercept)", "injurySpi", "time_factor60", "injurySpi:time_factor60")
  if (!identical(colnames(design), expected)) stop("Unexpected GSE16710 design columns")
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient GSE16710 design")
  rownames(design) <- meta$sample_id
  contrast_matrix <- cbind(
    spasticity_progression_interaction = c(0, 0, 0, 1),
    injury_at_day21 = c(0, 1, 0, 0),
    injury_at_day60 = c(0, 1, 0, 1),
    progression_within_spi = c(0, 0, 1, 1),
    time_within_sham = c(0, 0, 1, 0)
  )
  rownames(contrast_matrix) <- colnames(design)
  list(meta = meta, design = design, contrasts = contrast_matrix)
}

fit_limma_matrix <- function(mat, design_obj, robust = TRUE) {
  if (!identical(colnames(mat), rownames(design_obj$design))) stop("Expression/activity columns do not match design rows")
  fit <- limma::lmFit(mat, design_obj$design)
  fit <- limma::contrasts.fit(fit, design_obj$contrasts)
  fit <- limma::eBayes(fit, robust = robust)
  fit
}

extract_fit_long <- function(fit, feature_col = "feature_id") {
  out <- vector("list", ncol(fit$coefficients))
  for (j in seq_len(ncol(fit$coefficients))) {
    se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
    df <- rep_len(fit$df.total, nrow(fit$coefficients))
    crit <- stats::qt(0.975, df = df)
    p <- fit$p.value[, j]
    out[[j]] <- data.table(
      feature_id = rownames(fit$coefficients),
      contrast = colnames(fit$coefficients)[j],
      estimate = fit$coefficients[, j],
      standard_error = se,
      ci_low = fit$coefficients[, j] - crit * se,
      ci_high = fit$coefficients[, j] + crit * se,
      moderated_t = fit$t[, j],
      df_total = df,
      p_value = p,
      q_value_BH = stats::p.adjust(p, method = "BH"),
      B_statistic = fit$lods[, j]
    )
  }
  ans <- rbindlist(out)
  setnames(ans, "feature_id", feature_col)
  ans[]
}

collapse_probe_annotation <- function(map_dt) {
  map_dt[, .(
    gene_symbols = paste(sort(unique(gene_symbol)), collapse = ";"),
    entrez_ids = if ("entrez_id" %in% names(map_dt)) paste(sort(unique(entrez_id)), collapse = ";") else NA_character_,
    n_mapped_genes = uniqueN(gene_symbol)
  ), by = probe_id]
}

build_reactome_mapping <- function(frozen, available_genes, probe_map = NULL, fixed_terms = NULL) {
  ortho_rat <- frozen$ortholog_one_to_one[other_tax_id == 10116,
                                          .(human_symbol, rat_symbol = other_symbol, human_entrez, rat_entrez = other_entrez)]
  react <- merge(frozen$reactome[, .(term_id, term_name, human_symbol = gene_symbol)], ortho_rat,
                 by = "human_symbol", allow.cartesian = TRUE)
  react <- unique(react[, .(term_id, term_name, human_symbol, rat_symbol, human_entrez, rat_entrez)])
  coverage <- react[, .(
    source_human_genes = uniqueN(human_symbol),
    one_to_one_rat_genes = uniqueN(rat_symbol),
    mapped_expression_genes = uniqueN(rat_symbol[rat_symbol %in% available_genes])
  ), by = .(term_id, term_name)]
  coverage[, mapping_fraction_of_one_to_one := mapped_expression_genes / one_to_one_rat_genes]
  if (is.null(fixed_terms)) {
    coverage[, eligible := mapped_expression_genes >= 10L & mapped_expression_genes <= 500L]
  } else {
    coverage[, eligible := term_id %in% fixed_terms & mapped_expression_genes >= 10L]
  }
  term_ids <- if (is.null(fixed_terms)) coverage[eligible == TRUE, term_id] else intersect(fixed_terms, coverage[eligible == TRUE, term_id])
  pathways <- lapply(term_ids, function(id) sort(unique(react[term_id == id & rat_symbol %in% available_genes, rat_symbol])))
  names(pathways) <- term_ids
  if (!is.null(probe_map)) {
    probe_index <- lapply(pathways, function(g) sort(unique(probe_map[gene_symbol %in% g, probe_index])))
    probe_index <- probe_index[lengths(probe_index) > 0L]
    coverage[, mapped_probe_sets := lengths(probe_index)[match(term_id, names(probe_index))]]
  } else {
    probe_index <- NULL
  }
  list(long = react, coverage = coverage, pathways = pathways, probe_index = probe_index)
}

program_mean_activity <- function(gene_expr, pathways) {
  z <- zscore_rows(gene_expr)
  pathways <- lapply(pathways, function(g) intersect(g, rownames(z)))
  pathways <- pathways[lengths(pathways) >= 10L]
  score <- t(vapply(pathways, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(z))))
  colnames(score) <- colnames(z)
  score
}

run_camera_contrasts <- function(expr, index, design_obj, collection_label, contrasts = colnames(design_obj$contrasts)) {
  index <- index[lengths(index) > 0L]
  out <- vector("list", length(contrasts))
  for (k in seq_along(contrasts)) {
    nm <- contrasts[k]
    cam <- limma::camera(expr, index = index, design = design_obj$design,
                         contrast = design_obj$contrasts[, nm], inter.gene.cor = 0.01,
                         use.ranks = FALSE, allow.neg.cor = FALSE, sort = FALSE, directional = TRUE)
    dt <- as.data.table(cam, keep.rownames = "feature_id")
    rename_map <- c(NGenes = "n_genes", Direction = "direction", PValue = "p_value", FDR = "camera_reported_FDR")
    for (old in intersect(names(rename_map), names(dt))) setnames(dt, old, rename_map[[old]])
    if (!"p_value" %in% names(dt)) stop("camera result lacks PValue")
    dt[, `:=`(contrast = nm, collection = collection_label, q_value_BH = p.adjust(p_value, method = "BH"))]
    out[[k]] <- dt
  }
  rbindlist(out, fill = TRUE)
}

run_mroast_primary <- function(expr, index, design_obj) {
  set.seed(20260817)
  roast <- limma::mroast(expr, index = index, design = design_obj$design,
                         contrast = design_obj$contrasts[, "spasticity_progression_interaction"],
                         set.statistic = "mean", nrot = 9999, approx.zscore = TRUE,
                         adjust.method = "BH", midp = TRUE, sort = "none")
  as.data.table(roast, keep.rownames = "feature_id")
}

run_fgsea_primary <- function(gene_fit_long, pathways) {
  ranks <- gene_fit_long[contrast == "spasticity_progression_interaction",
                         setNames(moderated_t, gene_symbol)]
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  set.seed(20260817)
  fg <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500,
                               eps = 0, scoreType = "std", nPermSimple = 10000,
                               BPPARAM = BiocParallel::SerialParam())
  fg <- as.data.table(fg)
  if ("leadingEdge" %in% names(fg)) fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
  setnames(fg, "pathway", "feature_id")
  fg[]
}

build_ulm_networks <- function(frozen, available_genes) {
  ortho_rat <- frozen$ortholog_one_to_one[other_tax_id == 10116, .(human_symbol, rat_symbol = other_symbol)]

  prog <- merge(frozen$progeny_human[, .(source = as.character(pathway), human_symbol = as.character(gene_symbol), mor = weight)],
                ortho_rat, by = "human_symbol", allow.cartesian = TRUE)
  prog <- prog[rat_symbol %in% available_genes,
               .(mor = mean(mor), n_source_rows = .N), by = .(source, target = rat_symbol)]
  prog_count <- prog[, .(mapped_targets = uniqueN(target)), by = source]
  prog_keep <- prog_count[mapped_targets >= 10L & mapped_targets <= 500L, source]
  prog <- prog[source %in% prog_keep, .(source, target, mor)]

  doro <- frozen$dorothea_rat[, .(source = as.character(tf), target = as.character(target), mor = as.numeric(mor))]
  doro <- doro[target %in% available_genes, .(mor = mean(mor)), by = .(source, target)]
  doro_count <- doro[, .(mapped_targets = uniqueN(target)), by = source]
  doro_keep <- doro_count[mapped_targets >= 10L & mapped_targets <= 500L, source]
  doro <- doro[source %in% doro_keep]

  list(progeny = prog, progeny_coverage = prog_count,
       dorothea = doro, dorothea_coverage = doro_count)
}

run_ulm_activity <- function(gene_expr, network) {
  z <- zscore_rows(gene_expr)
  network <- network[target %in% rownames(z)]
  counts <- network[, .(mapped_targets = uniqueN(target)), by = source]
  keep <- counts[mapped_targets >= 10L & mapped_targets <= 500L, source]
  network <- network[source %in% keep]
  if (!nrow(network)) return(matrix(numeric(), nrow = 0L, ncol = ncol(z)))
  res <- decoupleR::run_ulm(z, network, .source = source, .target = target, .mor = mor,
                            sparse = FALSE, center = FALSE, na.rm = FALSE, minsize = 10L)
  res <- as.data.table(res)
  wide <- dcast(res, source ~ condition, value.var = "score")
  score <- as.matrix(wide[, -"source"])
  rownames(score) <- wide$source
  score <- score[, colnames(z), drop = FALSE]
  score
}

leave_one_out_effects <- function(activity, meta) {
  full_design <- make_design(meta)
  full_fit <- fit_limma_matrix(activity, full_design, robust = TRUE)
  full_est <- full_fit$coefficients[, "spasticity_progression_interaction"]
  out <- vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    keep <- seq_len(nrow(meta)) != i
    d <- make_design(meta[keep])
    f <- fit_limma_matrix(activity[, keep, drop = FALSE], d, robust = TRUE)
    out[[i]] <- data.table(
      feature_id = rownames(activity),
      omitted_sample_id = meta$sample_id[i],
      loo_estimate = f$coefficients[, "spasticity_progression_interaction"],
      loo_moderated_t = f$t[, "spasticity_progression_interaction"],
      loo_p_value = f$p.value[, "spasticity_progression_interaction"]
    )
  }
  loo <- rbindlist(out)
  summary <- loo[, .(
    full_estimate = full_est[match(feature_id[1L], names(full_est))],
    loo_min = min(loo_estimate, na.rm = TRUE),
    loo_max = max(loo_estimate, na.rm = TRUE),
    loo_direction_fraction = mean(safe_sign(loo_estimate) == safe_sign(full_est[match(feature_id[1L], names(full_est))]), na.rm = TRUE),
    loo_sign_flips = sum(safe_sign(loo_estimate) != safe_sign(full_est[match(feature_id[1L], names(full_est))]), na.rm = TRUE),
    max_relative_change = max(abs(loo_estimate - full_est[match(feature_id[1L], names(full_est))]) /
                                (abs(full_est[match(feature_id[1L], names(full_est))]) + 1e-12), na.rm = TRUE)
  ), by = feature_id]
  summary[, `:=`(
    passes_direction_fraction_0_80 = loo_direction_fraction >= 0.80,
    no_single_animal_sign_flip = loo_sign_flips == 0L
  )]
  list(long = loo, summary = summary)
}

ordinary_t <- function(y, design, contrast) {
  xtx_inv <- solve(crossprod(design))
  beta <- y %*% t(xtx_inv %*% t(design))
  fitted <- beta %*% t(design)
  residual <- y - fitted
  df <- ncol(y) - ncol(design)
  s2 <- rowSums(residual^2) / df
  estimate <- as.vector(beta %*% contrast)
  se <- sqrt(s2 * as.numeric(t(contrast) %*% xtx_inv %*% contrast))
  list(t = estimate / se, estimate = estimate, se = se, df = df)
}

freedman_lane <- function(activity, meta, nperm = 10000L) {
  meta2 <- copy(meta)
  meta2[, injury := factor(experimental_arm, levels = c("Sham", "Spi"))]
  meta2[, time_factor := factor(as.character(time_day), levels = c("21", "60"))]
  x0 <- model.matrix(~ injury + time_factor, data = meta2)
  x1 <- model.matrix(~ injury * time_factor, data = meta2)
  contrast <- c(0, 0, 0, 1)
  rownames(x0) <- rownames(x1) <- meta2$sample_id
  fit0 <- activity %*% t(solve(crossprod(x0)) %*% t(x0))
  fitted0 <- fit0 %*% t(x0)
  residual0 <- activity - fitted0
  observed <- ordinary_t(activity, x1, contrast)
  exceed <- integer(nrow(activity))
  valid <- integer(nrow(activity))
  null_tail_005 <- 0
  null_tail_001 <- 0
  null_total <- 0
  strata <- split(seq_len(nrow(meta2)), meta2$time_factor)
  set.seed(20260817)
  for (b in seq_len(nperm)) {
    perm <- seq_len(nrow(meta2))
    for (idx in strata) perm[idx] <- sample(idx, length(idx), replace = FALSE)
    yp <- fitted0 + residual0[, perm, drop = FALSE]
    tp <- ordinary_t(yp, x1, contrast)
    ok <- is.finite(tp$t) & is.finite(observed$t)
    exceed[ok] <- exceed[ok] + as.integer(abs(tp$t[ok]) >= abs(observed$t[ok]))
    valid[ok] <- valid[ok] + 1L
    p_theory <- 2 * stats::pt(-abs(tp$t[ok]), df = tp$df)
    null_tail_005 <- null_tail_005 + sum(p_theory <= 0.05)
    null_tail_001 <- null_tail_001 + sum(p_theory <= 0.01)
    null_total <- null_total + length(p_theory)
    if (b %% 1000L == 0L) cat("Freedman-Lane permutations:", b, "/", nperm, "\n")
  }
  result <- data.table(
    feature_id = rownames(activity),
    ordinary_t_observed = observed$t,
    empirical_p_value = (1 + exceed) / (1 + valid),
    valid_permutations = valid
  )
  calibration <- data.table(
    n_program_permutation_tests = null_total,
    observed_null_tail_p_le_0_05 = null_tail_005 / null_total,
    observed_null_tail_p_le_0_01 = null_tail_001 / null_total,
    expected_tail_0_05 = 0.05,
    expected_tail_0_01 = 0.01,
    anti_conservative_flag = (null_tail_005 / null_total) > 0.075 | (null_tail_001 / null_total) > 0.02
  )
  list(result = result, calibration = calibration)
}

# ----------------------------------------------------------------------------
# Frozen inputs, QC inclusion decision and model estimability
# ----------------------------------------------------------------------------

probe_expr_all <- readRDS(file.path(processed_dir, "rma_probe_expression.rds"))
gene_expr_current_all <- readRDS(file.path(processed_dir, "rma_gene_expression_current_bioconductor.rds"))
gene_expr_original_all <- readRDS(file.path(processed_dir, "rma_gene_expression_original_GPL1355_2016.rds"))
ledger <- fread(file.path(processed_dir, "sample_ledger_verified.tsv"), sep = "\t")
decision <- fread(file.path(adjudication_dir, "qc_inclusion_decision.tsv"), sep = "\t")
current_map <- fread(file.path(processed_dir, "probe_annotation_current_bioconductor.tsv"), sep = "\t")
original_map <- fread(file.path(processed_dir, "probe_annotation_original_GPL1355_2016.tsv"), sep = "\t")
frozen <- readRDS(file.path(gene_set_dir, "gene_set_resources_frozen_2026-08-18.rds"))

expected_decision_columns <- c(
  "sample_id", "primary_analysis_decision", "decision_stage", "technical_basis",
  "outcome_information_used", "adjudicated_on"
)
if (!identical(names(decision), expected_decision_columns)) stop("Unexpected frozen QC adjudication schema")
if (nrow(decision) != 27L || anyDuplicated(decision$sample_id)) stop("Frozen QC adjudication is not one row per array")
if (!setequal(ledger$sample_id, decision$sample_id)) stop("QC inclusion decision does not cover the frozen ledger exactly")
if (!all(decision$decision_stage == "outcome_blind_qc")) stop("Unexpected QC adjudication stage")
if (!all(toupper(as.character(decision$outcome_information_used)) == "FALSE")) {
  stop("Frozen QC adjudication records outcome information use")
}
if (!setequal(unique(decision$primary_analysis_decision),
              c("retain", "exclude_primary_include_all_sample_sensitivity"))) {
  stop("Unexpected frozen QC adjudication decision value")
}
ledger <- merge(ledger, decision, by = "sample_id", all.x = TRUE, sort = FALSE)
ledger[, sample_order_int__ := as.integer(sample_order)]
setorder(ledger, sample_order_int__)
ledger[, sample_order_int__ := NULL]
primary_samples <- ledger[primary_analysis_decision == "retain", sample_id]
excluded_samples <- ledger[primary_analysis_decision != "retain", sample_id]
if (!identical(excluded_samples, "GSM418682")) stop("Unexpected primary QC exclusion set")
if (!identical(colnames(probe_expr_all), ledger$sample_id)) stop("Probe expression columns differ from ledger")

meta_primary <- ledger[match(primary_samples, sample_id)]
meta_all <- ledger[match(colnames(probe_expr_all), sample_id)]
design_primary <- make_design(meta_primary)
design_all <- make_design(meta_all)

model_diag <- data.table(
  analysis = c("primary_outcome_blind_QC_set", "all_27_array_sensitivity"),
  n = c(nrow(meta_primary), nrow(meta_all)),
  design_columns = c(ncol(design_primary$design), ncol(design_all$design)),
  design_rank = c(qr(design_primary$design)$rank, qr(design_all$design)$rank),
  residual_df = c(nrow(meta_primary) - ncol(design_primary$design), nrow(meta_all) - ncol(design_all$design)),
  model_estimable = c(qr(design_primary$design)$rank == ncol(design_primary$design),
                      qr(design_all$design)$rank == ncol(design_all$design))
)
write_tsv(model_diag, file.path(diagnostic_dir, "model_estimability.tsv"))
write_tsv(as.data.table(design_primary$design, keep.rownames = "sample_id"), file.path(diagnostic_dir, "design_matrix_primary.tsv"))
write_tsv(as.data.table(design_primary$contrasts, keep.rownames = "coefficient"), file.path(diagnostic_dir, "contrast_matrix.tsv"))
write_tsv(ledger, file.path(diagnostic_dir, "analysis_sample_inclusion.tsv"))

# ----------------------------------------------------------------------------
# Probe- and gene-level limma models
# ----------------------------------------------------------------------------

cat("Fitting primary probe-level limma model on", length(primary_samples), "animals...\n")
probe_expr <- probe_expr_all[, primary_samples, drop = FALSE]
probe_fit <- fit_limma_matrix(probe_expr, design_primary, robust = TRUE)
probe_stats <- extract_fit_long(probe_fit, "probe_id")
probe_ann_collapsed <- collapse_probe_annotation(current_map)
probe_stats <- merge(probe_stats, probe_ann_collapsed, by = "probe_id", all.x = TRUE, sort = FALSE)
write_tsv(probe_stats, file.path(gene_result_dir, "limma_probe_statistics_all_contrasts.tsv.gz"))

rep_current <- fread(file.path(processed_dir, "gene_representative_probe_current_bioconductor.tsv"), sep = "\t")
gene_stats <- merge(rep_current[, .(gene_symbol, probe_id, entrez_id, probe_mean_expression)], probe_stats,
                    by = "probe_id", all.x = TRUE, allow.cartesian = TRUE)
write_tsv(gene_stats, file.path(gene_result_dir, "limma_gene_statistics_representative_probe_all_contrasts.tsv.gz"))

cat("Fitting all-27-array sensitivity model...\n")
probe_fit_all <- fit_limma_matrix(probe_expr_all, design_all, robust = TRUE)
probe_stats_all <- extract_fit_long(probe_fit_all, "probe_id")
probe_stats_all_primary <- probe_stats_all[contrast == "spasticity_progression_interaction"]
write_tsv(probe_stats_all_primary, file.path(gene_result_dir, "limma_probe_interaction_all27_sensitivity.tsv.gz"))

primary_gene_stats <- gene_stats[contrast == "spasticity_progression_interaction"]
write_tsv(primary_gene_stats, file.path(gene_result_dir, "limma_gene_interaction_primary.tsv.gz"))

# ----------------------------------------------------------------------------
# Reactome: camera primary, mroast/fgsea sensitivity, activity CIs
# ----------------------------------------------------------------------------

current_map[, probe_index := match(probe_id, rownames(probe_expr))]
current_map <- current_map[!is.na(probe_index)]
original_map[, probe_index := match(probe_id, rownames(probe_expr))]
original_map <- original_map[!is.na(probe_index)]

react_current <- build_reactome_mapping(frozen, rownames(gene_expr_current_all), current_map)
eligible_reactome <- react_current$coverage[eligible == TRUE, term_id]
if (!length(eligible_reactome)) stop("No eligible Reactome programs after GSE16710 mapping")
react_current$pathways <- react_current$pathways[eligible_reactome]
react_current$probe_index <- react_current$probe_index[eligible_reactome]
write_tsv(react_current$coverage, file.path(program_result_dir, "reactome_mapping_coverage_current.tsv"))
write_tsv(react_current$long, file.path(program_result_dir, "reactome_human_to_rat_one_to_one_mapping.tsv.gz"))

cat("Running Reactome camera for all prespecified contrasts...\n")
camera_current <- run_camera_contrasts(probe_expr, react_current$probe_index, design_primary, "Reactome_97")
write_tsv(camera_current, file.path(program_result_dir, "reactome_camera_primary_all_contrasts.tsv"))

cat("Running Reactome mroast sensitivity (9,999 rotations)...\n")
mroast_current <- run_mroast_primary(probe_expr, react_current$probe_index, design_primary)
write_tsv(mroast_current, file.path(program_result_dir, "reactome_mroast_interaction_sensitivity.tsv"))

cat("Computing Reactome sample-level activities and model-based CIs...\n")
react_activity_current <- program_mean_activity(gene_expr_current_all[, primary_samples, drop = FALSE], react_current$pathways)
react_activity_fit <- fit_limma_matrix(react_activity_current, design_primary, robust = TRUE)
react_activity_stats <- extract_fit_long(react_activity_fit, "feature_id")
write_tsv(react_activity_stats, file.path(program_result_dir, "reactome_mean_z_activity_all_contrasts.tsv"))
saveRDS(react_activity_current, file.path(processed_dir, "reactome_mean_z_activity_primary.rds"), compress = "xz")

cat("Running Reactome fgseaMultilevel sensitivity...\n")
fgsea_current <- run_fgsea_primary(primary_gene_stats[, .(gene_symbol, contrast, moderated_t)], react_current$pathways)
write_tsv(fgsea_current, file.path(program_result_dir, "reactome_fgsea_interaction_sensitivity.tsv"))

cat("Running original-GPL1355 annotation sensitivity...\n")
react_original <- build_reactome_mapping(frozen, rownames(gene_expr_original_all), original_map, fixed_terms = eligible_reactome)
write_tsv(react_original$coverage, file.path(program_result_dir, "reactome_mapping_coverage_original_GPL1355.tsv"))
original_index_fixed <- react_original$probe_index[intersect(eligible_reactome, names(react_original$probe_index))]
camera_original <- run_camera_contrasts(probe_expr, original_index_fixed, design_primary,
                                        "Reactome_97_original_GPL1355", "spasticity_progression_interaction")
write_tsv(camera_original, file.path(program_result_dir, "reactome_camera_original_annotation_sensitivity.tsv"))
react_activity_original <- program_mean_activity(gene_expr_original_all[, primary_samples, drop = FALSE],
                                                 react_original$pathways[intersect(eligible_reactome, names(react_original$pathways))])
react_activity_original_fit <- fit_limma_matrix(react_activity_original, design_primary, robust = TRUE)
react_activity_original_stats <- extract_fit_long(react_activity_original_fit, "feature_id")
write_tsv(react_activity_original_stats, file.path(program_result_dir, "reactome_activity_original_annotation_sensitivity.tsv"))

cat("Running all-27-array Reactome sensitivity...\n")
camera_all27 <- run_camera_contrasts(probe_expr_all, react_current$probe_index, design_all,
                                     "Reactome_97_all27", "spasticity_progression_interaction")
write_tsv(camera_all27, file.path(program_result_dir, "reactome_camera_all27_sensitivity.tsv"))
react_activity_all27 <- program_mean_activity(gene_expr_current_all, react_current$pathways)
react_activity_all27_fit <- fit_limma_matrix(react_activity_all27, design_all, robust = TRUE)
react_activity_all27_stats <- extract_fit_long(react_activity_all27_fit, "feature_id")
write_tsv(react_activity_all27_stats, file.path(program_result_dir, "reactome_activity_all27_sensitivity.tsv"))

cat("Computing Reactome leave-one-animal influence...\n")
react_loo <- leave_one_out_effects(react_activity_current, meta_primary)
write_tsv(react_loo$long, file.path(diagnostic_dir, "reactome_leave_one_animal_out_long.tsv.gz"))
write_tsv(react_loo$summary, file.path(diagnostic_dir, "reactome_leave_one_animal_out_summary.tsv"))

# ----------------------------------------------------------------------------
# PROGENy and DoRothEA activities (decoupleR ULM), with annotation/all-sample/LOO
# ----------------------------------------------------------------------------

cat("Computing PROGENy and DoRothEA activities using decoupleR ULM...\n")
net_current <- build_ulm_networks(frozen, rownames(gene_expr_current_all))
write_tsv(net_current$progeny_coverage, file.path(program_result_dir, "progeny_mapping_coverage_current.tsv"))
write_tsv(net_current$dorothea_coverage, file.path(program_result_dir, "dorothea_mapping_coverage_current.tsv"))
write_tsv(net_current$progeny, file.path(program_result_dir, "progeny_rat_network_used.tsv"))
write_tsv(net_current$dorothea, file.path(program_result_dir, "dorothea_rat_network_used.tsv"))

progeny_activity <- run_ulm_activity(gene_expr_current_all[, primary_samples, drop = FALSE], net_current$progeny)
dorothea_activity <- run_ulm_activity(gene_expr_current_all[, primary_samples, drop = FALSE], net_current$dorothea)
if (!nrow(progeny_activity) || !nrow(dorothea_activity)) stop("PROGENy or DoRothEA activity calculation returned no programs")
saveRDS(progeny_activity, file.path(processed_dir, "progeny_ulm_activity_primary.rds"), compress = "xz")
saveRDS(dorothea_activity, file.path(processed_dir, "dorothea_ulm_activity_primary.rds"), compress = "xz")

progeny_fit <- fit_limma_matrix(progeny_activity, design_primary, robust = TRUE)
dorothea_fit <- fit_limma_matrix(dorothea_activity, design_primary, robust = TRUE)
progeny_stats <- extract_fit_long(progeny_fit, "feature_id")
dorothea_stats <- extract_fit_long(dorothea_fit, "feature_id")
write_tsv(progeny_stats, file.path(program_result_dir, "progeny_ulm_activity_all_contrasts.tsv"))
write_tsv(dorothea_stats, file.path(program_result_dir, "dorothea_ulm_activity_all_contrasts.tsv"))

net_original <- build_ulm_networks(frozen, rownames(gene_expr_original_all))
progeny_activity_original <- run_ulm_activity(gene_expr_original_all[, primary_samples, drop = FALSE], net_original$progeny)
dorothea_activity_original <- run_ulm_activity(gene_expr_original_all[, primary_samples, drop = FALSE], net_original$dorothea)
progeny_original_stats <- extract_fit_long(fit_limma_matrix(progeny_activity_original, design_primary, robust = TRUE), "feature_id")
dorothea_original_stats <- extract_fit_long(fit_limma_matrix(dorothea_activity_original, design_primary, robust = TRUE), "feature_id")
write_tsv(progeny_original_stats, file.path(program_result_dir, "progeny_original_annotation_sensitivity.tsv"))
write_tsv(dorothea_original_stats, file.path(program_result_dir, "dorothea_original_annotation_sensitivity.tsv"))

progeny_activity_all27 <- run_ulm_activity(gene_expr_current_all, net_current$progeny)
dorothea_activity_all27 <- run_ulm_activity(gene_expr_current_all, net_current$dorothea)
progeny_all27_stats <- extract_fit_long(fit_limma_matrix(progeny_activity_all27, design_all, robust = TRUE), "feature_id")
dorothea_all27_stats <- extract_fit_long(fit_limma_matrix(dorothea_activity_all27, design_all, robust = TRUE), "feature_id")
write_tsv(progeny_all27_stats, file.path(program_result_dir, "progeny_all27_sensitivity.tsv"))
write_tsv(dorothea_all27_stats, file.path(program_result_dir, "dorothea_all27_sensitivity.tsv"))

progeny_loo <- leave_one_out_effects(progeny_activity, meta_primary)
dorothea_loo <- leave_one_out_effects(dorothea_activity, meta_primary)
write_tsv(progeny_loo$long, file.path(diagnostic_dir, "progeny_leave_one_animal_out_long.tsv.gz"))
write_tsv(progeny_loo$summary, file.path(diagnostic_dir, "progeny_leave_one_animal_out_summary.tsv"))
write_tsv(dorothea_loo$long, file.path(diagnostic_dir, "dorothea_leave_one_animal_out_long.tsv.gz"))
write_tsv(dorothea_loo$summary, file.path(diagnostic_dir, "dorothea_leave_one_animal_out_summary.tsv"))

# ----------------------------------------------------------------------------
# Prespecified Freedman-Lane residual permutation, preserving time strata
# ----------------------------------------------------------------------------

cat("Running 10,000 Freedman-Lane residual permutations across program activities...\n")
activity_combined <- rbind(
  react_activity_current,
  progeny_activity,
  dorothea_activity
)
family <- c(rep("Reactome", nrow(react_activity_current)), rep("PROGENy", nrow(progeny_activity)), rep("DoRothEA", nrow(dorothea_activity)))
original_id <- rownames(activity_combined)
rownames(activity_combined) <- paste(family, original_id, sep = "::")
fl <- freedman_lane(activity_combined, meta_primary, nperm = 10000L)
fl$result[, c("family", "feature_id_original") := tstrsplit(feature_id, "::", fixed = TRUE, keep = 1:2)]
fl$result[, empirical_q_value_BH := p.adjust(empirical_p_value, method = "BH"), by = family]
write_tsv(fl$result, file.path(diagnostic_dir, "program_Freedman_Lane_10000.tsv"))
write_tsv(fl$calibration, file.path(diagnostic_dir, "program_Freedman_Lane_null_calibration.tsv"))

# ----------------------------------------------------------------------------
# S1 gate: no story-based selection; every tested program retained in table
# ----------------------------------------------------------------------------

interaction_only <- function(x) x[contrast == "spasticity_progression_interaction"]

react_cam_int <- interaction_only(camera_current)
react_act_int <- interaction_only(react_activity_stats)
react_orig_int <- interaction_only(react_activity_original_stats)
react_all_int <- interaction_only(react_activity_all27_stats)
react_term <- unique(react_current$coverage[, .(feature_id = term_id, term_name)])
react_s1 <- merge(react_term, react_cam_int[, .(feature_id, camera_n_genes = n_genes, camera_direction = direction,
                                               camera_p_value = p_value, camera_q_value_BH = q_value_BH)],
                  by = "feature_id", all.y = TRUE)
react_s1 <- merge(react_s1, react_act_int[, .(feature_id, activity_estimate = estimate, activity_ci_low = ci_low,
                                             activity_ci_high = ci_high, activity_p_value = p_value,
                                             activity_q_value_BH = q_value_BH)], by = "feature_id", all.x = TRUE)
react_s1 <- merge(react_s1, react_orig_int[, .(feature_id, original_annotation_estimate = estimate)], by = "feature_id", all.x = TRUE)
react_s1 <- merge(react_s1, react_all_int[, .(feature_id, all27_estimate = estimate)], by = "feature_id", all.x = TRUE)
react_s1 <- merge(react_s1, react_loo$summary, by = "feature_id", all.x = TRUE)
react_s1 <- merge(react_s1, fl$result[family == "Reactome",
                                     .(feature_id = feature_id_original, permutation_p_value = empirical_p_value,
                                       permutation_q_value_BH = empirical_q_value_BH)], by = "feature_id", all.x = TRUE)
react_s1[, camera_sign := fifelse(camera_direction == "Up", 1L, fifelse(camera_direction == "Down", -1L, NA_integer_))]
react_s1[, activity_sign := safe_sign(activity_estimate)]
react_s1[, `:=`(
  primary_camera_pass = camera_q_value_BH <= 0.05,
  activity_ci_excludes_zero = (activity_ci_low > 0 & activity_ci_high > 0) | (activity_ci_low < 0 & activity_ci_high < 0),
  camera_activity_direction_agreement = camera_sign == activity_sign,
  original_annotation_direction_stable = safe_sign(original_annotation_estimate) == activity_sign,
  all27_direction_stable = safe_sign(all27_estimate) == activity_sign
)]
react_s1[, S1_pass := primary_camera_pass & activity_ci_excludes_zero & camera_activity_direction_agreement &
                         original_annotation_direction_stable & all27_direction_stable & no_single_animal_sign_flip]
react_s1[, collection := "Reactome"]

build_activity_s1 <- function(collection, primary_stats, original_stats, all_stats, loo_summary, perm_dt) {
  pri <- interaction_only(primary_stats)
  ori <- interaction_only(original_stats)
  alls <- interaction_only(all_stats)
  x <- merge(pri[, .(feature_id, activity_estimate = estimate, activity_standard_error = standard_error,
                     activity_ci_low = ci_low, activity_ci_high = ci_high, activity_p_value = p_value,
                     activity_q_value_BH = q_value_BH)],
             ori[, .(feature_id, original_annotation_estimate = estimate)], by = "feature_id", all.x = TRUE)
  x <- merge(x, alls[, .(feature_id, all27_estimate = estimate)], by = "feature_id", all.x = TRUE)
  x <- merge(x, loo_summary, by = "feature_id", all.x = TRUE)
  x <- merge(x, perm_dt[, .(feature_id = feature_id_original, permutation_p_value = empirical_p_value,
                            permutation_q_value_BH = empirical_q_value_BH)], by = "feature_id", all.x = TRUE)
  x[, activity_sign := safe_sign(activity_estimate)]
  x[, `:=`(
    primary_activity_pass = activity_q_value_BH <= 0.05,
    activity_ci_excludes_zero = (activity_ci_low > 0 & activity_ci_high > 0) | (activity_ci_low < 0 & activity_ci_high < 0),
    original_annotation_direction_stable = safe_sign(original_annotation_estimate) == activity_sign,
    all27_direction_stable = safe_sign(all27_estimate) == activity_sign
  )]
  x[, S1_pass := primary_activity_pass & activity_ci_excludes_zero & original_annotation_direction_stable &
                    all27_direction_stable & no_single_animal_sign_flip]
  x[, collection := collection]
  x[]
}

progeny_s1 <- build_activity_s1("PROGENy", progeny_stats, progeny_original_stats, progeny_all27_stats,
                                progeny_loo$summary, fl$result[family == "PROGENy"])
dorothea_s1 <- build_activity_s1("DoRothEA", dorothea_stats, dorothea_original_stats, dorothea_all27_stats,
                                 dorothea_loo$summary, fl$result[family == "DoRothEA"])

write_tsv(react_s1, file.path(program_result_dir, "S1_gate_Reactome.tsv"))
write_tsv(progeny_s1, file.path(program_result_dir, "S1_gate_PROGENy.tsv"))
write_tsv(dorothea_s1, file.path(program_result_dir, "S1_gate_DoRothEA.tsv"))

s1_summary <- rbindlist(list(
  react_s1[, .(collection, feature_id, feature_name = term_name, S1_pass, direction = activity_sign,
               primary_effect = activity_estimate, primary_ci_low = activity_ci_low, primary_ci_high = activity_ci_high,
               primary_q_value = camera_q_value_BH, permutation_q_value_BH, loo_direction_fraction,
               original_annotation_direction_stable, all27_direction_stable)],
  progeny_s1[, .(collection, feature_id, feature_name = feature_id, S1_pass, direction = activity_sign,
                 primary_effect = activity_estimate, primary_ci_low = activity_ci_low, primary_ci_high = activity_ci_high,
                 primary_q_value = activity_q_value_BH, permutation_q_value_BH, loo_direction_fraction,
                 original_annotation_direction_stable, all27_direction_stable)],
  dorothea_s1[, .(collection, feature_id, feature_name = feature_id, S1_pass, direction = activity_sign,
                  primary_effect = activity_estimate, primary_ci_low = activity_ci_low, primary_ci_high = activity_ci_high,
                  primary_q_value = activity_q_value_BH, permutation_q_value_BH, loo_direction_fraction,
                  original_annotation_direction_stable, all27_direction_stable)]
), fill = TRUE)
setorder(s1_summary, collection, -S1_pass, primary_q_value, feature_id)
write_tsv(s1_summary, file.path(program_result_dir, "S1_gate_all_collections.tsv"))
write_tsv(s1_summary, file.path(program_result_dir, "S1_GATE.tsv"))

# ----------------------------------------------------------------------------
# Diagnostics and non-story-based figures
# ----------------------------------------------------------------------------

gene_primary <- primary_gene_stats[!is.na(gene_symbol)]
gene_primary[, significance := q_value_BH <= 0.05]
p_volcano <- ggplot(gene_primary, aes(x = estimate, y = -log10(pmax(p_value, .Machine$double.xmin)), colour = significance)) +
  geom_point(alpha = 0.45, size = 0.8) + scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#B2182B")) +
  geom_vline(xintercept = 0, colour = "grey40", linetype = 2) + theme_bw(base_size = 11) +
  labs(x = "Interaction log2 expression contrast", y = "-log10(P)",
       title = "GSE16710 prespecified injury-by-time interaction", colour = "BH q <= 0.05")
ggsave(file.path(figure_dir, "GSE16710_gene_interaction_volcano.pdf"), p_volcano, width = 8, height = 6)
ggsave(file.path(figure_dir, "GSE16710_gene_interaction_volcano.png"), p_volcano, width = 8, height = 6, dpi = 300)

forest <- head(react_s1[order(camera_q_value_BH)], 25L)
forest[, label := paste0(term_name, " [", feature_id, "]")]
forest[, label := factor(label, levels = rev(label))]
p_forest <- ggplot(forest, aes(x = activity_estimate, y = label, colour = S1_pass)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = activity_ci_low, xmax = activity_ci_high), height = 0.15) +
  geom_point(size = 2) + theme_bw(base_size = 9) +
  labs(x = "Mean-z program interaction effect (95% CI)", y = NULL,
       title = "Reactome programs ranked only by prespecified camera q value")
ggsave(file.path(figure_dir, "GSE16710_Reactome_top25_statistical_forest.pdf"), p_forest, width = 11, height = 9)

perm_plot_dt <- merge(
  rbindlist(list(
    interaction_only(react_activity_stats)[, .(family = "Reactome", feature_id, nominal_p = p_value)],
    interaction_only(progeny_stats)[, .(family = "PROGENy", feature_id, nominal_p = p_value)],
    interaction_only(dorothea_stats)[, .(family = "DoRothEA", feature_id, nominal_p = p_value)]
  )),
  fl$result[, .(family, feature_id = feature_id_original, empirical_p = empirical_p_value)],
  by = c("family", "feature_id")
)
p_perm <- ggplot(perm_plot_dt, aes(x = -log10(pmax(nominal_p, 1e-12)), y = -log10(pmax(empirical_p, 1e-12)))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") + geom_point(alpha = 0.45, size = 0.8) +
  facet_wrap(~ family, scales = "free") + theme_bw(base_size = 10) +
  labs(x = "Nominal activity-model -log10(P)", y = "Freedman-Lane -log10(empirical P)",
       title = "Design-respecting label-null calibration")
ggsave(file.path(figure_dir, "GSE16710_program_permutation_calibration.pdf"), p_perm, width = 10, height = 4)

gene_counts <- gene_primary[, .(
  n_tested = .N,
  n_BH_q_le_0_05 = sum(q_value_BH <= 0.05, na.rm = TRUE),
  n_positive_q_le_0_05 = sum(q_value_BH <= 0.05 & estimate > 0, na.rm = TRUE),
  n_negative_q_le_0_05 = sum(q_value_BH <= 0.05 & estimate < 0, na.rm = TRUE)
)]
program_counts <- s1_summary[, .(n_tested = .N, n_S1_pass = sum(S1_pass, na.rm = TRUE)), by = collection]
write_tsv(gene_counts, file.path(gene_result_dir, "gene_interaction_summary_counts.tsv"))
write_tsv(program_counts, file.path(program_result_dir, "S1_summary_counts.tsv"))

run_summary <- list(
  dataset = "GSE16710",
  role = "spasticity_enriched_cross_disease_phenotype_anchor",
  primary_statistical_unit = "independent_rat_array",
  n_raw_arrays = ncol(probe_expr_all),
  n_primary_arrays = ncol(probe_expr),
  outcome_blind_qc_exclusions = excluded_samples,
  primary_group_counts = as.list(table(paste(meta_primary$experimental_arm, meta_primary$time_day, sep = "_"))),
  primary_model = "expression ~ injury * time_day",
  primary_estimand = "(Spi60-Spi21)-(Sham60-Sham21)",
  model_rank = qr(design_primary$design)$rank,
  residual_df = nrow(meta_primary) - ncol(design_primary$design),
  robust_empirical_bayes = TRUE,
  current_annotation = "rat2302.db_3.13.0 + org.Rn.eg.db_3.23.0",
  original_annotation_sensitivity = "GEO GPL1355 annotation dated 2016-08-09",
  reactome_release = "97",
  progeny_release = "1.34.0",
  dorothea_release = "1.23.0",
  n_S1_Reactome = program_counts[collection == "Reactome", n_S1_pass],
  n_S1_PROGENy = program_counts[collection == "PROGENy", n_S1_pass],
  n_S1_DoRothEA = program_counts[collection == "DoRothEA", n_S1_pass],
  inferential_ceiling = "associated_with_progression_of_experimental_spasticity; not_poststroke_spasticity_specific",
  completed_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
)
jsonlite::write_json(run_summary, file.path(result_dir, "GSE16710_analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_GSE16710_analysis.txt"))
source(file.path(project_root, "code", "R", "summarize_GSE16710_precision.R"), local = new.env(parent = globalenv()))

cat("GSE16710 phenotype-anchor modeling completed:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
