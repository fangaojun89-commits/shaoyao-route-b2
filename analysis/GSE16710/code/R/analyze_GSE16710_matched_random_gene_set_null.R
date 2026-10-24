#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 160)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
project_library <- normalizePath(file.path(project_root, "environment", "R-library", "4.6"), mustWork = TRUE)
.libPaths(c(project_library, .libPaths()))

required <- c("data.table", "limma", "matrixStats", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(matrixStats)
})

processed_dir <- file.path(project_root, "data", "processed", "GSE16710")
gene_set_path <- file.path(project_root, "data", "processed", "gene_sets", "gene_set_resources_frozen_2026-08-18.rds")
spec_path <- file.path(project_root, "protocol", "GSE16710_matched_random_gene_set_null_addendum_2026-08-31.md")
final_dir <- file.path(project_root, "results", "GSE16710", "diagnostics", "matched_random_gene_set_null_v1")
if (dir.exists(final_dir) || file.exists(final_dir)) stop("Final output already exists: ", final_dir)
staging_dir <- paste0(final_dir, "_STAGING_", Sys.getpid())
if (dir.exists(staging_dir) || file.exists(staging_dir)) stop("Staging path already exists: ", staging_dir)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(staging_dir, "matched_random_gene_set_null_run.log")
zz <- file(log_path, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink(type = "output")
  try(close(zz), silent = TRUE)
}, add = TRUE)

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), file = path, sep = "\t", quote = FALSE, na = "NA")
}

sha256_file <- function(path) {
  if (!file.exists(path)) stop("Cannot hash missing file: ", path)
  tolower(unname(tools::sha256sum(path)))
}

rank_bins <- function(values, ids, n_bins) {
  if (length(values) != length(ids) || anyDuplicated(ids)) stop("Invalid rank-bin inputs")
  ord <- order(values, ids, na.last = NA)
  if (length(ord) != length(values)) stop("Non-finite values entered rank bins")
  out <- integer(length(values))
  out[ord] <- pmin(n_bins, ceiling(seq_along(ord) * n_bins / length(ord)))
  names(out) <- ids
  out
}

make_design <- function(meta) {
  meta <- copy(meta)
  meta[, injury := factor(experimental_arm, levels = c("Sham", "Spi"))]
  meta[, time_factor := factor(as.character(time_day), levels = c("21", "60"))]
  design <- stats::model.matrix(~ injury * time_factor, data = meta)
  expected <- c("(Intercept)", "injurySpi", "time_factor60", "injurySpi:time_factor60")
  if (!identical(colnames(design), expected)) stop("Unexpected design columns")
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design")
  rownames(design) <- meta$sample_id
  design
}

cat("Matched random-gene-set calibration started:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
if (!file.exists(spec_path)) stop("Frozen addendum is missing")

n_random <- 10000L
n_expression_bins <- 20L
n_candidates <- 5L
sampling_seed <- 20260831L
family_reference_seed <- 20260832L
n_family_reference <- 10000L

expression_path <- file.path(processed_dir, "rma_gene_expression_current_bioconductor.rds")
ledger_path <- file.path(processed_dir, "sample_ledger_verified.tsv")
decision_path <- file.path(processed_dir, "qc_inclusion_decision.tsv")

expr_all <- readRDS(expression_path)
ledger <- fread(ledger_path, sep = "\t")
decision <- fread(decision_path, sep = "\t")
frozen <- readRDS(gene_set_path)

if (!is.matrix(expr_all) || anyDuplicated(rownames(expr_all)) || anyDuplicated(colnames(expr_all))) {
  stop("Expression object must be a uniquely named matrix")
}
if (!setequal(colnames(expr_all), ledger$sample_id)) stop("Expression columns and sample ledger differ")
if (!setequal(ledger$sample_id, decision$sample_id)) stop("QC decision does not cover the sample ledger")
if (!all(toupper(as.character(decision$outcome_information_used)) == "FALSE")) {
  stop("QC decision records outcome information use")
}

ledger <- merge(ledger, decision, by = "sample_id", all.x = TRUE, sort = FALSE)
ledger[, sample_order_int__ := as.integer(sample_order)]
setorder(ledger, sample_order_int__)
ledger[, sample_order_int__ := NULL]
primary_samples <- ledger[primary_analysis_decision == "retain", sample_id]
if (!identical(ledger[primary_analysis_decision != "retain", sample_id], "GSM418682")) {
  stop("Unexpected primary exclusion set")
}
expr <- expr_all[, primary_samples, drop = FALSE]
gene_mean <- rowMeans(expr)
gene_sd <- matrixStats::rowSds(expr)
gene_keep <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
expr <- expr[gene_keep, , drop = FALSE]
gene_mean <- gene_mean[gene_keep]
gene_sd <- gene_sd[gene_keep]
z_expr <- sweep(expr, 1L, gene_mean, "-")
z_expr <- sweep(z_expr, 1L, gene_sd, "/")
if (any(!is.finite(z_expr))) stop("Non-finite standardized expression")

ortho_rat <- frozen$ortholog_one_to_one[
  other_tax_id == 10116,
  .(human_symbol, rat_symbol = other_symbol)
]
react <- merge(
  frozen$reactome[, .(term_id, term_name, human_symbol = gene_symbol)],
  ortho_rat,
  by = "human_symbol",
  allow.cartesian = TRUE
)
react <- unique(react[rat_symbol %in% rownames(z_expr), .(term_id, term_name, rat_symbol)])
term_names <- unique(react[, .(term_id, term_name)])
if (any(term_names[, uniqueN(term_name), by = term_id]$V1 != 1L)) stop("Reactome term ID has multiple names")
pathways <- lapply(split(react$rat_symbol, react$term_id), function(x) sort(unique(x)))
pathways <- pathways[lengths(pathways) >= 10L & lengths(pathways) <= 500L]
pathways <- pathways[sort(names(pathways))]
if (length(pathways) != 1656L) stop("Unexpected eligible Reactome program count: ", length(pathways))

gene_bins <- rank_bins(gene_mean, names(gene_mean), n_expression_bins)
universe_by_bin <- split(names(gene_bins), factor(gene_bins, levels = seq_len(n_expression_bins)))
if (any(lengths(universe_by_bin) == 0L)) stop("Empty expression bin")

template_dt <- data.table(
  template_term_id = names(pathways),
  template_size = lengths(pathways),
  template_mean_expression = vapply(pathways, function(g) mean(gene_mean[g]), numeric(1))
)
template_dt <- merge(
  template_dt,
  term_names[, .(template_term_id = term_id, template_term_name = term_name)],
  by = "template_term_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(template_dt, template_term_id)
if (anyNA(template_dt$template_term_name)) stop("Missing Reactome term name")

real_activity <- t(vapply(pathways, function(g) colMeans(z_expr[g, , drop = FALSE]), numeric(ncol(z_expr))))
colnames(real_activity) <- colnames(z_expr)
rownames(real_activity) <- names(pathways)

set.seed(sampling_seed)
base_reps <- n_random %/% length(pathways)
remainder <- n_random - base_reps * length(pathways)
template_assignment <- rep(names(pathways), each = base_reps)
if (remainder > 0L) template_assignment <- c(template_assignment, sample(names(pathways), remainder, replace = FALSE))
template_assignment <- sample(template_assignment, length(template_assignment), replace = FALSE)
if (length(template_assignment) != n_random) stop("Template allocation length mismatch")

random_activity <- matrix(
  NA_real_,
  nrow = n_random,
  ncol = ncol(z_expr),
  dimnames = list(sprintf("RGS%05d", seq_len(n_random)), colnames(z_expr))
)
membership <- vector("list", n_random)
random_meta <- vector("list", n_random)

for (i in seq_len(n_random)) {
  random_set_id <- rownames(random_activity)[i]
  template_id <- template_assignment[i]
  template_genes <- pathways[[template_id]]
  template_bin_counts <- tabulate(gene_bins[template_genes], nbins = n_expression_bins)
  template_mean <- mean(gene_mean[template_genes])
  candidate_sets <- vector("list", n_candidates)
  candidate_delta <- rep(Inf, n_candidates)

  for (j in seq_len(n_candidates)) {
    selected <- unlist(lapply(seq_len(n_expression_bins), function(bin_id) {
      k <- template_bin_counts[bin_id]
      if (k == 0L) return(character())
      pool <- universe_by_bin[[bin_id]]
      if (length(pool) < k) stop("Expression bin cannot supply template count")
      sample(pool, k, replace = FALSE)
    }), use.names = FALSE)
    selected <- sort(selected)
    if (length(selected) != length(template_genes) || anyDuplicated(selected)) stop("Invalid random set")
    candidate_sets[[j]] <- selected
    candidate_delta[j] <- abs(mean(gene_mean[selected]) - template_mean)
  }

  chosen_index <- which.min(candidate_delta)
  selected <- candidate_sets[[chosen_index]]
  selected_bin_counts <- tabulate(gene_bins[selected], nbins = n_expression_bins)
  if (!identical(selected_bin_counts, template_bin_counts)) stop("Expression-bin composition mismatch")

  random_activity[i, ] <- colMeans(z_expr[selected, , drop = FALSE])
  overlap <- length(intersect(selected, template_genes))
  random_meta[[i]] <- data.table(
    random_set_id = random_set_id,
    template_term_id = template_id,
    template_size = length(template_genes),
    random_size = length(selected),
    template_mean_expression = template_mean,
    random_mean_expression = mean(gene_mean[selected]),
    absolute_mean_expression_difference = abs(mean(gene_mean[selected]) - template_mean),
    chosen_candidate_index = chosen_index,
    overlap_genes = overlap,
    jaccard_with_template = overlap / length(union(selected, template_genes)),
    exact_size_match = length(selected) == length(template_genes),
    exact_expression_bin_composition_match = identical(selected_bin_counts, template_bin_counts)
  )
  membership[[i]] <- data.table(
    random_set_id = random_set_id,
    template_term_id = template_id,
    gene_symbol = selected,
    expression_bin = gene_bins[selected],
    gene_mean_expression = gene_mean[selected]
  )
  if (i %% 1000L == 0L) cat("Constructed random sets:", i, "/", n_random, "\n")
}

random_meta <- rbindlist(random_meta)
membership <- rbindlist(membership)

meta_primary <- ledger[match(primary_samples, sample_id)]
design <- make_design(meta_primary)
if (!identical(colnames(real_activity), rownames(design)) || !identical(colnames(random_activity), rownames(design))) {
  stop("Activity columns and design rows differ")
}

combined_activity <- rbind(real_activity, random_activity)
rownames(combined_activity) <- c(paste0("REAL::", rownames(real_activity)), paste0("NULL::", rownames(random_activity)))
fit <- limma::lmFit(combined_activity, design)
fit <- limma::eBayes(fit, robust = TRUE)
coef_name <- "injurySpi:time_factor60"
coef_index <- match(coef_name, colnames(fit$coefficients))
if (is.na(coef_index)) stop("Interaction coefficient not found")
se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
model_dt <- data.table(
  combined_id = rownames(fit$coefficients),
  estimate = fit$coefficients[, coef_index],
  standard_error = se,
  moderated_t = fit$t[, coef_index],
  df_total = rep_len(fit$df.total, nrow(fit$coefficients)),
  p_value = fit$p.value[, coef_index]
)
model_dt[, absolute_moderated_t := abs(moderated_t)]

real_model <- model_dt[startsWith(combined_id, "REAL::")]
real_model[, template_term_id := sub("^REAL::", "", combined_id)]
real_model[, combined_id := NULL]
real_model[, q_value_BH_within_real_family := p.adjust(p_value, method = "BH")]
real_model <- merge(template_dt, real_model, by = "template_term_id", all.x = TRUE, sort = FALSE)
setorder(real_model, template_term_id)

null_model <- model_dt[startsWith(combined_id, "NULL::")]
null_model[, random_set_id := sub("^NULL::", "", combined_id)]
null_model[, combined_id := NULL]
null_model[, q_value_BH_within_random_family := p.adjust(p_value, method = "BH")]
random_results <- merge(random_meta, null_model, by = "random_set_id", all.x = TRUE, sort = FALSE)
random_results <- merge(
  random_results,
  real_model[, .(
    template_term_id,
    template_term_name,
    template_absolute_moderated_t = absolute_moderated_t,
    template_moderated_t = moderated_t,
    template_calibration_model_p_value = p_value
  )],
  by = "template_term_id",
  all.x = TRUE,
  sort = FALSE
)
random_results[, null_at_least_as_extreme_as_template := absolute_moderated_t >= template_absolute_moderated_t]
setorder(random_results, random_set_id)

template_summary <- random_results[, .(
  n_matched_random_sets = .N,
  mean_random_absolute_moderated_t = mean(absolute_moderated_t),
  median_random_absolute_moderated_t = median(absolute_moderated_t),
  random_at_least_as_extreme_fraction = mean(null_at_least_as_extreme_as_template),
  coarse_descriptive_empirical_tail = (1 + sum(null_at_least_as_extreme_as_template)) / (1 + .N),
  max_absolute_mean_expression_difference = max(absolute_mean_expression_difference)
), by = template_term_id]
template_summary <- merge(real_model, template_summary, by = "template_term_id", all.x = TRUE, sort = FALSE)
template_summary[, pathway_specific_inference_permitted := FALSE]
setorder(template_summary, template_term_id)

set.seed(family_reference_seed)
null_by_template <- split(random_results$absolute_moderated_t, random_results$template_term_id)
if (!setequal(names(null_by_template), real_model$template_term_id)) stop("Null template coverage differs from real family")
family_reference <- numeric(n_family_reference)
for (template_id in real_model$template_term_id) {
  values <- null_by_template[[template_id]]
  family_reference <- family_reference + sample(values, n_family_reference, replace = TRUE)
}
family_reference <- family_reference / nrow(real_model)
observed_family_mean <- mean(real_model$absolute_moderated_t)
family_reference_dt <- data.table(
  reference_draw = seq_len(n_family_reference),
  mean_absolute_moderated_t = family_reference
)
exploratory_family_tail <- (1 + sum(family_reference >= observed_family_mean)) / (1 + n_family_reference)

q_real <- stats::quantile(real_model$absolute_moderated_t, c(0.5, 0.95, 0.99, 1), names = FALSE, type = 8)
q_null <- stats::quantile(random_results$absolute_moderated_t, c(0.5, 0.95, 0.99, 1), names = FALSE, type = 8)
family_summary <- data.table(
  analysis_status = "PERFORMED_POST_OUTCOME_EXPLORATORY_ADJUNCT",
  primary_S1_state_unchanged = "COMPLETE_ZERO_PASS",
  eligible_reactome_programs = nrow(real_model),
  random_sets = nrow(random_results),
  random_sets_per_template_min = min(template_summary$n_matched_random_sets),
  random_sets_per_template_max = max(template_summary$n_matched_random_sets),
  median_absolute_mean_expression_difference = median(random_results$absolute_mean_expression_difference),
  p95_absolute_mean_expression_difference = as.numeric(stats::quantile(random_results$absolute_mean_expression_difference, 0.95, type = 8)),
  max_absolute_mean_expression_difference = max(random_results$absolute_mean_expression_difference),
  template_balanced_random_at_least_as_extreme_fraction = mean(template_summary$random_at_least_as_extreme_fraction),
  real_mean_absolute_moderated_t = observed_family_mean,
  random_mean_absolute_moderated_t = mean(random_results$absolute_moderated_t),
  real_median_absolute_moderated_t = q_real[1],
  random_median_absolute_moderated_t = q_null[1],
  real_p95_absolute_moderated_t = q_real[2],
  random_p95_absolute_moderated_t = q_null[2],
  real_p99_absolute_moderated_t = q_real[3],
  random_p99_absolute_moderated_t = q_null[3],
  real_max_absolute_moderated_t = q_real[4],
  random_max_absolute_moderated_t = q_null[4],
  real_nominal_p_le_0_05 = sum(real_model$p_value <= 0.05),
  random_nominal_p_le_0_05 = sum(random_results$p_value <= 0.05),
  real_BH_q_le_0_05_calibration_model = sum(real_model$q_value_BH_within_real_family <= 0.05),
  random_BH_q_le_0_05_calibration_model = sum(random_results$q_value_BH_within_random_family <= 0.05),
  family_reference_draws = n_family_reference,
  family_reference_median = median(family_reference),
  family_reference_p025 = as.numeric(stats::quantile(family_reference, 0.025, type = 8)),
  family_reference_p975 = as.numeric(stats::quantile(family_reference, 0.975, type = 8)),
  exploratory_family_reference_tail_area = exploratory_family_tail,
  confirmatory_p_value = NA_real_,
  candidate_gate_eligible = FALSE,
  pathway_specific_inference_permitted = FALSE
)

membership_counts <- membership[, .(membership_rows = .N, distinct_genes = uniqueN(gene_symbol)), by = random_set_id]
validation <- data.table(
  check_id = c(
    "primary_sample_count_26",
    "eligible_reactome_program_count_1656",
    "exactly_10000_random_sets",
    "template_allocation_balanced_6_or_7",
    "random_set_size_exact",
    "expression_bin_composition_exact",
    "membership_count_matches_random_size",
    "no_duplicate_gene_within_random_set",
    "all_model_statistics_finite",
    "real_family_BH_denom_1656",
    "calibration_gate_ineligible",
    "pathway_specific_inference_forbidden",
    "post_outcome_status_explicit",
    "primary_S1_state_unchanged"
  ),
  pass = c(
    length(primary_samples) == 26L,
    length(pathways) == 1656L,
    nrow(random_results) == 10000L,
    all(template_summary$n_matched_random_sets %in% c(6L, 7L)),
    all(random_results$exact_size_match),
    all(random_results$exact_expression_bin_composition_match),
    all(membership_counts$membership_rows == random_results$random_size[match(membership_counts$random_set_id, random_results$random_set_id)]),
    all(membership_counts$membership_rows == membership_counts$distinct_genes),
    all(is.finite(real_model$moderated_t)) && all(is.finite(random_results$moderated_t)),
    nrow(real_model) == 1656L,
    identical(family_summary$candidate_gate_eligible, FALSE),
    identical(family_summary$pathway_specific_inference_permitted, FALSE),
    identical(family_summary$analysis_status, "PERFORMED_POST_OUTCOME_EXPLORATORY_ADJUNCT"),
    identical(family_summary$primary_S1_state_unchanged, "COMPLETE_ZERO_PASS")
  )
)
validation[, observed := c(
  as.character(length(primary_samples)),
  as.character(length(pathways)),
  as.character(nrow(random_results)),
  paste(range(template_summary$n_matched_random_sets), collapse = "-"),
  as.character(all(random_results$exact_size_match)),
  as.character(all(random_results$exact_expression_bin_composition_match)),
  as.character(all(membership_counts$membership_rows == random_results$random_size[match(membership_counts$random_set_id, random_results$random_set_id)])),
  as.character(all(membership_counts$membership_rows == membership_counts$distinct_genes)),
  as.character(all(is.finite(real_model$moderated_t)) && all(is.finite(random_results$moderated_t))),
  as.character(nrow(real_model)),
  as.character(family_summary$candidate_gate_eligible),
  as.character(family_summary$pathway_specific_inference_permitted),
  family_summary$analysis_status,
  family_summary$primary_S1_state_unchanged
)]
validation[, expected := c(
  "26", "1656", "10000", "6-7", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "1656",
  "FALSE", "FALSE", "PERFORMED_POST_OUTCOME_EXPLORATORY_ADJUNCT", "COMPLETE_ZERO_PASS"
)]
validation[, status := fifelse(pass, "PASS", "FAIL")]
if (any(!validation$pass)) stop("Validation failed: ", paste(validation[!pass, check_id], collapse = ", "))

membership_path <- file.path(staging_dir, "matched_random_gene_set_membership_10000.tsv.gz")
random_results_path <- file.path(staging_dir, "matched_random_gene_set_null_10000.tsv.gz")
real_model_path <- file.path(staging_dir, "matched_random_gene_set_real_reactome_calibration_model.tsv.gz")
template_summary_path <- file.path(staging_dir, "matched_random_gene_set_template_summary.tsv.gz")
family_reference_path <- file.path(staging_dir, "matched_random_gene_set_family_reference_10000.tsv.gz")
family_summary_path <- file.path(staging_dir, "matched_random_gene_set_null_family_summary.tsv")
validation_path <- file.path(staging_dir, "matched_random_gene_set_null_validation.tsv")

write_tsv(membership, membership_path)
write_tsv(random_results, random_results_path)
write_tsv(real_model, real_model_path)
write_tsv(template_summary, template_summary_path)
write_tsv(family_reference_dt, family_reference_path)
write_tsv(family_summary, family_summary_path)
write_tsv(validation, validation_path)
writeLines(capture.output(sessionInfo()), file.path(staging_dir, "sessionInfo_matched_random_gene_set_null.txt"))

sink(type = "message")
sink(type = "output")
close(zz)

output_paths <- c(
  membership_path,
  random_results_path,
  real_model_path,
  template_summary_path,
  family_reference_path,
  family_summary_path,
  validation_path,
  log_path,
  file.path(staging_dir, "sessionInfo_matched_random_gene_set_null.txt")
)
output_manifest <- data.table(
  file_name = basename(output_paths),
  bytes = as.numeric(file.info(output_paths)$size),
  sha256 = vapply(output_paths, sha256_file, character(1))
)
write_tsv(output_manifest, file.path(staging_dir, "matched_random_gene_set_null_output_manifest.tsv"))

run_summary <- list(
  schema_version = "route_b2_gse16710_matched_random_gene_set_null_v1",
  status = "COMPLETE_POST_OUTCOME_EXPLORATORY_ADJUNCT",
  started_after_primary_outcome_known = TRUE,
  prospective_or_outcome_protected = FALSE,
  primary_S1_state_unchanged = "COMPLETE_ZERO_PASS",
  candidate_gate_eligible = FALSE,
  dataset = "GSE16710",
  organism = "Rattus norvegicus",
  tissue = "LCM sacral tail motor neurons",
  statistical_unit = "independent rat array",
  retained_primary_arrays = length(primary_samples),
  eligible_reactome_programs = nrow(real_model),
  random_sets = nrow(random_results),
  expression_bins = n_expression_bins,
  candidate_sets_per_draw = n_candidates,
  sampling_seed = sampling_seed,
  family_reference_seed = family_reference_seed,
  family_reference_draws = n_family_reference,
  calibration_statistic = "absolute robust-limma moderated t for injurySpi:time_factor60 on mean row-z activity",
  interpretation = "family-level exploratory calibration only; no pathway-specific inference and no change to primary camera or S1",
  frozen_addendum_relative_path = "protocol/GSE16710_matched_random_gene_set_null_addendum_2026-08-31.md",
  frozen_addendum_sha256 = sha256_file(spec_path),
  script_relative_path = "code/R/analyze_GSE16710_matched_random_gene_set_null.R",
  script_sha256 = sha256_file(script_path),
  inputs = list(
    expression = list(relative_path = "data/processed/GSE16710/rma_gene_expression_current_bioconductor.rds", sha256 = sha256_file(expression_path)),
    sample_ledger = list(relative_path = "data/processed/GSE16710/sample_ledger_verified.tsv", sha256 = sha256_file(ledger_path)),
    qc_decision = list(relative_path = "data/processed/GSE16710/qc_inclusion_decision.tsv", sha256 = sha256_file(decision_path)),
    frozen_gene_sets = list(relative_path = "data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds", sha256 = sha256_file(gene_set_path))
  ),
  family_summary = as.list(family_summary[1]),
  validation = list(passed = sum(validation$pass), total = nrow(validation)),
  output_manifest_sha256 = sha256_file(file.path(staging_dir, "matched_random_gene_set_null_output_manifest.tsv")),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "Asia/Shanghai")
)
jsonlite::write_json(
  run_summary,
  file.path(staging_dir, "matched_random_gene_set_null_COMPLETE.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16,
  na = "null"
)

cat("Validation:", sum(validation$pass), "/", nrow(validation), "PASS\n")
cat("Exploratory family-reference tail area:", exploratory_family_tail, "\n")
cat("Primary S1 state remains COMPLETE_ZERO_PASS\n")

if (!file.rename(staging_dir, final_dir)) stop("Failed to atomically publish the completed output directory")
cat("Published completed output:", final_dir, "\n")
