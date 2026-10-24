#!/usr/bin/env Rscript
# New equivalent reconstruction only. Never source the full historical entry,
# replace old outputs, infer historical execution, or tune against old tables.
suppressPackageStartupMessages({library(edgeR); library(limma); library(jsonlite); library(digest); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L, args[[3]] %in% c("preflight", "run"))
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
contract_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
mode <- args[[3]]
out <- gsub("\\\\", "/", args[[4]])
stopifnot(grepl("^R:/Shaoyao_storage/", out), !dir.exists(out),
          grepl("^R:/Shaoyao_storage/", normalizePath(tempdir(), winslash = "/")))
contract <- jsonlite::read_json(contract_path, simplifyVector = FALSE)
stopifnot(contract$scope == "NEW_EQUIVALENT_RECONSTRUCTION_NOT_HISTORICAL_AUTHENTICATION",
          contract$result_comparison_tolerance == 0,
          contract$seed == 20260817L, contract$nrot == 1999L,
          contract$min_genes == 10L, contract$max_genes == 500L)
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
sha <- function(p) digest::digest(file = p, algo = "sha256")
paths <- setNames(lapply(contract$files, function(f) file.path(root, f$relative_path)),
                  vapply(contract$files, `[[`, character(1), "key"))
check_inputs <- function() {
  lapply(contract$files, function(f) {
    p <- paths[[f$key]]
    stopifnot(file.exists(p), file.info(p)$size == f$size_bytes, sha(p) == f$sha256)
    list(key = f$key, relative_path = f$relative_path, size_bytes = f$size_bytes, sha256 = f$sha256)
  })
}
before <- check_inputs()
stopifnot(as.character(getRversion()) == contract$versions$R,
          as.character(packageVersion("edgeR")) == contract$versions$edgeR,
          as.character(packageVersion("limma")) == contract$versions$limma)
config <- yaml::read_yaml(paths$contrasts)
gs_config <- yaml::read_yaml(paths$gene_set_config)
stopifnot(config$global$random_seed == contract$seed,
          gs_config$testing$min_genes_after_mapping == contract$min_genes,
          gs_config$testing$max_genes_after_mapping == contract$max_genes,
          gs_config$primary_collections$reactome$release == 97L)
# Import ONLY two pure definitions from the pinned source AST; no top-level
# imports, directories, configurations, old run state, or model calls execute.
ast <- parse(paths$source)
definition <- function(name) {
  hit <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name(name)), as.list(ast))
  stopifnot(length(hit) == 1L, identical(hit[[1]][[3]][[1]], as.name("function")))
  hit[[1]]
}
eval(definition("write_tsv"))
eval(definition("collapse_expression_by_symbol"))
read_tsv <- function(p) utils::read.delim(p, check.names = FALSE, na.strings = c("", "NA"))
samples <- read_tsv(paths$samples)
ledger <- read_tsv(paths$sample_ledger)
ledger <- ledger[ledger$dataset_id == "GSE289929" & ledger$include_primary == "TRUE", , drop = FALSE]
stopifnot(nrow(samples) == 12L, !anyNA(samples$sample_id), !anyDuplicated(samples$sample_id),
          setequal(samples$sample_id, ledger$sample_id), !anyDuplicated(ledger$sample_id))
ledger <- ledger[match(samples$sample_id, ledger$sample_id), , drop = FALSE]
for (nm in intersect(names(ledger), names(samples)))
  stopifnot(identical(as.character(ledger[[nm]]), as.character(samples[[nm]])))
# Original normalize_sample_table retains row names of the full sample ledger;
# the exported TSV does not. Restore these cosmetic labels from that same ledger,
# not by copying any cached design coefficients or changing sample order.
rownames(samples) <- rownames(ledger)
stopifnot(all(samples$dataset_id == "GSE289929"), all(samples$exposure_duration_hour == 24),
          all(samples$tissue_or_cell_line == "PC12"), all(samples$batch == "A"))
cfg <- config$datasets$GSE289929$canonical_subset$gsm_sets
stopifnot(setequal(samples$sample_id[samples$dose_numeric == 0], unlist(cfg$vehicle)),
          setequal(samples$sample_id[samples$dose_numeric == 4], unlist(cfg$dose_4_ug_ml)),
          setequal(samples$sample_id[samples$dose_numeric == 20], unlist(cfg$dose_20_ug_ml)),
          setequal(samples$sample_id[samples$dose_numeric == 100], unlist(cfg$dose_100_ug_ml)))
samples$dose_group <- factor(samples$dose_group, levels = c("vehicle", "low", "medium", "high"))
stopifnot(!anyNA(samples$dose_group), identical(as.integer(table(samples$dose_group)), rep(3L, 4)),
          all(samples$ordered_dose == match(samples$dose_group, levels(samples$dose_group)) - 1L))
design_high <- model.matrix(~ 0 + dose_group, data = samples)
colnames(design_high) <- levels(samples$dose_group)
design_trend <- model.matrix(~ ordered_dose, data = samples)
contrast_high <- limma::makeContrasts(high - vehicle, levels = design_high)
stopifnot(qr(design_high)$rank == ncol(design_high), qr(design_trend)$rank == ncol(design_trend))
counts_table <- read_tsv(paths$counts)
annotation <- read_tsv(paths$annotation)
stopifnot(names(counts_table)[1] == "raw_gene_id", identical(names(counts_table)[-1], samples$sample_id),
          !anyDuplicated(counts_table$raw_gene_id), !anyNA(counts_table$raw_gene_id),
          identical(counts_table$raw_gene_id, annotation$raw_gene_id))
counts <- as.matrix(counts_table[, -1, drop = FALSE])
rownames(counts) <- counts_table$raw_gene_id
rownames(annotation) <- annotation$raw_gene_id
stopifnot(is.numeric(counts), all(is.finite(counts)), all(counts >= 0))
old <- readRDS(paths$object)
stopifnot(identical(samples$sample_id, old$samples$sample_id),
          identical(design_high, old$designs$high), identical(design_trend, old$designs$trend),
          identical(contrast_high, old$designs$contrast_high))
for (nm in names(annotation))
  stopifnot(identical(as.character(annotation[[nm]]), as.character(old$annotation[[nm]])))
y0 <- edgeR::DGEList(counts = counts, genes = annotation)
keep <- edgeR::filterByExpr(y0, design = design_high)
stopifnot(identical(keep, old$keep))
y <- y0[keep, , keep.lib.sizes = FALSE]
y <- edgeR::normLibSizes(y, method = "TMM")
tol <- contract$numeric_input_serialization_tolerance
stopifnot(isTRUE(all.equal(y$counts, old$y$counts, tolerance = tol)),
          isTRUE(all.equal(y$samples, old$y$samples, tolerance = tol)))
ann_keep <- annotation[keep, , drop = FALSE]
logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
human_collapsed <- collapse_expression_by_symbol(logcpm, ann_keep, "human_symbol")
stored_mapping <- read_tsv(paths$mapping)
for (nm in c("raw_gene_id", "gene_symbol", "mapped_symbol"))
  stopifnot(identical(human_collapsed$mapping[[nm]], stored_mapping[[nm]]))
mapping_max_abs <- max(abs(human_collapsed$mapping$mean_logCPM - stored_mapping$mean_logCPM))
stopifnot(isTRUE(all.equal(human_collapsed$mapping$mean_logCPM,
                         stored_mapping$mean_logCPM, tolerance = tol)))
# Synthetic tie/missing-symbol check of the imported source-defined selector.
synthetic <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4), nrow = 4, byrow = TRUE)
syn_annotation <- data.frame(raw_gene_id = c("b", "a", "c", "d"),
                             gene_symbol = c("x", "x", "y", "z"),
                             human_symbol = c("X", "X", NA, "Z"))
selected <- collapse_expression_by_symbol(synthetic, syn_annotation, "human_symbol")
stopifnot(identical(selected$selected_index, c(2L, 4L)))
tie <- collapse_expression_by_symbol(matrix(1, 4, 2), syn_annotation, "human_symbol")
stopifnot(identical(tie$selected_index, c(2L, 4L)))
y_human <- y[human_collapsed$selected_index, , keep.lib.sizes = FALSE]
rownames(y_human$counts) <- rownames(human_collapsed$expression)
y_human$genes <- data.frame(human_symbol = rownames(y_human))
y_human <- edgeR::normLibSizes(y_human, method = "TMM")
reactome_long <- read_tsv(paths$reactome)
stopifnot(all(reactome_long$release == 97L), all(reactome_long$tax_id == 9606L))
sets <- lapply(split(toupper(reactome_long$gene_symbol), reactome_long$term_id), unique)
mapped_sets <- lapply(sets, intersect, y = rownames(y_human))
mapped_sets <- mapped_sets[lengths(mapped_sets) >= contract$min_genes &
                             lengths(mapped_sets) <= contract$max_genes]
index <- limma::ids2indices(mapped_sets, identifiers = rownames(y_human), remove.empty = TRUE)
baselines <- list(high = read_tsv(paths$baseline_high), trend = read_tsv(paths$baseline_trend))
expected_columns <- c("term_id", "NGenes", "PropDown", "PropUp", "Direction", "PValue", "FDR", "PValue.Mixed", "FDR.Mixed")
for (b in baselines) {
  stopifnot(identical(names(b), expected_columns), identical(b$term_id, names(index)),
            !anyDuplicated(b$term_id), identical(as.integer(b$NGenes), as.integer(lengths(index))))
}
input_prov <- read_tsv(paths$input_provenance)
stopifnot(nrow(input_prov) == 1L, input_prov$reference == "RN6",
          unname(tools::md5sum(paths$raw_tar)) == input_prov$source_md5,
          file.info(paths$raw_tar)$size == input_prov$source_size_bytes)
metadata_qc <- read_tsv(paths$metadata_qc)
stopifnot(all(metadata_qc$status %in% c("PASS", "PASS_WITH_LIMITATION")),
          identical(metadata_qc$check_id[metadata_qc$status == "PASS_WITH_LIMITATION"], "exposure_duration"))
dir.create(out, recursive = TRUE)
save_json <- function(x, name) jsonlite::write_json(x, file.path(out, name), auto_unbox = TRUE,
                                                  pretty = TRUE, digits = NA, na = "null")
save_json(list(status = "PASS_CURRENT_NODE_INPUT_AND_METHOD_CHECKS_ONLY", mode = mode,
               historical_producer = "NOT_RECOVERED", sample_n_before = ncol(counts),
               sample_n_after = ncol(counts), samples_removed = 0L,
               input_genes = nrow(counts), filtered_genes = sum(keep),
               human_features = nrow(y_human), eligible_programs = length(index),
               sample_order = samples$sample_id, group_n = as.list(table(samples$dose_group)),
               source_ledger_row_names = rownames(samples),
               design_ranks = c(high = qr(design_high)$rank, trend = qr(design_trend)$rank),
               fractional_count_values = sum(counts != round(counts)), zero_count_values = sum(counts == 0),
               all_zero_input_features = sum(rowSums(counts) == 0),
               library_totals = as.list(colSums(counts)),
               mapping_max_abs_serialization_difference = mapping_max_abs,
               original_metadata_limitations = metadata_qc[metadata_qc$status == "PASS_WITH_LIMITATION", , drop = FALSE],
               selection_synthetic_checks_pass = TRUE,
               critical_errors = list(),
               warnings = c("One batch; culture-library replication is not independent product-lot replication.",
                            "Vehicle DMSO is not the water extraction solvent; labels were retained.",
                            "24 h exposure is locked from the study/ledger, not GEO sample characteristics.",
                            "No holdout training/test split or new efficacy/mechanism inference in this node.",
                            "Complete historical source/runtime/seed binding remains absent."),
               input_identities = before), "INPUT_VALIDATION.json")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
if (mode == "preflight") {
  check_inputs()
  message("PASS: PC12 input/preflight checks; no rotation models run.")
  quit(status = 0L)
}
# Source-defined common categorical dispersion object, shared by the two tests.
y_human <- edgeR::estimateDisp(y_human, design_high, robust = TRUE)
do.call(RNGkind, contract$rng_kind)
run_comparison <- function(label) {
  set.seed(contract$seed)
  if (label == "high") {
    fit <- limma::mroast(y_human, index = index, design = design_high,
                         contrast = contrast_high, nrot = contract$nrot,
                         set.statistic = "mean", gene.weights = NULL, sort = "none")
  } else {
    fit <- limma::mroast(y_human, index = index, design = design_trend,
                         contrast = 2L, nrot = contract$nrot,
                         set.statistic = "mean", gene.weights = NULL, sort = "none")
  }
  data.frame(term_id = rownames(fit), fit, row.names = NULL, check.names = FALSE)
}
comparisons <- list()
for (label in c("high", "trend")) {
  message("Running new equivalent PC12 rotation comparison: ", label)
  fresh <- run_comparison(label)
  filename <- if (label == "high") "Reactome97_mroast1999_high_vs_vehicle.tsv" else "Reactome97_mroast1999_ordered_trend.tsv"
  write_tsv(fresh, file.path(out, "output", filename))
  baseline <- baselines[[label]]
  observed <- read_tsv(file.path(out, "output", filename))
  stopifnot(identical(names(observed), names(baseline)), identical(observed$term_id, baseline$term_id))
  per_column <- lapply(names(baseline), function(nm) {
    a <- baseline[[nm]]; b <- observed[[nm]]
    equal <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    list(column = nm, unequal_n = sum(!equal),
         max_abs_difference = if (is.numeric(a)) max(abs(a - b), na.rm = TRUE) else NULL)
  })
  comparisons[[label]] <- list(rows = nrow(observed), columns = ncol(observed),
                               baseline_sha256 = sha(paths[[paste0("baseline_", label)]]),
                               fresh_sha256 = sha(file.path(out, "output", filename)),
                               bitwise_equal = identical(sha(paths[[paste0("baseline_", label)]]), sha(file.path(out, "output", filename))),
                               all_serialized_values_exact = all(vapply(per_column, function(x) x$unequal_n == 0, logical(1))),
                               per_column = per_column)
}
after <- check_inputs()
stopifnot(identical(before, after))
save_json(list(status = "NEW_EQUIVALENT_RECONSTRUCTION_COMPLETED_COMPARE_BEFORE_ADOPTION",
               historical_producer = "NOT_RECOVERED", started_at_utc = started,
               finished_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
               contract_sha256 = sha(contract_path), runner_sha256 = sha(paths$runner),
               seed = contract$seed, rng_kind = RNGkind(), seed_reset_each_comparison = TRUE,
               nrot = contract$nrot, statistic = "mean", gene_weights = NULL,
               sort = "none", common_dispersion_design = "~0+dose_group; robust=TRUE",
               inputs_unchanged = TRUE, comparison = comparisons,
               original_outputs_replaced = FALSE, original_history_authenticated = FALSE),
          "RECONSTRUCTION_COMPARISON.json")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("Completed both comparisons; old inputs/results untouched.")
