# One missing exploratory cell, with exact historical kernels and retained-input checks.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, args[[1L]] %in% c("prepare", "run"))
mode <- args[[1L]]
root <- normalizePath(args[[2L]], mustWork = TRUE, winslash = "/")
out <- args[[3L]]
options(stringsAsFactors = FALSE, width = 140)
stopifnot(as.character(getRversion()) == "4.6.1")
Sys.setlocale("LC_COLLATE", "C")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20260817L)
suppressPackageStartupMessages({library(data.table); library(limma); library(matrixStats); library(decoupleR)})
sha <- function(p) tolower(unname(tools::sha256sum(p)))
source_path <- file.path(root, "code/R/analyze_GSE16710_spasticity_anchor.R")
stopifnot(sha(source_path) == "37a0b7e1dd250b618ef1b507b6bd72a4ae0bb47193ac2271953a6f4bb026feaf")
functions <- c("safe_sign", "zscore_rows", "make_design", "fit_limma_matrix", "extract_fit_long",
               "build_reactome_mapping", "program_mean_activity", "run_camera_contrasts",
               "build_ulm_networks", "run_ulm_activity")
parsed <- parse(source_path)
selected <- list()
for (expr in parsed) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% functions &&
      is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
    selected[[as.character(expr[[2L]])]] <- expr
  }
}
stopifnot(setequal(names(selected), functions), length(selected) == length(functions))
for (expr in selected) eval(expr, envir = .GlobalEnv) # named definitions only, no analyzer top-level execution
processed <- file.path(root, "data/processed/GSE16710")
program <- file.path(root, "results/GSE16710/program_level")
comparison_names <- c("reactome_camera_original_annotation_sensitivity.tsv", "reactome_activity_original_annotation_sensitivity.tsv",
                      "progeny_original_annotation_sensitivity.tsv", "dorothea_original_annotation_sensitivity.tsv")
relative <- c("code/R/analyze_GSE16710_spasticity_anchor.R", "config/contrasts.yaml", "renv.lock",
              "protocol/statistical_analysis_plan.md", "protocol/GSE16710_ORIGINAL_ANNOTATION_ALL27_CROSSED_SENSITIVITY_20261014.md",
              "environment/run_r.ps1", paste0("data/processed/GSE16710/", c("rma_probe_expression.rds",
              "rma_gene_expression_current_bioconductor.rds", "rma_gene_expression_original_GPL1355_2016.rds",
              "probe_annotation_current_bioconductor.tsv", "probe_annotation_original_GPL1355_2016.tsv",
              "sample_ledger_verified.tsv", "qc_inclusion_decision.tsv")),
              "data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds",
              "results/GSE16710/program_level/reactome_mapping_coverage_current.tsv",
              paste0("results/GSE16710/program_level/", comparison_names),
              "results/GSE16710/program_level/S1_GATE.tsv")
paths <- file.path(root, relative)
input_identity <- data.table(relative_path = relative, bytes = file.info(paths)$size, sha256 = vapply(paths, sha, character(1)))
script_arg <- commandArgs(trailingOnly = FALSE)
helper <- normalizePath(sub("^--file=", "", script_arg[grep("^--file=", script_arg)]), mustWork = TRUE, winslash = "/")
ledger <- fread(file.path(processed, "sample_ledger_verified.tsv"))
qc <- fread(file.path(processed, "qc_inclusion_decision.tsv"))
stopifnot(nrow(ledger) == 27L, nrow(qc) == 27L, !anyDuplicated(ledger$sample_id), !anyDuplicated(qc$sample_id),
          setequal(ledger$sample_id, qc$sample_id), all(toupper(as.character(qc$outcome_information_used)) == "FALSE"))
ledger <- merge(ledger, qc, by = "sample_id", sort = FALSE)
ledger <- ledger[order(as.integer(sample_order))]
stopifnot(identical(ledger[primary_analysis_decision != "retain", sample_id], "GSM418682"))
primary <- ledger[primary_analysis_decision == "retain"]
design26 <- make_design(primary); design27 <- make_design(ledger)
needed <- c("limma", "data.table", "matrixStats", "decoupleR", "jsonlite", "statmod")
lock <- jsonlite::fromJSON(file.path(root, "renv.lock"), simplifyVector = FALSE)
versions <- vapply(needed, function(p) as.character(packageVersion(p)), character(1))
stopifnot(all(vapply(needed, function(p) versions[[p]] == lock$Packages[[p]]$Version, logical(1))))
registry_path <- file.path(out, "REGISTRY.json")
if (mode == "prepare") {
  if (file.exists(out) || !dir.exists(dirname(out))) stop("Use a fresh R output directory with existing parent")
  dir.create(out)
  registry <- list(status = "POST_OUTCOME_PLAN_FROZEN_BEFORE_NEW_CELL", utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    inputs = input_identity, source_sha256 = sha(source_path), helper_sha256 = sha(helper), helper_path = helper,
    sample_order_all27 = ledger$sample_id, sample_order_primary26 = primary$sample_id,
    groups_all27 = ledger[, .N, by = .(experimental_arm, time_day)],
    groups_primary26 = primary[, .N, by = .(experimental_arm, time_day)],
    design27 = unname(design27$design), design26 = unname(design26$design), contrasts = unname(design27$contrasts),
    package_versions = as.list(versions), R_version = as.character(getRversion()), locale = Sys.getlocale(),
    absolute_tolerance = 1e-8, relative_tolerance = 1e-7, seed = 20260817L, RNG_kind = RNGkind(),
    primary_changed = FALSE, new_cell = "all27_original_GPL1355", S1_gate_recomputed = FALSE,
    reference_refresh_or_raw_preprocessing = FALSE, new_GSE277064_fits = FALSE)
  jsonlite::write_json(registry, registry_path, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", digits = NA)
  cat("PREPARED ", registry_path, " SHA256=", sha(registry_path), "\n", sep = "")
  quit(status = 0L)
}
stopifnot(file.exists(registry_path), !file.exists(file.path(out, "RUN_STARTED.json")))
registry <- jsonlite::fromJSON(registry_path)
stopifnot(identical(registry$inputs$sha256, input_identity$sha256), registry$helper_sha256 == sha(helper),
          identical(registry$sample_order_all27, ledger$sample_id), identical(registry$sample_order_primary26, primary$sample_id))
jsonlite::write_json(list(status = "STARTED", utc = format(Sys.time(), tz = "UTC", usetz = TRUE), registry_sha256 = sha(registry_path),
                         invocation = commandArgs(), helper_sha256 = sha(helper)), file.path(out, "RUN_STARTED.json"), pretty = TRUE, auto_unbox = TRUE)
probe <- readRDS(file.path(processed, "rma_probe_expression.rds"))
current <- readRDS(file.path(processed, "rma_gene_expression_current_bioconductor.rds"))
original <- readRDS(file.path(processed, "rma_gene_expression_original_GPL1355_2016.rds"))
stopifnot(identical(colnames(probe), ledger$sample_id), identical(colnames(current), ledger$sample_id),
          identical(colnames(original), ledger$sample_id), all(is.finite(probe)), all(is.finite(original)))
frozen <- readRDS(file.path(root, "data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds"))
current_map <- fread(file.path(processed, "probe_annotation_current_bioconductor.tsv"))
original_map <- fread(file.path(processed, "probe_annotation_original_GPL1355_2016.tsv"))
for (map in c("current_map", "original_map")) {
  x <- get(map); x[, probe_index := match(probe_id, rownames(probe))]; assign(map, x[!is.na(probe_index)])
}
current_mapping <- build_reactome_mapping(frozen, rownames(current), current_map)
eligible <- current_mapping$coverage[eligible == TRUE, term_id]
retained_eligible <- fread(file.path(program, "reactome_mapping_coverage_current.tsv"))[eligible == TRUE, term_id]
stopifnot(identical(eligible, retained_eligible), length(eligible) == 1656L)
mapping <- build_reactome_mapping(frozen, rownames(original), original_map, fixed_terms = eligible)
index <- mapping$probe_index[intersect(eligible, names(mapping$probe_index))]
pathways <- mapping$pathways[intersect(eligible, names(mapping$pathways))]
networks <- build_ulm_networks(frozen, rownames(original))
calculate <- function(samples, design, camera_label) {
  x <- original[, samples, drop = FALSE]
  list(CAMERA = run_camera_contrasts(probe[, samples, drop = FALSE], index, design, camera_label, "spasticity_progression_interaction"),
       Reactome_mean_z = extract_fit_long(fit_limma_matrix(program_mean_activity(x, pathways), design, TRUE)),
       PROGENy_ULM = extract_fit_long(fit_limma_matrix(run_ulm_activity(x, networks$progeny), design, TRUE)),
       DoRothEA_ULM = extract_fit_long(fit_limma_matrix(run_ulm_activity(x, networks$dorothea), design, TRUE)))
}
cat("Reconstructing four original26 comparison tables before new-cell fitting\n")
reconstructed <- calculate(primary$sample_id, design26, "Reactome_97_original_GPL1355")
validation <- list()
for (i in seq_along(reconstructed)) {
  actual <- copy(reconstructed[[i]]); expected <- fread(file.path(program, comparison_names[[i]]))
  stopifnot(identical(names(actual), names(expected)), !anyDuplicated(actual[, .(feature_id, contrast)]),
            !anyDuplicated(expected[, .(feature_id, contrast)]))
  setorder(actual, feature_id, contrast); setorder(expected, feature_id, contrast)
  stopifnot(nrow(actual) == nrow(expected), identical(actual$feature_id, expected$feature_id), identical(actual$contrast, expected$contrast))
  for (field in names(actual)) {
    a <- actual[[field]]; b <- expected[[field]]
    if (is.numeric(a) && is.numeric(b)) {
      ok <- identical(is.na(a), is.na(b)) && all(abs(a - b) <= registry$absolute_tolerance + registry$relative_tolerance * abs(b), na.rm = TRUE)
      delta <- if (any(is.finite(a) & is.finite(b))) max(abs(a - b), na.rm = TRUE) else NA_real_
    } else {ok <- identical(as.character(a), as.character(b)); delta <- NA_real_}
    validation[[length(validation) + 1L]] <- data.table(layer = names(reconstructed)[[i]], field = field, rows = nrow(actual), pass = ok, maximum_absolute_difference = delta)
  }
}
checks <- rbindlist(validation)
fwrite(checks, file.path(out, "ORIGINAL26_COMPATIBILITY.tsv"), sep = "\t", na = "NA")
if (!all(checks$pass)) stop("Original26 compatibility mismatch; new cell not fitted")
cat("Original26 compatibility passed; fitting the missing original/all27 cell\n")
new <- calculate(ledger$sample_id, design27, "Reactome_97_original_GPL1355_all27")
new_files <- paste0(c("reactome_camera", "reactome_activity", "progeny", "dorothea"), "_original_all27_sensitivity.tsv")
summary <- list()
for (i in seq_along(new)) {
  x <- new[[i]]
  stopifnot(!anyDuplicated(x[, .(feature_id, contrast)]), all(is.finite(x$p_value)), all(is.finite(x$q_value_BH)))
  fwrite(x, file.path(out, new_files[[i]]), sep = "\t", quote = FALSE, na = "NA")
  summary[[i]] <- x[, .(tested = .N, q_le_005 = sum(q_value_BH <= 0.05), minimum_q = min(q_value_BH)), by = contrast][, layer := names(new)[[i]]]
}
family <- rbindlist(summary)
fwrite(family, file.path(out, "FAMILY_SUMMARY.tsv"), sep = "\t")
fwrite(mapping$coverage, file.path(out, "ORIGINAL_MAPPING_COVERAGE.tsv"), sep = "\t", na = "NA")
fwrite(as.data.table(design27$design, keep.rownames = "sample_id"), file.path(out, "DESIGN_ALL27.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "SESSION_INFO.txt"))
stopifnot(identical(input_identity$sha256, unname(vapply(paths, sha, character(1)))), registry$helper_sha256 == sha(helper))
output_files <- file.path(out, c(new_files, "FAMILY_SUMMARY.tsv", "ORIGINAL26_COMPATIBILITY.tsv", "ORIGINAL_MAPPING_COVERAGE.tsv", "DESIGN_ALL27.tsv", "SESSION_INFO.txt"))
receipt <- list(status = "COMPLETE_EXPLORATORY_ORIGINAL_ALL27_CELL", utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  registry_sha256 = sha(registry_path), helper_sha256 = sha(helper), inputs_unchanged = TRUE,
  original26_comparisons_passed = TRUE, comparison_field_rows = nrow(checks), native_family_summary = family,
  output_identities = data.table(path = output_files, bytes = file.info(output_files)$size, sha256 = vapply(output_files, sha, character(1))),
  source_derived_matrix_reproduction_only = TRUE, raw_RMA_recomputed = FALSE, new_cell_model_fitted = TRUE,
  primary_or_gate_changed = FALSE, independent_validation = FALSE, whole_multiverse_or_goal_complete = FALSE)
jsonlite::write_json(receipt, file.path(out, "COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", digits = NA)
cat("COMPLETE missing crossed cell; existing primary/gates unchanged\n")
print(family)


