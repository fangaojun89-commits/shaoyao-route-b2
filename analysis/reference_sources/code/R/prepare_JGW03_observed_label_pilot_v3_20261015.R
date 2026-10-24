#!/usr/bin/env Rscript
# Thin, one-cell, explicitly invoked pilot. No original driver or full null run.
# preflight: bind inputs/old endpoint identities/schema; no native model calls.
# observe: one original-label kernel call, persist it before any table bridge.

jgwp_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwp_pin <- function(path) list(bytes = unname(file.info(path)$size),
  sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
jgwp_json <- function(value, path) {
  jgwp_need(!file.exists(path), "New output only")
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE,
    digits = NA, na = "null", null = "null")
}
jgwp_read <- function(path, quote = "") utils::read.delim(path, quote = quote, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
jgwp_json_equal <- function(a, b) identical(
  jsonlite::toJSON(a, auto_unbox = TRUE, digits = NA, na = "null", null = "null"),
  jsonlite::toJSON(b, auto_unbox = TRUE, digits = NA, na = "null", null = "null"))
jgwp_entry_pin <- function() {
  argument <- grep("^--file=", commandArgs(), value = TRUE)
  jgwp_need(length(argument) == 1L, "Explicit Rscript entry file required")
  path <- normalizePath(sub("^--file=", "", argument), winslash = "/", mustWork = TRUE)
  c(list(path = path), jgwp_pin(path))
}
jgwp_sources <- function(root) {
  rel <- c(core = "code/R/JGW03_label_reference_core_20261015.R",
    admission = "code/R/JGW03_input_admission_v2_20261015.R",
    finite_mc = "code/R/JGW03_finite_population_mc_v2_20261015.R",
    native = "code/R/analyze_jgw03_three_cell_lines.R",
    equivalent_roast = "code/R/reconstruct_pc12_rotations_equivalent_20261013.R",
    ledger = "metadata/samples.tsv",
    manifest = "results/jgw03_three_cell_lines/provenance/output_manifest_sha256.tsv")
  expected <- c(
    "e3939f371f02b24cfdfefee0a9cd9e39bdf47c8920caa44f58434a93f7a2974e",
    "63a24a37bac60e1e1c99cd11e9a69c6a11ab8dbb70c44766d06f2e33b87e9c0c",
    "349a5de42e47c41c7211c9f64c0f012d3658c72d840f5b75e4deb8912d52db6b",
    "5df6378ed4dadfb4050d91e252a63af4f44d2be80d3070d889d1e0e007796775",
    "ea666b1e5aac06eea431e5f376148755bb47ae02c3b467ef56ba71e9b79d3577",
    "367894147ee329498e6491c9397d1fe146ca8be35b696110964c1365da3bdf66",
    "ace453263d4c7c524c34cb9fce92d5f5f7b9e8449f59b7ecb1fe3f51094c8074")
  answer <- lapply(seq_along(rel), function(i) {
    path <- file.path(root, rel[[i]])
    jgwp_need(file.exists(path) && identical(jgwp_pin(path)$sha256, expected[[i]]),
      paste("Source/metadata pin differs:", names(rel)[[i]]))
    c(list(role = names(rel)[[i]], relative_path = rel[[i]]), jgwp_pin(path))
  })
  names(answer) <- names(rel)
  answer
}
jgwp_runtime <- function() {
  packages <- c("edgeR", "limma", "DESeq2", "decoupleR", "fgsea", "BiocParallel",
    "digest", "jsonlite", "SummarizedExperiment", "S4Vectors")
  versions <- setNames(vapply(packages, function(x) {
    jgwp_need(requireNamespace(x, quietly = TRUE), paste("Existing package unavailable:", x))
    as.character(utils::packageVersion(x))
  }, character(1)), packages)
  jgwp_need(as.character(getRversion()) == "4.6.1", "Expected existing R 4.6.1")
  # Native fgsea calls use the registered default when BPPARAM is omitted.
  # This serial equivalent path is explicit, not historical-backend authentication.
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  jgwp_need(inherits(BiocParallel::bpparam(), "SerialParam") &&
    BiocParallel::bpnworkers(BiocParallel::bpparam()) == 1L, "Serial backend required")
  list(R = as.character(getRversion()), packages = as.list(versions),
    platform = R.version$platform, BLAS = unname(extSoftVersion()["BLAS"]),
    backend = class(BiocParallel::bpparam())[1L], workers = 1L,
    seed = 20260817L, RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"),
    historical_backend_authenticated = FALSE)
}
jgwp_expected_tables <- function(accession) {
  genes <- c(paste0("edgeR_QL_", c("high_vs_vehicle", "ordered_trend", "logdose_trend"),
                    "_all_genes.tsv.gz"),
    paste0("limma_voom_", c("high_vs_vehicle", "ordered_trend"), "_all_genes.tsv.gz"),
    paste0("DESeq2_rounded_", c("high_vs_vehicle", "ordered_trend"), "_all_genes.tsv.gz"))
  programs <- unlist(lapply(c("Reactome97", "QuickGO_SynGO_sensitivity"), function(x)
    paste0(x, c("_camera_high_vs_vehicle.tsv", "_camera_ordered_trend.tsv",
      "_camera_logdose_trend.tsv", "_score_high.tsv", "_score_trend.tsv",
      "_score_logtrend.tsv", "_fgsea_high_vs_vehicle.tsv", "_fgsea_ordered_trend.tsv",
      "_sample_scores.tsv", "_high_trend_dual.tsv"))), use.names = FALSE)
  if (accession == "GSE289929") programs <- c(programs,
    "Reactome97_mroast1999_high_vs_vehicle.tsv", "Reactome97_mroast1999_ordered_trend.tsv")
  activities <- unlist(lapply(c("PROGENy", "DoRothEA"), function(x)
    paste0(x, c("_high.tsv", "_trend.tsv", "_logtrend.tsv", "_sample_activity.tsv"))),
    use.names = FALSE)
  base <- file.path("results/jgw03_three_cell_lines", accession)
  c(file.path(base, "genes", genes), file.path(base, "programs", programs),
    file.path(base, "activities", activities))
}
jgwp_bound_records <- function(root, accession, admission) {
  own <- admission$jgw03_v2_expected_files(accession)
  own$size_bytes <- own$bytes
  own <- own[, c("relative_path", "size_bytes", "sha256")]
  manifest <- jgwp_read(file.path(root,
    "results/jgw03_three_cell_lines/provenance/output_manifest_sha256.tsv"), quote = '"')
  expected <- c(jgwp_expected_tables(accession), file.path(
    "data/processed/jgw03_three_cell_lines", accession, "analysis_objects.rds"))
  jgwp_need(!anyDuplicated(manifest$relative_path) &&
    all(expected %in% manifest$relative_path), "All retained endpoints and object must be manifest-bound")
  result <- rbind(own, manifest[match(expected, manifest$relative_path), names(own)])
  result <- result[!duplicated(result$relative_path), , drop = FALSE]
  rownames(result) <- NULL
  for (i in seq_len(nrow(result))) admission$jgw03_v2_check_file(
    file.path(root, result$relative_path[i]), result$size_bytes[i], result$sha256[i])
  result
}
jgwp_table_check <- function(fresh, old, key, endpoint) {
  # Exclude ONLY the three explicitly unvalidated historical QLF interval fields.
  excluded <- intersect(c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high"), names(old))
  retained_names <- names(old)[!names(old) %in% excluded]
  old <- old[, !names(old) %in% excluded, drop = FALSE]
  names(old) <- retained_names  # Preserve the original duplicate annotation columns.
  shape <- is.data.frame(fresh) && is.data.frame(old) && identical(names(fresh), names(old)) &&
    nrow(fresh) == nrow(old) && key %in% names(old) &&
    !anyNA(old[[key]]) && !anyDuplicated(old[[key]]) &&
    identical(as.character(fresh[[key]]), as.character(old[[key]]))
  if (!shape) return(data.frame(endpoint = endpoint, column = "SCHEMA_KEYS", n = nrow(old),
    pass = FALSE, compared_finite = 0L, unavailable = 0L, max_abs = NA_real_, max_log10P = NA_real_))
  checks <- lapply(seq_along(old), function(i) {
    a <- fresh[[i]]; b <- old[[i]]; nm <- names(old)[i]
    missing_equal <- identical(is.na(a), is.na(b)) && identical(is.nan(a), is.nan(b))
    numeric_pair <- is.numeric(a) && is.numeric(b)
    finite <- if (numeric_pair) is.finite(a) & is.finite(b) else rep(FALSE, length(b))
    probability <- nm %in% c("PValue", "FDR", "P.Value", "adj.P.Val", "pvalue", "padj",
      "pval", "PValue.Mixed", "FDR.Mixed", "FDR_high", "FDR_trend") ||
      grepl("(^|_)(PValue|FDR)$", nm)
    delta <- if (any(finite)) max(abs(a[finite] - b[finite])) else NA_real_
    log_delta <- NA_real_
    if (numeric_pair) {
      available <- !is.na(b)
      ok <- missing_equal && identical(is.infinite(a), is.infinite(b)) &&
        all(a[available & !finite] == b[available & !finite]) &&
        all(abs(a[finite] - b[finite]) <= 1e-10 + 1e-8 * abs(b[finite]))
      # No transport tolerance may silently change zeros or reported directions.
      if (nm %in% c("logFC", "effect", "t", "stat", "log2FoldChange", "ES", "NES",
                    "F", "CI_low", "CI_high") || grepl("_(effect|t|CI_low|CI_high)$", nm))
        ok <- ok && all(sign(a[finite]) == sign(b[finite]))
      if (probability) {
        positive <- finite & a > 0 & b > 0
        log_delta <- if (any(positive)) max(abs(log10(a[positive]) - log10(b[positive]))) else NA_real_
        ok <- ok && all(a[finite] >= 0 & a[finite] <= 1 & b[finite] >= 0 & b[finite] <= 1) &&
          identical(a == 0, b == 0) && all((a[finite] <= .05) == (b[finite] <= .05)) &&
          (!any(positive) || log_delta <= 1e-6)
      }
    } else ok <- missing_equal && identical(as.character(a[!is.na(a)]), as.character(b[!is.na(b)]))
    data.frame(endpoint = endpoint, column = paste0(nm, "[", i, "]"), n = length(b),
      pass = ok, compared_finite = sum(finite), unavailable = sum(is.na(b)),
      max_abs = delta, max_log10P = log_delta)
  })
  do.call(rbind, checks)
}
jgwp_fresh_tables <- function(fixed, result) {
  out <- result$native_exports
  stems <- c(high = "high_vs_vehicle", trend = "ordered_trend", logtrend = "logdose_trend")
  for (label in names(stems)) {
    value <- result$gene_edgeR[[label]]
    degrees <- result$native_edgeR_tests[[label]]$df.total
    if (length(degrees) == 1L) degrees <- rep(degrees, nrow(value))
    jgwp_need(length(degrees) == nrow(value), "Native QL df length differs")
    value$ql_df_total <- degrees
    out[[paste0("edgeR_QL_", stems[[label]], "_all_genes.tsv.gz")]] <- value
  }
  for (label in c("high", "trend")) {
    out[[paste0("limma_voom_", stems[[label]], "_all_genes.tsv.gz")]] <- result$gene_voom[[label]]$value
    out[[paste0("DESeq2_rounded_", stems[[label]], "_all_genes.tsv.gz")]] <- result$gene_DESeq2[[label]]$value
  }
  for (role in c("PROGENy", "DoRothEA")) {
    family <- if (role == "PROGENy") "PROGENy_top500_ULM" else "DoRothEA_AB_ULM"
    activity <- if (role == "PROGENy") fixed$progeny_activity else fixed$dorothea_activity
    out[[paste0(role, "_sample_activity.tsv")]] <- data.frame(feature = rownames(activity),
      activity, check.names = FALSE)
    for (label in c("high", "trend", "logtrend")) out[[paste0(role, "_", label, ".tsv")]] <-
      result$activities[[family]][[label]]
  }
  out
}
jgwp_roast_check <- function(result) {
  roast <- result$programs$Reactome97$roast
  if (is.null(roast)) return(list(applicable = FALSE))
  per_table <- lapply(roast, function(x) {
    check <- function(p, q) all(is.finite(p) & is.finite(q)) &&
      all(p >= 0 & p <= 1 & q >= 0 & q <= 1) &&
      all(abs(q - pmax(stats::p.adjust(p - .5 / 2000, "BH"), p)) <= 1e-12)
    list(n = nrow(x), nrot = 1999L, directional_midP_FDR = check(x$PValue, x$FDR),
      mixed_midP_FDR = check(x$PValue.Mixed, x$FDR.Mixed),
      raw_P_BH_is_not_native_midP_FDR = TRUE)
  })
  list(applicable = TRUE, tables = per_table,
    pass = all(vapply(per_table, function(x) x$directional_midP_FDR && x$mixed_midP_FDR, logical(1))))
}
jgwp_main <- function(args) {
  jgwp_need(length(args) %in% c(4L, 6L), "ROOT OUT ACCESSION preflight|observe [PREFLIGHT SHA] required")
  root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
  out <- gsub("\\\\", "/", args[[2]])
  accession <- args[[3]]; mode <- args[[4]]
  jgwp_need(mode %in% c("preflight", "observe") &&
    accession %in% c("GSE289929", "GSE227494", "GSE299063") && !file.exists(out) &&
    dir.exists(dirname(out)) &&
    startsWith(out, "R:/Shaoyao_storage/cache/JGW03_label_reference_20261015r01/"),
    "Explicit one-cell mode and fresh R-run leaf required")
  jgwp_need(startsWith(normalizePath(tempdir(), winslash = "/"), "R:/Shaoyao_storage/"),
    "TMPDIR must be on R before launching R")
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  runtime <- jgwp_runtime()
  entry <- jgwp_entry_pin()
  source_before <- jgwp_sources(root)
  admission <- new.env(parent = asNamespace("stats"))
  sys.source(file.path(root, source_before$admission$relative_path), envir = admission)
  kernel <- admission$jgw03_v2_bind_kernel(file.path(root, source_before$core$relative_path))
  records <- jgwp_bound_records(root, accession, admission)
  old_path <- file.path(root, "data/processed/jgw03_three_cell_lines", accession, "analysis_objects.rds")
  old <- readRDS(old_path)
  jgwp_need(is.list(old) && all(c("samples", "annotation", "keep", "y", "designs", "tests",
    "reactome", "syngo", "progeny", "dorothea") %in% names(old)), "Retained native schema differs")
  samples <- jgwp_read(file.path(root, "data/processed/jgw03_three_cell_lines", accession, "locked_samples.tsv"))
  admission$jgw03_v2_validate_metadata(samples, accession)
  ledger <- jgwp_read(file.path(root, source_before$ledger$relative_path))
  ledger <- ledger[ledger$dataset_id == accession & ledger$include_primary == "TRUE", , drop = FALSE]
  ledger <- ledger[order(as.integer(ledger$sample_order)), , drop = FALSE]
  jgwp_need(identical(samples$sample_id, ledger$sample_id) &&
    identical(samples$sample_id, old$samples$sample_id), "Independent ordered library metadata differs")
  # Restore cosmetic labels from independent ledger, not fitted coefficients.
  rownames(samples) <- rownames(ledger)
  for (field in intersect(names(ledger), names(samples))) jgwp_need(
    identical(is.na(ledger[[field]]) | as.character(ledger[[field]]) == "",
              is.na(samples[[field]]) | as.character(samples[[field]]) == "") &&
    all(as.character(ledger[[field]])[!is.na(ledger[[field]]) & as.character(ledger[[field]]) != ""] ==
        as.character(samples[[field]])[!is.na(samples[[field]]) & as.character(samples[[field]]) != ""]),
    "Original sample ledger literals differ")
  observed <- setNames(as.numeric(samples$ordered_dose), samples$sample_id)
  designs <- kernel$jgw03_designs(kernel$jgw03_relabel_model_metadata(samples, observed, accession))
  design_equal <- vapply(names(designs), function(nm) identical(designs[[nm]], old$designs[[nm]]), logical(1))
  jgwp_need(all(design_equal), "Original-label exact design metadata differs")
  registry <- kernel$jgw03_family_registry()
  registry <- registry[registry$accession == accession & registry$space == "four_group_MC", ]
  schema <- list(root_names = names(old), retained_gene_n = sum(old$keep),
    old_count_dim = dim(old$y$counts), old_designs = lapply(old$designs, dim),
    old_program_names = list(reactome = names(old$reactome), syngo = names(old$syngo)),
    ordered_sample_ids = samples$sample_id, original_ledger_row_names = rownames(samples))
  if (mode == "observe") {
    jgwp_need(length(args) == 6L, "Observe needs exact preflight receipt and SHA")
    preflight_path <- normalizePath(args[[5]], winslash = "/", mustWork = TRUE)
    jgwp_need(startsWith(preflight_path, "R:/Shaoyao_storage/cache/JGW03_label_reference_20261015r01/") &&
      identical(jgwp_pin(preflight_path)$sha256, args[[6]]), "Preflight receipt identity differs")
    pre <- jsonlite::read_json(preflight_path, simplifyVector = FALSE)
    jgwp_need(identical(pre$status, "PREFLIGHT_METADATA_IDENTITIES_NO_MODEL_CALLS") &&
      identical(pre$accession, accession) && jgwp_json_equal(pre$runtime, runtime) &&
      jgwp_json_equal(pre$sources, source_before) && jgwp_json_equal(pre$entry, entry),
      "Preflight runtime/source/entry/accession differs")
  }
  dir.create(out)
  utils::write.table(registry, file.path(out, "ENDPOINT_FAMILIES.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  common <- list(accession = accession, mode = mode, started_utc = started,
    sources = source_before, entry = entry, runtime = runtime, protected_inputs = records,
    schema = schema, observed_labels = as.list(observed), native_family_n = nrow(registry),
    real_null_allocations_fitted = 0L, historical_stochastic_backend_authenticated = FALSE,
    candidate_gate_changed = FALSE, full_reference_admitted = FALSE)
  capture.output(utils::sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  if (mode == "preflight") {
    jgwp_need(identical(records, jgwp_bound_records(root, accession, admission)) &&
      identical(source_before, jgwp_sources(root)), "Protected inputs changed during preflight")
    common$status <- "PREFLIGHT_METADATA_IDENTITIES_NO_MODEL_CALLS"
    common$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    jgwp_json(common, file.path(out, "PREFLIGHT.json"))
    return(0L)
  }
  # Everything before this point is metadata/hash/schema. No old native fit is reused.
  do.call(RNGkind, as.list(runtime$RNGkind)); set.seed(runtime$seed)
  stage <- proc.time()[["elapsed"]]
  bound <- admission$jgw03_v2_read_bound_inputs(root, accession)
  bound$samples <- samples
  rownames(bound$annotation) <- bound$annotation$raw_gene_id
  fixed <- kernel$jgw03_build_fixed_layer(bound$counts, bound$annotation, bound$samples,
    accession, bound$resources, file.path(root, source_before$native$relative_path), root)
  fixed_seconds <- proc.time()[["elapsed"]] - stage
  saveRDS(fixed, file.path(out, "FIXED_LAYER.rds"), compress = TRUE, version = 3L)
  fixed_bridge <- list(keep = identical(fixed$keep, old$keep),
    retained_gene_order = identical(rownames(fixed$y$counts), rownames(old$y$counts)),
    sample_order = identical(colnames(fixed$y$counts), colnames(old$y$counts)),
    retained_counts = identical(fixed$y$counts, old$y$counts),
    library_sizes = isTRUE(all.equal(fixed$y$samples$lib.size, old$y$samples$lib.size,
      tolerance = 1e-12, check.attributes = FALSE)),
    TMM_factors = isTRUE(all.equal(fixed$y$samples$norm.factors, old$y$samples$norm.factors,
      tolerance = 1e-12, check.attributes = FALSE)))
  jgwp_json(fixed_bridge, file.path(out, "FIXED_LAYER_BRIDGE.json"))
  jgwp_need(all(unlist(fixed_bridge)), "Fixed observed layer differs; no new native gene fit started")
  rm(old); invisible(gc())
  stage <- proc.time()[["elapsed"]]
  result <- kernel$jgw03_fit_one_allocation(fixed, observed,
    file.path(root, source_before$native$relative_path))
  fit_seconds <- proc.time()[["elapsed"]] - stage
  # Persist before reading/comparing any historical endpoint: transport/bridge
  # failure must not erase the completed original-label fit or trigger a refit.
  result_file <- file.path(out, "OBSERVED_NATIVE.rds")
  saveRDS(result, result_file, compress = TRUE, version = 3L)
  jgwp_json(c(common, list(status = "ORIGINAL_LABEL_NATIVE_FIT_SAVED_BRIDGE_PENDING",
    result = jgwp_pin(result_file), fixed_layer = jgwp_pin(file.path(out, "FIXED_LAYER.rds")),
    fixed_seconds = fixed_seconds, fit_seconds = fit_seconds)), file.path(out, "FIT_SAVED.json"))
  tables <- jgwp_fresh_tables(fixed, result)
  expected <- jgwp_expected_tables(accession)
  jgwp_need(setequal(names(tables), basename(expected)) && !anyDuplicated(names(tables)),
    "Full original endpoint plus auxiliary table coverage differs")
  export_dir <- file.path(out, "native_exports")
  dir.create(export_dir)
  for (name in names(tables)) {
    path <- file.path(export_dir, name)
    connection <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
    tryCatch(utils::write.table(tables[[name]], connection, sep = "\t", quote = FALSE,
      row.names = FALSE, na = "NA"), finally = close(connection))
  }
  checks <- lapply(expected, function(rel) {
    old_table <- jgwp_read(file.path(root, rel)); fresh <- tables[[basename(rel)]]
    key <- if ("raw_gene_id" %in% names(old_table)) "raw_gene_id" else
      if ("term_id" %in% names(old_table)) "term_id" else "feature"
    jgwp_table_check(fresh, old_table, key, basename(rel))
  })
  checks <- do.call(rbind, checks)
  utils::write.table(checks, file.path(out, "NATIVE_FIELD_BRIDGE.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA")
  roast <- jgwp_roast_check(result)
  # Native stochastic uncertainty is reported as saved by fgsea. No seed tuning,
  # alternative backend rerun, or imputation of unavailable log2err is performed.
  mc_diagnostics <- lapply(expected[grepl("_fgsea_", expected)], function(rel) {
    old_table <- jgwp_read(file.path(root, rel)); fresh <- tables[[basename(rel)]]
    data.frame(endpoint = basename(rel), term_id = old_table$term_id,
      old_P = old_table$pval, current_P = fresh$pval,
      old_log2err = old_table$log2err, current_log2err = fresh$log2err,
      old_NES = old_table$NES, current_NES = fresh$NES)
  })
  utils::write.table(do.call(rbind, mc_diagnostics), file.path(out, "FGSEA_NATIVE_MC_DIAGNOSTICS.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  after <- jgwp_bound_records(root, accession, admission)
  jgwp_need(identical(records, after) && identical(source_before, jgwp_sources(root)) &&
    identical(entry, jgwp_entry_pin()),
    "Protected input/source changed during original-label pilot")
  loaded <- readRDS(result_file)
  jgwp_need(identical(loaded$labels, observed) &&
    identical(digest::digest(jgwp_fresh_tables(fixed, loaded), "sha256", serializeVersion = 3L),
              digest::digest(tables, "sha256", serializeVersion = 3L)), "Saved table transport differs")
  passed <- all(checks$pass) && (!isTRUE(roast$applicable) || isTRUE(roast$pass))
  common$status <- if (passed) "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_ONLY" else
    "OBSERVED_FIT_RETAINED_BRIDGE_REVIEW_REQUIRED"
  common$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  common$fixed_seconds <- fixed_seconds; common$fit_seconds <- fit_seconds
  common$fixed_layer_bridge <- fixed_bridge
  common$RNG_after <- list(kind = RNGkind(), state_sha256 = digest::digest(.Random.seed, "sha256"))
  common$bridge <- list(table_n = length(expected), column_checks = nrow(checks),
    all_fields_pass = passed, failed_columns = sum(!checks$pass),
    floating_bound = "abs 1e-10 + rel 1e-8; positive P log10 delta <=1e-6",
    NA_NaN_zero_direction_and_native_0_05_masks = "exact",
    excluded_unvalidated_old_fields = c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high"))
  common$mroast_native_FDR <- roast
  common$DESeq2_diagnostics <- lapply(result$gene_DESeq2, function(x) list(
    pvalue_NA = sum(is.na(x$value$pvalue)), padj_NA = sum(is.na(x$value$padj)),
    native_default_independentFiltering_and_Cooks_retained = TRUE,
    minReplicatesForReplace = "Inf", fitType = attr(DESeq2::dispersionFunction(x$native_locals$dds), "fitType")))
  common$result <- jgwp_pin(result_file)
  common$protected_inputs_after <- after
  common$scientific_interpretation <- "Same-library exploratory conditional reference preparation, not biological null/FDR/candidate validation"
  jgwp_json(common, file.path(out, "OBSERVED_COMPLETE.json"))
  if (passed) 0L else 2L
}

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({ library(digest); library(jsonlite) })
  options(stringsAsFactors = FALSE, width = 140)
  exit_code <- tryCatch(jgwp_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    # Do not save/transcribe the raw condition. Successful FIT_SAVED remains if present.
    message("JGW03_PILOT_STOP: inspect retained new-node artifacts; no automatic refit")
    1L
  })
  quit(status = exit_code)
}




