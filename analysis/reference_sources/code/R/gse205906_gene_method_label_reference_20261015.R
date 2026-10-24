#!/usr/bin/env Rscript
# New method-specific conditional references; no historical model driver is run.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
library(data.table)
library(digest)
library(jsonlite)
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
pin <- function(path) list(path = normalizePath(path, mustWork = TRUE),
  bytes = unname(file.info(path)$size), sha256 = digest(file = path, algo = "sha256", serialize = FALSE))
verify <- function(value) {
  actual <- pin(value$path)
  need(actual$bytes == value$bytes && actual$sha256 == value$sha256, "Frozen file identity")
  invisible(actual)
}
jwrite <- function(value, path) {
  need(!file.exists(path), paste("Preserve existing output", basename(path)))
  write_json(value, path, pretty = TRUE, auto_unbox = TRUE, digits = 17, na = "null", null = "null")
}
twrite <- function(value, path) {
  need(!file.exists(path), paste("Preserve existing output", basename(path)))
  fwrite(value, path, sep = "\t", quote = TRUE, na = "NA")
}
groups <- c("sham", "SNI", "PLM")
contrasts <- cbind(local_disease = c(-1, 1, 0), intervention = c(0, -1, 1))
rownames(contrasts) <- groups
key_columns <- c("collection", "method", "contrast", "feature")
numeric_columns <- c("effect", "statistic", "p_value", "q_value")
versions <- c(edgeR = "4.10.1", limma = "3.68.4", DESeq2 = "1.52.0")
peak_memory_kb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line))
}

import_canonical <- function(path) {
  need(pin(path)$sha256 == "18eaab06d0bb980dae3373279c83d0e4e32f94dbee346ef4e2cc120e11ba3965",
    "Accepted canonical implementation SHA")
  expressions <- as.list(parse(path, keep.source = FALSE))
  selected <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    is.symbol(x[[2L]]) && as.character(x[[2L]]) == "canonical_spec" &&
    is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), expressions)
  need(length(selected) == 1L, "Only explicit canonical function declaration")
  scope <- new.env(parent = globalenv())
  eval(selected[[1L]], scope)
  scope$canonical_spec
}

fit_methods <- function(y, labels, canonical_spec) {
  started <- proc.time()[[3L]]
  spec <- canonical_spec(labels)
  need(identical(rownames(spec$design), colnames(y)), "Fixed count/sample order")
  v <- limma::voom(y, spec$design, plot = FALSE)
  fitted <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, spec$design), spec$contrasts), robust = TRUE)
  voom <- rbindlist(lapply(colnames(contrasts), function(name) {
    effect <- fitted$coefficients[, name] * spec$orientation[name]
    se <- fitted$stdev.unscaled[, name] * sqrt(fitted$s2.post)
    crit <- qt(.975, fitted$df.total)
    data.table(gene_id = rownames(fitted$coefficients), contrast = name,
      method = "limma_voom_robust_eBayes", logFC = effect,
      t = fitted$t[, name] * spec$orientation[name], p_value = fitted$p.value[, name],
      q_value = p.adjust(fitted$p.value[, name], "BH"), se = se,
      ci_low = effect - crit * se, ci_high = effect + crit * se, df_total = fitted$df.total)
  }))
  voom_seconds <- proc.time()[[3L]] - started
  canonical_group <- factor(colnames(spec$design)[max.col(spec$design)], levels = colnames(spec$design))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(y$counts),
    colData = data.frame(group = canonical_group, row.names = colnames(y)), design = ~ 0 + group)
  conditions <- new.env(parent = emptyenv())
  conditions$messages <- conditions$warnings <- character()
  dds <- withCallingHandlers(DESeq2::DESeq(dds, test = "Wald", fitType = "parametric",
    betaPrior = FALSE, minReplicatesForReplace = Inf, quiet = TRUE),
    message = function(m) {
      conditions$messages <- c(conditions$messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }, warning = function(w) {
      conditions$warnings <- c(conditions$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  actual_fit <- attr(DESeq2::dispersionFunction(dds), "fitType")
  token <- "a local regression fit was automatically substituted"
  substitution <- any(grepl(token, c(conditions$messages, conditions$warnings), fixed = TRUE))
  need(length(actual_fit) == 1L && actual_fit %in% c("parametric", "local") &&
    identical(actual_fit == "local", substitution), "Original DESeq2 fit-type/message contract")
  coefficient_names <- DESeq2::resultsNames(dds)
  expected_names <- paste0("group", rownames(spec$contrasts))
  need(all(expected_names %in% coefficient_names), "Canonical DESeq2 coefficient names")
  deseq <- rbindlist(lapply(colnames(contrasts), function(name) {
    vector <- setNames(numeric(length(coefficient_names)), coefficient_names)
    vector[expected_names] <- spec$contrasts[, name]
    result <- as.data.table(as.data.frame(DESeq2::results(dds, contrast = vector,
      alpha = .05, independentFiltering = FALSE, cooksCutoff = FALSE)), keep.rownames = "gene_id")
    result[, `:=`(log2FoldChange = log2FoldChange * spec$orientation[name],
      stat = stat * spec$orientation[name])]
    result[, `:=`(contrast = name, method = "DESeq2_Wald_sensitivity",
      requested_fit_type = "parametric", actual_fit_type_after_DESeq = actual_fit,
      ci_low = log2FoldChange - qnorm(.975) * lfcSE,
      ci_high = log2FoldChange + qnorm(.975) * lfcSE, p_value = pvalue, q_value = padj)]
    result
  }))
  result <- rbindlist(list(
    voom[, .(collection = "genes", method, contrast, feature = gene_id, effect = logFC,
      statistic = abs(t), p_value, q_value)],
    deseq[, .(collection = "genes", method, contrast, feature = gene_id, effect = log2FoldChange,
      statistic = abs(stat), p_value, q_value)]))
  setorderv(result, key_columns)
  metadata <- as.data.frame(S4Vectors::mcols(dds))
  finite_summary <- function(value) {
    selected <- value[is.finite(value)]
    list(n = length(value), nonfinite = sum(!is.finite(value)),
      minimum = if (length(selected)) min(selected) else NA_real_,
      maximum = if (length(selected)) max(selected) else NA_real_)
  }
  diagnostic <- list(actual_fit_type = actual_fit, automatic_local_substitution = substitution,
    message_n = length(conditions$messages), warning_n = length(conditions$warnings),
    condition_text_sha256 = digest(c(conditions$messages, conditions$warnings), algo = "sha256"),
    betaConv_available = "betaConv" %in% names(metadata),
    betaConv_false_n = if ("betaConv" %in% names(metadata)) sum(!metadata$betaConv, na.rm = TRUE) else NA_integer_,
    betaConv_missing_n = if ("betaConv" %in% names(metadata)) sum(is.na(metadata$betaConv)) else NA_integer_,
    betaIter = if ("betaIter" %in% names(metadata)) finite_summary(metadata$betaIter) else NULL,
    dispIter = if ("dispIter" %in% names(metadata)) finite_summary(metadata$dispIter) else NULL,
    dispOutlier_n = if ("dispOutlier" %in% names(metadata)) sum(metadata$dispOutlier, na.rm = TRUE) else NA_integer_,
    dispersion = finite_summary(DESeq2::dispersions(dds)),
    size_factors = unname(DESeq2::sizeFactors(dds)), voom_seconds = voom_seconds,
    deseq_seconds = proc.time()[[3L]] - started - voom_seconds)
  list(result = result, details = list(voom = voom, deseq = deseq), diagnostic = diagnostic,
    canonical_specification = spec, elapsed_seconds = proc.time()[[3L]] - started)
}

validate_result <- function(value, state, id) {
  result <- value$result
  need(identical(value$id, id) && identical(value$labels, state$labels[id, ]) &&
    identical(value$source_binding, state$source_binding), "Allocation ID/labels/source")
  need(identical(names(result), c(key_columns, numeric_columns)) &&
    identical(as.data.frame(result[, ..key_columns]), state$expected_keys), "Four complete ordered families")
  for (field in numeric_columns) need(all(is.finite(result[[field]])), "Nonfinite result prevents acceptance; no label omitted")
  need(all(result$statistic >= 0 & result$p_value > 0 & result$p_value <= 1 &
    result$q_value > 0 & result$q_value <= 1), "Statistic and P/q domain")
  for (index in split(seq_len(nrow(result)), paste(result$method, result$contrast))) need(
    max(abs(p.adjust(result$p_value[index], "BH") - result$q_value[index])) <= 1e-12,
    "Complete separate native BH family")
  invisible(TRUE)
}

observed_bridge <- function(fit, model) {
  fields <- list(voom = c("logFC", "t", "p_value", "q_value", "se", "ci_low", "ci_high", "df_total"),
    deseq = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "ci_low", "ci_high", "p_value", "q_value"))
  rows <- list()
  for (method in names(fields)) {
    expected <- fread(file.path(model, if (method == "voom") "gene_voom.tsv" else "gene_deseq.tsv"))
    actual <- copy(fit$details[[method]])
    setorderv(expected, c("contrast", "gene_id")); setorderv(actual, c("contrast", "gene_id"))
    need(identical(expected$gene_id, actual$gene_id) && identical(expected$contrast, actual$contrast) &&
      identical(expected$method, actual$method) && identical(expected$symbol, actual$symbol),
      "Native full gene/contrast/method/annotation domain")
    if (method == "deseq") for (field in c("requested_fit_type", "actual_fit_type_after_DESeq")) need(
      identical(actual[[field]], expected[[field]]), "Observed original DESeq2 fit-type metadata")
    for (field in fields[[method]]) {
      left <- actual[[field]]; right <- expected[[field]]
      need(identical(is.na(left), is.na(right)), "Original exact missingness contract")
      index <- !is.na(left)
      delta <- abs(left[index] - right[index])
      bound <- 1e-10 + 1e-8 * pmax(abs(left[index]), abs(right[index]))
      need(all(is.finite(left[index]) & is.finite(right[index])) && all(delta <= bound),
        paste("Original complete-field observed numeric contract", method, field))
      rows[[length(rows) + 1L]] <- data.table(method = method, field = field,
        comparisons = length(left), missing_n = sum(!index), maximum_absolute_difference = max(delta),
        maximum_fraction_of_fixed_bound = max(delta / bound),
        maximum_log10_difference = if (field %in% c("p_value", "q_value", "pvalue", "padj"))
          max(abs(log10(left[index]) - log10(right[index]))) else NA_real_)
    }
  }
  rbindlist(rows)
}

save_allocation <- function(value, run) {
  path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", value$id))
  need(!file.exists(path) && !file.exists(paste0(path, ".partial")), "Existing checkpoint retained")
  saveRDS(value, paste0(path, ".partial"), version = 3)
  need(file.rename(paste0(path, ".partial"), path), "Atomic allocation checkpoint")
  path
}

selftest <- function(run) {
  canonical <- import_canonical(file.path(run, "canonical_source.R"))
  set.seed(20260817)
  means <- exp(runif(800L, log(20), log(5000)))
  counts <- matrix(rnbinom(7200L, mu = rep(means, 9L) * rep(seq(.7, 1.3, length.out = 9L), each = 800L),
    size = rep(1 / (.1 + 1 / means), 9L)), 800L, 9L)
  rownames(counts) <- paste0("synthetic_", 1:800)
  colnames(counts) <- paste0("GSM", 6235292:6235300)
  y <- edgeR::normLibSizes(edgeR::DGEList(counts), method = "TMM")
  observed_labels <- setNames(rep(groups, each = 3L), colnames(counts))
  original <- fit_methods(y, observed_labels, canonical)
  checks <- list()
  for (name in colnames(contrasts)) {
    pair <- groups[contrasts[, name] != 0]
    labels <- observed_labels
    labels[observed_labels == pair[1L]] <- pair[2L]
    labels[observed_labels == pair[2L]] <- pair[1L]
    swapped <- fit_methods(y, labels, canonical)
    left <- original$result[contrast == name]; right <- swapped$result[contrast == name]
    need(identical(left[, ..key_columns], right[, ..key_columns]) &&
      identical(left$statistic, right$statistic) && identical(left$p_value, right$p_value) &&
      identical(left$q_value, right$q_value) && identical(left$effect, -right$effect),
      "Synthetic canonical two-sided symmetry/scientific direction")
    checks[[name]] <- list(feature_n = nrow(left), exact_statistic_P_q = TRUE, opposite_effect = TRUE)
  }
  list(status = "PASS_SYNTHETIC_MAPPING_AND_BOTH_METHODS", synthetic_only = TRUE,
    tests = checks, packages = as.list(versions), maximum_workers = 4L)
}

prepare <- function(run, canonical_run, model) {
  need(!file.exists(file.path(run, "PREPARED.json")), "Existing preparation retained")
  need(fromJSON(file.path(run, "SYNTHETIC_TESTS.json"))$status == "PASS_SYNTHETIC_MAPPING_AND_BOTH_METHODS",
    "Synthetic tests before real fitting")
  canonical_run <- normalizePath(canonical_run, mustWork = TRUE)
  model <- normalizePath(model, mustWork = TRUE)
  files <- list(canonical_prepared = file.path(canonical_run, "PREPARED.json"),
    canonical_state = file.path(canonical_run, "derived_inputs.rds"),
    canonical_completion = file.path(canonical_run, "ENUMERATION_COMPLETE.json"),
    canonical_filter_proof = file.path(canonical_run, "FILTER_INVARIANCE.tsv"),
    canonical_input_manifest = file.path(canonical_run, "INPUT_MANIFEST.tsv"),
    native_completion = file.path(model, "COMPLETED"), native_manifest = file.path(model, "file_manifest.json"),
    native_voom = file.path(model, "gene_voom.tsv"), native_deseq = file.path(model, "gene_deseq.tsv"),
    original_model_source = file.path(run, "model_source.R"), original_config = file.path(run, "config.json"),
    canonical_source = file.path(run, "canonical_source.R"), worker = file.path(run, "worker.R"),
    definition = file.path(run, "definition.md"), synthetic_tests = file.path(run, "SYNTHETIC_TESTS.json"),
    summary_source = file.path(run, "summary.R"))
  identities <- lapply(files, pin)
  known <- c(canonical_prepared = "0631f763b0decbd3d5345770185ba6b45b09fb399a9cd45ba4840b21ba2391b4",
    canonical_state = "7a4405a6c046cdf0ca6ba9b2b2692e0c79d6ae92e6b90eb279a1aa46cd41662f",
    canonical_completion = "3ae21310e0f06d23aa8b1163c8acc423551a62385a17fcb17bbfd108a8ad1015",
    native_completion = "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e",
    native_manifest = "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7",
    native_voom = "b9dd2ca38772e7f232f73f89ffaaad49f94623b40cec2904768c9553b02c0c3c",
    native_deseq = "4048a82a9d8029ce2b3b87b498185557f9ecd16b306fbee6567f1c9fd5e69909",
    original_model_source = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2",
    original_config = "ee6e22592b06ecdb31707234e97dc55b7b66cbad8dd11c364187af0624ece980")
  for (role in names(known)) need(identities[[role]]$sha256 == known[[role]], paste("Accepted input", role))
  config <- fromJSON(files$original_config, simplifyVector = FALSE)
  native_manifest <- fromJSON(files$native_manifest, simplifyVector = FALSE)
  need(config$analysis$DESeq2$requested_fit_type_before_fit == "parametric" &&
    identical(unlist(config$analysis$DESeq2$allowed_actual_fit_types_after_DESeq), c("parametric", "local")) &&
    config$analysis$DESeq2$independent_filtering == FALSE && config$analysis$DESeq2$cooks_cutoff == FALSE &&
    config$analysis$independent_numeric_revalidation$absolute_tolerance == 1e-10 &&
    config$analysis$independent_numeric_revalidation$relative_tolerance == 1e-8,
    "Original method contract not changed")
  upstream <- fread(files$canonical_input_manifest)
  for (i in seq_len(nrow(upstream))) verify(as.list(upstream[i]))
  original <- readRDS(files$canonical_state)
  need(identical(dim(original$y$counts), c(14217L, 9L)) && identical(dim(original$labels), c(1680L, 9L)) &&
    original$spaces$observed == 1L && all(colSums(original$y$counts) == original$y$samples$lib.size) &&
    all(original$y$counts >= 0 & original$y$counts == round(original$y$counts)), "Pinned retained integer counts")
  proof <- fread(files$canonical_filter_proof)
  need(nrow(proof) == 1680L && all(proof$same_mask) && all(proof$retained_gene_n == 14217L), "Full upstream filter proof")
  need(isTRUE(fromJSON(files$canonical_completion)$scoped_numerical_acceptance), "Accepted upstream node")
  input_manifest <- rbindlist(lapply(names(identities), function(role) {
    identity <- identities[[role]]
    data.table(role = role, path = identity$path, bytes = identity$bytes, sha256 = identity$sha256)
  }))
  twrite(input_manifest, file.path(run, "INPUT_MANIFEST.tsv"))
  state <- original[c("y", "labels", "spaces")]
  state$source_binding <- pin(file.path(run, "INPUT_MANIFEST.tsv"))$sha256
  canonical <- import_canonical(files$canonical_source)
  fit <- fit_methods(state$y, state$labels[1L, ], canonical)
  # Reuse the upstream GTF-validated retained annotation, not either new bridge target.
  annotation_path <- upstream$path[basename(upstream$path) == "gene_edger.tsv"]
  need(length(annotation_path) == 1L, "Bound upstream GTF-validated annotation table")
  annotation <- unique(fread(annotation_path)[, .(gene_id, symbol)])
  need(nrow(annotation) == 14217L && !anyDuplicated(annotation$gene_id), "Unique retained annotation")
  for (method in names(fit$details)) {
    index <- match(fit$details[[method]]$gene_id, annotation$gene_id)
    need(!anyNA(index), "Complete retained annotation join")
    fit$details[[method]][, symbol := annotation$symbol[index]]
  }
  bridge <- observed_bridge(fit, model)
  need(nrow(fit$result) == 56868L, "New four-family full result domain")
  state$expected_keys <- as.data.frame(fit$result[, ..key_columns])
  dir.create(file.path(run, "inputs"), showWarnings = FALSE)
  dir.create(file.path(run, "allocations"), showWarnings = FALSE)
  value <- fit[c("result", "diagnostic", "canonical_specification", "elapsed_seconds")]
  value$id <- 1L; value$labels <- state$labels[1L, ]; value$source_binding <- state$source_binding
  validate_result(value, state, 1L)
  observed_path <- save_allocation(value, run)
  twrite(bridge, file.path(run, "OBSERVED_COMPATIBILITY.tsv"))
  twrite(fit$result, file.path(run, "OBSERVED.tsv"))
  for (method in names(fit$details)) twrite(fit$details[[method]], file.path(run, paste0("OBSERVED_DETAILS_", method, ".tsv")))
  saveRDS(state, file.path(run, "inputs", "state.rds"), version = 3)
  capture.output(sessionInfo(), file = file.path(run, "sessionInfo.txt"))
  for (identity in identities) verify(identity)
  for (i in seq_len(nrow(upstream))) verify(as.list(upstream[i]))
  jwrite(list(status = "PREPARED_COMPLETE_NATIVE_OBSERVED_COMPATIBLE", completed_utc = utc(),
    derived = pin(file.path(run, "inputs", "state.rds")), observed = pin(file.path(run, "OBSERVED.tsv")),
    observed_checkpoint = pin(observed_path), worker = identities$worker, definition = identities$definition,
    input_manifest = pin(file.path(run, "INPUT_MANIFEST.tsv")), packages = as.list(versions),
    all_labels = 1680L, companion_labels = c(20L, 20L), retained_genes = 14217L, result_rows = 56868L,
    observed_seconds = fit$elapsed_seconds, observed_diagnostics = fit$diagnostic,
    process_peak_memory_kb = peak_memory_kb(),
    full_numeric_native_bridge = TRUE, historical_runtime_authenticated = FALSE,
    authenticated_exchangeability = FALSE, candidate_gate_eligible = FALSE, primary_results_modified = FALSE),
    file.path(run, "PREPARED.json"))
  cat("PREPARED_FULL_NATIVE_BRIDGE; observed seconds", fit$elapsed_seconds, "\n")
}

load_state <- function(run) {
  prepared <- fromJSON(file.path(run, "PREPARED.json"), simplifyVector = FALSE)
  for (identity in list(prepared$derived, prepared$observed, prepared$observed_checkpoint,
    prepared$worker, prepared$definition, prepared$input_manifest)) verify(identity)
  inputs <- fread(prepared$input_manifest$path)
  for (i in seq_len(nrow(inputs))) verify(as.list(inputs[i]))
  list(state = readRDS(prepared$derived$path), prepared = prepared)
}

run_worker <- function(run, worker_id, worker_n) {
  need(worker_n == 4L && worker_id %in% 1:4, "Four single-thread workers")
  loaded <- load_state(run); state <- loaded$state
  canonical <- import_canonical(file.path(run, "canonical_source.R"))
  start_path <- file.path(run, paste0("WORKER_", worker_id, "_STARTED.json"))
  if (!file.exists(start_path)) jwrite(list(status = "RUNNING", pid = Sys.getpid(), time_utc = utc(),
    command = commandArgs(), worker = loaded$prepared$worker, derived = loaded$prepared$derived), start_path)
  for (id in state$spaces$global[(state$spaces$global - 1L) %% worker_n == worker_id - 1L]) {
    path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", id))
    if (file.exists(path)) { validate_result(readRDS(path), state, id); next }
    fit <- fit_methods(state$y, state$labels[id, ], canonical)
    value <- fit[c("result", "diagnostic", "canonical_specification", "elapsed_seconds")]
    value$id <- id; value$labels <- state$labels[id, ]; value$source_binding <- state$source_binding
    # Save complete labelled output even when validation would stop; never silently omit a failed label.
    save_allocation(value, run)
    validate_result(value, state, id)
    if (id %% 40L < 4L) cat("completed", id, "of1680; seconds", value$elapsed_seconds,
      "; betaConv false", value$diagnostic$betaConv_false_n, "\n")
  }
  jwrite(list(status = "WORKER_COMPLETE", worker_id = worker_id, time_utc = utc(),
    worker = loaded$prepared$worker, derived = loaded$prepared$derived),
    file.path(run, paste0("WORKER_", worker_id, "_COMPLETE.json")))
  cat("WORKER_COMPLETE", worker_id, "\n")
}

args <- commandArgs(trailingOnly = TRUE)
need(length(args) >= 2L && args[1L] %in% c("selftest", "prepare", "run"), "Mode and explicit run directory required")
run <- normalizePath(args[2L], mustWork = TRUE)
for (package in names(versions)) need(as.character(packageVersion(package)) == versions[package], "Native package version")
set.seed(20260817)
data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
if (args[1L] == "selftest") {
  jwrite(selftest(run), file.path(run, "SYNTHETIC_TESTS.json")); cat("PASS_SYNTHETIC_ONLY\n")
} else if (args[1L] == "prepare") {
  need(length(args) == 4L, "Explicit upstream and native model directories")
  prepare(run, args[3L], args[4L])
} else {
  need(length(args) == 4L, "Explicit worker ID and count")
  run_worker(run, as.integer(args[3L]), as.integer(args[4L]))
}
