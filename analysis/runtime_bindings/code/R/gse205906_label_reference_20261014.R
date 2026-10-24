#!/usr/bin/env Rscript
# Explicit new calibration node. Never execute either historical model driver.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
pin <- function(path) list(path = normalizePath(path, mustWork = TRUE),
  bytes = unname(file.info(path)$size),
  sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
jwrite <- function(value, path) {
  need(!file.exists(path), paste("Existing output", basename(path)))
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, digits = 17,
    na = "null", null = "null")
}
twrite <- function(value, path) {
  need(!file.exists(path), paste("Existing output", basename(path)))
  data.table::fwrite(value, path, sep = "\t", quote = TRUE, na = "NA")
}
read_table <- function(path) data.table::fread(path, na.strings = "NA", check.names = FALSE)
close_num <- function(a, b, tolerance = 1e-8) length(a) == length(b) &&
  all(is.finite(a) & is.finite(b) & abs(a - b) <= tolerance * (1 + abs(b)))
groups <- c("sham", "SNI", "PLM")
contrast_matrix <- cbind(local_disease = c(-1, 1, 0), intervention = c(0, -1, 1))
rownames(contrast_matrix) <- groups
make_design <- function(group) {
  design <- model.matrix(~ 0 + factor(group, levels = groups))
  colnames(design) <- groups
  design
}
enumerate_labels <- function() {
  rows <- list()
  for (sham in combn(9L, 3L, simplify = FALSE)) {
    for (sni in combn(setdiff(1:9, sham), 3L, simplify = FALSE)) {
      group <- rep("PLM", 9L)
      group[sham] <- "sham"; group[sni] <- "SNI"
      rows[[length(rows) + 1L]] <- group
    }
  }
  do.call(rbind, rows)
}
check_labels <- function(labels) {
  observed <- rep(groups, each = 3L)
  identity <- which(apply(labels, 1, function(g) identical(unname(g), observed)))
  local <- which(apply(labels[, 7:9, drop = FALSE], 1, function(g) all(g == "PLM")))
  intervention <- which(apply(labels[, 1:3, drop = FALSE], 1, function(g) all(g == "sham")))
  need(nrow(labels) == 1680L && ncol(labels) == 9L && length(identity) == 1L &&
    !anyDuplicated(apply(labels, 1, paste, collapse = ":")) &&
    all(apply(labels, 1, function(g) all(table(factor(g, levels = groups)) == 3L))) &&
    length(local) == 20L && length(intervention) == 20L &&
    length(intersect(local, intervention)) == 1L, "Allocation space/observed identity")
  list(observed = identity, global = 1:1680, local_disease = local, intervention = intervention)
}
extract_helpers <- function(path) {
  need(pin(path)$sha256 == "68c01474c53840441977b9bba571f816ec7d2083677fbc78b112af46f3b64af1", "Helper source identity")
  allowed <- c("read_gtf_symbols", "plain_frame", "build_networks", "run_ulm", "analytic_ulm")
  expressions <- parse(path)
  scope <- new.env(parent = globalenv())
  for (name in allowed) {
    selected <- Filter(function(z) is.call(z) && identical(z[[1]], as.name("<-")) &&
      identical(z[[2]], as.name(name)) && is.call(z[[3]]) &&
      identical(z[[3]][[1]], as.name("function")), as.list(expressions))
    need(length(selected) == 1L, paste("Explicit helper function", name))
    eval(selected[[1]], envir = scope)
  }
  scope
}
self_test <- function() {
  labels <- enumerate_labels(); spaces <- check_labels(labels)
  toy <- matrix(round(exp((seq_len(4500) %% 53) / 15) * 12), 500, 9)
  y <- edgeR::DGEList(toy)
  keep <- edgeR::filterByExpr(y, make_design(labels[1, ]))
  for (i in spaces$global) {
    design <- make_design(labels[i, ])
    need(close_num(rowSums(qr.Q(qr(design))^2), rep(1/3, 9), 1e-12), "Synthetic leverage")
    need(identical(edgeR::filterByExpr(y, design), keep), "Synthetic filter invariance")
  }
  need(sum(c(2, 3, 4, 1) >= 2) / 4 == (1 + sum(c(3, 4, 1) >= 2)) / 4,
    "Exhaustive/plus-one observed counted once")
  need(sum(abs(c(-4, -2, 0, 2, 4)) >= 4) / 5 == 2/5, "Two-sided ties")
  list(status = "PASS_SYNTHETIC_ONLY", allocations = nrow(labels), pairwise_n = 20L,
    checks = c("exact_space_identity", "balanced_design_leverage", "filter_invariance_all_labels",
      "plus_one_without_double_observed", "two_sided_ties"))
}
verify_inputs <- function(table) {
  for (i in seq_len(nrow(table))) {
    actual <- pin(table$path[i])
    need(actual$bytes == table$bytes[i] && actual$sha256 == table$sha256[i],
      paste("Input identity", table$role[i]))
  }
}
load_fixed_inputs <- function(run) {
  formal <- "/root/autodl-tmp/GSE205906_full_read_counts_fasp_v18e_20261005/formal"
  model <- "/root/autodl-tmp/GSE205906_Ensembl101_checked_model_v5_20261006"
  public <- "/root/autodl-tmp/route_b2_public_analysis/public-input-stage-v1"
  recovery <- "/root/autodl-tmp/route_b2_gse205906_program_baseline_20261014r01"
  paths <- list(matrix = file.path(formal, "output/GSE205906_fixed9_gene_fragments_s0.tsv"),
    overlay = file.path(formal, "derived/NINE_LIBRARY_ANALYSIS_PREPARATION_OVERLAY_154.tsv"),
    GTF = "/root/autodl-tmp/GSE205906_full_read_reconstruction_156_20261003/references/derived/Rattus_norvegicus.Rnor_6.0.101.gtf",
    resources = file.path(public, "data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds"),
    recovery_scores = file.path(recovery, "restored_scores.rds"),
    helper = file.path(run, "baseline_helpers_source.R"),
    worker = file.path(run, "worker.R"), definition = file.path(run, "definition.md"))
  manifest_path <- file.path(model, "file_manifest.json")
  complete_path <- file.path(model, "COMPLETED")
  need(pin(manifest_path)$sha256 == "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7" &&
    pin(complete_path)$sha256 == "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e", "Accepted native node identity")
  original <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)$provenance$eleven_input_pins
  complete <- jsonlite::fromJSON(complete_path, simplifyVector = FALSE)
  expected <- c(matrix = original$matrix$sha256, overlay = original$overlay$sha256,
    GTF = original$GTF$sha256, resources = original$frozen_resource_RDS$sha256,
    recovery_scores = "cd7cc354abc7d64227b02ba45951ddf195c3234d8cb436c6a0b3f164a8c1532b")
  for (role in names(expected)) need(pin(paths[[role]])$sha256 == expected[[role]], paste("Frozen input", role))
  for (name in c("gene_edger.tsv", "activity_results.tsv", "unsigned_camera.tsv", "mapping_coverage.tsv", "design_out.tsv", "contrast_out.tsv")) {
    paths[[name]] <- file.path(model, name)
    need(pin(paths[[name]])$sha256 == complete$files[[name]]$sha256, paste("Native output", name))
  }
  paths$native_manifest <- manifest_path; paths$native_complete <- complete_path
  table <- data.table::rbindlist(lapply(names(paths), function(role) {
    value <- pin(paths[[role]])
    data.table::data.table(role = role, path = value$path, bytes = value$bytes, sha256 = value$sha256)
  }))
  list(paths = paths, manifest = table)
}
map_unsigned <- function(resource, symbols) {
  orth <- data.table::as.data.table(resource$ortholog_one_to_one)
  map <- unique(orth[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
    .(human_symbol = as.character(human_symbol), symbol = as.character(other_symbol))])
  need(!anyDuplicated(map$human_symbol) && !anyDuplicated(map$symbol), "Strict unsigned orthology")
  sources <- list(Reactome = resource$reactome, SynGO_sensitivity = resource$syngo_sensitivity)
  lapply(sources, function(source) {
    source <- data.table::as.data.table(source)
    mapped <- unique(merge(unique(source[, .(program_id = as.character(term_id),
      human_symbol = as.character(gene_symbol))]), map, by = "human_symbol", allow.cartesian = TRUE))
    members <- lapply(split(mapped$symbol, mapped$program_id), function(z) unique(intersect(z, symbols)))
    members <- members[lengths(members) >= 10L & lengths(members) <= 500L]
    lapply(members, function(z) which(symbols %in% z))
  })
}
fit_allocation <- function(state, labels, allocation_id) {
  started <- proc.time()[[3]]
  design <- make_design(labels)
  rownames(design) <- colnames(state$y)
  y <- edgeR::estimateDisp(state$y, design, robust = TRUE)
  qlf <- edgeR::glmQLFit(y, design, robust = TRUE)
  rows <- list()
  for (name in colnames(contrast_matrix)) {
    result <- edgeR::glmQLFTest(qlf, contrast = contrast_matrix[, name])$table
    rows[[length(rows) + 1L]] <- data.table::data.table(collection = "genes", method = "edgeR_robust_QL",
      contrast = name, feature = rownames(result), effect = result$logFC,
      statistic = result$F, p_value = result$PValue)
  }
  for (collection in names(state$scores)) {
    score <- state$scores[[collection]]
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(score, design), contrast_matrix), robust = TRUE)
    for (name in colnames(contrast_matrix)) {
      rows[[length(rows) + 1L]] <- data.table::data.table(collection = collection,
        method = "ULM_robust_limma", contrast = name, feature = rownames(score),
        effect = fit$coefficients[, name], statistic = abs(fit$t[, name]), p_value = fit$p.value[, name])
    }
  }
  for (collection in names(state$indices)) {
    for (name in colnames(contrast_matrix)) {
      result <- limma::camera(state$expression, state$indices[[collection]], design,
        contrast = contrast_matrix[, name], inter.gene.cor = .01)
      rows[[length(rows) + 1L]] <- data.table::data.table(collection = collection, method = "CAMERA_cor_0.01",
        contrast = name, feature = rownames(result), effect = ifelse(result$Direction == "Up", 1, -1),
        statistic = -log(result$PValue), p_value = result$PValue)
    }
  }
  result <- data.table::rbindlist(rows)
  result[, q_value := p.adjust(p_value, "BH"), by = .(collection, method, contrast)]
  need(!anyNA(result) && all(is.finite(result$effect)) && all(is.finite(result$statistic)) &&
    all(is.finite(result$p_value) & result$p_value > 0 & result$p_value <= 1) &&
    !anyDuplicated(result[, .(collection, method, contrast, feature)]), "Finite complete allocation inference")
  data.table::setorderv(result, c("collection", "method", "contrast", "feature"))
  list(id = allocation_id, labels = labels, result = result,
    source_binding = state$source_binding,
    elapsed_seconds = proc.time()[[3]] - started)
}
validate_allocation <- function(value, state, allocation_id) {
  need(is.list(value) && identical(value$id, allocation_id) &&
    identical(value$labels, state$labels[allocation_id, ]) &&
    identical(value$source_binding, state$source_binding), "Allocation run/source binding")
  result <- value$result
  keys <- c("collection", "method", "contrast", "feature")
  fields <- c(keys, "effect", "statistic", "p_value", "q_value")
  need(data.table::is.data.table(result) && identical(names(result), fields) &&
    identical(as.data.frame(result[, ..keys]), state$expected_keys) && !anyNA(result),
    "Allocation schema and complete ordered families")
  for (field in c("effect", "statistic", "p_value", "q_value"))
    need(all(is.finite(result[[field]])), "Allocation finite fields")
  need(all(result$statistic >= 0 & result$p_value > 0 & result$p_value <= 1 &
    result$q_value > 0 & result$q_value <= 1), "Allocation valid statistics/P/q")
  recalculated <- result[, .(q = p.adjust(p_value, "BH")), by = .(collection, method, contrast)]
  need(close_num(recalculated$q, result$q_value, 1e-12), "Allocation native BH family")
  invisible(TRUE)
}
compare_observed <- function(result, paths) {
  gene <- read_table(paths$gene_edger.tsv)
  activity <- read_table(paths$activity_results.tsv)
  camera <- read_table(paths$unsigned_camera.tsv)
  expected <- data.table::rbindlist(list(
    gene[, .(collection = "genes", method = "edgeR_robust_QL", contrast, feature = gene_id,
      effect = logFC, statistic = F, p_value = PValue, q_value)],
    activity[, .(collection, method = "ULM_robust_limma", contrast, feature = program_id,
      effect, statistic = abs(t), p_value, q_value = q_within_collection_contrast)],
    camera[, .(collection, method = "CAMERA_cor_0.01", contrast, feature = program_id,
      effect = ifelse(Direction == "Up", 1, -1), statistic = -log(p_value), p_value,
      q_value = q_within_collection_contrast)]))
  data.table::setorderv(expected, c("collection", "method", "contrast", "feature"))
  keys <- c("collection", "method", "contrast", "feature")
  need(identical(as.data.frame(result[, ..keys]), as.data.frame(expected[, ..keys])), "Observed complete native domains")
  for (field in c("effect", "statistic")) need(close_num(result[[field]], expected[[field]]), paste("Observed native", field))
  for (field in c("p_value", "q_value")) need(all(is.finite(log10(result[[field]]))) &&
    all(abs(log10(result[[field]]) - log10(expected[[field]])) <= 1e-6), paste("Observed native absolute log10", field))
  data.table::data.table(field = c("effect", "statistic", "p_value", "q_value"),
    maximum_absolute_difference = vapply(c("effect", "statistic", "p_value", "q_value"),
      function(field) max(abs(result[[field]] - expected[[field]])), numeric(1)))
}
prepare <- function(run) {
  need(!file.exists(file.path(run, "PREPARED.json")), "Prepared node already exists")
  input <- load_fixed_inputs(run); paths <- input$paths
  helper <- extract_helpers(paths$helper)
  meta <- read_table(paths$design_out.tsv); overlay <- read_table(paths$overlay)
  need(identical(meta$sample_id, paste0("GSM", 6235292:6235300)) &&
    identical(meta$group, rep(groups, each = 3L)) && identical(overlay$GEO_sample, meta$sample_id) &&
    identical(overlay$group, rep(c("JA", "SNI", "PLM"), each = 3L)), "Nine-library crosswalk")
  native_contrasts <- read_table(paths$contrast_out.tsv)
  need(identical(native_contrasts$coefficient, groups), "Contrast group order")
  for (name in colnames(contrast_matrix)) need(identical(as.numeric(native_contrasts[[name]]), unname(contrast_matrix[, name])), "Frozen contrast")
  table <- read_table(paths$matrix)
  need(names(table)[1] == "gene_id" && identical(names(table)[-1], overlay$label) &&
    nrow(table) == 32883L && !anyNA(table$gene_id) && !anyDuplicated(table$gene_id) &&
    identical(table$gene_id, sort(table$gene_id, method = "radix")), "Full count matrix schema/order")
  count <- as.matrix(table[, -1, with = FALSE]); rownames(count) <- table$gene_id
  colnames(count) <- meta$sample_id
  need(all(is.finite(count) & count >= 0 & count == round(count)), "Nonnegative integer counts")
  labels <- enumerate_labels(); colnames(labels) <- meta$sample_id
  spaces <- check_labels(labels)
  y0 <- edgeR::DGEList(count)
  cpm_keep <- rowSums(edgeR::cpm(y0, log = FALSE) >= 1) >= 3L
  keep <- edgeR::filterByExpr(y0, make_design(labels[spaces$observed, ])) & cpm_keep
  proof <- lapply(spaces$global, function(i) {
    design <- make_design(labels[i, ])
    candidate <- edgeR::filterByExpr(y0, design) & cpm_keep
    need(identical(candidate, keep), paste("Real filter mask invariance", i))
    need(close_num(rowSums(qr.Q(qr(design))^2), rep(1/3, 9), 1e-12), "Real leverage")
    data.table::data.table(allocation_id = i, retained_gene_n = sum(candidate), same_mask = TRUE)
  })
  need(sum(keep) == 14217L, "Native retained gene count")
  annotation <- helper$read_gtf_symbols(paths$GTF, rownames(count))
  native <- read_table(paths$gene_edger.tsv)
  for (name in colnames(contrast_matrix)) need(identical(native[contrast == name]$gene_id, rownames(count)[keep]) &&
    identical(native[contrast == name]$symbol, annotation$symbols[keep]), "Native gene/annotation alignment")
  y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  expression <- log2(sweep(y$counts + .5, 2, y$samples$lib.size * y$samples$norm.factors + 1, "/") * 1e6)
  need(close_num(as.vector(expression), as.vector(limma::voom(y, make_design(meta$group), plot = FALSE)$E), 1e-12), "Exact voom E expression")
  symbols <- annotation$symbols[keep]
  collapsed <- rowsum(expression, symbols, reorder = FALSE)
  collapsed <- collapsed / as.numeric(table(factor(symbols, levels = rownames(collapsed))))
  resource <- readRDS(paths$resources)
  networks <- helper$build_networks(resource, rownames(collapsed))
  coverage_native <- read_table(paths$mapping_coverage.tsv)
  coverage_new <- networks$coverage
  coverage_keys <- paste(coverage_new$collection, coverage_new$source)
  matched <- match(coverage_keys, paste(coverage_native$collection, coverage_native$source))
  need(!anyNA(matched), "Coverage domain match")
  for (field in c("source_feature_n", "tested_feature_n", "tested")) need(
    all(coverage_new[[field]] == coverage_native[[field]][matched]), paste("Native eligibility", field))
  scores <- lapply(networks$networks, function(network) helper$run_ulm(collapsed, network))
  recovered <- readRDS(paths$recovery_scores)
  all_scores <- do.call(rbind, scores)
  rownames(all_scores) <- unlist(lapply(names(scores), function(family) paste(family, rownames(scores[[family]]), sep = "::")))
  need(nrow(collapsed) == 14188L && nrow(all_scores) == 95L &&
    identical(dimnames(all_scores), dimnames(recovered$scores)) &&
    close_num(as.vector(all_scores), as.vector(recovered$scores), 1e-12), "Independent accepted score recovery bridge")
  state <- list(y = y, expression = expression, scores = scores,
    indices = map_unsigned(resource, symbols), labels = labels, spaces = spaces,
    source_binding = digest::digest(input$manifest, algo = "sha256"))
  observed <- fit_allocation(state, labels[spaces$observed, ], spaces$observed)
  comparison <- compare_observed(observed$result, paths)
  keys <- c("collection", "method", "contrast", "feature")
  state$expected_keys <- as.data.frame(observed$result[, ..keys])
  validate_allocation(observed, state, spaces$observed)
  verify_inputs(input$manifest)
  dir.create(file.path(run, "allocations"), showWarnings = FALSE)
  twrite(input$manifest, file.path(run, "INPUT_MANIFEST.tsv"))
  twrite(data.table::rbindlist(proof), file.path(run, "FILTER_INVARIANCE.tsv"))
  twrite(cbind(data.table::data.table(allocation_id = 1:1680), data.table::as.data.table(labels)), file.path(run, "ALLOCATIONS.tsv"))
  twrite(comparison, file.path(run, "OBSERVED_COMPATIBILITY.tsv"))
  twrite(observed$result, file.path(run, "OBSERVED.tsv"))
  saveRDS(state, file.path(run, "derived_inputs.rds"), version = 3)
  saveRDS(observed, file.path(run, "allocations", sprintf("allocation_%04d.rds", spaces$observed)), version = 3)
  capture.output(sessionInfo(), file = file.path(run, "sessionInfo.txt"))
  jwrite(list(status = "PREPARED_OBSERVED_NATIVE_COMPATIBLE", completed_utc = utc(),
    input_manifest = pin(file.path(run, "INPUT_MANIFEST.tsv")), derived = pin(file.path(run, "derived_inputs.rds")),
    observed = pin(file.path(run, "OBSERVED.tsv")), worker = pin(paths$worker), definition = pin(paths$definition),
    observed_checkpoint = pin(file.path(run, "allocations", sprintf("allocation_%04d.rds", spaces$observed))),
    library_n = 9L, raw_gene_n = 32883L, retained_gene_n = sum(keep), symbol_n = nrow(collapsed),
    signed_program_n = nrow(all_scores), all_labels_filter_verified = 1680L,
    observed_id = spaces$observed, observed_seconds = observed$elapsed_seconds,
    R = R.version.string, package_versions = setNames(lapply(c("edgeR", "limma", "decoupleR", "data.table", "jsonlite", "digest"),
      function(package) as.character(utils::packageVersion(package))), c("edgeR", "limma", "decoupleR", "data.table", "jsonlite", "digest")),
    authenticated_animal_batch_exchangeability = FALSE, primary_results_modified = FALSE,
    other_controls_complete = FALSE), file.path(run, "PREPARED.json"))
  cat("PREPARED; observed seconds:", observed$elapsed_seconds, "\n")
}
load_prepared <- function(run) {
  prepared <- jsonlite::fromJSON(file.path(run, "PREPARED.json"), simplifyVector = FALSE)
  for (item in list(prepared$derived, prepared$observed, prepared$worker, prepared$definition,
    prepared$input_manifest, prepared$observed_checkpoint))
    need(pin(item$path)$sha256 == item$sha256, "Frozen prepared stage identity")
  need(prepared$status == "PREPARED_OBSERVED_NATIVE_COMPATIBLE", "Prepared stage status")
  verify_inputs(read_table(file.path(run, "INPUT_MANIFEST.tsv")))
  readRDS(prepared$derived$path)
}
run_worker <- function(run, worker_id, worker_n) {
  need(worker_n >= 1 && worker_n <= 4 && worker_id >= 1 && worker_id <= worker_n, "Four-worker bound")
  state <- load_prepared(run)
  for (i in state$spaces$global[(state$spaces$global - 1L) %% worker_n == worker_id - 1L]) {
    path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", i))
    if (file.exists(path)) {
      result <- readRDS(path)
      validate_allocation(result, state, i)
      next
    }
    value <- fit_allocation(state, state$labels[i, ], i)
    validate_allocation(value, state, i)
    saveRDS(value, paste0(path, ".partial"), version = 3)
    need(file.rename(paste0(path, ".partial"), path), "Atomic allocation checkpoint")
    if (i %% 20L < worker_n) cat("completed", i, "seconds", value$elapsed_seconds, "\n")
  }
  jwrite(list(status = "WORKER_COMPLETE", worker_id = worker_id, worker_n = worker_n,
    completed_utc = utc(), worker = pin(file.path(run, "worker.R"))),
    file.path(run, paste0("WORKER_", worker_id, "_COMPLETE.json")))
}
args <- commandArgs(trailingOnly = TRUE)
need(length(args) >= 2L && args[1] %in% c("selftest", "prepare", "run"), "Mode and explicit run directory required")
mode <- args[1]; run <- normalizePath(args[2], mustWork = TRUE)
set.seed(20260817)
need(as.character(utils::packageVersion("edgeR")) == "4.10.1" &&
  as.character(utils::packageVersion("limma")) == "3.68.4", "Native inference versions")
if (mode == "selftest") {
  jwrite(self_test(), file.path(run, "SYNTHETIC_TESTS.json")); cat("PASS_SYNTHETIC_ONLY\n")
} else if (mode == "prepare") {
  need(jsonlite::fromJSON(file.path(run, "SYNTHETIC_TESTS.json"))$status == "PASS_SYNTHETIC_ONLY", "Synthetic tests required")
  prepare(run)
} else {
  need(length(args) == 4L, "Run requires worker ID/count")
  run_worker(run, as.integer(args[3]), as.integer(args[4]))
}
