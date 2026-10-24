#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
pin <- function(path) list(path = normalizePath(path, mustWork = TRUE),
  bytes = unname(file.info(path)$size),
  sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
verify <- function(value) {
  actual <- pin(value$path)
  need(actual$bytes == value$bytes && actual$sha256 == value$sha256, "Frozen file identity")
  invisible(actual)
}
jwrite <- function(value, path) {
  need(!file.exists(path), paste("Existing output", basename(path)))
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null")
}
twrite <- function(value, path) {
  need(!file.exists(path), paste("Existing output", basename(path)))
  data.table::fwrite(value, path, sep = "\t", quote = TRUE, na = "NA")
}
groups <- c("sham", "SNI", "PLM")
contrasts <- cbind(local_disease = c(-1, 1, 0), intervention = c(0, -1, 1))
rownames(contrasts) <- groups
key_columns <- c("collection", "method", "contrast", "feature")
numeric_columns <- c("effect", "statistic", "p_value", "q_value")
canonical_spec <- function(labels) {
  original <- model.matrix(~ 0 + factor(labels, levels = groups))
  colnames(original) <- groups
  rownames(original) <- names(labels)
  permutation <- order(vapply(seq_len(3L), function(j) which(original[, j] == 1)[1L], integer(1)))
  design <- original[, permutation, drop = FALSE]
  colnames(design) <- paste0("cell_", 1:3)
  vectors <- contrasts[permutation, , drop = FALSE]
  rownames(vectors) <- colnames(design)
  orientation <- vapply(seq_len(ncol(vectors)), function(j) sign(tail(vectors[vectors[, j] != 0, j], 1L)), numeric(1))
  names(orientation) <- colnames(vectors)
  vectors <- sweep(vectors, 2, orientation, "*")
  need(qr(design)$rank == 3L && all(rowSums(design) == 1) && all(colSums(design) == 3), "Canonical full rank/cells")
  for (name in colnames(contrasts)) need(all(vectors[, name] * orientation[name] == contrasts[permutation, name]),
    "Canonical scientific contrast mapping")
  list(design = design, contrasts = vectors, orientation = orientation, permutation = permutation)
}
all_labels <- function() {
  rows <- list()
  for (sham in combn(9L, 3L, simplify = FALSE)) {
    for (sni in combn(setdiff(1:9, sham), 3L, simplify = FALSE)) {
      labels <- rep("PLM", 9L)
      labels[sham] <- "sham"
      labels[sni] <- "SNI"
      rows[[length(rows) + 1L]] <- labels
    }
  }
  labels <- do.call(rbind, rows)
  colnames(labels) <- paste0("GSM", 6235292:6235300)
  labels
}
check_mapping <- function(labels) {
  signatures <- apply(labels, 1L, paste, collapse = ":")
  need(identical(dim(labels), c(1680L, 9L)) && !anyDuplicated(signatures), "Complete labelled space")
  for (i in 1:1680) {
    specification <- canonical_spec(labels[i, ])
    for (name in colnames(contrasts)) {
      pair <- groups[contrasts[, name] != 0]
      swapped <- labels[i, ]
      swapped[labels[i, ] == pair[1L]] <- pair[2L]
      swapped[labels[i, ] == pair[2L]] <- pair[1L]
      partner <- match(paste(swapped, collapse = ":"), signatures)
      need(!is.na(partner), "Swap closure")
      other <- canonical_spec(labels[partner, ])
      need(identical(specification$design, other$design) &&
        identical(specification$contrasts[, name], other$contrasts[, name]) &&
        specification$orientation[name] == -other$orientation[name], "Canonical swap identity/direction")
    }
  }
  need(sum(rowSums(labels[, 7:9, drop = FALSE] == "PLM") == 3L) == 20L &&
    sum(rowSums(labels[, 1:3, drop = FALSE] == "sham") == 3L) == 20L, "Companion spaces")
  invisible(TRUE)
}
fit_genes <- function(y, labels) {
  specification <- canonical_spec(labels)
  need(identical(rownames(specification$design), colnames(y)), "Fixed sample order")
  estimated <- edgeR::estimateDisp(y, specification$design, robust = TRUE)
  fitted <- edgeR::glmQLFit(estimated, specification$design, robust = TRUE)
  need(!is.null(fitted$average.ql.dispersion) && is.null(fitted$df.residual.zeros), "Original modern QL branch")
  rows <- lapply(colnames(contrasts), function(name) {
    result <- edgeR::glmQLFTest(fitted, contrast = specification$contrasts[, name])$table
    data.table::data.table(collection = "genes", method = "edgeR_robust_QL", contrast = name,
      feature = rownames(result), effect = result$logFC * specification$orientation[name],
      statistic = result$F, p_value = result$PValue, q_value = p.adjust(result$PValue, "BH"))
  })
  result <- data.table::rbindlist(rows)
  data.table::setorderv(result, key_columns)
  list(result = result, specification = specification)
}
validate_result <- function(result, state) {
  need(data.table::is.data.table(result) && identical(names(result), c(key_columns, numeric_columns)) &&
    identical(as.data.frame(result[, ..key_columns]), state$expected_keys) && !anyNA(result), "Full ordered result domain")
  for (field in numeric_columns) need(all(is.finite(result[[field]])), "Finite full inference")
  need(all((result$statistic >= 0 | (result$collection == "genes" & result$statistic >= -1e-10)) &
    result$p_value > 0 & result$p_value <= 1 & result$q_value > 0 & result$q_value <= 1), "Original statistic/P/q bounds")
  families <- split(seq_len(nrow(result)), paste(result$collection, result$method, result$contrast, sep = "::"))
  for (index in families) need(all(abs(p.adjust(result$p_value[index], "BH") - result$q_value[index]) <= 1e-12),
    "Complete separate native BH families")
  invisible(TRUE)
}
read_original <- function(state, manifest, id) {
  path <- file.path(state$original_run, "allocations", sprintf("allocation_%04d.rds", id))
  identity <- pin(path)
  need(identity$sha256 == manifest$sha256[id] && identity$bytes == manifest$bytes[id], "Original allocation SHA/bytes")
  value <- readRDS(path)
  need(identical(value$id, id) && identical(value$labels, state$labels[id, ]) &&
    identical(value$source_binding, state$original_source_binding), "Original ID/labels/source")
  validate_result(value$result, state)
  list(value = value, identity = identity)
}
assemble <- function(state, manifest, id) {
  started <- proc.time()[[3L]]
  original <- read_original(state, manifest, id)
  fresh <- fit_genes(state$y, state$labels[id, ])
  old_non_gene <- original$value$result[collection != "genes"]
  result <- data.table::rbindlist(list(old_non_gene, fresh$result))
  data.table::setorderv(result, key_columns)
  need(identical(result[collection != "genes"], old_non_gene), "Exact same-label non-gene field preservation")
  validate_result(result, state)
  old_gene <- original$value$result[collection == "genes"]
  new_gene <- result[collection == "genes"]
  impact <- data.table::rbindlist(lapply(colnames(contrasts), function(name) {
    old <- old_gene[contrast == name]
    new <- new_gene[contrast == name]
    need(identical(old$feature, new$feature), "Repair gene alignment")
    data.table::data.table(allocation_id = id, contrast = name, genes = nrow(new),
      maximum_absolute_effect_difference = max(abs(new$effect - old$effect)),
      maximum_scaled_effect_difference = max(abs(new$effect - old$effect) / (1 + pmax(abs(new$effect), abs(old$effect)))),
      maximum_absolute_F_difference = max(abs(new$statistic - old$statistic)),
      maximum_scaled_F_difference = max(abs(new$statistic - old$statistic) / (1 + pmax(abs(new$statistic), abs(old$statistic)))),
      maximum_absolute_log10_P_difference = max(abs(log10(new$p_value) - log10(old$p_value))),
      maximum_absolute_log10_q_difference = max(abs(log10(new$q_value) - log10(old$q_value))),
      old_BH_rejections = sum(old$q_value < .05), new_BH_rejections = sum(new$q_value < .05),
      BH_status_changed_n = sum((old$q_value < .05) != (new$q_value < .05)),
      minimum_new_F = min(new$statistic), negative_new_F_n = sum(new$statistic < 0))
  }))
  list(id = id, labels = state$labels[id, ], result = result,
    source_binding = state$source_binding, elapsed_seconds = proc.time()[[3L]] - started,
    original_allocation = original$identity, canonical_specification = fresh$specification,
    repair_impact = impact, new_independent_gene_refit = TRUE, same_label_non_gene_exact = TRUE)
}
observed_bridge <- function(value, original) {
  rows <- lapply(numeric_columns, function(field) {
    actual <- value$result[[field]]
    expected <- original$result[[field]]
    if (field %in% c("p_value", "q_value")) {
      difference <- max(abs(log10(actual) - log10(expected)))
      need(difference <= 1e-6, "Fixed observed absolute log10(P/q) compatibility")
    } else {
      difference <- max(abs(actual - expected) / (1 + abs(expected)))
      need(difference <= 1e-8, "Fixed observed effect/statistic compatibility")
    }
    data.table::data.table(field = field, maximum_absolute_difference = max(abs(actual - expected)),
      maximum_fixed_bridge_metric = difference)
  })
  data.table::rbindlist(rows)
}
save_allocation <- function(value, run) {
  path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", value$id))
  need(!file.exists(path) && !file.exists(paste0(path, ".partial")), "Existing allocation retained")
  saveRDS(value, paste0(path, ".partial"), version = 3)
  need(file.rename(paste0(path, ".partial"), path), "Atomic allocation checkpoint")
  path
}
selftest <- function() {
  labels <- all_labels()
  check_mapping(labels)
  counts <- matrix(round(exp(((seq_len(4500) * 37) %% 97) / 18) * 3), 500L, 9L)
  rownames(counts) <- paste0("gene_", 1:500)
  colnames(counts) <- colnames(labels)
  y <- edgeR::normLibSizes(edgeR::DGEList(counts), method = "TMM")
  for (partner_id in c(20L, 1281L)) {
    left <- fit_genes(y, labels[1L, ])$result
    right <- fit_genes(y, labels[partner_id, ])$result
    name <- if (partner_id == 20L) "intervention" else "local_disease"
    left <- left[contrast == name]
    right <- right[contrast == name]
    need(identical(left$statistic, right$statistic) && identical(left$p_value, right$p_value) &&
      identical(left$q_value, right$q_value) && all(abs(left$effect + right$effect) <= 1e-12),
      "Synthetic actual-fit symmetry and restored scientific direction")
  }
  list(status = "PASS_SYNTHETIC_CANONICAL_MAPPING_AND_FITS", allocation_n = 1680L,
    contrast_partner_checks = 3360L, companion_n = c(20L, 20L), synthetic_only = TRUE)
}
load_state <- function(run) {
  prepared <- jsonlite::fromJSON(file.path(run, "PREPARED.json"), simplifyVector = FALSE)
  for (item in list(prepared$derived, prepared$worker, prepared$definition, prepared$input_manifest)) verify(item)
  manifest <- data.table::fread(prepared$input_manifest$path)
  for (i in seq_len(nrow(manifest))) verify(as.list(manifest[i]))
  state <- readRDS(prepared$derived$path)
  list(state = state, original_manifest = data.table::fread(state$original_manifest$path), prepared = prepared)
}
prepare <- function(run, original_run) {
  need(!file.exists(file.path(run, "PREPARED.json")), "Existing preparation retained")
  need(jsonlite::fromJSON(file.path(run, "SYNTHETIC_TESTS.json"))$status ==
    "PASS_SYNTHETIC_CANONICAL_MAPPING_AND_FITS", "Synthetic tests required")
  original_run <- normalizePath(original_run, mustWork = TRUE)
  previous <- jsonlite::fromJSON(file.path(original_run, "PREPARED.json"), simplifyVector = FALSE)
  need(pin(file.path(original_run, "PREPARED.json"))$sha256 == "6e645e7c466bf9b5845e608b7f20a5c7fba3924183fbf4b243a0c48137e6dd1f" &&
    pin(file.path(original_run, "ALLOCATION_MANIFEST.tsv"))$sha256 == "20eab7838ff4368f6b03ef81c1d043483da9bc7be5fcb1a572062f01a8a62ea2" &&
    pin(file.path(original_run, "ENUMERATION_COMPLETE.json"))$sha256 == "3dafd55e27f2649b12d509baf2acb5317c36d764990697fb0733ee83eda5ad44", "Pinned completed original node")
  for (item in list(previous$derived, previous$observed, previous$worker, previous$definition,
    previous$input_manifest, previous$observed_checkpoint)) verify(item)
  inputs <- data.table::fread(previous$input_manifest$path)
  for (i in seq_len(nrow(inputs))) verify(as.list(inputs[i]))
  extras <- list(original_prepared = file.path(original_run, "PREPARED.json"),
    original_derived = previous$derived$path,
    original_allocation_manifest = file.path(original_run, "ALLOCATION_MANIFEST.tsv"),
    original_completion = file.path(original_run, "ENUMERATION_COMPLETE.json"),
    canonical_worker = file.path(run, "worker.R"), canonical_definition = file.path(run, "definition.md"),
    unchanged_summary_source = file.path(run, "summary.R"), unchanged_checker_source = file.path(run, "checker.R"))
  extra_rows <- data.table::rbindlist(lapply(names(extras), function(role) {
    value <- pin(extras[[role]])
    data.table::data.table(role = role, path = value$path, bytes = value$bytes, sha256 = value$sha256)
  }))
  inputs <- data.table::rbindlist(list(inputs, extra_rows))
  need(!anyDuplicated(inputs$path), "Unique input identity inventory")
  need(pin(file.path(run, "summary.R"))$sha256 == "b76bd6a9b7293976fe30328e7f8da578ff3a2647ff9e224c64c765b6f16bdffe" &&
    pin(file.path(run, "checker.R"))$sha256 == "5652292228d881d5e5e225dff68475a8f16ba728eb7f8101065df84c0d854bdb", "Unchanged summary and independent arithmetic")
  state <- readRDS(previous$derived$path)
  need(identical(state$labels, all_labels()), "Original fixed GSM order and full space")
  check_mapping(state$labels)
  state$original_source_binding <- state$source_binding
  state$original_run <- original_run
  state$original_manifest <- pin(file.path(original_run, "ALLOCATION_MANIFEST.tsv"))
  state$source_binding <- digest::digest(list(inputs = inputs, previous_source_binding = state$original_source_binding,
    numerical_method = "canonical_per_allocation_independent_refit_v1"), algo = "sha256")
  original_manifest <- data.table::fread(state$original_manifest$path)
  need(identical(original_manifest$allocation_id, 1:1680), "Complete original manifest")
  value <- assemble(state, original_manifest, state$spaces$observed)
  original_observed <- read_original(state, original_manifest, state$spaces$observed)$value
  bridge <- observed_bridge(value, original_observed)
  dir.create(file.path(run, "allocations"), showWarnings = FALSE)
  observed_path <- save_allocation(value, run)
  twrite(inputs, file.path(run, "INPUT_MANIFEST.tsv"))
  twrite(value$result, file.path(run, "OBSERVED.tsv"))
  twrite(bridge, file.path(run, "OBSERVED_COMPATIBILITY.tsv"))
  for (name in c("ALLOCATIONS.tsv", "FILTER_INVARIANCE.tsv")) {
    need(file.copy(file.path(original_run, name), file.path(run, name), overwrite = FALSE), "Preserved allocation/filter proof")
  }
  saveRDS(state, file.path(run, "derived_inputs.rds"), version = 3)
  capture.output(sessionInfo(), file = file.path(run, "sessionInfo.txt"))
  prepared <- previous
  prepared$status <- "PREPARED_OBSERVED_NATIVE_COMPATIBLE"
  prepared$derived <- pin(file.path(run, "derived_inputs.rds"))
  prepared$observed <- pin(file.path(run, "OBSERVED.tsv"))
  prepared$worker <- pin(file.path(run, "worker.R"))
  prepared$definition <- pin(file.path(run, "definition.md"))
  prepared$input_manifest <- pin(file.path(run, "INPUT_MANIFEST.tsv"))
  prepared$observed_checkpoint <- pin(observed_path)
  prepared$completed_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  prepared$observed_seconds <- value$elapsed_seconds
  prepared$original_preparation <- pin(file.path(original_run, "PREPARED.json"))
  prepared$canonical_gene_independent_refit <- TRUE
  prepared$same_label_non_gene_checkpoint_reuse <- TRUE
  for (i in seq_len(nrow(inputs))) verify(as.list(inputs[i]))
  jwrite(prepared, file.path(run, "PREPARED.json"))
  cat("PREPARED_CANONICAL_OBSERVED_COMPATIBLE; seconds", value$elapsed_seconds, "\n")
}
run_worker <- function(run, worker_id, worker_n) {
  need(worker_n == 4L && worker_id %in% 1:4, "Four single-thread workers")
  loaded <- load_state(run)
  state <- loaded$state
  started_path <- file.path(run, paste0("CANONICAL_WORKER_", worker_id, "_STARTED.json"))
  if (!file.exists(started_path)) jwrite(list(status = "CANONICAL_GENE_REFIT_STARTED", worker_id = worker_id,
    pid = Sys.getpid(), worker = loaded$prepared$worker, derived = loaded$prepared$derived,
    command_args = commandArgs(), time_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")), started_path)
  for (id in state$spaces$global[(state$spaces$global - 1L) %% worker_n == worker_id - 1L]) {
    path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", id))
    if (file.exists(path)) {
      value <- readRDS(path)
      need(identical(value$id, id) && identical(value$source_binding, state$source_binding) &&
        identical(value$labels, state$labels[id, ]), "Valid existing same-source checkpoint")
      validate_result(value$result, state)
      next
    }
    value <- assemble(state, loaded$original_manifest, id)
    save_allocation(value, run)
    if (id %% 100L < 4L) cat("completed", id, "of1680; seconds", value$elapsed_seconds, "\n")
  }
  jwrite(list(status = "CANONICAL_WORKER_COMPLETE", worker_id = worker_id,
    worker = loaded$prepared$worker, derived = loaded$prepared$derived,
    time_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
    file.path(run, paste0("CANONICAL_WORKER_", worker_id, "_COMPLETE.json")))
  cat("CANONICAL_WORKER_COMPLETE", worker_id, "\n")
}
args <- commandArgs(trailingOnly = TRUE)
need(length(args) >= 2L && args[1L] %in% c("selftest", "prepare", "run"), "Mode and explicit run required")
need(getRversion() == "4.6.1" && as.character(utils::packageVersion("edgeR")) == "4.10.1" &&
  as.character(utils::packageVersion("limma")) == "3.68.4", "Pinned native runtime")
data.table::setDTthreads(1L)
set.seed(20260817)
run <- normalizePath(args[2L], mustWork = TRUE)
if (args[1L] == "selftest") {
  jwrite(selftest(), file.path(run, "SYNTHETIC_TESTS.json"))
  cat("PASS_SYNTHETIC_CANONICAL_MAPPING_AND_FITS\n")
} else if (args[1L] == "prepare") {
  need(length(args) == 3L, "Prepare requires original completed node")
  prepare(run, args[3L])
} else {
  need(length(args) == 4L, "Run requires worker ID/count")
  run_worker(run, as.integer(args[3L]), as.integer(args[4L]))
}
