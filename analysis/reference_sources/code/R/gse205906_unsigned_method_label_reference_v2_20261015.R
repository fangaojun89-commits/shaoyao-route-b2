#!/usr/bin/env Rscript
# Full labelled unsigned-method reference v2. No historical driver is sourced.
# Preserve frozen v1 (0bac4621...). Canonical zero coefficients may carry -0
# after multiplying by orientation. Normalize only exact zeros on the ACTUAL
# mroast target vector before both hashing and calling: no numeric coefficient,
# hypothesis, RNG state, native parameter or nonzero value changes.
# CLI: selftest RUN; static-selftest RUN; prepare RUN CANONICAL_RUN NATIVE_MODEL;
#      run RUN WORKER_ID WORKER_N (1..4). Main freezes/reviews before real prepare/run.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)
library(data.table)
library(digest)
library(jsonlite)
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
pin <- function(path) list(path = normalizePath(path, mustWork = TRUE),
  bytes = unname(file.info(path)$size), sha256 = digest(file = path, algo = "sha256", serialize = FALSE))
verify <- function(expected) {
  actual <- pin(expected$path)
  need(actual$bytes == expected$bytes && actual$sha256 == expected$sha256, "Frozen input identity")
  invisible(actual)
}
jwrite <- function(value, path) {
  need(!file.exists(path), paste("Preserve existing output", basename(path)))
  write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
}
twrite <- function(value, path) {
  need(!file.exists(path), paste("Preserve existing output", basename(path)))
  transported <- copy(value)
  for (field in names(transported)) if (is.list(transported[[field]])) {
    transported[[field]] <- vapply(transported[[field]], function(z)
      as.character(toJSON(z, auto_unbox = FALSE, na = "null")), character(1))
  }
  fwrite(transported, path, sep = "\t", quote = TRUE, na = "NA")
}
save_new_rds <- function(value, path) {
  need(!file.exists(path) && !file.exists(paste0(path, ".partial")), "Preserve existing RDS/checkpoint")
  saveRDS(value, paste0(path, ".partial"), version = 3)
  need(file.rename(paste0(path, ".partial"), path), "Atomic new checkpoint")
  pin(path)
}
groups <- c("sham", "SNI", "PLM")
contrasts <- cbind(local_disease = c(-1, 1, 0), intervention = c(0, -1, 1))
rownames(contrasts) <- groups
collections <- c("Reactome", "SynGO_sensitivity")
key_columns <- c("collection", "method", "contrast", "feature")
numeric_columns <- c("effect", "statistic", "p_value", "q_value")
seed_base <- 20260817L
claim_role <- "EXPLORATORY_INJURED_SPINAL_CONTEXT_ONLY"
versions <- c(edgeR = "4.10.1", limma = "3.68.4", fgsea = "1.38.0",
  BiocParallel = "1.46.0", data.table = "1.18.4", digest = "0.6.39", jsonlite = "2.0.0")

peak_memory_kb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  value <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  if (length(value) == 1L) as.numeric(gsub("[^0-9]", "", value)) else NA_real_
}
backend_snapshot <- function() {
  backend <- BiocParallel::bpparam()
  list(class = class(backend), workers = BiocParallel::bpnworkers(backend),
    started = BiocParallel::bpisup(backend), RNGseed = BiocParallel::bpRNGseed(backend))
}
kernel_contract <- function() {
  roast <- get("mroast.default", envir = asNamespace("limma"))
  effects <- get(".roastEffects", envir = asNamespace("limma"))
  fg <- get("fgseaMultilevel", envir = asNamespace("fgsea"))
  structure <- function(fun) list(formals = formals(fun), body = body(fun))
  list(R = R.version.string, packages = as.list(versions), RNG_kind = RNGkind(),
    mroast_function = structure(roast), roast_effects_function = structure(effects),
    fgsea_function = structure(fg), backend = backend_snapshot(),
    mroast_nrot = 9999L, mroast_set_statistic = "mean", fgsea_eps = 0,
    minSize = 10L, maxSize = 500L,
    unchanged_implicit_defaults = TRUE, orientation_excluded_from_stochastic_hash = TRUE,
    exact_zero_mroast_target_coefficients_use_positive_zero_bits = TRUE)
}
runtime_snapshot <- function() {
  list(time_utc = utc(), R = R.version.string, packages = as.list(versions),
    library_paths = .libPaths(), platform = R.version$platform, locale = Sys.getlocale(),
    RNG_kind = RNGkind(), base_seed = seed_base, backend = backend_snapshot(),
    thread_environment = as.list(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"))),
    data_table_threads = getDTthreads(), command = commandArgs(), working_directory = getwd(),
    kernel_contract_sha256 = digest(kernel_contract(), algo = "sha256"),
    renv_restored_here = FALSE, Python_used = FALSE, Git_commit = NULL,
    historical_runtime_authenticated = FALSE, maximum_workers = 4L)
}
import_canonical <- function(path) {
  need(pin(path)$sha256 == "18eaab06d0bb980dae3373279c83d0e4e32f94dbee346ef4e2cc120e11ba3965", "Accepted canonical source SHA")
  expressions <- as.list(parse(path, keep.source = FALSE))
  selected <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    is.symbol(x[[2L]]) && as.character(x[[2L]]) == "canonical_spec" &&
    is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), expressions)
  need(length(selected) == 1L, "Import only the canonical function declaration, never its driver")
  scope <- new.env(parent = globalenv())
  eval(selected[[1L]], scope)
  scope$canonical_spec
}
check_space <- function(labels, spaces, canonical) {
  need(identical(dim(labels), c(1680L, 9L)) &&
    identical(colnames(labels), paste0("GSM", 6235292:6235300)), "Fixed nine GSM/allocation geometry")
  signatures <- apply(labels, 1L, paste, collapse = ":")
  need(!anyDuplicated(signatures) && all(rowSums(labels == "sham") == 3L) &&
    all(rowSums(labels == "SNI") == 3L) && all(rowSums(labels == "PLM") == 3L) &&
    identical(as.integer(spaces$global), seq_len(1680L)) && spaces$observed == 1L &&
    identical(unname(labels[1L, ]), rep(groups, each = 3L)), "Complete balanced space, observed once")
  companions <- list(local_disease = which(rowSums(labels[, 7:9, drop = FALSE] == "PLM") == 3L),
    intervention = which(rowSums(labels[, 1:3, drop = FALSE] == "sham") == 3L))
  for (name in names(companions)) need(length(companions[[name]]) == 20L &&
    identical(as.integer(companions[[name]]), as.integer(spaces[[name]])), "Full companion spaces")
  for (id in spaces$global) {
    left <- canonical(labels[id, ])
    for (name in colnames(contrasts)) {
      pair <- groups[contrasts[, name] != 0]
      swapped <- labels[id, ]
      swapped[labels[id, ] == pair[1L]] <- pair[2L]
      swapped[labels[id, ] == pair[2L]] <- pair[1L]
      partner <- match(paste(swapped, collapse = ":"), signatures)
      need(!is.na(partner), "Full label-swap closure")
      right <- canonical(labels[partner, ])
      need(identical(left$design, right$design) &&
        identical(left$contrasts[, name], right$contrasts[, name]) &&
        left$orientation[name] == -right$orientation[name], "Canonical target input/direction identity")
    }
  }
  invisible(TRUE)
}

# Rebuild the original memberships using the separately pinned orthology source.
# Keep duplicate-symbol gene rows: eligibility counts unique mapped symbols only.
map_resources <- function(resource, orthology, symbols) {
  need(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in% names(orthology)), "Original orthology schema")
  map <- unique(orthology[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
    .(human_symbol = as.character(human_symbol), rat_symbol = as.character(other_symbol))])
  need(!anyDuplicated(map$human_symbol) && !anyDuplicated(map$rat_symbol), "Strict one-to-one orthology")
  sources <- list(Reactome = resource$reactome, SynGO_sensitivity = resource$syngo_sensitivity)
  lapply(names(sources), function(collection) {
    x <- as.data.table(sources[[collection]])
    need(all(c("term_id", "term_name", "gene_symbol") %in% names(x)), "Original unsigned resource schema")
    x <- unique(x[, .(program_id = as.character(term_id), program_name = as.character(term_name),
      human_symbol = as.character(gene_symbol))])
    mapped <- unique(merge(x, map, by = "human_symbol", allow.cartesian = TRUE)[,
      .(program_id, program_name, symbol = rat_symbol)])
    members <- lapply(split(mapped$symbol, mapped$program_id), function(z) unique(intersect(z, symbols)))
    members <- members[lengths(members) >= 10L & lengths(members) <= 500L]
    list(collection = collection, members = members,
      indices = lapply(members, function(z) which(symbols %in% z)), names = unique(mapped[, .(program_id, program_name)]))
  }) |> setNames(names(sources))
}
add_flags <- function(value) {
  value[, `:=`(analysis_role = claim_role, candidate_gate_eligible = FALSE,
    can_generate_or_upgrade_candidate = FALSE, formal_sequential_inference = FALSE)]
  value
}
capture_call <- function(fun) {
  conditions <- new.env(parent = emptyenv())
  conditions$messages <- conditions$warnings <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(fun(),
    message = function(x) { conditions$messages <- c(conditions$messages, conditionMessage(x)); invokeRestart("muffleMessage") },
    warning = function(x) { conditions$warnings <- c(conditions$warnings, conditionMessage(x)); invokeRestart("muffleWarning") }),
    error = function(x) list(failed = TRUE, class = class(x),
      text_sha256 = digest(conditionMessage(x), algo = "sha256")))
  if (is.list(value) && isTRUE(value$failed) && !is.data.frame(value)) { error <- value; value <- NULL }
  list(value = value, diagnostic = list(message_n = length(conditions$messages),
    warning_n = length(conditions$warnings), condition_text_sha256 = digest(
      c(conditions$messages, conditions$warnings), algo = "sha256"), error = error))
}
stochastic_record <- function(method, collection, contrast, seed, arguments, contract) {
  set.seed(seed)
  initial <- .Random.seed
  payload <- list(method = method, collection = collection, contrast = contrast,
    actual_arguments = arguments, seed = seed, RNG_kind = RNGkind(), RNG_state = initial,
    backend = backend_snapshot(), kernel_contract = contract)
  list(method = method, collection = collection, contrast = contrast, seed = seed,
    full_stochastic_input_sha256 = digest(payload, algo = "sha256"),
    argument_sha256 = digest(arguments, algo = "sha256"), RNG_kind = RNGkind(),
    RNG_state = initial, backend = payload$backend,
    orientation_included = FALSE, component_sha256 = lapply(arguments, digest, algo = "sha256"))
}
native_columns <- function(method, decorated = TRUE) {
  columns <- if (method == "mroast") c("program_id", "NGenes", "PropDown", "PropUp", "Direction",
    "PValue", "FDR", "PValue.Mixed", "FDR.Mixed") else c("program_id", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
  if (decorated) columns <- c(columns, "collection", "contrast", "method", "p_value",
    "q_within_collection_contrast", "analysis_role", "candidate_gate_eligible",
    "can_generate_or_upgrade_candidate", "formal_sequential_inference")
  columns
}
finite_diagnostic <- function(value, resource, method, collection, contrast, ranks = NULL) {
  expected_programs <- names(resource$members)
  fields <- if (method == "mroast") c("NGenes", "PropDown", "PropUp", "PValue", "FDR", "PValue.Mixed", "FDR.Mixed",
    "p_value", "q_within_collection_contrast") else c("pval", "padj", "log2err", "ES", "NES", "size", "p_value", "q_within_collection_contrast")
  rows <- lapply(fields, function(field) data.table(field = field,
    values = if (is.null(value)) 0L else length(value[[field]]),
    missing_n = if (is.null(value)) NA_integer_ else sum(is.na(value[[field]])),
    nonfinite_n = if (is.null(value)) NA_integer_ else sum(!is.finite(value[[field]]))))
  probability_fields <- if (method == "mroast") c("PValue", "FDR", "PValue.Mixed", "FDR.Mixed",
    "p_value", "q_within_collection_contrast") else c("pval", "padj", "p_value", "q_within_collection_contrast")
  invalid_probability_n <- if (is.null(value)) NA_integer_ else sum(vapply(probability_fields,
    function(field) sum(!is.finite(value[[field]]) | value[[field]] <= 0 | value[[field]] > 1), integer(1)))
  text_fields <- c("program_id", "collection", "contrast", "method", "analysis_role",
    if (method == "mroast") "Direction")
  text_missing_n <- if (is.null(value)) NA_integer_ else sum(vapply(text_fields,
    function(field) sum(is.na(value[[field]]) | !nzchar(value[[field]])), integer(1)))
  leading_edge_missing_n <- if (method == "fgsea" && !is.null(value))
    sum(is.na(unlist(value$leadingEdge, use.names = FALSE))) else 0L
  expected_sizes <- if (method == "mroast") lengths(resource$indices) else
    vapply(resource$members, function(z) length(unique(intersect(z, names(ranks)))), integer(1))
  size_field <- if (method == "mroast") "NGenes" else "size"
  size_mismatch_n <- if (is.null(value) || !all(c("program_id", size_field) %in% names(value))) NA_integer_ else {
    expected_size <- expected_sizes[match(value$program_id, names(expected_sizes))]
    sum(is.na(expected_size) | !is.finite(value[[size_field]]) | value[[size_field]] != expected_size)
  }
  required <- native_columns(method)
  column_length_mismatch_n <- if (is.null(value)) length(required) else sum(vapply(required,
    function(field) length(value[[field]]) != nrow(value), logical(1)))
  direction_invalid_n <- if (method == "mroast" && !is.null(value))
    sum(is.na(value$Direction) | !value$Direction %in% c("Up", "Down")) else 0L
  prop_invalid_n <- if (method == "mroast" && !is.null(value)) sum(vapply(c("PropUp", "PropDown"),
    function(field) sum(!is.finite(value[[field]]) | value[[field]] < 0 | value[[field]] > 1), integer(1))) else 0L
  leading_edge_invalid_n <- if (method == "fgsea" && !is.null(value)) {
    if (!is.list(value$leadingEdge) || length(value$leadingEdge) != nrow(value)) nrow(value) + 1L else
      sum(!vapply(seq_len(nrow(value)), function(i) {
        edge <- value$leadingEdge[[i]]; id <- value$program_id[i]
        is.character(edge) && !anyNA(edge) && !anyDuplicated(edge) &&
          id %in% names(resource$members) && all(edge %in% intersect(resource$members[[id]], names(ranks)))
      }, logical(1)))
  } else 0L
  list(method = method, collection = collection, contrast = contrast,
    expected_program_n = length(expected_programs), returned_program_n = if (is.null(value)) 0L else nrow(value),
    absent_program_ids = if (is.null(value)) expected_programs else setdiff(expected_programs, value$program_id),
    unexpected_program_ids = if (is.null(value)) character() else setdiff(value$program_id, expected_programs),
    duplicated_program_n = if (is.null(value)) 0L else sum(duplicated(value$program_id)),
    missing_required_columns = setdiff(native_columns(method), names(value)),
    invalid_probability_n = invalid_probability_n, text_missing_n = text_missing_n,
    leading_edge_missing_n = leading_edge_missing_n, engine_size_mismatch_n = size_mismatch_n,
    required_column_length_mismatch_n = column_length_mismatch_n,
    direction_invalid_n = direction_invalid_n, proportion_invalid_n = prop_invalid_n,
    leading_edge_invalid_n = leading_edge_invalid_n,
    engine_size_rule = if (method == "mroast") "original_gene_row_indices_not_unique_symbol_eligibility" else
      "unique_original_members_intersect_actual_finite_rank_names", fields = rbindlist(rows))
}

fit_unsigned <- function(state, labels, canonical) {
  started <- proc.time()[[3L]]
  spec <- canonical(labels)
  need(identical(rownames(spec$design), colnames(state$y)), "Fixed sample order")
  v <- limma::voom(state$y, spec$design, plot = FALSE)
  contract <- kernel_contract()
  native <- list(mroast = list(), fgsea = list()); signed <- list()
  records <- diagnostics <- selectors <- list(); result_rows <- list()
  for (k in seq_along(collections)) {
    collection <- collections[k]; resource <- state$resources[[collection]]
    for (j in seq_len(ncol(contrasts))) {
      name <- colnames(contrasts)[j]; orientation <- spec$orientation[name]
      vector <- spec$contrasts[, name]; key <- paste(collection, name, sep = "::")
      # canonical_spec's zero * negative orientation is IEEE -0. Do not merely
      # hide its sign bit in a hash: make the actual kernel vector bit-identical.
      original_vector <- vector
      vector[vector == 0] <- 0
      need(all(vector == original_vector) && identical(names(vector), names(original_vector)),
        "Only exact-zero representation normalized, never any contrast value")
      rotation_seed <- seed_base + j + 100L * k
      roast_arguments <- list(E = v$E, indices = resource$indices, design = spec$design,
        contrast = vector, nrot = 9999L, set.statistic = "mean")
      roast_record <- stochastic_record("mroast", collection, name, rotation_seed, roast_arguments, contract)
      roast_call <- capture_call(function() as.data.table(limma::mroast(v$E, resource$indices,
        spec$design, contrast = vector, nrot = 9999L, set.statistic = "mean"), keep.rownames = "program_id"))
      roast <- roast_call$value
      if (!is.null(roast) && all(native_columns("mroast", FALSE) %in% names(roast))) {
        roast[, `:=`(collection = collection, contrast = name, method = "mroast_mean_9999_rotations", p_value = PValue)]
        roast[, q_within_collection_contrast := p.adjust(p_value, "BH")]
        roast <- add_flags(roast)
        result_rows[[length(result_rows) + 1L]] <- roast[, .(collection, method, contrast,
          feature = program_id, effect = (PropUp - PropDown) * orientation,
          statistic = -log(PValue), p_value, q_value = q_within_collection_contrast)]
      }
      native$mroast[[key]] <- roast
      records[[paste0(key, "::mroast")]] <- roast_record
      diagnostics[[paste0(key, "::mroast")]] <- list(call = roast_call$diagnostic,
        finite = finite_diagnostic(roast, resource, "mroast", collection, name))
      # Literal original weighted t/symbol selector, independently recomputed for this label/collection.
      weighted <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, spec$design),
        spec$contrasts[, name, drop = FALSE]), robust = TRUE)
      stat_table <- data.table(symbol = state$symbols, stat = as.numeric(weighted$t[, 1L]),
        gene_id = rownames(state$y), original_row = seq_along(state$symbols))
      stat_table <- stat_table[is.finite(stat) & !is.na(symbol) & nzchar(symbol)]
      stat_table <- stat_table[order(-abs(stat)), .SD[1L], by = symbol]
      ranks <- stat_table$stat; names(ranks) <- stat_table$symbol
      selectors[[key]] <- stat_table
      fg_seed <- seed_base + 1000L + j + 100L * k
      fg_arguments <- list(pathways = resource$members, stats = ranks,
        minSize = 10L, maxSize = 500L, eps = 0)
      fg_record <- stochastic_record("fgsea", collection, name, fg_seed, fg_arguments, contract)
      fg_call <- capture_call(function() as.data.table(fgsea::fgseaMultilevel(
        pathways = resource$members, stats = ranks, minSize = 10L, maxSize = 500L, eps = 0)))
      fg <- fg_call$value
      if (!is.null(fg) && "pathway" %in% names(fg)) setnames(fg, "pathway", "program_id")
      if (!is.null(fg) && all(native_columns("fgsea", FALSE) %in% names(fg))) {
        fg[, `:=`(collection = collection, contrast = name, method = "fgseaMultilevel",
          p_value = pval, q_within_collection_contrast = padj)]
        fg <- add_flags(fg)
        result_rows[[length(result_rows) + 1L]] <- fg[, .(collection, method, contrast,
          feature = program_id, effect = NES * orientation, statistic = abs(NES),
          p_value, q_value = q_within_collection_contrast)]
        signed[[key]] <- fg[, .(program_id, collection, contrast,
          derived_scientific_ES = ES * orientation, derived_scientific_NES = NES * orientation)]
      }
      native$fgsea[[key]] <- fg
      records[[paste0(key, "::fgsea")]] <- fg_record
      diagnostics[[paste0(key, "::fgsea")]] <- list(call = fg_call$diagnostic,
        finite = finite_diagnostic(fg, resource, "fgsea", collection, name, ranks))
    }
  }
  result <- if (length(result_rows)) rbindlist(result_rows) else data.table(collection = character(),
    method = character(), contrast = character(), feature = character(), effect = numeric(),
    statistic = numeric(), p_value = numeric(), q_value = numeric())
  setorderv(result, key_columns)
  list(result = result, native_payload = lapply(native, rbindlist),
    derived_signed_fgsea = rbindlist(signed), rank_selectors = selectors,
    stochastic_inputs = records, diagnostics = diagnostics,
    canonical_specification = spec, elapsed_seconds = proc.time()[[3L]] - started,
    kernel_contract_sha256 = digest(contract, algo = "sha256"),
    method_calls = c(mroast = 4L, fgsea = 4L), no_partition_fit_cache = TRUE,
    canonical_payload_Direction_and_leadingEdge_not_reversed = TRUE)
}
validate_result <- function(value, state, id) {
  result <- value$result
  need(identical(value$id, id) && identical(value$labels, state$labels[id, ]) &&
    identical(value$source_binding, state$source_binding), "Checkpoint ID/labels/source")
  need(identical(names(result), c(key_columns, numeric_columns)) &&
    identical(as.data.frame(result[, ..key_columns]), state$expected_keys), "Eight complete ordered native families")
  for (field in numeric_columns) need(all(is.finite(result[[field]])), "Nonfinite inference retained, full acceptance prevented")
  need(all(result$statistic >= 0 & result$p_value > 0 & result$p_value <= 1 &
    result$q_value > 0 & result$q_value <= 1), "Original statistic/P/q boundaries")
  families <- split(seq_len(nrow(result)), paste(result$collection, result$method, result$contrast))
  for (index in families) need(max(abs(p.adjust(result$p_value[index], "BH") - result$q_value[index])) <= 1e-12,
    "Project raw-P BH in each complete collection/method/contrast family")
  need(length(value$stochastic_inputs) == 8L && length(value$diagnostics) == 8L, "All stochastic call/diagnostic records")
  for (entry in value$diagnostics) {
    need(is.null(entry$call$error) && entry$finite$returned_program_n == entry$finite$expected_program_n &&
      !length(entry$finite$absent_program_ids) && !length(entry$finite$unexpected_program_ids) &&
      entry$finite$duplicated_program_n == 0L &&
      !length(entry$finite$missing_required_columns) &&
      entry$finite$invalid_probability_n == 0L &&
      entry$finite$text_missing_n == 0L && entry$finite$leading_edge_missing_n == 0L &&
      entry$finite$engine_size_mismatch_n == 0L &&
      entry$finite$required_column_length_mismatch_n == 0L &&
      entry$finite$direction_invalid_n == 0L && entry$finite$proportion_invalid_n == 0L &&
      entry$finite$leading_edge_invalid_n == 0L &&
      all(entry$finite$fields$values == entry$finite$expected_program_n) &&
      all(entry$finite$fields$missing_n == 0L & entry$finite$fields$nonfinite_n == 0L),
      "All native fields finite; no skipped programs, fill or native fallback")
  }
  # Resume validation reads the actual native payload, not merely its saved
  # diagnostic summary. No native method is fitted or stochastic call repeated.
  for (native_method in c("mroast", "fgsea")) {
    payload <- value$native_payload[[native_method]]
    expected_method <- if (native_method == "mroast") "mroast_mean_9999_rotations" else "fgseaMultilevel"
    need(is.data.table(payload) && all(native_columns(native_method) %in% names(payload)) &&
      nrow(payload) == 2L * sum(lengths(lapply(state$resources, `[[`, "members"))) &&
      all(payload$collection %in% collections & payload$contrast %in% colnames(contrasts)) &&
      all(payload$method == expected_method) && all(payload$analysis_role == claim_role) &&
      all(!payload$candidate_gate_eligible & !payload$can_generate_or_upgrade_candidate &
        !payload$formal_sequential_inference), "Actual complete native payload schema/flags")
    for (target_collection in collections) for (target_contrast in colnames(contrasts)) {
      key <- paste(target_collection, target_contrast, sep = "::")
      selected <- payload[collection == target_collection & contrast == target_contrast]
      selector <- value$rank_selectors[[key]]
      need(is.data.table(selector) && all(c("symbol", "stat", "gene_id", "original_row") %in% names(selector)) &&
        !anyNA(selector) && !anyDuplicated(selector$symbol) && all(is.finite(selector$stat)) &&
        all(selector$original_row %in% seq_along(state$symbols)) &&
        identical(selector$gene_id, rownames(state$y)[selector$original_row]) &&
        identical(selector$symbol, state$symbols[selector$original_row]), "Saved original-row symbol/rank selector alignment")
      ranks <- selector$stat; names(ranks) <- selector$symbol
      fresh <- finite_diagnostic(selected, state$resources[[target_collection]], native_method,
        target_collection, target_contrast, ranks)
      need(all(selected$p_value == if (native_method == "mroast") selected$PValue else selected$pval) &&
        all(selected$q_within_collection_contrast == if (native_method == "mroast")
          p.adjust(selected$PValue, "BH") else selected$padj), "Native/project aliases retain raw primary P and separate mid-P FDR")
      need(fresh$returned_program_n == fresh$expected_program_n &&
        !length(fresh$absent_program_ids) && !length(fresh$unexpected_program_ids) &&
        fresh$duplicated_program_n == 0L && !length(fresh$missing_required_columns) &&
        fresh$invalid_probability_n == 0L && fresh$text_missing_n == 0L &&
        fresh$leading_edge_missing_n == 0L && fresh$engine_size_mismatch_n == 0L &&
        fresh$required_column_length_mismatch_n == 0L && fresh$direction_invalid_n == 0L &&
        fresh$proportion_invalid_n == 0L && fresh$leading_edge_invalid_n == 0L &&
        all(fresh$fields$values == fresh$expected_program_n & fresh$fields$missing_n == 0L &
          fresh$fields$nonfinite_n == 0L), "Recomputed native payload finite/text/list/domain acceptance")
      orientation <- value$canonical_specification$orientation[target_contrast]
      linked <- if (native_method == "mroast") selected[, .(feature = program_id,
        effect = (PropUp - PropDown) * orientation, statistic = -log(PValue),
        p_value = PValue, q_value = q_within_collection_contrast)] else selected[, .(
          feature = program_id, effect = NES * orientation, statistic = abs(NES),
          p_value = pval, q_value = padj)]
      common <- result[collection == target_collection & contrast == target_contrast & method == expected_method]
      setorderv(linked, "feature"); setorderv(common, "feature")
      need(identical(linked$feature, common$feature), "Actual native/common program alignment")
      for (field in numeric_columns) need(all(linked[[field]] == common[[field]]),
        "Actual native payload deterministically links to outer fields; no reverse-native fabrication")
    }
  }
  invisible(TRUE)
}

# Full observed payload bridge. The original targets remain read-only. Numeric
# transport uses the frozen original abs+relative contract, not equality of RAM
# diagnostic maxima and later TSV-readback maxima. Ordered leadingEdge is text/list checked.
observed_bridge <- function(value, targets) {
  rows <- list()
  for (method in c("mroast", "fgsea")) {
    expected <- fread(targets[[method]], na.strings = "NA")
    actual <- copy(value$native_payload[[method]])
    keys <- c("collection", "contrast", "program_id")
    setorderv(expected, keys); setorderv(actual, keys)
    need(setequal(names(expected), names(actual)) && nrow(expected) == 3248L &&
      identical(as.data.frame(expected[, ..keys]), as.data.frame(actual[, ..keys])) &&
      !anyDuplicated(expected[, ..keys]), "Complete observed native payload schema/keys")
    for (field in names(expected)) {
      left <- actual[[field]]; right <- expected[[field]]
      if (is.list(left)) {
        # The checked production writer transports leadingEdge as a JSON array.
        need(field == "leadingEdge" && !anyNA(right), "Original ordered leadingEdge JSON transport")
        list_equal <- all(vapply(seq_along(left), function(i) identical(as.character(left[[i]]),
          as.character(fromJSON(right[[i]], simplifyVector = TRUE))), logical(1)))
        rows[[length(rows) + 1L]] <- data.table(method = method, field = field,
          comparisons = length(right), exact_NA_mask = !anyNA(unlist(left)), numeric = FALSE,
          maximum_absolute_difference = NA_real_, maximum_fraction_of_original_bound = NA_real_,
          accepted = list_equal && !anyNA(unlist(left)))
        next
      }
      mask <- identical(is.na(left), is.na(right))
      if (is.numeric(right)) {
        keep <- !is.na(right) & !is.na(left)
        delta <- abs(left[keep] - right[keep])
        bound <- 1e-10 + 1e-8 * pmax(abs(left[keep]), abs(right[keep]))
        accepted <- mask && all(is.finite(left[keep]) & is.finite(right[keep])) && all(delta <= bound)
        maximum <- if (length(delta)) max(delta) else 0
        bound_fraction <- if (length(delta)) max(delta / bound) else 0
      } else {
        accepted <- mask && identical(as.character(left), as.character(right))
        maximum <- bound_fraction <- NA_real_
      }
      rows[[length(rows) + 1L]] <- data.table(method = method, field = field,
        comparisons = length(right), exact_NA_mask = mask, numeric = is.numeric(right),
        maximum_absolute_difference = maximum, maximum_fraction_of_original_bound = bound_fraction,
        accepted = accepted)
    }
  }
  rbindlist(rows)
}
stamp <- function(fit, state, id) {
  fit$id <- id; fit$labels <- state$labels[id, ]; fit$source_binding <- state$source_binding
  fit$seed_base <- seed_base; fit$runtime <- runtime_snapshot()
  fit
}

selftest <- function(run, static_only = FALSE) {
  canonical <- import_canonical(file.path(run, "canonical_source.R"))
  labels <- setNames(rep(groups, each = 3L), paste0("synthetic_sample_", 1:9))
  tied <- list(Direction = "Down", PropUp = 0, PropDown = 0, leadingEdge = list(c("B", "A")))
  need(tied$Direction == "Down" && tied$PropUp - tied$PropDown == 0 &&
    identical(tied$leadingEdge, list(c("B", "A"))), "Canonical tied native payload must not be fabricated as reversed")
  for (name in colnames(contrasts)) {
    pair <- groups[contrasts[, name] != 0]; swapped <- labels
    swapped[labels == pair[1L]] <- pair[2L]; swapped[labels == pair[2L]] <- pair[1L]
    left <- canonical(labels); right <- canonical(swapped)
    need(identical(left$design, right$design) && identical(left$contrasts[, name], right$contrasts[, name]) &&
      left$orientation[name] == -right$orientation[name], "Synthetic coordinate-only coupling")
    actual_left_vector <- left$contrasts[, name]; actual_right_vector <- right$contrasts[, name]
    actual_left_vector[actual_left_vector == 0] <- 0
    actual_right_vector[actual_right_vector == 0] <- 0
    need(identical(actual_left_vector, actual_right_vector, num.eq = FALSE) &&
      identical(digest(actual_left_vector, algo = "sha256"), digest(actual_right_vector, algo = "sha256")),
      "Synthetic actual mroast vector has identical zero bits and full hash")
  }
  if (static_only) return(list(status = "PASS_UNSIGNED_STATIC_ALGEBRA_ONLY", synthetic_only = TRUE,
    method_calls = 0L, real_data_read = FALSE, no_scientific_acceptance = TRUE))
  set.seed(seed_base)
  means <- exp(runif(500L, log(20), log(4000)))
  counts <- matrix(rnbinom(4500L, mu = rep(means, 9L), size = 8), 500L, 9L)
  counts[1:24, 4:6] <- counts[1:24, 4:6] * 2L
  counts[25:48, 7:9] <- counts[25:48, 7:9] * 2L
  rownames(counts) <- paste0("synthetic_gene_", 1:500); colnames(counts) <- names(labels)
  symbols <- paste0("symbol_", 1:500); symbols[499:500] <- symbols[1:2]
  members <- list(positive = symbols[1:24], intervention = symbols[25:48],
    dispersed = symbols[seq(4, 490, by = 11)], background = symbols[100:140])
  resources <- setNames(lapply(collections, function(name) list(collection = name, members = members,
    indices = lapply(members, function(z) which(symbols %in% z)))), collections)
  state <- list(y = edgeR::normLibSizes(edgeR::DGEList(counts), method = "TMM"), symbols = symbols, resources = resources)
  original <- fit_unsigned(state, labels, canonical); checks <- list()
  state$labels <- matrix(labels, nrow = 1L, dimnames = list(NULL, names(labels)))
  state$source_binding <- "synthetic_only_no_scientific_source"
  state$expected_keys <- as.data.frame(original$result[, ..key_columns])
  validate_result(stamp(original, state, 1L), state, 1L)
  for (name in colnames(contrasts)) {
    pair <- groups[contrasts[, name] != 0]; swapped <- labels
    swapped[labels == pair[1L]] <- pair[2L]; swapped[labels == pair[2L]] <- pair[1L]
    other <- fit_unsigned(state, swapped, canonical)
    swapped_state <- state; swapped_state$labels[1L, ] <- swapped
    validate_result(stamp(other, swapped_state, 1L), swapped_state, 1L)
    for (target_collection in collections) for (method in c("mroast", "fgsea")) {
      key <- paste(target_collection, name, method, sep = "::")
      need(identical(original$stochastic_inputs[[key]], other$stochastic_inputs[[key]]), "Complete stochastic inputs, not merely equal seeds")
      a <- original$native_payload[[method]][contrast == name & collection == target_collection]
      b <- other$native_payload[[method]][contrast == name & collection == target_collection]
      need(identical(as.list(a), as.list(b)), "Synthetic canonical native payload P/mid-P/list/MC identity")
      checks[[key]] <- list(stochastic_input_exact = TRUE, complete_native_payload_exact = TRUE)
    }
  }
  list(status = "PASS_SYNTHETIC_UNSIGNED_WORKER_CANONICAL_PAYLOAD_AND_RNG", synthetic_only = TRUE,
    scientific_files_read = FALSE, kernel_calls = c(mroast = 12L, fgsea = 12L), checks = checks,
    nrot = 9999L, eps = 0, no_native_Direction_or_leadingEdge_reconstruction = TRUE,
    runtime = runtime_snapshot())
}

prepare <- function(run, canonical_run, model) {
  need(!file.exists(file.path(run, "PREPARED.json")), "Preserve prior preparation")
  synthetic <- fromJSON(file.path(run, "SYNTHETIC_TESTS.json"), simplifyVector = FALSE)
  need(synthetic$status == "PASS_SYNTHETIC_UNSIGNED_WORKER_CANONICAL_PAYLOAD_AND_RNG" &&
    synthetic$worker$sha256 == pin(file.path(run, "worker.R"))$sha256 &&
    synthetic$canonical_source$sha256 == pin(file.path(run, "canonical_source.R"))$sha256,
    "Worker-specific synthetic source binding and actual kernel tests required")
  canonical_run <- normalizePath(canonical_run, mustWork = TRUE); model <- normalizePath(model, mustWork = TRUE)
  files <- list(canonical_prepared = file.path(canonical_run, "PREPARED.json"),
    canonical_state = file.path(canonical_run, "derived_inputs.rds"),
    canonical_completion = file.path(canonical_run, "ENUMERATION_COMPLETE.json"),
    canonical_filter_proof = file.path(canonical_run, "FILTER_INVARIANCE.tsv"),
    canonical_input_manifest = file.path(canonical_run, "INPUT_MANIFEST.tsv"),
    native_completion = file.path(model, "COMPLETED"), native_manifest = file.path(model, "file_manifest.json"),
    native_mroast = file.path(model, "unsigned_mroast.tsv"), native_fgsea = file.path(model, "unsigned_fgsea.tsv"),
    model_source = file.path(run, "model_source.R"), config = file.path(run, "config.json"),
    canonical_source = file.path(run, "canonical_source.R"), worker = file.path(run, "worker.R"),
    definition = file.path(run, "definition.md"), synthetic_tests = file.path(run, "SYNTHETIC_TESTS.json"),
    summary_source = file.path(run, "summary.R"))
  expected <- c(canonical_prepared = "0631f763b0decbd3d5345770185ba6b45b09fb399a9cd45ba4840b21ba2391b4",
    canonical_state = "7a4405a6c046cdf0ca6ba9b2b2692e0c79d6ae92e6b90eb279a1aa46cd41662f",
    canonical_completion = "3ae21310e0f06d23aa8b1163c8acc423551a62385a17fcb17bbfd108a8ad1015",
    native_completion = "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e",
    native_manifest = "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7",
    native_mroast = "85e92b4bc52dd273903bfc2add70e4d1357e98003935785b1245aae56eac139f",
    native_fgsea = "8687d7f135a4c8464b0cffea2975cc908913b7b093780a211e115b3f417dd2a7",
    model_source = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2",
    config = "ee6e22592b06ecdb31707234e97dc55b7b66cbad8dd11c364187af0624ece980")
  for (role in names(expected)) need(pin(files[[role]])$sha256 == expected[role], paste("Accepted input", role))
  upstream <- fread(files$canonical_input_manifest)
  for (i in seq_len(nrow(upstream))) verify(as.list(upstream[i]))
  resource_path <- upstream$path[upstream$role == "resources"]
  annotation_path <- upstream$path[upstream$role == "gene_edger.tsv"]
  need(length(resource_path) == 1L && length(annotation_path) == 1L, "Unique accepted resource/annotation sources")
  files$resources <- resource_path
  files$orthology <- file.path(dirname(resource_path), "ncbi_human_mouse_rat_one_to_one.tsv")
  original_manifest <- fromJSON(files$native_manifest, simplifyVector = FALSE)
  need(pin(files$orthology)$sha256 == original_manifest$provenance$eleven_input_pins$orthology$sha256,
    "Original independent orthology authority, not embedded replacement")
  files$annotation <- annotation_path
  cfg <- fromJSON(files$config, simplifyVector = FALSE)
  need(cfg$analysis$mroast_rotations == 9999L && cfg$analysis$fgsea_eps == 0 &&
    cfg$analysis$program_size_min == 10L && cfg$analysis$program_size_max == 500L &&
    cfg$analysis$random_seed == seed_base && cfg$gate$analysis_role == claim_role,
    "Original unsigned method options/claim ceiling")
  original <- readRDS(files$canonical_state)
  canonical <- import_canonical(files$canonical_source)
  check_space(original$labels, original$spaces, canonical)
  need(identical(dim(original$y$counts), c(14217L, 9L)) &&
    all(original$y$counts >= 0 & original$y$counts == round(original$y$counts)) &&
    isTRUE(fromJSON(files$canonical_completion)$scoped_numerical_acceptance), "Accepted retained count state")
  proof <- fread(files$canonical_filter_proof)
  need(nrow(proof) == 1680L && all(proof$same_mask) && all(proof$retained_gene_n == 14217L), "Accepted full filter proof")
  ann <- fread(annotation_path)[contrast == "local_disease", .(gene_id, symbol)]
  ann <- ann[match(rownames(original$y), gene_id)]
  need(identical(ann$gene_id, rownames(original$y)) && !anyNA(ann) && !anyDuplicated(ann$gene_id), "GTF-validated retained annotation alignment")
  resources <- map_resources(readRDS(resource_path), fread(files$orthology), ann$symbol)
  for (collection in collections) need(identical(resources[[collection]]$indices, original$indices[[collection]]),
    "Original row indices independently reconstructed from separate orthology")
  need(identical(lengths(lapply(resources, `[[`, "members")), c(Reactome = 1546L, SynGO_sensitivity = 78L)), "All original unsigned programs")
  runtime_path <- file.path(run, "RUNTIME_BEFORE.json"); jwrite(runtime_snapshot(), runtime_path)
  files$runtime <- runtime_path
  identities <- lapply(files, pin)
  manifest <- rbindlist(lapply(names(files), function(role) data.table(role = role,
    path = identities[[role]]$path, bytes = identities[[role]]$bytes, sha256 = identities[[role]]$sha256)))
  twrite(manifest, file.path(run, "INPUT_MANIFEST.tsv"))
  state <- list(y = original$y, labels = original$labels, spaces = original$spaces,
    symbols = ann$symbol, resources = resources,
    source_binding = pin(file.path(run, "INPUT_MANIFEST.tsv"))$sha256)
  expected_keys <- rbindlist(lapply(collections, function(collection) rbindlist(lapply(
    c("mroast_mean_9999_rotations", "fgseaMultilevel"), function(method) rbindlist(lapply(
      colnames(contrasts), function(name) data.table(collection = collection, method = method,
        contrast = name, feature = names(resources[[collection]]$members))))))))
  setorderv(expected_keys, key_columns); state$expected_keys <- as.data.frame(expected_keys)
  dir.create(file.path(run, "inputs"), showWarnings = FALSE)
  dir.create(file.path(run, "allocations"), showWarnings = FALSE)
  state_pin <- save_new_rds(state, file.path(run, "inputs", "state.rds"))
  observed <- stamp(fit_unsigned(state, state$labels[1L, ], canonical), state, 1L)
  observed_pin <- save_new_rds(observed, file.path(run, "allocations", "allocation_0001.rds"))
  # Complete RDS is retained before validating finite fields or observed compatibility.
  need(all(observed$canonical_specification$orientation == 1), "Original observed orientation positive in both contrasts")
  validate_result(observed, state, 1L)
  need(nrow(observed$result) == 6496L, "Eight complete unsigned families")
  bridge <- observed_bridge(observed, list(mroast = files$native_mroast, fgsea = files$native_fgsea))
  twrite(bridge, file.path(run, "OBSERVED_COMPATIBILITY.tsv"))
  need(all(bridge$accepted), "All observed native numeric/text/list fields must bridge before launch")
  twrite(observed$result, file.path(run, "OBSERVED.tsv"))
  for (method in c("mroast", "fgsea")) twrite(observed$native_payload[[method]], file.path(run,
    paste0("OBSERVED_DETAILS_", method, ".tsv")))
  capture.output(sessionInfo(), file = file.path(run, "sessionInfo.txt"))
  for (identity in identities) verify(identity)
  for (i in seq_len(nrow(upstream))) verify(as.list(upstream[i]))
  jwrite(list(status = "PREPARED_FULL_UNSIGNED_NATIVE_PAYLOAD_BRIDGE", completed_utc = utc(),
    derived = state_pin, observed = pin(file.path(run, "OBSERVED.tsv")), observed_checkpoint = observed_pin,
    worker = identities$worker, definition = identities$definition,
    input_manifest = pin(file.path(run, "INPUT_MANIFEST.tsv")),
    native_bridge = pin(file.path(run, "OBSERVED_COMPATIBILITY.tsv")),
    all_labels = 1680L, companion_labels = c(20L, 20L), family_n = 8L, result_rows = 6496L,
    retained_genes = 14217L, program_n = c(Reactome = 1546L, SynGO_sensitivity = 78L),
    observed_seconds = observed$elapsed_seconds, process_peak_memory_kb = peak_memory_kb(),
    packages = as.list(versions), RNG_kind = RNGkind(), runtime = identities$runtime,
    canonical_payload = TRUE, restored_native_Direction_or_leadingEdge = FALSE,
    inner_MC_uncertainty_eliminated_by_outer_enumeration = FALSE,
    authenticated_exchangeability = FALSE, candidate_gate_eligible = FALSE,
    historical_runtime_authenticated = FALSE, primary_results_modified = FALSE), file.path(run, "PREPARED.json"))
  cat("PREPARED_FULL_UNSIGNED_NATIVE_PAYLOAD_BRIDGE; seconds", observed$elapsed_seconds, "\n")
}
load_state <- function(run) {
  prepared <- fromJSON(file.path(run, "PREPARED.json"), simplifyVector = FALSE)
  for (identity in list(prepared$derived, prepared$observed, prepared$observed_checkpoint,
    prepared$worker, prepared$definition, prepared$input_manifest, prepared$native_bridge, prepared$runtime)) verify(identity)
  manifest <- fread(prepared$input_manifest$path)
  for (i in seq_len(nrow(manifest))) verify(as.list(manifest[i]))
  upstream <- fread(manifest$path[manifest$role == "canonical_input_manifest"])
  for (i in seq_len(nrow(upstream))) verify(as.list(upstream[i]))
  state <- readRDS(prepared$derived$path)
  need(state$source_binding == prepared$input_manifest$sha256, "Exact manifest-file source binding")
  need(runtime_snapshot()$kernel_contract_sha256 == fromJSON(prepared$runtime$path)$kernel_contract_sha256,
    "Current kernel/default/RNG/backend identity matches preparation")
  list(state = state, prepared = prepared)
}
run_worker <- function(run, worker_id, worker_n) {
  need(worker_n %in% 1:4 && worker_id %in% seq_len(worker_n), "Measured-resource concurrency, at most four single-thread workers")
  loaded <- load_state(run); state <- loaded$state
  canonical <- import_canonical(file.path(run, "canonical_source.R"))
  started_path <- file.path(run, paste0("WORKER_", worker_id, "_STARTED.json"))
  if (!file.exists(started_path)) jwrite(list(status = "RUNNING", started_utc = utc(), pid = Sys.getpid(),
    command = commandArgs(), runtime = runtime_snapshot(), worker = loaded$prepared$worker,
    derived = loaded$prepared$derived), started_path)
  for (id in state$spaces$global[(state$spaces$global - 1L) %% worker_n == worker_id - 1L]) {
    path <- file.path(run, "allocations", sprintf("allocation_%04d.rds", id))
    if (file.exists(path)) { validate_result(readRDS(path), state, id); next }
    value <- stamp(fit_unsigned(state, state$labels[id, ], canonical), state, id)
    save_new_rds(value, path)
    validate_result(value, state, id)
    cat("completed", id, "of1680; seconds", value$elapsed_seconds, "\n")
  }
  # Verify source and input custody again before a completed worker receipt.
  load_state(run)
  jwrite(list(status = "WORKER_COMPLETE", worker_id = worker_id, worker_n = worker_n,
    completed_utc = utc(), runtime = runtime_snapshot(), worker = loaded$prepared$worker,
    derived = loaded$prepared$derived, primary_results_modified = FALSE),
    file.path(run, paste0("WORKER_", worker_id, "_COMPLETE.json")))
}

args <- commandArgs(trailingOnly = TRUE)
need(length(args) >= 2L && args[1L] %in% c("static-selftest", "selftest", "prepare", "run"), "Explicit mode and run path")
run <- normalizePath(args[2L], mustWork = TRUE)
need(as.character(getRversion()) == "4.6.1", "Pinned current R runtime")
for (package in names(versions)) need(as.character(packageVersion(package)) == versions[package], "Pinned current package version")
need(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")), "Original RNG kind")
setDTthreads(1L)
need(all(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")) == "1") &&
  Sys.getenv("MKL_NUM_THREADS", unset = "1") == "1", "Explicit single-thread BLAS/OpenMP environment")
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
need(is.null(BiocParallel::bpRNGseed(BiocParallel::bpparam())) &&
  BiocParallel::bpnworkers(BiocParallel::bpparam()) == 1L, "Original SerialParam RNG-null backend")
set.seed(seed_base)
if (args[1L] %in% c("static-selftest", "selftest")) {
  need(length(args) == 2L, "Synthetic mode explicit run only")
  static_only <- args[1L] == "static-selftest"
  receipt <- selftest(run, static_only)
  executed_source <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1L])
  receipt$worker <- pin(executed_source)
  receipt$canonical_source <- pin(file.path(run, "canonical_source.R"))
  jwrite(receipt, file.path(run, if (static_only) "STATIC_TESTS.json" else "SYNTHETIC_TESTS.json"))
  cat(receipt$status, "\n")
} else if (args[1L] == "prepare") {
  need(length(args) == 4L, "Explicit canonical and native model roots")
  prepare(run, args[3L], args[4L])
} else {
  need(length(args) == 4L, "Explicit worker ID and measured-resource worker count (1..4)")
  run_worker(run, as.integer(args[3L]), as.integer(args[4L]))
}


