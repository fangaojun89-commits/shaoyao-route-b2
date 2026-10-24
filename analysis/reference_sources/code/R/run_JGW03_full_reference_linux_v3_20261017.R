#!/usr/bin/env Rscript
# Versioned non-equivalent compatibility contract, 17 October 2026.
# Technical v3 successor only: canonicalize verified file records after JSON
# readback. The v2 scientific contract, kernels, comparisons and tails are unchanged.
# Functions only when sourced. A separate Linux execution identity; never a
# Windows checkpoint continuation. Original kernels, compact/tail and MC code
# are source-pinned imports, not reimplemented here. No dispatch or installation.
# CLI:
# CONFIG prepare ACCESSION
# CONFIG bridge ACCESSION 1|2
# CONFIG run ACCESSION SPACE FIRST LAST ACCEPTANCE_JSON ACCEPTANCE_SHA256
# CONFIG finish ACCESSION SPACE
# windows_receipts per cell: prepared, observed_complete,
# nonidentity_complete, nonidentity_saved (file records; extra records are pinned).
# Acceptance is a separately produced, pinned JSON, not created by this entry.
# V1 source/strict failure remain immutable. Only PC12 observed may reference
# its explicitly pinned r01 saved native; it is never fitted again by v2.

jgwl_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwl_digest <- function(x) digest::digest(x, "sha256", serializeVersion = 3L)
jgwl_pin <- function(path) {
  jgwl_need(file.exists(path) && !dir.exists(path), "Required regular file absent")
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
jgwl_record <- function(record) {
  jgwl_need(is.list(record) && is.character(record$path) && length(record$path) == 1L &&
    !is.na(record$path) && is.numeric(record$bytes) && length(record$bytes) == 1L &&
    is.finite(record$bytes) && record$bytes >= 0 && record$bytes == floor(record$bytes) &&
    is.character(record$sha256) && length(record$sha256) == 1L &&
    !is.na(record$sha256) && grepl("^[0-9a-f]{64}$", record$sha256), "Malformed file record")
  actual <- jgwl_pin(record$path)
  jgwl_need(actual$bytes == record$bytes && identical(actual$sha256, record$sha256),
    "Pinned file identity differs")
  actual
}
jgwl_records <- function(records) {
  jgwl_need(is.list(records) && length(records) > 0L, "Explicit file records required")
  lapply(records, jgwl_record)
}
jgwl_inside <- function(path, root) identical(path, root) || startsWith(path, paste0(root, "/"))
jgwl_entry <- function() {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  jgwl_need(length(arg) == 1L, "Explicit entry source required")
  jgwl_pin(sub("^--file=", "", arg))
}
jgwl_json <- function(value, path) {
  jgwl_need(!file.exists(path), "Immutable JSON already exists")
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE,
    digits = NA, na = "null", null = "null")
  jgwl_pin(path)
}
# Comparison-only UTF-8 views: no source/output mutation or model-locale change.
jgwl2_utf8 <- function(x) {
  x <- as.character(x)
  known <- !is.na(x) & Encoding(x) %in% c("UTF-8", "latin1")
  out <- x
  out[known] <- enc2utf8(x[known])
  other <- !is.na(x) & !known
  out[other] <- iconv(x[other], from = "UTF-8", to = "UTF-8", sub = NA_character_)
  jgwl_need(identical(is.na(out), is.na(x)), "Invalid UTF-8 comparison text")
  out
}
jgwl2_probability <- function(name) {
  name <- tolower(name)
  name %in% c("p", "q", "pvalue", "p.value", "pval", "padj", "adj.p.val", "fdr",
    "pvalue.mixed", "fdr.mixed", "qvalue", "q_value") ||
    grepl("(^|_)(pvalue|p\\.value|pval|padj|adj\\.p\\.val|fdr|qvalue|q_value)(_|$)", name)
}
jgwl2_views <- function(fresh, old, key) {
  jgwl_need(is.data.frame(fresh) && is.data.frame(old), "Native tables required")
  exclude <- c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high")
  subset <- function(x) {
    keep <- !names(x) %in% exclude; y <- x[, keep, drop = FALSE]
    names(y) <- names(x)[keep]
    y
  }
  a <- subset(fresh); b <- subset(old)
  jgwl_need(identical(jgwl2_utf8(names(a)), jgwl2_utf8(names(b))) &&
    nrow(a) == nrow(b) && key %in% names(a) && key %in% names(b), "Native schema/key column differs")
  ka <- jgwl2_utf8(a[[key]]); kb <- jgwl2_utf8(b[[key]])
  jgwl_need(!anyNA(ka) && !anyNA(kb) && all(nzchar(ka)) && all(nzchar(kb)) &&
    !anyDuplicated(ka) && !anyDuplicated(kb) && setequal(ka, kb), "Stable key bijection required")
  index <- match(kb, ka)
  a <- a[index, , drop = FALSE]; names(a) <- names(b)
  same_types <- identical(vapply(a, typeof, character(1L)), vapply(b, typeof, character(1L))) &&
    identical(lapply(a, class), lapply(b, class))
  jgwl_need(same_types, "Native field types/classes differ")
  convert <- function(x) {
    for (i in seq_along(x)) {
      if (is.character(x[[i]])) x[[i]] <- jgwl2_utf8(x[[i]])
      if (is.factor(x[[i]])) levels(x[[i]]) <- jgwl2_utf8(levels(x[[i]]))
    }
    names(x) <- jgwl2_utf8(names(x)); x
  }
  list(fresh = convert(a), old = convert(b), key = key, types_exact = same_types,
    row_bijection_exact = TRUE, raw_row_order_identical = identical(ka, kb))
}
jgwl2_semantic_table <- function(fresh, old, key, endpoint) {
  v <- jgwl2_views(fresh, old, key)
  rows <- lapply(seq_along(v$old), function(i) {
    a <- v$fresh[[i]]; b <- v$old[[i]]; name <- names(v$old)[i]
    na <- identical(is.na(a), is.na(b)); nan <- identical(is.nan(a), is.nan(b))
    inf <- zero <- sign <- probability <- text <- TRUE
    if (is.numeric(a) && is.numeric(b)) {
      inf <- identical(is.infinite(a) & a > 0, is.infinite(b) & b > 0) &&
        identical(is.infinite(a) & a < 0, is.infinite(b) & b < 0)
      zero <- identical(!is.na(a) & a == 0, !is.na(b) & b == 0)
      finite <- is.finite(a) & is.finite(b)
      sign <- identical(is.finite(a), is.finite(b)) && all(base::sign(a[finite]) == base::sign(b[finite]))
      if (jgwl2_probability(name)) probability <-
        all(is.na(a) | (is.finite(a) & a >= 0 & a <= 1)) &&
        all(is.na(b) | (is.finite(b) & b >= 0 & b <= 1)) &&
        identical(a < .05, b < .05) && identical(a <= .05, b <= .05)
    } else {
      text <- identical(as.character(a), as.character(b))
      if (is.factor(a)) text <- text && identical(levels(a), levels(b))
    }
    data.frame(endpoint = endpoint, column = paste0(name, "[", i, "]"), n = length(b),
      schema_types_keys = TRUE, NA_exact = na, NaN_exact = nan, signed_Inf_exact = inf,
      zero_exact = zero, finite_sign_exact = sign, probability_both_005_exact = probability,
      text_direction_exact = text, pass = na && nan && inf && zero && sign && probability && text,
      check.names = FALSE)
  })
  list(checks = do.call(rbind, rows), views = v)
}
jgwl2_diagnostics_valid <- function(x) {
  is.list(x) && length(x) == 2L && all(vapply(x, function(one) {
    flag <- one$betaConv
    fit <- one$fitType
    is.character(fit) && length(fit) == 1L && !is.na(fit) && nzchar(fit) &&
      is.list(flag) && identical(flag$status, "SAVED_NATIVE_FLAG_COUNTS") &&
      all(c("n", "false_n", "unavailable_n") %in% names(flag)) &&
      all(vapply(flag[c("n", "false_n", "unavailable_n")], function(n) is.numeric(n) &&
        length(n) == 1L && !is.na(n) && is.finite(n) && n >= 0 && n == floor(n), logical(1L))) &&
      flag$n > 0 && flag$false_n + flag$unavailable_n <= flag$n
  }, logical(1L)))
}
jgwl2_origin_records <- function(x) {
  required <- c("config", "prepared_receipt", "prepared", "fixed", "fit_saved",
    "bridge_complete", "process_exit", "process_command", "process_started", "native")
  root <- "/root/autodl-tmp/JGW03_full_reference_linux_20261017r01"
  jgwl_need(is.list(x) && identical(x$origin_run_root, root) &&
    all(required %in% names(x)), "Explicit original PC12 reuse records required")
  records <- jgwl_records(x[setdiff(names(x), "origin_run_root")])
  jgwl_need(all(vapply(records, function(r) jgwl_inside(r$path, root), logical(1L))),
    "Origin record escapes exact preserved r01 root")
  list(origin_run_root = root, records = records)
}
jgwl2_origin_links <- function(origin, metadata) {
  r <- origin$records
  eq <- function(a, b) isTRUE(a$bytes == b$bytes) && identical(a$sha256, b$sha256)
  fit <- metadata$fit_saved; bridge <- metadata$bridge_complete; process <- metadata$process_exit
  command <- metadata$process_command; started <- metadata$process_started
  integer_scalar <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x > 0 && x == floor(x)
  argv <- unlist(command$argv, use.names = FALSE)
  # The preserved EXIT has no entry/config hashes. Those are bound through the
  # original COMMAND and STARTED records, never retroactively added to EXIT.
  jgwl_need(eq(r$prepared, metadata$prepared_receipt$prepared) &&
    eq(r$native, fit$native) && eq(r$native, bridge$native) &&
    identical(fit$config_sha256, r$config$sha256) &&
    identical(fit$prepared_sha256, r$prepared$sha256) &&
    identical(bridge$prepared_sha256, r$prepared$sha256) &&
    identical(bridge$config_sha256, r$config$sha256) &&
    identical(fit$entry_sha256, "57e3643a06b5c17f66c46eb658d150665c9bc55a1b4cc025ceb3e21d6d2cdc8f") &&
    identical(bridge$entry_sha256, fit$entry_sha256) &&
    identical(bridge$status, "LINUX_NATIVE_BRIDGE_REVIEW_REQUIRED") &&
    identical(bridge$accession, "GSE289929") && bridge$allocation_id == 1 &&
    identical(fit$allocation_key, bridge$allocation_key) &&
    eq(command$config, r$config) && command$entry$bytes == 32584 &&
    identical(command$entry$sha256, fit$entry_sha256) &&
    identical(command$mode, "bridge") && identical(command$accession, "GSE289929") &&
    command$allocation == 1 && is.character(argv) && length(argv) == 7L &&
    identical(argv[2:7], c("--vanilla",
      paste0(origin$origin_run_root, "/project/code/R/run_JGW03_full_reference_linux_20261017.R"),
      r$config$path, "bridge", "GSE289929", "1")) &&
    identical(started$config_sha256, r$config$sha256) &&
    integer_scalar(started$actual_R_pid) && integer_scalar(started$start_ticks) &&
    integer_scalar(process$actual_R_pid) &&
    started$actual_R_pid == process$actual_R_pid && process$actual_R_exit_code == 2,
    "Preserved first Linux native/strict failure/actual exit-2 provenance differs")
  invisible(TRUE)
}
jgwl2_reuse_PC12 <- function(z, prepared, allocation) {
  origin <- z$before$original_PC12; r <- origin$records
  metadata <- lapply(r[c("config", "prepared_receipt", "fit_saved", "bridge_complete", "process_exit",
    "process_command", "process_started")],
    function(p) jsonlite::read_json(p$path, simplifyVector = FALSE))
  jgwl2_origin_links(origin, metadata)
  original_prepared <- readRDS(r$prepared$path)
  jgwl_need(identical(original_prepared$fixed$sha256, r$fixed$sha256) &&
    identical(original_prepared$binding$config$sha256, r$config$sha256) &&
    z$h$jgwp_json_equal(original_prepared$runtime, z$before$runtime) &&
    identical(original_prepared$registry$sha256, z$ref$registry$sha256),
    "Original first-fit fixed/registry/runtime binding differs")
  for (name in names(metadata$config$source_pins)) {
    old <- metadata$config$source_pins[[name]]; now <- z$cfg$source_pins[[name]]
    jgwl_need(is.list(now) && old$bytes == now$bytes && identical(old$sha256, now$sha256),
      "Original scientific/generation source changed")
  }
  value <- readRDS(r$native$path)
  jgwl_need(identical(value$accession, "GSE289929") && identical(value$labels, allocation),
    "Original Linux observed labels differ")
  session <- jgwl_record(metadata$fit_saved$session_info)
  execution <- metadata$fit_saved$execution
  execution$RNGkind <- unlist(execution$RNGkind, use.names = FALSE)
  list(value = value, saved = list(status = "LINUX_NATIVE_REUSED_SAVED_COMPATIBILITY_PENDING_V2",
    kind = "REUSED_ORIGINAL_LINUX_NATIVE", native = r$native, session_info = session,
    allocation_key = paste(as.integer(allocation), collapse = ""),
    prepared_sha256 = prepared$pin$sha256, config_sha256 = z$before$config$sha256,
    entry_sha256 = z$before$entry$sha256, execution = execution,
    fit_seconds = metadata$fit_saved$fit_seconds, new_native_calls = 0L,
    origin = list(records = r, actual_original_R_exit_code = 2L,
      actual_original_R_pid = metadata$process_exit$actual_R_pid,
      actual_original_R_start_ticks = metadata$process_started$start_ticks,
      source_generation_entry_sha256 = metadata$fit_saved$entry_sha256,
      new_linkage_entry = z$before$entry, original_strict_status = metadata$bridge_complete$status)))
}
jgwl_runtime <- function(h, expected) {
  actual <- h$jgwp_runtime()
  fields <- c("R", "packages", "platform", "BLAS", "backend", "workers", "seed", "RNGkind")
  jgwl_need(is.list(expected) && all(fields %in% names(expected)) &&
    identical(actual$R, "4.6.1") && grepl("linux", actual$platform, ignore.case = TRUE) &&
    setequal(names(expected$packages), names(actual$packages)) &&
    h$jgwp_json_equal(actual[fields], expected[fields]), "Explicit Linux runtime differs")
  threads <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"), unset = "")
  jgwl_need(all(threads == "1") &&
    inherits(BiocParallel::bpparam(), "SerialParam") &&
    BiocParallel::bpnworkers(BiocParallel::bpparam()) == 1L, "Actual serial backend/threads required")
  actual$library_paths <- normalizePath(.libPaths(), winslash = "/", mustWork = TRUE)
  if (!is.null(expected$library_paths)) jgwl_need(
    identical(actual$library_paths, unlist(expected$library_paths, use.names = FALSE)),
    "R library precedence differs")
  actual$thread_environment <- as.list(threads)
  actual$locale <- Sys.getlocale()
  actual
}
jgwl_open <- function(config_path, accession) {
  config_pin <- jgwl_pin(config_path)
  cfg <- jsonlite::read_json(config_pin$path, simplifyVector = FALSE)
  project <- normalizePath(cfg$project_root, winslash = "/", mustWork = TRUE)
  run <- normalizePath(cfg$run_root, winslash = "/", mustWork = TRUE)
  jgwl_need(identical(cfg$schema, "jgw03-separate-linux-reference-config-v2") &&
    identical(run, "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03") &&
    identical(project, paste0(run, "/project")) && identical(config_pin$path, paste0(run, "/CONFIG.json")) &&
    accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    setequal(names(cfg$cells), c("GSE289929", "GSE227494", "GSE299063")),
    "Exact separate Linux root/config/cell contract required")
  relative <- names(cfg$source_pins)
  jgwl_need(length(relative) > 0L && !anyDuplicated(relative) &&
    all(nzchar(relative)) && all(!grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))", relative)),
    "Project-relative source manifest required")
  sources <- lapply(seq_along(relative), function(i) {
    rec <- cfg$source_pins[[i]]; rec$path <- file.path(project, relative[i])
    pin <- jgwl_record(rec)
    jgwl_need(jgwl_inside(pin$path, project), "Source escapes staged project")
    pin
  }); names(sources) <- relative
  required <- c("code/R/prepare_JGW03_reference_first_nonidentity_20261015.R",
    "code/R/run_JGW03_full_reference_single_20261015.R",
    "code/R/prepare_JGW03_observed_label_pilot_v3_20261015.R",
    "code/R/JGW03_label_reference_core_20261015.R", "code/R/JGW03_input_admission_v2_20261015.R",
    "code/R/JGW03_reference_field_map_20261015.R", "code/R/JGW03_reference_compact_tail_20261015.R",
    "code/R/JGW03_finite_population_mc_v2_20261015.R", "code/R/analyze_jgw03_three_cell_lines.R",
    "code/R/run_JGW03_full_reference_linux_20261017.R",
    "code/R/run_JGW03_full_reference_linux_v3_20261017.R")
  jgwl_need(all(required %in% relative), "Required immutable source absent from manifest")
  jgwl_need(is.character(cfg$checker_source) && length(cfg$checker_source) == 1L &&
    cfg$checker_source %in% relative, "Independent checker must be source-pinned")
  jgwl_need(identical(sources[["code/R/run_JGW03_full_reference_linux_20261017.R"]]$sha256,
    "57e3643a06b5c17f66c46eb658d150665c9bc55a1b4cc025ceb3e21d6d2cdc8f"), "Preserved v1 entry differs")
  definition <- jgwl_record(cfg$compatibility_definition); approval <- jgwl_record(cfg$approval)
  jgwl_need(identical(definition$sha256,
    "5c06cb8e054c4b2e3a261c2250944c09939d6b6851a65043222ac2360de3aff8"), "Approved compatibility definition differs")
  origin <- jgwl2_origin_records(cfg$reused_first_PC12)
  entry <- jgwl_entry()
  jgwl_need(identical(entry, sources[[tail(required, 1L)]]), "Executing entry differs from manifest")
  # Both old entries are source-safe. Their main/bootstrap is NOT invoked.
  e <- new.env(parent = asNamespace("stats"))
  sys.source(sources[[required[1L]]]$path, envir = e)
  modules <- e$jgwx_modules(project)
  jgwl_need(identical(sources[[required[2L]]]$sha256,
    "f1de9c41003e0a96c6e08e3367a3e34f3a9ce898cce7a0be694bf9d4965ff4e5"),
    "Original full worker differs")
  sys.source(sources[[required[2L]]]$path, envir = e)
  runtime <- jgwl_runtime(modules$helper, cfg$expected_runtime)
  inputs <- jgwl_records(cfg$input_pins)
  jgwl_need(all(vapply(inputs, function(x) jgwl_inside(x$path, project), logical(1))),
    "Scientific inputs must belong to staged project")
  cell <- cfg$cells[[accession]]
  ref <- lapply(cell[c("registry", "windows_observed", "windows_nonidentity", "windows_fixed")], jgwl_record)
  jgwl_need(length(ref) == 4L && all(vapply(ref, function(x) jgwl_inside(x$path, run), logical(1))),
    "Explicit reference copies must belong to separate run")
  receipts <- jgwl_records(cell$windows_receipts)
  jgwl_need(all(c("prepared", "observed_complete", "nonidentity_complete", "nonidentity_saved") %in%
    names(receipts)) && all(vapply(receipts, function(x) jgwl_inside(x$path, run), logical(1))),
    "Windows provenance receipts required")
  metadata <- lapply(receipts, function(x) jsonlite::read_json(x$path, simplifyVector = FALSE))
  equal_pin <- function(left, right) left$bytes == right$bytes && identical(left$sha256, right$sha256)
  wp <- metadata$prepared; wo <- metadata$observed_complete
  wn <- metadata$nonidentity_complete; ws <- metadata$nonidentity_saved
  jgwl_need(identical(wo$accession, accession) && identical(wn$accession, accession) &&
    identical(wo$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_ONLY") && isTRUE(wo$bridge$all_fields_pass) &&
    identical(wn$status, "FIRST_NONIDENTITY_SAVED_AND_COMPACTED_ONLY") && wn$allocation_id == 2 &&
    equal_pin(ref$registry, wp$label_registry) && equal_pin(ref$windows_observed, wo$result) &&
    equal_pin(ref$windows_fixed, wp$observed_files$fixed) && equal_pin(ref$windows_nonidentity, wn$native) &&
    equal_pin(ref$windows_nonidentity, ws$native) &&
    identical(wn$prepared_sha256, wp$prepared$sha256) &&
    identical(ws$prepared_sha256, wp$prepared$sha256) && identical(ws$allocation_key, wn$allocation_key),
    "Windows reference-to-receipt association differs")
  # A serialized registry is imported byte-for-byte, never regenerated here.
  registry <- readRDS(ref$registry$path)
  keys <- e$jgwx_labels(registry)
  jgwl_need(identical(keys$four_group_MC[2L], wn$allocation_key), "Fixed Windows ID2 key differs")
  expected_inputs <- modules$admission$jgw03_v2_expected_files(accession)
  ledger_path <- file.path(project, "metadata/samples.tsv")
  mandatory <- c(file.path(project, expected_inputs$relative_path), ledger_path)
  jgwl_need(all(mandatory %in% vapply(inputs, `[[`, character(1L), "path")),
    "Eight frozen scientific inputs and independent sample ledger must be manifest-bound")
  before <- list(config = config_pin, sources = sources, inputs = inputs, references = ref,
    windows_receipts = receipts, entry = entry, runtime = runtime,
    compatibility_definition = definition, approval = approval, original_PC12 = origin)
  list(cfg = cfg, project = project, run = run, accession = accession, e = e, modules = modules,
    h = modules$helper, ref = ref, metadata = metadata, registry = registry, keys = keys,
    cell_dir = file.path(run, "intermediate", accession), before = before)
}
jgwl_after <- function(z) {
  records <- c(list(z$before$config, z$before$entry), z$before$sources, z$before$inputs,
    z$before$references, z$before$windows_receipts,
    list(z$before$compatibility_definition, z$before$approval), z$before$original_PC12$records)
  invisible(lapply(records, jgwl_record))
  jgwl_need(z$h$jgwp_json_equal(z$before$runtime,
    jgwl_runtime(z$h, z$cfg$expected_runtime)), "Linux runtime changed")
  invisible(TRUE)
}
jgwl_fixed_check <- function(fresh, old, admission) list(
  accession = identical(fresh$accession, old$accession), keep = identical(fresh$keep, old$keep),
  retained_gene_sample_order = identical(dimnames(fresh$y$counts), dimnames(old$y$counts)),
  retained_counts = identical(fresh$y$counts, old$y$counts),
  metadata = admission$jgw03_v2_table_equal(fresh$original_samples, old$original_samples, "log10_dose_plus_1"),
  metadata_row_names = identical(rownames(fresh$original_samples), rownames(old$original_samples)),
  annotation = admission$jgw03_v2_table_equal(fresh$annotation, old$annotation),
  library_sizes = isTRUE(all.equal(fresh$y$samples$lib.size, old$y$samples$lib.size,
    tolerance = 1e-12, check.attributes = FALSE)),
  TMM_factors = isTRUE(all.equal(fresh$y$samples$norm.factors, old$y$samples$norm.factors,
    tolerance = 1e-12, check.attributes = FALSE)))
jgwl_prepare <- function(z) {
  jgwl_need(!file.exists(z$cell_dir), "Preparation leaf already exists; retain without overwrite")
  jgwl_need(dir.create(z$cell_dir, recursive = TRUE), "Cannot exclusively create preparation leaf")
  jgwl_json(list(schema = "jgw03-linux-preparation-started-v2", binding = z$before),
    file.path(z$cell_dir, "PREPARATION_STARTED.json"))
  bound <- z$modules$admission$jgw03_v2_read_bound_inputs(z$project, z$accession)
  ledger <- z$h$jgwp_read(file.path(z$project, "metadata/samples.tsv"))
  ledger <- ledger[ledger$dataset_id == z$accession & ledger$include_primary == "TRUE", , drop = FALSE]
  ledger <- ledger[order(as.integer(ledger$sample_order)), , drop = FALSE]
  jgwl_need(identical(as.character(bound$samples$sample_id), as.character(ledger$sample_id)),
    "Independent sample ledger order differs")
  for (field in intersect(names(ledger), names(bound$samples))) {
    a <- as.character(ledger[[field]]); b <- as.character(bound$samples[[field]])
    missing_a <- is.na(a) | a == ""; missing_b <- is.na(b) | b == ""
    jgwl_need(identical(missing_a, missing_b) && identical(a[!missing_a], b[!missing_b]),
      "Independent sample ledger literal differs")
  }
  rownames(bound$samples) <- rownames(ledger)
  rownames(bound$annotation) <- bound$annotation$raw_gene_id
  observed <- setNames(as.numeric(bound$samples$ordered_dose), as.character(bound$samples$sample_id))
  jgwl_need(identical(colnames(z$registry$MC), names(observed)) &&
    identical(unname(as.integer(z$registry$MC[1L, ])), unname(as.integer(observed))),
    "Frozen observed registry/sample identity differs")
  do.call(RNGkind, as.list(z$before$runtime$RNGkind)); set.seed(z$before$runtime$seed)
  # This is the original label-independent fixed/scoring layer; no allocation
  # native kernel is called in prepare. No Windows fixed object is used as input.
  start <- proc.time()[["elapsed"]]
  fixed <- z$modules$kernel$jgw03_build_fixed_layer(bound$counts, bound$annotation, bound$samples,
    z$accession, bound$resources, z$modules$pins$native$path, z$project)
  fixed_seconds <- proc.time()[["elapsed"]] - start
  fixed_pin <- c(list(path = file.path(z$cell_dir, "FIXED_LAYER.rds")),
    z$e$jgwx_atomic_plain(fixed, file.path(z$cell_dir, "FIXED_LAYER.rds")))
  old <- readRDS(z$ref$windows_fixed$path)
  comparison <- jgwl_fixed_check(fixed, old, z$modules$admission)
  rm(old)
  if (z$accession == "GSE289929") {
    origin_fixed <- readRDS(z$before$original_PC12$records$fixed$path)
    # Only the path-dependent input-admission provenance differs between roots.
    # Every original scientific fixed-layer value/type/attribute must be exact.
    science <- setdiff(names(fixed), "input_admission_v2")
    comparison$original_Linux_fixed_exact <- identical(science,
      setdiff(names(origin_fixed), "input_admission_v2")) && identical(fixed[science], origin_fixed[science])
    rm(origin_fixed)
  }
  invisible(gc())
  check_pin <- jgwl_json(comparison, file.path(z$cell_dir, "FIXED_LAYER_BRIDGE.json"))
  jgwl_need(all(unlist(comparison)), "Fixed layer bridge requires review; no allocation fit performed")
  registry_path <- file.path(z$cell_dir, "FULL_LABEL_REGISTRY.rds")
  jgwl_need(!file.exists(registry_path) && file.copy(z$ref$registry$path, registry_path, overwrite = FALSE),
    "Registry byte-copy failed")
  registry_pin <- jgwl_record(c(list(path = registry_path), z$ref$registry[c("bytes", "sha256")]))
  session_path <- file.path(z$cell_dir, "sessionInfo.txt")
  utils::capture.output(utils::sessionInfo(), file = session_path)
  session_pin <- jgwl_pin(session_path)
  context <- list(schema = "jgw03-separate-linux-preparation-v2", accession = z$accession,
    root = z$project, run_root = z$run, binding = z$before, runtime = z$before$runtime,
    fixed = fixed_pin, fixed_bridge = check_pin, registry = registry_pin, session_info = session_pin,
    observed_labels = observed, sources = z$modules$pins, entry = z$before$entry,
    fixed_seconds = fixed_seconds, full_MC_nonidentity_draws = 10000L,
    full_companion_including_observed = 20L, full_MC_family_n = 89L, full_companion_family_n = 34L,
    native_allocation_calls = 0L, independent_bridge_and_resource_acceptance_required = TRUE)
  jgwl_after(z)
  prepared <- c(list(path = file.path(z$cell_dir, "PREPARED.rds")),
    z$e$jgwx_atomic_plain(context, file.path(z$cell_dir, "PREPARED.rds")))
  jgwl_json(list(status = "LINUX_FIXED_AND_ORIGINAL_REGISTRY_PREPARED_NO_ALLOCATION_FIT",
    accession = z$accession, prepared = prepared, fixed = fixed_pin, registry = registry_pin,
    session_info = session_pin,
    config_sha256 = z$before$config$sha256, entry_sha256 = z$before$entry$sha256),
    file.path(z$cell_dir, "PREPARED.json"))
  0L
}
jgwl_prepared <- function(z) {
  meta <- jsonlite::read_json(file.path(z$cell_dir, "PREPARED.json"), simplifyVector = FALSE)
  pin <- jgwl_record(meta$prepared); context <- readRDS(pin$path)
  jgwl_need(identical(context$schema, "jgw03-separate-linux-preparation-v2") &&
    identical(context$accession, z$accession) && identical(context$root, z$project) &&
    identical(context$run_root, z$run) && z$h$jgwp_json_equal(context$binding, z$before),
    "Prepared Linux input/source/runtime/config identity differs")
  jgwl_record(context$fixed); jgwl_record(context$fixed_bridge); jgwl_record(context$registry)
  jgwl_record(context$session_info)  # Provenance pin, not a cross-stage text equality test.
  fixed_bridge <- jsonlite::read_json(context$fixed_bridge$path, simplifyVector = TRUE)
  jgwl_need(all(unlist(fixed_bridge)), "Fixed bridge not accepted")
  list(context = context, pin = pin)
}
jgwl_bridge_record <- function(z, id, prepared) {
  path <- file.path(z$cell_dir, paste0("bridge_", id), "BRIDGE_COMPLETE.json")
  record <- jgwl_pin(path); meta <- jsonlite::read_json(path, simplifyVector = FALSE)
  key <- z$keys$four_group_MC[id]
  jgwl_need(identical(meta$status, "LINUX_NATIVE_COMPATIBILITY_PASS_V2") && meta$allocation_id == id &&
    identical(meta$allocation_key, key) && identical(meta$prepared_sha256, prepared$pin$sha256) &&
    identical(meta$config_sha256, z$before$config$sha256) &&
    identical(meta$entry_sha256, z$before$entry$sha256), "Accepted same-Linux bridge required")
  for (name in c("native", "compact", "checks", "semantic_checks", "aligned_strict_checks", "checks_saved", "session_info")) jgwl_record(meta[[name]])
  if (id == 1L) jgwl_record(meta$domain)
  if (z$accession == "GSE289929" && id == 1L) jgwl_need(
    identical(meta$native_generation_kind, "REUSED_ORIGINAL_LINUX_NATIVE") && meta$new_native_calls == 0 &&
    identical(jgwl_record(meta$native), jgwl_record(z$before$original_PC12$records$native)),
    "Observed PC12 lineage must remain the original exit-2 native")
  list(record = record, meta = meta)
}
jgwl_bridge <- function(z, id) {
  jgwl_need(id %in% c(1L, 2L), "Only fixed observed or ID2 bridge allowed")
  prepared <- jgwl_prepared(z); context <- prepared$context
  if (id == 2L) first <- jgwl_bridge_record(z, 1L, prepared)
  leaf <- file.path(z$cell_dir, paste0("bridge_", id))
  if (file.exists(file.path(leaf, "BRIDGE_COMPLETE.json"))) {
    jgwl_bridge_record(z, id, prepared); jgwl_after(z); return(0L)
  }
  if (!dir.exists(leaf)) jgwl_need(dir.create(leaf), "Cannot create bridge leaf")
  fixed <- readRDS(context$fixed$path)
  allocation <- setNames(as.numeric(z$registry$MC[id, ]), colnames(z$registry$MC))
  z$modules$kernel$jgw03_check_allocation(allocation, names(context$observed_labels))
  saved_path <- file.path(leaf, "FIT_SAVED.json"); native_path <- file.path(leaf, "NATIVE_CAPSULE.rds")
  if (file.exists(saved_path)) {
    saved <- jsonlite::read_json(saved_path, simplifyVector = FALSE)
    jgwl_need(identical(saved$prepared_sha256, prepared$pin$sha256) &&
      identical(saved$config_sha256, z$before$config$sha256) &&
      identical(saved$entry_sha256, z$before$entry$sha256) &&
      identical(saved$allocation_key, z$keys$four_group_MC[id]), "Retained Linux fit binding differs")
    jgwl_record(saved$native); jgwl_record(saved$session_info); value <- readRDS(saved$native$path)
    execution <- saved$execution
    execution$RNGkind <- unlist(execution$RNGkind, use.names = FALSE)
  } else if (z$accession == "GSE289929" && id == 1L) {
    jgwl_need(!file.exists(native_path) && !file.exists(paste0(native_path, ".partial")) &&
      !file.exists(file.path(leaf, "FIT_STARTED.json")), "Observed PC12 is reuse-only; no fresh fit allowed")
    reused <- jgwl2_reuse_PC12(z, prepared, allocation)
    value <- reused$value; saved <- reused$saved; execution <- saved$execution
    jgwl_json(saved, saved_path)
  } else {
    jgwl_need(!file.exists(native_path) && !file.exists(paste0(native_path, ".partial")) &&
      !file.exists(file.path(leaf, "FIT_STARTED.json")), "Incomplete bridge retained; no automatic refit")
    jgwl_json(list(status = "LINUX_BRIDGE_NATIVE_CALL_STARTED", prepared_sha256 = prepared$pin$sha256,
      allocation_key = z$keys$four_group_MC[id], entry = z$before$entry), file.path(leaf, "FIT_STARTED.json"))
    do.call(RNGkind, as.list(context$runtime$RNGkind)); set.seed(context$runtime$seed)
    started <- proc.time()[["elapsed"]]
    value <- z$modules$kernel$jgw03_fit_one_allocation(fixed, allocation, z$modules$pins$native$path)
    elapsed <- proc.time()[["elapsed"]] - started
    saveRDS(value, paste0(native_path, ".partial"), compress = "gzip", version = 3L)
    jgwl_need(file.rename(paste0(native_path, ".partial"), native_path), "Native capsule rename failed")
    session_path <- file.path(leaf, "sessionInfo.txt")
    jgwl_need(!file.exists(session_path), "Existing fit-session record retained")
    utils::capture.output(utils::sessionInfo(), file = session_path)
    execution <- list(seed = context$runtime$seed, RNGkind = context$runtime$RNGkind,
      backend = "SerialParam", workers = 1L, mroast_nrot = 1999L, fgsea_eps = 0,
      mroast_applicable = z$accession == "GSE289929", mroast_executed = z$accession == "GSE289929",
      actual_mroast_nrot = if (z$accession == "GSE289929") 1999L else NULL,
      provenance = paste("Separate Linux native bridge; prepared", prepared$pin$sha256),
      RNG_after_sha256 = jgwl_digest(.Random.seed), historical_backend_authenticated = FALSE)
    saved <- list(status = "LINUX_NATIVE_FIT_SAVED_BRIDGE_PENDING", native = jgwl_pin(native_path),
      allocation_key = z$keys$four_group_MC[id], prepared_sha256 = prepared$pin$sha256,
      config_sha256 = z$before$config$sha256, entry_sha256 = z$before$entry$sha256,
      fit_seconds = elapsed, execution = execution, session_info = jgwl_pin(session_path),
      kind = "NEW_LINUX_NATIVE_V2", new_native_calls = 1L,
      source_generation_entry_sha256 = z$before$entry$sha256)
    jgwl_json(saved, saved_path)
  }
  reused_PC12 <- z$accession == "GSE289929" && id == 1L
  if (reused_PC12) {
    jgwl_need(identical(saved$kind, "REUSED_ORIGINAL_LINUX_NATIVE") && saved$new_native_calls == 0 &&
      identical(jgwl_record(saved$native), jgwl_record(z$before$original_PC12$records$native)),
      "PC12 observed cannot be generated or relabelled as a new v2 fit")
  } else jgwl_need(identical(saved$kind, "NEW_LINUX_NATIVE_V2") && saved$new_native_calls == 1 &&
    identical(saved$native$path, native_path), "Other bridge must be its own new Linux native")
  jgwl_need(identical(value$accession, z$accession) && identical(value$labels, allocation),
    "Saved native accession/label identity differs")
  map <- z$h$jgw03_reference_field_map(z$modules$kernel$jgw03_family_registry())
  domain_path <- file.path(leaf, "COMPACT_DOMAIN.rds")
  domain <- if (id == 1L) z$h$jgw03c_make_domain(value, map) else readRDS(first$meta$domain$path)
  if (id == 1L) {
    if (file.exists(domain_path)) jgwl_need(identical(readRDS(domain_path), domain), "Retained Linux domain differs") else
      z$e$jgwx_atomic_plain(domain, domain_path)
    for (space in names(z$keys)) z$h$jgw03c_tail_start(
      z$h$jgw03c_extract(value, domain, execution), domain, space, z$keys[[space]])
  }
  compact <- z$h$jgw03c_extract(value, domain, execution)
  loaded <- readRDS(saved$native$path)
  transported <- z$h$jgw03c_extract(loaded, domain, execution)
  jgwl_need(identical(serialize(compact, NULL, version = 3L), serialize(transported, NULL, version = 3L)),
    "Linux native-to-compact transport differs")
  compact_path <- file.path(leaf, "COMPACT.rds")
  if (file.exists(compact_path)) jgwl_need(identical(readRDS(compact_path), compact), "Retained compact differs") else
    z$e$jgwx_atomic_plain(compact, compact_path)
  reference <- if (id == 1L) z$ref$windows_observed else z$ref$windows_nonidentity
  old <- readRDS(reference$path); old_fixed <- readRDS(z$ref$windows_fixed$path)
  jgwl_need(identical(old$accession, z$accession) && identical(old$labels, allocation),
    "Windows comparison capsule label differs")
  fresh <- z$h$jgwp_fresh_tables(fixed, value); prior <- z$h$jgwp_fresh_tables(old_fixed, old)
  jgwl_need(setequal(names(fresh), names(prior)) && !anyDuplicated(names(fresh)) &&
    !anyDuplicated(names(prior)), "Full endpoint table coverage differs")
  strict_original <- aligned_strict <- semantic <- list()
  for (name in names(prior)) {
    old_table <- prior[[name]]; table <- fresh[[name]]
    key <- if ("raw_gene_id" %in% names(old_table)) "raw_gene_id" else
      if ("term_id" %in% names(old_table)) "term_id" else "feature"
    # Keep the original raw-order Windows strict result, including its failures.
    strict_original[[name]] <- z$h$jgwp_table_check(table, old_table, key, name)
    compared <- jgwl2_semantic_table(table, old_table, key, name)
    semantic[[name]] <- compared$checks
    aligned_strict[[name]] <- z$h$jgwp_table_check(compared$views$fresh, compared$views$old, key, name)
  }
  reports <- list(strict_original = do.call(rbind, strict_original),
    strict_aligned_UTF8 = do.call(rbind, aligned_strict), semantic = do.call(rbind, semantic))
  report_paths <- file.path(leaf, c("NATIVE_FIELD_BRIDGE.tsv",
    "ALIGNED_UTF8_STRICT_FIELD_DIAGNOSTICS.tsv", "COMPATIBILITY_FIELDS.tsv"))
  names(report_paths) <- names(reports)
  checks_saved_path <- file.path(leaf, "CHECKS_SAVED.json")
  if (any(file.exists(report_paths))) {
    jgwl_need(all(file.exists(report_paths)) && file.exists(checks_saved_path),
      "Diagnostic write incomplete; preserve native without refit")
    prior_checks <- jsonlite::read_json(checks_saved_path, simplifyVector = FALSE)
    report_pins <- jgwl_records(prior_checks$files)
    jgwl_need(identical(vapply(report_pins, `[[`, character(1L), "path"), report_paths) &&
      identical(prior_checks$full_precision_checks_sha256, jgwl_digest(reports)) &&
      identical(prior_checks$native_sha256, saved$native$sha256) &&
      identical(prior_checks$prepared_sha256, prepared$pin$sha256), "Retained diagnostic coupling differs")
  } else {
    jgwl_need(!file.exists(checks_saved_path), "Orphaned diagnostic receipt retained")
    for (name in names(reports)) utils::write.table(reports[[name]], report_paths[[name]],
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    report_pins <- lapply(report_paths, jgwl_pin)
    jgwl_json(list(files = report_pins, full_precision_checks_sha256 = jgwl_digest(reports),
      native_sha256 = saved$native$sha256, prepared_sha256 = prepared$pin$sha256), checks_saved_path)
  }
  diagnostics <- z$e$jgwf_native_diagnostics(value)
  old_diagnostics <- z$e$jgwf_native_diagnostics(old)
  diagnostic_equal <- z$h$jgwp_json_equal(diagnostics, old_diagnostics)
  roast <- z$h$jgwp_roast_check(value)
  native_diagnostics_valid <- jgwl2_diagnostics_valid(diagnostics)
  passed <- all(reports$semantic$pass) && native_diagnostics_valid &&
    (!isTRUE(roast$applicable) || isTRUE(roast$pass))
  bridge <- list(status = if (passed) "LINUX_NATIVE_COMPATIBILITY_PASS_V2" else "LINUX_NATIVE_COMPATIBILITY_REVIEW_REQUIRED_V2",
    accession = z$accession, allocation_id = id, allocation_key = compact$allocation_key,
    prepared_sha256 = prepared$pin$sha256, config_sha256 = z$before$config$sha256,
    entry_sha256 = z$before$entry$sha256, native = saved$native, compact = jgwl_pin(compact_path),
    checks = report_pins$strict_original, semantic_checks = report_pins$semantic,
    aligned_strict_checks = report_pins$strict_aligned_UTF8,
    checks_saved = jgwl_pin(checks_saved_path), native_family_n = length(compact$families),
    session_info = saved$session_info,
    Windows_reference_only_never_used_as_Linux_result = reference,
    field_checks = nrow(reports$semantic), compatibility_all_fields_pass = passed,
    strict_original_order_all_fields_pass = all(reports$strict_original$pass),
    strict_aligned_UTF8_all_fields_pass = all(reports$strict_aligned_UTF8$pass),
    compatibility_is_platform_equivalence = FALSE, continuous_strict_bound_is_diagnostic_only = TRUE,
    fit_seconds = saved$fit_seconds, native_generation_kind = saved$kind,
    new_native_calls = saved$new_native_calls, origin = saved$origin,
    DESeq2_native_diagnostics = diagnostics, Windows_DESeq2_diagnostics = old_diagnostics,
    DESeq2_diagnostics_exact = diagnostic_equal, DESeq2_diagnostics_structurally_valid = native_diagnostics_valid,
    convergence_or_fitType_forced_to_match_Windows = FALSE, mroast_native_FDR = roast,
    native_field_bound = "abs 1e-10 + rel 1e-8; positive P log10 delta <=1e-6; masks/directions exact",
    tail_uses_no_tolerance = TRUE, original_Windows_result_or_gate_changed = FALSE,
    full_execution_admitted_here = FALSE)
  if (id == 1L) bridge$domain <- jgwl_pin(domain_path)
  jgwl_after(z)
  jgwl_json(bridge, file.path(leaf, "BRIDGE_COMPLETE.json"))
  if (passed) 0L else 2L
}
jgwl_boot <- function(z) {
  prepared <- jgwl_prepared(z)
  bridges <- lapply(1:2, function(id) jgwl_bridge_record(z, id, prepared))
  domain <- readRDS(bridges[[1L]]$meta$domain$path)
  observed <- readRDS(bridges[[1L]]$meta$compact$path)
  second <- readRDS(bridges[[2L]]$meta$compact$path)
  for (space in names(z$keys)) z$h$jgw03c_tail_start(observed, domain, space, z$keys[[space]])
  context <- prepared$context
  context$label_registry <- context$registry  # Name consumed by the unchanged prefix writer.
  context$observed_files <- list(fixed = context$fixed, domain = bridges[[1L]]$meta$domain,
    observed = bridges[[1L]]$meta$compact)
  identity <- list(schema = "jgw03-separate-linux-scientific-identity-v2", prepared = prepared$pin,
    config = z$before$config, entry = z$before$entry, sources = z$before$sources,
    inputs = z$before$inputs, runtime = z$before$runtime, registry = context$registry,
    domain_sha256 = domain$sha256, Linux_bridge_records = lapply(bridges, `[[`, "record"),
    observed = bridges[[1L]]$meta$compact, fixed = context$fixed)
  store <- file.path(z$run, "derived", z$accession)
  if (!dir.exists(store)) jgwl_need(dir.create(store, recursive = TRUE), "Cannot create separate Linux store")
  b <- list(e = z$e, modules = z$modules, h = z$h, context = context, registry = z$registry,
    keys = z$keys, domain = domain, observed = observed, store = store,
    identity = identity, identity_sha256 = jgwl_digest(identity), entry = z$before$entry,
    Linux_bridge_compacts = list(observed, second), Linux_bridge_records = bridges)
  # The only replaced worker reuse callback returns newly fitted Linux bridge
  # results. It never calls the original Windows first-pilot lookup.
  z$e$jgwf_reuse_first <- function(b, key, out) {
    index <- match(key, b$keys$four_group_MC[1:2])
    if (is.na(index)) return(NULL)
    value <- b$Linux_bridge_compacts[[index]]
    b$e$jgwf_compact_valid(value, b, key)
    list(compact = value, provenance = list(kind = "REUSED_SAME_LINUX_BRIDGE_NO_FIT",
      input = b$Linux_bridge_records[[index]]$meta$compact,
      completion = b$Linux_bridge_records[[index]]$record),
      native_diagnostics = value$saved_DESeq2_metadata, new_native_calls = 0L)
  }
  z$e$jgwf_compact_valid(observed, b, z$keys$four_group_MC[1L])
  z$e$jgwf_compact_valid(second, b, z$keys$four_group_MC[2L])
  list(b = b, prepared = prepared)
}
jgwl_accept <- function(z, file, sha256) {
  record <- jgwl_pin(file)
  jgwl_need(identical(record$sha256, sha256) && jgwl_inside(record$path, z$run),
    "Explicit independent acceptance pin/root differs")
  x <- jsonlite::read_json(record$path, simplifyVector = FALSE)
  jgwl_need(identical(x$schema, "jgw03-linux-full-execution-acceptance-v2") &&
    identical(x$status, "ACCEPTED_FOR_FULL_LINUX_REFERENCE_V2") &&
    identical(x$config_sha256, z$before$config$sha256) &&
    identical(x$entry_sha256, z$before$entry$sha256) &&
    x$original_PC12_native_actual_exit_code == 2 &&
    isTRUE(x$remaining_five_producer_actual_exits_zero) &&
    isTRUE(x$six_independent_checks_actual_exits_zero) &&
    isTRUE(x$independent_bridge_acceptance) && isTRUE(x$resources_admitted) &&
    x$full_MC_family_n == 89 && x$full_companion_family_n == 34 &&
    x$MC_nonidentity_per_cell == 10000 && x$companion_per_cell == 20 &&
    setequal(names(x$cells), names(z$cfg$cells)), "Independent v2 full-scope acceptance required")
  checker <- z$before$sources[[z$cfg$checker_source]]
  for (accession in names(x$cells)) {
    cell <- x$cells[[accession]]
    prep <- jgwl_record(cell$prepared)
    jgwl_need(identical(prep$path, file.path(z$run, "intermediate", accession, "PREPARED.rds")),
      "Acceptance preparation association differs")
    for (id in 1:2) {
      pin <- jgwl_record(cell$bridges[[as.character(id)]])
      jgwl_need(identical(pin$path, file.path(z$run, "intermediate", accession,
        paste0("bridge_", id), "BRIDGE_COMPLETE.json")), "Acceptance bridge association differs")
      bridge <- jsonlite::read_json(pin$path, simplifyVector = FALSE)
      jgwl_need(identical(bridge$status, "LINUX_NATIVE_COMPATIBILITY_PASS_V2") &&
        isTRUE(bridge$compatibility_all_fields_pass) &&
        identical(bridge$prepared_sha256, prep$sha256) &&
        identical(bridge$config_sha256, z$before$config$sha256) &&
        identical(bridge$entry_sha256, z$before$entry$sha256), "Acceptance bridge identity differs")
      check <- cell$checks[[as.character(id)]]
      receipt <- jgwl_record(check$receipt); process <- jgwl_record(check$process)
      cmeta <- jsonlite::read_json(receipt$path, simplifyVector = FALSE)
      pmeta <- jsonlite::read_json(process$path, simplifyVector = FALSE)
      jgwl_need(identical(cmeta$status, "PASS_SAVED_NATIVE_COMPATIBILITY_V2_NOT_EQUIVALENCE") &&
        identical(cmeta$accession, accession) && cmeta$allocation_id == id &&
        identical(cmeta$config_sha256, z$before$config$sha256) &&
        identical(cmeta$producer_entry_sha256, z$before$entry$sha256) &&
        identical(cmeta$prepared_sha256, prep$sha256) &&
        identical(cmeta$native_sha256, bridge$native$sha256) &&
        identical(cmeta$checker_sha256, checker$sha256) &&
        pmeta$actual_R_exit_code == 0 && identical(pmeta$entry_sha256, checker$sha256),
        "Actual independent checker identity/exit differs")
      if (!(accession == "GSE289929" && id == 1L)) {
        producer <- jgwl_record(cell$producer_processes[[as.character(id)]])
        terminal <- jsonlite::read_json(producer$path, simplifyVector = FALSE)
        jgwl_need(terminal$actual_R_exit_code == 0 &&
          identical(terminal$entry_sha256, z$before$entry$sha256),
          "Remaining new producer actual exit zero required")
      } else jgwl_need(identical(bridge$native_generation_kind, "REUSED_ORIGINAL_LINUX_NATIVE") &&
        bridge$new_native_calls == 0 && bridge$origin$actual_original_R_exit_code == 2,
        "First native origin exit two must remain disclosed")
    }
  }
  for (name in c("independent_check", "resource_receipt", "process_evidence")) jgwl_record(x[[name]])
  record
}
jgwl_execute <- function(z, mode, space, first = NULL, last = NULL, acceptance = NULL) {
  jgwl_need(space %in% c("four_group_MC", "high_companion_exact"), "Original complete space required")
  if (mode == "run") {
    jgwl_need(is.list(acceptance), "Explicit acceptance required before any run")
    acceptance <- jgwl_accept(z, acceptance$path, acceptance$sha256)
  }
  opened <- jgwl_boot(z); b <- opened$b; keys <- b$keys[[space]]
  if (mode == "run") {
    jgwl_need(!is.null(acceptance) && length(first) == 1L && length(last) == 1L &&
      is.finite(first) && is.finite(last) && first == floor(first) && last == floor(last) &&
      first >= 1L && last >= first && last <= length(keys), "Accepted explicit bounded range required")
    fixed <- readRDS(b$context$fixed$path)
    matrix <- if (space == "four_group_MC") b$registry$MC else b$registry$high_companion
    for (i in seq.int(first, last)) {
      labels <- setNames(as.numeric(matrix[i, ]), colnames(matrix))
      b$modules$kernel$jgw03_check_allocation(labels, b$domain$sample_ids)
      b$e$jgwf_get_draw(b, keys[i], labels, fixed, z$cell_dir, TRUE)
      if (i %% 10L == 0L) invisible(gc())
    }
    prefix <- b$e$jgwf_tail_prefix(b, space, last)
    jgwl_record(acceptance); jgwl_after(z)
    cat("SEPARATE_LINUX_BOUNDED_RANGE_COMPLETE", z$accession, space, first, last, "\n")
    return(0L)
  }
  jgwl_need(identical(mode, "finish"), "Saved-only finish mode required")
  state <- b$e$jgwf_tail_prefix(b, space, length(keys))$state
  summary <- b$h$jgw03c_tail_finish(state,
    list(MC = b$h$jgw03_finite_mc_interval, companion = b$h$jgw03_exact20_companion_interval),
    b$modules$pins$finite$sha256)
  target <- file.path(b$store, paste0(space, "_TAIL_SUMMARY.rds"))
  b$e$jgwx_atomic_plain(list(summary = summary, tail_checkpoint = state,
    scientific_identity_sha256 = b$identity_sha256, full_registry = b$context$registry,
    execution_entry = b$entry), target)
  jgwl_after(z)
  cat("SEPARATE_LINUX_FULL_SAVED_REFERENCE_FINISHED_NOT_INDEPENDENT_ACCEPTANCE", z$accession, space, "\n")
  0L
}
jgwl_main <- function(args) {
  jgwl_need(length(args) >= 3L, "CONFIG MODE ACCESSION required")
  mode <- args[2L]
  jgwl_need(mode %in% c("prepare", "bridge", "run", "finish"), "Known explicit mode required")
  z <- jgwl_open(args[1L], args[3L])
  if (mode == "prepare") { jgwl_need(length(args) == 3L, "Prepare has three arguments"); return(jgwl_prepare(z)) }
  if (mode == "bridge") {
    jgwl_need(length(args) == 4L && args[4L] %in% c("1", "2"), "Bridge needs fixed ID 1 or 2")
    return(jgwl_bridge(z, as.integer(args[4L])))
  }
  if (mode == "finish") {
    jgwl_need(length(args) == 4L, "Finish needs original space")
    return(jgwl_execute(z, mode, args[4L]))
  }
  jgwl_need(length(args) == 8L && all(grepl("^[0-9]+$", args[5:6])),
    "Run needs SPACE FIRST LAST ACCEPTANCE_JSON ACCEPTANCE_SHA256")
  permit <- list(path = args[7L], sha256 = args[8L])
  jgwl_execute(z, mode, args[4L], as.numeric(args[5L]), as.numeric(args[6L]), permit)
}
if (sys.nframe() == 0L) {
  options(stringsAsFactors = FALSE, width = 140)
  exit_code <- tryCatch(jgwl_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    # Raw error bodies are not persisted or echoed. Retained capsules/partials
    # and structured status evidence remain for bounded no-fit diagnosis.
    message("JGW03_LINUX_STOP_RETAINED_NO_AUTOMATIC_REFIT")
    1L
  })
  quit(status = exit_code)
}
