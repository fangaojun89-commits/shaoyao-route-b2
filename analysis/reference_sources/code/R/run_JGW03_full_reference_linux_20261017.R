#!/usr/bin/env Rscript
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
  jgwl_need(startsWith(run, "/root/autodl-tmp/") &&
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
    "code/R/run_JGW03_full_reference_linux_20261017.R")
  jgwl_need(all(required %in% relative), "Required immutable source absent from manifest")
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
    windows_receipts = receipts, entry = entry, runtime = runtime)
  list(cfg = cfg, project = project, run = run, accession = accession, e = e, modules = modules,
    h = modules$helper, ref = ref, metadata = metadata, registry = registry, keys = keys,
    cell_dir = file.path(run, "intermediate", accession), before = before)
}
jgwl_after <- function(z) {
  records <- c(list(z$before$config, z$before$entry), z$before$sources, z$before$inputs,
    z$before$references, z$before$windows_receipts)
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
  jgwl_json(list(schema = "jgw03-linux-preparation-started-v1", binding = z$before),
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
  rm(old); invisible(gc())
  check_pin <- jgwl_json(comparison, file.path(z$cell_dir, "FIXED_LAYER_BRIDGE.json"))
  jgwl_need(all(unlist(comparison)), "Fixed layer bridge requires review; no allocation fit performed")
  registry_path <- file.path(z$cell_dir, "FULL_LABEL_REGISTRY.rds")
  jgwl_need(!file.exists(registry_path) && file.copy(z$ref$registry$path, registry_path, overwrite = FALSE),
    "Registry byte-copy failed")
  registry_pin <- jgwl_record(c(list(path = registry_path), z$ref$registry[c("bytes", "sha256")]))
  session_path <- file.path(z$cell_dir, "sessionInfo.txt")
  utils::capture.output(utils::sessionInfo(), file = session_path)
  session_pin <- jgwl_pin(session_path)
  context <- list(schema = "jgw03-separate-linux-preparation-v1", accession = z$accession,
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
  jgwl_need(identical(context$schema, "jgw03-separate-linux-preparation-v1") &&
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
  jgwl_need(identical(meta$status, "LINUX_NATIVE_BRIDGE_PASS_ONLY") && meta$allocation_id == id &&
    identical(meta$allocation_key, key) && identical(meta$prepared_sha256, prepared$pin$sha256) &&
    identical(meta$config_sha256, z$before$config$sha256) &&
    identical(meta$entry_sha256, z$before$entry$sha256), "Accepted same-Linux bridge required")
  for (name in c("native", "compact", "checks", "checks_saved", "session_info")) jgwl_record(meta[[name]])
  if (id == 1L) jgwl_record(meta$domain)
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
      fit_seconds = elapsed, execution = execution, session_info = jgwl_pin(session_path))
    jgwl_json(saved, saved_path)
  }
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
  checks <- lapply(names(prior), function(name) {
    old_table <- prior[[name]]; table <- fresh[[name]]
    key <- if ("raw_gene_id" %in% names(old_table)) "raw_gene_id" else
      if ("term_id" %in% names(old_table)) "term_id" else "feature"
    result <- z$h$jgwp_table_check(table, old_table, key, name)
    excluded <- c("ql_approx_se", "ql_approx_CI_low", "ql_approx_CI_high")
    old_types <- vapply(old_table[!names(old_table) %in% excluded], typeof, character(1L))
    types <- vapply(table[!names(table) %in% excluded], typeof, character(1L))
    result$schema_types_exact <- identical(types, old_types)
    result$pass <- result$pass & result$schema_types_exact
    result
  })
  checks <- do.call(rbind, checks)
  check_path <- file.path(leaf, "NATIVE_FIELD_BRIDGE.tsv")
  checks_saved_path <- file.path(leaf, "CHECKS_SAVED.json")
  if (file.exists(check_path)) {
    # A text diagnostic may round its max-difference columns. Resume instead
    # authenticates its bytes AND the regenerated full-precision checks object;
    # native-field and exact-tail comparison rules are not changed.
    jgwl_need(file.exists(checks_saved_path), "Diagnostic write incomplete; retain saved native without refit")
    prior_checks <- jsonlite::read_json(checks_saved_path, simplifyVector = FALSE)
    check_pin <- jgwl_record(prior_checks$file)
    jgwl_need(identical(check_pin$path, check_path) &&
      identical(prior_checks$full_precision_checks_sha256, jgwl_digest(checks)) &&
      identical(prior_checks$native_sha256, saved$native$sha256) &&
      identical(prior_checks$prepared_sha256, prepared$pin$sha256), "Retained diagnostic/native coupling differs")
  } else {
    jgwl_need(!file.exists(checks_saved_path), "Orphaned diagnostic receipt retained")
    utils::write.table(checks, check_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    check_pin <- jgwl_pin(check_path)
    jgwl_json(list(file = check_pin, full_precision_checks_sha256 = jgwl_digest(checks),
      native_sha256 = saved$native$sha256, prepared_sha256 = prepared$pin$sha256), checks_saved_path)
  }
  diagnostics <- z$e$jgwf_native_diagnostics(value)
  old_diagnostics <- z$e$jgwf_native_diagnostics(old)
  diagnostic_equal <- z$h$jgwp_json_equal(diagnostics, old_diagnostics)
  roast <- z$h$jgwp_roast_check(value)
  passed <- all(checks$pass) && diagnostic_equal && (!isTRUE(roast$applicable) || isTRUE(roast$pass))
  bridge <- list(status = if (passed) "LINUX_NATIVE_BRIDGE_PASS_ONLY" else "LINUX_NATIVE_BRIDGE_REVIEW_REQUIRED",
    accession = z$accession, allocation_id = id, allocation_key = compact$allocation_key,
    prepared_sha256 = prepared$pin$sha256, config_sha256 = z$before$config$sha256,
    entry_sha256 = z$before$entry$sha256, native = saved$native, compact = jgwl_pin(compact_path),
    checks = check_pin, checks_saved = jgwl_pin(checks_saved_path), native_family_n = length(compact$families),
    session_info = saved$session_info,
    Windows_reference_only_never_used_as_Linux_result = reference,
    field_checks = nrow(checks), all_fields_pass = passed, fit_seconds = saved$fit_seconds,
    DESeq2_native_diagnostics = diagnostics, Windows_DESeq2_diagnostics = old_diagnostics,
    DESeq2_diagnostics_exact = diagnostic_equal, mroast_native_FDR = roast,
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
  identity <- list(schema = "jgw03-separate-linux-scientific-identity-v1", prepared = prepared$pin,
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
  jgwl_need(identical(x$schema, "jgw03-linux-full-execution-acceptance-v1") &&
    identical(x$status, "ACCEPTED_FOR_FULL_LINUX_REFERENCE") &&
    identical(x$config_sha256, z$before$config$sha256) &&
    identical(x$entry_sha256, z$before$entry$sha256) &&
    isTRUE(x$all_six_bridge_actual_exits_zero) && isTRUE(x$independent_bridge_acceptance) &&
    isTRUE(x$resources_admitted) && x$full_MC_family_n == 89 && x$full_companion_family_n == 34 &&
    x$MC_nonidentity_per_cell == 10000 && x$companion_per_cell == 20 &&
    setequal(names(x$cells), names(z$cfg$cells)), "Independent full-scope bridge/resource acceptance required")
  for (accession in names(x$cells)) {
    cell <- x$cells[[accession]]
    prep <- jgwl_record(cell$prepared)
    expected <- file.path(z$run, "intermediate", accession, "PREPARED.rds")
    jgwl_need(identical(prep$path, expected), "Acceptance preparation association differs")
    for (id in 1:2) {
      pin <- jgwl_record(cell$bridges[[as.character(id)]])
      jgwl_need(identical(pin$path, file.path(z$run, "intermediate", accession,
        paste0("bridge_", id), "BRIDGE_COMPLETE.json")), "Acceptance bridge association differs")
      bridge <- jsonlite::read_json(pin$path, simplifyVector = FALSE)
      jgwl_need(identical(bridge$status, "LINUX_NATIVE_BRIDGE_PASS_ONLY") &&
        identical(bridge$prepared_sha256, prep$sha256) &&
        identical(bridge$config_sha256, z$before$config$sha256) &&
        identical(bridge$entry_sha256, z$before$entry$sha256), "Acceptance bridge identity differs")
    }
  }
  for (name in c("independent_check", "resource_receipt", "process_evidence")) jgwl_record(x[[name]])
  record
}
jgwl_execute <- function(z, mode, space, first = NULL, last = NULL, acceptance = NULL) {
  jgwl_need(space %in% c("four_group_MC", "high_companion_exact"), "Original complete space required")
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
  permit <- jgwl_accept(z, args[7L], args[8L])
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
