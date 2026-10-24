#!/usr/bin/env Rscript
# Functions only when sourced. No dispatch/controller. Invoke an explicit one-
# process range only after the parent accepts the real resource measurements.
# Every MC family and full exact companion remain in scope; no faster kernel.
jgwf_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwf_digest <- function(x) digest::digest(x, "sha256", serializeVersion = 3L)
jgwf_boot <- function(root, out, accession) {
  path <- file.path(root, "code/R/prepare_JGW03_reference_first_nonidentity_20261015.R")
  jgwf_need(identical(digest::digest(file = path, algo = "sha256", serialize = FALSE),
    "99e4c614c0ddcbe6fd1970f6fc1fb4029731b812e14f92a6aafce43d63e9b82c"), "First-pilot entry source differs")
  e <- new.env(parent = asNamespace("stats")); sys.source(path, envir = e)
  prepared_json <- jsonlite::read_json(file.path(out, "PREPARED.json"), simplifyVector = FALSE)
  e$jgwx_assert_pin(file.path(out, "PREPARED.rds"), prepared_json$prepared)
  context <- readRDS(file.path(out, "PREPARED.rds"))
  modules <- e$jgwx_modules(root); h <- modules$helper
  jgwf_need(identical(context$schema, "jgw03-first-nonidentity-preparation-v1") &&
    identical(context$root, root) && identical(context$accession, accession) &&
    context$entry$sha256 %in% c(
      "99e4c614c0ddcbe6fd1970f6fc1fb4029731b812e14f92a6aafce43d63e9b82c",
      "49cf1d951602b2b4b78ce2dbe3a1aef3322c57cce03be1b78d509c9e1925ec5d") &&
    identical(context$sources, modules$pins) &&
    h$jgwp_json_equal(context$runtime, h$jgwp_runtime()), "Approved scientific preparation/runtime differs")
  e$jgwx_assert_pin(context$entry$path, context$entry[c("bytes", "sha256")])
  for (record in context$observed_files) e$jgwx_assert_pin(record$path, record[c("bytes", "sha256")])
  jgwf_need(identical(context$protected,
    h$jgwp_bound_records(root, accession, modules$admission)), "Protected scientific input differs")
  e$jgwx_assert_pin(file.path(out, "FULL_LABEL_REGISTRY.rds"), context$label_registry)
  registry <- readRDS(file.path(out, "FULL_LABEL_REGISTRY.rds"))
  keys <- e$jgwx_labels(registry)
  domain <- readRDS(context$observed_files$domain$path)
  observed <- readRDS(context$observed_files$observed$path)
  for (space in names(keys)) h$jgw03c_tail_start(observed, domain, space, keys[[space]])
  identity <- list(prepared_sha256 = prepared_json$prepared$sha256,
    registry = context$label_registry, domain_sha256 = domain$sha256,
    sources = context$sources, protected = context$protected,
    observed_files = context$observed_files, preparation_entry = context$entry, runtime = context$runtime,
    fixed_ram_sha256_is_checked_by_admission_before_every_new_fit = TRUE)
  store <- file.path(out, "full_reference_compacts_01")
  if (!dir.exists(store)) dir.create(store)
  list(e = e, modules = modules, h = h, context = context, registry = registry,
    keys = keys, domain = domain, observed = observed, store = store,
    identity = identity, identity_sha256 = jgwf_digest(identity), entry = h$jgwp_entry_pin())
}
jgwf_native_diagnostics <- function(value) {
  lapply(value$gene_DESeq2, function(x) {
    dds <- x$native_locals$dds
    flags <- S4Vectors::mcols(dds)$betaConv
    list(fitType = attr(DESeq2::dispersionFunction(dds), "fitType"),
      betaConv = if (is.null(flags)) list(status = "NOT_RETAINED", false_n = NULL, unavailable_n = NULL) else
        list(status = "SAVED_NATIVE_FLAG_COUNTS", n = length(flags),
          false_n = sum(!is.na(flags) & !flags), unavailable_n = sum(is.na(flags))),
      native_results_metadata = "Original res converted to data.frame; no new results() reconstruction")
  })
}
jgwf_compact_valid <- function(x, b, key) {
  jgwf_need(identical(x$schema, "jgw03-native-compact-v1") &&
    identical(x$accession, b$context$accession) && identical(x$domain_sha256, b$domain$sha256) &&
    identical(x$allocation_key, key) && identical(names(x$families), names(b$domain$families)) &&
    identical(x$definition, b$domain$definition), "Complete fixed compact domain/order differs")
  for (id in names(b$domain$families)) {
    spec <- b$domain$families[[id]]; family <- x$families[[id]]
    p <- family$native_numeric[[spec$map$native_P]]
    q <- family$native_numeric[[spec$map$native_q]]
    invalid <- function(v) sum(!is.na(v) & (!is.finite(v) | v < 0 | v > 1))
    jgwf_need(identical(names(family$native_numeric), spec$numeric_names) &&
      identical(names(family$native_text), spec$dynamic_names) &&
      identical(vapply(family$native_numeric, typeof, character(1)), spec$numeric_types) &&
      all(lengths(family$native_numeric) == length(spec$keys)) &&
      all(lengths(family$native_text) == length(spec$keys)) &&
      length(family$outer_statistic) == length(spec$keys) &&
      is.logical(family$outer_available) && !anyNA(family$outer_available) &&
      length(family$outer_available) == length(spec$keys) &&
      family$diagnostics$raw_P_invalid_n == 0L && family$diagnostics$raw_q_invalid_n == 0L,
      "Complete native vectors/probability domain differ; NA is retained separately")
    jgwf_need(is.numeric(p) && is.numeric(q) && invalid(p) == 0L && invalid(q) == 0L &&
      identical(family$outer_available, is.finite(p) & p >= 0 & p <= 1 & !is.na(family$outer_statistic)),
      "Native probability/availability domain differs; no imputation or missing-row deletion")
  }
  invisible(TRUE)
}
jgwf_read_archive <- function(b, key) {
  leaf <- file.path(b$store, key)
  file <- file.path(leaf, "COMPACT.rds"); meta <- file.path(leaf, "COMPLETE.json")
  if (!dir.exists(leaf)) return(NULL)
  jgwf_need(file.exists(file) && file.exists(meta), "Incomplete allocation retained; no automatic refit")
  receipt <- jsonlite::read_json(meta, simplifyVector = FALSE)
  b$e$jgwx_assert_pin(file, receipt$archive)
  value <- readRDS(file)
  jgwf_need(identical(value$schema, "jgw03-full-reference-allocation-v1") &&
    identical(value$scientific_identity_sha256, b$identity_sha256) &&
    identical(value$allocation_key, key) && identical(value$compact_sha256, jgwf_digest(value$compact)),
    "Checkpoint scientific/content identity differs")
  jgwf_compact_valid(value$compact, b, key)
  value
}
jgwf_reuse_first <- function(b, key, out) {
  if (key == b$observed$allocation_key) return(list(compact = b$observed,
    provenance = list(kind = "REUSED_ACCEPTED_OBSERVED_NO_FIT", input = b$context$observed_files$observed),
    native_diagnostics = b$observed$saved_DESeq2_metadata, new_native_calls = 0L))
  if (key != b$keys$four_group_MC[2L]) return(NULL)
  leaf <- file.path(out, "allocation_00002")
  meta <- jsonlite::read_json(file.path(leaf, "PILOT_COMPLETE.json"), simplifyVector = FALSE)
  file <- file.path(leaf, "COMPACT.rds")
  jgwf_need(identical(meta$status, "FIRST_NONIDENTITY_SAVED_AND_COMPACTED_ONLY") &&
    identical(meta$allocation_key, key) && identical(meta$domain_sha256, b$domain$sha256) &&
    identical(meta$prepared_sha256, b$identity$prepared_sha256) &&
    !isTRUE(meta$original_label_refitted), "Accepted first nonidentity checkpoint required")
  b$e$jgwx_assert_pin(file, meta$compact)
  # Its actual terminal process is an independent completion requirement.
  process_file <- paste0(out, "_pilot_01_process/PROCESS_EXIT.json")
  terminal <- jsonlite::read_json(process_file, simplifyVector = FALSE)
  jgwf_need(terminal$actual_R_exit_code == 0 && identical(terminal$mode, "pilot") &&
    identical(terminal$accession, b$context$accession) &&
    identical(normalizePath(terminal$output_leaf, winslash = "/", mustWork = TRUE), out) &&
    identical(terminal$entry_sha256, b$context$entry$sha256), "Actual same-cell first-pilot exit 0 required")
  value <- readRDS(file)
  list(compact = value, provenance = list(kind = "REUSED_ACCEPTED_MC_ID2_NO_FIT",
    input = c(list(path = file), meta$compact), completion = b$e$jgwx_pin(file.path(leaf, "PILOT_COMPLETE.json")),
    process = b$e$jgwx_pin(process_file)), native_diagnostics = value$saved_DESeq2_metadata,
    new_native_calls = 0L)
}
jgwf_get_draw <- function(b, key, allocation, fixed, out, allow_fit) {
  prior <- jgwf_read_archive(b, key)
  if (!is.null(prior)) return(prior)
  retained <- jgwf_reuse_first(b, key, out)
  if (is.null(retained)) jgwf_need(isTRUE(allow_fit), "Missing saved draw; finish never fits")
  leaf <- file.path(b$store, key)
  jgwf_need(dir.create(leaf, showWarnings = FALSE), "Allocation already owned/incomplete; no concurrent or repeated fit")
  b$e$jgwx_json(list(status = "ALLOCATION_OWNED_BEFORE_NATIVE_CALL", allocation_key = key,
    scientific_identity_sha256 = b$identity_sha256, execution_entry = b$entry,
    new_fit_intended = is.null(retained), utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(leaf, "STARTED.json"))
  if (is.null(retained)) {
    do.call(RNGkind, as.list(b$context$runtime$RNGkind)); set.seed(b$context$runtime$seed)
    start <- proc.time()[["elapsed"]]
    value <- b$modules$kernel$jgw03_fit_one_allocation(fixed, allocation, b$modules$pins$native$path)
    fit_seconds <- proc.time()[["elapsed"]] - start
    execution <- list(seed = b$context$runtime$seed, RNGkind = b$context$runtime$RNGkind,
      backend = "SerialParam", workers = 1L, mroast_nrot = 1999L, fgsea_eps = 0,
      mroast_applicable = b$context$accession == "GSE289929",
      mroast_executed = b$context$accession == "GSE289929",
      actual_mroast_nrot = if (b$context$accession == "GSE289929") 1999L else NULL,
      mroast_configuration_provenance = if (b$context$accession == "GSE289929")
        "Actual PC12 Reactome mroast" else "Unused source-config-only nrot=1999; no native mroast",
      provenance = paste("Actual single-worker full native call; scientific identity", b$identity_sha256),
      RNG_after_sha256 = jgwf_digest(.Random.seed), historical_backend_authenticated = FALSE)
    retained <- tryCatch({
      compact <- b$h$jgw03c_extract(value, b$domain, execution)
      jgwf_compact_valid(compact, b, key)
      list(compact = compact, provenance = list(kind = "NEW_UNCHANGED_NATIVE_KERNEL", fit_seconds = fit_seconds,
        input_fixed_fingerprint = fixed$input_admission_v2$fixed_ram_sha256),
        native_diagnostics = jgwf_native_diagnostics(value), new_native_calls = 1L)
      }, error = function(condition) {
        # Bounded exceptional capsule only, never every successful allocation.
        saveRDS(value, file.path(leaf, "POST_FIT_FAILURE_NATIVE_CAPSULE.rds"), compress = "gzip", version = 3L)
        stop("Post-fit extraction retained native capsule; stop without refit", call. = FALSE)
      })
  }
  jgwf_compact_valid(retained$compact, b, key)
  archive <- list(schema = "jgw03-full-reference-allocation-v1", scientific_identity_sha256 = b$identity_sha256,
    accession = b$context$accession, allocation_key = key, compact = retained$compact,
    compact_sha256 = jgwf_digest(retained$compact), provenance = retained$provenance,
    native_diagnostics = retained$native_diagnostics, new_native_calls = retained$new_native_calls,
    execution_entry = b$entry)
  completed <- tryCatch({
    pin <- b$e$jgwx_atomic_plain(archive, file.path(leaf, "COMPACT.rds"))
    b$e$jgwx_json(list(status = "IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT",
      allocation_key = key, scientific_identity_sha256 = b$identity_sha256,
      archive = pin, new_native_calls = retained$new_native_calls), file.path(leaf, "COMPLETE.json"))
    jgwf_read_archive(b, key)
  }, error = function(condition) {
    if (retained$new_native_calls == 1L) saveRDS(value,
      file.path(leaf, "POST_FIT_FAILURE_NATIVE_CAPSULE.rds"), compress = "gzip", version = 3L)
    stop("Post-fit checkpoint error retained; stop without automatic refit", call. = FALSE)
  })
  if (retained$new_native_calls == 1L) { rm(value); invisible(gc()) }
  completed
}
jgwf_tail_prefix <- function(b, space, last) {
  keys <- b$keys[[space]]
  state <- b$h$jgw03c_tail_start(b$observed, b$domain, space, keys)
  pins <- list()
  if (last >= 2L) for (i in seq.int(2L, last)) {
    archive <- jgwf_read_archive(b, keys[i])
    jgwf_need(!is.null(archive), "Prefix incomplete; no count/denominator shrink")
    state <- b$h$jgw03c_tail_update(state, archive$compact)
    pins[[as.character(i)]] <- list(allocation_key = keys[i],
      archive = b$e$jgwx_pin(file.path(b$store, keys[i], "COMPACT.rds")))
  }
  checkpoint <- list(schema = "jgw03-fixed-registry-tail-prefix-v1", space = space,
    through_allocation_id = last, scientific_identity_sha256 = b$identity_sha256,
    full_registry = b$context$label_registry, observed = b$context$observed_files$observed,
    completed_compact_pins = pins, state = state, provisional_P_computed = FALSE)
  target <- file.path(b$store, sprintf("%s_TAIL_PREFIX_%05d.rds", space, last))
  if (file.exists(target)) jgwf_need(identical(serialize(readRDS(target), NULL, version = 3L),
    serialize(checkpoint, NULL, version = 3L)), "Existing immutable tail prefix differs") else
      b$e$jgwx_atomic_plain(checkpoint, target)
  checkpoint
}
jgwf_check_after <- function(b, root, accession) {
  b$e$jgwx_assert_pin(b$context$entry$path, b$context$entry[c("bytes", "sha256")])
  for (record in b$context$observed_files) b$e$jgwx_assert_pin(record$path, record[c("bytes", "sha256")])
  jgwf_need(identical(b$context$protected,
    b$h$jgwp_bound_records(root, accession, b$modules$admission)) &&
    identical(b$context$sources, b$e$jgwx_modules(root)$pins) &&
    identical(b$entry, b$h$jgwp_entry_pin()), "Protected inputs/source/execution changed")
}
jgwf_main <- function(args) {
  jgwf_need(length(args) %in% c(5L, 7L), "ROOT PREPARED_CELL ACCESSION run SPACE FIRST LAST | finish SPACE")
  root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
  out <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
  accession <- args[[3L]]; mode <- args[[4L]]; space <- args[[5L]]
  jgwf_need(accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    space %in% c("four_group_MC", "high_companion_exact") &&
    startsWith(out, "R:/Shaoyao_storage/cache/JGW03_label_reference_20261015r01/") &&
    startsWith(normalizePath(tempdir(), winslash = "/"), "R:/Shaoyao_storage/"), "Known cell/space and R output/temp required")
  b <- jgwf_boot(root, out, accession)
  matrix <- if (space == "four_group_MC") b$registry$MC else b$registry$high_companion
  keys <- b$keys[[space]]
  if (mode == "run") {
    jgwf_need(length(args) == 7L && all(grepl("^[0-9]+$", args[6:7])), "Explicit integer bounded range required")
    ids <- as.integer(args[6:7])
    jgwf_need(all(is.finite(ids)) && ids[1L] >= 1L && ids[2L] >= ids[1L] && ids[2L] <= length(keys), "Range exceeds frozen space")
    fixed <- readRDS(b$context$observed_files$fixed$path)
    for (i in seq.int(ids[1L], ids[2L])) {
      labels <- setNames(as.numeric(matrix[i, ]), colnames(matrix))
      b$modules$kernel$jgw03_check_allocation(labels, b$domain$sample_ids)
      jgwf_get_draw(b, keys[i], labels, fixed, out, TRUE)
      if (i %% 10L == 0L) invisible(gc())
    }
    jgwf_tail_prefix(b, space, ids[2L])
    jgwf_check_after(b, root, accession)
    cat("BOUNDED_SINGLE_WORKER_RANGE_COMPLETE", accession, space, ids, "\n")
    return(0L)
  }
  jgwf_need(mode == "finish" && length(args) == 5L, "Explicit saved-only finish required")
  # Replay immutable full-registry checkpoints; never a smaller denominator.
  state <- jgwf_tail_prefix(b, space, length(keys))$state
  summary <- b$h$jgw03c_tail_finish(state,
    list(MC = b$h$jgw03_finite_mc_interval, companion = b$h$jgw03_exact20_companion_interval),
    b$modules$pins$finite$sha256)
  target <- file.path(b$store, paste0(space, "_TAIL_SUMMARY.rds"))
  b$e$jgwx_atomic_plain(list(summary = summary, tail_checkpoint = state,
    scientific_identity_sha256 = b$identity_sha256, full_registry = b$context$label_registry,
    execution_entry = b$entry), target)
  jgwf_check_after(b, root, accession)
  cat("FULL_REGISTRY_SAVED_ONLY_TAIL_FINISHED_NOT_INDEPENDENT_ACCEPTANCE", accession, space, "\n")
  0L
}
if (sys.nframe() == 0L) {
  options(stringsAsFactors = FALSE, width = 140)
  code <- tryCatch(jgwf_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_FULL_SINGLE_STOP_RETAINED_NO_AUTOMATIC_REFIT")
    1L
  })
  quit(status = code)
}
