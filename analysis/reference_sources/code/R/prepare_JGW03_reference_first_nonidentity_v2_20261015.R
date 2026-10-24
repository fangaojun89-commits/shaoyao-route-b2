#!/usr/bin/env Rscript
# Source safe. Prepare the COMPLETE frozen label spaces, then an explicitly
# invoked single-worker first-nonidentity pilot. No automatic full dispatch.
jgwx_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwx_pin <- function(path) list(bytes = unname(file.info(path)$size),
  sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
jgwx_assert_pin <- function(path, pin) {
  actual <- jgwx_pin(path)
  jgwx_need(file.exists(path) && is.list(pin) && is.numeric(pin$bytes) &&
    length(pin$bytes) == 1L && is.finite(pin$bytes) &&
    actual$bytes == pin$bytes && identical(actual$sha256, pin$sha256), "Bound file identity differs")
  invisible(TRUE)
}
jgwx_json <- function(value, path) {
  jgwx_need(!file.exists(path), "New JSON only")
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE,
    digits = NA, na = "null", null = "null")
}
jgwx_atomic_plain <- function(value, path) {
  partial <- paste0(path, ".partial")
  jgwx_need(!file.exists(path) && !file.exists(partial), "New plain checkpoint only")
  saveRDS(value, partial, compress = "gzip", version = 3L)
  jgwx_need(identical(serialize(value, NULL, version = 3L),
    serialize(readRDS(partial), NULL, version = 3L)), "Plain checkpoint transport differs")
  jgwx_need(file.rename(partial, path), "Checkpoint rename failed")
  jgwx_pin(path)
}
jgwx_modules <- function(root) {
  relative <- c(pilot = "code/R/prepare_JGW03_observed_label_pilot_v3_20261015.R",
    core = "code/R/JGW03_label_reference_core_20261015.R",
    admission = "code/R/JGW03_input_admission_v2_20261015.R",
    map = "code/R/JGW03_reference_field_map_20261015.R",
    compact = "code/R/JGW03_reference_compact_tail_20261015.R",
    finite = "code/R/JGW03_finite_population_mc_v2_20261015.R",
    native = "code/R/analyze_jgw03_three_cell_lines.R")
  expected <- c("fe01d32bf79239cc51a5b4b75140c3ad209b31501adc94f141623c9ecf991f38",
    "e3939f371f02b24cfdfefee0a9cd9e39bdf47c8920caa44f58434a93f7a2974e",
    "63a24a37bac60e1e1c99cd11e9a69c6a11ab8dbb70c44766d06f2e33b87e9c0c",
    "c349eab9c5f3718d95a06c0b51e44584f4adb60ecbfb2f2d75235af7bfd7da42",
    "1074cc5001899cf4f2500c72a52f6f218b1cbb763982c6a0a9b43bfeab82ce8d",
    "349a5de42e47c41c7211c9f64c0f012d3658c72d840f5b75e4deb8912d52db6b",
    "5df6378ed4dadfb4050d91e252a63af4f44d2be80d3070d889d1e0e007796775")
  pins <- lapply(seq_along(relative), function(i) {
    path <- file.path(root, relative[[i]])
    pin <- jgwx_pin(path)
    jgwx_need(identical(pin$sha256, expected[[i]]), "Approved helper source differs")
    c(list(path = path), pin)
  }); names(pins) <- names(relative)
  helper <- new.env(parent = asNamespace("stats"))
  for (name in setdiff(names(relative), c("native", "admission")))
    sys.source(pins[[name]]$path, envir = helper)
  admission <- new.env(parent = asNamespace("stats"))
  sys.source(pins$admission$path, envir = admission)
  kernel <- admission$jgw03_v2_bind_kernel(pins$core$path)
  list(helper = helper, admission = admission, kernel = kernel, pins = pins)
}
jgwx_labels <- function(registry) {
  jgwx_need(is.list(registry) && is.matrix(registry$MC) &&
    identical(dim(registry$MC), c(10001L, 12L)) &&
    is.matrix(registry$high_companion) && identical(dim(registry$high_companion), c(20L, 12L)),
    "Full frozen MC and exact companion registries required")
  key <- function(matrix) apply(matrix, 1L, paste, collapse = "")
  list(four_group_MC = unname(key(registry$MC)),
    high_companion_exact = unname(key(registry$high_companion)))
}
jgwx_bind_observed <- function(root, accession, readback_file, readback_sha, modules) {
  jgwx_need(identical(jgwx_pin(readback_file)$sha256, readback_sha), "Explicit compact readback SHA differs")
  readback <- jsonlite::read_json(readback_file, simplifyVector = FALSE)
  jgwx_need(identical(readback$accession, accession) &&
    grepl("^PASS_SAVED_.*COMPACT_EXACT_VECTOR_READBACK_ONLY$", readback$status) &&
    isTRUE(readback$all_vectors_exact) && isTRUE(readback$input_pins_unchanged) &&
    isTRUE(readback$domain_compact_roundtrip_raw_serialization_exact) && readback$new_fit_calls == 0,
    "Accepted saved-native compact readback required")
  leaf <- file.path("R:/Shaoyao_storage/cache/JGW03_label_reference_20261015r01",
    paste0(accession, if (accession == "GSE289929") "_observe_02" else "_observe_01"))
  completion <- file.path(leaf, "OBSERVED_COMPLETE.json")
  process <- paste0(leaf, "_process/PROCESS_EXIT.json")
  jgwx_assert_pin(completion, readback$before$completion)
  jgwx_assert_pin(process, readback$before$process)
  meta <- jsonlite::read_json(completion, simplifyVector = FALSE)
  terminal <- jsonlite::read_json(process, simplifyVector = FALSE)
  jgwx_need(identical(meta$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_ONLY") &&
    identical(meta$accession, accession) && isTRUE(meta$bridge$all_fields_pass) &&
    terminal$actual_R_exit_code == 0 && identical(terminal$entry_sha256, meta$entry$sha256),
    "Successful actual original-label bridge required")
  domain_file <- file.path(dirname(readback_file), "COMPACT_DOMAIN.rds")
  compact_file <- file.path(dirname(readback_file), "OBSERVED_COMPACT.rds")
  fixed_file <- file.path(leaf, "FIXED_LAYER.rds")
  saved_meta <- jsonlite::read_json(file.path(leaf, "FIT_SAVED.json"), simplifyVector = FALSE)
  files <- list(readback = c(list(path = readback_file), jgwx_pin(readback_file)),
    completion = c(list(path = completion), jgwx_pin(completion)),
    process = c(list(path = process), jgwx_pin(process)),
    domain = c(list(path = domain_file), readback$outputs$domain),
    observed = c(list(path = compact_file), readback$outputs$compact),
    fixed = c(list(path = fixed_file), saved_meta$fixed_layer))
  for (record in files) jgwx_assert_pin(record$path, record[c("bytes", "sha256")])
  # Full old-native capsule is NOT decoded here; only the accepted compact and
  # label-independent fixed layer are needed for future allocations.
  domain <- readRDS(domain_file); observed <- readRDS(compact_file); fixed <- readRDS(fixed_file)
  jgwx_need(identical(domain$accession, accession) && identical(observed$accession, accession) &&
    identical(fixed$accession, accession) && identical(names(fixed$observed_labels), domain$sample_ids) &&
    identical(fixed$observed_labels, domain$observed_labels) &&
    identical(fixed$input_admission_v2$fixed_ram_sha256,
      modules$admission$jgw03_v2_fingerprint(fixed)) &&
    identical(observed$domain_sha256, domain$sha256), "Fixed/compact/domain identity differs")
  list(files = files, meta = meta, domain = domain, observed = observed, fixed = fixed)
}
jgwx_prepare <- function(root, out, accession, readback_file, readback_sha) {
  modules <- jgwx_modules(root); h <- modules$helper
  runtime <- h$jgwp_runtime()
  bound <- jgwx_bind_observed(root, accession, readback_file, readback_sha, modules)
  jgwx_need(h$jgwp_json_equal(runtime, bound$meta$runtime), "Observed/first-pilot runtime differs")
  protected <- h$jgwp_bound_records(root, accession, modules$admission)
  jgwx_need(h$jgwp_json_equal(protected, bound$meta$protected_inputs_after), "Original protected records differ")
  registry <- modules$kernel$jgw03_allocation_registry(bound$domain$observed_labels)
  keys <- jgwx_labels(registry)
  # These starts validate ALL 10001/20 keys, identity and fixed nuisance roles;
  # they do not evaluate any nonidentity scientific value or make a tail P.
  starts <- lapply(names(keys), function(space) h$jgw03c_tail_start(
    bound$observed, bound$domain, space, keys[[space]])); names(starts) <- names(keys)
  dir.create(out)
  registry_pin <- jgwx_atomic_plain(registry, file.path(out, "FULL_LABEL_REGISTRY.rds"))
  jgwx_atomic_plain(starts, file.path(out, "TAIL_STARTS.rds"))
  for (space in names(keys)) {
    matrix <- if (space == "four_group_MC") registry$MC else registry$high_companion
    utils::write.table(data.frame(allocation_id = seq_len(nrow(matrix)), allocation_key = keys[[space]],
      matrix, check.names = FALSE), file.path(out, paste0(space, "_labels.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  field_map <- h$jgw03_reference_field_map(modules$kernel$jgw03_family_registry())
  utils::write.table(field_map, file.path(out, "FULL_123_FAMILY_MAP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  context <- list(schema = "jgw03-first-nonidentity-preparation-v1", root = root, accession = accession,
    sources = modules$pins, entry = h$jgwp_entry_pin(), runtime = runtime,
    observed_files = bound$files, protected = protected,
    label_registry = registry_pin, domain_sha256 = bound$domain$sha256,
    full_MC_nonidentity_draws = 10000L, full_companion_including_observed = 20L,
    full_global_method_family_n = 89L, full_companion_method_family_n = 34L,
    current_cell_method_family_n = length(bound$domain$families),
    pilot_only_allocation_id = 2L, full_reference_resource_admitted = FALSE,
    scientific_authority_is_external_not_created_by_this_object = TRUE)
  jgwx_atomic_plain(context, file.path(out, "PREPARED.rds"))
  for (record in bound$files) jgwx_assert_pin(record$path, record[c("bytes", "sha256")])
  jgwx_need(identical(protected, h$jgwp_bound_records(root, accession, modules$admission)), "Protected inputs changed")
  jgwx_json(c(context, list(status = "PREPARED_COMPLETE_SPACES_NO_NONIDENTITY_FIT",
    prepared = jgwx_pin(file.path(out, "PREPARED.rds")), new_native_fit_calls = 0L)), file.path(out, "PREPARED.json"))
  0L
}
jgwx_pilot <- function(root, out, accession) {
  context_path <- file.path(out, "PREPARED.rds")
  complete <- jsonlite::read_json(file.path(out, "PREPARED.json"), simplifyVector = FALSE)
  jgwx_assert_pin(context_path, complete$prepared)
  context <- readRDS(context_path)
  jgwx_need(identical(context$root, root) && identical(context$accession, accession) &&
    identical(context$schema, "jgw03-first-nonidentity-preparation-v1") &&
    identical(context$pilot_only_allocation_id, 2L), "Pinned first-allocation context required")
  modules <- jgwx_modules(root); h <- modules$helper
  jgwx_need(identical(context$sources, modules$pins) && identical(context$entry, h$jgwp_entry_pin()) &&
    h$jgwp_json_equal(context$runtime, h$jgwp_runtime()), "Source/entry/runtime changed")
  for (record in context$observed_files) jgwx_assert_pin(record$path, record[c("bytes", "sha256")])
  jgwx_need(identical(context$protected,
    h$jgwp_bound_records(root, accession, modules$admission)), "Protected input changed")
  registry_path <- file.path(out, "FULL_LABEL_REGISTRY.rds")
  jgwx_assert_pin(registry_path, context$label_registry)
  registry <- readRDS(registry_path)
  keys <- jgwx_labels(registry)
  domain <- readRDS(context$observed_files$domain$path)
  fixed <- readRDS(context$observed_files$fixed$path)
  allocation <- setNames(as.numeric(registry$MC[2L, ]), colnames(registry$MC))
  modules$kernel$jgw03_check_allocation(allocation, domain$sample_ids)
  jgwx_need(paste(allocation, collapse = "") == keys$four_group_MC[2L] &&
    keys$four_group_MC[2L] != keys$four_group_MC[1L], "Frozen first nonidentity label required")
  pilot <- file.path(out, "allocation_00002")
  if (!dir.exists(pilot)) dir.create(pilot)
  final_path <- file.path(pilot, "PILOT_COMPLETE.json")
  jgwx_need(!file.exists(final_path), "Pilot already complete; read its retained checkpoint, do not refit")
  native_path <- file.path(pilot, "NATIVE_CAPSULE.rds")
  saved_path <- file.path(pilot, "FIT_SAVED.json")
  reused <- file.exists(saved_path)
  if (reused) {
    saved <- jsonlite::read_json(saved_path, simplifyVector = FALSE)
    jgwx_need(identical(saved$allocation_key, keys$four_group_MC[2L]) &&
      identical(saved$prepared_sha256, complete$prepared$sha256), "Saved fit identity differs")
    jgwx_assert_pin(native_path, saved$native)
    value <- readRDS(native_path)
    execution <- saved$execution
    execution$RNGkind <- unlist(execution$RNGkind, use.names = FALSE)
    execution$seed <- as.integer(execution$seed)
    execution$workers <- as.integer(execution$workers)
    execution$mroast_nrot <- as.integer(execution$mroast_nrot)
    execution$fgsea_eps <- as.numeric(execution$fgsea_eps)
  } else {
    jgwx_need(!file.exists(native_path) && !file.exists(paste0(native_path, ".partial")) &&
      !file.exists(file.path(pilot, "FIT_STARTED.json")), "Incomplete fit retained; no automatic refit")
    jgwx_json(list(status = "SINGLE_FIRST_NONIDENTITY_NATIVE_CALL_START",
      allocation_key = keys$four_group_MC[2L], prepared_sha256 = complete$prepared$sha256,
      started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(pilot, "FIT_STARTED.json"))
    do.call(RNGkind, as.list(context$runtime$RNGkind)); set.seed(context$runtime$seed)
    started <- proc.time()[["elapsed"]]
    value <- modules$kernel$jgw03_fit_one_allocation(fixed, allocation, modules$pins$native$path)
    fit_seconds <- proc.time()[["elapsed"]] - started
    saveRDS(value, paste0(native_path, ".partial"), compress = "gzip", version = 3L)
    jgwx_need(file.rename(paste0(native_path, ".partial"), native_path), "Native capsule rename failed")
    execution <- list(seed = context$runtime$seed, RNGkind = context$runtime$RNGkind,
      backend = context$runtime$backend, workers = 1L, mroast_nrot = 1999L, fgsea_eps = 0,
      mroast_applicable = accession == "GSE289929",
      mroast_executed = accession == "GSE289929",
      actual_mroast_nrot = if (accession == "GSE289929") 1999L else NULL,
      mroast_configuration_provenance = if (accession == "GSE289929")
        "Actual PC12 Reactome native mroast; configured nrot=1999" else
        "Unused source-config-only nrot=1999; no native mroast called in this cell",
      provenance = paste("Original kernel actual single first nonidentity; prepared", complete$prepared$sha256),
      RNG_after_sha256 = digest::digest(.Random.seed, "sha256"), historical_backend_authenticated = FALSE)
    jgwx_json(list(status = "FIT_RETAINED_COMPACT_PENDING", allocation_key = keys$four_group_MC[2L],
      prepared_sha256 = complete$prepared$sha256, native = jgwx_pin(native_path),
      fit_seconds = fit_seconds, execution = execution), saved_path)
  }
  jgwx_need(identical(value$accession, accession) && identical(value$labels, allocation), "Saved allocation identity differs")
  compact <- h$jgw03c_extract(value, domain, execution)
  # Full S4/closure identity is not a serialization requirement. Verify every
  # transported native numeric/dynamic vector through the plain compact view.
  loaded <- readRDS(native_path)
  from_saved <- h$jgw03c_extract(loaded, domain, execution)
  jgwx_need(identical(serialize(compact, NULL, version = 3L),
    serialize(from_saved, NULL, version = 3L)), "Native capsule compact-view transport differs")
  compact_path <- file.path(pilot, "COMPACT.rds")
  if (file.exists(compact_path)) jgwx_need(identical(serialize(readRDS(compact_path), NULL, version = 3L),
    serialize(compact, NULL, version = 3L)), "Existing compact differs; no overwrite") else
      jgwx_atomic_plain(compact, compact_path)
  for (record in context$observed_files) jgwx_assert_pin(record$path, record[c("bytes", "sha256")])
  jgwx_need(identical(context$protected,
    h$jgwp_bound_records(root, accession, modules$admission)) && identical(context$sources, jgwx_modules(root)$pins),
    "Protected source/input changed during pilot")
  jgwx_json(list(status = "FIRST_NONIDENTITY_SAVED_AND_COMPACTED_ONLY", accession = accession,
    allocation_id = 2L, allocation_key = compact$allocation_key, domain_sha256 = domain$sha256,
    prepared_sha256 = complete$prepared$sha256, native = jgwx_pin(native_path), compact = jgwx_pin(compact_path),
    native_families = length(compact$families), reused_saved_native_without_fit = reused,
    new_native_kernel_calls_this_invocation = as.integer(!reused), original_label_refitted = FALSE,
    full_10000_MC_and_20_companion_not_completed = TRUE,
    full_reference_accepted = FALSE), final_path)
  0L
}
jgwx_main <- function(args) {
  jgwx_need(length(args) %in% c(4L, 6L), "ROOT RUN_CELL ACCESSION prepare [COMPACT_RECEIPT SHA] | pilot required")
  root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
  out <- gsub("\\\\", "/", args[[2L]])
  accession <- args[[3L]]; mode <- args[[4L]]
  jgwx_need(accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    startsWith(out, "R:/Shaoyao_storage/cache/JGW03_label_reference_20261015r01/") &&
    startsWith(normalizePath(tempdir(), winslash = "/"), "R:/Shaoyao_storage/"), "Known cell and R output/temp required")
  if (mode == "prepare") {
    jgwx_need(length(args) == 6L && !file.exists(out) && dir.exists(dirname(out)), "Fresh preparation leaf required")
    return(jgwx_prepare(root, out, accession, args[[5L]], args[[6L]]))
  }
  jgwx_need(mode == "pilot" && length(args) == 4L && dir.exists(out), "Explicit existing first-pilot leaf required")
  jgwx_pilot(root, out, accession)
}
if (sys.nframe() == 0L) {
  options(stringsAsFactors = FALSE, width = 140)
  exit_code <- tryCatch(jgwx_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_FIRST_NONIDENTITY_STOP_RETAINED_NO_AUTOMATIC_REFIT")
    1L
  })
  quit(status = exit_code)
}


