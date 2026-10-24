# Dedicated input/observed/job preparation for the admitted GSE128623 reference.
# No old driver is sourced; no fit occurs when declarations are imported.
# CLI: prepare RUN_ROOT | observe RUN_ROOT CELL | finalize RUN_ROOT ADMISSION_JSON

g128p_need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
g128p_sha <- function(path) digest::digest(file = path, algo = "sha256")
g128p_pin <- function(path) list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
  bytes = unname(file.info(path)$size), sha256 = g128p_sha(path))
g128p_cells <- function() c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
g128p_root <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  g128p_need(identical(root, "/root/autodl-tmp/gse128623_label_reference_20261015r01"),
    "Only the explicit new run root is accepted")
  root
}
g128p_json <- function(value, path) {
  g128p_need(!file.exists(path), "Refuse to overwrite existing metadata")
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17,
    na = "null", null = "null")
}
g128p_bundle_tables <- function(bundle) c(bundle$gene_camera$tables[c("edgeR", "voom", "DESeq2")],
  list(activity = bundle$activity$native, CAMERA = bundle$gene_camera$tables$CAMERA,
    mroast = bundle$unsigned$mroast, fgsea = bundle$unsigned$fgsea))
g128p_plain_tables <- function(value) {
  # Ignore only data.table's process-local self-reference pointer in a content
  # check, retaining all data-frame columns/keys/values and nested list order.
  if (inherits(value, "data.table")) return(as.data.frame(value))
  if (is.list(value) && !is.data.frame(value)) return(lapply(value, g128p_plain_tables))
  value
}
g128p_save <- function(value, path, native_bundle = FALSE) {
  g128p_need(!file.exists(path) && !file.exists(paste0(path, ".partial")), "Output/partial already exists")
  temporary <- paste0(path, ".partial")
  saveRDS(value, temporary, version = 3L, compress = "gzip")
  restored <- readRDS(temporary)
  if (native_bundle) {
    # S4 fitted objects contain closure environments: pointer identity is not
    # a meaningful RDS transport criterion. Check the retained result fields,
    # dependent count/E matrices, diagnostics and RNG records exactly instead.
    is_job <- all(c("context", "inputs", "observed_bundle", "allocations") %in% names(value))
    before_bundle <- if (is_job) value$observed_bundle else value
    after_bundle <- if (is_job) restored$observed_bundle else restored
    if (is_job) g128p_need(identical(value[c("context", "inputs", "allocations")],
      restored[c("context", "inputs", "allocations")]), "Job context/input/label transport changed")
    first <- g128p_bundle_tables(before_bundle); second <- g128p_bundle_tables(after_bundle)
    for (name in names(first)) g128p_need(identical(lapply(as.data.frame(first[[name]]), identity),
      lapply(as.data.frame(second[[name]]), identity)), "Observed native table transport changed")
    g128p_need(identical(before_bundle$gene_camera$native$dge$counts, after_bundle$gene_camera$native$dge$counts) &&
      identical(before_bundle$gene_camera$native$voom$E, after_bundle$gene_camera$native$voom$E) &&
      identical(before_bundle$integration_diagnostics, after_bundle$integration_diagnostics) &&
      identical(g128p_plain_tables(before_bundle$unsigned$prepared), g128p_plain_tables(after_bundle$unsigned$prepared)) &&
      identical(before_bundle$unsigned$rng_calls, after_bundle$unsigned$rng_calls) &&
      identical(before_bundle$preparation_evidence, after_bundle$preparation_evidence), "Observed dependent/RNG transport changed")
  } else g128p_need(identical(value, restored), "RDS roundtrip changed payload")
  g128p_need(file.link(temporary, path), "Exclusive atomic publication failed; retain partial")
  g128p_need(unlink(temporary) == 0L, "Own successful partial link cleanup failed")
  g128p_pin(path)
}
g128p_read_table <- function(path, ...) {
  # Base R decompresses the existing small TSV exports without adding R.utils.
  # fread still applies the same field types, NA handling and selection options.
  if (grepl("[.]gz$", path)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection))
    text <- paste(readLines(connection, warn = FALSE), collapse = "\n")
    return(as.data.frame(data.table::fread(text = text, na.strings = "NA",
      check.names = FALSE, ...), stringsAsFactors = FALSE))
  }
  as.data.frame(data.table::fread(path, na.strings = "NA", check.names = FALSE,
    ...), stringsAsFactors = FALSE)
}
g128p_check_record <- function(path, expected) {
  actual <- g128p_pin(path)
  g128p_need(identical(actual$sha256, expected$sha256) && actual$bytes == expected$bytes,
    paste("Retained identity mismatch:", basename(path)))
  actual
}
g128p_runtime <- function() {
  expected <- c(edgeR = "4.10.1", limma = "3.68.4", DESeq2 = "1.52.0", fgsea = "1.38.0",
    data.table = "1.18.4", openssl = "2.4.2", BiocParallel = "1.46.0", digest = "0.6.39", jsonlite = "2.0.0")
  g128p_need(identical(as.character(getRversion()), "4.6.1"), "R version changed")
  actual <- vapply(names(expected), function(p) as.character(utils::packageVersion(p)), character(1))
  g128p_need(identical(actual, expected), "Existing package version differs")
  g128p_need(all(vapply(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"),
    function(key) identical(Sys.getenv(key), "1"), logical(1))), "Single-thread launcher environment required")
  data.table::setDTthreads(1L)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  backend <- BiocParallel::bpparam()
  g128p_need(inherits(backend, "SerialParam") && is.null(BiocParallel::bpRNGseed(backend)), "Native backend differs")
  list(R = R.version.string, packages = as.list(actual), libraries = .libPaths(),
    backend_class = class(backend), backend_workers = BiocParallel::bpworkers(backend),
    backend_RNGseed = BiocParallel::bpRNGseed(backend), BLAS = extSoftVersion()["BLAS"],
    fgsea_formals = deparse(formals(fgsea::fgseaMultilevel)),
    mroast_formals = deparse(formals(limma::mroast)), command = commandArgs())
}
g128p_bind <- function(root) {
  directory <- file.path(root, "source")
  worker <- file.path(directory, "GSE128623_reference_shard_core_20261015.R")
  entry <- file.path(directory, "run_GSE128623_reference_shard_20261015.R")
  g128p_need(identical(g128p_sha(worker), "6f595b13ebaa00046ff2b204e810a819d7d3ae8639416d418d04a8e5373e5896") &&
    identical(g128p_sha(entry), "21a850894037b5e7bcb487909f15200f007ebecbdf149200b317b5670166c88f"), "Worker/entry pin changed")
  environment <- new.env(parent = baseenv())
  sys.source(worker, environment)
  paths <- setNames(file.path(directory, c("GSE128623_gene_camera_label_core_20261015.R",
    "GSE128623_activity_label_core_20261015.R", "GSE128623_unsigned_label_core_20261015.R",
    "GSE128623_label_bundle_core_v2_20261015.R", "GSE128623_reference_archive_core_20261015.R")),
    names(environment$g128w_source_pins()))
  binding <- environment$g128w_bind_sources(paths)
  identity <- c(binding$pins, worker = g128p_sha(worker), entry = g128p_sha(entry))
  list(binding = binding, worker_environment = environment, identity = identity,
    source_manifest_sha256 = digest::digest(identity, algo = "sha256"),
    paths = c(paths, worker = worker, entry = entry))
}
g128p_old_roots <- function() list(
  native = "/root/autodl-tmp/gse128623_real_run_20261007r04",
  scores = "/root/autodl-tmp/route_b2_gse128_program_baseline_20261014r01")
g128p_verify_manifest <- function(manifest) {
  g128p_need(!anyDuplicated(manifest$role) && !anyDuplicated(manifest$path), "Manifest duplicate role/path")
  for (i in seq_len(nrow(manifest))) g128p_check_record(manifest$path[[i]],
    list(bytes = manifest$bytes[[i]], sha256 = manifest$sha256[[i]]))
  invisible(TRUE)
}
g128p_load <- function(root) {
  path <- file.path(root, "PREPARED.rds")
  prepared <- readRDS(path)
  g128p_need(identical(prepared$status, "INPUTS_FROZEN_NO_MODELS") &&
    identical(prepared$run_root, root), "Wrong preparation object")
  g128p_need(identical(g128p_sha(file.path(root, "INPUT_MANIFEST.tsv")), prepared$input_manifest_sha256),
    "Input manifest drift")
  g128p_verify_manifest(prepared$manifest)
  prepared
}
g128p_native_paths <- function(model, cohort) setNames(file.path(model, "results", cohort, c(
  "gene_level/edgeR_primary.tsv.gz", "gene_level/limma_voom_sensitivity.tsv.gz",
  "gene_level/DESeq2_Wald_sensitivity.tsv.gz", "programs/activity_all_contrasts.tsv.gz",
  "programs/camera.tsv.gz", "programs/mroast_chronic.tsv.gz", "programs/fgsea_chronic.tsv.gz")),
  c("edgeR", "voom", "DESeq2", "activity", "CAMERA", "mroast", "fgsea"))
g128p_prepare <- function(root) {
  g128p_need(!file.exists(file.path(root, "PREPARED.rds")) &&
    !file.exists(file.path(root, "PREPARATION_STARTED.json")), "Prior preparation exists; preserve/inspect, do not overwrite")
  runtime <- g128p_runtime()
  bound <- g128p_bind(root)
  source <- file.path(root, "source")
  frozen <- c(original = "analyze_gse128623_exact30_linux_model_v2.R",
    proposal = "GSE128623_FULL_LABEL_REFERENCE_EXECUTION_PROPOSAL_20261015.md",
    allocation_definition = "GSE128623_LABEL_ALLOCATION_REGISTRY_20261015.md",
    contrasts = "contrasts.yaml", sap = "statistical_analysis_plan.md",
    native_config = "GSE128623_exact30_linux_model_v2.json")
  expected <- c(original = "e251c6ae9bb70cd944ebb1d781895878ee80963946f6e8ad94a52a83713add50",
    proposal = "15fec7caee23ac4eeb12af6a7363c5590aa14a001300a8686aa203be3b5096b4",
    allocation_definition = "bffcf90f2339f42f3138feb82ac564c9370db15209aea8de87289f726b161d58",
    contrasts = "98ffd869fb3703e38e3318f18887ad746d78a54c59e5a36e5181648f02b63492",
    sap = "c46c086cbba963436fc48cc88f77028713cb305f9e9bff327314fe6eb24e0b82",
    native_config = "dbff1a9799d6827314dbf8ac9e47ccbadeba13725edadf73bc8f0f72e07017af")
  frozen_paths <- setNames(file.path(root, "protocol", frozen), names(frozen))
  frozen_paths[["original"]] <- file.path(source, frozen[["original"]])
  g128p_need(identical(vapply(frozen_paths, g128p_sha, character(1)), expected), "Frozen definition/config/native source drift")
  label_root <- file.path(root, "registry")
  registry_path <- file.path(label_root, "PREPARED.json")
  g128p_need(identical(g128p_sha(registry_path), "63be688e2d548b2917479a9ecc6e1dba4b3ab37542ac8c4aca9244e78c3ef5d0"), "Registry receipt changed")
  registry <- jsonlite::read_json(registry_path, simplifyVector = FALSE)
  roots <- g128p_old_roots()
  recovery_path <- file.path(roots$scores, "RECOVERY_COMPLETE.json")
  recovery_manifest <- file.path(roots$scores, "INPUT_MANIFEST.tsv")
  g128p_need(identical(g128p_sha(recovery_path), "f88be11424933cbb13d6864dcf6a83e6d1ec8b018c1a762695f85f67433f0844") &&
    identical(g128p_sha(recovery_manifest), "537d1b5a30fea552f80f4d674ac6812253333c3740565ef480fe71444a1a6e36"), "Accepted recovery pins changed")
  recovery <- jsonlite::read_json(recovery_path, simplifyVector = FALSE)
  old_manifest <- g128p_read_table(recovery_manifest)
  g128p_verify_manifest(old_manifest)
  common <- setNames(file.path(roots$native, "analysis_inputs", c("sample_identity_design_ledger.tsv",
    "feature_annotation_from_existing_counts.tsv", "gene_set_resources_frozen_2026-08-18.rds", "S1_GATE_frozen.tsv")),
    c("ledger", "annotation", "resource", "s1_zero"))
  paths <- c(bound$paths, frozen_paths, preparer = file.path(source, "prepare_GSE128623_reference_observed_20261015.R"),
    registry = registry_path, recovery = recovery_path, recovery_manifest = recovery_manifest, common)
  cells <- list()
  for (processing in c("raw", "adapter")) {
    model <- file.path(roots$native, if (processing == "raw") "model_raw_copy_isolated_verified_20261013" else "model_adapter_copy_isolated_20261013")
    completion_path <- file.path(model, "MODEL_COMPLETE.json")
    completion <- jsonlite::read_json(completion_path, simplifyVector = FALSE)
    expected_sha <- if (processing == "raw") "5f92e23150275f678794597fe3f0a6c4fff361f0c13e975df27c1edd0b40b497" else "5ce97f6243b1adfe6e62fbbb62c7955d815810c7cfd78df8a6da32643a296548"
    g128p_need(identical(g128p_sha(completion_path), expected_sha), "Native completion changed")
    paths[paste0(processing, "_complete")] <- completion_path
    for (role in names(common)) g128p_check_record(common[[role]], completion$input_bindings[[role]])
    for (cohort in c("all30", "predefined29")) {
      old_cohort <- if (cohort == "all30") "primary_all30" else "complete29_sensitivity"
      count_stem <- if (cohort == "all30") "all30" else "complete29"
      id <- paste(processing, cohort, sep = "_")
      count_path <- file.path(roots$native, "derived/assembled", processing, paste0(count_stem, "_selected.tsv"))
      g128p_check_record(count_path, completion$input_bindings[[paste0(processing, "/", count_stem, "_selected.tsv")]])
      native_paths <- g128p_native_paths(model, old_cohort)
      for (path in native_paths) g128p_check_record(path, completion$outputs[[substring(path, nchar(file.path(model, "results")) + 2L)]])
      score <- recovery$cells[[paste(processing, old_cohort, sep = "_")]]$score
      g128p_check_record(score$path, score)
      library_path <- file.path(model, "results", old_cohort, "diagnostics/library_sizes_TMM.tsv")
      g128p_check_record(library_path, completion$outputs[[paste0(old_cohort, "/diagnostics/library_sizes_TMM.tsv")]])
      cells[[id]] <- list(processing = processing, cohort = cohort, score = score$path,
        counts = count_path, native = native_paths, library = library_path)
      paths <- c(paths, setNames(c(count_path, score$path, library_path, native_paths),
        paste0(id, "/", c("counts", "scores", "library", names(native_paths)))))
    }
  }
  for (cohort in c("all30", "predefined29")) {
    item <- registry$cohorts[[cohort]]
    spaces <- c(list(six_group = item$sampled_space), setNames(item$exact_day_control_spaces,
      paste0(names(item$exact_day_control_spaces), "_vs_control_exact")))
    for (space in names(spaces)) {
      path <- file.path(label_root, spaces[[space]]$file)
      g128p_need(identical(g128p_sha(path), spaces[[space]]$sha256), "Frozen allocation file changed")
      paths[paste0(cohort, "/", space)] <- path
    }
  }
  g128p_need(!anyDuplicated(paths) && !anyDuplicated(names(paths)), "Manifest roles/paths not unique")
  manifest <- do.call(rbind, lapply(names(paths), function(role) {
    record <- g128p_pin(paths[[role]])
    data.frame(role = role, path = record$path, bytes = record$bytes, sha256 = record$sha256)
  }))
  g128p_json(list(mode = "prepare", command = commandArgs(), pid = Sys.getpid(),
    runtime = runtime, models_run = FALSE), file.path(root, "PREPARATION_STARTED.json"))
  for (directory in c("inputs", "observed", "jobs")) g128p_need(dir.create(file.path(root, directory)), "Prepared output directory already exists")
  # Only the three named original mapping declarations are evaluated; no driver.
  original <- new.env(parent = baseenv())
  for (node in parse(frozen_paths[["original"]])) if (is.call(node) && identical(node[[1L]], as.name("<-")) &&
    is.symbol(node[[2L]]) && as.character(node[[2L]]) %in% c("v2a_stop", "v2a_assert", "v2a_program_members")) eval(node, original)
  members <- original$v2a_program_members(readRDS(common[["resource"]]))
  annotation_raw <- g128p_read_table(common[["annotation"]], colClasses = "character")
  annotation <- data.frame(gene_id = trimws(annotation_raw$GeneID), gene_symbol = trimws(annotation_raw$gene))
  ledger <- g128p_read_table(common[["ledger"]])
  for (id in names(cells)) {
    cell <- cells[[id]]
    recovered <- readRDS(cell$score)
    metadata <- ledger[if (cell$cohort == "all30") rep(TRUE, nrow(ledger)) else ledger$sample_id != "GSM3681230", , drop = FALSE]
    g128p_need(identical(as.character(metadata$sample_id), as.character(recovered$metadata$sample_id)) &&
      identical(as.character(metadata$group), as.character(recovered$metadata$group)), "Recovered metadata differs from fixed ledger")
    table <- g128p_read_table(cell$counts)
    counts <- as.matrix(table[, -1L, drop = FALSE]); rownames(counts) <- table$gene_id
    g128p_need(identical(names(table), c("gene_id", metadata$sample_id)) &&
      identical(annotation$gene_id, rownames(counts)), "Full count/annotation/sample order")
    genes <- g128p_read_table(cell$native[["edgeR"]], select = c("gene_id", "contrast"))
    gene_ids <- genes$gene_id[genes$contrast == "chronic_combined"]
    g128p_need(!anyDuplicated(gene_ids) && identical(rownames(counts)[rownames(counts) %in% gene_ids], gene_ids), "Original retained gene order")
    for (contrast in unique(genes$contrast)) g128p_need(identical(genes$gene_id[genes$contrast == contrast], gene_ids), "Original gene domain drift")
    inputs <- list(counts = counts, metadata = metadata, annotation = annotation,
      members = as.data.frame(members), scores = recovered$scores, registry = as.data.frame(recovered$registry))
    labels <- setNames(as.character(metadata$group), metadata$sample_id)
    gene_environment <- environment(bound$binding$kernels$gene_camera)
    gene_environment$g128lc_validate_inputs(counts, metadata, annotation, labels)
    g128p_need(identical(as.character(unlist(registry$cohorts[[cell$cohort]]$sample_ids)), metadata$sample_id) &&
      identical(colnames(inputs$scores), metadata$sample_id) && all(is.finite(inputs$scores)) &&
      identical(rownames(inputs$scores), paste(inputs$registry$collection, inputs$registry$program_id, sep = "::")) &&
      !anyDuplicated(rownames(inputs$scores)) && all(abs(rowMeans(inputs$scores)) <= 1e-8) &&
      all(abs(apply(inputs$scores, 1L, stats::sd) - 1) <= 1e-8), "Frozen score/registry/sample preprocessing identity")
    cells[[id]]$expected_gene_ids <- gene_ids
    cells[[id]]$inputs_file <- file.path(root, "inputs", paste0(id, ".rds"))
    cells[[id]]$inputs_pin <- g128p_save(inputs, cells[[id]]$inputs_file)
    cells[[id]]$inputs_object_sha256 <- digest::digest(inputs, algo = "sha256")
  }
  data.table::fwrite(manifest, file.path(root, "INPUT_MANIFEST.tsv"), sep = "\t", quote = FALSE)
  prepared <- list(status = "INPUTS_FROZEN_NO_MODELS", run_root = root, manifest = manifest,
    input_manifest_sha256 = g128p_sha(file.path(root, "INPUT_MANIFEST.tsv")),
    source_manifest_sha256 = bound$source_manifest_sha256, source_identity = bound$identity,
    registry = registry, cells = cells, runtime = runtime)
  g128p_verify_manifest(manifest)
  g128p_save(prepared, file.path(root, "PREPARED.rds"))
  g128p_json(list(status = prepared$status, input_manifest_sha256 = prepared$input_manifest_sha256,
    source_manifest_sha256 = prepared$source_manifest_sha256, cells = g128p_cells(),
    models_run = FALSE, full_run_admitted = FALSE), file.path(root, "INPUTS_PREPARED.json"))
  invisible(prepared)
}

g128p_component <- function(path, binding, compute) {
  if (file.exists(path)) {
    capsule <- readRDS(path)
    g128p_need(identical(capsule$binding, binding) &&
      is.raw(capsule$serialized_value) && identical(digest::digest(capsule$serialized_value,
        algo = "sha256", serialize = FALSE), capsule$serialized_value_sha256), "Completed component binding/digest mismatch")
    return(unserialize(capsule$serialized_value))
  }
  lease <- paste0(path, ".lease")
  g128p_need(dir.create(lease, showWarnings = FALSE), "Incomplete/owned component exists; inspect, never automatically refit")
  g128p_json(list(binding = binding, pid = Sys.getpid(), command = commandArgs()), file.path(lease, "STARTED.json"))
  start <- proc.time()[["elapsed"]]
  value <- compute()
  serialized_value <- serialize(value, NULL, version = 3L)
  capsule <- list(binding = binding, serialized_value = serialized_value,
    serialized_value_sha256 = digest::digest(serialized_value, algo = "sha256", serialize = FALSE),
    elapsed_seconds = proc.time()[["elapsed"]] - start, process_id = Sys.getpid(),
    process_high_water_RSS_bytes = g128p_peak())
  g128p_save(capsule, path)
  value
}
g128p_compare <- function(actual, expected, method) {
  actual <- as.data.frame(actual); expected <- as.data.frame(expected)
  g128p_need(identical(names(actual), names(expected)) && nrow(actual) == nrow(expected), paste("Native schema/row count:", method))
  keys <- if (method %in% c("edgeR", "voom", "DESeq2")) c("gene_id", "contrast") else c("collection", "program_id", "contrast")
  g128p_need(!anyDuplicated(actual[keys]) && identical(lapply(actual[keys], as.character), lapply(expected[keys], as.character)), paste("Native complete ordered keys:", method))
  rows <- lapply(names(actual), function(field) {
    a <- actual[[field]]; b <- expected[[field]]
    missing_equal <- identical(is.na(a), is.na(b))
    maximum <- 0; scaled <- 0; log_delta <- 0
    if (is.numeric(a)) {
      g128p_need(is.numeric(b) || all(is.na(b)), paste("Native numeric type:", method, field))
      finite <- is.finite(a) & is.finite(b)
      infinity_equal <- identical(which(is.infinite(a)), which(is.infinite(b))) &&
        identical(sign(a[is.infinite(a)]), sign(b[is.infinite(b)]))
      difference <- abs(a[finite] - b[finite]); maximum <- if (length(difference)) max(difference) else 0
      scaled <- if (length(difference)) max(difference / (1 + abs(b[finite]))) else 0
      matches <- missing_equal && infinity_equal && all(difference <= 1e-10 + 1e-8 * abs(b[finite]))
      if (field %in% c("PValue", "BH_q", "FDR", "PValue.Mixed", "FDR.Mixed", "pval", "padj")) {
        positive <- finite & a > 0 & b > 0
        log_delta <- if (any(positive)) max(abs(log10(a[positive]) - log10(b[positive]))) else 0
        matches <- matches && identical(which(!is.na(a) & a == 0), which(!is.na(b) & b == 0)) &&
          log_delta <= 1e-6 && identical(a[finite] < .05, b[finite] < .05) &&
          identical(a[finite] <= .05, b[finite] <= .05)
      }
      if (field %in% c("ci_low", "ci_high")) matches <- matches &&
        identical(a[finite] > 0, b[finite] > 0) && identical(a[finite] < 0, b[finite] < 0)
      if (field %in% c("logFC", "effect", "t", "WaldStatistic", "ES", "NES", "QLF_F")) matches <- matches &&
        identical(a[finite] == 0, b[finite] == 0) && identical(sign(a[finite]), sign(b[finite]))
    } else matches <- missing_equal && identical(as.character(a), as.character(b))
    data.frame(method = method, field = field, rows = length(a), missing_pattern_equal = missing_equal,
      max_absolute_difference = maximum, max_scaled_difference = scaled, max_log10_probability_difference = log_delta,
      compatible = matches, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
g128p_peak <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", line)) * 1024
}
g128p_observe <- function(root, id) {
  g128p_need(id %in% g128p_cells(), "Unknown fixed cell")
  runtime <- g128p_runtime(); prepared <- g128p_load(root); bound <- g128p_bind(root)
  cell <- prepared$cells[[id]]
  g128p_check_record(cell$inputs_file, cell$inputs_pin)
  inputs <- readRDS(cell$inputs_file)
  g128p_need(identical(digest::digest(inputs, algo = "sha256"), cell$inputs_object_sha256), "Input object changed")
  directory <- file.path(root, "observed", id)
  if (!dir.exists(directory)) g128p_need(dir.create(directory), "Cannot create observed leaf")
  binding <- list(cell_id = id, input_manifest_sha256 = prepared$input_manifest_sha256,
    source_manifest_sha256 = prepared$source_manifest_sha256, inputs_object_sha256 = cell$inputs_object_sha256,
    preparer_sha256 = g128p_sha(file.path(root, "source", "prepare_GSE128623_reference_observed_20261015.R")))
  labels <- setNames(as.character(inputs$metadata$group), inputs$metadata$sample_id)
  gene <- g128p_component(file.path(directory, "gene_camera.rds"), c(binding, component = "gene_camera"),
    function() do.call(bound$binding$kernels$gene_camera, c(inputs[c("counts", "metadata", "annotation", "members")], list(labels = labels))))
  g128p_need(identical(gene$retained_gene_ids, cell$expected_gene_ids), "Observed native mask changed")
  library <- g128p_read_table(cell$library)
  g128p_need(identical(library$sample_id, inputs$metadata$sample_id), "Historical TMM sample order")
  kept <- gene$native$dge$counts
  old_E <- t(log2(t(kept + .5) / (library$effective_library_size + 1) * 1e6))
  g128p_need(identical(dimnames(old_E), dimnames(gene$native$voom$E)) &&
    all(abs(old_E - gene$native$voom$E) <= 1e-12 * (1 + abs(old_E))), "Fixed-score preprocessing/TMM-E binding changed")
  wrappers <- list(gene_camera = function(...) gene,
    activity = function(...) { arguments <- list(...); g128p_component(file.path(directory, "activity.rds"),
      c(binding, component = "activity"), function() do.call(bound$binding$kernels$activity, arguments)) },
    unsigned = function(...) { arguments <- list(...); g128p_component(file.path(directory, "unsigned.rds"),
      c(binding, component = "unsigned"), function() do.call(bound$binding$kernels$unsigned, arguments)) })
  bundle_path <- file.path(directory, "OBSERVED_BUNDLE.rds")
  if (file.exists(bundle_path)) {
    bundle <- readRDS(bundle_path)
    g128p_need(identical(bundle$preparation_evidence, binding), "Existing observed bundle binding changed")
  } else {
    bundle <- do.call(bound$binding$fit, c(inputs, list(expected_retained_gene_ids = cell$expected_gene_ids,
      baseline_voom_E = gene$native$voom$E, labels = labels, cohort = cell$cohort, kernels = wrappers)))
    bundle$preparation_evidence <- binding
    g128p_save(bundle, bundle_path, native_bundle = TRUE)
    bundle <- readRDS(bundle_path)
  }
  packed <- bound$binding$pack(bundle)
  bound$binding$validate(packed)
  g128p_need(nrow(packed$families) == 67L, "Complete observed native family template differs")
  frames <- g128p_bundle_tables(bundle)
  comparisons <- do.call(rbind, lapply(names(frames), function(method) g128p_compare(frames[[method]],
    g128p_read_table(cell$native[[method]]), method)))
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20260819L)
  expected_initial <- get(".Random.seed", .GlobalEnv)
  do.call(RNGkind, as.list(old_kind))
  if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  calls <- bundle$unsigned$rng_calls
  g128p_need(length(calls) == 4L && all(vapply(calls, function(call) {
    identical(call$initial_state, expected_initial) && identical(call$seed, 20260819L) &&
      identical(call$rng_kind, c("Mersenne-Twister", "Inversion", "Rejection")) &&
      identical(call$final_rng_kind, call$rng_kind) && is.integer(call$final_state) &&
      length(call$final_state) == length(expected_initial) && "SerialParam" %in% call$backend$class &&
      call$backend$workers == 1L && is.null(call$backend$RNGseed)
  }, logical(1))), "Actual native inner RNG/backend trace differs")
  bridge_path <- file.path(directory, "OBSERVED_BRIDGE.rds")
  component_measurements <- lapply(c("gene_camera", "activity", "unsigned"), function(name) {
    saved <- readRDS(file.path(directory, paste0(name, ".rds")))
    saved[c("elapsed_seconds", "process_id", "process_high_water_RSS_bytes")]
  })
  names(component_measurements) <- c("gene_camera", "activity", "unsigned")
  measured_peaks <- c(g128p_peak(), vapply(component_measurements, function(record) record$process_high_water_RSS_bytes, numeric(1)))
  bridge <- list(binding = binding, comparisons = comparisons, compatible = all(comparisons$compatible),
    observed_bundle = g128p_pin(bundle_path), packed_rows = nrow(packed$keys), families = 67L,
    runtime = runtime, component_measurements = component_measurements,
    process_high_water_RSS_bytes = if (any(is.finite(measured_peaks))) max(measured_peaks, na.rm = TRUE) else NA_real_,
    RSS_definition = "maximum saved component-producing process VmHWM and current readback VmHWM; never reconstructed by refitting",
    comparison_layer = "saved full-precision observed RDS native fields versus original native TSV exports",
    RNG_calls = calls,
    source_score_preprocessing_maximum_absolute_difference = max(abs(old_E - gene$native$voom$E)),
    native_bound = "abs<=1e-10+1e-8*abs(expected); positive P/q log10<=1e-6; P/q zero and .05 boundaries, signed effect/statistic zero+sign, NA/text/CI-sign exact",
    actual_process_exit_confirmed_here = FALSE, full_run_admitted = FALSE)
  if (file.exists(bridge_path)) {
    saved <- readRDS(bridge_path)
    g128p_need(identical(saved$binding, binding) && identical(saved$comparisons, comparisons) &&
      identical(saved$observed_bundle, bridge$observed_bundle), "Existing observed bridge differs")
    bridge <- saved
  } else g128p_save(bridge, bridge_path)
  g128p_verify_manifest(prepared$manifest)
  g128p_need(identical(digest::digest(inputs, algo = "sha256"), cell$inputs_object_sha256), "Observed input object mutated")
  g128p_need(bridge$compatible, "Observed field comparison failed; completed fits retained; do not refit")
  receipt <- file.path(directory, "OBSERVED_COMPLETE.json")
  if (!file.exists(receipt)) g128p_json(list(status = "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED",
    cell_id = id, binding = binding, observed_bundle = g128p_pin(bundle_path), bridge = g128p_pin(bridge_path),
    family_n = 67L, field_n = nrow(comparisons), peak_memory_bytes = bridge$process_high_water_RSS_bytes,
    actual_process_exit_confirmed_here = FALSE, full_run_admitted = FALSE), receipt)
  invisible(bridge)
}

g128p_finalize <- function(root, admission_path) {
  prepared <- g128p_load(root); bound <- g128p_bind(root)
  admission_path <- normalizePath(admission_path, winslash = "/", mustWork = TRUE)
  g128p_need(startsWith(admission_path, paste0(root, "/")), "Admission must stay in the exact run")
  admission <- jsonlite::read_json(admission_path, simplifyVector = FALSE)
  g128p_need(identical(admission$status, "ALL_FOUR_OBSERVED_AND_RESOURCE_PASS") &&
    identical(admission$run_id, basename(root)) && identical(unlist(admission$cell_ids), g128p_cells()) &&
    identical(admission$input_manifest_sha256, prepared$input_manifest_sha256) &&
    identical(admission$source_manifest_sha256, prepared$source_manifest_sha256) &&
    isTRUE(admission$resource_acceptance) && identical(names(admission$observed_process_exit_codes), g128p_cells()) &&
    all(unlist(admission$observed_process_exit_codes) == 0L), "External all-four observed/process/resource admission required; never prefilled here")
  g128p_check_record(admission$resource_evidence$path, admission$resource_evidence)
  jobs <- list()
  for (id in g128p_cells()) {
    cell <- prepared$cells[[id]]
    receipt_path <- file.path(root, "observed", id, "OBSERVED_COMPLETE.json")
    g128p_check_record(receipt_path, admission$observed_receipts[[id]])
    receipt <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
    g128p_need(identical(receipt$status, "OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED"), "Observed bridge not complete")
    g128p_check_record(receipt$observed_bundle$path, receipt$observed_bundle)
    g128p_check_record(receipt$bridge$path, receipt$bridge)
    observed <- readRDS(receipt$observed_bundle$path)
    inputs <- readRDS(cell$inputs_file)
    item <- prepared$registry$cohorts[[cell$cohort]]
    spaces <- c(list(six_group = item$sampled_space), setNames(item$exact_day_control_spaces,
      paste0(names(item$exact_day_control_spaces), "_vs_control_exact")))
    for (space in names(spaces)) {
      table <- g128p_read_table(file.path(root, "registry", spaces[[space]]$file), colClasses = "character")
      samples <- as.character(inputs$metadata$sample_id)
      allocations <- as.matrix(table[, samples, drop = FALSE])
      rownames(allocations) <- table$allocation_id
      target <- if (space == "six_group") NULL else sub("_exact$", "", space)
      bound$worker_environment$g128w_validate_allocations(allocations, inputs$metadata, cell$cohort, target, FALSE)
      context <- list(run_id = basename(root), run_root = root, cell_id = id, processing = cell$processing,
        cohort = cell$cohort, reference_space = space, worker_id = "not_launched",
        input_manifest_sha256 = prepared$input_manifest_sha256, source_manifest_sha256 = prepared$source_manifest_sha256,
        allocation_matrix_sha256 = digest::digest(allocations, algo = "sha256"),
        inputs_object_sha256 = digest::digest(inputs, algo = "sha256"),
        observed_bundle_object_sha256 = digest::digest(observed, algo = "sha256"),
        observed_evidence = list(file = receipt$observed_bundle$path, sha256 = receipt$observed_bundle$sha256),
        worker_core_sha256 = unname(bound$identity[["worker"]]), entry_source_sha256 = unname(bound$identity[["entry"]]),
        admission = admission, admission_file = g128p_pin(admission_path))
      path <- file.path(root, "jobs", paste0(id, "_", space, ".rds"))
      jobs[[paste(id, space, sep = "/")]] <- g128p_save(list(context = context, inputs = inputs,
        observed_bundle = observed, allocations = allocations), path, native_bundle = TRUE)
    }
  }
  g128p_need(length(jobs) == 24L, "Full four-cell/six-space job coverage")
  g128p_verify_manifest(prepared$manifest)
  g128p_json(list(status = "JOBS_PREPARED_NOT_LAUNCHED", jobs = jobs,
    admission = g128p_pin(admission_path), worker_started = FALSE, whole_run_complete = FALSE), file.path(root, "JOBS_PREPARED.json"))
  invisible(jobs)
}

g128p_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  g128p_need(length(arguments) %in% 2:3, "prepare RUN_ROOT | observe RUN_ROOT CELL | finalize RUN_ROOT ADMISSION_JSON")
  root <- g128p_root(arguments[[2L]])
  invoked <- grep("^--file=", commandArgs(), value = TRUE)
  g128p_need(length(invoked) == 1L && identical(normalizePath(sub("^--file=", "", invoked),
    winslash = "/", mustWork = TRUE), file.path(root, "source", "prepare_GSE128623_reference_observed_20261015.R")),
    "Invoke the staged pinned preparer itself")
  if (identical(arguments[[1L]], "prepare") && length(arguments) == 2L) g128p_prepare(root) else
  if (identical(arguments[[1L]], "observe") && length(arguments) == 3L) g128p_observe(root, arguments[[3L]]) else
  if (identical(arguments[[1L]], "finalize") && length(arguments) == 3L) g128p_finalize(root, arguments[[3L]]) else
    stop("Unknown preparation mode", call. = FALSE)
  cat("Requested preparation stage returned; no shard worker launched.\n")
}
if (sys.nframe() == 0L) g128p_main()
