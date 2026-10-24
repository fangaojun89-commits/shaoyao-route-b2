# Fixed prospective scientific worker V4 absent-target writer repair; V3 serial backend retained. Loading defines functions only. Native
# parent owns authentication, runtime/resource guards, child retirement and the
# final COMPLETED marker; this source confers no permission or acceptance.
gse205906_ens101_worker_need_v1 <- function(ok) {
  if (!isTRUE(ok)) stop("ENS101_CHECKED_MODEL_WORKER_UNVERIFIED", call. = FALSE)
}

gse205906_ens101_worker_pin_v1 <- function(raw) {
  gse205906_ens101_worker_need_v1(is.raw(raw))
  list(bytes = length(raw), sha256 = digest::digest(raw, "sha256", serialize = FALSE))
}

gse205906_ens101_worker_inputs_v1 <- function() {
  p <- function(n, h, relative) list(bytes = n, sha256 = h, relative_path = relative)
  list(
    matrix = p(1475133L, "e98296fdf308dbaa82325a3ca6b0b041c6f2569e98be75a45757141b0a1f77a2", "output/GSE205906_fixed9_gene_fragments_s0.tsv"),
    GTF = p(280052998L, "84fcb6692d62e4d8c75617461ae5f471fe89e10b6982c82f349393b1a80045df", "references/derived/Rattus_norvegicus.Rnor_6.0.101.gtf"),
    overlay = p(10054L, "b0df6fa1354a83237bcfcd8bf8c323517deecd0e5680da8202489c4464936444", "derived/NINE_LIBRARY_ANALYSIS_PREPARATION_OVERLAY_154.tsv"),
    manifest = p(893848L, "751e43b3aba8c2245981dc0041fec273c455c6390e85ef7ec4514e56118bf3b0", "provenance/run_manifest.json"),
    Reactome = p(11599354L, "1ea55f41f5b0a4ea51f723fb82b940f8b13fce1cab13ea6825b91d207c72fdc0", "data/processed/gene_sets/reactome97_human_long.tsv"),
    SynGO_sensitivity = p(4835295L, "2a40eeb2bddcf556ea5ac8d6ef6a0cb1bd71f85ac8dddb1aad8165a8d075eb3a", "data/processed/gene_sets/quickgo_syngo_mammalian_direct_human.tsv"),
    PROGENy = p(609987L, "b78892ec4d56212b08eb33ea80551fcb018e89b3c698e3d804008537c96bbdb9", "data/processed/gene_sets/progeny_human_top500.tsv"),
    DoRothEA_AB = p(687990L, "e709ea4a9e0a1517cb00082bfbf6e8c4408d786694b6864f31ec3f57b458ecb8", "data/processed/gene_sets/dorothea_rat_AB_one_to_one.tsv"),
    orthology = p(1922349L, "e9dd0c9268e0044f615085e27292cbf977f5a648003ba90054cf3d01225cf998", "data/processed/gene_sets/ncbi_human_mouse_rat_one_to_one.tsv"),
    frozen_resource_RDS = p(1506178L, "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51", "data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds"),
    gene_set_manifest = p(32072L, "b3a9089a6e136c5271cf6ac3999d10304dc2316fafd10f5c170d9f43084bbce1", "metadata/gene_set_manifest.tsv"))
}

gse205906_ens101_worker_import_v1 <- function(raw, wanted = NULL) {
  need <- gse205906_ens101_worker_need_v1
  declarations <- as.list(parse(text = rawToChar(raw), keep.source = FALSE))
  if (!is.null(wanted)) declarations <- Filter(function(x) is.call(x) &&
    identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
    as.character(x[[2L]]) %in% wanted, declarations)
  e <- new.env(parent = baseenv())
  e$head <- utils::head; e$tail <- utils::tail; e$setNames <- stats::setNames
  found <- character()
  for (x in declarations) {
    need(is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
      is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")))
    found <- c(found, as.character(x[[2L]])); eval(x, e)
  }
  need(!anyDuplicated(found) && (is.null(wanted) || setequal(found, wanted)))
  lockEnvironment(e, bindings = TRUE); e
}

gse205906_ens101_worker_sources_v1 <- function(bundle, full_raw, baseline_raw) {
  need <- gse205906_ens101_worker_need_v1
  pin <- function(n, h) list(bytes = n, sha256 = h)
  extra <- list(
    consumer = pin(13213L, "faedb8d48c90f7c96790aadd43ba51339236ed04340e61b52a358592e6fc6872"),
    loader = pin(10884L, "45f21c67674b5cd2b424e9c145662910361baee1019cbc68d366a9aab1fcc4ba"),
    mapping = pin(18639L, "fc3ec2b20608cb13b358337f0fc439c51876d8ecc75f285f21e80a11edd6ee48"),
    resource_validator = pin(13034L, "e93aae6769a1173bb7842f213152f56695f5b7b5621b71106c54181496f994ae"),
    metadata_parser = pin(10769L, "1f0549994bb9aef084c8278fceb336f3bf436abcb589fe4c1755add018dbd61b"))
  need(is.list(bundle) && length(bundle) == 20L && identical(utils::tail(names(bundle), 5L), names(extra)))
  for (k in names(extra)) need(identical(gse205906_ens101_worker_pin_v1(bundle[[k]]), extra[[k]]))
  consumer <- gse205906_ens101_worker_import_v1(bundle$consumer)
  pins <- c(consumer$gse205906_ens101_consumer_pins_v1(), extra)
  need(identical(names(bundle), names(pins)))
  consumer$gse205906_ens101_consumer_sources_v1(bundle[seq_len(15L)])
  need(identical(gse205906_ens101_worker_pin_v1(full_raw), pin(10612L,
    "ee6e22592b06ecdb31707234e97dc55b7b66cbad8dd11c364187af0624ece980")) &&
    identical(gse205906_ens101_worker_pin_v1(baseline_raw), pin(40558L,
    "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c")))
  list(consumer = consumer, pins = pins,
    loader = gse205906_ens101_worker_import_v1(bundle$loader),
    prefit = gse205906_ens101_worker_import_v1(bundle$prefit),
    mapping = gse205906_ens101_worker_import_v1(bundle$mapping),
    config = gse205906_ens101_worker_import_v1(bundle$config_adapter),
    metadata = gse205906_ens101_worker_import_v1(bundle$metadata_parser,
      c("gse205906_ens101_nofit_need_v2", "gse205906_ens101_nofit_overlay_v2",
        "gse205906_ens101_nofit_metadata_consistency_v2")),
    resource = gse205906_ens101_worker_import_v1(bundle$resource_validator,
      c("gse205906_validate_resource_tables_v1", "gse205906_resource_label_adapter_v1")),
    decoder = gse205906_ens101_worker_import_v1(bundle$helpers,
      c("v2r4_stop", "v2r4_assert", "v2r4_as_raw", "v2r4_crc32_le", "v2r4_u32_le",
        "v2r4_xz_vli", "v2r4_validate_xz_container", "v2r4_gzip_drain",
        "v2r4_mem_decompress_exact", "v2r4_unserialize_full", "v2r4_read_rds_raw")))
}

gse205906_ens101_worker_finite_v1 <- function(x, depth = 0L) {
  need <- gse205906_ens101_worker_need_v1
  need(depth <= 100L)
  if (is.list(x)) {
    need(is.null(names(x)) || !anyDuplicated(names(x)))
    for (v in x) gse205906_ens101_worker_finite_v1(v, depth + 1L)
  } else if (is.numeric(x)) need(all(is.finite(x)))
  invisible(TRUE)
}

gse205906_ens101_worker_resource_manifest_v1 <- function(table, pins) {
  need <- gse205906_ens101_worker_need_v1
  keys <- c("relative_path", "size_bytes", "sha256", "analysis_eligible")
  need(is.data.frame(table) && !anyDuplicated(names(table)) && all(keys %in% names(table)) &&
    all(vapply(table[keys], function(x) is.character(x) && !anyNA(x), logical(1))) &&
    !anyDuplicated(table$relative_path))
  for (k in names(pins)[5:10]) {
    p <- pins[[k]]; row <- table[table$relative_path == p$relative_path, , drop = FALSE]
    need(nrow(row) == 1L && identical(row$analysis_eligible, "TRUE") &&
      grepl("^[1-9][0-9]*$", row$size_bytes) &&
      identical(as.numeric(row$size_bytes), as.numeric(p$bytes)) && identical(row$sha256, p$sha256))
  }
  invisible(TRUE)
}

gse205906_ens101_worker_resources_v1 <- function(raw, api) {
  need <- gse205906_ens101_worker_need_v1
  pins <- gse205906_ens101_worker_inputs_v1()
  manifest <- as.data.frame(data.table::fread(text = rawToChar(raw$gene_set_manifest),
    sep = "\t", colClasses = "character", showProgress = FALSE))
  gse205906_ens101_worker_resource_manifest_v1(manifest, pins)
  original <- api$decoder$v2r4_read_rds_raw(raw$frozen_resource_RDS, "FROZEN_PUBLIC_GENE_SETS", "gzip")
  adapted <- api$resource$gse205906_resource_label_adapter_v1(original)
  # All ten RDS tables retain original domain validation. Only FIVE exports are
  # independently supplied here; other five comparison inputs are label copies.
  exports <- adapted
  role_map <- c(Reactome = "reactome", SynGO_sensitivity = "syngo_sensitivity",
    PROGENy = "progeny_human", DoRothEA_AB = "dorothea_rat", orthology = "ortholog_one_to_one")
  for (k in names(role_map)) exports[[role_map[[k]]]] <- data.table::fread(
    text = rawToChar(raw[[k]]), sep = "\t", na.strings = c("NA", ""), showProgress = FALSE)
  exports$syngo_sensitivity <- unique(exports$syngo_sensitivity[, c("term_id", "term_name", "gene_symbol"), with = FALSE])
  audit <- api$resource$gse205906_validate_resource_tables_v1(original, exports)
  api$resource$gse205906_validate_resource_tables_v1(original, adapted)
  need(identical(audit$status, "PASS_FROZEN_RESOURCE_CONTENT_EXPORT_AGREEMENT_ONLY"))
  resources <- adapted[c("reactome", "syngo_sensitivity", "progeny_human", "dorothea_rat")]
  names_for_pins <- c("reactome", "syngo_sensitivity", "progeny_human", "dorothea_rat", "orthology")
  resource_pins <- lapply(pins[5:9], function(p) p[c("bytes", "sha256")]); names(resource_pins) <- names_for_pins
  list(resources = resources, orthology = adapted$ortholog_one_to_one, pins = resource_pins,
    summary = list(externally_supplied_export_n = 5L, RDS_only_comparison_role_n = 5L,
      RDS_domain_role_n = 10L, factor_labels_not_codes = TRUE,
      selected_export_agreement = TRUE, manifest_fixed_six_entities_agreement = TRUE,
      source_or_count_filtered_programs_replaced = FALSE))
}

gse205906_ens101_worker_mapping_cfg_v1 <- function(full) {
  a <- full$analysis
  list(schema = "GSE205906_ENS101_MEMBERSHIP_CONFIG/v2", analysis = list(
    program_size_min = as.integer(a$program_size_min), program_size_max = as.integer(a$program_size_max),
    orthology_authority = a$orthology_authority,
    program_domains = list(PROGENy = list(source_program_n = as.integer(a$program_domains$PROGENy$source_program_n),
      exact_program_ids = unlist(a$program_domains$PROGENy$exact_program_ids, use.names = FALSE)),
      DoRothEA_AB = list(frozen_source_program_upper_bound =
        as.integer(a$program_domains$DoRothEA_AB$frozen_source_program_upper_bound)))))
}

gse205906_ens101_worker_serial_backend_v3 <- function(full) {
  need <- gse205906_ens101_worker_need_v1
  binding <- full$runtime_binding
  need(is.list(binding) && identical(binding$threads, 1L) &&
    identical(binding$backend, "BiocParallel_SerialParam"))
  need(identical(as.list(RNGkind()), binding$RNG_kind))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_before <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  need(requireNamespace("BiocParallel", quietly = TRUE))
  # Consume the existing engineering declaration before every deferred fit.
  # Do not set RNGseed/RNGkind or alter the original per-method set.seed calls.
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  observed <- BiocParallel::bpparam()
  need(identical(class(observed), structure("SerialParam", package = "BiocParallel")) &&
    identical(BiocParallel::bpnworkers(observed), 1L) &&
    identical(BiocParallel::bpisup(observed), FALSE) &&
    is.null(BiocParallel::bpRNGseed(observed)))
  need(identical(as.list(RNGkind()), binding$RNG_kind) &&
    identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), had_seed) &&
    (!had_seed || identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed_before)))
  invisible(TRUE)
}

gse205906_ens101_worker_prepare_v3 <- function(raw, source_bundle, full_raw, baseline_raw, runtime_pin, api, diagnostic_state) {
  need <- gse205906_ens101_worker_need_v1
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "PREPARE_CONFIG")
  full <- api$config$gse205906_ens101_check_full_model_config_v1(full_raw, baseline_raw)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "PARALLEL_BACKEND")
  gse205906_ens101_worker_serial_backend_v3(full)
  view <- api$config$gse205906_ens101_prefit_view_v1(full, baseline_raw)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "PREPARE_METADATA")
  ledger <- api$metadata$gse205906_ens101_nofit_overlay_v2(raw$overlay)
  manifest <- jsonlite::fromJSON(rawToChar(raw$manifest), simplifyVector = FALSE)
  gse205906_ens101_worker_finite_v1(manifest)
  need(identical(manifest$schema, "GSE205906_FORMAL_RUN_MANIFEST_156/v1") &&
    identical(manifest$gene_order, "lexically_sorted_stable_gene_id") &&
    identical(unlist(manifest$sample_order, use.names = FALSE), ledger$label))
  api$metadata$gse205906_ens101_nofit_metadata_consistency_v2(manifest, ledger)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "INPUT_LOAD")
  converted <- api$loader$gse205906_ens101_input_loader_v2(raw$matrix, ledger,
    unlist(manifest$gene_ids, use.names = FALSE), source_bundle$prefit, gtf_bytes = raw$GTF)
  input_pins <- c(lapply(raw[1:4], gse205906_ens101_worker_pin_v1),
    list(gse205906_ens101_worker_pin_v1(full_raw), gse205906_ens101_worker_pin_v1(source_bundle$prefit), runtime_pin))
  names(input_pins) <- c("new_counts", "new_annotation", "new_metadata", "fixed_gene_manifest",
    "new_config", "new_prefit_SOURCE", "statistical_runtime_lock")
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "PREFIT")
  prefit <- api$prefit$gse205906_ens101_effect_independent_prefit_v1(converted$inputs, view,
    source_bundle[c("model", "helpers", "annotation_adapter")], input_pins)
  need(identical(as.integer(prefit$passing_phase_n), 10L) && identical(as.integer(prefit$failed_phase_n), 0L))
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "RESOURCES")
  resource <- gse205906_ens101_worker_resources_v1(raw, api)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "MAPPING")
  mapping <- api$mapping$gse205906_ens101_mapping_v2(prefit, resource$resources,
    resource$orthology, gse205906_ens101_worker_mapping_cfg_v1(full),
    source_bundle[c("model", "prefit")], resource$pins)
  need(identical(mapping$status, "PASS_TEN_MEMBERSHIP_ONLY") &&
    identical(as.integer(mapping$passing_phase_n), 10L) && identical(as.integer(mapping$failed_phase_n), 0L))
  rawfive <- converted$inputs[c("count_matrix", "gene_annotation", "sample_ledger")]
  rawfive$gene_resources <- resource$resources; rawfive$orthology_table <- resource$orthology
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "CONSUMER_LOAD")
  consumer <- api$consumer$load_gse205906_ens101_checked_consumer_v1(full_raw, baseline_raw,
    source_bundle[seq_len(15L)], prefit, mapping)
  core <- gse205906_ens101_worker_import_v1(source_bundle$core)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "CONSUMER_PREPARE")
  handle <- consumer$prepare(rawfive, core$gse205906_ens101_raw_input_identities_v1(rawfive))
  list(handle = handle, prefit = prefit, mapping = mapping, resource = resource$summary,
    transformation_sha256 = converted$transformation$transformation_object_sha256)
}

gse205906_ens101_worker_fd_map_v1 <- function(fds, output_fd) {
  need <- gse205906_ens101_worker_need_v1
  roles <- names(gse205906_ens101_worker_inputs_v1())
  integer_fd <- function(x) is.numeric(x) && !is.logical(x) && length(x) == 1L &&
    is.finite(x) && x == floor(x) && x >= 3 && x <= .Machine$integer.max
  need(is.list(fds) && identical(names(fds), roles) && all(vapply(fds, integer_fd, logical(1))) && integer_fd(output_fd))
  values <- c(unlist(fds, use.names = FALSE), output_fd); need(!anyDuplicated(values))
  invisible(TRUE)
}

gse205906_ens101_worker_read_fd_v1 <- function(fd, pin) {
  need <- gse205906_ens101_worker_need_v1
  con <- file(paste0("/proc/self/fd/", format(fd, scientific = FALSE, trim = TRUE)), "rb")
  on.exit(close(con), add = TRUE)
  remaining <- pin$bytes; parts <- list()
  while (remaining > 0) {
    x <- readBin(con, "raw", n = min(1024^2, remaining)); need(length(x) > 0L)
    parts[[length(parts) + 1L]] <- x; remaining <- remaining - length(x)
  }
  need(length(readBin(con, "raw", n = 1L)) == 0L)
  raw <- do.call(c, parts)
  need(identical(gse205906_ens101_worker_pin_v1(raw), pin[c("bytes", "sha256")]))
  raw
}

gse205906_ens101_worker_tables_v1 <- function() {
  c("gene_edger", "gene_voom", "gene_deseq", "unsigned_camera", "unsigned_mroast",
    "unsigned_fgsea", "activity_results", "mapping_coverage", "loao_long", "loao_summary",
    "design_out", "contrast_out", "filter_out", "annotation_metrics", "gate")
}

gse205906_ens101_worker_tsv_v1 <- function(table) {
  need <- gse205906_ens101_worker_need_v1
  need(is.data.frame(table) && ncol(table) > 0L && !anyDuplicated(names(table)))
  quote_text <- function(x) paste0('"', gsub('"', '""', enc2utf8(x), fixed = TRUE), '"')
  column <- function(x) {
    need(is.null(dim(x)))
    if (is.list(x)) {
      # fgsea leadingEdge remains a list in the complete RDS. TSV cells have
      # explicit JSON vectors, not silently collapsed/lost scientific columns.
      need(!is.object(x) && all(vapply(x, function(v) is.null(v) ||
        (is.atomic(v) && !is.object(v) && is.null(dim(v)) &&
          (is.character(v) || is.numeric(v) || is.logical(v))), logical(1))))
      for (v in x) if (is.numeric(v)) need(all(is.finite(v[!is.na(v)])))
      return(quote_text(vapply(x, function(v) as.character(jsonlite::toJSON(v,
        auto_unbox = FALSE, na = "null", null = "null", digits = NA)), character(1))))
    }
    need(is.atomic(x))
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      value <- quote_text(x)
    } else if (is.logical(x) || is.integer(x)) value <- as.character(x)
    else if (is.double(x)) {
      need(all(is.finite(x[!is.na(x)]))); value <- sprintf("%.17g", x)
    } else need(FALSE)
    value[is.na(x)] <- "NA"; value
  }
  values <- lapply(table, column)
  rows <- if (nrow(table)) vapply(seq_len(nrow(table)), function(i)
    paste(vapply(values, `[[`, character(1), i), collapse = "\t"), character(1)) else character()
  charToRaw(paste0(paste(c(paste(quote_text(names(table)), collapse = "\t"), rows), collapse = "\n"), "\n"))
}

gse205906_ens101_worker_write_plan_v1 <- function(checked, provenance) {
  need <- gse205906_ens101_worker_need_v1
  tables <- checked$payload$result$tables
  need(is.list(tables) && identical(names(tables), gse205906_ens101_worker_tables_v1()) &&
    is.list(checked$payload$output_sha256) && identical(names(checked$payload$output_sha256), names(tables)))
  files <- lapply(tables, gse205906_ens101_worker_tsv_v1); names(files) <- paste0(names(tables), ".tsv")
  files$checked_payload.rds <- serialize(checked, NULL, version = 3L)
  need(all(lengths(files) <= 256 * 1024^2) && sum(lengths(files)) <= 512 * 1024^2)
  hashes <- lapply(files, gse205906_ens101_worker_pin_v1)
  manifest <- list(schema = "GSE205906_ENS101_CHECKED_MODEL_FILE_MANIFEST/v1",
    status = "WORKER_FILES_PENDING_PARENT_ACCEPTANCE", files = hashes,
    table_object_sha256 = checked$payload$output_sha256, provenance = provenance,
    encoding = list(TSV = "UTF8_QUOTED_TEXT_17_SIGNIFICANT_DIGITS_NA_TOKEN",
      RDS = "UNCOMPRESSED_R_SERIALIZATION_VERSION_3",
      list_column_TSV = "QUOTED_JSON_VECTOR_CELLS_OR_NULL",
      table_column_classes = lapply(tables, function(t) lapply(t, class))),
    parent_completion_marker_written = FALSE, scientific_acceptance = FALSE, production_go = FALSE)
  manifest_raw <- charToRaw(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null", na = "null", digits = NA))
  need(length(manifest_raw) <= 128 * 1024L)
  files$file_manifest.json <- manifest_raw
  need(sum(lengths(files)) <= 512 * 1024^2)
  list(files = files, file_pins = lapply(files, gse205906_ens101_worker_pin_v1),
    total_bytes = sum(lengths(files)), manifest = manifest)
}

gse205906_ens101_worker_write_v1 <- function(plan, output_fd, backend = NULL) {
  # backend is an explicit RAM-only test seam; native entry never accepts one.
  need <- gse205906_ens101_worker_need_v1
  expected <- c(paste0(gse205906_ens101_worker_tables_v1(), ".tsv"), "checked_payload.rds", "file_manifest.json")
  need(is.list(plan) && is.list(plan$files) && identical(names(plan$files), expected) &&
    all(vapply(plan$files, is.raw, logical(1))) && all(lengths(plan$files) <= 256 * 1024^2) &&
    length(plan$files$file_manifest.json) <= 128 * 1024L && sum(lengths(plan$files)) <= 512 * 1024^2)
  if (is.null(backend)) {
    root <- paste0("/proc/self/fd/", format(output_fd, scientific = FALSE, trim = TRUE))
    backend <- list(list = function() list.files(root, all.files = TRUE, no.. = TRUE),
      write = function(name, raw) {
        path <- paste0(root, "/", name); need(!file.exists(path) && !nzchar(Sys.readlink(path)))
        con <- file(path, "wb"); on.exit(close(con), add = TRUE); writeBin(raw, con)
      }, rename = function(old, new) {
        need(!file.exists(paste0(root, "/", new)) && !nzchar(Sys.readlink(paste0(root, "/", new))))
        need(file.rename(paste0(root, "/", old), paste0(root, "/", new)))
      }, read = function(name, n) {
        con <- file(paste0(root, "/", name), "rb"); on.exit(close(con), add = TRUE)
        x <- readBin(con, "raw", n = n); need(length(readBin(con, "raw", n = 1L)) == 0L); x
      })
  }
  need(is.list(backend) && identical(names(backend), c("list", "write", "rename", "read")) &&
    all(vapply(backend, is.function, logical(1))) && length(backend$list()) == 0L)
  for (name in expected) {
    temporary <- paste0(name, ".tmp")
    backend$write(temporary, plan$files[[name]]); backend$rename(temporary, name)
    observed <- backend$read(name, length(plan$files[[name]]))
    need(identical(gse205906_ens101_worker_pin_v1(observed), plan$file_pins[[name]]))
  }
  need(setequal(backend$list(), expected) && length(backend$list()) == length(expected))
  invisible(TRUE)
}

# Linux readlink() can return NA for an absent target. The old nzchar(NA)
# test refused that valid state. Metadata warnings remain fixed failures.
gse205906_ens101_worker_absent_name_v4 <- function(root, name) {
  need <- gse205906_ens101_worker_need_v1
  expected <- c(paste0(gse205906_ens101_worker_tables_v1(), ".tsv"),
    "checked_payload.rds", "file_manifest.json")
  need(is.character(root) && !is.object(root) && length(root) == 1L &&
    !is.na(root) && grepl("^/proc/self/fd/[0-9]+$", root) &&
    is.character(name) && !is.object(name) && length(name) == 1L &&
    !is.na(name) && name %in% c(expected, paste0(expected, ".tmp")))
  metadata <- function(action) withCallingHandlers(action(),
    warning = function(ignored) need(FALSE))
  listed <- metadata(function() list.files(root, all.files = TRUE, no.. = TRUE))
  need(is.character(listed) && !is.object(listed) && !anyNA(listed) &&
    !anyDuplicated(listed) && !name %in% listed)
  path <- paste0(root, "/", name)
  link <- metadata(function() Sys.readlink(path))
  present <- metadata(function() file.exists(path))
  need(is.character(link) && !is.object(link) && length(link) == 1L &&
    is.logical(present) && !is.object(present) && length(present) == 1L &&
    identical(present, FALSE) && (is.na(link) || identical(link, "")))
  invisible(TRUE)
}

gse205906_ens101_worker_write_v4 <- function(plan, output_fd, backend = NULL) {
  # backend is an explicit RAM-only test seam; native entry never accepts one.
  need <- gse205906_ens101_worker_need_v1
  expected <- c(paste0(gse205906_ens101_worker_tables_v1(), ".tsv"), "checked_payload.rds", "file_manifest.json")
  need(is.list(plan) && is.list(plan$files) && identical(names(plan$files), expected) &&
    all(vapply(plan$files, is.raw, logical(1))) && all(lengths(plan$files) <= 256 * 1024^2) &&
    length(plan$files$file_manifest.json) <= 128 * 1024L && sum(lengths(plan$files)) <= 512 * 1024^2)
  if (is.null(backend)) {
    root <- paste0("/proc/self/fd/", format(output_fd, scientific = FALSE, trim = TRUE))
    backend <- list(list = function() list.files(root, all.files = TRUE, no.. = TRUE),
      write = function(name, raw) {
        gse205906_ens101_worker_absent_name_v4(root, name)
        path <- paste0(root, "/", name)
        con <- file(path, "wb"); on.exit(close(con), add = TRUE); writeBin(raw, con)
      }, rename = function(old, new) {
        gse205906_ens101_worker_absent_name_v4(root, new)
        need(file.rename(paste0(root, "/", old), paste0(root, "/", new)))
      }, read = function(name, n) {
        con <- file(paste0(root, "/", name), "rb"); on.exit(close(con), add = TRUE)
        x <- readBin(con, "raw", n = n); need(length(readBin(con, "raw", n = 1L)) == 0L); x
      })
  }
  need(is.list(backend) && identical(names(backend), c("list", "write", "rename", "read")) &&
    all(vapply(backend, is.function, logical(1))) && length(backend$list()) == 0L)
  for (name in expected) {
    temporary <- paste0(name, ".tmp")
    backend$write(temporary, plan$files[[name]]); backend$rename(temporary, name)
    observed <- backend$read(name, length(plan$files[[name]]))
    need(identical(gse205906_ens101_worker_pin_v1(observed), plan$file_pins[[name]]))
  }
  need(setequal(backend$list(), expected) && length(backend$list()) == length(expected))
  invisible(TRUE)
}

# V2 only adds a caller-owned, process-local finite diagnostic state.
# Failure projection never catches the worker error, releases a result or writes.
gse205906_ens101_diagnostic_sites_v2 <- function() {
  c("INITIAL", "FRAME_VALIDATE", "FD_BINDINGS", "SOURCE_BINDINGS", "INPUT_BYTES",
    "PREPARE_CONFIG", "PARALLEL_BACKEND", "PREPARE_METADATA", "INPUT_LOAD", "PREFIT", "RESOURCES",
    "MAPPING", "CONSUMER_LOAD", "CONSUMER_PREPARE", "EXECUTE",
    "WRITE_PLAN", "WRITE_FILES", "SUMMARY", "UNKNOWN")
}

gse205906_ens101_diagnostic_controller_sites_v2 <- function() {
  c("PREPARED", "SOURCES_PRE", "INPUT_PRE", "GENE_REFERENCE",
    "INPUT_AFTER_GENE_REFERENCE", "PROGRAM_REFERENCE", "INPUT_AFTER_PROGRAM_REFERENCE",
    "COMPUTE", "INPUT_AFTER_COMPUTE", "RESULT_SCHEMA", "GENE_COMPARE",
    "PROGRAM_COMPARE", "SOURCES_POST", "INPUT_POST", "FINALIZE", "COMPLETE",
    "NOT_AVAILABLE")
}

gse205906_ens101_diagnostic_pick_v2 <- function(value, choices, fallback) {
  if (is.character(value) && !is.object(value) && length(value) == 1L &&
      !is.na(value)) for (choice in choices) if (identical(value, choice)) return(choice)
  fallback
}

gse205906_ens101_new_diagnostic_state_v2 <- function() {
  state <- new.env(parent = emptyenv())
  state$failure_site <- "INITIAL"
  state$controller_receipt <- NULL
  state
}

gse205906_ens101_diagnostic_mark_v2 <- function(state, site) {
  gse205906_ens101_worker_need_v1(is.environment(state) &&
    identical(parent.env(state), emptyenv()) &&
    identical(sort(ls(state, all.names = TRUE)), c("controller_receipt", "failure_site")) &&
    identical(gse205906_ens101_diagnostic_pick_v2(site,
      gse205906_ens101_diagnostic_sites_v2(), "UNKNOWN"), site) &&
    !identical(site, "UNKNOWN"))
  state$failure_site <- site
  invisible(TRUE)
}

gse205906_ens101_failure_summary_v2 <- function(state) {
  # Only literal allowlist members may leave RAM. No coercion of foreign values,
  # conditionMessage(), class(), events, binding, scientific object or receipt body.
  site <- "UNKNOWN"
  controller <- "NOT_AVAILABLE"
  tryCatch({
    if (is.environment(state) && identical(parent.env(state), emptyenv()) &&
        identical(sort(ls(state, all.names = TRUE)), c("controller_receipt", "failure_site"))) {
      site <- gse205906_ens101_diagnostic_pick_v2(state$failure_site,
        gse205906_ens101_diagnostic_sites_v2(), "UNKNOWN")
      if (is.function(state$controller_receipt)) {
        observed <- state$controller_receipt()
        if (is.list(observed) && !is.object(observed) &&
            !anyDuplicated(names(observed))) controller <-
          gse205906_ens101_diagnostic_pick_v2(observed$last_stage,
            gse205906_ens101_diagnostic_controller_sites_v2(), "NOT_AVAILABLE")
      }
    }
  }, error = function(ignored) NULL)
  list(schema = "GSE205906_ENS101_CHECKED_MODEL_FINITE_FAILURE/v1",
    status = "FINITE_FAILURE_DIAGNOSTIC_ONLY_NOT_RESULT_OR_ACCEPTANCE",
    failure_site = site, controller_last_stage = controller)
}

gse205906_ens101_checked_model_RAM_worker_v4 <- function(frame, diagnostic_state) {
  need <- gse205906_ens101_worker_need_v1
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "FRAME_VALIDATE")
  need(is.list(frame) && identical(names(frame), c("schema", "source_bundle", "full_config_raw",
    "baseline_config_raw", "input_fds", "output_directory_fd", "statistical_runtime_pin")) &&
    identical(frame$schema, "GSE205906_ENS101_CHECKED_MODEL_RAM_FRAME/v1"))
  p <- frame$statistical_runtime_pin
  need(is.list(p) && identical(names(p), c("bytes", "sha256")) && is.numeric(p$bytes) &&
    !is.logical(p$bytes) && length(p$bytes) == 1L && is.finite(p$bytes) && p$bytes > 0 &&
    p$bytes == floor(p$bytes) && p$bytes <= 2^53 && is.character(p$sha256) &&
    length(p$sha256) == 1L && !is.na(p$sha256) && grepl("^[0-9a-f]{64}$", p$sha256))
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "FD_BINDINGS")
  gse205906_ens101_worker_fd_map_v1(frame$input_fds, frame$output_directory_fd)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "SOURCE_BINDINGS")
  api <- gse205906_ens101_worker_sources_v1(frame$source_bundle, frame$full_config_raw, frame$baseline_config_raw)
  input_pins <- gse205906_ens101_worker_inputs_v1()
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "INPUT_BYTES")
  raw <- lapply(names(input_pins), function(k) gse205906_ens101_worker_read_fd_v1(frame$input_fds[[k]], input_pins[[k]]))
  names(raw) <- names(input_pins)
  prepared <- gse205906_ens101_worker_prepare_v3(raw, frame$source_bundle, frame$full_config_raw,
    frame$baseline_config_raw, p, api, diagnostic_state)
  # Exactly the pinned controller: gene refit, program refit, then one core.
  need(is.function(prepared$handle$receipt))
  diagnostic_state$controller_receipt <- prepared$handle$receipt
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "EXECUTE")
  checked <- prepared$handle$execute()
  provenance <- list(source_pins = api$pins, eleven_input_pins = input_pins,
    declared_statistical_runtime_pin = p, current_prefit_object_sha256 = prepared$prefit$output_object_sha256,
    current_mapping_object_sha256 = prepared$mapping$output_object_sha256,
    transformation_object_sha256 = prepared$transformation_sha256,
    actual_file_runtime_authority_authenticated_by_worker = FALSE)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "WRITE_PLAN")
  plan <- gse205906_ens101_worker_write_plan_v1(checked, provenance)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "WRITE_FILES")
  gse205906_ens101_worker_write_v4(plan, frame$output_directory_fd)
  gse205906_ens101_diagnostic_mark_v2(diagnostic_state, "SUMMARY")
  summary <- list(schema = "GSE205906_ENS101_CHECKED_MODEL_RAM_WORKER_SUMMARY/v1",
    status = "WORKER_CHECKED_TABLES_WRITTEN_PENDING_PARENT_ACCEPTANCE",
    table_n = 15L, written_file_n = 17L, total_bytes = plan$total_bytes,
    file_pins = plan$file_pins, input_pins = input_pins, source_pins = api$pins,
    consumer_status = checked$receipt$status, consumer_last_stage = checked$receipt$last_stage,
    numerical_gene_table_n = 3L, numerical_program_table_n = 12L,
    source_input_gene_n = 32883L, source_input_library_n = 9L, original_metadata_column_n = 53L,
    removed_input_gene_n = 0L, removed_input_library_n = 0L,
    preprocessing_phase_n = 10L, mapping_phase_n = 10L, resource_scope = prepared$resource,
    observation_unit = "NINE_DEPOSITED_LIBRARIES_CONDITIONAL_INDEPENDENCE",
    animal_pool_batch_dose_resolved = FALSE,
    parent_completion_marker_written = FALSE, actual_file_runtime_authority_authenticated_here = FALSE,
    scientific_acceptance = FALSE, production_go = FALSE, candidate_gate_eligible = FALSE,
    can_generate_or_upgrade_candidate = FALSE, formal_sequential_inference = FALSE,
    raw_vectors_metadata_and_private_errors_exported_in_summary = FALSE)
  # receipt contains only the original closed stage schema; no input/result tables.
  need(nchar(jsonlite::toJSON(summary, auto_unbox = TRUE, null = "null", na = "null", digits = NA), type = "bytes") <= 32768L)
  summary
}
