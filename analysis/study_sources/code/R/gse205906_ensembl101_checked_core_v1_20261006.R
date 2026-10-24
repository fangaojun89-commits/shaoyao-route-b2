# Prospective in-memory checked-core successor; no file IO or launcher.
# Independent new references/checked consumer/native acceptance remain missing.
# prepare() is effect-free. Its compute() is deferred scientific code: do not
# invoke before the main-reviewed new references, consumer and actual caller.

gse205906_ens101_core_need_v1 <- function(ok, category) {
  if (!isTRUE(ok)) stop(category, call. = FALSE)
}

gse205906_ens101_core_hash_v1 <- function(object) {
  digest::digest(serialize(object, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_raw_input_identities_v1 <- function(inputs) {
  need <- gse205906_ens101_core_need_v1
  roles <- c("count_matrix", "gene_annotation", "sample_ledger", "gene_resources", "orthology_table")
  need(is.list(inputs) && identical(names(inputs), roles), "ENS101_CORE_RAW_FIVE_ROLE_ORDER")
  # Validate without canonicalizing or collapsing rows. The current prefit and
  # mapping bind serialize-v3 of these original objects, including their types.
  validate <- function(value, depth = 0L) {
    need(depth <= 12L, "ENS101_CORE_RAW_OBJECT_DEPTH")
    if (is.data.frame(value)) {
      need(identical(class(value), "data.frame") ||
             identical(class(value), c("data.table", "data.frame")), "ENS101_CORE_RAW_FRAME_CLASS")
      need(all(names(attributes(value)) %in% c("names", "row.names", "class", ".internal.selfref", "sorted")) &&
             !anyNA(names(value)) && !anyDuplicated(names(value)) && all(nzchar(names(value))),
           "ENS101_CORE_RAW_FRAME_ATTRIBUTES")
      for (column in value) validate(column, depth + 1L)
      return(invisible(TRUE))
    }
    if (is.matrix(value)) {
      need(is.numeric(value) && identical(class(value), c("matrix", "array")) &&
             all(names(attributes(value)) %in% c("dim", "dimnames")), "ENS101_CORE_RAW_MATRIX")
      return(invisible(TRUE))
    }
    if (is.list(value)) {
      need(!is.object(value) && !is.null(names(value)) && !anyNA(names(value)) &&
             !anyDuplicated(names(value)) && all(nzchar(names(value))) &&
             all(names(attributes(value)) %in% "names"), "ENS101_CORE_RAW_LIST")
      for (element in value) validate(element, depth + 1L)
      return(invisible(TRUE))
    }
    need((is.null(value) || is.character(value) || is.logical(value) ||
            is.integer(value) || is.double(value)) && !is.object(value) &&
           all(names(attributes(value)) %in% "names"), "ENS101_CORE_RAW_VALUE_TYPE")
    invisible(TRUE)
  }
  need(is.matrix(inputs$count_matrix) &&
         all(vapply(inputs[c("gene_annotation", "sample_ledger", "orthology_table")], is.data.frame, logical(1))) &&
         is.list(inputs$gene_resources), "ENS101_CORE_RAW_ROLE_TYPES")
  # The frozen annotation adapter indexes a plain data.frame, not data.table's
  # join syntax. The current GTF parser supplies this type. Do not silently
  # canonicalize an identity-bound annotation to make an unsupported type fit.
  need(identical(class(inputs$gene_annotation), "data.frame"), "ENS101_CORE_RAW_ANNOTATION_PLAIN_FRAME")
  for (value in inputs) validate(value)
  # A LIST, not the old consumer's vapply character vector. No source identity or
  # acceptance is created by this object hash function.
  lapply(inputs, gse205906_ens101_core_hash_v1)
}

gse205906_ens101_core_source_pins_v1 <- function() {
  list(
    old_core = list(bytes = 9530L, sha256 = "44754554efed7a51b0b2ec1bb68cb7310ac22367f1222edc7cf2d1a91cc90b28"),
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    helpers = list(bytes = 175131L, sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
    binding = list(bytes = 15697L, sha256 = "2aaa386df70c0031933c244eafa895d63c5dbd7a7be77277cfd4a57aceba538d"),
    annotation_adapter = list(bytes = 9943L, sha256 = "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470"),
    prefit = list(bytes = 18762L, sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
    config_adapter = list(bytes = 7852L, sha256 = "29ffe9705e675058c006d808305349275979efb58eeab2ea4236e5eb82fd339a"),
    program_contract = list(bytes = 12746L, sha256 = "bd86c229834a35317a701bfd53f727943d5238a44b99a94a549e3be8b5ba3adc")
  )
}

gse205906_ens101_derive_core_v1 <- function(old_core_text) {
  need <- gse205906_ens101_core_need_v1
  pin <- gse205906_ens101_core_source_pins_v1()$old_core
  need(is.character(old_core_text) && length(old_core_text) == 1L &&
         identical(digest::digest(charToRaw(old_core_text), "sha256", serialize = FALSE), pin$sha256),
       "ENS101_CORE_DERIVATION_FROZEN_SOURCE")
  text <- old_core_text
  amendments <- list()
  replace_once <- function(old, new, label) {
    positions <- gregexpr(old, text, fixed = TRUE)[[1L]]
    need(length(positions) == 1L && positions[[1L]] > 0L, "ENS101_CORE_DERIVATION_EXACT_ANCHOR")
    position <- positions[[1L]]
    text <<- paste0(substr(text, 1L, position - 1L), new,
                    substr(text, position + nchar(old), nchar(text)))
    amendments[[label]] <<- list(old_bytes = nchar(old, type = "bytes"), new_bytes = nchar(new, type = "bytes"),
      old_sha256 = digest::digest(charToRaw(old), "sha256", serialize = FALSE),
      new_sha256 = digest::digest(charToRaw(new), "sha256", serialize = FALSE))
  }
  replace_once("load_gse205906_bound_core_v2 <- function(project_root, program_contract)",
    "load_gse205906_ens101_core_engine_v1 <- function(source_bundle, cfg, program_contract, input_prepare)", "signature")
  replace_once(paste0('    "config/GSE205906_bulk_model_v2_r4.json" =\n',
      '      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"'),
    paste0('    "config/GSE205906_Ensembl101_full_model_v1_20261006.json" =\n',
      '      digest::digest(source_bundle$full_config, "sha256", serialize = FALSE)'), "new_config_pin")
  start <- regexpr("  capture_source <- function(relative, expected) {", text, fixed = TRUE)[[1L]]
  finish <- regexpr('  if (!requireNamespace("digest", quietly = TRUE) ||', text, fixed = TRUE)[[1L]]
  need(start > 0L && finish > start, "ENS101_CORE_CAPTURE_BOUNDARY")
  capture_old <- substr(text, start, finish - 1L)
  capture_new <- paste0(
    '  capture_source <- function(relative, expected) {\n',
    '    roles <- c("code/R/analyze_gse205906_bulk_model_v2_r4.R" = "model",\n',
    '      "code/R/gse205906_bulk_model_v2_r4_lib.R" = "helpers",\n',
    '      "config/GSE205906_Ensembl101_full_model_v1_20261006.json" = "full_config")\n',
    '    if (!relative %in% names(roles)) fail("ENS101_CORE_CAPTURE_ROLE")\n',
    '    payload <- source_bundle[[roles[[relative]]]]\n',
    '    if (!is.raw(payload) || length(payload) > 1000000L ||\n',
    '        !identical(digest::digest(payload, "sha256", serialize = FALSE), expected)) {\n',
    '      fail("CORE_SOURCE_IDENTITY_MISMATCH")\n',
    '    }\n    rawToChar(payload)\n  }\n')
  replace_once(capture_old, capture_new, "held_caller_source_bytes_not_path_reopen")
  replace_once("cfg <- program_contract$effective_config(jsonlite::fromJSON(source_text[[3L]], simplifyVector = FALSE))",
    "cfg <- program_contract$effective_config(cfg)", "explicit_verified_full_config")
  replace_once('orthology_table, domain = "SYNTHETIC")',
    'orthology_table, domain = "PRODUCTION")', "original_production_statistical_settings")
  replace_once('if (!identical(domain, "SYNTHETIC")) fail("BOUND_PUBLIC_EXECUTION_NOT_ACCEPTED")',
    'if (!identical(domain, "PRODUCTION")) fail("ENS101_CORE_PRODUCTION_SETTINGS_DOMAIN_ONLY")', "domain_not_authority")
  replace_once("    context <- new.env(parent = helpers)", paste0(
    "    prepared <- input_prepare(list(count_matrix=count_matrix, gene_annotation=gene_annotation,\n",
    "      sample_ledger=sample_ledger, gene_resources=gene_resources, orthology_table=orthology_table))\n",
    "    context <- new.env(parent = helpers)"), "raw_guard_before_internal_projection")
  replace_once("context$ledger <- helpers$v2r4_validate_ledger(data.table::copy(sample_ledger), cfg)",
    "context$ledger <- data.table::copy(prepared$raw_inputs$sample_ledger)", "original_53_columns_not_fabricated_animals")
  replace_once("context$assembled <- list(gene_annotation = data.table::copy(gene_annotation))",
    "context$assembled <- list(gene_annotation = data.table::copy(prepared$projection$computational_interface))",
    "original_adapter_computational_interface_only")
  replace_once('      source_amendment = "BOUND_MEMBERSHIP_AND_EXPLICIT_LOAO_FAILURES_V1",', paste0(
    '      source_amendment = "ENS101_RAW_BINDING_INTERNAL_ANNOTATION_PROJECTION_AND_CURRENT_MEMBERSHIP_V1",\n',
    '      raw_input_identities = prepared$raw_input_identities,\n',
    '      annotation_projection_object_sha256 = prepared$projection_object_sha256,\n',
    '      source_annotation_origin = "gene_id_gene_name_Ensembl101_NOT_RefSeq",\n',
    '      observation_unit = "NINE_DEPOSITED_LIBRARIES_CONDITIONAL_INDEPENDENCE",\n',
    '      legacy_omitted_animal_is_library_identifier_alias = TRUE,\n',
    '      independent_references_and_checked_consumer_implemented_here = FALSE,'), "explicit_non_acceptance_provenance")
  # No scientific-span, fit function, tolerance, RNG, contrast or threshold edit.
  list(text = text, amendments = amendments,
       derived_source_sha256 = digest::digest(charToRaw(text), "sha256", serialize = FALSE))
}

load_gse205906_ens101_checked_core_v1 <- function(config_bytes, legacy_config_bytes,
                                                source_bundle, prefit, mapping) {
  need <- gse205906_ens101_core_need_v1
  need(requireNamespace("digest", quietly = TRUE) && requireNamespace("data.table", quietly = TRUE),
       "ENS101_CORE_EXISTING_DEPENDENCIES_REQUIRED")
  pins <- gse205906_ens101_core_source_pins_v1()
  need(is.list(source_bundle) && identical(names(source_bundle), names(pins)), "ENS101_CORE_EXACT_SOURCE_ROLES")
  for (role in names(pins)) need(is.raw(source_bundle[[role]]) && length(source_bundle[[role]]) == pins[[role]]$bytes &&
    identical(digest::digest(source_bundle[[role]], "sha256", serialize = FALSE), pins[[role]]$sha256),
    "ENS101_CORE_SOURCE_BYTES_OR_SHA_DRIFT")
  import_functions <- function(bytes, names_needed) {
    environment <- new.env(parent = baseenv())
    parsed <- as.list(parse(text = rawToChar(bytes), keep.source = FALSE))
    selected <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% names_needed, parsed)
    found <- vapply(selected, function(x) as.character(x[[2L]]), character(1))
    need(setequal(found, names_needed) && !anyDuplicated(found), "ENS101_CORE_FUNCTION_COVERAGE")
    for (x in selected) {
      need(is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), "ENS101_CORE_NOT_FUNCTION_DECLARATION")
      eval(x, environment)
    }
    lockEnvironment(environment, bindings = TRUE)
    environment
  }
  config_api <- import_functions(source_bundle$config_adapter, c("gse205906_ens101_config_need_v1",
    "gse205906_ens101_config_hash_v1", "gse205906_ens101_legacy_config_v1",
    "gse205906_ens101_full_model_config_v1", "gse205906_ens101_check_full_model_config_v1", "gse205906_ens101_prefit_view_v1"))
  cfg <- config_api$gse205906_ens101_check_full_model_config_v1(config_bytes, legacy_config_bytes)
  cfg_identity <- gse205906_ens101_core_hash_v1(cfg)
  prefit_api <- import_functions(source_bundle$prefit, c("gse205906_ens101_prefit_need_v1",
    "gse205906_ens101_prefit_config_v1", "gse205906_ens101_prefit_ledger_v1"))
  view <- config_api$gse205906_ens101_prefit_view_v1(cfg, legacy_config_bytes)
  prefit_api$gse205906_ens101_prefit_config_v1(view)
  need(identical(prefit$config_object_sha256, gse205906_ens101_core_hash_v1(view)),
       "ENS101_CORE_CURRENT_PREFIT_CONFIG_VIEW_BINDING")
  adapter <- import_functions(source_bundle$annotation_adapter,
    c("gse205906_ens101_need_v1", "gse205906_ens101_character_v1", "gse205906_ens101_annotation_v1"))
  contract_api <- import_functions(source_bundle$program_contract,
    c("gse205906_ens101_program_hash_v2", "gse205906_ens101_program_contract_v2"))
  contract <- contract_api$gse205906_ens101_program_contract_v2(prefit, mapping,
    source_bundle[c("binding", "model", "helpers")], gse205906_ens101_raw_input_identities_v1)
  input_prepare <- function(raw_inputs) {
    identities <- gse205906_ens101_raw_input_identities_v1(raw_inputs)
    need(identical(gse205906_ens101_core_hash_v1(cfg), cfg_identity), "ENS101_CORE_CONFIG_CHANGED")
    prefit_api$gse205906_ens101_prefit_ledger_v1(raw_inputs$sample_ledger)
    counts <- raw_inputs$count_matrix
    need(identical(dim(counts), c(32883L, 9L)) &&
           identical(colnames(counts), paste0("GSM", 6235292:6235300)) &&
           identical(rownames(counts), prefit$counts_gene_ids) &&
           all(is.finite(counts)) && all(counts >= 0) && all(counts == floor(counts)) &&
           all(counts <= .Machine$integer.max), "ENS101_CORE_FIXED_INTEGER_COUNT_DOMAIN")
    projection <- adapter$gse205906_ens101_annotation_v1(raw_inputs$gene_annotation,
      prefit$counts_gene_ids, "FIXED_ENS101")
    need(identical(projection, prefit$annotation_projection), "ENS101_CORE_CURRENT_PREFIT_ANNOTATION_PROJECTION_BINDING")
    need(identical(identities, gse205906_ens101_raw_input_identities_v1(raw_inputs)), "ENS101_CORE_RAW_INPUT_CHANGED_DURING_PROJECTION")
    list(raw_inputs = raw_inputs, raw_input_identities = identities, projection = projection,
         projection_object_sha256 = gse205906_ens101_core_hash_v1(projection))
  }
  derived <- gse205906_ens101_derive_core_v1(rawToChar(source_bundle$old_core))
  engine_api <- import_functions(charToRaw(derived$text), "load_gse205906_ens101_core_engine_v1")
  engine <- engine_api$load_gse205906_ens101_core_engine_v1(
    c(source_bundle[c("model", "helpers")], list(full_config = config_bytes)), cfg, contract, input_prepare)
  prepare <- function(raw_inputs) {
    # Exact current raw object identities FIRST. No projected object is accepted
    # as raw evidence. No old public snapshot/wrapper is called.
    contract$check_inputs(raw_inputs, "PRODUCTION")
    serialized <- serialize(raw_inputs, NULL, version = 3L)
    snapshot <- unserialize(serialized)
    snapshot_identity <- gse205906_ens101_core_hash_v1(snapshot)
    prepared <- input_prepare(snapshot)
    check_snapshot <- function() {
      need(identical(gse205906_ens101_core_hash_v1(snapshot), snapshot_identity) &&
             identical(gse205906_ens101_core_hash_v1(prepared$projection), prepared$projection_object_sha256),
           "ENS101_CORE_CAPTURED_RAW_OR_PROJECTION_DRIFT")
      contract$check_inputs(snapshot, "PRODUCTION")
      invisible(TRUE)
    }
    inspect <- function() {
      check_snapshot()
      list(status = cfg$status, model_run = FALSE, production_go = FALSE,
        raw_input_identities = prepared$raw_input_identities,
        raw_snapshot_object_sha256 = snapshot_identity,
        annotation_projection_object_sha256 = prepared$projection_object_sha256,
        projection_computational_object_sha256 = gse205906_ens101_core_hash_v1(prepared$projection$computational_interface),
        raw_gene_rows = 32883L, count_library_n = 9L, raw_rows_removed = 0L,
        raw_attributes_preserved = TRUE, duplicate_symbol_gene_rows_collapsed = FALSE,
        legacy_omitted_animal_is_library_identifier_alias = TRUE,
        program_binding_sha256 = contract$binding_sha256,
        current_prefit_object_sha256 = contract$current_prefit_object_sha256,
        current_mapping_object_sha256 = contract$current_mapping_object_sha256,
        independent_references_implemented = FALSE, checked_consumer_implemented = FALSE,
        input_file_runtime_and_producer_authentication_supplied_here = FALSE)
    }
    compute <- function() {
      # This is the future in-memory core only, NOT dispatch/checked acceptance.
      # It changes RNG and calls original fits; SOURCE tests NEVER call it.
      check_snapshot()
      result <- do.call(engine$compute, c(snapshot, list(domain = "PRODUCTION")))
      check_snapshot()
      result
    }
    list(inspect = inspect, compute = compute)
  }
  list(status = cfg$status, prepare = prepare, config = cfg, prefit_view = view,
    source_pins = pins, engine_source_pins = engine$source_pins, source_spans = engine$source_spans,
    engine_derived_source_sha256 = derived$derived_source_sha256,
    interface_amendments = derived$amendments, packages = engine$packages,
    config_raw_sha256 = digest::digest(config_bytes, "sha256", serialize = FALSE),
    config_object_sha256 = cfg_identity, model_run = FALSE, production_go = FALSE,
    independent_references_implemented = FALSE, checked_consumer_implemented = FALSE,
    runtime_binary_library_custody_authenticated_here = FALSE,
    limitation = "Only SOURCE/config/raw-object/projection integration. Deferred core fits are not called by loading or prepare. Two new independent references, complete checked consumer, actual accepted input/runtime caller and separate model execution decision remain required.")
}
