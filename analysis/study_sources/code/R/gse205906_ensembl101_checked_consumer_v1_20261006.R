# Fixed in-memory consumer. No path IO, launcher or execution authority.
# Loading/prepare never call a statistical compute/refit. execute() is deferred
# science: caller authentication and the main's execution decision are separate.

gse205906_ens101_consumer_need_v1 <- function(ok, category) {
  if (!isTRUE(ok)) stop(category, call. = FALSE)
}

gse205906_ens101_consumer_hash_v1 <- function(value) {
  digest::digest(serialize(value, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_consumer_pins_v1 <- function() {
  pin <- function(n, h) list(bytes = n, sha256 = h)
  list(
    stage_machine = pin(15674L, "37bacea977472600ae9879d7e43448f953c5a0ddd93f0967e7586f8a1a5d3d49"),
    core = pin(17175L, "ac070364dca4369f9a9cb196550a8c6f3daced3114184f1572d5274ec03261e6"),
    gene_reference = pin(13197L, "4c6593f0ed6cdd6838699f8d844dbd7aa321baca990b2dac79abe9ee7c1f76d1"),
    program_reference = pin(9509L, "917c2ff7ec82615288c4a957a76beb7f67b040ebb4bc6d2c8e3169cb97f856ed"),
    old_core = pin(9530L, "44754554efed7a51b0b2ec1bb68cb7310ac22367f1222edc7cf2d1a91cc90b28"),
    model = pin(72490L, "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    helpers = pin(175131L, "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
    binding = pin(15697L, "2aaa386df70c0031933c244eafa895d63c5dbd7a7be77277cfd4a57aceba538d"),
    annotation_adapter = pin(9943L, "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470"),
    prefit = pin(18762L, "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
    config_adapter = pin(7852L, "29ffe9705e675058c006d808305349275979efb58eeab2ea4236e5eb82fd339a"),
    program_contract = pin(12746L, "bd86c229834a35317a701bfd53f727943d5238a44b99a94a549e3be8b5ba3adc"),
    old_gene_reference = pin(9590L, "4188175759fc7db35b18269bf54e92be0b217517226a676ed4f444aaabf8968f"),
    validator = pin(104432L, "c63296157a833611ef92b17f1cbb6b27ab3c23c1028cab79fd706aaab564225d"),
    old_program_reference = pin(25221L, "4e0f3cd999b327e06de1388f8fe92472a1af20bf62e9d71d7cea24ed3b44d493")
  )
}

gse205906_ens101_consumer_sources_v1 <- function(source_bundle) {
  need <- gse205906_ens101_consumer_need_v1
  pins <- gse205906_ens101_consumer_pins_v1()
  need(is.list(source_bundle) && identical(names(source_bundle), names(pins)),
       "ENS101_CONSUMER_SOURCE_ROLE_ORDER")
  for (role in names(pins)) {
    raw <- source_bundle[[role]]
    need(is.raw(raw) && length(raw) == pins[[role]]$bytes &&
           identical(digest::digest(raw, "sha256", serialize = FALSE), pins[[role]]$sha256),
         "ENS101_CONSUMER_COMPLETE_SOURCE_DRIFT")
  }
  invisible(TRUE)
}

gse205906_ens101_consumer_import_v1 <- function(raw, wanted = NULL) {
  need <- gse205906_ens101_consumer_need_v1
  declarations <- as.list(parse(text = rawToChar(raw), keep.source = FALSE))
  if (!is.null(wanted)) declarations <- Filter(function(x) is.call(x) &&
    identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
    as.character(x[[2L]]) %in% wanted, declarations)
  names_found <- character()
  private <- new.env(parent = baseenv())
  # The original source-controlled loaders use this explicit stats binding.
  private$setNames <- stats::setNames
  for (x in declarations) {
    need(is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
           is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")),
         "ENS101_CONSUMER_FUNCTION_DECLARATIONS_ONLY")
    names_found <- c(names_found, as.character(x[[2L]]))
    eval(x, private)
  }
  need(!anyDuplicated(names_found) && (is.null(wanted) || setequal(names_found, wanted)),
       "ENS101_CONSUMER_FUNCTION_COVERAGE")
  lockEnvironment(private, bindings = TRUE)
  private
}

gse205906_ens101_consumer_stage_v1 <- function(raw) {
  # Only the evidence label changes. The original controller ordering, single
  # execution, error rethrow and no-partial-release semantics are untouched.
  old <- 'c("PINNED_LOCAL_SYNTHETIC", "CONTROL_TEST_DOUBLE")'
  new <- 'c("PINNED_LOCAL_SYNTHETIC", "CONTROL_TEST_DOUBLE", "PINNED_ENS101_IN_MEMORY")'
  text <- rawToChar(raw)
  hits <- gregexpr(old, text, fixed = TRUE)[[1L]]
  gse205906_ens101_consumer_need_v1(length(hits) == 1L && hits[[1L]] > 0L,
                                  "ENS101_CONSUMER_STAGE_EXACT_SCOPE_ANCHOR")
  text <- paste0(substr(text, 1L, hits[[1L]] - 1L), new,
                 substr(text, hits[[1L]] + nchar(old), nchar(text)))
  gse205906_ens101_consumer_import_v1(charToRaw(text),
    "gse205906_checked_stage_machine_v1")$gse205906_checked_stage_machine_v1
}

gse205906_ens101_consumer_tables_v1 <- function() {
  c("gene_edger", "gene_voom", "gene_deseq", "unsigned_camera", "unsigned_mroast",
    "unsigned_fgsea", "activity_results", "mapping_coverage", "loao_long", "loao_summary",
    "design_out", "contrast_out", "filter_out", "annotation_metrics", "gate")
}

gse205906_ens101_consumer_handle_v1 <- function(raw_inputs, expected_identities,
    components, identities, source_check, stage_machine, evidence_class, source_pins) {
  # Explicit control seam for artificial tests; injected components never prove
  # real SOURCE identity. The fixed loader below supplies only pinned loaders.
  need <- gse205906_ens101_consumer_need_v1
  hash <- gse205906_ens101_consumer_hash_v1
  need(identical(names(components), c("core", "gene", "program")) &&
         all(vapply(components, function(x) is.list(x) && is.function(x$prepare), logical(1))) &&
         is.function(identities) && is.function(source_check) && is.function(stage_machine),
       "ENS101_CONSUMER_COMPONENT_INTERFACE")
  observed <- identities(raw_inputs)
  need(identical(observed, expected_identities), "ENS101_CONSUMER_RAW_IDENTITY_MISMATCH")
  snapshot <- unserialize(serialize(raw_inputs, NULL, version = 3L))
  snapshot_hash <- hash(snapshot)
  prepared <- lapply(components, function(x) x$prepare(snapshot))
  need(is.function(prepared$core$compute) && is.function(prepared$gene$refit) &&
         is.function(prepared$program$refit) &&
         all(vapply(prepared, function(x) is.function(x$inspect), logical(1))),
       "ENS101_CONSUMER_DEFERRED_INTERFACE")
  inspections <- lapply(prepared, function(x) x$inspect())
  common <- c("raw_input_identities", "raw_snapshot_object_sha256",
              "annotation_projection_object_sha256", "current_prefit_object_sha256",
              "current_mapping_object_sha256")
  for (x in inspections) need(all(common %in% names(x)) &&
    identical(x$raw_input_identities, observed) &&
    identical(x$raw_snapshot_object_sha256, snapshot_hash), "ENS101_CONSUMER_SNAPSHOT_BINDING")
  for (key in common) need(all(vapply(inspections[-1L], function(x)
    identical(x[[key]], inspections$core[[key]]), logical(1))), "ENS101_CONSUMER_COMPONENT_BINDING")
  check_input <- function() {
    need(identical(hash(snapshot), snapshot_hash) && identical(identities(snapshot), observed),
         "ENS101_CONSUMER_CAPTURED_INPUT_CHANGED")
    now <- lapply(prepared, function(x) x$inspect())
    need(identical(now, inspections), "ENS101_CONSUMER_PREPARED_INSPECTION_DRIFT")
    TRUE
  }
  tables <- gse205906_ens101_consumer_tables_v1()
  cfg <- components$core$config
  check_result <- function(result) {
    need(is.list(result) && identical(result$status, "IN_MEMORY_BOUND_SCIENCE_ONLY_NOT_PRODUCTION_ACCEPTANCE") &&
      identical(result$domain, "PRODUCTION") && identical(result$production_go, FALSE) &&
      is.list(result$tables) && identical(names(result$tables), tables), "ENS101_CONSUMER_RESULT_INTERFACE")
    need(identical(result$source_pins, components$core$engine_source_pins) &&
      identical(result$source_spans, components$core$source_spans) &&
      identical(result$packages, components$core$packages) &&
      identical(result$random_seed, cfg$analysis$random_seed) &&
      identical(result$raw_input_identities, observed) &&
      identical(result$annotation_projection_object_sha256, inspections$core$annotation_projection_object_sha256),
      "ENS101_CONSUMER_RESULT_PROVENANCE")
    need(is.list(result$program_binding) &&
      identical(result$program_binding$binding_sha256, inspections$core$program_binding_sha256) &&
      identical(result$program_binding$status, "PASS_BOUND_FAMILIES_AND_COMPLETE_REPORTING_ONLY") &&
      identical(result$program_binding$checked_gene_universes, 10L) &&
      identical(result$program_binding$checked_program_families, 22L), "ENS101_CONSUMER_COMPLETE_MEMBERSHIP")
    flags <- c("candidate_gate_eligible", "can_generate_or_upgrade_candidate", "formal_sequential_inference")
    for (table in result$tables) need(is.data.frame(table) && !anyDuplicated(names(table)) &&
      all(c("analysis_role", flags) %in% names(table)) && is.character(table$analysis_role) &&
      !anyNA(table$analysis_role) && all(table$analysis_role == cfg$gate$analysis_role) &&
      all(vapply(flags, function(key) identical(table[[key]], rep(FALSE, nrow(table))), logical(1))),
      "ENS101_CONSUMER_RESULT_FALSE_CLAIMS")
    TRUE
  }
  steps <- list(check_sources = source_check, check_input = check_input,
    gene_reference = function() prepared$gene$refit(),
    program_reference = function() prepared$program$refit(),
    compute = function() prepared$core$compute(),
    gene_compare = function(reference, result) {
      need(identical(reference$raw_input_identities, observed), "ENS101_CONSUMER_GENE_REFIT_RAW_BINDING")
      checked <- reference$compare(result)
      need(identical(checked$verified_tables, tables[1:3]), "ENS101_CONSUMER_GENE_TABLE_COVERAGE")
      checked
    },
    program_compare = function(reference, result) {
      need(identical(reference$raw_input_identities, observed), "ENS101_CONSUMER_PROGRAM_REFIT_RAW_BINDING")
      checked <- reference$compare(result)
      need(identical(checked$verified_tables, tables[-(1:3)]) &&
        identical(checked$program_binding, result$program_binding), "ENS101_CONSUMER_PROGRAM_TABLE_COVERAGE")
      checked
    }, check_result = check_result,
    finalize = function(result) list(result = result,
      output_sha256 = as.list(vapply(result$tables, function(x) hash(as.data.frame(x)), character(1))),
      observation_unit = "NINE_DEPOSITED_LIBRARIES_CONDITIONAL_INDEPENDENCE",
      raw_rows_removed = 0L, sample_columns_removed = 0L,
      original_metadata_column_n = ncol(snapshot$sample_ledger),
      input_file_runtime_and_producer_authentication_supplied_here = FALSE,
      limitation = "Complete in-memory comparisons, not biological replication, causal/mechanistic evidence, native custody, output admission or independent personnel audit. Original omitted_animal is a library-ID alias."))
  machine <- stage_machine(steps, binding = list(raw_input_identities = observed,
    raw_snapshot_object_sha256 = snapshot_hash, source_pins = source_pins,
    component_inspections = inspections, domain = "PRODUCTION",
    actual_file_runtime_and_model_authority = FALSE), evidence_class = evidence_class)
  list(execute = machine$execute, receipt = machine$receipt,
       inspect = function() { check_input(); machine$receipt() })
}

load_gse205906_ens101_checked_consumer_v1 <- function(config_bytes, legacy_config_bytes,
    source_bundle, prefit, mapping) {
  gse205906_ens101_consumer_sources_v1(source_bundle)
  import <- gse205906_ens101_consumer_import_v1
  core_api <- import(source_bundle$core)
  gene_api <- import(source_bundle$gene_reference)
  program_api <- import(source_bundle$program_reference)
  core_roles <- c("old_core", "model", "helpers", "binding", "annotation_adapter", "prefit", "config_adapter", "program_contract")
  common <- c(source_bundle[c("helpers", "model", "binding", "annotation_adapter", "prefit")],
    list(identity_adapter = source_bundle$core), source_bundle[c("config_adapter", "program_contract")])
  components <- list(
    core = core_api$load_gse205906_ens101_checked_core_v1(config_bytes, legacy_config_bytes,
      source_bundle[core_roles], prefit, mapping),
    gene = gene_api$load_gse205906_ens101_gene_reference_v1(config_bytes, legacy_config_bytes,
      c(list(old_reference = source_bundle$old_gene_reference, validator = source_bundle$validator), common), prefit, mapping),
    program = program_api$load_gse205906_ens101_program_reference_v1(config_bytes, legacy_config_bytes,
      c(list(old_reference = source_bundle$old_program_reference, reference_bootstrap = source_bundle$gene_reference), common), prefit, mapping))
  stage <- gse205906_ens101_consumer_stage_v1(source_bundle$stage_machine)
  pins <- gse205906_ens101_consumer_pins_v1()
  check_sources <- function() { gse205906_ens101_consumer_sources_v1(source_bundle); TRUE }
  prepare <- function(raw_inputs, expected_identities) {
    gse205906_ens101_consumer_handle_v1(raw_inputs, expected_identities, components,
      core_api$gse205906_ens101_raw_input_identities_v1, check_sources, stage,
      "PINNED_ENS101_IN_MEMORY", pins)
  }
  list(prepare = prepare, source_pins = pins, config = components$core$config,
    source_config_and_reference_loaders_bound = TRUE, production_go = FALSE,
    numerical_core_or_reference_run_by_loading = FALSE,
    actual_file_runtime_and_model_authority = FALSE)
}
