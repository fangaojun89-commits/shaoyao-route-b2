# SOURCE-only versioned configuration adapter. No file IO, fit or dispatch.
# Baseline statistical subtrees remain exact; new input provenance is separate.

gse205906_ens101_config_need_v1 <- function(ok, category) {
  if (!isTRUE(ok)) stop(category, call. = FALSE)
}

gse205906_ens101_config_hash_v1 <- function(object) {
  digest::digest(serialize(object, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_legacy_config_v1 <- function(legacy_config_bytes) {
  need <- gse205906_ens101_config_need_v1
  need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE),
       "ENS101_CONFIG_EXISTING_DEPENDENCIES_REQUIRED")
  need(is.raw(legacy_config_bytes) && length(legacy_config_bytes) == 40558L &&
         identical(digest::digest(legacy_config_bytes, "sha256", serialize = FALSE),
           "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"),
       "ENS101_CONFIG_FROZEN_BASELINE_BYTES_OR_SHA")
  jsonlite::fromJSON(rawToChar(legacy_config_bytes), simplifyVector = FALSE)
}

gse205906_ens101_full_model_config_v1 <- function(legacy_config_bytes) {
  baseline <- gse205906_ens101_legacy_config_v1(legacy_config_bytes)
  # Retain complete objects, including thresholds, methods, seed, CI, package
  # versions, synthetic controls and non-candidate flags. Legacy eligible
  # cardinalities are NOT current results; the reviewed contract overrides BOTH.
  retained <- c("design", "analysis", "software", "gate", "synthetic_contract")
  result <- list(
    schema = "GSE205906_ENS101_FULL_MODEL_CONFIG/v1",
    dataset_id = "GSE205906",
    model_version = "GSE205906_ENS101_CHECKED_CORE_V1_20261006",
    status = "NEEDS_REVIEW_INDEPENDENT_REFERENCES_AND_CHECKED_CONSUMER_NOT_IMPLEMENTED",
    reference = list(assembly = "Rnor_6.0", annotation = "Ensembl101",
      expected_gene_rows = 32883L, primary_strand_setting = 0L, count_unit = "gene_fragments"),
    design = baseline$design,
    analysis = baseline$analysis,
    software = baseline$software,
    gate = baseline$gate,
    synthetic_contract = baseline$synthetic_contract,
    input_snapshot = list(
      assembly = "Rnor_6.0", annotation_source = "Ensembl101",
      expected_gene_n_before_filter = 32883L,
      raw_annotation_columns = list("gene_id", "gene_name", "gene_biotype"),
      computational_annotation_columns = list("GeneID", "gene", "gene_biotype"),
      deterministic_annotation_projection = "gene_id_to_GeneID_gene_name_to_gene",
      source_attributes_preserved = TRUE, count_unit = "gene_fragments",
      selected_primary_strand_mode = 0L, strand_decision = "UNSTRANDED_CONDITIONAL_ON_ACTUAL_FROZEN_QC",
      minimum_nonblank_symbol_coverage = baseline$input_snapshot$minimum_nonblank_symbol_coverage,
      missing_symbol_policy = "ORIGINAL_GENE_ID_FALLBACK_AFTER_ORIGINAL_80_PERCENT_NONBLANK_RULE",
      source_gene_rows_removed = 0L, duplicate_symbol_gene_rows_collapsed = FALSE,
      historical_producer_species_reference_receipts_reusable = FALSE
    ),
    source_dependencies = list(
      baseline_config = list(bytes = 40558L,
        sha256 = "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"),
      model = list(bytes = 72490L,
        sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
      helpers = list(bytes = 175131L,
        sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
      annotation_adapter = list(bytes = 9943L,
        sha256 = "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470"),
      prefit = list(bytes = 18762L,
        sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
      mapping = list(bytes = 18639L,
        sha256 = "fc3ec2b20608cb13b358337f0fc439c51876d8ecc75f285f21e80a11edd6ee48"),
      program_contract = list(bytes = 12746L,
        sha256 = "bd86c229834a35317a701bfd53f727943d5238a44b99a94a549e3be8b5ba3adc")
    ),
    successor_changes = list(
      baseline_statistical_subtrees_exact = as.list(retained),
      replaced_domains = list("input_snapshot", "reference", "source_and_input_provenance", "namespace"),
      excluded_historical_entrypoint_contracts = list("governance", "external_runtime_authority",
        "qa_trust_root", "authority", "execution_contract", "required_output_payloads"),
      historical_eligible_counts_are_not_current = TRUE,
      current_eligible_counts_origin = "CURRENT_PREFIT_MAPPING_PROGRAM_CONTRACT_BOTH_ACTIVITY_FAMILIES",
      legacy_independent_rat_library_and_animals_strings_are_computational_labels_not_proof = TRUE
    ),
    input_provenance = list(
      status = "REQUIRES_ACTUAL_HELD_CALLER_BINDING_NOT_ESTABLISHED_HERE",
      accepted_counts_file_pin = NULL, raw_annotation_file_pin = NULL,
      metadata_file_pin = NULL, fixed_gene_manifest_pin = NULL,
      reference_FASTA_pin = NULL, reference_GTF_pin = NULL, reference_index_pins = NULL,
      fixed_nine_technical_acceptance_pin = NULL, actual_species_and_strand_pins = NULL,
      gene_resources_file_pin = NULL, orthology_file_pin = NULL,
      animals_pooling_batch_dose_resolved = FALSE,
      observation_unit = "NINE_DEPOSITED_LIBRARIES_CONDITIONAL_INDEPENDENCE",
      old_omitted_animal_column_is_library_identifier_alias = TRUE
    ),
    runtime_binding = list(
      expected_source_route = "EXISTING_LINUX_R_4_6_1_ROUTE_B2_PUBLIC_ANALYSIS_NOT_HISAT2_CONDA",
      original_runtime_receipt_reusable_without_fresh_verification = FALSE,
      actual_runtime_file_tree_pin = NULL, actual_session_and_backend_pin = NULL,
      threads = 1L, backend = "BiocParallel_SerialParam",
      RNG_kind = list("Mersenne-Twister", "Inversion", "Rejection"),
      current_binary_package_library_identity_authenticated_here = FALSE
    ),
    integration = list(
      independent_gene_reference_successor_implemented = FALSE,
      independent_program_reference_successor_implemented = FALSE,
      checked_consumer_successor_implemented = FALSE,
      production_caller_implemented = FALSE, model_dispatched = FALSE,
      production_go = FALSE, new_model_authorized_by_this_config = FALSE
    )
  )
  result
}

gse205906_ens101_check_full_model_config_v1 <- function(config_bytes, legacy_config_bytes) {
  need <- gse205906_ens101_config_need_v1
  need(is.raw(config_bytes) && length(config_bytes) > 0L && length(config_bytes) <= 65536L,
       "ENS101_CONFIG_RAW_BYTES_OR_CAPACITY")
  observed <- tryCatch(jsonlite::fromJSON(rawToChar(config_bytes), simplifyVector = FALSE),
    error = function(e) NULL)
  expected <- gse205906_ens101_full_model_config_v1(legacy_config_bytes)
  need(!is.null(observed) && identical(observed, expected), "ENS101_CONFIG_EXACT_SUCCESSOR_OBJECT_REQUIRED")
  expected
}

gse205906_ens101_prefit_view_v1 <- function(full_config, legacy_config_bytes) {
  need <- gse205906_ens101_config_need_v1
  need(identical(full_config, gse205906_ens101_full_model_config_v1(legacy_config_bytes)),
       "ENS101_PREFIT_VIEW_FULL_CONFIG_DRIFT")
  design_fields <- c("sample_ids_in_order", "groups_in_order", "group_levels", "coefficient_order",
                     "expected_rank", "contrasts")
  # Exact existing prefit schema; do not extend it with model/gate/runtime fields.
  list(
    schema = "GSE205906_ENS101_EFFECT_INDEPENDENT_PREFIT_CONFIG/v1",
    reference = list(assembly = "Rnor_6.0", annotation = "Ensembl101",
      expected_gene_rows = 32883L, primary_strand_setting = 0L, count_unit = "gene_fragments"),
    design = full_config$design[design_fields],
    analysis = full_config$analysis[c("minimum_cpm", "minimum_samples_at_cpm",
                                    "minimum_retained_genes_production", "random_seed")],
    input_snapshot = list(minimum_nonblank_symbol_coverage =
                           full_config$input_snapshot$minimum_nonblank_symbol_coverage)
  )
}
