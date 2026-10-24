# NEW Rnor_6.0 / Ensembl101 annotation interface. Inert base-R pure functions.
# No file IO, count ingestion/filtering/TMM, resource/membership computation,
# fit, score, model permission, config mutation or historical result reuse.
# Real input/SOURCE/runtime pin authentication belongs to the future reviewed
# new prefit caller; types or SOURCE_TOY labels do not establish that reality.

gse205906_ens101_need_v1 <- function(condition, category) {
  if (!isTRUE(condition)) {
    stop(category, call. = FALSE)
  }
}

gse205906_ens101_character_v1 <- function(values, nonblank = TRUE) {
  valid <- is.character(values) && is.null(dim(values)) &&
    !anyNA(values) && !anyDuplicated(values)
  if (valid && nonblank) {
    valid <- all(nzchar(values)) && identical(values, trimws(values))
  }
  valid
}

gse205906_ens101_annotation_v1 <- function(
    annotation, fixed_gene_ids, domain = "FIXED_ENS101") {
  need <- gse205906_ens101_need_v1
  need(identical(domain, "FIXED_ENS101") || identical(domain, "SOURCE_TOY"),
       "ENS101_INTERFACE_DOMAIN")
  need(is.data.frame(annotation) && !anyDuplicated(names(annotation)) &&
         all(c("gene_id", "gene_name", "gene_biotype") %in% names(annotation)),
       "ENS101_SOURCE_ATTRIBUTE_COLUMNS")
  need(!any(c("GeneID", "gene") %in% names(annotation)),
       "ENS101_AMBIGUOUS_ALREADY_PROJECTED_ATTRIBUTES")
  attributes <- annotation[c("gene_id", "gene_name", "gene_biotype")]
  need(all(vapply(attributes, function(values) {
    is.character(values) && is.null(dim(values)) && length(values) == nrow(annotation)
  }, logical(1))), "ENS101_ATTRIBUTE_TYPES")
  need(gse205906_ens101_character_v1(fixed_gene_ids) &&
         identical(annotation$gene_id, fixed_gene_ids) &&
         identical(nrow(annotation), length(fixed_gene_ids)),
       "ENS101_COMPLETE_GENE_ID_ORDER_IDENTITY")
  if (identical(domain, "FIXED_ENS101")) {
    need(length(fixed_gene_ids) == 32883L &&
           all(grepl("^ENSRNOG[0-9]+$", fixed_gene_ids)) &&
           identical(fixed_gene_ids, sort(fixed_gene_ids, method = "radix")),
         "ENS101_FIXED_32883_STABLE_ORDER_UNIVERSE")
  } else {
    need(length(fixed_gene_ids) >= 1L && length(fixed_gene_ids) <= 100L &&
           all(grepl("^TOY_ENS101_[0-9]+$", fixed_gene_ids)),
         "ENS101_TOY_ONLY_UNIVERSE")
  }

  # Exactly the old helper's trim/NA/empty semantics and threshold. Literal
  # strings such as "NA" are NOT silently recoded to missing. Duplicate symbols
  # remain gene rows: downstream activity mean / fgsea max-|t| policies differ.
  symbols <- trimws(annotation$gene_name)
  nonblank <- !is.na(symbols) & nzchar(symbols)
  coverage <- mean(nonblank)
  need(is.finite(coverage) && coverage >= 0.80,
       "ENS101_ORIGINAL_80_PERCENT_NONBLANK_SYMBOL_RULE")
  resolved_symbols <- symbols
  resolved_symbols[!nonblank] <- fixed_gene_ids[!nonblank]
  raw_annotation <- unserialize(serialize(annotation, NULL, version = 3L))
  computational_interface <- data.frame(
    GeneID = fixed_gene_ids,
    gene = annotation$gene_name,
    gene_biotype = annotation$gene_biotype,
    stringsAsFactors = FALSE
  )
  resolved_annotation <- data.frame(
    gene_id = fixed_gene_ids,
    symbol = resolved_symbols,
    used_gene_id_fallback = !nonblank,
    stringsAsFactors = FALSE
  )
  list(
    schema = "GSE205906_ENS101_ANNOTATION_INTERFACE/v1",
    domain = domain,
    source_reference = "Rnor_6.0_Ensembl101_NOT_RefSeq",
    raw_annotation = raw_annotation,
    computational_interface = computational_interface,
    resolved_annotation = resolved_annotation,
    transformation = c(
      gene_id = "GeneID", gene_name = "gene", gene_biotype = "gene_biotype"
    ),
    metrics = list(
      gene_id_source_column = "gene_id", symbol_source_column = "gene_name",
      biotype_source_column = "gene_biotype", total_gene_rows = length(fixed_gene_ids),
      nonblank_symbol_n = sum(nonblank), nonblank_symbol_coverage = coverage,
      unique_nonblank_symbol_n = length(unique(symbols[nonblank])),
      unique_nonblank_symbol_coverage = length(unique(symbols[nonblank])) / length(fixed_gene_ids),
      fallback_to_gene_id_n = sum(!nonblank), minimum_nonblank_symbol_coverage = 0.80
    ),
    rows_removed = 0L,
    orthology_invented = FALSE,
    actual_public_symbol_coverage_established = FALSE,
    model_run = FALSE
  )
}

gse205906_ens101_pin_v1 <- function(pin) {
  is.list(pin) && identical(sort(names(pin)), sort(c("bytes", "sha256"))) &&
    length(pin$bytes) == 1L && is.numeric(pin$bytes) && !is.logical(pin$bytes) &&
    is.finite(pin$bytes) && pin$bytes >= 0 && pin$bytes == floor(pin$bytes) &&
    pin$bytes <= 2^53 && is.character(pin$sha256) && length(pin$sha256) == 1L &&
    !is.na(pin$sha256) && grepl("^[0-9a-f]{64}$", pin$sha256)
}

gse205906_ens101_prefit_binding_v1 <- function(
    annotation_projection, sample_ids, phase_packets, input_pins) {
  need <- gse205906_ens101_need_v1
  need(is.list(annotation_projection) &&
         identical(annotation_projection$schema, "GSE205906_ENS101_ANNOTATION_INTERFACE/v1") &&
         annotation_projection$domain %in% c("FIXED_ENS101", "SOURCE_TOY"),
       "ENS101_PREFIT_PROJECTION_INTERFACE")
  # Recompute the pure projection from its retained raw source. A fabricated
  # metric/table or changed fallback cannot masquerade as this adapter output.
  checked <- gse205906_ens101_annotation_v1(
    annotation_projection$raw_annotation,
    annotation_projection$raw_annotation$gene_id, annotation_projection$domain
  )
  need(identical(annotation_projection, checked), "ENS101_PREFIT_PROJECTION_DRIFT")
  need(gse205906_ens101_character_v1(sample_ids) && length(sample_ids) == 9L,
       "ENS101_PREFIT_NINE_SAMPLE_IDENTIFIERS")
  if (identical(checked$domain, "FIXED_ENS101")) {
    need(identical(sample_ids, paste0("GSM", 6235292:6235300)),
         "ENS101_PREFIT_FIXED_SAMPLE_ORDER")
  } else {
    need(identical(sample_ids, paste0("TOY_LIB_", 1:9)),
         "ENS101_PREFIT_TOY_SAMPLE_ORDER")
  }
  input_roles <- c("new_counts", "new_annotation", "new_metadata", "new_config",
                   "new_prefit_SOURCE", "new_mapping_SOURCE", "statistical_runtime_lock")
  need(is.list(input_pins) && identical(names(input_pins), input_roles) &&
         all(vapply(input_pins, gse205906_ens101_pin_v1, logical(1))),
       "ENS101_PREFIT_CALLER_INPUT_PIN_SCHEMA")
  phases <- c("NONE", sample_ids)
  need(is.list(phase_packets) && identical(names(phase_packets), phases),
       "ENS101_PREFIT_COMPLETE_FULL_PLUS_NINE_PHASE_ORDER")
  collections <- c("Reactome", "SynGO_sensitivity", "PROGENy", "DoRothEA_AB")
  gene_ids <- checked$resolved_annotation$gene_id
  phase_views <- lapply(phases, function(phase) {
    packet <- phase_packets[[phase]]
    fields <- c("omitted_library", "sample_ids", "retained_gene_ids",
                "prefit_output_pin", "mapping_output_pin", "families")
    need(is.list(packet) && identical(names(packet), fields),
         "ENS101_PREFIT_PHASE_PACKET_SCHEMA")
    selected_samples <- if (identical(phase, "NONE")) sample_ids else sample_ids[sample_ids != phase]
    need(identical(packet$omitted_library, phase) &&
           identical(packet$sample_ids, selected_samples),
         "ENS101_PREFIT_OMISSION_SAMPLE_BINDING")
    retained <- packet$retained_gene_ids
    need(gse205906_ens101_character_v1(retained) && length(retained) > 0L &&
           all(retained %in% gene_ids) &&
           identical(retained, gene_ids[gene_ids %in% retained]),
         "ENS101_PREFIT_NEW_RETAINED_GENE_ORDER")
    need(gse205906_ens101_pin_v1(packet$prefit_output_pin) &&
           gse205906_ens101_pin_v1(packet$mapping_output_pin),
         "ENS101_PREFIT_PHASE_OUTPUT_PIN_SCHEMA")
    need(is.list(packet$families) && identical(names(packet$families), collections),
         "ENS101_PREFIT_COLLECTION_ORDER")
    for (collection in collections) {
      family <- packet$families[[collection]]
      need(is.list(family) && (length(family) == 0L ||
             gse205906_ens101_character_v1(names(family))),
           "ENS101_PREFIT_FAMILY_IDENTIFIERS")
      if (identical(collection, "PROGENy")) {
        need(length(family) <= 14L, "ENS101_PREFIT_PROGENY_SOURCE_DOMAIN_BOUND")
      }
      if (identical(collection, "DoRothEA_AB")) {
        need(length(family) <= 117L, "ENS101_PREFIT_DOROTHEA_SOURCE_DOMAIN_BOUND")
      }
      for (member in family) {
        need(is.list(member) && identical(names(member), c("n", "sha256")) &&
               is.numeric(member$n) && !is.logical(member$n) && length(member$n) == 1L &&
               is.finite(member$n) && member$n == floor(member$n) &&
               member$n >= 10 && member$n <= 500 &&
               gse205906_ens101_pin_v1(list(bytes = 0, sha256 = member$sha256)),
             "ENS101_PREFIT_ELIGIBLE_MEMBER_METADATA")
      }
    }
    list(
      omitted_library = phase,
      sample_ids = selected_samples,
      retained_gene_ids = retained,
      annotation = checked$resolved_annotation[match(retained, gene_ids), , drop = FALSE],
      prefit_output_pin = packet$prefit_output_pin,
      mapping_output_pin = packet$mapping_output_pin,
      families = packet$families,
      eligible_tested_program_n = lengths(packet$families),
      recomputation_required_for_this_phase = TRUE
    )
  })
  names(phase_views) <- phases
  list(
    schema = "GSE205906_ENS101_CALLER_METADATA_PREFIT_BINDING_INTERFACE/v1",
    domain = checked$domain,
    source_reference = checked$source_reference,
    input_pins = input_pins,
    phases = phase_views,
    source_program_domain_upper_bounds = c(PROGENy = 14L, DoRothEA_AB = 117L),
    membership_identity_is_caller_metadata_not_authenticated_here = TRUE,
    effect_independent_prefit_or_filtering_executed_here = FALSE,
    actual_public_symbol_coverage_established = FALSE,
    empty_family_requires_separate_reporting_review_before_fitting = TRUE,
    gene_symbol_rows_merged = FALSE,
    model_run = FALSE
  )
}
