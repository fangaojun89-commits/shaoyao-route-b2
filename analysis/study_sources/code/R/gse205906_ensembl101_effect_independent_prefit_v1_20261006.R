# NEW SOURCE candidate: effect-independent full + nine leave-library preprocessing.
# No file IO, model fit/dispersion/voom/score, resource mapping, deployment or
# authority. Raw arguments and old files remain unchanged. No TOY execution mode.

gse205906_ens101_prefit_need_v1 <- function(condition, category) {
  if (!isTRUE(condition)) stop(category, call. = FALSE)
}

gse205906_ens101_prefit_hash_v1 <- function(object) {
  digest::digest(serialize(object, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_prefit_sources_v1 <- function(source_bundle, cfg) {
  need <- gse205906_ens101_prefit_need_v1
  pins <- list(
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    helpers = list(bytes = 175131L, sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
    annotation_adapter = list(bytes = 9943L, sha256 = "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470")
  )
  need(is.list(source_bundle) && identical(names(source_bundle), names(pins)),
       "ENS101_PREFIT_EXACT_SOURCE_ROLES")
  for (role in names(pins)) {
    bytes <- source_bundle[[role]]
    need(is.raw(bytes) && length(bytes) == pins[[role]]$bytes &&
           identical(digest::digest(bytes, "sha256", serialize = FALSE), pins[[role]]$sha256),
         "ENS101_PREFIT_SOURCE_BYTES_OR_SHA_DRIFT")
  }
  # Authenticate complete SOURCE bytes first, then evaluate ONLY pinned function
  # declarations. The original model entrypoint is NEVER sourced or executed.
  declarations <- function(bytes, allowed, environment) {
    parsed <- parse(text = rawToChar(bytes), keep.source = FALSE)
    selected <- Filter(function(expr) {
      is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% allowed
    }, as.list(parsed))
    names_found <- vapply(selected, function(expr) as.character(expr[[2L]]), character(1))
    need(setequal(names_found, allowed) && !anyDuplicated(names_found),
         "ENS101_PREFIT_SOURCE_FUNCTION_COVERAGE")
    for (expr in selected) {
      need(is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function")),
           "ENS101_PREFIT_SOURCE_NOT_FUNCTION")
      eval(expr, envir = environment)
    }
    invisible(NULL)
  }
  helpers <- new.env(parent = baseenv())
  declarations(source_bundle$helpers,
               c("v2r4_stop", "v2r4_assert", "v2r4_design", "v2r4_validate_counts_original"), helpers)
  adapter <- new.env(parent = baseenv())
  declarations(source_bundle$annotation_adapter,
               c("gse205906_ens101_need_v1", "gse205906_ens101_character_v1",
                 "gse205906_ens101_annotation_v1"), adapter)
  model <- parse(text = rawToChar(source_bundle$model), keep.source = FALSE)
  selected <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("v2r4_fit_preprocess"))
  }, as.list(model))
  need(length(selected) == 1L, "ENS101_PREFIT_REFERENCE_FUNCTION_NOT_UNIQUE")
  declaration <- selected[[1L]][[3L]]
  need(is.call(declaration) && identical(declaration[[1L]], as.name("function")),
       "ENS101_PREFIT_REFERENCE_NOT_FUNCTION")
  statements <- as.list(declaration[[3L]])
  need(length(statements) == 11L &&
         identical(statements[[7L]][[3L]][[1L]], quote(edgeR::normLibSizes)) &&
         identical(statements[[8L]][[3L]][[1L]], quote(edgeR::estimateDisp)),
       "ENS101_PREFIT_REFERENCE_BOUNDARY_DRIFT")
  prefix <- as.call(c(list(as.name("{")), statements[2:7],
                     list(quote(list(keep = keep, keep_filter_by_expr = keep_fbe,
                                     keep_cpm = keep_cpm, y = y)))))
  reference_environment <- new.env(parent = helpers)
  reference_environment$cfg <- cfg
  reference <- eval(as.call(list(as.name("function"), declaration[[2L]], prefix)),
                    envir = reference_environment)
  lockEnvironment(helpers, bindings = TRUE)
  lockEnvironment(adapter, bindings = TRUE)
  lockEnvironment(reference_environment, bindings = TRUE)
  list(helpers = helpers, adapter = adapter, reference = reference,
       provenance = list(source_pins_authenticated_here = pins,
                         prefix_object_sha256 = gse205906_ens101_prefit_hash_v1(prefix),
                         retained_statement_n = 6L,
                         excluded_first_call = "edgeR::estimateDisp",
                         reference_kind = "SAME_PACKAGE_PINNED_SOURCE_PREFIX_NOT_EXTERNAL_VALIDATION"))
}

gse205906_ens101_prefit_config_v1 <- function(cfg) {
  need <- gse205906_ens101_prefit_need_v1
  need(is.list(cfg) && identical(names(cfg),
         c("schema", "reference", "design", "analysis", "input_snapshot")) &&
         identical(cfg$schema, "GSE205906_ENS101_EFFECT_INDEPENDENT_PREFIT_CONFIG/v1"),
       "ENS101_PREFIT_NEW_CONFIG_SCHEMA")
  need(identical(cfg$reference, list(assembly = "Rnor_6.0", annotation = "Ensembl101",
         expected_gene_rows = 32883L, primary_strand_setting = 0L, count_unit = "gene_fragments")),
       "ENS101_PREFIT_NEW_REFERENCE_DOMAIN")
  ids <- paste0("GSM", 6235292:6235300)
  groups <- rep(c("sham", "SNI", "PLM"), each = 3L)
  need(is.list(cfg$design) && identical(names(cfg$design),
         c("sample_ids_in_order", "groups_in_order", "group_levels", "coefficient_order",
           "expected_rank", "contrasts")), "ENS101_PREFIT_DESIGN_CONFIG_SCHEMA")
  need(identical(unlist(cfg$design$sample_ids_in_order, use.names = FALSE), ids) &&
         identical(unlist(cfg$design$groups_in_order, use.names = FALSE), groups) &&
         identical(unlist(cfg$design$group_levels, use.names = FALSE), c("sham", "SNI", "PLM")) &&
         identical(unlist(cfg$design$coefficient_order, use.names = FALSE), c("sham", "SNI", "PLM")) &&
         identical(as.numeric(cfg$design$expected_rank), 3), "ENS101_PREFIT_FIXED_DESIGN_DRIFT")
  need(is.list(cfg$design$contrasts) && identical(names(cfg$design$contrasts),
         c("local_disease", "intervention")), "ENS101_PREFIT_CONTRAST_SCHEMA")
  need(identical(as.numeric(unlist(cfg$design$contrasts$local_disease$vector)), c(-1, 1, 0)) &&
         identical(as.numeric(unlist(cfg$design$contrasts$intervention$vector)), c(0, -1, 1)),
       "ENS101_PREFIT_CONTRAST_DRIFT")
  need(is.list(cfg$analysis) && identical(names(cfg$analysis),
         c("minimum_cpm", "minimum_samples_at_cpm", "minimum_retained_genes_production", "random_seed")) &&
         identical(as.numeric(unlist(cfg$analysis, use.names = FALSE)), c(1, 3, 1000, 20260817)) &&
         !any(vapply(cfg$analysis, is.logical, logical(1))), "ENS101_PREFIT_FIXED_PARAMETERS_DRIFT")
  need(identical(cfg$input_snapshot, list(minimum_nonblank_symbol_coverage = 0.80)),
       "ENS101_PREFIT_ORIGINAL_SYMBOL_THRESHOLD_DRIFT")
  invisible(TRUE)
}

gse205906_ens101_prefit_ledger_v1 <- function(ledger) {
  need <- gse205906_ens101_prefit_need_v1
  fields <- strsplit(paste(
    "label group declared_biological_replicate GEO_sample BioSample SRA_sample SRA_experiment SRA_run",
    "BioProject SRA_study SRA_submission submission_alias library_name experiment_alias run_alias",
    "library_strategy library_source library_selection library_layout instrument total_spots total_bases",
    "run_count sample_members reads_per_spot biosample_submission_date biosample_publication_date",
    "run_published genotype tissue treatment animal_identifier animals_per_library physical_pool_membership",
    "independent_animal_crosswalk aliquot_identifier RNA_extraction_batch library_preparation_batch",
    "sequencing_batch_identifier original_FASTQ_files analysis_unit nominal_group_source dose_metadata_conflict",
    "nominal_PLM_group_dose_metadata_conflict paper_nominal_paeoniflorin_mg_per_ml",
    "paper_nominal_liquiritin_mg_per_ml administered_dose_verified author_reply_required",
    "full_run_QC_completed_in_supplied_packet ENA_generated_FASTQ_MD5 ENA_generated_FASTQ_bytes",
    "ENA_generated_FASTQ_HTTPS overlay_role"
  ), " ", fixed = TRUE)[[1L]]
  need(is.data.frame(ledger) && nrow(ledger) == 9L && identical(names(ledger), fields) &&
         all(vapply(ledger, function(x) is.character(x) && is.null(dim(x)) &&
                      length(x) == 9L && !anyNA(x), logical(1))),
       "ENS101_PREFIT_ORIGINAL_53_COLUMN_LEDGER_SCHEMA")
  expected <- list(
    label = paste0(rep(c("JA", "SNI", "PLM"), each = 3L), rep(1:3, 3L)),
    group = rep(c("JA", "SNI", "PLM"), each = 3L),
    declared_biological_replicate = as.character(rep(1:3, 3L)),
    GEO_sample = paste0("GSM", 6235292:6235300),
    BioSample = paste0("SAMN", 28987251 - 0:8),
    SRA_run = paste0("SRR", 19631207 - 0:8),
    SRA_experiment = paste0("SRX", 15681748 + 0:8),
    SRA_sample = paste0("SRS", 13379586 + 0:8)
  )
  need(all(vapply(names(expected), function(key) identical(ledger[[key]], expected[[key]]), logical(1))),
       "ENS101_PREFIT_FIXED_LIBRARY_ACCESSION_BIJECTION")
  need(all(ledger$library_layout == "PAIRED") && all(ledger$library_strategy == "RNA-Seq") &&
         all(ledger$library_source == "TRANSCRIPTOMIC"), "ENS101_PREFIT_LIBRARY_CLASS_DRIFT")
  expected
}

gse205906_ens101_prefit_phase_v1 <- function(counts, design, source) {
  need <- gse205906_ens101_prefit_need_v1
  if (!all(colSums(counts) > 0)) {
    return(list(status = "HELD_ZERO_INPUT_LIBRARY", retained_gene_n = NA_integer_,
                minimum_retained_genes = 1000L, normalized_counts_available = FALSE))
  }
  observed <- tryCatch({
    y0 <- edgeR::DGEList(counts = counts)
    fbe <- edgeR::filterByExpr(y0, design = design)
    cpm_mask <- rowSums(edgeR::cpm(y0, log = FALSE) >= 1) >= 3L
    list(y0 = y0, keep = fbe & cpm_mask, keep_filter_by_expr = fbe, keep_cpm = cpm_mask)
  }, error = function(e) NULL)
  if (is.null(observed)) {
    return(list(status = "UNVERIFIED_EDGER_FILTER", normalized_counts_available = FALSE))
  }
  mask_data <- observed[c("keep", "keep_filter_by_expr", "keep_cpm")]
  retained_n <- sum(observed$keep)
  masks <- c(mask_data, list(retained_gene_n = retained_n,
                            filter_by_expr_n = sum(observed$keep_filter_by_expr),
                            cpm_n = sum(observed$keep_cpm)))
  # No normalization or permissive threshold when a phase fails the original
  # minimum. Keep its actual mask/failure in the complete full + nine grid.
  if (retained_n < 1000L) {
    reference_failed <- tryCatch({ source$reference(counts, design, 1000L); FALSE },
      error = function(e) identical(conditionMessage(e),
                                   sprintf("TOO_FEW_FILTERED_GENES:%d", retained_n)))
    return(c(list(status = if (reference_failed) "HELD_MINIMUM_RETAINED" else "UNVERIFIED_REFERENCE_MINIMUM",
                  normalized_counts_available = FALSE, minimum_retained_genes = 1000L,
                  same_package_reference_minimum_rejection = reference_failed), masks))
  }
  result <- tryCatch({
    y <- edgeR::normLibSizes(observed$y0[observed$keep, , keep.lib.sizes = FALSE], method = "TMM")
    need(all(is.finite(y$samples$norm.factors)) && all(y$samples$norm.factors > 0) &&
           abs(sum(log(y$samples$norm.factors))) < 1e-10, "ENS101_PREFIT_TMM_INVALID")
    reference <- source$reference(counts, design, 1000L)
    need(all(vapply(names(mask_data), function(key) identical(mask_data[[key]], reference[[key]]), logical(1))) &&
           identical(y$counts, reference$y$counts) && identical(y$samples$lib.size, reference$y$samples$lib.size) &&
           max(abs(y$samples$norm.factors - reference$y$samples$norm.factors)) <= 1e-12,
         "ENS101_PREFIT_PINNED_REFERENCE_DISAGREEMENT")
    list(y = y)
  }, error = function(e) NULL)
  if (is.null(result)) {
    return(c(list(status = "UNVERIFIED_TMM_OR_PINNED_REFERENCE",
                  normalized_counts_available = FALSE), masks))
  }
  c(list(status = "PASS_DESIGN_FILTER_TMM_ONLY", normalized_counts_available = TRUE,
         y = result$y, frozen_prefix_agreement = TRUE, minimum_retained_genes = 1000L), masks)
}

gse205906_ens101_effect_independent_prefit_v1 <- function(inputs, cfg, source_bundle, input_pins) {
  need <- gse205906_ens101_prefit_need_v1
  need(requireNamespace("edgeR", quietly = TRUE) && requireNamespace("digest", quietly = TRUE),
       "ENS101_PREFIT_EXISTING_PACKAGES_REQUIRED_NO_AUTO_INSTALL")
  gse205906_ens101_prefit_config_v1(cfg)
  source <- gse205906_ens101_prefit_sources_v1(source_bundle, cfg)
  need(is.list(inputs) && identical(names(inputs),
         c("count_matrix", "gene_annotation", "sample_ledger", "fixed_gene_ids")),
       "ENS101_PREFIT_NEW_INPUT_ROLES")
  roles <- c("new_counts", "new_annotation", "new_metadata", "fixed_gene_manifest", "new_config",
             "new_prefit_SOURCE", "statistical_runtime_lock")
  pin_ok <- function(pin) is.list(pin) && identical(names(pin), c("bytes", "sha256")) &&
    is.numeric(pin$bytes) && !is.logical(pin$bytes) && length(pin$bytes) == 1L &&
    is.finite(pin$bytes) && pin$bytes >= 0 && pin$bytes <= 2^53 && pin$bytes == floor(pin$bytes) &&
    is.character(pin$sha256) && length(pin$sha256) == 1L && !is.na(pin$sha256) &&
    grepl("^[0-9a-f]{64}$", pin$sha256)
  need(is.list(input_pins) && identical(names(input_pins), roles) &&
         all(vapply(input_pins, pin_ok, logical(1))), "ENS101_PREFIT_DECLARED_INPUT_PIN_SCHEMA")
  before <- gse205906_ens101_prefit_hash_v1(list(inputs = inputs, cfg = cfg, pins = input_pins))
  ids <- paste0("GSM", 6235292:6235300)
  mapping <- gse205906_ens101_prefit_ledger_v1(inputs$sample_ledger)
  need(is.matrix(inputs$count_matrix) && identical(dim(inputs$count_matrix), c(32883L, 9L)) &&
         identical(rownames(inputs$count_matrix), inputs$fixed_gene_ids),
       "ENS101_PREFIT_NEW_32883_BY_NINE_GSM_INPUT")
  counts <- source$helpers$v2r4_validate_counts_original(inputs$count_matrix, ids)
  annotation <- source$adapter$gse205906_ens101_annotation_v1(
    inputs$gene_annotation, inputs$fixed_gene_ids, "FIXED_ENS101")
  design_contract <- source$helpers$v2r4_design(inputs$sample_ledger, cfg)
  groups <- design_contract$group
  phases <- c("NONE", ids)
  runs <- lapply(seq_along(phases), function(phase_index) {
    selected <- if (phase_index == 1L) seq_along(ids) else setdiff(seq_along(ids), phase_index - 1L)
    phase_groups <- factor(as.character(groups[selected]), levels = levels(groups))
    design <- stats::model.matrix(~ 0 + phase_groups)
    colnames(design) <- c("sham", "SNI", "PLM")
    rownames(design) <- ids[selected]
    storage.mode(design) <- "double"
    need(qr(design)$rank == 3L, "ENS101_PREFIT_PHASE_DESIGN_RANK")
    if (phase_index == 1L) need(
      identical(as.vector(design), as.vector(design_contract$design)) &&
        identical(dimnames(design), dimnames(design_contract$design)),
      "ENS101_PREFIT_FULL_DESIGN_DRIFT")
    phase_counts <- counts[, selected, drop = FALSE]
    computed <- gse205906_ens101_prefit_phase_v1(phase_counts, design, source)
    kept_ids <- if (is.null(computed$keep)) character() else rownames(phase_counts)[computed$keep]
    view <- list(
      omitted_library = phases[[phase_index]],
      omission = phases[[phase_index]],
      omitted_group = if (phase_index == 1L) "NONE" else as.character(groups[[phase_index - 1L]]),
      sample_ids = ids[selected], group_sizes = as.list(table(phase_groups)),
      design = design, rank = qr(design)$rank, nominal_residual_df = nrow(design) - qr(design)$rank,
      phase_input_object_sha256 = gse205906_ens101_prefit_hash_v1(phase_counts),
      original_gene_n = nrow(phase_counts), original_library_sizes = setNames(colSums(phase_counts), ids[selected]),
      preprocessing = computed, retained_gene_ids = kept_ids,
      retained_gene_n = computed$retained_gene_n,
      filter_by_expr_n = computed$filter_by_expr_n, cpm_n = computed$cpm_n,
      retained_annotation = annotation$resolved_annotation[match(kept_ids, inputs$fixed_gene_ids), , drop = FALSE],
      mask_object_sha256 = if (is.null(computed$keep)) NA_character_ else gse205906_ens101_prefit_hash_v1(computed$keep),
      membership_computed = FALSE, model_run = FALSE
    )
    if (isTRUE(computed$normalized_counts_available)) {
      view$filtered_library_sizes <- setNames(computed$y$samples$lib.size, ids[selected])
      view$tmm_factors <- setNames(computed$y$samples$norm.factors, ids[selected])
      view$retained_unique_symbol_n <- length(unique(view$retained_annotation$symbol))
      view$duplicate_symbol_excess <- length(kept_ids) - view$retained_unique_symbol_n
    }
    view$phase_output_object_sha256 <- gse205906_ens101_prefit_hash_v1(view)
    view
  })
  names(runs) <- phases
  after <- gse205906_ens101_prefit_hash_v1(list(inputs = inputs, cfg = cfg, pins = input_pins))
  need(identical(before, after), "ENS101_PREFIT_RAW_INPUT_OR_CONFIG_MUTATION")
  pass <- vapply(runs, function(run) identical(run$preprocessing$status, "PASS_DESIGN_FILTER_TMM_ONLY"), logical(1))
  result <- list(
    schema = "GSE205906_ENS101_EFFECT_INDEPENDENT_PREFIT/v1",
    status = if (all(pass)) "PASS_TEN_DESIGN_FILTER_TMM_ONLY" else "NEEDS_REVIEW_PHASE_PREFIT_FAILURES_RETAINED",
    input_reference = cfg$reference, sample_bijection = mapping, declared_input_file_pins = input_pins,
    declared_file_pins_authenticated_here = FALSE,
    input_object_identities = lapply(inputs, gse205906_ens101_prefit_hash_v1),
    config_object_sha256 = gse205906_ens101_prefit_hash_v1(cfg),
    counts_gene_ids = rownames(counts), annotation_projection = annotation,
    sample_ledger_preserved = unserialize(serialize(inputs$sample_ledger, NULL, version = 3L)),
    contrast_matrix = design_contract$contrast, runs = runs,
    passing_phase_n = sum(pass), failed_phase_n = sum(!pass),
    frozen_reference = source$provenance,
    parameters = cfg$analysis,
    runtime_observed = list(R_version = R.version.string,
                            edgeR_version = as.character(utils::packageVersion("edgeR")),
                            digest_version = as.character(utils::packageVersion("digest")),
                            locked_runtime_authenticated_here = FALSE),
    input_before_object_sha256 = before, input_after_object_sha256 = after,
    raw_genes_or_samples_removed = 0L, group_relabel_rule = "SOURCE_JA_TO_SHAM_ONLY",
    animal_pool_batch_metadata_resolved = FALSE, administered_dose_resolved = FALSE,
    membership_computed = FALSE, gene_dispersion_or_effect_fit_called = FALSE,
    program_scoring_or_testing_called = FALSE, full_activity_loao_completed = FALSE,
    production_go = FALSE, new_model_authorized = FALSE,
    random_seed = "NOT_APPLICABLE_DETERMINISTIC_PREFIT_ONLY",
    frozen_future_model_seed = 20260817L
  )
  result$output_object_sha256 <- gse205906_ens101_prefit_hash_v1(result)
  result
}
