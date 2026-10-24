# SOURCE candidate: membership/eligibility only, never scores, effects or tests.
# Duplicate-symbol gene rows stay intact. Only SET eligibility uses unique symbols.
# The genuine input/prefit/resource/runtime provenance must be authenticated by
# the caller; object consistency here is not a producer token or model permission.

gse205906_ens101_mapping_need_v2 <- function(ok, category) {
  if (!isTRUE(ok)) stop(category, call. = FALSE)
}

gse205906_ens101_mapping_hash_v2 <- function(object) {
  digest::digest(serialize(object, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_members_hash_v2 <- function(symbols) {
  digest::digest(charToRaw(enc2utf8(paste(symbols, collapse = "\n"))),
                 "sha256", serialize = FALSE)
}

gse205906_ens101_mapping_config_v2 <- function(cfg) {
  need <- gse205906_ens101_mapping_need_v2
  programs <- c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK", "NFkB",
                "PI3K", "p53", "TGFb", "TNFa", "Trail", "VEGF", "WNT")
  need(is.list(cfg) && identical(names(cfg), c("schema", "analysis")) &&
         identical(cfg$schema, "GSE205906_ENS101_MEMBERSHIP_CONFIG/v2"),
       "ENS101_MAPPING_NEW_CONFIG_SCHEMA")
  expected <- list(program_size_min = 10L, program_size_max = 500L,
    orthology_authority = "NCBI_Datasets_REST_v18.35.0_snapshot_2026-08-18_strict_one_to_one",
    program_domains = list(
      PROGENy = list(source_program_n = 14L, exact_program_ids = programs),
      DoRothEA_AB = list(frozen_source_program_upper_bound = 117L)))
  need(identical(cfg$analysis, expected), "ENS101_MAPPING_FIXED_RULES_OR_SOURCE_DOMAIN_DRIFT")
  invisible(TRUE)
}

gse205906_ens101_mapping_references_v2 <- function(source_bundle, cfg, rat_map, resources) {
  need <- gse205906_ens101_mapping_need_v2
  pins <- list(
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    prefit = list(bytes = 18762L, sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"))
  need(is.list(source_bundle) && identical(names(source_bundle), names(pins)),
       "ENS101_MAPPING_EXACT_SOURCE_ROLES")
  for (role in names(pins)) {
    bytes <- source_bundle[[role]]
    need(is.raw(bytes) && length(bytes) == pins[[role]]$bytes &&
           identical(digest::digest(bytes, "sha256", serialize = FALSE), pins[[role]]$sha256),
         "ENS101_MAPPING_SOURCE_BYTES_OR_SHA_DRIFT")
  }
  expressions <- as.list(parse(text = rawToChar(source_bundle$model), keep.source = FALSE))
  declaration <- function(name) {
    found <- Filter(function(expr) is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name(name)), expressions)
    need(length(found) == 1L && is.call(found[[1L]][[3L]]) &&
           identical(found[[1L]][[3L]][[1L]], as.name("function")),
         "ENS101_MAPPING_REFERENCE_FUNCTION")
    found[[1L]]
  }
  env <- new.env(parent = baseenv())
  env$cfg <- cfg
  env$rat_map <- rat_map
  env$gene_sets <- resources
  env$v2r4_assert <- function(ok, format, ...) if (!isTRUE(ok)) stop(sprintf(format, ...), call. = FALSE)
  unsigned <- declaration("v2r4_map_unsigned")
  eval(unsigned, env)
  activity <- declaration("v2r4_activity_domain")[[3L]]
  statements <- as.list(activity[[3L]])
  assigned <- function(expr, name) is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    identical(expr[[2L]], as.name(name))
  endpoints <- which(vapply(statements, function(expr) assigned(expr, "dcov"), logical(1)))
  need(length(endpoints) == 1L && endpoints > 3L && assigned(statements[[2L]], "expr") &&
         identical(statements[[endpoints + 1L]][[1L]], as.name("if")),
       "ENS101_MAPPING_REFERENCE_BOUNDARY")
  # Original coverage construction only: omit collapse/expression, expected
  # eligible-count enforcement, network score generation and all fitting.
  body <- as.call(c(list(as.name("{")), statements[3:endpoints],
    list(quote(list(pcov = pcov, dcov = dcov, prog = prog, doro = doro)))))
  env$cli <- list(execution_mode = "PUBLIC_MEMBERSHIP_DIAGNOSTIC")
  coverage <- eval(as.call(list(as.name("function"), pairlist(), body)), env)
  list(unsigned = function(table, collection, row_symbols) {
      env$retained_annotation <- data.table::data.table(symbol = row_symbols)
      env$v2r4_map_unsigned(table, collection)
    }, activity = function(unique_symbols) {
      env$expr <- matrix(numeric(), nrow = length(unique_symbols), ncol = 0L,
                         dimnames = list(unique_symbols, NULL))
      coverage()
    }, source_pins = pins,
    activity_prefix_sha256 = gse205906_ens101_mapping_hash_v2(body),
    unsigned_function_sha256 = gse205906_ens101_mapping_hash_v2(unsigned),
    excluded_after = "if (enforce_expected); expression collapse, scoring and fitting not evaluated")
}

gse205906_ens101_mapping_sets_v2 <- function(table, program_column, target_column, rat_map = NULL) {
  need <- gse205906_ens101_mapping_need_v2
  need(is.data.frame(table) && all(c(program_column, target_column) %in% names(table)),
       "ENS101_MAPPING_RESOURCE_SCHEMA")
  program <- table[[program_column]]
  target <- table[[target_column]]
  need(is.character(program) && is.character(target) && !anyNA(program) && !anyNA(target) &&
         all(nzchar(program)) && all(nzchar(target)), "ENS101_MAPPING_RESOURCE_IDS")
  if (!is.null(rat_map)) target <- rat_map$rat_symbol[match(target, rat_map$human_symbol)]
  ids <- sort(unique(program), method = "radix")
  sets <- lapply(ids, function(id) sort(unique(target[program == id & !is.na(target)]), method = "radix"))
  names(sets) <- ids
  sets
}

gse205906_ens101_mapping_coverage_v2 <- function(sets, unique_symbols, collection, minimum, maximum) {
  need <- gse205906_ens101_mapping_need_v2
  need(is.character(unique_symbols) && !anyNA(unique_symbols) && !anyDuplicated(unique_symbols) &&
         all(nzchar(unique_symbols)) && identical(unique_symbols, sort(unique_symbols, method = "radix")),
       "ENS101_MAPPING_ELIGIBILITY_SET_INVALID")
  # Always retain a typed zero-row schema; absence is not an effect or PASS fit.
  empty <- data.table::data.table(collection = character(), program_id = character(),
    source_feature_n = integer(), tested_feature_n = integer(), eligible = logical(),
    retained_members_sha256 = character())
  if (!length(sets)) return(empty)
  data.table::rbindlist(lapply(names(sets), function(id) {
    matched <- sort(intersect(sets[[id]], unique_symbols), method = "radix")
    data.table::data.table(collection = collection, program_id = id,
      source_feature_n = length(sets[[id]]), tested_feature_n = length(matched),
      eligible = length(matched) >= minimum && length(matched) <= maximum,
      retained_members_sha256 = gse205906_ens101_members_hash_v2(matched))
  }))
}

gse205906_ens101_mapping_v2 <- function(prefit, resources, orthology, cfg, source_bundle, resource_pins) {
  need <- gse205906_ens101_mapping_need_v2
  hash <- gse205906_ens101_mapping_hash_v2
  need(requireNamespace("data.table", quietly = TRUE) && requireNamespace("digest", quietly = TRUE),
       "ENS101_MAPPING_EXISTING_PACKAGES_REQUIRED_NO_INSTALL")
  gse205906_ens101_mapping_config_v2(cfg)
  input <- list(prefit = prefit, resources = resources, orthology = orthology, cfg = cfg,
                resource_pins = resource_pins)
  before <- hash(input)
  need(is.list(resources) && identical(names(resources),
         c("reactome", "syngo_sensitivity", "progeny_human", "dorothea_rat")),
       "ENS101_MAPPING_EXACT_RESOURCE_ROLES")
  pin_ok <- function(pin) is.list(pin) && identical(names(pin), c("bytes", "sha256")) &&
    is.numeric(pin$bytes) && !is.logical(pin$bytes) && length(pin$bytes) == 1L &&
    is.finite(pin$bytes) && pin$bytes >= 0 && pin$bytes <= 2^53 && pin$bytes == floor(pin$bytes) &&
    is.character(pin$sha256) && length(pin$sha256) == 1L && !is.na(pin$sha256) &&
    grepl("^[0-9a-f]{64}$", pin$sha256)
  need(is.list(resource_pins) && identical(names(resource_pins), c(names(resources), "orthology")) &&
         all(vapply(resource_pins, pin_ok, logical(1))), "ENS101_MAPPING_DECLARED_RESOURCE_PINS")
  need(is.data.frame(orthology) && all(c("group_method", "other_tax_id", "human_symbol", "other_symbol") %in% names(orthology)),
       "ENS101_MAPPING_ORTHOLOGY_SCHEMA")
  orth <- data.table::as.data.table(orthology)
  rat_map <- unique(orth[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
    .(human_symbol = as.character(human_symbol), rat_symbol = as.character(other_symbol))])
  need(!anyNA(rat_map) && all(nzchar(rat_map$human_symbol)) && all(nzchar(rat_map$rat_symbol)) &&
         !anyDuplicated(rat_map$human_symbol) && !anyDuplicated(rat_map$rat_symbol),
       "ENS101_MAPPING_NOT_STRICT_ONE_TO_ONE")
  sets <- list(
    Reactome = gse205906_ens101_mapping_sets_v2(resources$reactome, "term_id", "gene_symbol", rat_map),
    SynGO_sensitivity = gse205906_ens101_mapping_sets_v2(resources$syngo_sensitivity, "term_id", "gene_symbol", rat_map),
    PROGENy = gse205906_ens101_mapping_sets_v2(resources$progeny_human, "pathway", "gene_symbol", rat_map),
    DoRothEA_AB = gse205906_ens101_mapping_sets_v2(resources$dorothea_rat, "tf", "target"))
  need(all(resources$dorothea_rat$confidence %in% c("A", "B")), "ENS101_MAPPING_DOROTHEA_EXACT_AB")
  ref <- gse205906_ens101_mapping_references_v2(source_bundle, cfg, rat_map, resources)
  ids <- paste0("GSM", 6235292:6235300)
  phases <- c("NONE", ids)
  need(is.list(prefit) && identical(prefit$schema, "GSE205906_ENS101_EFFECT_INDEPENDENT_PREFIT/v1") &&
         identical(names(prefit$runs), phases) && is.character(prefit$counts_gene_ids) &&
         length(prefit$counts_gene_ids) == 32883L && !anyDuplicated(prefit$counts_gene_ids) &&
         all(grepl("^ENSRNOG[0-9]{11}$", prefit$counts_gene_ids)) &&
         identical(prefit$counts_gene_ids, sort(prefit$counts_gene_ids, method = "radix")),
       "ENS101_MAPPING_NEW_PREFIT_OR_FIXED_GENE_ORDER")
  need(identical(prefit$input_reference, list(assembly = "Rnor_6.0", annotation = "Ensembl101",
         expected_gene_rows = 32883L, primary_strand_setting = 0L, count_unit = "gene_fragments")),
       "ENS101_MAPPING_PREFIT_REFERENCE_DOMAIN")
  original_prefit <- prefit
  claimed <- original_prefit$output_object_sha256
  original_prefit$output_object_sha256 <- NULL
  need(identical(claimed, hash(original_prefit)), "ENS101_MAPPING_PREFIT_OBJECT_DRIFT")
  annotation <- prefit$annotation_projection$resolved_annotation
  need(is.data.frame(annotation) && identical(annotation$gene_id, prefit$counts_gene_ids) &&
         is.character(annotation$symbol) && !anyNA(annotation$symbol) && all(nzchar(annotation$symbol)),
       "ENS101_MAPPING_ANNOTATION_ROW_ASSOCIATION")
  compare <- function(observed, reference, key) {
    a <- observed[source_feature_n > 0, .(program_id, source_feature_n, tested_feature_n, eligible)]
    b <- data.table::data.table(program_id = as.character(reference[[key]]),
      source_feature_n = reference$source_feature_n, tested_feature_n = reference$tested_feature_n,
      eligible = reference$tested)
    data.table::setorder(a, program_id)
    data.table::setorder(b, program_id)
    need(identical(as.data.frame(a), as.data.frame(b)), "ENS101_MAPPING_PINNED_COVERAGE_DISAGREEMENT")
  }
  allowed_status <- c("PASS_DESIGN_FILTER_TMM_ONLY", "HELD_ZERO_INPUT_LIBRARY", "HELD_MINIMUM_RETAINED",
    "UNVERIFIED_EDGER_FILTER", "UNVERIFIED_REFERENCE_MINIMUM", "UNVERIFIED_TMM_OR_PINNED_REFERENCE")
  runs <- lapply(seq_along(phases), function(i) {
    run <- prefit$runs[[i]]
    selected_samples <- if (i == 1L) ids else ids[-(i - 1L)]
    need(is.list(run) && identical(run$omission, phases[[i]]) && identical(run$omitted_library, phases[[i]]) &&
           identical(run$sample_ids, selected_samples), "ENS101_MAPPING_PHASE_ASSOCIATION")
    claimed_phase <- run$phase_output_object_sha256
    phase_body <- run
    phase_body$phase_output_object_sha256 <- NULL
    need(identical(claimed_phase, hash(phase_body)), "ENS101_MAPPING_PHASE_OBJECT_DRIFT")
    stage_status <- run$preprocessing$status
    need(is.character(stage_status) && length(stage_status) == 1L && stage_status %in% allowed_status,
         "ENS101_MAPPING_PREFIT_STATUS_UNKNOWN")
    if (!identical(stage_status, "PASS_DESIGN_FILTER_TMM_ONLY")) {
      return(list(omission = phases[[i]], sample_ids = selected_samples,
        status = "NEEDS_REVIEW_PREFIT_PHASE_NOT_MAPPED", prefit_status = stage_status,
        gene_universe_sha256 = hash(run$retained_gene_ids), retained_gene_n = length(run$retained_gene_ids),
        coverage = NULL, families = NULL, eligible_program_n = NULL,
        failed_phase_not_substituted_by_empty_family = TRUE))
    }
    keep <- run$preprocessing$keep
    need(is.logical(keep) && length(keep) == 32883L && !anyNA(keep) &&
           isTRUE(run$preprocessing$normalized_counts_available) &&
           sum(keep) >= 1000L && is.numeric(run$retained_gene_n) && !is.logical(run$retained_gene_n) &&
           length(run$retained_gene_n) == 1L && is.finite(run$retained_gene_n) && run$retained_gene_n == sum(keep) &&
           is.numeric(run$preprocessing$retained_gene_n) && !is.logical(run$preprocessing$retained_gene_n) &&
           length(run$preprocessing$retained_gene_n) == 1L && is.finite(run$preprocessing$retained_gene_n) &&
           run$preprocessing$retained_gene_n == sum(keep) &&
           identical(run$retained_gene_ids, prefit$counts_gene_ids[keep]) &&
           identical(run$retained_annotation, annotation[keep, , drop = FALSE]) &&
           identical(run$mask_object_sha256, hash(keep)), "ENS101_MAPPING_ACTUAL_PHASE_MASK_OR_ANNOTATION")
    row_symbols <- run$retained_annotation$symbol
    unique_symbols <- sort(unique(row_symbols), method = "radix")
    coverages <- lapply(names(sets), function(collection) gse205906_ens101_mapping_coverage_v2(
      sets[[collection]], unique_symbols, collection, cfg$analysis$program_size_min, cfg$analysis$program_size_max))
    names(coverages) <- names(sets)
    families <- lapply(names(sets), function(collection) {
      selected <- coverages[[collection]][eligible == TRUE]$program_id
      members <- lapply(selected, function(id) sort(intersect(sets[[collection]][[id]], unique_symbols), method = "radix"))
      names(members) <- selected
      members
    })
    names(families) <- names(sets)
    unsigned_indices <- list()
    for (collection in c("Reactome", "SynGO_sensitivity")) {
      table <- if (collection == "Reactome") resources$reactome else resources$syngo_sensitivity
      reference <- ref$unsigned(table, collection, row_symbols)
      compare(coverages[[collection]], reference$coverage, "program_id")
      need(setequal(names(families[[collection]]), names(reference$members)), "ENS101_MAPPING_UNSIGNED_IDS")
      indices <- lapply(families[[collection]], function(members) which(row_symbols %in% members))
      for (id in names(indices)) need(setequal(families[[collection]][[id]], reference$members[[id]]) &&
        identical(indices[[id]], reference$indices[[id]]), "ENS101_MAPPING_UNSIGNED_ROWS_OR_MEMBERS")
      unsigned_indices[[collection]] <- indices
    }
    activity <- ref$activity(unique_symbols)
    compare(coverages$PROGENy, activity$pcov, "source")
    compare(coverages$DoRothEA_AB, activity$dcov, "source")
    family_n <- as.list(lengths(families))
    list(omission = phases[[i]], sample_ids = selected_samples, status = "PASS_MEMBERSHIP_ONLY",
      prefit_phase_object_sha256 = claimed_phase, gene_universe_sha256 = hash(run$retained_gene_ids),
      retained_gene_ids = run$retained_gene_ids, retained_gene_n = length(run$retained_gene_ids),
      retained_symbols = row_symbols, eligibility_symbols = unique_symbols,
      duplicate_symbol_excess = length(row_symbols) - length(unique_symbols),
      coverage = data.table::rbindlist(coverages), families = families,
      family_member_pins = lapply(families, function(family) lapply(family, function(members)
        list(matched_target_n = length(members), matched_targets_sha256 = gse205906_ens101_members_hash_v2(members)))),
      unsigned_gene_row_indices = unsigned_indices, eligible_program_n = family_n,
      empty_family_reporting = as.list(stats::setNames(vapply(families, function(family)
        if (length(family)) "NONEMPTY_MEMBERSHIP_NOT_AN_EFFECT" else "EMPTY_MEMBERSHIP_NO_TESTS_OR_ROBUSTNESS_CLAIM", character(1)), names(families))),
      frozen_coverage_agreement = TRUE, activity_scores_created = FALSE, effect_estimability_proven = FALSE)
  })
  names(runs) <- phases
  runs <- lapply(runs, function(run) {
    run$mapping_phase_object_sha256 <- hash(run)
    run
  })
  pass <- vapply(runs, function(run) identical(run$status, "PASS_MEMBERSHIP_ONLY"), logical(1))
  need(identical(before, hash(input)), "ENS101_MAPPING_INPUT_MUTATION")
  result <- list(schema = "GSE205906_ENS101_MEMBERSHIP/v2",
    status = if (all(pass)) "PASS_TEN_MEMBERSHIP_ONLY" else "NEEDS_REVIEW_PHASE_FAILURES_RETAINED",
    runs = runs, source_program_n = as.list(lengths(sets)), full_eligible_program_n = runs$NONE$eligible_program_n,
    source_zero_rat_mapping_ids = lapply(sets, function(set) names(set)[lengths(set) == 0L]),
    exact_phase_n = 10L, exact_collection_n = 4L, passing_phase_n = sum(pass), failed_phase_n = sum(!pass),
    declared_resource_file_pins = resource_pins, declared_file_pins_authenticated_here = FALSE,
    prefit_producer_provenance_authenticated_here = FALSE, prefit_output_object_sha256 = prefit$output_object_sha256,
    resource_object_sha256 = lapply(resources, hash), resource_bundle_object_sha256 = hash(resources),
    orthology_object_sha256 = hash(orthology),
    reference = list(source_pins = ref$source_pins, activity_prefix_sha256 = ref$activity_prefix_sha256,
      unsigned_function_sha256 = ref$unsigned_function_sha256, excluded_after = ref$excluded_after),
    input_before_object_sha256 = before, input_after_object_sha256 = hash(input),
    raw_gene_rows_removed = 0L, count_rows_collapsed_by_symbol = FALSE,
    full_eligible_cardinalities_are_observed_not_old_config_values = TRUE,
    program_scoring_or_testing_called = FALSE, model_fit_called = FALSE,
    full_activity_loao_completed = FALSE, production_go = FALSE, new_model_authorized = FALSE,
    empty_family_fit_reporting_and_actual_runtime_acceptance_still_required = TRUE,
    runtime_observed = list(R = R.version.string, data_table = as.character(utils::packageVersion("data.table")),
      digest = as.character(utils::packageVersion("digest")), locked_runtime_authenticated_here = FALSE))
  result$output_object_sha256 <- hash(result)
  result
}
