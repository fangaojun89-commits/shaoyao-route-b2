# V2 successor: optional day-companion scope; v1 remains unchanged.
# Definitions only. No scientific I/O, source(), RNG reset or invocation.
# The consuming worker must bind these explicit function references and all
# admitted inputs. This module neither reconstructs scores nor grants execution.

g128b_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

g128b_contract <- function() {
  list(
    core_sha256 = c(
      gene_camera = "dcb6b4ba4cdc33093df91bd3efcd55297a44339eb4be273c7aa943558ac6b613",
      activity = "12bb16ccfb85cc147180a721519fac47e414eeadfe1b315d997649599179c12b",
      unsigned = "f2efd9e8ba2160dcf01b66771c12345b8c71596a41937d554a745c2b06f0def6"
    ),
    groups = c("control", "D4", "D7", "D14", "D28", "D42"),
    collections = c("DoRothEA_AB", "PROGENy", "Reactome", "SynGO_assignedBy_proxy_sensitivity"),
    unsigned_seed = 20260819L, unsigned_nrot = 9999L,
    baseline_E_definition = "preserved observed native gene-core E from this same admitted run",
    baseline_E_comparison = "exact ordered numeric equality; no historical-export tolerance"
  )
}

g128b_fit_bundle <- function(counts, metadata, annotation, members, scores, registry,
                              expected_retained_gene_ids, baseline_voom_E, labels,
                              cohort, kernels, annotation_feature_id = "gene_id",
                              gene_symbol_column = "gene_symbol", target_contrast = NULL) {
  contract <- g128b_contract()
  g128b_assert(is.null(target_contrast) ||
                 (is.character(target_contrast) && length(target_contrast) == 1L &&
                    !is.na(target_contrast) && target_contrast %in% c(
                      "D4_vs_control", "D7_vs_control", "D14_vs_control",
                      "D28_vs_control", "D42_vs_control")),
               "Target is NULL for the full bundle or one of five fixed day-vs-control companions")
  g128b_assert(is.list(kernels) &&
                 identical(names(kernels), c("gene_camera", "activity", "unsigned")) &&
                 all(vapply(kernels, is.function, logical(1))),
               "Supply the three explicitly bound kernel function references")
  g128b_assert(is.character(cohort) && length(cohort) == 1L &&
                 !is.na(cohort) && cohort %in% c("all30", "predefined29"), "Unknown cohort")
  expected_n <- if (cohort == "all30") 30L else 29L
  g128b_assert(is.matrix(counts) && is.numeric(counts) && ncol(counts) == expected_n &&
                 !is.null(rownames(counts)) && !is.null(colnames(counts)) &&
                 !anyNA(rownames(counts)) && !anyNA(colnames(counts)) &&
                 all(nzchar(rownames(counts))) && all(nzchar(colnames(counts))) &&
                 !anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)),
               "Fixed-cohort counts require unique nonmissing row/sample keys")
  sample_ids <- colnames(counts)
  g128b_assert(is.data.frame(metadata) &&
                 all(c("sample_id", "group") %in% names(metadata)) &&
                 nrow(metadata) == expected_n &&
                 identical(as.character(metadata$sample_id), sample_ids), "Count/metadata sample order differs")
  g128b_assert(is.character(labels) && !anyNA(labels) &&
                 identical(names(labels), sample_ids), "Labels require unchanged named sample order")
  group_counts <- if (cohort == "all30") rep(5L, 6L) else c(5L, 4L, rep(5L, 4L))
  group <- factor(unname(labels), levels = contract$groups)
  g128b_assert(!anyNA(group) && identical(tabulate(group, nbins = 6L), group_counts),
               "Labels changed frozen cohort capacities")
  g128b_assert(is.character(expected_retained_gene_ids) &&
                 length(expected_retained_gene_ids) > 1L && !anyNA(expected_retained_gene_ids) &&
                 all(nzchar(expected_retained_gene_ids)) && !anyDuplicated(expected_retained_gene_ids) &&
                 all(expected_retained_gene_ids %in% rownames(counts)), "Invalid expected retained IDs")
  g128b_assert(is.matrix(baseline_voom_E) && is.numeric(baseline_voom_E) &&
                 all(is.finite(baseline_voom_E)) &&
                 identical(rownames(baseline_voom_E), expected_retained_gene_ids) &&
                 identical(colnames(baseline_voom_E), sample_ids), "Own baseline E keys/sample order differ")
  g128b_assert(is.matrix(scores) && is.numeric(scores) && all(is.finite(scores)) &&
                 nrow(scores) >= 2L && identical(colnames(scores), sample_ids),
               "Frozen score/sample order differs; do not reorder or approximate")
  g128b_assert(is.data.frame(registry) &&
                 all(c("collection", "program_id", "program_name") %in% names(registry)) &&
                 nrow(registry) == nrow(scores), "Frozen registry shape differs")
  for (field in c("collection", "program_id", "program_name")) {
    g128b_assert(is.character(registry[[field]]) && !anyNA(registry[[field]]) &&
                   all(nzchar(registry[[field]])), paste("Invalid registry field", field))
  }
  program_keys <- paste(registry$collection, registry$program_id, sep = "::")
  g128b_assert(identical(rownames(scores), program_keys) && !anyDuplicated(program_keys) &&
                 setequal(registry$collection, contract$collections), "Frozen score/registry key order differs")
  g128b_assert(all(abs(rowMeans(scores)) <= 1e-8) &&
                 all(abs(apply(scores, 1L, stats::sd) - 1) <= 1e-8),
               "Use the already fixed native standardized score layer without rescaling")
  g128b_assert(is.data.frame(annotation) && nrow(annotation) == nrow(counts) &&
                 all(c(annotation_feature_id, gene_symbol_column) %in% names(annotation)) &&
                 identical(trimws(as.character(annotation[[annotation_feature_id]])), rownames(counts)),
               "Full source annotation/count row order differs")

  # Exactly one gene/CAMERA call. No independent unsigned gene fit is allowed.
  gene_camera <- kernels$gene_camera(
    counts = counts, metadata = metadata, annotation = annotation, members = members,
    labels = labels, annotation_feature_id = annotation_feature_id,
    gene_symbol_column = gene_symbol_column
  )
  actual_counts <- gene_camera$native$dge$counts
  actual_E <- gene_camera$native$voom$E
  g128b_assert(identical(gene_camera$retained_gene_ids, expected_retained_gene_ids) &&
                 identical(rownames(actual_counts), expected_retained_gene_ids) &&
                 identical(colnames(actual_counts), sample_ids), "Current retained gene/sample domain changed")
  g128b_assert(is.matrix(actual_E) && identical(dimnames(actual_E), dimnames(baseline_voom_E)) &&
                 identical(as.numeric(actual_E), as.numeric(baseline_voom_E)),
               "Current native E differs from this run own observed baseline; cached scores cannot be reused")
  retained_annotation <- annotation[match(expected_retained_gene_ids, rownames(counts)), , drop = FALSE]
  edge_table <- as.data.frame(gene_camera$tables$edgeR)
  chronic <- edge_table[edge_table$contrast == "chronic_combined", , drop = FALSE]
  g128b_assert(identical(as.character(chronic$gene_id), expected_retained_gene_ids),
               "Current chronic edgeR table/domain differs")
  activity <- kernels$activity(scores = scores, registry = registry, labels = labels, cohort = cohort)
  g128b_assert(identical(as.numeric(activity$diagnostic$design), as.numeric(gene_camera$design)) &&
                 identical(colnames(activity$diagnostic$design), colnames(gene_camera$design)) &&
                 identical(activity$diagnostic$contrasts, gene_camera$contrasts),
               "Gene and activity designs/contrasts disagree")
  unsigned <- if (is.null(target_contrast)) kernels$unsigned(
    counts = actual_counts, voom = gene_camera$native$voom, design = gene_camera$design,
    contrasts = gene_camera$contrasts, members = members,
    annotation = retained_annotation, edge_chronic = chronic,
    seed = contract$unsigned_seed, nrot = contract$unsigned_nrot, minimum_members = 10L,
    gene_id_column = annotation_feature_id, symbol_column = gene_symbol_column
  ) else list(
    mroast = data.frame(), fgsea = data.frame(), prepared = NULL, rng_calls = list(),
    skip_reason = "NATIVE_UNSIGNED_IS_CHRONIC_ONLY_NOT_DEFINED_FOR_DAY_COMPANION",
    target_contrast = target_contrast
  )
  result <- list(
    gene_camera = gene_camera, activity = activity, unsigned = unsigned,
    integration_diagnostics = list(
      cohort = cohort, sample_ids = sample_ids, labels = labels,
      retained_gene_ids = expected_retained_gene_ids, program_keys = program_keys,
      gene_camera_call_n = 1L, activity_call_n = 1L,
      unsigned_call_n = if (is.null(target_contrast)) 1L else 0L,
      current_retained_ids_match = TRUE, current_native_E_exactly_matches_own_baseline = TRUE,
      baseline_E_maximum_absolute_difference = 0,
      frozen_score_layer_recomputed_or_approximated = FALSE,
      unsigned_received_current_gene_native_objects = is.null(target_contrast),
      separate_unsigned_gene_refit = FALSE,
      core_source_pin_contract = contract$core_sha256,
      historical_export_or_score_compatibility_certified = FALSE,
      real_execution_authorization_claimed = FALSE
    )
  )
  if (!is.null(target_contrast)) {
    result$integration_diagnostics$target_contrast <- target_contrast
    result$integration_diagnostics$target_scope <-
      "DAY_COMPANION_FULL_SEVEN_CONTRAST_GENE_AND_JOINT_ACTIVITY_FITS_UNSIGNED_SKIPPED"
    result$integration_diagnostics$archive_consumer_must_pack_target_contrast <- TRUE
  }
  result
}
