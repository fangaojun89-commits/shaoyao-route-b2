# Pure native activity kernel. No file IO, driver, RNG or real-data invocation.
# Input scores are the already pinned, standardized recovered program scores;
# never substitute ULM, a per-program/per-collection fit or unmoderated OLS.

g128a_require <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

g128a_specification <- function() {
  groups <- c("control", "D4", "D7", "D14", "D28", "D42")
  contrasts <- cbind(
    chronic_combined = c(-1, 0, 0, 0, 0.5, 0.5),
    D28_vs_control = c(-1, 0, 0, 0, 1, 0),
    D42_vs_control = c(-1, 0, 0, 0, 0, 1),
    acute_combined = c(-1, 0.5, 0.5, 0, 0, 0),
    D4_vs_control = c(-1, 1, 0, 0, 0, 0),
    D7_vs_control = c(-1, 0, 1, 0, 0, 0),
    D14_vs_control = c(-1, 0, 0, 1, 0, 0)
  )
  rownames(contrasts) <- groups
  list(groups = groups, contrasts = contrasts,
       collections = c("DoRothEA_AB", "PROGENy", "Reactome",
                       "SynGO_assignedBy_proxy_sensitivity"))
}

g128a_fit_activity <- function(scores, registry, labels, cohort) {
  g128a_require(requireNamespace("limma", quietly = TRUE), "Existing limma is required.")
  g128a_require(as.character(utils::packageVersion("limma")) == "3.68.4",
                "Native activity kernel requires reviewed limma 3.68.4.")
  specification <- g128a_specification()
  g128a_require(is.character(cohort) && length(cohort) == 1L &&
                  cohort %in% c("all30", "predefined29"), "Unknown fixed cohort.")
  expected_counts <- if (cohort == "all30") rep(5L, 6L) else c(5L, 4L, rep(5L, 4L))
  g128a_require(is.matrix(scores) && is.numeric(scores) && nrow(scores) >= 2L &&
                  ncol(scores) == sum(expected_counts) && all(is.finite(scores)),
                "Scores must be finite numeric program-by-fixed-cohort matrix.")
  g128a_require(!is.null(rownames(scores)) && !is.null(colnames(scores)) &&
                  !anyNA(rownames(scores)) && !anyNA(colnames(scores)) &&
                  all(nzchar(rownames(scores))) && all(nzchar(colnames(scores))) &&
                  !anyDuplicated(rownames(scores)) && !anyDuplicated(colnames(scores)),
                "Unique score row and sample identifiers required.")
  g128a_require(is.data.frame(registry) &&
                  all(c("collection", "program_id", "program_name") %in% names(registry)) &&
                  nrow(registry) == nrow(scores), "Registry schema/size mismatch.")
  registry <- as.data.frame(registry, stringsAsFactors = FALSE)
  for (field in c("collection", "program_id", "program_name")) {
    g128a_require(is.character(registry[[field]]) && !anyNA(registry[[field]]) &&
                    all(nzchar(registry[[field]])), paste("Invalid registry field:", field))
  }
  keys <- paste(registry$collection, registry$program_id, sep = "::")
  g128a_require(identical(rownames(scores), keys) && !anyDuplicated(keys) &&
                  setequal(registry$collection, specification$collections),
                "Fixed four-collection registry key/order mismatch.")
  g128a_require(is.character(labels) && !anyNA(labels) && !anyNA(names(labels)) &&
                  identical(names(labels), colnames(scores)),
                "Labels must be explicitly named in unchanged score-column order.")
  group <- factor(unname(labels), levels = specification$groups)
  g128a_require(!anyNA(group) &&
                  identical(tabulate(group, nbins = 6L), expected_counts),
                "Label reassignment changed fixed cohort group capacities.")
  # These are input invariants, never a new per-label filtering/scaling step.
  g128a_require(all(abs(rowMeans(scores)) <= 1e-8) &&
                  all(abs(apply(scores, 1L, stats::sd) - 1) <= 1e-8),
                "Use the pinned native standardized scores; do not rescale/filter here.")
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- specification$groups
  rownames(design) <- colnames(scores)
  g128a_require(qr(design)$rank == 6L, "Fixed design is rank deficient.")
  # One joint moderation across every source program and all seven contrasts.
  fit <- limma::lmFit(scores, design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, specification$contrasts), robust = TRUE)
  native <- do.call(rbind, lapply(seq_len(ncol(specification$contrasts)), function(index) {
    effect <- fit$coefficients[, index]
    standard_error <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
    critical <- stats::qt(0.975, df = fit$df.total)
    p_value <- fit$p.value[, index]
    data.frame(
      collection = registry$collection, program_id = registry$program_id,
      contrast = colnames(specification$contrasts)[[index]],
      effect = unname(effect), SE = unname(standard_error),
      ci_low = unname(effect - critical * standard_error),
      ci_high = unname(effect + critical * standard_error),
      PValue = unname(p_value),
      BH_q = unname(stats::ave(p_value, registry$collection,
                               FUN = function(values) stats::p.adjust(values, "BH"))),
      tested_program_count = length(effect), method = "weighted_voom_activity_descriptive",
      formal_sequential_inference = FALSE, candidate_generation_allowed = FALSE,
      s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT", stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(native) <- NULL
  signed_statistic <- fit$t
  dimnames(signed_statistic) <- list(keys, colnames(specification$contrasts))
  list(native = native, signed_statistic = signed_statistic,
       absolute_statistic = abs(signed_statistic), native_p = fit$p.value,
       native_q = matrix(native$BH_q, nrow = nrow(scores),
                         dimnames = dimnames(signed_statistic)),
       eligible_statistic = is.finite(signed_statistic),
       diagnostic = list(
         program_keys = keys, sample_ids = colnames(scores), group_counts = expected_counts,
         design = design, contrasts = specification$contrasts,
         residual_df = fit$df.residual, total_df = fit$df.total,
         prior_df = fit$df.prior, posterior_variance = fit$s2.post,
         coefficients = fit$coefficients, unscaled_se = fit$stdev.unscaled,
         robust = TRUE, trend = FALSE, joint_collections = TRUE,
         score_layer_recomputed = FALSE, source_filter_or_TMM_recomputed = FALSE,
         primary_models_replaced = FALSE, empirical_tail_or_q_computed = FALSE
       ))
}
