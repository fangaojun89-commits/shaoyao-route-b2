# Definitions only: assemble already-verified small summaries; no model or I/O.
# Caller must authenticate each upstream receipt/output and its native exit.
g128join_summaries <- function(parts, specification) {
  need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
  need(is.list(parts) && length(parts) == 24L && !is.null(names(parts)) &&
    !anyDuplicated(names(parts)) && nrow(specification) == 24L &&
    !anyDuplicated(specification$key) && setequal(names(parts), specification$key),
    'Exact24JobSummaryDomain')
  fields <- c('method', 'collection', 'contrast', 'family_id', 'n',
    'actual_BH_adjustment_n', 'cell_id', 'reference_space',
    'required_allocation_n', 'empirical_unavailable_n', 'empirical_q_le_05_n',
    'native_fraction_allocations_any_q_le_05')
  ordered <- lapply(seq_len(nrow(specification)), function(i) {
    spec <- specification[i, , drop = FALSE]
    x <- parts[[spec$key]]
    need(is.data.frame(x) && all(fields %in% names(x)) &&
      all(names(x) %in% c(fields, 'saved_output_QA')), 'SavedSummaryColumns')
    if ('saved_output_QA' %in% names(x)) need(all(x$saved_output_QA ==
      'PASS_SAVED_TABLES_AND_PINNED_MANIFEST_BINDING_ONLY'), 'SavedSummaryQAScope')
    need(nrow(x) == spec$family_n && !anyNA(x$family_id) &&
      !anyDuplicated(x$family_id) && all(x$cell_id == spec$cell_id) &&
      all(x$reference_space == spec$reference_space) &&
      all(x$required_allocation_n == spec$required_n), 'SavedSummaryJobBinding')
    for (name in c('n', 'actual_BH_adjustment_n', 'empirical_unavailable_n',
                   'empirical_q_le_05_n')) need(is.numeric(x[[name]]) &&
      all(is.finite(x[[name]]) & x[[name]] >= 0 & x[[name]] == floor(x[[name]])),
      'SavedSummaryIntegerCounts')
    need(all(x$actual_BH_adjustment_n + x$empirical_unavailable_n == x$n) &&
      all(x$empirical_q_le_05_n <= x$actual_BH_adjustment_n) &&
      all(is.finite(x$native_fraction_allocations_any_q_le_05) &
        x$native_fraction_allocations_any_q_le_05 >= 0 &
        x$native_fraction_allocations_any_q_le_05 <= 1), 'SavedSummaryCountAndFractionDomain')
    x[, fields, drop = FALSE]
  })
  answer <- do.call(rbind, ordered)
  rownames(answer) <- NULL
  need(nrow(answer) == 448L && sum(answer$reference_space == 'six_group') == 268L &&
    !anyDuplicated(answer[c('cell_id', 'reference_space', 'family_id')]),
    'Whole448FamilyInventory')
  answer
}
