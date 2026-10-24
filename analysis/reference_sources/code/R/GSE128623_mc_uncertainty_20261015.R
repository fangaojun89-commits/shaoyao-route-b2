# Functions only; pointwise outer-reference uncertainty, not a model or gate.
# The original empirical P and BH vector is never changed. This sidecar does
# not authenticate archive coverage, exchangeability, or the unobserved space.

g128mc_summarize <- function(tail, cohort) {
  need <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
  }
  need(is.character(cohort) && length(cohort) == 1L && !is.na(cohort) &&
         cohort %in% c("all30", "pre29"), "Expected all30 or pre29 cohort.")
  keys <- c("method", "collection", "contrast", "feature_id", "family_id")
  fields <- c(keys, "empirical_p", "empirical_BH_q", "hits_including_observed",
              "required_total_n", "missing_null_n", "observed_unavailable",
              "available_total_n", "reference_kind")
  need(is.data.frame(tail) && nrow(tail) > 0L && all(fields %in% names(tail)),
       "Expected a nonempty completed reference-tail table.")
  need(all(vapply(tail[keys], is.character, logical(1L))) &&
         !anyNA(tail[keys]) && !anyDuplicated(tail[keys]) &&
         all(nzchar(trimws(unlist(tail[keys], use.names = FALSE)))),
       "Invalid or duplicated feature keys.")
  need(is.character(tail$reference_kind) && !anyNA(tail$reference_kind) &&
         all(tail$reference_kind == "sampled"),
       "Only sampled six-group references need this outer Monte Carlo sidecar.")
  whole_counts <- function(value, upper) {
    is.numeric(value) && !anyNA(value) && all(is.finite(value)) &&
      all(value >= 0 & value <= upper & value == floor(value))
  }
  need(whole_counts(tail$required_total_n, 10001L) &&
         all(tail$required_total_n == 10001L),
       "The frozen reference requires observed once plus 10000 nonidentity draws.")
  need(whole_counts(tail$missing_null_n, 10000L) &&
         whole_counts(tail$available_total_n, 10001L) &&
         whole_counts(tail$hits_including_observed, 10001L) &&
         is.logical(tail$observed_unavailable) && !anyNA(tail$observed_unavailable),
       "Invalid required-draw, availability, or hit fields.")
  observed_valid <- !tail$observed_unavailable
  available <- 10001L - tail$missing_null_n - as.integer(!observed_valid)
  need(all(tail$available_total_n == available) &&
         all(tail$hits_including_observed <= available) &&
         all(tail$hits_including_observed[observed_valid] >= 1L) &&
         all(tail$hits_including_observed[!observed_valid] == 0L),
       "Observed-hit or availability arithmetic differs; never repair counts.")
  complete <- observed_valid & tail$missing_null_n == 0L
  for (field in c("empirical_p", "empirical_BH_q")) {
    value <- tail[[field]]
    need(is.numeric(value) && all(is.na(value) |
           (is.finite(value) & value >= 0 & value <= 1)) &&
           identical(is.na(value), !complete),
         paste("Invalid probability/completeness semantics:", field))
  }
  expected_p <- tail$hits_including_observed[complete] / 10001L
  need(identical(as.numeric(tail$empirical_p[complete]), as.numeric(expected_p)),
       "Empirical P differs from the frozen plus-one rule.")

  # M is retained as decimal text, not rounded to a double and called exact.
  # For both cohorts N=M-1 > 10^19. The collision union bound is 4.9995e-12;
  # round it UP to 5e-12 before subtracting it from the error budget.
  full_space <- c(all30 = "88832646059788350720", pre29 = "14805441009964725120")
  alpha <- 0.05
  coupling_bound <- 5e-12
  alpha_iid <- alpha - coupling_bound
  b <- rep(NA_real_, nrow(tail))
  b[complete] <- tail$hits_including_observed[complete] - 1L
  unique_b <- sort(unique(b[complete]))
  lower <- upper <- rep(NA_real_, length(unique_b))
  lower[unique_b == 0] <- 0
  selected <- which(unique_b > 0)
  lower[selected] <- stats::qbeta(alpha_iid / 2, unique_b[selected],
                                 10000L - unique_b[selected] + 1)
  upper[unique_b == 10000L] <- 1
  selected <- which(unique_b < 10000L)
  upper[selected] <- stats::qbeta(alpha_iid / 2, unique_b[selected] + 1,
                                 10000L - unique_b[selected], lower.tail = FALSE)
  need(all(is.finite(c(lower, upper))) && all(lower >= 0 & upper <= 1 & lower <= upper),
       "Numerical interval calculation failed.")
  indices <- match(b, unique_b)
  intervals <- tail[keys]
  intervals$nonidentity_hits <- b
  intervals$required_nonidentity_n <- 10000L
  intervals$theta_hat <- b / 10000L
  intervals$theta_lower <- lower[indices]
  intervals$theta_upper <- upper[indices]
  intervals$original_empirical_p <- tail$empirical_p
  intervals$original_empirical_BH_q <- tail$empirical_BH_q
  intervals$empirical_p_resolution <- 1 / 10001L
  intervals$mc_uncertainty_computed <- complete
  intervals$status <- ifelse(complete,
    "POINTWISE_INTERVAL_REQUIRES_TOTAL_INDICATOR_AND_SRSWOR",
    "NOT_COMPUTED_REQUIRED_STATISTIC_UNAVAILABLE")

  list(intervals = intervals, context = list(
    schema_version = "gse128623-outer-mc-uncertainty-v1",
    cohort = cohort, full_allocation_space_decimal = unname(full_space[[cohort]]),
    nonidentity_population_lower_bound_decimal = "10000000000000000000",
    required_nonidentity_n = 10000L, nominal_pointwise_confidence = 1 - alpha,
    binomial_interval_confidence = 1 - alpha_iid,
    collision_union_bound_before_upward_rounding = 4.9995e-12,
    conservative_coupling_error_budget = coupling_bound,
    method = "Clopper-Pearson equal-tailed with conservative SRSWOR coupling budget",
    target = "theta: fraction of all nonidentity assignments meeting the frozen upper-tail criterion",
    assumptions = c(
      "Uniform distinct nonidentity assignments, fixed B and no outcome-adaptive stopping",
      "Tail indicator defined for every assignment; sampled completeness does not prove this",
      "Native procedure, tie tolerance and inner RNG seeds fixed across assignments",
      "Pointwise repeated-reference-sampling coverage, not conditional on selection or no missingness"),
    exclusions = c(
      "Not simultaneous across features, families, cells, cohorts, or selected discoveries",
      "Not an interval for realized plus-one P, BH q, biological effects, or clinical benefit",
      "No inner mroast/fgsea RNG, biological exchangeability, animal/pooling/batch assurance",
      "No full-space exact P interval or acceptance of unobserved invalid statistics"),
    source_tail_unchanged = TRUE, original_p_and_BH_outputs_modified = FALSE,
    archive_coverage_authenticated_here = FALSE,
    cohort_binding_authenticated_here = FALSE,
    global_indicator_totality_verified = FALSE,
    exact_enumeration_applicability = "NOT_APPLICABLE_NO_OUTER_SAMPLING_ERROR",
    real_execution_authority_granted = FALSE))
}
