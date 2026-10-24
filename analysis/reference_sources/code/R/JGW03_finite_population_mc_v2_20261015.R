# Functions-only finite-population Monte Carlo uncertainty; no I/O or RNG calls.
# B | M ~ Hypergeometric(T, M, n), sampling distinct uniform nonidentity labels.
# The equal-tail confidence set concerns M, not biological uncertainty or FDR.
# Official phyper documentation defines the probability evaluator. The discrete
# confidence-set inversion and (1 + M)/(T + 1) transformation are project work.
# v2 retains v1 separately and adds a conservative floating equality guard.

jgw03_mc_integer_input <- function(x, name, scalar = FALSE) {
  valid <- (typeof(x) %in% c("integer", "double")) &&
    !is.object(x) && is.null(dim(x)) && length(x) > 0L &&
    all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 0) && all(x <= .Machine$integer.max - 1)
  if (scalar) valid <- valid && length(x) == 1L
  if (!valid) {
    stop(name, " must contain complete, finite, plain nonnegative integers",
         if (scalar) " as one scalar" else "", call. = FALSE)
  }
  as.integer(x)
}

jgw03_mc_confidence_input <- function(confidence) {
  valid <- (typeof(confidence) %in% c("integer", "double")) &&
    !is.object(confidence) && is.null(dim(confidence)) &&
    length(confidence) == 1L && is.finite(confidence) &&
    confidence > 0 && confidence < 1
  if (!valid) stop("confidence must be one plain numeric scalar in (0, 1)",
                   call. = FALSE)
  as.double(confidence)
}

jgw03_mc_invert_one <- function(b, n, population, log_cutoff) {
  if (n == population) return(c(M_lower = b, M_upper = b))
  # Include a tiny neighborhood below the ideal cutoff rather than losing
  # exact discrete equality to phyper/log rounding. This can only enlarge the
  # numerically accepted set; it is not exact-real or arbitrary-error proof.
  log_guard <- 8 * .Machine$double.eps * max(1, abs(log_cutoff))

  # Pr_M(B >= b) is nondecreasing in M. Search for its first accepted M.
  lower <- 0
  upper <- as.double(population)
  if (b == 0L) {
    M_lower <- 0L
  } else {
    while (lower < upper) {
      midpoint <- floor((lower + upper) / 2)
      tail_log <- stats::phyper(b - 1L, midpoint, population - midpoint, n,
                               lower.tail = FALSE, log.p = TRUE)
      if (tail_log >= log_cutoff - log_guard) upper <- midpoint else lower <- midpoint + 1
    }
    M_lower <- as.integer(lower)
  }

  # Pr_M(B <= b) is nonincreasing in M. Search for its last accepted M.
  lower <- 0
  upper <- as.double(population)
  if (b == n) {
    M_upper <- population
  } else {
    while (lower < upper) {
      midpoint <- ceiling((lower + upper) / 2)
      tail_log <- stats::phyper(b, midpoint, population - midpoint, n,
                               lower.tail = TRUE, log.p = TRUE)
      if (tail_log >= log_cutoff - log_guard) lower <- midpoint else upper <- midpoint - 1
    }
    M_upper <- as.integer(lower)
  }
  if (M_lower > M_upper) stop("Empty hypergeometric confidence set", call. = FALSE)
  c(M_lower = M_lower, M_upper = M_upper)
}

# tail_hits may be a vector; draws, population, and confidence must be scalars.
# Any unavailable/invalid count rejects the call: no NA-based denominator change.
# The default geometry is JGW03's 10,000 distinct draws from 369,599 nonidentity
# assignments; membership/uniformity/deterministic-statistic evidence is a caller
# responsibility, not authenticated by these arithmetic functions.
jgw03_finite_mc_interval <- function(tail_hits, draws = 10000L,
                                     nonidentity_population = 369599L,
                                     confidence = 0.95) {
  b <- jgw03_mc_integer_input(tail_hits, "tail_hits")
  n <- jgw03_mc_integer_input(draws, "draws", scalar = TRUE)
  population <- jgw03_mc_integer_input(nonidentity_population,
                                      "nonidentity_population", scalar = TRUE)
  confidence <- jgw03_mc_confidence_input(confidence)
  if (n < 1L || population < 1L || n > population || any(b > n)) {
    stop("Require 1 <= draws <= nonidentity_population and 0 <= tail_hits <= draws",
         call. = FALSE)
  }
  log_cutoff <- log1p(-confidence) - log(2)
  distinct_hits <- unique(b)
  bounds <- t(vapply(distinct_hits, jgw03_mc_invert_one, numeric(2),
                     n = n, population = population, log_cutoff = log_cutoff))
  bounds <- bounds[match(b, distinct_hits), , drop = FALSE]
  intervals <- data.frame(
    nonidentity_tail_hits_b = b,
    nonidentity_draws_n = n,
    nonidentity_population_T = population,
    full_space_including_observed = population + 1L,
    confidence = confidence,
    M_lower = as.integer(bounds[, "M_lower"]),
    M_upper = as.integer(bounds[, "M_upper"]),
    conditional_full_P_lower = (1 + bounds[, "M_lower"]) / (population + 1),
    conditional_full_P_upper = (1 + bounds[, "M_upper"]) / (population + 1),
    SAP_MC_P = (1 + as.double(b)) / (n + 1),
    sampled_nonidentity_tail_fraction = as.double(b) / n,
    interval_kind = if (n == population) "EXACT_CENSUS_POINT" else
      "POINTWISE_FINITE_POPULATION_MC_CONFIDENCE_SET",
    row.names = NULL, check.names = FALSE
  )
  list(
    schema = "jgw03-finite-population-mc-uncertainty-v2",
    intervals = intervals,
    scope = list(
      target = "full conditional P = (1 + M)/(T + 1), observed included once",
      sampling = "uniform distinct nonidentity assignments without replacement",
      tail_criterion = "frozen per-assignment statistic meeting the frozen tail criterion",
      confidence_set = "Pr_M(B>=b)>=alpha/2 and Pr_M(B<=b)>=alpha/2",
      probability_evaluation = "stats::phyper log tails; integer binary inversion",
      numerical_equality_guard_log = 8 * .Machine$double.eps *
        max(1, abs(log_cutoff)),
      guard_interpretation = paste(
        "include log tails within guard below cutoff; conservative numerical",
        "enlargement, not exact-real probability evaluation"
      ),
      biological_confidence_interval = FALSE,
      simultaneous_coverage = FALSE,
      inner_randomness_coverage = FALSE,
      FDR_control = FALSE,
      native_P_q_or_gate_modified = FALSE,
      BH_intervals_computed = FALSE,
      actual_allocation_membership_authenticated_here = FALSE
    )
  )
}

# A complete 20-allocation high-dose/control companion has observed and 19 other
# assignments. It has no outer-MC sampling uncertainty. Do not add another +1.
jgw03_exact20_companion_interval <- function(total_tail_hits, confidence = 0.95) {
  total <- jgw03_mc_integer_input(total_tail_hits, "total_tail_hits")
  if (any(total < 1L | total > 20L)) {
    stop("Complete companion total_tail_hits must be in 1..20, including observed",
         call. = FALSE)
  }
  answer <- jgw03_finite_mc_interval(total - 1L, draws = 19L,
                                    nonidentity_population = 19L,
                                    confidence = confidence)
  answer$intervals$total_tail_hits_including_observed <- total
  answer$intervals$exact_companion_P <- as.double(total) / 20
  answer$scope$companion <- "complete 20 assignments; exact total_tail_hits/20"
  answer
}
