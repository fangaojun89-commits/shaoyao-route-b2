# Allocation preparation only: no expression/count/model object is accepted.
# Usage: selftest FRESH_R_OUTPUT
#        prepare PROJECT_ROOT PINNED_LEDGER FRESH_R_OUTPUT
options(stringsAsFactors = FALSE, warn = 1)

assert_registry <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

allocation_keys <- function(labels) {
  apply(labels, 1L, paste, collapse = ",")
}

validate_allocations <- function(labels, observed, expected_rows) {
  assert_registry(is.matrix(labels) && is.integer(labels), "Integer label matrix required.")
  assert_registry(identical(dim(labels), c(as.integer(expected_rows), length(observed))),
                  "Allocation dimensions differ from contract.")
  assert_registry(!anyNA(labels) && all(labels %in% observed), "Invalid label code.")
  expected_counts <- tabulate(observed, nbins = max(observed))
  same_counts <- apply(labels, 1L, function(row) {
    identical(tabulate(row, nbins = max(observed)), expected_counts)
  })
  assert_registry(all(same_counts), "Allocation changed a group capacity.")
  keys <- allocation_keys(labels)
  assert_registry(!anyDuplicated(keys), "Duplicate label assignment.")
  assert_registry(identical(as.integer(labels[1L, ]), observed), "Observed assignment is not first.")
  assert_registry(sum(keys == paste(observed, collapse = ",")) == 1L,
                  "Observed assignment is not represented exactly once.")
  invisible(TRUE)
}

draw_allocations <- function(observed, draws, seed) {
  assert_registry(is.integer(observed) && !anyNA(observed) && all(observed > 0L),
                  "Positive integer observed labels required.")
  assert_registry(length(draws) == 1L && is.finite(draws) && draws >= 1L &&
                    draws == as.integer(draws), "Invalid draw count.")
  log_space <- lfactorial(length(observed)) - sum(lfactorial(tabulate(observed)))
  assert_registry(log(draws + 1) <= log_space + 1e-12,
                  "Requested more distinct nonidentity assignments than exist.")
  # Uniform index permutations induce uniform assignments of this fixed multiset.
  # Rejecting observed/seen assignments samples distinct nonidentity assignments.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  initial_rng_state <- .Random.seed
  labels <- matrix(NA_integer_, nrow = draws + 1L, ncol = length(observed))
  labels[1L, ] <- observed
  seen <- new.env(hash = TRUE, parent = emptyenv(), size = 2L * (draws + 1L))
  seen[[paste(observed, collapse = ",")]] <- TRUE
  accepted <- 0L
  attempts <- 0L
  attempt_limit <- max(10000L, draws * 100L)
  while (accepted < draws) {
    attempts <- attempts + 1L
    assert_registry(attempts <= attempt_limit, "Distinct-allocation sampling limit reached.")
    proposed <- observed[sample.int(length(observed), length(observed), replace = FALSE)]
    key <- paste(proposed, collapse = ",")
    if (!exists(key, envir = seen, inherits = FALSE)) {
      accepted <- accepted + 1L
      labels[accepted + 1L, ] <- proposed
      seen[[key]] <- TRUE
    }
  }
  validate_allocations(labels, observed, draws + 1L)
  list(labels = labels, seed = seed, initial_rng_state = initial_rng_state,
       final_rng_state = .Random.seed, attempts = attempts,
       rejected_observed_or_duplicate = attempts - draws)
}

enumerate_companion <- function(observed, control_code, day_code) {
  positions <- which(observed %in% c(control_code, day_code))
  assert_registry(control_code != day_code &&
                    all(c(control_code, day_code) %in% observed), "Invalid companion groups.")
  n_control <- sum(observed == control_code)
  control_subsets <- combn(positions, n_control)
  labels <- matrix(rep(observed, each = ncol(control_subsets)),
                   nrow = ncol(control_subsets), ncol = length(observed))
  labels[, positions] <- day_code
  for (row in seq_len(nrow(labels))) labels[row, control_subsets[, row]] <- control_code
  identity <- which(allocation_keys(labels) == paste(observed, collapse = ","))
  assert_registry(length(identity) == 1L, "Companion identity not unique.")
  labels <- labels[c(identity, setdiff(seq_len(nrow(labels)), identity)), , drop = FALSE]
  validate_allocations(labels, observed, choose(length(positions), n_control))
  fixed <- setdiff(seq_along(observed), positions)
  assert_registry(all(vapply(fixed, function(column) {
    all(labels[, column] == observed[column])
  }, logical(1L))), "Companion changed a nuisance-group label.")
  labels
}

read_registry_ledger <- function(path) {
  ledger <- read.delim(path, check.names = FALSE, colClasses = "character", quote = "",
                       comment.char = "", na.strings = character())
  expected <- c("sample_order", "sample_id", "biological_unit_id", "biosample",
                "sra_experiment", "run_accession", "sex", "group")
  assert_registry(identical(names(ledger), expected) && nrow(ledger) == 30L,
                  "Ledger schema or row count changed.")
  assert_registry(!anyNA(ledger) && all(nzchar(as.matrix(ledger))), "Missing ledger field.")
  assert_registry(identical(ledger$sample_order, as.character(seq_len(30L))),
                  "Original ledger order changed.")
  for (field in c("sample_id", "biosample", "sra_experiment", "run_accession")) {
    assert_registry(!anyDuplicated(ledger[[field]]), paste("Nonunique ledger field:", field))
  }
  assert_registry(identical(ledger$sample_id, ledger$biological_unit_id) &&
                    all(ledger$sex == "female"), "Ledger unit/sex differs from fixed metadata.")
  assert_registry(sum(ledger$sample_id == "GSM3681230" &
                        ledger$run_accession == "SRR8757507" & ledger$group == "D4") == 1L,
                  "Predefined exclusion identity differs.")
  ledger
}

create_registry_output <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  assert_registry(grepl("^[Rr]:/Shaoyao_storage/cache/[A-Za-z0-9_-]+$", normalized),
                  "Use a fresh, single named R:/Shaoyao_storage/cache directory.")
  assert_registry(!file.exists(normalized) && dir.exists(dirname(normalized)),
                  "Output exists or its explicit parent is missing; no overwrite/resume.")
  assert_registry(dir.create(normalized), "Could not create output directory.")
  normalized
}

write_registry_table <- function(labels, sample_ids, group_levels, path) {
  table_out <- data.frame(allocation_id = seq_len(nrow(labels)),
                          is_observed = seq_len(nrow(labels)) == 1L,
                          check.names = FALSE)
  label_text <- matrix(group_levels[labels], nrow = nrow(labels), ncol = ncol(labels))
  colnames(label_text) <- sample_ids
  table_out <- cbind(table_out, as.data.frame(label_text, check.names = FALSE))
  connection <- gzfile(path, open = "wt", encoding = "UTF-8", compression = 6L)
  tryCatch(write.table(table_out, connection, sep = "\t", quote = FALSE, row.names = FALSE,
                       na = "NA", eol = "\n"), finally = close(connection))
  # Read through the public file format, not an in-memory alias.
  read_connection <- gzfile(path, "rt")
  restored <- tryCatch(
    read.delim(read_connection, check.names = FALSE, colClasses = "character",
               quote = "", comment.char = ""),
    finally = close(read_connection)
  )
  assert_registry(identical(names(restored), names(table_out)) &&
                    identical(as.matrix(restored[, sample_ids, drop = FALSE]), label_text) &&
                    identical(restored$allocation_id, as.character(seq_len(nrow(labels)))) &&
                    identical(restored$is_observed, as.character(table_out$is_observed)),
                  "Allocation TSV transport mismatch.")
  list(file = basename(path), rows = nrow(labels), libraries = ncol(labels),
       sha256 = digest::digest(file = path, algo = "sha256"),
       assignment_order_sha256 = digest::digest(
         paste(c(paste(sample_ids, collapse = ","), allocation_keys(labels)), collapse = "\n"),
         algo = "sha256", serialize = FALSE))
}

registry_selftest <- function(output) {
  checks <- character()
  check <- function(name, condition) {
    assert_registry(condition, paste("Synthetic check failed:", name))
    checks <<- c(checks, name)
  }
  expect_error <- function(expression) inherits(try(force(expression), silent = TRUE), "try-error")
  toy <- as.integer(c(1, 1, 2, 2, 3, 3))
  draw <- draw_allocations(toy, 89L, 20260817L)
  check("finite toy space exhausts all 90 assignments", nrow(draw$labels) == 90L)
  check("observed first exactly once", sum(allocation_keys(draw$labels) == paste(toy, collapse = ",")) == 1L)
  check("sampling rejects duplicates", draw$attempts > 89L)
  check("seeded replay exact", identical(draw, draw_allocations(toy, 89L, 20260817L)))
  check("seed changes draw order", !identical(draw$labels, draw_allocations(toy, 89L, 20260818L)$labels))
  check("excess distinct draw request rejected", expect_error(draw_allocations(toy, 90L, 1L)))
  check("missing labels rejected", expect_error(draw_allocations(c(toy, NA_integer_), 1L, 1L)))
  check("noninteger labels rejected", expect_error(draw_allocations(as.numeric(toy), 1L, 1L)))
  check("zero draws rejected", expect_error(draw_allocations(toy, 0L, 1L)))
  companion <- enumerate_companion(toy, 1L, 2L)
  check("two plus two complete six assignments", nrow(companion) == 6L)
  check("nuisance labels unchanged", all(companion[, 5:6] == 3L))
  check("same companion group rejected", expect_error(enumerate_companion(toy, 1L, 1L)))
  check("absent companion group rejected", expect_error(enumerate_companion(toy, 1L, 4L)))
  asymmetric <- as.integer(c(rep(1, 5), rep(2, 4), rep(3, 5)))
  asymmetric_companion <- enumerate_companion(asymmetric, 1L, 2L)
  check("five plus four complete 126 assignments", nrow(asymmetric_companion) == 126L)
  balanced <- as.integer(c(rep(1, 5), rep(2, 5), rep(3, 5)))
  check("five plus five complete 252 assignments",
        nrow(enumerate_companion(balanced, 1L, 2L)) == 252L)
  damaged <- companion
  damaged[2L, ] <- damaged[1L, ]
  check("duplicate assignment rejected", expect_error(validate_allocations(damaged, toy, 6L)))
  damaged <- companion
  damaged[2L, 1L] <- 3L
  check("group capacity corruption rejected", expect_error(validate_allocations(damaged, toy, 6L)))
  check("wrong row count rejected", expect_error(validate_allocations(companion, toy, 5L)))
  pin <- write_registry_table(companion, paste0("synthetic_", seq_along(toy)),
                              c("control", "D4", "D7"), file.path(output, "toy.tsv.gz"))
  check("TSV codec verified", pin$rows == 6L && file.exists(file.path(output, pin$file)))
  check("no-overwrite directory rejected", expect_error(create_registry_output(output)))
  check("C drive output rejected", expect_error(create_registry_output("C:/not_allowed")))
  list(status = "SYNTHETIC_ALLOCATION_TESTS_PASS", checks = checks, count = length(checks),
       scientific_data_read = FALSE, models_fitted = 0L,
       toy_content_sha256 = pin$assignment_order_sha256)
}

prepare_registry <- function(project_root, ledger_path, output, script_path) {
  group_levels <- c("control", "D4", "D7", "D14", "D28", "D42")
  inputs <- c(ledger = ledger_path,
              config = file.path(project_root, "config", "contrasts.yaml"),
              sap = file.path(project_root, "protocol", "statistical_analysis_plan.md"),
              definition = file.path(project_root, "protocol", "GSE128623_LABEL_ALLOCATION_REGISTRY_20261015.md"),
              source = script_path)
  expected_pins <- c(
    ledger = "2f1463530d91ee9029c2c7ab5cba92b7c6eeb2ca74cb50767504012387a829c1",
    config = "98ffd869fb3703e38e3318f18887ad746d78a54c59e5a36e5181648f02b63492",
    sap = "c46c086cbba963436fc48cc88f77028713cb305f9e9bff327314fe6eb24e0b82"
  )
  pins <- vapply(inputs, function(path) digest::digest(file = path, algo = "sha256"), character(1L))
  assert_registry(identical(unname(pins[names(expected_pins)]), unname(expected_pins)),
                  "Pinned ledger/config/SAP changed.")
  ledger <- read_registry_ledger(ledger_path)
  assert_registry(identical(as.integer(table(factor(ledger$group, levels = group_levels))), rep(5L, 6L)),
                  "All30 group capacities differ.")
  cohorts <- list(all30 = ledger,
                  predefined29 = ledger[ledger$sample_id != "GSM3681230", , drop = FALSE])
  seeds <- c(all30 = 20260817L, predefined29 = 20260818L)
  exact_space_decimal <- c(all30 = "88832646059788350720", predefined29 = "14805441009964725120")
  receipts <- list()
  rng_states <- list()
  for (cohort_name in names(cohorts)) {
    cohort <- cohorts[[cohort_name]]
    observed <- as.integer(match(cohort$group, group_levels))
    expected_counts <- if (cohort_name == "all30") rep(5L, 6L) else c(5L, 4L, rep(5L, 4L))
    assert_registry(identical(tabulate(observed, nbins = 6L), expected_counts),
                    "Cohort group capacity differs.")
    metadata_file <- file.path(output, paste0(cohort_name, "_samples.tsv"))
    write.table(cohort, metadata_file, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
    sampled <- draw_allocations(observed, 10000L, seeds[[cohort_name]])
    global_pin <- write_registry_table(sampled$labels, cohort$sample_id, group_levels,
                                       file.path(output, paste0(cohort_name, "_six_group.tsv.gz")))
    companion_pins <- list()
    for (day in group_levels[-1L]) {
      labels <- enumerate_companion(observed, 1L, match(day, group_levels))
      companion_pins[[day]] <- write_registry_table(labels, cohort$sample_id, group_levels,
        file.path(output, paste0(cohort_name, "_", day, "_vs_control_exact.tsv.gz")))
    }
    receipts[[cohort_name]] <- list(
      n_libraries = nrow(cohort), group_counts = setNames(as.list(expected_counts), group_levels),
      original_sample_order = cohort$sample_order, sample_ids = cohort$sample_id,
      metadata_sha256 = digest::digest(file = metadata_file, algo = "sha256"),
      exact_six_group_assignment_space_decimal = exact_space_decimal[[cohort_name]],
      rng = list(kind = RNGkind(), seed = seeds[[cohort_name]],
                 attempts = sampled$attempts, rejected = sampled$rejected_observed_or_duplicate),
      sampled_nonidentity = 10000L, observed_assignment_id = 1L,
      sampled_space = global_pin, exact_day_control_spaces = companion_pins,
      reuse_processing = c("raw", "adapter"),
      monte_carlo_tail_policy = "(1 + nonidentity_tail_count) / 10001; ties included",
      exhaustive_tail_policy = "tail_count_including_observed / exact_space_size; ties included"
    )
    rng_states[[cohort_name]] <- sampled[c("seed", "initial_rng_state", "final_rng_state")]
  }
  saveRDS(rng_states, file.path(output, "RNG_STATES.rds"), version = 3L)
  pins_after <- vapply(inputs, function(path) digest::digest(file = path, algo = "sha256"), character(1L))
  assert_registry(identical(pins, pins_after), "Input/source definition changed during preparation.")
  list(status = "ALLOCATION_REGISTRY_PREPARED_ONLY", cohorts = receipts,
       input_manifest = data.frame(role = names(inputs), path = unname(inputs), sha256 = unname(pins)),
       rng_states_sha256 = digest::digest(file = file.path(output, "RNG_STATES.rds"), algo = "sha256"),
       scientific_matrix_or_rds_read = FALSE, models_fitted = 0L,
       empirical_p_or_q_computed = FALSE, family_coverage_completed = FALSE,
       animal_independence_authenticated = FALSE, batch_or_strata_absence_authenticated = FALSE,
       reference_interpretation = "conditional deposited-library label exchangeability only",
       chronic_or_acute_local_companion_prepared = FALSE,
       new_fit_authority_granted_by_this_receipt = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
assert_registry(length(args) %in% c(2L, 4L), "Use selftest OUTPUT or prepare PROJECT LEDGER OUTPUT.")
mode <- args[[1L]]
assert_registry((mode == "selftest" && length(args) == 2L) ||
                  (mode == "prepare" && length(args) == 4L), "Invalid registry mode/arguments.")
assert_registry(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE),
                  "Existing digest/jsonlite packages are required; no installation is attempted.")
script_flag <- grep("^--file=", commandArgs(), value = TRUE)
assert_registry(length(script_flag) == 1L, "Source path is unavailable.")
script_path <- normalizePath(sub("^--file=", "", script_flag), winslash = "/", mustWork = TRUE)
output <- create_registry_output(args[[length(args)]])
started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
receipt <- if (mode == "selftest") registry_selftest(output) else {
  prepare_registry(args[[2L]], args[[3L]], output, script_path)
}
receipt$started_utc <- started_utc
receipt$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
receipt$command <- commandArgs()
receipt$r_version <- R.version.string
receipt$package_versions <- list(digest = as.character(packageVersion("digest")),
                                  jsonlite = as.character(packageVersion("jsonlite")))
receipt$source_sha256 <- digest::digest(file = script_path, algo = "sha256")
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
jsonlite::write_json(receipt, file.path(output, if (mode == "selftest") "SYNTHETIC_TESTS.json" else "PREPARED.json"),
                      auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat(receipt$status, "\n")
