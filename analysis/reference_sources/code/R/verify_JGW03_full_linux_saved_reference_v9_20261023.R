#!/usr/bin/env Rscript
# Saved-only closed-eleven-to-fifteen handoff; five honest scientific actors. No native driver, fitting, RNG, tail_finish, new
# P/q adjustment or filesystem mutation outside one explicitly fresh QA leaf.
# CLI: CONFIG DISPATCH_EXIT_JSON INDEPENDENT_MATH_SOURCE FRESH_OUTPUT GATE_PY PYTHON PARENT_STARTED PARENT_EXIT
# The parent owns resource admission/deployment and the actual readback window.

jgwlc_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
jgwlc_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
jgwlc_digest <- function(value) digest::digest(value, "sha256", serializeVersion = 3L)
jgwlc_same_record <- function(left, right) {
  valid <- function(x) is.list(x) && is.character(x$path) && length(x$path) == 1L &&
    !is.na(x$path) && is.numeric(x$bytes) && length(x$bytes) == 1L &&
    is.finite(x$bytes) && x$bytes >= 0 && x$bytes == floor(x$bytes) &&
    is.character(x$sha256) && length(x$sha256) == 1L && !is.na(x$sha256) &&
    grepl("^[0-9a-f]{64}$", x$sha256)
  valid(left) && valid(right) && identical(left$path, right$path) &&
    left$bytes == right$bytes && identical(left$sha256, right$sha256)
}
jgwlc_archive_pointer <- function(pointer, key, record) {
  # The unchanged jgwx_pin used by the prefix writer has bytes/SHA only;
  # its enclosing key and the independently resolved store provide the path.
  jgwlc_need(identical(pointer$allocation_key, key) &&
    identical(names(pointer$archive), c("bytes", "sha256")) &&
    is.numeric(pointer$archive$bytes) && length(pointer$archive$bytes) == 1L &&
    is.finite(pointer$archive$bytes) && pointer$archive$bytes == record$bytes &&
    identical(pointer$archive$sha256, record$sha256), "Ordered full prefix archive pointer differs")
  invisible(TRUE)
}
jgwlc_constants <- function() list(
  root = "/root/autodl-tmp/JGW03_full_reference_linux_20261017r03",
  config = "eb1a72bbecd61878e61a073b025480a0683c6be76c135257ec2e7943f30bf01c",
  entry = "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d",
  driver = "8f5edfa505e44c6fbf36c605e8b582dfe23c73bb29b547645644003980fec4a0",
  math = "df8d33cf5bba9c4ddb32542b6c4c511d72fa3513dab52b1b48ef7a0ef1acef56",
  finite = "349a5de42e47c41c7211c9f64c0f012d3658c72d840f5b75e4deb8912d52db6b",
  tail = "1074cc5001899cf4f2500c72a52f6f218b1cbb763982c6a0a9b43bfeab82ce8d",
  cells = c("GSE289929", "GSE227494", "GSE299063"))

# Only these previously reviewed independent declarations are imported. The old
# Windows main/top-level code and all producer declarations are never evaluated.
jgwlc_math <- function(path) {
  jgwlc_need(identical(jgwlc_sha(path), jgwlc_constants()$math), "Independent math source differs")
  wanted <- c("jgwz_need", "jgwz_sha", "jgwz_digest", "jgwz_counts", "jgwz_invalid_p",
    "jgwz_validate_compact", "jgwz_increment", "jgwz_mc_endpoints", "jgwz_check_summary")
  environment <- new.env(parent = asNamespace("stats")); found <- character()
  for (expression in parse(path, keep.source = FALSE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% wanted) {
      jgwlc_need(is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function")),
        "Only pure function declarations may be imported")
      eval(expression, envir = environment); found <- c(found, as.character(expression[[2L]]))
    }
  }
  jgwlc_need(identical(found, wanted), "Independent declaration inventory/order differs")
  environment
}

jgwlc_plan <- function() {
  rounds <- list(); cells <- jgwlc_constants()$cells
  for (space in c("four_group_MC", "high_companion_exact")) {
    n <- if (space == "four_group_MC") 10001L else 20L
    for (first in seq.int(1L, n, 50L)) {
      # Frozen Python write_new uses json.dump(sort_keys=True), including each
      # nested job. Match that readback view; array/label order remains exact.
      rounds[[length(rounds) + 1L]] <- lapply(cells, function(accession)
        list(accession = accession, first = first, last = min(first + 49L, n),
          mode = "run", space = space))
    }
  }
  for (space in c("four_group_MC", "high_companion_exact")) {
    rounds[[length(rounds) + 1L]] <- lapply(cells, function(accession)
      list(accession = accession, first = NULL, last = NULL, mode = "finish", space = space))
  }
  rounds
}

jgwlc_terminal <- function(command, started, terminal, job, config_sha, entry_sha, argv, driver) {
  jgwlc_need(identical(command$job, job) && identical(started$job, job) &&
    identical(terminal$job, job) && identical(unlist(command$R_argv, use.names = FALSE), argv) &&
    identical(command$config_sha256, config_sha) && identical(started$config_sha256, config_sha) &&
    identical(command$entry_sha256, entry_sha) && identical(started$entry_sha256, entry_sha) &&
    identical(command$driver$sha256, driver) && is.numeric(command$assigned_CPU) &&
    length(command$assigned_CPU) == 1L && command$assigned_CPU >= 0 &&
    is.numeric(started$pid) && length(started$pid) == 1L && started$pid > 0 &&
    identical(started$pid, terminal$actual_R_pid) &&
    identical(started$start_ticks, terminal$start_ticks) &&
    identical(terminal$status, "ACTUAL_BOUNDED_EXIT_NOT_INDEPENDENT_NUMERICAL_ACCEPTANCE") &&
    identical(terminal$execution_kind, "AFFINITY_BOOTSTRAP_EXECV_SAME_PID") &&
    terminal$actual_native_exit_code == 0 && terminal$actual_R_exit_code == 0 &&
    isTRUE(terminal$memory_budget_met) && is.null(terminal$observer_failure_class_only) &&
    identical(terminal$raw_streams_persisted, FALSE) && identical(terminal$automatic_native_retry, FALSE),
    "Actual Linux COMMAND/STARTED/EXIT or native-zero identity differs")
  invisible(TRUE)
}

jgwlc_registry <- function(registry, domain, space) {
  n <- if (space == "four_group_MC") 10001L else 20L
  matrix <- if (space == "four_group_MC") registry$MC else registry$high_companion
  jgwlc_need(is.matrix(matrix) && is.numeric(matrix) && identical(dim(matrix), c(n, 12L)) &&
    identical(colnames(matrix), domain$sample_ids) && !anyNA(matrix) &&
    all(matrix %in% 0:3) && all(vapply(0:3, function(k)
      all(rowSums(matrix == k) == 3L), logical(1L))), "Full label geometry/sample order differs")
  keys <- unname(apply(matrix, 1L, paste0, collapse = ""))
  jgwlc_need(!anyDuplicated(keys) &&
    identical(keys[1L], paste(as.integer(domain$observed_labels), collapse = "")),
    "Distinct registry order/observed once differs")
  if (space == "high_companion_exact") {
    nuisance <- which(!domain$observed_labels %in% c(0, 3))
    jgwlc_need(all(vapply(seq_len(n), function(i) identical(as.integer(matrix[i, nuisance]),
      as.integer(domain$observed_labels[nuisance])), logical(1L))), "Companion nuisance roles changed")
  }
  keys
}

# Reconstruct the Linux producer identity from its preserved inputs, not by
# calling boot or copying identity from an allocation archive. File records in
# identity are canonical filesystem records (JSON integer bytes are not an RDS
# digest contract). The native fit capsules/fixed layer are never decoded here.
jgwlc_identity <- function(prepared_pin, context, bridges, domain, binding) {
  list(schema = "jgw03-separate-linux-scientific-identity-v2", prepared = prepared_pin,
    config = binding$config, entry = binding$entry, sources = binding$sources,
    inputs = binding$inputs, runtime = binding$runtime, registry = context$registry,
    domain_sha256 = domain$sha256, Linux_bridge_records = lapply(bridges, `[[`, "record"),
    observed = bridges[[1L]]$meta$compact, fixed = context$fixed)
}

# First producer owns a key across both spaces; exact reuse never changes it.
jgwlc_first_actor <- function(key, mc_keys, exact_keys, completed_rounds, old, successor, fit,
                              predecessor_rounds, eight_fit, eight_closed_rounds,
                              fifteen_fit, eleven_closed_rounds) {
  jgwlc_need(is.character(key) && length(key) == 1L && !is.na(key) &&
    !anyDuplicated(mc_keys) && !anyDuplicated(exact_keys) &&
    is.numeric(completed_rounds) && length(completed_rounds) == 1L &&
    is.finite(completed_rounds) && completed_rounds == floor(completed_rounds) &&
    completed_rounds == 66L && is.numeric(predecessor_rounds) && length(predecessor_rounds) == 1L &&
    is.finite(predecessor_rounds) && predecessor_rounds == floor(predecessor_rounds) &&
    predecessor_rounds >= 67L && predecessor_rounds <= 199L &&
    is.numeric(eight_closed_rounds) && length(eight_closed_rounds) == 1L &&
    is.finite(eight_closed_rounds) && eight_closed_rounds == floor(eight_closed_rounds) &&
    eight_closed_rounds > predecessor_rounds && eight_closed_rounds < eleven_closed_rounds &&
    is.numeric(eleven_closed_rounds) && length(eleven_closed_rounds) == 1L &&
    is.finite(eleven_closed_rounds) && eleven_closed_rounds == floor(eleven_closed_rounds) &&
    eleven_closed_rounds >= 93L && eleven_closed_rounds <= 200L,
    "Exact five-actor handoff key domain required")
  position <- match(key, mc_keys)
  if (!is.na(position)) {
    if (position <= 3300L) return(old)
    if (position <= 50L*predecessor_rounds) return(successor)
    if (position <= 50L*eight_closed_rounds) return(eight_fit)
    if (position <= 50L*eleven_closed_rounds) return(fit)
    return(fifteen_fit)
  }
  jgwlc_need(key %in% exact_keys, "Unknown frozen allocation key")
  old
}

jgwlc_main <- function(args) {
  jgwlc_need(length(args) == 8L, "CONFIG DISPATCH MATH FRESH_OUTPUT GATE_PY PYTHON PARENT_STARTED PARENT_EXIT required")
  constants <- jgwlc_constants(); root <- constants$root
  config_path <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
  dispatch_path <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
  out <- gsub("\\\\", "/", args[4L])
  jgwlc_need(identical(config_path, paste0(root, "/CONFIG.json")) &&
    startsWith(dispatch_path, paste0(root, "/provenance/")) &&
    basename(dispatch_path) == "DISPATCH_EXIT.json" &&
    startsWith(out, paste0(root, "/independent_")) && !grepl("(^|/)\\.\\.(/|$)", out) &&
    !file.exists(out) && dir.exists(dirname(out)), "Exact Linux input and fresh QA leaf required")
  checks <- 0; check <- function(ok, label) { jgwlc_need(ok, label); checks <<- checks + 1 }
  check(identical(as.character(getRversion()), "4.6.1"), "Readback R version differs")
  pins <- new.env(hash = TRUE, parent = emptyenv())
  bind <- function(path, expected = NULL) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    check(!identical(basename(path), "checked_payload.rds") && !dir.exists(path) &&
      identical(Sys.readlink(path), ""), "Forbidden payload or non-regular input")
    actual <- list(path = path, bytes = unname(file.info(path)$size), sha256 = jgwlc_sha(path))
    if (!is.null(expected)) {
      check(is.numeric(expected$bytes) && length(expected$bytes) == 1L &&
        actual$bytes == expected$bytes && identical(actual$sha256, expected$sha256), "File bytes/SHA differ")
      if (!is.null(expected$path)) check(identical(actual$path, expected$path), "Record path differs")
    }
    if (exists(path, pins, inherits = FALSE)) check(identical(get(path, pins), actual), "Input changed")
    assign(path, actual, pins); actual
  }
  json <- function(path, expected = NULL) { bind(path, expected); jsonlite::read_json(path, simplifyVector = FALSE) }
  load <- function(path, expected = NULL) { bind(path, expected); readRDS(path) }
  self <- grep("^--file=", commandArgs(), value = TRUE)
  check(length(self) == 1L, "Explicit source file required")
  self <- bind(sub("^--file=", "", self)); math_pin <- bind(args[3L])
  math <- jgwlc_math(math_pin$path)
  config_pin <- bind(config_path); check(identical(config_pin$sha256, constants$config), "Frozen CONFIG differs")
  cfg <- json(config_path)
  check(identical(cfg$schema, "jgw03-separate-linux-reference-config-v2") &&
    identical(cfg$run_root, root) && identical(cfg$project_root, paste0(root, "/project")) &&
    setequal(names(cfg$cells), constants$cells), "Complete separate Linux configuration required")
  sources <- lapply(seq_along(cfg$source_pins), function(i) bind(file.path(cfg$project_root,
    names(cfg$source_pins)[i]), cfg$source_pins[[i]])); names(sources) <- names(cfg$source_pins)
  entry <- sources[["code/R/run_JGW03_full_reference_linux_v3_20261017.R"]]
  check(identical(entry$sha256, constants$entry) &&
    identical(sources[["code/R/JGW03_reference_compact_tail_20261015.R"]]$sha256, constants$tail) &&
    identical(sources[["code/R/JGW03_finite_population_mc_v2_20261015.R"]]$sha256, constants$finite),
    "Frozen Linux entry/tail/finite MC sources differ")
  inputs <- lapply(cfg$input_pins, function(x) bind(x$path, x)) # Bytes only; no raw input decode.

  # Authenticate exact old/new execution union before decoding scientific RDS.
  dispatch <- json(dispatch_path)
  driver <- bind(dispatch$actual_driver$path, dispatch$actual_driver)
  check(identical(driver$sha256, "880b85ee7905d82b560c655b6dd0be6e2242ebffc2713e8b1996b69591988c0c"),
    "Final reviewed successor driver differs")
  gate_pin <- bind(args[5L])
  check(identical(gate_pin$sha256, "70c742ec7036bf9658f59073f51e901c9c92377b7b0d4ef2ddeadfee25cfac19"), "Independent five-actor execution gate differs")
  python_pin <- bind(args[6L])
  parent_start <- bind(args[7L]); parent_exit <- bind(args[8L])
  gate_text <- system2(python_pin$path, c("-I", "-B", shQuote(gate_pin$path),
    shQuote(driver$path), shQuote(dispatch_path), shQuote(parent_start$path), shQuote(parent_exit$path)),
    stdout = TRUE, stderr = FALSE)
  check(is.null(attr(gate_text, "status")), "Actual saved-only metadata gate exit nonzero")
  gate <- jsonlite::fromJSON(paste(gate_text, collapse = "\n"), simplifyVector = FALSE)
  check(gate$scientific_actor_epochs == 5L && gate$eleven_controller_attempts == 2L &&
    gate$eleven_first_attempt_closed_through == 92L &&
    gate$eleven_stopped_original_actual_parent_exit_code == 1L &&
    gate$eleven_continuation_actual_parent_exit_code == -15L &&
    gate$fifteen_actual_parent_exit_code == 0L,
    "Two historical eleven controllers and actual fifteen parent required")
  check(identical(gate$status, "PASS_EXACT_FIVE_ACTOR_EXECUTION_UNION_NOT_NUMERICAL_ACCEPTANCE") &&
    gate$actual_rounds == 204L && gate$legacy_rounds + gate$successor_rounds == 204L &&
    gate$legacy_processes + gate$successor_processes == 3615L-8L*gate$predecessor_closed_rounds-3L*gate$eight_worker_closed_rounds-4L*gate$eleven_worker_closed_rounds &&
    gate$legacy_rounds == 66L && gate$predecessor_closed_rounds >= 67L &&
    gate$predecessor_closed_rounds <= 199L &&
    gate$eight_worker_closed_rounds > gate$predecessor_closed_rounds &&
    gate$eight_worker_closed_rounds < gate$eleven_worker_closed_rounds &&
    gate$eleven_worker_closed_rounds >= 93L && gate$eleven_worker_closed_rounds <= 200L &&
    gate$incremental_rounds == gate$predecessor_closed_rounds-66L &&
    gate$eight_worker_rounds == gate$eight_worker_closed_rounds-gate$predecessor_closed_rounds &&
    gate$eleven_worker_rounds == gate$eleven_worker_closed_rounds-gate$eight_worker_closed_rounds &&
    gate$fifteen_worker_rounds == 204L-gate$eleven_worker_closed_rounds &&
    identical(gate$scientific_acceptance, FALSE) && jgwlc_same_record(gate$legacy_entry, entry),
    "Complete exact execution union required")
  invisible(lapply(gate$pins, function(record) bind(record$path, record)))
  successor_entry <- gate$incremental_entry; finish_entry <- gate$finish_entry
  fit_entry <- gate$fit_entry; merge_entry <- gate$merge_entry
  eight_fit_entry <- gate$eight_fit_entry; eight_merge_entry <- gate$eight_merge_entry
  fifteen_fit_entry <- gate$fifteen_fit_entry; fifteen_merge_entry <- gate$fifteen_merge_entry
  check(identical(fifteen_fit_entry$sha256, "7eabb266f29fb73fe172f08e090de8e1d9c1984b7ec9dec34104d0fa84d51450") &&
    identical(fifteen_merge_entry$sha256, "92945fe0dba456f6d2594e138a37115e7cfa8925738e1d640a5dde4683a76347"), "Reviewed fifteen-worker actors differ")
  check(identical(eight_fit_entry$sha256, "961b605d670812f2d9595dd03b3da60ad0fcd391becb82ef87ca6c058a8fd7db") &&
    identical(eight_merge_entry$sha256, "b48247c00e2a4bbb606d93fbe3bd0a542f75d31a73abc179eb65055bd0472056"),
    "Historical eight-worker actors differ")
  check(identical(fit_entry$sha256, "38fd274dcf96af76f3a4056dbd3021d396490da67da7064da8c8fe4c8357cd02") &&
    identical(merge_entry$sha256, "23450a49ce6a48ebf6c27d6de4835d4f3e5c27bcc6fdff90b3740ae08dfd9907"), "Reviewed distinct fit and saved-merge actors")
  check(identical(successor_entry$sha256, "452acfa10dbbaa29631ad41877feab120ba421dec6852bd42ab5180ef9b1b56c") &&
    identical(finish_entry$sha256, "69ba3dd31a7854462a22e03c0121d825ea8b228a0278d35d5862f1ac8a7d5613"),
    "Exact incremental/finish actors differ")
  summaries <- list(); cache <- new.env(hash = TRUE, parent = emptyenv())
  for (accession in constants$cells) {
    leaf <- file.path(root, "intermediate", accession)
    prepared <- json(file.path(leaf, "PREPARED.json")); prep_pin <- bind(prepared$prepared$path, prepared$prepared)
    context <- load(prep_pin$path, prep_pin)
    check(identical(context$schema, "jgw03-separate-linux-preparation-v2") &&
      identical(context$accession, accession) && identical(context$root, cfg$project_root) &&
      identical(context$run_root, root) && context$full_MC_nonidentity_draws == 10000L &&
      context$full_companion_including_observed == 20L && context$full_MC_family_n == 89L &&
      context$full_companion_family_n == 34L && identical(context$binding$config, config_pin) &&
      identical(context$binding$entry, entry) && identical(context$binding$sources, sources) &&
      identical(context$binding$inputs, inputs) && identical(context$binding$runtime, context$runtime),
      "Prepared Linux scientific binding/scope differs")
    bind(context$fixed$path, context$fixed) # Hash only. No fixed counts/S4 decode.
    bridge_fixed <- json(context$fixed_bridge$path, context$fixed_bridge)
    check(all(unlist(bridge_fixed)), "Saved fixed-layer bridge not accepted")
    bind(context$session_info$path, context$session_info)
    bridges <- lapply(1:2, function(id) {
      path <- file.path(leaf, paste0("bridge_", id), "BRIDGE_COMPLETE.json")
      record <- bind(path); meta <- json(path)
      check(identical(meta$status, "LINUX_NATIVE_COMPATIBILITY_PASS_V2") && meta$allocation_id == id &&
        identical(meta$prepared_sha256, prep_pin$sha256) && identical(meta$config_sha256, constants$config) &&
        identical(meta$entry_sha256, constants$entry) && isTRUE(meta$compatibility_all_fields_pass),
        "Own-Linux bridge identity differs")
      bind(meta$compact$path, meta$compact)
      list(record = record, meta = meta)
    })
    domain <- load(bridges[[1L]]$meta$domain$path, bridges[[1L]]$meta$domain)
    observed <- load(bridges[[1L]]$meta$compact$path, bridges[[1L]]$meta$compact)
    second <- load(bridges[[2L]]$meta$compact$path, bridges[[2L]]$meta$compact)
    registry <- load(context$registry$path, context$registry)
    check(context$registry$bytes == cfg$cells[[accession]]$registry$bytes &&
      identical(context$registry$sha256, cfg$cells[[accession]]$registry$sha256),
      "Frozen original registry was not byte-preserved")
    body <- domain; body$sha256 <- NULL; full_n <- if (accession == "GSE289929") 31L else 29L
    check(identical(domain$schema, "jgw03-native-compact-domain-v1") && identical(domain$accession, accession) &&
      identical(jgwlc_digest(body), domain$sha256) && identical(domain$observed_labels, context$observed_labels) &&
      length(domain$families) == full_n && length(domain$sample_ids) == 12L &&
      !anyDuplicated(domain$sample_ids) && identical(names(domain$observed_labels), domain$sample_ids) &&
      identical(domain$definition$empirical_BH_computed, FALSE) &&
      identical(domain$definition$native_P_q_changed, FALSE) &&
      identical(domain$definition$Mixed_outer_family_added, FALSE), "Own-Linux domain/family/BH scope differs")
    for (id in names(domain$families)) {
      map <- domain$families[[id]]$map
      check(identical(id, paste(map$method, map$collection, map$contrast, sep = "::")) &&
        identical(map$accession, accession) && identical(map$space, "four_group_MC") &&
        map$null_draw_n == 10000 && isTRUE(map$observed_once), "Original ordered family map differs")
    }
    binding <- context$binding
    identity <- jgwlc_identity(prep_pin, context, bridges, domain, binding)
    identity_sha <- jgwlc_digest(identity); store <- file.path(root, "derived", accession)
    cell_mc_keys <- jgwlc_registry(registry, domain, "four_group_MC")
    cell_exact_keys <- jgwlc_registry(registry, domain, "high_companion_exact")
    cell_actors <- lapply(unique(c(cell_mc_keys, cell_exact_keys)), function(key)
      jgwlc_first_actor(key, cell_mc_keys, cell_exact_keys, gate$legacy_rounds, entry, successor_entry,
        fit_entry, gate$predecessor_closed_rounds, eight_fit_entry, gate$eight_worker_closed_rounds,
        fifteen_fit_entry, gate$eleven_worker_closed_rounds))
    names(cell_actors) <- unique(c(cell_mc_keys, cell_exact_keys))
    cell_summary <- list()
    for (space in c("four_group_MC", "high_companion_exact")) {
      keys <- jgwlc_registry(registry, domain, space); n <- length(keys)
      check(identical(observed$allocation_key, keys[1L]) &&
        identical(bridges[[1L]]$meta$allocation_key, keys[1L]) &&
        identical(bridges[[2L]]$meta$allocation_key, jgwlc_registry(registry, domain, "four_group_MC")[2L]),
        "Own observed/ID2 label association differs")
      selected <- names(domain$families)
      if (space == "high_companion_exact") selected <- selected[vapply(domain$families,
        function(x) identical(x$map$contrast, "high_vs_vehicle"), logical(1L))]
      check(length(selected) == if (space == "four_group_MC") full_n else
        if (accession == "GSE289929") 12L else 11L, "Complete selected family set differs")
      prefix <- load(file.path(store, sprintf("%s_TAIL_PREFIX_%05d.rds", space, n)))
      product <- load(file.path(store, paste0(space, "_TAIL_SUMMARY.rds")))
      check(identical(prefix$schema, "jgw03-fixed-registry-tail-prefix-v1") &&
        identical(prefix$space, space) && prefix$through_allocation_id == n &&
        identical(prefix$scientific_identity_sha256, identity_sha) &&
        identical(prefix$full_registry, context$registry) && identical(prefix$observed, bridges[[1L]]$meta$compact) &&
        identical(prefix$provisional_P_computed, FALSE) &&
        identical(names(prefix$completed_compact_pins), as.character(seq.int(2L, n))) &&
        identical(product$tail_checkpoint, prefix$state) && identical(product$scientific_identity_sha256, identity_sha) &&
        identical(product$full_registry, context$registry) && jgwlc_same_record(product$execution_entry, finish_entry),
        "Full saved prefix/summary lineage differs")
      state <- prefix$state; summary <- product$summary
      check(identical(state$schema, "jgw03-tail-accumulator-v1") && identical(state$space, space) &&
        identical(state$domain_sha256, domain$sha256) && state$total_N == n &&
        state$processed_nonidentity == n - 1L && identical(state$allocation_keys, keys) &&
        identical(state$definition, domain$definition) && identical(names(state$families), selected),
        "Full tail state order/denominator differs")
      counters <- lapply(selected, function(id) {
        origin <- observed$families[[id]]
        check(all((origin$outer_statistic >= origin$outer_statistic)[origin$outer_available]),
          "Available observed must self-hit including Inf")
        list(hits = integer(length(origin$outer_statistic)), unknown = integer(length(origin$outer_statistic)))
      }); names(counters) <- selected
      totals <- c(rows = 0, numeric = 0, dynamic = 0); metadata_missing <- 0L
      for (i in seq_len(n)) {
        expected_actor <- cell_actors[[keys[i]]]
        archive_leaf <- file.path(store, keys[i]); meta <- json(file.path(archive_leaf, "COMPLETE.json"))
        owned <- json(file.path(archive_leaf, "STARTED.json"))
        archive_pin <- bind(file.path(archive_leaf, "COMPACT.rds"), meta$archive)
        archive <- load(archive_pin$path, archive_pin)
        check(identical(meta$status, "IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT") &&
          identical(meta$allocation_key, keys[i]) && identical(meta$scientific_identity_sha256, identity_sha) &&
          identical(owned$status, "ALLOCATION_OWNED_BEFORE_NATIVE_CALL") &&
          identical(owned$allocation_key, keys[i]) && identical(owned$scientific_identity_sha256, identity_sha) &&
          jgwlc_same_record(owned$execution_entry, expected_actor) &&
          identical(archive$schema, "jgw03-full-reference-allocation-v1") &&
          identical(archive$scientific_identity_sha256, identity_sha) && identical(archive$accession, accession) &&
          identical(archive$allocation_key, keys[i]) && jgwlc_same_record(archive$execution_entry, expected_actor) &&
          identical(archive$compact_sha256, jgwlc_digest(archive$compact)) &&
          archive$new_native_calls %in% 0:1 && archive$new_native_calls == meta$new_native_calls &&
          identical(owned$new_fit_intended, archive$new_native_calls == 1L),
          "Saved full allocation ownership/content identity differs")
        draw <- archive$compact; totals <- totals + math$jgwz_validate_compact(draw, domain, keys[i], check)
        if (i == 1L) check(identical(draw, observed) && archive$new_native_calls == 0L,
          "Linux observed must be exact zero-fit reuse") else {
          pointer <- prefix$completed_compact_pins[[as.character(i)]]
          jgwlc_archive_pointer(pointer, keys[i], archive_pin)
          if (identical(keys[i], second$allocation_key)) check(identical(draw, second) &&
            archive$new_native_calls == 0L, "Linux bridge ID2 must be exact zero-fit reuse")
          for (id in selected) {
            add <- math$jgwz_increment(observed$families[[id]], draw$families[[id]])
            counters[[id]]$hits <- counters[[id]]$hits + add$hits
            counters[[id]]$unknown <- counters[[id]]$unknown + add$unknown
          }
        }
        metadata_missing <- metadata_missing + as.integer(!identical(draw$execution_metadata$status,
          "CALLER_SUPPLIED_ACTUAL_EXECUTION_METADATA"))
        rm(archive, draw); if (i %% 25L == 0L) invisible(gc())
      }
      check(identical(summary$schema, "jgw03-tail-summary-v1") && identical(summary$space, space) &&
        identical(summary$domain_sha256, domain$sha256) && summary$processed_nonidentity == n - 1L &&
        identical(names(summary$families), selected) && identical(summary$finite_mc_source_sha256, constants$finite) &&
        summary$execution_metadata_missing_allocations == metadata_missing &&
        state$execution_metadata_missing_allocations == metadata_missing &&
        identical(summary$execution_metadata_coverage_complete, metadata_missing == 0L) &&
        identical(summary$execution_metadata_is_caller_supplied_not_independently_authenticated, TRUE) &&
        identical(summary$native_P_q_changed, FALSE) && identical(summary$empirical_BH_computed, FALSE) &&
        identical(summary$uncertainty_scope, "Pointwise conditional finite-population MC only; no inner-MC/biological/FDR coverage") &&
        identical(summary$reference_acceptance, "Arithmetic only, not independent native or scientific acceptance"),
        "Saved summary coverage/BH/uncertainty scope differs")
      for (id in selected) {
        origin <- observed$families[[id]]; saved <- state$families[[id]]; count <- counters[[id]]
        check(identical(saved$observed, origin$outer_statistic) &&
          identical(saved$observed_available, origin$outer_available) &&
          identical(saved$tail_hits_nonidentity, count$hits) && identical(saved$null_unavailable, count$unknown),
          "Independent inclusive full-space counters differ")
        math$jgwz_check_summary(summary$families[[id]], count$hits, count$unknown,
          origin$outer_available, space, cache, check)
      }
      cell_summary[[space]] <- list(allocation_n = n, selected_family_n = length(selected),
        archived_full_native_family_n_per_draw = full_n, feature_rows_readback = totals[["rows"]],
        numeric_values_readback = totals[["numeric"]], dynamic_values_readback = totals[["dynamic"]],
        execution_metadata_missing_allocations = metadata_missing,
        tail_hits_nonidentity = sum(vapply(counters, function(x) sum(x$hits), numeric(1))),
        null_unavailable_feature_draws = sum(vapply(counters, function(x) sum(x$unknown), numeric(1))))
      rm(prefix, product, state, summary, counters); invisible(gc())
    }
    summaries[[accession]] <- list(scientific_identity_sha256 = identity_sha,
      own_Linux_observed_and_ID2_reused_exactly = TRUE, spaces = cell_summary)
    rm(context, domain, observed, second, registry); invisible(gc())
  }
  # A complete second byte pass is intentional provenance, not a new native run.
  paths <- sort(ls(pins, all.names = TRUE)); manifest <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    record <- get(paths[i], pins)
    check(unname(file.info(paths[i])$size) == record$bytes && identical(jgwlc_sha(paths[i]), record$sha256),
      "Bound file/source changed during saved-only readback")
    manifest[[i]] <- as.data.frame(record, stringsAsFactors = FALSE)
  }
  check(!file.exists(out) && dir.create(out), "Fresh independent output was consumed")
  manifest_path <- file.path(out, "INPUT_PIN_MANIFEST.tsv.gz"); connection <- gzfile(manifest_path, "wt")
  tryCatch(utils::write.table(do.call(rbind, manifest), connection, sep = "\t", quote = TRUE,
    row.names = FALSE), finally = close(connection))
  receipt <- list(status = "PASS_INDEPENDENT_COMPLETE_LINUX_SAVED_REFERENCE_ONLY", checks = checks,
    UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), source = self, independent_math = math_pin,
    config = config_pin, dispatch = get(dispatch_path, pins), actual_rounds = 204L,
    execution_gate = gate_pin, parent_started = parent_start, parent_native_exit = parent_exit,
    scientific_actor_epochs = gate$scientific_actor_epochs,
    eleven_controller_attempts = gate$eleven_controller_attempts,
    eleven_first_attempt_closed_through = gate$eleven_first_attempt_closed_through,
    eleven_stopped_original_actual_parent_exit_code = gate$eleven_stopped_original_actual_parent_exit_code,
    eleven_continuation_actual_parent_exit_code = gate$eleven_continuation_actual_parent_exit_code,
    fifteen_actual_parent_exit_code = gate$fifteen_actual_parent_exit_code,
    legacy_rounds = gate$legacy_rounds, successor_rounds = gate$successor_rounds,
    legacy_processes = gate$legacy_processes, successor_processes = gate$successor_processes,
    old_controller_success_asserted = FALSE, actual_finish_entry = finish_entry,
    actual_R_process_exits_zero = gate$legacy_processes+gate$successor_processes,
    incremental_rounds = gate$incremental_rounds, eight_worker_rounds = gate$eight_worker_rounds,
    actual_fit_entry = fifteen_fit_entry, actual_saved_merge_entry = fifteen_merge_entry,
    historical_eleven_fit_entry = fit_entry, historical_eleven_merge_entry = merge_entry,
    historical_eight_fit_entry = eight_fit_entry, historical_eight_merge_entry = eight_merge_entry,
    eleven_worker_rounds = gate$eleven_worker_rounds, eight_worker_closed_rounds = gate$eight_worker_closed_rounds,
    fifteen_worker_rounds = gate$fifteen_worker_rounds, eleven_worker_closed_rounds = gate$eleven_worker_closed_rounds,
    predecessor_closed_rounds = gate$predecessor_closed_rounds,
    archive_actor_boundaries = list(legacy_through = 3300L,
      incremental_through = 50L*gate$predecessor_closed_rounds,
      eight_worker_through = 50L*gate$eight_worker_closed_rounds,
      eleven_worker_through = 50L*gate$eleven_worker_closed_rounds,
      exact = "Shared keys retain MC first-owner; exact-only keys legacy; merge is not a fit actor"),
    all_six_saved_only_finishes_zero = TRUE,
    cells = summaries, total_MC_families = 89L, total_companion_families = 34L,
    bound_files = length(paths), all_bound_pins_unchanged = TRUE,
    pin_manifest = list(path = manifest_path, bytes = unname(file.info(manifest_path)$size),
      sha256 = jgwlc_sha(manifest_path)),
    distinct_MC_hit_counts_with_independent_endpoint_checks = length(ls(cache, all.names = TRUE)),
    native_P_q_and_dynamic_fields = "Saved full vectors, types, NA/NaN/signed Inf/zero and diagnostic counts read back; not independently refitted",
    native_BH_recomputed = FALSE, outer_empirical_BH_computed = FALSE, BH_intervals_computed = FALSE,
    endpoint_oracle = "Pinned independent endpoint/neighbor hypergeometric log-tail checks with accepted conservative 8eps guard; no producer tail/interval import",
    scientific_payloads_decoded = "Only Linux PREPARED metadata, registry, domain, compact/prefix/summary RDS; no fixed/native S4 capsules or checked_payload",
    new_native_fit_calls = 0L, unavailable_values_imputed = FALSE, denominator_reduced = FALSE,
    finite_MC_scope = "Conditional fixed-per-assignment pointwise finite population; not inner-MC, biology, simultaneous or FDR coverage",
    full_native_algorithm_or_biological_acceptance = FALSE, original_records_removed = 0L)
  jsonlite::write_json(receipt, file.path(out, "INDEPENDENT_COMPLETE_LINUX_REFERENCE_READBACK.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
  cat(receipt$status, checks, "\n"); invisible(0L)
}

if (sys.nframe() == 0L) {
  code <- tryCatch(jgwlc_main(commandArgs(trailingOnly = TRUE)), error = function(condition) {
    message("JGW03_LINUX_SAVED_REFERENCE_STOP_NO_RAW_CONDITION_PERSISTED"); 1L
  })
  quit(status = code)
}


