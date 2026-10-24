# Prospective input-bound scientific execution, NOT a production/file-I/O gate.
# Control mechanics are separately testable with doubles; only the loader below
# binds them to the pinned real scientific components. No legacy entrypoint runs.

gse205906_checked_stage_machine_v1 <- function(steps, binding, evidence_class) {
  expected_steps <- c("check_sources", "check_input", "gene_reference", "program_reference",
    "compute", "gene_compare", "program_compare", "check_result", "finalize")
  if (!is.list(steps) || !identical(names(steps), expected_steps) ||
      !all(vapply(steps, is.function, logical(1))) ||
      length(evidence_class) != 1L || is.na(evidence_class) ||
      !evidence_class %in% c("PINNED_LOCAL_SYNTHETIC", "CONTROL_TEST_DOUBLE")) {
    stop("CHECKED_STAGE_CONTRACT_INVALID", call. = FALSE)
  }
  state <- new.env(parent = emptyenv())
  state$status <- "PREPARED"
  state$stage <- "PREPARED"
  state$events <- "PREPARED"
  state$error_class <- NULL
  state$output_released <- FALSE
  receipt <- function() list(status = state$status, last_stage = state$stage,
    events = state$events, error_class = state$error_class, output_released = state$output_released,
    evidence_class = evidence_class, binding = binding, production_go = FALSE,
    linux_native_acceptance = FALSE, independent_audit = FALSE,
    durable_single_use_enforcement = FALSE)
  execute <- function() {
    if (!identical(state$status, "PREPARED")) stop("CHECKED_HANDLE_ALREADY_EXECUTED", call. = FALSE)
    state$status <- "RUNNING"
    advance <- function(label, action) {
      state$stage <- label
      state$events <- c(state$events, label)
      action()
    }
    guard <- function(label, action) advance(label, function() {
      if (!isTRUE(action())) stop(paste0("CHECKED_GUARD_NOT_TRUE:", label), call. = FALSE)
      invisible(TRUE)
    })
    prepare_reference <- function(label, action) advance(label, function() {
      reference <- action()
      if (!is.list(reference) || !is.function(reference$compare)) {
        stop(paste0("CHECKED_REFERENCE_NOT_READY:", label), call. = FALSE)
      }
      reference
    })
    compare_reference <- function(label, expected_status, action) advance(label, function() {
      checked <- action()
      if (!is.list(checked) || !identical(checked$status, expected_status) ||
          !identical(checked$production_go, FALSE) || !identical(checked$independent_audit, FALSE) ||
          !identical(checked$full_result_validation, FALSE)) {
        stop(paste0("CHECKED_COMPARATOR_RECEIPT_INVALID:", label), call. = FALSE)
      }
      checked
    })
    tryCatch({
      guard("SOURCES_PRE", steps$check_sources)
      guard("INPUT_PRE", steps$check_input)
      gene <- prepare_reference("GENE_REFERENCE", steps$gene_reference)
      guard("INPUT_AFTER_GENE_REFERENCE", steps$check_input)
      program <- prepare_reference("PROGRAM_REFERENCE", steps$program_reference)
      guard("INPUT_AFTER_PROGRAM_REFERENCE", steps$check_input)
      # All expectations exist before the primary calculation can run.
      result <- advance("COMPUTE", steps$compute)
      guard("INPUT_AFTER_COMPUTE", steps$check_input)
      guard("RESULT_SCHEMA", function() steps$check_result(result))
      gene_check <- compare_reference("GENE_COMPARE", "PASS_GENE_NUMERIC_REFIT_ONLY",
        function() steps$gene_compare(gene, result))
      program_check <- compare_reference("PROGRAM_COMPARE", "PASS_NON_GENE_INPUT_DERIVED_REFERENCE_ONLY",
        function() steps$program_compare(program, result))
      guard("SOURCES_POST", steps$check_sources)
      guard("INPUT_POST", steps$check_input)
      finalized <- advance("FINALIZE", function() steps$finalize(result))
      if (!is.list(finalized)) stop("CHECKED_FINALIZED_PAYLOAD_INVALID", call. = FALSE)
      state$stage <- "COMPLETE"
      state$status <- "PASS_LOCAL_INPUT_BOUND_NUMERIC_ONLY"
      state$events <- c(state$events, "COMPLETE")
      state$output_released <- TRUE
      list(payload = finalized, numerical_checks = list(gene = gene_check, program = program_check),
        receipt = receipt())
    }, error = function(error) {
      state$status <- "FAILED_NO_RESULT_RELEASE"
      state$error_class <- class(error)[[1L]]
      state$output_released <- FALSE
      state$events <- c(state$events, "FAILED_NO_RESULT_RELEASE")
      # No tables or partial reference objects are returned via this interface.
      # Keep the original exception visible; receipt stores only its class/stage.
      stop(error)
    })
  }
  list(execute = execute, receipt = receipt)
}

load_gse205906_checked_execution_v1 <- function(project_root) {
  # Freeze the control function value at load time; do not resolve a mutable
  # global name when a later handle is prepared.
  stage_machine <- gse205906_checked_stage_machine_v1
  fail <- function(message) stop(message, call. = FALSE)
  require_true <- function(value, label) if (!isTRUE(value)) fail(label)
  pins <- c(
    "code/R/gse205906_in_memory_scientific_core_v1.R" =
      "fc18bafb9a79ac2cf7eaa073ece8f0a478480658404322120e716cf42eb3201b",
    "code/R/gse205906_in_memory_gene_validator_v1.R" =
      "4188175759fc7db35b18269bf54e92be0b217517226a676ed4f444aaabf8968f",
    "code/R/gse205906_in_memory_program_validator_v1.R" =
      "f91bee54e6228ec19312a3441144bc7bc8e647c9685add0224298d738e10f38c",
    "code/R/analyze_gse205906_bulk_model_v2_r4.R" =
      "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2",
    "code/R/validate_gse205906_bulk_model_v2_r4.R" =
      "c63296157a833611ef92b17f1cbb6b27ab3c23c1028cab79fd706aaab564225d",
    "code/R/gse205906_bulk_model_v2_r4_lib.R" =
      "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3",
    "config/GSE205906_bulk_model_v2_r4.json" =
      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c",
    "config/Route_B2_proportionate_governance_adopted_v1_20260908.json" =
      "6f38731652c497c11174da665e3531bd14e07479d5aa42c7f85f36743a728b4e"
  )
  capture <- function(relative) {
    path <- file.path(project_root, relative)
    info <- file.info(path)
    require_true(nrow(info) == 1L && !is.na(info$size) && !info$isdir &&
      info$size < 1000000 && !nzchar(Sys.readlink(path)), "CHECKED_SOURCE_FILE_INVALID")
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    raw <- readBin(con, "raw", n = as.integer(info$size) + 1L)
    require_true(length(raw) == info$size &&
      identical(digest::digest(raw, "sha256", serialize = FALSE), unname(pins[[relative]])),
      paste0("CHECKED_SOURCE_IDENTITY:", relative))
    raw
  }
  captured <- lapply(names(pins), capture)
  names(captured) <- names(pins)
  source_check <- function() {
    for (relative in names(pins)) capture(relative)
    invisible(TRUE)
  }
  modules <- new.env(parent = baseenv())
  for (relative in names(pins)[1:3]) {
    expressions <- parse(text = rawToChar(captured[[relative]]), keep.source = FALSE)
    require_true(length(expressions) == 1L, "CHECKED_MODULE_DECLARATION_COUNT")
    expression <- expressions[[1L]]
    require_true(is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function")),
      "CHECKED_MODULE_NOT_FUNCTION_DECLARATION")
    eval(expression, modules)
  }
  lockEnvironment(modules, bindings = TRUE)
  core <- modules$load_gse205906_in_memory_core(project_root)
  config <- core$config
  roles <- c("count_matrix", "gene_annotation", "sample_ledger", "gene_resources", "orthology_table")
  expected_tables <- c("gene_edger", "gene_voom", "gene_deseq", "unsigned_camera",
    "unsigned_mroast", "unsigned_fgsea", "activity_results", "mapping_coverage",
    "loao_long", "loao_summary", "design_out", "contrast_out", "filter_out", "annotation_metrics", "gate")
  # Canonical object encoding is for in-memory binding, NOT external file identity.
  # Reject arbitrary R objects that could contain executable closures or references.
  normalize <- function(value, depth = 0L) {
    require_true(depth <= 12L, "CHECKED_INPUT_NESTING")
    if (is.data.frame(value)) {
      require_true(identical(class(value), "data.frame") ||
        identical(class(value), c("data.table", "data.frame")), "CHECKED_INPUT_CLASS")
      require_true(all(names(attributes(value)) %in% c("names", "row.names", "class", ".internal.selfref", "sorted")),
        "CHECKED_INPUT_ATTRIBUTES")
      require_true(!anyDuplicated(names(value)) && !anyNA(names(value)) && all(nzchar(names(value))),
        "CHECKED_INPUT_COLUMN_NAMES")
      copy <- as.data.frame(data.table::copy(value))
      for (name in names(copy)) copy[[name]] <- normalize(copy[[name]], depth + 1L)
      return(copy)
    }
    if (is.matrix(value)) {
      require_true(is.numeric(value) && identical(class(value), c("matrix", "array")), "CHECKED_INPUT_MATRIX")
      require_true(all(names(attributes(value)) %in% c("dim", "dimnames")), "CHECKED_INPUT_ATTRIBUTES")
      return(value)
    }
    if (is.list(value)) {
      require_true(!is.object(value) && !is.null(names(value)) && !anyDuplicated(names(value)) &&
        !anyNA(names(value)) && all(nzchar(names(value))), "CHECKED_INPUT_LIST")
      require_true(all(names(attributes(value)) %in% "names"), "CHECKED_INPUT_ATTRIBUTES")
      return(lapply(value, normalize, depth = depth + 1L))
    }
    require_true((is.null(value) || is.character(value) || is.logical(value) || is.integer(value) || is.double(value)) &&
      !is.object(value), "CHECKED_INPUT_VALUE_TYPE")
    require_true(all(names(attributes(value)) %in% "names"), "CHECKED_INPUT_ATTRIBUTES")
    value
  }
  canonical_inputs <- function(inputs) {
    require_true(is.list(inputs) && !anyDuplicated(names(inputs)) && setequal(names(inputs), roles),
      "CHECKED_INPUT_ROLE_SET")
    inputs <- inputs[roles]
    require_true(is.matrix(inputs$count_matrix) && is.numeric(inputs$count_matrix) &&
      all(is.finite(inputs$count_matrix)) && all(inputs$count_matrix >= 0) &&
      all(inputs$count_matrix == floor(inputs$count_matrix)), "CHECKED_COUNTS_DOMAIN")
    require_true(identical(colnames(inputs$count_matrix), unlist(config$design$sample_ids_in_order, use.names = FALSE)) &&
      !is.null(rownames(inputs$count_matrix)) && !anyNA(rownames(inputs$count_matrix)) &&
      all(nzchar(rownames(inputs$count_matrix))) && !anyDuplicated(rownames(inputs$count_matrix)), "CHECKED_COUNT_KEYS")
    require_true(all(vapply(inputs[c("gene_annotation", "sample_ledger", "orthology_table")], is.data.frame, logical(1))) &&
      is.list(inputs$gene_resources), "CHECKED_INPUT_ROLE_TYPES")
    lapply(inputs, normalize)
  }
  hash_object <- function(value) digest::digest(serialize(value, NULL, version = 3), "sha256", serialize = FALSE)
  identities <- function(inputs) vapply(canonical_inputs(inputs), hash_object, character(1))
  prepare <- function(inputs, expected_input_sha256, domain = "SYNTHETIC") {
    # Real-file/runtime integration is not accepted; a domain string cannot grant it.
    require_true(identical(domain, "SYNTHETIC"), "CHECKED_EXECUTION_PRODUCTION_NOT_ACCEPTED")
    normalized <- canonical_inputs(inputs)
    require_true(is.character(expected_input_sha256) && identical(names(expected_input_sha256), roles) &&
      !anyNA(expected_input_sha256) && all(grepl("^[0-9a-f]{64}$", expected_input_sha256)),
      "CHECKED_INPUT_IDENTITY_MANIFEST")
    observed_hashes <- vapply(normalized, hash_object, character(1))
    require_true(identical(observed_hashes, expected_input_sha256), "CHECKED_INPUT_IDENTITY_MISMATCH")
    serialized <- serialize(normalized, NULL, version = 3)
    snapshot_hash <- digest::digest(serialized, "sha256", serialize = FALSE)
    # Bind one private value snapshot. Later changes to caller-owned tables do not
    # alter it, and none of the components reopen a study-data path.
    captured_inputs <- unserialize(serialized)
    check_input <- function() {
      require_true(identical(vapply(canonical_inputs(captured_inputs), hash_object, character(1)), observed_hashes) &&
        identical(digest::digest(serialized, "sha256", serialize = FALSE), snapshot_hash), "CHECKED_BOUND_INPUT_CHANGED")
      invisible(TRUE)
    }
    check_result <- function(result) {
      require_true(is.list(result) && identical(result$status, "IN_MEMORY_SCIENCE_ONLY_NOT_PRODUCTION_ACCEPTANCE") &&
        identical(result$domain, domain) && identical(result$production_go, FALSE) &&
        is.list(result$tables) && identical(names(result$tables), expected_tables), "CHECKED_RESULT_INTERFACE")
      require_true(identical(result$source_pins, core$source_pins) && identical(result$packages, core$packages) &&
        identical(result$source_spans, core$source_spans) && identical(result$random_seed, config$analysis$random_seed),
        "CHECKED_RESULT_PROVENANCE")
      for (name in expected_tables) {
        table <- result$tables[[name]]
        flags <- c("candidate_gate_eligible", "can_generate_or_upgrade_candidate", "formal_sequential_inference")
        require_true(is.data.frame(table) && !anyDuplicated(names(table)) &&
          all(c("analysis_role", flags) %in% names(table)) &&
          is.character(table$analysis_role) && !anyNA(table$analysis_role) &&
          all(table$analysis_role == config$gate$analysis_role) &&
          all(vapply(flags, function(field) identical(table[[field]], rep(FALSE, nrow(table))), logical(1))),
          paste0("CHECKED_RESULT_CLAIMS:", name))
      }
      invisible(TRUE)
    }
    steps <- list(
      check_sources = source_check, check_input = check_input,
      gene_reference = function() modules$prepare_gse205906_gene_validation(project_root,
        captured_inputs$count_matrix, captured_inputs$gene_annotation, captured_inputs$sample_ledger, domain),
      program_reference = function() do.call(modules$prepare_gse205906_program_validation,
        c(list(project_root = project_root), captured_inputs, list(domain = domain))),
      compute = function() do.call(core$compute, c(captured_inputs, list(domain = domain))),
      gene_compare = function(reference, result) {
        checked <- reference$compare(result)
        require_true(identical(checked$verified_tables, expected_tables[1:3]), "CHECKED_GENE_COVERAGE")
        checked
      },
      program_compare = function(reference, result) {
        checked <- reference$compare(result)
        require_true(identical(checked$verified_tables, expected_tables[-(1:3)]), "CHECKED_PROGRAM_COVERAGE")
        checked
      },
      check_result = check_result,
      finalize = function(result) list(result = result,
        output_sha256 = as.list(vapply(result$tables, function(table) hash_object(as.data.frame(table)), character(1))),
        limitations = c("No actual input-file, species, rights or runtime custody is certified.",
          "Object hashes are not file hashes or external authorization.",
          "All numerical checks share statistical packages; no independent personnel audit.",
          "In-memory handle consumption is not a durable protected single-use ledger."))
    )
    machine <- stage_machine(steps,
      binding = list(object_encoding = "CANONICAL_DATA_FRAMES_R_SERIALIZATION_V3",
        input_sha256 = as.list(observed_hashes), snapshot_sha256 = snapshot_hash,
        source_pins = as.list(pins), domain = domain, fixture_or_caller_manifest_not_external_authority = TRUE),
      evidence_class = "PINNED_LOCAL_SYNTHETIC")
    list(execute = machine$execute, receipt = machine$receipt)
  }
  list(input_identities = identities, prepare = prepare, source_pins = as.list(pins),
    production_go = FALSE, production_consumer_integrated = FALSE)
}
