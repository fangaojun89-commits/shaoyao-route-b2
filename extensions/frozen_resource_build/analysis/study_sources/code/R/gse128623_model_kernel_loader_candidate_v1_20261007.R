# Pure RAM definition loader, not an execution wrapper or model admission.
# Authenticates complete original model-v2 SOURCE bytes before parsing/eval.
# Only the known constants/functions and exact validator row-name assignment
# are imported. Both main definitions and the terminal CLI calls are excluded.
# Returned environments contain callable kernels; they must remain in RAM and
# be used only by a separately reviewed parent after genuine input/runtime gates.

g128_kernel_loader <- function(analyzer_raw, validator_raw) {
  need <- function(value, code) {
    if (!isTRUE(value)) {
      stop(paste0("GSE128623_KERNEL_LOADER:", code), call. = FALSE)
    }
  }
  pins <- list(
    analyzer = list(bytes = 43165L,
                    sha256 = "e251c6ae9bb70cd944ebb1d781895878ee80963946f6e8ad94a52a83713add50"),
    validator = list(bytes = 41328L,
                     sha256 = "87b7f3c26b088d352724fb7f5546fc31e8bec56c0c45e3aa490716240d17d08a")
  )
  authenticate <- function(raw, pin) {
    need(is.raw(raw) && is.null(attributes(raw)), "RAW_TYPE")
    need(length(raw) == pin$bytes, "SOURCE_SIZE")
    observed <- tolower(as.vector(as.character(openssl::sha256(raw)), "character"))
    need(identical(observed, pin$sha256), "SOURCE_SHA256")
    invisible(TRUE)
  }
  # Neither SOURCE is parsed, nor any environment populated, before both pins.
  authenticate(analyzer_raw, pins$analyzer)
  authenticate(validator_raw, pins$validator)

  analyzer_constants <- c(
    "V2A_CANDIDATE", "V2A_GROUPS", "V2A_CONTRASTS", "V2A_REQUIRED_PACKAGES",
    "V2A_EDGE_CI_STATUS", "V2A_ENRICHMENT_SEED", "V2A_RNG_KIND", "V2A_RESULT_FLAGS"
  )
  analyzer_functions <- c(
    "v2a_stop", "v2a_assert", "v2a_reset_enrichment_rng", "v2a_require_runtime",
    "v2a_file_identity", "v2a_assert_pin", "v2a_read_json", "v2a_exact_names",
    "v2a_parse_args", "v2a_validate_contract_gate", "v2a_validate_future_pin_set",
    "v2a_read_s1_zero", "v2a_validate_sample_ledger", "v2a_read_inputs",
    "v2a_design", "v2a_contrast_matrix", "v2a_filter", "v2a_edge_tables",
    "v2a_limma_tables", "v2a_deseq_tables", "v2a_collapse_symbol", "v2a_scale_rows",
    "v2a_program_members", "v2a_activity_tables", "v2a_camera_tables",
    "v2a_chronic_program_sensitivities", "v2a_gene_loao", "v2a_program_loao",
    "v2a_diagnostics", "v2a_write_table", "v2a_fit_branch", "v2a_compare_branches",
    "v2a_execute_model"
  )
  validator_constants <- c(
    "V2V_CANDIDATE", "V2V_GROUPS", "V2V_CONTRASTS", "V2V_VERSIONS",
    "V2V_EDGE_CI", "V2V_ENRICHMENT_SEED", "V2V_RNG_KIND"
  )
  validator_functions <- c(
    "v2v_stop", "v2v_require", "v2v_reset_enrichment_rng", "v2v_file_identity",
    "v2v_runtime", "v2v_args", "v2v_json", "v2v_contract_gate", "v2v_read_input",
    "v2v_s1_zero", "v2v_design", "v2v_branch_fit", "v2v_edge", "v2v_limma",
    "v2v_deseq", "v2v_numeric_equal", "v2v_compare_table", "v2v_compare_reconstruction",
    "v2v_members", "v2v_symbol_expression", "v2v_activity", "v2v_diagnostics",
    "v2v_camera", "v2v_ranking", "v2v_program_sensitivities",
    "v2v_gene_loao_keys_and_stats", "v2v_program_loao_stats",
    "v2v_validate_branch", "v2v_validate_output"
  )

  import_definitions <- function(raw, constants, functions, main_name, terminal,
                                 validator = FALSE) {
    expressions <- tryCatch(parse(text = rawToChar(raw), keep.source = FALSE),
                            error = function(error) {
                              stop("GSE128623_KERNEL_LOADER:SOURCE_PARSE", call. = FALSE)
                            })
    expected_names <- c(constants, functions, main_name)
    observed_names <- character()
    selected <- list()
    initialization_seen <- FALSE
    need(length(expressions) == length(expected_names) + 1L + as.integer(validator),
         "AST_COUNT")
    for (index in seq_along(expressions)) {
      node <- expressions[[index]]
      if (index == length(expressions)) {
        need(identical(node, terminal), "CLI_AST")
        next
      }
      if (validator && identical(node, quote(rownames(V2V_CONTRASTS) <- V2V_GROUPS))) {
        need(!initialization_seen && index == 4L, "VALIDATOR_INITIALIZATION_ORDER")
        initialization_seen <- TRUE
        selected[[length(selected) + 1L]] <- node
        next
      }
      need(is.call(node) && length(node) == 3L &&
             identical(node[[1L]], as.name("<-")) && is.symbol(node[[2L]]),
           "TOP_LEVEL_AST")
      name <- as.character(node[[2L]])
      observed_names <- c(observed_names, name)
      need(name %in% expected_names && !anyDuplicated(observed_names), "AST_NAME")
      is_function <- is.call(node[[3L]]) &&
        identical(node[[3L]][[1L]], as.name("function"))
      need(identical(is_function, name %in% c(functions, main_name)), "AST_DECLARATION")
      if (!identical(name, main_name)) selected[[length(selected) + 1L]] <- node
    }
    need(identical(observed_names, expected_names) &&
           identical(initialization_seen, validator), "AST_ORDER_OR_INITIALIZATION")
    target <- new.env(parent = baseenv())
    for (node in selected) eval(node, envir = target)
    need(identical(sort(ls(target, all.names = TRUE)), sort(c(constants, functions))),
         "IMPORTED_BINDINGS")
    for (name in functions) {
      definition <- get(name, envir = target, inherits = FALSE)
      need(is.function(definition) && identical(environment(definition), target),
           "FUNCTION_ENVIRONMENT")
    }
    lockEnvironment(target, bindings = TRUE)
    target
  }
  analyzer <- import_definitions(
    analyzer_raw, analyzer_constants, analyzer_functions, "v2a_main",
    quote(if (sys.nframe() == 0L) v2a_main())
  )
  validator <- import_definitions(
    validator_raw, validator_constants, validator_functions, "v2v_main",
    quote(if (sys.nframe() == 0L) v2v_main()), validator = TRUE
  )
  need(!identical(analyzer, validator) &&
         identical(parent.env(analyzer), baseenv()) &&
         identical(parent.env(validator), baseenv()), "SEPARATE_ENVIRONMENTS")
  need(identical(formals(analyzer$v2a_execute_model)$nrot, 9999L) &&
         identical(formals(analyzer$v2a_execute_model)$minimum_members, 10L) &&
         identical(formals(validator$v2v_validate_output)$nrot, 9999L) &&
         identical(formals(validator$v2v_validate_output)$minimum_members, 10L) &&
         identical(formals(analyzer$v2a_require_runtime)$test_only, FALSE) &&
         identical(formals(validator$v2v_runtime)$test_only, FALSE), "KERNEL_DEFAULTS")

  list(
    schema = "GSE128623_MODEL_KERNEL_DEFINITION_LOADER/v1",
    status = "CONTENT_KERNEL_DEFINITIONS_ONLY",
    analyzer = analyzer, validator = validator, source_pins = pins,
    production_go = FALSE, activation = FALSE, model_run = FALSE,
    runtime_validation_performed = FALSE, scientific_input_read = FALSE
  )
}
