# Prospective successor: explicit input-bound program family and failure reporting.
# Public execution remains refused; synthetic integration is the current scope.
# Prospective non-KMS workflow component, not a production launcher.
# It evaluates only pinned computational spans, never the legacy entrypoint.
# All study objects are explicit in-memory arguments; no study-data file I/O.
# Invocation sets the frozen RNG seed and advances the session RNG, as r4 did.
# Use an isolated R process; do not treat this component as an authorization gate.

load_gse205906_bound_core_v2 <- function(project_root, program_contract) {
  fail <- function(message) stop(message, call. = FALSE)
  pins <- c(
    "code/R/analyze_gse205906_bulk_model_v2_r4.R" =
      "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2",
    "code/R/gse205906_bulk_model_v2_r4_lib.R" =
      "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3",
    "config/GSE205906_bulk_model_v2_r4.json" =
      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"
  )
  capture_source <- function(relative, expected) {
    path <- file.path(project_root, relative)
    info <- file.info(path)
    if (nrow(info) != 1L || is.na(info$size) || info$isdir ||
        info$size > 1000000 || nzchar(Sys.readlink(path))) {
      fail("CORE_SOURCE_FILE_INVALID")
    }
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    payload <- readBin(con, "raw", n = as.integer(info$size) + 1L)
    observed <- digest::digest(payload, algo = "sha256", serialize = FALSE)
    if (length(payload) != info$size || !identical(observed, expected)) {
      fail("CORE_SOURCE_IDENTITY_MISMATCH")
    }
    rawToChar(payload)
  }
  if (!requireNamespace("digest", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    fail("CORE_CONFIG_DEPENDENCY_MISSING")
  }
  source_text <- Map(capture_source, names(pins), unname(pins))
  model <- program_contract$amend_model(source_text[[1L]])
  shared <- source_text[[2L]]
  cfg <- program_contract$effective_config(jsonlite::fromJSON(source_text[[3L]], simplifyVector = FALSE))
  if (!identical(as.character(getRversion()), cfg$software$R)) {
    fail("CORE_R_VERSION_MISMATCH")
  }
  packages <- c("data.table", "digest", "jsonlite", "edgeR", "limma",
                "DESeq2", "fgsea", "decoupleR")
  versions <- vapply(packages, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) fail("CORE_PACKAGE_MISSING")
    observed <- as.character(utils::packageVersion(package))
    if (!identical(observed, cfg$software[[package]])) {
      fail(paste0("CORE_PACKAGE_VERSION_MISMATCH:", package))
    }
    observed
  }, character(1))

  span_records <- list()
  extract_span <- function(text, start, end, label, origin) {
    starts <- gregexpr(start, text, fixed = TRUE)[[1L]]
    ends <- gregexpr(end, text, fixed = TRUE)[[1L]]
    if (length(starts) != 1L || length(ends) != 1L ||
        starts[[1L]] < 1L || ends[[1L]] <= starts[[1L]]) {
      fail(paste0("CORE_SPAN_AMBIGUOUS:", label))
    }
    payload <- substr(text, starts[[1L]], ends[[1L]] - 1L)
    # Records are local to this loader call, not global scientific state.
    span_records[[label]] <<- list(
      source = origin, start_anchor = start, end_anchor_exclusive = end,
      bytes = nchar(payload, type = "bytes"),
      sha256 = digest::digest(charToRaw(payload), "sha256", serialize = FALSE)
    )
    parse(text = payload, keep.source = FALSE)
  }
  helper_blocks <- list(
    extract_span(shared, "`%||%` <- function", "v2r4_child_request_contract <- function",
                 "basic_assertions", names(pins)[[2L]]),
    extract_span(shared, "v2r4_add_claim_flags <- function", "v2r4_normalize_exact_character_list <- function",
                 "claims_deseq_provenance", names(pins)[[2L]]),
    extract_span(shared, "v2r4_validate_ledger <- function", "v2r4_release_conditional_inclusion <- function",
                 "biological_unit_ledger", names(pins)[[2L]]),
    extract_span(shared, "v2r4_design <- function", "v2r4_output_manifest <- function",
                 "design_covariance_loao", names(pins)[[2L]])
  )
  compute_blocks <- list(
    extract_span(model, "design_contract <- v2r4_design(ledger, cfg)",
                 "gene_sets <- v2r4_read_rds_raw(", "gene_models", names(pins)[[1L]]),
    extract_span(model, "required_gene_sets <- c(", "parent <- dirname(output_dir)",
                 "programs_and_full_loao", names(pins)[[1L]]),
    extract_span(model, "design_out <- data.table::data.table(",
                 "runtime_bindings <- v2r4_add_claim_flags(", "result_tables", names(pins)[[1L]])
  )
  # Helpers must all be function declarations. No legacy bootstrap/control code
  # is executed. The source hashes, not a heuristic allowlist, bind the spans.
  for (block in helper_blocks) {
    for (expression in block) {
      if (!is.call(expression) || !identical(expression[[1L]], as.name("<-")) ||
          !is.call(expression[[3L]]) ||
          !identical(expression[[3L]][[1L]], as.name("function"))) {
        fail("CORE_HELPER_NOT_FUNCTION_DECLARATION")
      }
    }
  }
  helpers <- new.env(parent = baseenv())
  helpers$p.adjust <- stats::p.adjust
  for (block in helper_blocks) for (expression in block) eval(expression, helpers)
  lockEnvironment(helpers, bindings = TRUE)

  compute <- function(count_matrix, gene_annotation, sample_ledger, gene_resources,
                      orthology_table, domain = "SYNTHETIC") {
    if (!is.character(domain) || length(domain) != 1L || is.na(domain) ||
        !domain %in% c("SYNTHETIC", "PRODUCTION")) {
      fail("CORE_DOMAIN_INVALID")
    }
    if (!identical(domain, "SYNTHETIC")) fail("BOUND_PUBLIC_EXECUTION_NOT_ACCEPTED")
    program_contract$check_inputs(list(count_matrix=count_matrix, gene_annotation=gene_annotation,
      sample_ledger=sample_ledger, gene_resources=gene_resources, orthology_table=orthology_table), domain)
    # Domain selects frozen statistical settings, NOT execution permission.
    # Copies prevent data.table's by-reference behavior from altering callers.
    context <- new.env(parent = helpers)
    context$cfg <- cfg
    context$bound_programs <- program_contract$new_session()
    context$cli <- list(execution_mode = domain)
    context$analysis_role <- cfg$gate$analysis_role
    context$ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
    context$ledger <- helpers$v2r4_validate_ledger(data.table::copy(sample_ledger), cfg)
    context$counts <- helpers$v2r4_validate_counts_original(count_matrix, context$ids)
    if (is.null(rownames(context$counts)) || anyDuplicated(rownames(context$counts))) {
      fail("CORE_GENE_IDS_INVALID")
    }
    if (domain == "PRODUCTION" &&
        nrow(context$counts) != cfg$input_snapshot$expected_gene_n_before_filter) {
      fail("CORE_PRODUCTION_GENE_N_INVALID")
    }
    context$assembled <- list(gene_annotation = data.table::copy(gene_annotation))
    context$gene_sets <- data.table::copy(gene_resources)
    context$orthology <- data.table::as.data.table(data.table::copy(orthology_table))
    for (block in compute_blocks) for (expression in block) eval(expression, context)
    # The frozen flag helper expands a completely columnless data.table into
    # one flags-only row. That is not a tested program. Correct only this exact,
    # input-proven empty-family representation in the prospective API; frozen
    # source files and nonempty scientific values remain unchanged.
    empty_unsigned <- all(vapply(context$unsigned_resources,
      function(resource) length(resource$indices) == 0L, logical(1)))
    if (empty_unsigned) {
      flag_columns <- c("analysis_role", "candidate_gate_eligible",
        "can_generate_or_upgrade_candidate", "formal_sequential_inference")
      for (name in c("unsigned_camera", "unsigned_mroast", "unsigned_fgsea")) {
        value <- context[[name]]
        if (!identical(names(value), flag_columns) || nrow(value) != 1L ||
            !identical(value$analysis_role, cfg$gate$analysis_role) ||
            !all(vapply(value[, flag_columns[-1L], with = FALSE],
              function(flag) identical(flag, FALSE), logical(1)))) {
          fail("CORE_EMPTY_UNSIGNED_SENTINEL_SHAPE")
        }
        context[[name]] <- value[0L]
      }
    }
    outputs <- c("gene_edger", "gene_voom", "gene_deseq", "unsigned_camera",
                 "unsigned_mroast", "unsigned_fgsea", "activity_results",
                 "mapping_coverage", "loao_long", "loao_summary", "design_out",
                 "contrast_out", "filter_out", "annotation_metrics", "gate")
    list(
      status = "IN_MEMORY_BOUND_SCIENCE_ONLY_NOT_PRODUCTION_ACCEPTANCE",
      domain = domain, production_go = FALSE,
      program_binding = context$bound_programs$finish(),
      source_amendment = "BOUND_MEMBERSHIP_AND_EXPLICIT_LOAO_FAILURES_V1",
      tables = mget(outputs, envir = context, inherits = FALSE),
      deseq_fit_type = context$deseq_result$actual_fit_type,
      deseq_fit_messages = context$deseq_result$fit_messages,
      random_seed = cfg$analysis$random_seed,
      rng_kind_after = RNGkind(), source_pins = pins,
      source_spans = span_records, packages = versions,
      limitation = "No source-data authorization, species release, runtime custody, publication or independent validation is supplied by this compute-only API."
    )
  }
  list(compute = compute, source_pins = pins, source_spans = span_records,
       packages = versions, config = cfg)
}
