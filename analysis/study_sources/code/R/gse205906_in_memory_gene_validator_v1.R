# Separate gene-level refit implementation from the frozen validator.
# Not a full-result validator: program/activity/LOAO numerics remain unverified.
# Expected values are computed before, and without access to, reported results.

prepare_gse205906_gene_validation <- function(project_root, count_matrix,
                                             gene_annotation, sample_ledger,
                                             domain = "SYNTHETIC") {
  fail <- function(message) stop(message, call. = FALSE)
  if (!is.character(domain) || length(domain) != 1L || is.na(domain) ||
      !domain %in% c("SYNTHETIC", "PRODUCTION")) fail("VALIDATOR_DOMAIN_INVALID")
  pins <- c(
    "code/R/validate_gse205906_bulk_model_v2_r4.R" =
      "c63296157a833611ef92b17f1cbb6b27ab3c23c1028cab79fd706aaab564225d",
    "code/R/gse205906_bulk_model_v2_r4_lib.R" =
      "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3",
    "config/GSE205906_bulk_model_v2_r4.json" =
      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"
  )
  read_source <- function(relative, expected) {
    path <- file.path(project_root, relative)
    bytes <- file.info(path)$size
    if (length(bytes) != 1L || is.na(bytes) || bytes > 1000000 ||
        nzchar(Sys.readlink(path))) fail("VALIDATOR_SOURCE_INVALID")
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    raw <- readBin(con, "raw", n = as.integer(bytes) + 1L)
    if (length(raw) != bytes || !identical(digest::digest(raw, "sha256", serialize = FALSE), expected)) {
      fail("VALIDATOR_SOURCE_IDENTITY_MISMATCH")
    }
    rawToChar(raw)
  }
  sources <- Map(read_source, names(pins), unname(pins))
  validator <- sources[[1L]]
  shared <- sources[[2L]]
  cfg <- jsonlite::fromJSON(sources[[3L]], simplifyVector = FALSE)
  if (!identical(as.character(getRversion()), cfg$software$R)) fail("VALIDATOR_R_VERSION")
  packages <- c("data.table", "digest", "jsonlite", "edgeR", "limma", "DESeq2")
  versions <- vapply(packages, function(package) {
    observed <- as.character(utils::packageVersion(package))
    if (!identical(observed, cfg$software[[package]])) fail(paste0("VALIDATOR_PACKAGE_VERSION:", package))
    observed
  }, character(1))
  records <- new.env(parent = emptyenv())
  span <- function(text, start, end, label) {
    i <- gregexpr(start, text, fixed = TRUE)[[1L]]
    j <- gregexpr(end, text, fixed = TRUE)[[1L]]
    if (length(i) != 1L || length(j) != 1L || i < 1L || j <= i) {
      fail(paste0("VALIDATOR_SPAN_AMBIGUOUS:", label))
    }
    value <- substr(text, i, j - 1L)
    records[[label]] <- list(start_anchor = start, end_anchor_exclusive = end,
      bytes = nchar(value, type = "bytes"), sha256 = digest::digest(charToRaw(value), "sha256", serialize = FALSE))
    parse(text = value, keep.source = FALSE)
  }
  helpers <- new.env(parent = baseenv())
  helpers$p.adjust <- stats::p.adjust
  blocks <- list(
    span(shared, "`%||%` <- function", "v2r4_child_request_contract <- function", "basic"),
    span(shared, "v2r4_add_claim_flags <- function", "v2r4_normalize_exact_character_list <- function", "deseq_provenance"),
    span(shared, "v2r4_validate_ledger <- function", "v2r4_release_conditional_inclusion <- function", "ledger"),
    span(shared, "v2r4_design <- function", "v2r4_output_manifest <- function", "design_covariance"),
    span(validator, "assert_keyed_equal <- function", "\ncomplete <- v2r4_read_tsv_raw(", "comparison")
  )
  for (block in blocks) for (expression in block) {
    if (!is.call(expression) || !identical(expression[[1L]], as.name("<-")) ||
        !is.call(expression[[3L]]) || !identical(expression[[3L]][[1L]], as.name("function"))) {
      fail("VALIDATOR_HELPER_NOT_DECLARATION")
    }
    eval(expression, helpers)
  }
  lockEnvironment(helpers, bindings = TRUE)
  context <- new.env(parent = helpers)
  context$cfg <- cfg
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  ledger <- helpers$v2r4_validate_ledger(data.table::copy(sample_ledger), cfg)
  context$counts <- helpers$v2r4_validate_counts_original(count_matrix, ids)
  if (is.null(rownames(context$counts)) || anyDuplicated(rownames(context$counts))) fail("VALIDATOR_GENE_IDS")
  if (domain == "PRODUCTION" && nrow(context$counts) != cfg$input_snapshot$expected_gene_n_before_filter) {
    fail("VALIDATOR_PRODUCTION_GENE_N")
  }
  design <- helpers$v2r4_design(ledger, cfg)
  context$design <- design$design
  context$contrast_matrix <- design$contrast
  context$group <- design$group
  annotation <- helpers$v2r4_validate_annotation_gene(data.table::copy(gene_annotation),
    rownames(context$counts), cfg$input_snapshot$minimum_nonblank_symbol_coverage)$annotation
  # Input-side preprocessing, never derived from reported filtering or results.
  y0 <- edgeR::DGEList(counts = context$counts)
  fbe <- edgeR::filterByExpr(y0, design = design$design)
  cpm <- rowSums(edgeR::cpm(y0, log = FALSE) >= cfg$analysis$minimum_cpm) >= cfg$analysis$minimum_samples_at_cpm
  keep <- fbe & cpm
  minimum <- if (domain == "SYNTHETIC") cfg$synthetic_contract$minimum_retained_genes else cfg$analysis$minimum_retained_genes_production
  if (sum(keep) < minimum) fail("VALIDATOR_TOO_FEW_FILTERED_GENES")
  context$y0 <- y0
  context$keep <- keep
  context$retained_annotation <- annotation[match(rownames(context$counts)[keep], gene_id)]
  set.seed(as.integer(cfg$analysis$random_seed))
  refits <- list(
    span(validator, "y <- edgeR::normLibSizes(", "edge_reported <- payload_tsv(", "edgeR_refit"),
    span(validator, "voom_input <- limma::voom(", "assert_complete_numeric_family(\n  voom_reported,", "voom_refit"),
    span(validator, "deseq_messages_raw <- character()", "assert_complete_numeric_family(\n  deseq_reported,", "deseq_refit")
  )
  for (block in refits) for (expression in block) eval(expression, context)
  # No primary-core function or reported table is available to these refits.
  expected <- list(gene_edger = context$expected_edger,
                   gene_voom = context$voom_expected, gene_deseq = context$deseq_expected)
  expected$gene_edger[, symbol := context$retained_annotation$symbol[match(gene_id, context$retained_annotation$gene_id)]]
  expected$gene_edger[, method := "edgeR_QL_robust"]
  expected_deseq_fit <- context$deseq_actual_fit
  expected_deseq_messages <- context$deseq_messages
  input_identity <- digest::digest(serialize(list(counts = context$counts, annotation = gene_annotation,
    ledger = sample_ledger, domain = domain), NULL, version = 3), "sha256", serialize = FALSE)
  snapshot <- data.table::copy(expected)
  snapshot_hash <- digest::digest(serialize(snapshot, NULL, version = 3), "sha256", serialize = FALSE)

  compare <- function(result) {
    if (!identical(snapshot_hash, digest::digest(serialize(snapshot, NULL, version = 3), "sha256", serialize = FALSE))) {
      fail("VALIDATOR_REFERENCE_CHANGED")
    }
    if (!is.list(result) || !identical(result$domain, domain) || !identical(result$production_go, FALSE) ||
        !is.list(result$tables) || anyDuplicated(names(result$tables)) ||
        !all(names(snapshot) %in% names(result$tables))) fail("VALIDATOR_RESULT_INTERFACE")
    for (name in names(snapshot)) {
      observed <- data.table::copy(result$tables[[name]])
      flags <- c("candidate_gate_eligible", "can_generate_or_upgrade_candidate", "formal_sequential_inference")
      if (!all(c("analysis_role", flags) %in% names(observed)) ||
          !all(observed$analysis_role == cfg$gate$analysis_role) ||
          !all(vapply(flags, function(key) is.logical(observed[[key]]) && !anyNA(observed[[key]]) &&
                        all(!observed[[key]]), logical(1)))) fail("VALIDATOR_CLAIM_FLAGS")
      want <- snapshot[[name]]
      exact <- intersect(c("symbol", "method", "se_status", "requested_fit_type", "actual_fit_type_after_DESeq"), names(want))
      fields <- names(want)[vapply(want, is.numeric, logical(1))]
      # edgeR uses the original absolute 1e-8 tolerance; the two sensitivity
      # methods use the original absolute-plus-relative frozen tolerance.
      if (name == "gene_edger") {
        helpers$assert_complete_numeric_family(observed, want, c("gene_id", "contrast"), exact,
          character(), "EDGER_FAMILY", cfg)
        helpers$assert_keyed_equal(observed, want, c("gene_id", "contrast"),
          fields, "EDGER_NUMERIC", tolerance = 1e-8)
      } else {
        helpers$assert_complete_numeric_family(observed, want, c("gene_id", "contrast"), exact,
          fields, toupper(name), cfg)
      }
    }
    if (!identical(result$deseq_fit_type, expected_deseq_fit) ||
        !identical(result$deseq_fit_messages, expected_deseq_messages)) fail("VALIDATOR_DESEQ_FIT_PROVENANCE")
    list(status = "PASS_GENE_NUMERIC_REFIT_ONLY", production_go = FALSE,
      full_result_validation = FALSE, independent_audit = FALSE,
      verified_tables = names(snapshot), gene_contrast_rows = nrow(snapshot$gene_edger),
      unverified_numeric_tables = c("unsigned_camera", "unsigned_mroast", "unsigned_fgsea",
        "activity_results", "mapping_coverage", "loao_long", "loao_summary"),
      source_pins = as.list(pins), source_spans = as.list(records), packages = as.list(versions),
      input_object_identity = input_identity,
      limitation = "Separate gene refit code path with shared covariance helper/packages; not independent personnel or full program/LOAO numerical verification."
    )
  }
  list(compare = compare, source_pins = as.list(pins), source_spans = as.list(records),
       input_object_identity = input_identity, production_go = FALSE)
}
