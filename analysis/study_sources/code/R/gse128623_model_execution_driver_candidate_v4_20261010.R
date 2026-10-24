# Execution-driver component, not a native launcher or an authorization gate.
# Inputs are authenticated public SOURCE raw vectors and a parent's bound path
# list. A future native parent owns input/technical closure, one-way CI adoption,
# held-file/path custody, resource/owner lifecycle, logs and final publication.
# This function intentionally runs the unchanged kernels when explicitly called.
# It must not be called for this SOURCE/MOCK review or with real inputs on Windows.
# No CLI, caller GO flag, supplied function, test_only, nrot or member threshold.
# The original A/V CLI stops and statistical definitions remain unchanged.

g128_model_execution_driver <- function(analyzer_raw, validator_raw, loader_raw,
                                       member_guard_raw, comparison_guard_raw,
                                       paths) {
  need <- function(ok, code) {
    if (!isTRUE(ok)) stop(paste0("GSE128623_MODEL_DRIVER:", code), call. = FALSE)
  }
  source_pins <- list(
    loader = list(bytes = 6822L,
      sha256 = "7da5e5abdd028e3a5d237a8eb01763b5244f65252267ef6f1417cefef8988014"),
    member_guard = list(bytes = 6827L,
      sha256 = "eb4158c7849dac19e5e046b3b74e9dcf2724a1fb6929e791d20138bc15b2cf0e"),
    comparison_guard = list(bytes = 4559L,
      sha256 = "85633c1cd627744106087fd901142ee45a085cc51492be59906c73b07f46433c")
  )
  plain_sha <- function(raw) {
    tolower(as.vector(as.character(openssl::sha256(raw)), "character"))
  }
  authenticate <- function(raw, pin) {
    need(is.raw(raw) && is.null(attributes(raw)), "SOURCE_RAW_TYPE")
    need(length(raw) == pin$bytes && identical(plain_sha(raw), pin$sha256),
         "SOURCE_IDENTITY")
  }
  authenticate(loader_raw, source_pins$loader)
  authenticate(member_guard_raw, source_pins$member_guard)
  authenticate(comparison_guard_raw, source_pins$comparison_guard)
  import_function <- function(raw, name) {
    nodes <- tryCatch(parse(text = rawToChar(raw), keep.source = FALSE),
                      error = function(error) {
                        stop("GSE128623_MODEL_DRIVER:SOURCE_PARSE", call. = FALSE)
                      })
    need(length(nodes) == 1L && is.call(nodes[[1L]]) &&
           identical(nodes[[1L]][[1L]], as.name("<-")) &&
           identical(nodes[[1L]][[2L]], as.name(name)) &&
           is.call(nodes[[1L]][[3L]]) &&
           identical(nodes[[1L]][[3L]][[1L]], as.name("function")), "SOURCE_AST")
    module <- new.env(parent = baseenv())
    eval(nodes, module)
    lockEnvironment(module, bindings = TRUE)
    get(name, module, inherits = FALSE)
  }
  loader <- import_function(loader_raw, "g128_kernel_loader")
  member_guard <- import_function(member_guard_raw, "guard")
  comparison_guard <- import_function(comparison_guard_raw, "g128_comparison_guard")
  # The fixed loader authenticates both complete A/V files before AST import.
  kernels <- tryCatch(loader(analyzer_raw, validator_raw), error = function(error) {
    stop("GSE128623_MODEL_DRIVER:KERNEL_SOURCE", call. = FALSE)
  })
  analyzer <- kernels$analyzer
  validator <- kernels$validator

  path_names <- c("count_matrix", "feature_annotation", "sample_ledger",
                  "frozen_design", "gene_set_resource", "s1_zero", "staging",
                  "log_directory")
  need(is.list(paths) && !is.object(paths) &&
         identical(names(paths), path_names), "PATH_SCHEMA")
  for (path in paths) {
    need(is.character(path) && !is.object(path) && length(path) == 1L &&
           !is.na(path) && is.null(attributes(path)) && startsWith(path, "/") &&
           identical(path, trimws(path)) && !grepl("[\r\n\t]", path) &&
           !grepl("//|/\\.(/|$)|/\\.\\.(/|$)|/$", path), "PATH_VALUE")
  }
  need(!anyDuplicated(unlist(paths, use.names = FALSE)), "PATH_COLLISION")
  # Native parent must also enforce nofollow/held identity; these are only R-level
  # lexical/existence checks. No staging directory or log file is created here.
  staging_link <- Sys.readlink(paths$staging)
  need(!file.exists(paths$staging) &&
         (is.na(staging_link) || !nzchar(staging_link)),
       "STAGING_ALREADY_PRESENT")
  need(dir.exists(dirname(paths$staging)) && dir.exists(paths$log_directory),
       "OUTPUT_PARENT_REQUIRED")
  canonical_path <- function(path) {
    tryCatch(normalizePath(path, winslash = "/", mustWork = TRUE),
             error = function(error) {
               stop("GSE128623_MODEL_DRIVER:PATH_CANONICAL", call. = FALSE)
             })
  }
  canonical_log <- canonical_path(paths$log_directory)
  canonical_staging <- file.path(
    canonical_path(dirname(paths$staging)),
    basename(paths$staging)
  )
  inside <- function(path, root) identical(path, root) || startsWith(path, paste0(root, "/"))
  need(!inside(canonical_log, canonical_staging) &&
         !inside(canonical_staging, canonical_log), "LOG_STAGING_OVERLAP")
  for (name in path_names[1:6]) {
    info <- file.info(paths[[name]])
    need(nrow(info) == 1L && !is.na(info$isdir) && !info$isdir &&
           !inside(canonical_path(paths[[name]]),
                   canonical_staging), "INPUT_FILE_REQUIRED")
  }

  # No original exception message, scientific value or private path is exported.
  # The native parent captures this fixed phase code rather than raw R stderr.
  stage <- function(name, expression) {
    tryCatch(force(expression), error = function(error) {
      stop(paste0("GSE128623_MODEL_DRIVER:", name), call. = FALSE)
    })
  }
  file_identity <- function(path) {
    before <- file.info(path)
    need(nrow(before) == 1L && !is.na(before$isdir) && !before$isdir &&
           is.finite(before$size) && before$size >= 0, "FILE_IDENTITY")
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    sha <- plain_sha(connection)
    after <- file.info(path)
    fields <- c("size", "isdir", "mode", "mtime", "ctime")
    need(identical(before[fields], after[fields]), "FILE_DRIFT")
    list(bytes = unname(as.numeric(before$size)), sha256 = sha)
  }
  read_table <- function(path) {
    data.table::fread(path, sep = "\t", header = TRUE,
                      na.strings = c("NA", "NaN", ""))
  }
  branch_files <- c(
    "gene_level/edgeR_primary.tsv.gz", "gene_level/limma_voom_sensitivity.tsv.gz",
    "gene_level/DESeq2_Wald_sensitivity.tsv.gz", "diagnostics/library_sizes_TMM.tsv",
    "diagnostics/MDS_coordinates.tsv", "diagnostics/dispersion_BCV_QL.tsv.gz",
    "diagnostics/mean_variance.tsv.gz", "diagnostics/sample_influence.tsv",
    "diagnostics/prefilter_summary.tsv", "programs/camera.tsv.gz",
    "programs/mroast_chronic.tsv.gz", "programs/fgsea_chronic.tsv.gz",
    "programs/activity_all_contrasts.tsv.gz", "loao/gene_chronic_true_refit.tsv.gz",
    "loao/program_chronic_true_refit.tsv.gz"
  )
  expected_outputs <- c(file.path("primary_all30", branch_files),
                        file.path("complete29_sensitivity", branch_files),
                        "comparison/all30_vs_complete29_gene.tsv.gz",
                        "comparison/all30_vs_complete29_summary.tsv", "ANALYZER_COMPLETE.tsv")
  comparison_files <- file.path(paths$staging, c(
    "primary_all30/gene_level/edgeR_primary.tsv.gz",
    "complete29_sensitivity/gene_level/edgeR_primary.tsv.gz",
    "comparison/all30_vs_complete29_gene.tsv.gz",
    "comparison/all30_vs_complete29_summary.tsv"
  ))

  # Private closure, not a caller-injectable runner. Tests extract this exact AST
  # into an artificial environment; production has no mock/test bypass argument.
  execute_sequence <- function() {
    stage("ANALYZER_RUNTIME", analyzer$v2a_require_runtime(test_only = FALSE))
    stage("VALIDATOR_RUNTIME", validator$v2v_runtime(test_only = FALSE))
    stage("ANALYZER_S1_ZERO", analyzer$v2a_read_s1_zero(paths$s1_zero))
    stage("VALIDATOR_S1_ZERO", validator$v2v_s1_zero(paths$s1_zero))
    analyzer_input <- stage("ANALYZER_INPUT", analyzer$v2a_read_inputs(
      paths$count_matrix, paths$feature_annotation, paths$sample_ledger,
      paths$frozen_design, "gene_id", "GeneID", "gene"
    ))
    validator_input <- stage("VALIDATOR_INPUT", validator$v2v_read_input(
      paths$count_matrix, paths$feature_annotation, paths$sample_ledger,
      "gene_id", "GeneID", "gene"
    ))
    stage("INPUT_READER_RELATION", {
      need(identical(analyzer_input$counts, validator_input$counts) &&
             identical(analyzer_input$symbols, validator_input$symbols) &&
             identical(as.data.frame(analyzer_input$annotation),
                       as.data.frame(validator_input$annotation)) &&
             identical(as.data.frame(analyzer_input$ledger),
                       as.data.frame(validator_input$ledger)), "INPUT_READER_RELATION")
    })
    resource_pin <- stage("RESOURCE_IDENTITY", file_identity(paths$gene_set_resource))
    resources <- stage("RESOURCE_READ", readRDS(paths$gene_set_resource))
    stage("MEMBER_PREGUARD", {
      checked <- member_guard(resources)
      need(identical(checked$schema, "GSE128623_MODEL_MEMBER_CONTENT_GUARD/v1") &&
             identical(checked$status, "MEMBER_CONTENT_ONLY_PASS") &&
             identical(checked$resources, resources) &&
             identical(checked$input_modified, FALSE), "MEMBER_PREGUARD")
      need(identical(file_identity(paths$gene_set_resource), resource_pin), "RESOURCE_DRIFT")
    })
    written <- stage("ANALYZER_EXECUTION", analyzer$v2a_execute_model(
      analyzer_input, paths$gene_set_resource, paths$staging,
      nrot = 9999L, minimum_members = 10L
    ))
    stage("ANALYZER_OUTPUT_SET", {
      need(is.character(written) && identical(written, expected_outputs), "ANALYZER_OUTPUT_SET")
      need(identical(file_identity(paths$gene_set_resource), resource_pin), "RESOURCE_DRIFT")
    })
    comparison_pins <- stage("COMPARISON_BASELINE", lapply(comparison_files, file_identity))
    report <- stage("INDEPENDENT_VALIDATION", validator$v2v_validate_output(
      validator_input, paths$gene_set_resource, paths$staging,
      nrot = 9999L, minimum_members = 10L
    ))
    stage("VALIDATOR_REPORT", {
      fields <- c("candidate_id", "status", "primary_edge_key_count", "complete29_edge_key_count",
                  "primary_retained_feature_count", "complete29_retained_feature_count",
                  "edgeR_ci_status", "batch_adjustment_performed", "combat_performed",
                  "formal_sequential_inference", "candidate_generation_allowed", "s2_status")
      need(is.data.frame(report) && nrow(report) == 1L && identical(names(report), fields),
           "VALIDATOR_REPORT_SCHEMA")
      need(identical(report$candidate_id, "GSE128623_exact30_linux_model_v2") &&
             identical(report$status, "PASS_INDEPENDENT_RECONSTRUCTION") &&
             identical(report$edgeR_ci_status, "NOT_AVAILABLE_FROM_STANDARD_EDGER_QLF_OUTPUT") &&
             identical(report$s2_status, "NOT_IDENTIFIABLE_NO_S1_INPUT"), "VALIDATOR_REPORT_VALUE")
      for (name in fields[8:11]) need(identical(report[[name]], FALSE), "VALIDATOR_REPORT_FLAGS")
      for (name in fields[3:6]) {
        value <- report[[name]]
        need(typeof(value) %in% c("integer", "double") && !is.object(value) &&
               length(value) == 1L && is.finite(value) && value >= 2 && value == floor(value),
             "VALIDATOR_REPORT_COUNT")
      }
      need(report$primary_edge_key_count == 7 * report$primary_retained_feature_count &&
             report$complete29_edge_key_count == 7 * report$complete29_retained_feature_count,
           "VALIDATOR_REPORT_COUNT_RELATION")
      need(identical(file_identity(paths$gene_set_resource), resource_pin) &&
             identical(lapply(comparison_files, file_identity), comparison_pins), "POST_VALIDATOR_DRIFT")
    })
    stage("COMPARISON_CONTENT", {
      # These saved edge tables were just independently reconstructed/verified by
      # V. This is transitive content association, not a new independent fit or
      # authentication of V's private reconstructed objects. Native held-read and
      # output custody remain required even with these before/after whole hashes.
      primary <- as.data.frame(read_table(comparison_files[[1L]]))[
        c("gene_id", "contrast", "logFC", "BH_q")
      ]
      complete29 <- as.data.frame(read_table(comparison_files[[2L]]))[
        c("gene_id", "contrast", "logFC", "BH_q")
      ]
      checked <- comparison_guard(
        primary, complete29, read_table(comparison_files[[3L]]),
        read_table(comparison_files[[4L]]), nrow(analyzer_input$counts)
      )
      need(identical(checked$status, "PASS_BRANCH_COMPARISON_CONTENT_ONLY") &&
             identical(checked$comparison_numeric_checked, TRUE) &&
             identical(checked$summary_numeric_checked, TRUE) &&
             identical(checked$missingness_checked, TRUE), "COMPARISON_CONTENT")
      need(identical(lapply(comparison_files, file_identity), comparison_pins) &&
             identical(file_identity(paths$gene_set_resource), resource_pin), "FINAL_DRIFT")
    })
    list(
      schema = "GSE128623_MODEL_EXECUTION_DRIVER/v1",
      status = "ORIGINAL_MODEL_SEQUENCE_RETURNED_CONTENT_VALIDATION_ONLY",
      analyzer_returned = TRUE, independent_validator_returned = TRUE,
      comparison_content_checked = TRUE, model_run = TRUE,
      native_runtime_checks_returned = TRUE, s1_zero_checks_returned = TRUE,
      primary_sample_n = 30L, sensitivity_sample_n = 29L,
      rotation_count = 9999L, minimum_members = 10L,
      comparison_scope = "V_VALIDATED_SAVED_EDGE_TABLES_TRANSITIVE_ASSOCIATION",
      production_go = FALSE, scientific_acceptance = FALSE,
      closure_authenticated = FALSE, publication_authorized = FALSE
    )
  }
  execute_sequence()
}
