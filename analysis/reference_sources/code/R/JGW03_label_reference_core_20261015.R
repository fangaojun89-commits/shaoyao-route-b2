# Functions only; no I/O driver, real-data admission, observed bridge or fits
# execute when this file is sourced. The old report-only authority is unchanged.

jgw03_require <- function(value, label) {
  if (!isTRUE(value)) stop(label, call. = FALSE)
}

jgw03_label_contract <- function() {
  list(status = "PREPARATION_ONLY_REAL_EXECUTION_NOT_AUTHORIZED_HERE",
    native_source_sha256 = "5df6378ed4dadfb4050d91e252a63af4f44d2be80d3070d889d1e0e007796775",
    equivalent_rotation_source_sha256 = "ea666b1e5aac06eea431e5f376148755bb47ae02c3b467ef56ba71e9b79d3577",
    seed = 20260817L, rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    group_levels = c("vehicle", "low", "medium", "high"),
    libraries = 12L, group_capacity = rep(3L, 4L), full_assignment_n = 369600L,
    full_nonidentity_n = 369599L, MC_draws = 10000L,
    MC_rule = "(1+b)/(10000+1), observed separate and once",
    high_companion_n = 20L,
    high_companion_rule = "inclusive b/20, complete orbit includes observed once; high contrast only",
    mroast_nrot = 1999L, fgsea_eps = 0,
    native_q_rule = "q<=0.05 within the original collection/method/contrast family",
    source_only_permissions = "Definitions and invented fixtures are not permission for real fitting")
}

jgw03_check_allocation <- function(allocation, sample_ids) {
  jgw03_require(is.numeric(allocation) && length(allocation) == 12L &&
    !anyNA(allocation) && all(is.finite(allocation)) && all(allocation %in% 0:3) &&
    length(sample_ids) == 12L && !anyNA(sample_ids) && all(nzchar(sample_ids)) &&
    !anyDuplicated(sample_ids) && identical(names(allocation), sample_ids) &&
    identical(as.integer(tabulate(allocation + 1L, nbins = 4L)), rep(3L, 4L)),
    "Expected ordered, named 12-library allocation with four groups of three")
  invisible(TRUE)
}

jgw03_relabel_model_metadata <- function(samples, allocation, accession) {
  contract <- jgw03_label_contract()
  jgw03_require(accession %in% c("GSE289929", "GSE227494", "GSE299063") &&
    is.data.frame(samples) && nrow(samples) == 12L &&
    all(c("sample_id", "batch") %in% names(samples)) &&
    !anyNA(samples$batch) && length(unique(as.character(samples$batch))) == 1L,
    "One accession-specific 12-library single recorded-batch subset required")
  jgw03_check_allocation(allocation, as.character(samples$sample_id))
  doses <- if (accession == "GSE289929") c(0, 4, 20, 100) else c(0, 20, 100, 500)
  # Original metadata is copied, not rewritten. Only derived model-label columns
  # in this RAM object change; replicate numbers are not inferred paired blocks.
  model_metadata <- samples
  model_metadata$dose_group <- factor(contract$group_levels[allocation + 1L],
    levels = contract$group_levels)
  model_metadata$dose_numeric <- doses[allocation + 1L]
  model_metadata$ordered_dose <- as.numeric(model_metadata$dose_group) - 1
  model_metadata$log10_dose_plus_1 <- log10(model_metadata$dose_numeric + 1)
  model_metadata
}

jgw03_designs <- function(samples) {
  jgw03_require(is.factor(samples$dose_group) &&
    identical(levels(samples$dose_group), jgw03_label_contract()$group_levels),
    "Native dose factor levels required")
  high <- stats::model.matrix(~ 0 + dose_group, data = samples)
  colnames(high) <- levels(samples$dose_group)
  trend <- stats::model.matrix(~ ordered_dose, data = samples)
  logtrend <- stats::model.matrix(~ log10_dose_plus_1, data = samples)
  contrast_high <- limma::makeContrasts(high - vehicle, levels = high)
  jgw03_require(qr(high)$rank == 4L && qr(trend)$rank == 2L && qr(logtrend)$rank == 2L,
    "All native designs must retain full rank")
  list(high = high, trend = trend, logtrend = logtrend, contrast_high = contrast_high)
}

jgw03_allocation_registry <- function(observed) {
  contract <- jgw03_label_contract()
  jgw03_check_allocation(observed, names(observed))
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  do.call(RNGkind, as.list(contract$rng_kind))
  set.seed(contract$seed)
  initial_rng <- get(".Random.seed", envir = .GlobalEnv)
  allocations <- matrix(NA_integer_, contract$MC_draws + 1L, 12L,
    dimnames = list(as.character(seq_len(contract$MC_draws + 1L)), names(observed)))
  allocations[1L, ] <- as.integer(observed)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  identity_key <- paste(observed, collapse = "")
  assign(identity_key, TRUE, envir = seen)
  attempts <- rejected_identity <- rejected_duplicate <- 0L
  accepted <- 0L
  while (accepted < contract$MC_draws) {
    # Every labelled assignment has exactly (3!)^4 preimages under a uniform
    # permutation of these 12 positions. Rejection preserves uniform distinct
    # nonidentity draws; no balanced or favorable subset is selected.
    proposal <- as.integer(observed[sample.int(12L, 12L, replace = FALSE)])
    key <- paste(proposal, collapse = "")
    attempts <- attempts + 1L
    if (identical(key, identity_key)) {
      rejected_identity <- rejected_identity + 1L
    } else if (exists(key, envir = seen, inherits = FALSE)) {
      rejected_duplicate <- rejected_duplicate + 1L
    } else {
      accepted <- accepted + 1L
      allocations[accepted + 1L, ] <- proposal
      assign(key, TRUE, envir = seen)
    }
  }
  target <- which(observed %in% c(0L, 3L))
  companions <- vapply(utils::combn(target, 3L, simplify = FALSE), function(high_positions) {
    value <- as.integer(observed)
    value[target] <- 0L
    value[high_positions] <- 3L
    value
  }, integer(12L))
  companions <- t(companions)
  identity <- which(apply(companions, 1L, function(row) identical(row, as.integer(observed))))
  jgw03_require(length(identity) == 1L && nrow(companions) == 20L,
    "Full high/vehicle companion orbit must include identity exactly once")
  companions <- companions[c(identity, setdiff(seq_len(20L), identity)), , drop = FALSE]
  dimnames(companions) <- list(as.character(seq_len(20L)), names(observed))
  list(MC = allocations, high_companion = companions, attempts = attempts,
    rejected_identity = rejected_identity, rejected_duplicate = rejected_duplicate,
    seed = contract$seed, rng_kind = RNGkind(), initial_rng = initial_rng,
    final_rng = get(".Random.seed", envir = .GlobalEnv), full_space_n = 369600L,
    source_unit = "library observation; randomization and paired culture blocks not authenticated",
    MC_nonidentity_draw_n = 10000L, identity_row = 1L,
    MC_tail_convention = contract$MC_rule, companion_tail_convention = contract$high_companion_rule)
}

jgw03_family_registry <- function() {
  rows <- list()
  add <- function(accession, method, collection, contrasts) {
    rows[[length(rows) + 1L]] <<- data.frame(accession = accession, method = method,
      collection = collection, contrast = contrasts, space = "four_group_MC",
      null_draw_n = 10000L, observed_once = TRUE)
  }
  high <- "high_vs_vehicle"
  two <- c(high, "ordered_dose_per_step")
  three <- c(two, "log10_dose_plus_1")
  for (accession in c("GSE289929", "GSE227494", "GSE299063")) {
    add(accession, "edgeR_QL", "genes", three)
    add(accession, "limma_voom_trend_TRUE", "genes", two)
    add(accession, "DESeq2_native_IF_Cooks", "genes", two)
    for (collection in c("Reactome97", "QuickGO_SynGO_sensitivity")) {
      add(accession, "CAMERA", collection, three)
      add(accession, "mean_z_limma", collection, three)
      add(accession, "fgseaMultilevel", collection, two)
    }
    for (collection in c("PROGENy_top500_ULM", "DoRothEA_AB_ULM"))
      add(accession, "ULM_activity_limma", collection, three)
  }
  add("GSE289929", "mroast1999_equivalent_not_historical_authentication", "Reactome97", two)
  global <- do.call(rbind, rows)
  companions <- global[global$contrast == high, , drop = FALSE]
  companions$space <- "high_companion_exact"
  companions$null_draw_n <- 20L
  answer <- rbind(global, companions)
  rownames(answer) <- NULL
  jgw03_require(nrow(global) == 89L && nrow(companions) == 34L &&
    !anyDuplicated(answer[c("accession", "method", "collection", "contrast", "space")]),
    "Complete existing-native scope registry differs")
  answer
}

jgw03_native_definitions <- function(source_path) {
  contract <- jgw03_label_contract()
  jgw03_require(digest::digest(file = source_path, algo = "sha256", serialize = FALSE) ==
    contract$native_source_sha256, "Native canonical source pin differs")
  wanted <- c("stopf", "assert_true", "collapse_expression_by_symbol", "make_program_scores",
    "fit_activity_matrix", "fit_program_scores", "run_decoupler_activity",
    "limma_voom_table", "deseq_table", "run_camera_collection")
  syntax <- as.list(parse(source_path))
  environment <- new.env(parent = asNamespace("stats"))
  for (name in wanted) {
    definitions <- Filter(function(expression) is.call(expression) &&
      identical(expression[[1L]], as.name("<-")) && identical(expression[[2L]], as.name(name)), syntax)
    jgw03_require(length(definitions) == 1L && is.call(definitions[[1L]][[3L]]) &&
      identical(definitions[[1L]][[3L]][[1L]], as.name("function")),
      "One pure native function definition required")
    eval(definitions[[1L]], envir = environment)
  }
  # This callback captures native exports in RAM. It never writes to the old
  # result directory or invokes the whole driver or its unsafe resume path.
  environment$exports <- new.env(parent = emptyenv())
  environment$write_tsv <- function(value, path) {
    key <- basename(path)
    jgw03_require(!exists(key, envir = environment$exports, inherits = FALSE),
      "Repeated RAM-native export key")
    assign(key, value, envir = environment$exports)
    invisible(NULL)
  }
  environment$seed <- contract$seed
  environment
}

jgw03_capture_native_locals <- function(native_function, arguments, local_names) {
  # Retain fit objects without evaluating the kernel twice or altering any
  # instruction of the native body. These two selected source functions have
  # no explicit return(), so the added final capture is well-defined.
  jgw03_require(is.function(native_function) &&
    !"return" %in% all.names(body(native_function), functions = TRUE),
    "Capture supports the reviewed ordinary native function bodies only")
  original_body <- body(native_function)
  wrapped <- native_function
  body(wrapped) <- substitute({
    native_value <- ORIGINAL
    list(value = native_value,
      native_locals = mget(LOCALS, envir = environment(), inherits = FALSE))
  }, list(ORIGINAL = original_body, LOCALS = local_names))
  do.call(wrapped, arguments)
}

jgw03_build_fixed_layer <- function(counts, annotation, samples, accession,
                                     resources, native_source_path) {
  # Future approved execution only. No saved label-dependent gene fit/weights
  # or dispersion is accepted here; no scientific input is read by this module.
  native <- jgw03_native_definitions(native_source_path)
  observed <- setNames(as.numeric(samples$ordered_dose), as.character(samples$sample_id))
  model_metadata <- jgw03_relabel_model_metadata(samples, observed, accession)
  jgw03_require(is.matrix(counts) && is.numeric(counts) &&
    identical(colnames(counts), as.character(samples$sample_id)) &&
    !is.null(rownames(counts)) && !anyNA(rownames(counts)) &&
    !anyDuplicated(rownames(counts)) && all(is.finite(counts)) && all(counts >= 0) &&
    is.data.frame(annotation) &&
    all(c("raw_gene_id", "gene_symbol", "human_symbol") %in% names(annotation)) &&
    identical(rownames(counts), as.character(annotation$raw_gene_id)),
    "Full fractional expected-count/annotation/sample keys differ")
  designs <- jgw03_designs(model_metadata)
  y0 <- edgeR::DGEList(counts = counts, genes = annotation)
  keep <- edgeR::filterByExpr(y0, design = designs$high)
  jgw03_require(sum(keep) >= 1000L, "Native minimum expressed-gene population required")
  y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  retained_annotation <- annotation[keep, , drop = FALSE]
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  human <- native$collapse_expression_by_symbol(logcpm, retained_annotation, "human_symbol")
  species <- native$collapse_expression_by_symbol(logcpm, retained_annotation, "gene_symbol")
  y_human <- y[human$selected_index, , keep.lib.sizes = FALSE]
  rownames(y_human$counts) <- rownames(human$expression)
  y_human$genes <- data.frame(human_symbol = rownames(y_human))
  y_human <- edgeR::normLibSizes(y_human, method = "TMM")
  jgw03_require(identical(names(resources$collections),
    c("Reactome97", "QuickGO_SynGO_sensitivity")), "Both frozen collection definitions required")
  # The ULM score calculation is label-independent. Suppress only its following
  # dose-fit callback here, in this private imported-function environment.
  native_fit_activity <- native$fit_activity_matrix
  native$fit_activity_matrix <- function(...) NULL
  on.exit({ native$fit_activity_matrix <- native_fit_activity }, add = TRUE)
  progeny_expression <- if (accession == "GSE289929") human$expression else species$expression
  progeny <- native$run_decoupler_activity(progeny_expression, resources$progeny,
    "pathway", "gene_symbol", "weight", model_metadata, "PROGENy_top500_ULM", minsize = 50L)
  dorothea <- native$run_decoupler_activity(species$expression, resources$dorothea,
    "tf", "target", "mor", model_metadata, "DoRothEA_AB_ULM", minsize = 5L)
  list(accession = accession, original_samples = samples, observed_labels = observed,
    full_counts = counts, full_annotation = annotation, keep = keep, y = y,
    annotation = retained_annotation, expr_human = human$expression,
    expr_native = species$expression, human_mapping = human$mapping,
    native_mapping = species$mapping, y_human_unfitted = y_human,
    progeny_activity = progeny$activity, dorothea_activity = dorothea$activity,
    resources = resources, native_source_sha256 = jgw03_label_contract()$native_source_sha256,
    cached_gene_fit_or_weight = FALSE)
}

jgw03_fit_one_allocation <- function(fixed, allocation, native_source_path) {
  # Callable kernel definition, NOT a launcher or execution-authority check.
  # Its stochastic path is the native path; no unvalidated canonical swap or
  # historical-PC12-producer authentication is claimed by this preparation.
  jgw03_require(identical(fixed$native_source_sha256,
    jgw03_label_contract()$native_source_sha256) &&
    identical(fixed$cached_gene_fit_or_weight, FALSE),
    "A source-bound label-independent layer, not an old gene-fit cache, is required")
  native <- jgw03_native_definitions(native_source_path)
  samples <- jgw03_relabel_model_metadata(fixed$original_samples, allocation, fixed$accession)
  designs <- jgw03_designs(samples)
  tests <- fits <- list()
  for (label in c("high", "trend", "logtrend")) {
    y <- edgeR::estimateDisp(fixed$y, designs[[label]], robust = TRUE)
    fits[[label]] <- edgeR::glmQLFit(y, designs[[label]], robust = TRUE)
    tests[[label]] <- if (label == "high") edgeR::glmQLFTest(fits[[label]],
      contrast = designs$contrast_high) else edgeR::glmQLFTest(fits[[label]], coef = 2L)
  }
  labels <- c(high = "high_vs_vehicle", trend = "ordered_dose_per_step", logtrend = "log10_dose_plus_1")
  gene_edgeR <- lapply(names(tests), function(label) {
    tab <- edgeR::topTags(tests[[label]], n = Inf, sort.by = "none")$table
    jgw03_require(identical(rownames(tab), fixed$annotation$raw_gene_id), "Native edgeR gene order differs")
    # No ql_ci/QLF-reconstructed interval is imported or regenerated.
    answer <- cbind(fixed$annotation, contrast = labels[[label]], tab)
    answer$FDR <- stats::p.adjust(answer$PValue, method = "BH")
    answer
  })
  names(gene_edgeR) <- names(tests)
  voom <- list(
    high = jgw03_capture_native_locals(native$limma_voom_table,
      list(y = fixed$y, design = designs$high, contrast = designs$contrast_high,
        annotation = fixed$annotation, contrast_name = labels[["high"]]), c("v", "fit", "tab")),
    trend = jgw03_capture_native_locals(native$limma_voom_table,
      list(y = fixed$y, design = designs$trend, coef = 2L,
        annotation = fixed$annotation, contrast_name = labels[["trend"]]), c("v", "fit", "tab")))
  deseq <- lapply(c("high", "trend"), function(label) {
    jgw03_capture_native_locals(native$deseq_table,
      list(counts = fixed$full_counts, keep = fixed$keep, samples = samples,
        design_type = if (label == "high") "categorical" else "trend",
        annotation = fixed$full_annotation, contrast_name = labels[[label]]), c("dds", "res"))
  })
  names(deseq) <- c("high", "trend")
  common_human <- edgeR::estimateDisp(fixed$y_human_unfitted, designs$high, robust = TRUE)
  human_tests <- lapply(c("high", "trend"), function(label) {
    fit <- edgeR::glmQLFit(edgeR::estimateDisp(common_human, designs[[label]], robust = TRUE),
      designs[[label]], robust = TRUE)
    if (label == "high") edgeR::glmQLFTest(fit, contrast = designs$contrast_high) else
      edgeR::glmQLFTest(fit, coef = 2L)
  })
  names(human_tests) <- c("high", "trend")
  programs <- lapply(names(fixed$resources$collections), function(collection) {
    resource <- fixed$resources$collections[[collection]]
    native$run_camera_collection(fixed$accession, collection, common_human,
      fixed$expr_human, resource$sets, resource$term_meta, samples, designs,
      human_tests, result_dir = "RAM_NATIVE_EXPORTS_ONLY",
      run_roast = fixed$accession == "GSE289929" && collection == "Reactome97")
  })
  names(programs) <- names(fixed$resources$collections)
  activities <- list(
    PROGENy_top500_ULM = native$fit_activity_matrix(fixed$progeny_activity, samples, "PROGENy_top500_ULM"),
    DoRothEA_AB_ULM = native$fit_activity_matrix(fixed$dorothea_activity, samples, "DoRothEA_AB_ULM"))
  list(accession = fixed$accession, labels = allocation, model_metadata = samples,
    designs = designs, gene_edgeR = gene_edgeR, native_edgeR_fits = fits,
    native_edgeR_tests = tests, gene_voom = voom, gene_DESeq2 = deseq,
    programs = programs, activities = activities,
    native_exports = as.list(native$exports, all.names = TRUE),
    scope = "FUTURE_APPROVED_SAME_DATA_CONDITIONAL_REFERENCE_NOT_CANDIDATE_ADMISSION",
    native_source_sha256 = jgw03_label_contract()$native_source_sha256,
    copied_old_gene_fit_or_weight = FALSE, candidate_gate_eligible = FALSE,
    real_reference_admission_or_outer_acceptance_implemented = FALSE,
    stochastic_swap_canonicalization_accepted = FALSE,
    historical_PC12_roast_producer_authenticated = FALSE)
}
