# Input-derived reference for the twelve non-gene scientific-core tables.
# New orchestration; no primary-core fitting function or reported result is used.
# Shared input-contract helpers and package implementations are NOT independent
# personnel or independent mathematical implementations. No production authority.
# This API sets/advances RNG; invoke in a fresh process with recorded backend.

prepare_gse205906_program_validation <- function(project_root, count_matrix,
    gene_annotation, sample_ledger, gene_resources, orthology_table,
    domain = "SYNTHETIC") {
  fail <- function(message) stop(message, call. = FALSE)
  require_true <- function(value, message) if (!isTRUE(value)) fail(message)
  require_true(is.character(domain) && length(domain) == 1L && !is.na(domain) &&
    domain %in% c("SYNTHETIC", "PRODUCTION"), "PROGRAM_REFERENCE_DOMAIN")
  pins <- c(
    "code/R/gse205906_bulk_model_v2_r4_lib.R" =
      "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3",
    "config/GSE205906_bulk_model_v2_r4.json" =
      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"
  )
  capture <- function(relative, expected) {
    path <- file.path(project_root, relative)
    info <- file.info(path)
    require_true(nrow(info) == 1L && !is.na(info$size) && !info$isdir &&
      info$size < 1000000 && !nzchar(Sys.readlink(path)), "PROGRAM_REFERENCE_SOURCE")
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    payload <- readBin(con, "raw", n = as.integer(info$size) + 1L)
    require_true(length(payload) == info$size &&
      identical(digest::digest(payload, "sha256", serialize = FALSE), expected),
      "PROGRAM_REFERENCE_SOURCE_IDENTITY")
    rawToChar(payload)
  }
  sources <- Map(capture, names(pins), unname(pins))
  cfg <- jsonlite::fromJSON(sources[[2L]], simplifyVector = FALSE)
  require_true(identical(as.character(getRversion()), cfg$software$R), "PROGRAM_REFERENCE_R")
  packages <- c("digest", "jsonlite", "data.table", "edgeR", "limma", "fgsea", "decoupleR")
  versions <- vapply(packages, function(package) {
    observed <- as.character(utils::packageVersion(package))
    require_true(identical(observed, cfg$software[[package]]), paste0("PROGRAM_REFERENCE_PACKAGE:", package))
    observed
  }, character(1))
  # Only frozen input-validation function declarations are imported.
  helpers <- new.env(parent = baseenv())
  anchors <- list(
    c("`%||%` <- function", "v2r4_child_request_contract <- function"),
    c("v2r4_validate_ledger <- function", "v2r4_release_conditional_inclusion <- function"),
    c("v2r4_validate_counts_original <- function", "v2r4_edgeR_ql_covariance <- function")
  )
  for (anchor in anchors) {
    starts <- gregexpr(anchor[[1L]], sources[[1L]], fixed = TRUE)[[1L]]
    ends <- gregexpr(anchor[[2L]], sources[[1L]], fixed = TRUE)[[1L]]
    require_true(length(starts) == 1L && length(ends) == 1L && starts > 0 && ends > starts,
      "PROGRAM_REFERENCE_SOURCE_SPAN")
    block <- parse(text = substr(sources[[1L]], starts, ends - 1L))
    for (expression in block) {
      require_true(is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function")),
        "PROGRAM_REFERENCE_HELPER_NOT_DECLARATION")
      eval(expression, helpers)
    }
  }
  lockEnvironment(helpers, bindings = TRUE)
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  ledger <- helpers$v2r4_validate_ledger(data.table::copy(sample_ledger), cfg)
  counts <- helpers$v2r4_validate_counts_original(count_matrix, ids)
  require_true(!is.null(rownames(counts)) && !anyDuplicated(rownames(counts)), "PROGRAM_REFERENCE_GENE_IDS")
  if (domain == "PRODUCTION") require_true(
    nrow(counts) == cfg$input_snapshot$expected_gene_n_before_filter, "PROGRAM_REFERENCE_PRODUCTION_GENE_N")
  annotation_contract <- helpers$v2r4_validate_annotation_gene(data.table::copy(gene_annotation),
    rownames(counts), cfg$input_snapshot$minimum_nonblank_symbol_coverage)
  annotation <- annotation_contract$annotation
  group <- factor(unlist(cfg$design$groups_in_order, use.names = FALSE),
    levels = unlist(cfg$design$group_levels, use.names = FALSE))
  make_design <- function(indices) {
    selected_groups <- group[indices]
    value <- stats::model.matrix(~ 0 + selected_groups)
    colnames(value) <- unlist(cfg$design$coefficient_order, use.names = FALSE)
    rownames(value) <- ids[indices]
    require_true(qr(value)$rank == cfg$design$expected_rank, "PROGRAM_REFERENCE_DESIGN_RANK")
    value
  }
  design <- make_design(seq_along(ids))
  contrasts <- vapply(cfg$design$contrasts, function(item) unlist(item$vector, use.names = FALSE), numeric(3))
  rownames(contrasts) <- colnames(design)
  settings <- if (domain == "SYNTHETIC") cfg$synthetic_contract else cfg$analysis
  minimum <- if (domain == "SYNTHETIC") settings$minimum_retained_genes else settings$minimum_retained_genes_production
  domains <- settings$program_domains
  seed <- as.integer(cfg$analysis$random_seed)
  set.seed(seed)
  size_ok <- function(n) n >= cfg$analysis$program_size_min & n <= cfg$analysis$program_size_max
  add_flags <- function(value) {
    value <- data.table::as.data.table(value)
    value$analysis_role <- rep(cfg$gate$analysis_role, nrow(value))
    for (field in c("candidate_gate_eligible", "can_generate_or_upgrade_candidate", "formal_sequential_inference")) {
      value[[field]] <- rep(FALSE, nrow(value))
    }
    value
  }
  preprocess <- function(indices) {
    matrix_counts <- counts[, indices, drop = FALSE]
    model <- make_design(indices)
    initial <- edgeR::DGEList(counts = matrix_counts)
    retained_expression <- edgeR::filterByExpr(initial, design = model)
    retained_cpm <- rowSums(edgeR::cpm(initial, log = FALSE) >= cfg$analysis$minimum_cpm) >=
      cfg$analysis$minimum_samples_at_cpm
    retained <- retained_expression & retained_cpm
    require_true(sum(retained) >= minimum, "PROGRAM_REFERENCE_FILTERED_GENE_N")
    normalized <- edgeR::normLibSizes(initial[retained, , keep.lib.sizes = FALSE], method = "TMM")
    normalized <- edgeR::estimateDisp(normalized, model, robust = TRUE)
    # QL fitting is not needed for these outputs; voom consumes the same y.
    transformed <- limma::voom(normalized, design = model, plot = FALSE)
    list(v = transformed, design = model, annotation = annotation[retained],
      keep = retained, filter_by_expr = retained_expression, cpm = retained_cpm)
  }
  resources <- data.table::copy(gene_resources)
  require_true(is.list(resources) && all(c("reactome", "syngo_sensitivity", "progeny_human", "dorothea_rat") %in%
    names(resources)), "PROGRAM_REFERENCE_RESOURCES")
  orthology <- data.table::as.data.table(data.table::copy(orthology_table))
  require_true(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in% names(orthology)),
    "PROGRAM_REFERENCE_ORTHOLOGY_SCHEMA")
  map <- unique(orthology[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
    .(human_symbol = as.character(human_symbol), rat_symbol = as.character(other_symbol))])
  require_true(!anyNA(map) && !anyDuplicated(map$human_symbol) && !anyDuplicated(map$rat_symbol),
    "PROGRAM_REFERENCE_ORTHOLOGY_ONE_TO_ONE")
  unsigned_maps <- lapply(c(Reactome = "reactome", SynGO_sensitivity = "syngo_sensitivity"), function(name) {
    source <- data.table::as.data.table(resources[[name]])
    require_true(all(c("term_id", "term_name", "gene_symbol") %in% names(source)), "PROGRAM_REFERENCE_UNSIGNED_SCHEMA")
    source <- unique(source[, .(program_id = as.character(term_id), program_name = as.character(term_name),
      human_symbol = as.character(gene_symbol))])
    unique(merge(source, map, by = "human_symbol", allow.cartesian = TRUE)[,
      .(program_id, program_name, symbol = rat_symbol)])
  })
  progeny <- data.table::as.data.table(resources$progeny_human)
  require_true(all(c("pathway", "gene_symbol", "weight") %in% names(progeny)), "PROGRAM_REFERENCE_PROGENY_SCHEMA")
  progeny <- unique(merge(progeny[, .(source = as.character(pathway), human_symbol = as.character(gene_symbol),
    mor = as.numeric(weight))], map, by = "human_symbol", allow.cartesian = TRUE)[,
    .(source, target = rat_symbol, mor)])
  dorothea <- data.table::as.data.table(resources$dorothea_rat)
  require_true(all(c("tf", "target", "mor", "confidence") %in% names(dorothea)), "PROGRAM_REFERENCE_DOROTHEA_SCHEMA")
  dorothea <- unique(dorothea[, .(source = as.character(tf), target = as.character(target),
    mor = as.numeric(mor), confidence = as.character(confidence))])
  require_true(all(dorothea$confidence %in% c("A", "B")) && !anyNA(dorothea$confidence),
    "PROGRAM_REFERENCE_DOROTHEA_AB")
  require_true(data.table::uniqueN(progeny$source) == domains$PROGENy$source_program_n &&
    setequal(progeny$source, unlist(domains$PROGENy$exact_program_ids, use.names = FALSE)) &&
    data.table::uniqueN(dorothea$source) == domains$DoRothEA_AB$frozen_source_program_upper_bound,
    "PROGRAM_REFERENCE_SOURCE_DOMAINS")
  networks <- list(PROGENy = progeny, DoRothEA_AB = dorothea)
  base <- preprocess(seq_along(ids))
  expected <- list()
  unsigned_tables <- list(camera = list(), mroast = list(), fgsea = list())
  coverage_tables <- list()
  for (collection_index in seq_along(unsigned_maps)) {
    collection_name <- names(unsigned_maps)[[collection_index]]
    mapped <- unsigned_maps[[collection_index]]
    members <- lapply(split(mapped$symbol, mapped$program_id), function(x) unique(intersect(x, base$annotation$symbol)))
    members <- members[size_ok(lengths(members))]
    indices <- lapply(members, function(x) which(base$annotation$symbol %in% x))
    coverage <- mapped[, .(mapped_symbol_n = data.table::uniqueN(symbol)), by = .(program_id, program_name)]
    coverage$tested_feature_n <- vapply(coverage$program_id, function(id) {
      length(intersect(unique(mapped$symbol[mapped$program_id == id]), base$annotation$symbol))
    }, integer(1))
    coverage$collection <- collection_name
    coverage$source_feature_n <- coverage$mapped_symbol_n
    coverage$eligible_10_500 <- size_ok(coverage$tested_feature_n)
    coverage$tested <- coverage$program_id %in% names(members)
    coverage$confidence_domain <- "NOT_APPLICABLE"
    coverage$mapping_rule <- cfg$analysis$orthology_authority
    coverage_tables[[collection_name]] <- coverage
    if (!length(indices)) next
    for (j in seq_len(ncol(contrasts))) {
      camera <- data.table::as.data.table(limma::camera(base$v$E, indices, base$design,
        contrast = contrasts[, j], inter.gene.cor = cfg$analysis$camera_inter_gene_correlation),
        keep.rownames = "program_id")
      set.seed(seed + j + collection_index * 100L)
      mroast <- data.table::as.data.table(limma::mroast(base$v$E, indices, base$design,
        contrast = contrasts[, j], nrot = settings$mroast_rotations, set.statistic = "mean"),
        keep.rownames = "program_id")
      rank_fit <- limma::lmFit(base$v, base$design)
      rank_fit <- limma::eBayes(limma::contrasts.fit(rank_fit, contrasts[, j, drop = FALSE]), robust = TRUE)
      rank_table <- data.table::data.table(symbol = base$annotation$symbol, stat = as.numeric(rank_fit$t[, 1L]))
      rank_table <- rank_table[is.finite(stat) & !is.na(symbol) & nzchar(symbol)]
      rank_table <- rank_table[order(-abs(stat)), .SD[1L], by = symbol]
      ranks <- stats::setNames(rank_table$stat, rank_table$symbol)
      set.seed(seed + 1000L + j + collection_index * 100L)
      fgsea <- data.table::as.data.table(fgsea::fgseaMultilevel(members, ranks,
        minSize = cfg$analysis$program_size_min, maxSize = cfg$analysis$program_size_max,
        eps = cfg$analysis$fgsea_eps))
      data.table::setnames(fgsea, "pathway", "program_id")
      for (method in c("camera", "mroast", "fgsea")) {
        value <- switch(method, camera = camera, mroast = mroast, fgsea = fgsea)
        value$collection <- collection_name
        value$contrast <- colnames(contrasts)[[j]]
        value$method <- switch(method, camera = "camera_inter_gene_cor_0.01",
          mroast = paste0("mroast_mean_", settings$mroast_rotations, "_rotations"), fgsea = "fgseaMultilevel")
        value$p_value <- if (method == "fgsea") value$pval else value$PValue
        value$q_within_collection_contrast <- if (method == "fgsea") value$padj else stats::p.adjust(value$p_value, "BH")
        unsigned_tables[[method]][[length(unsigned_tables[[method]]) + 1L]] <- value
      }
    }
  }
  for (method in names(unsigned_tables)) expected[[paste0("unsigned_", method)]] <-
    add_flags(data.table::rbindlist(unsigned_tables[[method]], fill = TRUE))

  activity <- function(preprocessed, enforce_full) {
    # Average duplicate gene symbols before scoring, preserving first-seen order.
    symbols <- unique(preprocessed$annotation$symbol)
    expression <- t(vapply(symbols, function(symbol) {
      colMeans(preprocessed$v$E[preprocessed$annotation$symbol == symbol, , drop = FALSE])
    }, numeric(ncol(preprocessed$v$E))))
    rownames(expression) <- symbols
    colnames(expression) <- colnames(preprocessed$v$E)
    fits <- list()
    coverages <- list()
    scored_ids <- list()
    for (collection_name in names(networks)) {
      net <- networks[[collection_name]]
      coverage <- net[, .(source_feature_n = data.table::uniqueN(target),
        tested_feature_n = data.table::uniqueN(intersect(target, rownames(expression)))), by = source]
      coverage$collection <- collection_name
      coverage$eligible_10_500 <- size_ok(coverage$tested_feature_n)
      coverage$tested <- coverage$eligible_10_500
      coverage$confidence_domain <- if (collection_name == "PROGENy") "NOT_APPLICABLE" else
        paste(sort(unique(dorothea$confidence)), collapse = ",")
      coverage$mapping_rule <- if (collection_name == "PROGENy") cfg$analysis$orthology_authority else
        "frozen_rat_AB_strict_one_to_one_network"
      if (enforce_full) require_true(sum(coverage$tested) == domains[[collection_name]]$eligible_tested_program_n,
        "PROGRAM_REFERENCE_TESTED_DOMAIN")
      coverages[[collection_name]] <- coverage
      selected <- net[source %in% coverage$source[coverage$tested] & target %in% rownames(expression)]
      if (!nrow(selected)) {
        fits[[collection_name]] <- data.table::data.table()
        scored_ids[[collection_name]] <- character()
        next
      }
      scores <- data.table::as.data.table(decoupleR::run_ulm(expression, selected[, .(source, target, mor)],
        .source = source, .target = target, .mor = mor, minsize = cfg$analysis$program_size_min))
      wide <- data.table::dcast(scores, source ~ condition, value.var = "score")
      score_matrix <- as.matrix(wide[, colnames(expression), with = FALSE])
      rownames(score_matrix) <- wide$source
      scored_ids[[collection_name]] <- rownames(score_matrix)
      model <- limma::lmFit(score_matrix, preprocessed$design)
      model <- limma::eBayes(limma::contrasts.fit(model, contrasts), robust = TRUE)
      fits[[collection_name]] <- data.table::rbindlist(lapply(seq_len(ncol(contrasts)), function(j) {
        error <- model$stdev.unscaled[, j] * sqrt(model$s2.post)
        interval <- stats::qt(0.975, model$df.total) * error
        data.table::data.table(collection = collection_name, program_id = rownames(score_matrix),
          contrast = colnames(contrasts)[[j]], method = "decoupleR_ULM_then_limma",
          effect = model$coefficients[, j], se = error,
          ci_low = model$coefficients[, j] - interval, ci_high = model$coefficients[, j] + interval,
          t = model$t[, j], df_total = model$df.total, p_value = model$p.value[, j],
          q_within_collection_contrast = stats::p.adjust(model$p.value[, j], "BH"))
      }))
    }
    list(fits = fits, coverage = coverages, ids = scored_ids)
  }
  full_activity <- activity(base, TRUE)
  expected$activity_results <- add_flags(data.table::rbindlist(full_activity$fits, fill = TRUE))
  expected$mapping_coverage <- add_flags(data.table::rbindlist(c(coverage_tables, full_activity$coverage), fill = TRUE))
  omission_tables <- list()
  for (omission in seq_along(ids)) {
    refit <- activity(preprocess(setdiff(seq_along(ids), omission)), FALSE)
    for (collection_name in names(full_activity$ids)) {
      grid <- data.table::CJ(program_id = full_activity$ids[[collection_name]], contrast = colnames(contrasts), unique = TRUE)
      grid$collection <- collection_name
      grid$omitted_animal <- ids[[omission]]
      grid$omitted_group <- as.character(group[[omission]])
      grid$loao_method <- "ACTIVITY_FULL_PIPELINE_REFIT"
      fitted <- refit$fits[[collection_name]]
      if (nrow(fitted)) {
        grid <- merge(grid, fitted[, .(program_id, contrast, estimate = effect, se, ci_low, ci_high,
          p_value, q_within_collection_contrast)], by = c("program_id", "contrast"), all.x = TRUE, sort = FALSE)
      } else for (field in c("estimate", "se", "ci_low", "ci_high", "p_value", "q_within_collection_contrast")) {
        grid[[field]] <- rep(NA_real_, nrow(grid))
      }
      grid$estimable <- is.finite(grid$estimate)
      grid$nonestimable_reason <- ifelse(grid$estimable, "NONE", "PROGRAM_INELIGIBLE_OR_NOT_SCORED_AFTER_FULL_REFIT")
      omission_tables[[length(omission_tables) + 1L]] <- grid
    }
  }
  expected$loao_long <- add_flags(data.table::rbindlist(omission_tables, fill = TRUE))
  summary <- expected$loao_long[, {
    values <- estimate[estimable & is.finite(estimate)]
    list(n_omissions = .N, n_estimable = length(values),
      estimate_min = if (length(values)) min(values) else NA_real_,
      estimate_max = if (length(values)) max(values) else NA_real_,
      estimate_mean = if (length(values)) mean(values) else NA_real_,
      sign_consistency = if (length(values)) max(mean(values > 0), mean(values < 0)) else NA_real_)
  }, by = .(collection, program_id, contrast)]
  minimum_omissions <- if (domain == "SYNTHETIC") settings$minimum_estimable_omissions_per_program_contrast else
    settings$loao$minimum_estimable_omissions_per_program_contrast
  require_true(all(summary$n_omissions == length(ids)) && all(summary$n_estimable >= minimum_omissions),
    "PROGRAM_REFERENCE_LOAO_ESTIMABILITY")
  summary$all_steps_refit <- TRUE
  summary$loao_method <- "ACTIVITY_FULL_PIPELINE_REFIT"
  expected$loao_summary <- add_flags(summary)
  expected$design_out <- add_flags(data.table::data.table(sample_id = ids, group = as.character(group),
    data.table::as.data.table(design)))
  expected$contrast_out <- add_flags(data.table::data.table(coefficient = rownames(contrasts),
    data.table::as.data.table(contrasts)))
  expected$filter_out <- add_flags(data.table::data.table(dataset_id = cfg$dataset_id,
    input_gene_n = nrow(counts), retained_gene_n = sum(base$keep), filterByExpr_retained_n = sum(base$filter_by_expr),
    explicit_CPM_retained_n = sum(base$cpm), minimum_cpm = cfg$analysis$minimum_cpm,
    minimum_samples = cfg$analysis$minimum_samples_at_cpm, filtering = cfg$analysis$edgeR$filter,
    normalization = cfg$analysis$edgeR$normalization, primary_method = cfg$analysis$primary_method))
  expected$annotation_metrics <- add_flags(annotation_contract$metrics)
  expected$gate <- add_flags(data.table::data.table(dataset_id = cfg$dataset_id, gate = "E2", upstream_gate = "E1",
    eligible_input_n = 0L, evaluated_candidate_n = 0L, pass_n = NA_integer_, status = cfg$gate$machine_status))
  # Comparator tolerances fixed before this module's first numerical test.
  # Applied to floating outputs only; identifiers, flags, counts and NA types exact.
  tolerance <- list(absolute = 1e-10, relative = 1e-8)
  keys <- list(unsigned_camera = c("collection", "program_id", "contrast"),
    unsigned_mroast = c("collection", "program_id", "contrast"),
    unsigned_fgsea = c("collection", "program_id", "contrast"),
    activity_results = c("collection", "program_id", "contrast"),
    mapping_coverage = c("collection", "program_id", "source", "program_name"),
    loao_long = c("collection", "program_id", "contrast", "omitted_animal"),
    loao_summary = c("collection", "program_id", "contrast"), design_out = "sample_id",
    contrast_out = "coefficient", filter_out = "dataset_id", annotation_metrics = "gene_id_source_column", gate = "dataset_id")
  snapshot <- data.table::copy(expected)
  object_hash <- function(value) digest::digest(serialize(value, NULL, version = 3), "sha256", serialize = FALSE)
  snapshot_identity <- object_hash(snapshot)
  input_identity <- object_hash(list(counts = counts, annotation = gene_annotation, ledger = sample_ledger,
    resources = gene_resources, orthology = orthology_table, domain = domain))
  compare <- function(result) {
    require_true(identical(snapshot_identity, object_hash(snapshot)), "PROGRAM_REFERENCE_CHANGED")
    require_true(is.list(result) && identical(result$domain, domain) && identical(result$production_go, FALSE) &&
      is.list(result$tables) && !anyDuplicated(names(result$tables)) && all(names(snapshot) %in% names(result$tables)),
      "PROGRAM_REFERENCE_RESULT_INTERFACE")
    for (name in names(snapshot)) {
      want <- data.table::copy(snapshot[[name]])
      got <- result$tables[[name]]
      require_true(is.data.frame(got) && !anyDuplicated(names(got)) && setequal(names(got), names(want)) &&
        nrow(got) == nrow(want), paste0("PROGRAM_REFERENCE_SCHEMA:", name))
      got <- data.table::as.data.table(data.table::copy(got))
      key <- intersect(keys[[name]], names(want))
      # r4 legitimately emits zero-row, flags-only unsigned tables when neither
      # collection has an eligible program. Exact schema/empty-column checks
      # still apply; absent biological keys cannot be required of that result.
      if (nrow(want)) {
        require_true(length(key) > 0L && !anyDuplicated(got[, key, with = FALSE]) &&
          !anyDuplicated(want[, key, with = FALSE]), paste0("PROGRAM_REFERENCE_KEYS:", name))
        data.table::setorderv(want, key)
        data.table::setorderv(got, key)
      }
      for (field in names(want)) {
        observed <- got[[field]]
        reference <- want[[field]]
        label <- paste("PROGRAM_REFERENCE_VALUE", name, field, sep = ":")
        require_true(identical(is.na(observed), is.na(reference)), label)
        if (is.numeric(reference) && is.numeric(observed)) {
          require_true(identical(is.nan(observed), is.nan(reference)) &&
            identical(is.infinite(observed), is.infinite(reference)), label)
          finite <- is.finite(reference)
          # Discrete counts and design coefficients must never be tolerance-laundered.
          discrete <- is.integer(reference) || name %in% c("mapping_coverage", "design_out", "contrast_out", "filter_out", "gate") ||
            field %in% c("NGenes", "size", "n_omissions", "n_estimable")
          limit <- if (discrete) 0 else tolerance$absolute + tolerance$relative * abs(reference[finite])
          require_true(all(abs(observed[finite] - reference[finite]) <= limit) &&
            identical(sign(observed[is.infinite(reference)]), sign(reference[is.infinite(reference)])), label)
        } else require_true(identical(observed, reference), label)
      }
    }
    list(status = "PASS_NON_GENE_INPUT_DERIVED_REFERENCE_ONLY", verified_tables = names(snapshot),
      output_rows = as.list(vapply(snapshot, nrow, integer(1))), numeric_tolerance = tolerance,
      production_go = FALSE, full_result_validation = FALSE, independent_audit = FALSE,
      linux_native_acceptance = FALSE, input_object_identity = input_identity,
      source_pins = as.list(pins), packages = as.list(versions),
      limitation = "Twelve scientific tables only; use separate gene validator for three gene tables. Shared input helpers and statistical libraries permit common-mode errors. No input-file custody, runtime, production, independent personnel or biological validity certification.")
  }
  list(compare = compare, input_object_identity = input_identity, source_pins = as.list(pins),
    production_go = FALSE)
}
