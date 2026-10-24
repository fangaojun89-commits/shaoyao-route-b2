# Functions-only matched STAR consumer. No top-level IO/fit, dispatch or authority.
# The reviewed caller must authenticate held files, actual process/resource/QC
# receipts and a fresh cloud output leaf. This module never reads an old y,
# weight, score or label-null object. It returns complete tables/native objects;
# a thin caller owns immutable writing and independent acceptance.
g205m_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
g205m_hash <- function(x) digest::digest(x, "sha256", serializeVersion = 3L)
g205m_bytes_pin <- function(x) {
  g205m_need(is.raw(x), "Held raw bytes required")
  list(bytes = length(x), sha256 = digest::digest(x, "sha256", serialize = FALSE))
}
g205m_pin_ok <- function(x) {
  is.list(x) && is.character(x$path) && length(x$path) == 1L && !is.na(x$path) &&
    nzchar(x$path) && is.numeric(x$bytes) && !is.logical(x$bytes) && length(x$bytes) == 1L &&
    is.finite(x$bytes) && x$bytes >= 0 && x$bytes == floor(x$bytes) &&
    is.character(x$sha256) && length(x$sha256) == 1L && !is.na(x$sha256) &&
    grepl("^[0-9a-f]{64}$", x$sha256)
}
g205m_pins <- function() {
  list(
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    helpers = list(bytes = 175131L, sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
    config = list(bytes = 10612L, sha256 = "ee6e22592b06ecdb31707234e97dc55b7b66cbad8dd11c364187af0624ece980"),
    loader = list(bytes = 10884L, sha256 = "45f21c67674b5cd2b424e9c145662910361baee1019cbc68d366a9aab1fcc4ba"),
    prefit = list(bytes = 18762L, sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
    annotation = list(bytes = 9943L, sha256 = "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470"))
}
g205m_import <- function(bytes, wanted, environment) {
  expressions <- as.list(parse(text = rawToChar(bytes), keep.source = FALSE))
  selected <- Filter(function(x) is.call(x) && length(x) == 3L &&
    identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
    as.character(x[[2L]]) %in% wanted, expressions)
  found <- vapply(selected, function(x) as.character(x[[2L]]), character(1))
  g205m_need(setequal(found, wanted) && !anyDuplicated(found), "Exact function-only source declarations required")
  for (x in selected) {
    g205m_need(is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), "No top-level source execution")
    eval(x, environment)
  }
  invisible(found)
}
g205m_engine <- function(sources) {
  pins <- g205m_pins()
  g205m_need(is.list(sources) && identical(names(sources), names(pins)), "Exact six source roles required")
  for (role in names(pins)) {
    got <- g205m_bytes_pin(sources[[role]])
    g205m_need(got$bytes == pins[[role]]$bytes && identical(got$sha256, pins[[role]]$sha256), "Frozen source/config changed")
  }
  e <- new.env(parent = baseenv())
  # Only the actual non-base, unqualified dependencies in the imported bodies.
  # Do not rely on the caller having attached stats/utils or old driver state.
  e$p.adjust <- stats::p.adjust
  e$head <- utils::head
  e$tail <- utils::tail
  g205m_import(sources$helpers, c("%||%", "v2r4_stop", "v2r4_assert", "v2r4_design",
    "v2r4_validate_counts_original", "v2r4_validate_annotation_gene", "v2r4_assert_deseq_contract",
    "v2r4_collapse_condition_messages", "v2r4_assert_deseq_fit_provenance", "v2r4_add_claim_flags"), e)
  g205m_import(sources$model, c("v2r4_fit_preprocess", "v2r4_fit_voom", "v2r4_fit_deseq",
    "v2r4_map_unsigned", "v2r4_unsigned_tests", "v2r4_collapse_symbols", "v2r4_activity_domain", "v2r4_fit_activity"), e)
  g205m_import(sources$loader, c("gse205906_ens101_loader_need_v2", "gse205906_ens101_loader_hash_v2",
    "gse205906_ens101_gtf_attributes_v2", "gse205906_ens101_input_loader_v2"), e)
  g205m_import(sources$annotation, c("gse205906_ens101_need_v1", "gse205906_ens101_character_v1",
    "gse205906_ens101_annotation_v1"), e)
  e$cfg <- jsonlite::fromJSON(rawToChar(sources$config), simplifyVector = FALSE)
  e$cli <- list(execution_mode = "PRODUCTION")
  e$seed <- as.integer(e$cfg$analysis$random_seed)
  e$sources <- sources
  e$function_names <- sort(ls(e)[vapply(ls(e), function(k) is.function(e[[k]]), logical(1))])
  e
}
g205m_engine_hash <- function(e) {
  functions <- lapply(e$function_names, function(k) list(formals = formals(e[[k]]), body = body(e[[k]])))
  names(functions) <- e$function_names
  g205m_hash(list(sources = e$sources, cfg = e$cfg, cli = e$cli, seed = e$seed,
    functions = functions, resources = e$gene_sets, rat_map = e$rat_map))
}

# Validate dynamic STAR matrix SHA from the ACTUAL held producer COMPLETE,
# never from an invented future checksum. Filesystem custody remains caller-owned.
g205m_matrix_binding <- function(matrix_bytes, binding) {
  required <- c("schema", "matrix_path", "matrix_complete_bytes", "matrix_complete_file_pin",
    "producer_source_sha256", "held_files", "held_object_hashes", "hisat2_comparator_root")
  g205m_need(is.list(binding) && identical(names(binding), required) &&
    identical(binding$schema, "gse205906-matched-downstream-held-inputs-v1"), "Explicit held-input binding schema required")
  g205m_need(g205m_pin_ok(binding$matrix_complete_file_pin) &&
    is.character(binding$matrix_path) && length(binding$matrix_path) == 1L && !is.na(binding$matrix_path) &&
    identical(binding$matrix_path, "/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s0.tsv") &&
    identical(binding$producer_source_sha256, "7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d"),
    "Current reviewed STAR producer and explicit primary matrix path required")
  cp <- g205m_bytes_pin(binding$matrix_complete_bytes)
  g205m_need(cp$bytes == binding$matrix_complete_file_pin$bytes &&
    identical(cp$sha256, binding$matrix_complete_file_pin$sha256) &&
    identical(binding$matrix_complete_file_pin$path, paste0(dirname(binding$matrix_path), "/COMPLETE.json")),
    "Held matrix COMPLETE bytes/path differ from caller pin")
  complete <- jsonlite::fromJSON(rawToChar(binding$matrix_complete_bytes), simplifyVector = FALSE)
  g205m_need(identical(complete$status, "COMPLETE") &&
    identical(complete$source_sha256, binding$producer_source_sha256) &&
    identical(complete$metadata_sha256, "6dd7b53d40433baaacdeedf285204dee83f2d960c090e6a25f489f96d20fc5cb") &&
    identical(as.integer(unlist(complete$dimensions)), c(32883L, 9L)) && complete$primary_strand == 0 &&
    identical(complete$new_models_or_candidate_admission, FALSE) && length(complete$full_libraries) == 9L &&
    identical(complete$original_matrix_preserved$sha256, "e98296fdf308dbaa82325a3ca6b0b041c6f2569e98be75a45757141b0a1f77a2"),
    "Actual complete full-nine primary STAR matrix required")
  matches <- Filter(function(x) identical(x$path, binding$matrix_path), complete$small_products)
  g205m_need(length(matches) == 1L && g205m_pin_ok(matches[[1L]]), "Exactly one completed primary matrix record required")
  got <- g205m_bytes_pin(matrix_bytes); declared <- matches[[1L]]
  g205m_need(got$bytes == declared$bytes && identical(got$sha256, declared$sha256), "New matrix bytes differ from actual producer COMPLETE")
  roles <- c("metadata", "fixed_gene_manifest", "GTF", "gene_resources", "orthology",
    "reference_admission", "statistical_runtime", "hisat2_completion", "hisat2_manifest", "hisat2_matrix")
  g205m_need(is.list(binding$held_files) && identical(names(binding$held_files), roles) &&
    all(vapply(binding$held_files, g205m_pin_ok, logical(1))) &&
    is.character(binding$hisat2_comparator_root) && length(binding$hisat2_comparator_root) == 1L &&
    identical(binding$hisat2_comparator_root, "/root/autodl-tmp/GSE205906_Ensembl101_checked_model_v5_20261006"),
    "Explicit current resource/runtime and preserved comparator records required")
  fixed <- c(metadata = "b0df6fa1354a83237bcfcd8bf8c323517deecd0e5680da8202489c4464936444",
    fixed_gene_manifest = "751e43b3aba8c2245981dc0041fec273c455c6390e85ef7ec4514e56118bf3b0",
    GTF = "84fcb6692d62e4d8c75617461ae5f471fe89e10b6982c82f349393b1a80045df",
    gene_resources = "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51",
    orthology = "e9dd0c9268e0044f615085e27292cbf977f5a648003ba90054cf3d01225cf998",
    hisat2_completion = "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e",
    hisat2_manifest = "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7",
    hisat2_matrix = "e98296fdf308dbaa82325a3ca6b0b041c6f2569e98be75a45757141b0a1f77a2")
  for (role in names(fixed)) g205m_need(identical(binding$held_files[[role]]$sha256, unname(fixed[[role]])), "Frozen input/comparator record changed")
  list(matrix = declared, completion = binding$matrix_complete_file_pin,
    full_libraries = complete$full_libraries, input_file_custody_authenticated_here = FALSE)
}
g205m_prepare <- function(sources, matrix_bytes, sample_ledger, fixed_gene_ids,
                         gene_attributes, gene_resources, orthology_table, binding) {
  e <- g205m_engine(sources)
  matrix_binding <- g205m_matrix_binding(matrix_bytes, binding)
  objects <- list(sample_ledger = sample_ledger, fixed_gene_ids = fixed_gene_ids,
    gene_attributes = gene_attributes, gene_resources = gene_resources, orthology_table = orthology_table)
  g205m_need(identical(binding$held_object_hashes, lapply(objects, g205m_hash)), "Caller-bound held objects changed")
  converted <- e$gse205906_ens101_input_loader_v2(matrix_bytes, sample_ledger, fixed_gene_ids,
    sources$prefit, gene_attributes = gene_attributes)
  projected <- e$gse205906_ens101_annotation_v1(converted$inputs$gene_annotation, fixed_gene_ids)
  g205m_need(is.list(gene_resources) && all(c("reactome", "progeny_human", "dorothea_rat", "syngo_sensitivity") %in% names(gene_resources)) &&
    is.data.frame(orthology_table), "Bound frozen resource/orthology objects required")
  orthology <- data.table::as.data.table(orthology_table)
  g205m_need(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in% names(orthology)), "Strict external orthology columns missing")
  e$rat_map <- unique(orthology[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
    .(human_symbol = as.character(human_symbol), rat_symbol = as.character(other_symbol))])
  g205m_need(nrow(e$rat_map) > 0L && !anyNA(e$rat_map) && !anyDuplicated(e$rat_map$human_symbol) &&
    !anyDuplicated(e$rat_map$rat_symbol), "Strict rat orthology not one-to-one")
  e$gene_sets <- gene_resources
  design <- e$v2r4_design(converted$inputs$sample_ledger, e$cfg)
  counts <- e$v2r4_validate_counts_original(converted$inputs$count_matrix, rownames(design$design))
  annotation <- e$v2r4_validate_annotation_gene(projected$computational_interface, rownames(counts), .8)
  packet <- list(schema = "gse205906-matched-downstream-prepared-v1", counts = counts,
    annotation = annotation, design = design, binding = binding, matrix_binding = matrix_binding,
    source_pins = g205m_pins(), input_conversion = converted$transformation,
    engine_sha256 = g205m_engine_hash(e),
    gate = e$cfg$gate, sample_unit = "NINE_SUBMITTED_LIBRARIES_NOT_AUTHENTICATED_ANIMALS",
    input_file_custody_authenticated_here = FALSE, old_models_scores_or_nulls_reused = FALSE)
  packet$sha256 <- g205m_hash(packet)
  list(packet = packet, engine = e)
}

# Return-only capture: the original pinned body runs exactly once, with the
# original formals/environment. Only its answer plus selected local objects are
# retained. No model/statistic/RNG expression is replaced or called twice.
g205m_capture <- function(fun, locals, remove_condition_text = FALSE) {
  g205m_need(is.function(fun) && is.character(locals) && !anyNA(locals) && !anyDuplicated(locals), "Explicit original function/local names required")
  collect <- as.call(c(list(as.name("list")), stats::setNames(lapply(locals, as.name), locals)))
  body(fun) <- substitute({
    g205m_answer <- ORIGINAL
    g205m_conditions <- NULL
    if (REMOVE_CONDITIONS) {
      g205m_conditions <- list(raw_condition_text_saved = FALSE,
        collapsed_condition_sha256 = digest::digest(g205m_answer$fit_messages, "sha256", serializeVersion = 3L))
      g205m_answer$fit_messages <- NULL
      # dds@design can retain this call frame. Remove raw condition locals here,
      # not just from the outer result list, before returning the saved object.
      rm(list = intersect(c("fit_messages_raw", "fit_messages"), ls()))
    }
    list(value = g205m_answer, native = LOCALS, condition_fingerprint = g205m_conditions)
  }, list(ORIGINAL = body(fun), LOCALS = collect, REMOVE_CONDITIONS = remove_condition_text))
  fun
}
g205m_empty <- function(kind) {
  common <- data.table::data.table(collection = character(), program_id = character(), contrast = character(), method = character())
  if (kind == "activity") return(cbind(common, data.table::data.table(effect = double(), se = double(), ci_low = double(),
    ci_high = double(), t = double(), df_total = double(), p_value = double(), q_within_collection_contrast = double())))
  if (kind == "camera") return(cbind(common, data.table::data.table(NGenes = integer(), Direction = character(),
    PValue = double(), FDR = double(), p_value = double(), q_within_collection_contrast = double())))
  if (kind == "mroast") return(cbind(common, data.table::data.table(NGenes = integer(), PropDown = double(), PropUp = double(),
    Direction = character(), PValue = double(), FDR = double(), PValue.Mixed = double(), FDR.Mixed = double(),
    p_value = double(), q_within_collection_contrast = double())))
  g205m_need(kind == "fgsea", "Known typed-empty family required")
  cbind(common, data.table::data.table(pval = double(), padj = double(), log2err = double(), ES = double(), NES = double(),
    size = integer(), leadingEdge = I(list()), p_value = double(), q_within_collection_contrast = double()))
}
g205m_unsigned_schema <- function(value) {
  g205m_need(is.list(value) && identical(names(value), c("camera", "mroast", "fgsea")) &&
    all(vapply(value, is.data.frame, logical(1))), "Native unsigned return schema differs")
  # The original fgsea body only adds aliases when nrow(fg)>0. A legitimate
  # zero-row native return must remain zero rows, but needs the declared columns.
  if (nrow(value$fgsea) == 0L) value$fgsea <- g205m_empty("fgsea")
  value
}
g205m_slots <- function(tables) {
  specs <- rbind(data.frame(table = c("gene_edger", "gene_voom", "gene_deseq"), collection = "genes"),
    expand.grid(table = c("unsigned_camera", "unsigned_mroast", "unsigned_fgsea"),
      collection = c("Reactome", "SynGO_sensitivity"), stringsAsFactors = FALSE),
    data.frame(table = "activity_results", collection = c("PROGENy", "DoRothEA_AB")))
  out <- lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, ]; tab <- tables[[spec$table]]
    p_column <- if (spec$table == "gene_edger") "PValue" else "p_value"
    q_column <- if (spec$collection == "genes") "q_value" else "q_within_collection_contrast"
    g205m_need(is.data.frame(tab) && !anyDuplicated(names(tab)) &&
      all(c("contrast", p_column, q_column) %in% names(tab)) &&
      !anyNA(tab$contrast) && all(tab$contrast %in% c("local_disease", "intervention")), "Native table fields/contrast domain differ")
    for (column in c(p_column, q_column)) {
      value <- tab[[column]]
      g205m_need(is.numeric(value) && all(is.na(value) | (is.finite(value) & value >= 0 & value <= 1)),
        "Invalid native P/q; no imputation or silent row deletion")
    }
    do.call(rbind, lapply(c("local_disease", "intervention"), function(contrast_name) {
      keep <- tab$contrast == contrast_name
      if (spec$collection != "genes") keep <- keep & tab$collection == spec$collection
      p <- if (spec$table == "gene_edger") tab$PValue[keep] else tab$p_value[keep]
      q <- if (spec$collection == "genes") tab$q_value[keep] else tab$q_within_collection_contrast[keep]
      data.frame(table = spec$table, collection = spec$collection, contrast = contrast_name, rows = sum(keep),
        native_P_missing = sum(is.na(p)), native_q_missing = sum(is.na(q)),
        q_strict_lt_05 = sum(q < .05, na.rm = TRUE), q_exactly_05 = sum(q == .05, na.rm = TRUE),
        status = if (sum(keep)) "OWN_COMPLETE_ELIGIBLE_FAMILY" else "EMPTY_OWN_ELIGIBLE_FAMILY_NOT_A_PASS",
        stringsAsFactors = FALSE)
    }))
  })
  ans <- do.call(rbind, out)
  g205m_need(nrow(ans) == 22L && !anyDuplicated(ans[c("table", "collection", "contrast")]), "Full22 slots required")
  rownames(ans) <- NULL; ans
}
g205m_domains <- function(tables, slots, domains) {
  absent <- list()
  for (i in seq_len(nrow(slots))) {
    row <- slots[i, ]; tab <- tables[[row$table]]
    selected <- tab$contrast == row$contrast
    if (row$collection != "genes") selected <- selected & tab$collection == row$collection
    actual <- if (row$collection == "genes") tab$gene_id[selected] else tab$program_id[selected]
    expected <- domains[[row$table]][[row$collection]]
    g205m_need(is.character(actual) && !anyNA(actual) && !anyDuplicated(actual) &&
      all(actual %in% expected), "Unexpected/duplicate native domain keys")
    if (row$collection == "genes") g205m_need(identical(actual, expected), "Complete ordered native gene family required")
    else if (row$table != "unsigned_fgsea") g205m_need(setequal(actual, expected), "Complete native eligible program family required")
    missing <- setdiff(expected, actual)
    slots$eligible_input_program_or_gene_n[i] <- length(expected)
    slots$native_not_returned_n[i] <- length(missing)
    slots$expected_keys_sha256[i] <- g205m_hash(expected)
    slots$returned_keys_sha256[i] <- g205m_hash(actual)
    if (length(missing)) absent[[length(absent) + 1L]] <- data.frame(table = row$table,
      collection = row$collection, contrast = row$contrast, program_id = missing,
      reason = "NOT_RETURNED_BY_UNCHANGED_NATIVE_FGSEA_AFTER_FINITE_T_RANK_AND_SIZE_RULES", stringsAsFactors = FALSE)
    if (length(missing)) slots$status[i] <- "NATIVE_METHOD_DOMAIN_RESTRICTION_NOT_ELIGIBILITY_PASS"
  }
  list(slots = slots, native_not_returned = if (length(absent)) do.call(rbind, absent) else
    data.frame(table = character(), collection = character(), contrast = character(), program_id = character(), reason = character()))
}
g205m_runtime <- function(e) {
  g205m_need(identical(as.character(getRversion()), e$cfg$software$R), "Pinned R version required")
  packages <- c("edgeR", "limma", "DESeq2", "fgsea", "decoupleR", "data.table", "jsonlite", "digest")
  for (package in packages) g205m_need(requireNamespace(package, quietly = TRUE) &&
    identical(as.character(utils::packageVersion(package)), e$cfg$software[[package]]), "Pinned package version required")
  g205m_need(requireNamespace("BiocParallel", quietly = TRUE) &&
    identical(as.character(utils::packageVersion("BiocParallel")), "1.46.0"), "Pinned SerialParam backend required")
  g205m_need(all(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "R_DATATABLE_NUM_THREADS")) == "1"),
    "Caller must establish explicit single-thread runtime")
  BiocParallel::register(BiocParallel::SerialParam(workers = 1L, progressbar = FALSE), default = TRUE)
  do.call(RNGkind, list("Mersenne-Twister", "Inversion", "Rejection")); set.seed(e$seed)
  list(R = as.character(getRversion()), packages = stats::setNames(lapply(packages, function(x) as.character(utils::packageVersion(x))), packages),
    backend = "SerialParam", workers = 1L, seed = e$seed, RNGkind = RNGkind(),
    mroast_rotations = 9999L, fgsea_eps = 0, historical_runtime_authenticated = FALSE)
}
g205m_fit <- function(prepared) {
  p <- prepared$packet; e <- prepared$engine
  declared <- p$sha256; p$sha256 <- NULL
  g205m_need(identical(g205m_hash(p), declared) && identical(p$source_pins, g205m_pins()) &&
    identical(p$engine_sha256, g205m_engine_hash(e)), "Prepared input/source/resource identity changed")
  runtime <- g205m_runtime(e); started <- proc.time()[["elapsed"]]
  base <- e$v2r4_fit_preprocess(p$counts, p$design$design, e$cfg$analysis$minimum_retained_genes_production)
  ann <- p$annotation$annotation[match(rownames(base$y$counts), gene_id)]
  e$retained_annotation <- ann
  contrasts <- p$design$contrast; design <- p$design$design
  tests <- stats::setNames(lapply(seq_len(ncol(contrasts)), function(j) edgeR::glmQLFTest(base$ql, contrast = contrasts[, j])), colnames(contrasts))
  edge <- data.table::rbindlist(lapply(names(tests), function(name) {
    test <- tests[[name]]
    tab <- data.table::as.data.table(test$table, keep.rownames = "gene_id")
    g205m_need(identical(tab$gene_id, rownames(base$y$counts)), "Native QL gene order differs")
    tab[, `:=`(symbol = ann$symbol, contrast = name, method = "edgeR_QL_robust", FDR = stats::p.adjust(PValue, "BH"))]
    tab[, `:=`(q_value = FDR, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, se_df = NA_real_,
      se_status = "NOT_REPORTED_NO_VALIDATED_NATIVE_QL_INTERVAL", df_total = if (is.null(test$df.total)) rep(NA_real_, .N) else test$df.total)]
    tab
  }), fill = TRUE)
  voom <- g205m_capture(e$v2r4_fit_voom, "fit")(base$voom, design, contrasts, ann)
  deseq <- g205m_capture(e$v2r4_fit_deseq, "dds", remove_condition_text = TRUE)(base$y$counts, p$design$group, contrasts, ann)
  e$unsigned_resources <- list(Reactome = e$v2r4_map_unsigned(e$gene_sets$reactome, "Reactome"),
    SynGO_sensitivity = e$v2r4_map_unsigned(e$gene_sets$syngo_sensitivity, "SynGO_sensitivity"))
  unsigned <- lapply(e$unsigned_resources, function(resource) {
    if (!length(resource$indices)) return(list(camera = g205m_empty("camera"), mroast = g205m_empty("mroast"), fgsea = g205m_empty("fgsea")))
    g205m_unsigned_schema(e$v2r4_unsigned_tests(resource, base$voom, design, contrasts))
  })
  activity <- e$v2r4_activity_domain(base$voom, ann, enforce_expected = FALSE)
  activity_fits <- lapply(names(activity$scores), function(collection) {
    score <- activity$scores[[collection]]
    if (!nrow(score)) return(list(value = g205m_empty("activity"), native = list(fit = NULL)))
    g205m_capture(e$v2r4_fit_activity, "fit")(score, collection, design, contrasts)
  }); names(activity_fits) <- names(activity$scores)
  tables <- list(gene_edger = edge, gene_voom = voom$value, gene_deseq = deseq$value$table,
    unsigned_camera = data.table::rbindlist(lapply(unsigned, `[[`, "camera"), fill = TRUE),
    unsigned_mroast = data.table::rbindlist(lapply(unsigned, `[[`, "mroast"), fill = TRUE),
    unsigned_fgsea = data.table::rbindlist(lapply(unsigned, `[[`, "fgsea"), fill = TRUE),
    activity_results = data.table::rbindlist(lapply(activity_fits, `[[`, "value"), fill = TRUE))
  tables <- lapply(tables, e$v2r4_add_claim_flags, analysis_role = e$cfg$gate$analysis_role)
  coverage <- data.table::rbindlist(c(lapply(e$unsigned_resources, `[[`, "coverage"), list(activity$coverage)), fill = TRUE)
  unsigned_domains <- lapply(e$unsigned_resources, function(x) as.character(names(x$indices)))
  gene_domain <- list(genes = rownames(base$y$counts))
  domains <- list(gene_edger = gene_domain, gene_voom = gene_domain, gene_deseq = gene_domain,
    unsigned_camera = unsigned_domains, unsigned_mroast = unsigned_domains, unsigned_fgsea = unsigned_domains,
    activity_results = lapply(activity$scores, function(x) as.character(rownames(x))))
  checked_domains <- g205m_domains(tables, g205m_slots(tables), domains)
  dds <- deseq$native$dds; flags <- S4Vectors::mcols(dds)$betaConv
  diag <- list(actual_DESeq_fit_type = deseq$value$actual_fit_type,
    automatic_local_substitution = identical(deseq$value$actual_fit_type, "local"),
    condition_fingerprint = deseq$condition_fingerprint,
    betaConv_available = !is.null(flags), betaConv_false = if (is.null(flags)) NA_integer_ else sum(!flags, na.rm = TRUE),
    betaConv_missing = if (is.null(flags)) NA_integer_ else sum(is.na(flags)),
    source_input_genes = nrow(p$counts), own_retained_genes = sum(base$keep),
    own_filterByExpr_retained = sum(base$keep_filter_by_expr), own_CPM_retained = sum(base$keep_cpm),
    own_library_sizes = base$y$samples$lib.size, own_TMM = base$y$samples$norm.factors,
    original_label_only_no_null_or_LOAO = TRUE, native_mroast_FDR_not_raw_P_BH = TRUE)
  # Raw conditions are intentionally excluded from the return; caller must not
  # serialize a raw warning/error ledger. Scientific native models are retained.
  after <- prepared$packet; after$sha256 <- NULL
  g205m_need(identical(g205m_hash(after), declared) && identical(p$engine_sha256, g205m_engine_hash(e)),
    "Held input/source/resource changed during native calls")
  list(schema = "gse205906-matched-downstream-result-v1", tables = tables, slots = checked_domains$slots,
    native_not_returned_programs = checked_domains$native_not_returned,
    mapping_coverage = coverage, own_activity_scores = activity$scores,
    native = list(base_fit = base, QL_tests = tests, voom_fit = voom$native$fit, DESeq_fit = dds,
      activity_fits = lapply(activity_fits, `[[`, "native")), diagnostics = diag,
    runtime = c(runtime, list(RNG_after_sha256 = g205m_hash(.Random.seed))),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    binding = p$binding, matrix_binding = p$matrix_binding, prepared_sha256 = declared,
    statistical_source_pins = p$source_pins, gate = p$gate,
    candidate_admission = FALSE, independent_numeric_acceptance = FALSE,
    complete_file_custody_or_producer_QC_authenticated_here = FALSE,
    comparator_tables_read_here = FALSE, files_written_here = FALSE)
}
