# Read-only public-resource component for the prospective non-KMS consumer.
# No study counts, legacy launch, signatures, protected state or scientific fit.
# Exact captured gzip identity is primary; decoder checks are not a proof of
# standalone gzip member boundaries or a sandbox for arbitrary serialized R.

gse205906_validate_resource_tables_v1 <- function(resources, exports) {
  fail <- function(x) stop(x, call. = FALSE)
  roles <- c("reactome", "progeny_human", "progeny_mouse", "dorothea_human",
             "dorothea_mouse", "dorothea_rat", "syngo_sensitivity",
             "matrisome_human", "matrisome_mouse", "ortholog_one_to_one")
  if (!is.list(resources) || !identical(names(resources), roles) ||
      !is.list(exports) || !identical(names(exports), roles)) fail("RESOURCE_ROLE_ORDER")
  summaries <- list()
  for (role in roles) {
    a <- resources[[role]]
    b <- exports[[role]]
    if (!is.data.frame(a) || !is.data.frame(b) || !nrow(a) ||
        !identical(names(a), names(b)) || nrow(a) != nrow(b)) fail(paste0("RESOURCE_SCHEMA:", role))
    a <- as.data.frame(a)
    b <- as.data.frame(b)
    maximum_relative_difference <- 0
    for (column in names(a)) {
      x <- a[[column]]; y <- b[[column]]
      if (!identical(is.na(x), is.na(y))) fail(paste0("RESOURCE_MISSINGNESS:", role, ":", column))
      if (is.numeric(x)) {
        if (!is.numeric(y) || any(!is.finite(x[!is.na(x)])) ||
            any(!is.finite(y[!is.na(y)]))) fail(paste0("RESOURCE_NUMERIC_DOMAIN:", role, ":", column))
        keep <- !is.na(x)
        d <- abs(x[keep] - y[keep]); scale <- pmax(abs(x[keep]), abs(y[keep]))
        # TSV writer uses limited decimal precision. Fixed before first check:
        # relative 1e-10, zero absolute tolerance (including tiny p-values).
        if (any(d > 1e-10 * scale)) fail(paste0("RESOURCE_NUMERIC_MISMATCH:", role, ":", column))
        if (any(scale > 0)) maximum_relative_difference <- max(maximum_relative_difference, d[scale > 0] / scale[scale > 0])
      } else if ((!is.character(x) && !is.factor(x)) ||
                 !identical(as.character(x), as.character(y))) {
        fail(paste0("RESOURCE_TEXT_MISMATCH:", role, ":", column))
      }
    }
    if (role %in% c("reactome", "syngo_sensitivity")) {
      keys <- c("term_id", "gene_symbol")
      if (!all(c(keys, "term_name") %in% names(a)) ||
          anyNA(a[, c(keys, "term_name")]) || any(!nzchar(a$gene_symbol)) ||
          any(!nzchar(a$term_id)) || any(!nzchar(a$term_name))) fail(paste0("RESOURCE_UNSIGNED_DOMAIN:", role))
      program <- a$term_id
    } else if (startsWith(role, "progeny")) {
      keys <- c("pathway", "gene_symbol")
      if (!all(c(keys, "weight", "p.value", "rank_within_pathway") %in% names(a)) ||
          anyNA(a[, c(keys, "weight", "p.value", "rank_within_pathway")]) ||
          any(a$p.value < 0 | a$p.value > 1) ||
          length(unique(a$pathway)) != 14L || nrow(a) != 7000L ||
          any(table(a$pathway) != 500L) ||
          !all(vapply(split(a$rank_within_pathway, a$pathway), function(x) identical(sort(as.integer(x)), 1:500) && all(x == floor(x)), logical(1)))) fail(paste0("RESOURCE_PROGENY_DOMAIN:", role))
      program <- a$pathway
    } else if (startsWith(role, "dorothea")) {
      keys <- c("tf", "target")
      if (!all(c(keys, "mor", "confidence") %in% names(a)) ||
          anyNA(a[, c(keys, "mor", "confidence")]) ||
          !all(a$confidence %in% c("A", "B")) || !all(a$mor %in% c(-1, 1))) fail(paste0("RESOURCE_DOROTHEA_DOMAIN:", role))
      program <- a$tf
    } else if (startsWith(role, "matrisome")) {
      keys <- "source_gene_symbol"
      if (!all(c(keys, "category", "gene_symbol") %in% names(a))) fail(paste0("RESOURCE_MATRISOME_SCHEMA:", role))
      program <- a$category
    } else {
      keys <- c("human_symbol", "other_tax_id")
      if (!all(c(keys, "other_symbol", "group_method", "human_tax_id") %in% names(a)) ||
          anyNA(a[, c(keys, "other_symbol", "group_method", "human_tax_id")]) ||
          !all(a$group_method == "NCBI Ortholog") || !all(a$human_tax_id == 9606L) ||
          !all(a$other_tax_id %in% c(10090L, 10116L)) ||
          anyDuplicated(a[, keys]) || anyDuplicated(a[, c("other_symbol", "other_tax_id")])) fail("RESOURCE_ORTHOLOGY_DOMAIN")
      program <- a$other_tax_id
    }
    summaries[[role]] <- list(rows = nrow(a), columns = names(a),
      classes = lapply(a, class), factor_levels = lapply(a[vapply(a, is.factor, logical(1))], function(x) {
        lv <- levels(x)
        list(n_levels = length(lv), n_observed = length(unique(x)),
             n_unused = length(setdiff(lv, as.character(x))),
             levels_sha256 = digest::digest(serialize(lv, NULL, version = 3), "sha256", serialize = FALSE),
             levels_if_at_most_20 = if (length(lv) <= 20L) lv else NULL)
      }),
      missing_by_column = as.list(colSums(is.na(a))), duplicate_rows = sum(duplicated(a)),
      duplicate_keys = sum(duplicated(a[, keys, drop = FALSE])),
      unique_programs_or_species = length(unique(program)),
      numeric_max_relative_export_difference = maximum_relative_difference,
      original_rows_removed = 0L)
  }
  list(status = "PASS_FROZEN_RESOURCE_CONTENT_EXPORT_AGREEMENT_ONLY", tables = summaries,
       production_go = FALSE, independent_audit = FALSE,
       relative_tolerance = 1e-10, absolute_tolerance = 0,
       selected_or_count_filtered_programs = "NOT_COMPUTED")
}

gse205906_resource_label_adapter_v1 <- function(resources) {
  # In-memory representation only: preserve factor labels, never integer codes.
  # Original serialized resources and their unused levels remain unchanged.
  lapply(resources, function(table) {
    if (!is.data.frame(table)) stop("RESOURCE_ADAPTER_TABLE", call. = FALSE)
    result <- as.data.frame(data.table::copy(table))
    for (column in names(result)) {
      value <- result[[column]]
      if (is.factor(value)) {
        result[[column]] <- as.character(value)
      } else if (is.object(value) || !is.atomic(value)) {
        stop("RESOURCE_ADAPTER_NONPLAIN_COLUMN", call. = FALSE)
      }
    }
    result
  })
}

load_gse205906_resource_validation_v1 <- function(project_root) {
  fail <- function(x) stop(x, call. = FALSE)
  for (package in c("digest", "jsonlite", "data.table")) {
    if (!requireNamespace(package, quietly = TRUE)) fail(paste0("RESOURCE_DEPENDENCY:", package))
  }
  source_pins <- c(
    "code/R/gse205906_bulk_model_v2_r4_lib.R" = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3",
    "config/GSE205906_bulk_model_v2_r4.json" = "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c",
    "metadata/gene_set_manifest.tsv" = "b3a9089a6e136c5271cf6ac3999d10304dc2316fafd10f5c170d9f43084bbce1",
    "code/R/freeze_gene_sets.R" = "01b8f54303656fe063476ce272b30dd93aa668d023ecb6a6c3c947a7d908d6ea",
    "config/Route_B2_proportionate_governance_adopted_v1_20260908.json" = "6f38731652c497c11174da665e3531bd14e07479d5aa42c7f85f36743a728b4e"
  )
  prefix <- "data/processed/gene_sets/"
  paths <- paste0(prefix, c("reactome97_human_long.tsv", "progeny_human_top500.tsv",
    "progeny_mouse_top500.tsv", "dorothea_human_AB.tsv", "dorothea_mouse_AB.tsv",
    "dorothea_rat_AB_one_to_one.tsv", "quickgo_syngo_mammalian_direct_human.tsv",
    "matrisome_human.tsv", "matrisome_mouse.tsv", "ncbi_human_mouse_rat_one_to_one.tsv"))
  names(paths) <- c("reactome", "progeny_human", "progeny_mouse", "dorothea_human",
    "dorothea_mouse", "dorothea_rat", "syngo_sensitivity", "matrisome_human", "matrisome_mouse", "ortholog_one_to_one")
  rds_path <- paste0(prefix, "gene_set_resources_frozen_2026-08-18.rds")
  allowed <- c(names(source_pins), unname(paths), rds_path)
  capture <- function(relative, expected_sha, expected_bytes = NULL) {
    if (!relative %in% allowed) fail("RESOURCE_PATH_NOT_ALLOWLISTED")
    path <- file.path(project_root, relative)
    info <- file.info(path)
    if (nrow(info) != 1L || is.na(info$size) || info$isdir || info$size > 20000000 ||
        nzchar(Sys.readlink(path))) fail("RESOURCE_FILE_INVALID")
    con <- file(path, "rb"); on.exit(close(con), add = TRUE)
    raw <- readBin(con, "raw", n = as.integer(info$size) + 1L)
    if (length(raw) != info$size || (!is.null(expected_bytes) && length(raw) != expected_bytes) ||
        !identical(digest::digest(raw, "sha256", serialize = FALSE), expected_sha)) fail(paste0("RESOURCE_FILE_IDENTITY:", relative))
    raw
  }
  captured <- Map(capture, names(source_pins), unname(source_pins))
  names(captured) <- names(source_pins)
  cfg <- jsonlite::fromJSON(rawToChar(captured[[names(source_pins)[2L]]]), simplifyVector = FALSE)
  versions <- vapply(c("digest", "jsonlite", "data.table"), function(p) as.character(utils::packageVersion(p)), character(1))
  if (!identical(as.character(getRversion()), cfg$software$R) ||
      !all(vapply(names(versions), function(p) identical(versions[[p]], cfg$software[[p]]), logical(1)))) fail("RESOURCE_RUNTIME_VERSION")
  shared <- rawToChar(captured[[names(source_pins)[1L]]])
  helpers <- new.env(parent = baseenv())
  spans <- list()
  for (anchors in list(c("`%||%` <- function", "v2r4_child_request_contract <- function"),
                       c("v2r4_as_raw <- function", "v2r4_write_tsv <- function"))) {
    start <- gregexpr(anchors[1L], shared, fixed = TRUE)[[1L]]
    end <- gregexpr(anchors[2L], shared, fixed = TRUE)[[1L]]
    if (length(start) != 1L || length(end) != 1L || start < 1L || end <= start) fail("RESOURCE_DECODER_SPAN")
    text <- substr(shared, start, end - 1L)
    block <- parse(text = text, keep.source = FALSE)
    for (expr in block) {
      if (!is.call(expr) || !identical(expr[[1L]], as.name("<-")) ||
          !is.call(expr[[3L]]) || !identical(expr[[3L]][[1L]], as.name("function"))) fail("RESOURCE_HELPER_DECLARATION")
      eval(expr, helpers)
    }
    spans[[anchors[1L]]] <- digest::digest(charToRaw(text), "sha256", serialize = FALSE)
  }
  lockEnvironment(helpers, bindings = TRUE)
  manifest <- data.table::fread(text = rawToChar(captured[["metadata/gene_set_manifest.tsv"]]), sep = "\t", colClasses = "character", showProgress = FALSE)
  if (anyDuplicated(manifest$relative_path)) fail("RESOURCE_MANIFEST_DUPLICATE")
  provenance <- list()
  for (path in c(unname(paths), rds_path)) {
    row <- manifest[which(manifest$relative_path == path)]
    if (nrow(row) != 1L || row$analysis_eligible != "TRUE") fail("RESOURCE_MANIFEST_ROLE")
    captured[[path]] <- capture(path, row$sha256, as.numeric(row$size_bytes))
    provenance[[path]] <- list(bytes = length(captured[[path]]), sha256 = row$sha256)
  }
  if (!identical(provenance[[rds_path]]$sha256, "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51")) fail("RESOURCE_RDS_IDENTITY")
  # Every input is captured and hashed before this first object deserialization.
  resources <- helpers$v2r4_read_rds_raw(captured[[rds_path]], "FROZEN_PUBLIC_GENE_SETS", cfg$input_snapshot$rds_role_contracts$gene_set_object$compression)
  exports <- lapply(paths, function(path) data.table::fread(text = rawToChar(captured[[path]]), sep = "\t", na.strings = c("NA", ""), showProgress = FALSE))
  # Exact frozen producer projection, not replacement with the full annotation table.
  exports$syngo_sensitivity <- unique(exports$syngo_sensitivity[, c("term_id", "term_name", "gene_symbol"), with = FALSE])
  validator <- gse205906_validate_resource_tables_v1
  result <- validator(resources, exports)
  consumer_resources <- gse205906_resource_label_adapter_v1(resources)
  adapter_validation <- validator(resources, consumer_resources)
  result$adapter <- list(rule = "FACTOR_LABELS_TO_CHARACTER_IN_SEPARATE_PLAIN_DATA_FRAME_COPY",
    column_count = sum(vapply(resources, function(x) sum(vapply(x, is.factor, logical(1))), integer(1))),
    equivalence_status = adapter_validation$status,
    source_object_sha256 = digest::digest(serialize(resources, NULL, version = 3), "sha256", serialize = FALSE),
    consumer_object_sha256 = digest::digest(serialize(consumer_resources, NULL, version = 3), "sha256", serialize = FALSE),
    raw_or_original_object_modified = FALSE)
  result$source_pins <- source_pins
  result$input_provenance <- provenance
  result$decoder_span_sha256 <- spans
  result$packages <- as.list(versions)
  result$R <- as.character(getRversion())
  result$counts_read <- FALSE
  result$legacy_launcher_executed <- FALSE
  result$limitations <- c("Public resource content only; no species/count input release or model fit.",
    "Exports and RDS share a producer; agreement is not independent reconstruction from upstream databases.",
    "Exact payload identity protects this known object; deserialization is not safe for arbitrary untrusted RDS.",
    "No OS lease or hostile-administrator isolation; not an accepted production consumer.",
    "QuickGO assignedBy SynGO remains sensitivity-only and is not the SynGO portal bulk release.")
  list(resources = resources, exports = exports, consumer_resources = consumer_resources, validation = result)
}
