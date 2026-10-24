# Pure RAM entry: source/input consistency, lossless conversion and ten no-fit
# phases only. No file/launcher/network IO, model authority or automatic export.
gse205906_ens101_nofit_need_v2 <- function(ok) {
  if (!isTRUE(ok)) stop("ENS101_NOFIT_RAM_INPUT_OR_SOURCE_UNVERIFIED", call. = FALSE)
}

gse205906_ens101_nofit_pin_v2 <- function(raw) {
  gse205906_ens101_nofit_need_v2(is.raw(raw))
  list(bytes = length(raw), sha256 = digest::digest(raw, "sha256", serialize = FALSE))
}

gse205906_ens101_nofit_sources_v2 <- function(bundle) {
  need <- gse205906_ens101_nofit_need_v2
  pins <- list(
    loader = list(bytes = 10884L, sha256 = "45f21c67674b5cd2b424e9c145662910361baee1019cbc68d366a9aab1fcc4ba"),
    prefit = list(bytes = 18762L, sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
    config_adapter = list(bytes = 7852L, sha256 = "29ffe9705e675058c006d808305349275979efb58eeab2ea4236e5eb82fd339a"),
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    helpers = list(bytes = 175131L, sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
    annotation_adapter = list(bytes = 9943L, sha256 = "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470"))
  need(is.list(bundle) && identical(names(bundle), names(pins)))
  for (k in names(pins)) need(identical(gse205906_ens101_nofit_pin_v2(bundle[[k]]), pins[[k]]))
  import <- function(raw) {
    e <- new.env(parent = baseenv())
    e$head <- utils::head; e$tail <- utils::tail; e$setNames <- stats::setNames
    declarations <- as.list(parse(text = rawToChar(raw), keep.source = FALSE))
    for (x in declarations) {
      need(is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
             is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")))
      eval(x, e)
    }
    lockEnvironment(e, bindings = TRUE); e
  }
  list(loader = import(bundle$loader), prefit = import(bundle$prefit),
       config = import(bundle$config_adapter), pins = pins)
}

gse205906_ens101_nofit_overlay_v2 <- function(raw) {
  need <- gse205906_ens101_nofit_need_v2
  need(is.raw(raw) && length(raw) > 0L && length(raw) <= 1024^2 && !any(raw == as.raw(0)))
  text <- rawToChar(raw)
  need(!is.na(iconv(text, "UTF-8", "UTF-8", sub = NA)))
  text <- sub("^\ufeff", "", text)
  # textConnection otherwise materializes the final record terminator as an
  # extra empty physical line. Remove only one terminal LF/CRLF, not values.
  text <- sub("\r?\n$", "", text)
  con <- textConnection(text, open = "r", local = TRUE)
  on.exit(close(con), add = TRUE)
  # Match quoted TSV records, not physical lines. Character NA/blank/space and
  # embedded tab/newline/doubled quote values must remain original strings.
  result <- tryCatch(withCallingHandlers(
    utils::read.delim(con, header = TRUE, sep = "\t", quote = '"',
      colClasses = "character", na.strings = character(0), fill = FALSE,
      strip.white = FALSE, comment.char = "", check.names = FALSE,
      stringsAsFactors = FALSE, row.names = NULL, blank.lines.skip = FALSE,
      allowEscapes = FALSE, encoding = "UTF-8"),
    warning = function(w) stop("ENS101_NOFIT_RAM_INPUT_OR_SOURCE_UNVERIFIED", call. = FALSE)),
    error = function(e) NULL)
  need(is.data.frame(result) && nrow(result) == 9L && ncol(result) == 53L &&
    !anyDuplicated(names(result)) && all(vapply(result,
      function(x) is.character(x) && is.null(dim(x)) && length(x) == 9L && !anyNA(x),
      logical(1))))
  rownames(result) <- NULL
  result
}

gse205906_ens101_nofit_metadata_consistency_v2 <- function(manifest, ledger) {
  need <- gse205906_ens101_nofit_need_v2
  fields <- names(ledger)
  need(is.data.frame(ledger) && nrow(ledger) == 9L && length(fields) == 53L &&
    !anyDuplicated(fields))
  rows <- manifest$metadata_rows
  need(is.list(rows) && is.null(dim(rows)) && length(rows) == 9L)
  for (i in seq_len(9L)) {
    row <- rows[[i]]
    need(is.list(row) && !is.null(names(row)) && !anyDuplicated(names(row)) &&
      length(row) == 53L && setequal(names(row), fields))
    values <- vapply(row[fields], function(x) {
      need(is.character(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x))
      x
    }, character(1))
    need(identical(unname(values), unlist(ledger[i, fields, drop = FALSE], use.names = FALSE)))
  }
  invisible(TRUE)
}

gse205906_ens101_nofit_RAM_worker_v2 <- function(source_bundle, full_config_raw,
    baseline_config_raw, science_raw, input_pins) {
  need <- gse205906_ens101_nofit_need_v2
  need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE) &&
         requireNamespace("edgeR", quietly = TRUE))
  source <- gse205906_ens101_nofit_sources_v2(source_bundle)
  pin <- gse205906_ens101_nofit_pin_v2
  need(identical(pin(full_config_raw), list(bytes = 10612L,
    sha256 = "ee6e22592b06ecdb31707234e97dc55b7b66cbad8dd11c364187af0624ece980")))
  full <- source$config$gse205906_ens101_check_full_model_config_v1(full_config_raw, baseline_config_raw)
  cfg <- source$config$gse205906_ens101_prefit_view_v1(full, baseline_config_raw)
  roles <- c("new_counts", "new_annotation", "new_metadata", "fixed_gene_manifest", "new_config",
             "new_prefit_SOURCE", "statistical_runtime_lock")
  need(is.list(science_raw) && identical(names(science_raw), c("matrix", "GTF", "overlay", "manifest")))
  caps <- c(32 * 1024^2, 512 * 1024^2, 1024^2, 2 * 1024^2)
  need(all(vapply(seq_along(science_raw), function(i) is.raw(science_raw[[i]]) &&
    length(science_raw[[i]]) > 0L && length(science_raw[[i]]) <= caps[[i]], logical(1))))
  need(is.list(input_pins) && identical(names(input_pins), roles))
  valid_pin <- function(p) is.list(p) && identical(names(p), c("bytes", "sha256")) &&
    is.numeric(p$bytes) && !is.logical(p$bytes) && length(p$bytes) == 1L && is.finite(p$bytes) &&
    p$bytes > 0 && p$bytes <= 2^53 && p$bytes == floor(p$bytes) &&
    is.character(p$sha256) && length(p$sha256) == 1L && !is.na(p$sha256) && grepl("^[0-9a-f]{64}$", p$sha256)
  need(all(vapply(input_pins, valid_pin, logical(1))))
  observed <- c(lapply(science_raw, pin), list(pin(full_config_raw), pin(source_bundle$prefit)))
  for (i in seq_along(observed)) need(identical(as.numeric(input_pins[[i]]$bytes), as.numeric(observed[[i]]$bytes)) &&
    identical(input_pins[[i]]$sha256, observed[[i]]$sha256))
  # The seventh runtime pin is only a declaration: no runtime bytes are supplied.
  ledger <- gse205906_ens101_nofit_overlay_v2(science_raw$overlay)
  manifest <- jsonlite::fromJSON(rawToChar(science_raw$manifest), simplifyVector = FALSE)
  finite_unique <- function(x, depth = 0L) {
    need(depth <= 100L)
    if (is.list(x)) {
      need(is.null(names(x)) || !anyDuplicated(names(x)))
      for (v in x) finite_unique(v, depth + 1L)
    } else if (is.numeric(x)) need(all(is.finite(x)))
    invisible(TRUE)
  }
  finite_unique(manifest)
  need(is.list(manifest) && !is.null(names(manifest)) && !anyDuplicated(names(manifest)) &&
    identical(manifest$schema, "GSE205906_FORMAL_RUN_MANIFEST_156/v1") &&
    identical(manifest$gene_order, "lexically_sorted_stable_gene_id") &&
    identical(unlist(manifest$sample_order, use.names = FALSE), ledger$label))
  gse205906_ens101_nofit_metadata_consistency_v2(manifest, ledger)
  ids <- unlist(manifest$gene_ids, use.names = FALSE)
  converted <- source$loader$gse205906_ens101_input_loader_v2(science_raw$matrix, ledger, ids,
    source_bundle$prefit, gtf_bytes = science_raw$GTF)
  prefit <- source$prefit$gse205906_ens101_effect_independent_prefit_v1(converted$inputs, cfg,
    source_bundle[c("model", "helpers", "annotation_adapter")], input_pins)
  allowed <- c("PASS_DESIGN_FILTER_TMM_ONLY", "HELD_ZERO_INPUT_LIBRARY", "UNVERIFIED_EDGER_FILTER",
    "HELD_MINIMUM_RETAINED", "UNVERIFIED_REFERENCE_MINIMUM", "UNVERIFIED_TMM_OR_PINNED_REFERENCE")
  scalar <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    need(is.numeric(x) && length(x) == 1L && (is.na(x) || is.finite(x))); as.numeric(x)
  }
  limits <- function(x) {
    if (is.null(x) || length(x) == 0L) return(list(minimum = NA_real_, maximum = NA_real_))
    need(is.numeric(x) && all(is.finite(x))); list(minimum = min(x), maximum = max(x))
  }
  phases <- lapply(seq_along(prefit$runs), function(i) {
    r <- prefit$runs[[i]]; need(r$preprocessing$status %in% allowed)
    list(phase_index = i, status = r$preprocessing$status,
      library_n = length(r$sample_ids), design_rank = scalar(r$rank), nominal_residual_df = scalar(r$nominal_residual_df),
      retained_gene_n = scalar(r$retained_gene_n), filter_by_expr_n = scalar(r$filter_by_expr_n), cpm_n = scalar(r$cpm_n),
      retained_unique_symbol_n = scalar(r$retained_unique_symbol_n), duplicate_symbol_excess = scalar(r$duplicate_symbol_excess),
      original_library_sizes = limits(r$original_library_sizes), filtered_library_sizes = limits(r$filtered_library_sizes),
      tmm_factors = limits(r$tmm_factors), phase_object_sha256 = r$phase_output_object_sha256)
  })
  summary <- list(schema = "GSE205906_ENS101_NOFIT_RAM_WORKER_SUMMARY/v1", status = prefit$status,
    integration_status = "NEEDS_REVIEW_REAL_FD_RUNTIME_RESOURCE_CALLER_NOT_IMPLEMENTED",
    passing_phase_n = prefit$passing_phase_n, failed_phase_n = prefit$failed_phase_n, phases = phases,
    source_input_gene_n = 32883L, source_input_library_n = 9L, original_metadata_column_n = 53L,
    nonblank_symbol_coverage = prefit$annotation_projection$metrics$nonblank_symbol_coverage,
    observed_six_raw_byte_identities = stats::setNames(observed, roles[1:6]),
    declared_runtime_pin = input_pins$statistical_runtime_lock,
    prefit_object_sha256 = prefit$output_object_sha256,
    transformation_object_sha256 = converted$transformation$transformation_object_sha256,
    supplied_bytes_match_declared_pins = TRUE, actual_file_origin_authenticated_here = FALSE,
    manifest_metadata_rows_equal_all_overlay_values = TRUE,
    metadata_consistency_observed_on_supplied_RAM_only = TRUE,
    runtime_lock_authenticated_here = FALSE, Python_FD_resource_and_child_caller_implemented_here = FALSE,
    membership_computed = FALSE, effect_or_dispersion_fit_called = FALSE, model_run = FALSE,
    scientific_acceptance = FALSE, production_go = FALSE, new_model_authorized = FALSE,
    vectors_metadata_or_private_errors_exported_by_summary = FALSE)
  need(nchar(jsonlite::toJSON(summary, auto_unbox = TRUE, na = "null"), type = "bytes") <= 32768L)
  # Complete objects deliberately remain inside caller RAM. Do not serialize the
  # whole return to stdout; only summary is designed for an external safe channel.
  list(prefit = prefit, transformation = converted$transformation, summary = summary)
}
