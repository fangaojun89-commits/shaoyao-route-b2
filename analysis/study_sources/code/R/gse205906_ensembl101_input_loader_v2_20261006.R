# NEW pure in-memory input conversion, NOT a production caller or admission gate.
# No file/network IO, data selection, fitting, automatic database fallback or writes.
gse205906_ens101_loader_need_v2 <- function(condition, category) {
  if (!isTRUE(condition)) stop(category, call. = FALSE)
}

gse205906_ens101_loader_hash_v2 <- function(object) {
  digest::digest(serialize(object, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_gtf_attributes_v2 <- function(gtf_bytes) {
  need <- gse205906_ens101_loader_need_v2
  need(is.raw(gtf_bytes) && length(gtf_bytes) > 0L && length(gtf_bytes) <= 512 * 1024^2,
       "ENS101_LOADER_GTF_BYTES_OR_ENGINEERING_CAPACITY_UNVERIFIED")
  need(!any(gtf_bytes == as.raw(0L)), "ENS101_LOADER_GTF_NUL")
  text <- rawToChar(gtf_bytes)
  need(!is.na(iconv(text, from = "UTF-8", to = "UTF-8", sub = NA)),
       "ENS101_LOADER_GTF_UTF8")
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  if (length(lines) > 0L && identical(tail(lines, 1L), "")) lines <- head(lines, -1L)
  lines <- sub("\r$", "", lines)
  fields <- strsplit(lines[!startsWith(lines, "#")], "\t", fixed = TRUE)
  need(length(fields) > 0L && all(lengths(fields) == 9L), "ENS101_LOADER_GTF_NINE_FIELDS")
  gene_rows <- which(vapply(fields, function(row) identical(row[[3L]], "gene"), logical(1)))
  need(length(gene_rows) > 0L, "ENS101_LOADER_GTF_NO_GENE_RECORDS")
  pattern <- '([A-Za-z_][A-Za-z0-9_]*)[ \t]+"([^"]*)";'
  gene_attributes <- lapply(seq_along(gene_rows), function(ordinal) {
    i <- gene_rows[[ordinal]]
    row <- fields[[i]]
    raw_attributes <- row[[9L]]
    need(identical(trimws(gsub(pattern, "", raw_attributes, perl = TRUE)), ""),
         "ENS101_LOADER_GTF_ATTRIBUTE_GRAMMAR_UNVERIFIED")
    tokens <- regmatches(raw_attributes, gregexpr(pattern, raw_attributes, perl = TRUE))[[1L]]
    parsed <- lapply(tokens, function(token) regmatches(token, regexec(pattern, token, perl = TRUE))[[1L]])
    keys <- vapply(parsed, function(token) token[[2L]], character(1))
    values <- vapply(parsed, function(token) token[[3L]], character(1))
    need(!anyDuplicated(keys), "ENS101_LOADER_GTF_DUPLICATE_ATTRIBUTE_KEY")
    need(sum(keys == "gene_id") == 1L && sum(keys == "gene_biotype") == 1L &&
           sum(keys == "gene_name") <= 1L, "ENS101_LOADER_GTF_REQUIRED_GENE_ATTRIBUTES")
    id <- values[keys == "gene_id"]
    biotype <- values[keys == "gene_biotype"]
    need(grepl("^ENSRNOG[0-9]+$", id) && nzchar(biotype) && identical(biotype, trimws(biotype)),
         "ENS101_LOADER_GTF_GENE_ID_OR_BIOTYPE")
    c(
      gene_id = id, gene_name = if (sum(keys == "gene_name") == 0L) NA_character_ else values[keys == "gene_name"],
      gene_biotype = biotype, source_gene_record_ordinal = as.character(ordinal),
      source_noncomment_record_ordinal = as.character(i), seqname = row[[1L]], source = row[[2L]],
      feature = row[[3L]], start = row[[4L]], end = row[[5L]], score = row[[6L]], strand = row[[7L]],
      frame = row[[8L]], raw_gtf_attributes = raw_attributes
    )
  })
  result <- as.data.frame(do.call(rbind, gene_attributes), stringsAsFactors = FALSE)
  need(!anyDuplicated(result$gene_id), "ENS101_LOADER_GTF_DUPLICATE_GENE_RECORD")
  rownames(result) <- NULL
  result
}

gse205906_ens101_input_loader_v2 <- function(
    matrix_bytes, sample_ledger, fixed_gene_ids, prefit_source_bytes,
    gene_attributes = NULL, gtf_bytes = NULL) {
  need <- gse205906_ens101_loader_need_v2
  need(requireNamespace("digest", quietly = TRUE), "ENS101_LOADER_EXISTING_DIGEST_REQUIRED_NO_INSTALL")
  # Reuse ONLY the exact frozen pure metadata checker, never caller callbacks.
  need(is.raw(prefit_source_bytes) && length(prefit_source_bytes) == 18762L &&
         identical(digest::digest(prefit_source_bytes, "sha256", serialize = FALSE),
                   "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
       "ENS101_LOADER_PREFIT_SOURCE_IDENTITY")
  environment <- new.env(parent = baseenv())
  allowed <- c("gse205906_ens101_prefit_need_v1", "gse205906_ens101_prefit_ledger_v1")
  selected <- Filter(function(expr) is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
                       is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% allowed,
                     as.list(parse(text = rawToChar(prefit_source_bytes), keep.source = FALSE)))
  found <- vapply(selected, function(expr) as.character(expr[[2L]]), character(1))
  need(setequal(found, allowed) && !anyDuplicated(found), "ENS101_LOADER_METADATA_CHECKER_COVERAGE")
  for (expr in selected) {
    need(is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function")),
         "ENS101_LOADER_METADATA_CHECKER_DECLARATIONS")
    eval(expr, envir = environment)
  }
  lockEnvironment(environment, bindings = TRUE)
  mapping <- environment$gse205906_ens101_prefit_ledger_v1(sample_ledger)
  need(is.character(fixed_gene_ids) && is.null(dim(fixed_gene_ids)) &&
         length(fixed_gene_ids) == 32883L && !anyNA(fixed_gene_ids) && !anyDuplicated(fixed_gene_ids) &&
         all(grepl("^ENSRNOG[0-9]+$", fixed_gene_ids)) &&
         identical(fixed_gene_ids, sort(fixed_gene_ids, method = "radix")),
       "ENS101_LOADER_FIXED_COMPLETE_MANIFEST_ORDER")
  need(is.raw(matrix_bytes) && length(matrix_bytes) > 0L && length(matrix_bytes) <= 32 * 1024^2,
       "ENS101_LOADER_MATRIX_BYTES_OR_ENGINEERING_CAPACITY_UNVERIFIED")
  need(all(as.integer(matrix_bytes) %in% c(9L, 10L, 32:126)) &&
         identical(tail(matrix_bytes, 1L), charToRaw("\n")),
       "ENS101_LOADER_ORIGINAL_ASCII_LF_TSV")
  matrix_text <- rawToChar(matrix_bytes)
  lines <- strsplit(matrix_text, "\n", fixed = TRUE)[[1L]]
  if (identical(tail(lines, 1L), "")) lines <- head(lines, -1L)
  need(length(lines) == 32884L, "ENS101_LOADER_MATRIX_HEADER_PLUS_32883_ROWS")
  columns <- strsplit(lines, "\t", fixed = TRUE)
  need(all(lengths(columns) == 10L) &&
         identical(columns[[1L]], c("gene_id", mapping$label)),
       "ENS101_LOADER_ORIGINAL_MATRIX_HEADER_OR_FIELD_COUNT")
  data_rows <- columns[-1L]
  observed_ids <- vapply(data_rows, function(row) row[[1L]], character(1))
  need(identical(observed_ids, fixed_gene_ids), "ENS101_LOADER_NO_IMPLICIT_GENE_REORDER_OR_DROP")
  count_strings <- unlist(lapply(data_rows, function(row) row[-1L]), use.names = FALSE)
  need(all(grepl("^(0|[1-9][0-9]*)$", count_strings)),
       "ENS101_LOADER_CANONICAL_NONNEGATIVE_INTEGER_TEXT")
  numeric_counts <- as.numeric(count_strings)
  need(all(is.finite(numeric_counts)) && all(numeric_counts <= .Machine$integer.max),
       "ENS101_LOADER_ORIGINAL_R_INTEGER_BOUND")
  count_matrix <- matrix(as.integer(numeric_counts), nrow = 32883L, ncol = 9L, byrow = TRUE,
                         dimnames = list(fixed_gene_ids, mapping$GEO_sample))
  need(xor(is.null(gene_attributes), is.null(gtf_bytes)), "ENS101_LOADER_EXACTLY_ONE_BOUND_ANNOTATION_INPUT")
  annotation_kind <- if (is.null(gene_attributes)) "EXPLICIT_BOUND_PLAIN_GTF_UTF8_GENE_FEATURES" else "EXPLICIT_BOUND_GTF_GENE_ATTRIBUTE_TABLE"
  attributes <- if (is.null(gene_attributes)) gse205906_ens101_gtf_attributes_v2(gtf_bytes) else gene_attributes
  need(is.data.frame(attributes) && !anyDuplicated(names(attributes)) &&
         all(c("gene_id", "gene_name", "gene_biotype") %in% names(attributes)) &&
         !any(c("GeneID", "gene") %in% names(attributes)) && nrow(attributes) == 32883L,
       "ENS101_LOADER_COMPLETE_RAW_GTF_ATTRIBUTE_SCHEMA")
  need(all(vapply(attributes[c("gene_id", "gene_name", "gene_biotype")],
                  function(x) is.character(x) && is.null(dim(x)) && length(x) == 32883L, logical(1))) &&
         !anyNA(attributes$gene_id) && !anyDuplicated(attributes$gene_id) &&
         !anyNA(attributes$gene_biotype) && all(nzchar(attributes$gene_biotype)) &&
         identical(attributes$gene_biotype, trimws(attributes$gene_biotype)) &&
         setequal(attributes$gene_id, fixed_gene_ids), "ENS101_LOADER_ATTRIBUTE_IDENTITY_OR_BIOTYPE_CONFLICT")
  positions <- match(fixed_gene_ids, attributes$gene_id)
  # GTF coordinate order is not count/manifest order. This explicit complete
  # one-to-one key projection records the source row positions, never sorts the
  # count matrix or removes/merges any gene record.
  annotation <- attributes[positions, , drop = FALSE]
  rownames(annotation) <- NULL
  metadata <- unserialize(serialize(sample_ledger, NULL, version = 3L))
  inputs <- list(count_matrix = count_matrix, gene_annotation = annotation,
                 sample_ledger = metadata, fixed_gene_ids = fixed_gene_ids)
  label_mapping <- data.frame(original_matrix_label = mapping$label, GSM_sample = mapping$GEO_sample,
                              SRA_run = mapping$SRA_run, source_group = mapping$group, stringsAsFactors = FALSE)
  transformation <- list(
    schema = "GSE205906_ENS101_LOSSLESS_INPUT_CONVERSION/v1",
    original_matrix_header = columns[[1L]], column_mapping = label_mapping,
    rule = "RENAME_ONLY_FIXED_LABEL_TO_GSM_WITHOUT_COLUMN_REORDER",
    gene_row_order_unchanged = TRUE, sample_column_order_unchanged = TRUE,
    annotation_rule = "EXPLICIT_COMPLETE_ONE_TO_ONE_GTF_GENE_ID_JOIN_TO_MANIFEST_ORDER",
    annotation_source_row_positions = positions,
    annotation_input_kind = annotation_kind,
    input_gene_n = 32883L, output_gene_n = nrow(count_matrix),
    input_library_n = 9L, output_library_n = ncol(count_matrix),
    count_values_changed = FALSE, genes_removed_or_merged = 0L, samples_removed = 0L,
    source_gene_names_preserved_no_fallback = TRUE,
    original_matrix_bytes = list(bytes = length(matrix_bytes),
                                sha256 = digest::digest(matrix_bytes, "sha256", serialize = FALSE)),
    original_attribute_object_sha256 = gse205906_ens101_loader_hash_v2(attributes),
    original_metadata_object_sha256 = gse205906_ens101_loader_hash_v2(sample_ledger),
    original_manifest_object_sha256 = gse205906_ens101_loader_hash_v2(fixed_gene_ids),
    source_gtf_bytes = if (is.null(gtf_bytes)) NULL else
      list(bytes = length(gtf_bytes), sha256 = digest::digest(gtf_bytes, "sha256", serialize = FALSE)),
    prefit_source_bytes_authenticated_here = list(bytes = 18762L,
      sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
    external_file_origin_or_qc_admission_authenticated_here = FALSE,
    reference_and_statistical_runtime_authenticated_here = FALSE,
    animal_pool_batch_dose_resolved = FALSE, new_model_authorized = FALSE, model_run = FALSE
  )
  transformation$transformation_object_sha256 <- gse205906_ens101_loader_hash_v2(transformation)
  list(inputs = inputs, transformation = transformation,
       input_object_identities = lapply(inputs, gse205906_ens101_loader_hash_v2),
       status = "PASS_LOSSLESS_FIXED_SCHEMA_CONVERSION_ONLY_NOT_SCIENTIFIC_ADMISSION",
       missing_symbol_coverage_not_tested_here = TRUE, prefit_or_model_executed = FALSE)
}
