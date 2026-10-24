# Saved-only entry. Source defines functions; default CLI is a no-read dry-run.
# Compatibility entry: reuses the completed original HISAT preparation;
# --prepare-hisat is rejected. The original preparation entry remains intact.
# Rscript --vanilla THIS.R PROJECT_ROOT FRESH_OUTPUT --compare
#   --hisat-prepared OLD_LEAF --hisat-exit-receipt FILE --hisat-exit-sha SHA
#   --star-downstream STAR_LEAF --star-exit-receipt FILE --star-exit-sha SHA
# Existing frozen R/library only; never install, fit, source an old driver,
# read checked_payload.rds, TMM-normalize or invoke LOAO/null/RNG routines.
# Active modes require Linux, a fresh direct saved_* cloud leaf and threads=1.
# The future parent receipt is an explicit NOT-YET-PRODUCED contract:
# {schema:GSE205906_MATCHED_PARENT_PROCESS_EXIT/v1,run_root,actual_R_exit:0,
#  caller_source_sha256,complete_pin:{path,bytes,sha256}}.
# Parent must record an actually observed R terminal; this code creates none.

g205saved_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)

g205saved_pins <- function() list(
  core = "b7cac68a7ef3f5fc63ab727f7498b7b6d0c46b9bde18276affce3c87c4f3a71f",
  caller = "7cdafad1f838cd01ca8ae730762c059a78b2b3f515d47785fd09fe87f9e2275e",
  comparator = "eacfc026812bc40d25274415544fd2a2c10ad3cf5b7b46c5b678dd290520f599")

g205saved_pin <- function(path) {
  info <- file.info(path)
  g205saved_need(nrow(info) == 1L && !is.na(info$size) && !info$isdir &&
    identical(Sys.readlink(path), ""), "REGULAR_NONSYMLINK_SAVED_INPUT_REQUIRED")
  list(path = path, bytes = unname(info$size), sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

g205saved_raw <- function(path, sha, bytes = NULL) {
  pin <- g205saved_pin(path)
  g205saved_need(identical(pin$sha256, sha) && (is.null(bytes) || pin$bytes == bytes) &&
    pin$bytes <= 256 * 1024^2, "SAVED_INPUT_PIN_MISMATCH")
  raw <- readBin(path, "raw", n = pin$bytes)
  g205saved_need(length(raw) == pin$bytes && identical(pin, g205saved_pin(path)), "SAVED_INPUT_CHANGED_DURING_READ")
  list(raw = raw, pin = pin)
}

g205saved_functions <- function(raw, wanted, target) {
  expressions <- as.list(parse(text = rawToChar(raw), keep.source = FALSE))
  selected <- Filter(function(x) is.call(x) && length(x) == 3L && identical(x[[1L]], as.name("<-")) &&
    is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted, expressions)
  g205saved_need(setequal(vapply(selected, function(x) as.character(x[[2L]]), character(1)), wanted) &&
    length(selected) == length(wanted), "EXACT_FUNCTION_DECLARATIONS_REQUIRED")
  for (x in selected) {
    g205saved_need(is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), "FUNCTION_ONLY_IMPORT")
    eval(x, target)
  }
}

g205saved_modules <- function(project) {
  pins <- g205saved_pins()
  files <- c(core = "code/R/gse205906_matched_downstream_core_custody_v1_20261022.R",
    caller = "code/R/run_gse205906_matched_downstream_custody_v1_20261022.R",
    comparator = "code/R/gse205906_matched_saved_comparison_20261015.R")
  raw <- lapply(names(files), function(role) g205saved_raw(file.path(project, files[[role]]), pins[[role]]))
  names(raw) <- names(files)
  environments <- lapply(raw, function(held) {
    expressions <- as.list(parse(text = rawToChar(held$raw), keep.source = FALSE))
    wanted <- vapply(Filter(function(x) is.call(x) && length(x) == 3L && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), expressions),
      function(x) as.character(x[[2L]]), character(1))
    target <- new.env(parent = baseenv())
    g205saved_functions(held$raw, wanted, target); target
  })
  c(environments, list(module_pins = lapply(raw, `[[`, "pin")))
}

# Parse the frozen writer's quoted grammar without losing quotedness. In
# particular quoted "NA" is a character value; unquoted NA is missing. Tokens
# retain embedded tabs/newlines and doubled quotes. Roundtrip checks below bind
# every byte, column class and row, not only the displayed numeric fields.
g205saved_tsv <- function(raw, classes, writer) {
  g205saved_need(is.raw(raw) && length(raw) > 0L && !any(raw == as.raw(0L)), "TSV_RAW_REQUIRED")
  text <- rawToChar(raw)
  g205saved_need(!is.na(iconv(text, "UTF-8", "UTF-8", sub = NA)), "TSV_UTF8_REQUIRED")
  # Grammar delimiters are ASCII bytes. The complete text was validated above;
  # byte matching avoids the pathological multibyte whole-string PCRE path.
  # Restore UTF-8 labels, not contents, before per-field character operations.
  tokens <- regmatches(text, gregexpr('"(?:[^"]|"")*"|[^\\t\\n]+|[\\t\\n]', text,
    perl = TRUE, useBytes = TRUE))[[1L]]
  Encoding(tokens) <- "UTF-8"
  g205saved_need(length(tokens) %% 2L == 0L && identical(paste0(tokens, collapse = ""), text), "TSV_TOKEN_GRAMMAR")
  fields <- tokens[seq.int(1L, length(tokens), 2L)]
  separators <- tokens[seq.int(2L, length(tokens), 2L)]
  n <- length(classes)
  g205saved_need(n > 0L && !is.null(names(classes)) && !anyDuplicated(names(classes)) &&
    length(fields) %% n == 0L && identical(separators,
      rep(c(rep("\t", n - 1L), "\n"), length(fields) / n)), "TSV_RECORD_OR_COLUMN_STRUCTURE")
  quoted <- startsWith(fields, '"')
  g205saved_need(all(!quoted | endsWith(fields, '"')), "TSV_QUOTE_TERMINATOR")
  values <- fields
  values[quoted] <- gsub('""', '"', substring(fields[quoted], 2L, nchar(fields[quoted]) - 1L), fixed = TRUE)
  missing <- !quoted & fields == "NA"
  dim(values) <- c(n, length(values) / n); dim(missing) <- dim(values); dim(quoted) <- dim(values)
  g205saved_need(all(quoted[, 1L]) && identical(unname(values[, 1L]), names(classes)), "TSV_EXACT_HEADER_CLASSES")
  output <- lapply(seq_len(n), function(j) {
    cls <- unlist(classes[[j]], use.names = FALSE)
    g205saved_need(length(cls) == 1L && cls %in% c("character", "numeric", "integer", "logical", "list"),
      "UNSUPPORTED_MANIFEST_COLUMN_CLASS_NO_GUESSED_FACTOR_LEVELS")
    v <- unname(values[j, -1L]); na <- unname(missing[j, -1L]); q <- unname(quoted[j, -1L])
    if (cls == "character") {
      g205saved_need(all(q | na), "CHARACTER_TSV_VALUE_MUST_BE_QUOTED_OR_NA")
      v[na] <- NA_character_; return(v)
    }
    if (cls == "list") {
      g205saved_need(names(classes)[j] == "leadingEdge" && all(q) && !any(na), "ONLY_DECLARED_JSON_LEADING_EDGE_LIST_SUPPORTED")
      return(lapply(v, function(cell) {
        if (identical(cell, "null")) return(NULL)
        g205saved_need(startsWith(cell, "[") && endsWith(cell, "]"), "LEADING_EDGE_JSON_VECTOR_REQUIRED")
        z <- jsonlite::fromJSON(cell, simplifyVector = FALSE)
        g205saved_need(is.list(z) && is.null(names(z)) && all(vapply(z, function(x)
          is.null(x) || (is.character(x) && length(x) == 1L && !is.na(x)), logical(1))), "LEADING_EDGE_STRING_VECTOR_OR_NULL_ONLY")
        vapply(z, function(x) if (is.null(x)) NA_character_ else x, character(1))
      }))
    }
    g205saved_need(!any(q), "NUMERIC_LOGICAL_TSV_MUST_NOT_BE_QUOTED")
    if (cls == "logical") {
      g205saved_need(all(na | v %in% c("TRUE", "FALSE")), "LOGICAL_TSV_GRAMMAR")
      z <- v == "TRUE"; z[na] <- NA; return(z)
    }
    z <- suppressWarnings(as.numeric(v)); z[na] <- NA_real_
    g205saved_need(all(na | is.finite(z)) && all(!na | is.na(z)), "NUMERIC_TSV_GRAMMAR")
    if (cls == "integer") {
      g205saved_need(all(na | (z == floor(z) & abs(z) <= .Machine$integer.max)), "INTEGER_TSV_RANGE")
      return(as.integer(z))
    }
    z
  })
  names(output) <- names(classes)
  class(output) <- "data.frame"; attr(output, "row.names") <- .set_row_names(ncol(values) - 1L)
  g205saved_need(identical(writer(output), raw), "FROZEN_WRITER_FULL_TSV_BYTE_ROUNDTRIP_FAILED")
  output
}

g205saved_counts <- function(raw, ids, labels, samples) {
  g205saved_need(is.raw(raw) && all(as.integer(raw) %in% c(9L, 10L, 32:126)) &&
    identical(tail(raw, 1L), charToRaw("\n")), "OLD_MATRIX_ASCII_LF")
  rows <- strsplit(rawToChar(raw), "\n", fixed = TRUE)[[1L]]
  if (identical(tail(rows, 1L), "")) rows <- head(rows, -1L)
  fields <- strsplit(rows, "\t", fixed = TRUE)
  g205saved_need(length(fields) == length(ids) + 1L && all(lengths(fields) == length(samples) + 1L) &&
    identical(fields[[1L]], c("gene_id", labels)) && !anyDuplicated(ids) &&
    identical(vapply(fields[-1L], `[[`, character(1), 1L), ids), "COMPLETE_OLD_MATRIX_KEYS_AND_ORDER")
  strings <- unlist(lapply(fields[-1L], function(row) row[-1L]), use.names = FALSE)
  g205saved_need(all(grepl("^(0|[1-9][0-9]*)$", strings)), "CANONICAL_INTEGER_COUNT_STRINGS")
  value <- as.numeric(strings)
  g205saved_need(all(is.finite(value) & value <= .Machine$integer.max), "INTEGER_COUNTS_OVERFLOW")
  matrix(as.integer(value), nrow = length(ids), byrow = TRUE, dimnames = list(ids, samples))
}

g205saved_masks <- function(engine, counts, design) {
  # Execute ONLY these exact first four source assignments. The full frozen
  # preprocessing function is never called: no normLibSizes/TMM/dispersion,
  # QL/voom/DESeq, score or random-number operation occurs here.
  statements <- as.list(body(engine$v2r4_fit_preprocess))[-1L][seq_len(4L)]
  expected <- c("y0", "keep_fbe", "keep_cpm", "keep")
  g205saved_need(length(statements) == 4L && all(vapply(seq_len(4L), function(i) {
    x <- statements[[i]]
    is.call(x) && identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name(expected[i]))
  }, logical(1))), "FROZEN_FIRST_FOUR_MASK_ASSIGNMENTS_REQUIRED")
  local <- new.env(parent = environment(engine$v2r4_fit_preprocess))
  local$count_matrix <- counts; local$design_matrix <- design
  for (statement in statements) eval(statement, local)
  list(keep = local$keep, keep_filter_by_expr = local$keep_fbe, keep_cpm = local$keep_cpm,
    source_prefix_sha256 = digest::digest(statements, "sha256", serializeVersion = 3L),
    normalization_or_model_calls = 0L)
}

g205saved_check_design <- function(design, contrast, saved_design, saved_contrast) {
  g205saved_need(nrow(saved_design) == nrow(design$design) &&
    identical(saved_design$sample_id, rownames(design$design)) &&
    identical(saved_design$group, as.character(design$group)), "OLD_DESIGN_SAMPLE_GROUP_MISMATCH")
  for (name in colnames(design$design)) g205saved_need(is.numeric(saved_design[[name]]) &&
    identical(as.numeric(saved_design[[name]]), as.numeric(design$design[, name])), "OLD_DESIGN_COEFFICIENT_MISMATCH")
  g205saved_need(identical(saved_contrast$coefficient, rownames(contrast)), "OLD_CONTRAST_COEFFICIENT_ORDER")
  for (name in colnames(contrast)) g205saved_need(is.numeric(saved_contrast[[name]]) &&
    identical(as.numeric(saved_contrast[[name]]), as.numeric(contrast[, name])), "OLD_CONTRAST_VALUE_MISMATCH")
  TRUE
}

g205saved_check_masks <- function(counts, masks, tables, filter) {
  ids <- rownames(counts)
  for (key in c("keep", "keep_filter_by_expr", "keep_cpm")) g205saved_need(is.logical(masks[[key]]) &&
    length(masks[[key]]) == length(ids) && !anyNA(masks[[key]]), "COMPLETE_OLD_BOOLEAN_MASK_REQUIRED")
  g205saved_need(identical(unname(masks$keep), unname(masks$keep_filter_by_expr & masks$keep_cpm)), "OLD_MASK_AND_RELATION")
  g205saved_need(nrow(filter) == 1L && filter$input_gene_n == length(ids) &&
    filter$retained_gene_n == sum(masks$keep) && filter$filterByExpr_retained_n == sum(masks$keep_filter_by_expr) &&
    filter$explicit_CPM_retained_n == sum(masks$keep_cpm) && filter$minimum_cpm == 1 && filter$minimum_samples == 3,
    "OLD_MASK_COUNTS_DO_NOT_MATCH_SAVED_FILTER_SUMMARY")
  for (role in c("gene_edger", "gene_voom", "gene_deseq")) {
    tab <- tables[[role]]
    g205saved_need(setequal(unique(tab$contrast), c("local_disease", "intervention")), "TWO_OLD_GENE_CONTRASTS_REQUIRED")
    for (contrast in c("local_disease", "intervention")) g205saved_need(identical(
      tab$gene_id[tab$contrast == contrast], ids[masks$keep]), "OLD_RESTORED_MASK_DIFFERS_FROM_FULL_NATIVE_GENE_KEYS")
  }
  TRUE
}

g205saved_old <- function(project, api) {
  r <- api$caller; c <- api$core; constants <- r$g205run_constants()
  sources <- r$g205run_sources(project, c)
  engine <- c$g205m_engine(lapply(sources, `[[`, "raw"))
  g205saved_need(identical(as.character(getRversion()), engine$cfg$software$R) &&
    identical(as.character(utils::packageVersion("edgeR")), engine$cfg$software$edgeR), "FROZEN_R_EDGER_REQUIRED_FOR_MASK_ONLY")
  writer <- r$g205run_read(file.path(project, "code/R/gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R"), constants$writer_sha)
  parser <- r$g205run_read(file.path(project, "code/R/gse205906_ensembl101_nofit_RAM_worker_v2_20261006.R"), constants$parser_sha)
  adapter <- r$g205run_read(file.path(project, "code/R/gse205906_resource_object_validation_v1.R"), constants$adapter_sha)
  pure <- new.env(parent = baseenv())
  c$g205m_import(writer$raw, c("gse205906_ens101_worker_need_v1", "gse205906_ens101_worker_tsv_v1"), pure)
  c$g205m_import(parser$raw, c("gse205906_ens101_nofit_need_v2", "gse205906_ens101_nofit_overlay_v2", "gse205906_ens101_nofit_metadata_consistency_v2"), pure)
  c$g205m_import(adapter$raw, "gse205906_resource_label_adapter_v1", pure)
  file_manifest <- g205saved_raw(constants$hisat2_manifest, "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7", 14520)
  completed <- g205saved_raw(constants$hisat2_completion, "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e", 2137)
  manifest <- jsonlite::fromJSON(rawToChar(file_manifest$raw), simplifyVector = FALSE)
  encoding <- manifest$encoding
  g205saved_need(identical(manifest$schema, "GSE205906_ENS101_CHECKED_MODEL_FILE_MANIFEST/v1") &&
    identical(encoding$TSV, "UTF8_QUOTED_TEXT_17_SIGNIFICANT_DIGITS_NA_TOKEN") &&
    identical(encoding$list_column_TSV, "QUOTED_JSON_VECTOR_CELLS_OR_NULL"), "ORIGINAL_MANIFEST_ENCODING_REQUIRED")
  roles <- c("gene_edger", "gene_voom", "gene_deseq", "unsigned_camera", "unsigned_mroast", "unsigned_fgsea", "activity_results")
  extra <- c("mapping_coverage", "design_out", "contrast_out", "filter_out")
  held_tables <- lapply(c(roles, extra), function(role) {
    name <- paste0(role, ".tsv"); expected <- manifest$files[[name]]
    g205saved_need(!is.null(expected) && !is.null(encoding$table_column_classes[[role]]), "EXACT_OLD_TABLE_MANIFEST_RECORD_REQUIRED")
    held <- g205saved_raw(file.path(constants$hisat, name), expected$sha256, expected$bytes)
    list(table = g205saved_tsv(held$raw, encoding$table_column_classes[[role]], pure$gse205906_ens101_worker_tsv_v1),
      raw = held$raw, pin = held$pin)
  }); names(held_tables) <- c(roles, extra)
  matrix <- g205saved_raw(constants$hisat2_matrix, "e98296fdf308dbaa82325a3ca6b0b041c6f2569e98be75a45757141b0a1f77a2", 1475133)
  ledger_raw <- g205saved_raw(constants$metadata, "b0df6fa1354a83237bcfcd8bf8c323517deecd0e5680da8202489c4464936444", 10054)
  genes_raw <- g205saved_raw(constants$fixed_gene_manifest, "751e43b3aba8c2245981dc0041fec273c455c6390e85ef7ec4514e56118bf3b0")
  resource_raw <- g205saved_raw(constants$gene_resources, "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51", 1506178)
  meta_raw <- g205saved_raw(paste0(constants$run, "/metadata_preparation_01/METADATA_AND_RESERVOIR_BINDING.json"),
    "6dd7b53d40433baaacdeedf285204dee83f2d960c090e6a25f489f96d20fc5cb", 20661)
  ledger <- pure$gse205906_ens101_nofit_overlay_v2(ledger_raw$raw)
  genes <- jsonlite::fromJSON(rawToChar(genes_raw$raw), simplifyVector = FALSE)
  pure$gse205906_ens101_nofit_metadata_consistency_v2(genes, ledger)
  mapping_env <- new.env(parent = baseenv())
  c$g205m_import(sources$prefit$raw, c("gse205906_ens101_prefit_need_v1", "gse205906_ens101_prefit_ledger_v1"), mapping_env)
  mapping <- mapping_env$gse205906_ens101_prefit_ledger_v1(ledger)
  ids <- unlist(genes$gene_ids, use.names = FALSE)
  g205saved_need(length(ids) == 32883L && identical(ids, sort(ids, method = "radix")), "FULL_ORIGINAL_GENE_MANIFEST_REQUIRED")
  design <- engine$v2r4_design(ledger, engine$cfg)
  counts <- g205saved_counts(matrix$raw, ids, mapping$label, rownames(design$design))
  tables <- lapply(held_tables[roles], `[[`, "table")
  saved_design <- held_tables$design_out$table; saved_contrast <- held_tables$contrast_out$table
  g205saved_check_design(design, design$contrast, saved_design, saved_contrast)
  masks <- g205saved_masks(engine, counts, design$design)
  g205saved_check_masks(counts, masks, tables, held_tables$filter_out$table)
  resources <- pure$gse205906_resource_label_adapter_v1(readRDS(constants$gene_resources))
  fake_base <- c(masks[c("keep", "keep_filter_by_expr", "keep_cpm")], list(y = list(counts = counts[masks$keep, , drop = FALSE])))
  domains <- r$g205run_domains(list(counts = counts), list(native = list(base_fit = fake_base),
    slots = c$g205m_slots(tables), tables = tables, mapping_coverage = held_tables$mapping_coverage$table), resources)
  api$comparator$g205c_validate_branch(tables, domains)
  meta <- jsonlite::fromJSON(rawToChar(meta_raw$raw), simplifyVector = FALSE)
  sample_ids <- vapply(meta$libraries, `[[`, character(1), "GEO_sample")
  labels <- vapply(meta$libraries, `[[`, character(1), "label")
  runs <- vapply(meta$libraries, `[[`, character(1), "SRA_run")
  pairs <- stats::setNames(vapply(meta$libraries, function(x) as.numeric(x$registered_pairs), numeric(1)), sample_ids)
  g205saved_need(identical(sample_ids, colnames(counts)) && identical(labels, mapping$label) &&
    identical(runs, paste0("SRR", 19631207:19631199)) && all(is.finite(pairs) & pairs > 0 & pairs == floor(pairs)),
    "FULL_PAIR_DENOMINATOR_CROSSWALK_REQUIRED")
  tracked <- c(api$module_pins, lapply(sources, `[[`, "pin"), list(writer$pin, parser$pin, adapter$pin,
    file_manifest$pin, completed$pin, matrix$pin, ledger_raw$pin, genes_raw$pin, resource_raw$pin, meta_raw$pin),
    lapply(held_tables, `[[`, "pin"))
  list(schema = "GSE205906_HISAT_SAVED_COMPARISON_INPUTS/v1", counts = counts, tables = tables,
    source_domains = domains, source_universe = r$g205run_source_universe(resources), gene_ids = ids,
    input_pairs = pairs, design = design, masks = masks, input_pins = tracked,
    old_native_table_bytes_roundtrip = TRUE, old_mask_summary_and_two_contrast_gene_keys_match = TRUE,
    original_serialized_RDS_or_unused_factor_levels_authenticated_here = FALSE,
    old_historical_parent_exit_not_upgraded = TRUE, sample_unit = "SAME_NINE_SUBMITTED_LIBRARIES_NOT_AUTHENTICATED_ANIMALS",
    normalization_or_model_calls = 0L, source_tables = held_tables, writer = pure$gse205906_ens101_worker_tsv_v1)
}

g205saved_unchanged <- function(pins) {
  g205saved_need(all(vapply(pins, function(pin) identical(pin, g205saved_pin(pin$path)), logical(1))), "SAVED_SOURCE_OR_INPUT_CHANGED")
  TRUE
}

g205saved_rds <- function(value, path) {
  g205saved_need(!file.exists(path) && !file.exists(paste0(path, ".partial")), "NO_SAVED_OUTPUT_OVERWRITE")
  saveRDS(value, paste0(path, ".partial"), compress = "gzip", version = 3L)
  g205saved_need(file.rename(paste0(path, ".partial"), path), "SAVED_RDS_RENAME_FAILED")
  g205saved_pin(path)
}

g205saved_parent_exit <- function(path, sha, root, caller_sha) {
  held <- g205saved_raw(path, sha)
  value <- jsonlite::fromJSON(rawToChar(held$raw), simplifyVector = FALSE)
  g205saved_need(identical(value$schema, "GSE205906_MATCHED_PARENT_PROCESS_EXIT/v1") &&
    identical(value$run_root, root) && is.numeric(value$actual_R_exit) &&
    length(value$actual_R_exit) == 1L && value$actual_R_exit == 0 &&
    identical(value$caller_source_sha256, caller_sha) &&
    identical(value$complete_pin$path, file.path(root, "COMPLETE.json")), "EXPLICIT_ACTUAL_PARENT_R_EXIT_ZERO_REQUIRED")
  complete <- g205saved_raw(value$complete_pin$path, value$complete_pin$sha256, value$complete_pin$bytes)
  list(parent_pin = held$pin, complete_pin = complete$pin,
    complete = jsonlite::fromJSON(rawToChar(complete$raw), simplifyVector = FALSE))
}

g205saved_product <- function(complete, root, filename) {
  path <- file.path(root, filename)
  records <- Filter(function(x) identical(x$path, path), complete$products)
  g205saved_need(length(records) == 1L, "EXACT_SAVED_PRODUCT_PIN_REQUIRED")
  held <- g205saved_raw(path, records[[1L]]$sha256, records[[1L]]$bytes)
  held
}

# JSON may decode whole bytes as integer while file.info returns double.
# Compare byte values without accepting booleans or other coerced types.
g205saved_product_pin <- function(declared, actual) {
  valid_bytes <- function(x) is.numeric(x) && !is.object(x) &&
    length(x) == 1L && is.finite(x) && x >= 0 && x <= 2^53 && x == floor(x)
  g205saved_need(is.list(declared) && is.list(actual) &&
    is.character(declared$path) && length(declared$path) == 1L &&
    is.character(declared$sha256) && length(declared$sha256) == 1L &&
    identical(declared$path, actual$path) && identical(declared$sha256, actual$sha256) &&
    valid_bytes(declared$bytes) && valid_bytes(actual$bytes) &&
    declared$bytes == actual$bytes, "ACTUAL_STAR_CALLER_PRODUCT_PIN_CHANGED")
  actual
}

g205saved_compare <- function(options, api) {
  caller <- api$caller
  # Authenticate the already-finished HISAT preparation with its actual old actor,
  # never the new compatibility entry's SHA and never a new preparation run.
  legacy_entry <- g205saved_raw(file.path(options$project,
    "code/R/run_gse205906_matched_saved_comparison_20261015.R"),
    "62e821a18589a7f429f17d0b031d98a421a3bc075c3b5c97b86e97ec84e77595")
  old_exit <- g205saved_parent_exit(options$hisat_exit_receipt, options$hisat_exit_sha,
    options$hisat_prepared, legacy_entry$pin$sha256)
  new_exit <- g205saved_parent_exit(options$star_exit_receipt, options$star_exit_sha, options$star_downstream, g205saved_pins()$caller)
  g205saved_need(identical(old_exit$complete$status, "PREPARED_HISAT_SAVED_COMPARISON_INPUTS_NO_MODEL_FIT") &&
    identical(new_exit$complete$status, "COMPLETE_SAVED_NATIVE_NOT_INDEPENDENTLY_ACCEPTED") &&
    new_exit$complete$core_fit_invocations == 1 && identical(new_exit$complete$own_complete_source_domain_written, TRUE),
    "COMPLETE_OLD_PREPARATION_AND_NEW_SINGLE_CORE_RESULT_REQUIRED")
  expected_new_files <- c("SAVED_NATIVE_CHECKPOINT.rds", "STATISTICAL_RUNTIME.json",
    paste0(c("gene_edger", "gene_voom", "gene_deseq", "unsigned_camera", "unsigned_mroast", "unsigned_fgsea", "activity_results"), ".tsv"),
    "source_domains.tsv", "slots.tsv", "mapping_coverage.tsv", "native_not_returned_programs.tsv")
  product_paths <- vapply(new_exit$complete$products, `[[`, character(1), "path")
  g205saved_need(!anyDuplicated(product_paths) && setequal(product_paths, file.path(options$star_downstream, expected_new_files)),
    "ALL_ACTUAL_STAR_CALLER_PRODUCTS_REQUIRED")
  all_new_pins <- lapply(new_exit$complete$products, function(pin) {
    g205saved_product_pin(pin, g205saved_pin(pin$path))
  })
  old_held <- g205saved_product(old_exit$complete, options$hisat_prepared, "HISAT_SAVED_COMPARISON_INPUTS.rds")
  new_held <- g205saved_product(new_exit$complete, options$star_downstream, "SAVED_NATIVE_CHECKPOINT.rds")
  # Only the named new/old preparation checkpoints, never H/checked_payload.
  old <- readRDS(old_held$pin$path); checkpoint <- readRDS(new_held$pin$path)
  g205saved_need(identical(old$schema, "GSE205906_HISAT_SAVED_COMPARISON_INPUTS/v1") &&
    identical(checkpoint$schema, "GSE205906_MATCHED_SAVED_NATIVE_CHECKPOINT/v1") &&
    identical(checkpoint$statistical_fit_executed, TRUE), "EXPLICIT_SAVED_CHECKPOINT_SCHEMA")
  result <- checkpoint$result; packet <- checkpoint$prepared_packet
  g205saved_need(identical(result$schema, "gse205906-matched-downstream-result-v1") &&
    identical(packet$source_pins, api$core$g205m_pins()) &&
    identical(result$statistical_source_pins, api$core$g205m_pins()) &&
    identical(packet$binding$matrix_path, caller$g205run_constants()$matrix), "FROZEN_NEW_SOURCE_AND_MATRIX_BINDING_REQUIRED")
  writer_held <- caller$g205run_read(file.path(options$project, "code/R/gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R"), caller$g205run_constants()$writer_sha)
  pure <- new.env(parent = baseenv())
  api$core$g205m_import(writer_held$raw, c("gse205906_ens101_worker_need_v1", "gse205906_ens101_worker_tsv_v1"), pure)
  writer <- pure$gse205906_ens101_worker_tsv_v1
  new_domains <- g205saved_product(new_exit$complete, options$star_downstream, "source_domains.tsv")
  classes <- list(table = "character", collection = "character", contrast = "character", feature_id = "character",
    eligible = "logical", reason = "character", native_not_returned = "logical")
  domains <- g205saved_tsv(new_domains$raw, classes, writer)
  # Complete raw source universes must be present in each declared slot, not
  # merely returned/eligible rows. Both branches retain their own eligibility.
  for (branch in list(old$source_domains, domains)) for (i in seq_len(nrow(api$comparator$g205c_slots()))) {
    slot <- api$comparator$g205c_slots()[i, ]
    expected <- if (slot$collection == "genes") old$gene_ids else old$source_universe[[slot$collection]]
    actual <- branch$feature_id[branch$table == slot$table & branch$collection == slot$collection & branch$contrast == slot$contrast]
    g205saved_need(identical(actual, expected), "FULL_SOURCE_DOMAIN_ORDER_NOT_RETURNED_ONLY")
  }
  tracked <- c(list(legacy_entry$pin), api$module_pins, all_new_pins, list(old_exit$parent_pin, old_exit$complete_pin, new_exit$parent_pin, new_exit$complete_pin,
    old_held$pin, new_held$pin, new_domains$pin, writer_held$pin))
  for (role in names(old$tables)) {
    held <- g205saved_product(new_exit$complete, options$star_downstream, paste0(role, ".tsv"))
    g205saved_need(identical(writer(caller$g205run_tsv_view(result$tables[[role]])), held$raw), "NEW_NATIVE_TABLE_EXPORT_NOT_BYTE_IDENTICAL")
    tracked[[length(tracked) + 1L]] <- held$pin
  }
  new_pairs <- stats::setNames(vapply(checkpoint$custody$libraries, function(x) as.numeric(x$input_pairs), numeric(1)),
    vapply(checkpoint$custody$libraries, `[[`, character(1), "sample_id"))
  g205saved_need(identical(new_pairs, old$input_pairs), "SAME_NINE_FULL_PAIR_DENOMINATORS_REQUIRED")
  counts <- api$comparator$g205c_counts(old$counts, packet$counts, old$gene_ids, old$input_pairs)
  tables <- api$comparator$g205c_tables(old$tables, result$tables, old$source_domains, domains,
    topk_policies = NULL, caller_binding = list(input_pins = tracked, parent_exit_records_verified = TRUE))
  g205saved_unchanged(tracked)
  list(counts = counts, tables = tables, input_pins = tracked, writer = writer,
    normalization_or_model_calls = 0L, top_k = "PENDING_POLICY_BINDING", independent_numeric_acceptance = FALSE)
}

g205saved_options <- function(args) {
  g205saved_need(length(args) >= 2L, "EXPLICIT_PROJECT_OUTPUT_REQUIRED")
  out <- list(project = args[[1L]], output = args[[2L]], mode = "dry-run")
  flags <- c("--hisat-prepared", "--hisat-exit-receipt", "--hisat-exit-sha", "--star-downstream", "--star-exit-receipt", "--star-exit-sha")
  i <- 3L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--prepare-hisat", "--compare")) {
      g205saved_need(out$mode == "dry-run", "ONE_SAVED_MODE_ONLY")
      out$mode <- substring(key, 3L); i <- i + 1L
    } else {
      g205saved_need(key %in% flags && i < length(args), "UNKNOWN_SAVED_OPTION")
      name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
      g205saved_need(is.null(out[[name]]), "DUPLICATE_SAVED_OPTION")
      out[[name]] <- args[[i + 1L]]; i <- i + 2L
    }
  }
  g205saved_need(!grepl("(^|/|\\\\)[.][.](/|\\\\|$)", out$output), "NO_SAVED_OUTPUT_TRAVERSAL")
  if (out$mode == "compare") for (key in flags) {
    name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    g205saved_need(is.character(out[[name]]) && length(out[[name]]) == 1L && nzchar(out[[name]]), "EXPLICIT_OLD_NEW_LEAVES_AND_PARENT_EXIT_PINS_REQUIRED")
  }
  out
}

g205saved_main <- function(args) {
  options <- g205saved_options(args)
  if (options$mode == "dry-run") return(list(status = "DRY_RUN_NO_SCIENTIFIC_READ_WRITE_OR_COMPUTATION",
    mode = "dry-run", old_preparation = "REUSE_EXISTING_ONLY_NOT_RECOMPUTED", future_comparison = "FROZEN_SAVED_ONLY_COUNTS_AND_TABLES", top_k = "PENDING_POLICY_BINDING"))
  g205saved_need(options$mode == "compare", "COMPATIBILITY_ENTRY_REUSES_EXISTING_HISAT_PREPARATION_ONLY")
  g205saved_need(identical(unname(Sys.info()["sysname"]), "Linux") &&
    all(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "R_DATATABLE_NUM_THREADS")) == "1"), "EXISTING_SINGLE_THREAD_LINUX_RUNTIME_REQUIRED")
  cloud <- "/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01"
  g205saved_need(dirname(options$output) == cloud && startsWith(basename(options$output), "saved_") &&
    !file.exists(options$output) && dir.exists(cloud), "FRESH_EXPLICIT_SAVED_CLOUD_LEAF_REQUIRED")
  if (options$mode == "compare") g205saved_need(dirname(options$hisat_prepared) == cloud &&
    startsWith(basename(options$hisat_prepared), "saved_") && dirname(options$star_downstream) == cloud &&
    startsWith(basename(options$star_downstream), "downstream_"), "EXPLICIT_BOUND_CLOUD_SAVED_INPUT_LEAVES_REQUIRED")
  api <- g205saved_modules(options$project)
  self_pin <- g205saved_pin(file.path(options$project, "code/R/run_gse205906_matched_saved_comparison_custody_v2_20261022.R"))
  rng <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) digest::digest(get(".Random.seed", .GlobalEnv), "sha256") else NULL
  started <- proc.time()[["elapsed"]]; owned <- FALSE; phase <- "SAVED_INPUT_PREPARATION"
  tryCatch(withCallingHandlers({
    if (options$mode == "prepare-hisat") {
      result <- g205saved_old(options$project, api)
      result$input_pins <- c(list(self_pin), result$input_pins)
      writer <- result$writer; originals <- result$source_tables
      result$writer <- NULL; result$source_tables <- NULL
    } else {
      result <- g205saved_compare(options, api)
      result$input_pins <- c(list(self_pin), result$input_pins)
      writer <- result$writer; result$writer <- NULL
    }
    g205saved_unchanged(result$input_pins)
    after_rng <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) digest::digest(get(".Random.seed", .GlobalEnv), "sha256") else NULL
    g205saved_need(identical(rng, after_rng), "SAVED_ONLY_STAGE_CHANGED_RNG")
    phase <- "OWNED_FRESH_SAVED_EXPORT"
    g205saved_need(!file.exists(options$output) && dir.create(options$output, recursive = FALSE), "FRESH_SAVED_OUTPUT_CREATION_FAILED")
    owned <- TRUE
    products <- list()
    if (options$mode == "prepare-hisat") {
      products[[1L]] <- g205saved_rds(result, file.path(options$output, "HISAT_SAVED_COMPARISON_INPUTS.rds"))
      for (role in names(result$tables)) products[[length(products) + 1L]] <- api$caller$g205run_write_raw(
        originals[[role]]$raw, file.path(options$output, paste0(role, ".tsv")))
      exports <- list(source_domains = result$source_domains,
        full_counts_GSM = data.frame(gene_id = result$gene_ids, result$counts, check.names = FALSE),
        fixed_pairs = data.frame(sample_id = names(result$input_pairs), input_pairs = unname(result$input_pairs)))
      status <- "PREPARED_HISAT_SAVED_COMPARISON_INPUTS_NO_MODEL_FIT"
    } else {
      products[[1L]] <- g205saved_rds(result, file.path(options$output, "SAVED_COMPARISON.rds"))
      exports <- list(count_summary = result$counts$summary, count_MA = result$counts$MA,
        native_table_summary = result$tables$summary, full_source_union = result$tables$union)
      status <- "COMPLETE_SAVED_COMPARISON_NOT_INDEPENDENTLY_ACCEPTED"
    }
    for (name in names(exports)) products[[length(products) + 1L]] <- api$caller$g205run_write_raw(
      writer(api$caller$g205run_tsv_view(exports[[name]])), file.path(options$output, paste0(name, ".tsv")))
    g205saved_unchanged(result$input_pins)
    receipt <- list(schema = "GSE205906_MATCHED_SAVED_ENTRY/v1", status = status, mode = options$mode,
      input_pins = result$input_pins, products = products, elapsed_seconds = proc.time()[["elapsed"]] - started,
      old_counts_mask_recovery = if (options$mode == "prepare-hisat") "DGEList_FILTERBYEXPR_CPM_AND_ONLY" else "NOT_RECOMPUTED_IN_COMPARISON",
      model_or_TMM_or_dispersion_or_scores_or_null_calls = 0L, RNG_unchanged = TRUE,
      top_k = "PENDING_POLICY_BINDING", raw_condition_text_persisted = FALSE,
      same_nine_libraries_not_independent_estimates = TRUE, candidate_admission = FALSE,
      parent_must_record_actual_R_exit = TRUE, independent_numeric_acceptance = FALSE)
    api$caller$g205run_write_json(receipt, file.path(options$output, "COMPLETE.json"))
    list(status = status, output = options$output, models = 0L, independent_numeric_acceptance = FALSE)
  }, warning = function(condition) stop("SAVED_ONLY_WARNING_CLASS_REQUIRES_REVIEW", call. = FALSE),
    message = function(condition) invokeRestart("muffleMessage")), error = function(condition) {
      failure <- list(status = "FAILED", phase = phase, exception_class_only = class(condition)[[1L]], raw_condition_text_persisted = FALSE)
      if (owned) tryCatch(api$caller$g205run_write_json(failure, file.path(options$output, "FAILED.json")), error = function(ignored) NULL)
      failure
    })
}

if (sys.nframe() == 0L) {
  result <- tryCatch(g205saved_main(commandArgs(trailingOnly = TRUE)), error = function(condition)
    list(status = "FAILED_BEFORE_OWNED_OUTPUT", exception_class_only = class(condition)[[1L]], raw_condition_text_persisted = FALSE))
  if (requireNamespace("jsonlite", quietly = TRUE)) cat(as.character(jsonlite::toJSON(result, auto_unbox = TRUE)), "\n")
  else cat(result$status, "\n")
  quit(status = if (startsWith(result$status, "FAILED")) 1L else 0L, save = "no")
}





