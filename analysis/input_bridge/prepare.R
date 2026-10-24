# Portable input preparation only. Original readers remain byte-for-byte unchanged.
# Arguments: PUBLIC_REPOSITORY NEW_WORKSPACE OPTIONAL_EXISTING_R_LIBRARY FROZEN_GTF_GZIP
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
repo <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
work <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
if (nzchar(args[[3L]])) .libPaths(c(normalizePath(args[[3L]], mustWork = TRUE), .libPaths()))
required <- c(data.table = "1.18.4", openssl = "2.4.2", jsonlite = "2.0.0", digest = "0.6.39")
for (pkg in names(required)) stopifnot(requireNamespace(pkg, quietly = TRUE),
  identical(as.character(utils::packageVersion(pkg)), required[[pkg]]))
data.table::setDTthreads(1L)
had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
old_rng <- RNGkind()
contract <- jsonlite::fromJSON(file.path(repo, "analysis/input_bridge/CONTRACT.json"), simplifyVector = FALSE)
raw_file <- function(path) readBin(path, "raw", n = file.info(path)$size)
source_raw <- list()
for (record in contract$sources) {
  value <- raw_file(file.path(repo, record$path))
  stopifnot(length(value) == record$bytes,
    identical(digest::digest(value, "sha256", serialize = FALSE), record$sha256))
  source_raw[[record$role]] <- value
}
definitions <- function(raw, wanted = NULL) {
  nodes <- as.list(parse(text = rawToChar(raw), keep.source = FALSE))
  if (!is.null(wanted)) nodes <- Filter(function(x) is.call(x) &&
    identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted, nodes)
  target <- new.env(parent = baseenv())
  target$head <- utils::head; target$tail <- utils::tail; target$setNames <- stats::setNames
  names <- character()
  for (node in nodes) {
    stopifnot(is.call(node), identical(node[[1L]], as.name("<-")), is.symbol(node[[2L]]),
      is.call(node[[3L]]), identical(node[[3L]][[1L]], as.name("function")))
    names <- c(names, as.character(node[[2L]])); eval(node, target)
  }
  stopifnot(!anyDuplicated(names), is.null(wanted) || setequal(names, wanted))
  lockEnvironment(target, bindings = TRUE); target
}
g128 <- definitions(source_raw$g128_loader)$g128_kernel_loader(source_raw$g128_analyzer, source_raw$g128_validator)
inputs <- file.path(work, "original_inputs")
g128_root <- file.path(inputs, "gse128623_real_run_20261007r04")
ledger_path <- file.path(g128_root, "analysis_inputs/sample_identity_design_ledger.tsv")
annotation_path <- file.path(g128_root, "analysis_inputs/feature_annotation_from_existing_counts.tsv")
summaries <- list()
for (branch in c("raw", "adapter")) {
  matrix_path <- file.path(g128_root, "derived/assembled", branch, "all30_selected.tsv")
  a <- g128$analyzer$v2a_read_inputs(matrix_path, annotation_path, ledger_path, ledger_path,
    "gene_id", "GeneID", "gene")
  v <- g128$validator$v2v_read_input(matrix_path, annotation_path, ledger_path, "gene_id", "GeneID", "gene")
  stopifnot(identical(a$counts, v$counts), identical(a$symbols, v$symbols),
    identical(as.data.frame(a$annotation), as.data.frame(v$annotation)),
    identical(as.data.frame(a$ledger), as.data.frame(v$ledger)))
  design <- g128$analyzer$v2a_design(a$ledger$group)
  validator_design <- g128$validator$v2v_design(a$ledger$group)
  # The original functions name their local factors differently; compare the
  # scientific design values and axes, not the local-variable contrasts attribute.
  stopifnot(identical(as.vector(design), as.vector(validator_design)),
    identical(dimnames(design), dimnames(validator_design)), identical(dim(design), dim(validator_design)))
  small <- data.table::fread(file.path(g128_root, "derived/assembled", branch, "complete29_selected.tsv"))
  kept <- a$ledger$sample_id != "GSM3681230"
  stopifnot(identical(names(small)[-1L], colnames(a$counts)[kept]), identical(small$gene_id, rownames(a$counts)))
  values <- as.matrix(small[, -1L, with = FALSE]); storage.mode(values) <- "double"
  stopifnot(identical(unname(values), unname(a$counts[, kept, drop = FALSE])))
  output <- list(input = a, prespecified_29_columns = which(kept), design = design,
    contrast_matrix = g128$analyzer$v2a_contrast_matrix(), role = branch,
    specification = "UNFILTERED_INPUTS_NO_NORMALIZATION_OR_MODEL_FITTING")
  saveRDS(output, file.path(work, "prepared", paste0("G128_", branch, "_inputs.rds")), version = 3L)
  restored <- readRDS(file.path(work, "prepared", paste0("G128_", branch, "_inputs.rds")))
  stopifnot(identical(restored$input$counts, a$counts), identical(restored$input$symbols, a$symbols),
    identical(as.data.frame(restored$input$annotation), as.data.frame(a$annotation)),
    identical(as.data.frame(restored$input$ledger), as.data.frame(a$ledger)),
    identical(restored$design, design), identical(restored$prespecified_29_columns, which(kept)))
  summaries[[paste0("G128_", branch)]] <- list(genes = nrow(a$counts), libraries = ncol(a$counts),
    original_analyzer_validator_readers_agree = TRUE, exact_29_column_subset = TRUE,
    omitted_sample = "GSM3681230", design_rank = qr(design)$rank,
    symbols_missing = sum(is.na(a$symbols)))
  rm(a, v, small, values, output, restored); invisible(gc())
}
metadata_api <- definitions(source_raw$g205_metadata,
  c("gse205906_ens101_nofit_need_v2", "gse205906_ens101_nofit_overlay_v2"))
loader <- definitions(source_raw$g205_loader)
overlay <- raw_file(file.path(inputs,
  "GSE205906_full_read_counts_fasp_v18e_20261005/formal/derived/NINE_LIBRARY_ANALYSIS_PREPARATION_OVERLAY_154.tsv"))
ledger <- metadata_api$gse205906_ens101_nofit_overlay_v2(overlay)
attrs <- utils::read.delim(file.path(work, "derived/G205_gene_attributes.tsv"), sep = "\t", quote = '"',
  colClasses = "character", na.strings = "__ROUTE_B2_MISSING_GENE_NAME__", check.names = FALSE,
  stringsAsFactors = FALSE, comment.char = "", row.names = NULL)
# Compare every extracted field to the original GTF parser, in bounded chunks.
# Only its two source-record ordinals are offset to the full file coordinates.
gtf_con <- gzfile(args[[4L]], open = "rt", encoding = "UTF-8")
gene_offset <- 0L; record_offset <- 0L; annotation_chunks <- 0L
repeat {
  lines <- readLines(gtf_con, n = 10000L, warn = FALSE)
  if (!length(lines)) break
  records <- lines[!startsWith(lines, "#")]
  if (any(grepl("^[^\t]*\t[^\t]*\tgene\t", records))) {
    piece <- loader$gse205906_ens101_gtf_attributes_v2(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))))
    piece$source_gene_record_ordinal <- as.character(as.integer(piece$source_gene_record_ordinal) + gene_offset)
    piece$source_noncomment_record_ordinal <- as.character(as.integer(piece$source_noncomment_record_ordinal) + record_offset)
    expected <- attrs[gene_offset + seq_len(nrow(piece)), , drop = FALSE]
    rownames(expected) <- NULL
    stopifnot(identical(piece, expected))
    gene_offset <- gene_offset + nrow(piece)
    annotation_chunks <- annotation_chunks + 1L
  }
  record_offset <- record_offset + length(records)
}
close(gtf_con)
stopifnot(gene_offset == 32883L, gene_offset == nrow(attrs))
annotation_qa <- list(all_gene_fields_match_original_parser = TRUE,
  genes = gene_offset, noncomment_records = record_offset, fields = ncol(attrs),
  bounded_chunks_with_genes = annotation_chunks)
paths <- c(HISAT2_s0 = "GSE205906_full_read_counts_fasp_v18e_20261005/formal/output/GSE205906_fixed9_gene_fragments_s0.tsv",
  STAR_s0 = "GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s0.tsv",
  STAR_s1 = "GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s1.tsv",
  STAR_s2 = "GSE205906_matched_STAR_Ensembl101_20261015r01/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s2.tsv")
fixed_ids <- NULL
for (name in names(paths)) {
  path <- file.path(inputs, paths[[name]])
  original <- data.table::fread(path)
  if (is.null(fixed_ids)) fixed_ids <- original$gene_id
  stopifnot(identical(fixed_ids, original$gene_id))
  converted <- loader$gse205906_ens101_input_loader_v2(raw_file(path), ledger, fixed_ids,
    source_raw$g205_prefit, gene_attributes = attrs)
  matrix <- converted$inputs$count_matrix
  observed <- as.matrix(original[, -1L, with = FALSE]); storage.mode(observed) <- "integer"
  stopifnot(identical(unname(matrix), unname(observed)), identical(rownames(matrix), fixed_ids),
    identical(colnames(matrix), ledger$GEO_sample),
    identical(converted$transformation$count_values_changed, FALSE))
  saveRDS(converted, file.path(work, "prepared", paste0("G205_", name, "_inputs.rds")), version = 3L)
  restored <- readRDS(file.path(work, "prepared", paste0("G205_", name, "_inputs.rds")))
  stopifnot(identical(restored, converted))
  summaries[[paste0("G205_", name)]] <- list(genes = nrow(matrix), libraries = ncol(matrix),
    counts_unchanged = TRUE, sample_order_unchanged = TRUE, complete_annotation_join = TRUE,
    source = "Original input loader; fixed IDs recovered from the hash-authenticated complete matrix, not a recreated historical run manifest",
    role = if (name == "HISAT2_s0") "ORIGINAL" else if (name == "STAR_s0") "MATCHED_PRIMARY" else "STRAND_DIAGNOSTIC")
  rm(converted, original, observed, matrix, restored); invisible(gc())
}
stopifnot(identical(RNGkind(), old_rng), identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), had_seed),
  !had_seed || identical(get(".Random.seed", envir = .GlobalEnv), old_seed))
report <- list(status = "PASS_ORIGINAL_INPUT_READERS_WITHOUT_FITTING", branches = summaries,
  annotation_qa = annotation_qa,
  R_version = R.version.string, platform = R.version$platform, package_versions = as.list(required),
  scientific_model_runs = 0L, random_draws = 0L, original_files_modified = FALSE,
  historical_manifest_recreated = FALSE, biological_independence_authenticated = FALSE,
  complete_analysis_runtime_validated = FALSE)
jsonlite::write_json(report, file.path(work, "receipts/INPUT_READER_QA.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
writeLines(capture.output(utils::sessionInfo()), file.path(work, "receipts/R_SESSION.txt"))
cat("PASS_ORIGINAL_INPUT_READERS_WITHOUT_FITTING\n")
