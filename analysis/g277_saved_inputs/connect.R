# Connect retained preprocessing to original annotation/resource/input consumers.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
repo <- spec$repo; out <- spec$output
suppressPackageStartupMessages({library(data.table); library(limma); library(AnnotationDbi); library(org.Mm.eg.db)})
setDTthreads(1L)
stopifnot(as.character(getRversion()) == "4.6.1", as.character(packageVersion("org.Mm.eg.db")) == "3.23.0")
file_identity <- function(path) list(bytes = unname(file.info(path)$size),
  sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
matches_file <- function(path, expected) {
  actual <- file_identity(path)
  isTRUE(actual$bytes == expected$bytes) && identical(actual$sha256, expected$sha256)
}
object_hash <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)
for (name in names(spec$source_pins)) {
  expected <- spec$source_pins[[name]]
  stopifnot(matches_file(file.path(repo, name), expected))
}
contract <- jsonlite::fromJSON(file.path(repo, "data/g277_saved_inputs/PROVENANCE.json"), simplifyVector = FALSE)
stopifnot(matches_file(spec$prepared_input, contract$prepared_input_identity),
  matches_file(spec$resource_file, contract$resource_identity))
input <- readRDS(spec$prepared_input)
source_path <- file.path(repo, contract$original_source)
parsed <- parse(source_path, keep.source = TRUE)
lhs <- function(expressions) vapply(expressions, function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]])) as.character(x[[2L]]) else ""
}, character(1))
names_top <- lhs(parsed)
api <- new.env(parent = globalenv())
for (name in c("assert_true", "write_tsv", "make_reactome_resources")) {
  at <- which(names_top == name); stopifnot(length(at) == 1L)
  eval(parsed[[at]], api)
}
api$counts_all <- input$counts
api$result_root <- file.path(out, "annotation")
annotation_start <- which(names_top == "gene_ids")
annotation_end <- which(names_top == "gene_resources") - 1L
stopifnot(length(annotation_start) == 1L, length(annotation_end) == 1L)
for (i in seq.int(annotation_start, annotation_end)) eval(parsed[[i]], api)
annotation_file <- file.path(out, "annotation/audit/ensembl_to_mouse_annotation.tsv")
stopifnot(matches_file(annotation_file, contract$expected_annotation))
resources <- readRDS(spec$resource_file)
stopifnot(unique(resources$reactome$release) == "97", unique(resources$progeny_mouse$release) == "1.34.0",
  unique(resources$dorothea_mouse$release) == "1.23.0")

# Bind only the three original validation definitions, never fit functions.
core <- parse(file.path(repo, contract$gene_method_core))
core_names <- lhs(core)
consumer <- new.env(parent = globalenv())
for (name in c("gene_method_contract", "assert_gene_method", "validate_gene_method_inputs")) {
  at <- which(core_names == name); stopifnot(length(at) == 1L)
  eval(core[[at]], consumer)
}
details <- list(); natural_consumers <- 0L; rejected <- list()
for (key in names(contract$files)) {
  record <- contract$files[[key]]
  path <- file.path(repo, record$public_path)
  stopifnot(matches_file(path, record$identity))
  view <- readRDS(path)
  stopifnot(view$schema == contract$snapshot_schema,
    identical(view$masks$intersection, view$masks$filterByExpr & view$masks$cpm),
    identical(view$counts, input$counts[view$masks$intersection, as.character(view$sample_metadata$source_column_or_member), drop = FALSE]),
    identical(view$design, input$domains[[key]]$design), identical(view$contrasts, input$domains[[key]]$contrasts),
    identical(as.numeric(view$normalization$lib.size), as.numeric(colSums(view$counts))),
    all(is.finite(view$normalization$norm.factors)), all(view$normalization$norm.factors > 0),
    abs(sum(log(view$normalization$norm.factors))) < 1e-10,
    identical(dimnames(view$voom_checkpoint$E), dimnames(view$counts)),
    identical(dim(view$voom_checkpoint$weights), dim(view$counts)),
    all(is.finite(view$voom_checkpoint$E)), all(is.finite(view$voom_checkpoint$weights)),
    all(view$voom_checkpoint$weights > 0))
  selected_annotation <- api$annotation[match(rownames(view$counts), gene_id_versioned)]
  for (column in names(selected_annotation)) {
    stopifnot(identical(as.character(selected_annotation[[column]]), as.character(view$annotation[[column]])))
  }
  if (view$domain == "natural_course") {
    consumer$validate_gene_method_inputs(view$counts, as.character(view$sample_metadata$group),
      view$design, view$contrasts, view$voom_checkpoint)
    natural_consumers <- natural_consumers + 1L
    if (!length(rejected)) {
      for (kind in c("negative_count", "changed_design", "zero_weight")) {
        candidate <- view
        if (kind == "negative_count") candidate$counts[1L, 1L] <- -1
        if (kind == "changed_design") candidate$design[1L, 1L] <- candidate$design[1L, 1L] + 1
        if (kind == "zero_weight") candidate$voom_checkpoint$weights[1L, 1L] <- 0
        answer <- tryCatch({consumer$validate_gene_method_inputs(candidate$counts,
          as.character(candidate$sample_metadata$group), candidate$design, candidate$contrasts,
          candidate$voom_checkpoint); NULL}, error = identity)
        stopifnot(inherits(answer, "error")); rejected[[kind]] <- conditionMessage(answer)
      }
      rm(candidate)
    }
  }
  rr <- api$make_reactome_resources(api$annotation, as.data.table(resources$reactome),
    as.data.table(resources$ortholog_one_to_one), rownames(view$counts))
  stopifnot(identical(object_hash(rr$indices), record$expected$reactome_index_sha256),
    length(rr$indices) == record$expected$reactome_programs)
  # Membership coverage only: no run_ulm(), no expression collapse or scoring.
  measured_symbols <- unique(as.character(view$annotation$symbol))
  measured_symbols <- measured_symbols[!is.na(measured_symbols) & nzchar(measured_symbols)]
  pn <- as.data.table(resources$progeny_mouse)
  dn <- as.data.table(resources$dorothea_mouse)
  pc <- pn[, .(target_n = uniqueN(intersect(gene_symbol, measured_symbols))), by = pathway]
  dc <- dn[, .(target_n = uniqueN(intersect(target, measured_symbols))), by = tf]
  stopifnot(setequal(pc[target_n >= 10L, pathway], unlist(record$expected$progeny_program_ids)),
    setequal(dc[target_n >= 10L, tf], unlist(record$expected$dorothea_program_ids)))
  product <- list(native_preprocessing = view,
    reactome = list(indices = rr$indices, names = rr$names),
    resource_identity = contract$resource_identity,
    external_resource_role = "Use the exact separate frozen RDS for ULM networks; no scores or new fits here")
  output_path <- file.path(out, "prepared", paste0(key, "_native_inputs.rds"))
  saveRDS(product, output_path, version = 3L, compress = "gzip")
  details[[key]] <- list(genes = nrow(view$counts), libraries = ncol(view$counts),
    raw_count_projection_exact = TRUE, retained_mask_AND_exact = TRUE, saved_voom_E_weights_preserved = TRUE,
    annotation_fields_match_rebuilt_source = TRUE, original_Reactome_indices_match = TRUE,
    reactome_programs = length(rr$indices), progeny_programs = sum(pc$target_n >= 10L),
    dorothea_programs = sum(dc$target_n >= 10L))
  rm(view, product, rr); gc(FALSE)
}
stopifnot(natural_consumers == 3L)
writeLines(capture.output(sessionInfo()), file.path(out, "receipts/sessionInfo.txt"))
jsonlite::write_json(list(status = "PASS_ORIGINAL_ANNOTATION_RESOURCE_AND_SAVED_INPUT_CONSUMERS",
  annotation_bytes_exact = TRUE, domains = details, original_natural_input_consumers = natural_consumers,
  original_consumer_rejections = rejected, scientific_model_fits = 0L, scores_recomputed = 0L,
  filtering_or_TMM_or_voom_recomputed = FALSE, label_allocations = 0L,
  source_function_use = "Original mapIds/annotation block, make_reactome_resources, and natural gene-method validator only",
  package_versions = lapply(c("data.table", "limma", "AnnotationDbi", "org.Mm.eg.db", "jsonlite", "digest"),
    function(x) list(package = x, version = as.character(packageVersion(x)))),
  historical_runtime_authenticated = FALSE, whole_study_cold_restore = FALSE),
  file.path(out, "receipts/CONNECTION_QA.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
