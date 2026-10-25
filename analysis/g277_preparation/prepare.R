# Reuse only the pinned input block and pre-DGEList design prefix of the source.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
repo <- spec$repo; out <- spec$output
stopifnot(as.character(getRversion()) == "4.6.1")
for (name in names(spec$source_pins)) {
  path <- file.path(repo, name); pin <- spec$source_pins[[name]]
  stopifnot(file.info(path)$size == pin$bytes,
    identical(digest::digest(file = path, algo = "sha256", serialize = FALSE), pin$sha256))
}
suppressPackageStartupMessages({library(data.table); library(limma)})
source_path <- file.path(repo, spec$source)
parsed <- parse(source_path, keep.source = TRUE)
lhs <- vapply(parsed, function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]])) as.character(x[[2L]]) else ""
}, character(1))
api <- new.env(parent = globalenv())
api$raw_dir <- file.path(spec$project, "data/raw/GSE277064")
api$metadata_dir <- file.path(spec$project, "metadata")
api$result_root <- file.path(out, "unused_original_result_namespace")
for (name in c("assert_true", "natural_contrast_matrix", "exercise_contrast_matrix")) {
  at <- which(lhs == name); stopifnot(length(at) == 1L)
  eval(parsed[[at]], api)
}
refs <- attr(parsed, "srcref")
start <- which(lhs == "ledger")[[1L]]
end <- which(lhs == "mapping_audit")[[1L]] - 1L
stopifnot(grepl("ledger <- fread", paste(deparse(parsed[[start]]), collapse = " ")),
          grepl("Rotarod samples leaked", paste(deparse(parsed[[end]]), collapse = " ")))
for (i in seq.int(start, end)) eval(parsed[[i]], api)

# Build a separate, explicitly truncated function. No DGEList or later source
# expression is evaluated and the original source file is never changed.
definition <- parsed[[which(lhs == "run_domain_model")]][[3L]]
body_parts <- as.list(definition[[3L]])[-1L]
cut <- which(vapply(body_parts, function(x) is.call(x) &&
  identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("y_unfiltered")), logical(1)))
stopifnot(length(cut) == 1L)
prefix <- body_parts[seq_len(cut - 1L)]
stopifnot(!grepl("DGEList|estimateDisp|glmQLFit|voom\\(", paste(vapply(prefix,
  function(x) paste(deparse(x), collapse = " "), character(1)), collapse = "\n")))
definition[[3L]] <- as.call(c(list(as.name("{")), prefix,
  list(quote(list(samples = sample_meta, sample_columns = sample_columns,
                  design = design, contrasts = contrast_matrix)))))
admit_domain_design <- eval(definition, api)
domains <- list()
for (ct in c("v2a", "astrocyte", "microglia")) {
  domains[[paste0(ct, "_natural_course")]] <- admit_domain_design(ct, "natural_course")
}
for (ct in c("v2a", "astrocyte")) {
  domains[[paste0(ct, "_rotarod_boundary")]] <- admit_domain_design(ct, "rotarod_boundary")
}
all_uses <- unlist(lapply(domains, `[[`, "sample_columns"), use.names = FALSE)
stopifnot(length(all_uses) == 69L, length(unique(all_uses)) == 63L,
          sum(table(all_uses) == 2L) == 6L, !"microglia_rotarod_boundary" %in% names(domains))
for (name in names(domains)) {
  d <- domains[[name]]
  natural <- grepl("_natural_course$", name)
  expected <- if (natural) rbind(c(-1, 0, 0, .5, .5), c(-1, 0, 0, 1, 0),
    c(-1, 0, 0, 0, 1), c(-1, .5, .5, 0, 0), c(-1, 1, 0, 0, 0), c(-1, 0, 1, 0, 0)) else
    rbind(c(-1, 1, 0, 0), c(-1, 0, 1, 0), c(-1, 0, 0, 1))
  stopifnot(identical(unname(d$contrasts), unname(t(expected))),
    identical(as.numeric(colSums(d$design)), rep(3, ncol(d$design))),
    qr(d$design)$rank == ncol(d$design))
}
# Independently parse every supplied raw value using base R, outside the source reader.
fractional <- 0L
for (path in api$count_files) {
  tab <- utils::read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  values <- as.matrix(tab[, -1L, drop = FALSE]); storage.mode(values) <- "double"
  rownames(values) <- tab[[1L]]
  stopifnot(identical(values, api$counts_all[, colnames(values), drop = FALSE]))
  fractional <- fractional + sum(values != floor(values))
}
# Exercise original count assertions on small invented inputs, never on a model.
rejects <- list()
for (kind in c("duplicate_gene", "negative_count", "nonfinite_count")) {
  path <- file.path(out, "receipts", paste0("invented_", kind, ".tsv"))
  x <- data.frame(gene = c("g1", "g2"), sample = c(1, 2))
  if (kind == "duplicate_gene") x$gene[2L] <- "g1"
  if (kind == "negative_count") x$sample[2L] <- -1
  if (kind == "nonfinite_count") x$sample[2L] <- Inf
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE)
  answer <- tryCatch({api$read_count_table(path); NULL}, error = identity)
  stopifnot(inherits(answer, "error")); rejects[[kind]] <- conditionMessage(answer)
}
inputs <- list(counts = api$counts_all, ledger = api$ledger, domains = domains,
  source_sha256 = digest::digest(file = source_path, algo = "sha256", serialize = FALSE))
saveRDS(inputs, file.path(out, "prepared/GSE277064_unfitted_inputs.rds"), version = 3L, compress = "gzip")
for (name in names(domains)) {
  fwrite(domains[[name]]$samples, file.path(out, "prepared", paste0(name, "_samples.tsv")), sep = "\t", na = "NA")
}
writeLines(capture.output(sessionInfo()), file.path(out, "receipts/sessionInfo.txt"))
jsonlite::write_json(list(status = "PASS_ORIGINAL_INPUT_BLOCK_AND_PREFIT_DESIGN_PREFIX",
  count_rows = nrow(api$counts_all), unique_libraries = ncol(api$counts_all),
  natural_course_libraries = nrow(api$natural_ledger), trained_libraries = nrow(api$exercise_ledger),
  domain_models_not_fitted = length(domains), library_uses = length(all_uses), shared_D14_controls = 6L,
  fractional_count_cells = fractional, independent_base_R_count_readback = TRUE,
  original_input_block_lines = c(refs[[start]][1L], refs[[end]][3L]),
  original_prefit_expression_count = length(prefix), original_reader_rejections = rejects,
  R = as.character(getRversion()), package_versions = lapply(c("data.table", "limma", "jsonlite", "digest"),
    function(x) list(package = x, version = as.character(utils::packageVersion(x)))),
  scientific_model_fits = 0L, normalization_or_scores_recomputed = FALSE, random_allocations = 0L,
  annotation_resource_construction = FALSE, historical_execution_recreated = FALSE),
  file.path(out, "receipts/PREPARATION_QA.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
