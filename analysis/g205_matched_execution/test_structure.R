# Synthetic preparation and AST equivalence only; no study inputs or model fits.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
scratch <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
read_raw <- function(p) readBin(p, "raw", n = file.info(p)$size)
contract <- jsonlite::fromJSON(file.path(repo, "analysis/g205_matched_execution/CONTRACT.json"), simplifyVector = FALSE)
for (p in names(contract$sources)) {
  raw <- read_raw(file.path(repo, p)); expected <- contract$sources[[p]]
  stopifnot(length(raw) == expected$bytes, digest::digest(raw, "sha256", serialize = FALSE) == expected$sha256)
}
defs <- function(p) {
  env <- new.env(parent = baseenv())
  for (node in as.list(parse(p, keep.source = FALSE))) {
    stopifnot(is.call(node), identical(node[[1L]], as.name("<-")), is.call(node[[3L]]),
      identical(node[[3L]][[1L]], as.name("function")))
    eval(node, env)
  }
  env
}
core <- defs(file.path(repo, contract$auxiliary_roles$core))
adapter <- defs(file.path(repo, "analysis/g205_matched_execution/adapter.R"))
parse(file.path(repo, "analysis/g205_matched_execution/execute.R"), keep.source = FALSE) -> unused
# Compare the complete unchanged contiguous scientific preparation block.
block <- function(fun) {
  nodes <- as.list(body(fun))[-1L]
  lhs <- vapply(nodes, function(x) if (is.call(x) && identical(x[[1L]], as.name("<-"))) paste(deparse(x[[2L]]), collapse = "") else "", "")
  nodes[which(lhs == "projected"):which(lhs == "annotation")]
}
stopifnot(identical(block(core$g205m_prepare), block(adapter$g205p_prepare)))
fit_before <- core$g205m_function_descriptor(core$g205m_fit)
runtime_before <- core$g205m_function_descriptor(core$g205m_runtime)
had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
seed_before <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
rng <- RNGkind()
sources <- lapply(contract$engine_roles, function(p) read_raw(file.path(repo, p)))
ids <- sprintf("ENSRNOG%011d", seq_len(32883L))
samples <- paste0("GSM", 6235292:6235300)
counts <- matrix(rep(0:8, each = length(ids)), nrow = length(ids), dimnames = list(ids, samples))
converted <- list(inputs = list(count_matrix = counts,
  gene_annotation = data.frame(gene_id = ids, gene_name = paste0("SYNTHETIC", seq_along(ids)), gene_biotype = "protein_coding"),
  sample_ledger = data.frame(GEO_sample = samples), fixed_gene_ids = ids),
  transformation = list(count_values_changed = FALSE))
resources <- list(reactome = data.frame(synthetic = "a"), progeny_human = data.frame(synthetic = "a"),
  dorothea_rat = data.frame(synthetic = "a"), syngo_sensitivity = data.frame(synthetic = "a"))
orthology <- data.frame(group_method = "NCBI Ortholog", human_symbol = "SYNTHETIC_H", other_tax_id = "10116", other_symbol = "SYNTHETIC_R")
binding <- list(schema = "G205_MATCHED_PUBLIC_INPUT_BINDING/v1", prepared_input = list(test_only = TRUE),
  historical_native_custody_recreated = FALSE)
before <- serialize(list(converted, resources, orthology), NULL, version = 3L)
prepared <- adapter$g205p_prepare(core, sources, converted, resources, orthology, binding)
stopifnot(identical(prepared$packet$counts, counts), identical(rownames(prepared$packet$design$design), samples),
  identical(unname(prepared$packet$design$contrast[, "local_disease"]), c(-1, 1, 0)),
  identical(unname(prepared$packet$design$contrast[, "intervention"]), c(0, -1, 1)),
  identical(before, serialize(list(converted, resources, orthology), NULL, version = 3L)))
copy <- core$g205m_custody_working_copy(prepared)
stopifnot(identical(copy$packet$sha256, prepared$packet$sha256),
  identical(core$g205m_engine_hash(copy$engine), core$g205m_engine_hash(prepared$engine)))
rejects <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")
wrong <- converted; wrong$inputs$count_matrix <- counts[, 9:1]
stopifnot(rejects(adapter$g205p_prepare(core, sources, wrong, resources, orthology, binding)))
wrong <- converted; wrong$inputs$count_matrix[1L, 1L] <- -1L
stopifnot(rejects(adapter$g205p_prepare(core, sources, wrong, resources, orthology, binding)))
bad_orthology <- rbind(orthology, transform(orthology, other_symbol = "ANOTHER_R"))
stopifnot(rejects(adapter$g205p_prepare(core, sources, converted, resources, bad_orthology, binding)))
bad_binding <- binding; bad_binding$historical_native_custody_recreated <- TRUE
stopifnot(rejects(adapter$g205p_prepare(core, sources, converted, resources, orthology, bad_binding)))
fixture_path <- file.path(scratch, "synthetic_gzip_fixture.rds")
stopifnot(!file.exists(fixture_path))
toy <- list(a = 1:3, b = NA_real_)
saveRDS(toy, fixture_path, version = 3L, compress = "gzip")
stopifnot(identical(adapter$g205p_decode_gzip_rds(read_raw(fixture_path)), toy))
fixture <- as.raw(c(31,139,8,0,0,0,0,0,2,3))
stopifnot(suppressWarnings(rejects(adapter$g205p_decode_gzip_rds(fixture))), rejects(adapter$g205p_decode_gzip_rds(as.raw(1:5))))
pure <- new.env(parent = baseenv())
core$g205m_import(read_raw(file.path(repo, contract$auxiliary_roles$resource_adapter)), "gse205906_resource_label_adapter_v1", pure)
core$g205m_import(read_raw(file.path(repo, contract$auxiliary_roles$writer)),
  c("gse205906_ens101_worker_need_v1", "gse205906_ens101_worker_tsv_v1"), pure)
core$g205m_import(read_raw(file.path(repo, contract$auxiliary_roles$caller)), c("g205run_need", "g205run_tsv_view"), pure)
labels <- pure$gse205906_resource_label_adapter_v1(list(test = data.frame(name = factor(c("A", "B"), levels = c("B", "A", "C")))))
stopifnot(identical(labels$test$name, c("A", "B")))
empty <- core$g205m_empty("fgsea")
bytes <- pure$gse205906_ens101_worker_tsv_v1(pure$g205run_tsv_view(empty))
stopifnot(is.raw(bytes), length(bytes) > 0L, nrow(empty) == 0L, identical(class(empty$leadingEdge), "AsIs"))
stopifnot(identical(fit_before, core$g205m_function_descriptor(core$g205m_fit)),
  identical(runtime_before, core$g205m_function_descriptor(core$g205m_runtime)),
  identical(rng, RNGkind()), identical(had_seed, exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)),
  !had_seed || identical(seed_before, get(".Random.seed", envir = .GlobalEnv)))
cat("PASS_SYNTHETIC_PREPARATION_AST_AND_NEGATIVE_INPUT_TESTS_NO_FITS\n")
