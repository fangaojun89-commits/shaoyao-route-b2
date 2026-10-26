# Synthetic serialization test only: no biological input or statistical fitting.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !file.exists(args[[2L]]))
repo <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)
destination <- args[[2L]]
stopifnot(dir.create(destination, recursive = TRUE))
source(file.path(repo, 'analysis/study_sources/code/R/gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R'))
source(file.path(repo, 'analysis/g205_execution/write_result.R'))
example <- data.frame(text = c('quote" tab\t unicode\u03b2', NA_character_),
  effect = c(pi, NA_real_), n = c(2L, NA_integer_), yes = c(TRUE, FALSE),
  label = factor(c('A', 'B')), stringsAsFactors = FALSE)
example$leadingEdge <- list(c('GENE1', 'GENE2'), character())
tables <- setNames(rep(list(example), 15L), gse205906_ens101_worker_tables_v1())
checked <- list(payload = list(result = list(tables = tables),
  output_sha256 = lapply(tables, digest::digest, algo = 'sha256', serializeVersion = 3L)))
plan <- gse205906_ens101_worker_write_plan_v1(checked,
  list(scope = 'SYNTHETIC_SERIALIZATION_TEST_NO_MODEL', historical_custody = FALSE))
pin <- gse205906_ens101_worker_pin_v1
g205_write_result(plan, file.path(destination, 'output'), pin)
stopifnot(length(plan$files) == 17L,
  identical(readRDS(file.path(destination, 'output/checked_payload.rds')), checked))
observed <- read.delim(file.path(destination, 'output/gene_edger.tsv'),
  quote = '"', check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = 'UTF-8')
stopifnot(identical(observed$text, example$text), identical(observed$effect, example$effect),
  identical(observed$n, example$n), identical(observed$yes, example$yes),
  identical(observed$label, as.character(example$label)),
  identical(jsonlite::fromJSON(observed$leadingEdge[[1L]]), example$leadingEdge[[1L]]),
  length(jsonlite::fromJSON(observed$leadingEdge[[2L]])) == 0L)
rejected <- function(expr) inherits(tryCatch(force(expr), error = identity), 'error')
stopifnot(rejected(g205_write_result(plan, file.path(destination, 'output'), pin)))
bad <- plan
bad$file_pins[[1L]]$sha256 <- paste(rep('0', 64L), collapse = '')
stopifnot(rejected(g205_write_result(bad, file.path(destination, 'wrong_hash'), pin)))
bad <- plan
names(bad$files)[1L] <- names(bad$file_pins)[1L] <- '../outside.tsv'
stopifnot(rejected(g205_write_result(bad, file.path(destination, 'traversal'), pin)),
  !file.exists(file.path(destination, 'traversal')))
result <- list(status = 'PASS_SYNTHETIC_WRITER', scientific_fits = 0L,
  tests = c('original_15_table_serialization_17_files', 'all_output_hash_readback',
    'checked_RDS_roundtrip', 'TSV_text_numeric_NA_factor_JSON_roundtrip',
    'existing_destination_rejected', 'incorrect_hash_rejected', 'path_traversal_rejected'),
  output_bytes = plan$total_bytes)
jsonlite::write_json(result, file.path(destination, 'RESULT.json'), pretty = TRUE, auto_unbox = TRUE)
print(result$status)
