# Read-only rejection checks against the original contracts; no model execution.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- args[[1L]]; run <- args[[2L]]
spec <- jsonlite::fromJSON(file.path(run, "receipts/RUN.json"), simplifyVector = FALSE)
read <- function(p) readBin(p, "raw", n = file.info(p)$size)
worker <- new.env(parent = baseenv())
for (x in as.list(parse(text = rawToChar(read(file.path(repo, spec$worker)))))) {
  stopifnot(is.call(x), identical(x[[1L]], as.name("<-")), is.call(x[[3L]]),
    identical(x[[3L]][[1L]], as.name("function")))
  eval(x, worker)
}
role_map <- jsonlite::fromJSON(file.path(repo, "analysis/study_sources/G205_SOURCE_ROLES.json"), simplifyVector = FALSE)$roles
bundle <- lapply(role_map, function(x) read(file.path(repo, "analysis/study_sources", x$path)))
names(bundle) <- vapply(role_map, `[[`, "", "role")
full <- read(file.path(repo, spec$full_config)); base <- read(file.path(repo, spec$baseline_config))
api <- worker$gse205906_ens101_worker_sources_v1(bundle, full, base)
cfg <- api$config$gse205906_ens101_prefit_view_v1(
  api$config$gse205906_ens101_check_full_model_config_v1(full, base), base)
inputs <- readRDS(spec$prepared)$inputs
prefit <- readRDS(file.path(run, "derived/prefit.rds"))
expect_error <- function(f, pattern) {
  value <- tryCatch({f(); NULL}, error = function(e) conditionMessage(e))
  stopifnot(is.character(value), grepl(pattern, value, fixed = TRUE)); TRUE
}
bad <- bundle[c("model", "helpers", "annotation_adapter")]
bad$model[[1L]] <- as.raw(bitwXor(as.integer(bad$model[[1L]]), 1L))
source_rejected <- expect_error(function() api$prefit$gse205906_ens101_prefit_sources_v1(bad, cfg), "ENS101_PREFIT_SOURCE")
ledger <- inputs$sample_ledger
ledger$label[[1L]] <- "WRONG_LIBRARY"
ledger_rejected <- expect_error(function() api$prefit$gse205906_ens101_prefit_ledger_v1(ledger), "ENS101_PREFIT")
wrong <- inputs
wrong$count_matrix <- wrong$count_matrix[nrow(wrong$count_matrix):1L, , drop = FALSE]
gene_order_rejected <- expect_error(function() api$prefit$gse205906_ens101_effect_independent_prefit_v1(
  wrong, cfg, bundle[c("model", "helpers", "annotation_adapter")], prefit$declared_input_file_pins),
  "ENS101_PREFIT_NEW_32883_BY_NINE_GSM_INPUT")
result <- list(status = "PASS", altered_source_rejected = source_rejected,
  wrong_library_rejected = ledger_rejected, altered_gene_order_rejected = gene_order_rejected,
  mutated_files = 0L, fits = 0L)
cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), "\n")
