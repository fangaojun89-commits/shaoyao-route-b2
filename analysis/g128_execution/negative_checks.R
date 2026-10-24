# Checks original failure detection against a completed synthetic fixture.
# Does not refit a model and never edits the original synthetic outputs.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
run <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
qa <- jsonlite::fromJSON(file.path(run, "receipts/R_VALIDATION.json"))
stopifnot(qa$status == "PASS_FULL_SYNTHETIC_MODEL_AND_INDEPENDENT_RECONSTRUCTION")
contract <- jsonlite::fromJSON(file.path(repo, "analysis/g128_execution/CONTRACT.json"), simplifyVector = FALSE)
raws <- setNames(lapply(contract$sources, function(r) {
  path <- file.path(repo, r$path)
  value <- readBin(path, "raw", n = file.info(path)$size)
  stopifnot(length(value) == r$identity$bytes,
            tolower(as.character(openssl::sha256(value))) == r$identity$sha256)
  value
}), vapply(contract$sources, `[[`, "", "role"))
import <- function(role, name) {
  env <- new.env(parent = baseenv())
  eval(parse(text = rawToChar(raws[[role]]), keep.source = FALSE), env)
  get(name, env, inherits = FALSE)
}
loader <- import("loader", "g128_kernel_loader")
rejected <- function(expression) {
  error <- tryCatch({ force(expression); NULL }, error = function(e) conditionMessage(e))
  stopifnot(is.character(error), length(error) == 1L)
  error
}
changed <- raws$analyzer
changed[[1]] <- as.raw(bitwXor(as.integer(changed[[1]]), 1L))
source_rejection <- rejected(loader(changed, raws$validator))
a <- loader(raws$analyzer, raws$validator)$analyzer
root <- file.path(run, "negative_checks")
stopifnot(!file.exists(root), dir.create(root))
fixture <- file.path(run, "synthetic_inputs")
counts <- data.table::fread(file.path(fixture, "counts.tsv"))
reordered <- counts[, c(1L, 3L, 2L, 4:ncol(counts)), with = FALSE]
bad_path <- file.path(root, "wrong_sample_order.tsv")
data.table::fwrite(reordered, bad_path, sep = "\t")
order_rejection <- rejected(a$v2a_read_inputs(bad_path, file.path(fixture, "annotation.tsv"),
  file.path(fixture, "ledger.tsv"), file.path(fixture, "frozen_design.tsv"), "gene_id", "GeneID", "gene"))
stopifnot(grepl("COUNT_COLUMNS_NOT_EXACT_LEDGER_ORDER", order_rejection, fixed = TRUE))
resources <- readRDS(file.path(fixture, "resources.rds"))
extra <- data.table::copy(resources$progeny_mouse[1])
extra$weight <- extra$weight + 1
resources$progeny_mouse <- rbind(resources$progeny_mouse, extra)
weight_rejection <- rejected(import("member_guard", "guard")(resources))
stopifnot(grepl("WEIGHT_CONFLICT", weight_rejection, fixed = TRUE))
read <- function(path) data.table::fread(file.path(run, "model_output", path))
p <- as.data.frame(read("primary_all30/gene_level/edgeR_primary.tsv.gz"))[c("gene_id", "contrast", "logFC", "BH_q")]
q <- as.data.frame(read("complete29_sensitivity/gene_level/edgeR_primary.tsv.gz"))[c("gene_id", "contrast", "logFC", "BH_q")]
comparison <- read("comparison/all30_vs_complete29_gene.tsv.gz")
comparison$logFC_all30[[1]] <- comparison$logFC_all30[[1]] + 0.1
comparison_rejection <- rejected(import("comparison_guard", "g128_comparison_guard")(
  p, q, comparison, read("comparison/all30_vs_complete29_summary.tsv"), 360L))
stopifnot(grepl("BRANCH_COMPARISON_UNVERIFIED", comparison_rejection, fixed = TRUE))
result <- list(status = "PASS_FOUR_ORIGINAL_REJECTIONS_NO_MODEL_RERUN",
  altered_source = source_rejection, sample_order = order_rejection,
  conflicting_membership_weight = weight_rejection, altered_comparison_effect = comparison_rejection,
  original_outputs_modified = FALSE, real_data_inputs_read = FALSE)
jsonlite::write_json(result, file.path(root, "VALIDATION.json"), auto_unbox = TRUE, pretty = TRUE)
cat(result$status, "\n")
