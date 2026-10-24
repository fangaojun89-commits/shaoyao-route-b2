args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
contract <- jsonlite::fromJSON(file.path(spec$repo, "analysis/g128_execution/CONTRACT.json"),
                               simplifyVector = FALSE)
raws <- setNames(lapply(contract$sources, function(record) {
  path <- file.path(spec$repo, record$path)
  value <- readBin(path, "raw", n = file.info(path)$size)
  stopifnot(length(value) == record$identity$bytes,
            tolower(as.character(openssl::sha256(value))) == record$identity$sha256)
  value
}), vapply(contract$sources, `[[`, "", "role"))
import <- function(role, name) {
  env <- new.env(parent = baseenv())
  eval(parse(text = rawToChar(raws[[role]]), keep.source = FALSE), env)
  get(name, envir = env, inherits = FALSE)
}
kernels <- import("loader", "g128_kernel_loader")(raws$analyzer, raws$validator)
a <- kernels$analyzer
v <- import("comparison_helper", "g128_nonmutating_validator")(kernels$validator)
driver <- import("driver_factory", "g128_copy_isolated_model_driver")(
  raws$driver, raws$comparison_helper)
member_guard <- import("member_guard", "guard")
comparison_guard <- import("comparison_guard", "g128_comparison_guard")
data.table::setDTthreads(1L)
# Serial scheduling only; statistical algorithms, seeds and parameters are unchanged.
BiocParallel::register(BiocParallel::SerialParam())
a$v2a_require_runtime(test_only = spec$mode != "run")
v$v2v_runtime(test_only = spec$mode != "run")
stopifnot(identical(a$V2A_ENRICHMENT_SEED, 20260819L),
          identical(a$V2A_RNG_KIND, c("Mersenne-Twister", "Inversion", "Rejection")))
capture.output(sessionInfo(), file = file.path(spec$output, "receipts/sessionInfo.txt"))
read_pair <- function(paths) {
  ai <- a$v2a_read_inputs(paths$count_matrix, paths$feature_annotation, paths$sample_ledger,
                          paths$frozen_design, "gene_id", "GeneID", "gene")
  vi <- v$v2v_read_input(paths$count_matrix, paths$feature_annotation, paths$sample_ledger,
                         "gene_id", "GeneID", "gene")
  stopifnot(identical(ai$counts, vi$counts), identical(ai$symbols, vi$symbols),
            identical(as.data.frame(ai$ledger), as.data.frame(vi$ledger)),
            identical(as.data.frame(ai$annotation), as.data.frame(vi$annotation)))
  list(analyzer = ai, validator = vi)
}
if (spec$mode == "synthetic") {
  source(file.path(spec$repo, "analysis/g128_execution/synthetic.R"), local = TRUE)
  result <- run_synthetic(spec, a, v, member_guard, comparison_guard, read_pair)
} else {
  paths <- spec$paths[c("count_matrix", "feature_annotation", "sample_ledger", "frozen_design",
                        "gene_set_resource", "s1_zero", "staging", "log_directory")]
  if (spec$mode == "run") {
    result <- driver(raws$analyzer, raws$validator, raws$loader, raws$member_guard,
                     raws$comparison_guard, paths)
    result$portable_entry_scope <- "Original copy-isolated model sequence; not historical native custody or whole-study restoration"
  } else {
    inputs <- read_pair(paths)
    a$v2a_read_s1_zero(paths$s1_zero)
    v$v2v_s1_zero(paths$s1_zero)
    resources <- readRDS(paths$gene_set_resource)
    stopifnot(identical(member_guard(resources)$status, "MEMBER_CONTENT_ONLY_PASS"))
    am <- a$v2a_program_members(data.table::copy(resources))
    vm <- v$v2v_members(data.table::copy(resources))
    stopifnot(identical(as.data.frame(am), as.data.frame(vm)))
    result <- list(status = "PASS_ORIGINAL_MODEL_INPUT_RUNTIME_CHECK_NO_FITS",
                   genes = nrow(inputs$analyzer$counts), libraries = ncol(inputs$analyzer$counts),
                   pre_expression_filter_programs = nrow(unique(am[, .(collection, program_id)])),
                   scientific_model_runs = 0L, driver_factory_constructed = TRUE,
                   linux_execution_validated = FALSE)
  }
}
result$mode <- spec$mode
result$rotation_count <- 9999L
result$minimum_members <- 10L
result$seed <- 20260819L
result$system <- Sys.info()[["sysname"]]
result$whole_study_restoration_verified <- FALSE
namespaces <- loadedNamespaces()
data.table::fwrite(data.table::data.table(
  package = namespaces,
  version = vapply(namespaces, function(x) as.character(utils::packageVersion(x)), ""),
  path = vapply(namespaces, function(x) if (x == "base") system.file(package = "base") else
    getNamespaceInfo(asNamespace(x), "path"), "")),
  file.path(spec$output, "receipts/loaded_packages.tsv"), sep = "\t")
jsonlite::write_json(result, file.path(spec$output, "receipts/R_VALIDATION.json"),
                     auto_unbox = TRUE, pretty = TRUE, na = "null")
capture.output(sessionInfo(), file = file.path(spec$output, "receipts/sessionInfo.txt"))
