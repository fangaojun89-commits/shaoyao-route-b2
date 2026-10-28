# Resource content checks only; no counts, normalization, fitting or RNG calls.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- args[[1L]]; work <- args[[2L]]
data.table::setDTthreads(1L)
manifest <- jsonlite::fromJSON(file.path(repo, "MANIFEST.json"), simplifyVector = FALSE)
capture <- function(relative) {
  path <- file.path(repo, relative); raw <- readBin(path, "raw", n = file.info(path)$size)
  pin <- manifest$files[[relative]]
  stopifnot(length(raw) == pin$bytes, identical(digest::digest(raw, "sha256", serialize = FALSE), pin$sha256))
  raw
}
definitions <- function(raw, wanted) {
  nodes <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted,
    as.list(parse(text = rawToChar(raw), keep.source = FALSE)))
  stopifnot(setequal(vapply(nodes, function(x) as.character(x[[2L]]), ""), wanted), length(nodes) == length(wanted))
  target <- new.env(parent = baseenv())
  target$setNames <- stats::setNames
  for (node in nodes) {
    stopifnot(is.call(node[[3L]]), identical(node[[3L]][[1L]], as.name("function")))
    eval(node, target)
  }
  lockEnvironment(target, bindings = TRUE); target
}
prefix <- "analysis/study_sources/code/R/"
kernel <- definitions(capture(paste0(prefix, "gse128623_model_kernel_loader_candidate_v1_20261007.R")),
  "g128_kernel_loader")$g128_kernel_loader(
    capture(paste0(prefix, "analyze_gse128623_exact30_linux_model_v2.R")),
    capture(paste0(prefix, "validate_gse128623_exact30_linux_model_v2.R")))
guard <- definitions(capture(paste0(prefix, "gse128623_model_member_guard_candidate_v2_20261010.R")), "guard")$guard
api <- definitions(capture(paste0(prefix, "gse205906_resource_object_validation_v1.R")),
  c("gse205906_validate_resource_tables_v1", "gse205906_resource_label_adapter_v1"))
resource_path <- file.path(work, "data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds")
stopifnot(identical(digest::digest(file = resource_path, algo = "sha256"),
  "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51"))
resources <- readRDS(resource_path)
before <- serialize(resources, NULL, version = 3L)
guard_result <- guard(data.table::copy(resources))
stopifnot(identical(guard_result$status, "MEMBER_CONTENT_ONLY_PASS"))
a <- kernel$analyzer$v2a_program_members(data.table::copy(resources))
v <- kernel$validator$v2v_members(data.table::copy(resources))
stopifnot(identical(as.data.frame(a), as.data.frame(v)))
counts <- a[, .(member_rows = .N, programs = data.table::uniqueN(program_id)), by = collection]
files <- c(reactome = "reactome97_human_long.tsv", progeny_human = "progeny_human_top500.tsv",
  progeny_mouse = "progeny_mouse_top500.tsv", dorothea_human = "dorothea_human_AB.tsv",
  dorothea_mouse = "dorothea_mouse_AB.tsv", dorothea_rat = "dorothea_rat_AB_one_to_one.tsv",
  syngo_sensitivity = "quickgo_syngo_mammalian_direct_human.tsv", matrisome_human = "matrisome_human.tsv",
  matrisome_mouse = "matrisome_mouse.tsv", ortholog_one_to_one = "ncbi_human_mouse_rat_one_to_one.tsv")
exports <- lapply(files, function(name) data.table::fread(file.path(work, "data/processed/gene_sets", name),
  sep = "\t", na.strings = c("NA", ""), showProgress = FALSE))
exports$syngo_sensitivity <- unique(exports$syngo_sensitivity[, c("term_id", "term_name", "gene_symbol"), with = FALSE])
g205 <- api$gse205906_validate_resource_tables_v1(resources, exports)
adapted <- api$gse205906_resource_label_adapter_v1(resources)
adaptation <- api$gse205906_validate_resource_tables_v1(resources, adapted)
stopifnot(identical(before, serialize(resources, NULL, version = 3L)),
  identical(g205$status, "PASS_FROZEN_RESOURCE_CONTENT_EXPORT_AGREEMENT_ONLY"),
  identical(adaptation$status, g205$status))
jsonlite::write_json(list(status = "PASS_ORIGINAL_RESOURCE_CONSUMERS_WITHOUT_MODELS",
  G128_analyzer_validator_members_agree = TRUE, G128_member_inventory = as.data.frame(counts),
  G205_all_ten_table_export_validation = g205, G205_label_adapter_validation = adaptation$status,
  resource_object_unchanged = TRUE, counts_read = FALSE, model_fits = 0L,
  interpretation = "Full frozen resource domain before expression filtering; not eligible tested programs or independent biological validation"),
  file.path(work, "receipts/RESOURCE_CONSUMER_QA.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
