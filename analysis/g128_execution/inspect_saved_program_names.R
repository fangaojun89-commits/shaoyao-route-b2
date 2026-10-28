# Metadata/encoding inspection only; does not calculate or accept model results.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
destination <- args[[2L]]
requested_locale <- args[[3L]]
stopifnot(dir.exists(destination), identical(Sys.getlocale("LC_CTYPE"), requested_locale))
contract <- jsonlite::fromJSON(file.path(spec$repo, "analysis/g128_execution/CONTRACT.json"),
                              simplifyVector = FALSE)
raws <- setNames(lapply(contract$sources, function(record) {
  value <- readBin(file.path(spec$repo, record$path), "raw", n = record$identity$bytes)
  stopifnot(length(value) == record$identity$bytes,
            tolower(as.character(openssl::sha256(value))) == record$identity$sha256)
  value
}), vapply(contract$sources, `[[`, "", "role"))
module <- new.env(parent = baseenv())
eval(parse(text = rawToChar(raws$loader)), module)
validator <- module$g128_kernel_loader(raws$analyzer, raws$validator)$validator
data.table::setDTthreads(1L)
resources <- readRDS(spec$paths$gene_set_resource)
members <- validator$v2v_members(resources)
registry <- unique(members[, .(collection, program_id, program_name)])
stopifnot(!anyDuplicated(registry[, .(collection, program_id)]))
files <- c("primary_all30/loao/program_chronic_true_refit.tsv.gz",
           "complete29_sensitivity/loao/program_chronic_true_refit.tsv.gz")
records <- lapply(files, function(name) {
  observed <- data.table::fread(file.path(spec$paths$staging, name),
    select = c("collection", "program_id", "program_name"),
    colClasses = "character", na.strings = c("NA", "NaN", ""))
  observed <- unique(observed)
  stopifnot(!anyDuplicated(observed[, .(collection, program_id)]))
  expected <- merge(observed[, .(collection, program_id)], registry,
                    by = c("collection", "program_id"), all.x = TRUE, sort = FALSE)
  data.table::setorderv(observed, c("collection", "program_id"))
  data.table::setorderv(expected, c("collection", "program_id"))
  stopifnot(identical(observed$collection, expected$collection),
            identical(observed$program_id, expected$program_id),
            !anyNA(expected$program_name))
  left <- as.character(observed$program_name)
  right <- as.character(expected$program_name)
  equal <- vapply(seq_along(left), function(i) identical(left[[i]], right[[i]]), logical(1L))
  raw_equal <- vapply(seq_along(left), function(i)
    identical(charToRaw(left[[i]]), charToRaw(right[[i]])), logical(1L))
  # Tagging is a hypothetical interpretation on an in-memory copy only. It is
  # not applied to any saved table, original validator or acceptance decision.
  tagged <- left
  Encoding(tagged) <- "UTF-8"
  tagged_equal <- vapply(seq_along(left), function(i)
    identical(tagged[[i]], right[[i]]), logical(1L))
  selected <- head(which(!equal), 12L)
  list(file = name, program_count = length(left),
       exact_vector_predicate_used_by_original_validator = identical(left, right),
       exact_cell_mismatches = sum(!equal), raw_byte_mismatches = sum(!raw_equal),
       hypothetical_UTF8_tag_cell_mismatches = sum(!tagged_equal),
       observed_encoding_counts = as.list(table(Encoding(left))),
       expected_encoding_counts = as.list(table(Encoding(right))),
       first_mismatches = lapply(selected, function(i) list(
         collection = observed$collection[[i]], program_id = observed$program_id[[i]],
         observed_bytes_hex = paste(as.character(charToRaw(left[[i]])), collapse = ""),
         expected_bytes_hex = paste(as.character(charToRaw(right[[i]])), collapse = ""),
         observed_encoding = Encoding(left[[i]]), expected_encoding = Encoding(right[[i]]))))
})
jsonlite::write_json(list(status = "SAVED_PROGRAM_NAME_METADATA_INSPECTED_NO_MODEL_FITS",
  requested_locale = requested_locale, actual_locale = Sys.getlocale(),
  files = records, scientific_model_invocations = 0L,
  hypothetical_tag_applied_only_to_private_memory_copy = TRUE,
  source_or_saved_tables_modified = FALSE, acceptance_allowed = FALSE,
  original_hidden_failure_cause_proven = FALSE),
  file.path(destination, "RESULT.json"), auto_unbox = TRUE, pretty = TRUE)
capture.output(sessionInfo(), file = file.path(destination, "sessionInfo.txt"))
