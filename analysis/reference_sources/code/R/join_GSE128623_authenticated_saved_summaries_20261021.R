# CLI: MANIFEST_JSON MANIFEST_SHA256 NEW_OUTPUT. Saved TSV join only; no RDS.
# Manifest schema: gse128623-authenticated-saved-summary-join-manifest-v1;
# helper={path,bytes,sha256}; components=[{key,cell_id,reference_space,
# required_n,family_n,summary:{path,bytes,sha256}}]; upstream_evidence=[pins].
# Upstream acceptance is supplied/authenticated by the caller, not inferred here.
g128aj_need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)

g128aj_pin <- function(path) {
  g128aj_need(is.character(path) && length(path) == 1L && !is.na(path) &&
    (startsWith(path, "/") || grepl("^[A-Za-z]:", path)) && file.exists(path) &&
    !dir.exists(path), "ExistingAbsoluteRegularInput")
  link <- Sys.readlink(path)
  g128aj_need(is.na(link) || !nzchar(link), "NoSymlinkInput")
  list(path = path, bytes = unname(file.info(path)$size),
    sha256 = digest::digest(file = path, algo = "sha256"))
}

g128aj_check_pin <- function(record) {
  g128aj_need(is.list(record) && is.numeric(record$bytes) && length(record$bytes) == 1L &&
    is.finite(record$bytes) && record$bytes > 0 && record$bytes == floor(record$bytes) &&
    is.character(record$sha256) && length(record$sha256) == 1L &&
    !is.na(record$sha256) && grepl("^[a-f0-9]{64}$", record$sha256), "InputPinSchema")
  actual <- g128aj_pin(record$path)
  g128aj_need(identical(actual$sha256, record$sha256) &&
    actual$bytes == record$bytes, "InputBytesAndSHA256")
  actual
}

g128aj_specification <- function() {
  cells <- c("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
  spaces <- c("six_group", paste0(c("D4", "D7", "D14", "D28", "D42"), "_vs_control_exact"))
  x <- do.call(rbind, lapply(cells, function(cell) data.frame(cell_id = cell,
    reference_space = spaces, stringsAsFactors = FALSE)))
  x$key <- paste(x$cell_id, x$reference_space, sep = "/")
  x$required_n <- ifelse(x$reference_space == "six_group", 10001L,
    ifelse(grepl("predefined29$", x$cell_id) & x$reference_space == "D4_vs_control_exact", 126L, 252L))
  x$family_n <- ifelse(x$reference_space == "six_group", 67L, 9L)
  x
}

g128aj_main <- function(argv = commandArgs(trailingOnly = TRUE), source_path = NULL) {
  g128aj_need(length(argv) == 3L, "MANIFEST_JSON_SHA256_NEW_OUTPUTRequired")
  g128aj_need(requireNamespace("digest", quietly = TRUE) &&
    requireNamespace("jsonlite", quietly = TRUE), "ExistingDigestAndJsonliteRequiredNoInstall")
  g128aj_need(identical(as.character(getRversion()), "4.6.1"), "FrozenR461")
  if (is.null(source_path)) {
    files <- grep("^--file=", commandArgs(), value = TRUE)
    g128aj_need(length(files) == 1L, "ExplicitCLIEntrySource")
    source_path <- normalizePath(sub("^--file=", "", files), winslash = "/", mustWork = TRUE)
  }
  own <- g128aj_pin(source_path)
  manifest_pin <- g128aj_pin(argv[1L])
  g128aj_need(identical(manifest_pin$sha256, argv[2L]), "FrozenManifestSHA256")
  output <- argv[3L]
  g128aj_need((startsWith(output, "/") || grepl("^[A-Za-z]:", output)) &&
    !file.exists(output) && !dir.exists(output) && dir.exists(dirname(output)),
    "ExclusiveNewOutputExistingParent")
  manifest <- jsonlite::fromJSON(argv[1L], simplifyVector = FALSE)
  g128aj_need(identical(manifest$schema,
    "gse128623-authenticated-saved-summary-join-manifest-v1"), "FixedJoinManifestSchema")
  helper <- g128aj_check_pin(manifest$helper)
  g128aj_need(identical(helper$sha256,
    "0ef13e252d5a7c6859431a4dba4c0e26c3a665880df1a5e771c4bbe7187ed404"), "FrozenPureJoinHelper")
  components <- manifest$components
  g128aj_need(is.list(components) && length(components) == 24L &&
    is.list(manifest$upstream_evidence) && length(manifest$upstream_evidence) > 0L,
    "Exactly24ComponentsAndSuppliedUpstreamEvidence")
  text <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  count <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x == floor(x)
  g128aj_need(all(vapply(components, function(x) is.list(x) &&
    all(vapply(x[c("key", "cell_id", "reference_space")], text, logical(1L))) &&
    count(x$required_n) && count(x$family_n), logical(1L))), "FixedComponentFields")
  keys <- vapply(components, `[[`, character(1L), "key")
  expected <- g128aj_specification()
  g128aj_need(!anyDuplicated(keys) && setequal(keys, expected$key), "ExactFrozen24Keys")
  components <- components[match(expected$key, keys)]
  for (i in seq_along(components)) {
    x <- components[[i]]
    g128aj_need(identical(x$cell_id, expected$cell_id[i]) &&
      identical(x$reference_space, expected$reference_space[i]) &&
      x$required_n == expected$required_n[i] && x$family_n == expected$family_n[i],
      "FrozenCellSpaceAllocationAndFamilyDomain")
  }
  summary_records <- lapply(components, `[[`, "summary")
  summary_paths <- vapply(summary_records, `[[`, character(1L), "path")
  g128aj_need(!anyNA(summary_paths) && all(grepl("[.]tsv$", summary_paths)), "SavedSummaryTSVPaths")
  unique_paths <- unique(summary_paths)
  summaries <- lapply(unique_paths, function(path) {
    records <- summary_records[summary_paths == path]
    first <- g128aj_check_pin(records[[1L]])
    g128aj_need(all(vapply(records, function(x) is.numeric(x$bytes) && length(x$bytes) == 1L &&
      !is.na(x$bytes) && x$bytes == first$bytes && identical(x$sha256, first$sha256), logical(1L))),
      "SamePathMustHaveSameBytesAndSHA256")
    first
  })
  evidence <- lapply(manifest$upstream_evidence, g128aj_check_pin)
  all_inputs <- c(list(manifest_pin, helper), summaries, evidence)
  e <- new.env(parent = baseenv())
  sys.source(helper$path, envir = e)
  tables <- lapply(summaries, function(p) utils::read.delim(p$path,
    header = TRUE, sep = "\t", quote = "\"", comment.char = "", row.names = NULL,
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA"))
  parts <- stats::setNames(vector("list", 24L), expected$key)
  for (i in seq_along(tables)) {
    x <- tables[[i]]
    owners <- which(summary_paths == unique_paths[i])
    g128aj_need(!anyDuplicated(names(x)) && all(c("cell_id", "reference_space") %in% names(x)) &&
      !anyNA(x$cell_id) && !anyNA(x$reference_space), "UniqueColumnsAndPresentSavedJobKeys")
    keys_in_file <- paste(x$cell_id, x$reference_space, sep = "/")
    g128aj_need(setequal(unique(keys_in_file), expected$key[owners]) &&
      nrow(x) == sum(expected$family_n[owners]), "SharedTableExactAssignedJobDomainAndRows")
    for (j in owners) parts[[expected$key[j]]] <- x[keys_in_file == expected$key[j], , drop = FALSE]
  }
  answer <- e$g128join_summaries(parts, expected)
  recheck <- function() {
    invisible(lapply(c(list(own), all_inputs), g128aj_check_pin))
  }
  recheck()
  g128aj_need(dir.create(output, recursive = FALSE), "CreateExclusiveOutput")
  joined <- file.path(output, "JOINED_FAMILY_SUMMARY.tsv")
  utils::write.table(answer, joined, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  product <- g128aj_pin(normalizePath(joined, winslash = "/", mustWork = TRUE))
  recheck()
  receipt <- list(schema = "gse128623-authenticated-saved-summary-join-v1", scope = "JOIN_ONLY",
    status = "JOINED_SAVED_SUMMARIES_NOT_SCIENTIFIC_ACCEPTANCE", utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source = own, manifest = manifest_pin, helper = helper, components = components,
    authenticated_unique_summary_files = summaries, authenticated_upstream_evidence = evidence, joined_output = product,
    component_n = 24L, parsed_summary_file_n = length(tables), joined_family_n = nrow(answer), original_CLI_success = NULL,
    original_CLI_success_asserted = FALSE, scientific_acceptance = FALSE,
    upstream_scientific_closure_performed_here = FALSE, model_calls = 0L, RNG_calls = 0L,
    P_value_or_BH_recalculations = 0L, archive_or_scientific_RDS_reads = 0L,
    environment = list(R = R.version.string, digest = as.character(utils::packageVersion("digest")),
      jsonlite = as.character(utils::packageVersion("jsonlite")), sessionInfo = utils::capture.output(utils::sessionInfo())))
  jsonlite::write_json(receipt, file.path(output, "JOIN_RECEIPT.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(receipt)
}

if (sys.nframe() == 0L) g128aj_main()
