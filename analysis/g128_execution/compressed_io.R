# Operational prerequisite only. No counts, resources, models or inference.
g128_check_compressed_io <- function(receipts) {
  stopifnot(dir.exists(receipts))
  expected <- c(data.table = "1.18.4", R.methodsS3 = "1.8.2",
                R.oo = "1.27.1", R.utils = "2.13.0")
  result <- list(schema = "g128-compressed-table-io-v1",
                 status = "PREFIT_COMPRESSED_IO_CHECK_STARTED",
                 stage = "before_scientific_execution",
                 actual_locale = Sys.getlocale(),
                 native_UTF8 = isTRUE(l10n_info()[["UTF-8"]]),
                 scientific_model_invocations = 0L,
                 required_versions = as.list(expected), loaded = list())
  record <- file.path(receipts, "COMPRESSED_IO_CHECK.json")
  tryCatch({
    if (!isTRUE(result$native_UTF8))
      stop("UTF8_CTYPE_REQUIRED", call. = FALSE)
    for (name in names(expected)) {
      if (!requireNamespace(name, quietly = TRUE))
        stop(paste0("MISSING_PACKAGE:", name), call. = FALSE)
      actual <- as.character(utils::packageVersion(name))
      result$loaded[[name]] <- list(version = actual,
        namespace_path = getNamespaceInfo(asNamespace(name), "path"))
      if (!identical(actual, expected[[name]]))
        stop(paste0("PACKAGE_VERSION:", name), call. = FALSE)
    }
    path <- file.path(receipts, "GZIP_IO_PROBE.tsv.gz")
    if (file.exists(path)) stop("EXISTING_GZIP_PROBE", call. = FALSE)
    expected_table <- data.table::data.table(
      id = 1:3, value = c(1, NA_real_, -2.5), text = c("first", "second", "third_\u0394"))
    data.table::fwrite(expected_table, path, sep = "\t", quote = TRUE, na = "NA")
    magic <- readBin(path, "raw", n = 2L)
    if (!identical(magic, as.raw(c(0x1f, 0x8b))))
      stop("GZIP_MAGIC", call. = FALSE)
    observed <- data.table::fread(path, sep = "\t", na.strings = "NA")
    if (!identical(as.data.frame(observed), as.data.frame(expected_table)))
      stop("GZIP_ROUNDTRIP_CONTENT", call. = FALSE)
    result$status <- "PASS_EXACT_GZIP_WRITE_READ_PREFIT_CHECK"
    result$roundtrip_rows <- 3L
    result$probe_bytes <- unname(file.info(path)$size)
    jsonlite::write_json(result, record, auto_unbox = TRUE, pretty = TRUE)
    invisible(result)
  }, error = function(condition) {
    result$status <- "FAIL_COMPRESSED_IO_BEFORE_SCIENTIFIC_EXECUTION"
    result$condition <- conditionMessage(condition)
    jsonlite::write_json(result, record, auto_unbox = TRUE, pretty = TRUE)
    stop(paste0("G128_COMPRESSED_IO_PREFLIGHT:", conditionMessage(condition)), call. = FALSE)
  })
}
