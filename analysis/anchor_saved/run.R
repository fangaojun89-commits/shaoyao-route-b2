# Isolated driver for the unchanged retained-output validator; no model source.
stopifnot(as.character(getRversion()) == "4.6.1")
required <- c(data.table = "1.18.4", jsonlite = "2.0.0", R.utils = "2.13.0")
for (package in names(required)) {
  stopifnot(as.character(utils::packageVersion(package)) == required[[package]])
}
args <- commandArgs(trailingOnly = FALSE)
entry <- normalizePath(sub("^--file=", "", args[grep("^--file=", args)]), mustWork = TRUE)
project <- normalizePath(file.path(dirname(entry), "..", ".."), mustWork = TRUE)
receipts <- file.path(dirname(project), "receipts")
jsonlite::write_json(list(R = as.character(getRversion()), platform = R.version$platform,
  packages = as.list(required), locale = Sys.getlocale(), library_paths = .libPaths()),
  file.path(receipts, "runtime.json"), auto_unbox = TRUE, pretty = TRUE)
sys.source(file.path(dirname(entry), "validate_GSE16710_outputs.R"), envir = new.env(parent = globalenv()))
writeLines(capture.output(sessionInfo()), file.path(receipts, "sessionInfo.txt"), useBytes = TRUE)
