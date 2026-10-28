args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
contract <- jsonlite::fromJSON(file.path(args[[1L]], "analysis/resource_build/CONTRACT.json"))
stopifnot(identical(as.character(getRversion()), contract$R_version))
stopifnot(isTRUE(l10n_info()[["UTF-8"]]))
observed <- vapply(names(contract$packages), function(name) {
  stopifnot(requireNamespace(name, quietly = TRUE))
  as.character(utils::packageVersion(name))
}, character(1))
stopifnot(identical(unname(observed), unname(unlist(contract$packages))))
jsonlite::write_json(list(status = "PASS_EXISTING_BUILD_DEPENDENCY_VERSIONS",
  R_version = R.version.string, packages = as.list(observed),
  locale = Sys.getlocale(), encoding = l10n_info()),
  file.path(args[[2L]], "ENVIRONMENT.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(utils::sessionInfo()), file.path(args[[2L]], "PREFLIGHT_SESSION.txt"))
