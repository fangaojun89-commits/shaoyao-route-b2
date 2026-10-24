# Software acceptance only; no scientific input or model execution.
# Path-parameterized derivative of verify_g205_linux_core_20261023.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, Sys.info()[['sysname']] == 'Linux',
          as.character(getRversion()) == '4.6.1')
library_root <- normalizePath(args[[1]], mustWork = TRUE)
catalog_path <- normalizePath(args[[2]], mustWork = TRUE)
output <- args[[3]]
stopifnot(!file.exists(output), dir.exists(dirname(output)))
.libPaths(c(library_root, .Library), include.site = FALSE)
stopifnot(all(normalizePath(.libPaths()) %in%
              normalizePath(c(library_root, .Library))))
stopifnot(requireNamespace('jsonlite', quietly = TRUE))
stopifnot(identical(normalizePath(getNamespaceInfo(getNamespace('jsonlite'), 'path')),
                    normalizePath(file.path(library_root, 'jsonlite'))))
catalog <- jsonlite::fromJSON(catalog_path, simplifyVector = FALSE)
stopifnot(length(catalog$packages) == 70L)
versions <- lapply(catalog$packages, function(row) {
  path <- file.path(library_root, row$package, 'DESCRIPTION')
  d <- read.dcf(path)
  stopifnot(d[1, 'Package'] == row$package, d[1, 'Version'] == row$version)
  list(package = row$package, version = unname(d[1, 'Version']))
})
roots <- c('edgeR', 'limma', 'DESeq2', 'fgsea', 'decoupleR',
           'data.table', 'jsonlite', 'digest')
for (package in roots) stopifnot(requireNamespace(package, quietly = TRUE))
rebuilt <- vapply(versions, `[[`, character(1), 'package')
stopifnot(length(unique(rebuilt)) == 70L)
namespaces <- lapply(sort(loadedNamespaces()), function(package) {
  path <- if (package == 'base') file.path(.Library, 'base') else
    getNamespaceInfo(getNamespace(package), 'path')
  path <- normalizePath(path, mustWork = TRUE)
  if (package %in% rebuilt) {
    stopifnot(identical(path, normalizePath(file.path(library_root, package))))
  } else stopifnot(startsWith(path, paste0(normalizePath(.Library), '/')))
  list(package = package, path = path)
})
stopifnot(dir.create(output))
jsonlite::write_json(list(status = '70_INSTALLED_VERSIONS_AND_EIGHT_ROOT_LOADS_PASS',
  versions = versions, namespaces = namespaces, requested_roots = roots,
  R_version = as.character(getRversion()), scientific_execution = FALSE,
  existing_R_runtime_reused = TRUE, whole_study_cold_restore = FALSE),
  file.path(output, 'RESULT.json'), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output, 'sessionInfo.txt'))
cat('PASS: installed package versions and namespace origins\n')
