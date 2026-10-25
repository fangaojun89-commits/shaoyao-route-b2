# Load the explicit native library and parse scientific sources; do not evaluate them.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
spec <- jsonlite::fromJSON(args[[1L]], simplifyVector=FALSE)
contract <- jsonlite::fromJSON(spec$contract, simplifyVector=FALSE)
stopifnot(as.character(getRversion()) == contract$R)
versions <- character()
for (name in names(contract$packages)) {
  suppressPackageStartupMessages(library(name, character.only=TRUE))
  actual <- as.character(packageVersion(name))
  stopifnot(identical(actual, contract$packages[[name]]))
  versions[name] <- actual
}
data.table::setDTthreads(1L)
stopifnot(capabilities('cairo'), capabilities('png'))
for (rel in names(spec$staged)) {
  file <- file.path(spec$project, rel); pin <- spec$staged[[rel]]
  stopifnot(file.info(file)$size == pin$bytes,
    identical(digest::digest(file=file,algo='sha256',serialize=FALSE),pin$sha256))
}
parsed <- list()
for (rel in contract$sequence) {
  ex <- parse(file.path(spec$project, rel), keep.source=TRUE)
  parsed[[rel]] <- length(ex)
}
cfg <- yaml::read_yaml(file.path(spec$project,'config/contrasts.yaml'))
threshold <- yaml::read_yaml(file.path(spec$project,'config/thresholds.yaml'))
stopifnot(as.integer(cfg$global$random_seed) == contract$seed,
          threshold$expression_filtering$min_cpm == 1,
          Sys.getenv('GSE277064_ROAST_NROT') == '9999',
          Sys.getenv('GSE277064_RUN_FGSEA') == 'TRUE')
if (spec$mode != 'report-saved') {
  resource <- readRDS(file.path(spec$project,contract$resource_path))
  stopifnot(unique(resource$reactome$release) == '97',
    unique(resource$progeny_mouse$release) == '1.34.0',
    unique(resource$dorothea_mouse$release) == '1.23.0')
}
writeLines(capture.output(sessionInfo()),file.path(spec$receipts,'sessionInfo.txt'),useBytes=TRUE)
jsonlite::write_json(list(status='PASS_NATIVE_VERSIONS_SOURCE_PARSE_AND_STAGING',
  R=as.character(getRversion()), packages=as.list(versions), library_paths=.libPaths(),
  source_expressions_parsed=parsed, cairo=capabilities('cairo'),
  scientific_source_evaluated=FALSE, locale=Sys.getlocale()),
  file.path(spec$receipts,'PREFLIGHT.json'),pretty=TRUE,auto_unbox=TRUE)
