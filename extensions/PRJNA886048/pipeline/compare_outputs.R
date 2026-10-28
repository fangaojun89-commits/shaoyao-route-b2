# Compare regenerated exports against all frozen numeric/text table content.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
reference <- file.path(args[1], "data/output")
produced <- file.path(args[2], "data/output")
files <- list.files(reference, pattern = "\\.tsv$", recursive = TRUE)
stopifnot(length(files) > 10L)
for (name in files) {
  x <- read.delim(file.path(reference, name), check.names = FALSE, stringsAsFactors = FALSE)
  y <- read.delim(file.path(produced, name), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(names(x), names(y)), identical(dim(x), dim(y)))
  for (column in names(x)) {
    if (is.numeric(x[[column]]) && is.numeric(y[[column]])) {
      stopifnot(isTRUE(all.equal(x[[column]], y[[column]], tolerance = 1e-9, check.attributes = FALSE)))
    } else stopifnot(identical(x[[column]], y[[column]]))
  }
}
cat("PASS_FROZEN_TABLE_COMPARISON", length(files), "tables; tolerance 1e-9; no fitted models\n")
