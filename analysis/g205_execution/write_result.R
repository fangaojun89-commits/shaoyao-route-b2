# Ordinary directory output for the original serialization plan. No historical FD emulation.
g205_write_result <- function(plan, destination, pin) {
  stopifnot(!file.exists(destination), length(plan$files) > 0L,
    identical(names(plan$files), names(plan$file_pins)), !anyDuplicated(names(plan$files)),
    all(basename(names(plan$files)) == names(plan$files)),
    !any(names(plan$files) %in% c('', '.', '..')),
    all(vapply(plan$files, is.raw, logical(1))))
  stopifnot(dir.create(destination, recursive = FALSE))
  for (name in names(plan$files)) {
    target <- file.path(destination, name)
    stopifnot(!file.exists(target))
    writeBin(plan$files[[name]], target)
    actual <- readBin(target, 'raw', n = file.info(target)$size)
    stopifnot(identical(pin(actual), plan$file_pins[[name]]))
  }
  stopifnot(setequal(list.files(destination, all.files = TRUE, no.. = TRUE), names(plan$files)))
  invisible(TRUE)
}
