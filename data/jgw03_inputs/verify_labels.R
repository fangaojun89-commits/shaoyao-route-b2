# Read and validate existing labels only. No fitting or random-number generation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1]], mustWork = TRUE)
for (dataset in c('GSE227494', 'GSE289929', 'GSE299063')) {
  directory <- file.path(root, dataset)
  samples <- read.delim(file.path(directory, 'locked_samples.tsv'),
                        stringsAsFactors = FALSE, check.names = FALSE)
  registry <- readRDS(file.path(directory, 'FULL_LABEL_REGISTRY.rds'))
  expected <- as.integer(samples$ordered_dose)
  mc <- registry$MC
  companion <- registry$high_companion
  stopifnot(is.integer(mc), is.integer(companion),
            identical(dim(mc), c(10001L, 12L)),
            identical(dim(companion), c(20L, 12L)),
            identical(colnames(mc), samples$sample_id),
            identical(colnames(companion), samples$sample_id),
            identical(unname(mc[1L, ]), expected),
            identical(unname(companion[1L, ]), expected),
            registry$identity_row == 1L, registry$seed == 20260817L,
            registry$MC_nonidentity_draw_n == 10000L,
            identical(registry$rng_kind, c('Mersenne-Twister', 'Inversion', 'Rejection')),
            all(mc %in% 0:3), all(companion %in% 0:3))
  for (code in 0:3) {
    stopifnot(all(rowSums(mc == code) == 3L),
              all(rowSums(companion == code) == 3L))
  }
  stopifnot(!anyDuplicated(as.data.frame(mc)),
            !anyDuplicated(as.data.frame(companion)))
  fixed <- which(expected %in% c(1L, 2L))
  stopifnot(all(companion[, fixed, drop = FALSE] ==
              matrix(expected[fixed], nrow = 20L, ncol = length(fixed), byrow = TRUE)))
  chars <- unlist(registry[vapply(registry, is.character, logical(1))], use.names = FALSE)
  stopifnot(!any(grepl('Users[/\\\\]|autodl-tmp|seetacloud', chars)))
  cat(dataset, ': frozen 10001-row MC and 20-row companion labels PASS\n', sep = '')
}
