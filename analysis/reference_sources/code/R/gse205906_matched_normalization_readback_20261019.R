# Functions only. Takes a caller-held NEW matched checkpoint, never a file path.
# Reconstructs filtering, TMM and voom's expression transform without a model
# fit, RNG call, source import, file IO or acceptance-file creation. Voom
# precision weights, dispersion, mapping and external custody are not verified.

g205t_need <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
}

g205t_equal <- function(observed, expected, label) {
  g205t_need(is.numeric(observed) && is.numeric(expected) &&
    length(observed) == length(expected), paste0(label, ": numeric domain"))
  g205t_need(all(is.finite(observed)) && all(is.finite(expected)),
    paste0(label, ": finite values required"))
  g205t_need(all(abs(observed - expected) <= 1e-12 + 1e-10 * abs(expected)),
    paste0(label, ": values differ"))
  invisible(TRUE)
}

# A separate native normalization call; no reuse of the producer's function,
# saved factors or fitted statistical objects in deriving expected values.
g205t_readback <- function(checkpoint) {
  g205t_need(is.list(checkpoint) &&
    identical(checkpoint$schema, "GSE205906_MATCHED_SAVED_NATIVE_CHECKPOINT/v1") &&
    isTRUE(checkpoint$statistical_fit_executed), "New fitted checkpoint required")
  packet <- checkpoint$prepared_packet
  result <- checkpoint$result
  g205t_need(identical(result$schema, "gse205906-matched-downstream-result-v1") &&
    identical(result$binding, packet$binding) &&
    identical(result$prepared_sha256, packet$sha256) &&
    identical(result$candidate_admission, FALSE), "Prepared/result association differs")
  g205t_need(identical(as.character(getRversion()), result$runtime$R),
    "Saved R version required")
  for (package in c("edgeR", "limma")) {
    g205t_need(requireNamespace(package, quietly = TRUE) &&
      identical(as.character(utils::packageVersion(package)), result$runtime$packages[[package]]),
      "Saved normalization package version required")
  }

  counts <- packet$counts
  samples <- paste0("GSM", 6235292:6235300)
  g205t_need(is.matrix(counts) && is.numeric(counts) &&
    identical(dim(counts), c(32883L, 9L)) &&
    identical(colnames(counts), samples) &&
    is.character(rownames(counts)) && !anyNA(rownames(counts)) &&
    all(nzchar(rownames(counts))) && !anyDuplicated(rownames(counts)) &&
    all(is.finite(counts)) && all(counts >= 0 & counts == floor(counts)),
    "Full count/gene/sample domain differs")
  design <- packet$design$design
  expected_design <- diag(3L)[rep(seq_len(3L), each = 3L), , drop = FALSE]
  g205t_need(is.matrix(design) && is.numeric(design) &&
    identical(dim(design), c(9L, 3L)) && identical(rownames(design), samples) &&
    identical(colnames(design), c("sham", "SNI", "PLM")) &&
    identical(as.numeric(design), as.numeric(expected_design)),
    "Frozen no-intercept design differs")

  full <- edgeR::DGEList(counts = counts)
  keep_filter <- edgeR::filterByExpr(full, design = design)
  keep_cpm <- rowSums(edgeR::cpm(full, log = FALSE) >= 1) >= 3L
  keep <- keep_filter & keep_cpm
  base <- result$native$base_fit
  g205t_need(identical(base$keep_filter_by_expr, keep_filter) &&
    identical(base$keep_cpm, keep_cpm) && identical(base$keep, keep) &&
    sum(keep) >= 2L, "Own full-input filtering differs")
  retained <- counts[keep, , drop = FALSE]
  expected_y <- edgeR::DGEList(counts = retained)
  expected_y <- edgeR::normLibSizes(expected_y, method = "TMM")
  library_sizes <- unname(colSums(retained))
  factors <- unname(expected_y$samples$norm.factors)
  effective_sizes <- library_sizes * factors
  g205t_need(all(library_sizes > 0) && all(factors > 0) &&
    all(is.finite(effective_sizes)), "Positive normalization domain required")
  g205t_need(identical(base$y$counts, retained) &&
    identical(rownames(base$y$samples), samples), "Saved DGE count/sample keys differ")
  g205t_need(is.null(base$y$offset) && is.null(base$y$offset.prior),
    "Unexpected DGE offset in fixed no-offset pipeline")
  g205t_equal(base$y$samples$lib.size, library_sizes, "Reset library sizes")
  g205t_equal(base$y$samples$norm.factors, factors, "Independently reproduced TMM")
  g205t_equal(result$diagnostics$own_library_sizes, library_sizes, "Diagnostic library sizes")
  g205t_equal(result$diagnostics$own_TMM, factors, "Diagnostic TMM")

  voom <- base$voom
  g205t_need(identical(voom$design, design) &&
    is.matrix(voom$E) && is.numeric(voom$E) &&
    identical(dim(voom$E), dim(retained)) &&
    identical(dimnames(voom$E), dimnames(retained)) &&
    identical(rownames(voom$targets), samples) &&
    is.null(voom$other$offset.prior), "Saved voom domain/design/offset differs")
  g205t_equal(voom$targets$lib.size, effective_sizes, "Voom effective library sizes")
  g205t_equal(voom$targets$norm.factors, factors, "Voom target TMM")
  # Algebra from the pinned limma::voom implementation with default
  # normalize.method='none', no offsets, and 0.5 count/1 library priors.
  denominator <- matrix(effective_sizes + 1, nrow(retained), ncol(retained), byrow = TRUE)
  expected_expression <- log2((retained + 0.5) / denominator * 1e6)
  g205t_equal(as.numeric(voom$E), as.numeric(expected_expression), "Voom log-expression transform")

  list(status = "HELD_NEW_CHECKPOINT_NORMALIZATION_COMPONENT_VERIFIED",
    retained_genes = nrow(retained), libraries = ncol(retained),
    relative_tolerance = 1e-10, absolute_tolerance = 1e-12,
    own_filtering_and_reset_library_sizes_verified = TRUE,
    native_TMM_reproduction_verified = TRUE,
    voom_expression_transform_and_targets_verified = TRUE,
    voom_precision_weights_and_dispersion_verified = FALSE,
    model_fits = 0L, candidate_admission = FALSE,
    full_independent_numeric_acceptance = FALSE)
}
