# Source-faithful chronic mroast / fgsea kernel; definitions only, no I/O or fits
# on source(). The caller supplies already-filtered counts, its native voom
# object, a matching edgeR chronic table and the frozen mapped member table.
# Native source: analyze_gse128623_exact30_linux_model_v2.R, lines 471-557.
# SHA256 e251c6ae9bb70cd944ebb1d781895878ee80963946f6e8ad94a52a83713add50

g128_unsigned_spec <- function() {
  list(
    native_source_sha256 = "e251c6ae9bb70cd944ebb1d781895878ee80963946f6e8ad94a52a83713add50",
    seed = 20260819L,
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    nrot = 9999L, minimum_members = 10L, maximum_members = 500L,
    groups = c("control", "D4", "D7", "D14", "D28", "D42"),
    chronic = c(control = -1, D4 = 0, D7 = 0, D14 = 0, D28 = .5, D42 = .5),
    collections = c("Reactome", "SynGO_assignedBy_proxy_sensitivity"),
    packages = c(data.table = "1.18.4", limma = "3.68.4", fgsea = "1.38.0",
                 openssl = "2.4.2", BiocParallel = "1.46.0"),
    ranking_rule = "max_abs_signed_sqrt_QLF_F_per_symbol_then_score_desc_symbol_asc"
  )
}

g128_unsigned_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

g128_unsigned_require_runtime <- function() {
  spec <- g128_unsigned_spec()
  g128_unsigned_assert(as.character(getRversion()) == "4.6.1", "R 4.6.1 required")
  for (package in names(spec$packages)) {
    g128_unsigned_assert(requireNamespace(package, quietly = TRUE),
                         paste("Missing existing package:", package))
    g128_unsigned_assert(as.character(utils::packageVersion(package)) == spec$packages[[package]],
                         paste("Native package version mismatch:", package))
  }
  invisible(TRUE)
}

g128_unsigned_reset_rng <- function(seed, rng_kind) {
  RNGkind(kind = rng_kind[[1L]], normal.kind = rng_kind[[2L]],
          sample.kind = rng_kind[[3L]])
  set.seed(seed)
  invisible(TRUE)
}

# The native kernel resets RNG before each collection/method. Preserve that
# internal sequence, but restore the caller's RNG kind and seed on return/error.
# This deliberately changes only the external RNG side effect, not fit inputs.
g128_unsigned_preserve_rng <- function(code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(code)
}

g128_unsigned_project_mroast <- function(p_value, nrot = 9999L) {
  # Native project's extra directional BH projection. This does NOT replace
  # mroast's native FDR, PValue.Mixed or FDR.Mixed columns, including their NAs.
  pmax(p_value, stats::p.adjust(pmax(p_value - 0.5 / nrot, 0), "BH"))
}

g128_unsigned_ranks <- function(symbols, edge_chronic) {
  g128_unsigned_assert(length(symbols) == nrow(edge_chronic), "Rank input row mismatch")
  needed <- c("gene_id", "logFC", "QLF_F")
  g128_unsigned_assert(all(needed %in% names(edge_chronic)), "Missing native edgeR rank columns")
  rank_by_symbol <- data.table::data.table(
    symbol = symbols, gene_id = edge_chronic$gene_id,
    score = sign(edge_chronic$logFC) * sqrt(pmax(edge_chronic$QLF_F, 0))
  )
  rank_by_symbol <- rank_by_symbol[!is.na(symbol) & nzchar(symbol) & is.finite(score)]
  rank_by_symbol[, abs_score := abs(score)]
  data.table::setorder(rank_by_symbol, symbol, -abs_score, gene_id)
  rank_by_symbol <- unique(rank_by_symbol, by = "symbol")
  data.table::setorder(rank_by_symbol, -score, symbol)
  ranks <- rank_by_symbol$score
  names(ranks) <- rank_by_symbol$symbol
  payload <- paste0(paste(sprintf("%s\t%.17g", names(ranks), ranks),
                          collapse = "\n"), "\n")
  list(ranks = ranks, table = rank_by_symbol, payload = payload,
       sha256 = tolower(as.character(openssl::sha256(charToRaw(payload)))))
}

# Counts are retained/native-order counts, NOT a request to filter or fit them.
# Annotation columns are explicit and implement v2a_read_inputs' symbol trimming.
# Design row order is the caller's sample order; native model.matrix row labels
# need not equal GSM IDs. No canonicalization or sample/feature reordering occurs.
g128_unsigned_prepare <- function(counts, voom, design, contrasts, members,
                                  annotation, edge_chronic,
                                  gene_id_column = "gene_id",
                                  symbol_column = "gene_symbol",
                                  minimum_members = 10L) {
  spec <- g128_unsigned_spec()
  g128_unsigned_assert(is.matrix(counts) && is.numeric(counts) &&
                         all(is.finite(counts)) && all(counts >= 0),
                       "Retained finite nonnegative numeric counts required")
  g128_unsigned_assert(!is.null(rownames(counts)) && !anyDuplicated(rownames(counts)) &&
                         all(!is.na(rownames(counts)) & nzchar(trimws(rownames(counts)))) &&
                         !is.null(colnames(counts)) && !anyDuplicated(colnames(counts)) &&
                         all(!is.na(colnames(counts)) & nzchar(trimws(colnames(counts)))),
                       "Unique nonmissing nonblank count feature/sample names required")
  g128_unsigned_assert(is.list(voom) && is.matrix(voom$E) &&
                         all(is.finite(voom$E)) &&
                         identical(dim(voom$E), dim(counts)) &&
                         identical(dimnames(voom$E), dimnames(counts)),
                       "Native voom E/count domain and order must match")
  g128_unsigned_assert(is.matrix(design) && is.numeric(design) &&
                         nrow(design) == ncol(counts) && all(is.finite(design)) &&
                         identical(colnames(design), spec$groups) &&
                         all(design %in% c(0, 1)) && all(rowSums(design) == 1) &&
                         qr(design)$rank == ncol(design), "Invalid explicit six-group design")
  g128_unsigned_assert(is.matrix(contrasts) && is.numeric(contrasts) &&
                         identical(rownames(contrasts), colnames(design)) &&
                         sum(colnames(contrasts) == "chronic_combined") == 1L,
                       "Explicit chronic contrast/design coordinates required")
  contrast <- contrasts[, "chronic_combined"]
  g128_unsigned_assert(identical(as.numeric(contrast),
                                as.numeric(spec$chronic[colnames(design)])),
                       "Chronic scientific contrast changed")
  g128_unsigned_assert(length(minimum_members) == 1L && is.numeric(minimum_members) &&
                         is.finite(minimum_members) && minimum_members == 10L,
                       "Frozen minimum is 10 raw rows / subsequently 10 ranked symbols")
  g128_unsigned_assert(is.data.frame(annotation) && nrow(annotation) == nrow(counts) &&
                         all(c(gene_id_column, symbol_column) %in% names(annotation)),
                       "Explicit retained annotation required")
  annotation_ids <- trimws(as.character(annotation[[gene_id_column]]))
  g128_unsigned_assert(identical(annotation_ids, rownames(counts)),
                       "Annotation/count gene order mismatch")
  symbols <- trimws(as.character(annotation[[symbol_column]]))
  symbols[is.na(symbols) | !nzchar(symbols)] <- NA_character_
  g128_unsigned_assert(is.data.frame(edge_chronic) &&
                         identical(as.character(edge_chronic$gene_id), rownames(counts)),
                       "Native chronic edgeR/count gene order mismatch")
  if ("contrast" %in% names(edge_chronic)) {
    g128_unsigned_assert(all(edge_chronic$contrast == "chronic_combined"),
                         "Rank table must contain only chronic_combined")
  }
  needed <- c("collection", "program_id", "program_name", "gene_symbol", "test_class")
  g128_unsigned_assert(is.data.frame(members) && all(needed %in% names(members)),
                       "Frozen mapped member columns required")
  member_table <- data.table::copy(data.table::as.data.table(members))
  unsigned <- member_table[grepl("unsigned", test_class)]
  g128_unsigned_assert(all(unsigned$collection %in% spec$collections),
                       "Only frozen Reactome/SynGO unsigned collections are supported")
  symbol_to_row <- split(seq_along(symbols), symbols)
  ranking <- g128_unsigned_ranks(symbols, edge_chronic)
  collections <- vector("list", length(unique(unsigned$collection)))
  names(collections) <- unique(unsigned$collection)
  for (collection_name in names(collections)) {
    block <- unsigned[collection == collection_name]
    registry <- unique(block[, .(program_id, program_name)])
    index_list <- lapply(registry$program_id, function(id) {
      genes <- unique(block[program_id == id, gene_symbol])
      unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    })
    names(index_list) <- registry$program_id
    mapped_row_counts <- lengths(index_list)
    eligible <- mapped_row_counts >= minimum_members & mapped_row_counts <= 500L
    index_list <- index_list[eligible]
    paths <- lapply(index_list, function(index) unique(symbols[index]))
    paths <- lapply(paths, function(x) intersect(x, names(ranking$ranks)))
    paths <- paths[lengths(paths) >= minimum_members & lengths(paths) <= 500L]
    collections[[collection_name]] <- list(
      index = index_list, pathways = paths,
      registry = registry, mapped_row_counts = mapped_row_counts,
      eligible_mroast = eligible
    )
  }
  list(symbols = symbols, contrast = contrast, collections = collections,
       ranking = ranking, minimum_members = minimum_members)
}

g128_unsigned_chronic <- function(counts, voom, design, contrasts, members,
                                  annotation, edge_chronic, seed,
                                  nrot = 9999L, minimum_members = 10L,
                                  gene_id_column = "gene_id",
                                  symbol_column = "gene_symbol") {
  g128_unsigned_require_runtime()
  spec <- g128_unsigned_spec()
  g128_unsigned_assert(length(seed) == 1L && is.numeric(seed) && !is.na(seed) &&
                         seed == spec$seed, "Native enrichment seed is 20260819")
  g128_unsigned_assert(length(nrot) == 1L && is.numeric(nrot) && !is.na(nrot) &&
                         nrot == spec$nrot, "Frozen rotation count is 9999")
  seed <- as.integer(seed)
  nrot <- as.integer(nrot)
  # The original copy-isolated caller (not the helper) registers SerialParam().
  # Require the same caller state without mutating the process registry or
  # injecting additional fgsea arguments absent from the native source.
  backend <- BiocParallel::bpparam()
  g128_unsigned_assert(inherits(backend, "SerialParam") &&
                         is.null(BiocParallel::bpRNGseed(backend)),
                       "Native caller requires registered SerialParam with RNGseed NULL")
  g128_unsigned_assert(is.matrix(voom$weights) &&
                         identical(dim(voom$weights), dim(counts)) &&
                         (is.null(dimnames(voom$weights)) ||
                            identical(dimnames(voom$weights), dimnames(counts))) &&
                         all(is.finite(voom$weights)) && all(voom$weights > 0),
                       "Supplied native voom observation weights are required")
  prepared <- g128_unsigned_prepare(
    counts, voom, design, contrasts, members, annotation, edge_chronic,
    gene_id_column, symbol_column, minimum_members
  )
  g128_unsigned_preserve_rng({
    roast_rows <- list()
    fgsea_rows <- list()
    rng_calls <- list()
    cursor <- 1L
    for (collection_name in names(prepared$collections)) {
      index_list <- prepared$collections[[collection_name]]$index
      if (length(index_list) == 0L) next
      g128_unsigned_reset_rng(seed, spec$rng_kind)
      rng_initial <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      roast_native <- limma::mroast(
        voom, index = index_list, design = design, contrast = prepared$contrast,
        nrot = nrot, set.statistic = "mean", midp = TRUE, sort = "none"
      )
      rng_final <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      rng_calls[[length(rng_calls) + 1L]] <- list(
        collection = collection_name, method = "mroast", seed = seed,
        initial_state = rng_initial, final_state = rng_final,
        rng_kind = spec$rng_kind, final_rng_kind = RNGkind(),
        contrast = prepared$contrast, program_ids = names(index_list),
        ranking = list(rule = spec$ranking_rule, count = length(prepared$ranking$ranks),
                       sha256 = prepared$ranking$sha256),
        backend = list(class = class(backend), workers = BiocParallel::bpworkers(backend),
                       RNGseed = BiocParallel::bpRNGseed(backend))
      )
      roast <- data.table::as.data.table(roast_native, keep.rownames = "program_id")
      roast[, `:=`(
        collection = collection_name, contrast = "chronic_combined",
        BH_q = g128_unsigned_project_mroast(PValue, nrot),
        method = sprintf("mroast_mean_%d_rotations", nrot),
        package_version = as.character(utils::packageVersion("limma")),
        rng_seed = seed, rng_kind = spec$rng_kind[[1L]],
        rng_normal_kind = spec$rng_kind[[2L]], rng_sample_kind = spec$rng_kind[[3L]],
        rotation_count = nrot, set_statistic = "mean", midp = TRUE,
        sort_rule = "input_program_id_order", formal_sequential_inference = FALSE,
        candidate_generation_allowed = FALSE, s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
      )]
      roast_rows[[cursor]] <- roast
      paths <- prepared$collections[[collection_name]]$pathways
      if (length(paths) > 0L) {
        g128_unsigned_reset_rng(seed, spec$rng_kind)
        rng_initial <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        fg_native <- fgsea::fgseaMultilevel(
          pathways = paths, stats = prepared$ranking$ranks,
          minSize = minimum_members, maxSize = 500L,
          eps = 0, scoreType = "std", nPermSimple = 1000L
        )
        rng_final <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        rng_calls[[length(rng_calls) + 1L]] <- list(
          collection = collection_name, method = "fgseaMultilevel", seed = seed,
          initial_state = rng_initial, final_state = rng_final,
          rng_kind = spec$rng_kind, final_rng_kind = RNGkind(),
          contrast = prepared$contrast, program_ids = names(paths),
          ranking = list(rule = spec$ranking_rule, count = length(prepared$ranking$ranks),
                         sha256 = prepared$ranking$sha256),
          backend = list(class = class(backend), workers = BiocParallel::bpworkers(backend),
                         RNGseed = BiocParallel::bpRNGseed(backend))
        )
        fg <- data.table::as.data.table(fg_native)
        data.table::setnames(fg, "pathway", "program_id")
        fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ","),
                                   character(1L))]
        fg[, `:=`(
          collection = collection_name, contrast = "chronic_combined",
          BH_q = stats::p.adjust(pval, "BH"), method = "fgseaMultilevel_sensitivity",
          package_version = as.character(utils::packageVersion("fgsea")),
          rng_seed = seed, rng_kind = spec$rng_kind[[1L]],
          rng_normal_kind = spec$rng_kind[[2L]], rng_sample_kind = spec$rng_kind[[3L]],
          ranking_rule = spec$ranking_rule, ranking_count = length(prepared$ranking$ranks),
          ranking_sha256 = prepared$ranking$sha256, min_size = minimum_members,
          max_size = 500L, eps = 0, score_type = "std", n_perm_simple = 1000L,
          tie_rule = "score_desc_then_symbol_asc", formal_sequential_inference = FALSE,
          candidate_generation_allowed = FALSE, s2_status = "NOT_IDENTIFIABLE_NO_S1_INPUT"
        )]
        fgsea_rows[[cursor]] <- fg
      }
      cursor <- cursor + 1L
    }
    list(
      mroast = if (length(roast_rows)) data.table::rbindlist(roast_rows, fill = TRUE) else data.table::data.table(),
      fgsea = if (length(fgsea_rows)) data.table::rbindlist(fgsea_rows, fill = TRUE) else data.table::data.table(),
      prepared = prepared, rng_calls = rng_calls
    )
  })
}
