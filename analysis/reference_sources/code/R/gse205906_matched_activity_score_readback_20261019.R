# Functions only: caller-held NEW matched checkpoint and bound resources.
# Independent base-R mapping/aggregation and direct ULM arithmetic; no file IO,
# producer import, RNG, iterative model fitting or acceptance-file creation.

g205a_need <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
}

g205a_hash <- function(value) digest::digest(value, "sha256", serializeVersion = 3L)

g205a_equal <- function(actual, expected, label) {
  g205a_need(is.numeric(actual) && is.numeric(expected) &&
    identical(dim(actual), dim(expected)) && length(actual) == length(expected) &&
    identical(is.na(actual), is.na(expected)) && identical(is.nan(actual), is.nan(expected)) &&
    identical(is.finite(actual), is.finite(expected)), paste0(label, ": type/domain"))
  finite <- is.finite(expected)
  infinite <- !finite & !is.na(expected)
  g205a_need(all(actual[infinite] == expected[infinite]) &&
    all(abs(actual[finite] - expected[finite]) <= 1e-12 + 1e-10 * abs(expected[finite])),
    paste0(label, ": values"))
  invisible(TRUE)
}

g205a_collapse <- function(expression, annotation) {
  g205a_need(is.matrix(expression) && is.numeric(expression) &&
    all(is.finite(expression)) && !is.null(rownames(expression)) &&
    !is.null(colnames(expression)) && !anyDuplicated(rownames(expression)) &&
    !anyDuplicated(colnames(expression)) &&
    identical(as.character(annotation$gene_id), rownames(expression)), "Activity log-E/annotation keys")
  symbols <- as.character(annotation$symbol)
  symbols_kept <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  answer <- matrix(numeric(), length(symbols_kept), ncol(expression),
    dimnames = list(symbols_kept, colnames(expression)))
  for (i in seq_along(symbols_kept)) {
    rows <- which(!is.na(symbols) & symbols == symbols_kept[[i]])
    # Aggregate log-expression, not raw counts, precision weights or networks.
    answer[i, ] <- colMeans(expression[rows, , drop = FALSE])
  }
  answer
}

g205a_networks <- function(gene_resources, orthology_table) {
  orthology <- as.data.frame(orthology_table)
  g205a_need(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in%
    names(orthology)), "Activity orthology columns missing")
  selected <- which(orthology$group_method == "NCBI Ortholog" & orthology$other_tax_id == 10116L)
  rat_map <- unique(data.frame(human_symbol = as.character(orthology$human_symbol[selected]),
    rat_symbol = as.character(orthology$other_symbol[selected]), stringsAsFactors = FALSE))
  g205a_need(nrow(rat_map) > 0L && !anyNA(rat_map) &&
    !anyDuplicated(rat_map$human_symbol) && !anyDuplicated(rat_map$rat_symbol),
    "Activity rat orthology must be strict one-to-one")
  progeny <- as.data.frame(gene_resources$progeny_human)
  doro <- as.data.frame(gene_resources$dorothea_rat)
  g205a_need(all(c("pathway", "gene_symbol", "weight") %in% names(progeny)) &&
    all(c("tf", "target", "mor", "confidence") %in% names(doro)), "Activity resource columns missing")
  index <- match(as.character(progeny$gene_symbol), rat_map$human_symbol)
  mapped <- which(!is.na(index))
  prog <- unique(data.frame(source = as.character(progeny$pathway[mapped]),
    target = rat_map$rat_symbol[index[mapped]], mor = as.numeric(progeny$weight[mapped]),
    stringsAsFactors = FALSE))
  # The frozen DoRothEA input is rat A/B already; never remap or silently filter C.
  doro <- unique(data.frame(source = as.character(doro$tf), target = as.character(doro$target),
    mor = as.numeric(doro$mor), confidence = as.character(doro$confidence), stringsAsFactors = FALSE))
  g205a_need(!anyNA(doro$confidence) && all(doro$confidence %in% c("A", "B")),
    "Activity DoRothEA confidence must be A/B")
  expected_progeny <- c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK", "NFkB",
    "PI3K", "p53", "TGFb", "TNFa", "Trail", "VEGF", "WNT")
  g205a_need(setequal(unique(prog$source), expected_progeny) &&
    length(unique(doro$source)) == 117L, "Frozen activity source-program domain")
  list(PROGENy = prog, DoRothEA_AB = doro)
}

g205a_reconstruct <- function(expression, annotation, gene_resources, orthology_table) {
  collapsed <- g205a_collapse(expression, annotation)
  networks <- g205a_networks(gene_resources, orthology_table)
  scores <- coverage <- list()
  for (collection in names(networks)) {
    net <- networks[[collection]]
    g205a_need(!anyNA(net[c("source", "target")]) && all(nzchar(net$source)) &&
      all(nzchar(net$target)) && is.numeric(net$mor),
      "Activity network keys/weights invalid")
    ids <- unique(net$source)
    cov <- do.call(rbind, lapply(ids, function(id) {
      targets <- unique(net$target[net$source == id])
      tested <- length(intersect(targets, rownames(collapsed)))
      data.frame(source = id, source_feature_n = length(targets), tested_feature_n = tested,
        collection = collection, eligible_10_500 = tested >= 10L && tested <= 500L,
        tested = tested >= 10L && tested <= 500L,
        confidence_domain = if (collection == "PROGENy") "NOT_APPLICABLE" else
          paste(sort(unique(net$confidence)), collapse = ","),
        mapping_rule = if (collection == "PROGENy")
          "NCBI_Datasets_REST_v18.35.0_snapshot_2026-08-18_strict_one_to_one" else
          "frozen_rat_AB_strict_one_to_one_network", stringsAsFactors = FALSE)
    }))
    eligible <- sort(cov$source[cov$tested])
    answer <- matrix(numeric(), length(eligible), ncol(collapsed),
      dimnames = list(if (length(eligible)) eligible else NULL, colnames(collapsed)))
    for (i in seq_along(eligible)) {
      edges <- net[net$source == eligible[[i]] & net$target %in% rownames(collapsed), , drop = FALSE]
      g205a_need(all(is.finite(edges$mor)), "Native tested activity weights nonfinite")
      # decoupleR rejects repeated source/target pairs after confidence projection.
      g205a_need(!anyDuplicated(edges[c("source", "target")]), "Repeated native activity edges")
      weight <- numeric(nrow(collapsed))
      weight[match(edges$target, rownames(collapsed))] <- edges$mor
      g205a_need(length(weight) > 2L, "Activity expression universe too small")
      # Full expression-symbol universe, including non-target zero weights.
      # Direct installed-ULM equation, no run_ulm() or downstream limma refit.
      correlation <- as.numeric(stats::cor(weight, collapsed))
      degrees <- length(weight) - 2L
      answer[i, ] <- correlation * sqrt(degrees /
        ((1 - correlation + 1e-20) * (1 + correlation + 1e-20)))
    }
    scores[[collection]] <- answer
    coverage[[collection]] <- cov
  }
  list(scores = scores, coverage = do.call(rbind, coverage),
    expression_symbol_ids = rownames(collapsed))
}

g205a_readback <- function(checkpoint, gene_resources, orthology_table) {
  g205a_need(identical(checkpoint$schema, "GSE205906_MATCHED_SAVED_NATIVE_CHECKPOINT/v1") &&
    isTRUE(checkpoint$statistical_fit_executed), "New fitted checkpoint required")
  packet <- checkpoint$prepared_packet
  result <- checkpoint$result
  g205a_need(identical(result$schema, "gse205906-matched-downstream-result-v1") &&
    identical(result$binding, packet$binding) && identical(result$prepared_sha256, packet$sha256) &&
    identical(result$statistical_source_pins, packet$source_pins) &&
    identical(result$candidate_admission, FALSE), "Activity prepared/result association")
  g205a_need(is.list(gene_resources) &&
    all(c("progeny_human", "dorothea_rat") %in% names(gene_resources)) &&
    identical(g205a_hash(gene_resources), packet$binding$held_object_hashes$gene_resources) &&
    identical(g205a_hash(orthology_table), packet$binding$held_object_hashes$orthology_table),
    "Activity held resource/orthology binding")
  g205a_need(identical(as.character(getRversion()), result$runtime$R) &&
    identical(result$runtime$packages$decoupleR, "2.17.0") &&
    identical(as.character(utils::packageVersion("decoupleR")), "2.17.0"),
    "Saved activity runtime differs")
  counts <- packet$counts
  base <- result$native$base_fit
  annotation <- as.data.frame(packet$annotation$annotation)
  keep <- base$keep
  g205a_need(is.matrix(counts) && is.numeric(counts) && identical(dim(counts), c(32883L, 9L)) &&
    identical(colnames(counts), paste0("GSM", 6235292:6235300)) &&
    !is.null(rownames(counts)) && !anyNA(rownames(counts)) && all(nzchar(rownames(counts))) &&
    !anyDuplicated(rownames(counts)) && all(is.finite(counts)) &&
    all(counts >= 0 & counts == floor(counts)) &&
    all(c("gene_id", "symbol") %in% names(annotation)) &&
    identical(as.character(annotation$gene_id), rownames(counts)) && is.logical(keep) &&
    length(keep) == nrow(counts) && !anyNA(keep) &&
    identical(base$y$counts, counts[keep, , drop = FALSE]), "Activity retained count/annotation domain")
  retained <- annotation[keep, , drop = FALSE]
  expression <- base$voom$E
  g205a_need(identical(rownames(expression), rownames(base$y$counts)) &&
    identical(colnames(expression), colnames(counts)), "Activity voom log-E keys")
  expected <- g205a_reconstruct(expression, retained, gene_resources, orthology_table)
  saved <- result$own_activity_scores
  g205a_need(is.list(saved) && identical(names(saved), c("PROGENy", "DoRothEA_AB")),
    "Activity score collections")
  for (collection in names(saved)) {
    g205a_need(is.matrix(saved[[collection]]) &&
      identical(dimnames(saved[[collection]]), dimnames(expected$scores[[collection]])),
      "Activity score ordered keys")
    g205a_equal(saved[[collection]], expected$scores[[collection]], "Activity scores")
  }
  observed <- as.data.frame(result$mapping_coverage)
  wanted <- expected$coverage
  g205a_need("collection" %in% names(observed) && !anyNA(observed$collection) &&
    all(observed$collection %in% c("Reactome", "SynGO_sensitivity", "PROGENy", "DoRothEA_AB")),
    "Activity coverage collection domain")
  for (collection in names(saved)) {
    actual <- observed[observed$collection == collection, , drop = FALSE]
    target <- wanted[wanted$collection == collection, , drop = FALSE]
    g205a_need(all(names(target) %in% names(actual)) && !anyNA(actual$source) &&
      !anyDuplicated(actual$source) && setequal(as.character(actual$source), target$source),
      "Activity coverage source domain")
    actual <- actual[match(target$source, actual$source), , drop = FALSE]
    for (field in names(target)) {
      same <- if (is.numeric(target[[field]])) is.numeric(actual[[field]]) &&
        identical(as.numeric(actual[[field]]), as.numeric(target[[field]])) else
        identical(actual[[field]], target[[field]])
      g205a_need(same, paste0("Activity coverage field: ", field))
    }
  }
  list(status = "HELD_NEW_CHECKPOINT_ACTIVITY_SCORE_COMPONENT_VERIFIED",
    collapsed_expression_symbols = length(expected$expression_symbol_ids),
    program_counts = vapply(expected$scores, nrow, integer(1)),
    mapping_aggregation_eligibility_and_score_arithmetic_verified = TRUE,
    relative_tolerance = 1e-10, absolute_tolerance = 1e-12,
    downstream_fit_consumed_exact_score_authenticated = FALSE,
    upstream_normalization_or_external_custody_verified = FALSE,
    iterative_model_fits = 0L, candidate_admission = FALSE,
    full_independent_numeric_acceptance = FALSE)
}
