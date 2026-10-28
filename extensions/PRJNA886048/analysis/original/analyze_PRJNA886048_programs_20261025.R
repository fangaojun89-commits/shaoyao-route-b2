# Competitive membership tests use the existing gene fit; no refitting here.
# PROGENy/DoRothEA membership direction is not a signed activity estimate.
prjna886048_camera <- function(gene_fit, edges, min_size = 10L, max_size = 500L,
                              reference_terms = NULL) {
  stopifnot(all(c("collection", "term_id", "term_name", "gene_id") %in% names(edges)),
            min_size >= 2L, max_size >= min_size,
            identical(colnames(gene_fit$voom$E), rownames(gene_fit$design)),
            identical(rownames(gene_fit$voom$E), gene_fit$sensitivity_voom$gene_id),
            !anyDuplicated(edges[c("collection", "term_id", "gene_id")]))
  universe <- rownames(gene_fit$voom$E)
  effect <- setNames(gene_fit$sensitivity_voom$logFC, universe)
  coverage <- list()
  output <- list()
  if (is.null(reference_terms)) reference_terms <- unique(edges[c("collection", "term_id", "term_name")])
  stopifnot(all(c("collection", "term_id", "term_name") %in% names(reference_terms)),
            !anyDuplicated(reference_terms[c("collection", "term_id")]),
            all(paste(edges$collection, edges$term_id) %in%
                  paste(reference_terms$collection, reference_terms$term_id)))
  for (collection in sort(unique(reference_terms$collection))) {
    e <- edges[edges$collection == collection, , drop = FALSE]
    terms <- reference_terms[reference_terms$collection == collection, , drop = FALSE]
    sets <- split(e$gene_id, e$term_id)
    for (id in setdiff(terms$term_id, names(sets))) sets[id] <- list(character())
    mapped <- lapply(sets, function(x) which(universe %in% x))
    eligible <- lengths(mapped) >= min_size & lengths(mapped) <= max_size &
      lengths(mapped) < length(universe)
    cov <- data.frame(collection = collection, term_id = names(sets),
                      reference_members = lengths(sets), expressed_members = lengths(mapped),
                      tested = eligible, row.names = NULL)
    coverage[[collection]] <- cov
    if (!any(eligible)) next
    indices <- mapped[eligible]
    tested <- limma::camera(gene_fit$voom, index = indices, design = gene_fit$design,
                            contrast = gene_fit$contrast, inter.gene.cor = 0.01,
                            use.ranks = FALSE, sort = FALSE)
    stopifnot(identical(rownames(tested), names(indices)),
              all(is.finite(tested$PValue)), all(tested$PValue >= 0 & tested$PValue <= 1))
    labels <- vapply(names(indices), function(id) {
      values <- unique(terms$term_name[terms$term_id == id])
      stopifnot(length(values) == 1L)
      values
    }, character(1))
    output[[collection]] <- data.frame(
      collection = collection, term_id = names(indices), term_name = unname(labels),
      NGenes = tested$NGenes, Direction = tested$Direction,
      PValue = tested$PValue, FDR = p.adjust(tested$PValue, "BH"),
      mean_member_log2FC = vapply(indices, function(i) mean(effect[i]), numeric(1)),
      median_member_log2FC = vapply(indices, function(i) median(effect[i]), numeric(1)),
      background_genes = length(universe), row.names = NULL)
  }
  empty <- data.frame(collection = character(), term_id = character(), term_name = character(),
                       NGenes = integer(), Direction = character(), PValue = numeric(), FDR = numeric(),
                       mean_member_log2FC = numeric(), median_member_log2FC = numeric(),
                       background_genes = integer())
  list(camera = if (length(output)) do.call(rbind, output) else empty,
       coverage = do.call(rbind, coverage),
       method = "CAMERA on existing voom object; fixed inter-gene correlation 0.01; BH within collection",
       contrast = "baclofen minus vehicle",
       effect_summary = "Unweighted mean/median member gene log2FC; descriptive, not pathway activation",
       scope = "Unsigned membership comparison; signed footprint activities remain separate and pending")
}
