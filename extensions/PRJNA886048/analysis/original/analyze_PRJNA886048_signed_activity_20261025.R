# Signed footprint scores and independent-library contrasts; no gene refitting.
prjna886048_activity_contrast <- function(scores, samples) {
  setup <- prjna886048_design(samples)
  stopifnot(is.matrix(scores), identical(colnames(scores), samples$run), all(is.finite(scores)))
  fit <- limma::lmFit(scores, setup$design)
  fit <- limma::contrasts.fit(fit, setup$contrast)
  effect <- as.numeric(fit$coefficients[, 1L])
  se <- as.numeric(fit$stdev.unscaled[, 1L] * fit$sigma)
  df <- fit$df.residual
  total <- rowSums((scores - rowMeans(scores))^2)
  degenerate <- !is.finite(se) | (fit$sigma^2 * df <= 64 * .Machine$double.eps * pmax(1, total))
  t <- effect / se
  p <- 2 * pt(-abs(t), df)
  critical <- qt(0.975, df)
  low <- effect - critical * se
  high <- effect + critical * se
  se[degenerate] <- t[degenerate] <- p[degenerate] <- low[degenerate] <- high[degenerate] <- NA_real_
  data.frame(term_id = rownames(scores), estimate = effect, SE = se, t = t, df = df,
             PValue = p, FDR = p.adjust(p, "BH", n = nrow(scores)),
             ci_low = low, ci_high = high,
             mean_vehicle = rowMeans(scores[, setup$group == "Saline solution", drop = FALSE]),
             mean_baclofen = rowMeans(scores[, setup$group == "Baclofen", drop = FALSE]),
             status = ifelse(degenerate, "NUMERICALLY_DEGENERATE_NO_INFERENCE", "TESTED"),
             row.names = NULL)
}

prjna886048_signed_activity <- function(gene_fit, edges, samples, reference_terms) {
  prjna886048_design(samples)
  expression <- gene_fit$voom$E
  stopifnot(identical(colnames(expression), samples$run), !anyDuplicated(rownames(expression)),
            all(is.finite(expression)), all(c("collection", "term_id", "gene_id", "weight") %in% names(edges)),
            !anyDuplicated(edges[c("collection", "term_id", "gene_id")]),
            !anyDuplicated(reference_terms[c("collection", "term_id")]))
  all_scores <- list()
  all_results <- list()
  coverage <- list()
  for (collection in c("PROGENy", "DoRothEA")) {
    e <- edges[edges$collection == collection & edges$gene_id %in% rownames(expression), , drop = FALSE]
    e$weight <- as.numeric(e$weight)
    stopifnot(all(is.finite(e$weight)))
    e <- e[e$weight != 0, , drop = FALSE]
    terms <- reference_terms$term_id[reference_terms$collection == collection]
    stopifnot(all(e$term_id %in% terms))
    counts <- table(factor(e$term_id, levels = terms))
    keep <- counts >= 10L & counts <= 500L
    coverage[[collection]] <- data.frame(collection = collection, term_id = terms,
                                          nonzero_expressed_targets = as.integer(counts),
                                          eligible = as.logical(keep))
    if (!any(keep)) next
    e <- e[e$term_id %in% terms[keep], , drop = FALSE]
    network <- data.frame(source = e$term_id, target = e$gene_id, mor = e$weight)
    long <- decoupleR::run_ulm(mat = expression, network = network,
                               .source = "source", .target = "target", .mor = "mor",
                               minsize = 10L, center = FALSE)
    stopifnot(!anyDuplicated(long[c("source", "condition")]),
              nrow(long) == sum(keep) * ncol(expression), all(is.finite(long$score)),
              setequal(long$source, terms[keep]), setequal(long$condition, samples$run))
    scores <- matrix(NA_real_, nrow = sum(keep), ncol = ncol(expression),
                      dimnames = list(terms[keep], samples$run))
    scores[cbind(match(long$source, rownames(scores)), match(long$condition, colnames(scores)))] <- long$score
    stopifnot(!anyNA(scores))
    results <- prjna886048_activity_contrast(scores, samples)
    results$collection <- collection
    results$expressed_targets <- as.integer(counts[match(results$term_id, terms)])
    results$score_background_genes <- nrow(expression)
    all_scores[[collection]] <- scores
    all_results[[collection]] <- results
  }
  list(scores = all_scores, contrasts = if (length(all_results)) do.call(rbind, all_results) else data.frame(),
       coverage = do.call(rbind, coverage),
       method = "decoupleR ULM, center=FALSE, full filtered gene background; ordinary six-vs-six score contrast",
       inference = "Animal/library-level t test with 10 residual df and 95% CI; BH within collection; no gene-level ULM P used for treatment inference",
       units = "Difference in signed ULM score, not log2 fold change or measured activity",
       scope = "Retrospective pharmacological context; cross-source raw score magnitudes are not interchangeable")
}
