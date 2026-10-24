# Additive, content-only check after the independent v2 branch reconstruction.
# Does not read files, fit models, authenticate custody, or authorize publication.
g128_comparison_guard <- function(primary, complete29, observed_gene,
                                  observed_summary, max_features) {
  need <- function(ok) if (!isTRUE(ok)) stop("BRANCH_COMPARISON_UNVERIFIED", call. = FALSE)
  contrasts <- c("chronic_combined", "D28_vs_control", "D42_vs_control",
                 "acute_combined", "D4_vs_control", "D7_vs_control", "D14_vs_control")
  number <- function(x) typeof(x) %in% c("integer", "double") &&
    all(is.na(x) | is.finite(x))
  text <- function(x) is.character(x) && !anyNA(x) && all(nzchar(x)) &&
    all(x == trimws(x)) && !any(grepl('[\t\r\n"]', x))
  numeric_equal <- function(a, b) {
    if (length(a) != length(b) || !number(a) || !number(b)) return(FALSE)
    same_missing <- is.na(a) & is.na(b)
    finite <- is.finite(a) & is.finite(b)
    all(same_missing | (finite & abs(a - b) <= 2e-7 * pmax(1, abs(a), abs(b))))
  }
  need(length(max_features) == 1L && typeof(max_features) %in% c("integer", "double") &&
         is.finite(max_features) && max_features >= 2 && max_features <= .Machine$integer.max &&
         max_features == floor(max_features))
  edge <- function(x) {
    need(is.data.frame(x) && identical(names(x), c("gene_id", "contrast", "logFC", "BH_q")) &&
           nrow(x) > 0 && nrow(x) <= 7 * max_features)
    x <- as.data.frame(x)
    need(text(x$gene_id) && text(x$contrast) && number(x$logFC) && number(x$BH_q) &&
           all(is.na(x$BH_q) | (x$BH_q >= 0 & x$BH_q <= 1)) &&
           all(x$contrast %in% contrasts) && !anyDuplicated(x[c("gene_id", "contrast")]))
    need(all(table(factor(x$contrast, levels = contrasts)) == length(unique(x$gene_id))) &&
           all(table(x$gene_id) == 7L))
    x
  }
  p <- edge(primary); q <- edge(complete29)
  need(length(union(p$gene_id, q$gene_id)) <= max_features)
  names(p)[3:4] <- c("logFC_all30", "q_all30")
  names(q)[3:4] <- c("logFC_complete29", "q_complete29")
  expected <- merge(p, q, by = c("gene_id", "contrast"), all = TRUE, sort = TRUE)
  need(is.data.frame(observed_gene) && identical(names(observed_gene), names(expected)) &&
         nrow(observed_gene) == nrow(expected) &&
         !anyDuplicated(as.data.frame(observed_gene)[c("gene_id", "contrast")]))
  actual <- as.data.frame(observed_gene)
  need(text(actual$gene_id) && text(actual$contrast))
  actual <- actual[order(actual$gene_id, actual$contrast), , drop = FALSE]
  need(identical(actual$gene_id, expected$gene_id) && identical(actual$contrast, expected$contrast))
  for (name in names(expected)[3:6]) need(numeric_equal(actual[[name]], expected[[name]]))
  sums <- lapply(contrasts, function(name) {
    x <- expected[expected$contrast == name, , drop = FALSE]
    common <- is.finite(x$logFC_all30) & is.finite(x$logFC_complete29)
    # The old frozen analyzer's correlation needs at least one complete pair.
    need(any(common))
    data.frame(contrast = name, all30_tested = sum(is.finite(x$logFC_all30)),
      complete29_tested = sum(is.finite(x$logFC_complete29)), common_tested = sum(common),
      sign_concordance = mean(sign(x$logFC_all30) == sign(x$logFC_complete29), na.rm = TRUE),
      spearman_logFC = suppressWarnings(stats::cor(x$logFC_all30, x$logFC_complete29,
                                                  method = "spearman", use = "complete.obs")))
  })
  summary <- do.call(rbind, sums)
  need(is.data.frame(observed_summary) && identical(names(observed_summary), names(summary)) &&
         nrow(observed_summary) == 7L && text(observed_summary$contrast) &&
         !anyDuplicated(observed_summary$contrast) && setequal(observed_summary$contrast, contrasts))
  actual_summary <- as.data.frame(observed_summary)
  actual_summary <- actual_summary[match(contrasts, actual_summary$contrast), , drop = FALSE]
  for (name in names(summary)[2:4]) {
    need(number(actual_summary[[name]]) && identical(as.numeric(actual_summary[[name]]),
                                                   as.numeric(summary[[name]])))
  }
  for (name in names(summary)[5:6]) need(numeric_equal(actual_summary[[name]], summary[[name]]))
  list(status = "PASS_BRANCH_COMPARISON_CONTENT_ONLY", comparison_numeric_checked = TRUE,
       summary_numeric_checked = TRUE, missingness_checked = TRUE, model_run = FALSE,
       prior_independent_branch_validation_required = TRUE, closure_authenticated = FALSE,
       production_go = FALSE, publication_authorized = FALSE)
}
