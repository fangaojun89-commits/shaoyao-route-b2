# Functions-only metadata wiring. No I/O, model fitting, extraction or tail run.
# Pass the registry returned by the pinned jgw03_family_registry().
# Outer expressions are explicit runner-wiring recommendations, not evidence
# that a reference has run or that its numerical acceptance rules are frozen.
jgw03_reference_field_map <- function(registry) {
  need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
  columns <- c("accession", "method", "collection", "contrast", "space",
    "null_draw_n", "observed_once")
  need(is.data.frame(registry) && all(columns %in% names(registry)) &&
    nrow(registry) == 123L && !anyNA(registry[, columns]) &&
    !anyDuplicated(registry[, columns[1:5]]), "Complete unique 123-family registry required")
  methods <- c("edgeR_QL", "limma_voom_trend_TRUE", "DESeq2_native_IF_Cooks",
    "CAMERA", "mean_z_limma", "ULM_activity_limma", "fgseaMultilevel",
    "mroast1999_equivalent_not_historical_authentication")
  labels <- c(high_vs_vehicle = "high", ordered_dose_per_step = "trend",
    log10_dose_plus_1 = "logtrend")
  need(all(registry$method %in% methods) && all(registry$contrast %in% names(labels)) &&
    sum(registry$space == "four_group_MC") == 89L &&
    sum(registry$space == "high_companion_exact") == 34L &&
    all(registry$null_draw_n[registry$space == "four_group_MC"] == 10000L) &&
    all(registry$null_draw_n[registry$space == "high_companion_exact"] == 20L) &&
    all(registry$contrast[registry$space == "high_companion_exact"] == "high_vs_vehicle") &&
    all(registry$observed_once), "Frozen method/contrast/space geometry differs")
  answer <- lapply(seq_len(nrow(registry)), function(i) {
    row <- registry[i, , drop = FALSE]
    method <- as.character(row$method)
    label <- unname(labels[[as.character(row$contrast)]])
    collection <- as.character(row$collection)
    gene <- method %in% methods[1:3]
    need((gene && collection == "genes") ||
      (method %in% c("CAMERA", "mean_z_limma", "fgseaMultilevel") &&
       collection %in% c("Reactome97", "QuickGO_SynGO_sensitivity")) ||
      (method == "ULM_activity_limma" &&
       collection %in% c("PROGENy_top500_ULM", "DoRothEA_AB_ULM")) ||
      (method == methods[8] && row$accession == "GSE289929" && collection == "Reactome97"),
      "Method/collection/accession domain differs")
    need(!method %in% c(methods[2:3], methods[7:8]) || label %in% c("high", "trend"),
      "No additional log-dose sensitivity family may be invented")
    key <- if (gene) "raw_gene_id" else if (method %in%
      c("mean_z_limma", "ULM_activity_limma")) "feature" else "term_id"
    selector <- switch(method,
      edgeR_QL = sprintf("result$gene_edgeR$%s", label),
      limma_voom_trend_TRUE = sprintf("result$gene_voom$%s$value", label),
      DESeq2_native_IF_Cooks = sprintf("result$gene_DESeq2$%s$value", label),
      CAMERA = sprintf("result$programs$%s$%s", collection, label),
      mean_z_limma = sprintf("result$programs$%s$score_fit$%s", collection, label),
      ULM_activity_limma = sprintf("result$activities$%s$%s", collection, label),
      fgseaMultilevel = sprintf("result$programs$%s$fgsea_%s", collection, label),
      mroast1999_equivalent_not_historical_authentication =
        sprintf("result$programs$Reactome97$roast$%s", label))
    effect <- switch(method, edgeR_QL = "logFC", limma_voom_trend_TRUE = "logFC",
      DESeq2_native_IF_Cooks = "log2FoldChange", CAMERA = "NONE; Direction is not an effect size",
      mean_z_limma = "effect", ULM_activity_limma = "effect", fgseaMultilevel = "NES",
      mroast1999_equivalent_not_historical_authentication =
        "NONE; PropUp-PropDown is optional descriptive proportion difference only")
    statistic <- switch(method, edgeR_QL = "F", limma_voom_trend_TRUE = "t",
      DESeq2_native_IF_Cooks = "stat", CAMERA = "NONE",
      mean_z_limma = "t", ULM_activity_limma = "t", fgseaMultilevel = "NES",
      mroast1999_equivalent_not_historical_authentication = "NONE")
    outer <- switch(method, edgeR_QL = "F (unchanged, including native roundoff)",
      limma_voom_trend_TRUE = "abs(t)", DESeq2_native_IF_Cooks = "abs(stat)",
      CAMERA = "-log(PValue) (derived P-strength, not native statistic)",
      mean_z_limma = "abs(t)", ULM_activity_limma = "abs(t)", fgseaMultilevel = "abs(NES)",
      mroast1999_equivalent_not_historical_authentication =
        "-log(PValue) (directional P-strength; Mixed remains ancillary)")
    p <- switch(method, limma_voom_trend_TRUE = "P.Value",
      DESeq2_native_IF_Cooks = "pvalue", fgseaMultilevel = "pval", "PValue")
    q <- switch(method, limma_voom_trend_TRUE = "adj.P.Val",
      DESeq2_native_IF_Cooks = "padj", fgseaMultilevel = "padj", "FDR")
    q_rule <- switch(method,
      DESeq2_native_IF_Cooks = "Original DESeq2 results(alpha=.05) padj with default independentFiltering/Cooks; never overwrite",
      fgseaMultilevel = "Native fgsea padj over returned eligible pathways for this collection/contrast; never overwrite",
      mroast1999_equivalent_not_historical_authentication =
        "Native directional FDR=pmax(p.adjust(PValue-.5/2000,'BH'),PValue); Mixed separately same formula",
      "Original per-cell/method/collection/contrast stats::p.adjust(P,'BH') NA-aware default; no pooled family")
    file_stem <- c(high = "high_vs_vehicle", trend = "ordered_trend", logtrend = "logdose_trend")[[label]]
    relative <- switch(method,
      edgeR_QL = sprintf("genes/edgeR_QL_%s_all_genes.tsv.gz", file_stem),
      limma_voom_trend_TRUE = sprintf("genes/limma_voom_%s_all_genes.tsv.gz", file_stem),
      DESeq2_native_IF_Cooks = sprintf("genes/DESeq2_rounded_%s_all_genes.tsv.gz", file_stem),
      CAMERA = sprintf("programs/%s_camera_%s.tsv", collection, file_stem),
      mean_z_limma = sprintf("programs/%s_score_%s.tsv", collection, label),
      ULM_activity_limma = sprintf("activities/%s_%s.tsv",
        if (collection == "PROGENy_top500_ULM") "PROGENy" else "DoRothEA", label),
      fgseaMultilevel = sprintf("programs/%s_fgsea_%s.tsv", collection, file_stem),
      mroast1999_equivalent_not_historical_authentication =
        sprintf("programs/Reactome97_mroast1999_%s.tsv", file_stem))
    data.frame(row, kernel_label = label, native_key = key,
      native_order = "Preserve this table's key order; freeze observed keys then require identical keys/order for every allocation; no positional cross-method join",
      result_selector = selector, original_table = paste0(
        "results/jgw03_three_cell_lines/", row$accession, "/", relative),
      native_effect = effect, native_statistic = statistic,
      proposed_outer_statistic = outer, native_P = p, native_q = q,
      native_q_rule = q_rule, native_significance = "q<=0.05; NA remains unavailable",
      missing_policy = "Keep all fixed keys and all NA/NaN/+Inf/-Inf counts; no fill/drop or smaller reference denominator",
      outer_rule = if (row$space == "four_group_MC")
        "(1+b)/10001; b from 10000 distinct nonidentity labels; observed separate once" else
        "total_hits/20; complete high/vehicle companion contains observed once; no extra plus-one",
      stringsAsFactors = FALSE, check.names = FALSE)
  })
  answer <- do.call(rbind, answer)
  rownames(answer) <- NULL
  answer
}
