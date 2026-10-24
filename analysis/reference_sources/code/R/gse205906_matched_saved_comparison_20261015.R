# Functions-only, saved-table/count descriptions for the same nine libraries.
# No file IO, fitting, normalization, RNG, BH recalculation, nulls or admission.
# Caller authenticates both actual input/output manifests and supplies each
# branch's COMPLETE source-domain ledger; no future STAR hash/count is invented.
# Pipeline: bound held counts/tables/domains -> descriptions -> caller-owned
# fresh outputs. Existing frozen runtime/renv is reused, never restored here.
g205c_need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

g205c_ids <- function(x) {
  is.character(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
}

g205c_cor <- function(left, right, method) {
  g205c_need(is.numeric(left) && is.numeric(right) && length(left) == length(right) &&
    is.character(method) && length(method) == 1L && !is.na(method) &&
    method %in% c("pearson", "spearman", "kendall"), "Paired numeric correlation vectors/method required")
  keep <- is.finite(left) & is.finite(right)
  if (sum(keep) < 2L || length(unique(left[keep])) < 2L ||
      length(unique(right[keep])) < 2L) return(NA_real_)
  stats::cor(left[keep], right[keep], method = method)
}

g205c_policy <- function() {
  list(
    schema = "gse205906-matched-saved-comparison-policy-v1",
    contrasts = c("local_disease", "intervention"), q_rule = "STRICT_LT_0.05",
    q_scope = "EACH_ORIGINAL_COMPLETE_METHOD_COLLECTION_CONTRAST_FAMILY",
    topk_default = "PENDING_POLICY_BINDING",
    available_topk_policy = "EXISTING_BASELINE_P_TOP20_PAIRED_FINITE_ALL_BOUNDARY_TIES",
    topk_source = list(file = "code/R/gse205906_program_baselines_20261014_v2.R",
      function_name = "pair_metrics", sha256 = "68c01474c53840441977b9bba571f816ec7d2083677fbc78b112af46f3b64af1",
      scope = "EXISTING_PROGRAM_BASELINE_RULE_NOT_AUTOMATIC_ALL_METHOD_POLICY"),
    proposal_sha256 = "bb5207f1a0f616f0a2c2d9c4d7bbd5bdf32a17041b631a5a0c7a8dbdc3e8d622",
    SAP_topk_numeric_value_fixed = FALSE,
    repeated_libraries_not_independent_estimates = TRUE,
    effect_difference_uncertainty_computed = FALSE, native_q_recomputed = FALSE,
    candidate_or_gate_change = FALSE, biological_validation = FALSE,
    file_custody_authenticated_here = FALSE)
}

# Natural log1p count scale. M=right-left, A=(right+left)/2; no library-size
# adjustment. Exact same gene/sample ordering is required, never silently joined.
# input_pairs is a named full-pair denominator, not mapped reads or BAM records.
g205c_counts <- function(left, right, expected_gene_ids, input_pairs) {
  samples <- paste0("GSM", 6235292:6235300)
  g205c_need(g205c_ids(expected_gene_ids) && length(expected_gene_ids) > 0L,
    "Explicit complete caller-bound gene order required")
  for (value in list(left, right)) {
    g205c_need(is.matrix(value) && is.numeric(value) && !is.logical(value) &&
      identical(rownames(value), expected_gene_ids) && identical(colnames(value), samples) &&
      all(is.finite(value)) && all(value >= 0) && all(value == floor(value)) &&
      all(value <= 2^53), "Exact ordered nonnegative integer count matrices required")
  }
  g205c_need(is.numeric(input_pairs) && !is.logical(input_pairs) &&
    identical(names(input_pairs), samples) && length(input_pairs) == 9L &&
    all(is.finite(input_pairs)) && all(input_pairs > 0) &&
    all(input_pairs == floor(input_pairs)) && all(input_pairs <= 2^53),
    "Nine named complete integer pair denominators required")
  left_sum <- colSums(left); right_sum <- colSums(right)
  g205c_need(all(left_sum <= input_pairs) && all(right_sum <= input_pairs),
    "Assigned fragment sums exceed registered pair denominators")
  summaries <- list(); ma <- list()
  for (j in seq_along(samples)) {
    old <- left[, j]; new <- right[, j]
    old_log <- log1p(old); new_log <- log1p(new)
    scopes <- list(ALL_ANNOTATED_GENES = rep(TRUE, length(old)),
      UNION_POSITIVE_GENES = old > 0 | new > 0)
    for (scope in names(scopes)) {
      keep <- scopes[[scope]]
      summaries[[length(summaries) + 1L]] <- data.frame(
        sample_id = samples[j], scope = scope, n_genes = sum(keep),
        left_assigned_fragments = left_sum[j], right_assigned_fragments = right_sum[j],
        input_pairs = input_pairs[j], left_assigned_fraction = left_sum[j] / input_pairs[j],
        right_assigned_fraction = right_sum[j] / input_pairs[j],
        left_zero_n = sum(old == 0), right_zero_n = sum(new == 0),
        both_zero_n = sum(old == 0 & new == 0),
        left_only_zero_n = sum(old == 0 & new > 0),
        right_only_zero_n = sum(new == 0 & old > 0),
        spearman_log1p = g205c_cor(old_log[keep], new_log[keep], "spearman"),
        pearson_log1p = g205c_cor(old_log[keep], new_log[keep], "pearson"),
        stringsAsFactors = FALSE)
    }
    ma[[j]] <- data.frame(gene_id = expected_gene_ids, sample_id = samples[j],
      left_count = old, right_count = new,
      A_log1p = (old_log + new_log) / 2, M_log1p_right_minus_left = new_log - old_log,
      stringsAsFactors = FALSE)
  }
  list(summary = do.call(rbind, summaries), MA = do.call(rbind, ma),
    count_scale = "UNNORMALIZED_ASSIGNED_GENE_FRAGMENTS",
    MA_scale = "NATURAL_LOG1P_COUNTS", zero_counts_scope = "ALL_ANNOTATED_GENES_IN_EACH_SUMMARY_ROW",
    production_32883_gene_geometry = length(expected_gene_ids) == 32883L,
    expected_gene_order_supplied_by_caller = TRUE, file_custody_authenticated_here = FALSE)
}

g205c_slots <- function() {
  specs <- rbind(data.frame(table = c("gene_edger", "gene_voom", "gene_deseq"), collection = "genes"),
    expand.grid(table = c("unsigned_camera", "unsigned_mroast", "unsigned_fgsea"),
      collection = c("Reactome", "SynGO_sensitivity"), stringsAsFactors = FALSE),
    data.frame(table = "activity_results", collection = c("PROGENy", "DoRothEA_AB")))
  out <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    data.frame(table = specs$table[i], collection = specs$collection[i],
      contrast = c("local_disease", "intervention"), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

g205c_key <- function(table, collection, contrast) {
  if (!length(collection)) return(character())
  paste(table, collection, contrast, sep = "::")
}

# Domain rows: table, collection, contrast, feature_id, eligible (logical),
# reason (nonempty source-grounded reason). They include not-eligible source
# features, not just returned tests. Missing eligible outputs require explicit
# native_not_returned=TRUE and are permitted only for native fgsea restriction.
# This cannot certify an omitted source feature without the caller's resource
# manifest. Never infer ineligibility from a missing result or change old masks.
g205c_validate_branch <- function(tables, domains) {
  slots <- g205c_slots(); roles <- unique(slots$table)
  g205c_need(is.list(tables) && setequal(names(tables), roles) &&
    !anyDuplicated(names(tables)) && all(vapply(tables, is.data.frame, logical(1))),
    "All seven complete native table roles required, including typed empty tables")
  required <- c("table", "collection", "contrast", "feature_id", "eligible", "reason", "native_not_returned")
  g205c_need(is.data.frame(domains) && !anyDuplicated(names(domains)) && all(required %in% names(domains)),
    "Explicit complete source domain ledger required")
  domains <- as.data.frame(domains)
  for (field in c("table", "collection", "contrast", "feature_id", "reason")) {
    g205c_need(is.character(domains[[field]]) && !anyNA(domains[[field]]) &&
      all(nzchar(domains[[field]])), "Domain strings must be explicit and nonmissing")
  }
  for (field in c("eligible", "native_not_returned")) {
    g205c_need(is.logical(domains[[field]]) && !anyNA(domains[[field]]), "Domain flags must be nonmissing logicals")
  }
  domain_slot <- g205c_key(domains$table, domains$collection, domains$contrast)
  valid_slots <- g205c_key(slots$table, slots$collection, slots$contrast)
  g205c_need(all(domain_slot %in% valid_slots) &&
    !anyDuplicated(domains[c("table", "collection", "contrast", "feature_id")]),
    "Unexpected or duplicate domain keys")
  for (role in roles) {
    tab <- tables[[role]]
    required <- c("contrast", "method", if (startsWith(role, "gene_")) "gene_id" else c("collection", "program_id"))
    g205c_need(!anyDuplicated(names(tab)) && all(required %in% names(tab)), "Native table key/method columns missing")
    for (field in required) g205c_need(is.character(tab[[field]]) && !anyNA(tab[[field]]) && all(nzchar(tab[[field]])),
      "Native table key/method values invalid")
    collection <- if (startsWith(role, "gene_")) rep("genes", nrow(tab)) else tab$collection
    g205c_need(all(g205c_key(role, collection, tab$contrast) %in% valid_slots), "Native collection/contrast outside fixed scope")
  }
  out <- vector("list", nrow(slots)); names(out) <- valid_slots
  for (i in seq_len(nrow(slots))) {
    spec <- slots[i, ]; tab <- as.data.frame(tables[[spec$table]])
    keep <- tab$contrast == spec$contrast
    if (spec$collection != "genes") keep <- keep & tab$collection == spec$collection
    tab <- as.data.frame(tab[which(keep), , drop = FALSE])
    ledger <- as.data.frame(domains[domain_slot == valid_slots[i], , drop = FALSE])
    ids <- if (spec$collection == "genes") tab$gene_id else tab$program_id
    g205c_need(g205c_ids(ids) && all(ids %in% ledger$feature_id[ledger$eligible]), "Returned native keys are duplicated or outside own eligible domain")
    present <- ledger$feature_id %in% ids
    g205c_need(identical(ledger$native_not_returned, ledger$eligible & !present) &&
      (spec$table == "unsigned_fgsea" || !any(ledger$native_not_returned)),
      "Missing eligible outputs require explicit unchanged-native fgsea restriction")
    out[[i]] <- g205c_project(tab, ledger, spec$table)
  }
  out
}

# Fixed source-native columns; no CAMERA/mroast directional proxy is invented.
# Full original tables (including all mid-P/Mixed/list/native aliases/old SE)
# remain separately returned, unmodified. Old edgeR covariance intervals are
# not adopted as validated QL uncertainty; standardized QL SE/CI remain NA.
g205c_project <- function(tab, ledger, role) {
  gene <- startsWith(role, "gene_")
  fields <- switch(role,
    gene_edger = c("logFC", "F", "PValue", "q_value"),
    gene_voom = c("logFC", "t", "p_value", "q_value", "se", "ci_low", "ci_high", "df_total"),
    gene_deseq = c("log2FoldChange", "stat", "p_value", "q_value", "lfcSE", "ci_low", "ci_high"),
    activity_results = c("effect", "t", "p_value", "q_within_collection_contrast", "se", "ci_low", "ci_high", "df_total"),
    unsigned_fgsea = c("NES", "p_value", "q_within_collection_contrast"),
    c("Direction", "p_value", "q_within_collection_contrast"))
  g205c_need(all(fields %in% names(tab)), "Required source-native fields absent, including typed-empty schema")
  expected_method <- switch(role, gene_edger = "edgeR_QL_robust", gene_voom = "limma_voom_robust_eBayes",
    gene_deseq = "DESeq2_Wald_sensitivity", activity_results = "decoupleR_ULM_then_limma",
    unsigned_camera = "camera_inter_gene_cor_0.01", unsigned_mroast = "mroast_mean_9999_rotations",
    unsigned_fgsea = "fgseaMultilevel")
  g205c_need(all(tab$method == expected_method), "Native method must match fixed source convention")
  for (field in setdiff(fields, "Direction")) g205c_need(is.numeric(tab[[field]]) && !is.logical(tab[[field]]), "Native numeric field has invalid type")
  p_field <- if (role == "gene_edger") "PValue" else "p_value"
  q_field <- if (gene) "q_value" else "q_within_collection_contrast"
  for (field in c(p_field, q_field)) g205c_need(all(is.na(tab[[field]]) |
    (is.finite(tab[[field]]) & tab[[field]] >= 0 & tab[[field]] <= 1)), "Native P/q outside probability domain")
  direction_only <- role %in% c("unsigned_camera", "unsigned_mroast")
  if (direction_only) g205c_need(is.character(tab$Direction) &&
    all(is.na(tab$Direction) | tab$Direction %in% c("Up", "Down")), "Native direction domain differs")
  ids <- if (gene) tab$gene_id else tab$program_id
  matched <- match(ledger$feature_id, ids)
  get <- function(field) if (is.null(field)) rep(NA_real_, nrow(ledger)) else tab[[field]][matched]
  effect_field <- switch(role, gene_edger = "logFC", gene_voom = "logFC",
    gene_deseq = "log2FoldChange", activity_results = "effect", unsigned_fgsea = "NES", NULL)
  stat_field <- switch(role, gene_edger = "F", gene_voom = "t", gene_deseq = "stat", activity_results = "t", NULL)
  uncertainty <- role %in% c("gene_voom", "gene_deseq", "activity_results")
  out <- ledger
  out$native_present <- !is.na(matched)
  out$native_p <- get(p_field); out$native_q <- get(q_field)
  out$effect_or_NES <- get(effect_field); out$native_statistic <- get(stat_field)
  out$native_direction <- if (direction_only) tab$Direction[matched] else rep(NA_character_, nrow(out))
  out$native_se <- if (uncertainty) get(if (role == "gene_deseq") "lfcSE" else "se") else rep(NA_real_, nrow(out))
  out$native_ci_low <- if (uncertainty) get("ci_low") else rep(NA_real_, nrow(out))
  out$native_ci_high <- if (uncertainty) get("ci_high") else rep(NA_real_, nrow(out))
  out$native_df <- if (uncertainty && role != "gene_deseq") get("df_total") else rep(NA_real_, nrow(out))
  out$uncertainty_status <- rep(if (uncertainty) "OWN_NATIVE_MODEL_ONLY_PRESERVE_NA" else if (role == "gene_edger")
    "NO_VALIDATED_NATIVE_QL_INTERVAL_OLD_SOURCE_COLUMNS_RETAINED_NOT_ADOPTED" else "NO_NATIVE_EFFECT_INTERVAL_DEFINED", nrow(out))
  out$continuous_quantity <- rep(if (direction_only) "NONE_CATEGORICAL_DIRECTION_ONLY" else if (role == "unsigned_fgsea")
    "NATIVE_NES_NOT_GENE_OR_ACTIVITY_EFFECT" else if (role == "activity_results") "ULM_SCORE_CONTRAST" else "LOG2_FOLD_CHANGE", nrow(out))
  out$status <- ifelse(!out$eligible, "NOT_ELIGIBLE_IN_OWN_BRANCH", ifelse(!out$native_present,
    "NATIVE_NOT_RETURNED", ifelse(direction_only, ifelse(is.na(out$native_direction),
      "RETURNED_DIRECTION_NA", "RETURNED_CATEGORICAL_DIRECTION"), ifelse(is.finite(out$effect_or_NES),
        "ESTIMABLE_CONTINUOUS_QUANTITY", "RETURNED_NONFINITE_CONTINUOUS_QUANTITY"))))
  out
}

# Top-k policy is opt-in. SAP/config name top-k but give no numeric k. This
# exact existing baseline P-top20 policy uses paired finite P values and all
# boundary ties, NOT |effect|, P transforms, or independently truncated sets.
g205c_topk <- function(left_p, right_p, ids, policy) {
  allowed <- c("PENDING_POLICY_BINDING", g205c_policy()$available_topk_policy)
  g205c_need(is.character(policy) && length(policy) == 1L && !is.na(policy) && policy %in% allowed,
    "Explicit supported top-k policy required; no arbitrary k")
  g205c_need(is.numeric(left_p) && is.numeric(right_p) && length(left_p) == length(ids) &&
    length(right_p) == length(ids) && g205c_ids(ids), "Paired top-k keys/lengths invalid")
  for (value in list(left_p, right_p)) g205c_need(all(is.na(value) |
    (is.finite(value) & value >= 0 & value <= 1)), "Top-k inputs must be native probabilities or NA")
  base <- list(policy = policy, paired_finite_P_n = sum(is.finite(left_p) & is.finite(right_p)),
    k = NA_integer_, left_n = NA_integer_, right_n = NA_integer_, intersection_n = NA_integer_,
    union_n = NA_integer_, Jaccard = NA_real_, left_ids = character(), right_ids = character())
  if (policy == "PENDING_POLICY_BINDING") return(base)
  keep <- is.finite(left_p) & is.finite(right_p)
  old <- left_p[keep]; new <- right_p[keep]; ids <- ids[keep]
  k <- min(20L, length(ids))
  if (!k) return(utils::modifyList(base, list(k = 0L, left_n = 0L, right_n = 0L,
    intersection_n = 0L, union_n = 0L)))
  a <- old <= sort(old)[k]; b <- new <= sort(new)[k]
  utils::modifyList(base, list(k = k, left_n = sum(a), right_n = sum(b),
    intersection_n = sum(a & b), union_n = sum(a | b), Jaccard = sum(a & b) / sum(a | b),
    left_ids = ids[a], right_ids = ids[b]))
}

# tables: exact seven-role lists from retained HISAT and new STAR. domains:
# complete own eligibility/source ledgers described above. topk_policies: NULL
# (pending all22), or named choices for explicitly bound slot keys only.
# No source native table is modified; no inference is performed for differences.
g205c_tables <- function(left_tables, right_tables, left_domains, right_domains,
                         topk_policies = NULL, caller_binding = NULL) {
  left <- g205c_validate_branch(left_tables, left_domains)
  right <- g205c_validate_branch(right_tables, right_domains)
  if (!is.null(topk_policies)) g205c_need(is.character(topk_policies) &&
    !anyNA(topk_policies) && g205c_ids(names(topk_policies)) && all(names(topk_policies) %in% names(left)),
    "Top-k policies must name only explicitly chosen fixed slot keys")
  slots <- g205c_slots(); rows <- list(); summaries <- list(); top <- list()
  for (i in seq_along(left)) {
    key <- names(left)[i]; old <- left[[i]]; new <- right[[i]]
    ids <- union(old$feature_id, new$feature_id)
    a <- old[match(ids, old$feature_id), , drop = FALSE]
    b <- new[match(ids, new$feature_id), , drop = FALSE]
    a_known <- ids %in% old$feature_id; b_known <- ids %in% new$feature_id
    continuous <- !slots$table[i] %in% c("unsigned_camera", "unsigned_mroast")
    paired <- a_known & b_known & a$native_present & b$native_present &
      is.finite(a$effect_or_NES) & is.finite(b$effect_or_NES)
    paired[is.na(paired)] <- FALSE
    categorical <- a_known & b_known & a$native_present & b$native_present &
      !is.na(a$native_direction) & !is.na(b$native_direction)
    categorical[is.na(categorical)] <- FALSE
    nonzero_union <- paired & (a$effect_or_NES != 0 | b$effect_or_NES != 0)
    nonzero_both <- paired & a$effect_or_NES != 0 & b$effect_or_NES != 0
    delta <- rep(NA_real_, length(ids)); delta[paired] <- b$effect_or_NES[paired] - a$effect_or_NES[paired]
    old_hit <- a_known & is.finite(a$native_q) & a$native_q < .05
    new_hit <- b_known & is.finite(b$native_q) & b$native_q < .05
    z <- data.frame(table = rep(slots$table[i], length(ids)), collection = rep(slots$collection[i], length(ids)),
      contrast = rep(slots$contrast[i], length(ids)),
      feature_id = ids, left_domain_present = a_known, right_domain_present = b_known,
      paired_estimable_continuous = paired, paired_categorical_direction = categorical,
      right_minus_left = delta, stringsAsFactors = FALSE)
    fields <- c("eligible", "reason", "native_not_returned", "native_present", "status", "native_p", "native_q",
      "effect_or_NES", "native_statistic", "native_direction", "native_se", "native_ci_low", "native_ci_high",
      "native_df", "uncertainty_status", "continuous_quantity")
    for (field in fields) {
      z[[paste0("left_", field)]] <- a[[field]]
      z[[paste0("right_", field)]] <- b[[field]]
    }
    z$union_status <- ifelse(!a_known, "RIGHT_DOMAIN_ONLY", ifelse(!b_known, "LEFT_DOMAIN_ONLY",
      ifelse(a$eligible & !b$eligible, "LEFT_ELIGIBLE_ONLY", ifelse(!a$eligible & b$eligible,
        "RIGHT_ELIGIBLE_ONLY", ifelse(paired, "BOTH_ESTIMABLE_CONTINUOUS", ifelse(categorical,
          "BOTH_CATEGORICAL_DIRECTION", "BOTH_DOMAIN_SEE_NATIVE_AND_ELIGIBILITY_FLAGS"))))))
    same <- sign(a$effect_or_NES) == sign(b$effect_or_NES)
    rate <- function(mask, value) if (sum(mask)) mean(value[mask]) else NA_real_
    summaries[[i]] <- data.frame(table = slots$table[i], collection = slots$collection[i], contrast = slots$contrast[i],
      union_domain_n = length(ids), left_domain_n = nrow(old), right_domain_n = nrow(new),
      left_eligible_n = sum(old$eligible), right_eligible_n = sum(new$eligible),
      left_returned_n = sum(old$native_present), right_returned_n = sum(new$native_present),
      left_not_returned_n = sum(old$native_not_returned), right_not_returned_n = sum(new$native_not_returned),
      left_q_available_n = sum(is.finite(old$native_q)), right_q_available_n = sum(is.finite(new$native_q)),
      left_P_available_n = sum(is.finite(old$native_p)), right_P_available_n = sum(is.finite(new$native_p)),
      left_q_NA_n = sum(is.na(old$native_q) & old$native_present), right_q_NA_n = sum(is.na(new$native_q) & new$native_present),
      left_q_strict_lt_05_n = sum(old_hit), right_q_strict_lt_05_n = sum(new_hit),
      left_q_exactly_05_n = sum(old$native_q == .05, na.rm = TRUE), right_q_exactly_05_n = sum(new$native_q == .05, na.rm = TRUE),
      q_hit_intersection_n = sum(old_hit & new_hit), q_hit_union_n = sum(old_hit | new_hit),
      q_hit_Jaccard = if (sum(old_hit | new_hit)) sum(old_hit & new_hit) / sum(old_hit | new_hit) else NA_real_,
      paired_estimable_continuous_n = sum(paired), nonzero_union_n = sum(nonzero_union),
      paired_P_available_n = sum(is.finite(a$native_p) & is.finite(b$native_p)),
      paired_q_available_n = sum(is.finite(a$native_q) & is.finite(b$native_q)),
      paired_own_CI_finite_n = sum(is.finite(a$native_ci_low) & is.finite(a$native_ci_high) &
        is.finite(b$native_ci_low) & is.finite(b$native_ci_high)),
      both_nonzero_n = sum(nonzero_both), same_sign_nonzero_union_n = sum(same[nonzero_union]),
      same_sign_both_nonzero_n = sum(same[nonzero_both]),
      sign_agreement_nonzero_union = rate(nonzero_union, same), sign_agreement_both_nonzero = rate(nonzero_both, same),
      both_zero_n = sum(paired & a$effect_or_NES == 0 & b$effect_or_NES == 0, na.rm = TRUE),
      categorical_paired_n = sum(categorical), categorical_same_n = sum(a$native_direction[categorical] == b$native_direction[categorical]),
      categorical_agreement = rate(categorical, a$native_direction == b$native_direction),
      effect_or_NES_spearman = if (continuous) g205c_cor(a$effect_or_NES[paired], b$effect_or_NES[paired], "spearman") else NA_real_,
      effect_or_NES_Kendall = if (continuous) g205c_cor(a$effect_or_NES[paired], b$effect_or_NES[paired], "kendall") else NA_real_,
      median_right_minus_left = if (any(paired)) stats::median(delta[paired]) else NA_real_,
      median_abs_difference = if (any(paired)) stats::median(abs(delta[paired])) else NA_real_,
      q95_abs_difference = if (any(paired)) as.numeric(stats::quantile(abs(delta[paired]), .95, names = FALSE)) else NA_real_,
      stringsAsFactors = FALSE)
    policy <- if (!is.null(topk_policies) && key %in% names(topk_policies)) topk_policies[[key]] else "PENDING_POLICY_BINDING"
    top[[key]] <- g205c_topk(a$native_p, b$native_p, ids, policy)
    rows[[i]] <- z
  }
  list(schema = "gse205906-matched-saved-table-comparison-v1", summary = do.call(rbind, summaries),
    union = do.call(rbind, rows), topk = top, native_tables = list(left = left_tables, right = right_tables),
    policy = g205c_policy(), caller_binding = caller_binding,
    source_domain_completeness_authenticated_here = FALSE, real_comparison_executed_by_sourcing = FALSE)
}
