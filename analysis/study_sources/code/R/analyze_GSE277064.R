#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(file_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(yaml)
  library(decoupleR)
  library(fgsea)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

analysis_id <- "GSE277064_route_b2_edgeR_QL_v1.0.0"
random_seed <- 20260817L
set.seed(random_seed)
roast_nrot_run <- as.integer(Sys.getenv("GSE277064_ROAST_NROT", unset = "9999"))
run_fgsea_flag <- !identical(toupper(Sys.getenv("GSE277064_RUN_FGSEA", unset = "TRUE")), "FALSE")
assert_roast_nrot <- is.finite(roast_nrot_run) && roast_nrot_run >= 99L
if (!assert_roast_nrot) stop("GSE277064_ROAST_NROT must be an integer >=99", call. = FALSE)

raw_dir <- file.path(project_root, "data", "raw", "GSE277064")
gene_set_dir <- file.path(project_root, "data", "processed", "gene_sets")
metadata_dir <- file.path(project_root, "metadata")
result_root <- file.path(project_root, "results", "GSE277064")
dir_names <- c(
  "audit", "diagnostics", "figures", "gene_level", "objects", "programs",
  file.path("boundary_rotarod", "diagnostics"),
  file.path("boundary_rotarod", "figures"),
  file.path("boundary_rotarod", "gene_level"),
  file.path("boundary_rotarod", "objects"),
  file.path("boundary_rotarod", "programs")
)
for (d in dir_names) dir.create(file.path(result_root, d), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_root, "diagnostics", "analysis.log")
log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

cat("Analysis:", analysis_id, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

write_optional_program_tsv <- function(x, path, not_run_reason) {
  if (ncol(as.data.table(x)) == 0L) {
    write_tsv(data.table(status = "NOT_RUN", reason = not_run_reason), path)
  } else {
    write_tsv(x, path)
  }
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

safe_scale_rows <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1L, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1
  sweep(sweep(x, 1L, means, "-"), 1L, sds, "/")
}

break_rank_stat_ties <- function(stats) {
  stats <- stats[is.finite(stats)]
  original_ties <- sum(duplicated(stats))
  iteration <- 0L
  while (anyDuplicated(stats) && iteration < 12L) {
    iteration <- iteration + 1L
    tie_groups <- split(seq_along(stats), stats)
    tie_groups <- tie_groups[lengths(tie_groups) > 1L]
    for (indices in tie_groups) {
      ordered_indices <- indices[order(names(stats)[indices])]
      base_value <- stats[ordered_indices[1L]]
      step <- (2 ^ iteration) * .Machine$double.eps * max(1, abs(base_value))
      centered_rank <- seq_along(ordered_indices) - mean(seq_along(ordered_indices))
      stats[ordered_indices] <- base_value + centered_rank * step
    }
  }
  assert_true(!anyDuplicated(stats), "Unable to deterministically resolve preranked-statistic ties")
  attr(stats, "original_ties") <- original_ties
  stats
}

robust_z <- function(x) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - center) / spread
}

collapse_matrix_by_symbol <- function(x, symbols) {
  valid <- !is.na(symbols) & nzchar(symbols)
  x <- x[valid, , drop = FALSE]
  symbols <- symbols[valid]
  sums <- rowsum(x, group = symbols, reorder = FALSE, na.rm = TRUE)
  denominators <- as.numeric(table(factor(symbols, levels = rownames(sums))))
  sums / denominators
}

make_activity_matrix <- function(long_result, sample_ids) {
  dt <- as.data.table(long_result)
  assert_true(all(c("source", "condition", "score") %in% names(dt)), "Unexpected decoupleR output schema")
  wide <- dcast(dt, source ~ condition, value.var = "score")
  programs <- wide$source
  wide[, source := NULL]
  matrix_result <- as.matrix(wide)
  rownames(matrix_result) <- programs
  missing_samples <- setdiff(sample_ids, colnames(matrix_result))
  assert_true(length(missing_samples) == 0L, paste("Missing activity samples:", paste(missing_samples, collapse = ", ")))
  matrix_result[, sample_ids, drop = FALSE]
}

fit_activity_contrasts <- function(score_matrix, design, contrast_matrix, collection, method) {
  assert_true(identical(colnames(score_matrix), rownames(design)), "Activity matrix/design sample mismatch")
  fit <- lmFit(score_matrix, design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit, robust = TRUE)
  output <- rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(j) {
    se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
    df_total <- fit$df.total
    critical <- qt(0.975, df = df_total)
    data.table(
      collection = collection,
      method = method,
      program_id = rownames(score_matrix),
      contrast = colnames(contrast_matrix)[j],
      effect = fit$coefficients[, j],
      se = se,
      ci_low = fit$coefficients[, j] - critical * se,
      ci_high = fit$coefficients[, j] + critical * se,
      t = fit$t[, j],
      df_total = df_total,
      p_value = fit$p.value[, j]
    )
  }))
  output[, q_within_celltype_contrast := p.adjust(p_value, method = "BH"), by = .(collection, contrast)]
  output
}

fit_edgeR_contrasts <- function(fit, contrast_matrix, annotation, cell_type, domain) {
  rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(j) {
    qlf <- glmQLFTest(fit, contrast = contrast_matrix[, j])
    top_table <- topTags(qlf, n = Inf, sort.by = "none")$table
    if ("gene_id_versioned" %in% colnames(top_table)) {
      tab <- as.data.table(top_table)
    } else {
      tab <- as.data.table(top_table, keep.rownames = "gene_id_versioned")
    }
    tab <- tab[match(rownames(fit$counts), gene_id_versioned)]
    df_total <- qlf$df.total
    if (length(df_total) == 1L) df_total <- rep(df_total, nrow(tab))
    se <- ifelse(is.finite(tab$F) & tab$F > 0, abs(tab$logFC) / sqrt(tab$F), NA_real_)
    critical <- qt(0.975, df = df_total)
    tab[, signed_sqrt_F := sign(logFC) * sqrt(F)]
    tab[, ql_se := se]
    tab[, ci_low := logFC - critical * se]
    tab[, ci_high := logFC + critical * se]
    tab[, df_total := df_total]
    tab[, contrast := colnames(contrast_matrix)[j]]
    tab[, cell_type := cell_type]
    tab[, domain := domain]
    for (annotation_column in setdiff(names(annotation), names(tab))) {
      tab[, (annotation_column) := annotation[[annotation_column]][match(gene_id_versioned, annotation$gene_id_versioned)]]
    }
    tab
  }))
}

make_reactome_resources <- function(annotation, reactome, orthology, eligible_gene_ids) {
  annotation <- as.data.table(annotation)
  mouse_map <- orthology[other_tax_id == 10090, .(human_symbol, mouse_symbol = other_symbol)]
  mouse_map <- unique(mouse_map[!is.na(human_symbol) & !is.na(mouse_symbol)])
  mapped <- merge(
    reactome[, .(term_id, term_name, human_symbol = gene_symbol)],
    mouse_map,
    by = "human_symbol",
    allow.cartesian = TRUE
  )
  mapped <- unique(mapped[, .(term_id, term_name, mouse_symbol)])
  feature_table <- annotation[gene_id_versioned %in% eligible_gene_ids & !is.na(symbol) & nzchar(symbol)]
  indices <- split(feature_table$gene_id_versioned, feature_table$symbol)
  term_members <- split(mapped$mouse_symbol, mapped$term_id)
  term_names <- unique(mapped[, .(term_id, term_name)])
  index_list <- lapply(term_members, function(symbols) {
    unique(unlist(indices[intersect(names(indices), unique(symbols))], use.names = FALSE))
  })
  index_list <- lapply(index_list, function(ids) match(ids, eligible_gene_ids))
  sizes <- lengths(index_list)
  eligible <- sizes >= 10L & sizes <= 500L
  index_list <- index_list[eligible]
  term_names <- term_names[match(names(index_list), term_id)]
  mapping <- merge(
    reactome[, .(locked_human_n = uniqueN(gene_symbol)), by = .(term_id, term_name)],
    mapped[, .(mapped_mouse_n = uniqueN(mouse_symbol)), by = .(term_id, term_name)],
    by = c("term_id", "term_name"),
    all.x = TRUE
  )
  observed_sizes <- data.table(term_id = names(index_list), tested_feature_n = lengths(index_list))
  mapping <- merge(mapping, observed_sizes, by = "term_id", all.x = TRUE)
  mapping[, eligible_10_500 := !is.na(tested_feature_n) & tested_feature_n >= 10L & tested_feature_n <= 500L]
  list(indices = index_list, names = term_names, mapping = mapping, mapped_members = mapped)
}

make_reactome_scores <- function(log_expression, index_list) {
  z <- safe_scale_rows(log_expression)
  scores <- t(vapply(index_list, function(index) {
    colMeans(z[index, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(z))))
  rownames(scores) <- names(index_list)
  colnames(scores) <- colnames(log_expression)
  safe_scale_rows(scores)
}

run_reactome_tests <- function(v, design, contrast_matrix, reactome_resource, gene_stats,
                               roast_nrot, run_roast = TRUE, run_fgsea = TRUE,
                               sensitivity_contrasts = colnames(contrast_matrix)) {
  index_list <- reactome_resource$indices
  name_map <- reactome_resource$names
  camera_results <- list()
  roast_results <- list()
  fgsea_results <- list()

  for (j in seq_len(ncol(contrast_matrix))) {
    contrast_name <- colnames(contrast_matrix)[j]
    contrast_vector <- contrast_matrix[, j]
    camera_table <- as.data.table(camera(
      v,
      index = index_list,
      design = design,
      contrast = contrast_vector,
      inter.gene.cor = 0.01,
      sort = FALSE
    ), keep.rownames = "program_id")
    setnames(camera_table, c("NGenes", "Direction", "PValue", "FDR"),
             c("tested_feature_n", "direction", "p_value", "q_within_celltype_contrast"))
    camera_table[, contrast := contrast_name]
    camera_table[, method := "camera"]
    camera_results[[contrast_name]] <- camera_table

    if (run_roast && contrast_name %in% sensitivity_contrasts) {
      set.seed(random_seed + j)
      roast_table <- as.data.table(mroast(
        v,
        index = index_list,
        design = design,
        contrast = contrast_vector,
        set.statistic = "mean",
        nrot = roast_nrot,
        approx.zscore = TRUE,
        adjust.method = "BH",
        midp = TRUE,
        sort = "none"
      ), keep.rownames = "program_id")
      setnames(roast_table, c("NGenes", "Direction", "PValue", "FDR", "PValue.Mixed", "FDR.Mixed"),
               c("tested_feature_n", "direction", "p_value", "q_within_celltype_contrast",
                 "p_value_mixed", "q_mixed_within_celltype_contrast"))
      roast_table[, contrast := contrast_name]
      roast_table[, method := "mroast"]
      roast_results[[contrast_name]] <- roast_table
    }

    if (run_fgsea && contrast_name %in% sensitivity_contrasts) {
      stats_table <- gene_stats[contrast == contrast_name & !is.na(symbol) & nzchar(symbol),
                                .(rank_stat = mean(signed_sqrt_F, na.rm = TRUE)), by = symbol]
      stats_vector <- stats_table$rank_stat
      names(stats_vector) <- stats_table$symbol
      stats_vector <- break_rank_stat_ties(stats_vector)
      stats_vector <- sort(stats_vector[is.finite(stats_vector)], decreasing = TRUE)
      pathway_symbols <- lapply(index_list, function(index) {
        unique(na.omit(v$genes$symbol[index]))
      })
      set.seed(random_seed + 100L + j)
      fgsea_table <- as.data.table(fgseaMultilevel(
        pathways = pathway_symbols,
        stats = stats_vector,
        minSize = 10L,
        maxSize = 500L,
        eps = 0,
        nproc = 1
      ))
      fgsea_table[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
      setnames(fgsea_table, c("pathway", "pval", "padj", "size"),
               c("program_id", "p_value", "q_within_celltype_contrast", "tested_feature_n"))
      fgsea_table[, contrast := contrast_name]
      fgsea_table[, method := "fgseaMultilevel"]
      fgsea_results[[contrast_name]] <- fgsea_table
    }
  }

  add_names <- function(x) {
    if (length(x) == 0L) return(data.table())
    ans <- rbindlist(x, fill = TRUE)
    ans <- merge(ans, name_map, by.x = "program_id", by.y = "term_id", all.x = TRUE, sort = FALSE)
    ans[, collection := "Reactome"]
    ans
  }
  list(
    camera = add_names(camera_results),
    roast = add_names(roast_results),
    fgsea = add_names(fgsea_results)
  )
}

fit_loso_activity <- function(score_matrix, sample_meta, design_column, levels, contrast_vector,
                              collection, cell_type, domain) {
  full_design <- model.matrix(~ 0 + factor(sample_meta[[design_column]], levels = levels))
  colnames(full_design) <- levels
  full_fit <- eBayes(contrasts.fit(lmFit(score_matrix, full_design), contrast_vector), robust = TRUE)
  full_effect <- full_fit$coefficients[, 1L]
  leave_effects <- matrix(NA_real_, nrow = nrow(score_matrix), ncol = ncol(score_matrix),
                          dimnames = list(rownames(score_matrix), colnames(score_matrix)))
  for (i in seq_len(ncol(score_matrix))) {
    kept <- seq_len(ncol(score_matrix)) != i
    design_i <- model.matrix(~ 0 + factor(sample_meta[[design_column]][kept], levels = levels))
    colnames(design_i) <- levels
    assert_true(qr(design_i)$rank == ncol(design_i), "LOSO design became rank deficient")
    fit_i <- eBayes(contrasts.fit(lmFit(score_matrix[, kept, drop = FALSE], design_i), contrast_vector), robust = TRUE)
    leave_effects[, i] <- fit_i$coefficients[, 1L]
  }
  data.table(
    collection = collection,
    cell_type = cell_type,
    domain = domain,
    program_id = rownames(score_matrix),
    full_effect = full_effect,
    leave_one_min = apply(leave_effects, 1L, min, na.rm = TRUE),
    leave_one_max = apply(leave_effects, 1L, max, na.rm = TRUE),
    direction_fraction = rowMeans(sign(leave_effects) == sign(full_effect), na.rm = TRUE),
    sign_flip_any = rowSums(sign(leave_effects) != sign(full_effect), na.rm = TRUE) > 0,
    max_abs_delta = apply(abs(leave_effects - full_effect), 1L, max, na.rm = TRUE)
  )
}

plot_sample_qc <- function(qc_table, pca_table, cell_type, domain, output_dir) {
  palette <- c(control = "#3B6FB6", D3 = "#E69F00", D7 = "#D55E00", D14 = "#009E73", D28 = "#7A5195",
               D14_untrained = "#3B6FB6", H2 = "#E69F00", H4 = "#D55E00", H8 = "#009E73")
  p1 <- ggplot(qc_table, aes(x = group, y = library_size, color = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) + geom_point(size = 2.5) +
    scale_color_manual(values = palette) + theme_bw(base_size = 10) +
    labs(title = paste(cell_type, domain, "library size"), x = NULL, y = "Raw library size") +
    theme(legend.position = "none")
  ggsave(file.path(output_dir, paste0(cell_type, "_library_size.png")), p1, width = 6.5, height = 4.5, dpi = 180)

  p2 <- ggplot(pca_table, aes(x = PC1, y = PC2, color = group, label = source_column_or_member)) +
    geom_hline(yintercept = 0, color = "grey90") + geom_vline(xintercept = 0, color = "grey90") +
    geom_point(size = 3) + geom_text(size = 2.5, vjust = -0.8, check_overlap = TRUE) +
    scale_color_manual(values = palette) + theme_bw(base_size = 10) +
    labs(title = paste(cell_type, domain, "PCA"))
  ggsave(file.path(output_dir, paste0(cell_type, "_PCA.png")), p2, width = 7.5, height = 5.5, dpi = 180)
}

save_correlation_plot <- function(correlation_matrix, sample_meta, path, title) {
  pdf(path, width = 8, height = 7)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); dev.off()}, add = TRUE)
  par(mar = c(8, 8, 3, 2))
  image(seq_len(nrow(correlation_matrix)), seq_len(ncol(correlation_matrix)),
        t(correlation_matrix[nrow(correlation_matrix):1, , drop = FALSE]),
        axes = FALSE, col = hcl.colors(100, "Blue-Red 3"), zlim = c(min(correlation_matrix), 1),
        main = title, xlab = "", ylab = "")
  axis(1, at = seq_len(ncol(correlation_matrix)), labels = colnames(correlation_matrix), las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(nrow(correlation_matrix)), labels = rev(rownames(correlation_matrix)), las = 2, cex.axis = 0.65)
}

natural_contrast_matrix <- function(design) {
  expected <- c("control", "D3", "D7", "D14", "D28")
  assert_true(identical(colnames(design), expected), "Natural-course design levels differ from lock")
  makeContrasts(
    chronic_combined = 0.5 * D14 + 0.5 * D28 - control,
    D14_vs_control = D14 - control,
    D28_vs_control = D28 - control,
    early_combined = 0.5 * D3 + 0.5 * D7 - control,
    D3_vs_control = D3 - control,
    D7_vs_control = D7 - control,
    levels = design
  )
}

exercise_contrast_matrix <- function(design) {
  expected <- c("D14_untrained", "H2", "H4", "H8")
  assert_true(identical(colnames(design), expected), "Rotarod-boundary design levels differ from lock")
  makeContrasts(
    H2_vs_D14_untrained = H2 - D14_untrained,
    H4_vs_D14_untrained = H4 - D14_untrained,
    H8_vs_D14_untrained = H8 - D14_untrained,
    levels = design
  )
}

## -------------------------------------------------------------------------
## Frozen configuration and input identity audit
## -------------------------------------------------------------------------
contrasts_cfg <- read_yaml(file.path(project_root, "config", "contrasts.yaml"))
gene_sets_cfg <- read_yaml(file.path(project_root, "config", "gene_sets.yaml"))
thresholds_cfg <- read_yaml(file.path(project_root, "config", "thresholds.yaml"))
assert_true(identical(as.integer(contrasts_cfg$global$random_seed), random_seed), "Random seed differs from locked configuration")
assert_true(identical(gene_sets_cfg$freeze_status, "frozen_before_outcome_inference"), "Gene-set freeze is not locked")
assert_true(identical(as.integer(gene_sets_cfg$primary_collections$reactome$release), 97L), "Reactome release mismatch")
assert_true(identical(as.numeric(thresholds_cfg$expression_filtering$min_cpm), 1), "CPM threshold mismatch")

ledger <- fread(file.path(metadata_dir, "samples.tsv"), sep = "\t", quote = "", na.strings = c("", "NA"))
ledger <- ledger[dataset_id == "GSE277064"]
assert_true(nrow(ledger) == 63L, "Frozen sample ledger must contain 63 GSE277064 rows")
assert_true(uniqueN(ledger$sample_id) == 63L, "GSM identifiers are not unique")
assert_true(uniqueN(ledger$source_column_or_member) == 63L, "Count-matrix columns are not unique")
setorder(ledger, sample_order)

official_long <- fread(file.path(metadata_dir, "geo_series_metadata_long.tsv"), sep = "\t", quote = "")
official_titles <- official_long[dataset_id == "GSE277064" & field == "Sample_title",
                                 .(sample_id, official_title = value)]
assert_true(nrow(official_titles) == 63L && uniqueN(official_titles$sample_id) == 63L,
            "Official GEO title map is incomplete")

count_files <- c(
  v2a = file.path(raw_dir, "GSE277064_Chx10.txt.gz"),
  astrocyte = file.path(raw_dir, "GSE277064_Aldh1l1.txt.gz"),
  microglia = file.path(raw_dir, "GSE277064_Cx3cr1.txt.gz")
)
read_count_table <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  setnames(dt, 1L, "gene_id_versioned")
  assert_true(!anyDuplicated(dt$gene_id_versioned), paste("Duplicate Ensembl IDs in", basename(path)))
  mat <- as.matrix(dt[, -"gene_id_versioned"])
  storage.mode(mat) <- "numeric"
  # RSEM expected counts are valid non-negative count estimates and can be
  # fractional because multi-mapping reads are allocated probabilistically.
  # They must not be rounded before edgeR modeling.
  assert_true(all(is.finite(mat)) && all(mat >= 0), paste("Counts are not finite and non-negative in", basename(path)))
  rownames(mat) <- dt$gene_id_versioned
  storage.mode(mat) <- "double"
  mat
}
count_list <- lapply(count_files, read_count_table)
assert_true(all(vapply(count_list[-1], function(x) identical(rownames(x), rownames(count_list[[1]])), logical(1))),
            "Gene rows differ across official count matrices")
counts_all <- do.call(cbind, count_list)
assert_true(ncol(counts_all) == 63L && uniqueN(colnames(counts_all)) == 63L, "Official matrices do not contain 63 unique columns")
assert_true(setequal(colnames(counts_all), ledger$source_column_or_member), "Count columns and frozen ledger differ")

ledger <- merge(ledger, official_titles, by = "sample_id", all.x = TRUE, sort = FALSE)
setorder(ledger, sample_order)
cell_label <- c(v2a = "V2a", astrocyte = "Astrocyte", microglia = "Microglia")
ledger[, expected_official_title := ifelse(
  comparator_set %like% "_exercise$",
  sprintf("%s, %s h, rep%s", cell_label[cell_type], exercise_hour, replicate),
  ifelse(
    time_day == 0,
    sprintf("%s, Control, rep%s", cell_label[cell_type], replicate),
    sprintf("%s, Day %s, rep%s", cell_label[cell_type], time_day, replicate)
  )
)]
ledger[, count_column_present := source_column_or_member %chin% colnames(counts_all)]
ledger[, official_title_exact_match := official_title == expected_official_title]
assert_true(all(ledger$count_column_present), "At least one frozen count column is missing")
assert_true(all(ledger$official_title_exact_match), "Frozen GSM/count-column mapping disagrees with official GEO titles")

natural_ledger <- ledger[comparator_set %like% "_natural$"]
exercise_ledger <- ledger[comparator_set %like% "_exercise$"]
assert_true(nrow(natural_ledger) == 45L, "Natural-course population must contain exactly 45 samples")
assert_true(nrow(exercise_ledger) == 18L, "Rotarod-boundary population must contain exactly 18 trained samples")
natural_counts_check <- natural_ledger[, .N, by = .(cell_type, time_day)]
assert_true(nrow(natural_counts_check) == 15L && all(natural_counts_check$N == 3L), "Natural-course groups must all have n=3")
exercise_counts_check <- exercise_ledger[, .N, by = .(cell_type, exercise_hour)]
assert_true(nrow(exercise_counts_check) == 6L && all(exercise_counts_check$N == 3L), "Rotarod groups must all have n=3")
assert_true(!any(exercise_ledger$source_column_or_member %chin% natural_ledger$source_column_or_member),
            "Rotarod samples leaked into the natural-course population")

mapping_audit <- ledger[, .(
  sample_order, sample_id, source_column_or_member, biological_unit_id, cell_type,
  condition, time_day, exercise_hour, replicate, comparator_set, official_title,
  expected_official_title, count_column_present, official_title_exact_match
)]
write_tsv(mapping_audit, file.path(result_root, "audit", "sample_count_column_GSM_mapping.tsv"))
write_tsv(natural_counts_check, file.path(result_root, "audit", "natural_course_group_counts.tsv"))
write_tsv(exercise_counts_check, file.path(result_root, "audit", "rotarod_trained_group_counts.tsv"))

gene_ids <- rownames(counts_all)
gene_id_base <- sub("\\..*$", "", gene_ids)
gene_symbols <- unname(mapIds(org.Mm.eg.db, keys = gene_id_base, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
gene_entrez <- unname(mapIds(org.Mm.eg.db, keys = gene_id_base, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first"))
gene_names <- unname(mapIds(org.Mm.eg.db, keys = gene_id_base, keytype = "ENSEMBL", column = "GENENAME", multiVals = "first"))
annotation <- data.table(
  gene_id_versioned = gene_ids,
  gene_id_base = gene_id_base,
  symbol = as.character(gene_symbols),
  entrez_id = as.character(gene_entrez),
  gene_name = as.character(gene_names)
)
write_tsv(annotation, file.path(result_root, "audit", "ensembl_to_mouse_annotation.tsv"))

gene_resources <- readRDS(file.path(gene_set_dir, "gene_set_resources_frozen_2026-08-18.rds"))
reactome <- as.data.table(gene_resources$reactome)
progeny_mouse <- as.data.table(gene_resources$progeny_mouse)
dorothea_mouse <- as.data.table(gene_resources$dorothea_mouse)
orthology <- as.data.table(gene_resources$ortholog_one_to_one)
assert_true(unique(reactome$release) == "97", "Frozen Reactome data release is not 97")
assert_true(unique(progeny_mouse$release) == "1.34.0", "Frozen mouse PROGENy release mismatch")
assert_true(unique(dorothea_mouse$release) == "1.23.0", "Frozen mouse DoRothEA release mismatch")

## -------------------------------------------------------------------------
## Per-domain model runner. The natural course and rotarod boundary are never
## combined in a count model or multiplicity family.
## -------------------------------------------------------------------------
run_domain_model <- function(cell_type, domain = c("natural_course", "rotarod_boundary"), roast_nrot = 9999L) {
  domain <- match.arg(domain)
  cell_type_name <- cell_type
  if (domain == "natural_course") {
    sample_meta <- copy(natural_ledger[cell_type == cell_type_name])
    setorder(sample_meta, sample_order)
    assert_true(nrow(sample_meta) == 15L, paste(cell_type, "natural course does not contain 15 samples"))
    sample_meta[, group := fifelse(time_day == 0, "control", paste0("D", time_day))]
    levels <- c("control", "D3", "D7", "D14", "D28")
    output_base <- result_root
  } else {
    assert_true(cell_type %chin% c("v2a", "astrocyte"), "Rotarod boundary is available only for V2a and astrocytes")
    untrained <- copy(natural_ledger[cell_type == cell_type_name & time_day == 14])
    trained <- copy(exercise_ledger[cell_type == cell_type_name])
    sample_meta <- rbindlist(list(untrained, trained), fill = TRUE)
    setorder(sample_meta, sample_order)
    assert_true(nrow(sample_meta) == 12L, paste(cell_type, "rotarod boundary must contain 3 untrained and 9 trained samples"))
    sample_meta[, group := fifelse(is.na(exercise_hour), "D14_untrained", paste0("H", exercise_hour))]
    levels <- c("D14_untrained", "H2", "H4", "H8")
    output_base <- file.path(result_root, "boundary_rotarod")
  }

  sample_meta[, group := factor(group, levels = levels)]
  assert_true(all(table(sample_meta$group) == 3L), paste(cell_type, domain, "does not have n=3 per group"))
  sample_columns <- sample_meta$source_column_or_member
  counts <- counts_all[, sample_columns, drop = FALSE]
  design <- model.matrix(~ 0 + group, data = sample_meta)
  colnames(design) <- levels
  rownames(design) <- sample_columns
  assert_true(qr(design)$rank == ncol(design), paste(cell_type, domain, "design is rank deficient"))
  contrast_matrix <- if (domain == "natural_course") natural_contrast_matrix(design) else exercise_contrast_matrix(design)

  y_unfiltered <- DGEList(counts = counts, genes = annotation)
  keep_fbe <- filterByExpr(y_unfiltered, design = design)
  min_group_n <- min(table(sample_meta$group))
  keep_cpm <- rowSums(cpm(y_unfiltered) >= 1) >= min_group_n
  keep <- keep_fbe & keep_cpm
  assert_true(sum(keep) >= 5000L, paste(cell_type, domain, "retained unexpectedly few genes"))
  y <- y_unfiltered[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  v <- voom(y, design = design, plot = FALSE, normalize.method = "none")
  v$genes <- annotation[match(rownames(v$E), gene_id_versioned)]

  log_cpm <- cpm(y, log = TRUE, prior.count = 2)
  sample_cor <- cor(log_cpm, method = "pearson")
  pca <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
  pca_table <- data.table(
    source_column_or_member = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L]
  )
  pca_table <- merge(pca_table, sample_meta[, .(source_column_or_member, sample_id, group)],
                     by = "source_column_or_member", all.x = TRUE, sort = FALSE)
  pca_table[, group := as.character(group)]
  qc_table <- data.table(
    source_column_or_member = sample_columns,
    sample_id = sample_meta$sample_id,
    biological_unit_id = sample_meta$biological_unit_id,
    cell_type = cell_type,
    domain = domain,
    group = as.character(sample_meta$group),
    library_size = colSums(counts),
    detected_cpm1 = colSums(cpm(y_unfiltered) >= 1),
    TMM_norm_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    median_sample_correlation = apply(sample_cor, 2L, function(x) median(x[x < 0.999999], na.rm = TRUE))
  )
  qc_table[, library_robust_z := robust_z(log2(library_size))]
  qc_table[, detected_robust_z := robust_z(detected_cpm1)]
  qc_table[, correlation_robust_z := robust_z(median_sample_correlation)]
  qc_table[, objective_metric_flag := abs(library_robust_z) > 3.5 | abs(detected_robust_z) > 3.5 | correlation_robust_z < -3.5]
  qc_table[, exclusion := FALSE]
  qc_table[, exclusion_reason := "none; PCA position alone is not an exclusion criterion"]

  filter_summary <- data.table(
    cell_type = cell_type,
    domain = domain,
    n_samples = ncol(counts),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    residual_df = nrow(design) - qr(design)$rank,
    input_features = nrow(y_unfiltered),
    filterByExpr_features = sum(keep_fbe),
    cpm1_in_at_least_3_features = sum(keep_cpm),
    retained_intersection_features = sum(keep),
    retained_annotated_symbols = sum(!is.na(y$genes$symbol) & nzchar(y$genes$symbol)),
    common_dispersion = y$common.dispersion,
    biological_coefficient_of_variation = sqrt(y$common.dispersion),
    median_ql_prior_df = median(fit$df.prior, na.rm = TRUE),
    design_condition_number = kappa(design)
  )

  write_tsv(qc_table, file.path(output_base, "diagnostics", paste0(cell_type, "_sample_QC.tsv")))
  write_tsv(pca_table, file.path(output_base, "diagnostics", paste0(cell_type, "_PCA_coordinates.tsv")))
  write_tsv(as.data.table(sample_cor, keep.rownames = "sample"),
            file.path(output_base, "diagnostics", paste0(cell_type, "_sample_correlations.tsv")))
  write_tsv(filter_summary, file.path(output_base, "diagnostics", paste0(cell_type, "_model_summary.tsv")))
  write_tsv(as.data.table(design, keep.rownames = "sample"),
            file.path(output_base, "diagnostics", paste0(cell_type, "_design_matrix.tsv")))
  write_tsv(as.data.table(contrast_matrix, keep.rownames = "coefficient"),
            file.path(output_base, "diagnostics", paste0(cell_type, "_contrast_matrix.tsv")))

  plot_sample_qc(qc_table, pca_table, cell_type, domain, file.path(output_base, "figures"))
  save_correlation_plot(sample_cor, sample_meta,
                        file.path(output_base, "figures", paste0(cell_type, "_sample_correlation.pdf")),
                        paste(cell_type, domain, "sample correlation"))
  pdf(file.path(output_base, "figures", paste0(cell_type, "_BCV.pdf")), width = 6, height = 5)
  plotBCV(y, main = paste(cell_type, domain, "BCV"))
  dev.off()
  pdf(file.path(output_base, "figures", paste0(cell_type, "_QL_dispersion.pdf")), width = 6, height = 5)
  plotQLDisp(fit, main = paste(cell_type, domain, "QL dispersion"))
  dev.off()

  gene_stats <- fit_edgeR_contrasts(fit, contrast_matrix, annotation, cell_type, domain)
  for (contrast_name in unique(gene_stats$contrast)) {
    out <- gene_stats[contrast == contrast_name]
    write_tsv(out, file.path(output_base, "gene_level", cell_type, paste0(contrast_name, ".tsv.gz")))
  }
  gene_summary <- gene_stats[, .(
    tested_features = .N,
    FDR_0_05 = sum(FDR <= 0.05, na.rm = TRUE),
    FDR_0_05_up = sum(FDR <= 0.05 & logFC > 0, na.rm = TRUE),
    FDR_0_05_down = sum(FDR <= 0.05 & logFC < 0, na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ), by = .(cell_type, domain, contrast)]
  write_tsv(gene_summary, file.path(output_base, "gene_level", paste0(cell_type, "_contrast_summary.tsv")))

  reactome_resource <- make_reactome_resources(
    annotation = y$genes,
    reactome = reactome,
    orthology = orthology,
    eligible_gene_ids = rownames(y$counts)
  )
  write_tsv(reactome_resource$mapping,
            file.path(output_base, "programs", paste0(cell_type, "_Reactome_mapping_coverage.tsv")))
  reactome_scores <- make_reactome_scores(log_cpm, reactome_resource$indices)
  reactome_activity <- fit_activity_contrasts(
    reactome_scores, design, contrast_matrix, "Reactome", "standardized_mean_logCPM_activity"
  )
  reactome_tests <- run_reactome_tests(
    v, design, contrast_matrix, reactome_resource, gene_stats,
    roast_nrot = roast_nrot,
    run_roast = domain == "natural_course",
    run_fgsea = run_fgsea_flag && domain == "natural_course",
    sensitivity_contrasts = if (domain == "natural_course") "chronic_combined" else character()
  )

  expression_by_symbol <- collapse_matrix_by_symbol(v$E, v$genes$symbol)
  progeny_network <- progeny_mouse[, .(
    source = as.character(pathway),
    target = as.character(gene_symbol),
    mor = as.numeric(weight)
  )]
  progeny_long <- run_ulm(
    mat = expression_by_symbol,
    network = progeny_network,
    .source = source,
    .target = target,
    .mor = mor,
    minsize = 10L
  )
  progeny_scores <- make_activity_matrix(progeny_long, sample_columns)
  progeny_activity <- fit_activity_contrasts(progeny_scores, design, contrast_matrix, "PROGENy", "decoupleR_ULM")
  progeny_coverage <- progeny_network[, .(
    locked_target_n = uniqueN(target),
    tested_target_n = uniqueN(intersect(target, rownames(expression_by_symbol)))
  ), by = source]
  write_tsv(progeny_coverage, file.path(output_base, "programs", paste0(cell_type, "_PROGENy_mapping_coverage.tsv")))

  dorothea_network <- dorothea_mouse[, .(
    source = as.character(tf),
    target = as.character(target),
    mor = as.numeric(mor),
    confidence = as.character(confidence)
  )]
  dorothea_long <- run_ulm(
    mat = expression_by_symbol,
    network = dorothea_network[, .(source, target, mor)],
    .source = source,
    .target = target,
    .mor = mor,
    minsize = 10L
  )
  dorothea_scores <- make_activity_matrix(dorothea_long, sample_columns)
  dorothea_activity <- fit_activity_contrasts(dorothea_scores, design, contrast_matrix, "DoRothEA_AB", "decoupleR_ULM")
  dorothea_coverage <- dorothea_network[, .(
    confidence = paste(sort(unique(confidence)), collapse = ";"),
    locked_target_n = uniqueN(target),
    tested_target_n = uniqueN(intersect(target, rownames(expression_by_symbol)))
  ), by = source]
  write_tsv(dorothea_coverage, file.path(output_base, "programs", paste0(cell_type, "_DoRothEA_mapping_coverage.tsv")))

  activity_all <- rbindlist(list(reactome_activity, progeny_activity, dorothea_activity), fill = TRUE)
  activity_all[, cell_type := cell_type]
  activity_all[, domain := domain]
  write_tsv(activity_all, file.path(output_base, "programs", paste0(cell_type, "_activity_contrasts.tsv.gz")))
  write_tsv(reactome_tests$camera, file.path(output_base, "programs", paste0(cell_type, "_Reactome_camera.tsv.gz")))
  write_optional_program_tsv(
    reactome_tests$roast,
    file.path(output_base, "programs", paste0(cell_type, "_Reactome_mroast.tsv.gz")),
    "mroast is restricted to the natural-course primary chronic-combined contrast"
  )
  write_optional_program_tsv(
    reactome_tests$fgsea,
    file.path(output_base, "programs", paste0(cell_type, "_Reactome_fgseaMultilevel.tsv.gz")),
    "fgseaMultilevel is restricted to the natural-course primary chronic-combined contrast"
  )

  loso <- data.table()
  if (domain == "natural_course") {
    chronic_vector <- contrast_matrix[, "chronic_combined", drop = FALSE]
    loso <- rbindlist(list(
      fit_loso_activity(reactome_scores, sample_meta, "group", levels, chronic_vector, "Reactome", cell_type, domain),
      fit_loso_activity(progeny_scores, sample_meta, "group", levels, chronic_vector, "PROGENy", cell_type, domain),
      fit_loso_activity(dorothea_scores, sample_meta, "group", levels, chronic_vector, "DoRothEA_AB", cell_type, domain)
    ), fill = TRUE)
    write_tsv(loso, file.path(output_base, "programs", paste0(cell_type, "_chronic_activity_LOSO.tsv.gz")))
  }

  object <- list(
    analysis_id = analysis_id,
    cell_type = cell_type,
    domain = domain,
    sample_metadata = sample_meta,
    counts_unfiltered_dimensions = dim(y_unfiltered),
    keep_filterByExpr = keep_fbe,
    keep_cpm1_in_at_least_3 = keep_cpm,
    keep_primary_intersection = keep,
    dge = y,
    design = design,
    contrasts = contrast_matrix,
    ql_fit = fit,
    voom = v,
    reactome_indices = reactome_resource$indices,
    reactome_scores = reactome_scores,
    progeny_scores = progeny_scores,
    dorothea_scores = dorothea_scores,
    program_activity_results = activity_all,
    reactome_camera = reactome_tests$camera,
    reactome_mroast = reactome_tests$roast,
    reactome_fgsea = reactome_tests$fgsea,
    chronic_activity_loso = loso
  )
  saveRDS(object, file.path(output_base, "objects", paste0(cell_type, "_analysis.rds")), compress = "xz")

  list(
    filter_summary = filter_summary,
    sample_qc = qc_table,
    gene_summary = gene_summary,
    activities = activity_all,
    camera = reactome_tests$camera,
    roast = reactome_tests$roast,
    fgsea = reactome_tests$fgsea,
    loso = loso
  )
}

## Main natural-course models: 45 libraries only.
natural_results <- lapply(c("v2a", "astrocyte", "microglia"), function(ct) {
  cat("Running natural-course model:", ct, "\n")
  run_domain_model(ct, "natural_course", roast_nrot = roast_nrot_run)
})
names(natural_results) <- c("v2a", "astrocyte", "microglia")

## Separate boundary-domain models: 18 trained libraries plus the six matching
## cell-type-specific untrained D14 comparator libraries. These results are
## physically segregated and never enter S1-S3/E1-E4 gates.
boundary_results <- lapply(c("v2a", "astrocyte"), function(ct) {
  cat("Running rotarod-boundary model:", ct, "\n")
  run_domain_model(ct, "rotarod_boundary", roast_nrot = roast_nrot_run)
})
names(boundary_results) <- c("v2a", "astrocyte")

## -------------------------------------------------------------------------
## Across-cell-type multiplicity precomputation. This deliberately tests every
## eligible locked program because upstream S2 candidates are not yet present.
## Raw P values are retained so the sequential S2-only family can be recomputed
## without re-running GSE277064. These q values must not be called S3 evidence.
## -------------------------------------------------------------------------
activity_natural <- rbindlist(lapply(names(natural_results), function(ct) natural_results[[ct]]$activities), fill = TRUE)
camera_natural <- rbindlist(lapply(names(natural_results), function(ct) {
  x <- copy(natural_results[[ct]]$camera); x[, cell_type := ct]; x[, domain := "natural_course"]; x
}), fill = TRUE)
roast_natural <- rbindlist(lapply(names(natural_results), function(ct) {
  x <- copy(natural_results[[ct]]$roast); x[, cell_type := ct]; x[, domain := "natural_course"]; x
}), fill = TRUE)
fgsea_natural <- rbindlist(lapply(names(natural_results), function(ct) {
  x <- copy(natural_results[[ct]]$fgsea)
  if (nrow(x) == 0L) return(data.table())
  x[, cell_type := ct]
  x[, domain := "natural_course"]
  x
}), fill = TRUE)

chronic_contrasts <- c("chronic_combined", "D14_vs_control", "D28_vs_control")
activity_natural[, q_global_chronic_all_locked_precompute := NA_real_]
for (collection_name in unique(activity_natural$collection)) {
  idx <- activity_natural$collection == collection_name & activity_natural$contrast %chin% chronic_contrasts
  activity_natural$q_global_chronic_all_locked_precompute[idx] <- p.adjust(activity_natural$p_value[idx], method = "BH")
}
camera_natural[, q_global_chronic_all_locked_precompute := NA_real_]
idx_camera <- camera_natural$contrast %chin% chronic_contrasts
camera_natural$q_global_chronic_all_locked_precompute[idx_camera] <- p.adjust(camera_natural$p_value[idx_camera], method = "BH")
roast_natural[, q_global_chronic_all_locked_precompute := NA_real_]
idx_roast <- roast_natural$contrast %chin% chronic_contrasts
roast_natural$q_global_chronic_all_locked_precompute[idx_roast] <- p.adjust(roast_natural$p_value[idx_roast], method = "BH")
if (nrow(fgsea_natural) > 0L) {
  fgsea_natural[, q_global_chronic_all_locked_precompute := NA_real_]
  idx_fgsea <- fgsea_natural$contrast %chin% chronic_contrasts
  fgsea_natural$q_global_chronic_all_locked_precompute[idx_fgsea] <- p.adjust(fgsea_natural$p_value[idx_fgsea], method = "BH")
}

write_tsv(activity_natural, file.path(result_root, "programs", "all_celltypes_activity_contrasts.tsv.gz"))
write_tsv(camera_natural, file.path(result_root, "programs", "all_celltypes_Reactome_camera.tsv.gz"))
write_tsv(roast_natural, file.path(result_root, "programs", "all_celltypes_Reactome_mroast.tsv.gz"))
write_tsv(fgsea_natural, file.path(result_root, "programs", "all_celltypes_Reactome_fgseaMultilevel.tsv.gz"))

make_chronic_internal_summary <- function(activity_table, p_method_table = NULL, p_method = NULL) {
  components <- dcast(
    activity_table[contrast %chin% chronic_contrasts],
    collection + cell_type + program_id ~ contrast,
    value.var = c("effect", "ci_low", "ci_high", "p_value", "q_global_chronic_all_locked_precompute")
  )
  components[, D14_D28_opposite_signs := sign(effect_D14_vs_control) * sign(effect_D28_vs_control) < 0]
  components[, combined_ci_excludes_zero := ci_low_chronic_combined > 0 | ci_high_chronic_combined < 0]
  components[, combined_direction := sign(effect_chronic_combined)]
  components[, internal_support_all_locked_activity :=
               q_global_chronic_all_locked_precompute_chronic_combined <= 0.05 &
               combined_ci_excludes_zero & !D14_D28_opposite_signs]
  components[, upstream_disease_direction_checked := FALSE]
  components[, S3_status := "NOT_IDENTIFIABLE_NO_S2_INPUT"]
  if (!is.null(p_method_table)) {
    primary <- p_method_table[contrast == "chronic_combined",
                              .(collection = "Reactome", cell_type, program_id,
                                primary_method_direction = direction,
                                primary_method_p = p_value,
                                primary_method_q_all_locked = q_global_chronic_all_locked_precompute)]
    setnames(primary, "primary_method_direction", paste0(p_method, "_direction"))
    setnames(primary, "primary_method_p", paste0(p_method, "_p"))
    setnames(primary, "primary_method_q_all_locked", paste0(p_method, "_q_all_locked"))
    components <- merge(components, primary, by = c("collection", "cell_type", "program_id"), all.x = TRUE)
  }
  components
}

reactome_internal <- make_chronic_internal_summary(
  activity_natural[collection == "Reactome"], camera_natural, "camera"
)
other_internal <- make_chronic_internal_summary(activity_natural[collection != "Reactome"])
internal_summary <- rbindlist(list(reactome_internal, other_internal), fill = TRUE)
write_tsv(internal_summary, file.path(result_root, "programs", "chronic_internal_consistency_not_S3.tsv.gz"))

s3_gate <- data.table(
  gate = "S3",
  eligible_input_source = "S2",
  eligible_input_n = 0L,
  tested_candidate_n = 0L,
  passed_candidate_n = NA_integer_,
  status = "NOT_IDENTIFIABLE_NO_S2_INPUT",
  reason = "GSE16710 S1 freeze retained zero programs, so no S2 program can enter the internally specified and locally locked sequential GSE277064 family; this is an internal analysis lock, not a public preregistration",
  exploratory_all_program_results_may_upgrade_gate = FALSE,
  spasticity_phenotype_present_in_GSE277064 = FALSE
)
write_tsv(s3_gate, file.path(result_root, "programs", "S3_GATE.tsv"))

all_filter_summaries <- rbindlist(lapply(natural_results, function(x) x$filter_summary))
all_qc <- rbindlist(lapply(natural_results, function(x) x$sample_qc))
all_gene_summaries <- rbindlist(lapply(natural_results, function(x) x$gene_summary))
all_loso <- rbindlist(lapply(natural_results, function(x) x$loso), fill = TRUE)
write_tsv(all_filter_summaries, file.path(result_root, "diagnostics", "all_celltypes_model_summary.tsv"))
write_tsv(all_qc, file.path(result_root, "diagnostics", "all_celltypes_sample_QC.tsv"))
write_tsv(all_gene_summaries, file.path(result_root, "gene_level", "all_celltypes_contrast_summary.tsv"))
write_tsv(all_loso, file.path(result_root, "programs", "all_celltypes_chronic_activity_LOSO.tsv.gz"))

boundary_activity <- rbindlist(lapply(names(boundary_results), function(ct) boundary_results[[ct]]$activities), fill = TRUE)
write_tsv(boundary_activity, file.path(result_root, "boundary_rotarod", "programs", "all_celltypes_activity_contrasts.tsv.gz"))

## Compact machine- and human-readable completion summaries.
audit_summary <- data.table(
  check = c(
    "official_count_matrix_columns", "frozen_ledger_rows", "exact_GSM_title_mapping",
    "natural_course_primary_samples", "rotarod_trained_boundary_samples",
    "natural_course_group_size", "cell_types_modeled_separately",
    "rotarod_mixed_into_natural_course", "technical_exclusions"
  ),
  observed = c(
    ncol(counts_all), nrow(ledger), sum(mapping_audit$official_title_exact_match),
    nrow(natural_ledger), nrow(exercise_ledger), paste(sort(unique(natural_counts_check$N)), collapse = ","),
    paste(names(natural_results), collapse = ","), "no", sum(all_qc$exclusion)
  ),
  expected = c("63", "63", "63", "45", "18", "3", "v2a,astrocyte,microglia", "no", "outcome-independent only"),
  status = "PASS"
)
write_tsv(audit_summary, file.path(result_root, "audit", "analysis_audit_summary.tsv"))

program_counts <- activity_natural[, .(
  programs = uniqueN(program_id),
  chronic_combined_nominal_p_0_05 = sum(contrast == "chronic_combined" & p_value < 0.05, na.rm = TRUE),
  chronic_combined_within_q_0_05 = sum(contrast == "chronic_combined" & q_within_celltype_contrast <= 0.05, na.rm = TRUE),
  chronic_combined_global_all_locked_q_0_05 = sum(contrast == "chronic_combined" & q_global_chronic_all_locked_precompute <= 0.05, na.rm = TRUE)
), by = .(collection, cell_type)]
write_tsv(program_counts, file.path(result_root, "programs", "program_result_counts.tsv"))

summary_lines <- c(
  "# GSE277064 Route B2 analysis summary",
  "",
  paste0("**Analysis ID:** ", analysis_id),
  paste0("**Completed:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope and evidence ceiling",
  "",
  "The primary analysis contains exactly 45 natural-course RiboTag libraries: V2a interneurons, astrocytes and microglia, each modeled separately across control/D3/D7/D14/D28 with n=3 independent libraries per cell-type/time stratum. The 18 post-rotarod libraries were not included in these models. They were analyzed only in a physically separate boundary-domain branch against the cell-type-matched untrained D14 comparator.",
  "",
  "GSE277064 has no spasticity phenotype. These results could only localize a previously admitted S2 post-stroke program by cell type and time; they cannot validate post-stroke spasticity, clinical efficacy or a causal mechanism. The frozen upstream S1 set is empty, so S2 input is necessarily empty and S3 is formally NOT_IDENTIFIABLE_NO_S2_INPUT.",
  "",
  "## Analysis-plan provenance",
  "",
  "Terms such as prespecified and locked refer to an internally specified, locally frozen workflow. No public preregistration was made. This terminology clarification changes no sample-inclusion rule, contrast, method, threshold, numerical result or S3 gate decision.",
  "",
  "## Primary model",
  "",
  "Raw RSEM expected counts from the three official GEO matrices were filtered by the intersection of edgeR filterByExpr and CPM >=1 in at least three libraries, TMM-normalized, and fit with robust edgeR quasi-likelihood models using ~0+timepoint separately for each cell type. Prespecified chronic, component and stage-control contrasts were evaluated. Reactome used camera for all contrasts and mroast plus fgseaMultilevel sensitivity for the primary chronic-combined contrast; PROGENy and high-confidence DoRothEA A/B activities used decoupleR ULM and limma contrasts.",
  "",
  "## Key safeguards",
  "",
  "- All 63 count columns matched the frozen GSM ledger and official GEO sample titles exactly.",
  "- No technical sample was excluded; objective library/detection/correlation flags and PCA coordinates are retained.",
  "- Multiplicity-adjusted all-program results are exploratory precomputations only. With eligible_input_n=0, none can be upgraded to S3.",
  "- Leave-one-library activity-effect ranges are retained for the chronic combined contrast.",
  "",
  "## Model summary",
  "",
  paste(capture.output(print(all_filter_summaries)), collapse = "\n"),
  "",
  "## Program counts",
  "",
  paste(capture.output(print(program_counts)), collapse = "\n")
)
writeLines(summary_lines, file.path(result_root, "analysis_summary.md"), useBytes = TRUE)

session <- capture.output(sessionInfo())
writeLines(session, file.path(result_root, "diagnostics", "sessionInfo.txt"), useBytes = TRUE)
run_manifest <- data.table(
  manifest_key = c(
    "analysis_id", "script_sha256", "sample_ledger_sha256", "contrasts_yaml_sha256",
    "gene_sets_yaml_sha256", "frozen_gene_sets_rds_sha256", "random_seed", "roast_nrot",
    "primary_population", "boundary_population", "annotation_package", "annotation_package_version",
    "reference_reported_by_GEO", "counts_reported_by_GEO",
    "analysis_plan_provenance", "public_preregistration"
  ),
  value = c(
    analysis_id,
    digest(script_path, algo = "sha256", file = TRUE),
    digest(file.path(metadata_dir, "samples.tsv"), algo = "sha256", file = TRUE),
    digest(file.path(project_root, "config", "contrasts.yaml"), algo = "sha256", file = TRUE),
    digest(file.path(project_root, "config", "gene_sets.yaml"), algo = "sha256", file = TRUE),
    digest(file.path(gene_set_dir, "gene_set_resources_frozen_2026-08-18.rds"), algo = "sha256", file = TRUE),
    random_seed, roast_nrot_run, "45 natural-course libraries", "18 trained libraries plus 6 D14 untrained comparators",
    "org.Mm.eg.db", as.character(packageVersion("org.Mm.eg.db")), "GRCm39/m32", "RSEM 1.3.1 expected counts",
    "internally specified and locally locked", "none"
  )
)
write_tsv(run_manifest, file.path(result_root, "audit", "run_manifest.tsv"))

cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Outputs:", result_root, "\n")
