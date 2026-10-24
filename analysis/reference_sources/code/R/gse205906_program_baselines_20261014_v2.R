#!/usr/bin/env Rscript
# Fixed-input, normalization-only ULM recovery and supplemental baselines.
# Never source or execute the historical whole-model driver.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
pin <- function(path) {
  list(path = normalizePath(path, mustWork = TRUE),
       bytes = unname(file.info(path)$size),
       sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
jwrite <- function(value, path) {
  need(!file.exists(path), paste("Existing output:", basename(path)))
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE,
                       digits = 17, na = "null", null = "null")
}
twrite <- function(value, path) {
  need(!file.exists(path), paste("Existing output:", basename(path)))
  data.table::fwrite(value, path, sep = "\t", na = "NA", quote = TRUE)
}
read_table <- function(path) data.table::fread(path, na.strings = "NA", check.names = FALSE)
close_num <- function(a, b, tolerance = 1e-8) {
  length(a) == length(b) && length(a) > 0L &&
    all(is.finite(a) & is.finite(b) & abs(a - b) <= tolerance * (1 + abs(b)))
}
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
groups_fixed <- c("sham", "SNI", "PLM")
contrasts_fixed <- list(local_disease = c(-1, 1, 0), intervention = c(0, -1, 1))
program_key <- function(collection, id) paste(collection, id, sep = "::")
safe_max <- function(values) if (any(is.finite(values))) max(values[is.finite(values)]) else NA_real_
safe_cor <- function(a, b) {
  good <- is.finite(a) & is.finite(b)
  if (sum(good) < 2L || length(unique(a[good])) < 2L || length(unique(b[good])) < 2L) return(NA_real_)
  cor(a[good], b[good], method = "spearman")
}
pair_metrics <- function(left, right) {
  good <- is.finite(left) & is.finite(right)
  left <- left[good]; right <- right[good]
  n <- length(left); k <- min(20L, n)
  if (n == 0L) return(list(n = 0L, rho = NA_real_, top_k = 0L, left = 0L, right = 0L, overlap = 0L))
  a <- left <= sort(left)[k]; b <- right <= sort(right)[k]
  list(n = n, rho = safe_cor(left, right), top_k = k, left = sum(a), right = sum(b), overlap = sum(a & b))
}

cell_statistics <- function(scores, group, contrast) {
  need(is.matrix(scores) && all(is.finite(scores)) && ncol(scores) == length(group) &&
         !anyDuplicated(colnames(scores)) && setequal(unique(group), groups_fixed), "Score/design content")
  counts <- as.numeric(table(factor(group, levels = groups_fixed)))
  need(all(counts > 1L) && length(contrast) == length(groups_fixed), "Design rank/group counts")
  means <- sapply(groups_fixed, function(g) rowMeans(scores[, group == g, drop = FALSE]))
  variances <- sapply(groups_fixed, function(g) apply(scores[, group == g, drop = FALSE], 1, var))
  residual <- scores - means[, match(group, groups_fixed), drop = FALSE]
  sse <- rowSums(residual^2); df_residual <- ncol(scores) - length(groups_fixed)
  mse <- sse / df_residual; total <- rowSums((scores - rowMeans(scores))^2)
  invalid <- !is.finite(sse) | sse <= 64 * .Machine$double.eps * total | total == 0
  effect <- as.numeric(means %*% contrast); variance_factor <- sum(contrast^2 / counts)
  se <- sqrt(mse * variance_factor); t_value <- effect / se
  se[invalid] <- NA_real_; t_value[invalid] <- NA_real_
  p_value <- 2 * pt(-abs(t_value), df_residual); critical <- qt(.975, df_residual)
  restricted <- means - outer(effect, (contrast / counts) / variance_factor)
  sse_restricted <- rowSums((scores - restricted[, match(group, groups_fixed), drop = FALSE])^2)
  nested_f <- (sse_restricted - sse) / mse; nested_f[invalid] <- NA_real_
  ranks <- t(apply(scores, 1, rank, ties.method = "average"))
  rank_means <- sapply(groups_fixed, function(g) rowMeans(ranks[, group == g, drop = FALSE]))
  twice_rank_sums <- sapply(groups_fixed, function(g) rowSums(2 * ranks[, group == g, drop = FALSE]))
  common <- prod(counts)
  rank_effect <- as.numeric(twice_rank_sums %*% ((2 * contrast) * common / counts)) / (4 * common)
  leverage <- 1 / counts[match(group, groups_fixed)]
  studentized <- sweep(sweep(residual, 1, sqrt(mse), "/"), 2, sqrt(1 - leverage), "/")
  cooks <- sweep(studentized^2, 2, leverage / (length(groups_fixed) * (1 - leverage)), "*")
  studentized[invalid, ] <- NA_real_; cooks[invalid, ] <- NA_real_
  out <- data.table::data.table(
    program_key = rownames(scores), status = ifelse(invalid, "numerically_degenerate", "ok_nominal_working_model"),
    estimate = effect, SE = se, t = t_value, PValue = p_value,
    ci_low = effect - critical * se, ci_high = effect + critical * se,
    residual_df = df_residual, mse = mse, sse_full = sse, sse_restricted = sse_restricted,
    nested_F = nested_f, restricted_constraint = as.numeric(restricted %*% contrast),
    no_change_fixed_effect = 0, rank_effect = rank_effect,
    rank_tied_observations = apply(scores, 1, function(v) sum(duplicated(v) | duplicated(v, fromLast = TRUE))),
    variance_ratio = apply(variances, 1, function(v) if (min(v) > 0) max(v) / min(v) else NA_real_),
    max_abs_studentized = apply(abs(studentized), 1, safe_max), max_Cook = apply(cooks, 1, safe_max),
    max_Cook_sample = vapply(seq_len(nrow(scores)), function(i) {
      if (invalid[i]) NA_character_ else colnames(scores)[which.max(cooks[i, ])]
    }, character(1)))
  for (j in seq_along(groups_fixed)) {
    out[[paste0("mean_", groups_fixed[j])]] <- means[, j]
    out[[paste0("restricted_mean_", groups_fixed[j])]] <- restricted[, j]
    out[[paste0("rank_mean_", groups_fixed[j])]] <- rank_means[, j]
    out[[paste0("variance_", groups_fixed[j])]] <- variances[, j]
  }
  out
}

analytic_ulm <- function(expression, network) {
  sources <- sort(unique(network$source))
  weights <- matrix(0, nrow(expression), length(sources), dimnames = list(rownames(expression), sources))
  for (i in seq_len(nrow(network))) weights[network$target[i], network$source[i]] <- network$mor[i]
  r <- cor(weights, expression)
  r * sqrt((nrow(expression) - 2) / ((1 - r + 1e-20) * (1 + r + 1e-20)))
}
run_ulm <- function(expression, network) {
  long <- decoupleR::run_ulm(mat = expression, network = network[, c("source", "target", "mor"), with = FALSE],
                            .source = source, .target = target, .mor = mor, minsize = 10L)
  long <- data.table::as.data.table(long)
  wide <- data.table::dcast(long, source ~ condition, value.var = "score")
  ids <- wide$source; wide[, source := NULL]
  result <- as.matrix(wide); rownames(result) <- ids
  result[, colnames(expression), drop = FALSE]
}

self_test <- function() {
  group <- rep(groups_fixed, each = 3L)
  scores <- matrix(sin(seq_len(180L) / 7) + seq_len(180L) / 100, 20L, 9L,
                   dimnames = list(paste0("p", 1:20), paste0("s", 1:9)))
  design <- model.matrix(~ 0 + factor(group, levels = groups_fixed))
  for (contrast in contrasts_fixed) {
    result <- cell_statistics(scores, group, contrast)
    fit <- lm.fit(design, scores[1, ])
    v <- as.numeric(t(contrast) %*% solve(crossprod(design)) %*% contrast)
    need(close_num(result$estimate[1], sum(fit$coefficients * contrast)), "Synthetic QR effect")
    need(close_num(result$SE[1], sqrt(sum(fit$residuals^2) / 6 * v)), "Synthetic QR SE")
    need(close_num(result$nested_F, result$t^2) && all(abs(result$restricted_constraint) < 1e-12), "Synthetic restriction")
    need(close_num(cell_statistics(scores * 2 + 5, group, contrast)$estimate, 2 * result$estimate), "Native scale preserved")
    need(close_num(cell_statistics(scores * 2 + 5, group, contrast)$rank_effect, result$rank_effect), "Rank affine")
    need(close_num(cell_statistics(-scores, group, contrast)$rank_effect, -result$rank_effect), "Rank reversal")
  }
  tied <- rbind(c(1, 2, 9, 3, 4, 5, 6, 7, 8), rep(4, 9), rep(1:3, each = 3))
  colnames(tied) <- paste0("s", 1:9)
  result <- cell_statistics(tied, group, contrasts_fixed[[1]])
  need(result$rank_effect[1] == 0 && is.na(result$PValue[2]) && is.na(result$PValue[3]), "Ties/exact zero/degenerate")
  need(identical(rank(c(1, 1, 4)), c(1.5, 1.5, 3)), "Average ranks")
  bad <- scores; bad[1, 1] <- NA_real_
  need(inherits(tryCatch(cell_statistics(bad, group, contrasts_fixed[[1]]), error = function(e) e), "error"), "Nonfinite rejection")
  need(close_num(p.adjust(c(.01, .1, NA_real_), "BH", n = 3)[1:2], c(.03, .15)), "BH fixed denominator")
  toy <- round(exp(matrix(seq_len(1800) %% 41, 200, 9) / 20) * 30)
  y <- edgeR::normLibSizes(edgeR::DGEList(toy), method = "TMM")
  expression <- log2(sweep(y$counts + .5, 2, y$samples$lib.size * y$samples$norm.factors + 1, "/") * 1e6)
  need(close_num(as.vector(expression), as.vector(limma::voom(y, design = design, plot = FALSE)$E), 1e-12), "Synthetic voom E")
  rownames(expression) <- paste0("g", 1:200); colnames(expression) <- paste0("s", 1:9)
  network <- data.table::data.table(source = rep(c("A", "B"), each = 20),
                                  target = paste0("g", 1:40), mor = rep(c(-1, 1), 20))
  actual <- run_ulm(expression, network); expected <- analytic_ulm(expression, network)
  need(close_num(as.vector(actual), as.vector(expected), 1e-12), "Synthetic full-universe ULM")
  list(status = "PASS_SYNTHETIC_ONLY", checks = c("three_group_QR_effect_SE_df6", "single_constraint_F_t2",
    "original_score_scale", "rank_affine_reversal_ties_exactzero", "degenerate_unavailable_inference",
    "nonfinite_rejection", "fixed_BH", "voom_E", "whole_symbol_universe_ULM"))
}

read_gtf_symbols <- function(path, ids) {
  con <- file(path, "rt"); on.exit(close(con))
  parts <- list()
  extract <- function(text, field) {
    match <- regmatches(text, regexec(paste0('(?:^|;)[[:space:]]*', field, ' "([^"]*)"'), text, perl = TRUE))
    vapply(match, function(x) if (length(x) > 1L) x[2] else NA_character_, character(1))
  }
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    genes <- lines[!startsWith(lines, "#") & grepl("\tgene\t", lines, fixed = FALSE)]
    if (!length(genes)) next
    split <- strsplit(genes, "\t", fixed = TRUE)
    need(all(lengths(split) == 9L) && all(vapply(split, `[[`, character(1), 3L) == "gene"), "GTF gene record shape")
    attributes <- vapply(split, `[[`, character(1), 9L)
    parts[[length(parts) + 1L]] <- data.table::data.table(gene_id = extract(attributes, "gene_id"),
      gene_name = extract(attributes, "gene_name"), gene_biotype = extract(attributes, "gene_biotype"))
  }
  annotation <- data.table::rbindlist(parts)
  need(!anyNA(annotation$gene_id) && !anyDuplicated(annotation$gene_id) && setequal(ids, annotation$gene_id), "Complete GTF gene universe")
  annotation <- annotation[match(ids, annotation$gene_id)]
  symbols <- trimws(annotation$gene_name)
  missing <- is.na(symbols) | !nzchar(symbols)
  need(mean(!missing) >= .8, "Original annotation nonblank coverage rule")
  symbols[missing] <- ids[missing]
  list(symbols = symbols, fallback_n = sum(missing))
}
plain_frame <- function(value) {
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  value[] <- lapply(value, function(column) if (is.factor(column)) as.character(column) else column)
  data.table::as.data.table(value)
}
build_networks <- function(resource, expression_symbols) {
  orthology <- plain_frame(resource$ortholog_one_to_one)
  rat_map <- unique(orthology[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
                              .(human_symbol = as.character(human_symbol), rat_symbol = as.character(other_symbol))])
  need(!anyDuplicated(rat_map$human_symbol) && !anyDuplicated(rat_map$rat_symbol), "Strict rat orthology")
  progeny <- plain_frame(resource$progeny_human)
  progeny <- unique(merge(progeny[, .(source = as.character(pathway), human_symbol = as.character(gene_symbol),
                                    mor = as.numeric(weight))], rat_map, by = "human_symbol", allow.cartesian = TRUE)[,
                                    .(source, target = rat_symbol, mor)])
  doro <- plain_frame(resource$dorothea_rat)
  need(all(doro$confidence %in% c("A", "B")), "Frozen rat A/B confidence")
  doro <- unique(doro[, .(source = as.character(tf), target = as.character(target), mor = as.numeric(mor),
                        confidence = as.character(confidence))])
  networks <- list(PROGENy = progeny, DoRothEA_AB = doro)
  coverage <- lapply(names(networks), function(collection) {
    tab <- networks[[collection]][, .(source_feature_n = data.table::uniqueN(target),
       tested_feature_n = data.table::uniqueN(intersect(target, expression_symbols))), by = source]
    tab[, collection := collection]
    tab[, tested := tested_feature_n >= 10L & tested_feature_n <= 500L]
    tab
  })
  coverage <- data.table::rbindlist(coverage)
  network_names <- names(networks)
  networks <- lapply(network_names, function(collection) {
    target_collection <- collection
    eligible <- coverage[collection == target_collection & tested == TRUE]$source
    tab <- networks[[target_collection]][source %in% eligible & target %in% expression_symbols]
    need(!anyDuplicated(tab[, .(source, target)]) && all(is.finite(tab$mor)), "Network unique source-target weights")
    tab[, .(source, target, mor)]
  }) |> setNames(network_names)
  list(networks = networks, coverage = coverage)
}

build_inputs <- function(model, out) {
  manifest <- jsonlite::fromJSON(file.path(model, "file_manifest.json"), simplifyVector = FALSE)
  need(pin(file.path(model, "file_manifest.json"))$sha256 == "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7", "Original model manifest identity")
  need(pin(file.path(model, "COMPLETED"))$sha256 == "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e", "Original model completion identity")
  formal <- "/root/autodl-tmp/GSE205906_full_read_counts_fasp_v18e_20261005/formal"
  public <- "/root/autodl-tmp/route_b2_public_analysis/public-input-stage-v1"
  input <- manifest$provenance$eleven_input_pins
  paths <- list(matrix = file.path(formal, input$matrix$relative_path),
    GTF = "/root/autodl-tmp/GSE205906_full_read_reconstruction_156_20261003/references/derived/Rattus_norvegicus.Rnor_6.0.101.gtf",
    overlay = file.path(formal, input$overlay$relative_path), manifest = file.path(formal, input$manifest$relative_path))
  for (role in c("PROGENy", "DoRothEA_AB", "orthology", "frozen_resource_RDS", "gene_set_manifest")) {
    paths[[role]] <- file.path(public, input[[role]]$relative_path)
  }
  for (role in names(paths)) {
    actual <- pin(paths[[role]]); expected <- input[[role]]
    need(actual$bytes == expected$bytes && actual$sha256 == expected$sha256, paste("Frozen source identity", role))
  }
  completed <- jsonlite::fromJSON(file.path(model, "COMPLETED"), simplifyVector = FALSE)
  for (name in c("gene_edger.tsv", "activity_results.tsv", "mapping_coverage.tsv", "design_out.tsv",
                 "contrast_out.tsv", "filter_out.tsv", "annotation_metrics.tsv", "gate.tsv")) {
    actual <- pin(file.path(model, name)); expected <- completed$files[[name]]
    need(actual$bytes == expected$bytes && actual$sha256 == expected$sha256, paste("Retained product identity", name))
    paths[[name]] <- file.path(model, name)
  }
  paths$model_manifest <- file.path(model, "file_manifest.json"); paths$model_complete <- file.path(model, "COMPLETED")
  paths$worker <- file.path(out, "worker.R"); paths$definition <- file.path(out, "definition.md")
  table <- data.table::rbindlist(lapply(names(paths), function(role) {
    p <- pin(paths[[role]]); data.table::data.table(role = role, path = p$path, bytes = p$bytes, sha256 = p$sha256)
  }))
  list(paths = paths, table = table)
}
verify_inputs <- function(table) {
  for (i in seq_len(nrow(table))) {
    actual <- pin(table$path[i])
    need(actual$bytes == table$bytes[i] && actual$sha256 == table$sha256[i], paste("Input changed", table$role[i]))
  }
  invisible(TRUE)
}

recover <- function(model, out, continue_pre_score = FALSE) {
  if (continue_pre_score) {
    need(file.exists(file.path(out, "RECOVERY_STARTED.json")) &&
           !any(file.exists(file.path(out, c("RECOVERY_COMPLETE.json", "BASELINE_STARTED.json",
             "restored_scores.rds", "scores.tsv", "metadata.tsv", "retained_native.tsv",
             "new_mapping_coverage.tsv", "new_normalization.tsv")))),
         "Continuation requires previous pre-score stop with no scientific products")
    jwrite(list(started_utc = utc(), actual_worker = pin(file.path(out, "worker_v2.R")),
      change = "Empty nonfinite-vector comparison uses numeric values; original inputs and statistics unchanged",
      prior_process_exit = 1L, original_products_preserved = TRUE), file.path(out, "RECOVERY_CONTINUE_STARTED.json"))
  } else {
    need(!file.exists(file.path(out, "RECOVERY_STARTED.json")), "Prior recovery exists; do not restart")
    jwrite(list(started_utc = utc()), file.path(out, "RECOVERY_STARTED.json"))
  }
  input <- build_inputs(model, out); paths <- input$paths
  if (continue_pre_score) {
    previous <- read_table(file.path(out, "INPUT_MANIFEST.tsv")); verify_inputs(previous)
    need(identical(previous$role, input$table$role) && identical(previous$path, input$table$path) &&
           all(previous$bytes == input$table$bytes) && identical(previous$sha256, input$table$sha256),
         "Original input/source manifest identity")
    need(jsonlite::fromJSON(file.path(out, "SYNTHETIC_TESTS.json"))$status == "PASS_SYNTHETIC_ONLY", "Original synthetic record")
  } else {
    twrite(input$table, file.path(out, "INPUT_MANIFEST.tsv"))
    jwrite(self_test(), file.path(out, "SYNTHETIC_TESTS.json"))
  }
  meta <- read_table(paths$design_out.tsv); overlay <- read_table(paths$overlay)
  need(identical(meta$sample_id, paste0("GSM", 6235292:6235300)) && identical(meta$group, rep(groups_fixed, each = 3L)), "Fixed design")
  need(nrow(overlay) == 9L && identical(overlay$GEO_sample, meta$sample_id) &&
         identical(overlay$group, rep(c("JA", "SNI", "PLM"), each = 3L)) && !anyDuplicated(overlay$label), "Original JA/sham library crosswalk")
  meta[, count_label := overlay$label]
  meta[, leverage := 1 / 3]
  native_contrasts <- read_table(paths$contrast_out.tsv)
  need(identical(native_contrasts$coefficient, groups_fixed), "Contrast group order")
  for (name in names(contrasts_fixed)) need(identical(as.numeric(native_contrasts[[name]]), contrasts_fixed[[name]]), "Original contrast identity")
  tab <- read_table(paths$matrix)
  need(names(tab)[1] == "gene_id" && identical(names(tab)[-1], meta$count_label), "Count library alignment")
  ids <- tab$gene_id
  need(length(ids) == 32883L && !anyDuplicated(ids) && identical(ids, sort(ids, method = "radix")), "Complete fixed count gene universe")
  count <- as.matrix(tab[, -1, with = FALSE]); rownames(count) <- ids; colnames(count) <- meta$sample_id
  need(all(is.finite(count)) && all(count >= 0) && all(count == round(count)), "Nonnegative integer count content")
  annotation <- read_gtf_symbols(paths$GTF, ids)
  native_gene <- read_table(paths$gene_edger.tsv)
  gene_ids <- native_gene[contrast == "local_disease"]$gene_id
  need(length(gene_ids) == 14217L && !anyDuplicated(gene_ids) && identical(gene_ids, ids[ids %in% gene_ids]), "Native retained gene order")
  need(identical(native_gene[contrast == "intervention"]$gene_id, gene_ids), "Both native retained gene universes")
  selected <- match(gene_ids, ids)
  for (name in names(contrasts_fixed)) need(identical(native_gene[contrast == name]$symbol, annotation$symbols[selected]), "GTF-native symbol agreement")
  retained <- count[selected, , drop = FALSE]
  y <- edgeR::normLibSizes(edgeR::DGEList(retained), method = "TMM")
  effective <- y$samples$lib.size * y$samples$norm.factors
  expression <- log2(sweep(retained + .5, 2, effective + 1, "/") * 1e6)
  symbols <- annotation$symbols[selected]
  # The original collapse retains first-occurrence symbol order.
  collapsed <- rowsum(expression, symbols, reorder = FALSE)
  collapsed <- collapsed / as.numeric(table(factor(symbols, levels = rownames(collapsed))))
  resource <- readRDS(paths$frozen_resource_RDS)
  resource_roles <- c(PROGENy = "progeny_human", DoRothEA_AB = "dorothea_rat", orthology = "ortholog_one_to_one")
  for (role in names(resource_roles)) {
    original_frame <- as.data.frame(plain_frame(resource[[resource_roles[[role]]]]))
    character_fields <- names(original_frame)[vapply(original_frame, is.character, logical(1))]
    # TSV's automatic numeric-ID inference is not a change to RDS label semantics.
    source_frame <- as.data.frame(data.table::fread(paths[[role]],
      colClasses = list(character = character_fields), na.strings = "NA"))
    need(identical(names(original_frame), names(source_frame)) && nrow(original_frame) == nrow(source_frame),
         paste("RDS versus source TSV structure", role))
    for (field in names(original_frame)) {
      a <- original_frame[[field]]; b <- source_frame[[field]]
      need(identical(is.na(a), is.na(b)), paste("RDS versus source TSV missingness", role, field))
      if (is.numeric(a)) {
        need(is.numeric(b) && identical(is.finite(a), is.finite(b)), "Numeric resource types")
        finite <- is.finite(a)
        need(all(abs(a[finite] - b[finite]) <= 1e-10 * abs(a[finite])), "Resource numeric relative tolerance; zero exact")
        nonfinite <- !is.finite(a) & !is.na(a)
        need(all(a[nonfinite] == b[nonfinite]), "Resource nonfinite numeric identity")
      } else need(identical(a, b), paste("RDS versus source TSV exact labels", role, field))
    }
  }
  networks <- build_networks(resource, rownames(collapsed))
  native_coverage <- read_table(paths$mapping_coverage.tsv)[collection %in% names(networks$networks)]
  coverage_keys <- program_key(networks$coverage$collection, networks$coverage$source)
  native_coverage_keys <- program_key(native_coverage$collection, native_coverage$source)
  need(!anyDuplicated(native_coverage_keys) && setequal(coverage_keys, native_coverage_keys), "Complete source-domain coverage keys")
  native_coverage <- native_coverage[match(coverage_keys, native_coverage_keys)]
  for (field in c("source_feature_n", "tested_feature_n", "tested")) {
    need(all(networks$coverage[[field]] == native_coverage[[field]]), paste("Full native eligibility coverage", field))
  }
  scores <- lapply(networks$networks, function(network) run_ulm(collapsed, network))
  registry <- data.table::rbindlist(lapply(names(scores), function(collection) {
    data.table::data.table(collection = collection, program_id = rownames(scores[[collection]]))
  }))
  score <- do.call(rbind, scores); rownames(score) <- program_key(registry$collection, registry$program_id)
  native <- read_table(paths$activity_results.tsv)
  need(nrow(native) == 190L && nrow(score) == 95L && all(is.finite(score)) && identical(colnames(score), meta$sample_id), "Complete finite ULM scores")
  maximum <- list()
  for (name in names(contrasts_fixed)) {
    ref <- native[contrast == name]
    keys <- program_key(ref$collection, ref$program_id)
    need(length(keys) == nrow(score) && !anyDuplicated(keys) && setequal(keys, rownames(score)), "Every native program key")
    ref <- ref[match(rownames(score), keys)]
    means <- sapply(groups_fixed, function(g) rowMeans(score[, meta$group == g, drop = FALSE]))
    effect <- as.numeric(means %*% contrasts_fixed[[name]])
    need(close_num(effect, ref$effect), paste("Both native coefficient compatibility", name))
    maximum[[name]] <- max(abs(effect - ref$effect))
    for (family in unique(ref$collection)) {
      family_ref <- ref[collection == family]
      need(close_num(p.adjust(family_ref$p_value, "BH", n = nrow(family_ref)), family_ref$q_within_collection_contrast, 1e-12), "Retained native BH family")
    }
  }
  # Pin new score products before any supplemental baseline inference.
  restored <- list(scores = score, registry = registry, metadata = meta, coverage = networks$coverage,
                   normalization_role = "NEW_DETERMINISTIC_TMM_RECOVERY_NOT_HISTORICAL_CACHE")
  need(!file.exists(file.path(out, "restored_scores.rds")), "Existing score RDS")
  saveRDS(restored, file.path(out, "restored_scores.rds"), version = 3)
  exported <- cbind(registry, data.table::as.data.table(score))
  twrite(exported, file.path(out, "scores.tsv")); twrite(meta, file.path(out, "metadata.tsv"))
  twrite(native, file.path(out, "retained_native.tsv")); twrite(networks$coverage, file.path(out, "new_mapping_coverage.tsv"))
  normalization <- data.table::data.table(sample_id = meta$sample_id, group = meta$group,
    raw_library_size = colSums(count), retained_library_size = y$samples$lib.size,
    TMM_norm_factor = y$samples$norm.factors, effective_library_size = effective,
    role = "NEW_DETERMINISTIC_TMM_RECOVERY_NOT_HISTORICAL_CACHE")
  twrite(normalization, file.path(out, "new_normalization.tsv"))
  products <- lapply(file.path(out, c("restored_scores.rds", "scores.tsv", "metadata.tsv", "retained_native.tsv",
                                    "new_mapping_coverage.tsv", "new_normalization.tsv")), pin)
  names(products) <- vapply(products, function(x) basename(x$path), character(1))
  verify_inputs(input$table)
  jwrite(list(status = "RECOVERED_NEW_ULM_SCORES_BOTH_NATIVE_CONTRASTS_COMPATIBLE", completed_utc = utc(),
    input_manifest = pin(file.path(out, "INPUT_MANIFEST.tsv")), products = products,
    sample_n = 9L, input_gene_n = nrow(count), retained_gene_n = nrow(retained), symbol_n = nrow(collapsed),
    annotation_fallback_n = annotation$fallback_n, program_n = nrow(score), family_counts = as.list(table(registry$collection)),
    native_coefficient_maximum_difference = maximum, normalization_recomputed_only = TRUE,
    actual_recovery_worker = pin(file.path(out, "worker_v2.R")), continued_after_pre_score_comparator_stop = continue_pre_score,
    complete_RDS_source_TSV_semantics_checked = TRUE, all_source_program_eligibility_checked = TRUE,
    primary_gene_models_refit = FALSE, baseline_inference_run = FALSE, inputs_unchanged = TRUE,
    historical_sample_score_identity_authenticated = FALSE,
    runtime = c(list(R = R.version.string), setNames(lapply(c("edgeR", "limma", "decoupleR", "data.table", "jsonlite", "digest"),
      function(package) as.character(utils::packageVersion(package))), c("edgeR", "limma", "decoupleR", "data.table", "jsonlite", "digest")))),
    file.path(out, "RECOVERY_COMPLETE.json"))
  cat("PASS_RECOVERED_NEW_ULM_SCORES no primary gene refit\n")
}

run_baselines <- function(out) {
  need(!file.exists(file.path(out, "BASELINE_STARTED.json")), "Prior baseline stage exists; do not restart")
  recovery <- jsonlite::fromJSON(file.path(out, "RECOVERY_COMPLETE.json"), simplifyVector = FALSE)
  need(recovery$status == "RECOVERED_NEW_ULM_SCORES_BOTH_NATIVE_CONTRASTS_COMPATIBLE", "Recovery incomplete")
  inputs <- read_table(file.path(out, "INPUT_MANIFEST.tsv")); verify_inputs(inputs)
  for (product in recovery$products) {
    actual <- pin(product$path)
    need(actual$bytes == product$bytes && actual$sha256 == product$sha256, "Derived score product changed")
  }
  jwrite(list(started_utc = utc(), producer = pin(file.path(out, "worker_v2.R")),
              recovery = pin(file.path(out, "RECOVERY_COMPLETE.json"))), file.path(out, "BASELINE_STARTED.json"))
  restored <- readRDS(file.path(out, "restored_scores.rds")); native <- read_table(file.path(out, "retained_native.tsv"))
  score <- restored$scores; registry <- restored$registry; meta <- restored$metadata
  rows <- list(); families <- list()
  for (name in names(contrasts_fixed)) {
    result <- cell_statistics(score, meta$group, contrasts_fixed[[name]])
    result <- cbind(result, registry)
    result[, contrast := name]
    ref <- native[contrast == name]
    ref <- ref[match(result$program_key, program_key(ref$collection, ref$program_id))]
    for (field in c("effect", "se", "ci_low", "ci_high", "p_value", "q_within_collection_contrast")) result[[paste0("retained_", field)]] <- ref[[field]]
    result[, BH_q_nominal := p.adjust(PValue, "BH", n = .N), by = collection]
    result[, rank_inference_status := "DESCRIPTIVE_ONLY_NO_P_Q_CI"]
    result[, nominal_inference_status := "WORKING_HOMOSKEDASTIC_GAUSSIAN_NOT_PROVEN"]
    result[, score_scale := "ORIGINAL_ULM_T_LIKE_NOT_Z_STANDARDIZED"]
    for (family in unique(result$collection)) {
      tab <- result[collection == family]
      native_q <- tab$retained_q_within_collection_contrast; nominal_q <- tab$BH_q_nominal
      paired <- is.finite(tab$estimate) & is.finite(tab$rank_effect) & tab$estimate != 0 & tab$rank_effect != 0
      summaries <- list(collection = family, contrast = name, n = nrow(tab), sample_n = 9L, residual_df = 6L,
        invalid_n = sum(tab$status != "ok_nominal_working_model"),
        native_q_eligible_n = sum(is.finite(native_q)), OLS_q_eligible_n = sum(is.finite(nominal_q)),
        native_q05_n = sum(native_q < .05, na.rm = TRUE), OLS_q05_n = sum(nominal_q < .05, na.rm = TRUE),
        q05_overlap_n = sum(native_q < .05 & nominal_q < .05, na.rm = TRUE),
        native_vs_OLS_P = pair_metrics(tab$PValue, tab$retained_p_value),
        rank_direction_paired_n = sum(paired), rank_same_direction_n = sum(sign(tab$estimate[paired]) == sign(tab$rank_effect[paired])),
        rank_zero_n = sum(tab$rank_effect == 0), rank_missing_n = sum(!is.finite(tab$rank_effect)),
        rank_absolute_effect_rho = safe_cor(abs(tab$estimate), abs(tab$rank_effect)),
        native_effect_maximum_difference = max(abs(tab$estimate - tab$retained_effect)),
        variance_ratio = as.list(quantile(tab$variance_ratio, c(0, .5, .9, 1), na.rm = TRUE)),
        max_Cook = as.list(quantile(tab$max_Cook, c(0, .5, .9, 1), na.rm = TRUE)))
      families[[paste(family, name, sep = "/")]] <- summaries
    }
    rows[[name]] <- result
  }
  result <- data.table::rbindlist(rows)
  twrite(result, file.path(out, "BASELINES.tsv"))
  verify_inputs(inputs)
  jwrite(list(status = "BASELINES_COMPUTED_PENDING_INDEPENDENT_REVIEW", completed_utc = utc(),
    recovery = pin(file.path(out, "RECOVERY_COMPLETE.json")), result = pin(file.path(out, "BASELINES.tsv")),
    rows = nrow(result), family_cells = families, primary_gene_models_refit = FALSE,
    q_threshold_rule = "STRICT_LT_0.05_NATIVE_REPORTING_CONVENTION",
    gates_changed = FALSE, inputs_unchanged = TRUE, full_goal_complete = FALSE), file.path(out, "BASELINE_COMPLETE.json"))
  cat("PASS_BASELINES_COMPUTED pending independent numerical/report review\n")
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  need(length(args) >= 1L, "Mode required")
  if (args[1] == "self-test") print(self_test())
  else if (args[1] == "recover") {
    need(length(args) == 3L && dir.exists(args[3]), "recover MODEL_DIR NEW_OUTPUT_DIR")
    recover(args[2], args[3])
  } else if (args[1] == "continue-recover") {
    need(length(args) == 3L && dir.exists(args[3]), "continue-recover MODEL_DIR EXISTING_OUTPUT_DIR")
    recover(args[2], args[3], continue_pre_score = TRUE)
  } else if (args[1] == "run") {
    need(length(args) == 2L, "run OUTPUT_DIR")
    run_baselines(args[2])
  } else stop("Unknown mode", call. = FALSE)
}
