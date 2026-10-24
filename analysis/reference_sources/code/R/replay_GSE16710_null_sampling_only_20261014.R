# Reconstruct archived random-set sampling and family-reference draws only.
# Never source the historical analyzer, build activities, fit, adjust BH or gate.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1L]], mustWork = TRUE, winslash = "/")
out <- args[[2L]]
if (file.exists(out) || !dir.exists(dirname(out))) stop("Unused R output directory required")
if (Sys.setlocale("LC_COLLATE", "C") != "C") stop("Archived C collation unavailable")
if (as.character(getRversion()) != "4.6.1") stop("This replay specifies current/archived R4.6.1 only")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
suppressPackageStartupMessages(library(data.table))
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label)
  cat("PASS ", label, "\n", sep = "")
}
sha <- function(p) tolower(unname(tools::sha256sum(p)))
base <- file.path(root, "results/GSE16710/diagnostics/matched_random_gene_set_null_v1")
receipt_path <- file.path(base, "matched_random_gene_set_null_COMPLETE.json")
manifest_path <- file.path(base, "matched_random_gene_set_null_output_manifest.tsv")
check(sha(receipt_path) == "9901592fe81580e5ed814f05bfbf8dfa30e8c4a4cec080b640c754c065547ac4", "fixed original completion identity")
receipt <- jsonlite::fromJSON(receipt_path)
check(sha(manifest_path) == receipt$output_manifest_sha256, "fixed original output manifest")
source <- file.path(root, receipt$script_relative_path)
spec <- file.path(root, receipt$frozen_addendum_relative_path)
check(sha(source) == receipt$script_sha256 && sha(spec) == receipt$frozen_addendum_sha256, "original analyzer/spec identity; not executed")
paths <- c(receipt_path, manifest_path, source, spec)
roles <- c("retained_completion_metadata", "retained_output_manifest", "historical_analyzer_not_executed", "post_outcome_matching_specification")
for (item in receipt$inputs) {
  p <- file.path(root, item$relative_path)
  check(sha(p) == item$sha256, paste("fixed input", item$relative_path))
  paths <- c(paths, p); roles <- c(roles, "frozen_derived_input_or_linked_context")
}
mf <- fread(manifest_path)
for (i in seq_len(nrow(mf))) {
  p <- file.path(base, mf$file_name[[i]])
  check(file.info(p)$size == mf$bytes[[i]] && sha(p) == mf$sha256[[i]], paste("retained output identity", mf$file_name[[i]]))
  paths <- c(paths, p); roles <- c(roles, "retained_comparison_or_environment_not_raw")
}
before <- setNames(vapply(paths, sha, character(1)), paths)
check(identical(c(receipt$random_sets, receipt$expression_bins, receipt$candidate_sets_per_draw,
                  receipt$sampling_seed, receipt$family_reference_seed, receipt$family_reference_draws),
                c(10000L, 20L, 5L, 20260831L, 20260832L, 10000L)), "fixed matching parameters and seeds")

x <- readRDS(file.path(root, receipt$inputs$expression$relative_path))
ledger <- fread(file.path(root, receipt$inputs$sample_ledger$relative_path))
qc <- fread(file.path(root, receipt$inputs$qc_decision$relative_path))
check(is.matrix(x) && !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x)), "unique frozen expression IDs")
check(!anyDuplicated(ledger$sample_id) && !anyDuplicated(qc$sample_id) &&
      setequal(colnames(x), ledger$sample_id) && setequal(ledger$sample_id, qc$sample_id), "27-array matrix/metadata/QC coverage")
check(all(toupper(as.character(qc$outcome_information_used)) == "FALSE"), "stored QC field only; not historical blinding authentication")
meta <- merge(ledger, qc, by = "sample_id", sort = FALSE)
meta <- meta[order(as.integer(sample_order))]
check(identical(meta[primary_analysis_decision != "retain", sample_id], "GSM418682"), "fixed exclusion GSM418682")
samples <- meta[primary_analysis_decision == "retain", sample_id]
check(ncol(x) == 27L && length(samples) == 26L, "27 deposited arrays/26 retained arrays, not verified animal independence")
x <- x[, samples, drop = FALSE]
means <- rowMeans(x)
sds <- matrixStats::rowSds(x)
means <- means[is.finite(means) & is.finite(sds) & sds > 0]
ids <- names(means)
ord <- order(means, ids, na.last = NA)
bins <- setNames(integer(length(ids)), ids)
bins[ord] <- pmin(20L, ceiling(seq_along(ord) * 20L / length(ord)))
pool <- split(names(bins), factor(bins, levels = seq_len(20L)))
check(all(lengths(pool) > 0L), "twenty nonempty frozen-expression bins")
frozen <- readRDS(file.path(root, receipt$inputs$frozen_gene_sets$relative_path))
orth <- as.data.table(frozen$ortholog_one_to_one)[other_tax_id == 10116,
                                               .(human_symbol, rat_symbol = other_symbol)]
react <- merge(as.data.table(frozen$reactome)[, .(term_id, human_symbol = gene_symbol)],
               orth, by = "human_symbol", allow.cartesian = TRUE)
react <- unique(react[rat_symbol %in% ids, .(term_id, rat_symbol)])
templates <- lapply(split(react$rat_symbol, react$term_id), function(g) sort(unique(g)))
templates <- templates[lengths(templates) >= 10L & lengths(templates) <= 500L]
templates <- templates[sort(names(templates))]
check(length(templates) == 1656L, "1656 full eligible templates, no q selection")

saved <- fread(file.path(base, "matched_random_gene_set_null_10000.tsv.gz"))
members <- fread(file.path(base, "matched_random_gene_set_membership_10000.tsv.gz"))
real <- fread(file.path(base, "matched_random_gene_set_real_reactome_calibration_model.tsv.gz"))
references <- fread(file.path(base, "matched_random_gene_set_family_reference_10000.tsv.gz"))
setorder(saved, random_set_id)
expected_ids <- sprintf("RGS%05d", seq_len(10000L))
check(identical(saved$random_set_id, expected_ids) && !anyDuplicated(saved$random_set_id), "exact retained random-set IDs/order")
check(!anyDuplicated(members[, .(random_set_id, gene_symbol)]) && nrow(members) == 507966L,
      "507966 distinct within-set memberships")
saved_members <- split(members$gene_symbol, factor(members$random_set_id, levels = expected_ids))
check(setequal(real$template_term_id, names(templates)) && !anyDuplicated(real$template_term_id), "exact real-model template coverage only; no refit")
check(identical(references$reference_draw, seq_len(10000L)), "exact family-reference ID order")

# Tiny independent matching semantics check: exact ties select the first candidate.
check(which.min(c(0, 0, 1, 2, 3)) == 1L &&
      identical(tabulate(c(1L, 1L, 2L), nbins = 3L), c(2L, 1L, 0L)), "synthetic first-tie and exact-bin semantics")
dir.create(out, recursive = FALSE)
started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
set.seed(20260831L)
assignment <- rep(names(templates), each = 10000L %/% length(templates))
assignment <- c(assignment, sample(names(templates), 10000L - length(assignment), replace = FALSE))
assignment <- sample(assignment, length(assignment), replace = FALSE)
check(identical(assignment, saved$template_term_id), "all 10000 template assignments replayed exactly")
comparisons <- vector("list", 10000L)
for (i in seq_len(10000L)) {
  genes <- templates[[assignment[[i]]]]
  target_bins <- tabulate(bins[genes], nbins = 20L)
  target_mean <- mean(means[genes])
  candidates <- vector("list", 5L)
  deltas <- numeric(5L)
  for (j in seq_len(5L)) {
    picked <- unlist(lapply(seq_len(20L), function(b) {
      k <- target_bins[[b]]
      if (k == 0L) return(character())
      sample(pool[[b]], k, replace = FALSE)
    }), use.names = FALSE)
    candidates[[j]] <- sort(picked)
    deltas[[j]] <- abs(mean(means[candidates[[j]]]) - target_mean)
  }
  chosen <- which.min(deltas)
  picked <- candidates[[chosen]]
  # Use sorted-vector comparison rather than the stored TRUE flags.
  comparisons[[i]] <- data.table(
    random_set_id = expected_ids[[i]], template_term_id = assignment[[i]],
    selected_candidate = chosen, recorded_candidate = saved$chosen_candidate_index[[i]],
    candidate_index_match = chosen == saved$chosen_candidate_index[[i]],
    entire_member_vector_match = identical(picked, sort(saved_members[[i]])),
    exact_bin_vector_match = identical(tabulate(bins[picked], nbins = 20L), target_bins),
    size_match = length(picked) == saved$random_size[[i]] && length(picked) == length(genes),
    mean_expression_error = abs(mean(means[picked]) - saved$random_mean_expression[[i]]),
    delta_export_error = abs(deltas[[chosen]] - saved$absolute_mean_expression_difference[[i]]),
    best_of_five_only_expression_distance = deltas[[chosen]] == min(deltas))
  if (i %% 2000L == 0L) cat("SAMPLING_REPLAY ", i, "/10000; no activities/model fitted\n", sep = "")
}
cmp <- rbindlist(comparisons)
check(all(cmp$candidate_index_match & cmp$entire_member_vector_match & cmp$exact_bin_vector_match & cmp$size_match &
          cmp$best_of_five_only_expression_distance) && max(cmp$mean_expression_error) < 1e-12 && max(cmp$delta_export_error) < 1e-12,
      "10000/10000 candidate choices and entire memberships reproduced; exported means within 1e-12")
fwrite(cmp, file.path(out, "SAMPLING_REPLAY_COMPARISON.tsv"), sep = "\t")

# Reuse saved moderated statistics as read-only input, never regenerate them.
set.seed(20260832L)
by_template <- split(saved$absolute_moderated_t, saved$template_term_id)
family <- numeric(10000L)
for (term in real$template_term_id) family <- family + sample(by_template[[term]], 10000L, replace = TRUE)
family <- family / nrow(real)
ref_error <- abs(family - references$mean_absolute_moderated_t)
check(max(ref_error) < 1e-12, "10000/10000 family-reference draws reproduced within TSV tolerance")
observed <- mean(real$absolute_moderated_t)
tail <- (1 + sum(family >= observed)) / 10001
check(abs(tail - receipt$family_summary$exploratory_family_reference_tail_area) < 1e-12,
      "plus-one descriptive family tail reproduced; not a confirmatory P value")
fwrite(data.table(reference_draw = seq_len(10000L), reconstructed_mean_absolute_t = family,
                 recorded_mean_absolute_t = references$mean_absolute_moderated_t, absolute_export_error = ref_error),
       file.path(out, "FAMILY_REFERENCE_REPLAY_COMPARISON.tsv"), sep = "\t")
check(identical(before, setNames(vapply(paths, sha, character(1)), paths)), "all original inputs/source/output identities unchanged")
file_audit <- data.table(
  file_role = roles, expected_path = paths, actual_path_exists = paths,
  bytes = as.numeric(file.info(paths)$size), sha256 = unname(before),
  source_version_accession = "GSE16710 retained Aug31 post-outcome adjunct and its fixed inputs",
  download_date = "NOT_REAUTHENTICATED", licence_data_use = "NONCOMMERCIAL_RESEARCH; PUBLIC_REDISTRIBUTION_NOT_CLEARED",
  sample_metadata_match = "27 unique array IDs and QC cover expression columns; 26 retained; animal linkage not authenticated",
  inclusion_exclusion = "exclude GSM418682 only; no new exclusion", script_reference = receipt$script_relative_path,
  duplicate_stale_check = "exact receipt/manifest paths only; no latest/glob substitution",
  manifest_match = "EXACT_ORIGINAL_SHA_AND_BYTES", status = "NEEDS REVIEW",
  discrepancy = "original animal/historical custody/download/redistribution not proven by matching current bytes",
  required_action = "retain limited current-file/sampling-replay claim; do not upgrade gate or historical provenance")
fwrite(file_audit, file.path(out, "FILE_LEVEL_AUDIT.tsv"), sep = "\t")
session <- capture.output(sessionInfo())
writeLines(session, file.path(out, "SESSION_INFO.txt"))
outputs <- list.files(out, full.names = TRUE)
report <- list(
  status = "PASS_FIXED_SEED_MATCHING_AND_STORED_STATISTIC_REFERENCE_REPLAY_ONLY",
  overall_provenance_status = "NEEDS REVIEW", start_utc = started,
  end_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  input_files_bound_and_unchanged = length(paths), input_sha256 = as.list(before),
  producer_sha256 = sha(file.path(root, "code/R/replay_GSE16710_null_sampling_only_20261014.R")),
  R_version = as.character(getRversion()), RNG_kind = RNGkind(), collation = Sys.getlocale("LC_COLLATE"),
  sampling_seed = 20260831L, family_reference_seed = 20260832L,
  random_sets_replayed = 10000L, exact_template_assignments = 10000L,
  exact_candidate_indices = sum(cmp$candidate_index_match), exact_entire_member_vectors = sum(cmp$entire_member_vector_match),
  membership_rows_compared = nrow(members), candidate_sets_per_draw = 5L,
  max_mean_export_error = max(cmp$mean_expression_error), max_delta_export_error = max(cmp$delta_export_error),
  family_reference_draws_compared = length(family), max_reference_export_error = max(ref_error),
  exploratory_family_tail = tail, plus_one_numerator = 1 + sum(family >= observed), plus_one_denominator = 10001L,
  matching_rule = "first minimum of five absolute overall mean-expression distances, exact twenty-bin composition; no response/model statistic used for matching",
  original_model_statistics_recomputed = FALSE, model_fit_or_BH_or_gate_run = FALSE,
  primary_S1_state_unchanged = "COMPLETE_ZERO_PASS", biological_independence_authenticated = FALSE,
  prospective_outcome_protection_authenticated = FALSE, pathway_specific_inference_permitted = FALSE,
  full_label_null_or_full_multiverse_or_holdout_verified = FALSE, whole_goal_complete = FALSE,
  output_identities = lapply(outputs, function(p) list(path = p, bytes = as.numeric(file.info(p)$size), sha256 = sha(p))))
jsonlite::write_json(report, file.path(out, "SAMPLING_REPLAY_RECEIPT.json"), pretty = TRUE,
                     auto_unbox = TRUE, digits = 16, na = "null")
cat("COMPLETE fixed-seed memberships and family-reference replay; no scientific fit or primary change\n")
