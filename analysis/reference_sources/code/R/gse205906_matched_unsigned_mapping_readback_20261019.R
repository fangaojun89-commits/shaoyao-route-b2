# Functions only: independent base-R reconstruction from caller-held objects.
# No files, model fits, random calls or changes to supplied objects. The saved
# native checkpoint contains no member indices: NGenes/domain agreement cannot
# authenticate exactly which members were consumed by a native statistical call.

g205u_need <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
}

g205u_hash <- function(value) digest::digest(value, "sha256", serializeVersion = 3L)

g205u_reconstruct <- function(source, rat_map, annotation, collection) {
  source <- as.data.frame(source)
  g205u_need(all(c("term_id", "term_name", "gene_symbol") %in% names(source)),
    "Unsigned resource columns missing")
  source <- unique(data.frame(program_id = as.character(source$term_id),
    program_name = as.character(source$term_name),
    human_symbol = as.character(source$gene_symbol), stringsAsFactors = FALSE))
  g205u_need(!anyNA(source$program_id) && all(nzchar(source$program_id)),
    "Unsigned source program keys invalid")
  source_names <- unique(source[c("program_id", "program_name")])
  g205u_need(!anyDuplicated(source_names$program_id), "Ambiguous program display names")
  map_index <- match(source$human_symbol, rat_map$human_symbol)
  mapped <- source[!is.na(map_index), c("program_id", "program_name"), drop = FALSE]
  mapped$symbol <- rat_map$rat_symbol[map_index[!is.na(map_index)]]
  mapped <- unique(mapped)
  ids <- sort(unique(mapped$program_id))
  members <- indices <- list()
  coverage <- data.frame(program_id = character(), program_name = character(),
    mapped_symbol_n = integer(), tested_feature_n = integer(), stringsAsFactors = FALSE)
  for (id in ids) {
    rows <- mapped[mapped$program_id == id, , drop = FALSE]
    symbols <- unique(rows$symbol)
    retained <- intersect(symbols, as.character(annotation$symbol))
    coverage[nrow(coverage) + 1L, ] <- list(id, rows$program_name[[1L]],
      length(symbols), length(retained))
    # Eligibility uses unique symbols, whereas CAMERA/mroast use all gene rows.
    if (length(retained) >= 10L && length(retained) <= 500L) {
      members[[id]] <- retained
      indices[[id]] <- which(annotation$symbol %in% retained)
    }
  }
  coverage$collection <- rep(collection, nrow(coverage))
  coverage$source_feature_n <- coverage$mapped_symbol_n
  coverage$eligible_10_500 <- coverage$tested_feature_n >= 10L & coverage$tested_feature_n <= 500L
  coverage$tested <- coverage$program_id %in% names(members)
  coverage$confidence_domain <- rep("NOT_APPLICABLE", nrow(coverage))
  coverage$mapping_rule <- rep(
    "NCBI_Datasets_REST_v18.35.0_snapshot_2026-08-18_strict_one_to_one", nrow(coverage))
  list(collection = collection, source_program_ids = sort(unique(source$program_id)),
    members = members, indices = indices, coverage = coverage)
}

# This component verifies mapping consequences, not upstream source-file custody,
# normalization, statistical P values, or native consumption of exact indices.
g205u_readback <- function(checkpoint, gene_resources, orthology_table) {
  g205u_need(identical(checkpoint$schema, "GSE205906_MATCHED_SAVED_NATIVE_CHECKPOINT/v1") &&
    isTRUE(checkpoint$statistical_fit_executed), "New fitted checkpoint required")
  packet <- checkpoint$prepared_packet
  result <- checkpoint$result
  g205u_need(identical(result$schema, "gse205906-matched-downstream-result-v1") &&
    identical(result$binding, packet$binding) &&
    identical(result$prepared_sha256, packet$sha256) &&
    identical(result$statistical_source_pins, packet$source_pins) &&
    identical(result$candidate_admission, FALSE), "Prepared/result association differs")
  g205u_need(is.list(gene_resources) &&
    all(c("reactome", "syngo_sensitivity") %in% names(gene_resources)) &&
    identical(g205u_hash(gene_resources), packet$binding$held_object_hashes$gene_resources) &&
    identical(g205u_hash(orthology_table), packet$binding$held_object_hashes$orthology_table),
    "Held resource/orthology binding differs")

  counts <- packet$counts
  g205u_need(is.list(packet$annotation) && !is.null(packet$annotation$annotation),
    "Validated annotation packet missing")
  annotation <- as.data.frame(packet$annotation$annotation)
  base <- result$native$base_fit
  keep <- base$keep
  g205u_need(is.matrix(counts) && is.numeric(counts) &&
    identical(dim(counts), c(32883L, 9L)) &&
    identical(colnames(counts), paste0("GSM", 6235292:6235300)) &&
    !is.null(rownames(counts)) && !anyNA(rownames(counts)) &&
    all(nzchar(rownames(counts))) && !anyDuplicated(rownames(counts)) &&
    all(c("gene_id", "symbol") %in% names(annotation)) &&
    identical(as.character(annotation$gene_id), rownames(counts)) &&
    is.logical(keep) && length(keep) == nrow(counts) && !anyNA(keep) &&
    identical(base$y$counts, counts[keep, , drop = FALSE]), "Retained gene-row domain differs")
  retained_annotation <- annotation[keep, , drop = FALSE]
  orthology <- as.data.frame(orthology_table)
  g205u_need(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in%
    names(orthology)), "Orthology columns missing")
  selected <- which(orthology$group_method == "NCBI Ortholog" & orthology$other_tax_id == 10116L)
  rat_map <- unique(data.frame(human_symbol = as.character(orthology$human_symbol[selected]),
    rat_symbol = as.character(orthology$other_symbol[selected]), stringsAsFactors = FALSE))
  g205u_need(nrow(rat_map) > 0L && !anyNA(rat_map) &&
    !anyDuplicated(rat_map$human_symbol) && !anyDuplicated(rat_map$rat_symbol),
    "Selected rat mapping must be strict one-to-one")

  collections <- c(Reactome = "reactome", SynGO_sensitivity = "syngo_sensitivity")
  reconstructed <- lapply(names(collections), function(collection) {
    g205u_reconstruct(gene_resources[[collections[[collection]]]], rat_map,
      retained_annotation, collection)
  })
  names(reconstructed) <- names(collections)
  saved_coverage <- as.data.frame(result$mapping_coverage)
  g205u_need("collection" %in% names(saved_coverage) && !anyNA(saved_coverage$collection) &&
    all(saved_coverage$collection %in% c(names(collections), "PROGENy", "DoRothEA_AB")),
    "Coverage collection domain invalid")
  for (collection in names(collections)) {
    expected <- reconstructed[[collection]]$coverage
    observed <- saved_coverage[saved_coverage$collection == collection, , drop = FALSE]
    g205u_need(all(names(expected) %in% names(observed)) &&
      !anyNA(observed$program_id) && !anyDuplicated(observed$program_id) &&
      setequal(as.character(observed$program_id), expected$program_id),
      "Mapped coverage program domain differs")
    observed <- observed[match(expected$program_id, observed$program_id), , drop = FALSE]
    for (field in names(expected)) {
      actual <- observed[[field]]
      wanted <- expected[[field]]
      same <- if (is.numeric(wanted)) {
        is.numeric(actual) && !anyNA(actual) && identical(as.numeric(actual), as.numeric(wanted))
      } else identical(actual, wanted)
      g205u_need(same, paste0("Unsigned coverage field differs: ", field))
    }
  }
  contrasts <- c("local_disease", "intervention")
  for (table_name in c("unsigned_camera", "unsigned_mroast", "unsigned_fgsea")) {
    table <- as.data.frame(result$tables[[table_name]])
    g205u_need(all(c("collection", "contrast", "program_id") %in% names(table)) &&
      !anyNA(table[c("collection", "contrast", "program_id")]) &&
      all(table$collection %in% names(collections)) && all(table$contrast %in% contrasts) &&
      !anyDuplicated(table[c("collection", "contrast", "program_id")]),
      "Unsigned result key/domain invalid")
    for (collection in names(collections)) for (contrast in contrasts) {
      rows <- table[table$collection == collection & table$contrast == contrast, , drop = FALSE]
      expected_indices <- reconstructed[[collection]]$indices
      expected_ids <- names(expected_indices)
      g205u_need(all(rows$program_id %in% expected_ids), "Returned ineligible unsigned program")
      if (table_name != "unsigned_fgsea") {
        g205u_need(setequal(as.character(rows$program_id), expected_ids) &&
          "NGenes" %in% names(rows) && is.numeric(rows$NGenes) && !anyNA(rows$NGenes) &&
          identical(as.numeric(rows$NGenes),
            as.numeric(lengths(expected_indices[as.character(rows$program_id)]))),
          "Complete CAMERA/mroast domain or gene-row count differs")
      }
    }
  }
  list(status = "HELD_NEW_CHECKPOINT_UNSIGNED_MAPPING_COMPONENT_VERIFIED",
    reconstructed = reconstructed, retained_gene_rows = nrow(retained_annotation),
    own_mapping_and_coverage_verified = TRUE, camera_mroast_domain_and_NGenes_verified = TRUE,
    fgsea_returned_subset_verified = TRUE, exact_native_consumed_indices_authenticated = FALSE,
    statistical_P_values_verified = FALSE, external_source_file_custody_verified = FALSE,
    model_fits = 0L, candidate_admission = FALSE, full_independent_numeric_acceptance = FALSE)
}
