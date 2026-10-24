# Functions only: native result packing, exact RDS transport and streaming tails.
# This code runs no model and grants no scientific execution authority.

g128r_need <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

g128r_families <- function() {
  contrasts <- c("chronic_combined", "D28_vs_control", "D42_vs_control", "acute_combined",
                 "D4_vs_control", "D7_vs_control", "D14_vs_control")
  collections <- c("DoRothEA_AB", "PROGENy", "Reactome", "SynGO_assignedBy_proxy_sensitivity")
  rows <- list()
  for (method in c("edgeR", "voom", "DESeq2", "activity", "CAMERA", "mroast", "fgsea")) {
    groups <- if (method %in% c("edgeR", "voom", "DESeq2")) "GENE" else {
      if (method == "activity") collections else collections[3:4]
    }
    tests <- if (method %in% c("mroast", "fgsea")) "chronic_combined" else contrasts
    for (collection in groups) for (contrast in tests) {
      rows[[length(rows) + 1L]] <- data.frame(method = method, collection = collection,
                                             contrast = contrast, stringsAsFactors = FALSE)
    }
  }
  result <- do.call(rbind, rows)
  result$family_id <- with(result, paste(method, collection, contrast, sep = "::"))
  rownames(result) <- NULL
  result
}

g128r_project_q <- function(p_value, method) {
  if (method == "mroast") {
    return(pmax(p_value, stats::p.adjust(pmax(p_value - 0.5 / 9999L, 0), "BH")))
  }
  stats::p.adjust(p_value, "BH")
}

g128r_valid_statistic <- function(statistic, native_p, method) {
  # Only native -log10(P) transforms have a permitted positive infinity.
  is.finite(statistic) | (method %in% c("CAMERA", "mroast") &
    is.infinite(statistic) & statistic > 0 & !is.na(native_p) & native_p == 0)
}

g128r_validate_packed <- function(packed) {
  key_columns <- c("method", "collection", "contrast", "feature_id", "family_id")
  family_columns <- c("method", "collection", "contrast", "family_id")
  g128r_need(is.list(packed) && is.data.frame(packed$keys) &&
               identical(names(packed$keys), key_columns) &&
               all(vapply(packed$keys, is.character, logical(1))) &&
               !anyNA(packed$keys) && !anyDuplicated(packed$keys) &&
               all(nzchar(trimws(unlist(packed$keys, use.names = FALSE)))),
             "Invalid packed feature keys.")
  g128r_need(is.data.frame(packed$families) &&
               identical(names(packed$families), family_columns) &&
               all(vapply(packed$families, is.character, logical(1))) &&
               !anyNA(packed$families) && !anyDuplicated(packed$families$family_id),
             "Invalid packed family template.")
  full <- g128r_families()
  allowed <- full[match(packed$families$family_id, full$family_id), , drop = FALSE]
  rownames(allowed) <- rownames(packed$families) <- NULL
  g128r_need(nrow(packed$families) > 0L && identical(allowed, packed$families) &&
               all(packed$keys$family_id %in% packed$families$family_id) &&
               identical(with(packed$keys, paste(method, collection, contrast, sep = "::")),
                         packed$keys$family_id), "Packed family identity differs.")
  n <- nrow(packed$keys)
  g128r_need(is.double(packed$statistic) && is.double(packed$native_p) &&
               is.null(dim(packed$statistic)) && is.null(dim(packed$native_p)) &&
               length(packed$statistic) == n && length(packed$native_p) == n &&
               all(is.na(packed$native_p) |
                     (is.finite(packed$native_p) & packed$native_p >= 0 & packed$native_p <= 1)),
             "Packed statistic/P lengths, precision or probabilities differ.")
  invisible(TRUE)
}

g128r_pack_bundle <- function(bundle, target_contrast = NULL) {
  template <- g128r_families()
  tables <- c(bundle$gene_camera$tables[c("edgeR", "voom", "DESeq2")],
              list(activity = bundle$activity$native,
                   CAMERA = bundle$gene_camera$tables$CAMERA,
                   mroast = bundle$unsigned$mroast, fgsea = bundle$unsigned$fgsea))
  g128r_need(identical(names(tables), unique(template$method)), "Bundle method set/order differs.")
  rows <- list()
  statistics <- list()
  probabilities <- list()
  q_values <- list()
  mixed <- list()
  for (method in names(tables)) {
    tab <- as.data.frame(tables[[method]], stringsAsFactors = FALSE)
    if (nrow(tab) == 0L) next # Empty families are represented explicitly below.
    is_gene <- method %in% c("edgeR", "voom", "DESeq2")
    feature_column <- if (is_gene) "gene_id" else "program_id"
    p_column <- if (method == "fgsea") "pval" else "PValue"
    required <- c(feature_column, "contrast", p_column, "BH_q")
    if (!is_gene) required <- c(required, "collection")
    g128r_need(all(required %in% names(tab)), paste("Native table schema:", method))
    collection <- if (is_gene) rep("GENE", nrow(tab)) else as.character(tab$collection)
    feature <- as.character(tab[[feature_column]])
    contrast <- as.character(tab$contrast)
    g128r_need(!anyNA(c(collection, feature, contrast)) &&
                 all(nzchar(trimws(c(collection, feature, contrast)))), "Missing native result key.")
    family <- paste(method, collection, contrast, sep = "::")
    g128r_need(all(family %in% template$family_id), "Unexpected native family/contrast.")
    p_value <- as.numeric(tab[[p_column]])
    q_value <- as.numeric(tab$BH_q)
    g128r_need(is.numeric(tab[[p_column]]) && is.numeric(tab$BH_q), "Native P/q must be numeric.")
    g128r_need(all(is.na(p_value) | (is.finite(p_value) & p_value >= 0 & p_value <= 1)) &&
                 all(is.na(q_value) | (is.finite(q_value) & q_value >= 0 & q_value <= 1)),
               "Native P/q out of range; never repair or clip.")
    statistic_column <- switch(method, edgeR = "QLF_F", voom = "t",
                                DESeq2 = "WaldStatistic", fgsea = "NES", NULL)
    if (!is.null(statistic_column)) {
      g128r_need(statistic_column %in% names(tab) && is.numeric(tab[[statistic_column]]),
                 paste("Native statistic must be numeric:", method))
    }
    statistic <- switch(method,
      edgeR = as.numeric(tab$QLF_F), voom = abs(as.numeric(tab$t)),
      DESeq2 = abs(as.numeric(tab$WaldStatistic)),
      activity = {
        matrix_t <- bundle$activity$signed_statistic
        all_contrasts <- unique(template$contrast)
        g128r_need(is.matrix(matrix_t) && identical(colnames(matrix_t), all_contrasts),
                   "Activity moderated-t contrast order differs.")
        first <- which(contrast == all_contrasts[[1L]])
        expected_keys <- paste(collection[first], feature[first], sep = "::")
        g128r_need(identical(rownames(matrix_t), expected_keys) &&
                     length(matrix_t) == nrow(tab), "Activity moderated-t dimensions/keys differ.")
        for (test in all_contrasts) {
          selected <- which(contrast == test)
          g128r_need(identical(paste(collection[selected], feature[selected], sep = "::"), expected_keys),
                     "Activity native feature order differs between contrasts.")
        }
        g128r_need(identical(contrast, rep(all_contrasts, each = nrow(matrix_t))) &&
                     is.numeric(matrix_t), "Activity native table must match matrix contrast-major order.")
        abs(as.numeric(matrix_t))
      },
      CAMERA = -log10(p_value), mroast = -log10(p_value), fgsea = abs(as.numeric(tab$NES))
    )
    g128r_need(length(statistic) == nrow(tab), paste("Missing native statistic:", method))
    key_rows <- data.frame(method = method, collection = collection, contrast = contrast,
                            feature_id = feature, family_id = family, stringsAsFactors = FALSE)
    g128r_need(!anyDuplicated(key_rows), "Repeated feature within native family.")
    for (family_name in unique(family)) {
      indices <- which(family == family_name)
      expected_q <- g128r_project_q(p_value[indices], method)
      g128r_need(identical(q_value[indices], as.numeric(expected_q)),
                 paste("Native project BH differs from source rule:", family_name))
    }
    if (method == "mroast") {
      extra <- c("FDR", "PValue.Mixed", "FDR.Mixed")
      g128r_need(all(extra %in% names(tab)), "Native mixed/directional mroast fields missing.")
      g128r_need(all(vapply(tab[extra], is.numeric, logical(1))) &&
                   all(vapply(tab[extra], function(p) all(is.na(p) |
                     (is.finite(p) & p >= 0 & p <= 1)), logical(1))),
                 "Native auxiliary mroast probabilities out of range.")
      mixed[[method]] <- data.frame(key_rows,
        native_fdr = as.numeric(tab$FDR), mixed_p = as.numeric(tab$PValue.Mixed),
        mixed_fdr = as.numeric(tab$FDR.Mixed), stringsAsFactors = FALSE)
    }
    rows[[method]] <- key_rows
    statistics[[method]] <- statistic
    probabilities[[method]] <- p_value
    q_values[[method]] <- q_value
  }
  g128r_need(length(rows) > 0L, "Entire native bundle is empty.")
  keys <- do.call(rbind, rows)
  rownames(keys) <- NULL
  statistic <- unlist(statistics, use.names = FALSE)
  native_p <- unlist(probabilities, use.names = FALSE)
  native_q <- unlist(q_values, use.names = FALSE)
  g128r_need(!anyDuplicated(keys), "Bundle has duplicate keys.")
  if (!is.null(target_contrast)) {
    g128r_need(is.character(target_contrast) && length(target_contrast) == 1L &&
                 target_contrast %in% template$contrast && target_contrast != "chronic_combined" &&
                 target_contrast != "acute_combined", "Exact companion must target a fixed day/control contrast.")
    selected <- keys$contrast == target_contrast
    keys <- keys[selected, , drop = FALSE]
    statistic <- statistic[selected]
    native_p <- native_p[selected]
    native_q <- native_q[selected]
    template <- template[template$contrast == target_contrast, , drop = FALSE]
    mixed <- list() # No unsigned day contrasts were originally run.
    rownames(keys) <- NULL
    rownames(template) <- NULL
  }
  eligible <- g128r_valid_statistic(statistic, native_p, keys$method)
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(template)), function(index) {
    family <- template$family_id[[index]]
    selected <- which(keys$family_id == family)
    p <- native_p[selected]
    q <- native_q[selected]
    t <- statistic[selected]
    data.frame(template[index, , drop = FALSE], n = length(selected),
      available_statistic_n = sum(eligible[selected]),
      native_p_available_n = sum(!is.na(p)), native_q_available_n = sum(!is.na(q)),
      native_p_lt_05_n = sum(p < .05, na.rm = TRUE),
      native_q_le_05_n = sum(q <= .05, na.rm = TRUE),
      native_q_any_le_05 = any(q <= .05, na.rm = TRUE),
      native_p_zero_n = sum(p == 0, na.rm = TRUE),
      positive_infinite_statistic_n = sum(is.infinite(t) & t > 0),
      negative_edge_F_n = if (template$method[[index]] == "edgeR") sum(t < 0, na.rm = TRUE) else 0L,
      stringsAsFactors = FALSE)
  }))
  rownames(diagnostics) <- NULL
  packed <- list(keys = keys, statistic = as.numeric(statistic), native_p = as.numeric(native_p),
       families = template, native_diagnostics = diagnostics,
       mroast_auxiliary = if (length(mixed)) mixed[[1L]] else NULL,
       rng_calls = bundle$unsigned$rng_calls,
       scope = "native statistic/P packing only; no source/fit/biological acceptance")
  g128r_validate_packed(packed)
  packed
}

g128r_new_reference <- function(observed, null_draws, kind) {
  g128r_need(length(null_draws) == 1L && is.finite(null_draws) && null_draws >= 1L &&
               null_draws <= .Machine$integer.max - 1L &&
               null_draws == as.integer(null_draws), "Invalid required draw count.")
  g128r_need(identical(kind, "sampled") || identical(kind, "exhaustive"), "Unknown reference kind.")
  g128r_validate_packed(observed)
  valid <- g128r_valid_statistic(observed$statistic, observed$native_p, observed$keys$method)
  list(observed = observed, kind = kind, total = as.integer(null_draws + 1L),
       seen = c(TRUE, rep(FALSE, null_draws)), observed_valid = valid,
       hits_including_observed = as.integer(valid), missing_null_n = integer(length(valid)))
}

g128r_add_reference <- function(state, draw, allocation_id) {
  g128r_validate_packed(draw)
  g128r_need(length(allocation_id) == 1L && is.finite(allocation_id) &&
               allocation_id == as.integer(allocation_id) && allocation_id >= 2L &&
               allocation_id <= state$total && !state$seen[[allocation_id]],
             "Repeated/observed/out-of-range allocation ID.")
  g128r_need(identical(draw$keys, state$observed$keys) &&
               identical(draw$families, state$observed$families), "Draw feature/family domain changed.")
  valid <- g128r_valid_statistic(draw$statistic, draw$native_p, draw$keys$method)
  observed <- state$observed$statistic
  finite_observed <- is.finite(observed) & state$observed_valid
  extreme <- rep(FALSE, length(observed))
  selected <- which(finite_observed & valid)
  extreme[selected] <- draw$statistic[selected] >=
    observed[selected] - 1e-10 * (1 + abs(observed[selected]))
  selected <- which(state$observed_valid & is.infinite(observed) & observed > 0 & valid)
  extreme[selected] <- is.infinite(draw$statistic[selected]) & draw$statistic[selected] > 0
  state$hits_including_observed <- state$hits_including_observed + as.integer(extreme)
  state$missing_null_n <- state$missing_null_n + as.integer(!valid)
  state$seen[[allocation_id]] <- TRUE
  state
}

g128r_finish_reference <- function(state) {
  g128r_need(all(state$seen), "Required allocation coverage incomplete; no partial P/q.")
  complete <- state$observed_valid & state$missing_null_n == 0L
  p_value <- rep(NA_real_, length(complete))
  p_value[complete] <- state$hits_including_observed[complete] / state$total
  q_value <- rep(NA_real_, length(complete))
  for (family in state$observed$families$family_id) {
    indices <- which(state$observed$keys$family_id == family)
    q_value[indices] <- stats::p.adjust(p_value[indices], "BH")
  }
  data.frame(state$observed$keys, empirical_p = p_value, empirical_BH_q = q_value,
    hits_including_observed = state$hits_including_observed,
    required_total_n = state$total, missing_null_n = state$missing_null_n,
    observed_unavailable = !state$observed_valid,
    available_total_n = state$total - state$missing_null_n - as.integer(!state$observed_valid),
    reference_kind = state$kind, mc_uncertainty_computed = FALSE,
    stringsAsFactors = FALSE, check.names = FALSE)
}

g128r_write_archive <- function(packed, path, run_root, identity) {
  g128r_validate_packed(packed)
  root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
  g128r_need(grepl("^([Rr]:/Shaoyao_storage/cache/GSE128623_[A-Za-z0-9_-]+|/root/autodl-tmp/gse128623_label_reference_[A-Za-z0-9_-]+)$", root),
             "Use a specific owned GSE128623 R-cache/cloud run root.")
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  g128r_need(identical(parent, root) || startsWith(parent, paste0(root, "/")),
             "Archive is outside the explicit run root.")
  g128r_need(grepl("^[A-Za-z0-9_-]+[.]rds$", basename(path)), "Invalid archive filename.")
  final <- file.path(parent, basename(path))
  partial <- paste0(final, ".partial")
  g128r_need(!file.exists(final) && !file.exists(partial), "Archive/partial exists; no overwrite.")
  required <- c("run_id", "cell_id", "reference_space", "allocation_id", "input_manifest_sha256",
                "source_manifest_sha256", "allocation_row_sha256")
  g128r_need(is.list(identity) && all(required %in% names(identity)), "Incomplete archive coupling identity.")
  for (field in required[grepl("sha256$", required)]) {
    g128r_need(is.character(identity[[field]]) && length(identity[[field]]) == 1L &&
                 !is.na(identity[[field]]) && grepl("^[0-9a-f]{64}$", identity[[field]]),
               paste("Invalid archive coupling hash:", field))
  }
  for (field in c("run_id", "cell_id", "reference_space")) {
    g128r_need(is.character(identity[[field]]) && length(identity[[field]]) == 1L &&
                 !is.na(identity[[field]]) && nzchar(identity[[field]]), "Invalid archive run/cell/space identity.")
  }
  g128r_need(length(identity$allocation_id) == 1L && is.finite(identity$allocation_id) &&
               identity$allocation_id == as.integer(identity$allocation_id) && identity$allocation_id >= 1L,
             "Invalid archive allocation identity.")
  payload <- list(schema_version = "gse128623-statistic-native-p-archive-v1", identity = identity,
                   packed = packed, real_execution_authority_granted = FALSE)
  # Lossless compression avoids repeating large character key vectors on disk.
  # Numerical values remain 64-bit doubles and are checked exactly after readback.
  saveRDS(payload, partial, compress = "gzip", version = 3L)
  # Native S4 fits can retain closures: compare their complete serialized
  # contents, not process-local environment addresses after restoration.
  g128r_need(identical(serialize(readRDS(partial), NULL, version = 3L),
    serialize(payload, NULL, version = 3L)), "RDS payload roundtrip mismatch.")
  # Exclusive same-filesystem hardlink publication cannot replace an old file.
  g128r_need(file.link(partial, final), "Exclusive archive publication failed; preserve partial.")
  g128r_need(identical(serialize(readRDS(final), NULL, version = 3L),
    serialize(payload, NULL, version = 3L)), "Published archive differs.")
  cleanup_status <- unlink(partial, recursive = FALSE)
  list(file = final, bytes = unname(file.info(final)$size),
       sha256 = digest::digest(file = final, algo = "sha256"),
       own_partial_removed = identical(cleanup_status, 0L))
}
