# SOURCE-only new-input program binding; no old public-object reads or model run.
# Object consistency is not input custody, runtime acceptance or model authority.

gse205906_ens101_program_hash_v2 <- function(x) {
  digest::digest(serialize(x, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_program_contract_v2 <- function(prefit, mapping, source_bundle,
                                               input_identity_function) {
  need <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
  hash <- gse205906_ens101_program_hash_v2
  member_hash <- function(x) digest::digest(charToRaw(enc2utf8(paste(x, collapse = "\n"))),
                                            "sha256", serialize = FALSE)
  need(requireNamespace("digest", quietly = TRUE) && requireNamespace("data.table", quietly = TRUE),
       "ENS101_CONTRACT_EXISTING_DEPENDENCIES_REQUIRED")
  pins <- list(
    binding = list(bytes = 15697L, sha256 = "2aaa386df70c0031933c244eafa895d63c5dbd7a7be77277cfd4a57aceba538d"),
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    helpers = list(bytes = 175131L, sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"))
  need(is.list(source_bundle) && identical(names(source_bundle), names(pins)), "ENS101_CONTRACT_SOURCE_ROLES")
  for (role in names(pins)) need(is.raw(source_bundle[[role]]) &&
    length(source_bundle[[role]]) == pins[[role]]$bytes &&
    identical(digest::digest(source_bundle[[role]], "sha256", serialize = FALSE), pins[[role]]$sha256),
    "ENS101_CONTRACT_SOURCE_IDENTITY")
  sha_ok <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[0-9a-f]{64}$", x)
  checked_body <- function(x) {
    need(is.list(x) && sha_ok(x$output_object_sha256), "ENS101_CONTRACT_OBJECT_PIN_REQUIRED")
    declared <- x$output_object_sha256; x$output_object_sha256 <- NULL
    need(identical(declared, hash(x)), "ENS101_CONTRACT_OBJECT_DRIFT")
  }
  checked_body(prefit); checked_body(mapping)
  need(identical(prefit$schema, "GSE205906_ENS101_EFFECT_INDEPENDENT_PREFIT/v1") &&
       identical(prefit$status, "PASS_TEN_DESIGN_FILTER_TMM_ONLY") &&
       identical(mapping$schema, "GSE205906_ENS101_MEMBERSHIP/v2") &&
       identical(mapping$status, "PASS_TEN_MEMBERSHIP_ONLY") &&
       identical(prefit$production_go, FALSE) && identical(prefit$new_model_authorized, FALSE) &&
       identical(mapping$production_go, FALSE) && identical(mapping$new_model_authorized, FALSE) &&
       identical(mapping$prefit_output_object_sha256, prefit$output_object_sha256),
       "ENS101_CONTRACT_CURRENT_PHASES_NOT_READY")
  need(identical(prefit$input_reference, list(assembly = "Rnor_6.0", annotation = "Ensembl101",
    expected_gene_rows = 32883L, primary_strand_setting = 0L, count_unit = "gene_fragments")),
    "ENS101_CONTRACT_NEW_INPUT_DOMAIN")
  ids <- paste0("GSM", 6235292:6235300)
  phases <- c("NONE", ids)
  families <- c("Reactome", "SynGO_sensitivity", "PROGENy", "DoRothEA_AB")
  need(identical(names(prefit$runs), phases) && identical(names(mapping$runs), phases),
       "ENS101_CONTRACT_PHASE_KEYS")
  universe <- prefit$counts_gene_ids
  need(is.character(universe) && length(universe) == 32883L && !anyNA(universe) &&
       !anyDuplicated(universe) && all(grepl("^ENSRNOG[0-9]{11}$", universe)) &&
       identical(universe, sort(universe, method = "radix")), "ENS101_CONTRACT_NEW_GENE_UNIVERSE")
  packets <- lapply(seq_along(phases), function(i) {
    p <- prefit$runs[[i]]; m <- mapping$runs[[i]]
    samples <- if (i == 1L) ids else ids[-(i - 1L)]
    need(identical(p$omission, phases[[i]]) && identical(m$omission, phases[[i]]) &&
         identical(p$sample_ids, samples) && identical(m$sample_ids, samples) &&
         identical(p$preprocessing$status, "PASS_DESIGN_FILTER_TMM_ONLY") &&
         identical(m$status, "PASS_MEMBERSHIP_ONLY"), "ENS101_CONTRACT_PHASE_ASSOCIATION")
    pbody <- p; pbody$phase_output_object_sha256 <- NULL
    mbody <- m; mbody$mapping_phase_object_sha256 <- NULL
    need(identical(p$phase_output_object_sha256, hash(pbody)) &&
         identical(m$prefit_phase_object_sha256, p$phase_output_object_sha256) &&
         identical(m$mapping_phase_object_sha256, hash(mbody)),
         "ENS101_CONTRACT_PHASE_PREFIT_IDENTITY")
    genes <- p$retained_gene_ids
    need(is.character(genes) && !anyNA(genes) && !anyDuplicated(genes) && length(genes) >= 1000L &&
         identical(genes, universe[universe %in% genes]) && identical(m$retained_gene_ids, genes) &&
         identical(m$gene_universe_sha256, hash(genes)) &&
         identical(m$retained_symbols, p$retained_annotation$symbol) &&
         identical(m$eligibility_symbols, sort(unique(p$retained_annotation$symbol), method = "radix")),
         "ENS101_CONTRACT_RETAINED_GENES")
    need(identical(names(m$families), families) && identical(names(m$family_member_pins), families) &&
         identical(names(m$eligible_program_n), families), "ENS101_CONTRACT_FAMILIES")
    out <- lapply(families, function(collection) {
      f <- m$families[[collection]]; mp <- m$family_member_pins[[collection]]
      need(is.list(f) && is.list(mp) && identical(names(f), names(mp)) &&
           (length(f) == 0L || (!is.null(names(f)) && !anyNA(names(f)) && all(nzchar(names(f))) &&
             !anyDuplicated(names(f)))) && identical(m$eligible_program_n[[collection]], length(f)),
           "ENS101_CONTRACT_MEMBER_IDS")
      lapply(names(f), function(id) {
        symbols <- f[[id]]
        need(is.character(symbols) && !anyNA(symbols) && all(nzchar(symbols)) &&
             !anyDuplicated(symbols) && length(symbols) >= 10L && length(symbols) <= 500L &&
             identical(symbols, sort(symbols, method = "radix")) &&
             all(symbols %in% m$eligibility_symbols) &&
             identical(mp[[id]]$matched_target_n, length(symbols)) &&
             identical(mp[[id]]$matched_targets_sha256, member_hash(symbols)),
             "ENS101_CONTRACT_MEMBER_PIN")
        list(n = length(symbols), sha256 = member_hash(symbols))
      }) |> stats::setNames(names(f))
    }); names(out) <- families
    list(gene_ids = genes, families = out)
  }); names(packets) <- phases
  need(identical(mapping$full_eligible_program_n, as.list(lengths(packets$NONE$families))),
       "ENS101_CONTRACT_FULL_CARDINALITY_DRIFT")
  input_identity <- c(prefit$input_object_identities[c("count_matrix", "gene_annotation", "sample_ledger")],
    list(gene_resources = mapping$resource_bundle_object_sha256,
         orthology_table = mapping$orthology_object_sha256))
  need(identical(names(input_identity), c("count_matrix", "gene_annotation", "sample_ledger",
    "gene_resources", "orthology_table")) && all(vapply(input_identity, sha_ok, logical(1))) &&
    is.function(input_identity_function), "ENS101_CONTRACT_CURRENT_INPUT_IDENTITIES")
  spec <- list(domain = "PRODUCTION", sample_ids = ids, input_sha256 = input_identity,
    phases = packets, new_namespace = "RNOR6_ENS101_PROGRAM_CONTRACT_V2_20261006",
    current_prefit_object_sha256 = hash(prefit), current_mapping_object_sha256 = hash(mapping),
    declaration_is_provenance_authentication = FALSE)
  text <- rawToChar(source_bundle$binding)
  replace_once <- function(text, old, new) {
    hits <- gregexpr(old, text, fixed = TRUE)[[1L]]
    need(length(hits) == 1L && hits[[1L]] > 0L, "ENS101_CONTRACT_EXACT_SOURCE_ANCHOR")
    paste0(substr(text, 1L, hits[[1L]] - 1L), new, substr(text, hits[[1L]] + nchar(old), nchar(text)))
  }
  text <- replace_once(text,
    "    output[[branch]]$program_domains$DoRothEA_AB$eligible_tested_program_n <- length(spec$phases$NONE$families$DoRothEA_AB)",
    paste0("    output[[branch]]$program_domains$PROGENy$eligible_tested_program_n <- length(spec$phases$NONE$families$PROGENy)\n",
      "    output[[branch]]$program_domains$DoRothEA_AB$eligible_tested_program_n <- length(spec$phases$NONE$families$DoRothEA_AB)"))
  text <- replace_once(text, "program_id=names(spec$phases$NONE$families[[collection]])",
    "program_id=as.character(names(spec$phases$NONE$families[[collection]]))")
  text <- replace_once(text, "out$collection <- collection;out", "out$collection <- rep(collection,nrow(out));out")
  text <- replace_once(text,
    "robustness <<- list(all_programs_pass=all(x$robustness_pass),failed_program_contrast_n=sum(!x$robustness_pass))",
    paste0("robustness <<- list(all_programs_pass=if(nrow(x)) all(x$robustness_pass) else NA,",
      "failed_program_contrast_n=sum(!x$robustness_pass),tested_program_contrast_n=nrow(x),",
      "status=if(nrow(x)) 'ESTIMABILITY_COUNTS_ONLY' else 'NOT_APPLICABLE_EMPTY_FULL_ACTIVITY_FAMILIES')"))
  # Evaluate ONLY the first two pinned pure function declarations; no legacy
  # public/synthetic wrappers, old Sep9 paths, filtering or model are evaluated.
  declarations <- as.list(parse(text = text, keep.source = FALSE))[1:2]
  expected <- c("gse205906_amend_model_program_binding_v1", "gse205906_program_contract_v1")
  private <- new.env(parent = baseenv())
  for (i in 1:2) {
    x <- declarations[[i]]
    need(is.call(x) && identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name(expected[[i]])) &&
         is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), "ENS101_CONTRACT_PURE_DECLARATIONS")
    eval(x, private)
  }
  base <- private$gse205906_program_contract_v1(spec, input_identity_function)
  make_session <- base$new_session
  base$new_session <- function() {
    session <- make_session()
    empty_allowed <- function(collection = NULL) {
      if (is.null(collection)) need(all(lengths(packets$NONE$families[c("PROGENy", "DoRothEA_AB")]) == 0L),
                                    "ENS101_CONTRACT_NOT_EMPTY_FULL_ACTIVITY")
      else need(collection %in% c("PROGENy", "DoRothEA_AB") && length(packets$NONE$families[[collection]]) == 0L,
                "ENS101_CONTRACT_NOT_EMPTY_COLLECTION")
    }
    session$empty_activity <- function(collection) {
      empty_allowed(collection)
      data.table::data.table(collection = character(), program_id = character(), contrast = character(),
        method = character(), effect = double(), se = double(), ci_low = double(), ci_high = double(),
        t = double(), df_total = double(), p_value = double(), q_within_collection_contrast = double())
    }
    session$empty_loao_summary <- function(long, expected_ids) {
      empty_allowed()
      need(is.data.frame(long) && nrow(long) == 0L && identical(expected_ids, ids) &&
        all(c("collection", "program_id", "contrast", "omitted_animal", "estimable", "estimate") %in% names(long)),
        "ENS101_CONTRACT_EMPTY_LONG_SCHEMA")
      data.table::data.table(collection = character(), program_id = character(), contrast = character(),
        n_omissions = integer(), n_estimable = integer(), estimate_min = double(), estimate_max = double(),
        estimate_mean = double(), sign_consistency = double())
    }
    session
  }
  old_amend <- base$amend_model
  base$amend_model <- function(model_text) {
    need(is.character(model_text) && length(model_text) == 1L &&
      identical(digest::digest(charToRaw(model_text), "sha256", serialize = FALSE), pins$model$sha256),
      "ENS101_CONTRACT_MODEL_SOURCE_PIN")
    model_text <- old_amend(model_text)
    model_text <- replace_once(model_text,
      "v2r4_fit_activity <- function(score, collection, design_matrix, contrasts) {",
      paste0("v2r4_fit_activity <- function(score, collection, design_matrix, contrasts) {\n",
        "  if (nrow(score) == 0L) {\n",
        "    v2r4_assert(identical(colnames(score), rownames(design_matrix)), 'EMPTY_ACTIVITY_SAMPLE_ORDER')\n",
        "    return(bound_programs$empty_activity(collection))\n  }"))
    model_text <- replace_once(model_text, "expected_programs <- lapply(full_activity$scores, rownames)",
      "expected_programs <- lapply(full_activity$scores, function(score) as.character(rownames(score)))")
    replace_once(model_text, "loao_summary <- v2r4_recompute_loao_summary(loao_long, ids)",
      paste0("loao_summary <- if (nrow(loao_long)) v2r4_recompute_loao_summary(loao_long, ids) else ",
             "bound_programs$empty_loao_summary(loao_long, ids)"))
  }
  base$namespace <- spec$new_namespace
  base$current_prefit_object_sha256 <- spec$current_prefit_object_sha256
  base$current_mapping_object_sha256 <- spec$current_mapping_object_sha256
  base$source_pins <- pins
  base$full_family_reporting <- lapply(base$full_family_n, function(n)
    list(tested_n = n, status = if(n == 0L) "NOT_APPLICABLE_EMPTY_FAMILY" else "MEMBERSHIP_ONLY_NOT_EFFECT_ACCEPTANCE"))
  base$producer_provenance_authenticated_here <- FALSE
  base$new_model_authorized <- FALSE
  base
}
