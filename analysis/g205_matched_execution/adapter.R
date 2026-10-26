# Functions only. Public-input binding is separate from historical custody.
# No original function is patched, relaxed, or assigned a replacement body.
g205p_need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

g205p_decode_gzip_rds <- function(raw) {
  g205p_need(is.raw(raw) && length(raw) > 2L && identical(raw[1:2], as.raw(c(31, 139))), "Pinned gzip RDS bytes required")
  con <- gzcon(rawConnection(raw)); on.exit(close(con))
  readRDS(con)
}

g205p_prepare <- function(core, sources, converted, gene_resources, orthology_table, binding) {
  e <- core$g205m_engine(sources)
  g205m_need <- core$g205m_need
  inputs <- converted$inputs
  g205p_need(identical(names(inputs), c("count_matrix", "gene_annotation", "sample_ledger", "fixed_gene_ids")) &&
    identical(converted$transformation$count_values_changed, FALSE), "Original unfitted input bridge schema required")
  g205p_need(identical(binding$schema, "G205_MATCHED_PUBLIC_INPUT_BINDING/v1") &&
    identical(binding$historical_native_custody_recreated, FALSE), "Explicit new public-input provenance required")
  fixed_gene_ids <- inputs$fixed_gene_ids
  g205p_need(identical(rownames(inputs$count_matrix), fixed_gene_ids) &&
    identical(colnames(inputs$count_matrix), as.character(inputs$sample_ledger$GEO_sample)), "Prepared input axes differ")
  # Exact scientific preparation block from g205m_prepare, after its input loader.
  projected <- e$gse205906_ens101_annotation_v1(converted$inputs$gene_annotation, fixed_gene_ids)
  g205m_need(is.list(gene_resources) && all(c("reactome", "progeny_human", "dorothea_rat", "syngo_sensitivity") %in% names(gene_resources)) &&
    is.data.frame(orthology_table), "Bound frozen resource/orthology objects required")
  orthology <- data.table::as.data.table(orthology_table)
  g205m_need(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in% names(orthology)), "Strict external orthology columns missing")
  e$rat_map <- unique(orthology[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
    .(human_symbol = as.character(human_symbol), rat_symbol = as.character(other_symbol))])
  g205m_need(nrow(e$rat_map) > 0L && !anyNA(e$rat_map) && !anyDuplicated(e$rat_map$human_symbol) &&
    !anyDuplicated(e$rat_map$rat_symbol), "Strict rat orthology not one-to-one")
  e$gene_sets <- gene_resources
  design <- e$v2r4_design(converted$inputs$sample_ledger, e$cfg)
  counts <- e$v2r4_validate_counts_original(converted$inputs$count_matrix, rownames(design$design))
  annotation <- e$v2r4_validate_annotation_gene(projected$computational_interface, rownames(counts), .8)
  # New metadata describes the actual bridge input; no fake historical manifest.
  packet <- list(schema = "G205_MATCHED_PUBLIC_PREPARED/v1", counts = counts,
    annotation = annotation, design = design, binding = binding,
    matrix_binding = list(schema = "G205_MATCHED_PUBLIC_MATRIX_BINDING/v1", prepared_input = binding$prepared_input,
      role = "STAR_S0_MATCHED_PRIMARY", historical_producer_QC_authenticated_here = FALSE),
    source_pins = core$g205m_pins(), input_conversion = converted$transformation,
    engine_sha256 = core$g205m_engine_hash(e), gate = e$cfg$gate,
    sample_unit = "NINE_SUBMITTED_LIBRARIES_NOT_AUTHENTICATED_ANIMALS",
    input_file_custody_authenticated_here = FALSE, old_models_scores_or_nulls_reused = FALSE)
  packet$sha256 <- core$g205m_hash(packet)
  list(packet = packet, engine = e)
}

g205p_runtime_check <- function(cfg) {
  g205p_need(identical(as.character(getRversion()), cfg$software$R), "Pinned R required")
  versions <- cfg$software[setdiff(names(cfg$software), c("R", "PowerShell"))]
  versions$BiocParallel <- "1.46.0"
  observed <- lapply(names(versions), function(p) as.character(utils::packageVersion(p)))
  names(observed) <- names(versions)
  g205p_need(identical(observed, versions), "Pinned statistical packages required")
  g205p_need(all(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "R_DATATABLE_NUM_THREADS")) == "1"),
    "Explicit single-thread environment required")
  list(R = as.character(getRversion()), platform = R.version$platform, packages = observed,
    RNG_before = RNGkind(), historical_runtime_authenticated = FALSE)
}
