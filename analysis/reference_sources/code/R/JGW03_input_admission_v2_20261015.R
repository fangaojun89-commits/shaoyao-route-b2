# Functions only. Defines a versioned input adapter for the preserved v1 kernel.
# Sourcing this file reads no scientific inputs and runs no model or bridge.
# Calling its build/fit functions still requires separately approved execution.

jgw03_v2_need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

jgw03_v2_spec <- function(accession) {
  jgw03_v2_need(is.character(accession) && length(accession) == 1L &&
    !is.na(accession) && accession %in% c("GSE289929", "GSE227494", "GSE299063"),
    "One of the three frozen JGW03 accessions is required")
  switch(accession,
    GSE289929 = list(organism = "Rattus norvegicus", cell = "PC12", batch = "A",
      comparator = "PC12_GSE289929_JGW03", doses = c(0, 4, 20, 100),
      progeny = "progeny_human_top500.tsv", progeny_tax = 9606L,
      dorothea = "dorothea_rat_AB_one_to_one.tsv", dorothea_tax = 10116L),
    GSE227494 = list(organism = "Mus musculus", cell = "C2C12", batch = "Batch1",
      comparator = "C2C12_GSE227494_JGW03", doses = c(0, 20, 100, 500),
      progeny = "progeny_mouse_top500.tsv", progeny_tax = 10090L,
      dorothea = "dorothea_mouse_AB.tsv", dorothea_tax = 10090L),
    GSE299063 = list(organism = "Homo sapiens", cell = "HepG2", batch = "A",
      comparator = "HepG2_GSE299063_JGW03_batchA", doses = c(0, 20, 100, 500),
      progeny = "progeny_human_top500.tsv", progeny_tax = 9606L,
      dorothea = "dorothea_human_AB.tsv", dorothea_tax = 9606L))
}

jgw03_v2_validate_metadata <- function(samples, accession) {
  spec <- jgw03_v2_spec(accession)
  required <- c("dataset_id", "sample_id", "organism", "tissue_or_cell_line",
    "batch", "comparator_set", "exposure_duration_hour", "include_primary",
    "condition", "intervention_code", "formula_ratio", "solvent", "extraction_method",
    "dose_value", "dose_unit", "dose_level", "dose_group", "ordered_dose",
    "dose_numeric", "log10_dose_plus_1")
  jgw03_v2_need(is.data.frame(samples) && nrow(samples) == 12L &&
    all(required %in% names(samples)), "Complete 12-row original locked metadata required")
  exact <- c(dataset_id = accession, organism = spec$organism,
    tissue_or_cell_line = spec$cell, batch = spec$batch, comparator_set = spec$comparator)
  for (field in names(exact)) jgw03_v2_need(!anyNA(samples[[field]]) &&
    all(as.character(samples[[field]]) == exact[[field]]), paste("Original metadata differs:", field))
  ids <- as.character(samples$sample_id)
  jgw03_v2_need(!anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids), "Unique ordered sample IDs required")
  jgw03_v2_need(!anyNA(samples$include_primary) &&
    all(as.character(samples$include_primary) == "TRUE") &&
    is.numeric(samples$exposure_duration_hour) &&
    all(is.finite(samples$exposure_duration_hour)) && all(samples$exposure_duration_hour == 24),
    "Frozen inclusion and 24-hour exposure metadata required")
  groups <- c("vehicle", "low", "medium", "high")
  group <- as.character(samples$dose_group)
  ordinal <- match(group, groups) - 1L
  jgw03_v2_need(!anyNA(ordinal) && identical(tabulate(ordinal + 1L, 4L), rep(3L, 4L)),
    "Original four groups of three required")
  for (field in c("ordered_dose", "dose_numeric", "log10_dose_plus_1")) {
    jgw03_v2_need(is.numeric(samples[[field]]) && all(is.finite(samples[[field]])),
      paste("Finite original numeric dose field required:", field))
  }
  expected_dose <- spec$doses[ordinal + 1L]
  # Locked TSVs store 15 significant digits; preserve their values and accept
  # only the documented serialization precision, without overwriting originals.
  expected_log <- log10(expected_dose + 1)
  jgw03_v2_need(all(samples$ordered_dose == ordinal) &&
    all(samples$dose_numeric == expected_dose) &&
    all(abs(samples$log10_dose_plus_1 - expected_log) <= 1e-12 * (1 + abs(expected_log))),
    "Original dose group/ordinal/numeric/log-dose fields disagree")
  vehicle <- ordinal == 0L
  treatment <- !vehicle
  jgw03_v2_need(!anyNA(samples$condition) &&
    all(samples$condition[vehicle] == "vehicle") && all(samples$condition[treatment] == "treatment") &&
    !anyNA(samples$dose_level) && all(samples$dose_level[vehicle] == "control") &&
    all(samples$dose_level[treatment] == group[treatment]), "Original treatment/dose labels disagree")
  for (field in c("intervention_code", "formula_ratio", "solvent", "extraction_method", "dose_unit")) {
    value <- c(intervention_code = "JGW03", formula_ratio = "1:1", solvent = "water",
      extraction_method = "CEM", dose_unit = "ug/mL")[[field]]
    jgw03_v2_need(!anyNA(samples[[field]][treatment]) &&
      all(as.character(samples[[field]][treatment]) == value), paste("Frozen treatment metadata differs:", field))
  }
  jgw03_v2_need(is.numeric(samples$dose_value) &&
    all(is.finite(samples$dose_value[treatment])) &&
    all(samples$dose_value[treatment] == expected_dose[treatment]), "Original treatment dose differs")
  invisible(TRUE)
}

jgw03_v2_expected_files <- function(accession) {
  spec <- jgw03_v2_spec(accession)
  # Analysis outputs are pinned from the retained JGW03 output manifest;
  # resources use the separate gene-set manifest. No historical pin is replaced.
  own <- switch(accession,
    GSE289929 = c(6660, 1933217, 515087,
      "4a1e97a06c692c5750be6895deb66673dcce22d464def2fd7c9f42c1c7212e21",
      "40b640d7c8a8faf3c4bac7aadfa86a0a7101760c564adfa371453c4e3005b531",
      "63764607cc164ef85e845f5d1f263c0165e101e7e144e67b5c117d3cc5181efa"),
    GSE227494 = c(6009, 1486738, 613825,
      "e90434b1b1dcc38e61c958d9c8c18e5708847ef1d936db58750eee01b27bc583",
      "adbeed3211145e444e75ae10d3efd3ba99969912749feb5b4b4e392c8c76d075",
      "a3e19f4d59da5c8caee8da1d82695b184642af4f616da33b327346efa3b91eea"),
    GSE299063 = c(6531, 3187371, 725411,
      "18e414456468b3971e7ce0f9eaf2175bad0813d62892c0bc20236a4d07239bd6",
      "5d912686d05bb0e61f5ac3e3b299217c5dde099f8bc01011c8cc89b24d1f541b",
      "e74c1bc149fd56d525f534e3ff3124e65454d4618927d0e25935229a5af2cbb9"))
  resources <- data.frame(
    filename = c("progeny_human_top500.tsv", "progeny_mouse_top500.tsv",
      "dorothea_human_AB.tsv", "dorothea_mouse_AB.tsv", "dorothea_rat_AB_one_to_one.tsv",
      "reactome97_human_long.tsv", "quickgo_syngo_mammalian_direct_human.gmt",
      "ncbi_human_mouse_rat_one_to_one.tsv"),
    bytes = c(609987, 615522, 367202, 369008, 687990, 11599354, 47129, 1922349),
    sha256 = c("b78892ec4d56212b08eb33ea80551fcb018e89b3c698e3d804008537c96bbdb9",
      "97ac161e3a06b54da68fdb2479d444bc73a685f3794366975aa0c9760abaaf4f",
      "03914edc22ab31a653e4f407e2fcba6a0746939ad5ba9ee6348ad3ec9a9113ab",
      "5b594286c267f9410f7e9e1228fc07575c735b7152a5b8b157c7d4a8a37e2205",
      "e709ea4a9e0a1517cb00082bfbf6e8c4408d786694b6864f31ec3f57b458ecb8",
      "1ea55f41f5b0a4ea51f723fb82b940f8b13fce1cab13ea6825b91d207c72fdc0",
      "883dca3df6a7e5fc86bad1daa770aef167fa07f754afe5a7962e2b6a276f6448",
      "e9dd0c9268e0044f615085e27292cbf977f5a648003ba90054cf3d01225cf998"),
    stringsAsFactors = FALSE)
  chosen <- resources[match(c(spec$progeny, spec$dorothea, "reactome97_human_long.tsv",
    "quickgo_syngo_mammalian_direct_human.gmt", "ncbi_human_mouse_rat_one_to_one.tsv"), resources$filename), ]
  data.frame(role = c("samples", "annotation", "counts", "progeny", "dorothea", "reactome", "syngo", "orthology"),
    relative_path = c(file.path("data/processed/jgw03_three_cell_lines", accession,
      c("locked_samples.tsv", "gene_annotation.tsv", "locked_expected_counts.tsv.gz")),
      file.path("data/processed/gene_sets", chosen$filename)),
    bytes = c(as.numeric(own[1:3]), chosen$bytes), sha256 = c(own[4:6], chosen$sha256),
    stringsAsFactors = FALSE)
}

jgw03_v2_check_file <- function(path, bytes, sha256) {
  jgw03_v2_need(file.exists(path) && !dir.exists(path) && file.info(path)$size == bytes &&
    identical(digest::digest(file = path, algo = "sha256", serialize = FALSE), sha256),
    paste("Pinned input differs:", basename(path)))
  invisible(TRUE)
}

jgw03_v2_table_equal <- function(left, right, serialized_log_field = character()) {
  if (!is.data.frame(left) || !is.data.frame(right) ||
      !identical(names(left), names(right)) || nrow(left) != nrow(right)) return(FALSE)
  all(vapply(names(left), function(field) {
    a <- left[[field]]
    b <- right[[field]]
    if (!identical(is.na(a), is.na(b))) return(FALSE)
    keep <- !is.na(a)
    if (field %in% serialized_log_field && is.numeric(a) && is.numeric(b)) {
      return(all(abs(a[keep] - b[keep]) <= 1e-12 * (1 + abs(b[keep]))))
    }
    if (is.numeric(a) && is.numeric(b)) all(a[keep] == b[keep]) else
      identical(as.character(a[keep]), as.character(b[keep]))
  }, logical(1L)))
}

jgw03_v2_validate_resources <- function(resources, accession) {
  spec <- jgw03_v2_spec(accession)
  jgw03_v2_need(is.list(resources) && identical(names(resources), c("collections", "progeny", "dorothea")) &&
    identical(names(resources$collections), c("Reactome97", "QuickGO_SynGO_sensitivity")),
    "Exactly the original collections and activity resources are required")
  for (role in c("progeny", "dorothea")) {
    tab <- resources[[role]]
    required <- if (role == "progeny") c("pathway", "gene_symbol", "weight", "tax_id", "release") else
      c("tf", "target", "mor", "tax_id", "release", "confidence", "mapping")
    taxon <- if (role == "progeny") spec$progeny_tax else spec$dorothea_tax
    release <- if (role == "progeny") "1.34.0" else "1.23.0"
    jgw03_v2_need(is.data.frame(tab) && nrow(tab) > 0L && all(required %in% names(tab)) &&
      !anyNA(tab$tax_id) && all(tab$tax_id == taxon) &&
      !anyNA(tab$release) && all(as.character(tab$release) == release), paste("Resource species/release differs:", role))
    if (role == "dorothea") jgw03_v2_need(!anyNA(tab$confidence) && all(tab$confidence %in% c("A", "B")),
      "DoRothEA must retain the frozen A/B confidence set")
    if (role == "dorothea" && accession == "GSE289929") jgw03_v2_need(!anyNA(tab$mapping) &&
      all(tab$mapping == "NCBI_Datasets_REST_gene_group_strict_one_to_one_from_human"),
      "PC12 needs the frozen rat one-to-one DoRothEA mapping")
  }
  for (collection in resources$collections) {
    jgw03_v2_need(is.list(collection) && identical(names(collection), c("sets", "term_meta")) &&
      is.list(collection$sets) && length(collection$sets) > 0L &&
      !is.null(names(collection$sets)) && !anyNA(names(collection$sets)) &&
      !anyDuplicated(names(collection$sets)) && is.data.frame(collection$term_meta) &&
      all(c("term_id", "term_name") %in% names(collection$term_meta)) &&
      !anyNA(collection$term_meta$term_id) && !anyDuplicated(collection$term_meta$term_id) &&
      setequal(names(collection$sets), as.character(collection$term_meta$term_id)),
      "Collection keys and original term metadata differ")
  }
  invisible(TRUE)
}

jgw03_v2_read_bound_inputs <- function(project_root, accession) {
  # Defined now, but not invoked on real inputs in the preparation tests.
  records <- jgw03_v2_expected_files(accession)
  records$path <- file.path(normalizePath(project_root, winslash = "/", mustWork = TRUE), records$relative_path)
  for (i in seq_len(nrow(records))) jgw03_v2_check_file(records$path[i], records$bytes[i], records$sha256[i])
  paths <- setNames(records$path, records$role)
  read_table <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  samples <- read_table(paths[["samples"]])
  jgw03_v2_validate_metadata(samples, accession)
  annotation <- read_table(paths[["annotation"]])
  count_table <- read_table(paths[["counts"]])
  jgw03_v2_need(identical(names(count_table), c("raw_gene_id", as.character(samples$sample_id))) &&
    identical(as.character(count_table$raw_gene_id), as.character(annotation$raw_gene_id)) &&
    !anyNA(annotation$raw_gene_id) && !anyDuplicated(annotation$raw_gene_id), "Full count/annotation/sample order differs")
  counts <- as.matrix(count_table[, -1L, drop = FALSE])
  rownames(counts) <- count_table$raw_gene_id
  jgw03_v2_need(is.numeric(counts) && all(is.finite(counts)) && all(counts >= 0), "Finite nonnegative fractional counts required")
  reactome <- read_table(paths[["reactome"]])
  parts <- strsplit(readLines(paths[["syngo"]], warn = FALSE), "\t", fixed = TRUE)
  syngo_sets <- lapply(parts, function(z) unique(toupper(z[-c(1L, 2L)])))
  names(syngo_sets) <- vapply(parts, `[[`, character(1L), 1L)
  # Preserve the native first-column SynGO names, not the second-column GO IDs.
  resources <- list(collections = list(
    Reactome97 = list(sets = lapply(split(toupper(reactome$gene_symbol), reactome$term_id), unique),
      term_meta = unique(reactome[, c("term_id", "term_name")])),
    QuickGO_SynGO_sensitivity = list(sets = syngo_sets,
      term_meta = data.frame(term_id = names(syngo_sets), term_name = names(syngo_sets)))),
    progeny = read_table(paths[["progeny"]]), dorothea = read_table(paths[["dorothea"]]))
  jgw03_v2_validate_resources(resources, accession)
  for (i in seq_len(nrow(records))) jgw03_v2_check_file(records$path[i], records$bytes[i], records$sha256[i])
  list(counts = counts, annotation = annotation, samples = samples, resources = resources,
    source_records = records, accession = accession)
}

jgw03_v2_fingerprint <- function(fixed) {
  # Include the source-record identity too; exclude only this digest's own value.
  if (is.list(fixed$input_admission_v2)) fixed$input_admission_v2$fixed_ram_sha256 <- NULL
  digest::digest(fixed, algo = "sha256", serializeVersion = 3L)
}

jgw03_v2_bind_kernel <- function(core_path) {
  jgw03_v2_need(identical(digest::digest(file = core_path, algo = "sha256", serialize = FALSE),
    "e3939f371f02b24cfdfefee0a9cd9e39bdf47c8920caa44f58434a93f7a2974e"), "Preserved v1 kernel source differs")
  target <- new.env(parent = asNamespace("stats"))
  sys.source(core_path, envir = target)
  old_relabel <- target$jgw03_relabel_model_metadata
  old_build <- target$jgw03_build_fixed_layer
  old_fit <- target$jgw03_fit_one_allocation
  target$jgw03_relabel_model_metadata <- function(samples, allocation, accession) {
    jgw03_v2_validate_metadata(samples, accession)
    old_relabel(samples, allocation, accession)
  }
  target$jgw03_build_fixed_layer <- function(counts, annotation, samples, accession,
                                           resources, native_source_path, project_root) {
    jgw03_v2_validate_metadata(samples, accession)
    jgw03_v2_validate_resources(resources, accession)
    bound <- jgw03_v2_read_bound_inputs(project_root, accession)
    jgw03_v2_need(identical(dimnames(counts), dimnames(bound$counts)) &&
      is.numeric(counts) && all(is.finite(counts)) && all(counts == bound$counts) &&
      jgw03_v2_table_equal(annotation, bound$annotation) &&
      jgw03_v2_table_equal(samples, bound$samples, "log10_dose_plus_1"),
      "RAM inputs differ from exact frozen source files")
    jgw03_v2_need(identical(resources, bound$resources), "RAM resources differ from exact frozen source files")
    fixed <- old_build(counts, annotation, samples, accession, resources, native_source_path)
    fixed$input_admission_v2 <- list(schema_version = "jgw03-input-admission-v2",
      accession = accession, source_records = bound$source_records,
      input_identity_only_not_real_execution_authority = TRUE)
    fixed$input_admission_v2$fixed_ram_sha256 <- jgw03_v2_fingerprint(fixed)
    fixed
  }
  target$jgw03_fit_one_allocation <- function(fixed, allocation, native_source_path) {
    identity <- fixed$input_admission_v2
    jgw03_v2_need(is.list(identity) && identical(identity$schema_version, "jgw03-input-admission-v2") &&
      identical(identity$accession, fixed$accession) &&
      identical(identity$fixed_ram_sha256, jgw03_v2_fingerprint(fixed)),
      "A complete unmodified v2-bound fixed layer is required")
    result <- old_fit(fixed, allocation, native_source_path)
    jgw03_v2_need(is.list(result), "Original kernel must return its native result list")
    result$input_admission_v2 <- identity
    result
  }
  target
}
