#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(file_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
})

raw_root <- file.path(project_root, "data", "raw", "gene_sets")
processed_root <- file.path(project_root, "data", "processed", "gene_sets")
metadata_root <- file.path(project_root, "metadata")
diagnostic_root <- file.path(project_root, "results", "diagnostics")
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_root, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "")
}

write_gmt <- function(x, path, collection_column = "term_name", id_column = "term_id", gene_column = "gene_symbol") {
  stopifnot(all(c(collection_column, id_column, gene_column) %in% names(x)))
  split_key <- paste(x[[collection_column]], x[[id_column]], sep = "\t")
  genes <- split(x[[gene_column]], split_key)
  lines <- vapply(names(genes), function(key) {
    key_parts <- strsplit(key, "\t", fixed = TRUE)[[1]]
    members <- sort(unique(genes[[key]]))
    paste(c(key_parts[1], key_parts[2], members), collapse = "\t")
  }, character(1))
  writeLines(lines, path, useBytes = TRUE)
}

map_entrez_to_symbol <- function(ids, db) {
  values <- AnnotationDbi::mapIds(
    db,
    keys = as.character(ids),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(as.character(values))
}

map_symbol_to_entrez <- function(symbols, db) {
  values <- AnnotationDbi::mapIds(
    db,
    keys = as.character(symbols),
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  unname(as.character(values))
}

## -------------------------------------------------------------------------
## NCBI Datasets REST gene-group reports (API 18.35.0). Protein-coding gene
## metadata were retrieved separately for each species. We join on the NCBI
## Gene Group Identifier and retain groups containing exactly one gene in each
## species being compared. This is the sole gene-level cross-species map.
## -------------------------------------------------------------------------
orthology_root <- file.path(raw_root, "orthology", "ncbi_datasets_rest")
orthology_summary <- fread(file.path(orthology_root, "download_summary.tsv"))
stopifnot(
  nrow(orthology_summary) == 3L,
  uniqueN(orthology_summary$api_version) == 1L,
  unique(orthology_summary$api_version) == "18.35.0",
  identical(
    as.integer(orthology_summary$official_report_count),
    c(20598L, 26313L, 23246L)
  ),
  identical(
    as.integer(orthology_summary$empty_gene_reports),
    c(2L, 0L, 1L)
  )
)

read_gene_group_report <- function(label, expected_tax_id, expected_rows) {
  path <- file.path(
    orthology_root,
    sprintf("NCBI_%s_protein_coding_gene_groups_2026-08-18.tsv", label)
  )
  table <- fread(path, sep = "\t", quote = "", na.strings = c("", "NA"))
  stopifnot(identical(
    names(table),
    c(
      "NCBI GeneID", "Symbol", "Gene Group Identifier", "Gene Group Method",
      "Taxonomic ID", "Taxonomic Name", "Gene Type"
    )
  ))
  setnames(
    table,
    c(
      "NCBI GeneID", "Symbol", "Gene Group Identifier", "Gene Group Method",
      "Taxonomic ID", "Taxonomic Name", "Gene Type"
    ),
    c("entrez", "symbol", "group_id", "group_method", "tax_id", "taxon", "gene_type")
  )
  table[, `:=`(
    entrez = as.character(entrez),
    symbol = as.character(symbol),
    group_id = as.character(group_id),
    tax_id = as.integer(tax_id)
  )]
  stopifnot(
    nrow(table) == expected_rows,
    all(table$tax_id == expected_tax_id),
    all(table$gene_type == "PROTEIN_CODING"),
    !anyDuplicated(table$entrez)
  )
  table
}

orth_human_source <- read_gene_group_report("human", 9606L, 20596L)
orth_mouse_source <- read_gene_group_report("mouse", 10090L, 26313L)
orth_rat_source <- read_gene_group_report("rat", 10116L, 23245L)

make_one_to_one <- function(other_table, other_tax_id) {
  human <- orth_human_source[
    group_method == "NCBI Ortholog" & !is.na(group_id) & nzchar(group_id)
  ]
  other <- other_table[
    group_method == "NCBI Ortholog" & !is.na(group_id) & nzchar(group_id)
  ]
  human_groups <- human[, .(human_n = .N), by = group_id][human_n == 1L, group_id]
  other_groups <- other[, .(other_n = .N), by = group_id][other_n == 1L, group_id]
  eligible_groups <- intersect(human_groups, other_groups)
  normalized <- merge(
    human[group_id %chin% eligible_groups, .(
      group_id,
      group_method,
      human_entrez = entrez,
      human_symbol = symbol
    )],
    other[group_id %chin% eligible_groups, .(
      group_id,
      other_entrez = entrez,
      other_symbol = symbol
    )],
    by = "group_id",
    all = FALSE
  )
  normalized[, `:=`(human_tax_id = 9606L, other_tax_id = as.integer(other_tax_id))]
  setcolorder(normalized, c(
    "group_id", "group_method",
    "human_tax_id", "human_entrez", "human_symbol",
    "other_tax_id", "other_entrez", "other_symbol"
  ))
  setorder(normalized, group_id)
  stopifnot(
    !anyDuplicated(normalized$human_entrez),
    !anyDuplicated(normalized$other_entrez),
    !anyDuplicated(normalized$group_id)
  )
  normalized
}

orth_mouse <- make_one_to_one(orth_mouse_source, 10090L)
orth_rat <- make_one_to_one(orth_rat_source, 10116L)
orth_one_to_one <- rbindlist(list(orth_mouse, orth_rat), use.names = TRUE)
write_tsv(orth_one_to_one, file.path(processed_root, "ncbi_human_mouse_rat_one_to_one.tsv"))

## -------------------------------------------------------------------------
## Reactome release 97 (official current release observed 2026-08-17).
## Viral/non-human symbols are removed by intersection with HGNC symbols.
## -------------------------------------------------------------------------
reactome_file <- file.path(raw_root, "reactome", "release97", "ReactomePathways.gmt")
reactome_lines <- readLines(reactome_file, warn = FALSE)
reactome_parts <- strsplit(reactome_lines, "\t", fixed = TRUE)
reactome <- rbindlist(lapply(reactome_parts, function(parts) {
  if (length(parts) < 3L) return(NULL)
  data.table(term_name = parts[1], term_id = parts[2], gene_symbol = parts[-c(1, 2)])
}))
valid_human_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
reactome <- unique(reactome[grepl("^R-HSA-", term_id) & gene_symbol %chin% valid_human_symbols])
reactome[, `:=`(
  collection = "Reactome",
  release = "97",
  organism = "Homo sapiens",
  tax_id = 9606L
)]
setcolorder(reactome, c("collection", "release", "term_id", "term_name", "tax_id", "organism", "gene_symbol"))
setorder(reactome, term_id, gene_symbol)
write_tsv(reactome, file.path(processed_root, "reactome97_human_long.tsv"))
write_gmt(reactome, file.path(processed_root, "reactome97_human.gmt"))

## -------------------------------------------------------------------------
## PROGENy 1.34.0. Reconstruct getModel(top=500) deterministically from the
## packaged full model so no live OmniPath/Ensembl lookup can change results.
## -------------------------------------------------------------------------
suppressPackageStartupMessages(library(progeny))
data("model_human_full", package = "progeny")
data("model_mouse_full", package = "progeny")

top_model <- function(model, organism, tax_id) {
  table <- as.data.table(model)
  table[, rank_within_pathway := frank(p.value, ties.method = "first"), by = pathway]
  table <- table[rank_within_pathway <= 500L]
  table[, `:=`(
    collection = "PROGENy",
    release = as.character(packageVersion("progeny")),
    organism = organism,
    tax_id = as.integer(tax_id)
  )]
  setnames(table, "gene", "gene_symbol")
  setcolorder(table, c(
    "collection", "release", "pathway", "tax_id", "organism",
    "gene_symbol", "weight", "p.value", "rank_within_pathway"
  ))
  setorder(table, pathway, rank_within_pathway, gene_symbol)
  table
}

progeny_human <- top_model(model_human_full, "Homo sapiens", 9606L)
progeny_mouse <- top_model(model_mouse_full, "Mus musculus", 10090L)
write_tsv(progeny_human, file.path(processed_root, "progeny_human_top500.tsv"))
write_tsv(progeny_mouse, file.path(processed_root, "progeny_mouse_top500.tsv"))

## -------------------------------------------------------------------------
## DoRothEA 1.23.0, confidence A/B only. Rat is a fully transparent derived
## resource requiring one-to-one mapping of both TF and target from human.
## -------------------------------------------------------------------------
suppressPackageStartupMessages(library(dorothea))
data("dorothea_hs", package = "dorothea")
data("dorothea_mm", package = "dorothea")

format_dorothea <- function(x, organism, tax_id, mapping = "native") {
  table <- as.data.table(x)[confidence %chin% c("A", "B")]
  table[, `:=`(
    collection = "DoRothEA",
    release = as.character(packageVersion("dorothea")),
    organism = organism,
    tax_id = as.integer(tax_id),
    mapping = mapping
  )]
  setcolorder(table, c(
    "collection", "release", "confidence", "tax_id", "organism",
    "tf", "target", "mor", "mapping"
  ))
  unique(table)
}

dorothea_human <- format_dorothea(dorothea_hs, "Homo sapiens", 9606L)
dorothea_mouse <- format_dorothea(dorothea_mm, "Mus musculus", 10090L)
rat_tf_map <- orth_rat[, .(tf = human_symbol, rat_tf = other_symbol)]
rat_target_map <- orth_rat[, .(target = human_symbol, rat_target = other_symbol)]
dorothea_rat <- merge(
  merge(dorothea_human, rat_tf_map, by = "tf", all = FALSE),
  rat_target_map,
  by = "target",
  all = FALSE
)
dorothea_rat[, `:=`(
  tf = rat_tf,
  target = rat_target,
  organism = "Rattus norvegicus",
  tax_id = 10116L,
  mapping = "NCBI_Datasets_REST_gene_group_strict_one_to_one_from_human"
)]
dorothea_rat[, c("rat_tf", "rat_target") := NULL]
setcolorder(dorothea_rat, names(dorothea_human))
dorothea_rat <- unique(dorothea_rat)

write_tsv(dorothea_human, file.path(processed_root, "dorothea_human_AB.tsv"))
write_tsv(dorothea_mouse, file.path(processed_root, "dorothea_mouse_AB.tsv"))
write_tsv(dorothea_rat, file.path(processed_root, "dorothea_rat_AB_one_to_one.tsv"))

## -------------------------------------------------------------------------
## QuickGO snapshot of annotations explicitly assigned by SynGO. This is a
## conservative, direct-annotation, mammalian-only sensitivity resource and is
## not represented as the SynGO portal release 1.3 bulk dataset.
## -------------------------------------------------------------------------
syngo_file <- file.path(raw_root, "syngo", "QuickGO_SynGO_assignedBy_2026-08-17.tsv")
syngo <- fread(syngo_file, sep = "\t", quote = "", na.strings = c("", "NA"))
setnames(syngo, c(
  "GO NAME", "SYMBOL", "GO TERM", "GO ASPECT", "ECO ID",
  "GO EVIDENCE CODE", "REFERENCE", "TAXON ID", "TAXON NAME", "DATE"
), c(
  "term_name", "source_symbol", "term_id", "aspect", "eco_id",
  "go_evidence", "reference", "source_tax_id", "source_taxon", "annotation_date"
))
syngo <- syngo[source_tax_id %in% c(9606L, 10090L, 10116L)]
human_syngo <- syngo[source_tax_id == 9606L]
human_syngo[, gene_symbol := source_symbol]
mouse_syngo <- merge(
  syngo[source_tax_id == 10090L],
  orth_mouse[, .(source_symbol = other_symbol, gene_symbol = human_symbol)],
  by = "source_symbol",
  all = FALSE
)
rat_syngo <- merge(
  syngo[source_tax_id == 10116L],
  orth_rat[, .(source_symbol = other_symbol, gene_symbol = human_symbol)],
  by = "source_symbol",
  all = FALSE
)
syngo_mammalian <- rbindlist(list(human_syngo, mouse_syngo, rat_syngo), fill = TRUE, use.names = TRUE)
syngo_mammalian <- syngo_mammalian[gene_symbol %chin% valid_human_symbols]
syngo_direct <- unique(syngo_mammalian[, .(
  collection = "QuickGO_assignedBy_SynGO",
  release = "QuickGO_annotation_snapshot_2026-08-17",
  term_id,
  term_name,
  aspect,
  gene_symbol,
  source_tax_id,
  source_taxon,
  source_symbol,
  eco_id,
  go_evidence,
  reference,
  annotation_date,
  mapping = fifelse(source_tax_id == 9606L, "native_human", "NCBI_Datasets_REST_gene_group_strict_one_to_one")
)])
setorder(syngo_direct, term_id, gene_symbol, source_tax_id, reference)
write_tsv(syngo_direct, file.path(processed_root, "quickgo_syngo_mammalian_direct_human.tsv"))
syngo_gene_sets <- unique(syngo_direct[, .(term_id, term_name, gene_symbol)])
write_gmt(syngo_gene_sets, file.path(processed_root, "quickgo_syngo_mammalian_direct_human.gmt"))

## -------------------------------------------------------------------------
## Matrisome Portal official master lists (rev. 2014). These source workbooks
## expose a dedicated Gene Symbol column. This avoids treating the mixed
## synonym/UniProt/RefSeq identifier column in MatrisomeAnalyzeR as if every
## matching string were a canonical symbol. Exact current symbols are retained;
## obsolete symbols are updated only when an AnnotationDbi alias is unambiguous.
## -------------------------------------------------------------------------
matrisome_root <- file.path(raw_root, "matrisome", "official_masterlists")
mat_human_source <- fread(
  file.path(matrisome_root, "Hs_Matrisome_Masterlist_rev2014.csv"),
  skip = 1L,
  header = TRUE,
  check.names = FALSE,
  na.strings = c("", "NA", "None")
)
mat_mouse_source <- fread(
  file.path(matrisome_root, "Mm_Matrisome_Masterlist_rev2014.csv"),
  skip = 1L,
  header = TRUE,
  check.names = FALSE,
  na.strings = c("", "NA", "None")
)

resolve_current_symbols <- function(source_symbols, db) {
  result <- data.table(source_gene_symbol = as.character(source_symbols))
  valid_symbols <- keys(db, keytype = "SYMBOL")
  result[, `:=`(
    gene_symbol = fifelse(source_gene_symbol %chin% valid_symbols, source_gene_symbol, NA_character_),
    symbol_mapping = fifelse(source_gene_symbol %chin% valid_symbols, "exact_current_symbol", "unresolved")
  )]
  alias_keys <- keys(db, keytype = "ALIAS")
  candidates <- unique(result[is.na(gene_symbol) & source_gene_symbol %chin% alias_keys, source_gene_symbol])
  if (length(candidates) > 0L) {
    alias_hits <- as.data.table(AnnotationDbi::select(
      db,
      keys = candidates,
      keytype = "ALIAS",
      columns = "SYMBOL"
    ))
    alias_hits <- alias_hits[!is.na(SYMBOL), .(
      mapped_symbol = if (uniqueN(SYMBOL) == 1L) unique(SYMBOL) else NA_character_,
      mapping_n = uniqueN(SYMBOL)
    ), by = ALIAS]
    result[alias_hits, on = .(source_gene_symbol = ALIAS), `:=`(
      gene_symbol = i.mapped_symbol,
      symbol_mapping = fifelse(i.mapping_n == 1L, "unique_current_alias", "ambiguous_current_alias")
    )]
  }
  result[is.na(gene_symbol), gene_symbol := source_gene_symbol]
  result
}

format_matrisome <- function(source, organism, tax_id, db, expected_rows) {
  stopifnot(
    all(c("Matrisome Division", "Matrisome Category", "Gene Symbol") %in% names(source)),
    nrow(source) == expected_rows,
    uniqueN(source[["Gene Symbol"]]) == expected_rows
  )
  symbol_map <- resolve_current_symbols(source[["Gene Symbol"]], db)
  table <- cbind(
    source[, .(
      division = `Matrisome Division`,
      category = `Matrisome Category`,
      source_gene_symbol = `Gene Symbol`
    )],
    symbol_map[, .(gene_symbol, symbol_mapping)]
  )
  table[, `:=`(
    collection = "Matrisome_Portal_masterlist",
    release = "rev_2014_portal_snapshot_2026-08-18",
    tax_id = as.integer(tax_id),
    organism = organism
  )]
  setcolorder(table, c(
    "collection", "release", "tax_id", "organism", "division", "category",
    "source_gene_symbol", "gene_symbol", "symbol_mapping"
  ))
  unique(table)
}

mat_human <- format_matrisome(mat_human_source, "Homo sapiens", 9606L, org.Hs.eg.db, 1027L)
mat_mouse <- format_matrisome(mat_mouse_source, "Mus musculus", 10090L, org.Mm.eg.db, 1110L)
setorder(mat_human, division, category, gene_symbol)
setorder(mat_mouse, division, category, gene_symbol)
write_tsv(mat_human, file.path(processed_root, "matrisome_human.tsv"))
write_tsv(mat_mouse, file.path(processed_root, "matrisome_mouse.tsv"))

## Analysis-ready R object and diagnostic inventory.
resource_bundle <- list(
  reactome = reactome,
  progeny_human = progeny_human,
  progeny_mouse = progeny_mouse,
  dorothea_human = dorothea_human,
  dorothea_mouse = dorothea_mouse,
  dorothea_rat = dorothea_rat,
  syngo_sensitivity = syngo_gene_sets,
  matrisome_human = mat_human,
  matrisome_mouse = mat_mouse,
  ortholog_one_to_one = orth_one_to_one
)
saveRDS(resource_bundle, file.path(processed_root, "gene_set_resources_frozen_2026-08-18.rds"), version = 3)

inventory <- data.table(
  resource = c(
    "Reactome_human", "PROGENy_human", "PROGENy_mouse",
    "DoRothEA_human_AB", "DoRothEA_mouse_AB", "DoRothEA_rat_AB_derived",
    "QuickGO_SynGO_mammalian_direct_human_sensitivity",
    "Matrisome_human", "Matrisome_mouse",
    "NCBI_human_mouse_one_to_one", "NCBI_human_rat_one_to_one"
  ),
  release = c(
    "97", as.character(packageVersion("progeny")), as.character(packageVersion("progeny")),
    as.character(packageVersion("dorothea")), as.character(packageVersion("dorothea")), as.character(packageVersion("dorothea")),
    "QuickGO_snapshot_2026-08-17",
    "rev_2014_portal_snapshot_2026-08-18",
    "rev_2014_portal_snapshot_2026-08-18",
    "NCBI_Datasets_REST_v18.35.0_snapshot_2026-08-18",
    "NCBI_Datasets_REST_v18.35.0_snapshot_2026-08-18"
  ),
  records = c(
    nrow(reactome), nrow(progeny_human), nrow(progeny_mouse),
    nrow(dorothea_human), nrow(dorothea_mouse), nrow(dorothea_rat),
    nrow(syngo_gene_sets), nrow(mat_human), nrow(mat_mouse),
    nrow(orth_mouse), nrow(orth_rat)
  ),
  unique_features = c(
    uniqueN(reactome$gene_symbol), uniqueN(progeny_human$gene_symbol), uniqueN(progeny_mouse$gene_symbol),
    uniqueN(dorothea_human$target), uniqueN(dorothea_mouse$target), uniqueN(dorothea_rat$target),
    uniqueN(syngo_gene_sets$gene_symbol), uniqueN(mat_human$gene_symbol), uniqueN(mat_mouse$gene_symbol),
    uniqueN(orth_mouse$human_symbol), uniqueN(orth_rat$human_symbol)
  ),
  unique_programs = c(
    uniqueN(reactome$term_id), uniqueN(progeny_human$pathway), uniqueN(progeny_mouse$pathway),
    uniqueN(dorothea_human$tf), uniqueN(dorothea_mouse$tf), uniqueN(dorothea_rat$tf),
    uniqueN(syngo_gene_sets$term_id), uniqueN(mat_human$category), uniqueN(mat_mouse$category),
    NA_integer_, NA_integer_
  )
)
write_tsv(inventory, file.path(metadata_root, "gene_set_inventory.tsv"))

checks <- data.table(
  check_id = c(
    "reactome_release97_nonempty",
    "progeny_human_14x500",
    "progeny_mouse_14x500",
    "dorothea_human_AB_only",
    "dorothea_mouse_AB_only",
    "ncbi_rest_official_report_counts",
    "ncbi_rest_source_row_counts",
    "ncbi_rest_group_method",
    "ortholog_pairs_unique",
    "quickgo_syngo_export_complete_26133",
    "quickgo_syngo_sensitivity_nonempty",
    "matrisome_human_portal_1027",
    "matrisome_mouse_portal_1110",
    "matrisome_symbol_resolution"
  ),
  passed = c(
    nrow(reactome) > 0L,
    nrow(progeny_human) == 14L * 500L,
    nrow(progeny_mouse) == 14L * 500L,
    all(dorothea_human$confidence %chin% c("A", "B")),
    all(dorothea_mouse$confidence %chin% c("A", "B")),
    identical(
      as.integer(orthology_summary$official_report_count),
      c(20598L, 26313L, 23246L)
    ) && identical(
      as.integer(orthology_summary$empty_gene_reports),
      c(2L, 0L, 1L)
    ),
    identical(
      c(nrow(orth_human_source), nrow(orth_mouse_source), nrow(orth_rat_source)),
      c(20596L, 26313L, 23245L)
    ),
    all(vapply(
      list(orth_human_source, orth_mouse_source, orth_rat_source),
      function(x) all(x[!is.na(group_id), group_method] == "NCBI Ortholog"),
      logical(1)
    )),
    !anyDuplicated(orth_one_to_one[, .(human_symbol, other_tax_id)]) &&
      !anyDuplicated(orth_one_to_one[, .(other_symbol, other_tax_id)]),
    nrow(fread(syngo_file, select = 1L)) == 26133L,
    nrow(syngo_gene_sets) > 0L,
    nrow(mat_human_source) == 1027L && uniqueN(mat_human$source_gene_symbol) == 1027L,
    nrow(mat_mouse_source) == 1110L && uniqueN(mat_mouse$source_gene_symbol) == 1110L,
    mean(mat_human$symbol_mapping %chin% c("exact_current_symbol", "unique_current_alias")) >= 0.98 &&
      mean(mat_mouse$symbol_mapping %chin% c("exact_current_symbol", "unique_current_alias")) >= 0.98
  )
)
checks[, status := fifelse(passed, "PASS", "FAIL")]
write_tsv(checks, file.path(diagnostic_root, "gene_set_freeze_validation.tsv"))

if (any(!checks$passed)) {
  stop("Gene-set freeze validation failed: ", paste(checks[passed == FALSE, check_id], collapse = ", "))
}

cat("Gene-set resources frozen successfully.\n")
print(inventory)
print(checks[, .(check_id, status)])
