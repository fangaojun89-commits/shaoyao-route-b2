# Pure, additive admission of the five model-v2 member resources. This does not
# read RDS/files, run a model, authorize production, or modify the input. The
# returned resources are the original RAM objects, not physically immutable
# objects; callers must not export them as a safe diagnostic projection.
# Existing analyzer/validator statistical functions remain unchanged.

guard <- function(resources) {
  need <- function(value, code) {
    if (!isTRUE(value)) {
      stop(paste0("GSE128623_MODEL_MEMBER_GUARD:", code), call. = FALSE)
    }
  }
  roles <- c("reactome", "syngo_sensitivity", "ortholog_one_to_one",
             "progeny_mouse", "dorothea_mouse")
  need(is.list(resources) && !is.object(resources) && !is.data.frame(resources) &&
         !is.null(names(resources)) && !anyNA(names(resources)) &&
         all(nzchar(names(resources))) && !anyDuplicated(names(resources)) &&
         all(roles %in% names(resources)), "RESOURCE_SCHEMA")

  check_table <- function(value, columns, numeric_columns = character(),
                          factor_label_columns = character()) {
    ordinary_class <- identical(class(value), "data.frame") ||
      identical(class(value), c("data.table", "data.frame"))
    need(ordinary_class && is.data.frame(value) && nrow(value) > 0L &&
           !anyNA(names(value)) && !anyDuplicated(names(value)) &&
           all(columns %in% names(value)), "TABLE_SCHEMA")
    for (column in columns) {
      cells <- value[[column]]
      if (column %in% numeric_columns) {
        need((typeof(cells) %in% c("integer", "double")) &&
               !is.object(cells) && is.null(dim(cells)) &&
               !anyNA(cells) && all(is.finite(cells)), "NUMERIC_COLUMN")
      } else {
        # Only PROGENy labels may be factors. Compare their labels, never codes.
        # This local vector conversion leaves the original resource untouched.
        if (column %in% factor_label_columns && is.factor(cells)) {
          need(identical(class(cells), "factor"), "FACTOR_LABEL_CLASS")
          cells <- as.character(cells)
        }
        need(is.character(cells) && !is.object(cells) &&
               is.null(dim(cells)) && !anyNA(cells) &&
               all(nzchar(trimws(cells))), "CHARACTER_COLUMN")
      }
    }
    invisible(TRUE)
  }
  check_table(resources$reactome,
              c("release", "term_id", "term_name", "gene_symbol"))
  check_table(resources$syngo_sensitivity,
              c("term_id", "term_name", "gene_symbol"))
  check_table(resources$ortholog_one_to_one,
              c("other_tax_id", "human_symbol", "other_symbol"), "other_tax_id")
  check_table(resources$progeny_mouse,
              c("pathway", "gene_symbol", "weight"), "weight",
              factor_label_columns = c("pathway", "gene_symbol"))
  check_table(resources$dorothea_mouse, c("tf", "target", "mor"), "mor")
  need(identical(unique(resources$reactome$release), "97"), "REACTOME_RELEASE")
  taxa <- resources$ortholog_one_to_one$other_tax_id
  need(all(taxa > 0 & taxa == floor(taxa)), "TAXON_ID")

  # Work on ordinary new data frames; never use data.table by-reference updates.
  orthology <- resources$ortholog_one_to_one
  mouse_rows <- orthology$other_tax_id == 10090L
  mapping <- data.frame(
    human_symbol = orthology$human_symbol[mouse_rows],
    gene_symbol = orthology$other_symbol[mouse_rows], stringsAsFactors = FALSE
  )
  need(nrow(mapping) > 0L, "MOUSE_MAPPING_EMPTY")
  check_conflicts <- function(tab, keys, values, code) {
    if (nrow(tab) < 2L) return(invisible(TRUE))
    row_order <- do.call(order, c(unname(tab[keys]), list(method = "radix")))
    sorted <- tab[row_order, , drop = FALSE]
    same_key <- Reduce(`&`, lapply(keys, function(column) {
      sorted[[column]][-1L] == sorted[[column]][-nrow(sorted)]
    }))
    for (column in values) {
      differs <- sorted[[column]][-1L] != sorted[[column]][-nrow(sorted)]
      need(!any(same_key & differs), code)
    }
    invisible(TRUE)
  }
  check_conflicts(mapping, "human_symbol", "gene_symbol", "ORTHOLOGY_MAPPING_CONFLICT")
  check_conflicts(mapping, "gene_symbol", "human_symbol", "ORTHOLOGY_MAPPING_CONFLICT")

  human_members <- function(tab, id_column, name_column, collection) {
    source <- data.frame(
      human_symbol = tab$gene_symbol, program_id = tab[[id_column]],
      program_name = tab[[name_column]], stringsAsFactors = FALSE
    )
    # As in the kernels, symbols with no selected mouse mapping contribute no
    # mouse member. They remain untouched in the returned original resource.
    joined <- merge(source, mapping, by = "human_symbol", sort = FALSE)
    data.frame(
      collection = rep(collection, nrow(joined)), program_id = joined$program_id,
      program_name = joined$program_name, gene_symbol = joined$gene_symbol,
      weight = rep(1, nrow(joined)),
      test_class = rep("unsigned_and_activity", nrow(joined)),
      stringsAsFactors = FALSE
    )
  }
  progeny <- resources$progeny_mouse
  dorothea <- resources$dorothea_mouse
  candidates <- rbind(
    human_members(resources$reactome, "term_id", "term_name", "Reactome"),
    human_members(resources$syngo_sensitivity, "term_id", "term_name",
                  "SynGO_assignedBy_proxy_sensitivity"),
    data.frame(
      collection = rep("PROGENy", nrow(progeny)), program_id = progeny$pathway,
      program_name = progeny$pathway, gene_symbol = progeny$gene_symbol,
      weight = progeny$weight, test_class = rep("activity_only", nrow(progeny)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      collection = rep("DoRothEA_AB", nrow(dorothea)), program_id = dorothea$tf,
      program_name = dorothea$tf, gene_symbol = dorothea$target,
      weight = dorothea$mor, test_class = rep("activity_only", nrow(dorothea)),
      stringsAsFactors = FALSE
    )
  )
  need(nrow(candidates) > 0L, "MEMBERS_EMPTY")
  # Compare typed columns directly, without concatenated-key collisions or
  # selecting a first representative. No candidate row is discarded here.
  member_key <- c("collection", "program_id", "gene_symbol")
  check_conflicts(candidates, member_key, "weight", "WEIGHT_CONFLICT")
  check_conflicts(candidates, member_key, "program_name", "NAME_CONFLICT")
  check_conflicts(candidates, member_key, "test_class", "TEST_CLASS_CONFLICT")
  # A single program name is required independently of overlap between members.
  program_key <- c("collection", "program_id")
  check_conflicts(candidates, program_key, "program_name", "PROGRAM_NAME_CONFLICT")
  check_conflicts(candidates, program_key, "test_class", "PROGRAM_TEST_CLASS_CONFLICT")

  list(
    schema = "GSE128623_MODEL_MEMBER_CONTENT_GUARD/v1",
    status = "MEMBER_CONTENT_ONLY_PASS", resources = resources,
    input_modified = FALSE, model_run = FALSE,
    production_authorized = FALSE, scientific_acceptance = FALSE
  )
}
