# Versioned independent non-gene reference wiring. No core compute delegation.
# Deferred refit preserves the original numerical body/settings/comparator.
# SOURCE tests load/prepare only, with explicit MOCK prefit/mapping inputs.

gse205906_ens101_program_reference_derive_v1 <- function(old_text,config_sha256) {
  need <- function(ok) if (!isTRUE(ok)) stop("ENS101_PROGRAM_REFERENCE_EXACT_SOURCE_ANCHOR",call.=FALSE)
  need(is.character(old_text) && length(old_text)==1L &&
    identical(digest::digest(charToRaw(old_text),"sha256",serialize=FALSE),
      "4e0f3cd999b327e06de1388f8fe92472a1af20bf62e9d71d7cea24ed3b44d493"))
  text <- old_text; amendments <- list()
  replace <- function(old,new,label) {
    hits <- gregexpr(old,text,fixed=TRUE)[[1L]]
    need(length(hits)==1L && hits[[1L]]>0L)
    text <<- paste0(substr(text,1L,hits[[1L]]-1L),new,substr(text,hits[[1L]]+nchar(old),nchar(text)))
    amendments[[label]] <<- list(old_sha256=digest::digest(charToRaw(old),"sha256",serialize=FALSE),
      new_sha256=digest::digest(charToRaw(new),"sha256",serialize=FALSE))
  }
  start <- regexpr("prepare_gse205906_bound_program_validation_v2 <- function",text,fixed=TRUE)[[1L]]
  end <- regexpr('  fail <- function(message)',text,fixed=TRUE)[[1L]]
  replace(substr(text,start,end-1L),paste0(
    "prepare_gse205906_ens101_program_reference_engine_v1 <- function(source_bytes,cfg,count_matrix,\n",
    "    gene_annotation,sample_ledger,gene_resources,orthology_table,\n",
    "    domain=\"PRODUCTION\",program_contract,original_raw_inputs) {\n"),"explicit_source_cfg_raw_signature")
  replace('  require_true(identical(domain, "SYNTHETIC"), "BOUND_PUBLIC_EXECUTION_NOT_ACCEPTED")',
    '  require_true(identical(domain, "PRODUCTION"), "ENS101_PROGRAM_REFERENCE_PRODUCTION_SETTINGS_ONLY")',
    "production_settings_only")
  replace(paste0("  program_contract$check_inputs(list(count_matrix=count_matrix, gene_annotation=gene_annotation,\n",
    "    sample_ledger=sample_ledger, gene_resources=gene_resources, orthology_table=orthology_table), domain)"),
    "  program_contract$check_inputs(original_raw_inputs, domain)","raw_contract_identity_not_projected")
  replace(paste0('    "config/GSE205906_bulk_model_v2_r4.json" =\n',
    '      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"'),
    paste0('    "config/GSE205906_Ensembl101_full_model_v1_20261006.json" =\n      "',config_sha256,'"'),"new_config_raw_pin")
  a <- regexpr("  capture <- function(relative, expected) {",text,fixed=TRUE)[[1L]]
  b <- regexpr("  sources <- Map(capture, names(pins), unname(pins))",text,fixed=TRUE)[[1L]]
  replace(substr(text,a,b-1L),paste0('  capture <- function(relative, expected) {\n',
    '    roles <- c("code/R/gse205906_bulk_model_v2_r4_lib.R"="helpers",\n',
    '      "config/GSE205906_Ensembl101_full_model_v1_20261006.json"="full_config")\n',
    '    require_true(relative %in% names(roles),"ENS101_PROGRAM_REFERENCE_SOURCE_ROLE")\n',
    '    payload <- source_bytes[[roles[[relative]]]]\n',
    '    require_true(is.raw(payload) && identical(digest::digest(payload,"sha256",serialize=FALSE),expected),\n',
    '      "PROGRAM_REFERENCE_SOURCE_IDENTITY")\n    rawToChar(payload)\n  }\n'),"no_path_reopen")
  replace("  cfg <- program_contract$effective_config(jsonlite::fromJSON(sources[[2L]], simplifyVector = FALSE))",
    "  cfg <- program_contract$effective_config(cfg)","both_current_family_cardinalities")
  replace("  ledger <- helpers$v2r4_validate_ledger(data.table::copy(sample_ledger), cfg)",
    "  ledger <- data.table::copy(sample_ledger) # Original 53-column guard applied before projection.",
    "original_53_field_library_ledger")
  replace("        fits[[collection_name]] <- data.table::data.table()",
    "        fits[[collection_name]] <- bound_programs$empty_activity(collection_name)","typed_empty_activity")
  replace(paste0("      grid$collection <- collection_name\n",
    "      grid$omitted_animal <- ids[[omission]]\n",
    "      grid$omitted_group <- as.character(group[[omission]])\n",
    "      grid$loao_method <- \"ACTIVITY_FULL_PIPELINE_REFIT\""),paste0(
    "      grid$collection <- rep(collection_name,nrow(grid))\n",
    "      grid$omitted_animal <- rep(ids[[omission]],nrow(grid))\n",
    "      grid$omitted_group <- rep(as.character(group[[omission]]),nrow(grid))\n",
    "      grid$loao_method <- rep(\"ACTIVITY_FULL_PIPELINE_REFIT\",nrow(grid))"),"zero_grid_keeps_zero_rows")
  summary_start <- regexpr("  summary <- expected$loao_long[, {",text,fixed=TRUE)[[1L]]
  summary_end <- regexpr('  minimum_omissions <- if (domain == "SYNTHETIC")',text,fixed=TRUE)[[1L]]
  need(summary_start>0L && summary_end>summary_start)
  old_summary <- substr(text,summary_start,summary_end-1L)
  new_summary <- paste0("  if (nrow(expected$loao_long)) {\n",old_summary,
    "  } else summary <- bound_programs$empty_loao_summary(expected$loao_long,ids)\n")
  replace(old_summary,new_summary,"typed_empty_summary_nonempty_formula_exact")
  replace(paste0("  summary$all_steps_refit <- TRUE\n",
    "  summary$loao_method <- \"ACTIVITY_FULL_PIPELINE_REFIT\""),paste0(
    "  summary$all_steps_refit <- rep(TRUE,nrow(summary))\n",
    "  summary$loao_method <- rep(\"ACTIVITY_FULL_PIPELINE_REFIT\",nrow(summary))"),"zero_summary_keeps_zero_rows")
  # Full and omission preprocessing, 22 membership calls, scores/fits/seeds,
  # min7, two contrasts, twelve table names/keys and numerical limits unchanged.
  list(text=text,amendments=amendments,
    derived_source_sha256=digest::digest(charToRaw(text),"sha256",serialize=FALSE))
}

load_gse205906_ens101_program_reference_v1 <- function(config_bytes,legacy_config_bytes,
    source_bundle,prefit,mapping) {
  need <- function(ok,category) if (!isTRUE(ok)) stop(category,call.=FALSE)
  bootstrap_pin <- list(bytes=13197L,sha256="4c6593f0ed6cdd6838699f8d844dbd7aa321baca990b2dac79abe9ee7c1f76d1")
  need(requireNamespace("digest",quietly=TRUE) && is.list(source_bundle) &&
    is.raw(source_bundle$reference_bootstrap) && length(source_bundle$reference_bootstrap)==bootstrap_pin$bytes &&
    identical(digest::digest(source_bundle$reference_bootstrap,"sha256",serialize=FALSE),bootstrap_pin$sha256),
    "ENS101_PROGRAM_REFERENCE_BOOTSTRAP_SOURCE_PIN")
  wanted <- c("gse205906_ens101_reference_need_v1","gse205906_ens101_reference_hash_v1",
    "gse205906_ens101_reference_import_v1","gse205906_ens101_reference_common_pins_v1",
    "gse205906_ens101_reference_support_v1")
  declarations <- as.list(parse(text=rawToChar(source_bundle$reference_bootstrap),keep.source=FALSE))
  selected <- Filter(function(x) is.call(x) && identical(x[[1L]],as.name("<-")) &&
    is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted,declarations)
  need(setequal(vapply(selected,function(x)as.character(x[[2L]]),character(1)),wanted) &&
    length(selected)==length(wanted),"ENS101_PROGRAM_REFERENCE_PURE_BOOTSTRAP_COVERAGE")
  bootstrap <- new.env(parent=baseenv())
  for (x in selected) {
    need(is.call(x[[3L]]) && identical(x[[3L]][[1L]],as.name("function")),"ENS101_PROGRAM_REFERENCE_PURE_BOOTSTRAP_ONLY")
    eval(x,bootstrap)
  }
  lockEnvironment(bootstrap,bindings=TRUE)
  # Only pure source/config/input helpers are imported from the gene module;
  # neither its refit loader nor any core compute function is available here.
  pins <- c(list(old_reference=list(bytes=25221L,sha256="4e0f3cd999b327e06de1388f8fe92472a1af20bf62e9d71d7cea24ed3b44d493"),
    reference_bootstrap=bootstrap_pin),bootstrap$gse205906_ens101_reference_common_pins_v1())
  support <- bootstrap$gse205906_ens101_reference_support_v1(config_bytes,legacy_config_bytes,source_bundle,
    prefit,mapping,pins,c("digest","jsonlite","data.table","edgeR","limma","fgsea","decoupleR"))
  derived <- gse205906_ens101_program_reference_derive_v1(rawToChar(source_bundle$old_reference),support$config_raw_sha256)
  engine <- bootstrap$gse205906_ens101_reference_import_v1(charToRaw(derived$text),
    "prepare_gse205906_ens101_program_reference_engine_v1")
  prepare <- function(raw_inputs) {
    p <- support$prepare(raw_inputs)
    refit <- function() {
      p$guard()
      value <- engine$prepare_gse205906_ens101_program_reference_engine_v1(
        c(source_bundle["helpers"],list(full_config=config_bytes)),support$cfg,
        p$snapshot$count_matrix,p$projection$computational_interface,p$snapshot$sample_ledger,
        p$snapshot$gene_resources,p$snapshot$orthology_table,"PRODUCTION",support$contract,p$snapshot)
      p$guard()
      value$raw_input_identities <- p$raw_input_identities
      value$actual_file_producer_runtime_authentication_supplied_here <- FALSE
      value$legacy_omitted_animal_is_library_identifier_alias <- TRUE
      value
    }
    list(inspect=p$inspect,refit=refit)
  }
  list(status="NEEDS_REVIEW_NUMERICAL_PROGRAM_REFERENCE_NOT_RUN",prepare=prepare,
    source_pins=pins,source_amendments=derived$amendments,derived_source_sha256=derived$derived_source_sha256,
    packages=support$packages,config_raw_sha256=support$config_raw_sha256,
    numerical_reference_run=FALSE,model_run=FALSE,production_go=FALSE,checked_consumer_implemented=FALSE,
    limitation="Deferred independent original program/activity/full-LOAO reference body. No core compute/reported table is used to generate expectations. Shared input helpers/statistical libraries permit common-mode errors; SOURCE loading does not establish numerical comparison, independent personnel, real input/runtime custody or model acceptance.")
}
