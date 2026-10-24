# Versioned independent gene-reference wiring; old refit spans are retained.
# Loading/prepare/inspect are effect-free. Deferred refit performs science and
# advances RNG: it is NOT invoked by SOURCE tests and grants no model authority.

gse205906_ens101_reference_need_v1 <- function(ok, category) {
  if (!isTRUE(ok)) stop(category, call. = FALSE)
}

gse205906_ens101_reference_hash_v1 <- function(x) {
  digest::digest(serialize(x, NULL, version = 3L), "sha256", serialize = FALSE)
}

gse205906_ens101_reference_import_v1 <- function(bytes, wanted) {
  need <- gse205906_ens101_reference_need_v1
  declarations <- as.list(parse(text = rawToChar(bytes), keep.source = FALSE))
  selected <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted, declarations)
  found <- vapply(selected, function(x) as.character(x[[2L]]), character(1))
  need(setequal(found, wanted) && !anyDuplicated(found), "ENS101_REFERENCE_DECLARATION_COVERAGE")
  private <- new.env(parent = baseenv())
  for (x in selected) {
    need(is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), "ENS101_REFERENCE_FUNCTION_ONLY")
    eval(x, private)
  }
  lockEnvironment(private, bindings = TRUE)
  private
}

gse205906_ens101_reference_common_pins_v1 <- function() {
  list(
    helpers = list(bytes = 175131L, sha256 = "407bbf2acda7db749f33a4f679564e5cddd0c758f1d0ad9f45e870b2c25432c3"),
    model = list(bytes = 72490L, sha256 = "b8b274667dcb68e5c8a729e10f4df46d9e596371be4ee6a8a99efc7be1286cd2"),
    binding = list(bytes = 15697L, sha256 = "2aaa386df70c0031933c244eafa895d63c5dbd7a7be77277cfd4a57aceba538d"),
    annotation_adapter = list(bytes = 9943L, sha256 = "9f223ec6163a4ac6aca5e85edccac9594f2c5a631889c76bd81f6d560f187470"),
    prefit = list(bytes = 18762L, sha256 = "b394a9e301db4e17a3c0495ada7071337198e8aae15eccf82d87cfa7b32f9bc3"),
    identity_adapter = list(bytes = 17175L, sha256 = "ac070364dca4369f9a9cb196550a8c6f3daced3114184f1572d5274ec03261e6"),
    config_adapter = list(bytes = 7852L, sha256 = "29ffe9705e675058c006d808305349275979efb58eeab2ea4236e5eb82fd339a"),
    program_contract = list(bytes = 12746L, sha256 = "bd86c229834a35317a701bfd53f727943d5238a44b99a94a549e3be8b5ba3adc")
  )
}

gse205906_ens101_reference_support_v1 <- function(config_bytes, legacy_config_bytes,
    source_bundle, prefit, mapping, pins, packages) {
  need <- gse205906_ens101_reference_need_v1
  hash <- gse205906_ens101_reference_hash_v1
  import <- gse205906_ens101_reference_import_v1
  need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE) &&
       requireNamespace("data.table", quietly = TRUE), "ENS101_REFERENCE_EXISTING_DEPENDENCIES")
  need(is.list(source_bundle) && identical(names(source_bundle), names(pins)), "ENS101_REFERENCE_EXACT_SOURCE_ROLES")
  for (role in names(pins)) need(is.raw(source_bundle[[role]]) &&
    length(source_bundle[[role]]) == pins[[role]]$bytes &&
    identical(digest::digest(source_bundle[[role]], "sha256", serialize = FALSE), pins[[role]]$sha256),
    "ENS101_REFERENCE_SOURCE_BYTES_OR_SHA_DRIFT")
  need(is.raw(config_bytes) && length(config_bytes) == 10612L &&
    identical(digest::digest(config_bytes,"sha256",serialize=FALSE),
      "ee6e22592b06ecdb31707234e97dc55b7b66cbad8dd11c364187af0624ece980"),
    "ENS101_REFERENCE_CURRENT_FULL_CONFIG_RAW_PIN")
  config <- import(source_bundle$config_adapter, c("gse205906_ens101_config_need_v1",
    "gse205906_ens101_config_hash_v1", "gse205906_ens101_legacy_config_v1",
    "gse205906_ens101_full_model_config_v1", "gse205906_ens101_check_full_model_config_v1",
    "gse205906_ens101_prefit_view_v1"))
  cfg <- config$gse205906_ens101_check_full_model_config_v1(config_bytes,legacy_config_bytes)
  view <- config$gse205906_ens101_prefit_view_v1(cfg,legacy_config_bytes)
  pref <- import(source_bundle$prefit,c("gse205906_ens101_prefit_need_v1",
    "gse205906_ens101_prefit_config_v1","gse205906_ens101_prefit_ledger_v1"))
  pref$gse205906_ens101_prefit_config_v1(view)
  need(identical(prefit$config_object_sha256,hash(view)),"ENS101_REFERENCE_PREFIT_CONFIG_VIEW")
  identity <- import(source_bundle$identity_adapter,c("gse205906_ens101_core_need_v1",
    "gse205906_ens101_core_hash_v1","gse205906_ens101_raw_input_identities_v1"))
  # Only these three pure identity declarations are imported. Neither the new
  # core loader/compute nor any reported result is accessible to the refits.
  program <- import(source_bundle$program_contract,c("gse205906_ens101_program_hash_v2",
    "gse205906_ens101_program_contract_v2"))
  contract <- program$gse205906_ens101_program_contract_v2(prefit,mapping,
    source_bundle[c("binding","model","helpers")],identity$gse205906_ens101_raw_input_identities_v1)
  adapter <- import(source_bundle$annotation_adapter,c("gse205906_ens101_need_v1",
    "gse205906_ens101_character_v1","gse205906_ens101_annotation_v1"))
  need(identical(as.character(getRversion()),cfg$software$R),"ENS101_REFERENCE_R_VERSION")
  versions <- vapply(packages,function(package) {
    need(requireNamespace(package,quietly=TRUE),"ENS101_REFERENCE_PACKAGE_MISSING")
    version <- as.character(utils::packageVersion(package))
    need(identical(version,cfg$software[[package]]),"ENS101_REFERENCE_PACKAGE_VERSION")
    version
  },character(1))
  prepare <- function(raw_inputs) {
    contract$check_inputs(raw_inputs,"PRODUCTION") # raw identity BEFORE projection
    snapshot <- unserialize(serialize(raw_inputs,NULL,version=3L))
    snapshot_hash <- hash(snapshot)
    identities <- identity$gse205906_ens101_raw_input_identities_v1(snapshot)
    pref$gse205906_ens101_prefit_ledger_v1(snapshot$sample_ledger)
    counts <- snapshot$count_matrix
    need(identical(dim(counts),c(32883L,9L)) &&
      identical(colnames(counts),paste0("GSM",6235292:6235300)) &&
      identical(rownames(counts),prefit$counts_gene_ids) && all(is.finite(counts)) &&
      all(counts>=0) && all(counts==floor(counts)) && all(counts<=.Machine$integer.max),
      "ENS101_REFERENCE_FIXED_INTEGER_COUNTS")
    projection <- adapter$gse205906_ens101_annotation_v1(snapshot$gene_annotation,prefit$counts_gene_ids,"FIXED_ENS101")
    need(identical(projection,prefit$annotation_projection),"ENS101_REFERENCE_CURRENT_PREFIT_PROJECTION")
    projection_hash <- hash(projection)
    guard <- function() {
      need(identical(hash(snapshot),snapshot_hash) && identical(hash(projection),projection_hash),
        "ENS101_REFERENCE_CAPTURED_INPUT_OR_PROJECTION_DRIFT")
      contract$check_inputs(snapshot,"PRODUCTION")
      invisible(TRUE)
    }
    inspect <- function() {
      guard()
      list(status="SOURCE_INPUT_PREPARED_NUMERICAL_REFERENCE_NOT_RUN",raw_input_identities=identities,
        raw_snapshot_object_sha256=snapshot_hash,annotation_projection_object_sha256=projection_hash,
        current_prefit_object_sha256=contract$current_prefit_object_sha256,
        current_mapping_object_sha256=contract$current_mapping_object_sha256,
        count_gene_rows=32883L,library_n=9L,rows_removed=0L,symbol_gene_rows_collapsed=FALSE,
        observation_unit="NINE_DEPOSITED_LIBRARIES_CONDITIONAL_INDEPENDENCE",
        legacy_omitted_animal_is_library_identifier_alias=TRUE,
        numerical_reference_run=FALSE,model_run=FALSE,production_go=FALSE,
        actual_file_producer_runtime_authentication_supplied_here=FALSE)
    }
    list(snapshot=snapshot,projection=projection,guard=guard,inspect=inspect,raw_input_identities=identities)
  }
  list(cfg=cfg,contract=contract,prepare=prepare,packages=versions,pins=pins,
    config_raw_sha256=digest::digest(config_bytes,"sha256",serialize=FALSE),
    config_object_sha256=hash(cfg),prefit_view=view)
}

gse205906_ens101_derive_gene_reference_v1 <- function(old_text, config_sha256) {
  need <- gse205906_ens101_reference_need_v1
  need(is.character(old_text) && length(old_text)==1L &&
    identical(digest::digest(charToRaw(old_text),"sha256",serialize=FALSE),
      "4188175759fc7db35b18269bf54e92be0b217517226a676ed4f444aaabf8968f"),
    "ENS101_GENE_REFERENCE_FROZEN_SOURCE")
  text <- old_text; amendments <- list()
  replace <- function(old,new,label) {
    hits <- gregexpr(old,text,fixed=TRUE)[[1L]]
    need(length(hits)==1L && hits[[1L]]>0L,"ENS101_GENE_REFERENCE_EXACT_ANCHOR")
    text <<- paste0(substr(text,1L,hits[[1L]]-1L),new,substr(text,hits[[1L]]+nchar(old),nchar(text)))
    amendments[[label]] <<- list(old_sha256=digest::digest(charToRaw(old),"sha256",serialize=FALSE),
      new_sha256=digest::digest(charToRaw(new),"sha256",serialize=FALSE))
  }
  start <- regexpr("prepare_gse205906_gene_validation <- function",text,fixed=TRUE)[[1L]]
  end <- regexpr('  fail <- function(message)',text,fixed=TRUE)[[1L]]
  replace(substr(text,start,end-1L),paste0(
    "prepare_gse205906_ens101_gene_reference_engine_v1 <- function(source_bytes, cfg, count_matrix,\n",
    "    gene_annotation, sample_ledger, domain = \"PRODUCTION\") {\n"),"explicit_source_cfg_signature")
  replace('      !domain %in% c("SYNTHETIC", "PRODUCTION")) fail("VALIDATOR_DOMAIN_INVALID")',
    '      !identical(domain, "PRODUCTION")) fail("VALIDATOR_DOMAIN_INVALID")',"production_settings_only")
  replace(paste0('    "config/GSE205906_bulk_model_v2_r4.json" =\n',
    '      "d8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c"'),
    paste0('    "config/GSE205906_Ensembl101_full_model_v1_20261006.json" =\n      "',config_sha256,'"'),"new_config_raw_pin")
  a <- regexpr("  read_source <- function(relative, expected) {",text,fixed=TRUE)[[1L]]
  b <- regexpr("  sources <- Map(read_source, names(pins), unname(pins))",text,fixed=TRUE)[[1L]]
  replace(substr(text,a,b-1L),paste0('  read_source <- function(relative, expected) {\n',
    '    roles <- c("code/R/validate_gse205906_bulk_model_v2_r4.R"="validator",\n',
    '      "code/R/gse205906_bulk_model_v2_r4_lib.R"="helpers",\n',
    '      "config/GSE205906_Ensembl101_full_model_v1_20261006.json"="full_config")\n',
    '    if (!relative %in% names(roles)) fail("ENS101_GENE_REFERENCE_SOURCE_ROLE")\n',
    '    payload <- source_bytes[[roles[[relative]]]]\n',
    '    if (!is.raw(payload) || !identical(digest::digest(payload,"sha256",serialize=FALSE),expected))\n',
    '      fail("VALIDATOR_SOURCE_IDENTITY_MISMATCH")\n    rawToChar(payload)\n  }\n'),"no_path_reopen")
  replace("  cfg <- jsonlite::fromJSON(sources[[3L]], simplifyVector = FALSE)",
    "  # cfg is the authenticated full-config object supplied by the new loader.","explicit_full_config")
  replace("  ledger <- helpers$v2r4_validate_ledger(data.table::copy(sample_ledger), cfg)",
    "  ledger <- data.table::copy(sample_ledger) # Original 53-column input guard was applied before projection.","original_53_field_library_ledger")
  # Filtering/TMM/refits/seed/comparator code and tolerances are not edited.
  list(text=text,amendments=amendments,
    derived_source_sha256=digest::digest(charToRaw(text),"sha256",serialize=FALSE))
}

load_gse205906_ens101_gene_reference_v1 <- function(config_bytes,legacy_config_bytes,
    source_bundle,prefit,mapping) {
  pins <- c(list(
    old_reference=list(bytes=9590L,sha256="4188175759fc7db35b18269bf54e92be0b217517226a676ed4f444aaabf8968f"),
    validator=list(bytes=104432L,sha256="c63296157a833611ef92b17f1cbb6b27ab3c23c1028cab79fd706aaab564225d")),
    gse205906_ens101_reference_common_pins_v1())
  support <- gse205906_ens101_reference_support_v1(config_bytes,legacy_config_bytes,source_bundle,
    prefit,mapping,pins,c("data.table","digest","jsonlite","edgeR","limma","DESeq2"))
  derived <- gse205906_ens101_derive_gene_reference_v1(rawToChar(source_bundle$old_reference),support$config_raw_sha256)
  engine <- gse205906_ens101_reference_import_v1(charToRaw(derived$text),
    "prepare_gse205906_ens101_gene_reference_engine_v1")
  prepare <- function(raw_inputs) {
    p <- support$prepare(raw_inputs)
    refit <- function() {
      p$guard()
      value <- engine$prepare_gse205906_ens101_gene_reference_engine_v1(
        c(source_bundle[c("validator","helpers")],list(full_config=config_bytes)),support$cfg,
        p$snapshot$count_matrix,p$projection$computational_interface,p$snapshot$sample_ledger,"PRODUCTION")
      p$guard()
      value$raw_input_identities <- p$raw_input_identities
      value$actual_file_producer_runtime_authentication_supplied_here <- FALSE
      value$legacy_omitted_animal_is_library_identifier_alias <- TRUE
      value
    }
    list(inspect=p$inspect,refit=refit)
  }
  list(status="NEEDS_REVIEW_NUMERICAL_GENE_REFERENCE_NOT_RUN",prepare=prepare,
    source_pins=pins,source_amendments=derived$amendments,derived_source_sha256=derived$derived_source_sha256,
    packages=support$packages,config_raw_sha256=support$config_raw_sha256,
    numerical_reference_run=FALSE,model_run=FALSE,production_go=FALSE,checked_consumer_implemented=FALSE,
    limitation="Deferred separate original gene refits; no core compute/reported result is used to obtain expectations. Shared input helpers/covariance/packages permit common-mode errors and do not imply independent personnel, input/runtime authority or complete checked acceptance.")
}
