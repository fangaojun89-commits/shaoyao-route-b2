# Prospective program-family contract. Does not authorize public model execution.
# Public expected identities come from pinned pre-fit/mapping evidence, not effects.

gse205906_amend_model_program_binding_v1 <- function(text) {
  replace_once <- function(old,new) {
    matches <- gregexpr(old,text,fixed=TRUE)[[1L]]
    if(length(matches)!=1L || matches[[1L]]<1L) stop("BOUND_MODEL_AMENDMENT_ANCHOR",call.=FALSE)
    text <<- paste0(substr(text,1L,matches[[1L]]-1L),new,
      substr(text,matches[[1L]]+nchar(old),nchar(text)))
  }
  replace_once("  y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = \"TMM\")",
    paste0("  y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = \"TMM\")\n",
      "  bound_programs$check_genes(rownames(y$counts), colnames(y$counts))"))
  replace_once("v2r4_unsigned_tests <- function(resource, v, design_matrix, contrasts) {",
    paste0("bound_programs$check_unsigned(unsigned_resources, ids)\n\n",
      "v2r4_unsigned_tests <- function(resource, v, design_matrix, contrasts) {"))
  replace_once("  if (enforce_expected) {",
    paste0("  bound_programs$check_activity(prog, doro, pcov, dcov, rownames(expr), colnames(expr))\n",
      "  if (enforce_expected) {"))
  replace_once(paste0("v2r4_assert(all(loao_summary$n_omissions == 9L) && all(loao_summary$n_estimable >= min_estimable),\n",
      "          \"LOAO_INSUFFICIENT_ESTIMABLE_OMISSIONS\")"),
    paste0("v2r4_assert(all(loao_summary$n_omissions == 9L), \"LOAO_OMISSION_COUNT_INVALID\")\n",
      "loao_summary <- bound_programs$annotate_loao(loao_summary, min_estimable)"))
  text
}

gse205906_program_contract_v1 <- function(spec, input_identity_function) {
  amender <- gse205906_amend_model_program_binding_v1
  need <- function(x,label) if(!isTRUE(x)) stop(label,call.=FALSE)
  hash <- function(x) digest::digest(serialize(x,NULL,version=3),"sha256",serialize=FALSE)
  members_hash <- function(x) digest::digest(charToRaw(enc2utf8(paste(sort(unique(x),method="radix"),collapse="\n"))),"sha256",serialize=FALSE)
  spec <- unserialize(serialize(spec,NULL,version=3))
  immutable <- hash(spec)
  need(spec$domain %in% c("SYNTHETIC","PRODUCTION") && is.function(input_identity_function),"PROGRAM_BINDING_SPEC")
  ids <- spec$sample_ids
  families <- c("Reactome","SynGO_sensitivity","PROGENy","DoRothEA_AB")
  phases <- c("NONE",ids)
  need(length(ids)==9L && !anyDuplicated(ids) && identical(names(spec$phases),phases),"PROGRAM_BINDING_PHASES")
  for(phase in phases) {
    packet <- spec$phases[[phase]]
    need(is.character(packet$gene_ids) && !anyDuplicated(packet$gene_ids) &&
      identical(names(packet$families),families),"PROGRAM_BINDING_UNIVERSE")
    for(collection in families) {
      family <- packet$families[[collection]]
      need(is.list(family) && (length(family)==0L || (!is.null(names(family)) && !anyDuplicated(names(family)))) &&
        all(vapply(family,function(p) is.list(p) && length(p$n)==1L && is.finite(p$n) && p$n>=10 && p$n<=500 &&
          is.character(p$sha256) && length(p$sha256)==1L && grepl("^[0-9a-f]{64}$",p$sha256),logical(1))),"PROGRAM_BINDING_FAMILY_SPEC")
    }
  }
  check_integrity <- function() need(identical(hash(spec),immutable),"PROGRAM_BINDING_SPEC_MUTATION")
  check_inputs <- function(inputs,domain) {
    check_integrity()
    need(identical(domain,spec$domain) && identical(input_identity_function(inputs),spec$input_sha256),"PROGRAM_BINDING_INPUT_IDENTITY")
    invisible(TRUE)
  }
  effective_config <- function(cfg) {
    check_integrity()
    output <- unserialize(serialize(cfg,NULL,version=3))
    # Only the observed eligible cardinality changes; all filters, networks,
    # source cardinalities, contrasts, seed, methods and thresholds are retained.
    branch <- if(spec$domain=="SYNTHETIC") "synthetic_contract" else "analysis"
    output[[branch]]$program_domains$DoRothEA_AB$eligible_tested_program_n <- length(spec$phases$NONE$families$DoRothEA_AB)
    output
  }
  new_session <- function() {
    seen_genes <- character();seen_families <- character();annotated <- FALSE;robustness <- NULL
    phase_from_samples <- function(samples) {
      if(identical(samples,ids)) return("NONE")
      omitted <- setdiff(ids,samples)
      need(length(omitted)==1L && identical(samples,ids[ids!=omitted]),"PROGRAM_BINDING_SAMPLE_ORDER")
      omitted
    }
    check_genes <- function(gene_ids,samples) {
      check_integrity();phase <- phase_from_samples(samples)
      need(!phase %in% seen_genes && identical(gene_ids,spec$phases[[phase]]$gene_ids),"PROGRAM_BINDING_GENE_UNIVERSE")
      seen_genes <<- c(seen_genes,phase);invisible(TRUE)
    }
    check_family <- function(collection,members,samples) {
      check_integrity();phase <- phase_from_samples(samples);key <- paste(phase,collection,sep="|")
      need(collection %in% families && !key %in% seen_families && is.list(members),"PROGRAM_BINDING_FAMILY_CALL")
      expected <- spec$phases[[phase]]$families[[collection]]
      observed_ids <- if(length(members)) names(members) else character()
      expected_ids <- if(length(expected)) names(expected) else character()
      need((!length(members) || !is.null(names(members))) && !anyDuplicated(observed_ids) &&
        setequal(observed_ids,expected_ids),"PROGRAM_BINDING_PROGRAM_IDS")
      for(id in expected_ids) {
        targets <- members[[id]]
        need(is.character(targets) && !anyNA(targets) && all(nzchar(targets)) &&
          length(unique(targets))==expected[[id]]$n && identical(members_hash(targets),expected[[id]]$sha256),
          "PROGRAM_BINDING_MEMBERSHIP")
      }
      seen_families <<- c(seen_families,key);invisible(TRUE)
    }
    check_unsigned <- function(resources,samples) {
      need(identical(names(resources),c("Reactome","SynGO_sensitivity")),"PROGRAM_BINDING_UNSIGNED_COLLECTIONS")
      for(collection in names(resources)) check_family(collection,resources[[collection]]$members,samples)
      invisible(TRUE)
    }
    check_activity <- function(prog,doro,pcov,dcov,symbols,samples) {
      for(collection in c("PROGENy","DoRothEA_AB")) {
        net <- if(collection=="PROGENy") prog else doro
        coverage <- if(collection=="PROGENy") pcov else dcov
        selected <- as.character(coverage$source[coverage$tested])
        members <- lapply(selected,function(id) unique(intersect(net$target[net$source==id],symbols)))
        names(members) <- selected
        check_family(collection,members,samples)
      }
      invisible(TRUE)
    }
    annotate_loao <- function(summary,minimum) {
      check_integrity();need(!annotated && identical(as.integer(minimum),7L),"PROGRAM_BINDING_LOAO_THRESHOLD")
      keys <- c("collection","program_id","contrast")
      need(is.data.frame(summary) && all(c(keys,"n_omissions","n_estimable") %in% names(summary)),"PROGRAM_BINDING_LOAO_SCHEMA")
      x <- data.table::as.data.table(data.table::copy(summary))
      expected <- data.table::rbindlist(lapply(c("PROGENy","DoRothEA_AB"),function(collection) {
        out <- data.table::CJ(program_id=names(spec$phases$NONE$families[[collection]]),contrast=c("local_disease","intervention"),unique=TRUE)
        out$collection <- collection;out
      }))
      got <- x[,keys,with=FALSE];expected <- expected[,keys,with=FALSE]
      data.table::setorderv(got,keys);data.table::setorderv(expected,keys)
      need(!anyDuplicated(x[,keys,with=FALSE]) && identical(as.data.frame(got),as.data.frame(expected)),"PROGRAM_BINDING_LOAO_FULL_GRID")
      need(all(is.finite(x$n_estimable)) && all(x$n_estimable==floor(x$n_estimable)) &&
        all(x$n_estimable>=0 & x$n_estimable<=9) && all(x$n_omissions==9),"PROGRAM_BINDING_LOAO_COUNTS")
      available <- vapply(seq_len(nrow(x)),function(i) sum(vapply(ids,function(phase)
        x$program_id[[i]] %in% names(spec$phases[[phase]]$families[[x$collection[[i]]]]),logical(1))),integer(1))
      need(all(x$n_estimable<=available),"PROGRAM_BINDING_ESTIMABILITY_EXCEEDS_MEMBERSHIP")
      x$minimum_estimable_omissions <- rep(7L,nrow(x))
      x$robustness_pass <- x$n_estimable>=7L
      x$robustness_status <- ifelse(x$robustness_pass,"PASS_MINIMUM_ESTIMABLE_OMISSIONS","FAIL_INSUFFICIENT_ESTIMABLE_OMISSIONS")
      robustness <<- list(all_programs_pass=all(x$robustness_pass),failed_program_contrast_n=sum(!x$robustness_pass))
      annotated <<- TRUE;x
    }
    finish <- function() {
      expected_family_keys <- c(paste("NONE",families,sep="|"),
        unlist(lapply(ids,function(id) paste(id,c("PROGENy","DoRothEA_AB"),sep="|")),use.names=FALSE))
      need(setequal(seen_genes,phases) && length(seen_genes)==10L &&
        setequal(seen_families,expected_family_keys) && length(seen_families)==22L && annotated,"PROGRAM_BINDING_INCOMPLETE_SESSION")
      list(status="PASS_BOUND_FAMILIES_AND_COMPLETE_REPORTING_ONLY",binding_sha256=immutable,
        domain=spec$domain,full_family_n=as.list(lengths(spec$phases$NONE$families)),
        checked_gene_universes=10L,checked_program_families=22L,robustness=robustness,
        production_go=FALSE,independent_reproduction=FALSE,public_registration=FALSE)
    }
    list(check_genes=check_genes,check_family=check_family,check_unsigned=check_unsigned,
      check_activity=check_activity,annotate_loao=annotate_loao,finish=finish)
  }
  list(check_inputs=check_inputs,effective_config=effective_config,new_session=new_session,amend_model=amender,
    binding_sha256=immutable,domain=spec$domain,production_go=FALSE,
    public_execution_accepted=FALSE,full_family_n=as.list(lengths(spec$phases$NONE$families)))
}

gse205906_public_program_contract_v1 <- function(root,input_identity_function) {
  need <- function(x,label) if(!isTRUE(x)) stop(label,call.=FALSE)
  paths <- c(
    "results/diagnostics/GSE205906_public_mapping_v1_first_20260909/observed_program_universe_v1.json"="0ba8cb9d634d9b9c55da3098cb60e665024433c57d886a29e60c4814c031e1ff",
    "results/diagnostics/GSE205906_public_mapping_v1_first_20260909/worker_provisional.json"="5560ee1d72d91b15ecfe76ce30ec1dfcf11f311fc81d053874cc75e11c9d4cc5",
    "results/diagnostics/GSE205906_public_prefit_v1_first_20260909/worker_provisional.json"="7830f6af45ab008d85f3722a459841297f7dbadd9f99d716b142447a655ca95b")
  objects <- lapply(names(paths),function(path) {
    full <- file.path(root,path);size <- file.info(full)$size
    need(length(size)==1L && is.finite(size) && size>0 && size<15000000 && !nzchar(Sys.readlink(full)),"PUBLIC_BINDING_SOURCE")
    con <- file(full,"rb");on.exit(close(con),add=TRUE);raw <- readBin(con,"raw",as.integer(size)+1L)
    need(length(raw)==size && identical(digest::digest(raw,"sha256",serialize=FALSE),paths[[path]]),"PUBLIC_BINDING_SOURCE_IDENTITY")
    jsonlite::fromJSON(rawToChar(raw),simplifyVector=FALSE)
  })
  universe <- objects[[1L]];mapping <- objects[[2L]];prefit <- objects[[3L]]
  need(identical(universe$status,"OBSERVED_MAPPING_IDENTITY_NOT_MODEL_RELEASE") &&
    identical(universe$production_go,FALSE) && identical(mapping$input_sha256,prefit$input_sha256) &&
    identical(mapping$input_sha256,universe$input_sha256),"PUBLIC_BINDING_EVIDENCE_SCOPE")
  families <- c("Reactome","SynGO_sensitivity","PROGENy","DoRothEA_AB")
  ids <- unlist(prefit$count_validation$sample_ids,use.names=FALSE)
  phases <- lapply(seq_along(prefit$prefit$runs),function(i) {
    p <- prefit$prefit$runs[[i]];m <- mapping$mapping$runs[[i]]
    need(identical(p$omission,m$omission),"PUBLIC_BINDING_OMISSION_ORDER")
    cov <- m$coverage
    out <- lapply(families,function(collection) {
      selected <- Filter(function(x) identical(x$collection,collection) && isTRUE(x$eligible),cov)
      packet <- lapply(selected,function(x) list(n=as.integer(x$tested_feature_n),sha256=x$retained_members_sha256))
      names(packet) <- vapply(selected,`[[`,character(1),"program_id");packet
    });names(out) <- families
    list(gene_ids=unlist(p$retained_gene_ids,use.names=FALSE),families=out)
  });names(phases) <- c("NONE",ids)
  for(collection in families) {
    recorded <- universe$full_families[[collection]]$programs
    expected <- stats::setNames(lapply(recorded,function(x) list(n=as.integer(x$tested_feature_n),sha256=x$retained_members_sha256)),
      vapply(recorded,`[[`,character(1),"program_id"))
    observed <- phases$NONE$families[[collection]]
    need(identical(observed[sort(names(observed))],expected[sort(names(expected))]),"PUBLIC_BINDING_FULL_MANIFEST")
  }
  spec <- list(domain="PRODUCTION",sample_ids=ids,input_sha256=unlist(universe$input_sha256,use.names=TRUE),phases=phases)
  gse205906_program_contract_v1(spec,input_identity_function)
}

gse205906_synthetic_program_contract_v1 <- function(inputs,cfg,input_identity_function) {
  need <- function(x,label) if(!isTRUE(x)) stop(label,call.=FALSE)
  x <- inputs$count_matrix;ann <- inputs$gene_annotation
  ids <- unlist(cfg$design$sample_ids_in_order,use.names=FALSE)
  need(is.matrix(x) && nrow(x)>=20L && nrow(x)<=1000L && identical(colnames(x),ids) &&
    all(grepl("^SYNTHETIC_",rownames(x))) && all(is.finite(x)) && all(x>=0 & x==floor(x)),"BINDING_SYNTHETIC_FIXTURE_ONLY")
  need(identical(as.character(ann$GeneID),rownames(x)) && all(grepl("^SYN_",inputs$sample_ledger$biosample)),"BINDING_SYNTHETIC_LABELS")
  symbols <- trimws(as.character(ann$gene))
  need(!anyNA(symbols) && !anyDuplicated(symbols) && all(nzchar(symbols)),"BINDING_SYNTHETIC_SYMBOLS")
  groups <- factor(unlist(cfg$design$groups_in_order,use.names=FALSE),levels=unlist(cfg$design$group_levels,use.names=FALSE))
  map <- as.data.frame(inputs$orthology_table)
  map <- unique(map[map$group_method=="NCBI Ortholog" & map$other_tax_id==10116L,c("human_symbol","other_symbol")])
  need(!anyDuplicated(map$human_symbol) && !anyDuplicated(map$other_symbol),"BINDING_SYNTHETIC_ORTHOLOGY")
  resources <- inputs$gene_resources
  networks <- list()
  make_sets <- function(table,id,target,human) {
    table <- as.data.frame(table);program <- as.character(table[[id]]);gene <- as.character(table[[target]])
    if(human) gene <- as.character(map$other_symbol[match(gene,map$human_symbol)])
    out <- lapply(sort(unique(program)),function(p) unique(gene[program==p & !is.na(gene)]))
    names(out) <- sort(unique(program));out
  }
  networks$Reactome <- make_sets(resources$reactome,"term_id","gene_symbol",TRUE)
  networks$SynGO_sensitivity <- make_sets(resources$syngo_sensitivity,"term_id","gene_symbol",TRUE)
  networks$PROGENy <- make_sets(resources$progeny_human,"pathway","gene_symbol",TRUE)
  networks$DoRothEA_AB <- make_sets(resources$dorothea_rat,"tf","target",FALSE)
  phase_packets <- lapply(0:9,function(omit) {
    indices <- if(omit==0L) 1:9 else setdiff(1:9,omit)
    g <- groups[indices];d <- stats::model.matrix(~0+g)
    initial <- edgeR::DGEList(counts=x[,indices,drop=FALSE])
    retained <- edgeR::filterByExpr(initial,design=d) &
      rowSums(edgeR::cpm(initial,log=FALSE)>=cfg$analysis$minimum_cpm)>=cfg$analysis$minimum_samples_at_cpm
    need(sum(retained)>=cfg$synthetic_contract$minimum_retained_genes,"BINDING_SYNTHETIC_FILTER")
    families <- lapply(networks,function(net) {
      members <- lapply(net,function(z) sort(unique(intersect(z,symbols[retained])),method="radix"))
      members <- members[lengths(members)>=10L & lengths(members)<=500L]
      lapply(members,function(z) list(n=length(z),sha256=digest::digest(charToRaw(enc2utf8(paste(z,collapse="\n"))),"sha256",serialize=FALSE)))
    })
    list(gene_ids=rownames(x)[retained],families=families)
  });names(phase_packets) <- c("NONE",ids)
  spec <- list(domain="SYNTHETIC",sample_ids=ids,input_sha256=input_identity_function(inputs),phases=phase_packets)
  gse205906_program_contract_v1(spec,input_identity_function)
}
