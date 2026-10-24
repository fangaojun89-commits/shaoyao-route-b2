# Thin, explicitly invoked IO caller; source() defines functions only.
# Usage: Rscript --vanilla THIS.R PROJECT_ROOT PRODUCER_SOURCE FRESH_OUTPUT
#   [--admit | --execute] --python /absolute/python3
#   --quiet-window-confirmed --reserved-growth-bytes N --disk-budget-bytes N
# Default is dry-run: no scientific input reads, writes, import or native calls.
# --admit authenticates/converts inputs but does not fit. --execute calls the
# frozen core once. There is no resume/old driver/LOAO/null/top-k/comparison here.
# Environment is the existing frozen project R/library; no install or restore.

g205run_need <- function(ok, category) {
  if (!isTRUE(ok)) stop(category, call. = FALSE)
}

g205run_constants <- function() {
  n <- "/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01"
  f <- "/root/autodl-tmp/GSE205906_full_read_counts_fasp_v18e_20261005/formal"
  a <- "/root/autodl-tmp/GSE205906_full_read_reconstruction_156_20261003"
  p <- "/root/autodl-tmp/route_b2_public_analysis/public-input-stage-v1"
  h <- "/root/autodl-tmp/GSE205906_Ensembl101_checked_model_v5_20261006"
  list(run = n, hisat = h,
    producer_sha = "7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d",
    core_sha = "b7cac68a7ef3f5fc63ab727f7498b7b6d0c46b9bde18276affce3c87c4f3a71f",
    writer_sha = "3acd88e536a7808ad750c1d5aaa8d79ff6f1c84510f1bd9e4ed5d92f68363efc",
    parser_sha = "1f0549994bb9aef084c8278fceb336f3bf436abcb589fe4c1755add018dbd61b",
    adapter_sha = "e93aae6769a1173bb7842f213152f56695f5b7b5621b71106c54181496f994ae",
    matrix = paste0(n, "/matrix_01/GSE205906_STAR_fixed9_gene_fragments_s0.tsv"),
    metadata = paste0(f, "/derived/NINE_LIBRARY_ANALYSIS_PREPARATION_OVERLAY_154.tsv"),
    fixed_gene_manifest = paste0(f, "/provenance/run_manifest.json"),
    GTF = paste0(a, "/references/derived/Rattus_norvegicus.Rnor_6.0.101.gtf"),
    gene_resources = paste0(p, "/data/processed/gene_sets/gene_set_resources_frozen_2026-08-18.rds"),
    orthology = paste0(p, "/data/processed/gene_sets/ncbi_human_mouse_rat_one_to_one.tsv"),
    reference_admission = paste0(n, "/reference_admission_01/COMPLETE.json"),
    hisat2_completion = paste0(h, "/COMPLETED"), hisat2_manifest = paste0(h, "/file_manifest.json"),
    hisat2_matrix = paste0(f, "/output/GSE205906_fixed9_gene_fragments_s0.tsv"))
}

g205run_options <- function(args) {
  g205run_need(length(args) >= 3L, "EXPLICIT_PROJECT_PRODUCER_OUTPUT_REQUIRED")
  answer <- list(project = args[[1L]], producer = args[[2L]], output = args[[3L]],
    mode = "dry-run", python = NULL, quiet_window_confirmed = FALSE,
    reserved_growth_bytes = NULL, disk_budget_bytes = NULL)
  i <- 4L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--admit", "--execute")) {
      g205run_need(identical(answer$mode, "dry-run"), "ONE_MODE_ONLY")
      answer$mode <- substring(key, 3L); i <- i + 1L
    } else if (key == "--quiet-window-confirmed") {
      g205run_need(!answer$quiet_window_confirmed, "DUPLICATE_OPTION")
      answer$quiet_window_confirmed <- TRUE; i <- i + 1L
    } else {
      option_names <- c("--python" = "python", "--reserved-growth-bytes" = "reserved_growth_bytes",
        "--disk-budget-bytes" = "disk_budget_bytes")
      g205run_need(key %in% names(option_names) && i < length(args), "UNKNOWN_OR_INCOMPLETE_OPTION")
      target <- unname(option_names[[key]])
      g205run_need(is.null(answer[[target]]), "DUPLICATE_OPTION")
      value <- args[[i + 1L]]
      if (target != "python") {
        g205run_need(grepl("^[0-9]+$", value), "INTEGER_RESOURCE_BYTES_REQUIRED")
        value <- as.numeric(value)
        g205run_need(is.finite(value) && value <= 2^53 &&
          (target != "disk_budget_bytes" || value > 0), "INVALID_RESOURCE_BUDGET")
      }
      answer[[target]] <- value; i <- i + 2L
    }
  }
  g205run_need(!grepl("(^|/|\\\\)[.][.](/|\\\\|$)", answer$output), "NO_OUTPUT_TRAVERSAL")
  if (answer$mode != "dry-run") {
    g205run_need(identical(unname(Sys.info()["sysname"]), "Linux") &&
      is.character(answer$python) && startsWith(answer$python, "/") &&
      file.exists(answer$python) && answer$quiet_window_confirmed &&
      !is.null(answer$reserved_growth_bytes) && !is.null(answer$disk_budget_bytes),
      "EXPLICIT_EXISTING_LINUX_PYTHON_AND_PRODUCER_RESOURCE_CONTRACT_REQUIRED")
    prefix <- paste0(g205run_constants()$run, "/downstream_")
    g205run_need(startsWith(answer$output, prefix) &&
      dirname(answer$output) == g205run_constants()$run && !file.exists(answer$output) &&
      dir.exists(dirname(answer$output)), "FRESH_DIRECT_CLOUD_DOWNSTREAM_LEAF_REQUIRED")
  }
  answer
}

g205run_pin <- function(path) {
  info <- file.info(path)
  g205run_need(nrow(info) == 1L && !is.na(info$size) && !info$isdir &&
    identical(Sys.readlink(path), ""), "REGULAR_NONSYMLINK_INPUT_REQUIRED")
  list(path = path, bytes = unname(info$size),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

g205run_read <- function(path, expected = NULL, max_bytes = 512 * 1024^2) {
  before <- g205run_pin(path)
  g205run_need(before$bytes <= max_bytes && (is.null(expected) || identical(before$sha256, expected)),
    "INPUT_SIZE_OR_SHA_MISMATCH")
  value <- readBin(path, "raw", n = before$bytes)
  g205run_need(length(value) == before$bytes && identical(before, g205run_pin(path)), "INPUT_CHANGED_DURING_READ")
  list(raw = value, pin = before)
}

g205run_json <- function(raw) jsonlite::fromJSON(rawToChar(raw), simplifyVector = FALSE)

# Executed only after --admit/--execute. This bridge invokes ONLY already saved
# producer validation/read functions. Never stream/inflate/reference_dictionary,
# native commands, BAM/FASTQ scans, library(), matrix() or execute(). Small file
# SHA rechecks and raw/BAM/reference/index stat checks are deliberately distinct.
g205run_python <- function() paste(c(
  "import sys,json,hashlib,importlib.util,argparse,csv,io,contextlib",
  "from pathlib import Path",
  "captured_stderr=io.StringIO(); sys.stderr=captured_stderr",
  "try:",
  " source=Path(sys.argv[1]); expected=sys.argv[2]",
  " assert hashlib.sha256(source.read_bytes()).hexdigest()==expected",
  " spec=importlib.util.spec_from_file_location('held_star_producer',source)",
  " p=importlib.util.module_from_spec(spec); spec.loader.exec_module(p)",
  " pins={}",
  " def bind(path):",
  "  q=p.pin(path); pins[q['path']]=q; return q",
  " args=argparse.Namespace(quiet_window_confirmed=True,reserved_growth_bytes=int(sys.argv[3]),disk_budget_bytes=int(sys.argv[4]))",
  " resource=p.resources(args)",
  " answer={'schema':'GSE205906_STAR_SAVED_CUSTODY/v1','accepted':True,'resource':resource,'read_only':True}",
  " if sys.argv[5]=='all':",
  "  meta=p.metadata(); bind(p.META); bind(source)",
  "  admission=p.completed(p.RUN/'reference_admission_01',expected)",
  "  assert admission['plan_sha256']==p.PLAN_SHA and admission['authorization_sha256']==p.AUTH_SHA",
  "  assert admission['resource_before']['quiet_window_confirmed_by_caller'] is True",
  "  assert admission['events'] and all(x['actual_exit']==0 for x in admission['events'])",
  "  ref=admission['reference']; bind(p.RUN/'reference_admission_01/COMPLETE.json')",
  "  for q in admission['small_products']: assert bind(q['path'])==q",
  "  assert ref['reference_pins']==meta['expected_reference_pins_not_freshly_hashed'] and ref['exon_genes']==32883",
  "  assert [p.stat_file(Path(q['path'])) for q in ref['reference_stat']]==ref['reference_stat']",
  "  for tool in admission['tools'].values():",
  "   for q in [tool['binary'],*tool['linked_libraries']]: assert bind(q['path'])==q",
  "  index=p.completed(p.RUN/'index_01',expected)",
  "  assert index['reference_admission']==bind(p.RUN/'reference_admission_01/COMPLETE.json')",
  "  assert index['events'] and all(x['actual_exit']==0 for x in index['events'])",
  "  assert [p.stat_file(Path(q['path'])) for q in index['index_stat']]==index['index_stat']",
  "  bind(p.RUN/'index_01/COMPLETE.json')",
  "  complete=p.completed(p.RUN/'matrix_01',expected); bind(p.RUN/'matrix_01/COMPLETE.json')",
  "  assert complete['dimensions']==[32883,9] and complete['primary_strand']==0",
  "  assert complete['new_models_or_candidate_admission'] is False",
  "  assert complete['original_matrix_preserved']==bind(p.F/'output/GSE205906_fixed9_gene_fragments_s0.tsv')",
  "  assert complete['original_matrix_preserved']['sha256']==p.OLD_MATRIX_SHA",
  "  primary=p.RUN/'matrix_01/GSE205906_STAR_fixed9_gene_fragments_s0.tsv'",
  "  data=primary.read_bytes(); assert all(b in (9,10) or 32<=b<=126 for b in data) and data.endswith(b'\\n')",
  "  rows=list(csv.reader(io.StringIO(data.decode('ascii')),delimiter='\\t'))",
  "  assert rows[0]==['gene_id',*meta['original_sample_order']] and len(rows)==32884 and all(len(r)==10 for r in rows)",
  "  genes=[r[0] for r in rows[1:]]; assert len(set(genes))==32883 and genes==sorted(genes)",
  "  assert all(p.re.fullmatch(r'ENSRNOG[0-9]+',g) for g in genes)",
  "  full=[]; libraries=[]",
  "  for col,lib in enumerate(meta['libraries'],1):",
  "   folder=p.RUN/'libraries'/lib['SRA_run']; done=p.completed(folder,expected)",
  "   assert done['SRA_run']==lib['SRA_run'] and done['label']==lib['label'] and done['technical_reservoir_only'] is False",
  "   assert done['input_pairs']==lib['registered_pairs'] and done['native_counting_primary_strand']==0 and done['new_model_run'] is False",
  "   assert done['index_alignment_sort_not_overlapped'] is True and done['resource_before']['quiet_window_confirmed_by_caller'] is True",
  "   assert len(done['events'])>=8 and all(x['actual_exit']==0 for x in done['events'])",
  "   mates=[Path(x['path']) for x in lib['raw_mates']]; p.check_streams(folder,mates,expected)",
  "   stream={str(folder/'streams'/mate.name/f'pass{ordinal:02d}'/name) for mate in mates for ordinal in (1,2) for name in ('CRC.json','EXIT.json')}",
  "   counts={str(folder/f'counts_s{s}.tsv{suffix}') for s in (0,1,2) for suffix in ('','.summary')}",
  "   required=stream|counts|{str(folder/'QC.json'),str(folder/'alignment.Log.final.out')}",
  "   declared=[q['path'] for q in done['small_products']]; assert len(declared)==16 and len(set(declared))==16 and set(declared)==required",
  "   for q in done['small_products']: assert bind(q['path'])==q",
  "   qc=p.load(folder/'QC.json'); assert qc['input_pairs']==lib['registered_pairs'] and int(qc['STAR_Log_final']['Number of input reads'])==lib['registered_pairs']",
  "   assert qc['samtools_view_actual_exit']==0 and qc['equal_MAPQ_error_calibration_claimed'] is False",
  "   table=p.fc_table(folder/'counts_s0.tsv'); assert set(table)==set(genes)",
  "   assert all(p.re.fullmatch(r'0|[1-9][0-9]*',r[col]) and int(r[col])==table[r[0]] for r in rows[1:])",
  "   for s in ('0','1','2'):",
  "    q=qc['counts'][s]; assigned=qc['assignment'][s]['Assigned']",
  "    assert q['genes']==32883 and q['sum']==assigned and 0<=assigned<=lib['registered_pairs']",
  "    assert qc['assigned_fraction'][s]==assigned/lib['registered_pairs']",
  "   assert sum(table.values())==qc['counts']['0']['sum']",
  "   full.append(bind(folder/'COMPLETE.json'))",
  "   libraries.append({'sample_id':lib['GEO_sample'],'label':lib['label'],'SRA_run':lib['SRA_run'],'input_pairs':lib['registered_pairs'],'QC':bind(folder/'QC.json')})",
  "  assert complete['full_libraries']==full",
  "  assert set(q['path'] for q in complete['small_products'])=={str(p.RUN/f'matrix_01/GSE205906_STAR_fixed9_gene_fragments_s{s}.tsv') for s in (0,1,2)}",
  "  for q in complete['small_products']: assert bind(q['path'])==q",
  "  answer.update(libraries=libraries,pins=list(pins.values()),matrix_complete=bind(p.RUN/'matrix_01/COMPLETE.json'),matrix=bind(primary),full_nine=True,all_s0_cells_equal_FC=True,full_FASTQ_or_BAM_reread=False,index_or_FASTA_fresh_rehash=False)",
  " assert hashlib.sha256(source.read_bytes()).hexdigest()==expected",
  " err=captured_stderr.getvalue().encode('utf-8'); answer.update(stderr_bytes=len(err),stderr_sha256=hashlib.sha256(err).hexdigest(),stderr_text_saved=False)",
  " print(json.dumps(answer,allow_nan=False))",
  "except BaseException as exc:",
  " print(json.dumps({'accepted':False,'exception_class_only':type(exc).__name__})); sys.exit(1)"
), collapse = "\n")

g205run_bridge <- function(options, all = TRUE) {
  before <- g205run_pin(options$producer)
  g205run_need(identical(before$sha256, g205run_constants()$producer_sha), "PRODUCER_SHA_MISMATCH")
  # stderr stays RAM-only. The bridge itself emits safe JSON/classes only.
  text <- suppressWarnings(system2(options$python, c("-I", "-B", "-c", shQuote(g205run_python()),
    shQuote(options$producer), before$sha256, format(options$reserved_growth_bytes, scientific = FALSE),
    format(options$disk_budget_bytes, scientific = FALSE), if (all) "all" else "resource"),
    stdout = TRUE, stderr = TRUE))
  exit <- attr(text, "status"); if (is.null(exit)) exit <- 0L
  text_sha <- digest::digest(paste(text, collapse = "\n"), "sha256", serialize = FALSE)
  g205run_need(exit == 0L && length(text) == 1L && nchar(text, type = "bytes") <= 512 * 1024,
    "SAVED_PRODUCER_BRIDGE_EXIT_OR_SAFE_JSON_FAILED")
  answer <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  g205run_need(identical(answer$accepted, TRUE) && identical(answer$schema, "GSE205906_STAR_SAVED_CUSTODY/v1") &&
    identical(answer$read_only, TRUE) && identical(before, g205run_pin(options$producer)), "PRODUCER_BRIDGE_NOT_ACCEPTED")
  answer$actual_subprocess_exit <- exit
  answer$combined_safe_output_sha256 <- text_sha
  answer$raw_stdout_stderr_persisted <- FALSE
  answer
}

g205run_source_universe <- function(resources) {
  columns <- c(reactome = "term_id", syngo_sensitivity = "term_id", progeny_human = "pathway", dorothea_rat = "tf")
  out <- lapply(names(columns), function(role) {
    tab <- resources[[role]]
    g205run_need(is.data.frame(tab) && columns[[role]] %in% names(tab), "COMPLETE_RAW_RESOURCE_ID_COLUMN_REQUIRED")
    ids <- as.character(tab[[columns[[role]]]])
    g205run_need(!anyNA(ids) && all(nzchar(ids)), "MISSING_SOURCE_RESOURCE_ID")
    unique(ids)
  })
  names(out) <- c("Reactome", "SynGO_sensitivity", "PROGENy", "DoRothEA_AB")
  out
}

# No inference from significance or P/NA. Eligibility comes from the original
# masks and full raw-resource ID universes plus this fit's native coverage.
g205run_domains <- function(packet, result, resources) {
  ids <- rownames(packet$counts); base <- result$native$base_fit
  masks <- list(base$keep, base$keep_filter_by_expr, base$keep_cpm)
  g205run_need(all(vapply(masks, function(z) is.logical(z) && length(z) == length(ids) && !anyNA(z), logical(1))) &&
    identical(unname(base$keep), unname(base$keep_filter_by_expr & base$keep_cpm)), "ORIGINAL_FULL_GENE_MASKS_REQUIRED")
  retained <- ids[base$keep]
  g205run_need(identical(rownames(base$y$counts), retained), "RETAINED_GENE_ORDER_DIFFERS")
  gene_reason <- ifelse(base$keep, "ELIGIBLE_ORIGINAL_FILTER_BY_EXPR_AND_CPM",
    ifelse(!base$keep_filter_by_expr & !base$keep_cpm, "FAILED_FILTER_BY_EXPR_AND_CPM",
      ifelse(!base$keep_filter_by_expr, "FAILED_FILTER_BY_EXPR", "FAILED_CPM_AT_LEAST_THREE")))
  universe <- g205run_source_universe(resources)
  cov <- as.data.frame(result$mapping_coverage)
  g205run_need(all(c("collection", "program_id", "source", "tested", "tested_feature_n") %in% names(cov)),
    "NATIVE_MAPPING_COVERAGE_SCHEMA_REQUIRED")
  slots <- result$slots
  g205run_need(nrow(slots) == 22L && !anyDuplicated(slots[c("table", "collection", "contrast")]), "ALL_22_SLOTS_REQUIRED")
  rows <- lapply(seq_len(nrow(slots)), function(i) {
    slot <- slots[i, ]; tab <- as.data.frame(result$tables[[slot$table]])
    selected <- tab$contrast == slot$contrast
    if (slot$collection == "genes") {
      keys <- ids; eligible <- unname(base$keep); reason <- unname(gene_reason)
      returned <- as.character(tab$gene_id[selected])
      g205run_need(identical(returned, retained), "COMPLETE_NATIVE_RETAINED_GENE_DOMAIN_REQUIRED")
    } else {
      selected <- selected & tab$collection == slot$collection
      keys <- universe[[slot$collection]]
      z <- cov[which(cov$collection == slot$collection), , drop = FALSE]
      zkeys <- as.character(if (slot$collection %in% c("Reactome", "SynGO_sensitivity")) z$program_id else z$source)
      g205run_need(!anyNA(zkeys) && !anyDuplicated(zkeys) && all(zkeys %in% keys) &&
        is.logical(z$tested) && !anyNA(z$tested) && is.numeric(z$tested_feature_n) &&
        all(is.finite(z$tested_feature_n)), "NATIVE_COVERAGE_KEYS_FLAGS_REQUIRED")
      if (slot$collection %in% c("PROGENy", "DoRothEA_AB"))
        g205run_need(setequal(keys, zkeys), "COMPLETE_NATIVE_ACTIVITY_SOURCE_COVERAGE_REQUIRED")
      matched <- match(keys, zkeys)
      eligible <- rep(FALSE, length(keys)); known <- !is.na(matched)
      eligible[known] <- z$tested[matched[known]]
      reason <- rep("NO_STRICT_RAT_MAPPING_IN_FROZEN_RESOURCE", length(keys))
      reason[known] <- ifelse(eligible[known], "ELIGIBLE_ORIGINAL_MEMBER_COVERAGE",
        ifelse(z$tested_feature_n[matched[known]] < 10, "BELOW_ORIGINAL_MINIMUM_10_MEMBERS", "ABOVE_ORIGINAL_MAXIMUM_500_MEMBERS"))
      returned <- as.character(tab$program_id[selected])
    }
    g205run_need(!anyNA(returned) && !anyDuplicated(returned) && all(returned %in% keys[eligible]),
      "RETURNED_KEY_OUTSIDE_OWN_ELIGIBLE_SOURCE_DOMAIN")
    absent <- eligible & !(keys %in% returned)
    g205run_need(slot$table == "unsigned_fgsea" || !any(absent), "NON_FGSEA_ELIGIBLE_OUTPUT_MISSING")
    reason[absent] <- "NOT_RETURNED_BY_UNCHANGED_NATIVE_FGSEA_AFTER_FINITE_T_RANK_AND_SIZE_RULES"
    data.frame(table = rep(slot$table, length(keys)), collection = rep(slot$collection, length(keys)),
      contrast = rep(slot$contrast, length(keys)), feature_id = keys, eligible = eligible,
      reason = reason, native_not_returned = absent, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  g205run_need(!anyDuplicated(out[c("table", "collection", "contrast", "feature_id")]), "DUPLICATE_FULL_SOURCE_DOMAIN_KEY")
  out
}

g205run_write_raw <- function(raw, path) {
  g205run_need(is.raw(raw) && !file.exists(path) && !file.exists(paste0(path, ".partial")), "NO_OUTPUT_OVERWRITE")
  partial <- paste0(path, ".partial")
  con <- file(partial, "wb"); tryCatch(writeBin(raw, con), finally = close(con))
  g205run_need(file.rename(partial, path), "OWNED_OUTPUT_RENAME_FAILED")
  g205run_pin(path)
}

g205run_write_json <- function(value, path) {
  g205run_write_raw(charToRaw(paste0(as.character(jsonlite::toJSON(value,
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")), "\n")), path)
}

g205run_tsv_view <- function(table) {
  # Transport representation only. The frozen core's typed-empty leadingEdge
  # is an AsIs-wrapped plain list; the original writer rejects object lists.
  # Remove only this exact marker in a deep copy. No cell, NULL, empty vector,
  # NA, numeric value, native object or checkpoint is altered. All other
  # classes/objects/dimensions are left for the frozen writer's own checks.
  g205run_need(is.data.frame(table) && !anyDuplicated(names(table)), "TSV_TABLE_VIEW_REQUIRED")
  view <- data.table::copy(table)
  for (column in names(view)) {
    value <- view[[column]]
    if (identical(typeof(value), "list") && identical(class(value), "AsIs") &&
        is.null(dim(value)) && !base::isS4(value)) {
      class(value) <- NULL
      view[[column]] <- value
    }
  }
  view
}

g205run_save_checkpoint <- function(packet, result, custody, output) {
  # Preserve the completed fit BEFORE any derived-domain or post-fit check.
  # Later failures may prevent acceptance/export but must not discard native
  # computation or trigger an automatic second fit. Domains are separate TSV.
  checkpoint <- list(schema = "GSE205906_MATCHED_SAVED_NATIVE_CHECKPOINT/v1",
    prepared_packet = packet, result = result, source_domains = NULL,
    source_domains_deferred_to_separate_export = TRUE, custody = custody,
    statistical_fit_executed = !is.null(result), independent_numeric_acceptance = FALSE)
  path <- file.path(output, "SAVED_NATIVE_CHECKPOINT.rds")
  g205run_need(!file.exists(path) && !file.exists(paste0(path, ".partial")), "NO_CHECKPOINT_OVERWRITE")
  saveRDS(checkpoint, paste0(path, ".partial"), compress = "gzip", version = 3L)
  g205run_need(file.rename(paste0(path, ".partial"), path), "CHECKPOINT_RENAME_FAILED")
  g205run_pin(path)
}


# Isolated evidence only: this file is never an accepted checkpoint or COMPLETE.
g205run_candidate_callback <- function(output, custody) {
  state <- new.env(parent = emptyenv())
  state$invocations <- 0L; state$pin <- NULL
  callback <- function(candidate, original_packet, checks) {
    g205run_need(identical(state$invocations, 0L), "CANDIDATE_CALLBACK_ONCE_ONLY")
    state$invocations <- state$invocations + 1L
    required <- c("original_packet_unchanged", "original_engine_unchanged",
      "working_packet_unchanged", "working_engine_unchanged")
    g205run_need(is.list(checks) && identical(names(checks), required) &&
      all(vapply(checks, function(x) is.logical(x) && length(x) == 1L && !is.na(x), logical(1))),
      "FOUR_CUSTODY_FLAGS_REQUIRED")
    record <- list(schema = "GSE205906_QUARANTINED_NATIVE_CANDIDATE/v1",
      status = "QUARANTINE_ONLY", scientific_acceptance = FALSE,
      independent_numeric_acceptance = FALSE, formal_checkpoint = FALSE,
      candidate = candidate, original_packet = original_packet,
      integrity_checks = checks, field_diagnostics = candidate$integrity_field_diagnostics,
      custody = custody, callback_invocations = state$invocations,
      raw_condition_text_persisted = FALSE, automatic_refit = FALSE)
    path <- file.path(output, "QUARANTINED_NATIVE_CANDIDATE.rds")
    g205run_need(!file.exists(path) && !file.exists(paste0(path, ".partial")),
      "NO_QUARANTINE_OVERWRITE")
    saveRDS(record, paste0(path, ".partial"), compress = "gzip", version = 3L)
    g205run_need(file.rename(paste0(path, ".partial"), path), "QUARANTINE_RENAME_FAILED")
    state$pin <- g205run_pin(path)
    invisible(NULL)
  }
  list(callback = callback, state = state)
}

g205run_sources <- function(project, core) {
  paths <- c(model = "code/R/analyze_gse205906_bulk_model_v2_r4.R",
    helpers = "code/R/gse205906_bulk_model_v2_r4_lib.R",
    config = "config/GSE205906_Ensembl101_full_model_v1_20261006.json",
    loader = "code/R/gse205906_ensembl101_input_loader_v2_20261006.R",
    prefit = "code/R/gse205906_ensembl101_effect_independent_prefit_v1_20261006.R",
    annotation = "code/R/gse205906_ensembl101_annotation_adapter_v1_20261006.R")
  held <- lapply(names(paths), function(role) g205run_read(file.path(project, paths[[role]]), core$g205m_pins()[[role]]$sha256, 256 * 1024))
  names(held) <- names(paths)
  held
}

# Error text and call arguments are inspected in RAM only. Persist only a
# digest, class and code-function names; keep.source=FALSE cannot supply lines.
g205run_error_metadata <- function(condition, calls = sys.calls()) {
  function_name <- function(value) {
    if (is.symbol(value)) {
      name <- as.character(value)
      if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", name)) return(name)
    }
    if (is.call(value) && length(value) == 3L &&
        is.symbol(value[[1L]]) && as.character(value[[1L]]) %in% c("$", "::", ":::")) {
      left <- function_name(value[[2L]]); right <- function_name(value[[3L]])
      if (!is.null(left) && !is.null(right)) return(paste0(left, as.character(value[[1L]]), right))
    }
    NULL
  }
  head_name <- function(value) if (is.call(value) && length(value)) function_name(value[[1L]]) else NULL
  stack <- Filter(Negate(is.null), lapply(utils::tail(calls, 100L), head_name))
  list(exception_class_only = class(condition)[[1L]],
    condition_sha256 = digest::digest(conditionMessage(condition), "sha256", serialize = FALSE),
    code_location = list(condition_function = head_name(conditionCall(condition)),
      stack_function_names_no_arguments = unname(stack),
      source_line = NULL, source_line_status = "UNAVAILABLE_FUNCTION_ONLY_PARSE_KEEP_SOURCE_FALSE"),
    raw_condition_text_persisted = FALSE, call_arguments_persisted = FALSE)
}

g205run_execute <- function(options, core) {
  constants <- g205run_constants(); phase <- "SAVED_PRODUCER_CUSTODY"
  owned <- FALSE; started <- proc.time()[["elapsed"]]; conditions <- list()
  quarantine <- NULL
  checkpoint_pin <- NULL; core_fit_invocations <- 0L; error_metadata <- NULL
  diagnose_error <- function(condition) {
    error_metadata <<- tryCatch(g205run_error_metadata(condition, sys.calls()),
      error = function(ignored) list(exception_class_only = class(condition)[[1L]],
        diagnostic_collection_failure_class_only = class(ignored)[[1L]],
        raw_condition_text_persisted = FALSE, call_arguments_persisted = FALSE))
  }
  # Condition content is used in RAM for a fingerprint only, never printed or
  # written. This also prevents caller log capture from persisting raw messages.
  notice <- function(condition) {
    key <- paste(class(condition)[[1L]], digest::digest(conditionMessage(condition), "sha256", serialize = FALSE), sep = ":")
    conditions[[key]] <<- if (is.null(conditions[[key]])) 1L else conditions[[key]] + 1L
    restart <- if (inherits(condition, "warning")) "muffleWarning" else "muffleMessage"
    invokeRestart(restart)
  }
  tryCatch(withCallingHandlers({
    custody <- g205run_bridge(options)
    phase <- "FROZEN_SOURCE_AND_INPUT_BYTES"
    source <- g205run_sources(options$project, core)
    extras <- list(
      parser = g205run_read(file.path(options$project, "code/R/gse205906_ensembl101_nofit_RAM_worker_v2_20261006.R"), constants$parser_sha),
      adapter = g205run_read(file.path(options$project, "code/R/gse205906_resource_object_validation_v1.R"), constants$adapter_sha),
      writer = g205run_read(file.path(options$project, "code/R/gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R"), constants$writer_sha))
    pure <- new.env(parent = baseenv())
    core$g205m_import(extras$parser$raw, c("gse205906_ens101_nofit_need_v2", "gse205906_ens101_nofit_overlay_v2", "gse205906_ens101_nofit_metadata_consistency_v2"), pure)
    core$g205m_import(extras$adapter$raw, "gse205906_resource_label_adapter_v1", pure)
    core$g205m_import(extras$writer$raw, c("gse205906_ens101_worker_need_v1", "gse205906_ens101_worker_tsv_v1"), pure)
    fixed_sha <- c(metadata = "b0df6fa1354a83237bcfcd8bf8c323517deecd0e5680da8202489c4464936444",
      fixed_gene_manifest = "751e43b3aba8c2245981dc0041fec273c455c6390e85ef7ec4514e56118bf3b0",
      GTF = "84fcb6692d62e4d8c75617461ae5f471fe89e10b6982c82f349393b1a80045df",
      gene_resources = "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51",
      orthology = "e9dd0c9268e0044f615085e27292cbf977f5a648003ba90054cf3d01225cf998",
      hisat2_completion = "25c11e9824b72a45d6dc0e6625439ab73e2e19f8b0eeb524c41955c06b1d697e",
      hisat2_manifest = "428c0f9827b068673b7ecddad2e7fef3cecde09d10c3c35927fc8eb99d9743e7",
      hisat2_matrix = "e98296fdf308dbaa82325a3ca6b0b041c6f2569e98be75a45757141b0a1f77a2")
    held <- lapply(names(fixed_sha), function(role) g205run_read(constants[[role]], unname(fixed_sha[[role]])))
    names(held) <- names(fixed_sha)
    matrix <- g205run_read(constants$matrix, custody$matrix$sha256, 32 * 1024^2)
    completion <- g205run_read(paste0(dirname(constants$matrix), "/COMPLETE.json"), custody$matrix_complete$sha256, 1024^2)
    reference <- g205run_read(constants$reference_admission)
    phase <- "ORIGINAL_METADATA_AND_RESOURCE_REPRESENTATION"
    ledger <- pure$gse205906_ens101_nofit_overlay_v2(held$metadata$raw)
    manifest <- g205run_json(held$fixed_gene_manifest$raw)
    pure$gse205906_ens101_nofit_metadata_consistency_v2(manifest, ledger)
    gene_ids <- unlist(manifest$gene_ids, use.names = FALSE)
    engine <- core$g205m_engine(lapply(source, `[[`, "raw"))
    attributes <- engine$gse205906_ens101_gtf_attributes_v2(held$GTF$raw)
    # readRDS understands the frozen RDS's compression; raw unserialize() would
    # incorrectly assume an uncompressed serialization stream.
    resources <- pure$gse205906_resource_label_adapter_v1(readRDS(constants$gene_resources))
    con <- textConnection(rawToChar(held$orthology$raw), local = TRUE)
    orthology <- tryCatch(utils::read.delim(con, sep = "\t", quote = '"', colClasses = "character",
      na.strings = "NA", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE), finally = close(con))
    objects <- list(sample_ledger = ledger, fixed_gene_ids = gene_ids, gene_attributes = attributes,
      gene_resources = resources, orthology_table = orthology)
    phase <- "FRESH_OWNED_OUTPUT"
    g205run_need(!file.exists(options$output) && dir.create(options$output, recursive = FALSE), "FRESH_OUTPUT_CREATION_FAILED")
    owned <- TRUE
    runtime_record <- list(schema = "GSE205906_MATCHED_STATISTICAL_RUNTIME/v1", fit_started = FALSE,
      R = as.character(getRversion()), platform = R.version$platform,
      sessionInfo = capture.output(utils::sessionInfo()),
      environment = as.list(Sys.getenv(c("R_LIBS_USER", "R_LIBS_SITE", "LC_ALL", "LANG", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "R_DATATABLE_NUM_THREADS"))),
      command = commandArgs(), measured_producer_resource_contract = custody$resource,
      quiet_window_is_caller_confirmation_not_process_scheduler = TRUE)
    runtime_pin <- g205run_write_json(runtime_record, file.path(options$output, "STATISTICAL_RUNTIME.json"))
    roles <- c("metadata", "fixed_gene_manifest", "GTF", "gene_resources", "orthology", "reference_admission", "statistical_runtime", "hisat2_completion", "hisat2_manifest", "hisat2_matrix")
    pins <- lapply(held, `[[`, "pin"); pins$reference_admission <- reference$pin; pins$statistical_runtime <- runtime_pin
    binding <- list(schema = "gse205906-matched-downstream-held-inputs-v1", matrix_path = constants$matrix,
      matrix_complete_bytes = completion$raw, matrix_complete_file_pin = completion$pin,
      producer_source_sha256 = constants$producer_sha, held_files = pins[roles],
      held_object_hashes = lapply(objects, core$g205m_hash), hisat2_comparator_root = constants$hisat)
    phase <- "FROZEN_PREPARE_NO_FIT"
    prepared <- core$g205m_prepare(lapply(source, `[[`, "raw"), matrix$raw, ledger, gene_ids, attributes, resources, orthology, binding)
    result <- NULL; domains <- NULL
    if (options$mode == "execute") {
      phase <- "FRESH_MEASURED_RESOURCE_AND_SINGLE_NATIVE_FIT"
      resource_at_fit <- g205run_bridge(options, all = FALSE)
      phase <- "CORE_FIT_AND_IMMEDIATE_IMMUTABLE_CHECKPOINT"
      core_fit_invocations <- core_fit_invocations + 1L
      quarantine <- g205run_candidate_callback(options$output, custody)
      result <- core$g205m_fit(prepared, on_candidate = quarantine$callback)
      checkpoint_pin <- g205run_save_checkpoint(prepared$packet, result, custody, options$output)
      phase <- "COMPLETE_SOURCE_DOMAIN_LEDGER"
      domains <- g205run_domains(prepared$packet, result, resources)
    } else {
      resource_at_fit <- NULL
      checkpoint_pin <- g205run_save_checkpoint(prepared$packet, NULL, custody, options$output)
    }
    phase <- "INPUTS_UNCHANGED_AND_SAVED_EXPORT"
    tracked <- c(list(options$caller_source_pin, options$core_source_pin),
      lapply(source, `[[`, "pin"), lapply(extras, `[[`, "pin"), lapply(held, `[[`, "pin"),
      list(matrix$pin, completion$pin, reference$pin))
    g205run_need(all(vapply(tracked, function(pin) identical(pin, g205run_pin(pin$path)), logical(1))), "HELD_SOURCE_OR_INPUT_CHANGED")
    after_custody <- g205run_bridge(options)
    g205run_need(identical(custody$pins, after_custody$pins), "PRODUCER_PRODUCTS_CHANGED_DURING_CONSUMPTION")
    products <- list(checkpoint_pin, runtime_pin)
    if (!is.null(result)) {
      for (role in names(result$tables)) products[[length(products) + 1L]] <- g205run_write_raw(
        pure$gse205906_ens101_worker_tsv_v1(g205run_tsv_view(result$tables[[role]])), file.path(options$output, paste0(role, ".tsv")))
      for (role in c("source_domains", "slots", "mapping_coverage", "native_not_returned_programs")) {
        tab <- if (role == "source_domains") domains else result[[role]]
        products[[length(products) + 1L]] <- g205run_write_raw(pure$gse205906_ens101_worker_tsv_v1(g205run_tsv_view(tab)), file.path(options$output, paste0(role, ".tsv")))
      }
    }
    g205run_need(all(vapply(tracked, function(pin) identical(pin, g205run_pin(pin$path)), logical(1))), "INPUT_OR_SOURCE_CHANGED_DURING_EXPORT")
    receipt <- list(status = if (is.null(result)) "ADMITTED_NO_FIT" else "COMPLETE_SAVED_NATIVE_NOT_INDEPENDENTLY_ACCEPTED",
      schema = "GSE205906_MATCHED_DOWNSTREAM_CALLER/v1", phase = "SAVED_EXPORT_COMPLETE", mode = options$mode,
      quarantine_pin = if (is.null(quarantine)) NULL else quarantine$state$pin,
      candidate_callback_invocations = if (is.null(quarantine)) 0L else quarantine$state$invocations,
      quarantine_scientific_acceptance = FALSE,
      core_fit_invocations = core_fit_invocations, input_pins = tracked, producer_custody = custody,
      resource_at_fit = resource_at_fit, resource_after = after_custody$resource,
      input_sources_and_producer_products_unchanged = TRUE, products = products,
      original_tables_preserved = TRUE, raw_condition_text_persisted = FALSE, condition_fingerprints = conditions,
      TSV_transport = "COPY_ONLY_EXACT_AsIs_PLAIN_LIST_MARKER_REMOVED_BEFORE_FROZEN_WRITER",
      native_checkpoint_transport_modified = FALSE,
      own_complete_source_domain_written = !is.null(domains), source_domain_rows = if (is.null(domains)) NULL else nrow(domains),
      elapsed_seconds = proc.time()[["elapsed"]] - started, parent_must_record_actual_R_process_exit = TRUE,
      saved_comparison = "NOT_RUN_SEPARATE_SAVED_ONLY_STEP", top_k = "PENDING_POLICY_BINDING",
      LOAO_or_label_reference_calls = 0L, old_models_or_scores_reused = FALSE, independent_numeric_acceptance = FALSE,
      candidate_admission = FALSE, biological_validation = FALSE)
    g205run_write_json(receipt, file.path(options$output, if (is.null(result)) "ADMITTED.json" else "COMPLETE.json"))
    list(status = receipt$status, output = options$output, core_fit_invocations = core_fit_invocations,
      checkpoint_pin = checkpoint_pin, actual_R_exit_observed_here = FALSE)
  }, warning = notice, message = notice, error = diagnose_error), error = function(condition) {
    failure <- list(status = "FAILED", phase = phase, exception_class_only = class(condition)[[1L]],
      error_metadata = error_metadata,
      elapsed_seconds = proc.time()[["elapsed"]] - started, raw_condition_text_persisted = FALSE,
      core_fit_invocations = core_fit_invocations, retained_checkpoint_pin = checkpoint_pin,
      retained_quarantine_pin = if (is.null(quarantine)) NULL else quarantine$state$pin,
      candidate_callback_invocations = if (is.null(quarantine)) 0L else quarantine$state$invocations,
      quarantine_scientific_acceptance = FALSE,
      automatic_refit_or_checkpoint_overwrite = FALSE,
      actual_R_exit_observed_here = FALSE, successful_native_or_scientific_acceptance_inferred = FALSE)
    if (owned && !file.exists(file.path(options$output, "FAILED.json")))
      tryCatch(g205run_write_json(failure, file.path(options$output, "FAILED.json")), error = function(ignored) NULL)
    failure
  })
}

g205run_main <- function(args) {
  options <- g205run_options(args)
  if (options$mode == "dry-run") return(list(status = "DRY_RUN_NO_INPUT_READ_OR_FIT", options = options,
    matrix = g205run_constants()$matrix, core_function = "g205m_prepare_then_explicit_g205m_fit_once",
    no_output_created = TRUE, comparison = "SEPARATE_SAVED_ONLY_STEP", top_k = "PENDING_POLICY_BINDING"))
  g205run_need(requireNamespace("digest", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE), "EXISTING_PACKAGES_REQUIRED_NO_INSTALL")
  core_path <- file.path(options$project, "code/R/gse205906_matched_downstream_core_custody_v1_20261022.R")
  held <- g205run_read(core_path, g205run_constants()$core_sha)
  options$core_source_pin <- held$pin
  options$caller_source_pin <- g205run_pin(file.path(options$project, "code/R/run_gse205906_matched_downstream_custody_v1_20261022.R"))
  core <- new.env(parent = baseenv())
  # This particular SHA is the reviewed functions-only module, not an old driver.
  eval(parse(text = rawToChar(held$raw), keep.source = FALSE), core)
  g205run_execute(options, core)
}

if (sys.nframe() == 0L) {
  answer <- tryCatch(g205run_main(commandArgs(trailingOnly = TRUE)),
    error = function(condition) list(status = "FAILED_BEFORE_OWNED_OUTPUT", exception_class_only = class(condition)[[1L]], raw_condition_text_persisted = FALSE))
  if (requireNamespace("jsonlite", quietly = TRUE)) cat(as.character(jsonlite::toJSON(answer, auto_unbox = TRUE, null = "null")), "\n")
  else cat("MATCHED_DOWNSTREAM_SAFE_STATUS", answer$status, "\n")
  quit(status = if (startsWith(answer$status, "FAILED")) 1L else 0L, save = "no")
}



