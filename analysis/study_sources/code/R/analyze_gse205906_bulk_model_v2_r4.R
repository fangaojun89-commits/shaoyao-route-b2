#!/usr/bin/env Rscript

if (!identical(Sys.getenv("GSE205906_V2_R4_EXTERNAL_LAUNCH", unset = ""), "1") ||
    !nzchar(Sys.getenv("GSE205906_V2_R4_LAUNCH_CONTEXT", unset = "")) ||
    !nzchar(Sys.getenv("GSE205906_V2_R4_LAUNCH_CONTEXT_HMAC", unset = "")) ||
    !grepl("^[0-9a-f]{64}$", Sys.getenv("GSE205906_V2_R4_LAUNCHER_TOKEN", unset = ""))) {
  stop("V2_R4_EXTERNAL_LAUNCH_REQUIRED_BEFORE_LIBRARY_CAPTURE", call. = FALSE)
}

v2r4_bootstrap_capture <- function(path, label) {
  canonical <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (nzchar(Sys.readlink(canonical))) stop(label, "_SYMLINK_NOT_PERMITTED", call. = FALSE)
  info <- file.info(canonical)
  if (nrow(info) != 1L || is.na(info$size) || isTRUE(info$isdir)) {
    stop(label, "_NOT_REGULAR_FILE", call. = FALSE)
  }
  con <- file(canonical, open = "rb")
  on.exit(close(con), add = TRUE)
  remaining <- as.numeric(info$size)
  parts <- list()
  while (remaining > 0) {
    part <- readBin(con, what = "raw", n = as.integer(min(remaining, 64 * 1024^2)))
    if (!length(part)) stop(label, "_SHORT_READ", call. = FALSE)
    parts[[length(parts) + 1L]] <- part
    remaining <- remaining - length(part)
  }
  payload <- if (length(parts)) do.call(c, parts) else raw()
  if (as.numeric(length(payload)) != as.numeric(info$size)) {
    stop(label, "_BYTES_CHANGED_DURING_CAPTURE", call. = FALSE)
  }
  structure(list(
    raw = payload, canonical_path = canonical, bytes = as.numeric(length(payload)),
    sha256 = tolower(digest::digest(payload, algo = "sha256", serialize = FALSE)),
    md5 = tolower(digest::digest(payload, algo = "md5", serialize = FALSE)), label = label
  ), class = "v2r4_captured_bytes")
}

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
                             winslash = "/", mustWork = TRUE)
library_path <- normalizePath(file.path(dirname(script_path), "gse205906_bulk_model_v2_r4_lib.R"),
                              winslash = "/", mustWork = TRUE)
script_runtime_buffer <- v2r4_bootstrap_capture(script_path, "MODEL_SOURCE_RUNTIME")
library_runtime_buffer <- v2r4_bootstrap_capture(library_path, "SHARED_LIBRARY_RUNTIME")
library_expressions <- parse(text = enc2utf8(rawToChar(library_runtime_buffer$raw)),
                             keep.source = FALSE)
for (library_expression in library_expressions) {
  eval(library_expression, envir = .GlobalEnv)
}
rm(v2r4_bootstrap_capture, library_expressions, library_expression)

child_request <- v2r4_capture_child_request_args(
  commandArgs(trailingOnly = TRUE), "MODEL"
)
cli <- child_request$cli
project_root <- v2r4_project_root(script_path)
default_config <- file.path(project_root, "config", "GSE205906_bulk_model_v2_r4.json")
config_path <- cli$config

required_control_args <- c(
  "authority_manifest", "authority_signature", "authority_public_key", "snapshot_root",
  "lease_manifest", "lease_signature", "lease_public_key", "lease_request",
  "lease_identity_request", "peer_lease_identity_request",
  "state_root", "run_id", "holder_principal", "output_dir", "final_ack",
  "finalization_public_key", "validator_request", "runtime_post_receipt",
  "science_complete", "control_complete", "finalization_request",
  "final_challenge_request", "final_challenge_response",
  "model_consumer_receipt", "validator_consumer_receipt"
)
missing_control <- required_control_args[vapply(required_control_args, function(k) {
  is.null(cli[[k]]) || !nzchar(cli[[k]])
}, logical(1))]
if (length(missing_control)) {
  v2r4_stop("EXTERNAL_AUTHORITY_REQUIRED_BEFORE_ANY_DATA_READ:%s",
            paste(missing_control, collapse = ","))
}
for (k in c("authority_manifest", "authority_signature", "authority_public_key")) {
  if (!file.exists(cli[[k]])) {
    v2r4_stop("EXTERNAL_AUTHORITY_REQUIRED_BEFORE_ANY_DATA_READ:%s_MISSING", k)
  }
}

v2r4_require_packages(c(
  "data.table", "digest", "jsonlite", "openssl", "yaml", "edgeR", "limma",
  "DESeq2", "fgsea", "decoupleR"
))
config_path <- v2r4_norm_existing(config_path)
v2r4_assert(identical(config_path, v2r4_norm_existing(default_config)),
            "CONFIG_CANONICAL_PATH_NOT_LOCKED")
config_buffer <- v2r4_capture_file(config_path, "ANALYSIS_CONFIG_RUNTIME", project_root)
cfg <- v2r4_load_config_buffer(config_buffer)
analysis_role <- cfg$gate$analysis_role
launch_context <- v2r4_require_external_launch_context(cfg, "MODEL")
v2r4_validate_child_request_bindings(child_request, cfg, launch_context, "MODEL")

v2r4_assert(cli$execution_mode %in% c("SYNTHETIC", "PRODUCTION") &&
              identical(cli$execution_mode, launch_context$mode),
            "SIGNED_EXECUTION_MODE_INVALID")
required_execution <- c("attempt_root", "source_key_sha256", "lease_key_sha256")
missing_execution <- required_execution[vapply(required_execution, function(k) {
  is.null(cli[[k]]) || !nzchar(cli[[k]])
}, logical(1))]
v2r4_assert(!length(missing_execution), "EXECUTION_CONTROL_ARGUMENTS_MISSING:%s",
            paste(missing_execution, collapse = ","))
attempt_root <- v2r4_norm_existing(cli$attempt_root)
v2r4_require_outside(attempt_root, project_root, "EXTERNAL_ATTEMPT_ROOT")
v2r4_reject_symlink(attempt_root, "EXTERNAL_ATTEMPT_ROOT")
snapshot_root <- v2r4_norm_existing(cli$snapshot_root)
state_root <- v2r4_norm_existing(cli$state_root)
output_dir <- v2r4_norm_target(cli$output_dir)
lease_request_path <- v2r4_norm_target(cli$lease_request)
lease_identity_request_path <- v2r4_norm_target(cli$lease_identity_request)
peer_lease_identity_request_path <- v2r4_norm_target(cli$peer_lease_identity_request)
final_ack_path <- v2r4_norm_target(cli$final_ack)
finalization_public_key_path <- v2r4_norm_existing(cli$finalization_public_key)
validator_request_path <- v2r4_norm_existing(cli$validator_request)
runtime_post_receipt_path <- v2r4_norm_target(cli$runtime_post_receipt)
science_complete_path <- v2r4_norm_target(cli$science_complete)
control_complete_path <- v2r4_norm_target(cli$control_complete)
finalization_request_path <- v2r4_norm_target(cli$finalization_request)
model_consumer_receipt_path <- v2r4_norm_target(cli$model_consumer_receipt)
validator_consumer_receipt_path <- v2r4_norm_target(cli$validator_consumer_receipt)
external_authority_paths <- list(
  c(v2r4_norm_existing(cli$authority_manifest), "SOURCE_AUTHORITY"),
  c(v2r4_norm_existing(cli$authority_signature), "SOURCE_SIGNATURE"),
  c(v2r4_norm_existing(cli$authority_public_key), "SOURCE_PUBLIC_KEY"),
  c(snapshot_root, "SNAPSHOT_ROOT"),
  c(v2r4_norm_existing(cli$lease_public_key), "LEASE_PUBLIC_KEY")
)
for (pair in external_authority_paths) {
  v2r4_require_outside(pair[[1L]], project_root, pair[[2L]])
  v2r4_reject_symlink(pair[[1L]], pair[[2L]])
}
for (pair in list(
  c(state_root, "STATE_ROOT"),
  c(output_dir, "OUTPUT_ROOT"), c(lease_request_path, "LEASE_REQUEST"),
  c(lease_identity_request_path, "MODEL_LEASE_IDENTITY_REQUEST"),
  c(peer_lease_identity_request_path, "VALIDATOR_LEASE_IDENTITY_REQUEST"),
  c(v2r4_norm_target(cli$lease_manifest), "LEASE_MANIFEST"),
  c(v2r4_norm_target(cli$lease_signature), "LEASE_SIGNATURE"),
  c(final_ack_path, "FINAL_ACK"), c(finalization_public_key_path, "FINALIZATION_PUBLIC_KEY"),
  c(validator_request_path, "VALIDATOR_REQUEST"),
  c(runtime_post_receipt_path, "RUNTIME_POST_RECEIPT"),
  c(science_complete_path, "SCIENCE_COMPLETE"), c(control_complete_path, "CONTROL_COMPLETE"),
  c(finalization_request_path, "FINALIZATION_REQUEST"),
  c(model_consumer_receipt_path, "MODEL_CONSUMER_RECEIPT"),
  c(validator_consumer_receipt_path, "VALIDATOR_CONSUMER_RECEIPT"),
  c(v2r4_norm_existing(launch_context$guard_ready_receipt_path), "GUARD_READY_RECEIPT"),
  c(v2r4_norm_target(launch_context$guard_stop_request_path), "GUARD_STOP_REQUEST"),
  c(v2r4_norm_target(launch_context$guard_stopped_receipt_path), "GUARD_STOPPED_RECEIPT"),
  c(v2r4_norm_target(launch_context$runtime_final_receipt_path), "RUNTIME_FINAL_RECEIPT")
)) v2r4_require_within(pair[[1L]], attempt_root, pair[[2L]])
v2r4_assert(
  identical(lease_identity_request_path,
            v2r4_norm_target(launch_context$model_lease_identity_request_path)) &&
    identical(peer_lease_identity_request_path,
              v2r4_norm_target(launch_context$validator_lease_identity_request_path)) &&
    identical(lease_request_path,
              v2r4_norm_target(launch_context$combined_lease_request_path)),
  "MODEL_LEASE_REQUEST_PATH_CONTEXT_BINDING_INVALID"
)
v2r4_assert(!file.exists(output_dir), "OUTPUT_ALREADY_EXISTS")
source_key_sha256 <- tolower(cli$source_key_sha256)
lease_key_sha256 <- tolower(cli$lease_key_sha256)
expected_source_decision <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  cfg$authority$synthetic_source_decision
} else {
  cfg$authority$production_source_decision
}
if (identical(cli$execution_mode, "PRODUCTION")) {
  v2r4_assert(identical(source_key_sha256,
                        cfg$authority$production_source_public_key_sha256) &&
                identical(lease_key_sha256, cfg$lease$production_public_key_sha256),
              "PRODUCTION_SOURCE_OR_LEASE_KEY_PIN_INVALID")
}
v2r4_assert(!file.exists(lease_request_path), "LEASE_REQUEST_ALREADY_EXISTS")
v2r4_assert(!file.exists(lease_identity_request_path),
            "MODEL_LEASE_IDENTITY_REQUEST_ALREADY_EXISTS")
v2r4_assert(!file.exists(final_ack_path), "FINAL_ACK_ALREADY_EXISTS")
v2r4_assert(identical(source_key_sha256, launch_context$source_snapshot_public_key_sha256) &&
              identical(lease_key_sha256, launch_context$lease_public_key_sha256) &&
              identical(v2r4_capture_file(cli$authority_manifest,
                                          "SOURCE_MANIFEST_LAUNCH_BINDING")$sha256,
                        launch_context$source_snapshot_manifest_sha256) &&
              identical(finalization_public_key_path,
                        v2r4_norm_existing(launch_context$finalization_public_key_captured_path)) &&
              {
                captured_finalization_key <- v2r4_capture_file(
                  cli$finalization_public_key, "FINALIZATION_KEY_LAUNCH_BINDING", attempt_root
                )
                identical(as.numeric(captured_finalization_key$bytes),
                          as.numeric(launch_context$finalization_public_key_captured_bytes)) &&
                  identical(captured_finalization_key$sha256,
                            launch_context$finalization_public_key_captured_sha256) &&
                  identical(captured_finalization_key$sha256,
                            launch_context$finalization_public_key_sha256)
              },
            "EXTERNAL_LAUNCH_SOURCE_LEASE_FINALIZATION_BINDING_INVALID")

source_bundle <- v2r4_verify_detached_bundle(
  cli$authority_manifest, cli$authority_signature, cli$authority_public_key,
  source_key_sha256, "SOURCE_AUTHORITY"
)
auth <- v2r4_parse_authority_buffer(source_bundle$manifest, cfg, expected_source_decision)
expected_roles <- c(unlist(cfg$authority$required_dynamic_roles, use.names = FALSE),
                    unlist(cfg$authority$required_data_roles, use.names = FALSE))
v2r4_assert(setequal(auth$role, expected_roles), "SOURCE_AUTHORITY_ROLE_SET_NOT_EXACT")

guard_path <- v2r4_norm_existing(file.path(
  project_root, cfg$authority$dynamic_runtime_paths$os_guard_source
))
runtime_buffers <- list(
  model_source = script_runtime_buffer,
  validator_source = v2r4_capture_file(file.path(
    project_root, cfg$authority$dynamic_runtime_paths$validator_source
  ), "VALIDATOR_SOURCE_RUNTIME", project_root),
  shared_library_source = library_runtime_buffer,
  analysis_config = config_buffer,
  os_guard_source = v2r4_capture_file(guard_path, "OS_GUARD_RUNTIME", project_root)
)
runtime_bindings <- v2r4_bind_runtime_roles(auth, cfg, project_root, runtime_buffers)
v2r4_assert(identical(script_runtime_buffer$sha256, launch_context$model_source_sha256) &&
              identical(library_runtime_buffer$sha256, launch_context$library_source_sha256) &&
              identical(config_buffer$sha256, launch_context$config_sha256) &&
              identical(runtime_buffers$validator_source$sha256,
                        launch_context$validator_source_sha256) &&
              identical(runtime_buffers$os_guard_source$sha256,
                        launch_context$os_guard_sha256),
            "MODEL_RUNTIME_BYTES_NOT_BOUND_TO_EXTERNAL_LAUNCH")

run_id <- as.character(cli$run_id)
holder_principal <- as.character(cli$holder_principal)
v2r4_assert(grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", run_id), "RUN_ID_INVALID")
v2r4_assert(grepl("^[A-Za-z0-9][A-Za-z0-9._@-]{0,255}$", holder_principal),
            "HOLDER_PRINCIPAL_INVALID")
model_identity <- v2r4_process_identity(runtime_buffers$os_guard_source, Sys.getpid())
v2r4_assert_child_stop_key_absent()
guard_state <- v2r4_attach_existing_guard(
  runtime_buffers$os_guard_source, state_root, run_id, cfg,
  expected_owner_pid = as.integer(launch_context$launcher_pid),
  expected_owner_start_utc = launch_context$launcher_start_utc,
  expected_stop_hmac_key_sha256 = launch_context$guard_stop_hmac_key_sha256,
  expected_finalization_public_key_sha256 = launch_context$finalization_public_key_sha256,
  expected_finalization_public_key_captured_path =
    launch_context$finalization_public_key_captured_path,
  expected_finalization_public_key_captured_bytes =
    launch_context$finalization_public_key_captured_bytes,
  expected_finalization_public_key_captured_sha256 =
    launch_context$finalization_public_key_captured_sha256,
  expected_ready_path = launch_context$guard_ready_receipt_path,
  expected_stop_path = launch_context$guard_stop_request_path,
  expected_stopped_path = launch_context$guard_stopped_receipt_path
)
guard_state$initial_challenge <- v2r4_guard_challenge(
  guard_state, runtime_buffers$os_guard_source, cfg, Sys.getpid()
)
authorized_snapshot_tree_digest <- v2r4_authorized_snapshot_tree_digest(auth)
runtime_role_digest <- v2r4_runtime_role_digest(runtime_bindings)
model_identity_request_buffer <- v2r4_write_lease_identity_request(
  lease_identity_request_path, cfg, "MODEL", run_id, holder_principal,
  authorized_snapshot_tree_digest, guard_state, runtime_role_digest,
  launch_context, model_identity, child_request, attempt_root
)
lease_wait <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  cfg$execution_contract$lease_wait_seconds_synthetic
} else {
  cfg$execution_contract$lease_wait_seconds_production
}
v2r4_wait_for_files(c(peer_lease_identity_request_path, cli$lease_manifest,
                      cli$lease_signature, cli$lease_public_key),
                    lease_wait, "EXTERNAL_SIGNED_LEASE")
validator_request_buffer_prelease <- v2r4_capture_file(
  validator_request_path, "MODEL_PRELEASE_VALIDATOR_CHILD_REQUEST", attempt_root
)
validator_request_bundle_prelease <- v2r4_parse_child_request_buffer(
  validator_request_buffer_prelease, validator_request_buffer_prelease$sha256, "VALIDATOR"
)
v2r4_validate_child_request_bindings(
  validator_request_bundle_prelease, cfg, launch_context, "VALIDATOR"
)
validator_identity_request_buffer <- v2r4_capture_file(
  peer_lease_identity_request_path, "MODEL_CAPTURED_VALIDATOR_LEASE_IDENTITY_REQUEST",
  attempt_root
)
validator_identity_request <- v2r4_validate_lease_identity_request(
  validator_identity_request_buffer, cfg, "VALIDATOR",
  list(
    run_id = run_id,
    holder_principal = holder_principal,
    child_request_sha256 = validator_request_buffer_prelease$sha256,
    child_request_nonce = validator_request_bundle_prelease$request$request_nonce,
    peer_request_nonce = child_request$request$request_nonce,
    snapshot_tree_sha256 = authorized_snapshot_tree_digest,
    runtime_role_digest = runtime_role_digest,
    source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
    external_launch_context_sha256 = launch_context$context_sha256,
    runtime_pre_digest = launch_context$runtime_pre_digest,
    guard_owner_role = guard_state$owner_role,
    guard_owner_pid = guard_state$owner_pid,
    guard_owner_start_utc = guard_state$owner_start_utc,
    guard_pid = guard_state$guard_pid,
    guard_start_utc = guard_state$guard_start_utc,
    guard_nonce = guard_state$guard_nonce,
    guard_public_key_sha256 = guard_state$guard_public_key_sha256,
    guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
    guard_mutex_file_id = guard_state$mutex_file_id,
    guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
    finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
    finalization_public_key_captured_path =
      guard_state$finalization_public_key_captured_path,
    finalization_public_key_captured_bytes =
      guard_state$finalization_public_key_captured_bytes,
    finalization_public_key_captured_sha256 =
      guard_state$finalization_public_key_captured_sha256
  ), runtime_buffers$os_guard_source
)
lease_bundle <- v2r4_verify_detached_bundle(
  cli$lease_manifest, cli$lease_signature, cli$lease_public_key,
  lease_key_sha256, "LIVE_LEASE"
)
lease <- v2r4_validate_lease_buffer(
  lease_bundle$manifest, cfg, authorized_snapshot_tree_digest, run_id, holder_principal,
  runtime_buffers$os_guard_source, guard_state, launch_context, runtime_role_digest,
  state_root,
  expected_model_pid = Sys.getpid(),
  expected_model_start_utc = model_identity$start_utc,
  expected_model_request_sha256 = child_request$buffer$sha256,
  expected_model_request_nonce = child_request$request$request_nonce,
  expected_model_identity_request_sha256 = model_identity_request_buffer$sha256,
  expected_validator_pid = validator_identity_request$process_pid,
  expected_validator_start_utc = validator_identity_request$process_start_utc,
  expected_validator_request_sha256 = validator_request_buffer_prelease$sha256,
  expected_validator_request_nonce = validator_request_bundle_prelease$request$request_nonce,
  expected_validator_identity_request_sha256 = validator_identity_request_buffer$sha256
)

# No study object is opened before both live child identities are externally signed.
snapshot <- v2r4_capture_snapshot(auth, snapshot_root, defer_roles = "assembled_counts_rds")
producer_source_provenance <- v2r4_validate_producer_source_buffers(snapshot, cfg)
ledger <- v2r4_validate_ledger(v2r4_read_tsv_raw(v2r4_object(snapshot, "sample_ledger"), "SAMPLE_LEDGER"), cfg)
v2r4_validate_upstream_gate(
  v2r4_read_tsv_raw(v2r4_object(snapshot, "upstream_e1_gate"), "UPSTREAM_E1_GATE"), cfg
)
reference_manifest <- v2r4_read_tsv_raw(v2r4_object(snapshot, "reference_manifest"), "REFERENCE_MANIFEST")
v2r4_validate_reference_manifest(reference_manifest, cfg)
signed_fasta <- v2r4_role_row(auth, "reference_fasta_identity")
signed_gtf <- v2r4_role_row(auth, "reference_gtf_identity")
manifest_fasta <- reference_manifest[base::which(reference_manifest[["role"]] == "reference_fasta")]
manifest_gtf <- reference_manifest[base::which(reference_manifest[["role"]] == "gene_annotation_gtf")]
v2r4_assert(signed_fasta$materialization[[1L]] == "IDENTITY_ONLY" &&
              signed_gtf$materialization[[1L]] == "IDENTITY_ONLY" &&
              signed_fasta$sha256[[1L]] == cfg$input_snapshot$reference_fasta_sha256 &&
              signed_gtf$sha256[[1L]] == cfg$input_snapshot$gtf_sha256 &&
              as.numeric(signed_fasta$bytes[[1L]]) == as.numeric(manifest_fasta$bytes[[1L]]) &&
              as.numeric(signed_gtf$bytes[[1L]]) == as.numeric(manifest_gtf$bytes[[1L]]),
            "SIGNED_REFERENCE_IDENTITY_ROWS_INVALID")
v2r4_validate_strand(
  v2r4_read_tsv_raw(v2r4_object(snapshot, "strand_decision"), "STRAND_DECISION"),
  v2r4_read_tsv_raw(v2r4_object(snapshot, "strand_assignment_by_sample"), "STRAND_ASSIGNMENT"), cfg
)
cohort_species_buffer <- v2r4_object(snapshot, "cohort_species_gate")
sample_species_buffer <- v2r4_object(snapshot, "sample_species_table")
cohort_species <- v2r4_read_tsv_raw(cohort_species_buffer, "SPECIES_COHORT_GATE")
sample_species <- v2r4_read_tsv_raw(sample_species_buffer, "SPECIES_SAMPLE_GATE")
v2r4_validate_species(cohort_species, sample_species, cfg)
conditional_inclusion <- v2r4_release_conditional_inclusion(
  ledger, cohort_species, sample_species,
  cohort_species_buffer$sha256, sample_species_buffer$sha256, cfg
)
snapshot <- v2r4_complete_snapshot(snapshot)
v2r4_assert(identical(snapshot$tree_digest, authorized_snapshot_tree_digest),
            "SNAPSHOT_TREE_DIGEST_NOT_EXACTLY_PRELEASE_SIGNED")

assembled_buffer <- v2r4_object(snapshot, "assembled_counts_rds")
assembled <- v2r4_read_rds_raw(
  assembled_buffer, "ASSEMBLED_COUNTS_RDS",
  cfg$input_snapshot$rds_role_contracts$assembled_counts_rds$compression
)
required_assembled <- c("dataset_id", "counts", "gene_annotation", "samples",
                        "strand_decision", "assignment_qc")
v2r4_assert(is.list(assembled) && all(required_assembled %in% names(assembled)),
            "ASSEMBLED_RDS_SCHEMA_INVALID")
v2r4_assert(identical(assembled$dataset_id, cfg$dataset_id), "ASSEMBLED_RDS_DATASET_MISMATCH")
ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
counts <- v2r4_validate_counts_original(assembled$counts, ids)
v2r4_assert(!is.null(rownames(counts)) && !anyDuplicated(rownames(counts)),
            "ASSEMBLED_COUNTS_GENE_IDS_INVALID")
if (identical(cli$execution_mode, "PRODUCTION")) {
  v2r4_assert(nrow(counts) == cfg$input_snapshot$expected_gene_n_before_filter,
              "ASSEMBLED_COUNTS_GENE_N_INVALID")
}
v2r4_assert(identical(as.data.frame(assembled$samples), as.data.frame(ledger)),
            "ASSEMBLED_EMBEDDED_LEDGER_MISMATCH")

assembly_manifest <- v2r4_read_tsv_raw(
  v2r4_object(snapshot, "assembly_output_manifest"), "ASSEMBLY_OUTPUT_MANIFEST"
)
rds_idx <- base::which(assembly_manifest[["filename"]] == "gene_counts_primary.rds")
rds_row <- assembly_manifest[rds_idx]
v2r4_assert(nrow(rds_row) == 1L && as.numeric(rds_row$bytes[[1L]]) == assembled_buffer$bytes &&
              tolower(rds_row$md5[[1L]]) == assembled_buffer$md5,
            "ASSEMBLY_MANIFEST_RDS_BINDING_INVALID")
assembly_complete <- v2r4_read_tsv_raw(v2r4_object(snapshot, "assembly_complete"), "ASSEMBLY_COMPLETE")
v2r4_assert(nrow(assembly_complete) == 1L && assembly_complete$dataset_id[[1L]] == cfg$dataset_id &&
              assembly_complete$n_libraries[[1L]] == 9L && assembly_complete$status[[1L]] == "PASS" &&
              assembly_complete$selected_primary_strand_mode[[1L]] == 0L,
            "ASSEMBLY_COMPLETE_SEMANTICS_INVALID")

legacy_contrasts <- v2r4_read_yaml_raw(v2r4_object(snapshot, "contrasts_config"), "CONTRASTS_CONFIG")
legacy_thresholds <- v2r4_read_yaml_raw(v2r4_object(snapshot, "thresholds_config"), "THRESHOLDS_CONFIG")
gse_cfg <- legacy_contrasts$datasets$GSE205906
v2r4_assert(identical(gse_cfg$model, "~ 0 + group") &&
              identical(gse_cfg$contrasts$local_disease$expression, "SNI - sham") &&
              identical(gse_cfg$contrasts$intervention$expression, "PLM - SNI"),
            "INTERNAL_SAP_DESIGN_OR_CONTRAST_DRIFT")
v2r4_assert(as.numeric(legacy_thresholds$expression_filtering$min_cpm) == 1 &&
              as.numeric(legacy_thresholds$false_discovery_rate$primary_gene_set_q) == 0.05,
            "INTERNAL_SAP_THRESHOLD_DRIFT")

design_contract <- v2r4_design(ledger, cfg)
group <- design_contract$group
design <- design_contract$design
contrast_matrix <- design_contract$contrast
seed <- as.integer(cfg$analysis$random_seed)
set.seed(seed)

annotation_contract <- v2r4_validate_annotation_gene(
  assembled$gene_annotation, rownames(counts), cfg$input_snapshot$minimum_nonblank_symbol_coverage
)
annotation <- annotation_contract$annotation
annotation_metrics <- annotation_contract$metrics

min_retained <- if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$minimum_retained_genes else
  cfg$analysis$minimum_retained_genes_production

v2r4_fit_preprocess <- function(count_matrix, design_matrix, minimum_retained) {
  y0 <- edgeR::DGEList(counts = count_matrix)
  keep_fbe <- edgeR::filterByExpr(y0, design = design_matrix)
  keep_cpm <- rowSums(edgeR::cpm(y0, log = FALSE) >= cfg$analysis$minimum_cpm) >=
    cfg$analysis$minimum_samples_at_cpm
  keep <- keep_fbe & keep_cpm
  v2r4_assert(sum(keep) >= minimum_retained, "TOO_FEW_FILTERED_GENES:%d", sum(keep))
  y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- edgeR::estimateDisp(y, design_matrix, robust = TRUE)
  ql <- edgeR::glmQLFit(y, design_matrix, robust = TRUE)
  voom <- limma::voom(y, design = design_matrix, plot = FALSE)
  list(y = y, ql = ql, voom = voom, keep = keep,
       keep_filter_by_expr = keep_fbe, keep_cpm = keep_cpm)
}

base_fit <- v2r4_fit_preprocess(counts, design, min_retained)
retained_annotation <- annotation[match(rownames(base_fit$y$counts), gene_id)]

v2r4_fit_edger <- function(fit, contrasts, ann) {
  data.table::rbindlist(lapply(seq_len(ncol(contrasts)), function(j) {
    qlf <- edgeR::glmQLFTest(fit, contrast = contrasts[, j])
    tab <- data.table::as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table,
                                    keep.rownames = "gene_id")
    tab <- tab[match(rownames(fit$counts), gene_id)]
    covariance <- v2r4_edgeR_ql_covariance(fit, contrasts[, j])
    crit <- stats::qt(0.975, covariance$df)
    tab[, `:=`(
      symbol = ann$symbol[match(gene_id, ann$gene_id)],
      contrast = colnames(contrasts)[[j]], method = "edgeR_QL_robust",
      se = covariance$se,
      ci_low = logFC - crit * covariance$se,
      ci_high = logFC + crit * covariance$se,
      se_df = covariance$df,
      se_status = covariance$status,
      df_total = qlf$df.total,
      q_value = FDR
    )]
    tab
  }), fill = TRUE)
}

v2r4_fit_voom <- function(v, design_matrix, contrasts, ann) {
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design_matrix), contrasts), robust = TRUE)
  data.table::rbindlist(lapply(seq_len(ncol(contrasts)), function(j) {
    se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
    crit <- stats::qt(0.975, fit$df.total)
    data.table::data.table(
      gene_id = rownames(fit$coefficients),
      symbol = ann$symbol[match(rownames(fit$coefficients), ann$gene_id)],
      contrast = colnames(contrasts)[[j]], method = "limma_voom_robust_eBayes",
      logFC = fit$coefficients[, j], t = fit$t[, j], p_value = fit$p.value[, j],
      q_value = p.adjust(fit$p.value[, j], "BH"), se = se,
      ci_low = fit$coefficients[, j] - crit * se,
      ci_high = fit$coefficients[, j] + crit * se, df_total = fit$df.total
    )
  }))
}

v2r4_fit_deseq <- function(filtered_counts, groups, contrasts, ann) {
  fit_messages_raw <- character()
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(filtered_counts),
    colData = data.frame(group = groups, row.names = colnames(filtered_counts)),
    design = ~ 0 + group
  )
  dds <- withCallingHandlers(
    DESeq2::DESeq(dds, test = "Wald",
                  fitType = cfg$analysis$DESeq2$requested_fit_type_before_fit,
                  betaPrior = FALSE, minReplicatesForReplace = Inf, quiet = TRUE),
    message = function(m) {
      fit_messages_raw <<- c(fit_messages_raw, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      fit_messages_raw <<- c(fit_messages_raw, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  actual_fit <- attr(DESeq2::dispersionFunction(dds), "fitType") %||% "UNKNOWN"
  fit_messages <- v2r4_collapse_condition_messages(fit_messages_raw)
  v2r4_assert_deseq_fit_provenance(actual_fit, fit_messages, cfg, "MODEL_DESEQ2")
  result_names <- DESeq2::resultsNames(dds)
  expected_names <- paste0("group", rownames(contrasts))
  v2r4_assert(all(expected_names %in% result_names), "DESEQ2_COEFFICIENT_NAMES_INVALID")
  tab <- data.table::rbindlist(lapply(seq_len(ncol(contrasts)), function(j) {
    cvec <- numeric(length(result_names)); names(cvec) <- result_names
    cvec[expected_names] <- contrasts[, j]
    res <- DESeq2::results(dds, contrast = cvec, alpha = 0.05,
                          independentFiltering = FALSE, cooksCutoff = FALSE)
    x <- data.table::as.data.table(as.data.frame(res), keep.rownames = "gene_id")
    x[, `:=`(
      symbol = ann$symbol[match(gene_id, ann$gene_id)],
      contrast = colnames(contrasts)[[j]], method = "DESeq2_Wald_sensitivity",
      requested_fit_type = cfg$analysis$DESeq2$requested_fit_type_before_fit,
      actual_fit_type_after_DESeq = actual_fit,
      ci_low = log2FoldChange - stats::qnorm(0.975) * lfcSE,
      ci_high = log2FoldChange + stats::qnorm(0.975) * lfcSE,
      p_value = pvalue, q_value = padj
    )]
    x
  }), fill = TRUE)
  list(table = tab, actual_fit_type = actual_fit, fit_messages = fit_messages)
}

gene_edger <- v2r4_add_claim_flags(v2r4_fit_edger(base_fit$ql, contrast_matrix, retained_annotation), analysis_role)
gene_voom <- v2r4_add_claim_flags(v2r4_fit_voom(base_fit$voom, design, contrast_matrix, retained_annotation), analysis_role)
deseq_result <- v2r4_fit_deseq(base_fit$y$counts, group, contrast_matrix, retained_annotation)
gene_deseq <- v2r4_add_claim_flags(deseq_result$table, analysis_role)

gene_sets <- v2r4_read_rds_raw(
  v2r4_object(snapshot, "gene_set_object"), "GENE_SET_OBJECT",
  cfg$input_snapshot$rds_role_contracts$gene_set_object$compression
)
orthology <- v2r4_read_tsv_raw(v2r4_object(snapshot, "orthology_single_source"),
                               "ORTHOLOGY_SINGLE_SOURCE")
required_gene_sets <- c("reactome", "progeny_human", "dorothea_rat", "syngo_sensitivity")
v2r4_assert(is.list(gene_sets) && all(required_gene_sets %in% names(gene_sets)),
          "GENE_SET_OBJECT_SCHEMA_INVALID")
v2r4_assert(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in% names(orthology)),
          "ORTHOLOGY_SCHEMA_INVALID")
rat_map <- unique(orthology[group_method == "NCBI Ortholog" & other_tax_id == 10116L,
                            .(human_symbol = as.character(human_symbol),
                              rat_symbol = as.character(other_symbol))])
v2r4_assert(!anyDuplicated(rat_map$human_symbol) && !anyDuplicated(rat_map$rat_symbol),
          "ORTHOLOGY_NOT_STRICT_ONE_TO_ONE")

v2r4_map_unsigned <- function(table, collection) {
  target_collection <- as.character(collection)
  x <- data.table::as.data.table(table)
  v2r4_assert(all(c("term_id", "term_name", "gene_symbol") %in% names(x)),
            "UNSIGNED_SOURCE_SCHEMA_INVALID:%s", collection)
  x <- unique(x[, .(program_id = as.character(term_id), program_name = as.character(term_name),
                    human_symbol = as.character(gene_symbol))])
  mapped <- unique(merge(x, rat_map, by = "human_symbol", allow.cartesian = TRUE)[,
                         .(program_id, program_name, symbol = rat_symbol)])
  by_program <- split(mapped$symbol, mapped$program_id)
  eligible <- lapply(by_program, function(z) unique(intersect(z, retained_annotation$symbol)))
  sizes <- lengths(eligible)
  eligible <- eligible[sizes >= cfg$analysis$program_size_min & sizes <= cfg$analysis$program_size_max]
  names_table <- unique(mapped[, .(program_id, program_name)])
  coverage <- unique(mapped[, .(mapped_symbol_n = data.table::uniqueN(symbol)), by = .(program_id, program_name)])
  coverage[, tested_feature_n := vapply(program_id, function(p) {
    data.table::uniqueN(intersect(mapped[program_id == p]$symbol, retained_annotation$symbol))
  }, integer(1))]
  coverage[, `:=`(
    collection = target_collection,
    source_feature_n = mapped_symbol_n,
    eligible_10_500 = tested_feature_n >= cfg$analysis$program_size_min &
      tested_feature_n <= cfg$analysis$program_size_max,
    tested = program_id %in% names(eligible),
    confidence_domain = "NOT_APPLICABLE",
    mapping_rule = cfg$analysis$orthology_authority
  )]
  indices <- lapply(eligible, function(z) which(retained_annotation$symbol %in% z))
  list(collection = target_collection, members = eligible, indices = indices,
       names = names_table, coverage = coverage)
}

unsigned_resources <- list(
  Reactome = v2r4_map_unsigned(gene_sets$reactome, "Reactome"),
  SynGO_sensitivity = v2r4_map_unsigned(gene_sets$syngo_sensitivity, "SynGO_sensitivity")
)

v2r4_unsigned_tests <- function(resource, v, design_matrix, contrasts) {
  if (!length(resource$indices)) return(list(camera = data.table::data.table(),
                                              mroast = data.table::data.table(),
                                              fgsea = data.table::data.table()))
  cams <- list(); roasts <- list(); fgs <- list()
  nrot <- if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$mroast_rotations else cfg$analysis$mroast_rotations
  for (j in seq_len(ncol(contrasts))) {
    cname <- colnames(contrasts)[[j]]
    cam <- data.table::as.data.table(limma::camera(v$E, resource$indices, design_matrix,
                                                   contrast = contrasts[, j],
                                                   inter.gene.cor = cfg$analysis$camera_inter_gene_correlation),
                                     keep.rownames = "program_id")
    cam[, `:=`(collection = resource$collection, contrast = cname,
               method = "camera_inter_gene_cor_0.01", p_value = PValue)]
    cam[, q_within_collection_contrast := p.adjust(p_value, "BH")]
    cams[[j]] <- cam
    set.seed(seed + j + match(resource$collection, names(unsigned_resources)) * 100L)
    roast <- data.table::as.data.table(limma::mroast(v$E, resource$indices, design_matrix,
                                                     contrast = contrasts[, j], nrot = nrot,
                                                     set.statistic = "mean"),
                                       keep.rownames = "program_id")
    roast[, `:=`(collection = resource$collection, contrast = cname,
                 method = paste0("mroast_mean_", nrot, "_rotations"), p_value = PValue)]
    roast[, q_within_collection_contrast := p.adjust(p_value, "BH")]
    roasts[[j]] <- roast
    lf <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design_matrix), contrasts[, j, drop = FALSE]),
                        robust = TRUE)
    stat_dt <- data.table::data.table(symbol = retained_annotation$symbol,
                                      stat = as.numeric(lf$t[, 1L]))
    stat_dt <- stat_dt[is.finite(stat) & !is.na(symbol) & nzchar(symbol)]
    stat_dt <- stat_dt[order(-abs(stat)), .SD[1L], by = symbol]
    ranks <- stat_dt$stat; names(ranks) <- stat_dt$symbol
    set.seed(seed + 1000L + j + match(resource$collection, names(unsigned_resources)) * 100L)
    fg <- data.table::as.data.table(fgsea::fgseaMultilevel(
      pathways = resource$members, stats = ranks,
      minSize = cfg$analysis$program_size_min, maxSize = cfg$analysis$program_size_max,
      eps = cfg$analysis$fgsea_eps
    ))
    if (nrow(fg)) {
      data.table::setnames(fg, "pathway", "program_id")
      fg[, `:=`(collection = resource$collection, contrast = cname,
                method = "fgseaMultilevel", p_value = pval,
                q_within_collection_contrast = padj)]
    }
    fgs[[j]] <- fg
  }
  list(camera = data.table::rbindlist(cams, fill = TRUE),
       mroast = data.table::rbindlist(roasts, fill = TRUE),
       fgsea = data.table::rbindlist(fgs, fill = TRUE))
}

unsigned_all <- lapply(unsigned_resources, v2r4_unsigned_tests,
                       v = base_fit$voom, design_matrix = design, contrasts = contrast_matrix)
unsigned_camera <- v2r4_add_claim_flags(data.table::rbindlist(lapply(unsigned_all, `[[`, "camera"), fill = TRUE), analysis_role)
unsigned_mroast <- v2r4_add_claim_flags(data.table::rbindlist(lapply(unsigned_all, `[[`, "mroast"), fill = TRUE), analysis_role)
unsigned_fgsea <- v2r4_add_claim_flags(data.table::rbindlist(lapply(unsigned_all, `[[`, "fgsea"), fill = TRUE), analysis_role)

v2r4_collapse_symbols <- function(expr, ann) {
  dt <- data.table::data.table(symbol = ann$symbol, data.table::as.data.table(expr))
  dt <- dt[!is.na(symbol) & nzchar(symbol), lapply(.SD, mean), by = symbol]
  rn <- dt$symbol; dt[, symbol := NULL]
  out <- as.matrix(dt); rownames(out) <- rn; out
}

v2r4_activity_domain <- function(v, ann, enforce_expected = TRUE) {
  expr <- v2r4_collapse_symbols(v$E, ann)
  prog_h <- data.table::as.data.table(gene_sets$progeny_human)
  v2r4_assert(all(c("pathway", "gene_symbol", "weight") %in% names(prog_h)),
            "PROGENY_SOURCE_SCHEMA_INVALID")
  prog <- unique(merge(
    prog_h[, .(source = as.character(pathway), human_symbol = as.character(gene_symbol),
               mor = as.numeric(weight))], rat_map, by = "human_symbol", allow.cartesian = TRUE
  )[, .(source, target = rat_symbol, mor)])
  doro_raw <- data.table::as.data.table(gene_sets$dorothea_rat)
  v2r4_assert(all(c("tf", "target", "mor", "confidence") %in% names(doro_raw)),
            "DOROTHEA_SOURCE_SCHEMA_INVALID")
  doro <- unique(doro_raw[, .(source = as.character(tf), target = as.character(target),
                              mor = as.numeric(mor), confidence = as.character(confidence))])
  v2r4_assert(all(doro$confidence %in% c("A", "B")), "DOROTHEA_CONFIDENCE_NOT_EXACT_AB")
  domains <- if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$program_domains else cfg$analysis$program_domains
  v2r4_assert(data.table::uniqueN(prog$source) == domains$PROGENy$source_program_n,
            "PROGENY_SOURCE_PROGRAM_COUNT_INVALID")
  v2r4_assert(setequal(unique(as.character(prog$source)),
                     unlist(domains$PROGENy$exact_program_ids, use.names = FALSE)),
            "PROGENY_SOURCE_PROGRAM_ID_SET_INVALID")
  v2r4_assert(data.table::uniqueN(doro$source) == domains$DoRothEA_AB$frozen_source_program_upper_bound,
            "DOROTHEA_SOURCE_PROGRAM_UPPER_BOUND_INVALID")
  coverage_one <- function(net, collection, confidence_domain) {
    target_collection <- as.character(collection)
    target_confidence_domain <- as.character(confidence_domain)
    z <- net[, .(source_feature_n = data.table::uniqueN(target),
                 tested_feature_n = data.table::uniqueN(intersect(target, rownames(expr)))), by = source]
    z[, `:=`(
      collection = target_collection,
      eligible_10_500 = tested_feature_n >= cfg$analysis$program_size_min &
        tested_feature_n <= cfg$analysis$program_size_max,
      tested = tested_feature_n >= cfg$analysis$program_size_min &
        tested_feature_n <= cfg$analysis$program_size_max,
      confidence_domain = target_confidence_domain,
      mapping_rule = if (target_collection == "PROGENy") cfg$analysis$orthology_authority else
        "frozen_rat_AB_strict_one_to_one_network"
    )]
    z
  }
  pcov <- coverage_one(prog, "PROGENy", "NOT_APPLICABLE")
  dcov <- coverage_one(doro, "DoRothEA_AB", paste(sort(unique(doro$confidence)), collapse = ","))
  if (enforce_expected) {
    v2r4_assert(sum(pcov$tested) == domains$PROGENy$eligible_tested_program_n,
              "PROGENY_ELIGIBLE_TESTED_PROGRAM_COUNT_INVALID")
    v2r4_assert(sum(dcov$tested) == domains$DoRothEA_AB$eligible_tested_program_n,
              "DOROTHEA_ELIGIBLE_TESTED_PROGRAM_COUNT_INVALID")
  }
  prog_keep <- pcov[tested == TRUE]$source
  doro_keep <- dcov[tested == TRUE]$source
  prog_test <- prog[source %in% prog_keep & target %in% rownames(expr)]
  doro_test <- doro[source %in% doro_keep & target %in% rownames(expr)]
  v2r4_assert(all(prog_test[, data.table::uniqueN(target), by = source]$V1 <= cfg$analysis$program_size_max),
            "PROGENY_MAXSIZE_FILTER_FAILED")
  v2r4_assert(all(doro_test[, data.table::uniqueN(target), by = source]$V1 <= cfg$analysis$program_size_max),
            "DOROTHEA_MAXSIZE_FILTER_FAILED")
  run_one <- function(net) {
    if (!nrow(net)) return(matrix(numeric(), 0L, ncol(expr), dimnames = list(NULL, colnames(expr))))
    long <- decoupleR::run_ulm(mat = expr, network = net[, .(source, target, mor)],
                               .source = source, .target = target, .mor = mor,
                               minsize = cfg$analysis$program_size_min)
    long <- data.table::as.data.table(long)
    v2r4_assert(all(c("source", "condition", "score") %in% names(long)),
              "DECOUPLER_ULM_OUTPUT_SCHEMA_INVALID")
    wide <- data.table::dcast(long, source ~ condition, value.var = "score")
    rn <- wide$source; wide[, source := NULL]
    mat <- as.matrix(wide); rownames(mat) <- rn
    mat[, colnames(expr), drop = FALSE]
  }
  list(scores = list(PROGENy = run_one(prog_test), DoRothEA_AB = run_one(doro_test)),
       coverage = data.table::rbindlist(list(pcov, dcov), fill = TRUE))
}

v2r4_fit_activity <- function(score, collection, design_matrix, contrasts) {
  v2r4_assert(nrow(score) > 0L && identical(colnames(score), rownames(design_matrix)),
            "ACTIVITY_MATRIX_INVALID:%s", collection)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(score, design_matrix), contrasts), robust = TRUE)
  out <- data.table::rbindlist(lapply(seq_len(ncol(contrasts)), function(j) {
    se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
    crit <- stats::qt(0.975, fit$df.total)
    data.table::data.table(
      collection = collection, program_id = rownames(score),
      contrast = colnames(contrasts)[[j]], method = "decoupleR_ULM_then_limma",
      effect = fit$coefficients[, j], se = se,
      ci_low = fit$coefficients[, j] - crit * se,
      ci_high = fit$coefficients[, j] + crit * se,
      t = fit$t[, j], df_total = fit$df.total, p_value = fit$p.value[, j]
    )
  }))
  out[, q_within_collection_contrast := p.adjust(p_value, "BH"), by = .(collection, contrast)]
  out
}

full_activity <- v2r4_activity_domain(base_fit$voom, retained_annotation, enforce_expected = TRUE)
activity_results <- data.table::rbindlist(lapply(names(full_activity$scores), function(nm) {
  v2r4_fit_activity(full_activity$scores[[nm]], nm, design, contrast_matrix)
}), fill = TRUE)
activity_results <- v2r4_add_claim_flags(activity_results, analysis_role)

unsigned_coverage <- data.table::rbindlist(lapply(unsigned_resources, `[[`, "coverage"), fill = TRUE)
mapping_coverage <- data.table::rbindlist(list(unsigned_coverage, full_activity$coverage), fill = TRUE)
mapping_coverage <- v2r4_add_claim_flags(mapping_coverage, analysis_role)

## Full-pipeline activity LOAO: every omission recomputes filtering, TMM,
## voom, symbol aggregation, 10-500 eligibility, ULM scores, and limma contrasts.
loao_rows <- list()
expected_programs <- lapply(full_activity$scores, rownames)
for (omit_i in seq_along(ids)) {
  keep_samples <- seq_along(ids) != omit_i
  sub_ids <- ids[keep_samples]
  sub_groups <- factor(as.character(group[keep_samples]), levels = levels(group))
  sub_design <- stats::model.matrix(~ 0 + sub_groups)
  colnames(sub_design) <- colnames(design)
  rownames(sub_design) <- sub_ids
  v2r4_assert(qr(sub_design)$rank == 3L, "LOAO_DESIGN_RANK_INVALID:%s", ids[[omit_i]])
  sub_fit <- v2r4_fit_preprocess(counts[, keep_samples, drop = FALSE], sub_design,
                               if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$minimum_retained_genes else
                                 cfg$analysis$minimum_retained_genes_production)
  sub_ann <- annotation[match(rownames(sub_fit$y$counts), gene_id)]
  sub_activity <- v2r4_activity_domain(sub_fit$voom, sub_ann, enforce_expected = FALSE)
  for (collection in names(expected_programs)) {
    collection_value <- as.character(collection)
    available <- rownames(sub_activity$scores[[collection]])
    fit_table <- if (length(available)) {
      v2r4_fit_activity(sub_activity$scores[[collection]], collection, sub_design, contrast_matrix)
    } else data.table::data.table()
    grid <- data.table::CJ(program_id = expected_programs[[collection]],
                           contrast = colnames(contrast_matrix), unique = TRUE)
    grid[, `:=`(collection = collection_value, omitted_animal = ids[[omit_i]],
                omitted_group = as.character(group[[omit_i]]),
                loao_method = "ACTIVITY_FULL_PIPELINE_REFIT")]
    if (nrow(fit_table)) {
      grid <- merge(grid, fit_table[, .(program_id, contrast, estimate = effect, se, ci_low, ci_high,
                                        p_value, q_within_collection_contrast)],
                    by = c("program_id", "contrast"), all.x = TRUE, sort = FALSE)
    } else {
      grid[, `:=`(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                  p_value = NA_real_, q_within_collection_contrast = NA_real_)]
    }
    grid[, `:=`(estimable = is.finite(estimate),
                nonestimable_reason = ifelse(is.finite(estimate), "NONE", "PROGRAM_INELIGIBLE_OR_NOT_SCORED_AFTER_FULL_REFIT"))]
    loao_rows[[length(loao_rows) + 1L]] <- grid
  }
}
loao_long <- data.table::rbindlist(loao_rows, fill = TRUE)
loao_long <- v2r4_add_claim_flags(loao_long, analysis_role)
loao_summary <- v2r4_recompute_loao_summary(loao_long, ids)
loao_summary[, `:=`(all_steps_refit = TRUE, loao_method = "ACTIVITY_FULL_PIPELINE_REFIT")]
min_estimable <- if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$minimum_estimable_omissions_per_program_contrast else
  cfg$analysis$loao$minimum_estimable_omissions_per_program_contrast
v2r4_assert(all(loao_summary$n_omissions == 9L) && all(loao_summary$n_estimable >= min_estimable),
          "LOAO_INSUFFICIENT_ESTIMABLE_OMISSIONS")
loao_summary <- v2r4_add_claim_flags(loao_summary, analysis_role)

parent <- dirname(output_dir)
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
v2r4_assert(!file.exists(output_dir), "OUTPUT_ALREADY_EXISTS")
staging <- file.path(parent, paste0(basename(output_dir), ".staging_v2r4_", Sys.getpid(), "_",
                                   sprintf("%08x", sample.int(.Machine$integer.max, 1L))))
v2r4_assert(!file.exists(staging), "STAGING_ALREADY_EXISTS")
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
for (d in c("audit", "diagnostics", "gene_level", "programs", "gates")) {
  dir.create(file.path(staging, d), recursive = TRUE, showWarnings = FALSE)
}

design_out <- data.table::data.table(sample_id = rownames(design), group = as.character(group),
                                     data.table::as.data.table(design))
contrast_out <- data.table::data.table(coefficient = rownames(contrast_matrix),
                                       data.table::as.data.table(contrast_matrix))
filter_out <- data.table::data.table(
  dataset_id = cfg$dataset_id, input_gene_n = nrow(counts), retained_gene_n = sum(base_fit$keep),
  filterByExpr_retained_n = sum(base_fit$keep_filter_by_expr),
  explicit_CPM_retained_n = sum(base_fit$keep_cpm), minimum_cpm = cfg$analysis$minimum_cpm,
  minimum_samples = cfg$analysis$minimum_samples_at_cpm,
  filtering = cfg$analysis$edgeR$filter, normalization = cfg$analysis$edgeR$normalization,
  primary_method = cfg$analysis$primary_method
)
gate <- data.table::data.table(
  dataset_id = cfg$dataset_id, gate = "E2", upstream_gate = "E1", eligible_input_n = 0L,
  evaluated_candidate_n = 0L, pass_n = NA_integer_, status = cfg$gate$machine_status,
  analysis_role = analysis_role, candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE, formal_sequential_inference = FALSE
)
design_out <- v2r4_add_claim_flags(design_out, analysis_role)
contrast_out <- v2r4_add_claim_flags(contrast_out, analysis_role)
filter_out <- v2r4_add_claim_flags(filter_out, analysis_role)
annotation_metrics <- v2r4_add_claim_flags(annotation_metrics, analysis_role)
runtime_bindings <- v2r4_add_claim_flags(runtime_bindings, analysis_role)
verification <- v2r4_add_claim_flags(snapshot$verification, analysis_role)
lease_validation <- v2r4_add_claim_flags(data.table::data.table(
  dataset_id = cfg$dataset_id,
  model_version = cfg$model_version,
  lease_id = lease$lease_id[[1L]],
  holder_principal = lease$holder_principal[[1L]],
  run_id = lease$run_id[[1L]],
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  model_pid = as.integer(lease$model_pid[[1L]]),
  model_start_utc = lease$model_start_utc[[1L]],
  model_request_sha256 = lease$model_request_sha256[[1L]],
  model_request_nonce = lease$model_request_nonce[[1L]],
  model_identity_request_sha256 = lease$model_identity_request_sha256[[1L]],
  validator_pid = as.integer(lease$validator_pid[[1L]]),
  validator_start_utc = lease$validator_start_utc[[1L]],
  validator_request_sha256 = lease$validator_request_sha256[[1L]],
  validator_request_nonce = lease$validator_request_nonce[[1L]],
  validator_identity_request_sha256 = lease$validator_identity_request_sha256[[1L]],
  guard_owner_role = lease$guard_owner_role[[1L]],
  guard_owner_pid = as.integer(lease$guard_owner_pid[[1L]]),
  guard_owner_start_utc = lease$guard_owner_start_utc[[1L]],
  snapshot_tree_sha256 = lease$snapshot_tree_sha256[[1L]],
  mutex_relative_path = lease$mutex_relative_path[[1L]],
  guard_pid = as.integer(lease$guard_pid[[1L]]),
  guard_start_utc = lease$guard_start_utc[[1L]],
  guard_nonce = lease$guard_nonce[[1L]],
  guard_public_key_sha256 = lease$guard_public_key_sha256[[1L]],
  guard_ready_receipt_sha256 = lease$guard_ready_receipt_sha256[[1L]],
  guard_challenge_request_sha256 = lease$guard_challenge_request_sha256[[1L]],
  guard_challenge_response_sha256 = lease$guard_challenge_response_sha256[[1L]],
  guard_mutex_file_id = lease$guard_mutex_file_id[[1L]],
  guard_stop_hmac_key_sha256 = lease$guard_stop_hmac_key_sha256[[1L]],
  finalization_public_key_sha256 = lease$finalization_public_key_sha256[[1L]],
  finalization_public_key_captured_path =
    lease$finalization_public_key_captured_path[[1L]],
  finalization_public_key_captured_bytes =
    as.integer(lease$finalization_public_key_captured_bytes[[1L]]),
  finalization_public_key_captured_sha256 =
    lease$finalization_public_key_captured_sha256[[1L]],
  status = "PASS_SIGNED_LIVE_LEASE_AND_HELD_OS_MUTEX"
), analysis_role)
observed_group_counts <- table(group)
group_count_text <- paste(names(observed_group_counts), as.integer(observed_group_counts),
                          sep = "=", collapse = ";")
run_manifest <- data.table::data.table(
  dataset_id = cfg$dataset_id, model_version = cfg$model_version,
  implementation_status = cfg$status,
  execution_mode = if (identical(cli$execution_mode, "SYNTHETIC")) "SYNTHETIC_ISOLATED_ONLY" else "PRODUCTION_AUTHORITY_GATED",
  source_authority_manifest_sha256 = source_bundle$manifest$sha256,
  source_authority_public_key_sha256 = source_bundle$public_key$sha256,
  source_authority_decision = expected_source_decision,
  external_source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
  external_launch_context_sha256 = launch_context$context_sha256,
  external_launcher_sha256 = launch_context$bootstrap_launcher_sha256,
  runtime_pre_digest = launch_context$runtime_pre_digest,
  runtime_semantic_digest = launch_context$runtime_semantic_digest,
  snapshot_tree_sha256 = snapshot$tree_digest,
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  lease_id = lease$lease_id[[1L]],
  holder_principal = holder_principal,
  run_id = run_id,
  launcher_pid = as.integer(launch_context$launcher_pid),
  launcher_start_utc = launch_context$launcher_start_utc,
  guard_owner_role = guard_state$owner_role,
  guard_owner_pid = guard_state$owner_pid,
  guard_owner_start_utc = guard_state$owner_start_utc,
  guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
  finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
  finalization_public_key_captured_path = guard_state$finalization_public_key_captured_path,
  finalization_public_key_captured_bytes = guard_state$finalization_public_key_captured_bytes,
  finalization_public_key_captured_sha256 = guard_state$finalization_public_key_captured_sha256,
  model_pid = as.integer(model_identity$pid),
  model_start_utc = model_identity$start_utc,
  model_request_sha256 = child_request$buffer$sha256,
  model_request_nonce = child_request$request$request_nonce,
  model_identity_request_sha256 = model_identity_request_buffer$sha256,
  validator_pid = as.integer(lease$validator_pid[[1L]]),
  validator_start_utc = lease$validator_start_utc[[1L]],
  validator_request_sha256 = validator_request_buffer_prelease$sha256,
  validator_request_nonce = validator_request_bundle_prelease$request$request_nonce,
  validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  mutex_relative_path = cfg$lease$mutex_relative_path,
  guard_pid = guard_state$guard_pid,
  guard_start_utc = guard_state$guard_start_utc,
  guard_nonce = guard_state$guard_nonce,
  guard_public_key_sha256 = guard_state$guard_public_key_sha256,
  guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
  producer_source_provenance_status = producer_source_provenance$status,
  producer_source_digest = producer_source_provenance$digest,
  statistical_unit = cfg$design$biological_unit,
  independent_rat_libraries = as.integer(nrow(ledger)),
  unique_biosamples = as.integer(data.table::uniqueN(ledger$biosample)),
  unique_runs = as.integer(data.table::uniqueN(ledger$run_accession)),
  animal_independence_basis = "GEO_ENA_ONE_TO_ONE_BIOLOGICAL_REPLICATE_METADATA",
  repository_inferred_animal_n = as.integer(nrow(ledger)), group_counts = group_count_text,
  sample_include_primary_marker = conditional_inclusion$sample_include_primary_marker,
  sample_metadata_status = conditional_inclusion$sample_metadata_status,
  conditional_inclusion_release_status =
    conditional_inclusion$conditional_inclusion_release_status,
  conditional_inclusion_n = conditional_inclusion$conditional_inclusion_n,
  cohort_species_gate_sha256 = conditional_inclusion$cohort_species_gate_sha256,
  sample_species_table_sha256 = conditional_inclusion$sample_species_table_sha256,
  plm_metadata_conflict_policy = conditional_inclusion$plm_metadata_conflict_policy,
  design_formula = cfg$design$formula,
  design_rank = qr(design)$rank, coefficient_order = paste(colnames(design), collapse = ","),
  local_disease_expression = cfg$design$contrasts$local_disease$expression,
  intervention_expression = cfg$design$contrasts$intervention$expression,
  primary_method = cfg$analysis$primary_method,
  primary_filtering = cfg$analysis$edgeR$filter, normalization = cfg$analysis$edgeR$normalization,
  random_seed = seed, orthology_authority = cfg$analysis$orthology_authority,
  orthology_single_source_role = cfg$analysis$orthology_single_source_role,
  assembly = cfg$input_snapshot$assembly, assembly_accession = cfg$input_snapshot$assembly_accession,
  selected_primary_strand_mode = 0L,
  species_gate_status = "PASS_BEFORE_COUNTS_RDS_CAPTURE_OR_SEMANTIC_PARSE",
  gene_id_source_column = annotation_metrics$gene_id_source_column[[1L]],
  symbol_source_column = annotation_metrics$symbol_source_column[[1L]],
  biotype_source_column = annotation_metrics$biotype_source_column[[1L]],
  nonblank_symbol_coverage = annotation_metrics$nonblank_symbol_coverage[[1L]],
  unique_nonblank_symbol_coverage = annotation_metrics$unique_nonblank_symbol_coverage[[1L]],
  fallback_to_gene_id_n = annotation_metrics$fallback_to_gene_id_n[[1L]],
  requested_DESeq2_fit_type = cfg$analysis$DESeq2$requested_fit_type_before_fit,
  actual_DESeq2_fit_type_after_DESeq = deseq_result$actual_fit_type,
  DESeq2_fit_messages = deseq_result$fit_messages,
  progeny_source_program_n = data.table::uniqueN(data.table::as.data.table(gene_sets$progeny_human)$pathway),
  progeny_tested_program_n = sum(full_activity$coverage$collection == "PROGENy" & full_activity$coverage$tested),
  dorothea_source_program_n = data.table::uniqueN(data.table::as.data.table(gene_sets$dorothea_rat)$tf),
  dorothea_tested_program_n = sum(full_activity$coverage$collection == "DoRothEA_AB" & full_activity$coverage$tested),
  dorothea_confidence_domain = paste(sort(unique(data.table::as.data.table(gene_sets$dorothea_rat)$confidence)), collapse = ","),
  loao_method = "ACTIVITY_FULL_PIPELINE_REFIT",
  gate_status = cfg$gate$machine_status,
  analysis_role = analysis_role,
  candidate_gate_eligible = FALSE, can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)

v2r4_write_tsv(verification, file.path(staging, "audit", "input_snapshot_verification.tsv"))
v2r4_write_tsv(lease_validation, file.path(staging, "audit", "lease_validation.tsv"))
v2r4_write_tsv(run_manifest, file.path(staging, "audit", "run_manifest.tsv"))
v2r4_write_tsv(runtime_bindings, file.path(staging, "audit", "runtime_bindings.tsv"))
v2r4_write_tsv(annotation_metrics, file.path(staging, "diagnostics", "annotation_mapping.tsv"))
v2r4_write_tsv(design_out, file.path(staging, "diagnostics", "design_matrix.tsv"))
v2r4_write_tsv(contrast_out, file.path(staging, "diagnostics", "contrast_matrix.tsv"))
v2r4_write_tsv(filter_out, file.path(staging, "diagnostics", "filter_summary.tsv"))
v2r4_write_tsv(gene_edger, file.path(staging, "gene_level", "edgeR_QL_all_contrasts.tsv.gz"))
v2r4_write_tsv(gene_voom, file.path(staging, "gene_level", "limma_voom_all_contrasts.tsv.gz"))
v2r4_write_tsv(gene_deseq, file.path(staging, "gene_level", "DESeq2_all_contrasts.tsv.gz"))
v2r4_write_tsv(gate, file.path(staging, "gates", "E2_GATE.tsv"))
v2r4_write_tsv(activity_results, file.path(staging, "programs", "activity_contrasts.tsv.gz"))
v2r4_write_tsv(loao_long, file.path(staging, "programs", "activity_full_pipeline_LOAO_long.tsv.gz"))
v2r4_write_tsv(loao_summary, file.path(staging, "programs", "activity_full_pipeline_LOAO_summary.tsv.gz"))
v2r4_write_tsv(mapping_coverage, file.path(staging, "programs", "program_mapping_coverage.tsv"))
v2r4_write_tsv(unsigned_camera, file.path(staging, "programs", "unsigned_camera.tsv.gz"))
v2r4_write_tsv(unsigned_fgsea, file.path(staging, "programs", "unsigned_fgseaMultilevel.tsv.gz"))
v2r4_write_tsv(unsigned_mroast, file.path(staging, "programs", "unsigned_mroast.tsv.gz"))

required_payloads <- unlist(cfg$required_output_payloads, use.names = FALSE)
out_manifest <- v2r4_output_manifest(staging, required_payloads)
v2r4_write_tsv(out_manifest, file.path(staging, "OUTPUT_MANIFEST.tsv"))
output_manifest_buffer <- v2r4_capture_file(file.path(staging, "OUTPUT_MANIFEST.tsv"),
                                            "MODEL_OUTPUT_MANIFEST", staging)
complete <- data.table::data.table(
  dataset_id = cfg$dataset_id, model_version = cfg$model_version,
  implementation_status = cfg$status,
  computation_status = if (identical(cli$execution_mode, "SYNTHETIC")) "PASS_SYNTHETIC_ISOLATED_ONLY" else
    "COMPUTATION_COMPLETE_AWAITING_INDEPENDENT_REDTEAM_NOT_AUTHORIZED",
  run_id = run_id,
  model_pid = as.integer(Sys.getpid()),
  snapshot_tree_sha256 = snapshot$tree_digest,
  output_manifest_sha256 = output_manifest_buffer$sha256,
  gate_status = cfg$gate$machine_status,
  analysis_role = analysis_role,
  candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)
v2r4_write_tsv(complete, file.path(staging, "COMPLETE.tsv"))
v2r4_assert(file.rename(staging, output_dir), "ATOMIC_OUTPUT_PUBLICATION_FAILED")
published_output_manifest <- v2r4_capture_file(file.path(output_dir, "OUTPUT_MANIFEST.tsv"),
                                                "PUBLISHED_OUTPUT_MANIFEST", output_dir)
published_complete <- v2r4_capture_file(file.path(output_dir, "COMPLETE.tsv"),
                                         "PUBLISHED_COMPLETE", output_dir)
v2r4_assert(identical(published_output_manifest$sha256, output_manifest_buffer$sha256),
            "PUBLISHED_OUTPUT_MANIFEST_IDENTITY_DRIFT")

final_wait <- cfg$execution_contract$signed_final_ack_wait_seconds
v2r4_assert(identical(final_wait, 900L) &&
              identical(launch_context$final_ack_wait_seconds, 900L),
            "MODEL_FINAL_CONTROL_WAIT_MUST_BE_CAPTURED_EXACT_900")
v2r4_wait_for_files(c(
  validator_request_path, science_complete_path, finalization_request_path,
  runtime_post_receipt_path, cli$final_challenge_request,
  cli$final_challenge_response, control_complete_path
), final_wait, "MODEL_TWO_PHASE_FINAL_CONTROL")

validator_request_buffer <- v2r4_capture_file(
  validator_request_path, "MODEL_CAPTURED_VALIDATOR_REQUEST", attempt_root
)
validator_request_bundle <- v2r4_parse_child_request_buffer(
  validator_request_buffer, validator_request_buffer$sha256, "VALIDATOR"
)
v2r4_assert(
  identical(validator_request_buffer$sha256, validator_request_buffer_prelease$sha256) &&
    identical(validator_request_bundle$request$request_nonce,
              validator_request_bundle_prelease$request$request_nonce),
  "VALIDATOR_CHILD_REQUEST_DRIFT_AFTER_SIGNED_LEASE"
)
v2r4_validate_child_request_bindings(
  validator_request_bundle, cfg, launch_context, "VALIDATOR"
)
v2r4_assert(
  identical(child_request$request$peer_request_nonce,
            validator_request_bundle$request$request_nonce) &&
    identical(validator_request_bundle$request$peer_request_nonce,
              child_request$request$request_nonce),
  "MODEL_VALIDATOR_REQUEST_NONCE_CROSS_BINDING_INVALID"
)
finalization_request_buffer <- v2r4_capture_file(
  finalization_request_path, "MODEL_FINALIZATION_REQUEST", attempt_root
)
science <- v2r4_validate_science_complete(
  science_complete_path, cfg,
  list(
    run_id = run_id,
    validator_pid = as.integer(lease$validator_pid[[1L]]),
    validator_start_utc = lease$validator_start_utc[[1L]],
    validator_request_sha256 = validator_request_buffer$sha256,
    validator_request_nonce = validator_request_bundle$request$request_nonce,
    validation_report_sha256 = v2r4_capture_file(
      launch_context$validation_report_path, "MODEL_VALIDATION_REPORT", attempt_root
    )$sha256,
    output_manifest_sha256 = published_output_manifest$sha256,
    complete_sha256 = published_complete$sha256,
    finalization_request_sha256 = finalization_request_buffer$sha256
  ), attempt_root
)
runtime_post <- v2r4_validate_runtime_post_receipt(
  runtime_post_receipt_path, cfg,
  list(
    run_id = run_id,
    runtime_pre_digest = launch_context$runtime_pre_digest,
    runtime_post_digest = launch_context$runtime_pre_digest,
    runtime_pre_post_equal = TRUE,
    runtime_semantic_digest = launch_context$runtime_semantic_digest,
    execution_identity_manifest_sha256 =
      launch_context$execution_identity_manifest_sha256,
    analysis_config_sha256 = config_buffer$sha256,
    science_complete_sha256 = science$buffer$sha256,
    validator_pid = as.integer(lease$validator_pid[[1L]]),
    validator_start_utc = lease$validator_start_utc[[1L]],
    validator_request_sha256 = validator_request_buffer$sha256,
    validator_request_nonce = validator_request_bundle$request$request_nonce,
    wait_seconds = 900L
  ), attempt_root
)
challenge_nonce_expected <- sub(
  "^final_challenge_([0-9a-f]{64})\\.request\\.json$", "\\1",
  basename(cli$final_challenge_request)
)
v2r4_assert(grepl("^[0-9a-f]{64}$", challenge_nonce_expected),
            "MODEL_FINAL_CHALLENGE_PATH_NONCE_INVALID")
challenge <- v2r4_verify_final_challenge(
  cli$final_challenge_request, cli$final_challenge_response, cfg, guard_state,
  list(
    run_id = run_id,
    lease_id = lease$lease_id[[1L]],
    lease_manifest_sha256 = lease_bundle$manifest$sha256,
    source_authority_payload_sha256 =
      launch_context$source_authorization_payload_sha256,
    launcher_pid = as.integer(launch_context$launcher_pid),
    launcher_start_utc = launch_context$launcher_start_utc,
    guard_owner_role = guard_state$owner_role,
    guard_owner_pid = guard_state$owner_pid,
    guard_owner_start_utc = guard_state$owner_start_utc,
    guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
    finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
    finalization_public_key_captured_path =
      guard_state$finalization_public_key_captured_path,
    finalization_public_key_captured_bytes =
      guard_state$finalization_public_key_captured_bytes,
    finalization_public_key_captured_sha256 =
      guard_state$finalization_public_key_captured_sha256,
    model_pid = as.integer(Sys.getpid()),
    model_start_utc = model_identity$start_utc,
    validator_pid = as.integer(lease$validator_pid[[1L]]),
    validator_start_utc = lease$validator_start_utc[[1L]],
    guard_pid = guard_state$guard_pid,
    guard_start_utc = guard_state$guard_start_utc,
    guard_nonce = guard_state$guard_nonce,
    guard_public_key_sha256 = guard_state$guard_public_key_sha256,
    guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
    guard_mutex_file_id = guard_state$mutex_file_id,
    model_request_sha256 = child_request$buffer$sha256,
    model_request_nonce = child_request$request$request_nonce,
    model_identity_request_sha256 = model_identity_request_buffer$sha256,
    validator_request_sha256 = validator_request_buffer$sha256,
    validator_request_nonce = validator_request_bundle$request$request_nonce,
    validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
    science_complete_sha256 = science$buffer$sha256,
    finalization_request_sha256 = finalization_request_buffer$sha256,
    validation_report_sha256 = science$value$validation_report_sha256,
    output_manifest_sha256 = published_output_manifest$sha256,
    complete_sha256 = published_complete$sha256,
    snapshot_tree_sha256 = snapshot$tree_digest,
    runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
    runtime_pre_digest = launch_context$runtime_pre_digest,
    runtime_post_digest = launch_context$runtime_pre_digest,
    runtime_pre_post_equal = TRUE,
    challenge_nonce = challenge_nonce_expected
  ), attempt_root
)
control <- v2r4_validate_control_complete(
  control_complete_path, cfg,
  list(
    run_id = run_id,
    lease_id = lease$lease_id[[1L]],
    lease_manifest_sha256 = lease_bundle$manifest$sha256,
    launcher_pid = as.integer(launch_context$launcher_pid),
    launcher_start_utc = launch_context$launcher_start_utc,
    guard_owner_role = guard_state$owner_role,
    guard_owner_pid = guard_state$owner_pid,
    guard_owner_start_utc = guard_state$owner_start_utc,
    guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
    finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
    finalization_public_key_captured_path =
      guard_state$finalization_public_key_captured_path,
    finalization_public_key_captured_bytes =
      guard_state$finalization_public_key_captured_bytes,
    finalization_public_key_captured_sha256 =
      guard_state$finalization_public_key_captured_sha256,
    model_pid = as.integer(Sys.getpid()),
    model_start_utc = model_identity$start_utc,
    validator_pid = as.integer(lease$validator_pid[[1L]]),
    validator_start_utc = lease$validator_start_utc[[1L]],
    guard_pid = guard_state$guard_pid,
    guard_start_utc = guard_state$guard_start_utc,
    guard_nonce = guard_state$guard_nonce,
    guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
    guard_mutex_file_id = guard_state$mutex_file_id,
    model_request_sha256 = child_request$buffer$sha256,
    model_request_nonce = child_request$request$request_nonce,
    model_identity_request_sha256 = model_identity_request_buffer$sha256,
    validator_request_sha256 = validator_request_buffer$sha256,
    validator_request_nonce = validator_request_bundle$request$request_nonce,
    validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
    science_complete_sha256 = science$buffer$sha256,
    finalization_request_sha256 = finalization_request_buffer$sha256,
    validation_report_sha256 = science$value$validation_report_sha256,
    output_manifest_sha256 = published_output_manifest$sha256,
    complete_sha256 = published_complete$sha256,
    runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
    runtime_pre_digest = launch_context$runtime_pre_digest,
    runtime_post_digest = launch_context$runtime_pre_digest,
    runtime_pre_post_equal = TRUE,
    final_challenge_request_sha256 = challenge$request$sha256,
    final_challenge_response_sha256 = challenge$response$sha256,
    final_challenge_payload_sha256 = challenge$payload_sha256,
    final_challenge_nonce = challenge$request_value$challenge_nonce,
    final_challenge_issued_utc = challenge$request_value$issued_utc,
    final_challenge_responded_utc = challenge$response_value$responded_utc,
    lease_issued_utc = lease$issued_utc[[1L]],
    lease_expires_utc = lease$expires_utc[[1L]]
  ), attempt_root
)

ack_expected <- list(
  holder_principal = holder_principal,
  run_id = run_id,
  lease_id = lease$lease_id[[1L]],
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
  model_source_sha256 = script_runtime_buffer$sha256,
  validator_source_sha256 = runtime_buffers$validator_source$sha256,
  shared_library_sha256 = library_runtime_buffer$sha256,
  analysis_config_sha256 = config_buffer$sha256,
  external_launcher_sha256 = launch_context$bootstrap_launcher_sha256,
  os_guard_sha256 = runtime_buffers$os_guard_source$sha256,
  launcher_pid = as.integer(launch_context$launcher_pid),
  launcher_start_utc = launch_context$launcher_start_utc,
  guard_owner_role = guard_state$owner_role,
  guard_owner_pid = guard_state$owner_pid,
  guard_owner_start_utc = guard_state$owner_start_utc,
  guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
  finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
  finalization_public_key_captured_path = guard_state$finalization_public_key_captured_path,
  finalization_public_key_captured_bytes = guard_state$finalization_public_key_captured_bytes,
  finalization_public_key_captured_sha256 = guard_state$finalization_public_key_captured_sha256,
  guard_pid = guard_state$guard_pid,
  guard_start_utc = guard_state$guard_start_utc,
  guard_nonce = guard_state$guard_nonce,
  guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
  output_manifest_sha256 = published_output_manifest$sha256,
  complete_sha256 = published_complete$sha256,
  validation_report_sha256 = science$value$validation_report_sha256,
  validation_request_sha256 = finalization_request_buffer$sha256,
  snapshot_tree_sha256 = snapshot$tree_digest,
  runtime_pre_digest = launch_context$runtime_pre_digest,
  runtime_post_digest = runtime_post$value$runtime_post_digest,
  runtime_pre_post_equal = TRUE,
  runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
  model_request_sha256 = child_request$buffer$sha256,
  model_request_nonce = child_request$request$request_nonce,
  model_identity_request_sha256 = model_identity_request_buffer$sha256,
  validator_request_sha256 = validator_request_buffer$sha256,
  validator_request_nonce = validator_request_bundle$request$request_nonce,
  validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  science_complete_sha256 = science$buffer$sha256,
  control_complete_sha256 = control$buffer$sha256,
  finalization_request_sha256 = finalization_request_buffer$sha256,
  final_challenge_request_sha256 = challenge$request$sha256,
  final_challenge_response_sha256 = challenge$response$sha256,
  final_challenge_payload_sha256 = challenge$payload_sha256,
  final_challenge_nonce = challenge$request_value$challenge_nonce,
  final_challenge_issued_utc = challenge$request_value$issued_utc,
  final_challenge_responded_utc = challenge$response_value$responded_utc,
  lease_issued_utc = lease$issued_utc[[1L]],
  lease_expires_utc = lease$expires_utc[[1L]],
  model_pid = as.integer(model_identity$pid),
  model_start_utc = model_identity$start_utc,
  validator_pid = as.integer(lease$validator_pid[[1L]]),
  validator_start_utc = lease$validator_start_utc[[1L]],
  guard_mutex_file_id = guard_state$mutex_file_id,
  named_mutex_unavailable = TRUE,
  guard_mutex_held = TRUE,
  e2_gate_status = cfg$gate$machine_status,
  candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)
ack <- v2r4_verify_signed_final_ack(
  final_ack_path, finalization_public_key_path,
  launch_context$finalization_public_key_captured_sha256, cfg, ack_expected,
  final_wait
)
v2r4_write_consumer_receipt(
  model_consumer_receipt_path, cfg, "MODEL", ack,
  ack_expected, attempt_root
)
message(sprintf("PASS %s %s", cfg$model_version,
                if (identical(cli$execution_mode, "SYNTHETIC")) "SYNTHETIC_ISOLATED_ONLY" else cfg$status))
