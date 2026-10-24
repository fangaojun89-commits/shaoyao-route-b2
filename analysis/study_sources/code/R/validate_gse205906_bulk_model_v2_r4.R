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
script_runtime_buffer <- v2r4_bootstrap_capture(script_path, "VALIDATOR_SOURCE_RUNTIME")
library_runtime_buffer <- v2r4_bootstrap_capture(library_path, "SHARED_LIBRARY_RUNTIME")
library_expressions <- parse(text = enc2utf8(rawToChar(library_runtime_buffer$raw)),
                             keep.source = FALSE)
for (library_expression in library_expressions) {
  eval(library_expression, envir = .GlobalEnv)
}
rm(v2r4_bootstrap_capture, library_expressions, library_expression)

child_request <- v2r4_capture_child_request_args(
  commandArgs(trailingOnly = TRUE), "VALIDATOR"
)
cli <- child_request$cli
project_root <- v2r4_project_root(script_path)
default_config <- file.path(project_root, "config", "GSE205906_bulk_model_v2_r4.json")
config_path <- cli$config

required_control_args <- c(
  "authority_manifest", "authority_signature", "authority_public_key", "snapshot_root",
  "lease_manifest", "lease_signature", "lease_public_key", "lease_request",
  "lease_identity_request", "peer_lease_identity_request", "state_root",
  "run_id", "holder_principal", "output_authority_manifest", "output_authority_signature",
  "output_authority_public_key", "output_dir", "report", "finalization_request",
  "model_request", "runtime_post_receipt", "science_complete", "control_complete",
  "final_challenge_request", "final_challenge_response", "final_ack",
  "finalization_public_key", "model_consumer_receipt", "validator_consumer_receipt"
)
missing_control <- required_control_args[vapply(required_control_args, function(k) {
  is.null(cli[[k]]) || !nzchar(cli[[k]])
}, logical(1))]
if (length(missing_control)) {
  v2r4_stop("EXTERNAL_SOURCE_LEASE_AND_OUTPUT_AUTHORITY_REQUIRED_BEFORE_ANY_OUTPUT_READ:%s",
            paste(missing_control, collapse = ","))
}
for (k in c("authority_manifest", "authority_signature", "authority_public_key")) {
  if (!file.exists(cli[[k]])) {
    v2r4_stop("EXTERNAL_SOURCE_LEASE_AND_OUTPUT_AUTHORITY_REQUIRED_BEFORE_ANY_OUTPUT_READ:%s_MISSING", k)
  }
}

v2r4_require_packages(c("data.table", "digest", "jsonlite", "openssl", "yaml", "edgeR",
                        "limma", "DESeq2"))
config_path <- v2r4_norm_existing(config_path)
v2r4_assert(identical(config_path, v2r4_norm_existing(default_config)),
            "CONFIG_CANONICAL_PATH_NOT_LOCKED")
config_buffer <- v2r4_capture_file(config_path, "ANALYSIS_CONFIG_RUNTIME", project_root)
cfg <- v2r4_load_config_buffer(config_buffer)
analysis_role <- cfg$gate$analysis_role
launch_context <- v2r4_require_external_launch_context(cfg, "VALIDATOR")
v2r4_validate_child_request_bindings(child_request, cfg, launch_context, "VALIDATOR")

v2r4_assert(cli$execution_mode %in% c("SYNTHETIC", "PRODUCTION") &&
              identical(cli$execution_mode, launch_context$mode),
            "SIGNED_EXECUTION_MODE_INVALID")
required_execution <- c("attempt_root", "source_key_sha256", "lease_key_sha256",
                        "output_key_sha256")
missing_execution <- required_execution[vapply(required_execution, function(k) {
  is.null(cli[[k]]) || !nzchar(cli[[k]])
}, logical(1))]
v2r4_assert(!length(missing_execution), "EXECUTION_CONTROL_ARGUMENTS_MISSING:%s",
            paste(missing_execution, collapse = ","))
attempt_root <- v2r4_norm_existing(cli$attempt_root)
v2r4_require_outside(attempt_root, project_root, "EXTERNAL_ATTEMPT_ROOT")
snapshot_root <- v2r4_norm_existing(cli$snapshot_root)
state_root <- v2r4_norm_existing(cli$state_root)
output_dir <- v2r4_norm_target(cli$output_dir)
report_path <- v2r4_norm_target(cli$report)
finalization_request_path <- v2r4_norm_target(cli$finalization_request)
model_request_path <- v2r4_norm_existing(cli$model_request)
lease_request_path <- v2r4_norm_target(cli$lease_request)
lease_identity_request_path <- v2r4_norm_target(cli$lease_identity_request)
peer_lease_identity_request_path <- v2r4_norm_target(cli$peer_lease_identity_request)
runtime_post_receipt_path <- v2r4_norm_target(cli$runtime_post_receipt)
science_complete_path <- v2r4_norm_target(cli$science_complete)
control_complete_path <- v2r4_norm_target(cli$control_complete)
pre_final_challenge_ready_path <- v2r4_norm_target(file.path(
  dirname(control_complete_path), "validator_pre_final_challenge_ready.json"
))
final_ack_path <- v2r4_norm_target(cli$final_ack)
finalization_public_key_path <- v2r4_norm_existing(cli$finalization_public_key)
output_public_key_path <- v2r4_norm_existing(cli$output_authority_public_key)
model_consumer_receipt_path <- v2r4_norm_target(cli$model_consumer_receipt)
validator_consumer_receipt_path <- v2r4_norm_target(cli$validator_consumer_receipt)
runtime_final_receipt_path <- v2r4_norm_target(launch_context$runtime_final_receipt_path)
guard_ready_receipt_path <- v2r4_norm_existing(launch_context$guard_ready_receipt_path)
guard_stop_request_path <- v2r4_norm_target(launch_context$guard_stop_request_path)
guard_stopped_receipt_path <- v2r4_norm_target(launch_context$guard_stopped_receipt_path)
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
  c(v2r4_norm_target(cli$lease_manifest), "LEASE_MANIFEST"),
  c(v2r4_norm_target(cli$lease_signature), "LEASE_SIGNATURE"),
  c(v2r4_norm_target(cli$output_authority_manifest), "OUTPUT_AUTHORITY"),
  c(v2r4_norm_target(cli$output_authority_signature), "OUTPUT_SIGNATURE"),
  c(output_public_key_path, "OUTPUT_PUBLIC_KEY"),
  c(state_root, "STATE_ROOT"), c(output_dir, "OUTPUT_ROOT"), c(report_path, "REPORT"),
  c(lease_request_path, "COMBINED_LEASE_REQUEST"),
  c(lease_identity_request_path, "VALIDATOR_LEASE_IDENTITY_REQUEST"),
  c(peer_lease_identity_request_path, "MODEL_LEASE_IDENTITY_REQUEST"),
  c(finalization_request_path, "FINALIZATION_REQUEST"), c(model_request_path, "MODEL_REQUEST"),
  c(runtime_post_receipt_path, "RUNTIME_POST_RECEIPT"),
  c(science_complete_path, "SCIENCE_COMPLETE"), c(control_complete_path, "CONTROL_COMPLETE"),
  c(pre_final_challenge_ready_path, "VALIDATOR_PRE_FINAL_CHALLENGE_READY"),
  c(final_ack_path, "FINAL_ACK"), c(finalization_public_key_path, "FINALIZATION_PUBLIC_KEY"),
  c(model_consumer_receipt_path, "MODEL_CONSUMER_RECEIPT"),
  c(validator_consumer_receipt_path, "VALIDATOR_CONSUMER_RECEIPT"),
  c(runtime_final_receipt_path, "RUNTIME_FINAL_RECEIPT"),
  c(guard_ready_receipt_path, "GUARD_READY_RECEIPT"),
  c(guard_stop_request_path, "GUARD_STOP_REQUEST"),
  c(guard_stopped_receipt_path, "GUARD_STOPPED_RECEIPT")
)) v2r4_require_within(pair[[1L]], attempt_root, pair[[2L]])
v2r4_assert(
  identical(lease_identity_request_path,
            v2r4_norm_target(launch_context$validator_lease_identity_request_path)) &&
    identical(peer_lease_identity_request_path,
              v2r4_norm_target(launch_context$model_lease_identity_request_path)) &&
    identical(lease_request_path,
              v2r4_norm_target(launch_context$combined_lease_request_path)),
  "VALIDATOR_LEASE_REQUEST_PATH_CONTEXT_BINDING_INVALID"
)
source_key_sha256 <- tolower(cli$source_key_sha256)
lease_key_sha256 <- tolower(cli$lease_key_sha256)
output_key_sha256 <- tolower(cli$output_key_sha256)
expected_source_decision <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  cfg$authority$synthetic_source_decision
} else cfg$authority$production_source_decision
expected_output_decision <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  "SYNTHETIC_GSE205906_BULK_MODEL_V2_R4_OUTPUT_MARKERS"
} else cfg$output_authority$decision
if (identical(cli$execution_mode, "PRODUCTION")) {
  v2r4_assert(identical(source_key_sha256,
                        cfg$authority$production_source_public_key_sha256) &&
                identical(lease_key_sha256, cfg$lease$production_public_key_sha256) &&
                identical(output_key_sha256,
                          cfg$output_authority$production_public_key_sha256),
              "PRODUCTION_SOURCE_LEASE_OR_OUTPUT_KEY_PIN_INVALID")
}
v2r4_assert(!file.exists(report_path), "VALIDATION_REPORT_ALREADY_EXISTS")
v2r4_assert(!file.exists(lease_identity_request_path),
            "VALIDATOR_LEASE_IDENTITY_REQUEST_ALREADY_EXISTS")
v2r4_assert(!file.exists(finalization_request_path), "FINALIZATION_REQUEST_ALREADY_EXISTS")
v2r4_assert(!file.exists(science_complete_path) && !file.exists(control_complete_path) &&
              !file.exists(pre_final_challenge_ready_path) &&
              !file.exists(final_ack_path) && !file.exists(validator_consumer_receipt_path),
            "VALIDATOR_FINAL_CONTROL_TARGET_ALREADY_EXISTS")
v2r4_assert(identical(source_key_sha256, launch_context$source_snapshot_public_key_sha256) &&
              identical(lease_key_sha256, launch_context$lease_public_key_sha256) &&
              identical(output_key_sha256,
                        launch_context$output_public_key_captured_sha256) &&
              identical(output_public_key_path,
                        v2r4_norm_existing(launch_context$output_public_key_captured_path)) &&
              {
                captured_output_key <- v2r4_capture_file(
                  output_public_key_path, "VALIDATOR_OUTPUT_PUBLIC_KEY_CAPTURED_COPY",
                  attempt_root
                )
                identical(as.numeric(captured_output_key$bytes),
                          as.numeric(launch_context$output_public_key_captured_bytes)) &&
                  identical(captured_output_key$sha256,
                            launch_context$output_public_key_captured_sha256) &&
                  identical(captured_output_key$sha256,
                            launch_context$output_public_key_sha256)
              } &&
              identical(finalization_public_key_path,
                        v2r4_norm_existing(
                          launch_context$finalization_public_key_captured_path
                        )) &&
              {
                captured_finalization_key <- v2r4_capture_file(
                  finalization_public_key_path,
                  "VALIDATOR_FINALIZATION_PUBLIC_KEY_CAPTURED_COPY",
                  attempt_root
                )
                identical(as.numeric(captured_finalization_key$bytes),
                          as.numeric(launch_context$finalization_public_key_captured_bytes)) &&
                  identical(captured_finalization_key$sha256,
                            launch_context$finalization_public_key_captured_sha256) &&
                  identical(captured_finalization_key$sha256,
                            launch_context$finalization_public_key_sha256)
              } &&
              identical(v2r4_capture_file(cli$authority_manifest,
                                          "VALIDATOR_SOURCE_MANIFEST_LAUNCH_BINDING")$sha256,
                        launch_context$source_snapshot_manifest_sha256),
            "VALIDATOR_EXTERNAL_LAUNCH_SOURCE_LEASE_FINALIZATION_BINDING_INVALID")

checks <- list()
add_check <- function(id, detail) {
  checks[[length(checks) + 1L]] <<- data.table::data.table(
    check_id = id, status = "PASS", detail = as.character(detail),
    analysis_role = analysis_role, candidate_gate_eligible = FALSE,
    can_generate_or_upgrade_candidate = FALSE, formal_sequential_inference = FALSE
  )
}

source_bundle <- v2r4_verify_detached_bundle(
  cli$authority_manifest, cli$authority_signature, cli$authority_public_key,
  source_key_sha256, "SOURCE_AUTHORITY"
)
source_auth <- v2r4_parse_authority_buffer(source_bundle$manifest, cfg, expected_source_decision)
expected_roles <- c(unlist(cfg$authority$required_dynamic_roles, use.names = FALSE),
                    unlist(cfg$authority$required_data_roles, use.names = FALSE))
v2r4_assert(setequal(source_auth$role, expected_roles), "SOURCE_AUTHORITY_ROLE_SET_NOT_EXACT")
add_check("SOURCE_AUTHORITY_SAME_BUFFER_SIGNATURE_PARSE", source_bundle$manifest$sha256)

guard_path <- v2r4_norm_existing(file.path(
  project_root, cfg$authority$dynamic_runtime_paths$os_guard_source
))
runtime_buffers <- list(
  model_source = v2r4_capture_file(file.path(
    project_root, cfg$authority$dynamic_runtime_paths$model_source
  ), "MODEL_SOURCE_RUNTIME", project_root),
  validator_source = script_runtime_buffer,
  shared_library_source = library_runtime_buffer,
  analysis_config = config_buffer,
  os_guard_source = v2r4_capture_file(guard_path, "OS_GUARD_RUNTIME", project_root)
)
runtime_bindings <- v2r4_bind_runtime_roles(source_auth, cfg, project_root, runtime_buffers)
v2r4_assert(identical(script_runtime_buffer$sha256,
                        launch_context$validator_source_sha256) &&
              identical(library_runtime_buffer$sha256,
                        launch_context$library_source_sha256) &&
              identical(config_buffer$sha256, launch_context$config_sha256) &&
              identical(runtime_buffers$model_source$sha256,
                        launch_context$model_source_sha256) &&
              identical(runtime_buffers$os_guard_source$sha256,
                        launch_context$os_guard_sha256),
            "VALIDATOR_RUNTIME_BYTES_NOT_BOUND_TO_EXTERNAL_LAUNCH")
add_check("DYNAMIC_RUNTIME_ROLE_BINDING", sprintf("%d exact runtime roles", nrow(runtime_bindings)))

validator_identity <- v2r4_process_identity(runtime_buffers$os_guard_source, Sys.getpid())
v2r4_assert_child_stop_key_absent()
guard_state <- v2r4_attach_existing_guard(
  runtime_buffers$os_guard_source, state_root, as.character(cli$run_id), cfg,
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
  expected_ready_path = guard_ready_receipt_path,
  expected_stop_path = guard_stop_request_path,
  expected_stopped_path = guard_stopped_receipt_path
)
guard_state$initial_challenge <- v2r4_guard_challenge(
  guard_state, runtime_buffers$os_guard_source, cfg, Sys.getpid()
)
authorized_snapshot_tree_digest <- v2r4_authorized_snapshot_tree_digest(source_auth)
runtime_role_digest <- v2r4_runtime_role_digest(runtime_bindings)
validator_identity_request_buffer <- v2r4_write_lease_identity_request(
  lease_identity_request_path, cfg, "VALIDATOR", as.character(cli$run_id),
  as.character(cli$holder_principal), authorized_snapshot_tree_digest, guard_state,
  runtime_role_digest, launch_context, validator_identity, child_request, attempt_root
)
lease_wait <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  cfg$execution_contract$lease_wait_seconds_synthetic
} else {
  cfg$execution_contract$lease_wait_seconds_production
}
v2r4_wait_for_files(
  c(peer_lease_identity_request_path, cli$lease_manifest, cli$lease_signature,
    cli$lease_public_key), lease_wait, "VALIDATOR_EXTERNAL_SIGNED_LEASE"
)
model_request_buffer_prelease <- v2r4_capture_file(
  model_request_path, "VALIDATOR_PRELEASE_MODEL_CHILD_REQUEST", attempt_root
)
model_request_bundle_prelease <- v2r4_parse_child_request_buffer(
  model_request_buffer_prelease, model_request_buffer_prelease$sha256, "MODEL"
)
v2r4_validate_child_request_bindings(
  model_request_bundle_prelease, cfg, launch_context, "MODEL"
)
model_identity_request_buffer <- v2r4_capture_file(
  peer_lease_identity_request_path, "VALIDATOR_CAPTURED_MODEL_LEASE_IDENTITY_REQUEST",
  attempt_root
)
model_identity_request <- v2r4_validate_lease_identity_request(
  model_identity_request_buffer, cfg, "MODEL",
  list(
    run_id = as.character(cli$run_id),
    holder_principal = as.character(cli$holder_principal),
    child_request_sha256 = model_request_buffer_prelease$sha256,
    child_request_nonce = model_request_bundle_prelease$request$request_nonce,
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
guard_state$initial_challenge <- list(
  request = structure(list(sha256 = model_identity_request$guard_challenge_request_sha256),
                      class = "v2r4_lease_bound_identity"),
  response = structure(list(sha256 = model_identity_request$guard_challenge_response_sha256),
                       class = "v2r4_lease_bound_identity")
)
lease <- v2r4_validate_lease_buffer(
  lease_bundle$manifest, cfg, authorized_snapshot_tree_digest, as.character(cli$run_id),
  as.character(cli$holder_principal), runtime_buffers$os_guard_source, guard_state,
  launch_context, runtime_role_digest, state_root,
  expected_model_pid = model_identity_request$process_pid,
  expected_model_start_utc = model_identity_request$process_start_utc,
  expected_model_request_sha256 = model_request_buffer_prelease$sha256,
  expected_model_request_nonce = model_request_bundle_prelease$request$request_nonce,
  expected_model_identity_request_sha256 = model_identity_request_buffer$sha256,
  expected_validator_pid = Sys.getpid(),
  expected_validator_start_utc = validator_identity$start_utc,
  expected_validator_request_sha256 = child_request$buffer$sha256,
  expected_validator_request_nonce = child_request$request$request_nonce,
  expected_validator_identity_request_sha256 = validator_identity_request_buffer$sha256
)
validator_guard_challenge <- v2r4_guard_challenge(
  guard_state, runtime_buffers$os_guard_source, cfg, Sys.getpid()
)
add_check("SIGNED_LIVE_LEASE_AND_OS_MUTEX", lease$lease_id[[1L]])

# The validator remains data-blind until the model's externally attested outputs exist.
output_wait <- cfg$execution_contract$signed_final_ack_wait_seconds
v2r4_assert(identical(output_wait, 900L),
            "VALIDATOR_OUTPUT_WAIT_MUST_BE_CAPTURED_EXACT_900")
v2r4_wait_for_files(
  c(file.path(output_dir, "OUTPUT_MANIFEST.tsv"), file.path(output_dir, "COMPLETE.tsv"),
    cli$output_authority_manifest, cli$output_authority_signature,
    cli$output_authority_public_key),
  output_wait, "VALIDATOR_MODEL_OUTPUT_AND_AUTHORITY"
)
output_dir <- v2r4_norm_existing(output_dir)

snapshot <- v2r4_capture_snapshot(source_auth, snapshot_root, defer_roles = "assembled_counts_rds")
producer_source_provenance <- v2r4_validate_producer_source_buffers(snapshot, cfg)
add_check("PRODUCER_SOURCE_PROVENANCE",
          "4 signed producer sources matched before assembled-count capture")
ledger <- v2r4_validate_ledger(
  v2r4_read_tsv_raw(v2r4_object(snapshot, "sample_ledger"), "SAMPLE_LEDGER"), cfg
)
add_check(
  "LEDGER_INDEPENDENT_LIBRARY_UNIT",
  "nine independent rat libraries, treated as one library per animal based on repository metadata."
)
v2r4_validate_upstream_gate(
  v2r4_read_tsv_raw(v2r4_object(snapshot, "upstream_e1_gate"), "UPSTREAM_E1_GATE"), cfg
)
add_check("E1_EMPTY_RECOMPUTATION", "E1 eligible_input_n=0")
reference_manifest <- v2r4_read_tsv_raw(
  v2r4_object(snapshot, "reference_manifest"), "REFERENCE_MANIFEST"
)
v2r4_validate_reference_manifest(reference_manifest, cfg)
signed_fasta <- v2r4_role_row(source_auth, "reference_fasta_identity")
signed_gtf <- v2r4_role_row(source_auth, "reference_gtf_identity")
manifest_fasta <- reference_manifest[base::which(reference_manifest[["role"]] == "reference_fasta")]
manifest_gtf <- reference_manifest[base::which(reference_manifest[["role"]] == "gene_annotation_gtf")]
v2r4_assert(signed_fasta$materialization[[1L]] == "IDENTITY_ONLY" &&
              signed_gtf$materialization[[1L]] == "IDENTITY_ONLY" &&
              signed_fasta$sha256[[1L]] == cfg$input_snapshot$reference_fasta_sha256 &&
              signed_gtf$sha256[[1L]] == cfg$input_snapshot$gtf_sha256 &&
              as.numeric(signed_fasta$bytes[[1L]]) == as.numeric(manifest_fasta$bytes[[1L]]) &&
              as.numeric(signed_gtf$bytes[[1L]]) == as.numeric(manifest_gtf$bytes[[1L]]),
            "SIGNED_REFERENCE_IDENTITY_ROWS_INVALID")
add_check("REFERENCE_ASSEMBLY_TRUST_ROOT", "mRatBN7.2 RefSeq identities independently matched")
v2r4_validate_strand(
  v2r4_read_tsv_raw(v2r4_object(snapshot, "strand_decision"), "STRAND_DECISION"),
  v2r4_read_tsv_raw(v2r4_object(snapshot, "strand_assignment_by_sample"), "STRAND_ASSIGNMENT"), cfg
)
add_check("STRAND_RECOMPUTATION", "UNSTRANDED selected without labels")
cohort_species_buffer <- v2r4_object(snapshot, "cohort_species_gate")
sample_species_buffer <- v2r4_object(snapshot, "sample_species_table")
cohort_species <- v2r4_read_tsv_raw(cohort_species_buffer, "SPECIES_COHORT_GATE")
sample_species <- v2r4_read_tsv_raw(sample_species_buffer, "SPECIES_SAMPLE_GATE")
v2r4_validate_species(cohort_species, sample_species, cfg)
conditional_inclusion <- v2r4_release_conditional_inclusion(
  ledger, cohort_species, sample_species,
  cohort_species_buffer$sha256, sample_species_buffer$sha256, cfg
)
add_check("SPECIES_RECOMPUTATION_BEFORE_COUNTS_CAPTURE", "9/9 Rattus norvegicus screen PASS")
add_check("CONDITIONAL_INCLUSION_RELEASE",
          "9 CONDITIONAL libraries released only after exact sample and cohort species rules PASS")
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
v2r4_assert(is.list(assembled) && all(required_assembled %in% names(assembled)) &&
              identical(assembled$dataset_id, cfg$dataset_id), "ASSEMBLED_RDS_SCHEMA_INVALID")
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
annotation_contract <- v2r4_validate_annotation_gene(
  assembled$gene_annotation, rownames(counts), cfg$input_snapshot$minimum_nonblank_symbol_coverage
)
annotation <- annotation_contract$annotation
annotation_metrics_expected <- annotation_contract$metrics
assembly_manifest <- v2r4_read_tsv_raw(
  v2r4_object(snapshot, "assembly_output_manifest"), "ASSEMBLY_OUTPUT_MANIFEST"
)
rds_row <- assembly_manifest[base::which(assembly_manifest[["filename"]] == "gene_counts_primary.rds")]
v2r4_assert(nrow(rds_row) == 1L && as.numeric(rds_row$bytes[[1L]]) == assembled_buffer$bytes &&
              tolower(rds_row$md5[[1L]]) == assembled_buffer$md5,
            "ASSEMBLY_MANIFEST_RDS_BINDING_INVALID")
assembly_complete <- v2r4_read_tsv_raw(v2r4_object(snapshot, "assembly_complete"), "ASSEMBLY_COMPLETE")
v2r4_assert(nrow(assembly_complete) == 1L && assembly_complete$dataset_id[[1L]] == cfg$dataset_id &&
              assembly_complete$n_libraries[[1L]] == 9L && assembly_complete$status[[1L]] == "PASS" &&
              assembly_complete$selected_primary_strand_mode[[1L]] == 0L,
            "ASSEMBLY_COMPLETE_SEMANTICS_INVALID")
add_check("COUNTS_AND_GENE_INTERFACE", "raw count domain and frozen gene symbol column independently validated")

legacy_contrasts <- v2r4_read_yaml_raw(v2r4_object(snapshot, "contrasts_config"), "CONTRASTS_CONFIG")
legacy_thresholds <- v2r4_read_yaml_raw(v2r4_object(snapshot, "thresholds_config"), "THRESHOLDS_CONFIG")
gse_cfg <- legacy_contrasts$datasets$GSE205906
v2r4_assert(identical(gse_cfg$model, "~ 0 + group") &&
              identical(gse_cfg$contrasts$local_disease$expression, "SNI - sham") &&
              identical(gse_cfg$contrasts$intervention$expression, "PLM - SNI") &&
              as.numeric(legacy_thresholds$expression_filtering$min_cpm) == 1 &&
              as.numeric(legacy_thresholds$false_discovery_rate$primary_gene_set_q) == 0.05,
            "INTERNAL_SAP_DESIGN_CONTRAST_OR_THRESHOLD_DRIFT")
design_contract <- v2r4_design(ledger, cfg)
design <- design_contract$design
contrast_matrix <- design_contract$contrast
group <- design_contract$group
add_check("DESIGN_AND_CONTRAST_RECOMPUTATION", "rank 3; SNI-sham and PLM-SNI exact")

output_bundle <- v2r4_verify_detached_bundle(
  cli$output_authority_manifest, cli$output_authority_signature,
  cli$output_authority_public_key, output_key_sha256, "OUTPUT_AUTHORITY"
)
output_auth <- v2r4_parse_output_authority(
  output_bundle$manifest, cfg, expected_output_decision, cli$execution_mode,
  as.character(cli$run_id)
)
authorized_output <- v2r4_capture_authorized_outputs(
  output_auth, output_dir, unlist(cfg$required_output_payloads, use.names = FALSE)
)
add_check("OUTPUT_AUTHORITY_EXACT_MARKERS_SAME_BUFFER", output_bundle$manifest$sha256)

payload_tsv <- function(relative_path, label = basename(relative_path),
                        character_columns = character()) {
  v2r4_read_tsv_raw(authorized_output$payloads[[relative_path]], label,
                    gzip = endsWith(relative_path, ".gz"),
                    character_columns = character_columns)
}

assert_keyed_equal <- function(observed, expected, keys, fields, label, tolerance = 1e-10) {
  observed <- data.table::copy(data.table::as.data.table(observed))
  expected <- data.table::copy(data.table::as.data.table(expected))
  v2r4_assert(all(c(keys, fields) %in% names(observed)) && all(c(keys, fields) %in% names(expected)),
              "%s_COMPARE_SCHEMA_INVALID", label)
  data.table::setorderv(observed, keys)
  data.table::setorderv(expected, keys)
  v2r4_assert(identical(observed[, ..keys], expected[, ..keys]), "%s_KEY_SET_MISMATCH", label)
  for (field in fields) {
    got <- observed[[field]]
    want <- expected[[field]]
    if (is.numeric(got) || is.integer(got) || is.numeric(want) || is.integer(want)) {
      got <- as.numeric(got)
      want <- as.numeric(want)
      same <- (is.na(got) & is.na(want)) |
        (is.finite(got) & is.finite(want) & abs(got - want) <= tolerance)
    } else {
      same <- (is.na(got) & is.na(want)) | as.character(got) == as.character(want)
    }
    v2r4_assert(all(same), "%s_FIELD_MISMATCH:%s", label, field)
  }
  invisible(TRUE)
}

assert_complete_numeric_family <- function(observed, expected, keys, exact_fields,
                                           numeric_fields, label, cfg) {
  observed <- data.table::copy(data.table::as.data.table(observed))
  expected <- data.table::copy(data.table::as.data.table(expected))
  required <- c(keys, exact_fields, numeric_fields)
  v2r4_assert(all(required %in% names(observed)) && all(required %in% names(expected)),
              "%s_COMPLETE_NUMERIC_SCHEMA_INVALID", label)
  v2r4_assert(!anyDuplicated(observed[, ..keys]) && !anyDuplicated(expected[, ..keys]),
              "%s_DUPLICATE_TESTED_FAMILY_KEY", label)
  data.table::setorderv(observed, keys)
  data.table::setorderv(expected, keys)
  for (key in keys) {
    v2r4_assert(identical(as.character(observed[[key]]), as.character(expected[[key]])),
                "%s_COMPLETE_TESTED_FAMILY_KEY_MISMATCH:%s", label, key)
  }
  for (field in exact_fields) {
    same <- (is.na(observed[[field]]) & is.na(expected[[field]])) |
      (!is.na(observed[[field]]) & !is.na(expected[[field]]) &
         as.character(observed[[field]]) == as.character(expected[[field]]))
    v2r4_assert(all(same), "%s_EXACT_PROVENANCE_FIELD_MISMATCH:%s", label, field)
  }
  contract <- cfg$analysis$independent_numeric_revalidation
  abs_tol <- as.numeric(contract$absolute_tolerance)
  rel_tol <- as.numeric(contract$relative_tolerance)
  for (field in numeric_fields) {
    got <- as.numeric(observed[[field]])
    want <- as.numeric(expected[[field]])
    v2r4_assert(identical(is.na(got), is.na(want)),
                "%s_NA_PATTERN_MISMATCH:%s", label, field)
    finite <- !is.na(want)
    v2r4_assert(all(is.finite(got[finite])) && all(is.finite(want[finite])),
                "%s_NONFINITE_NUMERIC_DOMAIN:%s", label, field)
    if (any(finite)) {
      bound <- abs_tol + rel_tol * pmax(abs(got[finite]), abs(want[finite]))
      v2r4_assert(all(abs(got[finite] - want[finite]) <= bound),
                  "%s_NUMERIC_MISMATCH:%s", label, field)
    }
  }
  invisible(TRUE)
}

complete <- v2r4_read_tsv_raw(authorized_output$markers$complete, "AUTHORIZED_COMPLETE")
required_complete <- c(
  "dataset_id", "model_version", "implementation_status", "computation_status",
  "run_id", "model_pid", "snapshot_tree_sha256", "output_manifest_sha256",
  "gate_status", "analysis_role", "candidate_gate_eligible",
  "can_generate_or_upgrade_candidate", "formal_sequential_inference"
)
v2r4_assert(nrow(complete) == 1L && all(required_complete %in% names(complete)) &&
              complete$dataset_id[[1L]] == cfg$dataset_id &&
              complete$model_version[[1L]] == cfg$model_version &&
              complete$implementation_status[[1L]] == cfg$status &&
              complete$run_id[[1L]] == cli$run_id &&
              as.integer(complete$model_pid[[1L]]) == as.integer(lease$model_pid[[1L]]) &&
              complete$snapshot_tree_sha256[[1L]] == snapshot$tree_digest &&
              complete$output_manifest_sha256[[1L]] == authorized_output$markers$output_manifest$sha256 &&
              complete$gate_status[[1L]] == cfg$gate$machine_status,
            "COMPLETE_TRUST_OR_GATE_SEMANTICS_INVALID")
v2r4_validate_uniform_flags(complete, cfg, "COMPLETE")
add_check("COMPLETE_TRUST_BINDING", "source tree, live model and exact output manifest bound")

run <- payload_tsv(
  "audit/run_manifest.tsv", "RUN_MANIFEST",
  character_columns = c("launcher_start_utc", "guard_owner_start_utc",
                        "model_start_utc", "validator_start_utc", "guard_start_utc")
)
required_run <- c(
  "dataset_id", "model_version", "implementation_status", "execution_mode",
  "source_authority_manifest_sha256",
  "source_authority_public_key_sha256", "source_authority_decision",
  "external_source_authority_payload_sha256", "external_launch_context_sha256",
  "external_launcher_sha256", "runtime_pre_digest", "runtime_semantic_digest",
  "snapshot_tree_sha256",
  "lease_manifest_sha256", "lease_id", "holder_principal", "run_id", "launcher_pid",
  "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
  "guard_owner_start_utc", "guard_stop_hmac_key_sha256",
  "finalization_public_key_sha256", "finalization_public_key_captured_path",
  "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256",
  "model_pid",
  "model_start_utc", "model_request_sha256", "model_request_nonce",
  "model_identity_request_sha256", "validator_pid", "validator_start_utc",
  "validator_request_sha256", "validator_request_nonce",
  "validator_identity_request_sha256", "mutex_relative_path", "guard_pid",
  "guard_start_utc", "guard_nonce",
  "guard_public_key_sha256", "guard_ready_receipt_sha256",
  "producer_source_provenance_status", "producer_source_digest",
  "statistical_unit", "independent_rat_libraries", "unique_biosamples", "unique_runs",
  "animal_independence_basis", "repository_inferred_animal_n", "group_counts",
  "sample_include_primary_marker", "sample_metadata_status",
  "conditional_inclusion_release_status", "conditional_inclusion_n",
  "cohort_species_gate_sha256", "sample_species_table_sha256",
  "plm_metadata_conflict_policy",
  "design_formula", "design_rank", "coefficient_order", "local_disease_expression",
  "intervention_expression", "primary_method", "primary_filtering", "normalization",
  "random_seed", "orthology_authority", "orthology_single_source_role", "assembly",
  "assembly_accession", "selected_primary_strand_mode", "species_gate_status",
  "gene_id_source_column", "symbol_source_column", "biotype_source_column",
  "nonblank_symbol_coverage", "unique_nonblank_symbol_coverage",
  "fallback_to_gene_id_n", "requested_DESeq2_fit_type", "actual_DESeq2_fit_type_after_DESeq",
  "DESeq2_fit_messages",
  "progeny_source_program_n", "progeny_tested_program_n", "dorothea_source_program_n",
  "dorothea_tested_program_n", "dorothea_confidence_domain", "loao_method", "gate_status",
  "analysis_role", "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
  "formal_sequential_inference"
)
v2r4_assert(nrow(run) == 1L && identical(names(run), required_run),
            "RUN_MANIFEST_EXACT_SCHEMA_INVALID")
v2r4_assert_utc_fixed6(run$model_start_utc[[1L]], "RUN_MANIFEST_MODEL_START")
expected_execution_mode <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  "SYNTHETIC_ISOLATED_ONLY"
} else {
  "PRODUCTION_AUTHORITY_GATED"
}
observed_group_counts <- table(group)
expected_group_count_text <- paste(names(observed_group_counts),
                                   as.integer(observed_group_counts),
                                   sep = "=", collapse = ";")
v2r4_assert(run$dataset_id[[1L]] == cfg$dataset_id && run$model_version[[1L]] == cfg$model_version &&
              run$implementation_status[[1L]] == cfg$status &&
              run$execution_mode[[1L]] == expected_execution_mode &&
              run$source_authority_manifest_sha256[[1L]] == source_bundle$manifest$sha256 &&
              run$source_authority_public_key_sha256[[1L]] == source_bundle$public_key$sha256 &&
              run$source_authority_decision[[1L]] == expected_source_decision &&
              run$external_source_authority_payload_sha256[[1L]] ==
                launch_context$source_authorization_payload_sha256 &&
              run$external_launch_context_sha256[[1L]] == launch_context$context_sha256 &&
              run$external_launcher_sha256[[1L]] == launch_context$bootstrap_launcher_sha256 &&
              run$runtime_pre_digest[[1L]] == launch_context$runtime_pre_digest &&
              run$runtime_semantic_digest[[1L]] == launch_context$runtime_semantic_digest &&
              run$snapshot_tree_sha256[[1L]] == snapshot$tree_digest &&
              run$lease_manifest_sha256[[1L]] == lease_bundle$manifest$sha256 &&
              run$lease_id[[1L]] == lease$lease_id[[1L]] && run$holder_principal[[1L]] == cli$holder_principal &&
              run$run_id[[1L]] == cli$run_id && as.integer(run$model_pid[[1L]]) == as.integer(lease$model_pid[[1L]]) &&
              identical(run$model_start_utc[[1L]], lease$model_start_utc[[1L]]) &&
              as.integer(run$validator_pid[[1L]]) == as.integer(lease$validator_pid[[1L]]) &&
              identical(run$validator_start_utc[[1L]], lease$validator_start_utc[[1L]]) &&
              run$model_request_sha256[[1L]] == lease$model_request_sha256[[1L]] &&
              run$model_request_nonce[[1L]] == lease$model_request_nonce[[1L]] &&
              run$model_identity_request_sha256[[1L]] ==
                lease$model_identity_request_sha256[[1L]] &&
              run$validator_request_sha256[[1L]] == lease$validator_request_sha256[[1L]] &&
              run$validator_request_nonce[[1L]] == lease$validator_request_nonce[[1L]] &&
              run$validator_identity_request_sha256[[1L]] ==
                lease$validator_identity_request_sha256[[1L]] &&
              run$mutex_relative_path[[1L]] == cfg$lease$mutex_relative_path &&
              as.integer(run$guard_pid[[1L]]) == guard_state$guard_pid &&
              identical(run$guard_start_utc[[1L]], guard_state$guard_start_utc) &&
              run$guard_nonce[[1L]] == guard_state$guard_nonce &&
              run$guard_public_key_sha256[[1L]] == guard_state$guard_public_key_sha256 &&
              run$guard_ready_receipt_sha256[[1L]] == guard_state$ready_receipt_sha256 &&
              run$finalization_public_key_sha256[[1L]] ==
                guard_state$finalization_public_key_sha256 &&
              v2r4_norm_existing(run$finalization_public_key_captured_path[[1L]]) ==
                guard_state$finalization_public_key_captured_path &&
              as.integer(run$finalization_public_key_captured_bytes[[1L]]) ==
                as.integer(guard_state$finalization_public_key_captured_bytes) &&
              run$finalization_public_key_captured_sha256[[1L]] ==
                 guard_state$finalization_public_key_captured_sha256,
            "RUN_MANIFEST_EXTERNAL_TRUST_BINDING_INVALID")
v2r4_assert(run$producer_source_provenance_status[[1L]] ==
              producer_source_provenance$status &&
              run$producer_source_digest[[1L]] == producer_source_provenance$digest,
            "RUN_MANIFEST_PRODUCER_SOURCE_PROVENANCE_INVALID")
v2r4_assert(run$statistical_unit[[1L]] == "independent_rat_library" &&
              as.integer(run$independent_rat_libraries[[1L]]) == nrow(ledger) &&
              as.integer(run$unique_biosamples[[1L]]) == data.table::uniqueN(ledger$biosample) &&
              as.integer(run$unique_runs[[1L]]) == data.table::uniqueN(ledger$run_accession) &&
              run$animal_independence_basis[[1L]] ==
                "GEO_ENA_ONE_TO_ONE_BIOLOGICAL_REPLICATE_METADATA" &&
              as.integer(run$repository_inferred_animal_n[[1L]]) == nrow(ledger) &&
              run$group_counts[[1L]] == expected_group_count_text &&
              run$design_formula[[1L]] == "~ 0 + group" && run$design_rank[[1L]] == 3L &&
              run$coefficient_order[[1L]] == "sham,SNI,PLM" &&
              run$local_disease_expression[[1L]] == "SNI - sham" &&
              run$intervention_expression[[1L]] == "PLM - SNI",
            "RUN_MANIFEST_DESIGN_OR_CONTRAST_INVALID")
v2r4_assert(
  run$sample_include_primary_marker[[1L]] ==
    conditional_inclusion$sample_include_primary_marker &&
    run$sample_metadata_status[[1L]] == conditional_inclusion$sample_metadata_status &&
    run$conditional_inclusion_release_status[[1L]] ==
      conditional_inclusion$conditional_inclusion_release_status &&
    as.integer(run$conditional_inclusion_n[[1L]]) == conditional_inclusion$conditional_inclusion_n &&
    run$cohort_species_gate_sha256[[1L]] == conditional_inclusion$cohort_species_gate_sha256 &&
    run$sample_species_table_sha256[[1L]] == conditional_inclusion$sample_species_table_sha256 &&
    run$plm_metadata_conflict_policy[[1L]] ==
      conditional_inclusion$plm_metadata_conflict_policy,
  "RUN_MANIFEST_CONDITIONAL_INCLUSION_RELEASE_INVALID"
)
v2r4_assert(run$primary_method[[1L]] == "edgeR_QL_robust" &&
              run$primary_filtering[[1L]] == cfg$analysis$edgeR$filter &&
              run$normalization[[1L]] == "TMM" &&
              run$random_seed[[1L]] == cfg$analysis$random_seed &&
              run$orthology_authority[[1L]] == cfg$analysis$orthology_authority &&
              run$orthology_single_source_role[[1L]] == "orthology_single_source",
            "RUN_MANIFEST_METHOD_SEED_OR_ORTHOLOGY_INVALID")
v2r4_assert(run$assembly[[1L]] == cfg$input_snapshot$assembly &&
              run$assembly_accession[[1L]] == cfg$input_snapshot$assembly_accession &&
              run$selected_primary_strand_mode[[1L]] == 0L &&
              run$species_gate_status[[1L]] == "PASS_BEFORE_COUNTS_RDS_CAPTURE_OR_SEMANTIC_PARSE" &&
              run$gene_id_source_column[[1L]] == "GeneID" &&
              run$symbol_source_column[[1L]] == "gene" &&
              run$biotype_source_column[[1L]] == "gene_biotype" &&
              run$requested_DESeq2_fit_type[[1L]] == cfg$analysis$DESeq2$requested_fit_type_before_fit,
            "RUN_MANIFEST_REFERENCE_SYMBOL_OR_DESEQ_INVALID")
v2r4_assert_deseq_fit_provenance(
  run$actual_DESeq2_fit_type_after_DESeq[[1L]], run$DESeq2_fit_messages[[1L]], cfg,
  "RUN_MANIFEST_DESEQ2"
)
v2r4_assert(run$loao_method[[1L]] == "ACTIVITY_FULL_PIPELINE_REFIT" &&
              run$gate_status[[1L]] == cfg$gate$machine_status,
            "RUN_MANIFEST_LOAO_OR_GATE_INVALID")
v2r4_validate_uniform_flags(run, cfg, "RUN_MANIFEST")
add_check("RUN_MANIFEST_SEMANTICS", "external trust, design, methods and gate independently matched")

runtime_reported <- payload_tsv("audit/runtime_bindings.tsv", "RUNTIME_BINDINGS")
assert_keyed_equal(runtime_reported, runtime_bindings, "role",
                   c("runtime_relative_path", "bytes", "sha256", "status"),
                   "RUNTIME_BINDINGS")
v2r4_validate_uniform_flags(runtime_reported, cfg, "RUNTIME_BINDINGS")
lease_reported <- payload_tsv(
  "audit/lease_validation.tsv", "LEASE_VALIDATION",
  character_columns = c("model_start_utc", "validator_start_utc",
                        "guard_owner_start_utc", "guard_start_utc")
)
v2r4_assert(nrow(lease_reported) == 1L, "LEASE_REPORT_NOT_EXTERNALLY_REPRODUCED")
v2r4_assert_utc_fixed6(c(lease_reported$model_start_utc[[1L]],
                         lease_reported$validator_start_utc[[1L]],
                         lease_reported$guard_owner_start_utc[[1L]],
                         lease_reported$guard_start_utc[[1L]]),
                       "LEASE_REPORT_PROCESS_STARTS")
v2r4_assert(nrow(lease_reported) == 1L && lease_reported$lease_id[[1L]] == lease$lease_id[[1L]] &&
              lease_reported$holder_principal[[1L]] == lease$holder_principal[[1L]] &&
              lease_reported$run_id[[1L]] == lease$run_id[[1L]] &&
              as.integer(lease_reported$model_pid[[1L]]) == as.integer(lease$model_pid[[1L]]) &&
              identical(lease_reported$model_start_utc[[1L]], lease$model_start_utc[[1L]]) &&
              lease_reported$lease_manifest_sha256[[1L]] == lease_bundle$manifest$sha256 &&
              lease_reported$model_request_sha256[[1L]] == lease$model_request_sha256[[1L]] &&
              lease_reported$model_request_nonce[[1L]] == lease$model_request_nonce[[1L]] &&
              lease_reported$model_identity_request_sha256[[1L]] ==
                lease$model_identity_request_sha256[[1L]] &&
              as.integer(lease_reported$validator_pid[[1L]]) ==
                as.integer(lease$validator_pid[[1L]]) &&
              identical(lease_reported$validator_start_utc[[1L]],
                        lease$validator_start_utc[[1L]]) &&
              lease_reported$validator_request_sha256[[1L]] ==
                lease$validator_request_sha256[[1L]] &&
              lease_reported$validator_request_nonce[[1L]] ==
                lease$validator_request_nonce[[1L]] &&
              lease_reported$validator_identity_request_sha256[[1L]] ==
                lease$validator_identity_request_sha256[[1L]] &&
              lease_reported$snapshot_tree_sha256[[1L]] == snapshot$tree_digest &&
              lease_reported$mutex_relative_path[[1L]] == cfg$lease$mutex_relative_path &&
              as.integer(lease_reported$guard_pid[[1L]]) == guard_state$guard_pid &&
              identical(lease_reported$guard_start_utc[[1L]], guard_state$guard_start_utc) &&
              lease_reported$guard_nonce[[1L]] == guard_state$guard_nonce &&
              lease_reported$guard_public_key_sha256[[1L]] ==
                guard_state$guard_public_key_sha256 &&
              lease_reported$guard_ready_receipt_sha256[[1L]] ==
                guard_state$ready_receipt_sha256 &&
              lease_reported$guard_challenge_request_sha256[[1L]] ==
                lease$guard_challenge_request_sha256[[1L]] &&
              lease_reported$guard_challenge_response_sha256[[1L]] ==
                lease$guard_challenge_response_sha256[[1L]] &&
              lease_reported$guard_mutex_file_id[[1L]] == guard_state$mutex_file_id &&
              lease_reported$finalization_public_key_sha256[[1L]] ==
                lease$finalization_public_key_sha256[[1L]] &&
              v2r4_norm_existing(
                lease_reported$finalization_public_key_captured_path[[1L]]
              ) == guard_state$finalization_public_key_captured_path &&
              as.integer(lease_reported$finalization_public_key_captured_bytes[[1L]]) ==
                as.integer(lease$finalization_public_key_captured_bytes[[1L]]) &&
              lease_reported$finalization_public_key_captured_sha256[[1L]] ==
                lease$finalization_public_key_captured_sha256[[1L]],
            "LEASE_REPORT_NOT_EXTERNALLY_REPRODUCED")
v2r4_validate_uniform_flags(lease_reported, cfg, "LEASE_REPORT")
add_check("MODEL_REPORTED_RUNTIME_AND_LEASE_RECONCILIATION", "reported audit equals independent observations")

design_reported <- payload_tsv("diagnostics/design_matrix.tsv", "DESIGN_MATRIX")
expected_design <- data.table::data.table(sample_id = rownames(design), group = as.character(group),
                                          data.table::as.data.table(design))
assert_keyed_equal(design_reported, expected_design, "sample_id",
                   c("group", colnames(design)), "DESIGN_MATRIX", tolerance = 0)
contrast_reported <- payload_tsv("diagnostics/contrast_matrix.tsv", "CONTRAST_MATRIX")
expected_contrast <- data.table::data.table(coefficient = rownames(contrast_matrix),
                                            data.table::as.data.table(contrast_matrix))
assert_keyed_equal(contrast_reported, expected_contrast, "coefficient",
                   colnames(contrast_matrix), "CONTRAST_MATRIX", tolerance = 0)
v2r4_validate_uniform_flags(design_reported, cfg, "DESIGN_MATRIX")
v2r4_validate_uniform_flags(contrast_reported, cfg, "CONTRAST_MATRIX")

annotation_reported <- payload_tsv("diagnostics/annotation_mapping.tsv", "ANNOTATION_MAPPING")
assert_keyed_equal(annotation_reported, annotation_metrics_expected,
                   c("gene_id_source_column", "symbol_source_column", "biotype_source_column"),
                   c("total_gene_rows", "nonblank_symbol_n", "nonblank_symbol_coverage",
                     "unique_nonblank_symbol_n", "unique_nonblank_symbol_coverage",
                     "fallback_to_gene_id_n"), "ANNOTATION_MAPPING")
v2r4_validate_uniform_flags(annotation_reported, cfg, "ANNOTATION_MAPPING")
v2r4_assert(annotation_reported$symbol_source_column[[1L]] == "gene", "ANNOTATION_SOURCE_NOT_GENE")
add_check("ANNOTATION_GENE_SOURCE_AND_COVERAGE", "frozen gene column and coverage reproduced from signed RDS")

y0 <- edgeR::DGEList(counts = counts)
keep_filter_by_expr <- edgeR::filterByExpr(y0, design = design)
keep_cpm <- rowSums(edgeR::cpm(y0, log = FALSE) >= cfg$analysis$minimum_cpm) >=
  cfg$analysis$minimum_samples_at_cpm
keep <- keep_filter_by_expr & keep_cpm
minimum_retained <- if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$minimum_retained_genes else
  cfg$analysis$minimum_retained_genes_production
v2r4_assert(sum(keep) >= minimum_retained, "VALIDATOR_TOO_FEW_FILTERED_GENES")
filter_reported <- payload_tsv("diagnostics/filter_summary.tsv", "FILTER_SUMMARY")
v2r4_assert(nrow(filter_reported) == 1L && filter_reported$input_gene_n[[1L]] == nrow(counts) &&
              filter_reported$retained_gene_n[[1L]] == sum(keep) &&
              filter_reported$filterByExpr_retained_n[[1L]] == sum(keep_filter_by_expr) &&
              filter_reported$explicit_CPM_retained_n[[1L]] == sum(keep_cpm) &&
              filter_reported$normalization[[1L]] == "TMM" &&
              filter_reported$primary_method[[1L]] == "edgeR_QL_robust",
            "FILTER_SUMMARY_RECOMPUTATION_MISMATCH")
v2r4_validate_uniform_flags(filter_reported, cfg, "FILTER_SUMMARY")
retained_annotation <- annotation[match(rownames(counts)[keep], gene_id)]
add_check("FILTER_RECOMPUTATION", sprintf("%d genes independently retained", sum(keep)))

y <- edgeR::normLibSizes(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
y <- edgeR::estimateDisp(y, design, robust = TRUE)
ql_fit <- edgeR::glmQLFit(y, design, robust = TRUE)
expected_edger <- data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(j) {
  qlf <- edgeR::glmQLFTest(ql_fit, contrast = contrast_matrix[, j])
  tab <- data.table::as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table,
                                  keep.rownames = "gene_id")
  tab <- tab[match(rownames(ql_fit$counts), gene_id)]
  covariance <- v2r4_edgeR_ql_covariance(ql_fit, contrast_matrix[, j])
  critical <- stats::qt(0.975, covariance$df)
  tab[, `:=`(
    contrast = colnames(contrast_matrix)[[j]],
    se = covariance$se,
    ci_low = logFC - critical * covariance$se,
    ci_high = logFC + critical * covariance$se,
    se_df = covariance$df,
    se_status = covariance$status,
    df_total = qlf$df.total,
    q_value = FDR
  )]
  tab
}), fill = TRUE)
edge_reported <- payload_tsv("gene_level/edgeR_QL_all_contrasts.tsv.gz", "EDGER_QL")
expected_gene_keys <- data.table::CJ(gene_id = rownames(y$counts),
                                     contrast = colnames(contrast_matrix), unique = TRUE)
observed_gene_keys <- edge_reported[order(gene_id, contrast), .(
  gene_id = as.character(gene_id), contrast = as.character(contrast)
)]
expected_gene_keys <- expected_gene_keys[order(gene_id, contrast)]
v2r4_assert(identical(observed_gene_keys[["gene_id"]],
                      as.character(expected_gene_keys[["gene_id"]])) &&
                identical(observed_gene_keys[["contrast"]],
                          as.character(expected_gene_keys[["contrast"]])),
              "EDGER_GENE_CONTRAST_FAMILY_INVALID")
assert_keyed_equal(edge_reported, expected_edger, c("gene_id", "contrast"),
                   c("logFC", "F", "PValue", "FDR", "q_value", "se", "ci_low", "ci_high",
                     "se_df", "se_status", "df_total"),
                   "EDGER_PRIMARY_REFIT", tolerance = 1e-8)
v2r4_assert(all(edge_reported$method == "edgeR_QL_robust"), "EDGER_METHOD_LABEL_INVALID")
estimable_se <- edge_reported$se_status == "ESTIMABLE_QL_FISHER_COVARIANCE"
v2r4_assert(all(is.finite(edge_reported$se[estimable_se]) & edge_reported$se[estimable_se] > 0) &&
              all(is.finite(edge_reported$ci_low[estimable_se])) &&
              all(is.finite(edge_reported$ci_high[estimable_se])) &&
              all(edge_reported$ci_low[estimable_se] < edge_reported$ci_high[estimable_se]),
            "EDGER_ESTIMABLE_SE_CI_DOMAIN_INVALID")
if (any(!estimable_se)) {
  v2r4_assert(all(is.na(edge_reported$se[!estimable_se])) &&
                all(is.na(edge_reported$ci_low[!estimable_se])) &&
                all(is.na(edge_reported$ci_high[!estimable_se])),
              "EDGER_NONESTIMABLE_SE_CI_MUST_BE_NA")
}
v2r4_validate_uniform_flags(edge_reported, cfg, "EDGER_PRIMARY")
add_check("EDGER_QL_INDEPENDENT_REFIT_AND_COVARIANCE", "primary tests and covariance SE/CI reproduced")

voom_reported <- payload_tsv("gene_level/limma_voom_all_contrasts.tsv.gz", "LIMMA_VOOM")
deseq_reported <- payload_tsv("gene_level/DESeq2_all_contrasts.tsv.gz", "DESEQ2")
set.seed(as.integer(cfg$analysis$random_seed))
voom_input <- limma::voom(y, design = design, plot = FALSE)
voom_fit <- limma::eBayes(
  limma::contrasts.fit(limma::lmFit(voom_input, design), contrast_matrix),
  robust = TRUE
)
voom_expected <- data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(j) {
  se <- voom_fit$stdev.unscaled[, j] * sqrt(voom_fit$s2.post)
  critical <- stats::qt(0.975, voom_fit$df.total)
  data.table::data.table(
    gene_id = rownames(voom_fit$coefficients),
    symbol = retained_annotation$symbol[
      match(rownames(voom_fit$coefficients), retained_annotation$gene_id)
    ],
    contrast = colnames(contrast_matrix)[[j]],
    method = "limma_voom_robust_eBayes",
    logFC = voom_fit$coefficients[, j],
    t = voom_fit$t[, j],
    p_value = voom_fit$p.value[, j],
    q_value = stats::p.adjust(voom_fit$p.value[, j], method = "BH"),
    se = se,
    ci_low = voom_fit$coefficients[, j] - critical * se,
    ci_high = voom_fit$coefficients[, j] + critical * se,
    df_total = voom_fit$df.total
  )
}))
assert_complete_numeric_family(
  voom_reported, voom_expected, c("gene_id", "contrast"),
  c("symbol", "method"),
  c("logFC", "t", "p_value", "q_value", "se", "ci_low", "ci_high", "df_total"),
  "LIMMA_VOOM_COMPLETE_INDEPENDENT_REFIT", cfg
)
for (contrast_name in colnames(contrast_matrix)) {
  idx <- base::which(voom_reported[["contrast"]] == contrast_name)
  expected_q <- stats::p.adjust(voom_reported$p_value[idx], method = "BH")
  bound <- cfg$analysis$independent_numeric_revalidation$absolute_tolerance +
    cfg$analysis$independent_numeric_revalidation$relative_tolerance *
      pmax(abs(expected_q), abs(voom_reported$q_value[idx]))
  v2r4_assert(identical(is.na(voom_reported$q_value[idx]), is.na(expected_q)) &&
                all(abs(voom_reported$q_value[idx] - expected_q) <= bound, na.rm = TRUE),
              "LIMMA_VOOM_CONTRAST_BH_RECOMPUTATION_MISMATCH:%s", contrast_name)
}
v2r4_validate_uniform_flags(voom_reported, cfg, "SENSITIVITY_VOOM")

deseq_messages_raw <- character()
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(y$counts),
  colData = data.frame(group = group, row.names = colnames(y$counts)),
  design = ~ 0 + group
)
dds <- withCallingHandlers(
  DESeq2::DESeq(
    dds, test = "Wald", fitType = cfg$analysis$DESeq2$requested_fit_type_before_fit,
    betaPrior = FALSE, minReplicatesForReplace = Inf, quiet = TRUE
  ),
  message = function(m) {
    deseq_messages_raw <<- c(deseq_messages_raw, conditionMessage(m))
    invokeRestart("muffleMessage")
  },
  warning = function(w) {
    deseq_messages_raw <<- c(deseq_messages_raw, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
deseq_actual_fit <- attr(DESeq2::dispersionFunction(dds), "fitType") %||% "UNKNOWN"
deseq_messages <- v2r4_collapse_condition_messages(deseq_messages_raw)
v2r4_assert_deseq_fit_provenance(deseq_actual_fit, deseq_messages, cfg,
                                 "VALIDATOR_DESEQ2")
result_names <- DESeq2::resultsNames(dds)
expected_names <- paste0("group", rownames(contrast_matrix))
v2r4_assert(all(expected_names %in% result_names),
            "VALIDATOR_DESEQ2_COEFFICIENT_NAMES_INVALID")
deseq_expected <- data.table::rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(j) {
  contrast_vector <- numeric(length(result_names))
  names(contrast_vector) <- result_names
  contrast_vector[expected_names] <- contrast_matrix[, j]
  result <- DESeq2::results(
    dds, contrast = contrast_vector, alpha = 0.05,
    independentFiltering = FALSE, cooksCutoff = FALSE
  )
  value <- data.table::as.data.table(as.data.frame(result), keep.rownames = "gene_id")
  value[, `:=`(
    symbol = retained_annotation$symbol[match(gene_id, retained_annotation$gene_id)],
    contrast = colnames(contrast_matrix)[[j]],
    method = "DESeq2_Wald_sensitivity",
    requested_fit_type = cfg$analysis$DESeq2$requested_fit_type_before_fit,
    actual_fit_type_after_DESeq = deseq_actual_fit,
    ci_low = log2FoldChange - stats::qnorm(0.975) * lfcSE,
    ci_high = log2FoldChange + stats::qnorm(0.975) * lfcSE,
    p_value = pvalue,
    q_value = padj
  )]
  value
}), fill = TRUE)
assert_complete_numeric_family(
  deseq_reported, deseq_expected, c("gene_id", "contrast"),
  c("symbol", "method", "requested_fit_type", "actual_fit_type_after_DESeq"),
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
    "ci_low", "ci_high", "p_value", "q_value"),
  "DESEQ2_COMPLETE_INDEPENDENT_REFIT", cfg
)
for (contrast_name in colnames(contrast_matrix)) {
  idx <- base::which(deseq_reported[["contrast"]] == contrast_name)
  expected_q <- stats::p.adjust(deseq_reported$pvalue[idx], method = "BH")
  v2r4_assert(identical(is.na(deseq_reported$padj[idx]), is.na(expected_q)),
              "DESEQ2_CONTRAST_BH_NA_PATTERN_MISMATCH:%s", contrast_name)
  finite <- !is.na(expected_q)
  if (any(finite)) {
    bound <- cfg$analysis$independent_numeric_revalidation$absolute_tolerance +
      cfg$analysis$independent_numeric_revalidation$relative_tolerance *
        pmax(abs(expected_q[finite]), abs(deseq_reported$padj[idx][finite]))
    v2r4_assert(all(abs(deseq_reported$padj[idx][finite] - expected_q[finite]) <= bound),
                "DESEQ2_CONTRAST_BH_RECOMPUTATION_MISMATCH:%s", contrast_name)
  }
}
v2r4_assert(identical(deseq_actual_fit,
                        run$actual_DESeq2_fit_type_after_DESeq[[1L]]) &&
              identical(deseq_messages, run$DESeq2_fit_messages[[1L]]) &&
              all(deseq_reported$actual_fit_type_after_DESeq == deseq_actual_fit),
            "DESEQ2_INDEPENDENT_FIT_PROVENANCE_MISMATCH")
v2r4_validate_uniform_flags(deseq_reported, cfg, "SENSITIVITY_DESEQ2")
add_check(
  "VOOM_AND_DESEQ2_COMPLETE_INDEPENDENT_REFIT",
  sprintf("%d retained genes x %d exact contrasts; full numerics, BH and fit provenance",
          nrow(y$counts), ncol(contrast_matrix))
)

gate <- payload_tsv("gates/E2_GATE.tsv", "E2_GATE")
v2r4_assert(nrow(gate) == 1L && gate$dataset_id[[1L]] == cfg$dataset_id &&
              gate$gate[[1L]] == "E2" && gate$upstream_gate[[1L]] == "E1" &&
              gate$eligible_input_n[[1L]] == 0L && gate$evaluated_candidate_n[[1L]] == 0L &&
              is.na(gate$pass_n[[1L]]) && gate$status[[1L]] == "E2_NOT_IDENTIFIABLE_NO_E1_INPUT",
            "E2_GATE_SEMANTICS_INVALID")
v2r4_validate_uniform_flags(gate, cfg, "E2_GATE")
add_check("E2_STRUCTURAL_NONIDENTIFIABILITY", "E2_NOT_IDENTIFIABLE_NO_E1_INPUT; no candidate eligibility")

gene_sets <- v2r4_read_rds_raw(
  v2r4_object(snapshot, "gene_set_object"), "GENE_SET_OBJECT",
  cfg$input_snapshot$rds_role_contracts$gene_set_object$compression
)
orthology <- v2r4_read_tsv_raw(v2r4_object(snapshot, "orthology_single_source"),
                               "ORTHOLOGY_SINGLE_SOURCE")
v2r4_assert(is.list(gene_sets) && all(c(
  "reactome", "progeny_human", "dorothea_rat", "syngo_sensitivity"
) %in% names(gene_sets)), "GENE_SET_OBJECT_SCHEMA_INVALID")
v2r4_assert(all(c("group_method", "human_symbol", "other_tax_id", "other_symbol") %in% names(orthology)),
            "ORTHOLOGY_SCHEMA_INVALID")
orthology_idx <- base::which(orthology[["group_method"]] == "NCBI Ortholog" &
                              orthology[["other_tax_id"]] == 10116L)
rat_map <- unique(data.table::data.table(
  human_symbol = as.character(orthology$human_symbol[orthology_idx]),
  rat_symbol = as.character(orthology$other_symbol[orthology_idx])
))
v2r4_assert(!anyDuplicated(rat_map$human_symbol) && !anyDuplicated(rat_map$rat_symbol),
            "ORTHOLOGY_NOT_STRICT_ONE_TO_ONE")
retained_symbols <- unique(retained_annotation$symbol)

rebuild_unsigned <- function(source_table, collection_name) {
  source_table <- data.table::as.data.table(source_table)
  v2r4_assert(all(c("term_id", "term_name", "gene_symbol") %in% names(source_table)),
              "UNSIGNED_SOURCE_SCHEMA_INVALID:%s", collection_name)
  source_unique <- unique(data.table::data.table(
    program_id = as.character(source_table$term_id),
    program_name = as.character(source_table$term_name),
    human_symbol = as.character(source_table$gene_symbol)
  ))
  mapped <- unique(merge(source_unique, rat_map, by = "human_symbol", allow.cartesian = TRUE))
  mapped[, symbol := rat_symbol]
  coverage <- mapped[, .(mapped_symbol_n = data.table::uniqueN(symbol)),
                     by = .(program_id, program_name)]
  coverage[, tested_feature_n := vapply(program_id, function(program_value) {
    idx <- base::which(mapped[["program_id"]] == program_value)
    data.table::uniqueN(intersect(mapped$symbol[idx], retained_symbols))
  }, integer(1))]
  coverage[, `:=`(
    collection = collection_name,
    source_feature_n = mapped_symbol_n,
    eligible_10_500 = tested_feature_n >= cfg$analysis$program_size_min &
      tested_feature_n <= cfg$analysis$program_size_max,
    tested = tested_feature_n >= cfg$analysis$program_size_min &
      tested_feature_n <= cfg$analysis$program_size_max,
    confidence_domain = "NOT_APPLICABLE",
    mapping_rule = cfg$analysis$orthology_authority
  )]
  list(coverage = coverage, tested = coverage[tested == TRUE]$program_id)
}

unsigned_expected <- list(
  Reactome = rebuild_unsigned(gene_sets$reactome, "Reactome"),
  SynGO_sensitivity = rebuild_unsigned(gene_sets$syngo_sensitivity, "SynGO_sensitivity")
)

progeny_human <- data.table::as.data.table(gene_sets$progeny_human)
v2r4_assert(all(c("pathway", "gene_symbol", "weight") %in% names(progeny_human)),
            "PROGENY_SOURCE_SCHEMA_INVALID")
progeny_net <- unique(merge(data.table::data.table(
  source = as.character(progeny_human$pathway),
  human_symbol = as.character(progeny_human$gene_symbol),
  mor = as.numeric(progeny_human$weight)
), rat_map, by = "human_symbol", allow.cartesian = TRUE)[,
  .(source, target = rat_symbol, mor)])
dorothea_raw <- data.table::as.data.table(gene_sets$dorothea_rat)
v2r4_assert(all(c("tf", "target", "mor", "confidence") %in% names(dorothea_raw)),
            "DOROTHEA_SOURCE_SCHEMA_INVALID")
dorothea_net <- unique(data.table::data.table(
  source = as.character(dorothea_raw$tf),
  target = as.character(dorothea_raw$target),
  mor = as.numeric(dorothea_raw$mor),
  confidence = as.character(dorothea_raw$confidence)
))
v2r4_assert(all(dorothea_net$confidence %in% c("A", "B")), "DOROTHEA_CONFIDENCE_NOT_EXACT_AB")

rebuild_activity_coverage <- function(network, collection_name, confidence_value) {
  out <- network[, .(
    source_feature_n = data.table::uniqueN(target),
    tested_feature_n = data.table::uniqueN(intersect(target, retained_symbols))
  ), by = source]
  out[, `:=`(
    collection = collection_name,
    eligible_10_500 = tested_feature_n >= cfg$analysis$program_size_min &
      tested_feature_n <= cfg$analysis$program_size_max,
    tested = tested_feature_n >= cfg$analysis$program_size_min &
      tested_feature_n <= cfg$analysis$program_size_max,
    confidence_domain = confidence_value,
    mapping_rule = if (collection_name == "PROGENy") cfg$analysis$orthology_authority else
      "frozen_rat_AB_strict_one_to_one_network"
  )]
  out
}

progeny_coverage <- rebuild_activity_coverage(progeny_net, "PROGENy", "NOT_APPLICABLE")
dorothea_coverage <- rebuild_activity_coverage(
  dorothea_net, "DoRothEA_AB", paste(sort(unique(dorothea_net$confidence)), collapse = ",")
)
domains <- if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$program_domains else cfg$analysis$program_domains
v2r4_assert(data.table::uniqueN(progeny_net$source) == domains$PROGENy$source_program_n &&
              setequal(unique(progeny_net$source), unlist(domains$PROGENy$exact_program_ids, use.names = FALSE)) &&
              sum(progeny_coverage$tested) == domains$PROGENy$eligible_tested_program_n,
            "PROGENY_SOURCE_OR_TESTED_DOMAIN_INVALID")
v2r4_assert(data.table::uniqueN(dorothea_net$source) == domains$DoRothEA_AB$frozen_source_program_upper_bound &&
              sum(dorothea_coverage$tested) == domains$DoRothEA_AB$eligible_tested_program_n &&
              paste(sort(unique(dorothea_net$confidence)), collapse = ",") == "A,B",
            "DOROTHEA_SOURCE_OR_TESTED_DOMAIN_INVALID")

mapping_expected <- data.table::rbindlist(list(
  unsigned_expected$Reactome$coverage,
  unsigned_expected$SynGO_sensitivity$coverage,
  progeny_coverage,
  dorothea_coverage
), fill = TRUE)
mapping_reported <- payload_tsv("programs/program_mapping_coverage.tsv", "PROGRAM_MAPPING_COVERAGE")
normalize_mapping <- function(x) {
  x <- data.table::copy(data.table::as.data.table(x))
  if (!"program_id" %in% names(x)) x[, program_id := NA_character_]
  if (!"source" %in% names(x)) x[, source := NA_character_]
  x[, canonical_program_id := ifelse(!is.na(program_id) & nzchar(program_id), program_id, source)]
  x[, .(
    collection = as.character(collection),
    canonical_program_id = as.character(canonical_program_id),
    source_feature_n = as.integer(source_feature_n),
    tested_feature_n = as.integer(tested_feature_n),
    eligible_10_500 = as.logical(eligible_10_500),
    tested = as.logical(tested),
    confidence_domain = as.character(confidence_domain),
    mapping_rule = as.character(mapping_rule)
  )]
}
mapping_expected_normalized <- normalize_mapping(mapping_expected)
mapping_reported_normalized <- normalize_mapping(mapping_reported)
assert_keyed_equal(mapping_reported_normalized, mapping_expected_normalized,
                   c("collection", "canonical_program_id"),
                   c("source_feature_n", "tested_feature_n", "eligible_10_500", "tested",
                     "confidence_domain", "mapping_rule"),
                   "PROGRAM_MAPPING_COVERAGE", tolerance = 0)
v2r4_validate_uniform_flags(mapping_reported, cfg, "PROGRAM_MAPPING_COVERAGE")
add_check("PROGRAM_FAMILY_RECONSTRUCTION", "unsigned and signed tested families rebuilt from signed inputs")

tested_families <- list(
  Reactome = sort(unsigned_expected$Reactome$tested),
  SynGO_sensitivity = sort(unsigned_expected$SynGO_sensitivity$tested),
  PROGENy = sort(progeny_coverage[tested == TRUE]$source),
  DoRothEA_AB = sort(dorothea_coverage[tested == TRUE]$source)
)

validate_program_table <- function(table_value, collections, expected_method, label,
                                  p_column = "p_value", q_column = "q_within_collection_contrast") {
  table_value <- data.table::as.data.table(table_value)
  v2r4_assert(all(c("collection", "program_id", "contrast", "method", p_column, q_column) %in%
                    names(table_value)), "%s_SCHEMA_INVALID", label)
  expected_families <- tested_families[collections]
  v2r4_assert_exact_program_grid(table_value, expected_families, collections,
                                 colnames(contrast_matrix), label)
  v2r4_assert(all(table_value$method == expected_method) &&
                all(is.finite(table_value[[p_column]]) & table_value[[p_column]] >= 0 &
                      table_value[[p_column]] <= 1) &&
                all(is.finite(table_value[[q_column]]) & table_value[[q_column]] >= 0 &
                      table_value[[q_column]] <= 1), "%s_METHOD_OR_NUMERIC_DOMAIN_INVALID", label)
  bh_input <- data.table::copy(table_value)
  bh_input[, tested_for_bh__ := TRUE]
  v2r4_recompute_bh(bh_input, p_col = p_column, q_col = q_column, tested_col = "tested_for_bh__")
  v2r4_validate_uniform_flags(table_value, cfg, label)
  invisible(TRUE)
}

activity <- payload_tsv("programs/activity_contrasts.tsv.gz", "ACTIVITY_CONTRASTS")
validate_program_table(activity, c("PROGENy", "DoRothEA_AB"),
                       "decoupleR_ULM_then_limma", "ACTIVITY_CONTRASTS")
unsigned_camera <- payload_tsv("programs/unsigned_camera.tsv.gz", "UNSIGNED_CAMERA")
validate_program_table(unsigned_camera, c("Reactome", "SynGO_sensitivity"),
                       "camera_inter_gene_cor_0.01", "UNSIGNED_CAMERA")
unsigned_mroast <- payload_tsv("programs/unsigned_mroast.tsv.gz", "UNSIGNED_MROAST")
expected_mroast_method <- paste0("mroast_mean_",
  if (identical(cli$execution_mode, "SYNTHETIC")) cfg$synthetic_contract$mroast_rotations else cfg$analysis$mroast_rotations,
  "_rotations")
validate_program_table(unsigned_mroast, c("Reactome", "SynGO_sensitivity"),
                       expected_mroast_method, "UNSIGNED_MROAST")
unsigned_fgsea <- payload_tsv("programs/unsigned_fgseaMultilevel.tsv.gz", "UNSIGNED_FGSEA")
validate_program_table(unsigned_fgsea, c("Reactome", "SynGO_sensitivity"),
                       "fgseaMultilevel", "UNSIGNED_FGSEA")
add_check("PROGRAM_MULTIPLICITY_RECOMPUTATION", "all collection-by-contrast BH families independently reproduced")

v2r4_assert(run$progeny_source_program_n[[1L]] == domains$PROGENy$source_program_n &&
              run$progeny_tested_program_n[[1L]] == length(tested_families$PROGENy) &&
              run$dorothea_source_program_n[[1L]] == domains$DoRothEA_AB$frozen_source_program_upper_bound &&
              run$dorothea_tested_program_n[[1L]] == length(tested_families$DoRothEA_AB) &&
              run$dorothea_confidence_domain[[1L]] == "A,B",
            "RUN_MANIFEST_PROGRAM_DOMAIN_MISMATCH")

loao_long <- payload_tsv("programs/activity_full_pipeline_LOAO_long.tsv.gz", "LOAO_LONG")
expected_loao_grid <- data.table::rbindlist(lapply(c("PROGENy", "DoRothEA_AB"),
                                                   function(collection_name) {
  data.table::CJ(
    collection = collection_name,
    program_id = tested_families[[collection_name]],
    contrast = colnames(contrast_matrix),
    omitted_animal = ids,
    unique = TRUE
  )
}))
v2r4_assert(identical(
  loao_long[order(collection, program_id, contrast, omitted_animal),
            .(collection, program_id, contrast, omitted_animal)],
  expected_loao_grid[order(collection, program_id, contrast, omitted_animal)]
), "LOAO_LONG_EXACT_NINE_ANIMAL_GRID_INVALID")
expected_omitted_group <- setNames(as.character(group), ids)
v2r4_assert(all(loao_long$omitted_group == expected_omitted_group[loao_long$omitted_animal]) &&
              all(loao_long$loao_method == "ACTIVITY_FULL_PIPELINE_REFIT") &&
              all(loao_long$estimable == is.finite(loao_long$estimate)) &&
              all((loao_long$estimable & loao_long$nonestimable_reason == "NONE") |
                    (!loao_long$estimable & loao_long$nonestimable_reason ==
                       "PROGRAM_INELIGIBLE_OR_NOT_SCORED_AFTER_FULL_REFIT")),
            "LOAO_LONG_ESTIMABILITY_OR_OMISSION_SEMANTICS_INVALID")
loao_bh <- data.table::copy(loao_long)
loao_bh[, expected_q__ := stats::p.adjust(p_value, method = "BH"),
        by = .(collection, contrast, omitted_animal)]
same_loao_q <- (is.na(loao_bh$q_within_collection_contrast) & is.na(loao_bh$expected_q__)) |
  (is.finite(loao_bh$q_within_collection_contrast) & is.finite(loao_bh$expected_q__) &
     abs(loao_bh$q_within_collection_contrast - loao_bh$expected_q__) <= 1e-12)
v2r4_assert(all(same_loao_q), "LOAO_WITHIN_OMISSION_BH_MISMATCH")
loao_summary_expected <- v2r4_recompute_loao_summary(loao_long, ids)
loao_summary <- payload_tsv("programs/activity_full_pipeline_LOAO_summary.tsv.gz", "LOAO_SUMMARY")
assert_keyed_equal(loao_summary, loao_summary_expected,
                   c("collection", "program_id", "contrast"),
                   c("n_omissions", "n_estimable", "estimate_min", "estimate_max",
                     "estimate_mean", "sign_consistency"),
                   "LOAO_SUMMARY_RECOMPUTATION", tolerance = 1e-12)
minimum_estimable <- if (identical(cli$execution_mode, "SYNTHETIC")) {
  cfg$synthetic_contract$minimum_estimable_omissions_per_program_contrast
} else cfg$analysis$loao$minimum_estimable_omissions_per_program_contrast
v2r4_assert(all(loao_summary$n_omissions == 9L) &&
              all(loao_summary$n_estimable >= minimum_estimable) &&
              all(loao_summary$all_steps_refit %in% TRUE) &&
              all(loao_summary$loao_method == "ACTIVITY_FULL_PIPELINE_REFIT"),
            "LOAO_SUMMARY_FULL_REFIT_OR_ESTIMABILITY_INVALID")
v2r4_validate_uniform_flags(loao_long, cfg, "LOAO_LONG")
v2r4_validate_uniform_flags(loao_summary, cfg, "LOAO_SUMMARY")
add_check("LOAO_NINE_ANIMAL_KEYED_RECOMPUTATION", "long grid and every summary statistic independently rebuilt")

verification_reported <- payload_tsv("audit/input_snapshot_verification.tsv", "INPUT_SNAPSHOT_VERIFICATION")
assert_keyed_equal(verification_reported, snapshot$verification, "role",
                   c("materialization", "relative_path", "expected_bytes", "observed_bytes",
                     "expected_sha256", "observed_sha256", "status"),
                   "INPUT_SNAPSHOT_VERIFICATION", tolerance = 0)
v2r4_validate_uniform_flags(verification_reported, cfg, "INPUT_SNAPSHOT_VERIFICATION")
add_check("SNAPSHOT_TREE_AND_CAPTURE_LEDGER", snapshot$tree_digest)

all_payloads <- unlist(cfg$required_output_payloads, use.names = FALSE)
for (relative_path in all_payloads) {
  table_value <- payload_tsv(relative_path, paste0("UNIFORM_FLAGS_", basename(relative_path)))
  v2r4_validate_uniform_flags(table_value, cfg, paste0("PAYLOAD_", basename(relative_path)))
}
add_check("ALL_PAYLOAD_CANDIDATE_FLAGS_FALSE", sprintf("%d/%d payloads", length(all_payloads), length(all_payloads)))

v2r4_validate_lease_buffer(
  lease_bundle$manifest, cfg, snapshot$tree_digest, as.character(cli$run_id),
  as.character(cli$holder_principal), runtime_buffers$os_guard_source, guard_state,
  launch_context, runtime_role_digest, state_root,
  expected_model_pid = model_identity_request$process_pid,
  expected_model_start_utc = model_identity_request$process_start_utc,
  expected_model_request_sha256 = model_request_buffer_prelease$sha256,
  expected_model_request_nonce = model_request_bundle_prelease$request$request_nonce,
  expected_model_identity_request_sha256 = model_identity_request_buffer$sha256,
  expected_validator_pid = Sys.getpid(),
  expected_validator_start_utc = validator_identity$start_utc,
  expected_validator_request_sha256 = child_request$buffer$sha256,
  expected_validator_request_nonce = child_request$request$request_nonce,
  expected_validator_identity_request_sha256 = validator_identity_request_buffer$sha256
)
final_validator_guard_challenge <- v2r4_guard_challenge(
  guard_state, runtime_buffers$os_guard_source, cfg, Sys.getpid()
)
add_check("FINAL_LIVE_LEASE_RECHECK",
          "model and guard PID/start, fixed mutex and fresh signed challenge remain live")

report <- data.table::rbindlist(checks, fill = TRUE)
v2r4_assert(all(report$status == "PASS") &&
              all(report$candidate_gate_eligible %in% FALSE) &&
              all(report$can_generate_or_upgrade_candidate %in% FALSE) &&
              all(report$formal_sequential_inference %in% FALSE),
            "VALIDATION_REPORT_INTERNAL_STATUS_INVALID")
v2r4_write_tsv(report, report_path)
report_buffer <- v2r4_capture_file(report_path, "FINAL_VALIDATION_REPORT", attempt_root)

model_request_buffer <- v2r4_capture_file(
  model_request_path, "VALIDATOR_CAPTURED_MODEL_REQUEST", attempt_root
)
model_request_bundle <- v2r4_parse_child_request_buffer(
  model_request_buffer, model_request_buffer$sha256, "MODEL"
)
v2r4_assert(
  identical(model_request_buffer$sha256, model_request_buffer_prelease$sha256) &&
    identical(model_request_bundle$request$request_nonce,
              model_request_bundle_prelease$request$request_nonce),
  "MODEL_CHILD_REQUEST_DRIFT_AFTER_SIGNED_LEASE"
)
v2r4_validate_child_request_bindings(model_request_bundle, cfg, launch_context, "MODEL")
v2r4_assert(
  identical(child_request$request$peer_request_nonce,
            model_request_bundle$request$request_nonce) &&
    identical(model_request_bundle$request$peer_request_nonce,
              child_request$request$request_nonce),
  "VALIDATOR_MODEL_REQUEST_NONCE_CROSS_BINDING_INVALID"
)
validator_identity <- v2r4_process_identity(
  runtime_buffers$os_guard_source, Sys.getpid()
)
v2r4_assert_process_identity_exact(
  validator_identity, Sys.getpid(), validator_identity$start_utc
)
v2r4_assert(
  identical(as.integer(lease$validator_pid[[1L]]), as.integer(Sys.getpid())) &&
    identical(lease$validator_start_utc[[1L]], validator_identity$start_utc) &&
    identical(lease$validator_request_sha256[[1L]], child_request$buffer$sha256) &&
    identical(lease$validator_request_nonce[[1L]], child_request$request$request_nonce) &&
    identical(lease$validator_identity_request_sha256[[1L]],
              validator_identity_request_buffer$sha256),
  "VALIDATOR_RUNTIME_IDENTITY_NOT_EXACTLY_SIGNED_IN_LEASE"
)
launcher_identity <- v2r4_process_identity(
  runtime_buffers$os_guard_source, as.integer(launch_context$launcher_pid)
)
model_identity <- v2r4_process_identity(
  runtime_buffers$os_guard_source, as.integer(lease$model_pid[[1L]])
)
guard_identity <- v2r4_process_identity(
  runtime_buffers$os_guard_source, guard_state$guard_pid
)
v2r4_assert_process_identity_exact(
  launcher_identity, as.integer(launch_context$launcher_pid),
  launch_context$launcher_start_utc
)
v2r4_assert_process_identity_exact(
  model_identity, as.integer(lease$model_pid[[1L]]), lease$model_start_utc[[1L]]
)
v2r4_assert_process_identity_exact(
  guard_identity, guard_state$guard_pid, guard_state$guard_start_utc
)
finalization_request <- list(
  schema_version = cfg$finalization$finalization_request_schema_version,
  status = "VALIDATOR_INDEPENDENT_RECOMPUTATION_PASS_AWAITING_EXTERNAL_FINAL_ACK",
  dataset_id = cfg$dataset_id,
  model_version = cfg$model_version,
  run_id = as.character(cli$run_id),
  lease_id = lease$lease_id[[1L]],
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  holder_principal = as.character(cli$holder_principal),
  source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
  model_source_sha256 = runtime_buffers$model_source$sha256,
  validator_source_sha256 = script_runtime_buffer$sha256,
  shared_library_sha256 = library_runtime_buffer$sha256,
  analysis_config_sha256 = config_buffer$sha256,
  external_launcher_sha256 = launch_context$bootstrap_launcher_sha256,
  os_guard_sha256 = runtime_buffers$os_guard_source$sha256,
  launcher_pid = as.integer(launch_context$launcher_pid),
  launcher_start_utc = launch_context$launcher_start_utc,
  guard_owner_role = guard_state$owner_role,
  guard_owner_pid = as.integer(guard_state$owner_pid),
  guard_owner_start_utc = guard_state$owner_start_utc,
  guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
  finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
  finalization_public_key_captured_path =
    guard_state$finalization_public_key_captured_path,
  finalization_public_key_captured_bytes =
    guard_state$finalization_public_key_captured_bytes,
  finalization_public_key_captured_sha256 =
    guard_state$finalization_public_key_captured_sha256,
  model_pid = as.integer(lease$model_pid[[1L]]),
  model_start_utc = lease$model_start_utc[[1L]],
  validator_pid = as.integer(Sys.getpid()),
  validator_start_utc = validator_identity$start_utc,
  guard_pid = guard_state$guard_pid,
  guard_start_utc = guard_state$guard_start_utc,
  guard_nonce = guard_state$guard_nonce,
  guard_public_key_sha256 = guard_state$guard_public_key_sha256,
  guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
  guard_mutex_file_id = guard_state$mutex_file_id,
  model_request_sha256 = model_request_buffer$sha256,
  model_request_nonce = model_request_bundle$request$request_nonce,
  model_identity_request_sha256 = model_identity_request_buffer$sha256,
  validator_request_sha256 = child_request$buffer$sha256,
  validator_request_nonce = child_request$request$request_nonce,
  validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  output_manifest_sha256 = authorized_output$markers$output_manifest$sha256,
  complete_sha256 = authorized_output$markers$complete$sha256,
  validation_report_sha256 = report_buffer$sha256,
  snapshot_tree_sha256 = snapshot$tree_digest,
  runtime_pre_digest = launch_context$runtime_pre_digest,
  lease_issued_utc = lease$issued_utc[[1L]],
  lease_expires_utc = lease$expires_utc[[1L]],
  validator_guard_challenge_request_sha256 = final_validator_guard_challenge$request$sha256,
  validator_guard_challenge_response_sha256 = final_validator_guard_challenge$response$sha256,
  guard_ready_receipt_path = guard_state$ready_path,
  guard_challenge_directory = guard_state$challenge_directory,
  e2_gate_status = cfg$gate$machine_status,
  candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)
v2r4_assert(identical(names(finalization_request), v2r4_finalization_request_fields()),
            "FINALIZATION_REQUEST_FIELD_ORDER_INVALID")
finalization_request_buffer <- v2r4_write_json_atomic(
  finalization_request, finalization_request_path, attempt_root
)
science_complete <- list(
  schema_version = cfg$finalization$science_complete_schema_version,
  status = "SCIENCE_RECOMPUTATION_COMPLETE_PROCESS_HELD",
  dataset_id = cfg$dataset_id,
  model_version = cfg$model_version,
  run_id = as.character(cli$run_id),
  validator_pid = as.integer(Sys.getpid()),
  validator_start_utc = validator_identity$start_utc,
  validator_request_sha256 = child_request$buffer$sha256,
  validator_request_nonce = child_request$request$request_nonce,
  validation_report_sha256 = report_buffer$sha256,
  output_manifest_sha256 = authorized_output$markers$output_manifest$sha256,
  complete_sha256 = authorized_output$markers$complete$sha256,
  finalization_request_sha256 = finalization_request_buffer$sha256,
  issued_utc = v2r4_now_fixed6(),
  e2_gate_status = cfg$gate$machine_status,
  candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)
v2r4_assert(identical(names(science_complete), v2r4_science_complete_fields()),
            "SCIENCE_COMPLETE_FIELD_ORDER_INVALID")
science_complete_buffer <- v2r4_write_json_atomic(
  science_complete, science_complete_path, attempt_root
)

final_wait <- cfg$execution_contract$signed_final_ack_wait_seconds
v2r4_assert(identical(final_wait, 900L) &&
              identical(launch_context$final_ack_wait_seconds, 900L),
            "VALIDATOR_FINAL_CONTROL_WAIT_MUST_BE_CAPTURED_EXACT_900")
v2r4_wait_for_files(
  runtime_post_receipt_path, final_wait,
  "VALIDATOR_PRE_FINAL_CHALLENGE_RUNTIME_POST"
)
runtime_post <- v2r4_validate_runtime_post_receipt(
  runtime_post_receipt_path, cfg,
  list(
    run_id = as.character(cli$run_id),
    runtime_pre_digest = launch_context$runtime_pre_digest,
    runtime_post_digest = launch_context$runtime_pre_digest,
    runtime_pre_post_equal = TRUE,
    runtime_semantic_digest = launch_context$runtime_semantic_digest,
    execution_identity_manifest_sha256 =
      launch_context$execution_identity_manifest_sha256,
    analysis_config_sha256 = config_buffer$sha256,
    science_complete_sha256 = science_complete_buffer$sha256,
    validator_pid = as.integer(Sys.getpid()),
    validator_start_utc = validator_identity$start_utc,
    validator_request_sha256 = child_request$buffer$sha256,
    validator_request_nonce = child_request$request$request_nonce,
    wait_seconds = 900L
  ), attempt_root
)
v2r4_validate_lease_buffer(
  lease_bundle$manifest, cfg, snapshot$tree_digest, as.character(cli$run_id),
  as.character(cli$holder_principal), runtime_buffers$os_guard_source, guard_state,
  launch_context, runtime_role_digest, state_root,
  expected_model_pid = model_identity_request$process_pid,
  expected_model_start_utc = model_identity_request$process_start_utc,
  expected_model_request_sha256 = model_request_buffer_prelease$sha256,
  expected_model_request_nonce = model_request_bundle_prelease$request$request_nonce,
  expected_model_identity_request_sha256 = model_identity_request_buffer$sha256,
  expected_validator_pid = Sys.getpid(),
  expected_validator_start_utc = validator_identity$start_utc,
  expected_validator_request_sha256 = child_request$buffer$sha256,
  expected_validator_request_nonce = child_request$request$request_nonce,
  expected_validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  now_utc = Sys.time()
)
for (identity_check in list(
  list(v2r4_process_identity(runtime_buffers$os_guard_source,
                             as.integer(launch_context$launcher_pid)),
       as.integer(launch_context$launcher_pid), launch_context$launcher_start_utc),
  list(v2r4_process_identity(runtime_buffers$os_guard_source,
                             as.integer(lease$model_pid[[1L]])),
       as.integer(lease$model_pid[[1L]]), lease$model_start_utc[[1L]]),
  list(v2r4_process_identity(runtime_buffers$os_guard_source,
                             as.integer(Sys.getpid())),
       as.integer(Sys.getpid()), validator_identity$start_utc),
  list(v2r4_process_identity(runtime_buffers$os_guard_source,
                             guard_state$guard_pid),
       guard_state$guard_pid, guard_state$guard_start_utc)
)) {
  v2r4_assert_process_identity_exact(
    identity_check[[1L]], identity_check[[2L]], identity_check[[3L]]
  )
}
v2r4_assert(!file.exists(cli$final_challenge_request) &&
              !file.exists(cli$final_challenge_response) &&
              !file.exists(control_complete_path),
            "VALIDATOR_CHALLENGE_OR_CONTROL_EXISTED_BEFORE_PRE_FINAL_READY")
pre_final_ready_fields <- c(
  "schema_version", "status", "authorization_effect", "dataset_id", "model_version",
  "run_id", "ready_nonce", "lease_id", "lease_manifest_sha256",
  "source_authority_payload_sha256", "external_launch_context_sha256",
  "runtime_post_receipt_sha256", "runtime_pre_digest", "runtime_post_digest",
  "runtime_pre_post_equal", "science_complete_sha256", "finalization_request_sha256",
  "validation_report_sha256", "output_manifest_sha256", "complete_sha256",
  "snapshot_tree_sha256", "launcher_pid", "launcher_start_utc", "model_pid",
  "model_start_utc", "model_request_sha256", "model_request_nonce",
  "model_identity_request_sha256", "validator_pid", "validator_start_utc",
  "validator_request_sha256", "validator_request_nonce",
  "validator_identity_request_sha256", "guard_pid", "guard_start_utc", "guard_nonce",
  "guard_public_key_sha256", "guard_ready_receipt_sha256", "guard_mutex_file_id",
  "emitter_role", "emitter_pid", "emitter_start_utc", "emitter_source_sha256",
  "issued_utc", "expires_utc", "e2_gate_status", "candidate_gate_eligible",
  "can_generate_or_upgrade_candidate", "formal_sequential_inference"
)
pre_final_ready_time <- Sys.time()
pre_final_ready_issued <- format(
  pre_final_ready_time, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"
)
pre_final_ready_expires <- format(
  pre_final_ready_time + 120, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"
)
pre_final_ready_nonce <- paste(
  sprintf("%02x", as.integer(openssl::rand_bytes(32))), collapse = ""
)
v2r4_assert(
  as.numeric(v2r4_parse_fixed6_time(runtime_post$value$checked_utc,
                                    "VALIDATOR_PRE_FINAL_RUNTIME_POST_CHECKED")) <=
    as.numeric(pre_final_ready_time) &&
    as.numeric(v2r4_parse_fixed6_time(science_complete$issued_utc,
                                      "VALIDATOR_PRE_FINAL_SCIENCE_ISSUED")) <=
    as.numeric(pre_final_ready_time),
  "VALIDATOR_PRE_FINAL_READY_INPUT_TIME_ORDER_INVALID"
)
pre_final_ready <- list(
  schema_version = "GSE205906_bulk_model_v2_r4_validator_pre_final_challenge_ready_v1",
  status = "VALIDATOR_STABLE_PRECHECKS_COMPLETE_AWAITING_ONE_TIME_CHALLENGE",
  authorization_effect = "UNTRUSTED_SCHEDULING_HINT_NOT_AUTHORIZATION",
  dataset_id = cfg$dataset_id,
  model_version = cfg$model_version,
  run_id = as.character(cli$run_id),
  ready_nonce = pre_final_ready_nonce,
  lease_id = lease$lease_id[[1L]],
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
  external_launch_context_sha256 = launch_context$context_sha256,
  runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
  runtime_pre_digest = launch_context$runtime_pre_digest,
  runtime_post_digest = runtime_post$value$runtime_post_digest,
  runtime_pre_post_equal = TRUE,
  science_complete_sha256 = science_complete_buffer$sha256,
  finalization_request_sha256 = finalization_request_buffer$sha256,
  validation_report_sha256 = report_buffer$sha256,
  output_manifest_sha256 = authorized_output$markers$output_manifest$sha256,
  complete_sha256 = authorized_output$markers$complete$sha256,
  snapshot_tree_sha256 = snapshot$tree_digest,
  launcher_pid = as.integer(launch_context$launcher_pid),
  launcher_start_utc = launch_context$launcher_start_utc,
  model_pid = as.integer(lease$model_pid[[1L]]),
  model_start_utc = lease$model_start_utc[[1L]],
  model_request_sha256 = model_request_buffer$sha256,
  model_request_nonce = model_request_bundle$request$request_nonce,
  model_identity_request_sha256 = model_identity_request_buffer$sha256,
  validator_pid = as.integer(Sys.getpid()),
  validator_start_utc = validator_identity$start_utc,
  validator_request_sha256 = child_request$buffer$sha256,
  validator_request_nonce = child_request$request$request_nonce,
  validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  guard_pid = guard_state$guard_pid,
  guard_start_utc = guard_state$guard_start_utc,
  guard_nonce = guard_state$guard_nonce,
  guard_public_key_sha256 = guard_state$guard_public_key_sha256,
  guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
  guard_mutex_file_id = guard_state$mutex_file_id,
  emitter_role = "VALIDATOR",
  emitter_pid = as.integer(Sys.getpid()),
  emitter_start_utc = validator_identity$start_utc,
  emitter_source_sha256 = script_runtime_buffer$sha256,
  issued_utc = pre_final_ready_issued,
  expires_utc = pre_final_ready_expires,
  e2_gate_status = cfg$gate$machine_status,
  candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)
v2r4_assert(identical(names(pre_final_ready), pre_final_ready_fields),
            "VALIDATOR_PRE_FINAL_READY_FIELD_ORDER_INVALID")
pre_final_ready_buffer <- v2r4_write_json_atomic(
  pre_final_ready, pre_final_challenge_ready_path, attempt_root
)
v2r4_assert(identical(pre_final_ready_buffer$canonical_path,
                      pre_final_challenge_ready_path),
            "VALIDATOR_PRE_FINAL_READY_ATOMIC_PATH_INVALID")

v2r4_wait_for_files(
  cli$final_challenge_request, final_wait,
  "VALIDATOR_FINAL_CHALLENGE_REQUEST_AFTER_READY"
)
v2r4_wait_for_files(
  cli$final_challenge_response,
  cfg$finalization$final_challenge_max_age_seconds,
  "VALIDATOR_FINAL_CHALLENGE_RESPONSE_AFTER_READY"
)
challenge_nonce_expected <- sub(
  "^final_challenge_([0-9a-f]{64})\\.request\\.json$", "\\1",
  basename(cli$final_challenge_request)
)
v2r4_assert(grepl("^[0-9a-f]{64}$", challenge_nonce_expected),
            "VALIDATOR_FINAL_CHALLENGE_PATH_NONCE_INVALID")
challenge <- v2r4_verify_final_challenge(
  cli$final_challenge_request, cli$final_challenge_response, cfg, guard_state,
  list(
    run_id = as.character(cli$run_id),
    lease_id = lease$lease_id[[1L]],
    lease_manifest_sha256 = lease_bundle$manifest$sha256,
    source_authority_payload_sha256 =
      launch_context$source_authorization_payload_sha256,
    launcher_pid = as.integer(launch_context$launcher_pid),
    launcher_start_utc = launch_context$launcher_start_utc,
    guard_owner_role = guard_state$owner_role,
    guard_owner_pid = as.integer(guard_state$owner_pid),
    guard_owner_start_utc = guard_state$owner_start_utc,
    guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
    finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
    finalization_public_key_captured_path =
      guard_state$finalization_public_key_captured_path,
    finalization_public_key_captured_bytes =
      guard_state$finalization_public_key_captured_bytes,
    finalization_public_key_captured_sha256 =
      guard_state$finalization_public_key_captured_sha256,
    model_pid = as.integer(lease$model_pid[[1L]]),
    model_start_utc = lease$model_start_utc[[1L]],
    validator_pid = as.integer(Sys.getpid()),
    validator_start_utc = validator_identity$start_utc,
    guard_pid = guard_state$guard_pid,
    guard_start_utc = guard_state$guard_start_utc,
    guard_nonce = guard_state$guard_nonce,
    guard_public_key_sha256 = guard_state$guard_public_key_sha256,
    guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
    guard_mutex_file_id = guard_state$mutex_file_id,
    model_request_sha256 = model_request_buffer$sha256,
    model_request_nonce = model_request_bundle$request$request_nonce,
    model_identity_request_sha256 = model_identity_request_buffer$sha256,
    validator_request_sha256 = child_request$buffer$sha256,
    validator_request_nonce = child_request$request$request_nonce,
    validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
    science_complete_sha256 = science_complete_buffer$sha256,
    finalization_request_sha256 = finalization_request_buffer$sha256,
    validation_report_sha256 = report_buffer$sha256,
    output_manifest_sha256 = authorized_output$markers$output_manifest$sha256,
    complete_sha256 = authorized_output$markers$complete$sha256,
    snapshot_tree_sha256 = snapshot$tree_digest,
    runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
    runtime_pre_digest = launch_context$runtime_pre_digest,
    runtime_post_digest = launch_context$runtime_pre_digest,
    runtime_pre_post_equal = TRUE,
    challenge_nonce = challenge_nonce_expected
  ), attempt_root
)
control_issued_time <- Sys.time()
challenge_issued_time <- v2r4_parse_fixed6_time(
  challenge$request_value$issued_utc, "VALIDATOR_CONTROL_CHALLENGE_ISSUED"
)
challenge_responded_time <- v2r4_parse_fixed6_time(
  challenge$response_value$responded_utc, "VALIDATOR_CONTROL_CHALLENGE_RESPONDED"
)
challenge_expires_time <- v2r4_parse_fixed6_time(
  challenge$request_value$expires_utc, "VALIDATOR_CONTROL_CHALLENGE_EXPIRES"
)
v2r4_assert(
  as.numeric(v2r4_parse_fixed6_time(pre_final_ready_issued,
                                    "VALIDATOR_PRE_FINAL_READY_ISSUED")) <=
    as.numeric(challenge_issued_time) &&
    as.numeric(challenge_responded_time) <= as.numeric(control_issued_time) &&
    as.numeric(control_issued_time) < as.numeric(challenge_expires_time),
  "VALIDATOR_CHALLENGE_RESPONSE_CONTROL_TIME_ORDER_INVALID"
)
control_issued_utc <- format(
  control_issued_time, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"
)
control_complete <- list(
  schema_version = cfg$finalization$control_complete_schema_version,
  status = "FINAL_CONTROL_VERIFIED_VALIDATOR_ALIVE_WAITING_FOR_ACK",
  dataset_id = cfg$dataset_id,
  model_version = cfg$model_version,
  run_id = as.character(cli$run_id),
  lease_id = lease$lease_id[[1L]],
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  launcher_pid = as.integer(launch_context$launcher_pid),
  launcher_start_utc = launch_context$launcher_start_utc,
  guard_owner_role = guard_state$owner_role,
  guard_owner_pid = as.integer(guard_state$owner_pid),
  guard_owner_start_utc = guard_state$owner_start_utc,
  guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
  finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
  finalization_public_key_captured_path =
    guard_state$finalization_public_key_captured_path,
  finalization_public_key_captured_bytes =
    guard_state$finalization_public_key_captured_bytes,
  finalization_public_key_captured_sha256 =
    guard_state$finalization_public_key_captured_sha256,
  model_pid = as.integer(lease$model_pid[[1L]]),
  model_start_utc = lease$model_start_utc[[1L]],
  validator_pid = as.integer(Sys.getpid()),
  validator_start_utc = validator_identity$start_utc,
  guard_pid = guard_state$guard_pid,
  guard_start_utc = guard_state$guard_start_utc,
  guard_nonce = guard_state$guard_nonce,
  guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
  guard_mutex_file_id = guard_state$mutex_file_id,
  model_request_sha256 = model_request_buffer$sha256,
  model_request_nonce = model_request_bundle$request$request_nonce,
  model_identity_request_sha256 = model_identity_request_buffer$sha256,
  validator_request_sha256 = child_request$buffer$sha256,
  validator_request_nonce = child_request$request$request_nonce,
  validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  science_complete_sha256 = science_complete_buffer$sha256,
  finalization_request_sha256 = finalization_request_buffer$sha256,
  validation_report_sha256 = report_buffer$sha256,
  output_manifest_sha256 = authorized_output$markers$output_manifest$sha256,
  complete_sha256 = authorized_output$markers$complete$sha256,
  runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
  runtime_pre_digest = launch_context$runtime_pre_digest,
  runtime_post_digest = runtime_post$value$runtime_post_digest,
  runtime_pre_post_equal = TRUE,
  final_challenge_request_sha256 = challenge$request$sha256,
  final_challenge_response_sha256 = challenge$response$sha256,
  final_challenge_payload_sha256 = challenge$payload_sha256,
  final_challenge_nonce = challenge$request_value$challenge_nonce,
  final_challenge_issued_utc = challenge$request_value$issued_utc,
  final_challenge_responded_utc = challenge$response_value$responded_utc,
  lease_issued_utc = lease$issued_utc[[1L]],
  lease_expires_utc = lease$expires_utc[[1L]],
  issued_utc = control_issued_utc,
  e2_gate_status = cfg$gate$machine_status,
  candidate_gate_eligible = FALSE,
  can_generate_or_upgrade_candidate = FALSE,
  formal_sequential_inference = FALSE
)
v2r4_assert(identical(names(control_complete), v2r4_control_complete_fields()),
            "CONTROL_COMPLETE_FIELD_ORDER_INVALID")
control_complete_buffer <- v2r4_write_json_atomic(
  control_complete, control_complete_path, attempt_root
)

ack_expected <- list(
  holder_principal = as.character(cli$holder_principal),
  run_id = as.character(cli$run_id),
  lease_id = lease$lease_id[[1L]],
  lease_manifest_sha256 = lease_bundle$manifest$sha256,
  source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
  model_source_sha256 = runtime_buffers$model_source$sha256,
  validator_source_sha256 = script_runtime_buffer$sha256,
  shared_library_sha256 = library_runtime_buffer$sha256,
  analysis_config_sha256 = config_buffer$sha256,
  external_launcher_sha256 = launch_context$bootstrap_launcher_sha256,
  os_guard_sha256 = runtime_buffers$os_guard_source$sha256,
  launcher_pid = as.integer(launch_context$launcher_pid),
  launcher_start_utc = launch_context$launcher_start_utc,
  guard_owner_role = guard_state$owner_role,
  guard_owner_pid = as.integer(guard_state$owner_pid),
  guard_owner_start_utc = guard_state$owner_start_utc,
  guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
  finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
  finalization_public_key_captured_path =
    guard_state$finalization_public_key_captured_path,
  finalization_public_key_captured_bytes =
    guard_state$finalization_public_key_captured_bytes,
  finalization_public_key_captured_sha256 =
    guard_state$finalization_public_key_captured_sha256,
  guard_pid = guard_state$guard_pid,
  guard_start_utc = guard_state$guard_start_utc,
  guard_nonce = guard_state$guard_nonce,
  guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
  output_manifest_sha256 = authorized_output$markers$output_manifest$sha256,
  complete_sha256 = authorized_output$markers$complete$sha256,
  validation_report_sha256 = report_buffer$sha256,
  validation_request_sha256 = finalization_request_buffer$sha256,
  snapshot_tree_sha256 = snapshot$tree_digest,
  runtime_pre_digest = launch_context$runtime_pre_digest,
  runtime_post_digest = runtime_post$value$runtime_post_digest,
  runtime_pre_post_equal = TRUE,
  runtime_post_receipt_sha256 = runtime_post$buffer$sha256,
  model_request_sha256 = model_request_buffer$sha256,
  model_request_nonce = model_request_bundle$request$request_nonce,
  model_identity_request_sha256 = model_identity_request_buffer$sha256,
  validator_request_sha256 = child_request$buffer$sha256,
  validator_request_nonce = child_request$request$request_nonce,
  validator_identity_request_sha256 = validator_identity_request_buffer$sha256,
  science_complete_sha256 = science_complete_buffer$sha256,
  control_complete_sha256 = control_complete_buffer$sha256,
  finalization_request_sha256 = finalization_request_buffer$sha256,
  final_challenge_request_sha256 = challenge$request$sha256,
  final_challenge_response_sha256 = challenge$response$sha256,
  final_challenge_payload_sha256 = challenge$payload_sha256,
  final_challenge_nonce = challenge$request_value$challenge_nonce,
  final_challenge_issued_utc = challenge$request_value$issued_utc,
  final_challenge_responded_utc = challenge$response_value$responded_utc,
  lease_issued_utc = lease$issued_utc[[1L]],
  lease_expires_utc = lease$expires_utc[[1L]],
  model_pid = as.integer(lease$model_pid[[1L]]),
  model_start_utc = lease$model_start_utc[[1L]],
  validator_pid = as.integer(Sys.getpid()),
  validator_start_utc = validator_identity$start_utc,
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
  launch_context$finalization_public_key_captured_sha256, cfg, ack_expected, final_wait
)
v2r4_write_consumer_receipt(
  validator_consumer_receipt_path, cfg, "VALIDATOR", ack, ack_expected, attempt_root
)
message(sprintf("PASS %s VALIDATION_CHECKS=%d STATUS=%s",
                cfg$model_version, nrow(report), cfg$status))
