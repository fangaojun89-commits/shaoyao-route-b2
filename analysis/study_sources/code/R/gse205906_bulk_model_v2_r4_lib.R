`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

v2r4_stop <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

v2r4_assert <- function(ok, fmt, ...) {
  if (!isTRUE(ok)) v2r4_stop(fmt, ...)
  invisible(TRUE)
}

v2r4_child_request_contract <- function(expected_role) {
  v2r4_assert(is.character(expected_role) && length(expected_role) == 1L &&
                !is.na(expected_role) && expected_role %in% c("MODEL", "VALIDATOR"),
              "CHILD_REQUEST_EXPECTED_ROLE_INVALID")
  common_fields <- c(
    "schema_version", "transport", "role", "purpose", "run_id",
    "request_nonce", "peer_request_nonce", "launch_context_sha256", "launcher_token_sha256",
    "runtime_pre_digest", "source_authorization_payload_sha256"
  )
  if (identical(expected_role, "MODEL")) {
    cli_fields <- c(
      "execution_mode", "attempt_root", "config", "authority_manifest",
      "authority_signature", "authority_public_key", "source_key_sha256",
      "snapshot_root", "lease_manifest", "lease_signature", "lease_public_key",
      "lease_key_sha256", "lease_request", "lease_identity_request",
      "peer_lease_identity_request", "state_root", "holder_principal",
      "output_dir", "final_ack", "finalization_public_key", "validator_request",
      "runtime_post_receipt", "science_complete", "control_complete",
      "finalization_request", "final_challenge_request", "final_challenge_response",
      "model_consumer_receipt", "validator_consumer_receipt"
    )
    schema_version <- "GSE205906_bulk_model_v2_r4_model_request_v1"
    purpose <- "EXECUTE_GSE205906_BULK_MODEL_V2_R4_SIGNED_MODE"
  } else {
    cli_fields <- c(
      "execution_mode", "attempt_root", "config", "authority_manifest",
      "authority_signature", "authority_public_key", "source_key_sha256",
      "snapshot_root", "lease_manifest", "lease_signature", "lease_public_key",
      "lease_key_sha256", "lease_request", "lease_identity_request",
      "peer_lease_identity_request", "state_root", "holder_principal",
      "output_authority_manifest", "output_authority_signature",
      "output_authority_public_key", "output_key_sha256", "output_dir", "report",
      "finalization_request", "model_request", "runtime_post_receipt",
      "science_complete", "control_complete", "final_challenge_request",
      "final_challenge_response", "final_ack", "finalization_public_key",
      "model_consumer_receipt", "validator_consumer_receipt"
    )
    schema_version <- "GSE205906_bulk_model_v2_r4_validator_request_v1"
    purpose <- "VALIDATE_GSE205906_BULK_MODEL_V2_R4_SIGNED_MODE"
  }
  list(
    schema_version = schema_version,
    transport = "CAPTURED_UTF8_NO_BOM_EXACT_JSON_SHA256_POSITIONAL_V1",
    role = expected_role,
    purpose = purpose,
    maximum_bytes = 65536L,
    exact_fields = c(common_fields, cli_fields),
    cli_fields = c("run_id", cli_fields),
    logical_fields = character()
  )
}

v2r4_child_request_utf8 <- function(x, label) {
  bytes <- v2r4_as_raw(x)
  v2r4_assert(length(bytes) > 1L, "%s_EMPTY_OR_TRUNCATED", label)
  v2r4_assert(!(length(bytes) >= 3L &&
                  identical(bytes[1:3], as.raw(c(0xef, 0xbb, 0xbf)))),
              "%s_UTF8_BOM_FORBIDDEN", label)
  v2r4_assert(identical(bytes[[length(bytes)]], as.raw(0x0a)) &&
                !any(bytes == as.raw(0x0d)) && !any(bytes == as.raw(0x00)),
              "%s_CANONICAL_TEXT_BYTES_INVALID", label)
  text <- rawToChar(bytes)
  validated <- iconv(text, from = "UTF-8", to = "UTF-8", sub = NA_character_)
  v2r4_assert(length(validated) == 1L && !is.na(validated) &&
                identical(charToRaw(enc2utf8(validated)), bytes),
              "%s_UTF8_ROUNDTRIP_INVALID", label)
  validated
}

v2r4_parse_child_request_buffer <- function(request_buffer, expected_sha256,
                                             expected_role) {
  v2r4_require_packages(c("digest", "jsonlite"))
  contract <- v2r4_child_request_contract(expected_role)
  v2r4_assert(inherits(request_buffer, "v2r4_captured_bytes"),
              "CHILD_REQUEST_CAPTURED_BUFFER_REQUIRED")
  v2r4_assert(is.character(expected_sha256) && length(expected_sha256) == 1L &&
                !is.na(expected_sha256) &&
                grepl("^[0-9a-f]{64}$", expected_sha256),
              "CHILD_REQUEST_EXPECTED_SHA256_INVALID")
  v2r4_assert(request_buffer$bytes > 0 &&
                request_buffer$bytes <= contract$maximum_bytes,
              "CHILD_REQUEST_SIZE_INVALID")
  v2r4_assert(identical(request_buffer$sha256, expected_sha256),
              "CHILD_REQUEST_SHA256_MISMATCH")
  request_text <- v2r4_child_request_utf8(request_buffer, "CHILD_REQUEST")
  request <- tryCatch(
    jsonlite::fromJSON(request_text, simplifyVector = FALSE),
    error = function(e) v2r4_stop("CHILD_REQUEST_JSON_INVALID:%s", conditionMessage(e))
  )
  v2r4_assert(is.list(request) && !is.null(names(request)) &&
                !anyDuplicated(names(request)) &&
                identical(names(request), contract$exact_fields),
              "CHILD_REQUEST_EXACT_FIELD_ORDER_INVALID")
  string_fields <- setdiff(contract$exact_fields, contract$logical_fields)
  valid_strings <- vapply(request[string_fields], function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(value)
  }, logical(1))
  v2r4_assert(all(valid_strings), "CHILD_REQUEST_STRING_FIELD_TYPE_OR_EMPTY_INVALID")
  v2r4_assert(identical(request$execution_mode, "SYNTHETIC") ||
                identical(request$execution_mode, "PRODUCTION"),
              "CHILD_REQUEST_EXECUTION_MODE_INVALID")
  v2r4_assert(identical(request$schema_version, contract$schema_version) &&
                identical(request$transport, contract$transport) &&
                identical(request$role, contract$role),
              "CHILD_REQUEST_ROLE_OR_SCHEMA_INVALID")
  v2r4_assert(identical(request$purpose, contract$purpose),
              "CHILD_REQUEST_PURPOSE_INVALID")
  v2r4_assert(grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", request$run_id),
              "CHILD_REQUEST_RUN_ID_LEXICAL_INVALID")
  for (field in c(
    "request_nonce", "peer_request_nonce", "launch_context_sha256", "launcher_token_sha256",
    "runtime_pre_digest", "source_authorization_payload_sha256"
  )) {
    v2r4_assert(grepl("^[0-9a-f]{64}$", request[[field]]),
                "CHILD_REQUEST_HASH_OR_NONCE_INVALID:%s", field)
  }
  cli <- request[contract$cli_fields]
  list(request = request, cli = cli, buffer = request_buffer, contract = contract)
}

v2r4_capture_child_request_args <- function(args, expected_role) {
  v2r4_assert(is.character(args) && length(args) == 2L && !anyNA(args) &&
                all(nzchar(args)),
              "CHILD_REQUEST_CLI_MUST_BE_EXACT_PATH_AND_BARE_SHA256")
  v2r4_assert(grepl("^[0-9a-f]{64}$", args[[2L]]),
              "CHILD_REQUEST_CLI_SHA256_NOT_BARE_LOWERCASE")
  request_buffer <- v2r4_capture_file(args[[1L]], "CHILD_REQUEST_RUNTIME")
  v2r4_parse_child_request_buffer(request_buffer, args[[2L]], expected_role)
}

v2r4_validate_child_request_bindings <- function(bundle, cfg, launch_context,
                                                  expected_role) {
  contract <- v2r4_child_request_contract(expected_role)
  v2r4_assert(is.list(bundle) && is.list(bundle$request) &&
                identical(bundle$contract$exact_fields, contract$exact_fields),
              "CHILD_REQUEST_BUNDLE_INVALID")
  request <- bundle$request
  transport <- cfg$external_runtime_authority$child_request_transport
  role_config <- if (identical(expected_role, "MODEL")) transport$model else transport$validator
  v2r4_assert(identical(request$schema_version, role_config$schema_version) &&
                identical(request$transport, transport$transport) &&
                identical(request$role, role_config$role) &&
                identical(request$purpose, role_config$purpose),
              "CHILD_REQUEST_CONFIG_BINDING_INVALID")
  launcher_token <- Sys.getenv("GSE205906_V2_R4_LAUNCHER_TOKEN", unset = "")
  v2r4_assert(grepl("^[0-9a-f]{64}$", launcher_token),
              "CHILD_REQUEST_LAUNCHER_TOKEN_UNAVAILABLE")
  expected_nonce <- tolower(digest::digest(charToRaw(paste(
    launcher_token, launch_context$context_sha256, expected_role, request$run_id,
    sep = "\n"
  )), algo = "sha256", serialize = FALSE))
  peer_role <- if (identical(expected_role, "MODEL")) "VALIDATOR" else "MODEL"
  expected_peer_nonce <- tolower(digest::digest(charToRaw(paste(
    launcher_token, launch_context$context_sha256, peer_role, request$run_id,
    sep = "\n"
  )), algo = "sha256", serialize = FALSE))
  v2r4_assert(identical(request$run_id, launch_context$run_id),
              "CHILD_REQUEST_RUN_BINDING_INVALID")
  v2r4_assert(identical(request$request_nonce, expected_nonce),
              "CHILD_REQUEST_REPLAY_BINDING_INVALID")
  v2r4_assert(identical(request$peer_request_nonce, expected_peer_nonce) &&
                !identical(request$request_nonce, request$peer_request_nonce),
              "CHILD_REQUEST_PEER_NONCE_BINDING_INVALID")
  v2r4_assert(identical(request$launch_context_sha256,
                        launch_context$context_sha256),
              "CHILD_REQUEST_CONTEXT_BINDING_INVALID")
  v2r4_assert(identical(request$launcher_token_sha256,
                        launch_context$launcher_token_sha256),
              "CHILD_REQUEST_TOKEN_BINDING_INVALID")
  v2r4_assert(identical(request$runtime_pre_digest,
                        launch_context$runtime_pre_digest) &&
                identical(request$source_authorization_payload_sha256,
                          launch_context$source_authorization_payload_sha256),
              "CHILD_REQUEST_DIGEST_BINDING_INVALID")
  v2r4_assert(identical(request$execution_mode, launch_context$mode) &&
                identical(request$attempt_root, launch_context$attempt_root),
              "CHILD_REQUEST_SIGNED_MODE_OR_ATTEMPT_BINDING_INVALID")
  v2r4_assert(identical(request$holder_principal,
                        launch_context$holder_principal) &&
                identical(request$state_root, launch_context$state_root) &&
                identical(request$output_dir, launch_context$output_dir),
              "CHILD_REQUEST_CONTROL_PATH_BINDING_INVALID")
  if (identical(expected_role, "MODEL")) {
    v2r4_assert(identical(request$final_ack,
                          launch_context$final_ack_envelope_path),
                "CHILD_REQUEST_MODEL_FINAL_ACK_BINDING_INVALID")
    v2r4_assert(identical(request$validator_request,
                          launch_context$validator_request_path),
                "CHILD_REQUEST_MODEL_PEER_PATH_BINDING_INVALID")
  } else {
    v2r4_assert(identical(request$report,
                          launch_context$validation_report_path) &&
                  identical(request$finalization_request,
                            launch_context$finalization_request_path),
                "CHILD_REQUEST_VALIDATOR_OUTPUT_BINDING_INVALID")
    v2r4_assert(identical(request$model_request,
                          launch_context$model_request_path) &&
                  identical(request$final_ack,
                            launch_context$final_ack_envelope_path),
                "CHILD_REQUEST_VALIDATOR_PEER_OR_ACK_PATH_BINDING_INVALID")
  }
  control_path_fields <- c(
    runtime_post_receipt = "runtime_post_receipt_path",
    science_complete = "science_complete_path",
    control_complete = "control_complete_path",
    finalization_request = "finalization_request_path",
    final_challenge_request = "final_challenge_request_path",
    final_challenge_response = "final_challenge_response_path",
    model_consumer_receipt = "model_consumer_receipt_path",
    validator_consumer_receipt = "validator_consumer_receipt_path"
  )
  v2r4_assert(all(vapply(names(control_path_fields), function(field) {
    identical(request[[field]], launch_context[[control_path_fields[[field]]]])
  }, logical(1))), "CHILD_REQUEST_FINAL_CONTROL_PATH_BINDING_INVALID")
  invisible(TRUE)
}

v2r4_script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  v2r4_assert(length(arg) == 1L, "SCRIPT_PATH_NOT_IDENTIFIABLE")
  normalizePath(sub("^--file=", "", arg), winslash = "/", mustWork = TRUE)
}

v2r4_project_root <- function(script_path = v2r4_script_path()) {
  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
}

v2r4_norm_existing <- function(path) {
  v2r4_assert(length(path) == 1L && !is.na(path) && nzchar(path), "EMPTY_PATH")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

v2r4_norm_target <- function(path) {
  v2r4_assert(length(path) == 1L && !is.na(path) && nzchar(path), "EMPTY_PATH")
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  paste0(parent, "/", basename(path))
}

v2r4_is_within <- function(path, root, allow_equal = FALSE) {
  p <- tolower(gsub("/+$", "", path))
  r <- tolower(gsub("/+$", "", root))
  (identical(p, r) && allow_equal) || startsWith(p, paste0(r, "/"))
}

v2r4_require_outside <- function(path, project_root, label) {
  v2r4_assert(!v2r4_is_within(path, project_root, allow_equal = TRUE),
              "%s_MUST_BE_OUTSIDE_PROJECT:%s", label, path)
}

v2r4_require_within <- function(path, root, label, allow_equal = FALSE) {
  v2r4_assert(v2r4_is_within(path, root, allow_equal = allow_equal),
              "%s_ESCAPES_REQUIRED_ROOT:%s", label, path)
}

v2r4_reject_symlink <- function(path, label) {
  v2r4_assert(!nzchar(Sys.readlink(path)), "%s_SYMLINK_NOT_PERMITTED:%s", label, path)
  invisible(TRUE)
}

v2r4_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  v2r4_assert(!length(missing), "MISSING_R_PACKAGES:%s", paste(missing, collapse = ","))
}

v2r4_hmac_sha256 <- function(key, bytes) {
  v2r4_assert(length(key) == 1L && grepl("^[0-9a-f]{64}$", key),
              "LAUNCHER_TOKEN_LEXICAL_INVALID")
  v2r4_assert(is.raw(bytes), "LAUNCH_CONTEXT_RAW_REQUIRED")
  tolower(digest::hmac(key = key, object = bytes, algo = "sha256", serialize = FALSE))
}

v2r4_require_external_launch_context <- function(cfg, expected_role) {
  v2r4_require_packages(c("digest", "jsonlite"))
  token <- Sys.getenv("GSE205906_V2_R4_LAUNCHER_TOKEN", unset = "")
  context_path <- Sys.getenv("GSE205906_V2_R4_LAUNCH_CONTEXT", unset = "")
  mac_path <- Sys.getenv("GSE205906_V2_R4_LAUNCH_CONTEXT_HMAC", unset = "")
  v2r4_assert(identical(Sys.getenv("GSE205906_V2_R4_EXTERNAL_LAUNCH", unset = ""), "1") &&
                grepl("^[0-9a-f]{64}$", token) && nzchar(context_path) && nzchar(mac_path),
              "V2_R4_EXTERNAL_LAUNCH_REQUIRED_BEFORE_CANDIDATE_EXECUTION")
  context_buffer <- v2r4_capture_file(context_path, "EXTERNAL_LAUNCH_CONTEXT")
  mac_buffer <- v2r4_capture_file(mac_path, "EXTERNAL_LAUNCH_CONTEXT_HMAC")
  observed_mac <- trimws(v2r4_parse_text_raw(mac_buffer, "EXTERNAL_LAUNCH_CONTEXT_HMAC"))
  v2r4_assert(grepl("^[0-9a-f]{64}$", observed_mac) &&
                identical(observed_mac, v2r4_hmac_sha256(token, context_buffer$raw)),
              "EXTERNAL_LAUNCH_CONTEXT_HMAC_INVALID")
  context <- v2r4_read_json_raw(context_buffer, "EXTERNAL_LAUNCH_CONTEXT")
  required <- c(
    "schema_version", "mode", "project_root", "bundle_root", "attempt_root", "state_root",
    "run_id", "holder_principal", "output_dir", "nonce_ledger_root",
    "validation_report_path",
    "finalization_request_path", "final_ack_envelope_path", "runtime_post_receipt_path",
    "science_complete_path", "control_complete_path", "final_challenge_request_path",
    "final_challenge_response_path", "model_consumer_receipt_path",
    "validator_consumer_receipt_path", "runtime_final_receipt_path",
    "guard_ready_receipt_path", "guard_stop_request_path", "guard_stopped_receipt_path",
    "guard_stop_hmac_key_sha256", "model_request_path", "validator_request_path",
    "model_lease_identity_request_path", "validator_lease_identity_request_path",
    "combined_lease_request_path",
    "source_authorization_envelope_path",
    "source_public_key_path", "source_authorization_payload_sha256",
    "source_authorization_envelope_sha256", "source_public_key_sha256",
    "source_snapshot_manifest_sha256", "source_snapshot_public_key_sha256",
    "lease_public_key_sha256", "output_public_key_source_path",
    "output_public_key_source_bytes", "output_public_key_sha256",
    "output_public_key_captured_path", "output_public_key_captured_bytes",
    "output_public_key_captured_sha256", "finalization_public_key_source_path",
    "finalization_public_key_source_bytes", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "captured_manifest_path",
    "captured_manifest_sha256", "allowlist_sha256", "bootstrap_launcher_sha256",
    "production_control_manifest_sha256", "sealed_production_runner_sha256",
    "execution_identity_manifest_path", "execution_identity_manifest_sha256",
    "runtime_file_manifest_path", "runtime_file_manifest_sha256", "runtime_lock_path",
    "runtime_lock_sha256", "package_sbom_path", "package_sbom_sha256", "runtime_pre_digest",
    "runtime_semantic_digest", "pwsh_path", "pwsh_sha256", "run_r_path", "run_r_sha256",
    "rscript_path", "rscript_sha256", "model_source_sha256", "validator_source_sha256",
    "library_source_sha256", "config_path", "config_sha256", "os_guard_path",
    "os_guard_sha256", "finalizer_sha256", "final_ack_wait_seconds",
    "launcher_pid", "launcher_start_utc", "launcher_token_sha256", "global_mutex_name",
    "created_utc", "permitted_child_roles"
  )
  v2r4_assert(identical(sort(names(context)), sort(required)),
              "EXTERNAL_LAUNCH_CONTEXT_SCHEMA_INVALID")
  v2r4_assert(identical(context$schema_version,
                        "GSE205906_bulk_model_v2_r4_external_launch_context_v1") &&
                context$mode %in% c("SYNTHETIC", "PRODUCTION") &&
                identical(context$global_mutex_name,
                          cfg$external_runtime_authority$global_cooperative_mutex) &&
                identical(context$launcher_token_sha256,
                          tolower(digest::digest(charToRaw(token), algo = "sha256",
                                                 serialize = FALSE))) &&
                expected_role %in% unlist(context$permitted_child_roles, use.names = FALSE) &&
                identical(context$final_ack_wait_seconds, 900L),
              "EXTERNAL_LAUNCH_CONTEXT_FIXED_BINDING_INVALID")
  v2r4_assert(isTRUE(cfg$governance$production_execution_authorized) &&
                isTRUE(cfg$governance$production_validation_authorized) &&
                isTRUE(cfg$governance$production_entrypoint_enabled) &&
                isTRUE(cfg$governance$production_validation_entrypoint_enabled) &&
                identical(cfg$governance$production_static_authorization, FALSE) &&
                isTRUE(cfg$governance$authorization_is_per_run_only),
              "V2_R4_PRODUCTION_ENTRYPOINT_CONTRACT_INVALID")
  v2r4_assert_utc_fixed6(c(context$launcher_start_utc, context$created_utc),
                          "EXTERNAL_LAUNCH_CONTEXT_TIMES")
  for (field in c(
    "source_authorization_payload_sha256", "source_authorization_envelope_sha256",
    "source_public_key_sha256", "source_snapshot_manifest_sha256",
    "source_snapshot_public_key_sha256", "lease_public_key_sha256",
    "output_public_key_sha256", "output_public_key_captured_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_sha256",
    "captured_manifest_sha256", "allowlist_sha256",
    "bootstrap_launcher_sha256", "production_control_manifest_sha256",
    "sealed_production_runner_sha256", "execution_identity_manifest_sha256",
    "runtime_file_manifest_sha256", "runtime_lock_sha256", "package_sbom_sha256",
    "runtime_pre_digest", "runtime_semantic_digest", "pwsh_sha256", "run_r_sha256",
    "rscript_sha256", "model_source_sha256", "validator_source_sha256",
    "library_source_sha256", "config_sha256", "os_guard_sha256", "finalizer_sha256",
    "guard_stop_hmac_key_sha256"
  )) {
    v2r4_assert(grepl("^[0-9a-f]{64}$", context[[field]]),
                "EXTERNAL_LAUNCH_CONTEXT_HASH_INVALID:%s", field)
  }
  v2r4_assert(is.numeric(context$finalization_public_key_source_bytes) &&
                length(context$finalization_public_key_source_bytes) == 1L &&
                !is.na(context$finalization_public_key_source_bytes) &&
                context$finalization_public_key_source_bytes > 0 &&
                is.numeric(context$finalization_public_key_captured_bytes) &&
                length(context$finalization_public_key_captured_bytes) == 1L &&
                !is.na(context$finalization_public_key_captured_bytes) &&
                identical(as.numeric(context$finalization_public_key_captured_bytes),
                          as.numeric(context$finalization_public_key_source_bytes)) &&
                identical(context$finalization_public_key_captured_sha256,
                          context$finalization_public_key_sha256),
              "EXTERNAL_LAUNCH_FINALIZATION_KEY_CAPTURE_IDENTITY_INVALID")
  v2r4_assert(is.numeric(context$output_public_key_source_bytes) &&
                length(context$output_public_key_source_bytes) == 1L &&
                !is.na(context$output_public_key_source_bytes) &&
                context$output_public_key_source_bytes > 0 &&
                is.numeric(context$output_public_key_captured_bytes) &&
                length(context$output_public_key_captured_bytes) == 1L &&
                !is.na(context$output_public_key_captured_bytes) &&
                identical(as.numeric(context$output_public_key_captured_bytes),
                          as.numeric(context$output_public_key_source_bytes)) &&
                identical(context$output_public_key_captured_sha256,
                          context$output_public_key_sha256) &&
                !identical(context$output_public_key_sha256,
                           context$finalization_public_key_sha256),
              "EXTERNAL_LAUNCH_OUTPUT_KEY_CAPTURE_OR_SEPARATION_INVALID")
  output_source_key_path <- v2r4_norm_existing(context$output_public_key_source_path)
  output_captured_key_path <- v2r4_norm_existing(context$output_public_key_captured_path)
  v2r4_require_within(output_captured_key_path,
                      v2r4_norm_existing(context$attempt_root),
                      "OUTPUT_PUBLIC_KEY_CAPTURED_COPY")
  v2r4_assert(!identical(output_source_key_path, output_captured_key_path),
              "OUTPUT_PUBLIC_KEY_SOURCE_AND_COPY_PATH_COLLISION")
  source_key_path <- v2r4_norm_existing(context$finalization_public_key_source_path)
  captured_key_path <- v2r4_norm_existing(context$finalization_public_key_captured_path)
  v2r4_require_within(captured_key_path, v2r4_norm_existing(context$attempt_root),
                      "FINALIZATION_PUBLIC_KEY_CAPTURED_COPY")
  v2r4_assert(!identical(source_key_path, captured_key_path),
              "FINALIZATION_PUBLIC_KEY_SOURCE_AND_COPY_PATH_COLLISION")
  attempt_root <- v2r4_norm_existing(context$attempt_root)
  nonce_ledger_root <- v2r4_norm_existing(context$nonce_ledger_root)
  if (identical(context$mode, "PRODUCTION")) {
    v2r4_assert(!v2r4_is_within(nonce_ledger_root, attempt_root, allow_equal = TRUE) &&
                  !v2r4_is_within(attempt_root, nonce_ledger_root, allow_equal = TRUE),
                "PRODUCTION_NONCE_LEDGER_MUST_BE_PERSISTENT_OUTSIDE_ATTEMPT")
    v2r4_require_outside(nonce_ledger_root,
                         v2r4_norm_existing(context$project_root),
                         "PRODUCTION_NONCE_LEDGER")
  } else {
    v2r4_require_within(nonce_ledger_root, attempt_root,
                        "SYNTHETIC_NONCE_LEDGER", allow_equal = FALSE)
  }
  context$context_sha256 <- context_buffer$sha256
  context$context_hmac_sha256 <- observed_mac
  context
}

v2r4_capture_file <- function(path, label, required_root = NULL) {
  canonical <- v2r4_norm_existing(path)
  if (!is.null(required_root)) {
    v2r4_require_within(canonical, v2r4_norm_existing(required_root), label)
  }
  v2r4_reject_symlink(canonical, label)
  info <- file.info(canonical)
  v2r4_assert(nrow(info) == 1L && !is.na(info$size) && !isTRUE(info$isdir),
              "%s_NOT_REGULAR_FILE:%s", label, canonical)
  con <- file(canonical, open = "rb")
  on.exit(close(con), add = TRUE)
  remaining <- as.numeric(info$size)
  parts <- list()
  while (remaining > 0) {
    part <- readBin(con, what = "raw", n = as.integer(min(remaining, 64 * 1024^2)))
    v2r4_assert(length(part) > 0L, "%s_SHORT_READ:%s", label, canonical)
    parts[[length(parts) + 1L]] <- part
    remaining <- remaining - length(part)
  }
  payload <- if (length(parts)) do.call(c, parts) else raw()
  v2r4_assert(as.numeric(length(payload)) == as.numeric(info$size),
              "%s_BYTES_CHANGED_DURING_CAPTURE:%s", label, canonical)
  structure(list(
    raw = payload,
    canonical_path = canonical,
    bytes = as.numeric(length(payload)),
    sha256 = tolower(digest::digest(payload, algo = "sha256", serialize = FALSE)),
    md5 = tolower(digest::digest(payload, algo = "md5", serialize = FALSE)),
    label = label
  ), class = "v2r4_captured_bytes")
}

v2r4_as_raw <- function(x) {
  if (inherits(x, "v2r4_captured_bytes")) return(x$raw)
  v2r4_assert(is.raw(x), "IMMUTABLE_BUFFER_REQUIRED")
  x
}

v2r4_sha256_raw <- function(x) {
  tolower(digest::digest(v2r4_as_raw(x), algo = "sha256", serialize = FALSE))
}

v2r4_md5_raw <- function(x) {
  tolower(digest::digest(v2r4_as_raw(x), algo = "md5", serialize = FALSE))
}

v2r4_parse_text_raw <- function(x, label) {
  bytes <- v2r4_as_raw(x)
  v2r4_assert(!any(bytes == as.raw(0L)), "%s_CONTAINS_NUL", label)
  enc2utf8(rawToChar(bytes))
}

v2r4_read_tsv_raw <- function(x, label = "TSV_BUFFER", gzip = FALSE,
                             character_columns = character()) {
  bytes <- v2r4_as_raw(x)
  if (gzip) {
    bytes <- tryCatch(memDecompress(bytes, type = "gzip"),
                      error = function(e) v2r4_stop("%s_GZIP_INVALID:%s", label, conditionMessage(e)))
  }
  character_columns <- as.character(character_columns)
  v2r4_assert(!anyNA(character_columns) && !anyDuplicated(character_columns),
              "%s_CHARACTER_COLUMN_SPEC_INVALID", label)
  column_classes <- if (length(character_columns)) list(character = character_columns) else NULL
  data.table::fread(text = v2r4_parse_text_raw(bytes, label), sep = "\t",
                    na.strings = c("NA", ""), showProgress = FALSE,
                    colClasses = column_classes)
}

v2r4_read_json_raw <- function(x, label = "JSON_BUFFER") {
  jsonlite::fromJSON(v2r4_parse_text_raw(x, label), simplifyVector = FALSE)
}

v2r4_read_yaml_raw <- function(x, label = "YAML_BUFFER") {
  yaml::yaml.load(v2r4_parse_text_raw(x, label), eval.expr = FALSE)
}

v2r4_crc32_le <- function(bytes) {
  v2r4_assert(is.raw(bytes), "CRC32_RAW_REQUIRED")
  hex <- tolower(digest::digest(bytes, algo = "crc32", serialize = FALSE))
  v2r4_assert(grepl("^[0-9a-f]{8}$", hex), "CRC32_RESULT_INVALID")
  big_endian <- as.raw(strtoi(substring(hex, c(1L, 3L, 5L, 7L),
                                        c(2L, 4L, 6L, 8L)), base = 16L))
  rev(big_endian)
}

v2r4_u32_le <- function(bytes, label) {
  v2r4_assert(is.raw(bytes) && length(bytes) == 4L, "%s_U32_WIDTH_INVALID", label)
  sum(as.numeric(as.integer(bytes)) * 256^(0:3))
}

v2r4_xz_vli <- function(bytes, position, limit, label, field) {
  value <- 0
  multiplier <- 1
  for (width in seq_len(9L)) {
    v2r4_assert(position <= limit, "%s_XZ_%s_VLI_TRUNCATED", label, field)
    octet <- as.integer(bytes[[position]])
    position <- position + 1L
    payload <- bitwAnd(octet, 0x7F)
    value <- value + as.numeric(payload) * multiplier
    v2r4_assert(is.finite(value) && value <= 2^53 - 1,
                "%s_XZ_%s_VLI_EXCEEDS_EXACT_INTEGER_DOMAIN", label, field)
    if (bitwAnd(octet, 0x80) == 0L) {
      v2r4_assert(width == 1L || payload != 0L,
                  "%s_XZ_%s_VLI_NONCANONICAL", label, field)
      return(list(value = value, next_position = position, width = width))
    }
    v2r4_assert(width < 9L, "%s_XZ_%s_VLI_TOO_LONG", label, field)
    multiplier <- multiplier * 128
  }
  v2r4_stop("%s_XZ_%s_VLI_INVALID", label, field)
}

v2r4_validate_xz_container <- function(bytes, label) {
  magic <- as.raw(c(0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00))
  flags <- as.raw(c(0x00, 0x01))
  n_bytes <- as.numeric(length(bytes))
  v2r4_assert(n_bytes >= 36, "%s_XZ_LENGTH_INVALID", label)
  v2r4_assert(identical(bytes[seq_len(6L)], magic), "%s_XZ_MAGIC_INVALID", label)
  v2r4_assert(identical(bytes[7:8], flags), "%s_XZ_STREAM_FLAGS_INVALID", label)
  v2r4_assert(identical(bytes[9:12], v2r4_crc32_le(bytes[7:8])),
              "%s_XZ_HEADER_CRC_INVALID", label)

  footer_start <- length(bytes) - 11L
  footer <- bytes[footer_start:length(bytes)]
  v2r4_assert(identical(footer[11:12], as.raw(c(0x59, 0x5A))),
              "%s_XZ_FOOTER_MAGIC_INVALID", label)
  v2r4_assert(identical(footer[9:10], flags), "%s_XZ_FOOTER_FLAGS_INVALID", label)
  v2r4_assert(identical(footer[1:4], v2r4_crc32_le(footer[5:10])),
              "%s_XZ_FOOTER_CRC_INVALID", label)
  backward_size <- v2r4_u32_le(footer[5:8], paste0(label, "_XZ_BACKWARD_SIZE"))
  index_size <- 4 * (backward_size + 1)
  v2r4_assert(is.finite(index_size) && index_size >= 8 &&
                index_size == floor(index_size) && index_size <= n_bytes - 24,
              "%s_XZ_BACKWARD_SIZE_INVALID", label)
  index_start <- n_bytes - index_size - 11
  index_end <- n_bytes - 12
  v2r4_assert(index_start >= 13 && index_start <= index_end,
              "%s_XZ_INDEX_SPAN_INVALID", label)
  index_bytes <- bytes[index_start:index_end]
  index_content <- index_bytes[seq_len(length(index_bytes) - 4L)]
  index_crc <- index_bytes[(length(index_bytes) - 3L):length(index_bytes)]
  v2r4_assert(identical(index_crc, v2r4_crc32_le(index_content)),
              "%s_XZ_INDEX_CRC_INVALID", label)
  v2r4_assert(length(index_content) >= 2L && index_content[[1L]] == as.raw(0x00),
              "%s_XZ_INDEX_INDICATOR_INVALID", label)

  cursor <- 2L
  count_vli <- v2r4_xz_vli(index_content, cursor, length(index_content), label,
                           "INDEX_RECORD_COUNT")
  record_count <- count_vli$value
  cursor <- count_vli$next_position
  v2r4_assert(record_count >= 1 && record_count <= .Machine$integer.max &&
                record_count <= floor((length(index_content) - cursor + 1L) / 2L),
              "%s_XZ_INDEX_RECORD_COUNT_INVALID", label)
  record_count <- as.integer(record_count)
  unpadded_sizes <- numeric(record_count)
  uncompressed_sizes <- numeric(record_count)
  for (record_index in seq_len(record_count)) {
    unpadded <- v2r4_xz_vli(index_content, cursor, length(index_content), label,
                            sprintf("RECORD_%d_UNPADDED_SIZE", record_index))
    cursor <- unpadded$next_position
    uncompressed <- v2r4_xz_vli(index_content, cursor, length(index_content), label,
                                sprintf("RECORD_%d_UNCOMPRESSED_SIZE", record_index))
    cursor <- uncompressed$next_position
    v2r4_assert(unpadded$value > 4 && uncompressed$value > 0,
                "%s_XZ_INDEX_RECORD_DOMAIN_INVALID:%d", label, record_index)
    unpadded_sizes[[record_index]] <- unpadded$value
    uncompressed_sizes[[record_index]] <- uncompressed$value
  }
  padding_n <- length(index_content) - cursor + 1L
  v2r4_assert(padding_n >= 0L && padding_n <= 3L,
              "%s_XZ_INDEX_PADDING_WIDTH_INVALID", label)
  if (padding_n > 0L) {
    v2r4_assert(all(index_content[cursor:length(index_content)] == as.raw(0x00)),
                "%s_XZ_INDEX_PADDING_INVALID", label)
  }

  padded_sizes <- ceiling(unpadded_sizes / 4) * 4
  v2r4_assert(sum(padded_sizes) == index_start - 13,
              "%s_XZ_BLOCK_SPAN_INVALID", label)
  block_cursor <- 13
  for (record_index in seq_len(record_count)) {
    padding_n <- padded_sizes[[record_index]] - unpadded_sizes[[record_index]]
    if (padding_n > 0) {
      padding_start <- block_cursor + unpadded_sizes[[record_index]] - 4
      padding_end <- padding_start + padding_n - 1
      v2r4_assert(padding_start >= block_cursor && padding_end < index_start &&
                    all(bytes[padding_start:padding_end] == as.raw(0x00)),
                  "%s_XZ_BLOCK_PADDING_INVALID:%d", label, record_index)
    }
    block_cursor <- block_cursor + padded_sizes[[record_index]]
  }
  v2r4_assert(block_cursor == index_start, "%s_XZ_BLOCK_CURSOR_INVALID", label)
  total_uncompressed <- sum(uncompressed_sizes)
  v2r4_assert(is.finite(total_uncompressed) && total_uncompressed <= 2^53 - 1,
              "%s_XZ_UNCOMPRESSED_SIZE_DOMAIN_INVALID", label)
  list(record_count = record_count, total_uncompressed = total_uncompressed,
       index_start = index_start, index_size = index_size)
}

v2r4_gzip_drain <- function(bytes, label) {
  state <- new.env(parent = emptyenv())
  state$con <- NULL
  on.exit({
    if (!is.null(state$con)) try(suppressWarnings(close(state$con)), silent = TRUE)
  }, add = TRUE)
  drained <- tryCatch(
    withCallingHandlers({
      state$con <- rawConnection(bytes, open = "rb")
      state$con <- gzcon(state$con, level = 6, allowNonCompressed = FALSE, text = FALSE)
      chunks <- list()
      repeat {
        chunk <- readBin(state$con, what = "raw", n = 1024L * 1024L)
        if (!length(chunk)) break
        chunks[[length(chunks) + 1L]] <- chunk
      }
      close(state$con)
      state$con <- NULL
      if (length(chunks)) do.call(c, chunks) else raw()
    }, warning = function(w) stop(simpleError(conditionMessage(w)))),
    error = identity
  )
  if (inherits(drained, "error")) {
    v2r4_stop("%s_GZIP_DRAIN_FAILED:%s", label, conditionMessage(drained))
  }
  drained
}

v2r4_mem_decompress_exact <- function(bytes, compression, label) {
  decompressed <- tryCatch(
    withCallingHandlers(
      memDecompress(bytes, type = compression),
      warning = function(w) stop(simpleError(conditionMessage(w)))
    ),
    error = identity
  )
  if (inherits(decompressed, "error")) {
    v2r4_stop("%s_%s_DECOMPRESSION_FAILED:%s", label, toupper(compression),
              conditionMessage(decompressed))
  }
  decompressed
}

v2r4_unserialize_full <- function(decompressed, label) {
  con <- rawConnection(decompressed, open = "rb")
  on.exit(close(con), add = TRUE)
  value <- tryCatch(unserialize(con), error = identity)
  if (inherits(value, "error")) {
    v2r4_stop("%s_UNSERIALIZE_FAILED:%s", label, conditionMessage(value))
  }
  consumed <- seek(con, where = NA, origin = "start")
  v2r4_assert(identical(as.numeric(consumed), as.numeric(length(decompressed))),
              "%s_SERIALIZATION_TRAILING_BYTES:consumed=%s:total=%s",
              label, as.character(consumed), as.character(length(decompressed)))
  value
}

v2r4_read_rds_raw <- function(x, label = "RDS_BUFFER", compression) {
  bytes <- v2r4_as_raw(x)
  compression <- as.character(compression)
  v2r4_assert(length(compression) == 1L && !is.na(compression) &&
                compression %in% c("xz", "gzip"),
              "%s_COMPRESSION_CONTRACT_INVALID", label)
  if (identical(compression, "xz")) {
    container <- v2r4_validate_xz_container(bytes, label)
    decompressed <- v2r4_mem_decompress_exact(bytes, "xz", label)
    v2r4_assert(as.numeric(length(decompressed)) == container$total_uncompressed,
                "%s_XZ_UNCOMPRESSED_SIZE_MISMATCH:index=%s:observed=%s", label,
                as.character(container$total_uncompressed), as.character(length(decompressed)))
  } else {
    v2r4_assert(length(bytes) >= 18L, "%s_GZIP_LENGTH_INVALID", label)
    v2r4_assert(identical(bytes[1:3], as.raw(c(0x1F, 0x8B, 0x08))),
                "%s_GZIP_MAGIC_INVALID", label)
    v2r4_assert(bytes[[4L]] == as.raw(0x00), "%s_GZIP_FLAGS_INVALID", label)
    drained <- v2r4_gzip_drain(bytes, label)
    v2r4_assert(as.numeric(length(drained)) < 2^32,
                "%s_GZIP_ISIZE_DOMAIN_INVALID", label)
    trailer <- bytes[(length(bytes) - 7L):length(bytes)]
    v2r4_assert(identical(trailer[1:4], v2r4_crc32_le(drained)),
                "%s_GZIP_TRAILER_CRC_INVALID", label)
    isize <- v2r4_u32_le(trailer[5:8], paste0(label, "_GZIP_ISIZE"))
    v2r4_assert(isize == as.numeric(length(drained)),
                "%s_GZIP_ISIZE_INVALID:trailer=%s:observed=%s", label,
                as.character(isize), as.character(length(drained)))
    decompressed <- v2r4_mem_decompress_exact(bytes, "gzip", label)
    v2r4_assert(identical(decompressed, drained),
                "%s_GZIP_DECODER_DISAGREEMENT", label)
  }
  v2r4_unserialize_full(decompressed, label)
}

v2r4_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

v2r4_add_claim_flags <- function(x, analysis_role) {
  out <- data.table::as.data.table(x)
  role_value <- as.character(analysis_role)
  out[, `:=`(
    analysis_role = role_value,
    candidate_gate_eligible = FALSE,
    can_generate_or_upgrade_candidate = FALSE,
    formal_sequential_inference = FALSE
  )]
  out
}

v2r4_assert_deseq_contract <- function(cfg) {
  contract <- cfg$analysis$DESeq2
  allowed_actual <- unlist(contract$allowed_actual_fit_types_after_DESeq,
                           use.names = FALSE)
  v2r4_assert(identical(contract$requested_fit_type_before_fit, "parametric") &&
                identical(allowed_actual, c("parametric", "local")) &&
                identical(contract$automatic_local_substitution_message_token,
                          "a local regression fit was automatically substituted") &&
                identical(contract$manual_gene_wise_fallback_allowed, FALSE) &&
                identical(contract$fit_failure_actual_type_policy,
                          "NO_ACTUAL_FIT_TYPE_ON_FAILURE") &&
                isTRUE(contract$capture_actual_fit_type_after_DESeq) &&
                isTRUE(contract$capture_fallback_messages) &&
                identical(contract$independent_filtering, FALSE) &&
                identical(contract$cooks_cutoff, FALSE),
              "CONFIG_DESEQ2_FIT_PROVENANCE_CONTRACT_INVALID")
  invisible(TRUE)
}

v2r4_collapse_condition_messages <- function(messages) {
  values <- as.character(messages)
  values <- trimws(gsub(" +", " ", gsub("[\t\r\n]+", " ", values)))
  values <- unique(values[!is.na(values) & nzchar(values)])
  if (length(values)) paste(values, collapse = " | ") else "NONE"
}

v2r4_assert_deseq_fit_provenance <- function(actual_fit_type, fit_messages, cfg,
                                             label = "DESEQ2") {
  v2r4_assert_deseq_contract(cfg)
  actual <- as.character(actual_fit_type)
  messages <- as.character(fit_messages)
  allowed <- unlist(cfg$analysis$DESeq2$allowed_actual_fit_types_after_DESeq,
                    use.names = FALSE)
  token <- cfg$analysis$DESeq2$automatic_local_substitution_message_token
  v2r4_assert(length(actual) == 1L && !is.na(actual) && actual %in% allowed,
              "%s_ACTUAL_FIT_TYPE_NOT_ALLOWED", label)
  v2r4_assert(length(messages) == 1L && !is.na(messages) && nzchar(messages),
              "%s_FIT_MESSAGES_INVALID", label)
  reports_local_substitution <- grepl(token, messages, fixed = TRUE)
  if (identical(actual, "local")) {
    v2r4_assert(reports_local_substitution,
                "%s_LOCAL_SUBSTITUTION_MESSAGE_REQUIRED", label)
  } else {
    v2r4_assert(!reports_local_substitution,
                "%s_PARAMETRIC_SUBSTITUTION_MESSAGE_CONTRADICTION", label)
  }
  invisible(TRUE)
}

v2r4_normalize_exact_character_list <- function(x, expected, label) {
  v2r4_assert(is.character(expected) && length(expected) > 0L &&
                all(!is.na(expected)) && all(nzchar(expected)) &&
                !anyDuplicated(expected),
              "%s_EXPECTED_CHARACTER_CONTRACT_INVALID", label)
  v2r4_assert(is.list(x) && is.null(names(x)) &&
                length(x) == length(expected),
              "%s_JSON_LIST_SHAPE_INVALID", label)
  scalar_character <- vapply(
    x,
    function(value) {
      is.character(value) && length(value) == 1L &&
        !is.na(value) && nzchar(value)
    },
    logical(1)
  )
  v2r4_assert(all(scalar_character), "%s_JSON_LIST_ELEMENT_INVALID", label)
  normalized <- unlist(x, use.names = FALSE)
  v2r4_assert(!anyDuplicated(normalized), "%s_JSON_LIST_DUPLICATE_INVALID", label)
  v2r4_assert(identical(normalized, expected),
              "%s_VALUE_OR_ORDER_INVALID", label)
  normalized
}

v2r4_final_ack_binding_fields <- function() {
  c(
    "run_id", "lease_id", "lease_manifest_sha256", "source_authority_payload_sha256",
    "model_source_sha256",
    "validator_source_sha256", "shared_library_sha256", "analysis_config_sha256",
    "external_launcher_sha256", "os_guard_sha256", "launcher_pid", "launcher_start_utc",
    "guard_owner_role", "guard_owner_pid", "guard_owner_start_utc",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "guard_pid", "guard_start_utc",
    "guard_nonce", "guard_ready_receipt_sha256", "output_manifest_sha256",
    "complete_sha256", "validation_report_sha256", "validation_request_sha256",
    "snapshot_tree_sha256", "runtime_pre_digest", "runtime_post_digest",
    "runtime_pre_post_equal", "runtime_post_receipt_sha256",
    "model_request_sha256", "model_request_nonce", "model_identity_request_sha256",
    "validator_request_sha256", "validator_request_nonce",
    "validator_identity_request_sha256", "science_complete_sha256",
    "control_complete_sha256",
    "finalization_request_sha256", "final_challenge_request_sha256",
    "final_challenge_response_sha256", "final_challenge_payload_sha256",
    "final_challenge_nonce", "final_challenge_issued_utc",
    "final_challenge_responded_utc", "lease_issued_utc", "lease_expires_utc",
    "model_pid", "model_start_utc", "validator_pid", "validator_start_utc",
    "guard_mutex_file_id", "named_mutex_unavailable", "guard_mutex_held",
    "e2_gate_status", "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_guard_ready_fields <- function() {
  c(
    "schema_version", "status", "run_id", "owner_role", "owner_pid",
    "owner_start_utc", "guard_pid", "guard_start_utc", "guard_nonce",
    "guard_public_key_pem_b64", "guard_public_key_sha256", "stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_path",
    "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256",
    "mutex_file_id", "mutex_path", "challenge_directory", "stop_request_path",
    "stopped_receipt_path", "source_sha256"
  )
}

v2r4_guard_orderly_stop_payload_fields <- function() {
  c(
    "schema_version", "decision", "purpose", "status", "run_id", "lease_id",
    "source_authority_payload_sha256", "launch_context_sha256", "launcher_pid",
    "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_pid", "guard_start_utc", "guard_nonce",
    "guard_public_key_sha256", "guard_ready_receipt_sha256", "guard_mutex_file_id",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "runtime_pre_digest", "runtime_post_digest",
    "runtime_final_digest", "runtime_pre_post_final_equal", "runtime_final_receipt_sha256",
    "final_ack_envelope_sha256", "final_ack_payload_sha256",
    "model_consumer_receipt_sha256", "validator_consumer_receipt_sha256",
    "shutdown_nonce", "issued_utc", "expires_utc", "e2_gate_status",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_guard_abort_stop_payload_fields <- function() {
  c(
    "schema_version", "decision", "purpose", "status", "run_id",
    "source_authority_payload_sha256", "launch_context_sha256", "launcher_pid",
    "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_pid", "guard_start_utc", "guard_nonce",
    "guard_public_key_sha256", "guard_ready_receipt_sha256", "guard_mutex_file_id",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "shutdown_nonce", "failure_phase",
    "failure_token_sha256", "issued_utc", "expires_utc", "e2_gate_status",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_guard_orderly_stopped_payload_fields <- function() {
  c(
    "schema_version", "status", "request_sha256", "request_payload_sha256",
    "request_hmac_sha256", "decision", "purpose", "run_id", "lease_id",
    "launcher_pid", "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_pid", "guard_start_utc", "guard_nonce",
    "guard_public_key_sha256", "guard_ready_receipt_sha256", "guard_mutex_file_id",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "runtime_pre_digest", "runtime_post_digest",
    "runtime_final_digest", "runtime_pre_post_final_equal", "runtime_final_receipt_sha256",
    "final_ack_envelope_sha256", "final_ack_payload_sha256",
    "model_consumer_receipt_sha256", "validator_consumer_receipt_sha256",
    "shutdown_nonce", "mutex_released", "released_utc", "e2_gate_status",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_guard_abort_stopped_payload_fields <- function() {
  c(
    "schema_version", "status", "request_sha256", "request_payload_sha256",
    "request_hmac_sha256", "decision", "purpose", "run_id", "launcher_pid",
    "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_pid", "guard_start_utc", "guard_nonce",
    "guard_public_key_sha256", "guard_ready_receipt_sha256", "guard_mutex_file_id",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "shutdown_nonce", "failure_phase",
    "failure_token_sha256", "mutex_released", "released_utc", "e2_gate_status",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_load_config_buffer <- function(config_buffer) {
  cfg <- v2r4_read_json_raw(config_buffer, "ANALYSIS_CONFIG")
  v2r4_assert(identical(cfg$schema_version, 7L), "CONFIG_SCHEMA_VERSION_MISMATCH")
  v2r4_assert(identical(cfg$dataset_id, "GSE205906"), "CONFIG_DATASET_MISMATCH")
  v2r4_assert(identical(cfg$model_version, "GSE205906_BULK_MODEL_V2_R4_2026-08-22"),
              "CONFIG_MODEL_VERSION_MISMATCH")
  v2r4_assert(identical(cfg$status,
                        "PRODUCTION_CAPABLE_EXTERNAL_SINGLE_USE_AUTHORIZATION_REQUIRED"),
              "CONFIG_STATUS_EXCEEDS_ALLOWED_CEILING")
  v2r4_assert(isTRUE(cfg$governance$production_execution_authorized) &&
                isTRUE(cfg$governance$production_validation_authorized) &&
                identical(cfg$governance$production_static_authorization, FALSE) &&
                isTRUE(cfg$governance$authorization_is_per_run_only),
              "CONFIG_PRODUCTION_SINGLE_USE_GOVERNANCE_INVALID")
  red_team <- cfg$v2_r2_independent_red_team_trust_root
  v2r4_assert(is.list(red_team) && length(red_team) == 3L &&
                identical(vapply(red_team, `[[`, character(1), "relative_path"), c(
                  "results/diagnostics/GSE205906_bulk_model_v2_r2_independent_red_team_2026-08-21.md",
                  "results/diagnostics/GSE205906_bulk_model_v2_r2_independent_red_team_2026-08-21.tsv",
                  "results/diagnostics/GSE205906_bulk_model_v2_r2_independent_red_team_2026-08-21_SHA256.tsv"
                )) && identical(vapply(red_team, `[[`, integer(1), "bytes"),
                                        c(19772L, 12718L, 334L)) &&
                identical(vapply(red_team, `[[`, character(1), "sha256"), c(
                  "35717d353746a0f56f0515a331eb4a7a9b51438a6c3e51b7e601285049905636",
                  "daedcb6d5bb6752eb9bd727b941b69e012e7d3214b8615665bb6f5e54887c357",
                  "90e7b15eb10d6ea2bd7b86d444d25e96ee5736b14e1ff9e1814bba5d3f1252fe"
                )), "CONFIG_V2_R2_INDEPENDENT_RED_TEAM_TRUST_ROOT_INVALID")
  ext <- cfg$external_runtime_authority
  expected_runtime_rehash_phases <- c(
    "PRE_CHILD_UNDER_GLOBAL_LEASE",
    "POST_CHILD_UNDER_GLOBAL_LEASE"
  )
  runtime_rehash_phases <- v2r4_normalize_exact_character_list(
    ext$runtime_rehash_phases,
    expected_runtime_rehash_phases,
    "CONFIG_EXTERNAL_RUNTIME_REHASH_PHASES"
  )
  request_transport <- ext$child_request_transport
  request_transport_fields <- c(
    "schema_version", "transport", "argc_after_script", "request_maximum_bytes",
    "encoding", "expected_sha_argument", "capture_policy",
    "launcher_file_hold_policy", "child_environment_exact_allowlist",
    "model", "validator"
  )
  request_role_fields <- c(
    "schema_version", "role", "purpose", "exact_fields", "field_type_policy"
  )
  v2r4_assert(is.list(request_transport) &&
                identical(names(request_transport), request_transport_fields) &&
                is.list(request_transport$model) &&
                identical(names(request_transport$model), request_role_fields) &&
                is.list(request_transport$validator) &&
                identical(names(request_transport$validator), request_role_fields),
              "CONFIG_CHILD_REQUEST_SCHEMA_INVALID")
  model_request_contract <- v2r4_child_request_contract("MODEL")
  validator_request_contract <- v2r4_child_request_contract("VALIDATOR")
  model_request_fields <- v2r4_normalize_exact_character_list(
    request_transport$model$exact_fields,
    model_request_contract$exact_fields,
    "CONFIG_MODEL_REQUEST_FIELDS"
  )
  validator_request_fields <- v2r4_normalize_exact_character_list(
    request_transport$validator$exact_fields,
    validator_request_contract$exact_fields,
    "CONFIG_VALIDATOR_REQUEST_FIELDS"
  )
  expected_child_environment <- c(
    "SystemRoot", "WINDIR", "ComSpec", "PATH", "PATHEXT", "TEMP", "TMP",
    "USERPROFILE", "HOMEDRIVE", "HOMEPATH", "APPDATA", "LOCALAPPDATA",
    "PROGRAMDATA", "ProgramFiles", "ProgramFiles(x86)", "CommonProgramFiles",
    "CommonProgramFiles(x86)", "NUMBER_OF_PROCESSORS", "PROCESSOR_ARCHITECTURE", "OS",
    "GSE205906_V2_R4_EXTERNAL_LAUNCH", "GSE205906_V2_R4_LAUNCHER_TOKEN",
    "GSE205906_V2_R4_LAUNCH_CONTEXT", "GSE205906_V2_R4_LAUNCH_CONTEXT_HMAC",
    "GSE205906_V2_R4_PWSH_PATH", "GSE205906_V2_R4_PWSH_SHA256"
  )
  child_environment <- v2r4_normalize_exact_character_list(
    request_transport$child_environment_exact_allowlist,
    expected_child_environment, "CONFIG_CHILD_ENVIRONMENT_ALLOWLIST"
  )
  v2r4_assert(identical(request_transport$schema_version, 1L) &&
                identical(request_transport$transport,
                          model_request_contract$transport) &&
                identical(request_transport$argc_after_script, 2L) &&
                identical(request_transport$request_maximum_bytes,
                          model_request_contract$maximum_bytes) &&
                identical(request_transport$encoding, "UTF8_NO_BOM_FINAL_LF") &&
                identical(request_transport$expected_sha_argument,
                          "BARE_LOWERCASE_SHA256") &&
                identical(request_transport$capture_policy,
                          paste0("SINGLE_CAPTURE_SIZE_THEN_SHA256_EXACT_SCHEMA_PARSE_",
                                 "SAME_RAW_NO_REOPEN")) &&
                identical(request_transport$launcher_file_hold_policy,
                           paste0("HELD_DENY_WRITE_DELETE_FILE_ID_BYTES_SHA256_",
                                  "PRE_POST_CHILD")) &&
                identical(child_environment, expected_child_environment) &&
                identical(request_transport$model$schema_version,
                          model_request_contract$schema_version) &&
                identical(request_transport$model$role, model_request_contract$role) &&
                identical(request_transport$model$purpose,
                          model_request_contract$purpose) &&
                identical(request_transport$model$field_type_policy,
                          "ALL_FIELDS_EXACT_NONEMPTY_STRING") &&
                identical(model_request_fields, model_request_contract$exact_fields) &&
                identical(request_transport$validator$schema_version,
                          validator_request_contract$schema_version) &&
                identical(request_transport$validator$role,
                          validator_request_contract$role) &&
                identical(request_transport$validator$purpose,
                          validator_request_contract$purpose) &&
                identical(request_transport$validator$field_type_policy,
                          "ALL_FIELDS_EXACT_NONEMPTY_STRING") &&
                identical(validator_request_fields,
                          validator_request_contract$exact_fields),
              "CONFIG_CHILD_REQUEST_TRANSPORT_CONTRACT_INVALID")
  v2r4_assert(identical(ext$decision,
                        "AUTHORIZE_GSE205906_BULK_MODEL_V2_R4_SOURCE_RUNTIME_SNAPSHOT") &&
                identical(ext$purpose, "SOURCE_RUNTIME_PRE_MODEL_EXECUTION") &&
                identical(ext$production_public_key_sha256,
                          "9b9193a9fa0ce53507b07e1ad8a207fb7a4751ace52258aaef96eba1481abbb3") &&
                identical(ext$path_execution_identity_policy,
                          "HELD_DENY_WRITE_DELETE_FILE_ID_BYTES_SHA256_PRE_POST") &&
                identical(runtime_rehash_phases, expected_runtime_rehash_phases) &&
                identical(ext$production_external_assets_installed, TRUE) &&
                isTRUE(ext$launcher_token_required_before_candidate_library_load),
              "CONFIG_EXTERNAL_RUNTIME_AUTHORITY_CONTRACT_INVALID")
  key_pins <- unlist(ext$key_separation[c(
    "source_snapshot_public_key_sha256", "run_authorization_public_key_sha256",
    "live_lease_public_key_sha256", "output_marker_public_key_sha256",
    "final_ack_public_key_sha256"
  )], use.names = FALSE)
  v2r4_assert(length(key_pins) == 5L && all(grepl("^[0-9a-f]{64}$", key_pins)) &&
                length(unique(key_pins)) == 5L &&
                isTRUE(ext$key_separation$all_five_pairwise_distinct_required) &&
                isTRUE(ext$key_separation$private_keys_outside_project_required) &&
                identical(key_pins[[1L]], cfg$authority$production_source_public_key_sha256) &&
                identical(key_pins[[2L]], ext$production_public_key_sha256) &&
                identical(key_pins[[3L]], cfg$lease$production_public_key_sha256) &&
                identical(key_pins[[4L]], cfg$output_authority$production_public_key_sha256) &&
                identical(key_pins[[5L]], cfg$finalization$production_public_key_sha256),
              "CONFIG_FIVE_KEY_SEPARATION_INVALID")
  lifecycle <- ext$guard_lifecycle
  orderly_stop_fields <- v2r4_normalize_exact_character_list(
    lifecycle$orderly_stop_payload_exact_fields,
    v2r4_guard_orderly_stop_payload_fields(), "CONFIG_GUARD_ORDERLY_STOP_FIELDS"
  )
  abort_stop_fields <- v2r4_normalize_exact_character_list(
    lifecycle$abort_stop_payload_exact_fields,
    v2r4_guard_abort_stop_payload_fields(), "CONFIG_GUARD_ABORT_STOP_FIELDS"
  )
  orderly_stopped_fields <- v2r4_normalize_exact_character_list(
    lifecycle$orderly_stopped_payload_exact_fields,
    v2r4_guard_orderly_stopped_payload_fields(), "CONFIG_GUARD_ORDERLY_STOPPED_FIELDS"
  )
  abort_stopped_fields <- v2r4_normalize_exact_character_list(
    lifecycle$abort_stopped_payload_exact_fields,
    v2r4_guard_abort_stopped_payload_fields(), "CONFIG_GUARD_ABORT_STOPPED_FIELDS"
  )
  v2r4_assert(identical(lifecycle$owner_role, "EXTERNAL_LAUNCHER") &&
                identical(lifecycle$launch_policy,
                          "LAUNCHER_DIRECT_PROCESSSTARTINFO_ARGUMENTLIST_PROCESS_HANDLE_HELD") &&
                identical(lifecycle$owner_identity_policy,
                          paste0("ACTUAL_LAUNCHER_PID_START_CONTEXT_HMAC_SIGNED_LEASE_",
                                 "AND_INDEPENDENT_OS_PROBE")) &&
                identical(lifecycle$stop_hmac_key_bytes, 32L) &&
                identical(lifecycle$stop_hmac_key_environment_variable,
                          "GSE205906_V2_R4_GUARD_STOP_HMAC_KEY") &&
                identical(lifecycle$stop_hmac_key_distribution,
                          paste0("LAUNCHER_MEMORY_AND_GUARD_SPECIFIC_PROCESSSTARTINFO_",
                                 "ENVIRONMENT_ONLY")) &&
                identical(lifecycle$stop_hmac_key_forbidden_locations,
                           "DISK_ARGV_MODEL_VALIDATOR_PROJECT") &&
                identical(lifecycle$finalization_public_key_capture_policy,
                          "EXTERNAL_SOURCE_CAPTURED_ONCE_THEN_CREATE_NEW_EXACT_COPY_UNDER_ATTEMPT") &&
                identical(lifecycle$finalization_public_key_guard_consumption,
                          "ATTEMPT_CAPTURED_COPY_ONLY_NO_EXTERNAL_SOURCE_REOPEN") &&
                identical(lifecycle$finalization_public_key_copy_hold_policy,
                          "HELD_DENY_WRITE_DELETE_FILE_ID_BYTES_SHA256_UNTIL_GUARD_EXIT") &&
                identical(lifecycle$stop_wait_seconds, 900L) &&
                isTRUE(lifecycle$stop_nonce_single_use_create_new_required) &&
                isTRUE(lifecycle$stop_request_captured_bytes_hmac_required) &&
                isTRUE(lifecycle$stop_request_launcher_hold_deny_write_delete_until_guard_exit) &&
                isTRUE(lifecycle$stopped_receipt_guard_ephemeral_signature_required) &&
                isTRUE(lifecycle$guard_orderly_stop_final_ack_signature_and_consumer_binding_revalidation_required) &&
                isTRUE(lifecycle$guard_unlock_before_stopped_receipt_publish) &&
                isTRUE(lifecycle$launcher_reacquires_same_file_id_exclusive_until_success) &&
                identical(orderly_stop_fields, v2r4_guard_orderly_stop_payload_fields()) &&
                identical(abort_stop_fields, v2r4_guard_abort_stop_payload_fields()) &&
                identical(orderly_stopped_fields, v2r4_guard_orderly_stopped_payload_fields()) &&
                identical(abort_stopped_fields, v2r4_guard_abort_stopped_payload_fields()),
              "CONFIG_LAUNCHER_OWNED_GUARD_LIFECYCLE_INVALID")
  lease_identity_fields <- v2r4_normalize_exact_character_list(
    cfg$lease$identity_request_exact_fields,
    v2r4_lease_identity_request_fields(), "CONFIG_LEASE_IDENTITY_REQUEST_FIELDS"
  )
  combined_lease_fields <- v2r4_normalize_exact_character_list(
    cfg$lease$combined_request_exact_fields,
    v2r4_combined_lease_request_fields(), "CONFIG_COMBINED_LEASE_REQUEST_FIELDS"
  )
  signed_lease_fields <- v2r4_normalize_exact_character_list(
    cfg$lease$signed_lease_exact_fields,
    v2r4_signed_lease_fields(), "CONFIG_SIGNED_LEASE_FIELDS"
  )
  v2r4_assert(identical(cfg$lease$schema_version, 5L) &&
                identical(cfg$lease$purpose, "LIVE_LAUNCHER_OWNED_GUARD_SNAPSHOT_LEASE") &&
                identical(cfg$lease$identity_request_schema_version,
                          "GSE205906_bulk_model_v2_r4_lease_identity_request_v1") &&
                identical(cfg$lease$identity_request_status,
                          "LIVE_CHILD_IDENTITY_PUBLISHED_AWAITING_COMBINED_SIGNED_LEASE") &&
                identical(cfg$lease$combined_request_schema_version,
                          "GSE205906_bulk_model_v2_r4_combined_lease_request_v1") &&
                identical(cfg$lease$combined_request_status,
                          paste0("MODEL_VALIDATOR_GUARD_IDENTITIES_VERIFIED_",
                                 "AWAITING_EXTERNAL_SIGNED_LEASE")) &&
                identical(lease_identity_fields, v2r4_lease_identity_request_fields()) &&
                identical(combined_lease_fields, v2r4_combined_lease_request_fields()) &&
                identical(signed_lease_fields, v2r4_signed_lease_fields()) &&
                identical(cfg$lease$guard_owner_role, "EXTERNAL_LAUNCHER") &&
                isTRUE(cfg$lease$launcher_must_hold_guard_mutex_until_final_runtime_rehash_and_authenticated_stop) &&
                isTRUE(cfg$lease$model_and_validator_must_remain_live_until_signed_final_ack_consumption) &&
                isTRUE(cfg$lease$validator_must_recheck_live_launcher_model_validator_guard_pid_start_and_mutex) &&
                isTRUE(cfg$lease$guard_must_start_before_model_and_validator) &&
                isTRUE(cfg$lease$model_and_validator_identity_requests_required_before_combined_lease_request) &&
                isTRUE(cfg$lease$model_and_validator_must_validate_same_captured_signed_lease_before_data_read) &&
                identical(cfg$lease$replacement_guard_permitted, FALSE),
              "CONFIG_V2_R4_LIVE_LEASE_CONTRACT_INVALID")
  finalization_bindings <- v2r4_normalize_exact_character_list(
    cfg$finalization$required_bindings,
    v2r4_final_ack_binding_fields(),
    "CONFIG_FINAL_ACK_REQUIRED_BINDINGS"
  )
  v2r4_assert(identical(cfg$finalization$schema_version, 3L) &&
                identical(cfg$finalization$decision,
                         "ATTEST_GSE205906_BULK_MODEL_V2_R4_FINAL_VALIDATED_RUN") &&
                identical(cfg$finalization$purpose,
                          "POST_VALIDATION_FINAL_ACKNOWLEDGEMENT") &&
                identical(cfg$finalization$production_public_key_sha256,
                          "52532de691dcb71b9f912eaaf2ee80ed80ab754ef320bc989c760265cc71caa5") &&
                identical(cfg$finalization$unsigned_ack_permitted, FALSE) &&
                identical(cfg$finalization$production_external_finalizer_installed, TRUE) &&
                isTRUE(cfg$finalization$runtime_pre_post_equality_required) &&
                identical(cfg$finalization$runtime_post_receipt_schema_version,
                          "GSE205906_bulk_model_v2_r4_runtime_post_receipt_v1") &&
                identical(cfg$finalization$finalization_request_schema_version,
                          "GSE205906_bulk_model_v2_r4_finalization_request_v1") &&
                identical(cfg$finalization$science_complete_schema_version,
                          "GSE205906_bulk_model_v2_r4_science_complete_v1") &&
                identical(cfg$finalization$control_complete_schema_version,
                          "GSE205906_bulk_model_v2_r4_control_complete_v1") &&
                identical(cfg$finalization$final_challenge_request_schema_version,
                          "GSE205906_bulk_model_v2_r4_final_challenge_v1") &&
                identical(cfg$finalization$final_challenge_response_payload_schema_version,
                          "GSE205906_bulk_model_v2_r4_final_challenge_response_payload_v1") &&
                identical(cfg$finalization$final_challenge_response_envelope_schema_version,
                          "GSE205906_bulk_model_v2_r4_final_challenge_response_envelope_v1") &&
                identical(cfg$finalization$consumer_receipt_schema_version,
                          "GSE205906_bulk_model_v2_r4_final_ack_consumer_receipt_v1") &&
                identical(cfg$finalization$final_ack_payload_schema_version,
                          "GSE205906_bulk_model_v2_r4_final_ack_payload_v1") &&
                identical(cfg$finalization$final_ack_envelope_schema_version,
                          "GSE205906_bulk_model_v2_r4_final_ack_envelope_v1") &&
                identical(cfg$finalization$guard_stop_request_envelope_schema_version,
                          "GSE205906_bulk_model_v2_r4_guard_stop_request_envelope_v1") &&
                identical(cfg$finalization$guard_stopped_envelope_schema_version,
                          "GSE205906_bulk_model_v2_r4_guard_stopped_envelope_v1") &&
                identical(cfg$finalization$wait_seconds_source,
                          "SAME_CAPTURED_ANALYSIS_CONFIG_EXECUTION_CONTRACT_EXACT_900") &&
                identical(cfg$finalization$wait_clock,
                          paste0("LOCAL_SYSTEM_DIAGNOSTICS_STOPWATCH_WITH_UTC_",
                                 "ISSUED_DEADLINE_BINDING")) &&
                isTRUE(cfg$finalization$launcher_holds_global_mutex_and_execution_handles_through_guard_stopped_receipt) &&
                isTRUE(cfg$finalization$consumer_receipt_copies_every_final_ack_binding) &&
                isTRUE(cfg$finalization$final_challenge_nonce_single_use_required) &&
                identical(cfg$finalization$final_challenge_max_age_seconds, 10L) &&
                isTRUE(cfg$finalization$finalizer_requires_live_launcher_model_validator_guard_pid_start) &&
                isTRUE(cfg$finalization$finalizer_requires_named_mutex_unavailable) &&
                isTRUE(cfg$finalization$finalizer_requires_same_guard_mutex_file_id_and_holder_signature) &&
                identical(finalization_bindings, v2r4_final_ack_binding_fields()),
              "CONFIG_V2_R4_FINALIZATION_CONTRACT_INVALID")
  v2r4_assert(
    identical(cfg$execution_contract$signed_final_ack_wait_seconds, 900L) &&
      identical(unlist(cfg$execution_contract$execution_modes, use.names = FALSE),
                c("SYNTHETIC", "PRODUCTION")) &&
      isTRUE(cfg$execution_contract$production_nonce_create_new_required_before_guard_or_child) &&
      isTRUE(cfg$execution_contract$output_and_final_ack_signers_must_be_distinct),
    "CONFIG_EXECUTION_CONTRACT_INVALID"
  )
  numeric_contract <- cfg$analysis$independent_numeric_revalidation
  v2r4_assert(isTRUE(numeric_contract$required) &&
                identical(numeric_contract$scope,
                          "COMPLETE_RETAINED_GENE_BY_EXACT_CONTRAST_FAMILY") &&
                identical(numeric_contract$absolute_tolerance, 1e-10) &&
                identical(numeric_contract$relative_tolerance, 1e-8) &&
                isTRUE(numeric_contract$na_pattern_must_be_exact) &&
                isTRUE(numeric_contract$tolerance_may_not_be_tuned_after_results),
              "CONFIG_INDEPENDENT_NUMERIC_REVALIDATION_CONTRACT_INVALID")
  v2r4_assert_deseq_contract(cfg)
  v2r4_assert(
    identical(cfg$input_snapshot$required_sample_include_primary_marker, "CONDITIONAL") &&
      identical(cfg$input_snapshot$required_sample_metadata_status,
                "derived_from_official_GEO_series_matrix") &&
      identical(cfg$input_snapshot$conditional_inclusion_release_gate,
                "ALL_9_SAMPLE_AND_COHORT_SPECIES_RULES_PASS_BEFORE_COUNTS_CAPTURE") &&
      identical(cfg$input_snapshot$conditional_inclusion_release_status,
                paste0("RELEASED_AFTER_",
                       "ALL_9_SAMPLE_AND_COHORT_SPECIES_RULES_PASS_BEFORE_COUNTS_CAPTURE")) &&
      identical(cfg$input_snapshot$plm_metadata_conflict_policy,
                paste0("GROUP_CONTRAST_ONLY_NO_ABSOLUTE_DOSE_DOSE_RESPONSE_",
                       "COMPONENT_RATIO_OR_SYNERGY")),
    "CONFIG_CONDITIONAL_INCLUSION_CONTRACT_INVALID"
  )
  v2r4_assert(
    identical(cfg$input_snapshot$species_diagnostic_method,
              "same_deterministic_prefix_dual_reference_unique_Subjunc") &&
      identical(as.integer(cfg$input_snapshot$species_fragments_per_library), 500000L) &&
      identical(as.numeric(cfg$input_snapshot$species_minimum_rat_unique_mapping_rate), 0.60) &&
      identical(as.numeric(cfg$input_snapshot$species_maximum_human_unique_mapping_rate), 0.20) &&
      identical(as.numeric(cfg$input_snapshot$species_minimum_rat_minus_human_rate), 0.10) &&
      identical(as.numeric(cfg$input_snapshot$species_minimum_rat_to_human_rate_ratio), 1.10) &&
      identical(as.numeric(cfg$input_snapshot$species_minimum_cohort_median_rate_difference), 0.20),
    "CONFIG_SPECIES_THRESHOLD_CONTRACT_INVALID"
  )
  producer_roles <- c("align_count_source", "assembly_source",
                      "species_sample_source", "species_cohort_source")
  producer_paths <- c(
    align_count_source = "code/R/align_count_one_sample.R",
    assembly_source = "code/R/assemble_rnaseq_counts.R",
    species_sample_source = "code/R/species_screen_one_sample.R",
    species_cohort_source = "code/R/assemble_species_screen.R"
  )
  producer_bytes <- c(
    align_count_source = 9837, assembly_source = 14780,
    species_sample_source = 8553, species_cohort_source = 5589
  )
  producer_sha256 <- c(
    align_count_source = "6875cbd2a1fc2be03ce752324644228824c0591a2108e74e0b771457c0020701",
    assembly_source = "38be1eb124703f7c79168b86ee1d7044378311d15e4723fa10fb217a142d523f",
    species_sample_source = "5ffefdac0b2af26bed105e6b21b436b621ca8a146cf6b5c3a9b89b77f7422dab",
    species_cohort_source = "d219c7e20dfb49fc80ace0267be776893072ab1bafd0161ad1611732aedf5f95"
  )
  producer_contract <- cfg$input_snapshot$producer_source_identities
  v2r4_assert(identical(names(producer_contract), producer_roles),
              "CONFIG_PRODUCER_SOURCE_ROLE_ORDER_INVALID")
  for (role in producer_roles) {
    entry <- producer_contract[[role]]
    v2r4_assert(identical(names(entry), c("relative_path", "bytes", "sha256")) &&
                  identical(entry$relative_path, unname(producer_paths[[role]])) &&
                  identical(as.numeric(entry$bytes), unname(producer_bytes[[role]])) &&
                  identical(entry$sha256, unname(producer_sha256[[role]])),
                paste0("CONFIG_PRODUCER_SOURCE_IDENTITY_INVALID:", role))
  }
  v2r4_assert(identical(cfg$input_snapshot$featurecounts_gene_id_output_column, "GeneID") &&
                identical(cfg$input_snapshot$featurecounts_symbol_attribute, "gene") &&
                identical(cfg$input_snapshot$featurecounts_biotype_attribute, "gene_biotype") &&
                identical(unlist(cfg$input_snapshot$featurecounts_required_annotation_columns,
                                 use.names = FALSE), c("GeneID", "gene", "gene_biotype")) &&
                !isTRUE(cfg$input_snapshot$gene_name_fallback_permitted),
              "CONFIG_FROZEN_FEATURECOUNTS_INTERFACE_INVALID")
  parser_contract <- cfg$input_snapshot$rds_parser_runtime_contract
  v2r4_assert(identical(parser_contract$R_exact, "4.6.1") &&
                identical(parser_contract$digest_exact, "0.6.39") &&
                identical(parser_contract$R_exact, cfg$software$R) &&
                identical(parser_contract$digest_exact, cfg$software$digest) &&
                isTRUE(parser_contract$captured_raw_only) &&
                isTRUE(parser_contract$unknown_compression_forbidden) &&
                isTRUE(parser_contract$uncompressed_rds_forbidden) &&
                isTRUE(parser_contract$serialization_full_consumption_required),
              "CONFIG_RDS_PARSER_RUNTIME_CONTRACT_INVALID")
  rds_contracts <- cfg$input_snapshot$rds_role_contracts
  v2r4_assert(setequal(names(rds_contracts), c("assembled_counts_rds", "gene_set_object")) &&
                identical(rds_contracts$assembled_counts_rds$compression, "xz") &&
                identical(rds_contracts$assembled_counts_rds$magic_hex, "fd377a585a00") &&
                identical(rds_contracts$assembled_counts_rds$container_contract,
                          "XZ_SINGLE_STREAM_NO_STREAM_PADDING_NO_CONCAT") &&
                identical(rds_contracts$assembled_counts_rds$stream_flags_hex, "0001") &&
                identical(rds_contracts$assembled_counts_rds$check_type, "CRC32") &&
                isTRUE(rds_contracts$assembled_counts_rds$header_footer_index_crc_required) &&
                isTRUE(rds_contracts$assembled_counts_rds$canonical_vli_index_and_block_span_required) &&
                isTRUE(rds_contracts$assembled_counts_rds$index_uncompressed_size_total_required) &&
                identical(rds_contracts$assembled_counts_rds$exact_memDecompress_type, "xz") &&
                identical(rds_contracts$assembled_counts_rds$writer_source_relative_path,
                          "code/R/assemble_rnaseq_counts.R") &&
                identical(rds_contracts$assembled_counts_rds$writer_source_sha256,
                          "38be1eb124703f7c79168b86ee1d7044378311d15e4723fa10fb217a142d523f") &&
                identical(rds_contracts$assembled_counts_rds$writer_call_contract,
                          "saveRDS_compress_xz") &&
                identical(rds_contracts$gene_set_object$compression, "gzip") &&
                identical(rds_contracts$gene_set_object$magic_hex, "1f8b08") &&
                identical(rds_contracts$gene_set_object$container_contract,
                          "RFC1952_FLG0_LAYERED_WITH_SIGNED_CAPTURED_IDENTITY") &&
                as.integer(rds_contracts$gene_set_object$flg_exact) == 0L &&
                identical(rds_contracts$gene_set_object$authority_identity_contract,
                          "EXACT_CAPTURED_BYTES_SIZE_SHA256_BEFORE_PARSE") &&
                isTRUE(rds_contracts$gene_set_object$gzcon_full_drain_warnings_errors_fatal) &&
                isTRUE(rds_contracts$gene_set_object$trailer_crc32_and_isize_required) &&
                identical(rds_contracts$gene_set_object$exact_memDecompress_type, "gzip") &&
                isTRUE(rds_contracts$gene_set_object$decoder_outputs_must_be_identical) &&
                identical(rds_contracts$gene_set_object$member_boundary_observable, FALSE) &&
                identical(rds_contracts$gene_set_object$member_boundary_claim,
                          "NOT_OBSERVABLE_FROM_GZCON_RAWCONNECTION") &&
                identical(rds_contracts$gene_set_object$concat_append_exclusion_claim,
                          paste0("AUTHORITY_IDENTITY_PRIMARY_CONTAINER_CHECKS_ADDITIONAL_",
                                 "NOT_STANDALONE_BOUNDARY_PROOF")) &&
                identical(rds_contracts$gene_set_object$writer_source_relative_path,
                          "code/R/freeze_gene_sets.R") &&
                identical(rds_contracts$gene_set_object$writer_source_sha256,
                          "01b8f54303656fe063476ce272b30dd93aa668d023ecb6a6c3c947a7d908d6ea") &&
                identical(rds_contracts$gene_set_object$writer_call_contract,
                          "saveRDS_version_3_default_compress_TRUE_gzip"),
              "CONFIG_RDS_ROLE_COMPRESSION_CONTRACT_INVALID")
  cfg
}

v2r4_verify_detached_bundle <- function(manifest_path, signature_path, public_key_path,
                                       expected_public_key_sha256, label) {
  v2r4_require_packages(c("digest", "openssl"))
  manifest <- v2r4_capture_file(manifest_path, paste0(label, "_MANIFEST"))
  signature <- v2r4_capture_file(signature_path, paste0(label, "_SIGNATURE"))
  public_key <- v2r4_capture_file(public_key_path, paste0(label, "_PUBLIC_KEY"))
  v2r4_assert(grepl("^[0-9a-f]{64}$", expected_public_key_sha256) &&
                identical(public_key$sha256, tolower(expected_public_key_sha256)),
              "%s_PUBLIC_KEY_PIN_MISMATCH", label)
  key <- tryCatch(openssl::read_pubkey(v2r4_parse_text_raw(public_key, paste0(label, "_PUBLIC_KEY"))),
                  error = function(e) v2r4_stop("%s_PUBLIC_KEY_INVALID:%s", label, conditionMessage(e)))
  verified <- tryCatch(
    openssl::signature_verify(manifest$raw, signature$raw, hash = openssl::sha256, pubkey = key),
    error = function(e) FALSE
  )
  v2r4_assert(isTRUE(verified), "%s_SIGNATURE_INVALID", label)
  list(manifest = manifest, signature = signature, public_key = public_key)
}

v2r4_parse_authority_buffer <- function(manifest_buffer, cfg, expected_decision,
                                       allowed_materialization = c("CONTENT_OBJECT", "IDENTITY_ONLY")) {
  auth <- v2r4_read_tsv_raw(manifest_buffer, "SOURCE_AUTHORITY")
  required <- unlist(cfg$authority$required_authority_columns, use.names = FALSE)
  v2r4_assert(identical(names(auth), required), "AUTHORITY_SCHEMA_NOT_EXACT")
  v2r4_assert(nrow(auth) > 0L && !anyDuplicated(auth$role) &&
                !anyDuplicated(auth$relative_path), "AUTHORITY_ROWS_NOT_UNIQUE")
  v2r4_assert(all(auth$schema_version == cfg$authority$schema_version) &&
                all(auth$dataset_id == cfg$dataset_id) &&
                all(auth$model_version == cfg$model_version) &&
                all(auth$decision == expected_decision), "AUTHORITY_HEADER_SEMANTICS_INVALID")
  v2r4_assert(all(auth$materialization %in% allowed_materialization),
              "AUTHORITY_MATERIALIZATION_INVALID")
  v2r4_assert(all(grepl("^[0-9a-f]{64}$", auth$sha256)) &&
                all(is.finite(auth$bytes) & auth$bytes >= 0 & auth$bytes == floor(auth$bytes)),
              "AUTHORITY_IDENTITY_DOMAIN_INVALID")
  content_idx <- base::which(auth$materialization == "CONTENT_OBJECT")
  if (length(content_idx)) {
    expected_rel <- sprintf("objects/sha256/%s/%s", auth$sha256[content_idx],
                            basename(auth$relative_path[content_idx]))
    v2r4_assert(identical(auth$relative_path[content_idx], expected_rel),
                "AUTHORITY_CONTENT_ADDRESS_PATH_INVALID")
  }
  auth
}

v2r4_role_row <- function(authority, role) {
  role_value <- as.character(role)
  idx <- base::which(as.character(authority[["role"]]) == role_value)
  v2r4_assert(length(idx) == 1L, "AUTHORITY_ROLE_NOT_UNIQUE:%s", role_value)
  authority[idx]
}

v2r4_bind_runtime_roles <- function(authority, cfg, project_root, runtime_buffers) {
  expected_map <- cfg$authority$dynamic_runtime_paths
  required_roles <- names(expected_map)
  v2r4_assert(setequal(required_roles, unlist(cfg$authority$required_dynamic_roles, use.names = FALSE)),
              "CONFIG_DYNAMIC_ROLE_MAP_INVALID")
  v2r4_assert(setequal(names(runtime_buffers), required_roles), "RUNTIME_BUFFER_ROLE_SET_INVALID")
  rows <- lapply(required_roles, function(role_name) {
    row <- v2r4_role_row(authority, role_name)
    expected_rel <- expected_map[[role_name]]
    expected_abs <- v2r4_norm_existing(file.path(project_root, expected_rel))
    buffer <- runtime_buffers[[role_name]]
    v2r4_assert(inherits(buffer, "v2r4_captured_bytes") &&
                  identical(buffer$canonical_path, expected_abs),
                "RUNTIME_CANONICAL_PATH_MISMATCH:%s", role_name)
    v2r4_assert(identical(row$runtime_relative_path[[1L]], expected_rel) &&
                  row$materialization[[1L]] == "CONTENT_OBJECT" &&
                  as.numeric(row$bytes[[1L]]) == buffer$bytes &&
                  identical(row$sha256[[1L]], buffer$sha256),
                "RUNTIME_ROLE_BINDING_MISMATCH:%s", role_name)
    data.table::data.table(role = role_name, runtime_relative_path = expected_rel,
                           bytes = buffer$bytes, sha256 = buffer$sha256,
                           status = "PASS_RUNTIME_EXACT")
  })
  data.table::rbindlist(rows)
}

v2r4_snapshot_tree_identity <- function(content) {
  canonical_rows <- data.table::copy(content)[order(relative_path), .(
    relative_path, bytes = as.character(as.numeric(bytes)), sha256
  )]
  tree_text <- paste0(paste(canonical_rows$relative_path, canonical_rows$bytes,
                            canonical_rows$sha256, sep = "\t", collapse = "\n"), "\n")
  list(
    bytes = charToRaw(enc2utf8(tree_text)),
    sha256 = tolower(digest::digest(charToRaw(enc2utf8(tree_text)), algo = "sha256",
                                     serialize = FALSE))
  )
}

v2r4_capture_snapshot <- function(authority, snapshot_root, defer_roles = character()) {
  root <- v2r4_norm_existing(snapshot_root)
  v2r4_reject_symlink(root, "SNAPSHOT_ROOT")
  objects_root <- v2r4_norm_existing(file.path(root, "objects", "sha256"))
  content <- authority[base::which(authority[["materialization"]] == "CONTENT_OBJECT")]
  identity_only <- authority[base::which(authority[["materialization"]] == "IDENTITY_ONLY")]
  v2r4_assert(nrow(content) > 0L, "SNAPSHOT_CONTENT_SET_EMPTY")
  v2r4_assert(all(defer_roles %in% content$role), "SNAPSHOT_DEFER_ROLE_INVALID")
  listed <- list.files(objects_root, recursive = TRUE, full.names = TRUE,
                       all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  listed <- listed[file.info(listed)$isdir %in% FALSE]
  listed_rel <- substring(gsub("\\\\", "/", listed), nchar(gsub("\\\\", "/", root)) + 2L)
  v2r4_assert(setequal(listed_rel, content$relative_path), "SNAPSHOT_PHYSICAL_OBJECT_SET_NOT_EXACT")
  capture_rows <- content[base::which(!(content[["role"]] %in% defer_roles))]
  deferred_rows <- content[base::which(content[["role"]] %in% defer_roles)]
  captured <- lapply(seq_len(nrow(capture_rows)), function(i) {
    row <- capture_rows[i]
    path <- file.path(root, row$relative_path[[1L]])
    buffer <- v2r4_capture_file(path, paste0("SNAPSHOT_OBJECT_", row$role[[1L]]), root)
    v2r4_assert(buffer$bytes == as.numeric(row$bytes[[1L]]) &&
                  buffer$sha256 == row$sha256[[1L]],
                "SNAPSHOT_OBJECT_IDENTITY_MISMATCH:%s", row$role[[1L]])
    buffer
  })
  names(captured) <- capture_rows$role
  listed_after <- list.files(objects_root, recursive = TRUE, full.names = TRUE,
                             all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  listed_after <- listed_after[file.info(listed_after)$isdir %in% FALSE]
  listed_after_rel <- substring(gsub("\\\\", "/", listed_after),
                                nchar(gsub("\\\\", "/", root)) + 2L)
  v2r4_assert(identical(sort(listed_rel), sort(listed_after_rel)),
              "SNAPSHOT_OBJECT_SET_CHANGED_DURING_CAPTURE")
  verification <- data.table::rbindlist(list(
    capture_rows[, .(role, materialization, relative_path, expected_bytes = bytes,
                observed_bytes = vapply(captured, `[[`, numeric(1), "bytes"),
                expected_sha256 = sha256,
                observed_sha256 = vapply(captured, `[[`, character(1), "sha256"),
                status = "PASS_CAPTURED_ONCE")],
    deferred_rows[, .(role, materialization, relative_path, expected_bytes = bytes,
                      observed_bytes = NA_real_, expected_sha256 = sha256,
                      observed_sha256 = NA_character_, status = "DEFERRED_NOT_READ_BEFORE_SPECIES_GATE")],
    identity_only[, .(role, materialization, relative_path, expected_bytes = bytes,
                      observed_bytes = NA_real_, expected_sha256 = sha256,
                      observed_sha256 = NA_character_, status = "PASS_SIGNED_IDENTITY_ONLY")]
  ), fill = TRUE)
  structure(list(root = root, authority = authority, objects = captured,
                 tree_digest = if (nrow(deferred_rows)) NA_character_ else v2r4_snapshot_tree_identity(content)$sha256,
                 verification = verification,
                 canonical_tree_bytes = if (nrow(deferred_rows)) raw() else v2r4_snapshot_tree_identity(content)$bytes,
                 deferred_roles = deferred_rows$role,
                 initial_physical_paths = sort(listed_rel)),
            class = "v2r4_snapshot")
}

v2r4_complete_snapshot <- function(snapshot) {
  v2r4_assert(inherits(snapshot, "v2r4_snapshot"), "SNAPSHOT_OBJECT_REQUIRED")
  if (!length(snapshot$deferred_roles)) return(snapshot)
  content <- snapshot$authority[base::which(snapshot$authority[["materialization"]] == "CONTENT_OBJECT")]
  deferred <- content[base::which(content[["role"]] %in% snapshot$deferred_roles)]
  additions <- lapply(seq_len(nrow(deferred)), function(i) {
    row <- deferred[i]
    buffer <- v2r4_capture_file(file.path(snapshot$root, row$relative_path[[1L]]),
                                paste0("SNAPSHOT_OBJECT_", row$role[[1L]]), snapshot$root)
    v2r4_assert(buffer$bytes == as.numeric(row$bytes[[1L]]) && buffer$sha256 == row$sha256[[1L]],
                "SNAPSHOT_OBJECT_IDENTITY_MISMATCH:%s", row$role[[1L]])
    buffer
  })
  names(additions) <- deferred$role
  snapshot$objects <- c(snapshot$objects, additions)
  objects_root <- v2r4_norm_existing(file.path(snapshot$root, "objects", "sha256"))
  listed_after <- list.files(objects_root, recursive = TRUE, full.names = TRUE,
                             all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  listed_after <- listed_after[file.info(listed_after)$isdir %in% FALSE]
  listed_after_rel <- substring(gsub("\\\\", "/", listed_after),
                                nchar(gsub("\\\\", "/", snapshot$root)) + 2L)
  v2r4_assert(identical(snapshot$initial_physical_paths, sort(listed_after_rel)),
              "SNAPSHOT_OBJECT_SET_CHANGED_BEFORE_COMPLETION")
  for (role_name in names(additions)) {
    idx <- base::which(snapshot$verification[["role"]] == role_name)
    snapshot$verification$observed_bytes[idx] <- additions[[role_name]]$bytes
    snapshot$verification$observed_sha256[idx] <- additions[[role_name]]$sha256
    snapshot$verification$status[idx] <- "PASS_CAPTURED_ONCE_AFTER_SPECIES_GATE"
  }
  tree <- v2r4_snapshot_tree_identity(content)
  snapshot$tree_digest <- tree$sha256
  snapshot$canonical_tree_bytes <- tree$bytes
  snapshot$deferred_roles <- character()
  snapshot
}

v2r4_object <- function(snapshot, role) {
  v2r4_assert(inherits(snapshot, "v2r4_snapshot"), "SNAPSHOT_OBJECT_REQUIRED")
  role_value <- as.character(role)
  v2r4_assert(role_value %in% names(snapshot$objects), "SNAPSHOT_ROLE_NOT_CAPTURED:%s", role_value)
  snapshot$objects[[role_value]]
}

v2r4_validate_producer_source_buffers <- function(snapshot, cfg) {
  roles <- c("align_count_source", "assembly_source",
             "species_sample_source", "species_cohort_source")
  contract <- cfg$input_snapshot$producer_source_identities
  v2r4_assert(identical(names(contract), roles),
              "PRODUCER_SOURCE_CONTRACT_ROLE_ORDER_INVALID")
  lines <- vapply(roles, function(role) {
    buffer <- v2r4_object(snapshot, role)
    entry <- contract[[role]]
    v2r4_assert(identical(buffer$bytes, as.numeric(entry$bytes)) &&
                  identical(buffer$sha256, entry$sha256) &&
                  identical(basename(buffer$canonical_path), basename(entry$relative_path)),
                paste0("PRODUCER_SOURCE_BUFFER_IDENTITY_INVALID:", role))
    paste(role, format(buffer$bytes, scientific = FALSE, trim = TRUE),
          buffer$sha256, sep = "\t")
  }, character(1))
  text <- paste0(paste(lines, collapse = "\n"), "\n")
  list(
    status = "PASS_4_SIGNED_CONTENT_OBJECTS_BEFORE_COUNTS_CAPTURE",
    digest = tolower(digest::digest(charToRaw(enc2utf8(text)), algo = "sha256",
                                    serialize = FALSE))
  )
}

v2r4_pwsh <- function() {
  executable <- Sys.getenv("GSE205906_V2_R4_PWSH_PATH", unset = "")
  expected_sha256 <- Sys.getenv("GSE205906_V2_R4_PWSH_SHA256", unset = "")
  v2r4_assert(nzchar(executable) && grepl("^[0-9a-f]{64}$", expected_sha256),
              "EXTERNAL_LAUNCH_PINNED_PWSH_REQUIRED")
  buffer <- v2r4_capture_file(executable, "EXTERNAL_LAUNCH_PINNED_PWSH")
  v2r4_assert(identical(buffer$sha256, expected_sha256),
              "EXTERNAL_LAUNCH_PINNED_PWSH_IDENTITY_MISMATCH")
  buffer$canonical_path
}

v2r4_process_identity <- function(guard_buffer, pid) {
  v2r4_assert(inherits(guard_buffer, "v2r4_captured_bytes"),
              "OS_GUARD_CAPTURED_BUFFER_REQUIRED")
  result <- suppressWarnings(system2(v2r4_pwsh(), c(
    "-NoProfile", "-NonInteractive", "-File", shQuote(guard_buffer$canonical_path),
    "-SourceSha256", guard_buffer$sha256,
    "-Command", "ProcessInfo", "-TargetPid", as.character(as.integer(pid))
  ), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status") %||% 0L
  v2r4_assert(as.integer(status) == 0L && length(result) == 1L,
              "LEASE_PROCESS_NOT_LIVE:%s", paste(result, collapse = " | "))
  jsonlite::fromJSON(result[[1L]], simplifyVector = TRUE)
}

v2r4_assert_utc_fixed6 <- function(values, label) {
  values <- as.character(values)
  pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\\.[0-9]{6}Z$"
  v2r4_assert(length(values) > 0L && !anyNA(values) && all(grepl(pattern, values)),
              "%s_UTC_FIXED6_FORMAT_INVALID", label)
  invisible(TRUE)
}

v2r4_assert_process_identity_exact <- function(observed, expected_pid, expected_start_utc) {
  v2r4_assert(is.list(observed) && isTRUE(observed$live) &&
                identical(as.integer(observed$pid), as.integer(expected_pid)) &&
                identical(observed$start_utc, as.character(expected_start_utc)),
              "LEASE_MODEL_PROCESS_IDENTITY_MISMATCH")
  invisible(TRUE)
}

v2r4_probe_mutex <- function(guard_buffer, allowed_root, mutex_path) {
  v2r4_assert(inherits(guard_buffer, "v2r4_captured_bytes"),
              "OS_GUARD_CAPTURED_BUFFER_REQUIRED")
  result <- suppressWarnings(system2(v2r4_pwsh(), c(
    "-NoProfile", "-NonInteractive", "-File", shQuote(guard_buffer$canonical_path),
    "-SourceSha256", guard_buffer$sha256,
    "-Command", "Probe", "-AllowedRoot", shQuote(allowed_root),
    "-LockPath", shQuote(mutex_path)
  ), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status") %||% 0L
  v2r4_assert(as.integer(status) == 0L && identical(result, "HELD"),
              "LEASE_MUTEX_NOT_HELD:%s", paste(result, collapse = " | "))
  invisible(TRUE)
}

v2r4_wait_for_files <- function(paths, timeout_seconds, label) {
  deadline <- Sys.time() + as.numeric(timeout_seconds)
  repeat {
    if (all(file.exists(paths))) return(invisible(TRUE))
    if (Sys.time() >= deadline) {
      v2r4_stop("%s_TIMEOUT_MISSING:%s", label, paste(basename(paths[!file.exists(paths)]), collapse = ","))
    }
    Sys.sleep(0.10)
  }
}

v2r4_write_json_atomic <- function(value, path, required_root) {
  root <- v2r4_norm_existing(required_root)
  target <- v2r4_norm_target(path)
  v2r4_require_within(target, root, "ATOMIC_JSON_TARGET")
  v2r4_assert(!file.exists(target), "ATOMIC_JSON_TARGET_ALREADY_EXISTS:%s", basename(target))
  bytes <- charToRaw(paste0(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null",
                                             digits = NA), "\n"))
  temporary <- file.path(dirname(target), paste0(".", basename(target), ".", Sys.getpid(),
                                                 ".", paste(sample(c(0:9, letters[1:6]), 32,
                                                                    replace = TRUE), collapse = ""),
                                                 ".tmp"))
  con <- file(temporary, open = "wb")
  on.exit({
    try(close(con), silent = TRUE)
    if (file.exists(temporary)) unlink(temporary, force = TRUE)
  }, add = TRUE)
  writeBin(bytes, con)
  close(con)
  con <- NULL
  v2r4_assert(file.rename(temporary, target), "ATOMIC_JSON_PUBLICATION_FAILED:%s", basename(target))
  v2r4_capture_file(target, "ATOMIC_JSON_PUBLISHED", root)
}

v2r4_now_fixed6 <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

v2r4_guard_challenge <- function(guard_state, guard_buffer, cfg,
                                 challenger_pid = Sys.getpid()) {
  v2r4_assert(is.list(guard_state) && inherits(guard_buffer, "v2r4_captured_bytes"),
              "GUARD_CHALLENGE_STATE_INVALID")
  challenger <- v2r4_process_identity(guard_buffer, challenger_pid)
  challenge_nonce <- paste(sprintf("%02x", as.integer(openssl::rand_bytes(32L))), collapse = "")
  stem <- paste0("challenge_", challenger_pid, "_", challenge_nonce)
  request_path <- file.path(guard_state$challenge_directory, paste0(stem, ".request.json"))
  response_path <- file.path(guard_state$challenge_directory, paste0(stem, ".response.json"))
  request <- list(
    schema_version = "GSE205906_bulk_model_v2_r4_guard_challenge_v1",
    run_id = guard_state$run_id,
    guard_pid = as.integer(guard_state$guard_pid),
    guard_start_utc = guard_state$guard_start_utc,
    guard_nonce = guard_state$guard_nonce,
    challenger_pid = as.integer(challenger$pid),
    challenger_start_utc = challenger$start_utc,
    challenge_nonce = challenge_nonce,
    issued_utc = v2r4_now_fixed6()
  )
  request_buffer <- v2r4_write_json_atomic(request, request_path, guard_state$root)
  v2r4_wait_for_files(response_path, as.numeric(cfg$lease$challenge_max_age_seconds),
                      "GUARD_CHALLENGE_RESPONSE")
  response_buffer <- v2r4_capture_file(response_path, "GUARD_CHALLENGE_RESPONSE", guard_state$root)
  envelope <- v2r4_read_json_raw(response_buffer, "GUARD_CHALLENGE_RESPONSE")
  v2r4_assert(identical(names(envelope),
                        c("schema_version", "signed_payload_b64", "signature_b64")) &&
                identical(envelope$schema_version,
                          "GSE205906_bulk_model_v2_r4_guard_response_envelope_v1"),
              "GUARD_CHALLENGE_ENVELOPE_INVALID")
  payload_raw <- tryCatch(jsonlite::base64_dec(envelope$signed_payload_b64),
                          error = function(e) raw())
  signature_raw <- tryCatch(jsonlite::base64_dec(envelope$signature_b64),
                            error = function(e) raw())
  v2r4_assert(length(payload_raw) > 0L && length(signature_raw) > 0L,
              "GUARD_CHALLENGE_BASE64_INVALID")
  public_pem_raw <- tryCatch(jsonlite::base64_dec(guard_state$guard_public_key_pem_b64),
                             error = function(e) raw())
  v2r4_assert(length(public_pem_raw) > 0L &&
                identical(v2r4_sha256_raw(public_pem_raw), guard_state$guard_public_key_sha256),
              "GUARD_CHALLENGE_PUBLIC_KEY_IDENTITY_INVALID")
  public_key <- tryCatch(openssl::read_pubkey(rawToChar(public_pem_raw)),
                         error = function(e) NULL)
  verified <- !is.null(public_key) && tryCatch(
    openssl::signature_verify(payload_raw, signature_raw, hash = openssl::sha256,
                              pubkey = public_key), error = function(e) FALSE
  )
  v2r4_assert(isTRUE(verified), "GUARD_CHALLENGE_SIGNATURE_INVALID")
  payload <- tryCatch(jsonlite::fromJSON(rawToChar(payload_raw), simplifyVector = TRUE),
                      error = function(e) NULL)
  required <- c("schema_version", "status", "run_id", "request_sha256", "challenge_nonce",
                "guard_pid", "guard_start_utc", "guard_nonce", "mutex_file_id",
                "challenger_pid", "challenger_start_utc", "responded_utc")
  v2r4_assert(is.list(payload) && identical(sort(names(payload)), sort(required)),
              "GUARD_CHALLENGE_PAYLOAD_SCHEMA_INVALID")
  v2r4_assert_utc_fixed6(payload$responded_utc, "GUARD_CHALLENGE_RESPONSE_TIME")
  response_time <- as.POSIXct(payload$responded_utc, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
  v2r4_assert(identical(payload$schema_version,
                        "GSE205906_bulk_model_v2_r4_guard_response_payload_v1") &&
                identical(payload$status, "LIVE_MUTEX_HELD_CHALLENGE_SIGNED") &&
                identical(payload$run_id, guard_state$run_id) &&
                identical(payload$request_sha256, request_buffer$sha256) &&
                identical(payload$challenge_nonce, challenge_nonce) &&
                identical(as.integer(payload$guard_pid), as.integer(guard_state$guard_pid)) &&
                identical(payload$guard_start_utc, guard_state$guard_start_utc) &&
                identical(payload$guard_nonce, guard_state$guard_nonce) &&
                identical(payload$mutex_file_id, guard_state$mutex_file_id) &&
                identical(as.integer(payload$challenger_pid), as.integer(challenger$pid)) &&
                identical(payload$challenger_start_utc, challenger$start_utc) &&
                is.finite(as.numeric(response_time)) &&
                abs(as.numeric(Sys.time()) - as.numeric(response_time)) <=
                  as.numeric(cfg$lease$challenge_max_age_seconds),
              "GUARD_CHALLENGE_RESPONSE_BINDING_INVALID")
  list(request = request_buffer, response = response_buffer,
       response_payload_sha256 = v2r4_sha256_raw(payload_raw), payload = payload)
}

v2r4_start_mutex_guard <- function(...) {
  v2r4_stop("GUARD_MUST_BE_STARTED_AND_OWNED_BY_CAPTURED_EXTERNAL_LAUNCHER")
}

v2r4_stop_mutex_guard <- function(...) {
  v2r4_stop("MODEL_OR_VALIDATOR_MAY_NOT_STOP_LAUNCHER_OWNED_GUARD")
}

v2r4_attach_existing_guard <- function(guard_buffer, state_root, run_id, cfg,
                                       expected_owner_pid, expected_owner_start_utc,
                                       expected_stop_hmac_key_sha256,
                                       expected_finalization_public_key_sha256,
                                       expected_finalization_public_key_captured_path,
                                       expected_finalization_public_key_captured_bytes,
                                       expected_finalization_public_key_captured_sha256,
                                       expected_ready_path, expected_stop_path,
                                       expected_stopped_path) {
  state_root <- v2r4_norm_existing(state_root)
  state_dir <- v2r4_norm_existing(file.path(state_root, "state"))
  ready_path <- v2r4_norm_existing(expected_ready_path)
  challenge_directory <- v2r4_norm_existing(file.path(
    state_dir, paste0("guard_challenges_", run_id)
  ))
  ready_buffer <- v2r4_capture_file(ready_path, "EXISTING_GUARD_READY_RECEIPT", state_root)
  ready <- v2r4_read_json_raw(ready_buffer, "EXISTING_GUARD_READY_RECEIPT")
  captured_key_path <- v2r4_norm_existing(expected_finalization_public_key_captured_path)
  captured_key <- v2r4_capture_file(
    captured_key_path, "EXISTING_GUARD_FINALIZATION_PUBLIC_KEY_COPY",
    dirname(state_root)
  )
  required <- v2r4_guard_ready_fields()
  v2r4_assert(identical(names(ready), required) &&
                identical(ready$schema_version, "GSE205906_bulk_model_v2_r4_guard_ready_v1") &&
                identical(ready$status, cfg$lease$guard_ready_status) &&
                identical(ready$run_id, run_id) &&
                identical(ready$owner_role, "EXTERNAL_LAUNCHER") &&
                identical(as.integer(ready$owner_pid), as.integer(expected_owner_pid)) &&
                identical(ready$owner_start_utc, expected_owner_start_utc) &&
                identical(ready$stop_hmac_key_sha256, expected_stop_hmac_key_sha256) &&
                identical(ready$finalization_public_key_sha256,
                          expected_finalization_public_key_sha256) &&
                identical(v2r4_norm_existing(ready$finalization_public_key_captured_path),
                          captured_key_path) &&
                identical(as.numeric(ready$finalization_public_key_captured_bytes),
                          as.numeric(expected_finalization_public_key_captured_bytes)) &&
                identical(ready$finalization_public_key_captured_sha256,
                          expected_finalization_public_key_captured_sha256) &&
                identical(as.numeric(captured_key$bytes),
                          as.numeric(expected_finalization_public_key_captured_bytes)) &&
                identical(captured_key$sha256,
                          expected_finalization_public_key_captured_sha256) &&
                identical(captured_key$sha256, expected_finalization_public_key_sha256) &&
                identical(ready$source_sha256, guard_buffer$sha256) &&
                identical(v2r4_norm_existing(ready$challenge_directory), challenge_directory) &&
                identical(v2r4_norm_target(ready$stop_request_path),
                          v2r4_norm_target(expected_stop_path)) &&
                identical(v2r4_norm_target(ready$stopped_receipt_path),
                          v2r4_norm_target(expected_stopped_path)),
              "EXISTING_GUARD_READY_BINDING_INVALID")
  v2r4_assert_utc_fixed6(c(ready$owner_start_utc, ready$guard_start_utc),
                         "EXISTING_GUARD_READY_TIMES")
  owner <- v2r4_process_identity(guard_buffer, expected_owner_pid)
  guard <- v2r4_process_identity(guard_buffer, as.integer(ready$guard_pid))
  v2r4_assert_process_identity_exact(owner, expected_owner_pid, expected_owner_start_utc)
  v2r4_assert_process_identity_exact(guard, as.integer(ready$guard_pid), ready$guard_start_utc)
  mutex_path <- v2r4_norm_existing(ready$mutex_path)
  v2r4_probe_mutex(guard_buffer, state_root, mutex_path)
  list(root = state_root, mutex_path = mutex_path, ready_path = ready_path,
       ready_receipt_sha256 = ready_buffer$sha256,
       challenge_directory = challenge_directory, owner_role = ready$owner_role,
       owner_pid = as.integer(ready$owner_pid), owner_start_utc = ready$owner_start_utc,
       stop_path = v2r4_norm_target(ready$stop_request_path),
        stopped_receipt_path = v2r4_norm_target(ready$stopped_receipt_path),
        stop_hmac_key_sha256 = ready$stop_hmac_key_sha256,
        finalization_public_key_sha256 = ready$finalization_public_key_sha256,
        finalization_public_key_captured_path = captured_key_path,
        finalization_public_key_captured_bytes = as.integer(captured_key$bytes),
        finalization_public_key_captured_sha256 = captured_key$sha256,
       guard_pid = as.integer(ready$guard_pid), guard_start_utc = ready$guard_start_utc,
       guard_nonce = ready$guard_nonce,
       guard_public_key_pem_b64 = ready$guard_public_key_pem_b64,
       guard_public_key_sha256 = ready$guard_public_key_sha256,
       mutex_file_id = ready$mutex_file_id, run_id = run_id)
}

v2r4_lease_identity_request_fields <- function() {
  c(
    "schema_version", "status", "role", "dataset_id", "model_version", "run_id",
    "holder_principal", "process_pid", "process_start_utc", "child_request_sha256",
    "child_request_nonce", "peer_request_nonce", "snapshot_tree_sha256",
    "runtime_role_digest", "source_authority_payload_sha256",
    "external_launch_context_sha256", "runtime_pre_digest", "guard_owner_role",
    "guard_owner_pid", "guard_owner_start_utc", "guard_pid", "guard_start_utc",
    "guard_nonce", "guard_public_key_sha256", "guard_ready_receipt_sha256",
    "guard_challenge_request_sha256", "guard_challenge_response_sha256",
    "guard_mutex_file_id", "guard_stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_path",
    "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256",
    "issued_utc",
    "e2_gate_status", "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_combined_lease_request_fields <- function() {
  c(
    "schema_version", "dataset_id", "model_version", "decision", "purpose",
    "request_status", "holder_principal", "run_id", "model_pid", "model_start_utc",
    "model_request_sha256", "model_request_nonce", "model_identity_request_sha256",
    "validator_pid", "validator_start_utc", "validator_request_sha256",
    "validator_request_nonce", "validator_identity_request_sha256",
    "snapshot_tree_sha256", "mutex_relative_path", "runtime_role_digest",
    "source_authority_payload_sha256", "external_launch_context_sha256",
    "runtime_pre_digest", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_pid", "guard_start_utc", "guard_nonce",
    "guard_public_key_sha256", "guard_ready_receipt_sha256",
    "guard_challenge_request_sha256", "guard_challenge_response_sha256",
    "guard_mutex_file_id", "guard_stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_path",
    "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256"
  )
}

v2r4_signed_lease_fields <- function() {
  c(
    "schema_version", "dataset_id", "model_version", "decision", "purpose",
    "lease_id", "status", "holder_principal", "issued_utc", "expires_utc", "run_id",
    "model_pid", "model_start_utc", "model_request_sha256", "model_request_nonce",
    "model_identity_request_sha256", "validator_pid", "validator_start_utc",
    "validator_request_sha256", "validator_request_nonce",
    "validator_identity_request_sha256", "snapshot_tree_sha256", "mutex_relative_path",
    "runtime_role_digest", "source_authority_payload_sha256",
    "external_launch_context_sha256", "runtime_pre_digest", "guard_owner_role",
    "guard_owner_pid", "guard_owner_start_utc", "guard_pid", "guard_start_utc",
    "guard_nonce", "guard_public_key_sha256", "guard_ready_receipt_sha256",
    "guard_challenge_request_sha256", "guard_challenge_response_sha256",
    "guard_mutex_file_id", "guard_stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_path",
    "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256"
  )
}

v2r4_assert_child_stop_key_absent <- function() {
  secret_name <- "GSE205906_V2_R4_GUARD_STOP_HMAC_KEY"
  secret_value <- Sys.getenv(secret_name, unset = NA_character_)
  argv <- commandArgs(FALSE)
  v2r4_assert(length(secret_value) == 1L && is.na(secret_value),
              "GUARD_STOP_HMAC_KEY_PRESENT_IN_MODEL_OR_VALIDATOR_ENVIRONMENT")
  v2r4_assert(!any(grepl(secret_name, argv, fixed = TRUE)),
              "GUARD_STOP_HMAC_KEY_NAME_PRESENT_IN_MODEL_OR_VALIDATOR_ARGV")
  invisible(TRUE)
}

v2r4_authorized_snapshot_tree_digest <- function(authority) {
  content <- authority[base::which(authority[["materialization"]] == "CONTENT_OBJECT")]
  v2r4_assert(nrow(content) > 0L && !anyDuplicated(content[["relative_path"]]),
              "AUTHORIZED_SNAPSHOT_TREE_DOMAIN_INVALID")
  tree <- v2r4_snapshot_tree_identity(content)
  v2r4_assert(grepl("^[0-9a-f]{64}$", tree$sha256),
              "AUTHORIZED_SNAPSHOT_TREE_DIGEST_INVALID")
  tree$sha256
}

v2r4_runtime_role_digest <- function(runtime_bindings) {
  runtime_bindings <- data.table::as.data.table(runtime_bindings)
  v2r4_assert(identical(names(runtime_bindings)[1:4],
                        c("role", "runtime_relative_path", "bytes", "sha256")) &&
                nrow(runtime_bindings) > 0L && !anyDuplicated(runtime_bindings$role) &&
                all(grepl("^[0-9a-f]{64}$", runtime_bindings$sha256)),
              "RUNTIME_ROLE_DIGEST_INPUT_INVALID")
  tolower(digest::digest(
    charToRaw(paste(runtime_bindings$role, runtime_bindings$sha256,
                    sep = "\t", collapse = "\n")),
    algo = "sha256", serialize = FALSE
  ))
}

v2r4_write_lease_identity_request <- function(path, cfg, role, run_id,
                                               holder_principal,
                                               snapshot_tree_digest, guard_state,
                                               runtime_role_digest, launch_context,
                                               process_identity, child_request,
                                               root) {
  v2r4_assert(role %in% c("MODEL", "VALIDATOR"),
              "LEASE_IDENTITY_REQUEST_ROLE_INVALID")
  v2r4_assert(is.list(child_request) && inherits(child_request$buffer,
                                                  "v2r4_captured_bytes"),
              "LEASE_IDENTITY_REQUEST_CHILD_REQUEST_INVALID")
  value <- list(
    schema_version = cfg$lease$identity_request_schema_version,
    status = cfg$lease$identity_request_status,
    role = role,
    dataset_id = cfg$dataset_id,
    model_version = cfg$model_version,
    run_id = run_id,
    holder_principal = holder_principal,
    process_pid = as.integer(process_identity$pid),
    process_start_utc = process_identity$start_utc,
    child_request_sha256 = child_request$buffer$sha256,
    child_request_nonce = child_request$request$request_nonce,
    peer_request_nonce = child_request$request$peer_request_nonce,
    snapshot_tree_sha256 = snapshot_tree_digest,
    runtime_role_digest = runtime_role_digest,
    source_authority_payload_sha256 = launch_context$source_authorization_payload_sha256,
    external_launch_context_sha256 = launch_context$context_sha256,
    runtime_pre_digest = launch_context$runtime_pre_digest,
    guard_owner_role = guard_state$owner_role,
    guard_owner_pid = as.integer(guard_state$owner_pid),
    guard_owner_start_utc = guard_state$owner_start_utc,
    guard_pid = as.integer(guard_state$guard_pid),
    guard_start_utc = guard_state$guard_start_utc,
    guard_nonce = guard_state$guard_nonce,
    guard_public_key_sha256 = guard_state$guard_public_key_sha256,
    guard_ready_receipt_sha256 = guard_state$ready_receipt_sha256,
    guard_challenge_request_sha256 = guard_state$initial_challenge$request$sha256,
    guard_challenge_response_sha256 = guard_state$initial_challenge$response$sha256,
    guard_mutex_file_id = guard_state$mutex_file_id,
    guard_stop_hmac_key_sha256 = guard_state$stop_hmac_key_sha256,
    finalization_public_key_sha256 = guard_state$finalization_public_key_sha256,
    finalization_public_key_captured_path = guard_state$finalization_public_key_captured_path,
    finalization_public_key_captured_bytes =
      as.integer(guard_state$finalization_public_key_captured_bytes),
    finalization_public_key_captured_sha256 =
      guard_state$finalization_public_key_captured_sha256,
    issued_utc = v2r4_now_fixed6(),
    e2_gate_status = cfg$gate$machine_status,
    candidate_gate_eligible = FALSE,
    can_generate_or_upgrade_candidate = FALSE,
    formal_sequential_inference = FALSE
  )
  v2r4_assert(identical(names(value), v2r4_lease_identity_request_fields()),
              "LEASE_IDENTITY_REQUEST_FIELD_ORDER_INVALID")
  v2r4_write_json_atomic(value, path, root)
}

v2r4_validate_lease_identity_request <- function(buffer, cfg, expected_role,
                                                   expected, guard_buffer,
                                                   now_utc = Sys.time()) {
  value <- v2r4_read_json_raw(buffer, paste0(expected_role, "_LEASE_IDENTITY_REQUEST"))
  v2r4_assert(identical(names(value), v2r4_lease_identity_request_fields()) &&
                identical(value$schema_version, cfg$lease$identity_request_schema_version) &&
                identical(value$status, cfg$lease$identity_request_status) &&
                identical(value$role, expected_role) &&
                identical(value$dataset_id, cfg$dataset_id) &&
                identical(value$model_version, cfg$model_version) &&
                identical(value$e2_gate_status, cfg$gate$machine_status) &&
                identical(value$candidate_gate_eligible, FALSE) &&
                identical(value$can_generate_or_upgrade_candidate, FALSE) &&
                identical(value$formal_sequential_inference, FALSE),
              "LEASE_IDENTITY_REQUEST_SCHEMA_OR_SEMANTICS_INVALID:%s", expected_role)
  v2r4_assert_utc_fixed6(c(value$process_start_utc, value$guard_owner_start_utc,
                           value$guard_start_utc, value$issued_utc),
                         paste0(expected_role, "_LEASE_IDENTITY_REQUEST_TIMES"))
  hash_fields <- c(
    "child_request_sha256", "child_request_nonce", "peer_request_nonce",
    "snapshot_tree_sha256", "runtime_role_digest", "source_authority_payload_sha256",
    "external_launch_context_sha256", "runtime_pre_digest", "guard_nonce",
    "guard_public_key_sha256", "guard_ready_receipt_sha256",
    "guard_challenge_request_sha256", "guard_challenge_response_sha256",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_sha256"
  )
  v2r4_assert(all(vapply(value[hash_fields], function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[0-9a-f]{64}$", x)
  }, logical(1L))), "LEASE_IDENTITY_REQUEST_HASH_INVALID:%s", expected_role)
  v2r4_assert(is.character(value$guard_mutex_file_id) &&
                length(value$guard_mutex_file_id) == 1L &&
                grepl("^[0-9a-f]{8}:[0-9a-f]{16}$", value$guard_mutex_file_id),
              "LEASE_IDENTITY_REQUEST_MUTEX_FILE_ID_INVALID:%s", expected_role)
  v2r4_assert(is.character(value$finalization_public_key_captured_path) &&
                length(value$finalization_public_key_captured_path) == 1L &&
                nzchar(value$finalization_public_key_captured_path) &&
                identical(as.integer(value$finalization_public_key_captured_bytes),
                          as.integer(expected$finalization_public_key_captured_bytes)),
              "LEASE_IDENTITY_REQUEST_FINALIZATION_KEY_COPY_INVALID:%s", expected_role)
  issued <- as.POSIXct(value$issued_utc, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
  v2r4_assert(is.finite(as.numeric(issued)) &&
                abs(as.numeric(as.POSIXct(now_utc, tz = "UTC")) - as.numeric(issued)) <= 900,
              "LEASE_IDENTITY_REQUEST_NOT_FRESH:%s", expected_role)
  v2r4_assert_named_expected(value, expected,
                             paste0(expected_role, "_LEASE_IDENTITY_REQUEST"))
  identity <- v2r4_process_identity(guard_buffer, as.integer(value$process_pid))
  v2r4_assert_process_identity_exact(identity, as.integer(value$process_pid),
                                     value$process_start_utc)
  value
}

v2r4_write_lease_request <- function(...) {
  v2r4_stop("LEASE_V4_SINGLE_MODEL_REQUEST_FORBIDDEN_USE_DUAL_IDENTITY_COMBINED_REQUEST")
}

v2r4_validate_lease_buffer <- function(lease_buffer, cfg, snapshot_tree_digest,
                                        run_id, holder_principal, guard_buffer,
                                        guard_state, launch_context, runtime_role_digest,
                                        mutex_root, expected_model_pid,
                                        expected_model_start_utc,
                                        expected_model_request_sha256,
                                        expected_model_request_nonce,
                                        expected_model_identity_request_sha256,
                                        expected_validator_pid,
                                        expected_validator_start_utc,
                                        expected_validator_request_sha256,
                                        expected_validator_request_nonce,
                                        expected_validator_identity_request_sha256,
                                        now_utc = Sys.time()) {
  lease <- v2r4_read_tsv_raw(
    lease_buffer, "SIGNED_LEASE",
    character_columns = c("issued_utc", "expires_utc", "model_start_utc",
                          "validator_start_utc", "guard_owner_start_utc",
                          "guard_start_utc")
  )
  required <- v2r4_signed_lease_fields()
  v2r4_assert(identical(names(lease), required) && nrow(lease) == 1L,
              "LEASE_SCHEMA_NOT_EXACT")
  v2r4_assert(lease$schema_version[[1L]] == cfg$lease$schema_version &&
                lease$dataset_id[[1L]] == cfg$dataset_id &&
                lease$model_version[[1L]] == cfg$model_version &&
                lease$decision[[1L]] == cfg$lease$decision &&
                lease$purpose[[1L]] == cfg$lease$purpose &&
                lease$status[[1L]] == "ACTIVE_EXCLUSIVE_IMMUTABLE_SNAPSHOT",
              "LEASE_HEADER_SEMANTICS_INVALID")
  v2r4_assert(lease$holder_principal[[1L]] == holder_principal &&
                lease$run_id[[1L]] == run_id &&
                lease$snapshot_tree_sha256[[1L]] == snapshot_tree_digest &&
                lease$runtime_role_digest[[1L]] == runtime_role_digest &&
                lease$source_authority_payload_sha256[[1L]] ==
                  launch_context$source_authorization_payload_sha256 &&
                lease$external_launch_context_sha256[[1L]] == launch_context$context_sha256 &&
                lease$runtime_pre_digest[[1L]] == launch_context$runtime_pre_digest,
              "LEASE_HOLDER_RUN_OR_TREE_MISMATCH")
  v2r4_assert(grepl("^[0-9a-f]{64}$", lease$lease_id[[1L]]),
              "LEASE_ID_INVALID")
  v2r4_assert_utc_fixed6(
    c(lease$issued_utc[[1L]], lease$expires_utc[[1L]], lease$model_start_utc[[1L]],
      lease$validator_start_utc[[1L]], lease$guard_owner_start_utc[[1L]],
      lease$guard_start_utc[[1L]]),
    "LEASE_TIMESTAMPS"
  )
  issued <- as.POSIXct(lease$issued_utc[[1L]], format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
  expires <- as.POSIXct(lease$expires_utc[[1L]], format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
  now_num <- as.numeric(as.POSIXct(now_utc, tz = "UTC"))
  v2r4_assert(is.finite(as.numeric(issued)) && is.finite(as.numeric(expires)) &&
                as.numeric(issued) < as.numeric(expires) &&
                as.numeric(expires) - as.numeric(issued) <=
                  as.numeric(cfg$lease$maximum_duration_seconds) &&
                as.numeric(issued) <= now_num && now_num < as.numeric(expires),
              "LEASE_TIME_WINDOW_INVALID")
  hash_fields <- c(
    "model_request_sha256", "model_request_nonce", "model_identity_request_sha256",
    "validator_request_sha256", "validator_request_nonce",
    "validator_identity_request_sha256", "snapshot_tree_sha256", "runtime_role_digest",
    "source_authority_payload_sha256", "external_launch_context_sha256",
    "runtime_pre_digest", "guard_nonce", "guard_public_key_sha256",
    "guard_ready_receipt_sha256", "guard_challenge_request_sha256",
    "guard_challenge_response_sha256", "guard_stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_sha256"
  )
  v2r4_assert(all(vapply(hash_fields, function(field) {
    grepl("^[0-9a-f]{64}$", lease[[field]][[1L]])
  }, logical(1))), "LEASE_HASH_OR_NONCE_DOMAIN_INVALID")
  model_pid <- as.integer(lease$model_pid[[1L]])
  v2r4_assert(identical(model_pid, as.integer(expected_model_pid)) &&
                identical(lease$model_start_utc[[1L]], expected_model_start_utc) &&
                identical(lease$model_request_sha256[[1L]],
                          expected_model_request_sha256) &&
                identical(lease$model_request_nonce[[1L]], expected_model_request_nonce) &&
                identical(lease$model_identity_request_sha256[[1L]],
                          expected_model_identity_request_sha256),
              "LEASE_MODEL_SIGNED_IDENTITY_BINDING_INVALID")
  observed <- v2r4_process_identity(guard_buffer, model_pid)
  v2r4_assert_process_identity_exact(observed, model_pid, lease$model_start_utc[[1L]])
  validator_pid <- as.integer(lease$validator_pid[[1L]])
  v2r4_assert(identical(validator_pid, as.integer(expected_validator_pid)) &&
                identical(lease$validator_start_utc[[1L]], expected_validator_start_utc) &&
                identical(lease$validator_request_sha256[[1L]],
                          expected_validator_request_sha256) &&
                identical(lease$validator_request_nonce[[1L]],
                          expected_validator_request_nonce) &&
                identical(lease$validator_identity_request_sha256[[1L]],
                          expected_validator_identity_request_sha256),
              "LEASE_VALIDATOR_SIGNED_IDENTITY_BINDING_INVALID")
  validator_observed <- v2r4_process_identity(guard_buffer, validator_pid)
  v2r4_assert_process_identity_exact(validator_observed, validator_pid,
                                     lease$validator_start_utc[[1L]])
  owner_pid <- as.integer(lease$guard_owner_pid[[1L]])
  v2r4_assert(identical(lease$guard_owner_role[[1L]], "EXTERNAL_LAUNCHER") &&
                identical(owner_pid, as.integer(launch_context$launcher_pid)) &&
                identical(lease$guard_owner_start_utc[[1L]],
                          launch_context$launcher_start_utc) &&
                identical(owner_pid, as.integer(guard_state$owner_pid)) &&
                identical(lease$guard_owner_start_utc[[1L]], guard_state$owner_start_utc),
              "LEASE_GUARD_OWNER_BINDING_INVALID")
  owner_observed <- v2r4_process_identity(guard_buffer, owner_pid)
  v2r4_assert_process_identity_exact(owner_observed, owner_pid,
                                     lease$guard_owner_start_utc[[1L]])
  v2r4_assert(identical(lease$mutex_relative_path[[1L]], cfg$lease$mutex_relative_path),
              "LEASE_MUTEX_IDENTIFIER_MISMATCH")
  v2r4_assert(identical(as.integer(lease$guard_pid[[1L]]), guard_state$guard_pid) &&
                identical(lease$guard_start_utc[[1L]], guard_state$guard_start_utc) &&
                identical(lease$guard_nonce[[1L]], guard_state$guard_nonce) &&
                identical(lease$guard_public_key_sha256[[1L]],
                          guard_state$guard_public_key_sha256) &&
                identical(lease$guard_ready_receipt_sha256[[1L]],
                          guard_state$ready_receipt_sha256) &&
                identical(lease$guard_challenge_request_sha256[[1L]],
                          guard_state$initial_challenge$request$sha256) &&
                identical(lease$guard_challenge_response_sha256[[1L]],
                           guard_state$initial_challenge$response$sha256) &&
                identical(lease$guard_mutex_file_id[[1L]], guard_state$mutex_file_id) &&
                identical(lease$guard_stop_hmac_key_sha256[[1L]],
                          guard_state$stop_hmac_key_sha256) &&
                identical(lease$finalization_public_key_sha256[[1L]],
                          guard_state$finalization_public_key_sha256) &&
                identical(lease$finalization_public_key_captured_path[[1L]],
                          guard_state$finalization_public_key_captured_path) &&
                identical(as.integer(lease$finalization_public_key_captured_bytes[[1L]]),
                          as.integer(guard_state$finalization_public_key_captured_bytes)) &&
                identical(lease$finalization_public_key_captured_sha256[[1L]],
                          guard_state$finalization_public_key_captured_sha256),
              "LEASE_GUARD_IDENTITY_OR_CHALLENGE_BINDING_INVALID")
  guard_observed <- v2r4_process_identity(guard_buffer, guard_state$guard_pid)
  v2r4_assert_process_identity_exact(guard_observed, guard_state$guard_pid,
                                     guard_state$guard_start_utc)
  mutex_root <- v2r4_norm_existing(mutex_root)
  mutex_path <- v2r4_norm_existing(file.path(mutex_root, lease$mutex_relative_path[[1L]]))
  v2r4_probe_mutex(guard_buffer, mutex_root, mutex_path)
  lease
}

v2r4_assert_named_expected <- function(value, expected, label) {
  v2r4_assert(is.list(value) && is.list(expected) && !is.null(names(expected)) &&
                !anyDuplicated(names(expected)), "%s_EXPECTED_SCHEMA_INVALID", label)
  ok <- vapply(names(expected), function(field) {
    observed <- value[[field]]
    target <- expected[[field]]
    if (is.logical(target)) identical(observed, target) else
      identical(as.character(observed), as.character(target))
  }, logical(1))
  v2r4_assert(all(ok), "%s_EXPECTED_BINDING_INVALID:%s", label,
               paste(names(expected)[!ok], collapse = ","))
  invisible(TRUE)
}

v2r4_parse_fixed6_time <- function(x, label) {
  v2r4_assert_utc_fixed6(x, label)
  parsed <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
  v2r4_assert(length(parsed) == 1L && is.finite(as.numeric(parsed)),
              "%s_PARSE_INVALID", label)
  parsed
}

v2r4_finalization_request_fields <- function() {
  c(
    "schema_version", "status", "dataset_id", "model_version", "run_id", "lease_id",
    "lease_manifest_sha256",
    "holder_principal", "source_authority_payload_sha256", "model_source_sha256",
    "validator_source_sha256", "shared_library_sha256", "analysis_config_sha256",
    "external_launcher_sha256", "os_guard_sha256", "launcher_pid", "launcher_start_utc",
    "guard_owner_role", "guard_owner_pid", "guard_owner_start_utc",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "model_pid", "model_start_utc",
    "validator_pid", "validator_start_utc", "guard_pid", "guard_start_utc",
    "guard_nonce", "guard_public_key_sha256", "guard_ready_receipt_sha256",
    "guard_mutex_file_id", "model_request_sha256", "model_request_nonce",
    "model_identity_request_sha256", "validator_request_sha256",
    "validator_request_nonce", "validator_identity_request_sha256", "output_manifest_sha256",
    "complete_sha256", "validation_report_sha256", "snapshot_tree_sha256",
    "runtime_pre_digest", "lease_issued_utc", "lease_expires_utc",
    "validator_guard_challenge_request_sha256",
    "validator_guard_challenge_response_sha256", "guard_ready_receipt_path",
    "guard_challenge_directory", "e2_gate_status", "candidate_gate_eligible",
    "can_generate_or_upgrade_candidate", "formal_sequential_inference"
  )
}

v2r4_science_complete_fields <- function() {
  c(
    "schema_version", "status", "dataset_id", "model_version", "run_id",
    "validator_pid", "validator_start_utc", "validator_request_sha256",
    "validator_request_nonce", "validation_report_sha256", "output_manifest_sha256",
    "complete_sha256", "finalization_request_sha256", "issued_utc", "e2_gate_status",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_runtime_post_receipt_fields <- function() {
  c(
    "schema_version", "status", "run_id", "runtime_pre_digest", "runtime_post_digest",
    "runtime_pre_post_equal", "runtime_semantic_digest",
    "execution_identity_manifest_sha256", "analysis_config_sha256",
    "science_complete_sha256", "validator_pid", "validator_start_utc",
    "validator_request_sha256", "validator_request_nonce", "wait_seconds",
    "wait_issued_utc", "final_ack_deadline_utc", "checked_utc",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_final_challenge_request_fields <- function() {
  c(
    "schema_version", "status", "run_id", "lease_id", "lease_manifest_sha256",
    "source_authority_payload_sha256",
    "launcher_pid", "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_path",
    "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256",
    "model_pid", "model_start_utc", "validator_pid", "validator_start_utc",
    "guard_pid", "guard_start_utc", "guard_nonce", "guard_public_key_sha256",
    "guard_ready_receipt_sha256", "guard_mutex_file_id", "model_request_sha256",
    "model_request_nonce", "model_identity_request_sha256", "validator_request_sha256",
    "validator_request_nonce", "validator_identity_request_sha256",
    "science_complete_sha256", "finalization_request_sha256",
    "validation_report_sha256", "output_manifest_sha256", "complete_sha256",
    "snapshot_tree_sha256", "runtime_post_receipt_sha256", "runtime_pre_digest",
    "runtime_post_digest", "runtime_pre_post_equal", "challenge_nonce", "issued_utc",
    "expires_utc", "e2_gate_status", "candidate_gate_eligible",
    "can_generate_or_upgrade_candidate", "formal_sequential_inference"
  )
}

v2r4_final_challenge_response_fields <- function() {
  c(
    "schema_version", "status", "request_sha256", "run_id", "lease_id",
    "lease_manifest_sha256",
    "source_authority_payload_sha256", "launcher_pid", "launcher_start_utc",
    "guard_owner_role", "guard_owner_pid", "guard_owner_start_utc",
    "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
    "finalization_public_key_captured_path", "finalization_public_key_captured_bytes",
    "finalization_public_key_captured_sha256", "model_pid", "model_start_utc",
    "validator_pid", "validator_start_utc", "guard_pid", "guard_start_utc",
    "guard_nonce", "guard_public_key_sha256", "guard_ready_receipt_sha256",
    "guard_mutex_file_id", "model_request_sha256", "model_request_nonce",
    "model_identity_request_sha256", "validator_request_sha256",
    "validator_request_nonce", "validator_identity_request_sha256",
    "science_complete_sha256",
    "finalization_request_sha256", "validation_report_sha256",
    "output_manifest_sha256", "complete_sha256", "snapshot_tree_sha256",
    "runtime_post_receipt_sha256", "runtime_pre_digest", "runtime_post_digest",
    "runtime_pre_post_equal", "challenge_nonce", "challenge_issued_utc",
    "challenge_expires_utc", "mutex_file_id_current", "responded_utc",
    "e2_gate_status", "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_control_complete_fields <- function() {
  c(
    "schema_version", "status", "dataset_id", "model_version", "run_id", "lease_id",
    "lease_manifest_sha256",
    "launcher_pid", "launcher_start_utc", "guard_owner_role", "guard_owner_pid",
    "guard_owner_start_utc", "guard_stop_hmac_key_sha256",
    "finalization_public_key_sha256", "finalization_public_key_captured_path",
    "finalization_public_key_captured_bytes", "finalization_public_key_captured_sha256",
    "model_pid", "model_start_utc", "validator_pid", "validator_start_utc",
    "guard_pid", "guard_start_utc", "guard_nonce", "guard_ready_receipt_sha256",
    "guard_mutex_file_id", "model_request_sha256", "model_request_nonce",
    "model_identity_request_sha256", "validator_request_sha256",
    "validator_request_nonce", "validator_identity_request_sha256",
    "science_complete_sha256",
    "finalization_request_sha256", "validation_report_sha256",
    "output_manifest_sha256", "complete_sha256", "runtime_post_receipt_sha256",
    "runtime_pre_digest", "runtime_post_digest", "runtime_pre_post_equal",
    "final_challenge_request_sha256", "final_challenge_response_sha256",
    "final_challenge_payload_sha256", "final_challenge_nonce",
    "final_challenge_issued_utc", "final_challenge_responded_utc", "lease_issued_utc",
    "lease_expires_utc", "issued_utc", "e2_gate_status", "candidate_gate_eligible",
    "can_generate_or_upgrade_candidate", "formal_sequential_inference"
  )
}

v2r4_runtime_final_receipt_fields <- function() {
  c(
    "schema_version", "status", "run_id", "runtime_pre_digest", "runtime_post_digest",
    "runtime_final_digest", "runtime_pre_post_final_equal", "final_ack_envelope_sha256",
    "final_ack_payload_sha256", "model_consumer_receipt_sha256",
    "validator_consumer_receipt_sha256", "checked_utc", "e2_gate_status",
    "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
    "formal_sequential_inference"
  )
}

v2r4_consumer_receipt_fields <- function() {
  c(
    "schema_version", "status", "consumer_role", "consumer_pid",
    "consumer_start_utc", "final_ack_envelope_sha256",
    "final_ack_payload_sha256", "holder_principal",
    v2r4_final_ack_binding_fields(), "issued_utc"
  )
}

v2r4_validate_science_complete <- function(path, cfg, expected, root) {
  buffer <- v2r4_capture_file(path, "SCIENCE_COMPLETE", root)
  value <- v2r4_read_json_raw(buffer, "SCIENCE_COMPLETE")
  v2r4_assert(identical(names(value), v2r4_science_complete_fields()) &&
                identical(value$schema_version, cfg$finalization$science_complete_schema_version) &&
                identical(value$status, "SCIENCE_RECOMPUTATION_COMPLETE_PROCESS_HELD") &&
                identical(value$dataset_id, cfg$dataset_id) &&
                identical(value$model_version, cfg$model_version) &&
                identical(value$e2_gate_status, cfg$gate$machine_status) &&
                identical(value$candidate_gate_eligible, FALSE) &&
                identical(value$can_generate_or_upgrade_candidate, FALSE) &&
                identical(value$formal_sequential_inference, FALSE),
              "SCIENCE_COMPLETE_SCHEMA_OR_SEMANTICS_INVALID")
  v2r4_assert_utc_fixed6(c(value$validator_start_utc, value$issued_utc),
                          "SCIENCE_COMPLETE_TIMES")
  v2r4_assert_named_expected(value, expected, "SCIENCE_COMPLETE")
  list(buffer = buffer, value = value)
}

v2r4_validate_runtime_post_receipt <- function(path, cfg, expected, root) {
  buffer <- v2r4_capture_file(path, "RUNTIME_POST_RECEIPT", root)
  value <- v2r4_read_json_raw(buffer, "RUNTIME_POST_RECEIPT")
  v2r4_assert(identical(names(value), v2r4_runtime_post_receipt_fields()) &&
                identical(value$schema_version,
                          cfg$finalization$runtime_post_receipt_schema_version) &&
                identical(value$status,
                          "PASS_RUNTIME_PRE_POST_EQUAL_BEFORE_FINALIZATION_UNDER_GLOBAL_LEASE") &&
                identical(value$runtime_pre_post_equal, TRUE) &&
                identical(value$runtime_post_digest, value$runtime_pre_digest) &&
                identical(value$wait_seconds, 900L) &&
                identical(value$candidate_gate_eligible, FALSE) &&
                identical(value$can_generate_or_upgrade_candidate, FALSE) &&
                identical(value$formal_sequential_inference, FALSE),
              "RUNTIME_POST_RECEIPT_SCHEMA_OR_EQUALITY_INVALID")
  issued <- v2r4_parse_fixed6_time(value$wait_issued_utc, "RUNTIME_POST_WAIT_ISSUED")
  deadline <- v2r4_parse_fixed6_time(value$final_ack_deadline_utc,
                                      "RUNTIME_POST_WAIT_DEADLINE")
  v2r4_assert_utc_fixed6(c(value$validator_start_utc, value$checked_utc),
                          "RUNTIME_POST_RECEIPT_TIMES")
  v2r4_assert(abs(as.numeric(deadline) - as.numeric(issued) - 900) <= 1e-6,
              "RUNTIME_POST_RECEIPT_DEADLINE_NOT_EXACT_900")
  v2r4_assert_named_expected(value, expected, "RUNTIME_POST_RECEIPT")
  list(buffer = buffer, value = value)
}

v2r4_verify_final_challenge <- function(request_path, response_path, cfg, guard_state,
                                         expected, root, now_utc = Sys.time()) {
  request_buffer <- v2r4_capture_file(request_path, "FINAL_CHALLENGE_REQUEST", root)
  response_buffer <- v2r4_capture_file(response_path, "FINAL_CHALLENGE_RESPONSE", root)
  request <- v2r4_read_json_raw(request_buffer, "FINAL_CHALLENGE_REQUEST")
  v2r4_assert(identical(names(request), v2r4_final_challenge_request_fields()) &&
                identical(request$schema_version,
                          cfg$finalization$final_challenge_request_schema_version) &&
                identical(request$status, "ONE_TIME_FINALIZATION_CHALLENGE") &&
                identical(request$runtime_pre_post_equal, TRUE) &&
                identical(request$runtime_post_digest, request$runtime_pre_digest) &&
                identical(request$candidate_gate_eligible, FALSE) &&
                identical(request$can_generate_or_upgrade_candidate, FALSE) &&
                identical(request$formal_sequential_inference, FALSE),
              "FINAL_CHALLENGE_REQUEST_SCHEMA_OR_SEMANTICS_INVALID")
  issued <- v2r4_parse_fixed6_time(request$issued_utc, "FINAL_CHALLENGE_ISSUED")
  expires <- v2r4_parse_fixed6_time(request$expires_utc, "FINAL_CHALLENGE_EXPIRES")
  now_num <- as.numeric(as.POSIXct(now_utc, tz = "UTC"))
  v2r4_assert(as.numeric(issued) <= now_num && now_num < as.numeric(expires) &&
                as.numeric(expires) - as.numeric(issued) <=
                  as.numeric(cfg$finalization$final_challenge_max_age_seconds),
              "FINAL_CHALLENGE_REQUEST_NOT_CURRENT")
  v2r4_assert_named_expected(request, expected, "FINAL_CHALLENGE_REQUEST")
  envelope <- v2r4_read_json_raw(response_buffer, "FINAL_CHALLENGE_RESPONSE")
  v2r4_assert(identical(names(envelope),
                         c("schema_version", "signed_payload_b64", "signature_b64")) &&
                identical(envelope$schema_version,
                          cfg$finalization$final_challenge_response_envelope_schema_version),
              "FINAL_CHALLENGE_RESPONSE_ENVELOPE_INVALID")
  payload_raw <- tryCatch(jsonlite::base64_dec(envelope$signed_payload_b64),
                          error = function(e) raw())
  signature_raw <- tryCatch(jsonlite::base64_dec(envelope$signature_b64),
                            error = function(e) raw())
  public_raw <- tryCatch(jsonlite::base64_dec(guard_state$guard_public_key_pem_b64),
                         error = function(e) raw())
  v2r4_assert(length(payload_raw) > 0L && length(signature_raw) > 0L &&
                length(public_raw) > 0L &&
                identical(v2r4_sha256_raw(public_raw), guard_state$guard_public_key_sha256),
              "FINAL_CHALLENGE_RESPONSE_KEY_OR_BASE64_INVALID")
  key <- tryCatch(openssl::read_pubkey(rawToChar(public_raw)), error = function(e) NULL)
  verified <- !is.null(key) && tryCatch(openssl::signature_verify(
    payload_raw, signature_raw, hash = openssl::sha256, pubkey = key
  ), error = function(e) FALSE)
  v2r4_assert(isTRUE(verified), "FINAL_CHALLENGE_RESPONSE_SIGNATURE_INVALID")
  payload <- tryCatch(jsonlite::fromJSON(rawToChar(payload_raw), simplifyVector = TRUE),
                      error = function(e) NULL)
  v2r4_assert(is.list(payload) && identical(names(payload),
                                            v2r4_final_challenge_response_fields()) &&
                identical(payload$schema_version,
                          cfg$finalization$final_challenge_response_payload_schema_version) &&
                identical(payload$status,
                          "LIVE_ORIGINAL_LAUNCHER_OWNED_GUARD_MUTEX_HELD_FINAL_CHALLENGE_SIGNED") &&
                identical(payload$request_sha256, request_buffer$sha256) &&
                identical(payload$mutex_file_id_current, guard_state$mutex_file_id) &&
                identical(payload$guard_mutex_file_id, guard_state$mutex_file_id) &&
                identical(payload$runtime_pre_post_equal, TRUE) &&
                identical(payload$runtime_post_digest, payload$runtime_pre_digest) &&
                identical(payload$candidate_gate_eligible, FALSE) &&
                identical(payload$can_generate_or_upgrade_candidate, FALSE) &&
                identical(payload$formal_sequential_inference, FALSE),
              "FINAL_CHALLENGE_RESPONSE_SCHEMA_OR_HOLDER_INVALID")
  echo_fields <- setdiff(v2r4_final_challenge_request_fields(),
                         c("schema_version", "status", "issued_utc", "expires_utc"))
  v2r4_assert(all(vapply(echo_fields, function(field) {
    if (is.logical(request[[field]])) identical(payload[[field]], request[[field]]) else
      identical(as.character(payload[[field]]), as.character(request[[field]]))
  }, logical(1))) &&
                identical(payload$challenge_issued_utc, request$issued_utc) &&
                identical(payload$challenge_expires_utc, request$expires_utc),
              "FINAL_CHALLENGE_RESPONSE_ECHO_BINDING_INVALID")
  responded <- v2r4_parse_fixed6_time(payload$responded_utc,
                                      "FINAL_CHALLENGE_RESPONDED")
  response_now_num <- as.numeric(as.POSIXct(Sys.time(), tz = "UTC"))
  v2r4_assert(as.numeric(responded) >= as.numeric(issued) &&
                as.numeric(responded) < as.numeric(expires) &&
                as.numeric(responded) <= response_now_num &&
                response_now_num < as.numeric(expires) &&
                response_now_num - as.numeric(responded) <=
                  as.numeric(cfg$finalization$final_challenge_max_age_seconds),
              "FINAL_CHALLENGE_RESPONSE_NOT_FRESH")
  list(request = request_buffer, response = response_buffer, request_value = request,
       response_value = payload, payload_sha256 = v2r4_sha256_raw(payload_raw))
}

v2r4_validate_control_complete <- function(path, cfg, expected, root) {
  buffer <- v2r4_capture_file(path, "CONTROL_COMPLETE", root)
  value <- v2r4_read_json_raw(buffer, "CONTROL_COMPLETE")
  v2r4_assert(identical(names(value), v2r4_control_complete_fields()) &&
                identical(value$schema_version, cfg$finalization$control_complete_schema_version) &&
                identical(value$status,
                          "FINAL_CONTROL_VERIFIED_VALIDATOR_ALIVE_WAITING_FOR_ACK") &&
                identical(value$dataset_id, cfg$dataset_id) &&
                identical(value$model_version, cfg$model_version) &&
                identical(value$runtime_pre_post_equal, TRUE) &&
                identical(value$runtime_post_digest, value$runtime_pre_digest) &&
                identical(value$e2_gate_status, cfg$gate$machine_status) &&
                identical(value$candidate_gate_eligible, FALSE) &&
                identical(value$can_generate_or_upgrade_candidate, FALSE) &&
                identical(value$formal_sequential_inference, FALSE),
              "CONTROL_COMPLETE_SCHEMA_OR_SEMANTICS_INVALID")
  v2r4_assert_utc_fixed6(c(
    value$launcher_start_utc, value$guard_owner_start_utc,
    value$model_start_utc, value$validator_start_utc, value$guard_start_utc,
    value$final_challenge_issued_utc, value$final_challenge_responded_utc,
    value$lease_issued_utc, value$lease_expires_utc, value$issued_utc
  ), "CONTROL_COMPLETE_TIMES")
  control_challenge_issued <- v2r4_parse_fixed6_time(
    value$final_challenge_issued_utc, "CONTROL_COMPLETE_CHALLENGE_ISSUED"
  )
  control_challenge_responded <- v2r4_parse_fixed6_time(
    value$final_challenge_responded_utc, "CONTROL_COMPLETE_CHALLENGE_RESPONDED"
  )
  control_issued <- v2r4_parse_fixed6_time(
    value$issued_utc, "CONTROL_COMPLETE_ISSUED"
  )
  v2r4_assert(as.numeric(control_challenge_issued) <=
                as.numeric(control_challenge_responded) &&
                as.numeric(control_challenge_responded) <= as.numeric(control_issued),
              "CONTROL_COMPLETE_CHALLENGE_TIME_ORDER_INVALID")
  v2r4_assert_named_expected(value, expected, "CONTROL_COMPLETE")
  list(buffer = buffer, value = value)
}

v2r4_write_consumer_receipt <- function(path, cfg, role, ack, expected, root) {
  v2r4_assert(role %in% c("MODEL", "VALIDATOR"), "CONSUMER_RECEIPT_ROLE_INVALID")
  identity <- list(pid = as.integer(Sys.getpid()), start_utc = if (identical(role, "MODEL")) {
    expected$model_start_utc
  } else expected$validator_start_utc)
  value <- c(list(
    schema_version = cfg$finalization$consumer_receipt_schema_version,
    status = "SIGNED_FINAL_ACK_INDEPENDENTLY_VERIFIED",
    consumer_role = role,
    consumer_pid = identity$pid,
    consumer_start_utc = identity$start_utc,
    final_ack_envelope_sha256 = ack$envelope$sha256,
    final_ack_payload_sha256 = ack$payload_sha256
  ), expected[c("holder_principal", v2r4_final_ack_binding_fields())], list(
    issued_utc = v2r4_now_fixed6()
  ))
  v2r4_assert(identical(names(value), v2r4_consumer_receipt_fields()),
              "CONSUMER_RECEIPT_FIELD_ORDER_INVALID")
  v2r4_write_json_atomic(value, path, root)
}

v2r4_verify_signed_final_ack <- function(envelope_path, public_key_path, expected_key_sha256,
                                         cfg, expected, timeout_seconds) {
  v2r4_assert(identical(as.integer(timeout_seconds), 900L),
              "SIGNED_FINAL_ACK_TIMEOUT_MUST_BE_EXACT_900")
  v2r4_wait_for_files(c(envelope_path, public_key_path), timeout_seconds,
                      "EXTERNAL_SIGNED_FINAL_ACK")
  envelope_buffer <- v2r4_capture_file(envelope_path, "SIGNED_FINAL_ACK_ENVELOPE")
  public_buffer <- v2r4_capture_file(public_key_path, "SIGNED_FINAL_ACK_PUBLIC_KEY")
  v2r4_assert(grepl("^[0-9a-f]{64}$", expected_key_sha256) &&
                identical(public_buffer$sha256, expected_key_sha256),
              "SIGNED_FINAL_ACK_PUBLIC_KEY_PIN_MISMATCH")
  envelope <- v2r4_read_json_raw(envelope_buffer, "SIGNED_FINAL_ACK_ENVELOPE")
  v2r4_assert(identical(names(envelope),
                        c("schema_version", "signed_payload_b64", "signature_b64")) &&
                identical(envelope$schema_version,
                          cfg$finalization$final_ack_envelope_schema_version),
              "SIGNED_FINAL_ACK_ENVELOPE_SCHEMA_INVALID")
  payload_raw <- tryCatch(jsonlite::base64_dec(envelope$signed_payload_b64),
                          error = function(e) raw())
  signature_raw <- tryCatch(jsonlite::base64_dec(envelope$signature_b64),
                            error = function(e) raw())
  key <- tryCatch(openssl::read_pubkey(v2r4_parse_text_raw(public_buffer,
                                                           "SIGNED_FINAL_ACK_PUBLIC_KEY")),
                  error = function(e) NULL)
  verified <- !is.null(key) && length(payload_raw) > 0L && length(signature_raw) > 0L &&
    tryCatch(openssl::signature_verify(payload_raw, signature_raw, hash = openssl::sha256,
                                       pubkey = key), error = function(e) FALSE)
  v2r4_assert(isTRUE(verified), "SIGNED_FINAL_ACK_SIGNATURE_INVALID")
  payload <- tryCatch(jsonlite::fromJSON(rawToChar(payload_raw), simplifyVector = TRUE),
                      error = function(e) NULL)
  required <- c(
    "schema_version", "dataset_id", "model_version", "decision", "purpose", "status",
    "holder_principal", v2r4_final_ack_binding_fields(), "issued_utc"
  )
  v2r4_assert(is.list(payload) && identical(names(payload), required),
              "SIGNED_FINAL_ACK_PAYLOAD_SCHEMA_INVALID")
  v2r4_assert_utc_fixed6(c(
    payload$issued_utc, payload$launcher_start_utc, payload$guard_owner_start_utc,
    payload$guard_start_utc, payload$model_start_utc,
    payload$validator_start_utc, payload$final_challenge_issued_utc,
    payload$final_challenge_responded_utc, payload$lease_issued_utc,
    payload$lease_expires_utc
  ),
                          "SIGNED_FINAL_ACK_TIMES")
  v2r4_assert(identical(names(expected),
                         c("holder_principal", v2r4_final_ack_binding_fields())),
               "SIGNED_FINAL_ACK_EXPECTED_BINDING_SCHEMA_INVALID")
  v2r4_assert(all(vapply(names(expected), function(name) {
    if (is.logical(expected[[name]])) {
      identical(payload[[name]], expected[[name]])
    } else {
      identical(as.character(payload[[name]]), as.character(expected[[name]]))
    }
  }, logical(1))), "SIGNED_FINAL_ACK_EXPECTED_BINDING_MISMATCH")
  v2r4_assert(identical(payload$schema_version,
                        cfg$finalization$final_ack_payload_schema_version) &&
                identical(payload$dataset_id, cfg$dataset_id) &&
                identical(payload$model_version, cfg$model_version) &&
                identical(payload$decision, cfg$finalization$decision) &&
                identical(payload$purpose, cfg$finalization$purpose) &&
                 identical(payload$status, "FINAL_VALIDATION_ATTESTED_NO_CANDIDATE_GATE") &&
                identical(payload$runtime_pre_post_equal, TRUE) &&
                identical(payload$runtime_post_digest, payload$runtime_pre_digest) &&
                identical(payload$named_mutex_unavailable, TRUE) &&
                identical(payload$guard_mutex_held, TRUE) &&
                identical(payload$e2_gate_status, cfg$gate$machine_status) &&
                identical(payload$candidate_gate_eligible, FALSE) &&
                identical(payload$can_generate_or_upgrade_candidate, FALSE) &&
                identical(payload$formal_sequential_inference, FALSE) &&
                all(vapply(c(
                  "lease_manifest_sha256", "source_authority_payload_sha256",
                  "model_source_sha256",
                  "validator_source_sha256", "shared_library_sha256",
                   "analysis_config_sha256", "external_launcher_sha256", "os_guard_sha256",
                   "guard_stop_hmac_key_sha256", "finalization_public_key_sha256",
                   "finalization_public_key_captured_sha256",
                  "guard_nonce", "guard_ready_receipt_sha256", "output_manifest_sha256",
                  "complete_sha256", "validation_report_sha256", "validation_request_sha256",
                  "snapshot_tree_sha256", "runtime_pre_digest", "runtime_post_digest",
                  "runtime_post_receipt_sha256", "model_request_sha256",
                  "model_request_nonce", "model_identity_request_sha256",
                  "validator_request_sha256", "validator_request_nonce",
                  "validator_identity_request_sha256", "science_complete_sha256",
                  "control_complete_sha256", "finalization_request_sha256",
                  "final_challenge_request_sha256", "final_challenge_response_sha256",
                  "final_challenge_payload_sha256", "final_challenge_nonce"
                ), function(field) grepl("^[0-9a-f]{64}$", payload[[field]]), logical(1))),
              "SIGNED_FINAL_ACK_SEMANTIC_INVALID")
  list(envelope = envelope_buffer, public_key = public_buffer, payload = payload,
       payload_sha256 = v2r4_sha256_raw(payload_raw))
}

v2r4_validate_species <- function(cohort, sample_gate, cfg) {
  cohort <- data.table::as.data.table(cohort)
  sample_gate <- data.table::as.data.table(sample_gate)
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  v2r4_assert(nrow(cohort) == 1L && all(c(
    "dataset_id", "n_libraries", "all_sample_level_rules_pass",
    "median_rat_minus_human_rate", "threshold_median_rat_minus_human_rate", "status"
  ) %in% names(cohort)), "SPECIES_COHORT_SCHEMA_INVALID")
  required <- c("dataset_id", "sample_id", "diagnostic_method", "mapping_unit",
                "total_reads_or_pairs",
                "rat_unique_mapping_rate", "human_unique_mapping_rate", "rat_minus_human_rate",
                "rat_to_human_rate_ratio", "threshold_rat_rate", "threshold_rate_difference",
                "threshold_rate_ratio", "status")
  v2r4_assert(nrow(sample_gate) == 9L && all(required %in% names(sample_gate)) &&
                identical(sample_gate$sample_id, ids), "SPECIES_SAMPLE_SCHEMA_OR_ORDER_INVALID")
  finite_fields <- c("rat_unique_mapping_rate", "human_unique_mapping_rate",
                     "rat_minus_human_rate", "threshold_rat_rate",
                     "threshold_rate_difference", "threshold_rate_ratio")
  v2r4_assert(all(vapply(sample_gate[, ..finite_fields], function(value) {
    all(is.finite(value))
  }, logical(1))) &&
                all((is.finite(sample_gate$rat_to_human_rate_ratio) &
                       sample_gate$rat_to_human_rate_ratio >= 0) |
                      (is.infinite(sample_gate$rat_to_human_rate_ratio) &
                         sample_gate$rat_to_human_rate_ratio > 0)),
              "SPECIES_SAMPLE_NUMERIC_DOMAIN_INVALID")
  v2r4_assert(all(sample_gate$rat_unique_mapping_rate >= 0 &
                    sample_gate$rat_unique_mapping_rate <= 1) &&
                all(sample_gate$human_unique_mapping_rate >= 0 &
                      sample_gate$human_unique_mapping_rate <= 1) &&
                all(sample_gate$rat_minus_human_rate >= -1 &
                      sample_gate$rat_minus_human_rate <= 1),
              "SPECIES_SAMPLE_RATE_RANGE_INVALID")
  v2r4_assert(all(sample_gate$dataset_id == cfg$dataset_id) &&
                all(sample_gate$diagnostic_method == cfg$input_snapshot$species_diagnostic_method) &&
                all(sample_gate$mapping_unit == "fragments") &&
                all(sample_gate$total_reads_or_pairs ==
                      cfg$input_snapshot$species_fragments_per_library),
              "SPECIES_SAMPLE_IDENTITY_OR_METHOD_INVALID")
  v2r4_assert(
    all(sample_gate$threshold_rat_rate ==
          cfg$input_snapshot$species_minimum_rat_unique_mapping_rate) &&
      all(sample_gate$threshold_rate_difference ==
            cfg$input_snapshot$species_minimum_rat_minus_human_rate) &&
      all(sample_gate$threshold_rate_ratio ==
            cfg$input_snapshot$species_minimum_rat_to_human_rate_ratio),
    "SPECIES_SAMPLE_FROZEN_THRESHOLD_DRIFT"
  )
  recomputed_difference <- sample_gate$rat_unique_mapping_rate -
    sample_gate$human_unique_mapping_rate
  recomputed_ratio <- ifelse(
    sample_gate$human_unique_mapping_rate == 0,
    Inf,
    sample_gate$rat_unique_mapping_rate / sample_gate$human_unique_mapping_rate
  )
  v2r4_assert(all(abs(sample_gate$rat_minus_human_rate - recomputed_difference) <= 1e-12),
              "SPECIES_SAMPLE_DIFFERENCE_RECOMPUTATION_MISMATCH")
  v2r4_assert(all((is.finite(recomputed_ratio) &
                     abs(sample_gate$rat_to_human_rate_ratio - recomputed_ratio) <= 1e-12) |
                    (is.infinite(recomputed_ratio) & recomputed_ratio > 0 &
                       is.infinite(sample_gate$rat_to_human_rate_ratio) &
                       sample_gate$rat_to_human_rate_ratio > 0)),
              "SPECIES_SAMPLE_RATIO_RECOMPUTATION_MISMATCH")
  recomputed_sample_pass <-
    sample_gate$rat_unique_mapping_rate >=
      cfg$input_snapshot$species_minimum_rat_unique_mapping_rate &
    sample_gate$human_unique_mapping_rate <=
      cfg$input_snapshot$species_maximum_human_unique_mapping_rate &
    recomputed_difference >= cfg$input_snapshot$species_minimum_rat_minus_human_rate &
    recomputed_ratio >= cfg$input_snapshot$species_minimum_rat_to_human_rate_ratio
  v2r4_assert(all(recomputed_sample_pass) && all(sample_gate$status == "PASS"),
              "SPECIES_SAMPLE_RECOMPUTED_RULE_FAILURE")
  observed <- stats::median(recomputed_difference)
  recomputed_cohort_pass <- all(recomputed_sample_pass) &&
    observed >= cfg$input_snapshot$species_minimum_cohort_median_rate_difference
  v2r4_assert(cohort$dataset_id[[1L]] == cfg$dataset_id &&
                as.integer(cohort$n_libraries[[1L]]) == 9L &&
                isTRUE(cohort$all_sample_level_rules_pass[[1L]]) &&
                identical(cohort$status[[1L]], "PASS") && recomputed_cohort_pass,
              "SPECIES_COHORT_NOT_RECOMPUTED_PASS")
  v2r4_assert(abs(observed - cohort$median_rat_minus_human_rate[[1L]]) <= 1e-12 &&
                cohort$threshold_median_rat_minus_human_rate[[1L]] ==
                  cfg$input_snapshot$species_minimum_cohort_median_rate_difference,
              "SPECIES_COHORT_MEDIAN_OR_THRESHOLD_MISMATCH")
  invisible(TRUE)
}

v2r4_validate_reference_manifest <- function(x, cfg) {
  x <- data.table::as.data.table(x)
  required <- c("relative_path", "role", "bytes", "md5", "sha256")
  v2r4_assert(all(required %in% names(x)), "REFERENCE_MANIFEST_SCHEMA_INVALID")
  fasta <- x[base::which(x[["role"]] == "reference_fasta")]
  gtf <- x[base::which(x[["role"]] == "gene_annotation_gtf")]
  v2r4_assert(nrow(fasta) == 1L && nrow(gtf) == 1L, "REFERENCE_FASTA_GTF_ROWS_INVALID")
  v2r4_assert(fasta$md5[[1L]] == cfg$input_snapshot$reference_fasta_md5 &&
                fasta$sha256[[1L]] == cfg$input_snapshot$reference_fasta_sha256 &&
                gtf$md5[[1L]] == cfg$input_snapshot$gtf_md5 &&
                gtf$sha256[[1L]] == cfg$input_snapshot$gtf_sha256,
              "REFERENCE_IDENTITY_MISMATCH")
  invisible(TRUE)
}

v2r4_validate_strand <- function(decision, assignment, cfg) {
  decision <- data.table::as.data.table(decision)
  assignment <- data.table::as.data.table(assignment)
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  v2r4_assert(nrow(decision) == 1L && nrow(assignment) == 9L &&
                identical(assignment$sample_id, ids), "STRAND_ROW_COUNT_OR_ORDER_INVALID")
  v2r4_assert(all(c("assigned_0", "assigned_1", "assigned_2") %in% names(assignment)),
              "STRAND_ASSIGNMENT_SCHEMA_INVALID")
  dominant_mode <- ifelse(assignment$assigned_1 >= assignment$assigned_2, 1L, 2L)
  dominant <- ifelse(dominant_mode == 1L, assignment$assigned_1, assignment$assigned_2)
  minor <- ifelse(dominant_mode == 1L, assignment$assigned_2, assignment$assigned_1)
  candidate <- as.integer(names(which.max(table(dominant_mode))))
  same_fraction <- mean(dominant_mode == candidate)
  median_dominant_unstranded <- stats::median(dominant / assignment$assigned_0)
  median_minor_dominant <- stats::median(minor / dominant)
  selected <- if (same_fraction >= 0.8 && median_dominant_unstranded >= 0.8 &&
                  median_minor_dominant <= 0.2) candidate else 0L
  v2r4_assert(selected == cfg$input_snapshot$selected_primary_strand_mode &&
                decision$selected_primary_strand_mode[[1L]] == selected &&
                decision$decision[[1L]] == cfg$input_snapshot$strand_decision &&
                !isTRUE(decision$labels_used_for_decision[[1L]]) && decision$status[[1L]] == "PASS",
              "STRAND_DECISION_RECOMPUTATION_MISMATCH")
  v2r4_assert(abs(decision$dominant_mode_fraction[[1L]] - same_fraction) < 1e-12 &&
                abs(decision$median_dominant_over_unstranded[[1L]] - median_dominant_unstranded) < 1e-12 &&
                abs(decision$median_minor_over_dominant[[1L]] - median_minor_dominant) < 1e-12,
              "STRAND_DECISION_METRIC_MISMATCH")
  invisible(TRUE)
}

v2r4_validate_upstream_gate <- function(x, cfg) {
  x <- data.table::as.data.table(x)
  v2r4_assert(all(c("gate", "eligible_input_n", "status") %in% names(x)),
              "UPSTREAM_GATE_SCHEMA_INVALID")
  e1 <- x[base::which(x[["gate"]] == "E1")]
  v2r4_assert(nrow(e1) == 1L && as.integer(e1$eligible_input_n[[1L]]) == 0L &&
                e1$status[[1L]] == "NOT_IDENTIFIABLE_NO_S2_OR_J2_INPUT",
              "UPSTREAM_E1_NOT_EXACTLY_EMPTY")
  invisible(TRUE)
}

v2r4_validate_ledger <- function(x, cfg) {
  x <- data.table::as.data.table(x)
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  groups <- unlist(cfg$design$groups_in_order, use.names = FALSE)
  required <- c("dataset_id", "sample_order", "sample_id", "biological_unit_id",
                "statistical_unit", "organism", "experimental_arm", "condition",
                "replicate", "library_layout", "biosample", "sra_experiment",
                "run_accession", "include_primary", "metadata_status")
  v2r4_assert(nrow(x) == 9L && all(required %in% names(x)), "LEDGER_SCHEMA_OR_N_INVALID")
  v2r4_assert(identical(x$sample_id, ids) && identical(as.integer(x$sample_order), 1:9),
              "LEDGER_SAMPLE_ORDER_INVALID")
  observed_groups <- ifelse(x$experimental_arm == "JA", "sham", x$experimental_arm)
  v2r4_assert(identical(observed_groups, groups), "LEDGER_GROUP_MAPPING_INVALID")
  v2r4_assert(!anyDuplicated(x$sample_id) && !anyDuplicated(x$biological_unit_id) &&
                !anyDuplicated(x$biosample) && !anyDuplicated(x$sra_experiment) &&
                !anyDuplicated(x$run_accession) && identical(x$sample_id, x$biological_unit_id),
              "LEDGER_UNIT_NOT_ONE_TO_ONE")
  v2r4_assert(all(x$statistical_unit == "independent_rat_library") &&
                all(x$organism == "Rattus norvegicus") && all(x$library_layout == "PAIRED"),
              "LEDGER_UNIT_SPECIES_LAYOUT_INVALID")
  v2r4_assert(identical(as.character(x$include_primary),
                        rep(cfg$input_snapshot$required_sample_include_primary_marker, 9L)) &&
                identical(as.character(x$metadata_status),
                          rep(cfg$input_snapshot$required_sample_metadata_status, 9L)),
              "LEDGER_CONDITIONAL_MARKER_OR_METADATA_STATUS_INVALID")
  x
}

v2r4_release_conditional_inclusion <- function(ledger, cohort, sample_gate,
                                                cohort_sha256, sample_gate_sha256,
                                                cfg) {
  ledger <- data.table::as.data.table(ledger)
  cohort <- data.table::as.data.table(cohort)
  sample_gate <- data.table::as.data.table(sample_gate)
  v2r4_assert(grepl("^[0-9a-f]{64}$", cohort_sha256) &&
                grepl("^[0-9a-f]{64}$", sample_gate_sha256),
              "CONDITIONAL_INCLUSION_SPECIES_EVIDENCE_SHA_INVALID")
  v2r4_validate_species(cohort, sample_gate, cfg)
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  v2r4_assert(nrow(ledger) == 9L && identical(ledger$sample_id, ids) &&
                identical(as.character(ledger$include_primary), rep("CONDITIONAL", 9L)) &&
                identical(sample_gate$sample_id, ids),
              "CONDITIONAL_INCLUSION_LEDGER_OR_SPECIES_ORDER_INVALID")
  list(
    sample_include_primary_marker = cfg$input_snapshot$required_sample_include_primary_marker,
    sample_metadata_status = cfg$input_snapshot$required_sample_metadata_status,
    conditional_inclusion_release_status =
      cfg$input_snapshot$conditional_inclusion_release_status,
    conditional_inclusion_n = as.integer(nrow(ledger)),
    cohort_species_gate_sha256 = cohort_sha256,
    sample_species_table_sha256 = sample_gate_sha256,
    plm_metadata_conflict_policy = cfg$input_snapshot$plm_metadata_conflict_policy
  )
}

v2r4_design <- function(samples, cfg) {
  ids <- unlist(cfg$design$sample_ids_in_order, use.names = FALSE)
  groups <- factor(unlist(cfg$design$groups_in_order, use.names = FALSE),
                   levels = unlist(cfg$design$group_levels, use.names = FALSE))
  names(groups) <- ids
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- unlist(cfg$design$coefficient_order, use.names = FALSE)
  rownames(design) <- ids
  storage.mode(design) <- "double"
  v2r4_assert(identical(dim(design), c(9L, 3L)) && qr(design)$rank == cfg$design$expected_rank,
              "DESIGN_DIM_OR_RANK_INVALID")
  contrast <- cbind(
    local_disease = unlist(cfg$design$contrasts$local_disease$vector, use.names = FALSE),
    intervention = unlist(cfg$design$contrasts$intervention$vector, use.names = FALSE)
  )
  rownames(contrast) <- colnames(design)
  storage.mode(contrast) <- "double"
  v2r4_assert(identical(as.vector(contrast[, "local_disease"]), c(-1, 1, 0)) &&
                identical(as.vector(contrast[, "intervention"]), c(0, -1, 1)),
              "CONTRAST_DIRECTION_INVALID")
  list(group = groups, design = design, contrast = contrast)
}

v2r4_validate_counts_original <- function(counts, expected_sample_ids) {
  v2r4_assert(is.matrix(counts) && is.numeric(counts), "ASSEMBLED_COUNTS_NOT_NUMERIC_MATRIX")
  v2r4_assert(ncol(counts) == length(expected_sample_ids) &&
                identical(colnames(counts), expected_sample_ids),
              "ASSEMBLED_COUNTS_SAMPLE_ORDER_INVALID")
  v2r4_assert(all(is.finite(counts)), "ASSEMBLED_COUNTS_NONFINITE")
  v2r4_assert(all(counts >= 0), "ASSEMBLED_COUNTS_NEGATIVE")
  v2r4_assert(all(counts == floor(counts)), "ASSEMBLED_COUNTS_FRACTIONAL")
  v2r4_assert(all(counts <= .Machine$integer.max), "ASSEMBLED_COUNTS_INTEGER_OVERFLOW")
  out <- counts
  storage.mode(out) <- "integer"
  out
}

v2r4_validate_annotation_gene <- function(annotation, gene_ids, minimum_coverage = 0.80) {
  annotation <- data.table::as.data.table(annotation)
  v2r4_assert(all(c("GeneID", "gene", "gene_biotype") %in% names(annotation)),
              "GENE_ANNOTATION_FROZEN_FEATURECOUNTS_COLUMNS_MISSING")
  v2r4_assert(!("gene_name" %in% names(annotation)),
              "GENE_ANNOTATION_NONCONTRACT_GENE_NAME_PRESENT")
  v2r4_assert(nrow(annotation) == length(gene_ids), "GENE_ANNOTATION_ROW_COUNT_INVALID")
  v2r4_assert(identical(as.character(annotation[["GeneID"]]), as.character(gene_ids)) &&
                !anyDuplicated(annotation[["GeneID"]]),
              "GENE_ANNOTATION_GENEID_ORDER_OR_IDENTITY_INVALID")
  symbols <- trimws(as.character(annotation[["gene"]]))
  nonblank <- !is.na(symbols) & nzchar(symbols)
  coverage <- mean(nonblank)
  v2r4_assert(length(minimum_coverage) == 1L && is.finite(minimum_coverage) &&
                minimum_coverage >= 0 && minimum_coverage <= 1,
              "GENE_ANNOTATION_MINIMUM_COVERAGE_INVALID")
  v2r4_assert(is.finite(coverage) && coverage >= minimum_coverage,
              "GENE_ANNOTATION_SYMBOL_COVERAGE_TOO_LOW")
  fallback <- sum(!nonblank)
  symbols[!nonblank] <- gene_ids[!nonblank]
  list(
    annotation = data.table::data.table(gene_id = gene_ids, symbol = symbols),
    metrics = data.table::data.table(
      gene_id_source_column = "GeneID", symbol_source_column = "gene",
      biotype_source_column = "gene_biotype", total_gene_rows = length(gene_ids),
      nonblank_symbol_n = sum(nonblank), nonblank_symbol_coverage = coverage,
      unique_nonblank_symbol_n = data.table::uniqueN(symbols[nonblank]),
      unique_nonblank_symbol_coverage = data.table::uniqueN(symbols[nonblank]) / length(gene_ids),
      fallback_to_gene_id_n = fallback
    )
  )
}

v2r4_edgeR_ql_covariance <- function(fit, contrast) {
  contrast <- as.numeric(contrast)
  n_genes <- nrow(fit$coefficients)
  empty <- list(se = rep(NA_real_, n_genes), df = rep(NA_real_, n_genes),
                status = rep("QL_COVARIANCE_INPUT_UNAVAILABLE", n_genes))
  required <- c("design", "fitted.values", "dispersion", "s2.post")
  if (!all(required %in% names(fit)) || length(contrast) != ncol(fit$design)) return(empty)
  design <- as.matrix(fit$design)
  mu <- as.matrix(fit$fitted.values)
  if (!identical(dim(mu), c(n_genes, nrow(design)))) return(empty)
  dispersion <- fit$dispersion
  prior_weights <- fit$weights
  df_resid <- fit$df.residual.zeros %||% fit$df.residual.adj %||% fit$df.residual
  df_total_cap <- sum(fit$df.residual)
  out <- empty
  for (i in seq_len(n_genes)) {
    phi <- if (is.matrix(dispersion)) dispersion[i, ] else if (length(dispersion) == 1L) {
      rep(dispersion, nrow(design))
    } else rep(dispersion[[i]], nrow(design))
    pweights <- if (is.null(prior_weights)) rep(1, nrow(design)) else if (is.matrix(prior_weights)) {
      prior_weights[i, ]
    } else if (length(prior_weights) == 1L) rep(prior_weights, nrow(design)) else prior_weights
    qscale <- fit$s2.post[[i]]
    if (!all(is.finite(c(mu[i, ], phi, pweights, qscale))) || qscale <= 0 ||
        any(mu[i, ] < 0) || any(phi < 0) || any(pweights < 0)) {
      out$status[[i]] <- "QL_COVARIANCE_NONFINITE_OR_NONPOSITIVE_INPUT"
      next
    }
    fisher_w <- pweights * mu[i, ] / (1 + phi * mu[i, ])
    information <- crossprod(design, design * fisher_w)
    if (qr(information)$rank < ncol(design)) {
      out$status[[i]] <- "QL_COVARIANCE_INFORMATION_SINGULAR"
      next
    }
    covariance <- tryCatch(solve(information), error = function(e) NULL)
    if (is.null(covariance)) {
      out$status[[i]] <- "QL_COVARIANCE_INFORMATION_NOT_INVERTIBLE"
      next
    }
    variance_log2 <- as.numeric(crossprod(contrast, covariance %*% contrast)) * qscale / log(2)^2
    if (!is.finite(variance_log2) || variance_log2 <= 0) {
      out$status[[i]] <- "QL_COVARIANCE_VARIANCE_NOT_POSITIVE_FINITE"
      next
    }
    df_i <- min(as.numeric(fit$df.prior[[min(i, length(fit$df.prior))]]) +
                  as.numeric(df_resid[[min(i, length(df_resid))]]), df_total_cap)
    if (!is.finite(df_i) || df_i <= 0) {
      out$status[[i]] <- "QL_COVARIANCE_DF_NOT_POSITIVE_FINITE"
      next
    }
    out$se[[i]] <- sqrt(variance_log2)
    out$df[[i]] <- df_i
    out$status[[i]] <- "ESTIMABLE_QL_FISHER_COVARIANCE"
  }
  out
}

v2r4_recompute_bh <- function(x, collection_col = "collection", contrast_col = "contrast",
                             p_col = "p_value", q_col = "q_value", tested_col = "tested") {
  x <- data.table::copy(data.table::as.data.table(x))
  x[, expected_q__ := NA_real_]
  eligible <- base::which(as.logical(x[[tested_col]]) & is.finite(x[[p_col]]))
  if (length(eligible)) {
    parts <- split(eligible, interaction(x[[collection_col]][eligible], x[[contrast_col]][eligible],
                                         drop = TRUE, lex.order = TRUE))
    for (idx in parts) x$expected_q__[idx] <- stats::p.adjust(x[[p_col]][idx], method = "BH")
  }
  observed <- x[[q_col]]
  same <- (is.na(observed) & is.na(x$expected_q__)) |
    (is.finite(observed) & is.finite(x$expected_q__) & abs(observed - x$expected_q__) <= 1e-12)
  v2r4_assert(all(same), "COLLECTION_CONTRAST_BH_RECOMPUTATION_MISMATCH")
  x[, expected_q__ := NULL]
  invisible(TRUE)
}

v2r4_assert_exact_program_grid <- function(x, tested_families, collections, contrasts, label) {
  x <- data.table::as.data.table(x)
  required <- c("collection", "program_id", "contrast")
  v2r4_assert(all(required %in% names(x)), "%s_PROGRAM_GRID_SCHEMA_INVALID", label)
  collection_values <- as.character(collections)
  contrast_values <- as.character(contrasts)
  v2r4_assert(setequal(names(tested_families), collection_values),
              "%s_EXPECTED_FAMILY_COLLECTION_SET_INVALID", label)
  expected <- data.table::rbindlist(lapply(collection_values, function(collection_value) {
    programs <- as.character(tested_families[[collection_value]])
    v2r4_assert(length(programs) > 0L && !anyDuplicated(programs),
                "%s_EXPECTED_PROGRAM_FAMILY_INVALID:%s", label, collection_value)
    data.table::CJ(collection = collection_value, program_id = programs,
                   contrast = contrast_values, unique = TRUE)
  }))
  observed <- data.table::copy(x[, .(
    collection = as.character(collection), program_id = as.character(program_id),
    contrast = as.character(contrast)
  )])
  data.table::setorder(observed, collection, program_id, contrast)
  data.table::setorder(expected, collection, program_id, contrast)
  v2r4_assert(identical(observed, expected), "%s_TESTED_FAMILY_OR_CONTRAST_SET_INVALID", label)
  invisible(TRUE)
}

v2r4_recompute_loao_summary <- function(long, expected_ids) {
  long <- data.table::as.data.table(long)
  keys <- c("collection", "program_id", "contrast")
  required <- c(keys, "omitted_animal", "estimable", "estimate")
  v2r4_assert(all(required %in% names(long)), "LOAO_LONG_SCHEMA_INVALID")
  per_key <- long[, .(
    n_omissions = .N,
    n_unique_omissions = data.table::uniqueN(omitted_animal),
    omission_set = paste(sort(unique(omitted_animal)), collapse = ","),
    n_estimable = sum(as.logical(estimable) & is.finite(estimate)),
    estimate_min = if (any(as.logical(estimable) & is.finite(estimate))) min(estimate[as.logical(estimable) & is.finite(estimate)]) else NA_real_,
    estimate_max = if (any(as.logical(estimable) & is.finite(estimate))) max(estimate[as.logical(estimable) & is.finite(estimate)]) else NA_real_,
    estimate_mean = if (any(as.logical(estimable) & is.finite(estimate))) mean(estimate[as.logical(estimable) & is.finite(estimate)]) else NA_real_,
    sign_consistency = if (any(as.logical(estimable) & is.finite(estimate))) {
      vals <- estimate[as.logical(estimable) & is.finite(estimate)]
      max(mean(vals > 0), mean(vals < 0))
    } else NA_real_
  ), by = keys]
  expected_set <- paste(sort(expected_ids), collapse = ",")
  v2r4_assert(all(per_key$n_omissions == length(expected_ids)) &&
                all(per_key$n_unique_omissions == length(expected_ids)) &&
                all(per_key$omission_set == expected_set), "LOAO_NOT_EXACT_NINE_ANIMALS_PER_KEY")
  per_key[, c("n_unique_omissions", "omission_set") := NULL]
  data.table::setorderv(per_key, keys)
  per_key
}

v2r4_output_manifest <- function(root, required_payloads) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  rel <- substring(gsub("\\\\", "/", files), nchar(gsub("\\\\", "/", root)) + 2L)
  keep <- !rel %in% c("OUTPUT_MANIFEST.tsv", "COMPLETE.tsv")
  files <- files[keep]
  rel <- rel[keep]
  v2r4_assert(setequal(rel, required_payloads), "OUTPUT_PAYLOAD_SET_MISMATCH")
  ord <- order(rel)
  identities <- lapply(files[ord], v2r4_capture_file, label = "OUTPUT_PAYLOAD_FOR_MANIFEST")
  data.table::data.table(
    relative_path = rel[ord],
    bytes = vapply(identities, `[[`, numeric(1), "bytes"),
    sha256 = vapply(identities, `[[`, character(1), "sha256")
  )
}

v2r4_parse_output_authority <- function(manifest_buffer, cfg, expected_decision,
                                         expected_mode, expected_run_id) {
  auth <- v2r4_read_tsv_raw(manifest_buffer, "OUTPUT_AUTHORITY")
  required <- unlist(cfg$output_authority$required_columns, use.names = FALSE)
  v2r4_assert(identical(names(auth), required) && nrow(auth) == 2L &&
                !anyDuplicated(auth$role) && !anyDuplicated(auth$relative_path),
              "OUTPUT_AUTHORITY_SCHEMA_OR_UNIQUENESS_INVALID")
  v2r4_assert(all(auth$schema_version == cfg$output_authority$schema_version) &&
                all(auth$dataset_id == cfg$dataset_id) &&
                all(auth$model_version == cfg$model_version) &&
                all(auth$mode == expected_mode) &&
                all(auth$run_id == expected_run_id) &&
                all(auth$decision == expected_decision) &&
                all(auth$materialization == "OUTPUT_FILE"), "OUTPUT_AUTHORITY_HEADER_INVALID")
  expected <- data.table::data.table(
    role = c("output_manifest", "complete"),
    relative_path = c("OUTPUT_MANIFEST.tsv", "COMPLETE.tsv")
  )
  got <- auth[match(expected$role, auth$role), .(role, relative_path)]
  v2r4_assert(identical(got, expected), "OUTPUT_AUTHORITY_MARKER_PATHS_NOT_EXACT")
  v2r4_assert(all(grepl("^[0-9a-f]{64}$", auth$sha256)) &&
                all(is.finite(auth$bytes) & auth$bytes >= 0 & auth$bytes == floor(auth$bytes)),
              "OUTPUT_AUTHORITY_IDENTITY_DOMAIN_INVALID")
  auth
}

v2r4_capture_authorized_outputs <- function(output_authority, output_root, required_payloads) {
  root <- v2r4_norm_existing(output_root)
  marker_buffers <- lapply(seq_len(nrow(output_authority)), function(i) {
    row <- output_authority[i]
    buffer <- v2r4_capture_file(file.path(root, row$relative_path[[1L]]),
                                paste0("AUTHORIZED_OUTPUT_", row$role[[1L]]), root)
    v2r4_assert(buffer$bytes == as.numeric(row$bytes[[1L]]) &&
                  buffer$sha256 == row$sha256[[1L]],
                "AUTHORIZED_OUTPUT_MARKER_IDENTITY_MISMATCH:%s", row$role[[1L]])
    buffer
  })
  names(marker_buffers) <- output_authority$role
  manifest <- v2r4_read_tsv_raw(marker_buffers$output_manifest, "AUTHORIZED_OUTPUT_MANIFEST")
  v2r4_assert(identical(names(manifest), c("relative_path", "bytes", "sha256")) &&
                !anyDuplicated(manifest$relative_path) && setequal(manifest$relative_path, required_payloads),
              "OUTPUT_MANIFEST_SCHEMA_OR_PAYLOAD_SET_INVALID")
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  rel <- substring(gsub("\\\\", "/", files), nchar(gsub("\\\\", "/", root)) + 2L)
  v2r4_assert(setequal(rel, c(required_payloads, "OUTPUT_MANIFEST.tsv", "COMPLETE.tsv")),
              "OUTPUT_DIRECTORY_FILE_SET_INVALID")
  payloads <- lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i]
    buffer <- v2r4_capture_file(file.path(root, row$relative_path[[1L]]),
                                paste0("OUTPUT_PAYLOAD_", row$relative_path[[1L]]), root)
    v2r4_assert(buffer$bytes == as.numeric(row$bytes[[1L]]) &&
                  buffer$sha256 == row$sha256[[1L]],
                "OUTPUT_PAYLOAD_IDENTITY_MISMATCH:%s", row$relative_path[[1L]])
    buffer
  })
  names(payloads) <- manifest$relative_path
  list(markers = marker_buffers, manifest = manifest, payloads = payloads)
}

v2r4_validate_uniform_flags <- function(x, cfg, label) {
  required <- c("analysis_role", "candidate_gate_eligible", "can_generate_or_upgrade_candidate",
                "formal_sequential_inference")
  v2r4_assert(all(required %in% names(x)), "%s_FLAG_COLUMNS_MISSING", label)
  v2r4_assert(all(x$analysis_role == cfg$gate$analysis_role) &&
                all(x$candidate_gate_eligible %in% FALSE) &&
                all(x$can_generate_or_upgrade_candidate %in% FALSE) &&
                all(x$formal_sequential_inference %in% FALSE),
              "%s_CANDIDATE_FLAGS_INVALID", label)
  invisible(TRUE)
}
