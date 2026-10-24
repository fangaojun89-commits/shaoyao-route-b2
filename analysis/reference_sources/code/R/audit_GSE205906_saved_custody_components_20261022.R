# Saved-only component audit. No iterative fit or normalization dispatcher.
# Usage: Rscript --vanilla THIS.R PROJECT MODULE_DIRECTORY PARENT_SHA FRESH_OUTPUT
# No arguments: plan only. The genuine successful parent must already exist.
main <- function(args) {
  if (!length(args)) {
    cat("PLAN_ONLY_SAVED_COMPONENT_AUDIT_NO_READS_NO_FITS\n")
    return(invisible(NULL))
  }
  stopifnot(length(args) == 4L)
  project <- args[[1L]]; modules <- args[[2L]]; parent_sha <- args[[3L]]; output <- args[[4L]]
  root <- "/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01"
  run <- file.path(root, "downstream_custody_01")
  parent <- file.path(root, "provenance/matched_downstream_custody_parent_20261022_r01/PARENT_NATIVE_EXIT.json")
  stopifnot(dirname(output) == root, startsWith(basename(output), "saved_"), !file.exists(output))
  pin <- function(path) {
    info <- file.info(path)
    stopifnot(nrow(info) == 1L, !is.na(info$size), !info$isdir, identical(Sys.readlink(path), ""))
    list(path = path, bytes = unname(info$size),
      sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }
  need_pin <- function(value) {
    actual <- pin(value$path)
    stopifnot(identical(value$path, actual$path), value$bytes == actual$bytes,
      identical(value$sha256, actual$sha256))
    invisible(value)
  }
  json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
  self_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  stopifnot(length(self_file) == 1L)
  self_pin <- pin(self_file)
  parent_pin <- pin(parent)
  stopifnot(identical(parent_pin$sha256, parent_sha))
  terminal <- json(parent)
  stopifnot(identical(terminal$schema, "GSE205906_MATCHED_PARENT_PROCESS_EXIT/v1"),
    is.numeric(terminal$actual_R_exit), length(terminal$actual_R_exit) == 1L,
    is.finite(terminal$actual_R_exit), terminal$actual_R_exit == 0L, identical(terminal$run_root, run),
    identical(terminal$caller_source_sha256, "7cdafad1f838cd01ca8ae730762c059a78b2b3f515d47785fd09fe87f9e2275e"),
    is.null(terminal$observer_failure_class_only))
  need_pin(terminal$complete_pin)
  need_pin(terminal$started)
  stopifnot(identical(terminal$started$path, file.path(dirname(parent), "STARTED.json")))
  started <- json(terminal$started$path)
  stopifnot(identical(started$schema, "GSE205906_MATCHED_ACTUAL_R_STARTED/v1"),
    identical(started$actual_R_pid, terminal$actual_R_pid),
    identical(started$actual_parent_pid, terminal$actual_parent_pid),
    identical(started$start_ticks, terminal$start_ticks), identical(started$argv, terminal$argv),
    identical(started$source$sha256, terminal$caller_source_sha256))
  need_pin(started$source)
  stopifnot(identical(terminal$complete_pin$path, file.path(run, "COMPLETE.json")))
  completed <- json(terminal$complete_pin$path)
  stopifnot(identical(completed$status, "COMPLETE_SAVED_NATIVE_NOT_INDEPENDENTLY_ACCEPTED"))
  invisible(lapply(completed$products, need_pin))
  checkpoint_path <- file.path(run, "SAVED_NATIVE_CHECKPOINT.rds")
  checkpoint_pin <- pin(checkpoint_path)
  stopifnot(any(vapply(completed$products, function(value)
    identical(value$path, checkpoint_pin$path) && value$bytes == checkpoint_pin$bytes &&
      identical(value$sha256, checkpoint_pin$sha256), logical(1))))
  caller_path <- file.path(project, "code/R/run_gse205906_matched_downstream_custody_v1_20261022.R")
  stopifnot(identical(pin(caller_path)$sha256, terminal$caller_source_sha256))
  caller <- new.env(parent = baseenv()); sys.source(caller_path, caller)
  constants <- caller$g205run_constants()
  resources_pin <- pin(constants$gene_resources); orthology_pin <- pin(constants$orthology)
  stopifnot(identical(resources_pin$sha256, "0a8526c5f84435a025ec4bdb07be021ce050eda39777f30cbf40bb1ec2daba51"),
    identical(orthology_pin$sha256, "e9dd0c9268e0044f615085e27292cbf977f5a648003ba90054cf3d01225cf998"))
  adapter <- caller$g205run_read(file.path(project, "code/R/gse205906_resource_object_validation_v1.R"), constants$adapter_sha)
  core_path <- file.path(project, "code/R/gse205906_matched_downstream_core_custody_v1_20261022.R")
  core_pin <- pin(core_path)
  stopifnot(identical(core_pin$sha256, constants$core_sha))
  core <- new.env(parent = baseenv()); sys.source(core_path, core)
  pure <- new.env(parent = baseenv())
  core$g205m_import(adapter$raw, "gse205906_resource_label_adapter_v1", pure)
  resources <- pure$gse205906_resource_label_adapter_v1(readRDS(constants$gene_resources))
  orthology <- utils::read.delim(constants$orthology, sep = "\t", quote = '"', colClasses = "character",
    na.strings = "NA", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  checkpoint <- readRDS(checkpoint_path)
  specs <- c(native_numeric = "5a82c2cdccf93dbcd5916eead2c4cf6e4cb0566a094e8c35ad0b96ea549826da",
    normalization = "c785b505f9a2cfa71653bf6267afbe51c3bdaf3c8a7dca59c2a5314db5690f9b",
    unsigned_mapping = "628a67ff470240bd245d8d05f290cff51e5994dc1d040d9763043e555075da9c",
    activity_score = "ad0d6ea1f4ab0c97ea5c6effbe9a1e92d4294c4fa1c34dfd9f23f35bae966d1e")
  loaded <- lapply(names(specs), function(role) {
    path <- file.path(modules, paste0("gse205906_matched_", role, "_readback_20261019.R"))
    stopifnot(identical(pin(path)$sha256, unname(specs[[role]])))
    env <- new.env(parent = baseenv()); sys.source(path, env)
    list(env = env, pin = pin(path))
  })
  rng_pin <- function() if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    digest::digest(get(".Random.seed", .GlobalEnv), "sha256", serializeVersion = 3L) else NULL
  # Saved DESeq2 numeric-contrast extraction initializes an absent seed through
  # its native entry. With an existing seed it consumes no random draws.
  # Establish an audit-only seed when absent, require exact preservation during
  # all components, then restore the caller's original absent/present state.
  rng_existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  original_rng <- if (rng_existed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (rng_existed) assign(".Random.seed", original_rng, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!rng_existed) set.seed(1L)
  before <- rng_pin(); phase <- "SAVED_COMPONENT_READBACK"
  stopifnot(dir.create(output, recursive = FALSE))
  tryCatch({
    results <- list(
      native = loaded[[1L]]$env$g205n_readback(checkpoint),
      normalization = loaded[[2L]]$env$g205t_readback(checkpoint),
      unsigned = loaded[[3L]]$env$g205u_readback(checkpoint, resources, orthology),
      activity = loaded[[4L]]$env$g205a_readback(checkpoint, resources, orthology))
    stopifnot(identical(before, rng_pin()))
    invisible(lapply(c(list(self_pin, parent_pin, terminal$complete_pin, terminal$started,
      started$source, core_pin, adapter$pin, checkpoint_pin, resources_pin, orthology_pin),
      lapply(loaded, `[[`, "pin")), need_pin))
    receipt <- list(status = "FOUR_SAVED_COMPONENTS_VERIFIED_NOT_FULL_CHAIN_ACCEPTANCE",
      audit_source = self_pin, parent = parent_pin, checkpoint = checkpoint_pin, components = results,
      modules = lapply(loaded, `[[`, "pin"), iterative_fits = 0L, RNG_unchanged = TRUE,
      RNG_audit_seed_initialized_if_absent = !rng_existed,
      RNG_original_state_restored_on_exit = TRUE,
      full_independent_numeric_acceptance = FALSE, candidate_admission = FALSE,
      limitations = c("Unsigned P values not replayed", "Voom weights and dispersion not independently reconstructed",
        "Exact native consumption of reconstructed members and scores not authenticated"))
    jsonlite::write_json(receipt, file.path(output, "COMPONENTS_VERIFIED.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  }, error = function(condition) {
    jsonlite::write_json(list(status = "FAILED_COMPONENT_AUDIT_NO_REFIT", phase = phase,
      exception_class = class(condition)[[1L]], condition_sha256 = digest::digest(conditionMessage(condition), "sha256", serialize = FALSE),
      raw_error_saved = FALSE), file.path(output, "FAILED.json"), auto_unbox = TRUE, pretty = TRUE)
    stop("SAVED_COMPONENT_AUDIT_FAILED_SEE_REDACTED_RECEIPT", call. = FALSE)
  })
}
main(commandArgs(trailingOnly = TRUE))


