# Explicit scientific execution; called only by g205_model.py --execute.
# The original preparation script defines spec, worker, handle and capture below.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
entry_spec <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
prepare_path <- file.path(entry_spec$repo, 'analysis/g205_preparation/prepare.R')
prepare_raw <- readBin(prepare_path, 'raw', n = file.info(prepare_path)$size)
prepare_pin <- entry_spec$source_pins[['analysis/g205_preparation/prepare.R']]
stopifnot(length(prepare_raw) == prepare_pin$bytes,
  identical(digest::digest(prepare_raw, 'sha256', serialize = FALSE), prepare_pin$sha256))
eval(parse(text = rawToChar(prepare_raw), keep.source = FALSE), envir = .GlobalEnv)
writer_raw <- capture('analysis/g205_execution/write_result.R')
eval(parse(text = rawToChar(writer_raw), keep.source = FALSE), envir = .GlobalEnv)
checked <- handle$execute()
provenance <- list(scope = 'NEW_PATH_PARAMETERIZED_PRIMARY_EXECUTION',
  historical_native_custody_recreated = FALSE,
  input_preparation = spec$prepared_identity, sources = spec$source_pins,
  current_runtime = runtime,
  limitation = 'Original deterministic preparation and checked consumer; new outputs, not historical completion markers or original Linux process custody.')
plan <- worker$gse205906_ens101_worker_write_plan_v1(checked, provenance)
g205_write_result(plan, file.path(out, 'results'), worker$gse205906_ens101_worker_pin_v1)
qa <- list(status = 'PASS_ORIGINAL_CONSUMER_AND_OUTPUT_READBACK',
  scientific_execution = TRUE, original_consumer_receipt = handle$inspect(),
  table_names = names(checked$payload$result$tables),
  file_pins = plan$file_pins, total_bytes = plan$total_bytes,
  historical_native_custody_recreated = FALSE,
  numerical_equivalence_to_historical_results_established = FALSE)
jsonlite::write_json(qa, file.path(out, 'receipts/EXECUTION_QA.json'), auto_unbox = TRUE,
  pretty = TRUE, null = 'null', digits = NA)
writeLines(capture.output(sessionInfo()), file.path(out, 'receipts/sessionInfo.txt'), useBytes = TRUE)
