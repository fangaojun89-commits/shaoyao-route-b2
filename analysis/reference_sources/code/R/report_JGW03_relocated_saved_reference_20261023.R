# Report-only relocation adapter. Original provenance bytes and numeric helpers remain unchanged.
# CLI: RELOCATION_JSON FRESH_OUTPUT. No model/RNG/null-tail or P/q recalculation.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
stopifnot(identical(Sys.getlocale('LC_CTYPE'),'C')) # Required by accepted RDS serialization identity.
sha <- function(p) digest::digest(file=p,algo='sha256',serialize=FALSE)
relocation_path <- normalizePath(args[1],mustWork=TRUE)
relocation <- jsonlite::read_json(relocation_path,simplifyVector=FALSE)
stopifnot(identical(relocation$schema,'jgw-saved-relocation-v1'))
anchor <- relocation$accepted_export
stopifnot(sha(anchor)=='f87400c832bd8b77210076fd6e6231eba823e069dde79d07f6c447216794cfe7')
accepted_export <- jsonlite::read_json(anchor,simplifyVector=FALSE)
stopifnot(identical(accepted_export$status,'PASS_SEPARATE_LINUX_SAVED_SUMMARY_DESCRIPTIVE_EXPORT_ONLY'))
originals <- names(accepted_export$inputs)
stopifnot(setequal(names(relocation$files),originals), !anyDuplicated(names(relocation$files)))
paths <- vapply(relocation$files,function(x) normalizePath(x,mustWork=TRUE),character(1))
stopifnot(!anyDuplicated(paths))
check_inputs <- function() {
  for (p in originals) {
    record <- accepted_export$inputs[[p]]
    stopifnot(identical(record$path,p),file.info(paths[[p]])$size==record$bytes,
              identical(sha(paths[[p]]),record$sha256))
  }
  stopifnot(sha(anchor)=='f87400c832bd8b77210076fd6e6231eba823e069dde79d07f6c447216794cfe7')
  TRUE
}
stopifnot(check_inputs())
helper <- new.env(parent=globalenv())
source(paths[[accepted_export$source$path]],local=helper)
need <- helper$jgwlrep_need
stopifnot(identical(accepted_export$source$sha256,
                  '1a7069210ad804b02cfbf88bcd0bd71cbbba6499f3c768cbc912123615ce49a4'))
binding_key <- originals[endsWith(originals,'/REPORT_BINDING.json')]
stopifnot(length(binding_key)==1L)
binding <- jsonlite::read_json(paths[[binding_key]],simplifyVector=FALSE)
stopifnot(identical(binding$schema,'jgw03-complete-linux-saved-report-binding-v1'),
          identical(names(binding$cells),helper$jgwlrep_cells()))
resolve <- function(record) {
  stopifnot(record$path %in% originals)
  original <- accepted_export$inputs[[record$path]]
  stopifnot(identical(record$sha256,original$sha256),record$bytes==original$bytes)
  paths[[record$path]]
}
readback <- jsonlite::read_json(resolve(binding$accepted_readback),simplifyVector=FALSE)
terminal <- jsonlite::read_json(resolve(binding$actual_readback_parent_exit),simplifyVector=FALSE)
stopifnot(identical(readback$status,'PASS_INDEPENDENT_COMPLETE_LINUX_SAVED_REFERENCE_ONLY'),
          identical(readback$source$sha256,'fd1bb54a1798cf51431dcc49044e3d6af039ed7b6b72549b6c7cbccf6c29e38f'),
          identical(readback$pin_manifest$path,binding$accepted_manifest$path),
          identical(readback$pin_manifest$sha256,binding$accepted_manifest$sha256),
          readback$pin_manifest$bytes==binding$accepted_manifest$bytes,
          terminal$actual_native_exit_code==0,
          identical(terminal$output,dirname(binding$accepted_readback$path)),
          readback$total_MC_families==89,readback$total_companion_families==34,
          readback$new_native_fit_calls==0,isTRUE(readback$all_bound_pins_unchanged))
connection <- gzfile(resolve(binding$accepted_manifest),'rt')
manifest <- tryCatch(read.delim(connection,quote='"',check.names=FALSE,stringsAsFactors=FALSE),
                     finally=close(connection))
stopifnot(identical(names(manifest),c('path','bytes','sha256')),
          nrow(manifest)==readback$bound_files,!anyNA(manifest),!anyDuplicated(manifest$path))
results <- list()
for (accession in helper$jgwlrep_cells()) {
  records <- binding$cells[[accession]]
  stopifnot(setequal(names(records),c('domain','four_group_MC','high_companion_exact')))
  for (record in records) {
    resolve(record)
    # Match the unchanged historical path, not the new storage location.
    helper$jgwlrep_manifest_match(record,manifest)
  }
  domain <- readRDS(resolve(records$domain))
  for (space in c('four_group_MC','high_companion_exact')) {
    results[[paste(accession,space)]] <- helper$jgwlrep_cell(domain,readRDS(resolve(records[[space]])),
                                                          accession,space,readback$cells[[accession]])
  }
}
families <- do.call(rbind,lapply(results,`[[`,'families'))
features <- do.call(rbind,lapply(results,`[[`,'features'))
rownames(families) <- rownames(features) <- NULL
stopifnot(nrow(families)==123L,nrow(features)==484723L)
out <- file.path(normalizePath(dirname(args[2]),mustWork=TRUE),basename(args[2]))
stopifnot(!file.exists(out),!basename(out) %in% c('.','..'),check_inputs(),dir.create(out))
helper$jgwlrep_write_tsv(families,file.path(out,'FAMILY_REFERENCE_SUMMARY.tsv'))
connection <- gzfile(file.path(out,'FEATURE_SAVED_REFERENCE.tsv.gz'),'wt')
tryCatch(helper$jgwlrep_write_tsv(features,connection),finally=close(connection))
stopifnot(check_inputs())
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
jsonlite::write_json(list(status='RELOCATED_SAVED_EXPORT_COMPLETE_COMPARISON_PENDING',
  original_export_sha256=sha(anchor),relocation_sha256=sha(relocation_path),
  family_n=nrow(families),feature_rows=nrow(features),original_input_bytes_unchanged=TRUE,
  scientific_refits=0L,RNG_calls=0L,saved_value_helper_unchanged=TRUE,
  path_mapping=relocation$files),file.path(out,'RESULT.json'),pretty=TRUE,auto_unbox=TRUE)
cat('RELOCATED_EXPORT_COMPLETE: 123 families; 484723 feature rows\n')
