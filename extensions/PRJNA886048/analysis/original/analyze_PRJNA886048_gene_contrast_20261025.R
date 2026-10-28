# New-source gene comparison. No existing source or candidate gate is modified.
# This module is sourced by a driver after quantification/QC acceptance.

prjna886048_design <- function(samples) {
  required <- c("run", "biosample", "treatment")
  stopifnot(all(required %in% names(samples)), nrow(samples) == 12L,
            !anyDuplicated(samples$run), !anyDuplicated(samples$biosample))
  group <- factor(samples$treatment, levels = c("Saline solution", "Baclofen"))
  stopifnot(!anyNA(group), identical(as.integer(table(group)), c(6L, 6L)))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("vehicle", "baclofen")
  rownames(design) <- samples$run
  list(group = group, design = design, contrast = c(vehicle = -1, baclofen = 1))
}

prjna886048_import <- function(samples, quant_root, tx2gene_file) {
  prjna886048_design(samples)
  map <- read.delim(tx2gene_file, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("transcript_id_versioned", "gene_id") %in% names(map)),
            !anyDuplicated(map$transcript_id_versioned), !anyNA(map$gene_id))
  paths <- setNames(file.path(quant_root, samples$run, "abundance.tsv"), samples$run)
  stopifnot(all(file.exists(paths)))
  observed_ids <- NULL
  for (path in paths) {
    x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("target_id", "length", "eff_length", "est_counts", "tpm") %in% names(x)),
              !anyDuplicated(x$target_id), all(x$target_id %in% map$transcript_id_versioned),
              all(is.finite(as.matrix(x[c("length", "eff_length", "est_counts", "tpm")]))),
              all(x$length > 0), all(x$eff_length > 0),
              all(x$est_counts >= 0), all(x$tpm >= 0), sum(x$est_counts) > 0)
    ids <- x$target_id
    if (is.null(observed_ids)) observed_ids <- ids else stopifnot(identical(ids, observed_ids))
  }
  omitted <- setdiff(map$transcript_id_versioned, observed_ids)
  # Native quantification may omit targets shorter than the fixed k=31.
  # Any other omission needs explicit investigation, not silent gene loss.
  if (length(omitted)) {
    stopifnot("length" %in% names(map),
              all(map$length[match(omitted, map$transcript_id_versioned)] < 31))
  }
  pair <- map[, c("transcript_id_versioned", "gene_id")]
  importer <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  raw <- tximport::tximport(paths, type = "kallisto", tx2gene = pair,
                           importer = importer, countsFromAbundance = "no")
  adjusted <- tximport::tximport(paths, type = "kallisto", tx2gene = pair,
                                importer = importer, countsFromAbundance = "lengthScaledTPM")
  stopifnot(identical(colnames(adjusted$counts), samples$run),
            identical(rownames(raw$counts), rownames(adjusted$counts)),
            isTRUE(all.equal(colSums(raw$counts), colSums(adjusted$counts), tolerance = 1e-8)))
  list(raw = raw, length_adjusted = adjusted, omitted_short_transcripts = omitted,
       inputs = data.frame(run = samples$run, file = unname(paths),
                           sha256 = vapply(paths, digest::digest, character(1),
                                           algo = "sha256", file = TRUE), row.names = NULL))
}

prjna886048_fit_genes <- function(counts, samples) {
  setup <- prjna886048_design(samples)
  stopifnot(is.matrix(counts), is.numeric(counts), ncol(counts) == 12L,
            identical(colnames(counts), samples$run),
            !is.null(rownames(counts)), !anyDuplicated(rownames(counts)),
            all(is.finite(counts)), all(counts >= 0), all(colSums(counts) > 0))
  d <- edgeR::DGEList(counts = counts, group = setup$group)
  d <- edgeR::calcNormFactors(d, method = "TMM")
  keep <- edgeR::filterByExpr(d, design = setup$design) &
    rowSums(edgeR::cpm(d) >= 1) >= 6L
  stopifnot(sum(keep) >= 20L)
  retained <- edgeR::calcNormFactors(d[keep, , keep.lib.sizes = FALSE], method = "TMM")
  retained <- edgeR::estimateDisp(retained, setup$design, robust = TRUE)
  ql <- edgeR::glmQLFit(retained, setup$design, robust = TRUE)
  test <- edgeR::glmQLFTest(ql, contrast = setup$contrast)
  primary <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  primary <- data.frame(gene_id = rownames(primary), primary, row.names = NULL)
  v <- limma::voom(retained, setup$design, plot = FALSE)
  fit <- limma::contrasts.fit(limma::lmFit(v, setup$design), setup$contrast)
  fit <- limma::eBayes(fit, robust = TRUE)
  secondary <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "none", confint = TRUE)
  secondary <- data.frame(gene_id = rownames(secondary), secondary, row.names = NULL)
  stopifnot(identical(primary$gene_id, secondary$gene_id),
            isTRUE(all.equal(primary$FDR, p.adjust(primary$PValue, "BH"))))
  list(primary_edgeR = primary, sensitivity_voom = secondary, dge = retained,
       voom = v, gene_filter = data.frame(gene_id = rownames(counts), retained = keep),
       design = setup$design, contrast = setup$contrast,
       input_definition = "tximport lengthScaledTPM count-scale estimates; no second length offset",
       contrast_definition = "baclofen minus vehicle; independent six-versus-six groups",
       scope = "Retrospective source-wide pharmacological context; no S1 candidate admission")
}

prjna886048_save_genes <- function(result, output) {
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  checkpoint <- file.path(output, "gene_model_objects.rds")
  if (file.exists(checkpoint)) {
    stopifnot(identical(readRDS(checkpoint), result))
  } else {
    # Preserve completed fits before exporting tables, so an export failure
    # can resume from the existing object without fitting again.
    temporary <- paste0(checkpoint, ".partial")
    stopifnot(!file.exists(temporary))
    saveRDS(result, temporary)
    stopifnot(identical(readRDS(temporary), result), file.rename(temporary, checkpoint))
  }
  for (key in c("primary_edgeR", "sensitivity_voom", "gene_filter")) {
    write.table(result[[key]], file.path(output, paste0(key, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  invisible(output)
}
