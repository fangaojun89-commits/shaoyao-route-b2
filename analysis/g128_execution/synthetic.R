# Generated counts and memberships only. Accession-shaped labels satisfy the
# original frozen 30/29 exclusion interface; they do not represent observations.
run_synthetic <- function(spec, a, v, member_guard, comparison_guard, read_pair) {
  root <- file.path(spec$output, "synthetic_inputs")
  stopifnot(dir.create(root))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20261024L)
  n <- 360L
  genes <- sprintf("SYNTH_G%04d", seq_len(n))
  human <- sprintf("SYNTH_H%04d", seq_len(n))
  mouse <- sprintf("SynthM%04d", seq_len(n))
  ledger <- data.table::data.table(
    sample_order = 1:30, sample_id = paste0("GSM", 3681209:3681238),
    biological_unit_id = paste0("SYNTH_UNIT_", 1:30), biosample = paste0("SYNTH_BIO_", 1:30),
    sra_experiment = paste0("SYNTH_EXP_", 1:30), run_accession = paste0("SRR", 8757486:8757515),
    sex = "female", group = rep(c("D42", "D28", "D14", "D7", "D4", "control"), each = 5L))
  base <- exp(runif(n, log(40), log(600)))
  mu <- outer(base, exp(seq(-0.25, 0.25, length.out = 30L)))
  chronic <- ledger$group %in% c("D28", "D42")
  mu[1:25, chronic] <- mu[1:25, chronic] * 1.8
  mu[26:50, chronic] <- mu[26:50, chronic] / 1.8
  counts <- matrix(rnbinom(n * 30L, mu = as.vector(mu), size = 12), nrow = n)
  counts[(n-9L):n, ] <- 0L
  colnames(counts) <- ledger$sample_id
  annotation <- data.table::data.table(GeneID = genes, gene = mouse)
  resources <- list(
    reactome = data.table::data.table(release = "97", term_id = rep(c("SYNTH_R1", "SYNTH_R2"), each = 20),
      term_name = rep(c("Synthetic program A", "Synthetic program B"), each = 20), gene_symbol = human[1:40]),
    syngo_sensitivity = data.table::data.table(term_id = "SYNTH_S1", term_name = "Synthetic synaptic set",
      gene_symbol = human[41:60]),
    ortholog_one_to_one = data.table::data.table(other_tax_id = 10090L, human_symbol = human, other_symbol = mouse),
    progeny_mouse = data.table::data.table(pathway = "SYNTH_P1", gene_symbol = mouse[61:80],
      weight = rep(c(-0.8, 1.2), 10)),
    dorothea_mouse = data.table::data.table(tf = "SYNTH_TF1", target = mouse[81:100], mor = rep(c(-1, 1), 10)))
  paths <- list(count_matrix = file.path(root, "counts.tsv"), feature_annotation = file.path(root, "annotation.tsv"),
                sample_ledger = file.path(root, "ledger.tsv"), frozen_design = file.path(root, "frozen_design.tsv"),
                gene_set_resource = file.path(root, "resources.rds"), staging = spec$paths$staging)
  data.table::fwrite(data.table::data.table(gene_id = genes, counts), paths$count_matrix, sep = "\t")
  data.table::fwrite(annotation, paths$feature_annotation, sep = "\t")
  data.table::fwrite(ledger, paths$sample_ledger, sep = "\t")
  stopifnot(file.copy(paths$sample_ledger, paths$frozen_design),
            identical(member_guard(resources)$status, "MEMBER_CONTENT_ONLY_PASS"))
  saveRDS(resources, paths$gene_set_resource)
  input_pins <- lapply(paths[names(paths) != "staging"], a$v2a_file_identity)
  inputs <- read_pair(paths)
  written <- a$v2a_execute_model(inputs$analyzer, paths$gene_set_resource, paths$staging,
                                 nrot = 9999L, minimum_members = 10L)
  report <- v$v2v_validate_output(inputs$validator, paths$gene_set_resource, paths$staging,
                                  nrot = 9999L, minimum_members = 10L)
  stopifnot(length(written) == 33L, identical(report$status, "PASS_INDEPENDENT_RECONSTRUCTION"))
  read <- function(relative) data.table::fread(file.path(paths$staging, relative))
  primary <- read("primary_all30/gene_level/edgeR_primary.tsv.gz")
  sensitivity <- read("complete29_sensitivity/gene_level/edgeR_primary.tsv.gz")
  compared <- comparison_guard(
    as.data.frame(primary)[c("gene_id", "contrast", "logFC", "BH_q")],
    as.data.frame(sensitivity)[c("gene_id", "contrast", "logFC", "BH_q")],
    read("comparison/all30_vs_complete29_gene.tsv.gz"),
    read("comparison/all30_vs_complete29_summary.tsv"), n)
  stopifnot(identical(compared$status, "PASS_BRANCH_COMPARISON_CONTENT_ONLY"),
            all(!primary$candidate_generation_allowed),
            all(!sensitivity$candidate_generation_allowed),
            !any(genes[(n-9L):n] %in% primary$gene_id))
  stopifnot(identical(input_pins, lapply(paths[names(paths) != "staging"], a$v2a_file_identity)))
  list(status = "PASS_FULL_SYNTHETIC_MODEL_AND_INDEPENDENT_RECONSTRUCTION",
       synthetic_genes = n, deliberately_zero_genes = 10L, synthetic_programs = 5L,
       branches = c(30L, 29L), contrasts_per_branch = 7L,
       analyzer_output_files = length(written), validation_report = as.list(report[1]),
       branch_comparison = compared$status, real_scientific_inputs_read = FALSE,
       synthetic_inputs_unchanged = TRUE,
       real_data_model_runs = 0L, original_linux_driver_executed = FALSE,
       scope = "Full original statistical kernels on synthetic fixture with copy-isolated validator; not a real-data reproduction")
}
