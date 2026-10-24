#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(file_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})

result_root <- file.path(project_root, "results", "GSE277064")
figure_root <- file.path(result_root, "figures", "programs")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

activity <- fread(file.path(result_root, "programs", "all_celltypes_activity_contrasts.tsv.gz"))
counts <- fread(file.path(result_root, "programs", "exploratory_consensus_counts_not_S3.tsv"))
boundary <- fread(file.path(result_root, "boundary_rotarod", "programs", "all_celltypes_activity_contrasts.tsv.gz"))

cell_labels <- c(v2a = "V2a", astrocyte = "Astrocyte", microglia = "Microglia")
contrast_labels <- c(
  D3_vs_control = "D3", D7_vs_control = "D7", D14_vs_control = "D14",
  D28_vs_control = "D28", chronic_combined = "D14/D28 mean"
)

save_both <- function(plot, filename, width, height) {
  ggsave(file.path(figure_root, paste0(filename, ".png")), plot, width = width, height = height, dpi = 240)
  ggsave(file.path(figure_root, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
}

progeny_time <- activity[
  collection == "PROGENy" & contrast %chin% names(contrast_labels)
]
progeny_time[, cell_type_label := factor(cell_labels[cell_type], levels = unname(cell_labels))]
progeny_time[, contrast_label := factor(contrast_labels[contrast], levels = unname(contrast_labels))]
progeny_time[, pathway := factor(program_id, levels = rev(sort(unique(program_id))))]
progeny_time[, display_effect := pmax(pmin(effect, 2.5), -2.5)]
progeny_time[, marker := fifelse(q_within_celltype_contrast <= 0.05, "•", "")]

p_time <- ggplot(progeny_time, aes(contrast_label, pathway, fill = display_effect)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = marker), size = 2.6) +
  facet_wrap(~ cell_type_label, nrow = 1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
                       limits = c(-2.5, 2.5), name = "Activity\neffect") +
  coord_fixed(ratio = 0.55) +
  theme_bw(base_size = 9) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), strip.background = element_rect(fill = "grey95")) +
  labs(title = "PROGENy activity across the natural post-stroke course",
       subtitle = "Exploratory all-program analysis; dot denotes within-cell-type/contrast BH q<=0.05; not S3")
save_both(p_time, "PROGENy_natural_timecourse_not_S3", 10.5, 6.6)

progeny_chronic <- activity[collection == "PROGENy" & contrast == "chronic_combined"]
progeny_chronic[, cell_type_label := factor(cell_labels[cell_type], levels = unname(cell_labels))]
progeny_chronic[, pathway := factor(program_id, levels = rev(sort(unique(program_id))))]
progeny_chronic[, display_effect := pmax(pmin(effect, 2.5), -2.5)]
progeny_chronic[, marker := fifelse(q_global_chronic_all_locked_precompute <= 0.05, "•", "")]
p_chronic <- ggplot(progeny_chronic, aes(cell_type_label, pathway, fill = display_effect)) +
  geom_tile(color = "white", linewidth = 0.3) + geom_text(aes(label = marker), size = 3) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
                       limits = c(-2.5, 2.5), name = "Chronic\neffect") +
  coord_fixed(ratio = 0.8) + theme_bw(base_size = 9) +
  theme(axis.title = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Chronic-combined PROGENy activity by RiboTag cell type",
       subtitle = "Dot denotes all-locked global BH q<=0.05; exploratory only; S3 not identifiable")
save_both(p_chronic, "PROGENy_chronic_by_celltype_not_S3", 6.8, 6.7)

counts[, cell_type_label := factor(cell_labels[cell_type], levels = unname(cell_labels))]
counts[, collection_label := factor(collection, levels = c("Reactome", "PROGENy", "DoRothEA_AB"))]
p_counts <- ggplot(counts, aes(cell_type_label, strict_all_locked_exploratory_support_n,
                               fill = collection_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = strict_all_locked_exploratory_support_n),
            position = position_dodge(width = 0.8), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c(Reactome = "#4477AA", PROGENy = "#EE6677", DoRothEA_AB = "#228833"),
                    name = "Collection") +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank()) +
  labs(y = "Strict all-locked exploratory programs",
       title = "All-program chronic consistency counts",
       subtitle = "Descriptive only: eligible S2 input=0; S3 pass count=NA")
save_both(p_counts, "exploratory_consensus_counts_not_S3", 7.6, 5.2)

boundary_progeny <- boundary[collection == "PROGENy"]
boundary_labels <- c(
  H2_vs_D14_untrained = "2 h", H4_vs_D14_untrained = "4 h", H8_vs_D14_untrained = "8 h"
)
boundary_progeny[, cell_type_label := factor(cell_labels[cell_type], levels = c("V2a", "Astrocyte"))]
boundary_progeny[, contrast_label := factor(boundary_labels[contrast], levels = unname(boundary_labels))]
boundary_progeny[, pathway := factor(program_id, levels = rev(sort(unique(program_id))))]
boundary_progeny[, display_effect := pmax(pmin(effect, 2.5), -2.5)]
boundary_progeny[, marker := fifelse(q_within_celltype_contrast <= 0.05, "•", "")]
p_boundary <- ggplot(boundary_progeny, aes(contrast_label, pathway, fill = display_effect)) +
  geom_tile(color = "white", linewidth = 0.25) + geom_text(aes(label = marker), size = 2.8) +
  facet_wrap(~ cell_type_label, nrow = 1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
                       limits = c(-2.5, 2.5), name = "Activity\neffect") +
  coord_fixed(ratio = 0.7) + theme_bw(base_size = 9) +
  theme(axis.title = element_blank(), panel.grid = element_blank()) +
  labs(title = "Post-rotarod PROGENy activity relative to untrained D14",
       subtitle = "Boundary-domain recovery response only; dot denotes within-contrast BH q<=0.05; never an S/E gate")
save_both(p_boundary, "PROGENy_rotarod_boundary_only", 8.4, 6.5)

manifest_path <- file.path(result_root, "audit", "run_manifest.tsv")
manifest <- fread(manifest_path)
manifest <- manifest[manifest_key != "program_plot_script_sha256"]
manifest <- rbind(manifest, data.table(
  manifest_key = "program_plot_script_sha256",
  value = digest(script_path, algo = "sha256", file = TRUE)
))
fwrite(manifest, manifest_path, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote four GSE277064 exploratory program figures in PNG and PDF.\n")

