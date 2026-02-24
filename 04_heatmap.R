library(pheatmap)
library(DESeq2)
library(dplyr)
library(readxl)

groups <- list(
  "gf_24h" = list(
    cols = c("GF1_5_24_C1", "GF1_5_24_C2", "GF1_5_24_C3",
             "GF2_5_24_A1", "GF2_5_24_A2", "GF2_5_24_A3",
             "GF3_5_24_B1", "GF3_5_24_B2", "GF3_5_24_B3",
             "GF4_5_24_C1", "GF4_5_24_C2", "GF4_5_24_C3"),
    cond = rep(c("Untreated", "SG1B", "SG1C"), each = 4)
  ),
  "gf_6h" = list(
    cols = c("GF1_5_6h_C1", "GF1_5_6h_C2", "GF1_5_6h_C3",
             "GF2_5_6h_A1", "GF2_5_6h_A2", "GF2_5_6h_A3",
             "GF3_5_6h_B1", "GF3_5_6h_B2", "GF3_5_6h_B3",
             "GF4_5_6h_C1", "GF4_5_6h_C2", "GF4_5_6h_C3"),
    cond = rep(c("Untreated", "SG1B", "SG1C"), each = 4)
  ),
  "spf_24h" = list(
    cols = c("SPF1_5_24_A1", "SPF1_5_24_A2", "SPF1_5_24_A3",
             "SPF2_5_24_B1", "SPF2_5_24_B2", "SPF2_5_24_B3"),
    cond = rep(c("Untreated", "SG1B", "SG1C"), each = 2)
  ),
  "spf_6h" = list(
    cols = c("SPF1_5_6h_A1", "SPF1_5_6h_A2", "SPF1_5_6h_A3",
             "SPF2_5_6h_B1", "SPF2_5_6h_B2", "SPF2_5_6h_B3"),
    cond = rep(c("Untreated", "SG1B", "SG1C"), each = 2)
  )
)

df <- read_excel("Experiment_5.xlsx")

sig_genes <- readRDS("spf_24h/rds/res_Untreated_vs_SG1C.rds")
sig_genes <- subset(as.data.frame(sig_genes), padj < 0.05 & abs(log2FoldChange) > 1)
sig_genes <- sig_genes[order(sig_genes$padj), ]
top_genes <- head(sig_genes$gene, 50)

cat(sprintf("Using top %d significant genes from SPF 24h (Untreated vs SG1C - pathogenic)\n", length(top_genes)))

all_vst <- list()
all_conditions <- c()

for (group_name in names(groups)) {
  vst_data <- readRDS(file.path(group_name, "rds", "vst_counts.rds"))
  sample_cols <- groups[[group_name]]$cols
  
  vst_subset <- vst_data[, c("gene_name", sample_cols)]
  vst_subset <- vst_subset %>% filter(gene_name %in% top_genes)
  
  gene_names <- vst_subset$gene_name
  expr_matrix <- as.matrix(vst_subset[, sample_cols])
  rownames(expr_matrix) <- gene_names
  
  all_vst[[group_name]] <- expr_matrix
  all_conditions <- c(all_conditions, groups[[group_name]]$cond)
}

combined_matrix <- do.call(cbind, all_vst)

annotation_col <- data.frame(
  Condition = factor(all_conditions, levels = c("Untreated", "SG1B", "SG1C")),
  Time = rep(c("24h", "6h", "24h", "6h"), c(12, 12, 6, 6)),
  Microbiome = rep(c("GF", "GF", "SPF", "SPF"), c(12, 12, 6, 6))
)

rownames(annotation_col) <- colnames(combined_matrix)

annotation_colors <- list(
  Condition = c(Untreated = "#808080", SG1B = "#2E8B57", SG1C = "#DC143C"),
  Time = c("24h" = "#1E90FF", "6h" = "#FFA500"),
  Microbiome = c(GF = "#9370DB", SPF = "#20B2AA")
)

dir.create("heatmaps", showWarnings = FALSE)

pheatmap(combined_matrix,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = 45,
         main = "Top 50 DE Genes: SPF 24h (Untreated vs SG1C) across all groups",
         filename = "heatmaps/heatmap_SPF24h_Untreated_vs_SG1C_top50.png",
         width = 14,
         height = 10)

cat("\nHeatmap saved to: heatmaps/heatmap_all_groups_top50.png\n")

pheatmap(combined_matrix[1:25, ],
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         fontsize_row = 9,
         fontsize_col = 8,
         angle_col = 45,
         main = "Top 25 DE Genes: SPF 24h (Untreated vs SG1C) across all groups",
         filename = "heatmaps/heatmap_SPF24h_Untreated_vs_SG1C_top25.png",
         width = 14,
         height = 8)

cat("Heatmap saved to: heatmaps/heatmap_all_groups_top25.png\n")

cat("\n=== Heatmap generation complete! ===\n")
cat("\nSummary:\n")
cat("- Genes selected from: SPF 24h (Untreated vs SG1C - pathogenic strain)\n")
cat("- Shows how these genes behave in:\n")
cat("  * Germ Free (GF) mice at 24 hours\n")
cat("  * Germ Free (GF) mice at 6 hours\n")
cat("  * Specific Pathogen Free (SPF) mice at 24 hours\n")
cat("  * Specific Pathogen Free (SPF) mice at 6 hours\n")
