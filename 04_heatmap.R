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

contrasts <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")

all_conditions <- c()
for (group_name in names(groups)) {
  all_conditions <- c(all_conditions, groups[[group_name]]$cond)
}

annotation_col <- data.frame(
  Condition = factor(all_conditions, levels = c("Untreated", "SG1B", "SG1C")),
  Time = rep(c("24h", "6h", "24h", "6h"), c(12, 12, 6, 6)),
  Microbiome = rep(c("GF", "GF", "SPF", "SPF"), c(12, 12, 6, 6))
)

annotation_colors <- list(
  Condition = c(Untreated = "#808080", SG1B = "#2E8B57", SG1C = "#DC143C"),
  Time = c("24h" = "#1E90FF", "6h" = "#FFA500"),
  Microbiome = c(GF = "#9370DB", SPF = "#20B2AA")
)

all_vst <- list()
for (group_name in names(groups)) {
  vst_data <- readRDS(file.path(group_name, "rds", "vst_counts.rds"))
  sample_cols <- groups[[group_name]]$cols
  vst_subset <- vst_data[, c("gene_name", sample_cols)]
  all_vst[[group_name]] <- vst_subset
}

generate_heatmap <- function(source_group, contrast, top_n, genes, total_sig) {
  n_genes <- length(genes)
  
  if (n_genes < 2) {
    cat(sprintf("    Only %d gene(s) available, skipping (need >=2 for clustering)\n", n_genes))
    return(NULL)
  }
  
  cat(sprintf("    Generating heatmap with %d genes\n", n_genes))
  
  combined_matrix <- NULL
  col_names_all <- character()
  
  for (group_name in names(groups)) {
    vst_subset <- all_vst[[group_name]] %>% filter(gene_name %in% genes)
    
    if (nrow(vst_subset) == 0) {
      next
    }
    
    sample_cols <- groups[[group_name]]$cols
    gene_names <- vst_subset$gene_name
    expr_matrix <- as.matrix(vst_subset[, sample_cols])
    rownames(expr_matrix) <- gene_names
    
    if (is.null(combined_matrix)) {
      combined_matrix <- expr_matrix
    } else {
      combined_matrix <- cbind(combined_matrix, expr_matrix)
    }
    col_names_all <- c(col_names_all, sample_cols)
  }
  
  if (is.null(combined_matrix) || nrow(combined_matrix) == 0) {
    cat("    No genes found in VST data, skipping...\n")
    return(NULL)
  }
  
  col_annot <- annotation_col[1:ncol(combined_matrix), , drop = FALSE]
  rownames(col_annot) <- colnames(combined_matrix)
  
  if (n_genes == total_sig) {
    title_suffix <- sprintf("(%d genes)", total_sig)
  } else {
    title_suffix <- sprintf("top %d of %d genes", n_genes, total_sig)
  }
  
  output_file <- sprintf("%s/results/heatmap_%s_top%d.png", source_group, contrast, n_genes)
  
  pheatmap(combined_matrix,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "ward.D2",
           annotation_col = col_annot,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 8,
           angle_col = 45,
           main = sprintf("DE Genes: %s (%s) - %s across all groups", source_group, contrast, title_suffix),
           filename = output_file,
           width = 14,
           height = 10)
  
  cat(sprintf("    Saved: %s\n", output_file))
  
  return(output_file)
}

cat("\n=== Generating heatmaps for all comparisons ===\n\n")

for (group_name in names(groups)) {
  cat(sprintf("--- %s ---\n", group_name))
  
  dir.create(file.path(group_name, "results"), showWarnings = FALSE)
  
  for (contrast in contrasts) {
    res_file <- file.path(group_name, "rds", paste0("res_", contrast, ".rds"))
    
    if (!file.exists(res_file)) {
      cat(sprintf("  %s: File not found, skipping\n", contrast))
      next
    }
    
    sig_genes <- readRDS(res_file)
    sig_genes <- subset(as.data.frame(sig_genes), padj < 0.05 & abs(log2FoldChange) > 1)
    sig_genes <- sig_genes[order(sig_genes$padj), ]
    
    n_sig <- nrow(sig_genes)
    
    if (n_sig == 0) {
      cat(sprintf("  %s: 0 significant genes - skipping\n", contrast))
      next
    }
    
    cat(sprintf("  %s: %d significant genes (padj<0.05, |log2FC|>1)\n", contrast, n_sig))
    
    if (n_sig >= 50) {
      top_genes_50 <- head(sig_genes$gene, 50)
      top_genes_25 <- head(sig_genes$gene, 25)
      generate_heatmap(group_name, contrast, 50, top_genes_50, n_sig)
      generate_heatmap(group_name, contrast, 25, top_genes_25, n_sig)
    } else if (n_sig >= 25) {
      top_genes_25 <- head(sig_genes$gene, 25)
      generate_heatmap(group_name, contrast, 25, top_genes_25, n_sig)
    } else {
      top_genes_n <- head(sig_genes$gene, n_sig)
      generate_heatmap(group_name, contrast, n_sig, top_genes_n, n_sig)
    }
  }
  cat("\n")
}

if (dir.exists("heatmaps")) {
  unlink("heatmaps", recursive = TRUE)
  cat("\nDeleted: heatmaps/ folder\n")
}

cat("\n=== Heatmap generation complete! ===\n")
