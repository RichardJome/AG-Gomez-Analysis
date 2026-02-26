# 03_volcano_plots.R ------------------------------------------------------------------------
# Volcano Plot Generator for Experiment 5
# Author: Richard Jome
# Date: 2025
# Purpose: Create volcano plots for differential expression results

# Load Libraries -------------------------------------------------------------------
library(EnhancedVolcano)

# Function: Create Volcano Plots --------------------------------------------------
#' Create volcano plots for DESeq2 results
#' @param group_name Name of the experimental group
#' @param input_dir Directory containing DESeq2 results (RDS files)
#' @param output_dir Output directory for volcano plot PNGs
#' @return Saves volcano plots to PNG files
create_volcano_plots <- function(group_name, input_dir, output_dir) {
  
  cat(sprintf("\n=== Creating volcano plots for: %s ===\n", group_name))
  
  res1 <- readRDS(file.path(input_dir, "rds", "res_Untreated_vs_SG1B.rds"))
  res2 <- readRDS(file.path(input_dir, "rds", "res_Untreated_vs_SG1C.rds"))
  res3 <- readRDS(file.path(input_dir, "rds", "res_SG1B_vs_SG1C.rds"))
  
  df1 <- as.data.frame(res1)
  df2 <- as.data.frame(res2)
  df3 <- as.data.frame(res3)
  
  # Calculate y-axis limits
  ylim1_val <- max(-log10(min(df1$padj[df1$padj > 0 & df1$padj < 1], na.rm = TRUE)) + 1, 5)
  ylim2_val <- max(-log10(min(df2$padj[df2$padj > 0 & df2$padj < 1], na.rm = TRUE)) + 1, 5)
  ylim3_val <- max(-log10(min(df3$padj[df3$padj > 0 & df3$padj < 1], na.rm = TRUE)) + 1, 5)
  
  ylim1 <- c(0, ylim1_val)
  ylim2 <- c(0, ylim2_val)
  ylim3 <- c(0, ylim3_val)
  
  # Extract significant genes for labeling
  sig_df1 <- df1[df1$padj < 0.05 & abs(df1$log2FoldChange) > 1 & !is.na(df1$gene) & !is.na(df1$padj) & !is.na(df1$log2FoldChange), ]
  sig_df1 <- sig_df1[order(sig_df1$padj), ]
  sig_genes1 <- sig_df1$gene
  
  sig_df2 <- df2[df2$padj < 0.05 & abs(df2$log2FoldChange) > 1 & !is.na(df2$gene) & !is.na(df2$padj) & !is.na(df2$log2FoldChange), ]
  sig_df2 <- sig_df2[order(sig_df2$padj), ]
  sig_genes2 <- sig_df2$gene
  
  sig_df3 <- df3[df3$padj < 0.05 & abs(df3$log2FoldChange) > 1 & !is.na(df3$gene) & !is.na(df3$padj) & !is.na(df3$log2FoldChange), ]
  sig_df3 <- sig_df3[order(sig_df3$padj), ]
  sig_genes3 <- sig_df3$gene
  
  # Volcano plot 1: Untreated vs SG1B
  p1 <- EnhancedVolcano(df1,
    lab = df1$gene,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = sig_genes1,
    xlim = c(-10, 10),
    ylim = ylim1,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'p-value'),
    title = sprintf('Untreated vs SG1B - %s', group_name),
    subtitle = 'Differentially expressed genes',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    vline = c(-1, 1),
    vlineType = "longdash",
    vlineCol = "black",
    vlineWidth = 0.5)
  
  # Volcano plot 2: Untreated vs SG1C
  p2 <- EnhancedVolcano(df2,
    lab = df2$gene,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = sig_genes2,
    xlim = c(-10, 10),
    ylim = ylim2,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'p-value'),
    title = sprintf('Untreated vs SG1C - %s', group_name),
    subtitle = 'Differentially expressed genes',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    vline = c(-1, 1),
    vlineType = "longdash",
    vlineCol = "black",
    vlineWidth = 0.5)
  
  # Volcano plot 3: SG1B vs SG1C
  p3 <- EnhancedVolcano(df3,
    lab = df3$gene,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = sig_genes3,
    xlim = c(-10, 10),
    ylim = ylim3,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'p-value'),
    title = sprintf('SG1B vs SG1C - %s', group_name),
    subtitle = 'Differentially expressed genes',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    vline = c(-1, 1),
    vlineType = "longdash",
    vlineCol = "black",
    vlineWidth = 0.5)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(output_dir, "volcano_Untreated_vs_SG1B.png"), p1, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "volcano_Untreated_vs_SG1C.png"), p2, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "volcano_SG1B_vs_SG1C.png"), p3, width = 10, height = 8, dpi = 300)
  
  cat(sprintf("Volcano plots saved to: %s/\n", output_dir))
  cat(sprintf("  - Untreated vs SG1B: %d significant genes\n", length(sig_genes1)))
  cat(sprintf("  - Untreated vs SG1C: %d significant genes\n", length(sig_genes2)))
  cat(sprintf("  - SG1B vs SG1C: %d significant genes\n", length(sig_genes3)))
  
  # Clean up
  rm(df1, df2, df3, sig_df1, sig_df2, sig_df3)
  gc()
  
  return(invisible(NULL))
}

# Create Volcano Plots ------------------------------------------------------------
create_volcano_plots("GF 24h", "gf_24h", "gf_24h/results")
create_volcano_plots("GF 6h", "gf_6h", "gf_6h/results")
create_volcano_plots("SPF 24h", "spf_24h", "spf_24h/results")
create_volcano_plots("SPF 6h", "spf_6h", "spf_6h/results")

cat("\n=== All volcano plots complete! ===\n")

# Save Session Info ---------------------------------------------------------------
sink("session_info_03.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_03.txt\n")
