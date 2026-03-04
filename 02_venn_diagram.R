# 02_venn_diagram.R ------------------------------------------------------------------------
# Venn Diagram Generator for Experiment 5
# Author: Richard Jome
# Date: 2025
# Purpose: Create Venn diagrams showing overlap of differentially expressed genes

# Load Libraries -------------------------------------------------------------------
library(ggVennDiagram)
library(ggplot2)

# Function: Create Venn Diagram ----------------------------------------------------
#' Create Venn diagram showing overlap of DEGs between comparisons
#' @param group_name Name of the experimental group
#' @param input_dir Directory containing DESeq2 results (RDS files)
#' @param output_dir Output directory for Venn diagram PNG
#' @return Saves Venn diagram to PNG file
create_venn <- function(group_name, input_dir, output_dir) {
  
  cat(sprintf("\n=== Creating Venn diagram for: %s ===\n", group_name))
  
  res1 <- readRDS(file.path(input_dir, "rds", "res_Untreated_vs_SG1B.rds"))
  res2 <- readRDS(file.path(input_dir, "rds", "res_Untreated_vs_SG1C.rds"))
  res3 <- readRDS(file.path(input_dir, "rds", "res_SG1B_vs_SG1C.rds"))
  
  # Standard thresholds: padj < 0.05, |log2FC| > 1.0
  sig1 <- subset(as.data.frame(res1), padj < 0.05 & abs(log2FoldChange) > 1.0)
  sig2 <- subset(as.data.frame(res2), padj < 0.05 & abs(log2FoldChange) > 1.0)
  sig3 <- subset(as.data.frame(res3), padj < 0.05 & abs(log2FoldChange) > 1.0)
  
  genes_A <- sig1$gene[!is.na(sig1$gene)]
  genes_B <- sig2$gene[!is.na(sig2$gene)]
  genes_C <- sig3$gene[!is.na(sig3$gene)]
  
  venn_list <- list(
    "A" = genes_A,
    "B" = genes_B,
    "C" = genes_C
  )
  
  p <- ggVennDiagram(venn_list, 
                     category.names = c("Untreated\nvs SG1B", 
                                        "Untreated\nvs SG1C", 
                                        "SG1B\nvs SG1C"))
  
  p <- p + 
    ggtitle(sprintf("Differentially Expressed Genes\n%s", group_name)) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, "venn_diagram.png"), p, dpi = 300, bg = "white", width = 10, height = 8)
  
  cat(sprintf("Venn diagram saved to: %s/venn_diagram.png\n", output_dir))
  cat(sprintf("A (Untreated vs SG1B): %d genes\n", length(genes_A)))
  cat(sprintf("B (Untreated vs SG1C): %d genes\n", length(genes_B)))
  cat(sprintf("C (SG1B vs SG1C): %d genes\n", length(genes_C)))
  
  # Clean up
  rm(res1, res2, res3, sig1, sig2, sig3)
  gc()
  
  return(invisible(NULL))
}

# Create Venn Diagrams ------------------------------------------------------------
create_venn("Germ Free 24h", "gf_24h", "gf_24h/results")
create_venn("Germ Free 6h", "gf_6h", "gf_6h/results")
create_venn("Specific Pathogen Free 24h", "spf_24h", "spf_24h/results")
create_venn("Specific Pathogen Free 6h", "spf_6h", "spf_6h/results")

cat("\n=== All Venn diagrams complete! ===\n")

# Save Session Info ---------------------------------------------------------------
sink("session_info_02.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_02.txt\n")
