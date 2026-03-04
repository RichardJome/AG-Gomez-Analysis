# 05_pathway_analysis.R ------------------------------------------------------------------------
# Pathway Enrichment Analysis for Experiment 5
# Author: Richard Jome
# Date: 2025
# Purpose: GO and KEGG pathway enrichment analysis (UP/DOWN separated)

# Load Libraries -------------------------------------------------------------------
library(clusterProfiler)
library(ggplot2)
library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Define Groups --------------------------------------------------------------------
groups <- list(
  "gf_24h" = "Germ Free 24h",
  "gf_6h" = "Germ Free 6h",
  "spf_24h" = "SPF 24h",
  "spf_6h" = "SPF 6h"
)

comparisons <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")

# Function: Convert ENSEMBL to Entrez IDs -----------------------------------------
#' Convert Ensembl gene IDs to Entrez IDs for pathway analysis
#' @param ensembl_ids Vector of Ensembl gene IDs
#' @return Vector of Entrez gene IDs
convert_ensembl_to_entrez <- function(ensembl_ids) {
  mapped <- tryCatch({
    mapIds(
      org.Mm.eg.db,
      keys = ensembl_ids,
      keytype = "ENSEMBL",
      column = "ENTREZID",
      multiVals = "first"
    )
  }, error = function(e) NULL)
  
  if (is.null(mapped)) {
    return(character(0))
  }
  
  entrez_ids <- na.omit(as.character(mapped))
  return(entrez_ids[entrez_ids != ""])
}

# Function: Run Pathway Analysis for Gene Set -------------------------------------
#' Run GO and KEGG enrichment analysis
#' @param entrez_ids Vector of Entrez gene IDs
#' @param title_base Base title for plots
#' @param results_dir Output directory for plots
#' @param prefix Filename prefix
#' @param min_genes Minimum number of genes required
#' @return List of enrichment results
run_pathway_for_genes <- function(entrez_ids, title_base, results_dir, prefix, min_genes = 3) {
  
  if (length(entrez_ids) < min_genes) {
    cat(sprintf("      %s: Only %d genes, skipping\n", prefix, length(entrez_ids)))
    return(NULL)
  }
  
  cat(sprintf("      Running %s analysis with %d genes...\n", prefix, length(entrez_ids)))
  
  # GO enrichment
  ego <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  if (!is.null(ego) && nrow(ego) > 0) {
    go_file <- file.path(results_dir, paste0("pathway_go_", prefix, ".png"))
    png(go_file, width = 14, height = 12, units = "in", res = 300)
    print(dotplot(ego, showCategory = 15, font.size = 10) + 
            ggtitle(paste("GO Enrichment -", title_base)) +
            theme(plot.title = element_text(size = 12, face = "bold"),
                  axis.text.y = element_text(size = 9)))
    dev.off()
  }
  
  # KEGG enrichment
  ekegg <- enrichKEGG(
    gene       = entrez_ids,
    organism   = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    kegg_file <- file.path(results_dir, paste0("pathway_kegg_", prefix, ".png"))
    png(kegg_file, width = 14, height = 12, units = "in", res = 300)
    print(dotplot(ekegg, showCategory = 15, font.size = 10) + 
            ggtitle(paste("KEGG Enrichment -", title_base)) +
            theme(plot.title = element_text(size = 12, face = "bold"),
                  axis.text.y = element_text(size = 9)))
    dev.off()
  }
  
  return(list(go = ego, kegg = ekegg))
}

# Function: Run Pathway Analysis for Group -----------------------------------------
#' Run pathway analysis for all comparisons in a group
#' @param group_name Directory name for the group
#' @param group_label Display label for the group
#' @return Saves pathway plots to results directory
run_pathway_analysis <- function(group_name, group_label) {
  
  cat(sprintf("\n=== Processing: %s ===\n", group_label))
  
  results_dir <- file.path(group_name, "results")
  rds_dir <- file.path(group_name, "rds")
  
  dir.create(results_dir, showWarnings = FALSE)
  
  for (comp in comparisons) {
    res_file <- file.path(rds_dir, paste0("res_", comp, ".rds"))
    
    if (!file.exists(res_file)) {
      cat(sprintf("  %s: File not found\n", comp))
      next
    }
    
    res <- readRDS(res_file)
    res_df <- as.data.frame(res)
    
    # Standard thresholds: padj < 0.05, |log2FC| > 1.0
    sig_all <- res_df[res_df$padj < 0.05 & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= 1.0, ]
    
    sig_up <- sig_all[sig_all$log2FoldChange > 0, ]
    sig_down <- sig_all[sig_all$log2FoldChange < 0, ]
    
    n_up <- nrow(sig_up)
    n_down <- nrow(sig_down)
    n_total <- nrow(sig_all)
    
    cat(sprintf("  %s: %d significant (%d up, %d down)\n", comp, n_total, n_up, n_down))
    
    title_base <- paste(group_label, "-", comp)
    
    # Upregulated genes
    if (n_up >= 3) {
      cat("    Processing UPREGULATED genes...\n")
      up_genes <- rownames(sig_up)
      entrez_up <- convert_ensembl_to_entrez(up_genes)
      if (length(entrez_up) >= 3) {
        run_pathway_for_genes(entrez_up, paste(title_base, "(UP)"), results_dir, 
                              paste0(comp, "_UP"), min_genes = 3)
      }
    }
    
    # Downregulated genes
    if (n_down >= 3) {
      cat("    Processing DOWNREGULATED genes...\n")
      down_genes <- rownames(sig_down)
      entrez_down <- convert_ensembl_to_entrez(down_genes)
      if (length(entrez_down) >= 3) {
        run_pathway_for_genes(entrez_down, paste(title_base, "(DOWN)"), results_dir, 
                              paste0(comp, "_DOWN"), min_genes = 3)
      }
    }
    
    cat(sprintf("    Completed: %s\n", comp))
  }
}

# Main: Run Pathway Analysis ------------------------------------------------------
cat("=== PATHWAY ANALYSIS (UP/DOWN SEPARATED - PNG OUTPUT) ===\n")

for (group_name in names(groups)) {
  run_pathway_analysis(group_name, groups[[group_name]])
}

cat("\n=== COMPLETE ===\n")

# Save Session Info ---------------------------------------------------------------
sink("session_info_05.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_05.txt\n")
