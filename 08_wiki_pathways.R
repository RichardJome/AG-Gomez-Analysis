# 08_wiki_pathways.R ------------------------------------------------------------------------
# WikiPathways Analysis for Experiment 5
# Analyzes WikiPathways enrichment for UP/DOWN genes separately

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

# Parameters (standard thresholds)
fdr_threshold <- 0.05   # padj < 0.05 (standard)
log2fc_threshold <- 1.0  # |log2FC| > 1.0 (standard 2-fold change)

# Define groups
groups <- list(
  "gf_24h" = "Germ Free 24h",
  "gf_6h" = "Germ Free 6h", 
  "spf_24h" = "SPF 24h",
  "spf_6h" = "SPF 6h"
)

comparisons <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")

# Function to convert Ensembl to Entrez ---------------------------------------------
convert_ensembl_to_entrez <- function(genes) {
  result <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  return(result$ENTREZID)
}

# Function: Run WikiPathways analysis -----------------------------------------------
run_wiki_pathway_analysis <- function(group_name, group_label) {
  
  cat(sprintf("\n=== Processing: %s ===\n", group_label))
  
  rds_dir <- file.path(group_name, "rds")
  results_dir <- file.path(group_name, "results")
  dir.create(results_dir, showWarnings = FALSE)
  
  for (comp in comparisons) {
    res_file <- file.path(rds_dir, paste0("res_", comp, ".rds"))
    
    if (!file.exists(res_file)) {
      cat(sprintf("  %s: File not found\n", comp))
      next
    }
    
    res <- readRDS(res_file)
    res_df <- as.data.frame(res)
    
    sig_all <- res_df[res_df$padj < fdr_threshold & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= log2fc_threshold, ]
    
    sig_up <- sig_all[sig_all$log2FoldChange > 0, ]
    sig_down <- sig_all[sig_all$log2FoldChange < 0, ]
    
    n_up <- nrow(sig_up)
    n_down <- nrow(sig_down)
    n_total <- nrow(sig_all)
    
    cat(sprintf("  %s: %d significant (%d up, %d down)\n", comp, n_total, n_up, n_down))
    
    title_base <- paste(group_label, "-", comp)
    
    # Upregulated genes
    if (n_up >= 3) {
      cat("    Processing UPREGULATED genes (WikiPathways)...\n")
      up_genes <- rownames(sig_up)
      entrez_up <- convert_ensembl_to_entrez(up_genes)
      if (length(entrez_up) >= 3) {
        tryCatch({
          wp_result <- enrichWP(entrez_up, organism = "Mus musculus", pvalueCutoff = 0.05)
          if (!is.null(wp_result) && nrow(wp_result) > 0) {
            wp_result@result$Description <- substr(wp_result@result$Description, 1, 50)
            p <- dotplot(wp_result, showCategory = 12, font.size = 10) + 
              ggtitle(paste(title_base, "(UP) - WikiPathways")) +
              theme(plot.title = element_text(size = 11, face = "bold"),
                    axis.text.y = element_text(size = 9))
            ggsave(file.path(results_dir, paste0("wiki_", comp, "_UP.png")), p, width = 12, height = 10, dpi = 300)
            cat(sprintf("    Saved: wiki_%s_UP.png\n", comp))
          }
        }, error = function(e) {
          cat(sprintf("    Error in UP: %s\n", e$message))
        })
      }
    }
    
    # Downregulated genes
    if (n_down >= 3) {
      cat("    Processing DOWNREGULATED genes (WikiPathways)...\n")
      down_genes <- rownames(sig_down)
      entrez_down <- convert_ensembl_to_entrez(down_genes)
      if (length(entrez_down) >= 3) {
        tryCatch({
          wp_result <- enrichWP(entrez_down, organism = "Mus musculus", pvalueCutoff = 0.05)
          if (!is.null(wp_result) && nrow(wp_result) > 0) {
            wp_result@result$Description <- substr(wp_result@result$Description, 1, 50)
            p <- dotplot(wp_result, showCategory = 12, font.size = 10) + 
              ggtitle(paste(title_base, "(DOWN) - WikiPathways")) +
              theme(plot.title = element_text(size = 11, face = "bold"),
                    axis.text.y = element_text(size = 9))
            ggsave(file.path(results_dir, paste0("wiki_", comp, "_DOWN.png")), p, width = 12, height = 10, dpi = 300)
            cat(sprintf("    Saved: wiki_%s_DOWN.png\n", comp))
          }
        }, error = function(e) {
          cat(sprintf("    Error in DOWN: %s\n", e$message))
        })
      }
    }
  }
  cat(sprintf("  Completed: %s\n", group_label))
}

# Main: Run WikiPathways Analysis --------------------------------------------------
cat("=== WIKIPATHWAYS ANALYSIS ===\n")

for (group_name in names(groups)) {
  run_wiki_pathway_analysis(group_name, groups[[group_name]])
}

cat("\n=== WIKIPATHWAYS COMPLETE ===\n")

sink("session_info_08.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_08.txt\n")
