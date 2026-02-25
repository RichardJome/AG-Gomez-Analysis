library(clusterProfiler)
library(ggplot2)
library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)

groups <- list(
  "gf_24h" = "Germ Free 24h",
  "gf_6h" = "Germ Free 6h",
  "spf_24h" = "SPF 24h",
  "spf_6h" = "SPF 6h"
)

comparisons <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")

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
    
    sig_all <- res_df[res_df$padj < 0.05 & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= 1, ]
    
    if (nrow(sig_all) < 3) {
      cat(sprintf("  %s: Only %d significant genes, skipping\n", comp, nrow(sig_all)))
      next
    }
    
    sig_genes <- rownames(sig_all)
    cat(sprintf("  %s: %d significant genes\n", comp, length(sig_genes)))
    
    entrez_ids <- convert_ensembl_to_entrez(sig_genes)
    cat(sprintf("    -> Converted to Entrez IDs: %d\n", length(entrez_ids)))
    
    if (length(entrez_ids) < 3) {
      cat(sprintf("    -> Not enough mapped genes, skipping\n"))
      next
    }
    
    title_base <- paste(group_label, "-", comp)
    
    cat("    Running GO enrichment...\n")
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
      cat("    Saving GO plot...\n")
      go_file <- file.path(results_dir, paste0("pathway_go_", comp, ".png"))
      png(go_file, width = 14, height = 10, units = "in", res = 300)
      print(dotplot(ego, showCategory = 20) + ggtitle(paste("GO Enrichment -", title_base)))
      dev.off()
    }
    
    cat("    Running KEGG enrichment...\n")
    ekegg <- enrichKEGG(
      gene       = entrez_ids,
      organism   = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      cat("    Saving KEGG plot...\n")
      kegg_file <- file.path(results_dir, paste0("pathway_kegg_", comp, ".png"))
      png(kegg_file, width = 14, height = 10, units = "in", res = 300)
      print(dotplot(ekegg, showCategory = 20) + ggtitle(paste("KEGG Enrichment -", title_base)))
      dev.off()
    }
    
    cat(sprintf("    Completed: %s\n", comp))
  }
}

cat("=== PATHWAY ANALYSIS (Separate Comparisons - PNG OUTPUT) ===\n")

for (group_name in names(groups)) {
  run_pathway_analysis(group_name, groups[[group_name]])
}

cat("\n=== COMPLETE ===\n")
