# 09_gene_histogram.R ------------------------------------------------------------------------
# Histogram of UP/DOWN/COMMON genes for Experiment 5
# Shows genes shared and unique between GF vs SPF mice and conditions

library(dplyr)
library(ggplot2)

# Parameters ------------------------------------------------------------------------
fdr_threshold <- 0.2
log2fc_threshold <- 0.5

# Function to get significant genes -------------------------------------------------
get_sig_genes <- function(group_name, comparison, direction = "all") {
  res_file <- file.path(group_name, "rds", paste0("res_", comparison, ".rds"))
  
  if (!file.exists(res_file)) {
    return(character(0))
  }
  
  res <- readRDS(res_file)
  res_df <- as.data.frame(res)
  
  sig <- res_df[res_df$padj < fdr_threshold & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= log2fc_threshold, ]
  
  if (direction == "up") {
    sig <- sig[sig$log2FoldChange > 0, ]
  } else if (direction == "down") {
    sig <- sig[sig$log2FoldChange < 0, ]
  }
  
  genes <- sig$gene[!is.na(sig$gene)]
  return(genes)
}

# Function to create histogram ------------------------------------------------------
create_gene_histogram <- function() {
  
  cat("\n=== Creating Gene Distribution Histogram ===\n")
  
  # Get significant genes for each group and comparison
  # For Untreated vs SG1C (pathogenic) as main comparison
  
  # GF 24h
  gf24_up <- get_sig_genes("gf_24h", "Untreated_vs_SG1C", "up")
  gf24_down <- get_sig_genes("gf_24h", "Untreated_vs_SG1C", "down")
  
  # SPF 24h
  spf24_up <- get_sig_genes("spf_24h", "Untreated_vs_SG1C", "up")
  spf24_down <- get_sig_genes("spf_24h", "Untreated_vs_SG1C", "down")
  
  # Calculate overlaps
  up_common <- intersect(gf24_up, spf24_up)
  down_common <- intersect(gf24_down, spf24_down)
  
  gf24_only_up <- setdiff(gf24_up, spf24_up)
  spf24_only_up <- setdiff(spf24_up, gf24_up)
  
  gf24_only_down <- setdiff(gf24_down, spf24_down)
  spf24_only_down <- setdiff(spf24_down, gf24_down)
  
  cat(sprintf("\nUntreated vs SG1C (24h):\n"))
  cat(sprintf("  GF Upregulated: %d\n", length(gf24_up)))
  cat(sprintf("  SPF Upregulated: %d\n", length(spf24_up)))
  cat(sprintf("  Common UP: %d\n", length(up_common)))
  cat(sprintf("  GF Downregulated: %d\n", length(gf24_down)))
  cat(sprintf("  SPF Downregulated: %d\n", length(spf24_down)))
  cat(sprintf("  Common DOWN: %d\n", length(down_common)))
  
  # Create data for bar plot
  plot_data <- data.frame(
    Category = c("GF Only UP", "SPF Only UP", "Common UP", "GF Only DOWN", "SPF Only DOWN", "Common DOWN"),
    Count = c(length(gf24_only_up), length(spf24_only_up), length(up_common),
              length(gf24_only_down), length(spf24_only_down), length(down_common)),
    Type = c("Upregulated", "Upregulated", "Upregulated", "Downregulated", "Downregulated", "Downregulated")
  )
  
  # Bar plot
  p1 <- ggplot(plot_data, aes(x = Category, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Count), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")) +
    labs(title = "Gene Distribution: GF vs SPF (24h, Untreated vs SG1C)",
         x = "Category", y = "Number of Genes", fill = "Direction") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("gene_distribution_histogram.png", p1, width = 12, height = 8, dpi = 300)
  cat("\nSaved: gene_distribution_histogram.png\n")
  
  # Also create comparison table
  cat("\n=== DETAILED GENE COUNTS ===\n")
  
  # All comparisons for both time points
  all_results <- data.frame()
  
  for (group in c("gf_24h", "gf_6h", "spf_24h", "spf_6h")) {
    for (comp in c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")) {
      up <- get_sig_genes(group, comp, "up")
      down <- get_sig_genes(group, comp, "down")
      all <- c(up, down)
      
      all_results <- rbind(all_results, data.frame(
        Group = group,
        Comparison = comp,
        UP = length(up),
        DOWN = length(down),
        TOTAL = length(all)
      ))
    }
  }
  
  print(all_results)
  
  # Save table
  write.csv(all_results, "gene_counts_summary.csv", row.names = FALSE)
  cat("\nSaved: gene_counts_summary.csv\n")
  
  # Create grouped bar chart
  p2 <- ggplot(all_results, aes(x = Comparison, y = TOTAL, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = TOTAL), vjust = -0.5, size = 3) +
    labs(title = "Total Significant Genes by Group and Comparison",
         x = "Comparison", y = "Number of Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("gene_counts_by_group.png", p2, width = 12, height = 8, dpi = 300)
  cat("Saved: gene_counts_by_group.png\n")
  
  return(invisible(all_results))
}

# Run ------------------------------------------------------------------------
create_gene_histogram()

cat("\n=== HISTOGRAM COMPLETE ===\n")
