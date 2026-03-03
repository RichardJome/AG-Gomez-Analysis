# 06_summary_files.R ------------------------------------------------------------------------
# Summary File Generator for Experiment 5
# Author: Richard Jome
# Date: 2025
# Purpose: Generate summary.txt files for each group with statistics and file listings

# Load Libraries -------------------------------------------------------------------
library(DESeq2)

# Define Groups --------------------------------------------------------------------
groups <- list(
  "gf_24h" = "Germ Free 24h",
  "gf_6h" = "Germ Free 6h",
  "spf_24h" = "SPF 24h",
  "spf_6h" = "SPF 6h"
)

comparisons <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")

# Function: Generate Summary for Group ---------------------------------------------
#' Generate summary.txt file for a group
#' @param group_name Directory name for the group
#' @param group_label Display label for the group
#' @return Data frame of statistics
generate_summary <- function(group_name, group_label) {
  
  cat(sprintf("\n=== Generating summary: %s ===\n", group_label))
  
  results_dir <- file.path(group_name, "results")
  rds_dir <- file.path(group_name, "rds")
  
  dir.create(results_dir, showWarnings = FALSE)
  
  # Build statistics section
  summary_lines <- c(
    paste("===", toupper(group_label), "SUMMARY ==="),
    "",
    "STATISTICS (padj < 0.2, |log2FC| > 0.5):",
    "-----------------------------------"
  )
  
  stats_data <- data.frame(
    comparison = character(),
    significant = integer(),
    upregulated = integer(),
    downregulated = integer(),
    stringsAsFactors = FALSE
  )
  
  for (comp in comparisons) {
    res_file <- file.path(rds_dir, paste0("res_", comp, ".rds"))
    
    if (!file.exists(res_file)) {
      summary_lines <- c(summary_lines, sprintf("%s: No results file", comp))
      next
    }
    
    res <- readRDS(res_file)
    res_df <- as.data.frame(res)
    
    sig_all <- res_df[res_df$padj < 0.2 & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= 0.5, ]
    sig_up <- sig_all[sig_all$log2FoldChange > 0, ]
    sig_down <- sig_all[sig_all$log2FoldChange < 0, ]
    
    n_sig <- nrow(sig_all)
    n_up <- nrow(sig_up)
    n_down <- nrow(sig_down)
    
    summary_lines <- c(summary_lines, sprintf(
      "%s: %d sig (%d up, %d down)",
      comp, n_sig, n_up, n_down
    ))
    
    stats_data <- rbind(stats_data, data.frame(
      comparison = comp,
      significant = n_sig,
      upregulated = n_up,
      downregulated = n_down,
      stringsAsFactors = FALSE
    ))
  }
  
  # List generated files
  summary_lines <- c(summary_lines, "", "GENERATED FILES:", "---------------")
  
  files <- list.files(results_dir, pattern = "\\.(png|txt)$")
  if (length(files) == 0) {
    summary_lines <- c(summary_lines, "(No output files yet)")
  } else {
    summary_files <- sort(files)
    for (f in summary_files) {
      summary_lines <- c(summary_lines, f)
    }
  }
  
  # Pathway analysis status
  summary_lines <- c(summary_lines, "", "PATHWAY ANALYSIS STATUS:", "----------------------")
  
  for (comp in comparisons) {
    res_file <- file.path(rds_dir, paste0("res_", comp, ".rds"))
    if (!file.exists(res_file)) {
      next
    }
    
    res <- readRDS(res_file)
    res_df <- as.data.frame(res)
    sig_all <- res_df[res_df$padj < 0.2 & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= 0.5, ]
    sig_up <- sig_all[sig_all$log2FoldChange > 0, ]
    sig_down <- sig_all[sig_all$log2FoldChange < 0, ]
    
    n_total <- nrow(sig_all)
    n_up <- nrow(sig_up)
    n_down <- nrow(sig_down)
    
    status_parts <- c()
    if (n_up >= 3) status_parts <- c(status_parts, "UP")
    if (n_down >= 3) status_parts <- c(status_parts, "DOWN")
    
    if (length(status_parts) > 0) {
      summary_lines <- c(summary_lines, sprintf("  %s: %s", comp, paste(status_parts, collapse = ", ")))
    } else {
      summary_lines <- c(summary_lines, sprintf("  %s: <3 genes, skipped", comp))
    }
  }
  
  summary_text <- paste(summary_lines, collapse = "\n")
  cat(summary_text)
  cat("\n")
  
  summary_file <- file.path(results_dir, "summary.txt")
  writeLines(summary_text, summary_file)
  cat(sprintf("Summary saved to: %s\n", summary_file))
  
  return(invisible(stats_data))
}

# Main: Generate All Summaries -----------------------------------------------------
cat("=== GENERATING SUMMARY FILES ===\n")

all_stats <- data.frame()

for (group_name in names(groups)) {
  stats <- generate_summary(group_name, groups[[group_name]])
  if (!is.null(stats)) {
    stats$group <- group_name
    all_stats <- rbind(all_stats, stats)
  }
}

cat("\n=== ALL STATISTICS COMBINED ===\n")
print(all_stats)

cat("\n=== COMPLETE ===\n")

# Save Session Info ---------------------------------------------------------------
sink("session_info_06.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_06.txt\n")
