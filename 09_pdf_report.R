# 07_pdf_report.R ------------------------------------------------------------------------
# PDF Report Generator for Experiment 5
# Author: Richard Jome
# Date: 2025
# Purpose: Generate comprehensive PDF report with all plots and summaries

# Load Libraries -------------------------------------------------------------------
library(png)
library(grid)
library(ggplot2)
library(ggVennDiagram)

# Parameters (standard thresholds)
fdr_threshold <- 0.05   # padj < 0.05 (standard)
log2fc_threshold <- 1.0 # |log2FC| > 1 (standard 2-fold change)

# Define Groups (order for report) ------------------------------------------------
groups_order <- list(
  "gf_6h" = "Germ Free 6h",
  "gf_24h" = "Germ Free 24h",
  "spf_6h" = "SPF 6h",
  "spf_24h" = "SPF 24h"
)

comparisons <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")
groups_order <- list(
  "gf_6h" = "Germ Free 6h",
  "gf_24h" = "Germ Free 24h",
  "spf_6h" = "SPF 6h",
  "spf_24h" = "SPF 24h"
)

comparisons <- c("Untreated_vs_SG1B", "Untreated_vs_SG1C", "SG1B_vs_SG1C")

# Function: Get Statistics --------------------------------------------------------
#' Extract statistics from DESeq2 results
#' @param group_name Directory name for the group
#' @return Data frame of statistics
get_stats <- function(group_name) {
  rds_dir <- file.path(group_name, "rds")
  stats_df <- data.frame()
  
  for (comp in comparisons) {
    res_file <- file.path(rds_dir, paste0("res_", comp, ".rds"))
    if (!file.exists(res_file)) next
    
    res <- readRDS(res_file)
    res_df <- as.data.frame(res)
    
    sig_all <- res_df[res_df$padj < fdr_threshold & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= log2fc_threshold, ]
    sig_up <- sig_all[sig_all$log2FoldChange > 0, ]
    sig_down <- sig_all[sig_all$log2FoldChange < 0, ]
    
    stats_df <- rbind(stats_df, data.frame(
      Comparison = comp,
      Significant = nrow(sig_all),
      Upregulated = nrow(sig_up),
      Downregulated = nrow(sig_down)
    ))
  }
  return(stats_df)
}

# Page Counter ---------------------------------------------------------------------
page_counter <- 0

# Function: Add PNG Page ----------------------------------------------------------
#' Add a PNG image to the PDF on its own page
#' @param img_path Path to PNG file
#' @param title Optional title for the page
#' @return NULL (adds page to PDF)
add_png_page <- function(img_path, title = "") {
  page_counter <<- page_counter + 1
  
  if (!file.exists(img_path)) {
    grid.newpage()
    grid.text("File not found", y = 0.5, 
              gp = gpar(fontsize = 16))
    grid.text(paste("Page", page_counter), y = 0.03, 
              gp = gpar(fontsize = 10, col = "gray50"))
    return()
  }
  
  img <- readPNG(img_path)
  
  grid.newpage()
  
  # Title at top
  if (title != "") {
    grid.text(title, y = 0.97, 
              gp = gpar(fontsize = 14, fontface = "bold"))
  }
  
  # Figure centered
  grid.raster(img, y = 0.50, height = 0.84, interpolate = FALSE)
  
  # Footer with page number
  grid.text(paste("Page", page_counter), y = 0.03, 
            gp = gpar(fontsize = 10, col = "gray50"))
}

# Function: Add Title Page --------------------------------------------------------
add_title_page <- function() {
  page_counter <<- page_counter + 1
  
  grid.newpage()
  grid.text("Experiment 5 Analysis", y = 0.70, 
            gp = gpar(fontsize = 28, fontface = "bold"))
  grid.text("Differential Gene Expression Analysis", y = 0.55, 
            gp = gpar(fontsize = 16))
  grid.text("Performed by Richard Jome", y = 0.45, 
            gp = gpar(fontsize = 14))
  grid.text(paste("Date:", Sys.Date()), y = 0.35, 
            gp = gpar(fontsize = 12, col = "gray50"))
  grid.text(paste("Page", page_counter), y = 0.03, 
            gp = gpar(fontsize = 10, col = "gray50"))
}

# Function: Add Summary Page -------------------------------------------------------
add_summary_page <- function() {
  page_counter <<- page_counter + 1
  
  grid.newpage()
  
  # Title
  grid.text("Summary Statistics", y = 0.97, 
            gp = gpar(fontsize = 18, fontface = "bold"))
  
  # Table header
  grid.text("Group", x = 0.08, y = 0.88, 
            gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("Comparison", x = 0.32, y = 0.88, 
            gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("Sig", x = 0.58, y = 0.88, 
            gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("Up", x = 0.68, y = 0.88, 
            gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("Down", x = 0.80, y = 0.88, 
            gp = gpar(fontsize = 12, fontface = "bold"))
  
  # Horizontal line
  grid.lines(c(0.05, 0.95), c(0.85, 0.85), 
             gp = gpar(col = "black", lwd = 1))
  
  # Get all stats
  all_stats <- data.frame()
  for (grp in names(groups_order)) {
    stats <- get_stats(grp)
    stats$Group <- groups_order[[grp]]
    all_stats <- rbind(all_stats, stats)
  }
  
  all_stats <- all_stats[, c("Group", "Comparison", "Significant", "Upregulated", "Downregulated")]
  
  # Table rows
  y_pos <- 0.80
  for (i in 1:nrow(all_stats)) {
    grid.text(all_stats$Group[i], x = 0.08, y = y_pos, 
              gp = gpar(fontsize = 10), just = "left")
    grid.text(all_stats$Comparison[i], x = 0.32, y = y_pos, 
              gp = gpar(fontsize = 10), just = "left")
    grid.text(as.character(all_stats$Significant[i]), x = 0.58, y = y_pos, 
              gp = gpar(fontsize = 10), just = "center")
    grid.text(as.character(all_stats$Upregulated[i]), x = 0.68, y = y_pos, 
              gp = gpar(fontsize = 10), just = "center")
    grid.text(as.character(all_stats$Downregulated[i]), x = 0.80, y = y_pos, 
              gp = gpar(fontsize = 10), just = "center")
    
    y_pos <- y_pos - 0.045
    
    if (i %% 3 == 0) {
      y_pos <- y_pos - 0.01
    }
  }
  
  # Parameters note
  grid.text(sprintf("Parameters: padj < %.2f, |log2FC| > %.1f", fdr_threshold, log2fc_threshold), y = 0.08, 
            gp = gpar(fontsize = 10, col = "gray50"))
  
  # Page number
  grid.text(paste("Page", page_counter), y = 0.03, 
            gp = gpar(fontsize = 10, col = "gray50"))
}

# Function: Add Group Section -----------------------------------------------------
#' Add all plots for a group to the PDF
#' @param group_name Directory name for the group
#' @param group_label Display label for the group
#' @return NULL (adds pages to PDF)
add_group_section <- function(group_name, group_label) {
  cat(sprintf("Processing: %s\n", group_label))
  
  results_dir <- file.path(group_name, "results")
  
  # Section title page
  page_counter <<- page_counter + 1
  grid.newpage()
  grid.text(group_label, y = 0.50, 
            gp = gpar(fontsize = 28, fontface = "bold"))
  grid.text(paste("Page", page_counter), y = 0.03, 
            gp = gpar(fontsize = 10, col = "gray50"))
  
  # Volcano plots
  for (comp in comparisons) {
    res_file <- file.path(group_name, "rds", paste0("res_", comp, ".rds"))
    volcano_file <- file.path(results_dir, paste0("volcano_", comp, ".png"))
    
    if (file.exists(volcano_file)) {
      title_str <- paste(group_label, "- Volcano:", gsub("_", " vs ", comp))
      add_png_page(volcano_file, title_str)
    } else if (file.exists(res_file)) {
      res <- readRDS(res_file)
      res_df <- as.data.frame(res)
      sig_count <- sum(res_df$padj < fdr_threshold & !is.na(res_df$padj) & abs(res_df$log2FoldChange) >= log2fc_threshold)
      
      page_counter <<- page_counter + 1
      grid.newpage()
      
      title_str <- paste(group_label, "- Volcano:", gsub("_", " vs ", comp), 
                        "(n=", sig_count, ")")
      grid.text(title_str, y = 0.97, 
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      par(mar = c(4, 4, 2, 1))
      with(res_df, plot(log2FoldChange, -log10(padj), 
                       main = title_str,
                       pch = 16, cex = 0.5,
                       col = ifelse(padj < fdr_threshold & abs(log2FoldChange) >= log2fc_threshold, 
                                    ifelse(log2FoldChange > 0, "red", "blue"), "gray80"),
                       xlab = "log2FoldChange", ylab = "-log10(padj)"))
      abline(h = -log10(fdr_threshold), col = "gray40", lty = 2)
      abline(v = c(-log2fc_threshold, log2fc_threshold), col = "gray40", lty = 2)
      
      grid.text(paste("Page", page_counter), y = 0.03, 
                gp = gpar(fontsize = 10, col = "gray50"))
    } else {
      page_counter <<- page_counter + 1
      grid.newpage()
      grid.text(paste("No data for", comp), y = 0.50, 
                gp = gpar(fontsize = 14))
      grid.text(paste("Page", page_counter), y = 0.03, 
                gp = gpar(fontsize = 10, col = "gray50"))
    }
  }
  
  # Venn diagram
  venn_file <- file.path(results_dir, "venn_diagram.png")
  if (file.exists(venn_file)) {
    add_png_page(venn_file, paste(group_label, "- Venn Diagram"))
  }
  
  # Heatmaps
  heatmap_files <- list.files(results_dir, pattern = "heatmap.*\\.png$", full.names = TRUE)
  heatmap_files <- sort(heatmap_files)
  
  for (hfile in heatmap_files) {
    fname <- basename(hfile)
    title_str <- paste(group_label, "-", gsub("_", " ", fname))
    add_png_page(hfile, title_str)
  }
  
  # Pathway plots (GO/KEGG)
  pathway_files <- list.files(results_dir, pattern = "pathway_.*\\.png$", full.names = TRUE)
  pathway_files <- sort(pathway_files)
  
  for (pfile in pathway_files) {
    fname <- basename(pfile)
    fname_clean <- gsub("pathway_", "", gsub("\\.png$", "", fname))
    fname_clean <- gsub("_", " ", fname_clean)
    title_str <- paste(group_label, "- GO/KEGG:", fname_clean)
    add_png_page(pfile, title_str)
  }
  
  # WikiPathways plots
  wiki_files <- list.files(results_dir, pattern = "wiki_.*\\.png$", full.names = TRUE)
  wiki_files <- sort(wiki_files)
  
  for (wfile in wiki_files) {
    fname <- basename(wfile)
    fname_clean <- gsub("wiki_", "", gsub("\\.png$", "", fname))
    fname_clean <- gsub("_", " ", fname_clean)
    title_str <- paste(group_label, "- WikiPathways:", fname_clean)
    add_png_page(wfile, title_str)
  }
}

# Main: Generate PDF ---------------------------------------------------------------
cat("=== GENERATING PDF REPORT ===\n\n")

pdf("Experiment_5_Report.pdf", width = 11, height = 8.5, onefile = TRUE)

add_title_page()
add_summary_page()

for (grp in names(groups_order)) {
  add_group_section(grp, groups_order[[grp]])
}

# Add histogram and summary plots at the end
page_counter <<- page_counter + 1
grid.newpage()
grid.text("Gene Distribution Analysis", y = 0.97, 
          gp = gpar(fontsize = 18, fontface = "bold"))

# Gene distribution histogram
hist_file <- "gene_distribution_histogram.png"
if (file.exists(hist_file)) {
  add_png_page(hist_file, "Gene Distribution: UP/DOWN/COMMON")
}

# Gene counts by group
counts_file <- "gene_counts_by_group.png"
if (file.exists(counts_file)) {
  add_png_page(counts_file, "Gene Counts by Group and Comparison")
}

dev.off()

cat("\n=== PDF REPORT GENERATED ===\n")
cat("Output: Experiment_5_Report.pdf\n")
cat(sprintf("Total pages: %d\n", page_counter))

# Save Session Info ---------------------------------------------------------------
sink("session_info_07.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_07.txt\n")
