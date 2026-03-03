# 01_de_analysis.R ------------------------------------------------------------------------
# Differential Expression Analysis for Experiment 5
# Corrected sample mapping based on Word file:
#   _1 = Untreated
#   _2 = SG1B (commensal-like)
#   _3 = SG1C (pathogenic-like)
# Filter: -log10(padj) > 0.7 (padj < 0.2), |log2FC| > 0.5

library(DESeq2)
library(dplyr)
library(readxl)
library(tidyr)

# Read Data -------------------------------------------------------------------------
df <- read_excel("Experiment_5.xlsx")

# Parameters ------------------------------------------------------------------------
params <- list(
  fdr_threshold = 0.2,        # -log10(padj) > 0.7 = padj < 0.2
  log2fc_threshold = 0.5      # |log2FC| > 0.5
)

cat("=== DE Analysis Parameters ===\n")
cat(sprintf("FDR threshold: %.2f (-log10(padj) > 0.7)\n", params$fdr_threshold))
cat(sprintf("|log2FC| threshold: %.1f\n", params$log2fc_threshold))
cat("\n")

# Function to get sample columns for a group ---------------------------------------
get_sample_cols <- function(df, time, mouse_type) {
  all_cols <- colnames(df)[-c(1:7)]  # Exclude metadata columns
  
  # Filter by time (6h or 24h) and mouse type (GF or SPF)
  if (time == "6h") {
    cols <- all_cols[grepl("6h", all_cols)]
  } else if (time == "24h") {
    cols <- all_cols[grepl("24", all_cols)]
  }
  
  if (mouse_type == "GF") {
    cols <- cols[grepl("GF", cols)]
  } else if (mouse_type == "SPF") {
    cols <- cols[grepl("SPF", cols)]
  }
  
  return(cols)
}

# Function to categorize sample by condition -------------------------------------
get_condition <- function(col_name) {
  if (grepl("_[^_]*1$", col_name)) return("Untreated")  # Ends with 1 (e.g., _C1, _A1)
  if (grepl("_[^_]*2$", col_name)) return("SG1B")       # Ends with 2
  if (grepl("_[^_]*3$", col_name)) return("SG1C")       # Ends with 3
  return(NA)
}

# Function: Run DE Analysis -----------------------------------------------------
run_de_analysis <- function(group_name, sample_cols, output_dir, params) {
  
  cat(sprintf("\n=== Running DE analysis for: %s ===\n", group_name))
  
  # Get relevant columns
  counts <- df %>% select(all_of(c("gene_id", "gene_name", sample_cols)))
  
  # Categorize samples by condition
  conditions <- sapply(sample_cols, get_condition)
  n_untreated <- sum(conditions == "Untreated")
  n_sg1b <- sum(conditions == "SG1B")
  n_sg1c <- sum(conditions == "SG1C")
  n_replicates <- min(n_untreated, n_sg1b, n_sg1c)
  
  # Create coldata
  coldata <- data.frame(
    sample = sample_cols,
    condition = conditions,
    replicate = rep(1:n_replicates, each = 3, length.out = length(sample_cols))
  )
  coldata$condition <- factor(coldata$condition, levels = c("Untreated", "SG1B", "SG1C"))
  
  cat(sprintf("Conditions: Untreated=%d, SG1B=%d, SG1C=%d\n", 
              sum(coldata$condition == "Untreated"),
              sum(coldata$condition == "SG1B"),
              sum(coldata$condition == "SG1C")))
  
  # Create count matrix
  count_matrix <- counts %>% select(-gene_id, -gene_name) %>% as.matrix()
  rownames(count_matrix) <- counts$gene_id
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
  
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  cat(sprintf("After filtering: %d genes retained\n", nrow(dds)))
  
  dds <- DESeq(dds)
  
  # Get contrasts
  res_U_vs_B <- results(dds, contrast = c("condition", "SG1B", "Untreated"))
  res_U_vs_B$gene <- counts$gene_name[match(rownames(res_U_vs_B), counts$gene_id)]
  
  res_U_vs_C <- results(dds, contrast = c("condition", "SG1C", "Untreated"))
  res_U_vs_C$gene <- counts$gene_name[match(rownames(res_U_vs_C), counts$gene_id)]
  
  res_B_vs_C <- results(dds, contrast = c("condition", "SG1B", "SG1C"))
  res_B_vs_C$gene <- counts$gene_name[match(rownames(res_B_vs_C), counts$gene_id)]
  
  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  
  # Save DESeq2 results
  saveRDS(res_U_vs_B, file.path(output_dir, "rds", "res_Untreated_vs_SG1B.rds"))
  saveRDS(res_U_vs_C, file.path(output_dir, "rds", "res_Untreated_vs_SG1C.rds"))
  saveRDS(res_B_vs_C, file.path(output_dir, "rds", "res_SG1B_vs_SG1C.rds"))
  saveRDS(dds, file.path(output_dir, "rds", "dds_object.rds"))
  
  # VST normalized counts
  vst_counts <- assay(vst(dds, blind = TRUE))
  vst_df <- as.data.frame(vst_counts)
  vst_df$gene_id <- rownames(vst_counts)
  vst_df$gene_name <- counts$gene_name[match(vst_df$gene_id, counts$gene_id)]
  saveRDS(vst_df, file.path(output_dir, "rds", "vst_counts.rds"))
  
  # Summary statistics with new thresholds
  fdr <- params$fdr_threshold
  lfc <- params$log2fc_threshold
  
  sig1 <- subset(as.data.frame(res_U_vs_B), padj < fdr & abs(log2FoldChange) > lfc)
  sig2 <- subset(as.data.frame(res_U_vs_C), padj < fdr & abs(log2FoldChange) > lfc)
  sig3 <- subset(as.data.frame(res_B_vs_C), padj < fdr & abs(log2FoldChange) > lfc)
  
  cat(sprintf("Untreated vs SG1B: %d significant genes (%d up, %d down)\n", 
              nrow(sig1), sum(sig1$log2FoldChange > 0), sum(sig1$log2FoldChange < 0)))
  cat(sprintf("Untreated vs SG1C: %d significant genes (%d up, %d down)\n", 
              nrow(sig2), sum(sig2$log2FoldChange > 0), sum(sig2$log2FoldChange < 0)))
  cat(sprintf("SG1B vs SG1C: %d significant genes (%d up, %d down)\n", 
              nrow(sig3), sum(sig3$log2FoldChange > 0), sum(sig3$log2FoldChange < 0)))
  
  rm(count_matrix, sig1, sig2, sig3)
  gc()
  
  return(invisible(NULL))
}

# Define groups ------------------------------------------------------------------
# GF 24h
gf_24h_cols <- get_sample_cols(df, "24h", "GF")
# GF 6h
gf_6h_cols <- get_sample_cols(df, "6h", "GF")
# SPF 24h
spf_24h_cols <- get_sample_cols(df, "24h", "SPF")
# SPF 6h
spf_6h_cols <- get_sample_cols(df, "6h", "SPF")

cat("Sample columns identified:\n")
cat(sprintf("GF 24h: %d samples\n", length(gf_24h_cols)))
cat(sprintf("GF 6h: %d samples\n", length(gf_6h_cols)))
cat(sprintf("SPF 24h: %d samples\n", length(spf_24h_cols)))
cat(sprintf("SPF 6h: %d samples\n", length(spf_6h_cols)))

# Run Analysis --------------------------------------------------------------------
run_de_analysis("GF 24h", gf_24h_cols, "gf_24h", params)
run_de_analysis("GF 6h", gf_6h_cols, "gf_6h", params)
run_de_analysis("SPF 24h", spf_24h_cols, "spf_24h", params)
run_de_analysis("SPF 6h", spf_6h_cols, "spf_6h", params)

cat("\n=== All DE analyses complete! ===\n")

# Save Session Info ---------------------------------------------------------------
sink("session_info_01.txt")
sessionInfo()
sink()

cat("Session info saved to session_info_01.txt\n")
