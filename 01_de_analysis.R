library(DESeq2)
library(dplyr)
library(readxl)
library(tidyr)

df <- read_excel("Experiment_5.xlsx")

run_de_analysis <- function(group_name, sample_cols, output_dir) {
  
  cat(sprintf("\n=== Running DE analysis for: %s ===\n", group_name))
  
  counts <- df %>% select(all_of(c("gene_id", "gene_name", sample_cols)))
  
  n_replicates <- length(sample_cols) / 3
  
  coldata <- data.frame(
    sample = sample_cols,
    condition = rep(c("Untreated", "SG1B", "SG1C"), each = n_replicates),
    replicate = rep(1:n_replicates, times = 3)
  )
  coldata$condition <- factor(coldata$condition, levels = c("Untreated", "SG1B", "SG1C"))
  
  count_matrix <- counts %>% select(-gene_id, -gene_name) %>% as.matrix()
  rownames(count_matrix) <- counts$gene_id
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  
  res_U_vs_B <- results(dds, contrast = c("condition", "SG1B", "Untreated"))
  res_U_vs_B$gene <- counts$gene_name[match(rownames(res_U_vs_B), counts$gene_id)]
  
  res_U_vs_C <- results(dds, contrast = c("condition", "SG1C", "Untreated"))
  res_U_vs_C$gene <- counts$gene_name[match(rownames(res_U_vs_C), counts$gene_id)]
  
  res_B_vs_C <- results(dds, contrast = c("condition", "SG1B", "SG1C"))
  res_B_vs_C$gene <- counts$gene_name[match(rownames(res_B_vs_C), counts$gene_id)]
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(res_U_vs_B, file.path(output_dir, "rds", "res_Untreated_vs_SG1B.rds"))
  saveRDS(res_U_vs_C, file.path(output_dir, "rds", "res_Untreated_vs_SG1C.rds"))
  saveRDS(res_B_vs_C, file.path(output_dir, "rds", "res_SG1B_vs_SG1C.rds"))
  
  sig1 <- subset(as.data.frame(res_U_vs_B), padj < 0.05 & abs(log2FoldChange) > 1)
  sig2 <- subset(as.data.frame(res_U_vs_C), padj < 0.05 & abs(log2FoldChange) > 1)
  sig3 <- subset(as.data.frame(res_B_vs_C), padj < 0.05 & abs(log2FoldChange) > 1)
  
  cat(sprintf("Untreated vs SG1B: %d significant genes\n", nrow(sig1)))
  cat(sprintf("Untreated vs SG1C: %d significant genes\n", nrow(sig2)))
  cat(sprintf("SG1B vs SG1C: %d significant genes\n", nrow(sig3)))
  cat(sprintf("Results saved to: %s/\n", output_dir))
  
  return(invisible(NULL))
}

gf_24h_cols <- c("GF1_5_24_C1", "GF1_5_24_C2", "GF1_5_24_C3",
                  "GF2_5_24_A1", "GF2_5_24_A2", "GF2_5_24_A3",
                  "GF3_5_24_B1", "GF3_5_24_B2", "GF3_5_24_B3",
                  "GF4_5_24_C1", "GF4_5_24_C2", "GF4_5_24_C3")

gf_6h_cols <- c("GF1_5_6h_C1", "GF1_5_6h_C2", "GF1_5_6h_C3",
                 "GF2_5_6h_A1", "GF2_5_6h_A2", "GF2_5_6h_A3",
                 "GF3_5_6h_B1", "GF3_5_6h_B2", "GF3_5_6h_B3",
                 "GF4_5_6h_C1", "GF4_5_6h_C2", "GF4_5_6h_C3")

spf_24h_cols <- c("SPF1_5_24_A1", "SPF1_5_24_A2", "SPF1_5_24_A3",
                   "SPF2_5_24_B1", "SPF2_5_24_B2", "SPF2_5_24_B3")

spf_6h_cols <- c("SPF1_5_6h_A1", "SPF1_5_6h_A2", "SPF1_5_6h_A3",
                  "SPF2_5_6h_B1", "SPF2_5_6h_B2", "SPF2_5_6h_B3")

run_de_analysis("GF 24h", gf_24h_cols, "gf_24h")
run_de_analysis("GF 6h", gf_6h_cols, "gf_6h")
run_de_analysis("SPF 24h", spf_24h_cols, "spf_24h")
run_de_analysis("SPF 6h", spf_6h_cols, "spf_6h")

cat("\n=== All DE analyses complete! ===\n")
