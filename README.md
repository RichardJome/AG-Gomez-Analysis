# AG Gomez Analysis - Differential Gene Expression

## Project Overview
Differential gene expression analysis of epidermal organoids from mice treated with different Staphylococcus epidermidis strains.

## Author
**Richard Jome**

## Experiment: Experiment 5

### Mice
- **Germ Free (GF)**: n=4 biological replicates
- **Specific Pathogen Free (SPF)**: n=2 biological replicates

### Conditions
| Condition | Description |
|-----------|-------------|
| Untreated | Control (no bacteria) |
| SG1B | S. epidermidis SG1B (commensal-like, beneficial) |
| SG1C | S. epidermidis SG1C (pathogenic-like, harmful) |

### Time Points
- 6 hours (6h)
- 24 hours (24h)

---

## Analysis Groups

| Group | Mice | Time | Replicates |
|-------|------|------|-------------|
| GF 24h | 4 | 24h | 4 |
| GF 6h | 4 | 6h | 4 |
| SPF 24h | 2 | 24h | 2 |
| SPF 6h | 2 | 6h | 2 |

---

## Folder Structure

```
AG Gomez Analysis/
├── gf_24h/
│   ├── rds/
│   │   ├── res_Untreated_vs_SG1B.rds
│   │   ├── res_Untreated_vs_SG1C.rds
│   │   └── res_SG1B_vs_SG1C.rds
│   └── results/
│       ├── summary.txt
│       ├── venn_diagram.png
│       ├── volcano_Untreated_vs_SG1B.png
│       ├── volcano_Untreated_vs_SG1C.png
│       ├── volcano_SG1B_vs_SG1C.png
│       ├── heatmap_*.png
│       ├── pathway_go_*_UP.png
│       ├── pathway_go_*_DOWN.png
│       ├── pathway_kegg_*_UP.png
│       └── pathway_kegg_*_DOWN.png
│
├── gf_6h/
│   ├── rds/
│   └── results/
│
├── spf_24h/
│   ├── rds/
│   └── results/
│
├── spf_6h/
│   ├── rds/
│   └── results/
│
├── 01_de_analysis.R
├── 02_venn_diagram.R
├── 03_volcano_plots.R
├── 04_heatmap.R
├── 05_pathway_analysis.R
├── 06_summary_files.R
├── 07_pdf_report.R
├── README.md
├── Experiment_5.xlsx
└── Experiment_5_Report.pdf
```

---

## Comparisons

For each group, three DESeq2 contrasts:

1. **Untreated vs SG1B** - Commensal response
2. **Untreated vs SG1C** - Pathogenic response  
3. **SG1B vs SG1C** - Commensal vs Pathogenic difference

---

## Statistical Methods

- **Tool**: DESeq2
- **Design**: ~ condition
- **Significance**: padj < 0.05, |log2FC| > 1

---

## Results Summary

### Total Significant Genes (padj < 0.05, |log2FC| > 1 - STANDARD THRESHOLDS)

| Group | U vs SG1B | U vs SG1C | SG1B vs SG1C |
|-------|------------|------------|---------------|
| GF 24h | 1,331 | 1,777 | 0 |
| GF 6h | 0 | 0 | 0 |
| SPF 24h | 491 | 305 | 0 |
| SPF 6h | 1 | 0 | 1 |

### Upregulated / Downregulated Breakdown

| Group | Comparison | Significant | Upregulated | Downregulated |
|-------|------------|-------------|-------------|----------------|
| GF 24h | U vs SG1B | 1,331 | 681 | 650 |
| GF 24h | U vs SG1C | 1,777 | 992 | 785 |
| GF 24h | SG1B vs SG1C | 0 | 0 | 0 |
| GF 6h | U vs SG1B | 0 | 0 | 0 |
| GF 6h | U vs SG1C | 0 | 0 | 0 |
| GF 6h | SG1B vs SG1C | 0 | 0 | 0 |
| SPF 24h | U vs SG1B | 491 | 122 | 369 |
| SPF 24h | U vs SG1C | 305 | 102 | 203 |
| SPF 24h | SG1B vs SG1C | 0 | 0 | 0 |
| SPF 6h | U vs SG1B | 1 | 0 | 1 |
| SPF 6h | U vs SG1C | 0 | 0 | 0 |
| SPF 6h | SG1B vs SG1C | 1 | 0 | 1 |

---

## Volcano Plot Interpretation

- **LEFT side**: Genes higher in Untreated (control)
- **RIGHT side**: Genes higher in treatment (SG1B/SG1C)
- **Colors**: Grey → Green → Blue → Red (increasing significance)

---

## Run Analysis

```bash
# 1. Run differential expression
Rscript 01_de_analysis.R

# 2. Create Venn diagrams
Rscript 02_venn_diagram.R

# 3. Generate volcano plots
Rscript 03_volcano_plots.R

# 4. Generate heatmaps
Rscript 04_heatmap.R

# 5. Pathway enrichment analysis (UP/DOWN separated)
Rscript 05_pathway_analysis.R

# 6. Generate summary files
Rscript 06_summary_files.R

# 7. Generate PDF report
Rscript 07_pdf_report.R
```

---

## R Script Best Practices

All R scripts follow these conventions:

- **Headers**: Script name, author, date, purpose
- **Section dividers**: `# Section Name -------------------------------------------------------------------`
- **Function docs**: `@param` and `@return` roxygen-style comments
- **Memory cleanup**: `rm()` + `gc()` after large operations
- **Session info**: Each script saves `session_info_XX.txt` for reproducibility

Example:
```r
# 01_script_name.R ------------------------------------------------------------------------
# Description of script
# Author: Name
# Date: 2025

# Load Libraries -------------------------------------------------------------------
library(dplyr)

# Function: Do Something -----------------------------------------------------------
#' Description of function
#' @param input Description of input parameter
#' @return Description of return value
do_something <- function(input) {
  # code
  rm(obj)
  gc()
}

# Main ----------------------------------------------------------------------------
do_something()

# Save Session Info ---------------------------------------------------------------
sink("session_info_XX.txt")
sessionInfo()
sink()
```

---

## Pathway Analysis

- Generates GO and KEGG enrichment plots for each comparison
- **Separates upregulated and downregulated genes**:
  - `pathway_go_<comparison>_UP.png` - Upregulated genes only
  - `pathway_go_<comparison>_DOWN.png` - Downregulated genes only
- Only runs for comparisons with ≥3 significant genes in each direction
- Each group folder has a `summary.txt` with statistics and file listing
