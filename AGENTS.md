# Analysis Context - Experiment 5 (Richard & Omar)

## OpenCode Workflow
- Use **Tab** to switch to **Plan agent** for read-only analysis before making changes
- Use **@explore** for fast codebase searching - faster than manual grep
- Use **@general** for parallel multi-step tasks
- After code changes, run lint/typecheck if available

## Tool Usage
- **Bash**: Run Rscript commands (`Rscript 01_de_analysis.R`)
- **Read/Edit**: File operations
- **grep/glob**: Finding files and code patterns

## Project Structure
```
D:\AG Gomez Analysis/
├── 01_de_analysis.R        # Main DE analysis (corrected sample mapping)
├── 02_venn_diagram.R      # Venn diagram generator
├── 03_volcano_plots.R     # Volcano plot generator
├── 04_heatmap.R           # Heatmap generator (genes across all groups)
├── 05_pathway_analysis.R  # Pathway enrichment (GO/KEGG - UP/DOWN separated)
├── 06_summary_files.R     # Generate summary.txt for each group
├── 08_wiki_pathways.R     # WikiPathways enrichment analysis
├── 09_gene_histogram.R    # Gene distribution histogram (UP/DOWN/COMMON)
├── 09_pdf_report.R        # Generate PDF report with all plots
├── AGENTS.md              # This file
├── README.md              # Project documentation
├── gf_24h/               # GF 24h analysis results
│   ├── rds/              # DESeq2 results and VST counts
│   └── results/          # Plots (venn, volcano, heatmap, pathway)
├── gf_6h/
├── spf_24h/
├── spf_6h/
└── data/                 # Input data (gitignored)
```

## R Code Standards
- Use `<-` for assignment (not `=`)
- Comment with `#` for section headers
- Variable names: lowercase with underscores
- Function names: lowercase with underscores
- Load packages explicitly with `library()` at script top
- Use tidyverse style for dplyr operations
- Include session info at end of each script
- Clean up memory with `rm()` + `gc()` after large operations

## Project Overview
- **Project**: Differential gene expression analysis of epidermal organoids
- **Organism**: Mouse (Mus musculus)
- **Tissue**: Epidermal organoids

## R Dependencies
```r
install.packages(c("BiocManager", "dplyr", "tidyr", "readxl", "ggplot2", "ggVennDiagram", "pheatmap"))
BiocManager::install(c("DESeq2", "EnhancedVolcano", "clusterProfiler", "org.Mm.eg.db"))
```

## Experiment Design

### Mice
- **Germ Free (GF)**: 4 biological replicates
- **Specific Pathogen Free (SPF)**: 2 biological replicates

### Conditions
1. **Untreated** - Control (no bacteria)
2. **SG1B** - S. epidermidis SG1B (commensal-like)
3. **SG1C** - S. epidermidis SG1C (pathogenic-like)

### Time Points
- 6 hours (6h)
- 24 hours (24h)

## Sample Mapping (CRITICAL - CORRECTED!)
- **_1** = Untreated (condition 1)
- **_2** = SG1B (condition 2)
- **_3** = SG1C (condition 3)

Example: GF1_5_24_C1 = GF mouse 1, 24h, replicate C, condition 1 (Untreated)

## Statistical Parameters
- **Tool**: DESeq2
- **Design**: ~ condition
- **padj cutoff**: 0.05 (standard)
- **|log2FC| cutoff**: 1.0 (standard 2-fold change)

## Analysis Groups (4 total)
1. **gf_24h** - Germ Free mice, 24 hours
2. **gf_6h** - Germ Free mice, 6 hours
3. **spf_24h** - Specific Pathogen Free mice, 24 hours
4. **spf_6h** - Specific Pathogen Free mice, 6 hours

## Comparisons (3 per group)
1. `Untreated vs SG1B` - Commensal response
2. `Untreated vs SG1C` - Pathogenic response
3. `SG1B vs SG1C` - Commensal vs Pathogenic difference

## Results Summary

| Group | Comparison | Significant | Up | Down |
|-------|------------|-------------|-----|------|
| GF 24h | U vs B | 1,331 | 681 | 650 |
| GF 24h | U vs C | 1,777 | 992 | 785 |
| GF 24h | SG1B vs SG1C | 0 | 0 | 0 |
| GF 6h | All comparisons | 0 | 0 | 0 |
| SPF 24h | U vs B | 491 | 122 | 369 |
| SPF 24h | U vs C | 305 | 102 | 203 |
| SPF 24h | SG1B vs SG1C | 0 | 0 | 0 |
| SPF 6h | U vs B | 1 | 0 | 1 |
| SPF 6h | U vs C | 0 | 0 | 0 |
| SPF 6h | SG1B vs SG1C | 1 | 0 | 1 |

## Output Structure
```
gf_24h/
├── rds/
│   ├── dds_object.rds
│   ├── vst_counts.rds
│   ├── res_Untreated_vs_SG1B.rds
│   ├── res_Untreated_vs_SG1C.rds
│   └── res_SG1B_vs_SG1C.rds
└── results/
    ├── summary.txt
    ├── venn_diagram.png
    ├── volcano_*.png
    ├── heatmap_*.png
    ├── pathway_go_*_UP.png
    ├── pathway_go_*_DOWN.png
    ├── pathway_kegg_*_UP.png
    ├── pathway_kegg_*_DOWN.png
    └── wiki_*_UP.png / wiki_*_DOWN.png
```

## Scripts Status
- [x] 01_de_analysis.R - DE analysis (corrected mapping)
- [x] 02_venn_diagram.R - Venn diagrams
- [x] 03_volcano_plots.R - Volcano plots
- [x] 04_heatmap.R - Clustered heatmaps
- [x] 05_pathway_analysis.R - GO/KEGG pathway enrichment
- [x] 06_summary_files.R - Summary files
- [x] 08_wiki_pathways.R - WikiPathways analysis
- [x] 09_gene_histogram.R - Gene distribution histograms
- [x] 09_pdf_report.R - PDF report with all plots

## How to Run

```bash
# 1. Run differential expression analysis
Rscript 01_de_analysis.R

# 2. Create Venn diagrams
Rscript 02_venn_diagram.R

# 3. Generate volcano plots
Rscript 03_volcano_plots.R

# 4. Generate heatmaps
Rscript 04_heatmap.R

# 5. Pathway enrichment analysis (GO/KEGG)
Rscript 05_pathway_analysis.R

# 6. Generate summary files
Rscript 06_summary_files.R

# 7. WikiPathways analysis
Rscript 08_wiki_pathways.R

# 8. Gene distribution histograms
Rscript 09_gene_histogram.R

# 9. Generate PDF report
Rscript 09_pdf_report.R
```

## Key Outputs
- **Experiment_5_Report.pdf** - Full PDF report with all plots
- **gene_distribution_histogram.png** - UP/DOWN/COMMON gene distribution
- **gene_counts_summary.csv** - Complete gene counts table
- **gene_counts_by_group.png** - Bar chart of genes by group

## Notes
- Sample mapping: _1=Untreated, _2=SG1B, _3=SG1C
- **Standard thresholds**: padj < 0.05, |log2FC| > 1 (updated from lenient thresholds)
- Upregulated and downregulated genes are analyzed separately in pathway enrichment
