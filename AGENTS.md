# Analysis Context - Experiment 5 (Richard)

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
├── 01_de_analysis.R        # Main DE analysis script
├── 02_venn_diagram.R       # Venn diagram generator
├── 03_volcano_plots.R      # Volcano plot generator
├── 04_heatmap.R            # Heatmap generator (genes across all groups)
├── 05_pathway_analysis.R   # Pathway enrichment (UP/DOWN separated)
├── 06_summary_files.R      # Generate summary.txt for each group
├── 07_pdf_report.R        # Generate PDF report with all plots
├── AGENTS.md               # This file
├── gf_24h/                 # GF 24h analysis results
│   ├── rds/                # DESeq2 results and VST counts
│   └── results/            # Plots (venn, volcano, heatmap)
├── gf_6h/
├── spf_24h/
├── spf_6h/
└── data/                   # Input data (gitignored)
```

## R Code Standards
- Use `<-` for assignment (not `=`)
- Comment with `#` for section headers
- Variable names: lowercase with underscores (e.g., `result_df`)
- Function names: lowercase with underscores (e.g., `read_count_data`)
- Load packages explicitly with `library()` at script top
- Use tidyverse style for dplyr operations

## Project Overview
- **Project**: Differential gene expression analysis of epidermal organoids
- **Organism**: Mouse (Mus musculus)
- **Tissue**: Epidermal organoids

## R Dependencies

All packages available via Bioconductor or CRAN:

```r
install.packages(c("BiocManager", "dplyr", "tidyr", "readxl", "ggplot2", "ggVennDiagram"))
BiocManager::install(c("DESeq2", "EnhancedVolcano", "clusterProfiler", "org.Mm.eg.db"))
```

| Package | Purpose |
|---------|---------|
| DESeq2 | Differential expression analysis |
| dplyr, tidyr | Data manipulation |
| readxl | Read Excel files |
| ggplot2 | Base plotting |
| ggVennDiagram | Venn diagrams |
| EnhancedVolcano | Volcano plots |
| clusterProfiler | Pathway enrichment (GO/KEGG) |
| org.Mm.eg.db | Mouse gene annotation |

## How to Run

```bash
# 1. Run differential expression analysis
Rscript 01_de_analysis.R

# 2. Create Venn diagrams
Rscript 02_venn_diagram.R

# 3. Generate volcano plots
Rscript 03_volcano_plots.R
```

Outputs are saved to respective group folders (gf_24h/, gf_6h/, spf_24h/, spf_6h/).

## What NOT to Modify

- Do not change sample column names in 01_de_analysis.R without updating sample mapping
- Do not commit large data files (.xlsx, .rds, .png) - they are gitignored
- Do not modify DESeq2 design formula without re-running all analyses

## Experiment Design

### Mice
- **Germ Free (GF)**: 4 biological replicates
- **Specific Pathogen Free (SPF)**: 2 biological replicates

### Conditions
1. **Untreated** - Control (no bacteria)
2. **SG1B** - S. epidermidis SG1B (commensal-like, beneficial)
3. **SG1C** - S. epidermidis SG1C (pathogenic-like, harmful)

### Time Points
- 6 hours (6h)
- 24 hours (24h)

## Sample Mapping

| Group | Sample Columns | Mice | Replicates |
|-------|---------------|------|-------------|
| GF 24h | GF*_5_24_* | GF1, GF2, GF3, GF4 | 4 |
| GF 6h | GF*_5_6h_* | GF1, GF2, GF3, GF4 | 4 |
| SPF 24h | SPF*_5_24_* | SPF1, SPF2 | 2 |
| SPF 6h | SPF*_5_6h_* | SPF1, SPF2 | 2 |

## Analysis Groups (4 total)
1. **gf_24h** - Germ Free mice, 24 hours
2. **gf_6h** - Germ Free mice, 6 hours
3. **spf_24h** - Specific Pathogen Free mice, 24 hours
4. **spf_6h** - Specific Pathogen Free mice, 6 hours

## Comparisons (3 per group)

For each group, three DESeq2 contrasts:
1. `Untreated vs SG1B` - Commensal response
2. `Untreated vs SG1C` - Pathogenic response
3. `SG1B vs SG1C` - Commensal vs Pathogenic difference

## Statistical Parameters
- **Tool**: DESeq2
- **Design**: ~ condition
- **padj cutoff**: 0.05
- **|log2FC| cutoff**: 1.0

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
    ├── summary.txt                    # Statistics and file listing
    ├── venn_diagram.png
    ├── volcano_Untreated_vs_SG1B.png
    ├── volcano_Untreated_vs_SG1C.png
    ├── volcano_SG1B_vs_SG1C.png
    ├── heatmap_*.png (varies by available genes)
    ├── pathway_go_*_UP.png           # GO enrichment - Upregulated
    ├── pathway_go_*_DOWN.png         # GO enrichment - Downregulated
    ├── pathway_kegg_*_UP.png        # KEGG enrichment - Upregulated
    └── pathway_kegg_*_DOWN.png      # KEGG enrichment - Downregulated

gf_6h/   (same structure)
spf_24h/ (same structure)
spf_6h/  (same structure)
```

Note: 
- Pathway plots are only generated for comparisons with ≥3 significant genes
- Pathway analysis separates UP and DOWN genes (best practice)
- Each group folder has a summary.txt with statistics and file listing

## Key Biological Findings (from initial GF 24h analysis)

### Barrier Genes (DOWNREGULATED in treatment - LEFT side of volcano)
- **Flg** (Filaggrin) - Skin barrier protein
- **Dsg1a** (Desmoglein-1) - Epidermal adhesion
- **Klk5** (Kallikrein 5) - Protease for skin desquamation
- **Aqp5** (Aquaporin 5) - Water channel
- **Svep1** - Lymphatic development
- **Nell1** - Neural EGFL-like protein

### Interpretation
- Bacterial treatment downregulates skin barrier genes
- Both commensal and pathogenic strains trigger barrier remodeling
- Suggests skin is actively responding to bacterial colonization

## Scripts Status
- [x] 01_de_analysis.R - Runs DE analysis for all 4 groups
- [x] 02_venn_diagram.R - Creates Venn diagrams for all 4 groups
- [x] 03_volcano_plots.R - Generates volcano plots for all 4 groups
- [x] 04_heatmap.R - Generates clustered heatmaps comparing all groups
- [x] 05_pathway_analysis.R - Pathway enrichment (UP/DOWN separated)
- [x] 06_summary_files.R - Generates summary.txt for each group
- [x] 07_pdf_report.R - Generates PDF report with all plots

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

# 5. Pathway enrichment analysis (UP/DOWN separated)
Rscript 05_pathway_analysis.R

# 6. Generate summary files
Rscript 06_summary_files.R

# 7. Generate PDF report
Rscript 07_pdf_report.R
```

## Heatmap Output
- Heatmaps saved to each group's `results/` folder
- For each comparison, shows significant genes across all 4 groups
- Dynamic sizing: shows available genes (skip if <2 genes)
- Z-score row scaling for relative expression visualization
- Dendrograms on both rows and columns for clustering

## Notes for Omar
- Volcano plot interpretation: LEFT = higher in Untreated, RIGHT = higher in treatment
- All significant genes are labeled (sorted by p-value, most significant first)
- X-axis fixed at -10 to +10 for consistency
- Color scheme: grey → green → blue → red (by significance)
