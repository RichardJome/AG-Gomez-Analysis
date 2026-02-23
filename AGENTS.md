# Analysis Context - Experiment 5 (Omar)

## Project Overview
- **Project**: Differential gene expression analysis of epidermal organoids
- **Organism**: Mouse (Mus musculus)
- **Tissue**: Epidermal organoids

## R Dependencies

All packages available via Bioconductor or CRAN:

```r
# Core analysis
install.packages("BiocManager")
BiocManager::install("DESeq2")

# Data manipulation
install.packages(c("dplyr", "tidyr", "readxl"))

# Visualization
install.packages(c("ggplot2", "ggVennDiagram"))
BiocManager::install("EnhancedVolcano")
```

### Package List
| Package | Purpose |
|---------|---------|
| DESeq2 | Differential expression analysis |
| dplyr | Data manipulation |
| tidyr | Data tidying |
| readxl | Read Excel files |
| ggplot2 | Base plotting |
| ggVennDiagram | Venn diagrams |
| EnhancedVolcano | Volcano plots |

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
personal webpage/
├── gf_24h/
│   ├── rds/ (3 RDS files)
│   └── results/ (1 venn + 3 volcano plots)
├── gf_6h/
│   ├── rds/
│   └── results/
├── spf_24h/
│   ├── rds/
│   └── results/
├── spf_6h/
│   ├── rds/
│   └── results/
├── 01_de_analysis.R
├── 02_venn_diagram.R
├── 03_volcano_plots.R
└── analysis_readme.txt
```

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
- [x] AGENTS.md - Updated with instructions and dependencies

## Notes for Omar
- Volcano plot interpretation: LEFT = higher in Untreated, RIGHT = higher in treatment
- All significant genes are labeled (sorted by p-value, most significant first)
- X-axis fixed at -10 to +10 for consistency
- Color scheme: grey → green → blue → red (by significance)
