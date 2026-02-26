# AG Gomez Analysis - Differential Gene Expression

## Project Overview
Differential gene expression analysis of epidermal organoids from mice treated with different Staphylococcus epidermidis strains.

## Experiment: Experiment 5 (Omar)

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

### Total Significant Genes (padj < 0.05, |log2FC| > 1)

| Group | U vs SG1B | U vs SG1C | SG1B vs SG1C |
|-------|------------|------------|---------------|
| GF 24h | 19 | 50 | 10 |
| GF 6h | 6 | 1 | 0 |
| SPF 24h | 0 | 165 | 0 |
| SPF 6h | 22 | 152 | 2 |

### Upregulated / Downregulated Breakdown

| Group | Comparison | Significant | Upregulated | Downregulated |
|-------|------------|-------------|-------------|----------------|
| GF 24h | U vs SG1B | 19 | 8 | 11 |
| GF 24h | U vs SG1C | 50 | 22 | 28 |
| GF 24h | SG1B vs SG1C | 10 | 0 | 10 |
| GF 6h | U vs SG1B | 6 | 3 | 3 |
| GF 6h | U vs SG1C | 1 | 1 | 0 |
| GF 6h | SG1B vs SG1C | 0 | 0 | 0 |
| SPF 24h | U vs SG1B | 0 | 0 | 0 |
| SPF 24h | U vs SG1C | 165 | 102 | 63 |
| SPF 24h | SG1B vs SG1C | 0 | 0 | 0 |
| SPF 6h | U vs SG1B | 22 | 22 | 0 |
| SPF 6h | U vs SG1C | 152 | 142 | 10 |
| SPF 6h | SG1B vs SG1C | 2 | 0 | 2 |

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

## Pathway Analysis

- Generates GO and KEGG enrichment plots for each comparison
- **Separates upregulated and downregulated genes**:
  - `pathway_go_<comparison>_UP.png` - Upregulated genes only
  - `pathway_go_<comparison>_DOWN.png` - Downregulated genes only
- Only runs for comparisons with ≥3 significant genes in each direction
- Each group folder has a `summary.txt` with statistics and file listing

---

## Key Biological Findings

### Barrier Genes (Downregulated by bacteria - LEFT side of volcano)
- **Flg** (Filaggrin) - Skin barrier protein
- **Dsg1a** (Desmoglein-1) - Epidermal adhesion
- **Klk5** (Kallikrein 5) - Protease for skin desquamation
- **Aqp5** (Aquaporin 5) - Water channel

### Interpretation
- Bacterial treatment downregulates skin barrier genes
- Both commensal and pathogenic strains trigger barrier remodeling
- Suggests skin is actively responding to bacterial colonization
