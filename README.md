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
│       ├── venn_diagram.png
│       ├── volcano_Untreated_vs_SG1B.png
│       ├── volcano_Untreated_vs_SG1C.png
│       ├── volcano_SG1B_vs_SG1C.png
│       ├── heatmap_*.png
│       ├── pathway_go_*.png
│       └── pathway_kegg_*.png
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
├── README.md
└── Experiment_5.xlsx
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

| Group | U vs SG1B | U vs SG1C | SG1B vs SG1C |
|-------|------------|------------|---------------|
| GF 24h | 19 | 50 | 10 |
| GF 6h | 6 | 1 | 0 |
| SPF 24h | 0 | 165 | 0 |
| SPF 6h | 22 | 152 | 2 |

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

# 5. Pathway enrichment analysis (GO/KEGG)
Rscript 05_pathway_analysis.R
```

---

## Pathway Analysis

- Generates GO and KEGG enrichment plots for each comparison
- Only runs for comparisons with ≥3 significant genes
- Output: `pathway_go_<comparison>.png` and `pathway_kegg_<comparison>.png`
- Titles include group and comparison name (e.g., "GO Enrichment - SPF 24h - Untreated vs SG1C")

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
