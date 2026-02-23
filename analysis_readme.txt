================================================================================
              DIFFERENTIAL GENE EXPRESSION ANALYSIS - README
================================================================================

EXPERIMENT: Experiment 5 (Omar)
================================================================================

OBJECTIVE:
Analyze gene expression changes in epidermal organoids from mice treated with
different Staphylococcus epidermidis strains (commensal vs pathogenic).

================================================================================

EXPERIMENTAL DESIGN:
================================================================================

MICE:
- Germ Free (GF) mice: n=4 biological replicates
- Specific Pathogen Free (SPF) mice: n=2 biological replicates

CONDITIONS:
- Untreated: Control group (no bacteria)
- SG1B: S. epidermidis SG1B (commensal-like strain - beneficial)
- SG1C: S. epidermidis SG1C (pathogenic-like strain - harmful)

TIME POINTS:
- 6 hours post-treatment
- 24 hours post-treatment

================================================================================

ANALYSIS GROUPS:
================================================================================

| Group     | Mice | Time | Replicates | Description                    |
|-----------|------|------|------------|--------------------------------|
| GF 24h    | 4    | 24h  | 4          | Germ Free mice at 24 hours     |
| GF 6h     | 4    | 6h   | 4          | Germ Free mice at 6 hours     |
| SPF 24h   | 2    | 24h  | 2          | SPF mice at 24 hours          |
| SPF 6h    | 2    | 6h   | 2          | SPF mice at 6 hours           |

================================================================================

FOLDER STRUCTURE:
================================================================================

personal webpage/
├── gf_24h/
│   ├── rds/
│   │   ├── res_Untreated_vs_SG1B.rds
│   │   ├── res_Untreated_vs_SG1C.rds
│   │   └── res_SG1B_vs_SG1C.rds
│   └── results/
│       ├── venn_diagram.png
│       ├── volcano_Untreated_vs_SG1B.png
│       ├── volcano_Untreated_vs_SG1C.png
│       └── volcano_SG1B_vs_SG1C.png
│
├── gf_6h/
│   ├── rds/
│   │   ├── res_Untreated_vs_SG1B.rds
│   │   ├── res_Untreated_vs_SG1C.rds
│   │   └── res_SG1B_vs_SG1C.rds
│   └── results/
│       ├── venn_diagram.png
│       ├── volcano_Untreated_vs_SG1B.png
│       ├── volcano_Untreated_vs_SG1C.png
│       └── volcano_SG1B_vs_SG1C.png
│
├── spf_24h/
│   ├── rds/
│   │   ├── res_Untreated_vs_SG1B.rds
│   │   ├── res_Untreated_vs_SG1C.rds
│   │   └── res_SG1B_vs_SG1C.rds
│   └── results/
│       ├── venn_diagram.png
│       ├── volcano_Untreated_vs_SG1B.png
│       ├── volcano_Untreated_vs_SG1C.png
│       └── volcano_SG1B_vs_SG1C.png
│
├── spf_6h/
│   ├── rds/
│   │   ├── res_Untreated_vs_SG1B.rds
│   │   ├── res_Untreated_vs_SG1C.rds
│   │   └── res_SG1B_vs_SG1C.rds
│   └── results/
│       ├── venn_diagram.png
│       ├── volcano_Untreated_vs_SG1B.png
│       ├── volcano_Untreated_vs_SG1C.png
│       └── volcano_SG1B_vs_SG1C.png
│
├── 01_de_analysis.R
├── 02_venn_diagram.R
├── 03_volcano_plots.R
└── analysis_readme.txt (this file)

================================================================================

COMPARISONS PER GROUP:
================================================================================

For EACH group (GF 24h, GF 6h, SPF 24h, SPF 6h), three comparisons are made:

1. Untreated vs SG1B
   - Genes responding to COMMENSAL bacteria vs no treatment
   - Identifies host response to beneficial bacteria

2. Untreated vs SG1C
   - Genes responding to PATHOGENIC bacteria vs no treatment
   - Identifies host response to harmful bacteria

3. SG1B vs SG1C
   - DIFFERENCE between commensal vs pathogenic response
   - Identifies genes that distinguish beneficial from harmful bacteria

================================================================================

STATISTICAL METHODS:
================================================================================

- Tool: DESeq2 (R package)
- Design: ~ condition
- Significance criteria:
  * Adjusted p-value (padj) < 0.05
  * |log2FoldChange| > 1 (at least 2-fold change)

================================================================================

RDS FILES CONTAIN:
================================================================================

Each RDS file contains a DESeq2 results object with:
- log2FoldChange: Log2 fold change between conditions
  * Positive = higher in second condition (e.g., SG1B > Untreated)
  * Negative = lower in second condition
- padj: Adjusted p-value (Benjamini-Hochberg FDR)
- gene: Gene name (mapped from gene_id)

================================================================================

HOW TO RE-RUN ANALYSIS:
================================================================================

1. Run differential expression analysis:
   Rscript 01_de_analysis.R

2. Create Venn diagrams:
   Rscript 02_venn_diagram.R

3. Generate volcano plots:
   Rscript 03_volcano_plots.R

================================================================================

OUTPUT FILES:
================================================================================

VENN DIAGRAM:
- Shows overlap of differentially expressed genes between the 3 comparisons
- Circle A: Untreated vs SG1B
- Circle B: Untreated vs SG1C  
- Circle C: SG1B vs SG1C

VOLCANO PLOTS:
- X-axis: log2 fold change
- Y-axis: -log10 adjusted p-value
- Colors indicate significance:
  * Grey: Not significant
  * Green: Significant but |log2FC| < 1
  * Blue: Significant with |log2FC| > 1 (below FC cutoff)
  * Red: Significant with |log2FC| > 1 (above FC cutoff)
- Labels show gene names for significant genes

================================================================================

INTERPRETATION:
================================================================================

LEFT SIDE of volcano plot:
- Genes with negative log2FoldChange
- Higher expression in UNTREATED (control)
- These genes are DOWNREGULATED by bacterial treatment

RIGHT SIDE of volcano plot:
- Genes with positive log2FoldChange
- Higher expression in TREATMENT (SG1B or SG1C)
- These genes are UPREGULATED by bacterial treatment

================================================================================

BIOLOGICAL INTERPRETATION:
================================================================================

Barrier Genes (typically downregulated by treatment):
- Filaggrin (Flg), Desmoglein-1 (Dsg1), KLK5, AQP5
- These are skin differentiation/barrier markers
- Changes indicate skin barrier remodeling in response to bacteria

Immune/Inflammatory Genes (may be upregulated):
- Cytokines, chemokines, antimicrobial peptides
- Defense response to bacterial colonization

================================================================================

CREATED: $(date)
LAST UPDATED: $(date)
================================================================================
