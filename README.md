# schizoplacenta

A reproducible pipeline for analyzing placenta transcriptome data to investigate sex differences and neurodevelopmental outcomes (cognition, language, and motor skills). This repository accompanies the study “coming soon...” (available soon).

Project structure
-----------------
.
├── preprocessing
│   ├── 1_run_qc.sh          # FASTQ quality control using FastQC & MultiQC
│   ├── 2_run_kallisto.sh    # Transcript quantification with Kallisto
│   └── README.md            # Preprocessing instructions and dependencies
├── 1_Sex.R                  # Sex-differential expression analysis & visualization
├── 2_0_Cog.R                # Global cognition outcome analysis
├── 2_1_Cog_15fet_tissue.R   # Cognition-specific analysis in fetal samples
├── 2_2_Cog_15mat_tissue.R   # Cognition-specific analysis in maternal samples
├── 3_0_Lang.R               # Global language outcome analysis
├── 3_1_Lang_15fet_tissue.R  # Language-specific analysis in fetal samples
├── 3_2_Lang_15mat_tissue.R  # Language-specific analysis in maternal samples
├── 4_0_Motor.R              # Global motor outcome analysis
├── 4_1_Motor_15fet_tissue.R # Motor-specific analysis in fetal samples
├── 4_2_Motor_15mat_tissue.R # Motor-specific analysis in maternal samples
├── 5_Intersection.R         # Intersection of DEGs across analyses
├── 6_Cog_Heatmaps.R         # Heatmaps for cognitive outcome DEGs
├── 7_Cog_WGCNA.R            # Weighted gene co-expression network analysis (WGCNA) for cognition
└── README.md                # This file

Prerequisites
-------------
- Unix-like environment (Linux or macOS)
- R (>= 4.0) with the following packages:
  - data.table, dplyr, tximport, ensembldb, EnsDb.Hsapiens.v86, edgeR, limma, sva
  - cowplot, ggplot2, plotly, ComplexHeatmap, circlize, RColorBrewer, patchwork
  - GSEABase, GSVA, gprofiler2, clusterProfiler, msigdbr, enrichplot
- Kallisto (>= 0.46.0)
- FastQC & MultiQC (for QC)
- (Optional) Basespace CLI for downloading FASTQ from Illumina BaseSpace

`Note: GSEA, GSVA, gprofiler2, msigdbr should have a network available to proceed with the multiple queries.`

Usage
-----
1. Preprocessing:
   ```bash
   cd preprocessing
   ./1_run_qc.sh        # Generate raw QC reports
   ./2_run_kallisto.sh  # Build index and run pseudoalignment
   ```
2. Differential expression & downstream analysis:
   - Edit the working directory (`setwd(...)`) at the top of each `.R` script or launch R in the project root.
   - Run the R scripts in numerical order:
     ```bash
     Rscript 1_Sex.R
     Rscript 2_0_Cog.R
     # ... through 7_Cog_WGCNA.R
     ```
   - Outputs (plots, tables) will be saved under `results_paper/`.

Analysis overview
-----------------
- **1_Sex.R**: Sex-based differential expression analysis, PCA, volcano plots, heatmaps, and enrichment (GSEA, GO).
- **2*_*.R**, **3*_*.R**, **4*_*.R**: Parallel pipelines for cognition, language, and motor outcomes, each in global, fetal, and maternal contexts.
- **5_Intersection.R**: Identifies overlapping DEGs across outcome categories.
- **6_Cog_Heatmaps.R**: Generates detailed heatmaps for cognition-related DEGs.
- **7_Cog_WGCNA.R**: Constructs co-expression modules, correlates modules with cognitive metrics, and performs network visualization.

Citation
--------
Please cite “coming soon". 

License
-------
Need to decide the license
