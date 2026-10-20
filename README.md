# schizoplacenta

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17333992.svg)](https://doi.org/10.5281/zenodo.17333992)

Placental transcriptome (bulk RNA-seq, maternal-facing and fetal-facing biopsies) from 15 pregnancies with third-trimester SARS-CoV-2 infection, paired with the children's Bayley-III outcomes (cognitive, language, motor) at 6–24 months. Code for the manuscript *Placental transcriptomic dysregulation following maternal SARS-CoV-2 infection links immune pathways to early neurodevelopment and schizophrenia risk* (Molecular Psychiatry, in revision).

## Layout

```
.
├── preprocessing/            FASTQ QC (FastQC/MultiQC) and kallisto quantification
├── data/                     small derived inputs shipped with the repo
│   ├── external_gene_sets.tsv                       published placental gene sets (Ursini 2023, Cilleros-Portet 2025, Gao 2023)
│   └── GSE171381_key_gene_donor_celltype_pseudobulk.tsv   donor x cell-type pseudobulk from GEO GSE171381
├── R/
│   ├── 00_setup.R                     packages, paths, sample metadata, counts (tximport, cached), gene sets
│   ├── 01_differential_expression.R   primary model, sensitivity models, sex model, gene-level tables
│   ├── 02a_export_rankings.R          rankings (moderated t) and GMTs for GSEA
│   ├── 02b_gsea_prerank.py            pre-ranked GSEA (gseapy), chunked and resumable
│   ├── 02_gene_set_tests.R            GSEA results, over-representation of nominal DEGs, limma::fry non-directional test
│   ├── 03_main_figures.R              Figures 2–5
│   ├── 04_supplementary_figures.R     Figures S4–S22
│   └── 05_figure1_sex.R               Figure 1 (sex-biased expression, adjusted sex model)
├── run_all.R                 entry point (see Usage)
├── 1_Sex.R                   sex-differential expression and Figure 1 (first submission; the covariate-adjusted re-fit is in R/01)
├── legacy/                   first-submission outcome scripts, superseded (kept for provenance, see legacy/README.md)
├── renv.lock                 package versions
└── results_revision/         created by run_all.R (tables/, figures/, rds/); not tracked
```

## Statistical model

One covariate-adjusted, pregnancy-paired limma-voom model per Bayley-III domain, fitted to all 30 biopsies:

```
log2 expression ~ face * outcome + maternal_age_z + offspring_sex
```

`face` = maternal-facing (reference) or fetal-facing biopsy; `outcome` = DD (mean composite score < 90) or TD; pregnancy is the blocking factor (`duplicateCorrelation`), voom weights are re-estimated with the consensus correlation, `eBayes(robust = TRUE)`, Benjamini–Hochberg within contrast. Contrasts: fetal-facing DD−TD, maternal-facing DD−TD, face × outcome. All libraries were prepared and sequenced in one batch, so batch is not a model term. Sensitivity models: unadjusted paired, RIN complete-case, closest-to-24-month assessment, continuous score, leave-one-out (participant RN01). Pathway-level inference uses pre-ranked GSEA (the method of the first submission), implemented with gseapy on the moderated t statistics of each contrast (`python3 R/02b_gsea_prerank.py <root>`; `pip install gseapy`), against MSigDB C2 (all, canonical pathways), C7 and the published placental sets, with BH FDR within collection; over-representation of the nominal fetal-facing DEGs is tested across GO/HPO (from `msigdb_C5.tsv`), KEGG, Reactome, PID and BioCarta; `limma::fry` gives a non-directional P for the published sets. `Rscript run_all.R genesets` re-runs this step from the cached differential-expression models.

## Inputs

Place in the project root (or point `SCHIZO_ROOT` to the folder that holds them):

* `0_TableCode.tsv` — sample/participant metadata (one row per Bayley-III visit)
* kallisto output folders listed in `0_TableCode.tsv` (`<sample>/abundance.h5`)
* `c2.all.v7.1.symbols.gmt`, `c2.cp.v7.1.symbols.gmt` (MSigDB v7.1)
* optional: `msigdb_C7.tsv`, `msigdb_C5.tsv` exported with `msigdbr` (`gs_name, gs_subcat, gene_symbol`) to add immunologic signatures and HPO terms

Raw reads are deposited at ENA (accession to be added on acceptance).

## Usage

```bash
Rscript run_all.R                       # models + gene-set tests + all figures (≈15 min)
Rscript run_all.R figures               # redraw figures from cached models (results_revision/rds)
SCHIZO_ROOT=/path/to/data Rscript run_all.R
```

Preprocessing (once): `cd preprocessing && ./1_run_qc.sh && ./2_run_kallisto.sh`.

## Requirements

R ≥ 4.3 with data.table, dplyr, edgeR, limma, tximport, ensembldb, EnsDb.Hsapiens.v86, ggplot2, ggrepel, patchwork, cowplot, gridExtra, ComplexHeatmap, circlize, RColorBrewer, stringr, scales (`renv::restore()` reproduces the exact versions). No network access is needed once the GMT files are present.

## Outputs

`results_revision/tables/` — `DE_primary_<domain>_<contrast>.tsv`, `DE_sensitivity_*.tsv`, `DE_sex_adjusted.tsv`, `GSEA_*.tsv`, `ORA_*.tsv`, `FRY_*.tsv`, `S_*` summaries (model diagnostics, effect stability, gene-set manifest, candidate genes across models).
`results_revision/figures/` — Figures 2–5 (PNG/PDF/TIFF 400 dpi) and S4–S22 (PNG/PDF).

## Citation

de Farias AA et al. *schizoplacenta*. Zenodo. https://doi.org/10.5281/zenodo.17333992
