#!/usr/bin/env Rscript
# schizoplacenta — revised analysis. Run from the project root:
#   Rscript run_all.R            # everything
#   Rscript run_all.R genesets   # re-run gene-set tests + figures from cached DE models (results_revision/rds)
#   Rscript run_all.R figures    # only re-draw figures from cached models and gene-set results
# Set SCHIZO_ROOT to point to the data folder if it differs from the working directory.
step <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "all"
source("R/00_setup.R")
if (step == "all") {
  source("R/01_differential_expression.R")
  source("R/02_gene_set_tests.R")
} else if (step == "genesets") {
  m <- readRDS(file.path(RDS, "de_models.rds")); list2env(m, envir = globalenv())
  source("R/02_gene_set_tests.R")
} else {
  m <- readRDS(file.path(RDS, "de_models.rds")); list2env(m, envir = globalenv())
  gsea <- readRDS(file.path(RDS, "gene_set_results.rds")); ora <- readRDS(file.path(RDS, "ora_results.rds")); fr <- readRDS(file.path(RDS, "fry_results.rds"))
}
if (step != "all") rawkey <- fread(file.path(TAB, "S_fetal_log2CPM_previously_highlighted_genes.tsv"))
source("R/03_main_figures.R")
source("R/04_supplementary_figures.R")
source("R/05_figure1_sex.R")
