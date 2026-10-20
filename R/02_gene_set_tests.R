# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 02_gene_set_tests.R: pathway-level analyses
#   1. pre-ranked GSEA (the method of the first submission), run by R/02b_gsea_prerank.py (gseapy) on the
#      moderated t statistics exported by 02a_export_rankings.R: MSigDB C2 (all, CP), C7, published placental sets
#   2. over-representation (hypergeometric) of the nominal fetal-facing DEGs, g:Profiler-like sources (Figs. S4/S7/S10)
#   3. non-directional self-contained test (limma::fry, mixed statistic) for the published placental sets
# Source via run_all.R from the project root.
# =============================================================================
source("R/02a_export_rankings.R")
gsea_file <- file.path(TAB, "GSEA_prerank_all.tsv")
if (!file.exists(gsea_file) || file.mtime(gsea_file) < file.mtime(file.path(RDS, "gsea_rankings.tsv"))) {
  cat("Running pre-ranked GSEA (python3 + gseapy); this takes ~1-2 h on 2 cores ...\n")
  status <- system2("python3", c(file.path("R", "02b_gsea_prerank.py"), shQuote(ROOT), "1000"))
  if (status != 0) stop("GSEA step failed; install gseapy (pip install gseapy) and re-run")
}
gsea <- fread(gsea_file)
gsea[, domain := DOM[domain]]
gsea_sum <- gsea[, .(sets_FDR05 = sum(FDR < 0.05), up = sum(FDR < 0.05 & NES > 0), down = sum(FDR < 0.05 & NES < 0), tested = .N), by = .(domain, model, collection, contrast)]
wt(gsea_sum, "S_GSEA_summary.tsv")
for (cl in unique(gsea$collection)) wt(gsea[collection == cl & model == "primary"][order(domain, contrast, FDR)], sprintf("GSEA_%s_primary.tsv", cl))
wt(gsea[collection == "external"][order(domain, model, contrast, FDR)], "GSEA_external_sets_all_models.tsv")

# ---- 2. over-representation of nominal fetal-facing DEGs (up / down in DD), hypergeometric, BH within source
ora <- rbindlist(lapply(names(DOM), function(d) {
  t <- prim_tabs[domain == d & contrast == "fetal"]; universe <- t$gene
  deg <- list(`higher in DD` = t[P.Value < 0.05 & logFC >= 1, gene], `lower in DD` = t[P.Value < 0.05 & logFC <= -1, gene])
  rbindlist(lapply(names(ora_sets), function(src) {
    sets <- lapply(ora_sets[[src]], intersect, universe); sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
    if (!length(sets)) return(NULL)
    rbindlist(lapply(names(deg), function(dr) {
      q <- deg[[dr]]; N <- length(universe); K <- length(q)
      r <- rbindlist(lapply(names(sets), function(s) { m <- length(sets[[s]]); k <- length(intersect(sets[[s]], q))
        data.table(term = s, term_size = m, overlap = k, p = if (k == 0) 1 else phyper(k - 1, m, N - m, K, lower.tail = FALSE), genes = paste(intersect(sets[[s]], q), collapse = ";")) }))
      r[, padj := p.adjust(p, "BH")][, `:=`(domain = DOM[d], direction = dr, source = src, query_size = K)]
    }))
  }))
}))
setcolorder(ora, c("domain", "direction", "source", "term", "term_size", "overlap", "query_size", "p", "padj", "genes"))
wt(ora[order(domain, direction, padj)], "ORA_fetal_nominal_DEGs_all_sources.tsv")
wt(ora[padj < 0.05][order(domain, direction, padj)], "ORA_fetal_nominal_DEGs_significant.tsv")

# ---- 3. non-directional test for the published sets (limma::fry; PValue.Mixed ignores the sign of the changes)
run_fry <- function(obj, sets, label) {
  idx <- ids2indices(sets, rownames(obj$v)); idx <- idx[lengths(idx) >= 5]
  rbindlist(lapply(colnames(obj$cm), function(cc) { r <- fry(obj$v, idx, obj$design, contrast = obj$cm[, cc]); as.data.table(r, keep.rownames = "gene_set")[, `:=`(contrast = cc, collection = label)] }))
}
fr <- rbindlist(lapply(names(DOM), function(d) rbind(run_fry(primary[[d]], extsets_ok, "external")[, `:=`(domain = DOM[d], model = "primary")],
                                                     run_fry(sens[[paste0(d, "_LOO_RN01")]], extsets_ok, "external")[, `:=`(domain = DOM[d], model = "LOO_RN01")])))
setcolorder(fr, c("domain", "model", "collection", "contrast", "gene_set", "NGenes", "Direction", "PValue", "FDR", "PValue.Mixed", "FDR.Mixed"))
wt(fr, "FRY_external_sets_directional_and_mixed.tsv")
gsea[fr, on = c("domain", "model", "collection", "contrast", "gene_set"), `:=`(PValue_mixed = i.PValue.Mixed, FDR_mixed = i.FDR.Mixed)]
wt(gsea[collection == "external"][order(domain, model, contrast, FDR)], "GSEA_external_sets_all_models.tsv")

# ---- leading-edge genes of the top fetal-facing C2:CP sets
lead <- gsea[model == "primary" & collection == "C2CP" & contrast == "fetal"][order(domain, FDR, -abs(NES))][, head(.SD, 15), by = domain][, .(domain, gene_set, NES, FDR, NGenes, leading_edge)]
wt(lead, "S_GSEA_C2CP_fetal_top_sets_leading_genes.tsv")
saveRDS(gsea, file.path(RDS, "gene_set_results.rds")); saveRDS(ora, file.path(RDS, "ora_results.rds")); saveRDS(fr, file.path(RDS, "fry_results.rds"))
