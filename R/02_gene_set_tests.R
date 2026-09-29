# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 02_gene_set_tests.R: competitive gene-set tests (limma::camera) on MSigDB and published placental sets
# Source via run_all.R from the project root (folder that contains 0_TableCode.tsv).
# =============================================================================
run_camera <- function(obj, sets, label) {
  idx <- ids2indices(sets, rownames(obj$v)); idx <- idx[lengths(idx) >= 5]
  rbindlist(lapply(colnames(obj$cm), function(cc) {
    r <- camera(obj$v, idx, obj$design, contrast = obj$cm[, cc], inter.gene.cor = 0.01)
    as.data.table(r, keep.rownames = "gene_set")[, `:=`(contrast = cc, collection = label)]
  }))
}
cam <- list()
for (d in names(DOM)) {
  for (cl in names(gs)) cam[[paste(d, cl)]] <- run_camera(primary[[d]], gs[[cl]], cl)[, `:=`(domain = DOM[d], model = "primary")]
  cam[[paste(d, "EXT")]] <- run_camera(primary[[d]], extsets_ok, "external")[, `:=`(domain = DOM[d], model = "primary")]
  cam[[paste(d, "EXT LOO")]] <- run_camera(sens[[paste0(d, "_LOO_RN01")]], extsets_ok, "external")[, `:=`(domain = DOM[d], model = "LOO_RN01")]
  cam[[paste(d, "C2CP LOO")]] <- run_camera(sens[[paste0(d, "_LOO_RN01")]], gs$C2CP, "C2CP")[, `:=`(domain = DOM[d], model = "LOO_RN01")]
  cam[[paste(d, "C2CP cont")]] <- run_camera(sens[[paste0(d, "_continuous")]], gs$C2CP, "C2CP")[, `:=`(domain = DOM[d], model = "continuous")]
  cam[[paste(d, "EXT cont")]] <- run_camera(sens[[paste0(d, "_continuous")]], extsets_ok, "external")[, `:=`(domain = DOM[d], model = "continuous")]
}
cam <- rbindlist(cam)
setcolorder(cam, c("domain", "model", "collection", "contrast", "gene_set", "NGenes", "Direction", "PValue", "FDR"))
wt(cam, "GSEA_camera_all.tsv")
for (cl in unique(cam$collection)) wt(cam[collection == cl & model == "primary"][order(domain, contrast, FDR)], sprintf("GSEA_camera_%s_primary.tsv", cl))
wt(cam[collection == "external"][order(domain, model, contrast, FDR)], "GSEA_camera_external_sets_all_models.tsv")
cam_sum <- cam[, .(sets_FDR05 = sum(FDR < 0.05), up = sum(FDR < 0.05 & Direction == "Up"), down = sum(FDR < 0.05 & Direction == "Down"), tested = .N), by = .(domain, model, collection, contrast)]
wt(cam_sum, "S_camera_summary.tsv")

# leading-edge style annotation: top genes of top sets (fetal, C2CP)
lead <- rbindlist(lapply(names(DOM), function(d) {
  t <- prim_tabs[domain == d & contrast == "fetal"]; top <- cam[domain == DOM[d] & model == "primary" & collection == "C2CP" & contrast == "fetal"][order(FDR)][1:15]
  rbindlist(lapply(seq_len(nrow(top)), function(i) {
    g <- intersect(gs$C2CP[[top$gene_set[i]]], t$gene); tt <- t[gene %in% g][order(if (top$Direction[i] == "Down") t else -t)][1:min(10, .N)]
    data.table(domain = DOM[d], gene_set = top$gene_set[i], Direction = top$Direction[i], FDR = top$FDR[i], n_in_set = length(g), leading_genes = paste(tt$gene, collapse = ";"))
  }))
}))
wt(lead, "S_camera_C2CP_fetal_top_sets_leading_genes.tsv")

saveRDS(cam, file.path(RDS, "gene_set_results.rds"))
