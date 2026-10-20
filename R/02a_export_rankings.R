# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 02a_export_rankings.R: gene rankings (moderated t, DD - TD) for pre-ranked GSEA
#   primary model: fetal-facing, maternal-facing; sensitivity: leave-one-out (RN01) and continuous score (fetal-facing)
# Source via run_all.R from the project root.
# =============================================================================
rk <- data.table(gene = rownames(logcpm))
for (d in names(DOM)) {
  for (cc in c("fetal", "maternal")) { t <- prim_tabs[domain == d & contrast == cc]; rk[[paste(d, "primary", cc, sep = "|")]] <- t$t[match(rk$gene, t$gene)] }
  t <- sens_tabs[domain == d & model == "LOO_RN01" & contrast == "fetal"]; rk[[paste(d, "LOO_RN01", "fetal", sep = "|")]] <- t$t[match(rk$gene, t$gene)]
  t <- sens_tabs[domain == d & model == "continuous" & contrast == "fetal"]; rk[[paste(d, "continuous", "fetal", sep = "|")]] <- -t$t[match(rk$gene, t$gene)]  # higher score = TD; flip so positive = DD-like
}
rkf <- file.path(RDS, "gsea_rankings.tsv")
if (!file.exists(rkf) || !isTRUE(all.equal(as.data.frame(fread(rkf)), as.data.frame(rk), tolerance = 1e-8, check.attributes = FALSE))) fwrite(rk, rkf, sep = "\t")  # keep the mtime when unchanged, so GSEA is not re-run
# gene-set collections as GMT for the python step
write_gmt <- function(sets, f) writeLines(vapply(names(sets), function(n) paste(c(n, "na", unique(sets[[n]])), collapse = "\t"), ""), f)
write_gmt(extsets_ok, file.path(RDS, "external_sets.gmt"))
if (!is.null(gs$C7)) write_gmt(gs$C7, file.path(RDS, "C7.gmt"))
cat("Rankings and GMTs written to", RDS, "\n")
