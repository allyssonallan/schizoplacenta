# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 03_main_figures.R: Figures 2-5
# Source via run_all.R from the project root (folder that contains 0_TableCode.tsv).
# =============================================================================
save3 <- function(p, fn, w, h) { ggsave(file.path(FIG, paste0(fn, ".png")), p, width = w, height = h, dpi = 400, bg = "white"); ggsave(file.path(FIG, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white"); ggsave(file.path(FIG, paste0(fn, ".tiff")), p, width = w, height = h, dpi = 400, bg = "white", compression = "lzw") }
theme_pub <- function(b = 10) theme_classic(base_size = b) + theme(plot.title = element_text(face = "bold", size = rel(1)), axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), strip.background = element_blank(), strip.text = element_text(face = "bold"))
pal_grp <- c(DD = "#E84A5F", TD = "#00A8E8"); pal_face <- c(maternal = "#8C564B", fetal = "#17BECF")
pca <- prcomp(t(logcpm)); pcv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcdf <- data.table(pca$x[, 1:2], sample = colnames(logcpm)); pcdf <- merge(pcdf, meta, by = "sample")
clean_set <- function(x, width = 34, trunc = 58) str_wrap(str_trunc(str_to_sentence(gsub("_", " ", sub("^(REACTOME|KEGG|WP|BIOCARTA|PID|NABA|SA|SIG|ST)_", "", x))), trunc), width)

domain_figure <- function(d, fn, letter_title) {
  t <- prim_tabs[domain == d]; tf <- t[contrast == "fetal"]; grp <- paste0(d, "_grp")
  pa <- ggplot(pcdf, aes(PC1, PC2)) + geom_point(aes(shape = face, colour = .data[[grp]]), size = 3, stroke = 1.1) + scale_shape_manual(values = c(maternal = 1, fetal = 16), name = "Placental face", labels = c("maternal-facing", "fetal-facing")) +
    scale_colour_manual(values = pal_grp, name = "Group") + stat_ellipse(aes(colour = .data[[grp]]), linetype = 2, linewidth = 0.4) + labs(x = paste0("PC1 (", pcv[1], "%)"), y = paste0("PC2 (", pcv[2], "%)")) + theme_pub()
  tf[, cat := fifelse(adj.P.Val < 0.05, "FDR < 0.05", fifelse(P.Value < 0.05 & abs(logFC) >= 1, "P < 0.05, |log2FC| >= 1", "n.s."))]
  lab <- tf[order(P.Value)][1:10]
  pb <- ggplot(tf, aes(logFC, -log10(P.Value))) + geom_point(aes(colour = cat), size = 0.8, alpha = 0.7) + scale_colour_manual(values = c("FDR < 0.05" = "#C43C39", "P < 0.05, |log2FC| >= 1" = "#2C7FB8", "n.s." = "grey75"), name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") + geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") + geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = Inf, min.segment.length = 0) +
    labs(x = "log2 fold change (DD vs TD), fetal-facing", y = expression(-log[10]~italic(P))) + theme_pub() + theme(legend.position = "bottom", legend.text = element_text(size = 9))
  tab <- tf[order(P.Value)][1:15, .(gene, logFC = sprintf("%.2f", logFC), P = formatC(P.Value, format = "e", digits = 1), FDR = sprintf("%.2f", adj.P.Val))]
  pc <- ggplot() + annotation_custom(tableGrob(tab, rows = NULL, theme = ttheme_minimal(base_size = 9, padding = unit(c(3.5, 2.5), "mm")))) + theme_void()
  # heatmap: top 100 fetal-facing genes, all 30 samples
  hg <- tf[order(P.Value)][1:100, gene]; mat <- t(scale(t(logcpm[hg, ])))
  ord <- meta[order(face, get(grp), get(d))]; mat <- mat[, ord$sample]
  ha <- HeatmapAnnotation(Face = ord$face, Group = ord[[grp]], col = list(Face = pal_face, Group = pal_grp), annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9), simple_anno_size = unit(3, "mm"), annotation_legend_param = list(Face = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8), labels = c("maternal-facing   ", "fetal-facing   ")), Group = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))))
  ht <- Heatmap(mat, name = "Z-score", col = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")), cluster_columns = FALSE, column_split = factor(ifelse(ord$face == "fetal", "fetal-facing", "maternal-facing"), levels = c("maternal-facing", "fetal-facing")), show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), top_annotation = ha, row_title = "Top 100 fetal-facing DD-vs-TD genes", row_title_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 9), heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)))
  pd <- ggplot() + annotation_custom(grid.grabExpr(draw(ht, merge_legends = TRUE, padding = unit(c(1, 1, 1, 1), "mm")), width = unit(7.6, "in"), height = unit(7.6, "in"))) + theme_void()
  # enrichment: pre-ranked GSEA (MSigDB C2 canonical pathways), fetal & maternal-facing, top 12 by FDR (then |NES|) each
  ce <- gsea[domain == DOM[d] & model == "primary" & collection == "C2CP" & contrast %in% c("fetal", "maternal")]
  ce <- ce[order(contrast, FDR, -abs(NES))][, head(.SD, 12), by = contrast]
  ce[, label := clean_set(gene_set)]; ce[, contrast := factor(ifelse(contrast == "fetal", "Fetal-facing DD vs TD", "Maternal-facing DD vs TD"), levels = c("Fetal-facing DD vs TD", "Maternal-facing DD vs TD"))]
  ce[, sig := ifelse(FDR < 0.05, "FDR < 0.05", "FDR >= 0.05")]
  pe <- ggplot(ce, aes(NES, reorder(paste(contrast, label), NES))) + geom_segment(aes(x = 0, xend = NES, yend = reorder(paste(contrast, label), NES), colour = Direction), linewidth = 0.5) + geom_point(aes(size = NGenes, colour = Direction, alpha = sig)) +
    facet_wrap(~contrast, scales = "free_y", ncol = 1) + scale_y_discrete(labels = function(x) sub("^[^ ]+-facing DD vs TD ", "", x)) + scale_colour_manual(values = c(Up = "#C43C39", Down = "#2C7FB8"), name = "Direction in DD") +
    scale_alpha_manual(values = c(`FDR < 0.05` = 1, `FDR >= 0.05` = 0.35), name = NULL) + scale_size_continuous(range = c(1.5, 4.5), name = "Genes") + geom_vline(xintercept = 0, colour = "grey60") + labs(x = "Normalized enrichment score (GSEA)", y = NULL) + theme_pub(10) + theme(axis.text.y = element_text(size = 8.5, lineheight = 0.85), legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"), panel.spacing = unit(6, "mm"))
  top <- (pa + pb + pc) + plot_layout(widths = c(1, 1, 0.9))
  fig <- (top / ((pd | pe) + plot_layout(widths = c(1.2, 1)))) + plot_layout(heights = c(0.9, 1.6)) + plot_annotation(tag_levels = "a", title = letter_title, theme = theme(plot.title = element_text(size = 11, face = "bold")))
  save3(fig, fn, 14, 13)
}
domain_figure("cog", "Fig2_Cognitive", "Cognitive domain - adjusted, pregnancy-paired model (n = 15 pregnancies, 30 biopsies)")
domain_figure("lang", "Fig3_Language", "Language domain - adjusted, pregnancy-paired model (n = 15 pregnancies, 30 biopsies)")
domain_figure("mot", "Fig4_Motor", "Motor domain - adjusted, pregnancy-paired model (n = 15 pregnancies, 30 biopsies)")

# ---- Fig 5: intersection of nominal fetal-facing DEGs (top 25 up / 25 down by P among P<0.05 & |log2FC|>=1)
top50 <- rbindlist(lapply(names(DOM), function(d) { t <- prim_tabs[domain == d & contrast == "fetal" & P.Value < 0.05 & abs(logFC) >= 1]
  rbind(t[logFC > 0][order(P.Value)][1:min(25, .N)][, dir := "Up"], t[logFC < 0][order(P.Value)][1:min(25, .N)][, dir := "Down"])[, domain := DOM[d]] }))
wt(top50[, .(domain, direction = dir, gene, logFC, P.Value, adj.P.Val)], "S_top25_up_down_fetal_per_domain.tsv")
sets5 <- split(top50$gene, paste(top50$domain, top50$dir, sep = "_"))
allg <- unique(top50$gene); memb <- sapply(sets5, function(s) allg %in% s); rownames(memb) <- allg
combo <- apply(memb, 1, function(r) paste(names(sets5)[r], collapse = "&"))
cdt <- data.table(gene = allg, combo = combo)[, .(n = .N, genes = list(gene)), by = combo][order(-n)]
cdt[, x := seq_len(.N)]
bars <- ggplot(cdt, aes(x, n)) + geom_col(fill = "grey80", colour = "grey30", width = 0.85) + geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  geom_text(data = cdt[, .(gene = unlist(genes), y = seq_len(n) - 0.5), by = x], aes(x, y, label = gene), size = 2.6) + labs(y = "Intersection size", x = NULL) + theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dots <- rbindlist(lapply(seq_len(nrow(cdt)), function(i) data.table(x = i, set = names(sets5), on = names(sets5) %in% strsplit(cdt$combo[i], "&")[[1]])))
dots[, set := factor(set, levels = rev(names(sets5)))]
pm <- ggplot(dots, aes(x, set)) + geom_point(aes(colour = on), size = 2.5) + scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), guide = "none") + geom_line(data = dots[on == TRUE], aes(group = x), linewidth = 0.6) + labs(x = NULL, y = NULL) + theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ss <- data.table(set = factor(names(sets5), levels = rev(names(sets5))), n = lengths(sets5))
ps <- ggplot(ss, aes(n, set)) + geom_col(fill = "grey40", width = 0.7) + scale_x_reverse() + labs(x = "Set size", y = NULL) + theme_pub() + theme(axis.text.y = element_blank())
fig5 <- (plot_spacer() + bars + plot_layout(widths = c(0.25, 1))) / (ps + pm + plot_layout(widths = c(0.25, 1))) + plot_layout(heights = c(2, 1))
save3(fig5, "Fig5_intersection", 11, 9)
