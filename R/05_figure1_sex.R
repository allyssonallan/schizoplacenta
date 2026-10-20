# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 05_figure1_sex.R: Figure 1 (sex-biased expression, male vs female placentas), same layout as Figures 2-4
#   model: ~ sex + face + maternal age (z), pregnancy-blocked (01_differential_expression.R -> DE_sex_adjusted.tsv)
#   (a) PCA by sex and face, (b) volcano, (c) top 10 DEGs by adjusted p, (d) heatmap of all FDR<0.05 genes,
#   (e) over-representation of the sex-biased DEGs (higher in males / higher in females) across the ORA sources + table of significant terms
# Source via run_all.R from the project root (after 03_main_figures.R, which defines save3/theme_pub/pcdf).
# =============================================================================
sex_tab <- fread(file.path(TAB, "DE_sex_adjusted.tsv"))
pal_sex <- c(F = "#E84A5F", M = "#00A8E8")
# (a) PCA
pa <- ggplot(pcdf, aes(PC1, PC2)) + geom_point(aes(shape = face, colour = sex), size = 3, stroke = 1.1) + scale_shape_manual(values = c(maternal = 1, fetal = 16), name = "Placental face", labels = c("maternal-facing", "fetal-facing")) +
  scale_colour_manual(values = pal_sex, name = "Offspring sex", labels = c(F = "Female", M = "Male")) + stat_ellipse(aes(colour = sex), linetype = 2, linewidth = 0.4) + labs(x = paste0("PC1 (", pcv[1], "%)"), y = paste0("PC2 (", pcv[2], "%)")) + theme_pub()
# (b) volcano (male vs female), y = -log10 adjusted p as in the first submission
sex_tab[, cat := fifelse(adj.P.Val < 0.05 & abs(logFC) >= 2, "FDR < 0.05, |log2FC| >= 2", fifelse(adj.P.Val < 0.05, "FDR < 0.05", "n.s."))]
lab <- sex_tab[order(adj.P.Val)][1:10]
pb <- ggplot(sex_tab, aes(logFC, -log10(adj.P.Val))) + geom_point(aes(colour = cat), size = 0.8, alpha = 0.7) + scale_colour_manual(values = c(`FDR < 0.05, |log2FC| >= 2` = "black", `FDR < 0.05` = "#C43C39", n.s. = "grey75"), name = NULL) +
  geom_vline(xintercept = c(-2, 2), linetype = 2, colour = "grey60") + geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") + geom_text_repel(data = lab, aes(label = gene), size = 3, fontface = "italic", max.overlaps = Inf, min.segment.length = 0) +
  labs(x = "log2 fold change (male vs female)", y = expression(-log[10]~"adjusted "~italic(P))) + theme_pub() + theme(legend.position = "bottom", legend.text = element_text(size = 9))
# (c) top 10 by adjusted p
tab <- sex_tab[order(adj.P.Val)][1:10, .(gene, logFC = sprintf("%.2f", logFC), P = formatC(P.Value, format = "e", digits = 1), FDR = formatC(adj.P.Val, format = "e", digits = 1))]
pc <- ggplot() + annotation_custom(tableGrob(tab, rows = NULL, theme = ttheme_minimal(base_size = 9, padding = unit(c(3.5, 2.5), "mm")))) + theme_void()
# (d) heatmap of all sex-biased DEGs (FDR < 0.05), columns ordered by face then sex
hg <- sex_tab[adj.P.Val < 0.05, gene]; mat <- t(scale(t(logcpm[hg, ])))
ord <- meta[order(face, sex, participant)]; mat <- mat[, ord$sample]
ha <- HeatmapAnnotation(Face = ord$face, Sex = ord$sex, col = list(Face = pal_face, Sex = pal_sex), annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9), simple_anno_size = unit(3, "mm"),
                        annotation_legend_param = list(Face = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8), labels = c("maternal-facing   ", "fetal-facing   ")), Sex = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8), labels = c("Female", "Male"))))
ht <- Heatmap(mat, name = "Z-score", col = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")), cluster_columns = FALSE, column_split = factor(ifelse(ord$face == "fetal", "fetal-facing", "maternal-facing"), levels = c("maternal-facing", "fetal-facing")),
              show_row_names = length(hg) <= 70, row_names_gp = gpar(fontsize = 6, fontface = "italic"), column_names_gp = gpar(fontsize = 7), top_annotation = ha, row_title = sprintf("%d sex-biased genes (FDR < 0.05)", length(hg)), row_title_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)))
pd <- ggplot() + annotation_custom(grid.grabExpr(draw(ht, merge_legends = TRUE, padding = unit(c(1, 1, 1, 1), "mm")), width = unit(7.6, "in"), height = unit(7.6, "in"))) + theme_void()
# (e) over-representation of the sex-biased DEGs (higher in males / higher in females), same sources as Figs S4/S7/S10
universe <- sex_tab$gene
deg <- list(`higher in males` = sex_tab[adj.P.Val < 0.05 & logFC > 0, gene], `higher in females` = sex_tab[adj.P.Val < 0.05 & logFC < 0, gene])
ora_sex <- rbindlist(lapply(names(ora_sets), function(src) {
  sets <- lapply(ora_sets[[src]], intersect, universe); sets <- sets[lengths(sets) >= 3 & lengths(sets) <= 500]  # min 3 expressed genes: small Y-linked / histone terms
  rbindlist(lapply(names(deg), function(dr) { q <- deg[[dr]]; N <- length(universe); K <- length(q)
    r <- rbindlist(lapply(names(sets), function(s) { m <- length(sets[[s]]); k <- length(intersect(sets[[s]], q))
      data.table(term = s, term_size = m, overlap = k, p = if (k == 0) 1 else phyper(k - 1, m, N - m, K, lower.tail = FALSE), genes = paste(intersect(sets[[s]], q), collapse = ";")) }))
    r[, padj := p.adjust(p, "BH")][, `:=`(direction = dr, source = src, query_size = K)] }))
}))
wt(ora_sex[order(direction, padj)], "ORA_sex_biased_DEGs_all_sources.tsv")
x <- ora_sex[overlap > 0]; x[, source := factor(source, levels = c("GO:BP", "GO:CC", "GO:MF", "HP", "KEGG", "REAC", "PID", "BIOCARTA"))]; x[, direction := factor(direction, levels = names(deg))]
set.seed(1); x[, xpos := as.numeric(source) + runif(.N, -0.35, 0.35)]; x[, sig := padj < 0.05]
labx <- x[padj < 0.05][order(padj)][, head(.SD, 3), by = .(direction, source)]; labx[, term_lab := str_trunc(str_to_sentence(gsub("_", " ", sub("^(GOBP|GOCC|GOMF|HP|KEGG|REACTOME|PID|BIOCARTA)_", "", term))), 40)]
pal <- c(`GO:BP` = "#FF9900", `GO:CC` = "#109618", `GO:MF` = "#DC3912", HP = "#990099", KEGG = "#DD4477", REAC = "#3366CC", PID = "#66AA00", BIOCARTA = "#B82E2E")
pe1 <- ggplot(x, aes(xpos, -log10(p), colour = source)) + geom_point(aes(size = overlap, alpha = sig)) + scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.25), guide = "none") +
  geom_text_repel(data = labx, aes(label = term_lab), size = 2.6, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4, show.legend = FALSE) +
  facet_wrap(~direction, ncol = 2, labeller = labeller(direction = function(v) paste("Sex-biased DEGs", v))) + scale_x_continuous(breaks = seq_along(levels(x$source)), labels = levels(x$source)) + scale_colour_manual(values = pal, guide = "none") + scale_size_continuous(range = c(1, 4), name = "DEGs in term") +
  labs(x = NULL, y = expression(-log[10]~italic(P)~"(solid: adjusted "~italic(P)~"< 0.05)")) + theme_pub(10) + theme(legend.position = "bottom", axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1))
tabe <- ora_sex[padj < 0.05][order(direction, padj)][, .(direction, source, term = str_trunc(gsub("_", " ", sub("^(GOBP|GOCC|GOMF|HP|KEGG|REACTOME|PID|BIOCARTA)_", "", term)), 48), size = term_size, DEGs = overlap, FDR = formatC(padj, format = "e", digits = 1))][1:min(.N, 14)]
pe2 <- if (nrow(tabe)) ggplot() + annotation_custom(tableGrob(tabe, rows = NULL, theme = ttheme_minimal(base_size = 7.5, padding = unit(c(2.5, 1.8), "mm")))) + theme_void() else ggplot() + theme_void()
pe <- (pe1 / pe2) + plot_layout(heights = c(1.3, 1))
top <- (pa + pb + pc) + plot_layout(widths = c(1, 1, 0.9))
fig1 <- (top / ((pd | pe) + plot_layout(widths = c(1.2, 1)))) + plot_layout(heights = c(0.9, 1.6)) + plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e", "")), title = "Sex-biased placental expression - adjusted, pregnancy-paired model (n = 15 pregnancies, 30 biopsies)", theme = theme(plot.title = element_text(size = 11, face = "bold")))
save3(fig1, "Fig1_Sex", 14, 13)
cat("Figure 1 written\n")
