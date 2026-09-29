# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 04_supplementary_figures.R: Figures S4-S22
# Source via run_all.R from the project root (folder that contains 0_TableCode.tsv).
# =============================================================================
XLAB <- expression("signed "-log[10]~"FDR (competitive gene-set test)")
clean <- function(x) str_trunc(str_to_sentence(gsub("_", " ", x)), 62)
lolli2 <- function(ce, fn, title, n = 15, w = 9, h = 14) {
  ce <- ce[order(contrast, FDR)]
  ce <- rbind(ce[Direction == "Up", head(.SD, n), by = contrast], ce[Direction == "Down", head(.SD, n), by = contrast])
  ce[, label := clean(gene_set)]; ce[, signed := ifelse(Direction == "Up", 1, -1) * -log10(pmax(FDR, 1e-30))]
  ce[, contrast := factor(ifelse(contrast == "fetal", "Fetal-facing DD vs TD", "Maternal-facing DD vs TD"), levels = c("Fetal-facing DD vs TD", "Maternal-facing DD vs TD"))]
  ce[, key := paste(contrast, label)]
  p <- ggplot(ce, aes(signed, reorder(key, signed))) + geom_segment(aes(x = 0, xend = signed, yend = reorder(key, signed), colour = Direction)) + geom_point(aes(size = NGenes, colour = Direction)) +
    facet_wrap(~contrast, scales = "free_y", ncol = 1) + scale_y_discrete(labels = function(x) sub("^[^ ]+-facing DD vs TD ", "", x)) +
    scale_colour_manual(values = c(Up = "#C43C39", Down = "#2C7FB8"), name = "Direction in DD") + scale_size_continuous(range = c(1.5, 5), name = "Genes tested") +
    geom_vline(xintercept = 0, colour = "grey60") + geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = 2, colour = "grey75") +
    labs(x = XLAB, y = NULL, title = str_wrap(title, 75)) + theme_pub(11) + theme(axis.text.y = element_text(size = 9), legend.position = "bottom", legend.box = "vertical")
  ggsave(file.path(FIG, paste0(fn, ".png")), p, width = w, height = h, dpi = 300, bg = "white"); ggsave(file.path(FIG, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}
imm <- "INTERFERON|IFN|COVID|SARS|VIRUS|VIRAL|INFLUENZA|INFECTION|INFLAMM|CYTOKINE|INNATE|TLR|NFKB|ANTIVIRAL|ANTIGEN|MHC|COMPLEMENT|CHEMOKINE|INTERLEUKIN|IL6|IL2|IL12|TNF|IMMUN"
for (d in c("Cognitive", "Language", "Motor")) {
  base <- cam[domain == d & model == "primary" & contrast %in% c("fetal", "maternal")]
  s <- switch(d, Cognitive = c("S4", "S5", "S6"), Language = c("S7", "S8", "S9"), Motor = c("S10", "S11", "S12"))
  lolli2(base[collection == "C2ALL"], paste0("Fig", s[1], "_", d, "_C2all"), paste(d, "domain - MSigDB C2 curated gene sets, top 15 per direction"))
  lolli2(base[collection == "C2CP"],  paste0("Fig", s[2], "_", d, "_C2CP"),  paste(d, "domain - MSigDB C2 canonical pathways, top 15 per direction"))
  if (!is.null(gs$C7)) lolli2(base[collection == "C7"], paste0("Fig", s[3], "_", d, "_C7"), paste(d, "domain - MSigDB C7 immunologic signatures, top 15 per direction"))
  else lolli2(base[collection == "C2ALL" & grepl(imm, gene_set)], paste0("Fig", s[3], "_", d, "_immune"), paste(d, "domain - immune, interferon and infection-related curated gene sets (MSigDB C2)"))
}

# face-specific supplementary figure (S13-S15): fetal vs maternal-facing contrasts and interaction
face_fig <- function(d, fn) {
  t <- prim_tabs[domain == d]; w <- dcast(t[, .(gene, contrast, logFC, P.Value)], gene ~ contrast, value.var = c("logFC", "P.Value"))
  vol <- function(cc, ttl) { x <- t[contrast == cc]; x[, sig := P.Value < 0.05 & abs(logFC) >= 1]; ggplot(x, aes(logFC, -log10(P.Value))) + geom_point(aes(colour = sig), size = 0.6, alpha = 0.6) + scale_colour_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey75"), guide = "none") +
    geom_text_repel(data = x[order(P.Value)][1:8], aes(label = gene), size = 2.8, max.overlaps = Inf) + geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") + geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") + labs(title = ttl, x = "log2 fold change (DD vs TD)", y = expression(-log[10]~italic(P))) + theme_pub(11) }
  pa <- vol("fetal", "Fetal-facing DD vs TD"); pb <- vol("maternal", "Maternal-facing DD vs TD"); pc <- vol("interaction", "Face x outcome interaction")
  r <- cor(w$logFC_fetal, w$logFC_maternal, method = "spearman")
  pd <- ggplot(w, aes(logFC_maternal, logFC_fetal)) + geom_point(size = 0.5, alpha = 0.3, colour = "#264653") + geom_abline(linetype = 2, colour = "#C43C39") + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    labs(title = sprintf("Concordance of DD-TD effects between faces (Spearman rho = %.2f)", r), x = "Maternal-facing log2FC", y = "Fetal-facing log2FC") + theme_pub(11)
  fig <- (pa | pb | pc) / pd + plot_annotation(tag_levels = "a", title = paste(DOM[d], "domain - face-specific contrasts from the joint adjusted model"))
  ggsave(file.path(FIG, paste0(fn, ".png")), fig, width = 12, height = 8, dpi = 300, bg = "white"); ggsave(file.path(FIG, paste0(fn, ".pdf")), fig, width = 12, height = 8, device = cairo_pdf)
}
face_fig("cog", "FigS13_Cognitive_faces"); face_fig("lang", "FigS14_Language_faces"); face_fig("mot", "FigS15_Motor_faces")

# S20: published placental gene sets, primary and leave-one-out models
ce <- cam[collection == "external" & contrast %in% c("fetal", "maternal") & model %in% c("primary", "LOO_RN01")]
ce[, signed := ifelse(Direction == "Up", 1, -1) * -log10(pmax(FDR, 1e-30))]
ce[, contrast := factor(ifelse(contrast == "fetal", "Fetal-facing", "Maternal-facing"), levels = c("Fetal-facing", "Maternal-facing"))]
ce[, model := factor(ifelse(model == "primary", "All 15 pregnancies", "Without RN01"), levels = c("All 15 pregnancies", "Without RN01"))]
short <- c("SCZ TWAS Ursini 2023" = "SCZ placental TWAS (Ursini 2023)", "SCZ prioritized Ursini 2023" = "SCZ placenta-prioritized (Ursini 2023)",
           "SCZ GWAS placenta-implicated (Ursini 2023 S1; 56)" = "SCZ GWAS placenta-implicated (Ursini 2023)", "PlacPRS1 schizophrenia genes (Ursini 2021; 66)" = "PlacPRS1 genes (Ursini 2021)",
           "SCZ CpG annotation Cilleros 2025" = "SCZ CpG-annotated (Cilleros-Portet 2025)", "SCZ eQTM Cilleros 2025" = "SCZ eQTM targets (Cilleros-Portet 2025)",
           "ADHD TWAS Ursini 2023" = "ADHD placental TWAS (Ursini 2023)", "COVID up Gao 2023" = "Up in COVID placenta (Gao 2023)", "COVID down Gao 2023" = "Down in COVID placenta (Gao 2023)")
ce[, gene_set := ifelse(gene_set %in% names(short), short[gene_set], gene_set)]
ce[, gene_set := factor(gene_set, levels = rev(sort(unique(gene_set))))]
p20 <- ggplot(ce, aes(signed, gene_set, colour = Direction, shape = model)) + geom_vline(xintercept = 0, colour = "grey60") + geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = 2, colour = "grey75") +
  geom_point(aes(size = NGenes), position = position_dodge(width = 0.6), alpha = 0.9) + facet_grid(domain ~ contrast) + scale_colour_manual(values = c(Up = "#C43C39", Down = "#2C7FB8"), name = "Direction in DD") + scale_shape_manual(values = c(16, 1), name = "Model") + scale_size_continuous(range = c(1.5, 5), name = "Genes tested") +
  labs(x = XLAB, y = NULL, title = "Published placental gene sets in the DD-vs-TD contrasts") + theme_pub(11) + theme(axis.text.y = element_text(size = 9.5), legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(FIG, "FigS20_external_gene_sets.png"), p20, width = 11, height = 11, dpi = 300, bg = "white"); ggsave(file.path(FIG, "FigS20_external_gene_sets.pdf"), p20, width = 11, height = 11, device = cairo_pdf)

# S22: leave-one-out (RN01) and continuous-score sensitivity at the pathway level (C2:CP fetal)
cc <- cam[collection == "C2CP" & contrast == "fetal" & model %in% c("primary", "LOO_RN01", "continuous")]
cc[, signed := ifelse(Direction == "Up", 1, -1) * -log10(pmax(FDR, 1e-30))]
ww <- dcast(cc, domain + gene_set ~ model, value.var = "signed")
ww[, `:=`(continuous = -continuous)]  # continuous model: higher score = TD, flip sign so that positive = higher in lower scores (DD-like)
rho <- ww[, .(rho_LOO = cor(primary, LOO_RN01, method = "spearman"), rho_cont = cor(primary, continuous, method = "spearman")), by = domain]
wt(rho, "S_camera_C2CP_fetal_sensitivity_concordance.tsv")
wl <- melt(ww, id.vars = c("domain", "gene_set", "primary"), variable.name = "sensitivity", value.name = "signed_sens")
wl[, sensitivity := factor(sensitivity, labels = c("Without RN01", "Continuous score (sign flipped)"))]
p22 <- ggplot(wl, aes(primary, signed_sens)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") + geom_abline(linetype = 2, colour = "#C43C39") + geom_point(size = 0.8, alpha = 0.5, colour = "#264653") +
  facet_grid(sensitivity ~ domain) + labs(x = expression("primary model: signed "-log[10]~"FDR"), y = expression("sensitivity model: signed "-log[10]~"FDR"), title = "Pathway-level stability of fetal-facing DD-vs-TD signatures (MSigDB C2 canonical pathways)") + theme_pub(11)
ggsave(file.path(FIG, "FigS22_pathway_sensitivity.png"), p22, width = 10, height = 6.5, dpi = 300, bg = "white"); ggsave(file.path(FIG, "FigS22_pathway_sensitivity.pdf"), p22, width = 10, height = 6.5, device = cairo_pdf)

# S19: previously highlighted genes - fetal-facing expression by participant
rk <- melt(rawkey[, c("participant", "cog", intersect(c("NTRK2", "GRIN2B", "HLA-DRB5", "SCGN", "TNFSF14", "BTNL9", "MEGF10", "DGCR8", "MAPK3", "FURIN"), names(rawkey))), with = FALSE], id.vars = c("participant", "cog"), variable.name = "gene", value.name = "log2CPM")
rk[, grp := ifelse(cog < 90, "DD", "TD")]
p19 <- ggplot(rk, aes(cog, log2CPM)) + geom_point(aes(colour = grp), size = 2.4) + geom_text_repel(data = rk[participant == "RN01"], aes(label = participant), size = 3) + facet_wrap(~gene, scales = "free_y", ncol = 5) + scale_colour_manual(values = pal_grp, name = "Cognitive group") +
  labs(x = "Mean Bayley-III cognitive composite", y = "Fetal-facing log2 CPM", title = "Genes highlighted in the previous version: fetal-facing expression per pregnancy") + theme_pub(11)
ggsave(file.path(FIG, "FigS19_highlighted_genes_per_pregnancy.png"), p19, width = 11, height = 5.5, dpi = 300, bg = "white"); ggsave(file.path(FIG, "FigS19_highlighted_genes_per_pregnancy.pdf"), p19, width = 11, height = 5.5, device = cairo_pdf)

# S21: cell-type localization (GSE171381 donor-balanced pseudobulk) for the top fetal-facing genes and previously highlighted genes, if the table exists
ctf <- file.path(ROOT, "data", "GSE171381_key_gene_donor_celltype_pseudobulk.tsv")  # donor x cell-type pseudobulk extracted from GEO GSE171381
if (file.exists(ctf)) {
  ct <- fread(ctf)[condition == "Control"]
  ct <- ct[, .(mean_cp10k = mean(pseudobulk_cp10k), detected = mean(detected_percent), donors = uniqueN(sample_id)), by = .(gene, cell_type)]
  wt(ct, "S_GSE171381_control_donor_balanced_expression.tsv")
  ct[, cell_type := factor(cell_type, levels = sort(unique(cell_type)))]
  p21 <- ggplot(ct, aes(cell_type, gene, fill = log2(mean_cp10k + 0.01))) + geom_tile() + geom_point(aes(size = detected), shape = 21, fill = NA, colour = "black", stroke = 0.3) + scale_fill_gradient(low = "white", high = "#B2182B", name = "log2 mean CP10K") + scale_size_continuous(range = c(0, 3), name = "% cells detected") +
    labs(x = "Cell type (GSE171381, control term placentas, donor-balanced)", y = NULL, title = "Placental cell-type expression of top fetal-facing genes and previously highlighted genes") + theme_pub(10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(FIG, "FigS21_celltype_expression.png"), p21, width = 11, height = 8, dpi = 300, bg = "white"); ggsave(file.path(FIG, "FigS21_celltype_expression.pdf"), p21, width = 11, height = 8, device = cairo_pdf)
}

writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
cat("Done. Outputs in", OUT, "\n")