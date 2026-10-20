# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 04_supplementary_figures.R: Figures S4-S15 and S19 (same slots as the first submission) + S20-S22 (new) + a rebuttal-only figure
# Source via run_all.R from the project root.
# =============================================================================
clean <- function(x) str_trunc(str_to_sentence(gsub("_", " ", x)), 62)
sv <- function(p, fn, w, h) { ggsave(file.path(FIG, paste0(fn, ".png")), p, width = w, height = h, dpi = 300, bg = "white"); ggsave(file.path(FIG, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf) }

# ---- GSEA lollipops (S5/S8/S11 C2:CP; S6/S9/S12 C7): top n sets higher and lower in DD per contrast, by FDR then |NES|
lolli_gsea <- function(ce, fn, title, n = 15, w = 9, h = 14) {
  ce <- ce[order(contrast, FDR, -abs(NES))]
  ce <- rbind(ce[NES > 0, head(.SD, n), by = contrast], ce[NES < 0, head(.SD, n), by = contrast])
  ce[, label := clean(gene_set)]; ce[, sig := ifelse(FDR < 0.05, "FDR < 0.05", "FDR >= 0.05")]
  ce[, contrast := factor(ifelse(contrast == "fetal", "Fetal-facing DD vs TD", "Maternal-facing DD vs TD"), levels = c("Fetal-facing DD vs TD", "Maternal-facing DD vs TD"))]
  ce[, key := paste(contrast, label)]
  p <- ggplot(ce, aes(NES, reorder(key, NES))) + geom_segment(aes(x = 0, xend = NES, yend = reorder(key, NES), colour = Direction)) + geom_point(aes(size = NGenes, colour = Direction, alpha = sig)) +
    facet_wrap(~contrast, scales = "free_y", ncol = 1) + scale_y_discrete(labels = function(x) sub("^[^ ]+-facing DD vs TD ", "", x)) +
    scale_colour_manual(values = c(Up = "#C43C39", Down = "#2C7FB8"), name = "Direction in DD") + scale_alpha_manual(values = c(`FDR < 0.05` = 1, `FDR >= 0.05` = 0.35), name = NULL) + scale_size_continuous(range = c(1.5, 5), name = "Genes in set") +
    geom_vline(xintercept = 0, colour = "grey60") + labs(x = "Normalized enrichment score (pre-ranked GSEA; positive = higher in DD)", y = NULL, title = str_wrap(title, 75)) + theme_pub(11) + theme(axis.text.y = element_text(size = 9), legend.position = "bottom", legend.box = "vertical")
  sv(p, fn, w, h)
}
# ---- ORA Manhattan-like plot (S4/S7/S10): nominal fetal-facing DEGs higher / lower in DD across sources
ora_plot <- function(d, fn) {
  x <- ora[domain == DOM[d] & overlap > 0]; x[, source := factor(source, levels = c("GO:BP", "GO:CC", "GO:MF", "HP", "KEGG", "REAC", "PID", "BIOCARTA"))]
  x[, direction := factor(direction, levels = c("higher in DD", "lower in DD"))]
  set.seed(1); x[, xpos := as.numeric(source) + runif(.N, -0.35, 0.35)]
  lab <- x[padj < 0.05][order(padj)][, head(.SD, 3), by = .(direction, source)]; lab[, term_lab := str_trunc(str_to_sentence(gsub("_", " ", sub("^(GOBP|GOCC|GOMF|HP|KEGG|REACTOME|WP|PID|BIOCARTA)_", "", term))), 45)]
  pal <- c(`GO:BP` = "#FF9900", `GO:CC` = "#109618", `GO:MF` = "#DC3912", HP = "#990099", KEGG = "#DD4477", REAC = "#3366CC", WP = "#0099C6", PID = "#66AA00", BIOCARTA = "#B82E2E")
  x[, sig := padj < 0.05]
  p <- ggplot(x, aes(xpos, -log10(p), colour = source)) + geom_point(aes(size = overlap, alpha = sig)) + scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.25), guide = "none") +
    geom_text_repel(data = lab, aes(label = term_lab), size = 2.7, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4, show.legend = FALSE) +
    facet_wrap(~direction, ncol = 1, labeller = labeller(direction = function(v) paste("Nominal fetal-facing DEGs", v))) + scale_x_continuous(breaks = seq_along(levels(x$source)), labels = levels(x$source)) + scale_colour_manual(values = pal, guide = "none") + scale_size_continuous(range = c(1, 5), name = "DEGs in term") +
    labs(x = NULL, y = expression(-log[10]~italic(P)~"(hypergeometric; solid and labelled: adjusted "~italic(P)~"< 0.05 within source)"), title = str_wrap(paste(DOM[d], "domain: over-representation of nominal fetal-facing DEGs (p < 0.05, |log2FC| >= 1) across GO, HPO, KEGG, Reactome, PID and BioCarta"), 95)) + theme_pub(11) + theme(legend.position = "bottom")
  sv(p, fn, 11, 10)
}
for (d in names(DOM)) {
  D <- DOM[d]; s <- switch(D, Cognitive = c("S4", "S5", "S6"), Language = c("S7", "S8", "S9"), Motor = c("S10", "S11", "S12"))
  base <- gsea[domain == D & model == "primary" & contrast %in% c("fetal", "maternal")]
  ora_plot(d, paste0("Fig", s[1], "_", D, "_ORA"))
  lolli_gsea(base[collection == "C2CP"], paste0("Fig", s[2], "_", D, "_GSEA_C2CP"), paste(D, "domain - GSEA, MSigDB C2 canonical pathways, top 15 per direction"))
  if ("C7" %in% base$collection) lolli_gsea(base[collection == "C7"], paste0("Fig", s[3], "_", D, "_GSEA_C7"), paste(D, "domain - GSEA, MSigDB C7 immunologic signatures, top 15 per direction"))
}

# ---- S13-S15: face-specific contrasts of the joint model
face_fig <- function(d, fn) {
  t <- prim_tabs[domain == d]; w <- dcast(t[, .(gene, contrast, logFC, P.Value)], gene ~ contrast, value.var = c("logFC", "P.Value"))
  vol <- function(cc, ttl) { x <- t[contrast == cc]; x[, sig := P.Value < 0.05 & abs(logFC) >= 1]; ggplot(x, aes(logFC, -log10(P.Value))) + geom_point(aes(colour = sig), size = 0.6, alpha = 0.6) + scale_colour_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey75"), guide = "none") +
    geom_text_repel(data = x[order(P.Value)][1:8], aes(label = gene), size = 2.8, max.overlaps = Inf) + geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") + geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") + labs(title = ttl, x = "log2 fold change (DD vs TD)", y = expression(-log[10]~italic(P))) + theme_pub(11) }
  pa <- vol("fetal", "Fetal-facing DD vs TD"); pb <- vol("maternal", "Maternal-facing DD vs TD"); pc <- vol("interaction", "Face x outcome interaction")
  r <- cor(w$logFC_fetal, w$logFC_maternal, method = "spearman")
  pd <- ggplot(w, aes(logFC_maternal, logFC_fetal)) + geom_point(size = 0.5, alpha = 0.3, colour = "#264653") + geom_abline(linetype = 2, colour = "#C43C39") + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    labs(title = sprintf("Concordance of DD-TD effects between faces (Spearman rho = %.2f)", r), x = "Maternal-facing log2FC", y = "Fetal-facing log2FC") + theme_pub(11)
  fig <- (pa | pb | pc) / pd + plot_annotation(tag_levels = "a", title = paste(DOM[d], "domain - face-specific contrasts from the joint adjusted model"))
  sv(fig, fn, 12, 8)
}
face_fig("cog", "FigS13_Cognitive_faces"); face_fig("lang", "FigS14_Language_faces"); face_fig("mot", "FigS15_Motor_faces")

# ---- S19 (same slot as the first submission): fetal- versus maternal-facing DD-vs-TD effects, cognitive domain, candidate genes labelled
CAND <- c("BTNL9", "HLA-DRB5", "GRIN2B", "GRIN2A", "CACNA1F", "KCNQ3", "EGR2", "TECTA", "IL7", "TREM2", "ABT1", "CRHBP", "SERINC2", "DKK2", "ABCA7", "CSMD1", "MTHFD2", "NTRK2", "SCGN", "TNFSF14")
t <- prim_tabs[domain == "cog"]; w <- dcast(t[, .(gene, contrast, logFC, P.Value)], gene ~ contrast, value.var = c("logFC", "P.Value"))
w[, cls := fifelse(P.Value_fetal < 0.05 & abs(logFC_fetal) >= 1 & !(P.Value_maternal < 0.05 & abs(logFC_maternal) >= 1), "fetal-facing only",
          fifelse(P.Value_maternal < 0.05 & abs(logFC_maternal) >= 1 & !(P.Value_fetal < 0.05 & abs(logFC_fetal) >= 1), "maternal-facing only",
          fifelse(P.Value_fetal < 0.05 & abs(logFC_fetal) >= 1, "both faces", "n.s.")))]
w[, cls := factor(cls, levels = c("fetal-facing only", "maternal-facing only", "both faces", "n.s."))]
p19 <- ggplot(w[order(-as.integer(cls))], aes(logFC_maternal, logFC_fetal)) + geom_abline(linetype = 2, colour = "grey60") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = cls, alpha = cls), size = 1) + scale_colour_manual(values = c(`fetal-facing only` = "#17BECF", `maternal-facing only` = "#8C564B", `both faces` = "#7B3294", n.s. = "grey80"), name = "Nominal DEG (p < 0.05, |log2FC| >= 1)") + scale_alpha_manual(values = c(0.9, 0.9, 0.9, 0.3), guide = "none") +
  geom_point(data = w[gene %in% CAND], colour = "black", shape = 21, size = 2.6, fill = NA) + geom_text_repel(data = w[gene %in% CAND], aes(label = gene), size = 3, fontface = "italic", max.overlaps = Inf, min.segment.length = 0, box.padding = 0.5) +
  labs(x = "Maternal-facing log2 fold change (DD vs TD)", y = "Fetal-facing log2 fold change (DD vs TD)", title = str_wrap("Cognitive domain: fetal-facing versus maternal-facing DD-vs-TD effects (covariate-adjusted, pregnancy-paired model)", 90)) + theme_pub(11) + theme(legend.position = "bottom")
sv(p19, "FigS19_cognitive_fetal_vs_maternal_effects", 10, 9)
wt(w[gene %in% CAND][order(P.Value_fetal)], "S_candidate_genes_fetal_vs_maternal_cog.tsv")

# ---- S20 (new): published placental gene sets: (a) GSEA NES, primary and leave-one-out; (b) non-directional fry mixed P
short <- c("SCZ TWAS Ursini 2023" = "SCZ placental TWAS (Ursini 2023)", "SCZ prioritized Ursini 2023" = "SCZ placenta-prioritized (Ursini 2023)",
           "SCZ GWAS placenta-implicated (Ursini 2023 S1; 56)" = "SCZ GWAS placenta-implicated (Ursini 2023)", "PlacPRS1 schizophrenia genes (Ursini 2021; 66)" = "PlacPRS1 genes (Ursini 2021)",
           "SCZ CpG annotation Cilleros 2025" = "SCZ CpG-annotated (Cilleros-Portet 2025)", "SCZ eQTM Cilleros 2025" = "SCZ eQTM targets (Cilleros-Portet 2025)",
           "ADHD TWAS Ursini 2023" = "ADHD placental TWAS (Ursini 2023)", "COVID up Gao 2023" = "Up in COVID placenta (Gao 2023)", "COVID down Gao 2023" = "Down in COVID placenta (Gao 2023)")
prep_ext <- function(x) {
  x <- copy(x)[contrast %in% c("fetal", "maternal") & model %in% c("primary", "LOO_RN01")]
  x[, contrast := factor(ifelse(contrast == "fetal", "Fetal-facing", "Maternal-facing"), levels = c("Fetal-facing", "Maternal-facing"))]
  x[, model := factor(ifelse(model == "primary", "All 15 pregnancies", "Without RN01"), levels = c("All 15 pregnancies", "Without RN01"))]
  x[, gene_set := ifelse(gene_set %in% names(short), short[gene_set], gene_set)]; x[, gene_set := factor(gene_set, levels = rev(sort(unique(gene_set))))]
}
ce <- prep_ext(gsea[collection == "external"]); ce[, sig := ifelse(FDR < 0.05, "FDR < 0.05", "FDR >= 0.05")]
p20a <- ggplot(ce, aes(NES, gene_set, colour = Direction, shape = model)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(aes(size = NGenes, alpha = sig), position = position_dodge(width = 0.6)) + facet_grid(domain ~ contrast) +
  scale_colour_manual(values = c(Up = "#C43C39", Down = "#2C7FB8"), name = "Direction in DD") + scale_shape_manual(values = c(16, 1), name = "Model") + scale_alpha_manual(values = c(`FDR < 0.05` = 1, `FDR >= 0.05` = 0.35), name = NULL) + scale_size_continuous(range = c(1.5, 5), name = "Genes tested") +
  labs(x = "Normalized enrichment score (pre-ranked GSEA; positive = higher in DD)", y = NULL, title = "Directional: pre-ranked GSEA") + theme_pub(11) + theme(axis.text.y = element_text(size = 9.5), legend.position = "bottom", legend.box = "vertical")
fe <- prep_ext(fr[collection == "external"])
p20b <- ggplot(fe, aes(-log10(PValue.Mixed), gene_set, shape = model)) + geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey75") + geom_point(aes(size = NGenes), position = position_dodge(width = 0.6), colour = "#264653", alpha = 0.9) + facet_grid(domain ~ contrast) +
  scale_shape_manual(values = c(16, 1), name = "Model") + scale_size_continuous(range = c(1.5, 5), name = "Genes tested") + labs(x = expression(-log[10]~italic(P)~"(non-directional rotation test, mixed statistic)"), y = NULL, title = "Non-directional: any change in either direction (limma fry, mixed statistic)") + theme_pub(11) + theme(axis.text.y = element_text(size = 9.5), legend.position = "bottom", legend.box = "vertical")
sv((p20a / p20b) + plot_annotation(tag_levels = "a", title = "Published placental gene sets in the DD-vs-TD contrasts"), "FigS20_external_gene_sets", 11, 16)

# ---- S21 (new): cell-type localization (GSE171381 donor-balanced pseudobulk) of the top fetal-facing genes and retained candidate genes
ctf <- file.path(ROOT, "data", "GSE171381_key_gene_donor_celltype_pseudobulk.tsv")
if (file.exists(ctf)) {
  topfet <- unique(c(prim_tabs[contrast == "fetal"][order(P.Value), head(.SD, 15), by = domain]$gene, "BTNL9", "HLA-DRB5", "GRIN2B", "GRIN2A", "CACNA1F"))
  ct <- fread(ctf)[condition == "Control" & gene %in% topfet]
  ct <- ct[, .(mean_cp10k = mean(pseudobulk_cp10k), detected = mean(detected_percent), donors = uniqueN(sample_id)), by = .(gene, cell_type)]
  wt(ct, "S_GSE171381_control_donor_balanced_expression.tsv")
  ct[, cell_type := factor(cell_type, levels = sort(unique(cell_type)))]
  p21 <- ggplot(ct, aes(cell_type, gene, fill = log2(mean_cp10k + 0.01))) + geom_tile() + geom_point(aes(size = detected), shape = 21, fill = NA, colour = "black", stroke = 0.3) + scale_fill_gradient(low = "white", high = "#B2182B", name = "log2 mean CP10K") + scale_size_continuous(range = c(0, 3), name = "% cells detected") +
    labs(x = "Cell type (GSE171381, control term placentas, donor-balanced)", y = NULL, title = "Placental cell-type expression of the top fetal-facing DD-vs-TD genes and candidate genes") + theme_pub(10) + theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text.y = element_text(face = "italic"))
  sv(p21, "FigS21_celltype_expression", 11, 7.5)
}

# ---- S22 (new): pathway-level stability (GSEA NES, C2:CP fetal): primary vs leave-one-out and vs continuous score
cc <- gsea[collection == "C2CP" & contrast == "fetal" & model %in% c("primary", "LOO_RN01", "continuous")]
ww <- dcast(cc, domain + gene_set ~ model, value.var = "NES")
rho <- ww[, .(rho_LOO = cor(primary, LOO_RN01, method = "spearman", use = "complete.obs"), rho_cont = cor(primary, continuous, method = "spearman", use = "complete.obs")), by = domain]
wt(rho, "S_GSEA_C2CP_fetal_sensitivity_concordance.tsv")
wl <- melt(ww, id.vars = c("domain", "gene_set", "primary"), variable.name = "sensitivity", value.name = "NES_sens")
wl[, sensitivity := factor(sensitivity, levels = c("LOO_RN01", "continuous"), labels = c("Without RN01", "Continuous score"))]
p22 <- ggplot(wl, aes(primary, NES_sens)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") + geom_abline(linetype = 2, colour = "#C43C39") + geom_point(size = 0.8, alpha = 0.5, colour = "#264653") +
  facet_grid(sensitivity ~ domain) + labs(x = "Primary model: NES", y = "Sensitivity model: NES", title = str_wrap("Pathway-level stability of the fetal-facing DD-vs-TD signatures (pre-ranked GSEA, MSigDB C2 canonical pathways)", 90)) + theme_pub(11)
sv(p22, "FigS22_pathway_sensitivity", 10, 6.5)

# ---- Rebuttal-only figure: genes highlighted in the first submission, fetal-facing expression per pregnancy
rk <- melt(rawkey[, c("participant", "cog", intersect(c("NTRK2", "GRIN2B", "HLA-DRB5", "SCGN", "TNFSF14", "BTNL9", "MEGF10", "DGCR8", "MAPK3", "FURIN"), names(rawkey))), with = FALSE], id.vars = c("participant", "cog"), variable.name = "gene", value.name = "log2CPM")
rk[, grp := ifelse(cog < 90, "DD", "TD")]
prb <- ggplot(rk, aes(cog, log2CPM)) + geom_point(aes(colour = grp), size = 2.4) + geom_text_repel(data = rk[participant == "RN01"], aes(label = participant), size = 3) + facet_wrap(~gene, scales = "free_y", ncol = 5) + scale_colour_manual(values = pal_grp, name = "Cognitive group") +
  labs(x = "Mean Bayley-III cognitive composite", y = "Fetal-facing log2 CPM", title = "Genes highlighted in the first submission: fetal-facing expression per pregnancy") + theme_pub(11)
sv(prb, "Rebuttal_highlighted_genes_per_pregnancy", 11, 5.5)

writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
cat("Done. Outputs in", OUT, "\n")
