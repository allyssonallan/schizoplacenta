### Load libraries ###

library(dplyr)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(sva)
library(limma)
library(DT)
library(gt)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplotify)
library(ggplot2)
library(grid)
library(gridExtra)
library(patchwork)
library(gtable)

### set env ###
getwd()
setwd(here::here("Documents/Dani/Transcriptoma/"))

### load meta and pseudoalignment data ###
targets <- data.table::fread("0_TableCode.tsv")
# Metadata validation
required_cols <- c("mascgenRN", "masctransM1", "masctranscRN1", "masctransM", "masctranscRN")
missing_cols <- setdiff(required_cols, names(targets))
if (length(missing_cols) > 0) stop(paste("Metadata is missing required columns:", paste(missing_cols, collapse=", ")))
# Validate path and sample columns
path_cols <- required_cols[-1]
for (col in path_cols) {
  if (!is.character(targets[[col]])) stop(paste("Column", col, "must be character"))
  if (any(is.na(targets[[col]]) | targets[[col]] == "")) stop(paste("Column", col, "contains empty or NA values"))
}
# Validate sex column if present
if ("sex" %in% names(targets)) {
  invalid_sex <- setdiff(unique(targets$sex), c("M", "F"))
  if (length(invalid_sex) > 0) stop(paste("Invalid values in 'sex' column:", paste(invalid_sex, collapse=", ")))
}
# Validate numeric metadata columns if present
num_cols <- intersect(c("QuantCRN", "QuantLRN", "QuantQRN"), names(targets))
for (col in num_cols) {
  if (!is.numeric(targets[[col]])) stop(paste("Column", col, "must be numeric"))
  if (any(is.na(targets[[col]]))) stop(paste("Column", col, "contains NA values"))
}
targetsnd <- targets[!duplicated(targets$mascgenRN), ]
pathM <- file.path(targetsnd$masctransM1, "abundance.h5")
pathRN <- file.path(targetsnd$masctranscRN1, "abundance.h5")

Tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(
  c(pathM, pathRN),
  type = "kallisto",
  tx2gene = Tx,
  txOut = FALSE,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE
)

### proceed with normalization ###

sampleLabels <- c(targetsnd$masctransM, targetsnd$masctranscRN)
myDGEList <- DGEList(Txi_gene$counts)
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm > 1) >= 5
myDGEList.filtered <- myDGEList[keepers,]
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log = TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- tidyr::pivot_longer(
  log2.cpm.filtered.norm.df,
  cols = m_S24:f_S11,
  names_to = "samples",
  values_to = "expression"
)

p <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x = samples, y = expression, fill = samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 95,
               size = 10,
               color = "black",
               show.legend = FALSE) +
  labs(y = "log2 expression", x = "sample",
       title = "Log2 Counts per Million (CPM)",
       subtitle = "filtered, TMM normalized") +
  theme_bw() +
  coord_flip()

p

#### SEX CATEGORY from here to the end ####
#### pca plot for sex category ####

group <- c(targetsnd$sex,targetsnd$sex)
group <- factor(group)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale. = F, retx = T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pca.res.df <- as_tibble(pca.res$x)

pca.plot <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, label = sampleLabels, color = group) +
  geom_point(size = 4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title = "") +
  coord_fixed() +
  theme_bw()

p1 <- pca.plot
ggplotly(pca.plot)

### split F and M ###

log2.cpm.filtered.norm.df <- log2.cpm.filtered.norm.df[!duplicated(names(log2.cpm.filtered.norm.df))]
male <- c(targetsnd$masctranscRN[targetsnd$sex == "M"], targetsnd$masctransM[targetsnd$sex == "M"])
female <- c(targetsnd$masctranscRN[targetsnd$sex == "F"], targetsnd$masctransM[targetsnd$sex == "F"])

# Compute group averages safely without eval/parse
mydata.df <- log2.cpm.filtered.norm.df %>%
  mutate(
    male.AVG   = rowMeans(select(., all_of(male))),
    female.AVG = rowMeans(select(., all_of(female))),
    LogFC      = male.AVG - female.AVG
  ) %>%
  mutate_if(is.numeric, round, 2)

### 18,972 genes table ###

datatable(mydata.df,
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))


### DEG analysis ###

group <- c(targetsnd$sex,targetsnd$sex)
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)

fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(M - F, levels = design)

### calculate the logFC and plot the adjusted p-values ###
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, sort.by="logFC")

myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

top_upregulated_genes <- myTopHits.df[order(-myTopHits.df$logFC), ][1:25, ]  # Top 25 upregulated
top_downregulated_genes <- myTopHits.df[order(myTopHits.df$logFC), ][1:25, ] # Top 25 downregulated
top50_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
write.table(top50_genes, file = "1_Sex_top50_25up25down", sep = "\t", row.names = FALSE, col.names = TRUE)


df_clean <- na.omit(top_upregulated_genes)
df_clean <- df_clean[c("geneID", "logFC", "adj.P.Val")]
deg_table_grob <- tableGrob(df_clean, rows = NULL, theme = ttheme_minimal(base_size = 8))
n_cols <- ncol(deg_table_grob)
deg_table_grob <- gtable_add_grob(deg_table_grob,
                                  grobs = segmentsGrob(x0 = unit(0, "npc"),
                                                       y0 = unit(0, "npc"),
                                                       x1 = unit(1, "npc"),
                                                       y1 = unit(0, "npc"),
                                                       gp = gpar(lwd = 2)),
                                  t = 1, l = 1, r = n_cols)
n_rows <- nrow(deg_table_grob)
deg_table_grob <- gtable_add_grob(deg_table_grob,
                                  grobs = segmentsGrob(x0 = unit(0, "npc"),
                                                       y0 = unit(0, "npc"),
                                                       x1 = unit(1, "npc"),
                                                       y1 = unit(0, "npc"),
                                                       gp = gpar(lwd = 2)),
                                  t = n_rows, l = 1, r = n_cols)
p_table <- as.ggplot(deg_table_grob)
print(p_table)

vplot <- ggplot(myTopHits.df) +
  aes(
    y = -log10(adj.P.Val),
    x = logFC,
    text = paste("Symbol:", geneID)
  ) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "longdash",
    colour = "grey",
    size = 1
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "longdash",
    colour = "grey",
    size = 1
  ) +
  geom_vline(
    xintercept = -1,
    linetype = "longdash",
    colour = "grey",
    size = 1
  ) +
#  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
#  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(
    title = "",
    subtitle = ""
  ) +
  theme_bw()

p2 <- vplot
ggplotly(vplot)

### DEG table plot ###

results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 1e-3, lfc = 1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] != 0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in samples',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = c(2:11), digits = 2)

write.table(diffGenes.df, file = "2_Sex_DEG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

### heatmap

clustRows <- hclust(as.dist(1 - cor(t(diffGenes), method = "pearson")), method = "complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1 - cor(diffGenes, method = "spearman")), method = "complete")
module.assign <- cutree(clustRows, k = 2)
module.color <- rainbow(length(unique(module.assign)), start = 0.1, end = 0.9)
module.color <- module.color[as.vector(module.assign)]


# Prepare the heatmap data
diffGenes_scaled <- t(scale(t(diffGenes)))
row_dend <- as.dendrogram(clustRows)
col_dend <- as.dendrogram(clustColumns)
cluster_labels <- ifelse(module.assign == 1, "1", "2")
annotation_df <- data.frame(Gene = rownames(diffGenes), Cluster = cluster_labels)
row_annotation <- rowAnnotation(Cluster = cluster_labels, col = list(Cluster = c("1" = "red", "2" = "blue")))

heatmap <- Heatmap(diffGenes_scaled,
                   name = "Expression",
                   cluster_rows = row_dend,
                   cluster_columns = col_dend,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   right_annotation = row_annotation,
                   heatmap_legend_param = list(title = "Z-score"),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 10))


heatmap_grob <- grid.grabExpr(draw(heatmap))
heatmap_ggplot <- as.ggplot(heatmap_grob)
print(heatmap_ggplot)

p3 <- draw(heatmap, heatmap_legend_side = "right")


### Fig1Sex ####

combined <- ((p1 + p2 + p_table) / heatmap_ggplot) +
  plot_annotation(
    title = "Sex category",
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(size = 14),
      plot.tag = element_text(size = 14),
      plot.margin = margin(5, 5, 5, 5)  # Reduce margin to tighten the layout
    )
  ) +
  plot_layout(
    widths = c(0.5, 0.5, 0.5),
    heights = c(1, 1)
#    guides = "collect"  # Align legends
  )

combined <- combined | p7 +
  plot_layout(ncol = 1)

combined = combined +
  plot_annotation(title = "Sex DEGs category",
                  subtitle = "Fetal and maternal samples (n = 30)",
                  tag_levels = 'a')
combined
ggsave("results_paper/4_Sex.png", combined, width = 600, height = 300, units = "mm", dpi = 400)
ggsave("results_paper/4_Sex.tiff", combined, width = 600, height = 300, units = "mm", dpi = 400)



# Suppose p1, p2, p_table, and heatmap_ggplot are already created

# Two-by-two layout:
# First row: p1, p2
# Second row: p_table, heatmap_ggplot
combined <- (p1 | p2) /
  (p_table | heatmap_ggplot) | p8 +
  plot_annotation(
    title = "",
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(size = 14),
      plot.tag = element_text(size = 14),
      plot.margin = margin(5, 5, 5, 5)
    )
  ) +
  plot_layout(
  heights = c(0.8, 1.2),
  widths = c(1, 1),
  )

combined

# Save the combined figure
ggsave("results_paper/Fig1_Sex.png",
       combined,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)
ggsave("results_paper/Fig1_Sex.tiff",
       combined,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)

### Enrichment analysis ###
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

modulePick <- 2
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_up <- hclust(as.dist(1 - cor(t(myModule_up), method = "pearson")), method = "complete")
gost.res_up <- gost(rownames(myModule_up), organism = "hsapiens", correction_method = "fdr")
p4 <- gostplot(gost.res_up, interactive = F, capped = T)
publish_gostplot(p4)

modulePick <- 1
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_down <- hclust(as.dist(1 - cor(t(myModule_down), method = "pearson")), method = "complete")
gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
p5 <- gostplot(gost.res_down, interactive = F, capped = T)
p5
publish_gostplot(p5, highlight_terms = c("HP:0001450", "GO:0032452",
                                         "REAC:R-HSA-3214842", "HP:0000027"))

multigp = gost(list("up-regulated" = row.names(myModule_up),
                    "down-regulated" = row.names(myModule_down)))
p6 <- gostplot(multigp, interactive = F, capped = T)

p7 <- publish_gostplot(p6, highlight_terms = c("HP:0001450", "GO:0032452",
                                         "REAC:R-HSA-3214842", "HP:0000027"))

# Save the combined figure
ggsave("results_paper/Fig2_SexGostplot.png",
       p7,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)
ggsave("results_paper/Fig2_SexGostplot.tiff",
       p7,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)

### GSEA ###
hs_gsea_c2 <- msigdbr(species = "Homo sapiens",
                      category = "C2") %>% dplyr::select(gs_name, gene_symbol)
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

myGSEA.res_c2 <- GSEA(mydata.gsea, TERM2GENE = hs_gsea_c2, verbose = FALSE)
myGSEA.df_2 <- as_tibble(myGSEA.res_c2@result)

datatable(myGSEA.df_2,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in Sex differences',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>% formatRound(columns=c(3:10), digits=2)

hs_gsea_c7 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C7") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature


mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
myGSEA.res_c7 <- GSEA(mydata.gsea, TERM2GENE = hs_gsea_c7, verbose = FALSE)
myGSEA.df_7 <- as_tibble(myGSEA.res_c7@result)

# view results as an interactive table
datatable(myGSEA.df_7,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Immuno signatures enriched in Sex differences',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = c(3:10), digits = 2)

### load libraries ###
library(ggplot2)
library(enrichplot)
library(ggraph)
library(clusterProfiler)

# Sort for the top upregulated (highest NES) and downregulated (lowest NES) gene sets
top_upregulated <- myGSEA.df_2[order(-myGSEA.df_2$NES), ][1:25, ]  # Top 25 upregulated
top_downregulated <- myGSEA.df_2[order(myGSEA.df_2$NES), ][1:25, ] # Top 25 downregulated
top50_c2 <- rbind(top_upregulated, top_downregulated)

library(dplyr)
library(viridis)  # for viridis color scales

top50_c2 <- top50_c2 %>%
  mutate(logP = -log10(pvalue))

p8 <- ggplot(top50_c2, aes(x = reorder(Description, NES), y = NES)) +
  geom_segment(aes(xend = Description, y = 0, yend = NES, color = logP), size = 1) +
  geom_point(aes(size = setSize, color = logP), alpha = 0.7) +
  coord_flip() +
  scale_color_viridis(option = "D", direction = -1) +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = "Gene Set",
       y = "Normalized Enrichment Score (NES)",
       color = "-log10(p-value)",
       size = "Set size",
       title = "GSEA results for curated gene sets - c2") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))

# Save the combined figure
ggsave("results_paper/Supp1_Sex_GSEA_c2.png",
       p8,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)
ggsave("results_paper/Supp1_Sex_GSEA_c2.tiff",
       p8,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)


# Sort for the top upregulated (highest NES) and downregulated (lowest NES) gene sets
top_upregulated <- myGSEA.df_7[order(-myGSEA.df_7$NES), ][1:25, ]  # Top 25 upregulated
top_downregulated <- myGSEA.df_7[order(myGSEA.df_7$NES), ][1:25, ] # Top 25 downregulated
top50_immune <- rbind(top_upregulated, top_downregulated)
top50_immune <- top50_immune %>%
  mutate(logP = -log10(pvalue))

p9 <- ggplot(top50_immune, aes(x = reorder(Description, NES), y = NES)) +
  geom_segment(aes(xend = Description, y = 0, yend = NES, color = logP), size = 1) +
  geom_point(aes(size = setSize, color = logP), alpha = 0.7) +
  coord_flip() +
  scale_color_viridis(option = "D", direction = -1) +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = "Gene Set",
       y = "Normalized Enrichment Score (NES)",
       color = "-log10(p-value)",
       size = "Set size",
       title = "GSEA results for immune curated gene sets - c7") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))

# Save the combined figure
ggsave("results_paper/Supp2_Sex_GSEA_immune_c7.png",
       p9,
       width = 500,
       height = 300,
       units = "mm",
       dpi = 400)
ggsave("results_paper/Supp2_Sex_GSEA_immune_c7.tiff",
       p9,
       width = 500,
       height = 300,
       units = "mm",
       dpi = 400)

