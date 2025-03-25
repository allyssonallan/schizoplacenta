library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
sampleLabels <- c(targetsnd$masctransM, targetsnd$masctranscRN)
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log = TRUE)

### PCA plot

library(tidyverse)
library(DT)
library(gt)
library(plotly)

#categorical
# group <- c(targetsnd$QualiQRN,targetsnd$QualiQRN)
# group <- factor(group)
#continuous
targetsnd$motor_group <- ifelse(targetsnd$QuantQRNmean >= 90,
                                "typical",
                                "altered")
targetsnd$motor_group <- factor(targetsnd$motor_group)
group <- c(targetsnd$motor_group, targetsnd$motor_group)
group <- factor(group)
groupcont <- c(targetsnd$QuantQRNmean, targetsnd$QuantQRNmean)
groupcont <- factor(groupcont)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  coord_fixed() +
  theme_bw()

p1 <- pca.plot
ggplotly(pca.plot)


### Filtered and normalized counts per million (CPM)

log2.cpm.filtered.norm.df <- log2.cpm.filtered.norm.df[!duplicated(names(log2.cpm.filtered.norm.df))]
#continuous data
healthy <- targetsnd %>%
  dplyr::filter(QuantQRNmean >= 90) %>%
  pivot_longer(
    cols = c("masctranscRN", "masctransM"),
    names_to = "type",
    values_to = "QuantQRN_value"
  )

altered <- targetsnd %>%
  dplyr::filter(QuantQRNmean < 90) %>%
  pivot_longer(
    cols = c("masctranscRN", "masctransM"),
    names_to = "type",
    values_to = "QuantQRN_value"
  )

healthy_concatenated <- paste(healthy$QuantQRN_value, collapse = " + ")
altered_concatenated <- paste(altered$QuantQRN_value, collapse = " + ")

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = eval(parse(text = paste0("(", healthy_concatenated, ")/length(healthy)"))),
                    altered.AVG = eval(parse(text = paste0("(", altered_concatenated, ")/length(altered)"))),
                    LogFC = (healthy.AVG - altered.AVG)) %>%
  mutate_if(is.numeric, round, 2)

datatable(mydata.df,
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))



### Volcano plot

library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(knitr)
#categorical

design <- model.matrix(~ groupcont)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
ebFit <- eBayes(fit)
myTopHits <- topTable(ebFit, adjust ="BH", coef=2, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID", number=4000)


top_upregulated_genes <- myTopHits.df[order(-myTopHits.df$logFC), ][1:25, ]  # Top 25 upregulated
top_downregulated_genes <- myTopHits.df[order(myTopHits.df$logFC), ][1:25, ] # Top 25 downregulated
top50_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
write.table(top50_genes, file = "Supp_9_Mot_top50_up25down25.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

df_clean <- na.omit(top50_genes)
df_clean <- df_clean[c("geneID", "logFC", "adj.P.Val")]
df_clean <- df_clean[df_clean$adj.P.Val <= 0.05,]
df_clean <- df_clean[order(-df_clean$logFC), ]
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
  theme_bw()

p2 <- vplot
ggplotly(vplot)

### DEG table

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=1e-10, lfc=3)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in samples',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:31), digits=2)

write.table(diffGenes.df, file = "10_Mot_DEG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

### heatmap

library(tidyverse)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)

myheatcolors <- rev(brewer.pal(name = "RdBu", n = 11))
clustRows <- hclust(as.dist(1 - cor(t(diffGenes), method = "pearson")), method =
                      "complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1 - cor(diffGenes, method = "spearman")), method =
                         "complete")
module.assign <- cutree(clustRows, k = 2)
module.color <- rainbow(length(unique(module.assign)), start = 0.1, end =
                          0.9)
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
                   show_row_names = FALSE,
                   show_column_names = TRUE,
                   right_annotation = row_annotation,
                   heatmap_legend_param = list(title = "Z-score"),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 10))

heatmap_grob <- grid.grabExpr(draw(heatmap))
heatmap_ggplot <- as.ggplot(heatmap_grob)
print(heatmap_ggplot)

p3 <- draw(heatmap, heatmap_legend_side = "right")

### GO - enriched in low motor

modulePick <- 1
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete")

library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
gost.res_up <- gost(rownames(myModule_up), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

### GO - enriched in normal motor category

modulePick <- 2
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete")
gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

p5 <- gostplot(gost.res_down, interactive = F, capped = T)

publish_gostplot(p5)

multigp = gost(list("altered" = row.names(myModule_up),
                    "typical" = row.names(myModule_down)))

p6 <- gostplot(multigp, interactive = F, capped = T)
p6

multi_gostres2 <- gost(query = list("altered" = rownames(myModule_down),
                                    "typical" = rownames(myModule_up)),
                       organism = "hsapiens",
                       correction_method = "fdr",
                       sources = c("KEGG", "REAC", "MIRNA","HP", "WP"),
                       multi_query = TRUE)
p6 <- gostplot(multi_gostres2, interactive = F, capped = T)
p6

publish_gosttable(multi_gostres2,
                  highlight_terms = multi_gostres2$result[1:100,],
                  use_colors = TRUE,
                  show_columns = c("source", "term_name", "term_size"),
                  filename = "Supp_Motor_top100.pdf")

p7 <- publish_gostplot(p6, highlight_terms = c(
  "MIRNA:hsa-miR-16-5p",
  "MIRNA:hsa-miR-16-5p",
  "MIRNA:hsa-miR-155-5p",
  "MIRNA:hsa-miR-30a-5p",
  "REAC:R-HSA-2262752",
  "REAC:R-HSA-8953897",
  "MIRNA:hsa-miR-92a-3p",
  "MIRNA:hsa-miR-21-5p",
  "MIRNA:hsa-miR-17-5p",
  "WP:WP2005",
  "WP:WP2004",
  "REAC:R-HSA-9824446",
  "REAC:R-HSA-5663205",
  "KEGG:05014",
  "HP:0002060",
  "HP:0000240",
  "HP:0100547",
  "HP:0011804",
  "HP:0002977",
  "HP:0003808",
  "HP:0001263",
  "HP:0003011",
  "KEGG:05012",
  "KEGG:05016",
  "HP:0012759"
))
p7

### FigMot###

library(patchwork)

combined <- ((p1 + p2 + p_table) / heatmap_ggplot) +
  plot_annotation(
    title = "Motor category",
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
  plot_annotation(title = "Motor DEGs category",
                  subtitle = "Fetal and maternal samples (n = 30)",
                  tag_levels = 'a')
combined
ggsave("results_paper/Fig4_Mot.png", combined, width = 600, height = 400, units = "mm", dpi = 400)
ggsave("results_paper/Fig4_Mot.tiff", combined, width = 600, height = 400, units = "mm", dpi = 400)
### general table

hs_gsea_c2 <- msigdbr(species = "Homo sapiens",
                      category = "C2") %>% dplyr::select(gs_name, gene_symbol)
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res_c2 <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df_2 <- as_tibble(myGSEA.res_c2@result)

## load libraries ###
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
       title = "GSEA results for curated gene sets - c2",
       subtitle = "Motor domain") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))
p8

# Save the combined figure
ggsave("results_paper/Supp7_Mot_GSEA_c2.png",
       p8,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)
ggsave("results_paper/Supp7_Mot_GSEA_c2.tiff",
       p8,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400)

### Immuno table

hs_gsea_c7 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C7") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res_c7 <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c7, verbose=FALSE)
myGSEA.df_7 <- as_tibble(myGSEA.res_c7@result)

# view results as an interactive table
datatable(myGSEA.df_7,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Immuno signatures enriched in Language B-III',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

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
       title = "GSEA results for immune curated gene sets - c7",
       subtile = "Motor domain") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))
p9
# Save the combined figure
ggsave("results_paper/Supp8_Lang_GSEA_immune_c7.png",
       p9,
       width = 500,
       height = 300,
       units = "mm",
       dpi = 400)
ggsave("results_paper/Supp8_Lang_GSEA_immune_c7.tiff",
       p9,
       width = 500,
       height = 300,
       units = "mm",
       dpi = 400)

### load libraries ###
library(ggplot2)
library(enrichplot)
library(ggraph)
library(clusterProfiler)

# Sort for the top upregulated (highest ES) and downregulated (lowest NES) gene sets
top_upregulated <- myGSEA.df_7[order(-myGSEA.df_7$NES), ][1:25, ]  # Top 25 upregulated
top_downregulated <- myGSEA.df_7[order(myGSEA.df_7$NES), ][1:25, ] # Top 25 downregulated
top50_immune <- rbind(top_upregulated, top_downregulated)

ggplot(top50_immune, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(x = "Gene Set", y = "Normalized Enrichment Score (NES)", fill = "p.adjust",
       title = "GSEA Results for Immune curated gene sets") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
