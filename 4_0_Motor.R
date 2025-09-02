#!/usr/bin/env Rscript

#### imports ####
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
library(tidyverse)
library(DT)
library(gt)
library(plotly)

#### normalization and pca ####
sampleLabels <- c(targetsnd$masctransM, targetsnd$masctranscRN)
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log = TRUE)
targetsnd$motor_group <- ifelse(targetsnd$QuantQRNmean >= 90,
                                "TD",
                                "DD")
targetsnd$motor_group <- factor(targetsnd$motor_group)
group <- c(targetsnd$motor_group, targetsnd$motor_group)
group <- factor(group)
groupcont <- c(targetsnd$QuantQRNmean, targetsnd$QuantQRNmean)
groupcont <- factor(groupcont)

pca.res <- prcomp(t(log2.cpm.filtered.norm),
                  scale. = F,
                  retx = T)
pc.var <- pca.res$sdev ^ 2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var / sum(pc.var) * 100, 1)
pca.res.df <- as_tibble(pca.res$x)
compartment <- factor(
  c(rep("maternal", nrow(targetsnd)), rep("fetal", nrow(targetsnd))),
  levels = c("maternal", "fetal")
)
pca.res.df$sample      <- sampleLabels
pca.res.df$group       <- group
pca.res.df$compartment <- compartment

pca.plot <- ggplot(pca.res.df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(shape = compartment, color = group),
    size = 3.5, stroke = 1.1, alpha = 0.95
  ) +
  scale_shape_manual(values = c(maternal = 1, fetal = 16)) +
  stat_ellipse(aes(color = group), linetype = "dashed", linewidth = 0.4) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  coord_fixed() +
  theme_bw()

p1 <- pca.plot
ggplotly(pca.plot)

#### healthy altered ####

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



#### Volcano plot ####

library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(knitr)

design <- model.matrix(~ groupcont)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
ebFit <- eBayes(fit)
myTopHits <- topTable(ebFit, adjust.method = "BH", coef = 2, number = Inf, sort.by = "logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")


top_upregulated_genes <- myTopHits.df[order(-myTopHits.df$logFC), ][1:25, ]
top_downregulated_genes <- myTopHits.df[order(myTopHits.df$logFC), ][1:25, ]
top50_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
write.table(top50_genes, file = "Supp9_Mot_top50_up25down25.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

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

#### DEG table ####

results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.05)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] != 0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(
  diffGenes.df,
  extensions = c('KeyTable', "FixedHeader"),
  caption = 'Table 1: DEGs in samples',
  options = list(
    keys = TRUE,
    searchHighlight = TRUE,
    pageLength = 10,
    lengthMenu = c("10", "25", "50", "100")
  )
) %>%
  formatRound(columns = c(2:31), digits = 2)

write.table(diffGenes.df, file = "results_paper/Supp10_Mot_DEG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

#### heatmap ####

library(tidyverse)
library(gplots)
library(RColorBrewer)
library(magick)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(patchwork)
library(cowplot)
library(limma)

results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 1e-3)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] != 0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

clustRows <- hclust(as.dist(1 - cor(t(diffGenes), method = "pearson")), method = "complete")
clustColumns <- hclust(as.dist(1 - cor(diffGenes, method = "spearman")), method = "complete")
col_hc <- clustColumns
row_hc <- clustRows
k <- 2
row_cluster_id  <- cutree(row_hc, k = k)
row_cluster_lab <- factor(paste0("C", row_cluster_id), levels = paste0("C", 1:k))
cluster_pal <- c(C1 = "#F39C12", C2 = "#2E86C1")
comp_pal    <- c(maternal = "#8C564B", fetal = "#17BECF")
cog_pal     <- c(DD = "#E84A5F", TD = "#00A8E8")
row_ha <- rowAnnotation(
  Cluster = row_cluster_lab,
  col = list(Cluster = cluster_pal),
  annotation_legend_param = list(title = "Cluster")
)
nM <- nrow(targetsnd)
nF <- nrow(targetsnd)
compartment <- factor(c(rep("maternal", nM), rep("fetal", nF)), levels = c("maternal","fetal"))
cog <- factor(c(targetsnd$motor_group, targetsnd$motor_group), levels = c("DD","TD"))
col_ha <- HeatmapAnnotation(
  Compartment = compartment,
  Motor   = cog,
  col = list(
    Compartment = comp_pal,
    Motor   = cog_pal
  ),
  annotation_name_side = "left",
  gp = gpar(col = NA)
)
col_fun <- colorRamp2(c(-4, 0, 4), brewer.pal(11, "RdBu")[c(11, 6, 1)])
diffGenes_scaled <- t(scale(t(diffGenes)))
ht <- Heatmap(
  diffGenes_scaled,
  name = "Z-score",
  col  = col_fun,
  cluster_rows    = as.dendrogram(row_hc),
  cluster_columns = as.dendrogram(col_hc),
  show_row_names  = FALSE,
  show_column_names = TRUE,
  column_names_rot  = 45,
  top_annotation    = col_ha,
  right_annotation  = row_ha,
  row_split         = k,
  row_title         = paste0(" "),
  heatmap_legend_param = list(title = "Z-score"),
  border = TRUE
)

heatmap_grob <- grid.grabExpr(draw(ht))
heatmap_ggplot <- as.ggplot(heatmap_grob)
print(heatmap_ggplot)

p3 <- draw(ht, heatmap_legend_side = "right", merge_legends = TRUE)
p3

#### gostplot ####
module.assign <- cutree(clustRows, k = 2)

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

modulePick <- 2
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete")
gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

p5 <- gostplot(gost.res_down, interactive = F, capped = T)

publish_gostplot(p5)

multigp = gost(list("DD" = row.names(myModule_up),
                    "TD" = row.names(myModule_down)))

p6 <- gostplot(multigp, interactive = F, capped = T)
p6

multi_gostres2 <- gost(query = list("DD" = rownames(myModule_down),
                                    "TD" = rownames(myModule_up)),
                       organism = "hsapiens",
                       correction_method = "fdr",
                       sources = c("KEGG", "REAC", "HP", "WP"),
                       multi_query = TRUE)
p6 <- gostplot(multi_gostres2, interactive = F, capped = T)
p6

publish_gosttable(multi_gostres2,
                  highlight_terms = multi_gostres2$result[1:100,],
                  use_colors = TRUE,
                  show_columns = c("source", "term_name", "term_size"),
                  filename = "results_paper/Supp11_Motor_top100.pdf")

p7 <- publish_gostplot(p6, highlight_terms = c(
  "REAC:R-HSA-2262752",
  "REAC:R-HSA-8953897",
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

#### Fig4 ####

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
#### general table ####

hs_gsea_c2 <- msigdbr(species = "Homo sapiens",
                      category = "C2") %>% dplyr::select(gs_name, gene_symbol)
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
myGSEA.res_c2 <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df_2 <- as_tibble(myGSEA.res_c2@result)

#### plot c2 ####
library(ggplot2)
library(enrichplot)
library(ggraph)
library(clusterProfiler)
library(dplyr)
library(viridis)

top_up <- myGSEA.df_2 %>% arrange(desc(NES)) %>% slice_head(n = 25)
top_dn <- myGSEA.df_2 %>% arrange(NES)       %>% slice_head(n = 25)
top50_c2 <- bind_rows(top_up, top_dn) %>%
  mutate(
    FDR = p.adjust,
    mlogFDR = -log10(pmax(FDR, 1e-300)),
    # leading edge count from "A/B/C" string:
    leading_n = if_else(is.na(core_enrichment) | core_enrichment == "",
                        0L, str_count(core_enrichment, "/") + 1L),
    lead_pct = leading_n / setSize,
    direction = if_else(NES >= 0, "Upregulated", "Downregulated"),
    Description_clean = str_replace_all(Description, "_", " "),
    Description_wrap  = str_wrap(Description_clean, width = 45)
  ) %>%
  group_by(direction) %>%
  arrange(direction, desc(abs(NES))) %>%
  mutate(order_y = row_number()) %>%
  ungroup()

# choose a sparse set to label
to_label <- top50_c2 %>%
  group_by(direction) %>%
  slice_max(order_by = abs(NES), n = 6) %>%
  pull(Description_wrap)

p8 <- ggplot(top50_c2,
             aes(x = NES, y = fct_reorder(Description_wrap, order_y))) +
  # lollipop
  geom_segment(aes(x = 0, xend = NES,
                   yend = fct_reorder(Description_wrap, order_y),
                   color = mlogFDR)) +
  geom_point(aes(size = lead_pct, color = mlogFDR), alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "solid", color = "grey40") +
  geom_vline(xintercept = c(-1.5, 1.5), linewidth = 0.3, linetype = "dashed", color = "grey70") +
  scale_color_viridis_c(name = expression(-log[10]("FDR")), option = "D", direction = -1) +
  scale_size_continuous(name = "Leading-edge fraction *") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "GSEA results for curated gene sets - c2",
    subtitle = "Motor domain"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(size = 11),
    legend.position = "right"
  )

p8 <- p8 +
  labs(
    caption = paste0(
      "n = ", nrow(top50_c2), " gene sets; vertical lines at NES = ±1.5\n",
      "*Leading-edge fraction = proportion of genes in the set that actually drive the enrichment score.\n",
      "For example, 60% means that 60% of the genes in that pathway are enriched at the top (or bottom) of the ranked list."
    )
  )

p8

# Save the combined figure
ggsave("results_paper/Supp12_Motor_GSEA_c2.png",
       p8,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400,
       bg = "white")
ggsave("results_paper/Supp12_Motor_GSEA_c2.tiff",
       p8,
       width = 300,
       height = 300,
       units = "mm",
       dpi = 400,
       bg = "white")

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
          caption = 'Immuno signatures enriched in Motor B-III',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

top_up <- myGSEA.df_7 %>% arrange(desc(NES)) %>% slice_head(n = 25)
top_dn <- myGSEA.df_7 %>% arrange(NES)       %>% slice_head(n = 25)
top50_c7 <- bind_rows(top_up, top_dn) %>%
  mutate(
    FDR = p.adjust,
    mlogFDR = -log10(pmax(FDR, 1e-300)),
    # leading edge count from "A/B/C" string:
    leading_n = if_else(is.na(core_enrichment) | core_enrichment == "",
                        0L, str_count(core_enrichment, "/") + 1L),
    lead_pct = leading_n / setSize,
    direction = if_else(NES >= 0, "Upregulated", "Downregulated"),
    Description_clean = str_replace_all(Description, "_", " "),
    Description_wrap  = str_wrap(Description_clean, width = 45)
  ) %>%
  group_by(direction) %>%
  arrange(direction, desc(abs(NES))) %>%
  mutate(order_y = row_number()) %>%
  ungroup()

# choose a sparse set to label
to_label <- top50_c7 %>%
  group_by(direction) %>%
  slice_max(order_by = abs(NES), n = 6) %>%
  pull(Description_wrap)

p9 <- ggplot(top50_c7,
             aes(x = NES, y = fct_reorder(Description_wrap, order_y))) +
  # lollipop
  geom_segment(aes(x = 0, xend = NES,
                   yend = fct_reorder(Description_wrap, order_y),
                   color = mlogFDR)) +
  geom_point(aes(size = lead_pct, color = mlogFDR), alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "solid", color = "grey40") +
  geom_vline(xintercept = c(-1.5, 1.5), linewidth = 0.3, linetype = "dashed", color = "grey70") +
  scale_color_viridis_c(name = expression(-log[10]("FDR")), option = "D", direction = -1) +
  scale_size_continuous(name = "Leading-edge fraction *") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "GSEA results for curated gene sets - c7",
    subtitle = "Motor domain"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(size = 11),
    legend.position = "right"
  )

p9 <- p9 +
  labs(
    caption = paste0(
      "n = ", nrow(top50_c2), " gene sets; vertical lines at NES = ±1.5\n",
      "*Leading-edge fraction = proportion of genes in the set that actually drive the enrichment score.\n",
      "For example, 60% means that 60% of the genes in that pathway are enriched at the top (or bottom) of the ranked list."
    )
  )

p9
# Save the combined figure
ggsave("results_paper/Supp13_Mot_GSEA_immune_c7.png",
       p9,
       width = 500,
       height = 300,
       units = "mm",
       dpi = 400,
       bg = "white")
ggsave("results_paper/Supp13_Mot_GSEA_immune_c7.tiff",
       p9,
       width = 500,
       height = 300,
       units = "mm",
       dpi = 400,
       bg = "white")
