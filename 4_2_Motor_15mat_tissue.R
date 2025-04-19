
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
#sampleLabels <- c(targetsnd$masctransM, targetsnd$masctranscRN)
sampleLabels <- targetsnd$masctransM1
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

### PCA plot

library(tidyverse)
library(DT)
library(gt)
library(plotly)

#categorical
# group <- c(targetsnd$QualiQRN,targetsnd$QualiQRN)
# group <- factor(group)
#continuous

## samples transformation
targets <- read_tsv("0_TableCode.tsv")
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
# Validate numeric metadata columns if present
num_cols <- intersect(c("QuantCRN", "QuantLRN", "QuantQRN"), names(targets))
for (col in num_cols) {
  if (!is.numeric(targets[[col]])) stop(paste("Column", col, "must be numeric"))
  if (any(is.na(targets[[col]]))) stop(paste("Column", col, "contains NA values"))
}
targetsmean <- targets %>%
  group_by(masctransM1) %>%
  mutate(QuantCRNmean = mean(QuantCRN, na.rm = TRUE),
         QuantLRNmean = mean(QuantLRN, na.rm = TRUE),
         QuantQRNmean = mean(QuantQRN, na.rm = TRUE))

targetsnd <- targetsmean[!duplicated(targetsmean$masctransM1), ]

targetsnd$motor_group <- ifelse(targetsnd$QuantQRNmean >= 90,
                                "typical",
                                "altered")
targetsnd$motor_group <- factor(targetsnd$motor_group)
group <- targetsnd$motor_group
group <- factor(group)
groupcont <- targetsnd$QuantQRNmean
groupcont <- factor(groupcont)


## process transcripts
path <- file.path(sampleLabels, "abundance.h5")
Tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

myDGEList <- DGEList(Txi_gene$counts)



cpm <- cpm(myDGEList)
keepers <- rowSums(cpm > 1) >= 5
myDGEList.filtered <- myDGEList[keepers,]
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log = TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = AAL_materno_S24:SLS_materno_S27,
                                                names_to = "samples",
                                                values_to = "expression")
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
       subtitle = "filtered, TMM normalized",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

plot(p)


pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, label = sampleLabels, color = group) +
  geom_point(size = 4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  coord_fixed() +
  theme_bw()

p1_m <- pca.plot
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

healthy_concatenated <- healthy[!duplicated(healthy$masctransM1),]
healthy_concatenated <- paste(healthy_concatenated$masctransM1, collapse = " + ")
altered_concatenated <- altered[!duplicated(altered$masctransM1),]
altered_concatenated <- paste(altered_concatenated$masctransM1, collapse = " + ")


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


write.table(diffGenes.df, file = "All_Maternal_Mot_DEG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

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


top_upregulated_genes <- myTopHits.df[order(-myTopHits.df$logFC), ][1:10, ]  # Top 25 upregulated
top_downregulated_genes <- myTopHits.df[order(myTopHits.df$logFC), ][1:10, ] # Top 25 downregulated
top50_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
write.table(top50_genes, file = "top50Motor_maternal.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


df_clean <- na.omit(top50_genes)
df_clean <- df_clean[c("geneID", "logFC", "P.Value", "adj.P.Val")]
#df_clean <- df_clean[df_clean$adj.P.Val <= 0.999,]
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
p_table_m <- as.ggplot(deg_table_grob)
print(p_table_m)

vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="grey", size=1) +
  theme_bw()

p2_m <- vplot
ggplotly(vplot)

p1_f1 <- p1_f + ggtitle("Fetal N = 15")

p1_m <- p1 + ggtitle("Maternal N = 15")

# with 1_Motor_15plac.R

combined <- (p1_f1 + p2_f + p_table_f) / (p1_m + p2_m + p_table_m) +
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

combined = combined +
  plot_annotation(title = "Motor DEGs category",
                  tag_levels = 'a')
combined
ggsave("results_paper/Fig4_Mot.png", combined, width = 600, height = 400, units = "mm", dpi = 400)
ggsave("results_paper/Fig4_Mot.tiff", combined, width = 600, height = 400, units = "mm", dpi = 400)

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


### Fig2Motor ####

library(patchwork)

combined <- (p1 + p2) / heatmap_ggplot +
  plot_annotation(
    title = "Expression data based on motor group",
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(size = 14),
      plot.tag = element_text(size = 14),
      plot.margin = margin(5, 5, 5, 5)  # Reduce margin to tighten the layout
    )
  ) +
  plot_layout(
    widths = c(1, 1),  # Ensure equal widths for panels a and b
    heights = c(0.5, 1)  # Adjust heights to ensure c is fully visible
    #    guides = "collect"  # Align legends
  )

combined <- combined + plot_layout(ncol = 1, heights = c(0.4, 1))
ggsave("results_paper/Fig6Motor.png", combined, width = 170, height = 200, units = "mm", dpi = 300)

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

multigp = gost(list("up-regulated" = row.names(myModule_up),
                    "down-regulated" = row.names(myModule_down)))

p6 <- gostplot(multigp, interactive = F, capped = T)
p6

p7 <- publish_gostplot(p6)

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

# view results as an interactive table
datatable(myGSEA.df_2,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in Motor B-III',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>% formatRound(columns=c(2:10), digits=2)

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

### load libraries ###
library(ggplot2)
library(enrichplot)
library(ggraph)
library(clusterProfiler)

# Sort for the top upregulated (highest NES) and downregulated (lowest NES) gene sets
top_upregulated <- myGSEA.df_2[order(-myGSEA.df_2$NES), ][1:25, ]  # Top 25 upregulated
top_downregulated <- myGSEA.df_2[order(myGSEA.df_2$NES), ][1:25, ] # Top 25 downregulated
top50_immune <- rbind(top_upregulated, top_downregulated)

ggplot(top50_immune, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(x = "Gene Set", y = "Normalized Enrichment Score (NES)", fill = "p.adjust",
       title = "GSEA Results for curated gene sets") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))



# Sort for the top upregulated (highest NES) and downregulated (lowest NES) gene sets
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


## Schizo Table DEG

schizogenes <- c("CD14", "CENPM", "CRTC2", "DDX28", "ESAM", "GALNT10", "GLTP", "NRGN", "NUDT1", "RNF38",
                     "RP5_1148A213", "RPS17", "SLC25A3", "TMBIM6", "DGCR8", "TAC3", "WNT7", "IL-2R",
                     "LINC02050", "WDR6", "XXcos-LUCA16.1", "RP4-555D20.2", "RP4-555D20.1", "LINC02484",
                     "RP3-407E4.3", "RP5-1148A21.3", "RP4-778K6.3", "MIR4688", "VSIG2", "ESAM-AS1", "TMTC1",
                     "RILPL2", "GPR135", "CTD-2574D22.3", "CIAO2B", "RP11-818O24.3", "MIR138-1", "CUL7",
                     "ATP5MGP8", "ATP5F1B", "P2RX4", "CYP19A1", "ATP13A1", "IRF3", "ATF4", "HSPE1",
                     "STAB1", "BAP1", "FXR1", "GAB1", "TRIM8", "TRIM8-DT", "EIF5", "FURIN", "MAPK3", "SREBP2")

filtered_diffGenes.df <- diffGenes.df[diffGenes.df$geneID %in% schizogenes, ]

datatable(filtered_diffGenes.df,
          extensions = c('KeyTable', 'FixedHeader'),
          caption = 'Filtered Schizo DEGs in samples',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c(10, 25, 50, 100))) %>%
  formatRound(columns = 2:31, digits = 2)

write.table(filtered_diffGenes.df, file = "Motor_Schizo_DEGs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
