### Communication (or language domain) ###
library(tidyverse)
library(DT)
library(gt)
library(plotly)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(ensembldb)
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
# Validate numeric metadata columns if present
num_cols <- intersect(c("QuantCRN", "QuantLRN", "QuantQRN"), names(targets))
for (col in num_cols) {
  if (!is.numeric(targets[[col]])) stop(paste("Column", col, "must be numeric"))
  if (any(is.na(targets[[col]]))) stop(paste("Column", col, "contains NA values"))
}

targetsmean <- targets %>% # getting the mean of the follow-up visit's scores
  group_by(masctranscRN) %>%
  mutate(
    QuantCRNmean = mean(QuantCRN, na.rm = TRUE),
    QuantLRNmean = mean(QuantLRN, na.rm = TRUE),
    QuantQRNmean = mean(QuantQRN, na.rm = TRUE)
  )

targetsnd <- targetsmean[!duplicated(targetsmean$mascgenRN), ]

sampleLabels <- targetsnd$masctranscRN
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log = TRUE)

#targetsnd <- targets[!duplicated(targets$mascgenRN), ]
pathM <- file.path(targetsnd$masctransM1, "abundance.h5")
pathRN <- file.path(targetsnd$masctranscRN1, "abundance.h5")

targetsnd$language_group <- ifelse(targetsnd$QuantLRNmean >= 90,
                                   "typical",
                                   "altered")

targetsnd$language_group <- factor(targetsnd$language_group)
group <- targetsnd$language_group
group <- factor(group)
groupcont <- targetsnd$QuantLRNmean
groupcont <- factor(groupcont)

path <- file.path(sampleLabels, "abundance.h5")
Tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(
  pathRN,
  type = "kallisto",
  tx2gene = Tx,
  txOut = FALSE,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE
)

myDGEList <- DGEList(Txi_gene$counts)
sampleLabels <- targetsnd$masctranscRN

### proceed with normalization ###

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
  cols = f_S8:f_S11,
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

pca.res <- prcomp(t(log2.cpm.filtered.norm),
                  scale. = F,
                  retx = T)
pc.var <- pca.res$sdev ^ 2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var / sum(pc.var) * 100, 1)
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(
    x = PC1,
    y = PC2,
    label = sampleLabels,
    color = group
  ) +
  geom_point(size = 4) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%", ")")) +
  ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
  coord_fixed() +
  theme_bw()

p1_f <- pca.plot
ggplotly(pca.plot)

log2.cpm.filtered.norm.df <- log2.cpm.filtered.norm.df[!duplicated(names(log2.cpm.filtered.norm.df))]

#continuous data
healthy <- targetsnd %>%
  dplyr::filter(QuantLRNmean >= 90) %>%
  pivot_longer(cols = "masctranscRN",
               names_to = "type",
               values_to = "QuantLRN_value")

disease <- targetsnd %>%
  dplyr::filter(QuantLRNmean < 90) %>%
  pivot_longer(cols = "masctranscRN",
               names_to = "type",
               values_to = "QuantLRN_value")

healthy_concatenated <- healthy[!duplicated(healthy$masctranscRN1),]
healthy_concatenated <- paste(healthy_concatenated$QuantLRN_value, collapse = " + ")
disease_concatenated <- disease[!duplicated(disease$masctranscRN1),]
disease_concatenated <- paste(disease_concatenated$QuantLRN_value, collapse = " + " )

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = eval(parse(text = paste0("(", healthy_concatenated, ")/length(healthy)"))),
                    disease.AVG = eval(parse(text = paste0("(", disease_concatenated, ")/length(disease)"))),
                    LogFC = (disease.AVG - healthy.AVG)) %>%
  mutate_if(is.numeric, round, 2)

datatable(mydata.df,
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))


### Volcano plot
library(gt)
library(DT)
library(plotly)
library(knitr)
library(limma)
library(edgeR)
library(tibble)
library(dplyr)

design <- model.matrix(~ groupcont)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
ebFit <- eBayes(fit)
myTopHits <- topTable(ebFit, adjust = "BH", coef = 2, number = 40000, sort.by = "logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

top_upregulated_genes <- myTopHits.df[order(-myTopHits.df$logFC), ][1:10, ]  # Top 25 upregulated
top_downregulated_genes <- myTopHits.df[order(myTopHits.df$logFC), ][1:10, ] # Top 25 downregulated
top50_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
write.table(top50_genes, file = "Supp_4_Lang_top50_up25down25.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

df_clean <- na.omit(top50_genes)
df_clean <- df_clean[c("geneID", "logFC", "P.Value", "adj.P.Val")]
#df_clean <- df_clean[df_clean$adj.P.Val <= 0.05,]
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
p_table_f <- as.ggplot(deg_table_grob)
print(p_table_f)


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

p2_f <- vplot
ggplotly(vplot)

### DEG table

results <- decideTests(
  ebFit,
  method = "global",
  adjust.method = "BH",
  p.value = 1e-10,
  lfc = 3
)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[, 1] != 0, ]
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

write.table(diffGenes.df, file = "3_1_Lang_Fet_Tissue.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
