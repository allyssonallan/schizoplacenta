library(tidyverse)
library(DT)
library(gt)
library(plotly)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(ggplot2)
library(reshape2)
library(gplots)
library(rhdf5)

getwd()
setwd(here::here("Documents/Dani/Transcriptoma/"))

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
  group_by(masctransM) %>%
  mutate(QuantCRNmean = mean(QuantCRN, na.rm = TRUE),
         QuantLRNmean = mean(QuantLRN, na.rm = TRUE),
         QuantQRNmean = mean(QuantQRN, na.rm = TRUE))
targetsnd <- targetsmean[!duplicated(targetsmean$mascgenRN), ]
sampleLabels <- c(targetsnd$masctransM1, targetsnd$masctranscRN1)

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
sampleLabels <- c(targetsnd$masctransM, targetsnd$masctranscRN)



cpm <- cpm(myDGEList)
keepers <- rowSums(cpm > 1) >= 5
myDGEList.filtered <- myDGEList[keepers,]
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log = TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = m_S24:f_S11,
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


# Load required libraries
library(ggstatsplot)
library(patchwork) # For combining plots
library(gplots)
library(data.table)

# Combine maternal and fetal samples for pairwise correlation
#### Heatmap to help visualize clusters ####
mycolors <- colorRampPalette(c("yellow", "red"))

maternal_samples <- sampleLabels[1:15]
fetal_samples <- sampleLabels[16:30]
combined_samples <- c(maternal_samples, fetal_samples)
cor_matrix_maternal <- cor(log2.cpm.filtered.norm.df[, maternal_samples])
cor_matrix_fetal <- cor(log2.cpm.filtered.norm.df[, fetal_samples])
cor_matrix_combined <- cor(log2.cpm.filtered.norm.df[, combined_samples])
cor_matrix_maternal_dt <- as.data.table(cor_matrix_maternal, keep.rownames = "Var1")
cor_matrix_fetal_dt <- as.data.table(cor_matrix_fetal, keep.rownames = "Var1")
cor_matrix_combined_dt <- as.data.table(cor_matrix_combined, keep.rownames = "Var1")
maternal_long <- melt(cor_matrix_maternal_dt, id.vars = "Var1", variable.name = "Var2", value.name = "Correlation")
fetal_long <- melt(cor_matrix_fetal_dt, id.vars = "Var1", variable.name = "Var2", value.name = "Correlation")
combined_long <- melt(cor_matrix_combined_dt, id.vars = "Var1", variable.name = "Var2", value.name = "Correlation")


# Plot using ggcorrmat
combined_heatmap <- ggcorrmat(
  data = as.data.frame(cor_matrix_combined),  # Maternal correlation matrix as data frame
  cor.vars = colnames(cor_matrix_combined),  # Use all columns for correlation
  cor.vars.names = combined_samples,         # Friendly names for samples                            # Non-parametric correlations
  title = "Cognitive: Maternal and fetal",        # Title for the plot
  subtitle = "Full correlation matrix",      # Subtitle for context
  colors = c("yellow", "orange", "red"),               # Yellow to red gradient
  ggcorrplot.args = list(
    lab_col = "black",                       # Correlation coefficient text color
    lab_size = 4,                            # Correlation coefficient text size
    tl.srt = 45                              # Rotate axis labels
  )
)

maternal_heatmap <- ggcorrmat(
  data = as.data.frame(cor_matrix_maternal),  # Maternal correlation matrix as data frame
  cor.vars = colnames(cor_matrix_maternal),  # Use all columns for correlation
  cor.vars.names = maternal_samples,         # Friendly names for samples
  title = "Cognitive: Maternal side",        # Title for the plot
  subtitle = "Full correlation matrix",      # Subtitle for context
  colors = c("yellow", "orange", "red"),     # Yellow to red gradient
  ggcorrplot.args = list(
    lab_col = "black",                       # Correlation coefficient text color
    lab_size = 4,                            # Correlation coefficient text size
    tl.srt = 45                              # Rotate axis labels
  )
)

fetal_heatmap <- ggcorrmat(
  data = as.data.frame(cor_matrix_fetal),    # Fetal correlation matrix as data frame
  cor.vars = colnames(cor_matrix_fetal),     # Use all columns for correlation
  cor.vars.names = fetal_samples,            # Friendly names for samples
  title = "Cognitive: Fetal side",           # Title for the plot
  subtitle = "Full correlation matrix",      # Subtitle for context
  colors = c("yellow", "orange", "red"),     # Yellow to red gradient
  ggcorrplot.args = list(
    lab_col = "black",                       # Correlation coefficient text color
    lab_size = 4,                            # Correlation coefficient text size
    tl.srt = 45                              # Rotate axis labels
  )
)

# Combine all plots in a grid layout
final_plot <- (combined_heatmap / (maternal_heatmap + fetal_heatmap)) +
  plot_layout(heights = c(1, 1), guides = "collect") # Adjust relative heights of rows

# Display the final plot
print(final_plot)

combined_heatmap <- ggplot(combined_long, aes(x = Column, y = Row, fill = Value)) +
  geom_tile() +
  heatmap_colors +
  labs(title = "Fetal and maternal") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Plot the maternal heatmap
maternal_heatmap <- ggplot(maternal_long, aes(x = Column, y = Row, fill = Value)) +
  geom_tile() +
  heatmap_colors +
  labs(title = "Maternal") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Plot the fetal heatmap
fetal_heatmap <- ggplot(fetal_long, aes(x = Column, y = Row, fill = Value)) +
  geom_tile() +
  heatmap_colors +
  labs(title = "Fetal") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Combine all plots in a grid layout
final_plot <- (combined_heatmap / (maternal_heatmap + fetal_heatmap)) +
  plot_layout(heights = c(1, 1), guides = "collect") # Adjust relative heights of rows

print(final_plot)


##########

library(gplots)

# Define the color palette
mycolors <- colorRampPalette(c("yellow", "red"))
colnames(log2.cpm.filtered.norm) <- sampleLabels
cor_matrix <- cor(log2.cpm.filtered.norm)

couples <- data.frame(f = c("f_S8","f_S16","f_S9","f_S4","f_S1","f_S12","f_S10",
                             "f_S5","f_S6","f_S15","f_S3","f_S2","f_S7","f_S14",
                             "f_S11"),
                      m = c("m_S24","m_S32","m_S25","m_S20","m_S17","m_S28",
                             "m_S26","m_S21","m_S22","m_S31","m_S19","m_S18",
                             "m_S23","m_S30","m_S27"))
mycolors <- rainbow(nrow(couples))
col_vector <- setNames(rep(mycolors, each=2), c(rbind(couples$f, couples$m)))

# Create the heatmap
heatmap.2(cor_matrix, trace = "none", col = mycolors,
          dendrogram="row",
          Colv="NA",
          RowSideColors = col_vector[rownames(cor_matrix)],
          main = "Gene expression TMM Normalized",
          labRow = colnames(log2.cpm.filtered.norm),
          labCol = colnames(log2.cpm.filtered.norm),
          margins = c(8, 8))


targetsmean <- targets %>%
  group_by(masctranscRN) %>%
  mutate(
    QuantCRNmean = mean(QuantCRN, na.rm = TRUE),
    QuantLRNmean = mean(QuantLRN, na.rm = TRUE),
    QuantQRNmean = mean(QuantQRN, na.rm = TRUE)
  )
targetsnd <- targetsmean[!duplicated(targetsmean$mascgenRN), ]
#continuous
targetsnd$cognitive_group <- ifelse(targetsnd$QuantCRN >= 90,
                                    "normal_cognitive",
                                    "low_cognitive")
targetsnd$cognitive_group <- factor(targetsnd$cognitive_group)
group <- c(targetsnd$cognitive_group, targetsnd$cognitive_group)
group <- factor(group)
groupcont <- c(targetsnd$QuantCRNmean, targetsnd$QuantCRNmean)
groupcont <- factor(groupcont)
sampleLabels <- c(targetsnd$masctransM, targetsnd$masctranscRN)


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

p1 <- pca.plot
ggplotly(pca.plot)


### Filtered and normalized counts per million (CPM)

log2.cpm.filtered.norm.df <- log2.cpm.filtered.norm.df[!duplicated(names(log2.cpm.filtered.norm.df))]

#continuous data
healthy <- targetsnd %>%
  dplyr::filter(QuantCRN >= 90) %>%
  tidyr::pivot_longer(
    cols = c("masctranscRN", "masctransM"),
    names_to = "type",
    values_to = "QuantCRN_value"
  )

altered <- targetsnd %>%
  dplyr::filter(QuantCRN < 90) %>%
  tidyr::pivot_longer(
    cols = c("masctranscRN", "masctransM"),
    names_to = "type",
    values_to = "QuantCRN_value"
  )


healthy_concatenated <- paste(healthy$QuantCRN_value, collapse = " + ")
altered_concatenated <- paste(altered$QuantCRN_value, collapse = " + ")

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = eval(parse(text = paste0("(", healthy_concatenated, ")/length(healthy)"))),
                    altered.AVG = eval(parse(text = paste0("(", altered_concatenated, ")/length(altered)"))),
                    LogFC = (healthy.AVG - altered.AVG)) %>%
  mutate_if(is.numeric, round, 2)

extract_samples <- function(concatenated_string, pattern) {
  stringr::str_extract_all(concatenated_string, pattern) %>%
    unlist()
}

fetal_pattern <- "f_[^ +]+"
maternal_pattern <- "m_[^ +]+"
healthy_fetal <- extract_samples(healthy_concatenated, fetal_pattern)
healthy_maternal <- extract_samples(healthy_concatenated, maternal_pattern)
altered_fetal <- extract_samples(altered_concatenated, fetal_pattern)
altered_maternal <- extract_samples(altered_concatenated, maternal_pattern)
cat("Healthy Fetal Samples: ", paste(healthy_fetal, collapse = ", "), "\n")
cat("Healthy Maternal Samples: ", paste(healthy_maternal, collapse = ", "), "\n")
cat("Altered Fetal Samples: ", paste(altered_fetal, collapse = ", "), "\n")
cat("Altered Maternal Samples: ", paste(altered_maternal, collapse = ", "), "\n")

mydata.df.mat <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = eval(parse(text = paste0("(", healthy_maternal, ")/length(healthy_maternal)"))),
                    altered.AVG = eval(parse(text = paste0("(", altered_maternal, ")/length(altered_maternal)"))),
                    LogFC = (healthy.AVG - altered.AVG)) %>%
  mutate_if(is.numeric, round, 2)

mydata.df.fet <- mutate(log2.cpm.filtered.norm.df,
                        healthy.AVG = eval(parse(text = paste0("(", healthy_fetal, ")/length(healthy_fetal)"))),
                        altered.AVG = eval(parse(text = paste0("(", altered_fetal, ")/length(altered_fetal)"))),
                        LogFC = (healthy.AVG - altered.AVG)) %>%
  mutate_if(is.numeric, round, 2)

library(dplyr)
library(ggstatsplot)

mydata.df.mat <- mydata.df.mat %>%
  mutate(Group = "Maternal") %>%
  dplyr::select(geneID, LogFC, Group)

mydata.df.fet <- mydata.df.fet %>%
  mutate(Group = "Fetal") %>%
  dplyr::select(geneID, LogFC, Group)

combined_data <- bind_rows(mydata.df.mat, mydata.df.fet)

comparison_data <- combined_data %>%
  pivot_wider(names_from = Group, values_from = LogFC) %>%
  drop_na()  # Remove genes without values in both groups

comparison_data <- comparison_data %>%
  mutate(Difference = Maternal - Fetal)

library(ggstatsplot)

comparison_data <- comparison_data %>%
  dplyr::mutate(label_condition = ifelse(Maternal > 0.9 & Fetal > 0.9, geneID, NA))

comparison_data_down <- comparison_data %>%
  dplyr::mutate(label_condition = ifelse(Maternal < -0.8 & Fetal < -0.8, geneID, NA))

ggstatsplot::ggscatterstats(
  data = comparison_data_down,
  x = Maternal,
  y = Fetal,
  label.var = geneID,
  label.expression = c(Maternal < -0.8 | Maternal > 0.66),
  title = "Cognitive Domain",
  xlab = "LogFC (Maternal)",
  ylab = "LogFC (Fetal)",
  ggtheme = ggplot2::theme_minimal()
)

library(presto)

results <- wilcoxauc(
  data = myTopHits.df,
  group_by = "group",
  features = "logFC"
)

library(ggstatsplot)

# Assuming 'log2.cpm.filtered.norm' is your data frame and 'side' is the column indicating maternal or fetal
# Add a 'side' column to your data if not present
log2.cpm.filtered.norm$side <- rep(c("Maternal", "Fetal"), each = 15)

# Perform the comparison using ggbetweenstats
ggbetweenstats(
  data = log2.cpm.filtered.norm,
  x = side,  # Grouping variable
  y = logFC,  # The variable to compare, replace with the appropriate column
  title = "Comparison of Gene Expression between Maternal and Fetal Sides"
)


top_upregulated_genes_mat <- mydata.df.mat[order(-mydata.df.mat$LogFC),][1:500, ]
top_upregulated_genes_fet <- mydata.df.mat[order(-mydata.df.fet$LogFC),][1:500, ]

log2.cpm.filtered.norm.df.mat <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% top_upregulated_genes_mat$geneID)
log2.cpm.filtered.norm.df.fet <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% top_upregulated_genes_fet$geneID)



library(WGCNA)
library(tibble)
library(dynamicTreeCut)

combined_genes <- unique(c(top_upregulated_genes_mat$geneID, top_upregulated_genes_fet$geneID))
log2.cpm.filtered.norm.df.combined <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% combined_genes)


# Set gene IDs as row names and transpose the data
datExpr <- log2.cpm.filtered.norm.df.combined %>%
  tibble::column_to_rownames("geneID") %>%  # Set geneID as row names
  t() %>%                                   # Transpose the data
  as.data.frame()

# Assuming you have variables 'number_of_maternal_samples' and 'number_of_fetal_samples'
sampleInfo <- data.frame(
  SampleID = rownames(datExpr),
  TissueType = c(rep("Maternal", 15), rep("Fetal", 15))
)
rownames(sampleInfo) <- sampleInfo$SampleID

# Define covariates (e.g., batch, tissue type)
# Replace `covariates` with the actual variables you want to adjust for
covariates <- data.frame(
  Batch = c(rep("Batch1", 15), rep("Batch2", 15)),  # Example batch info
  TissueType = c(rep("Maternal", 15), rep("Fetal", 15))
)
rownames(covariates) <- rownames(datExpr)

# Fit a linear model for each gene and extract residuals
adjusted_expr <- apply(datExpr, 2, function(x) {
  model <- lm(x ~ ., data = covariates)
  resid(model)  # Extract residuals
})

# Convert residual matrix back to data frame
datExpr <- as.data.frame(adjusted_expr)
rownames(datExpr) <- rownames(covariates)

# Proceed with WGCNA analysis
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Check for good samples and genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Hierarchical clustering of samples
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub = "", xlab = "")

# Choose a set of soft-thresholding powers
powers <- c(1:20)

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results to choose the power
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale-Free Topology Model Fit,signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.90, col="red")  # Reference line indicating R^2 cut-off

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")

# Select the soft-thresholding power
softPower <- sft$powerEstimate

# Construct the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

# Calculate the TOM and its dissimilarity
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Cluster the genes using the dissimilarity TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot the gene dendrogram
plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based Dissimilarity", labels = FALSE, hang = 0.04)

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)

# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the eigengene dendrogram
plot(METree, main = "Clustering of Module Eigengenes", xlab = "", sub = "")

# Ensure TissueType is numeric
sampleInfo$TissueTypeNumeric <- as.numeric(as.factor(sampleInfo$TissueType))

# Number of samples
nSamples <- nrow(datExpr)

# Correlate module eigengenes with TissueType
moduleTraitCor <- cor(MEs, sampleInfo$TissueTypeNumeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
genesInModule <- colnames(datExpr)
moduleGenes <- dynamicColors
genesInModule <- colnames(datExpr)[moduleGenes]
textMatrix <- paste(signif(moduleTraitCor*10e13, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Plot the module-trait relationships
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor*10e13,
               xLabels = "TissueType",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.025,0.025),
               main = "Module-TissueType Relationships (10e13)")

dev.off()



sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub = "", xlab = "")

# For example, extract genes from a module of inter
module <- "grey"  # Replace with your module color of interest
moduleGenes <- dynamicColors == module
genesInModule <- colnames(datExpr)[moduleGenes]
cat(genesInModule, sep = "\n")
cat(genesInModule, sep = ",")
# AC012501.3,AC068533.7,AC093724.2,AL109923.1,BCL2L2-PABPN1,C16orf45,CSPG4P11,
# DHRS9,EFCAB8,FAAH2,GOLGA8O,IGHV1-46,IGHV3-74,IGSF22,NPIPB6,PMS2,PPIAP22,PPP1R17,
# RP11-497H16.7,RP11-498E11.2,RP11-529K1.3,RP11-566K11.2,RP11-848G14.2,
# RP3-461P17.10,RPL7P57,SEC13P1,SLC12A1,SP5,SPSB4,SYT10,TNFRSF4,TPSAB1,
# USP32P1,ZNF275,ZNF526
length(genesInModule)

