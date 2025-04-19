library(WGCNA)
library(tibble)
library(dynamicTreeCut)
library(dplyr)
library(ellipse)
library(Rmisc)
library(cluster)
library(flashClust)
library(car)
library(gProfileR)
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


#### QC process

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



top_upregulated_genes_mat <- mydata.df.mat[order(-mydata.df.mat$LogFC),][1:500, ]
top_upregulated_genes_fet <- mydata.df.mat[order(-mydata.df.fet$LogFC),][1:500, ]

log2.cpm.filtered.norm.df.mat <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% top_upregulated_genes_mat$geneID)
log2.cpm.filtered.norm.df.fet <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% top_upregulated_genes_fet$geneID)

# Combine gene lists and filter the expression data.
combined_genes <- unique(c(top_upregulated_genes_mat$geneID, top_upregulated_genes_fet$geneID))
log2.cpm.filtered.norm.df.combined <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% combined_genes)

# Prepare the expression matrix: set gene IDs as row names and transpose so samples are rows.
datExpr <- log2.cpm.filtered.norm.df.combined %>%
  column_to_rownames("geneID") %>%
  t() %>%
  as.data.frame()

# Define sample information.
sampleInfo <- data.frame(
  SampleID = rownames(datExpr),
  TissueType = c(rep("Maternal", 15), rep("Fetal", 15)),
  Batch = c(rep("Batch1", 15), rep("Batch2", 15))
)
rownames(sampleInfo) <- sampleInfo$SampleID

# For downstream processing, define p (genes x samples) and pd (sample information).
p <- t(datExpr)   # now rows are genes, columns are samples
pd <- sampleInfo

### PART 2. Surrogate Variable Extraction via PCA

# Create an intercept-only model matrix from the phenotype data.
mod0 <- model.matrix(~ 1, data = pd)

# Set the method to use PCA (since we don't have a housekeeping gene list).
preMethod <- "pca"
# Manually set the number of surrogate variables (adjust as needed).
n.sv <- 2

# Extract surrogate variables via PCA on the entire expression matrix.
# (Here we assume that the major unwanted variation can be captured by the top n.sv PCs.)
PCAcomp <- prcomp(t(p), scale = TRUE, center = TRUE)$x[, 1:n.sv, drop = FALSE]
colnames(PCAcomp) <- paste0("pca", 1:n.sv)

# Append the PCA surrogate variables to the phenotype data.
pd <- cbind(pd, PCAcomp)

### PART 3. Empirical Bayes Adjustment

# Define which covariates to retain and which to remove.
# In this case, we retain TissueType and Batch.
retCov <- c("Batch")
# Here we set the unwanted covariates based on the PCA components.
uwv <- if(n.sv > 0) pd[, paste0("pca", 1:n.sv)] else NULL
retCovMat <- pd[, retCov, drop = FALSE]

# Run the empirical Bayes linear model adjustment.
# Note: The empiricalBayesLM function expects data with samples in rows,
# so we transpose p (which is genes x samples) to have samples as rows.
eBLM <- empiricalBayesLM(data = t(p), removedCovariates = uwv, retainedCovariates = retCovMat,
                         weights = NULL, weightType = "empirical", stopOnSmallWeights = TRUE,
                         tol = 1e-4, maxIterations = 1000, scaleMeanToSamples = NULL,
                         robustPriors = TRUE, automaticWeights = "bicov", aw.maxPOutliers = 0.1,
                         verbose = 5, indent = 3)
# The adjusted expression data (with unwanted variation removed) is:
input_matrix <- eBLM$adjustedData
dim(input_matrix)

### PART 4. Outlier Detection

library(ggplot2)
library(ellipse)

# PCA
pr_obj <- prcomp(input_matrix, center = TRUE, scale. = TRUE)

# Create a data frame for plotting
pca_df <- data.frame(
  PC1 = scale(pr_obj$x[,1], scale = TRUE),
  PC2 = scale(pr_obj$x[,2], scale = TRUE),
  TissueType = pd$TissueType
)

# Calculate explained variance
expl_var <- round((pr_obj$sdev^2 / sum(pr_obj$sdev^2))[1:2] * 100, 2)

# Use a full 2x2 correlation matrix
cor_mat <- cor(pca_df[, c("PC1", "PC2")])
el_matrix <- ellipse(cor_mat, level = 0.95, center = c(0, 0), radius = sqrt(qchisq(0.95, 2)), npoints = 10000)
el_df <- as.data.frame(el_matrix)
centre = colMeans(pca_df[, c("PC1", "PC2")])
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TissueType)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("blue", "orange")) +
  xlim(-5, 5) + ylim(-5, 5) +
  labs(
    x = paste0("PC1 (", expl_var[1], "%)"),
    y = paste0("PC2 (", expl_var[2], "%)"),
    title = "TissueType"
  ) +
#  geom_path(data = el_df, aes(x = x, y = y), inherit.aes = FALSE, alpha = 0.3) +
  annotate("text", x = -2.5, y = 2, label = "CI 95%", size = 5) +
  theme_minimal(base_size = 14)

p1 <- pca_plot
print(pca_plot)

library(ggplot2)
library(ggdendro)
library(flashClust)
library(WGCNA)   # for numbers2colors and plotDendroAndColors if needed
library(ellipse) # in case you need ellipse functions elsewhere

# Compute Cook's distance
fit_lm <- lm(pr_obj$x[,1] ~ pr_obj$x[,2])
cooksd <- cooks.distance(fit_lm)
df_cook <- data.frame(
  Subject = seq_along(cooksd),
  Cook = cooksd,
  SubjectName = names(cooksd)
)
mean_cook <- mean(cooksd, na.rm = TRUE)
threshold <- 4 * mean_cook
ci95_val <- CI(cooksd, ci = 0.95)[1]

p_cook <- ggplot(df_cook, aes(x = Subject, y = Cook)) +
  geom_point(shape = 8, size = 3) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_hline(yintercept = ci95_val, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = ifelse(Cook > ci95_val, SubjectName, "")),
            vjust = -0.5, color = "darkred", size = 3) +
  geom_text(data = data.frame(x = 5.2, y = threshold + 0.07, label = "4*mean"),
            aes(x = x, y = y, label = label), color = "red", size = 3) +
  geom_text(data = data.frame(x = 3, y = ci95_val + 0.05, label = "CI 95%"),
            aes(x = x, y = y, label = label), color = "darkred", size = 3) +
  ylim(0, 2.5) +
  labs(title = "Influential Obs by Cook's distance",
       x = "Subjects", y = "Cook's distance") +
  theme_minimal(base_size = 14)
p2 <- p_cook

ggsave("cooks_distance_plot.png", p_cook, width = 8, height = 6)

# Calculate inter-array correlations (IAC)
IAC <- cor(t(input_matrix), use = "p")
IAC_vals <- IAC[upper.tri(IAC)]
df_IAC <- data.frame(IAC = IAC_vals)
mean_IAC <- mean(IAC_vals)
p_hist <- ggplot(df_IAC, aes(x = IAC)) +
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  labs(title = "Histogram of IAC",
       subtitle = paste("Mean =", format(mean_IAC, digits = 3)),
       x = "Correlation", y = "Count") +
  theme_minimal(base_size = 14)

p3 <- p_hist

ggsave("IAC_histogram.png", p_hist, width = 8, height = 6)

# Cluster dendrogram of IAC using hclust
cluster1 <- hclust(as.dist(1 - IAC), method = "average")
dendro_data1 <- dendro_data(as.dendrogram(cluster1), type = "rectangle")
p_dendro <- ggplot() +
  geom_segment(data = segment(dendro_data1),
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro_data1$labels,
            aes(x = x, y = y, label = label),
            hjust = 1, angle = 90, size = 3) +
  scale_y_reverse() +
  labs(title = "Cluster Dendrogram of IAC") +
  theme_minimal(base_size = 14)

p4 <- p_dendro

ggsave("IAC_dendrogram.png", p_dendro, width = 8, height = 6)

# SD distance of average IAC plot
meanIAC <- apply(IAC, 2, mean)
sdCorr <- sd(meanIAC)
numbersd <- (meanIAC - mean(meanIAC)) / sdCorr
df_sd <- data.frame(
  Subject = seq_along(numbersd),
  SDdistance = numbersd,
  SubjectName = names(numbersd)
)
numSD <- 2
numSD3 <- 3

p_sd <- ggplot(df_sd, aes(x = Subject, y = SDdistance)) +
  geom_point() +
  geom_hline(yintercept = -numSD, linetype = "dashed", color = "darkred", size = 1) +
  geom_hline(yintercept = -numSD3, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(label = ifelse(SDdistance < -numSD, SubjectName, "")),
            vjust = 1, color = "darkred", size = 3) +
  ylim(-7, 1) +
  labs(title = "SD distance of average IAC", x = "Subjects", y = "# SD") +
  theme_minimal(base_size = 14)

p5 <- p_sd

ggsave("SD_distance_plot.png", p_sd, width = 8, height = 6)

# Identify outlier samples from Cook's distance and IAC
cooksd_1 <- names(cooksd)[cooksd > 1]
outliers_3sd <- names(numbersd)[numbersd < -numSD3]
removeOut <- union(cooksd_1, outliers_3sd)
good.Samples <- setdiff(rownames(pd), removeOut)
input_matrix <- input_matrix[good.Samples, ]
pd <- pd[good.Samples, ]
if (!identical(rownames(input_matrix), rownames(pd))) stop("Row names do not match")

# Build sample tree using flashClust and create a ggplot dendrogram with batch colors
sampleTree <- flashClust(dist(input_matrix), method = "average")
# Convert TissueType and Batch to colors
dx <- numbers2colors(as.numeric(as.factor(pd$TissueType)), signed = FALSE)
batch_colors <- numbers2colors(as.numeric(as.factor(pd$Batch)), signed = FALSE)
# Create dendrogram data for sampleTree
dendro_data2 <- dendro_data(as.dendrogram(sampleTree), type = "rectangle")
# Merge batch colors with dendrogram labels. Assuming names(batch_colors) match sample labels.
dendro_labels <- dendro_data2$labels
dendro_labels$BatchColor <- batch_colors[match(dendro_labels$label, names(batch_colors))]

p_sampleTree <- ggplot() +
  geom_segment(data = segment(dendro_data2),
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro_labels,
            aes(x = x, y = y, label = label, color = BatchColor),
            hjust = 1, angle = 90, size = 3) +
  scale_color_identity() +
  scale_y_reverse() +
  labs(title = "Subjects Dendrogram with Batch Colors") +
  theme_minimal(base_size = 14)
p6 <- p_sampleTree

ggsave("sample_tree_dendrogram.png", p_sampleTree, width = 8, height = 6)


### PART 5. Network Construction with WGCNA

# Set network parameters.
CPU <- 10
modSize <- 5
networkType <- "signed hybrid"
corType <- "pearson"
corType.sft <- "cor"
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
enableWGCNAThreads(nThreads = 2)
sft <- pickSoftThreshold(input_matrix, dataIsExpr = TRUE, powerVector = powers,
                         RsquaredCut = 0.8, networkType = networkType, corFnc = corType.sft,
                         blockSize = 1, verbose = 5)

#sft <- pickSoftThreshold(input_matrix, dataIsExpr = TRUE, powerVector = powers,
#                         RsquaredCut = 0.8, networkType = networkType, corFnc = corType.sft,
#                         verbose = 5)
beta_value <- sft$powerEstimate
par(mfrow = c(1,2), mar = c(5,5,5,2), lwd = 96/72, ps = 12, las = 1)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "red", lwd = 3)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = "red")

myCor <- function(x, y = NULL, use = "pairwise.complete.obs", method = "pearson", weights.x = NULL, weights.y = NULL, cosine = FALSE) {
  stats::cor(x, y, use = use, method = method)
}

cor <- WGCNA::cor #nstead of it the cor is stats::cor not WGCNA ones

# Construct the network and detect modules.
net <- blockwiseModules(input_matrix, power = beta_value, minModuleSize = modSize,
                        corType = "pearson", networkType = networkType,
                        reassignThreshold = 0, mergeCutHeight = 0.05,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        TOMType = "signed", saveTOMs = FALSE, maxBlockSize = 4000,
                        verbose = 3)

moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(input_matrix, moduleColors, softPower = beta_value)$eigengenes
MEList <- moduleEigengenes(input_matrix, moduleColors, softPower = beta_value)
MEs <- MEList$eigengenes
varExp <- MEList$varExplained

if(is.null(varExp) || length(varExp) == 0) {
  # If varExp is not available, you can either compute it manually or omit it.
  moduleData <- as.data.frame(table(moduleColors))
  colnames(moduleData) <- c("Module", "Size")
} else {
  moduleData <- data.frame(table(moduleColors), varExp = varExp)
}

rownames(moduleData) <- paste0("ME", 0:(nrow(moduleData)-1))
par(mar = c(5,5,5,2), lwd = 96/72, ps = 12, las = 1)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors",
                    main = paste("Genes =", length(moduleColors), "Beta =", beta_value, "\nModules =",
                                 length(unique(moduleColors)) - 1, "Grey =", sum(moduleColors == "grey")),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, cex.main = 1.5)

# Display module eigengene relationships.
plotEigengeneNetworks(MEs, setLabels = NULL, letterSubPlots = FALSE, Letters = NULL,
                      excludeGrey = TRUE, greyLabel = "ME0", plotDendrograms = TRUE,
                      plotHeatmaps = TRUE, setMargins = TRUE, cex.text = 0.8)
title("Eigengenes Dendrogram and Heatmap")

# Compute intramodular connectivity.
IMk <- intramodularConnectivity.fromExpr(datExpr = input_matrix, colors = moduleColors,
                                         corFnc = corType.sft, corOptions = "use = 'p'",
                                         distFnc = "dist", distOptions = "method = 'euclidean'",
                                         networkType = networkType, power = if(networkType == "distance") 1 else 6,
                                         scaleByMax = FALSE, ignoreColors = if(is.numeric(moduleColors)) 0 else "grey",
                                         getWholeNetworkConnectivity = TRUE)
IMksc <- intramodularConnectivity.fromExpr(datExpr = input_matrix, colors = moduleColors,
                                           corFnc = corType.sft, corOptions = "use = 'p'",
                                           distFnc = "dist", distOptions = "method = 'euclidean'",
                                           networkType = networkType, power = if(networkType == "distance") 1 else 6,
                                           scaleByMax = TRUE, ignoreColors = if(is.numeric(moduleColors)) 0 else "grey",
                                           getWholeNetworkConnectivity = FALSE)
# Create a mapping object for gene names.
map <- data.frame(Gene = colnames(datExpr))
imkx <- data.frame(map, moduleColors, IMk, IMksc)

# assign(paste0("MEs.", networkType, ".analysis"), MEs)
# assign(paste0("IMk.", networkType, ".analysis"), imkx)
# assign(paste0("net.", networkType, ".analysis"), net)
# assign(paste0("sft.", networkType, ".analysis"), sft)
# assign(paste0("moduleData.", networkType, ".analysis"), moduleData)
# save(list = c(paste0("MEs.", networkType, ".analysis"),
#               paste0("IMk.", networkType, ".analysis"),
#               paste0("net.", networkType, ".analysis"),
#               paste0("sft.", networkType, ".analysis"),
#               paste0("moduleData.", networkType, ".analysis")),
#      file = paste0("WGCNA.", networkType, ".analysis.RData"))
# dev.off()

### PART 6. Functional Enrichment Analyses

## Gene Ontology Enrichment with g:Profiler

# Get module labels (exclude grey).


#mod <- mod[mod == "grey"]
# Create a list of genes per module.
query.list <- sapply(mod, function(x) imkx$Gene[imkx$moduleColors == x])
background <- imkx$Gene
correction <- "fdr"
GOresults <- list()

mod <- sort(unique(imkx$moduleColors))
grey_genes <- query.list$grey
readr::write_tsv(grey_genes, "results_paper/grey.tsv")

# [1] "ABT1"          "AC007192.4"    "AC011330.5"    "AC023283.1"    "ACE2"          "AP000350.6"    "ASTN1"
#[8] "ATP7A"         "CXCL8"         "LEMD1"         "MUC20P1"       "PGA4"          "PKD1L1"        "RP11-1319K7.1"
#[15] "RP11-195E2.2"  "RP11-673D15.9" "RPL13P12"      "USP32P1"       "UTS2"          "WASH1"


mod <- mod[mod != "grey"]
# Install and load gprofiler2 if not already installed
# install.packages("gprofiler2")
library(gprofiler2)

# Define your parameters (adjust these values as needed)
onto <- c("GO:BP", "GO:MF", "GO:CC", "REAC", "TF", "CORUM", "HP", "HPA")
maxPvalue <- 1        # user threshold (p-value cutoff)
minSize <- 20         # minimum gene set size
maxSize <- 2000       # maximum gene set size
minOverlap <- 0       # minimum overlap between your query and annotated sets
correction <- "fdr"   # correction method

# Assume:
# query.list is a list of character vectors where each element is a set of genes (one per module)
# background is a character vector of all genes (e.g., rownames(imkx))
# networkType is defined (e.g., "signed hybrid") and used in file naming

GOresults <- gost(query = query.list,
                  organism = "hsapiens",
                  ordered_query = FALSE,
                  multi_query = TRUE,
                  significant = TRUE,
                  exclude_iea = FALSE,
                  correction_method = correction,
                  domain_scope = "annotated",
                  user_threshold = maxPvalue,
                  sources = onto,
                  custom_bg = background)
p <- gostplot(GOresults, interactive = T, capped = F)
p

p <- p + scale_y_continuous(limits = c(0, 4))
p

combined_results <- do.call(rbind, lapply(GOresults, function(x) x$result))

top10_terms <- p$data %>%
  arrange(p_value) %>%
  head(10) %>%
  pull(term_id)

pp <- publish_gostplot(p, width = NA, height = NA, filename = NULL, highlight_terms = top10_terms)
pp

# If the column negative_log10_of_adjusted_p_value is missing, create it:
if(!"negative_log10_of_adjusted_p_value" %in% colnames(combined_results)){
  combined_results$negative_log10_of_adjusted_p_value <- -log10(combined_results$p_value)
}

# It also helps to have a query column. If missing, add a dummy one:
if(!"query" %in% colnames(combined_results)){
  combined_results$query <- "combined"
}

# Create a meta element; include at least the query and organism.
combined_meta <- list(query = unique(combined_results$query), organism = "hsapiens")

# Combine into one object.
GOcombined <- list(meta = combined_meta, result = combined_results)

grey_genes <- imkx$Gene[imkx$moduleColors == "grey"]
cat(paste(grey_genes, collapse = ", "))


# # Optionally, assign and save the GO results:
# assign(paste0("GOprofileR.", correction, "0.05", ".min", minSize, ".max", maxSize,
#               ".", networkType, ".analysis"), GOresults)
# save(list = paste0("GOprofileR.", correction, "0.05", ".min", minSize, ".max", maxSize,
#                    ".", networkType, ".analysis"),
#      file = paste0("GOprofileR.", correction, "0.05", ".min", minSize, ".max", maxSize,
#                    ".", networkType, ".analysis.RData"))

## Hypergeometric Enrichment for PGC Schizophrenia Risk Genes

# Load your PGC gene lists.

# Define your schizophrenia gene set.
schizogenes <- c("CD14", "CENPM", "CRTC2", "DDX28", "ESAM", "GALNT10", "GLTP", "NRGN", "NUDT1", "RNF38",
                 "RP5_1148A213", "RPS17", "SLC25A3", "TMBIM6", "DGCR8", "TAC3", "WNT7", "IL-2R",
                 "LINC02050", "WDR6", "XXcos-LUCA16.1", "RP4-555D20.2", "RP4-555D20.1", "LINC02484",
                 "RP3-407E4.3", "RP5-1148A21.3", "RP4-778K6.3", "MIR4688", "VSIG2", "ESAM-AS1", "TMTC1",
                 "RILPL2", "GPR135", "CTD-2574D22.3", "CIAO2B", "RP11-818O24.3", "MIR138-1", "CUL7",
                 "ATP5MGP8", "ATP5F1B", "P2RX4", "CYP19A1", "ATP13A1", "IRF3", "ATF4", "HSPE1",
                 "STAB1", "BAP1", "FXR1", "GAB1", "TRIM8", "TRIM8-DT", "EIF5", "FURIN", "MAPK3", "SREBP2")

# First, create a list of gene sets for each module.
module_list <- lapply(unique(imkx$moduleColors), function(mod) {
  imkx$Gene[imkx$moduleColors == mod]
})
names(module_list) <- unique(imkx$moduleColors)

# Now, for each module, compute the intersection with your schizophrenia genes.
module_matches <- lapply(module_list, function(genes) {
  intersect(genes, genes_binary$geneID) # from 5_intersection
})

# To see how many matches occur in each module, you can print the counts:
match_counts <- sapply(module_matches, length)
print(match_counts)

# To print the matching genes for each module, separated by commas:
for(mod in names(module_matches)) {
  cat(mod, ":", paste(module_matches[[mod]], collapse = ", "), "\n")
}

# Set up the target gene set as a named list.
target <- "Schizophrenia"
ll.target <- list(schizogenes)
names(ll.target) <- target

# Prepare your enrichment data.
imk.enrich <- imkx
# If the gene names are stored as rownames, use them (or change this if you already have a Gene column)
imk.enrich$Gene <- rownames(imk.enrich)
all.genes <- unique(imk.enrich$Gene)

# Get the modules (excluding the "grey" module).
mod <- sort(unique(imk.enrich$moduleColors))
mod <- mod[mod == "grey"]
# For each module, collect the unique gene symbols.
ll.modules <- sapply(mod, function(i) unique(imk.enrich$Gene[imk.enrich$moduleColors == i]), simplify = FALSE)
ll.modules.size <- sapply(ll.modules, length)

# Determine the overlap (hits) between each target and each module.
hit.genes <- lapply(ll.target, function(g) lapply(ll.modules, function(h) intersect(g, h)))

# Build the enrichment matrix. For each target, create a data frame with the number of hits per module.
enrich.matrix <- lapply(ll.target, function(g) {
  data.frame(hit = sapply(ll.modules, function(h) length(intersect(g, h))),
             stringsAsFactors = FALSE)
})
# For each target, count how many of its genes are present in the background.
hit.pop <- sapply(ll.target, function(g) length(intersect(g, all.genes)))
enrich.matrix <- lapply(1:length(enrich.matrix), function(i) {
  cbind(enrich.matrix[[i]],
        hit.pop = rep(hit.pop[i], length(ll.modules)))
})
# Calculate the number of genes in the background that are not in the target.
fail.pop <- sapply(ll.target, function(g) length(all.genes) - length(intersect(g, all.genes)))
enrich.matrix <- lapply(1:length(enrich.matrix), function(i) {
  cbind(enrich.matrix[[i]],
        fail.pop = rep(fail.pop[i], length(ll.modules)))
})
# Add the module (drawn) sizes.
enrich.matrix <- lapply(1:length(enrich.matrix), function(i) {
  cbind(enrich.matrix[[i]],
        drawn = ll.modules.size)
})

# Compute the hypergeometric p-values for each module.
enrich.matrix.pval <- sapply(1:length(enrich.matrix), function(i)
  sapply(1:nrow(enrich.matrix[[i]]), function(r)
    phyper(enrich.matrix[[i]][r, "hit"] - 1,
           enrich.matrix[[i]][r, "hit.pop"],
           enrich.matrix[[i]][r, "fail.pop"],
           enrich.matrix[[i]][r, "drawn"],
           lower.tail = FALSE)))
enrich.matrix.pval <- data.frame(enrich.matrix.pval)
rownames(enrich.matrix.pval) <- names(ll.modules)
colnames(enrich.matrix.pval) <- names(ll.target)

# Append the computed p-values and adjust for multiple testing.
enrich.matrix <- lapply(1:length(enrich.matrix), function(i) {
  df <- enrich.matrix[[i]]
  df$hyper.geo.pval <- enrich.matrix.pval[, i]
  df$hyper.geo.fdr <- p.adjust(df$hyper.geo.pval, method = "fdr")
  df$hyper.geo.bonf <- p.adjust(df$hyper.geo.pval, method = "bonferroni")
  df$hyper.geo.log10.pval <- -log10(df$hyper.geo.pval)
  df
})
names(enrich.matrix) <- names(ll.target)


# assign(paste("enrich.matrices", target, "analysis", sep = "."), enrich.matrix)
# assign(paste("enrich.pval", target, "analysis", sep = "."), enrich.matrix.pval)
# assign(paste("modules", "analysis", sep = "."), ll.modules)
# assign(paste("hit", "analysis", sep = "."), hit.genes)
# save(list = c(paste("enrich.matrices", target, "analysis", sep = "."),
#               paste("enrich.pval", target, "analysis", sep = "."),
#               paste("modules", "analysis", sep = "."),
#               paste("hit", "analysis", sep = ".")),
#      file = paste0(target, "enrich.", networkType, ".analysis.RData"))

