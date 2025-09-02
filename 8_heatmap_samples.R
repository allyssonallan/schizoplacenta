# Load libraries
library(ComplexHeatmap)
library(circlize)
library(grid)

# Your groups (same as before)
lang_1 <- c("f_S1", "m_S17", "f_S15", "m_S31", "f_S14", "m_S30", "f_S11", "m_S27")
lang_2 <- c("f_S8", "m_S24", "f_S16", "m_S32", "f_S9", "m_S25", "f_S4", "m_S20", "f_S12", "m_S28",
            "f_S10", "m_S26", "f_S5", "m_S21", "f_S6", "m_S22", "f_S3", "m_S19", "f_S2", "m_S18", "f_S7", "m_S23")
cog_1 <- c("f_S1", "m_S17", "f_S12", "m_S28", "f_S10", "m_S26", "f_S5", "m_S21", "f_S6", "m_S22",
           "f_S15", "m_S31", "f_S3", "m_S19", "f_S7", "m_S23", "f_S11", "m_S27")
cog_2 <- c("f_S8", "m_S24", "f_S16", "m_S32", "f_S9", "m_S25", "f_S4", "m_S20", "f_S2", "m_S18", "f_S14", "m_S30")
mot_1 <- c("f_S16", "m_S32", "f_S1", "m_S17", "f_S10", "m_S26", "f_S5", "m_S21", "f_S6", "m_S22",
           "f_S15", "m_S31", "f_S7", "m_S23", "f_S11", "m_S27")
mot_2 <- c("f_S8", "m_S24", "f_S9", "m_S25", "f_S4", "m_S20", "f_S12", "m_S28", "f_S3", "m_S19",
           "f_S2", "m_S18", "f_S14", "m_S30")

sex_1 <- c("f_S9", "f_S5", "f_S6", "f_S3", "f_S11", "m_S25", "m_S21", "m_S22", "m_S19", "m_S27")
sex_2 <- c("f_S8", "f_S16", "f_S4", "f_S1", "f_S12", "f_S10", "f_S15", "f_S2", "f_S7", "f_S14",
           "m_S24", "m_S32", "m_S20", "m_S17", "m_S28", "m_S26", "m_S31", "m_S18", "m_S23", "m_S30")

# Combine and manually sort by sex and ID
extract_numeric_id <- function(x) as.numeric(sub(".*_S", "", x))
samples_f <- unique(c(cog_1, cog_2, lang_1, lang_2, mot_1, mot_2))
samples_f <- samples_f[grepl("^f_", samples_f)]
samples_f <- samples_f[order(extract_numeric_id(samples_f))]

samples_m <- unique(c(cog_1, cog_2, lang_1, lang_2, mot_1, mot_2))
samples_m <- samples_m[grepl("^m_", samples_m)]
samples_m <- samples_m[order(extract_numeric_id(samples_m))]

samples <- c(samples_f, samples_m)

# Status function
get_status <- function(s, g1, g2) {
  if (s %in% g1) return(1)
  if (s %in% g2) return(2)
  return(0)
}

# Build matrix
data_mat <- t(sapply(samples, function(s) {
  c(Sex = ifelse(s %in% sex_1, "Male", ifelse(s %in% sex_2, "Female", "Unknown")),
    Cognitive = get_status(s, cog_1, cog_2),
    Language = get_status(s, lang_1, lang_2),
    Motor = get_status(s, mot_1, mot_2))
}))
rownames(data_mat) <- samples

# Convert domain columns to numeric matrix
domain_matrix <- apply(data_mat[, c("Cognitive", "Language", "Motor")], 2, as.numeric)

# Set color mappings
sex_col_fun <- structure(c("skyblue", "firebrick", "gray90"),
                         names = c("Male", "Female", "Unknown"))

domain_col_fun <- structure(c("gray90", "grey90", "grey30"),
                            names = c("0", "1", "2"))

col = list(Sex = sex_col_fun)
col = sex_col_fun

ht1 <- Heatmap(
  matrix = matrix(data_mat[, "Sex"], ncol = 1, dimnames = list(rownames(data_mat), "Sex")),
  col = sex_col_fun,
  column_names_side = "top",
  name = "Sex",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  width = unit(1.5, "cm")
)

ht2 <- Heatmap(
  matrix = domain_matrix,
  col = domain_col_fun,
  name = "Status",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = samples,
  row_names_gp = gpar(fontsize = 12),
  column_title = "Domain",
  column_names_side = "top",
  width = unit(3, "cm"),
  heatmap_legend_param = list(
    at = c(1, 2),
    labels = c("Healthy", "Altered")
  )
)

ht1 + ht2
