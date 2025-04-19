# Load necessary libraries
library(ggplot2)
library(ComplexUpset)
library(tibble)
library(dplyr)
library(tidyr)

### set env ###
getwd()
setwd(here::here("Documents/Dani/Transcriptoma/"))

# Read the data
sex = read.csv("results_paper/tables/Supp_2_Sex_DEG.tsv", sep = "\t")
cog = read.csv("results_paper/tables/Supp_4_Cog_DEG.tsv", sep = "\t")
lang = read.csv("results_paper/tables/Supp_7_Lang_DEG.tsv", sep = "\t")
mot = read.csv("results_paper/tables/Supp_10_Motor_DEG.tsv", sep = "\t")
cog50 = read.csv("results_paper/tables/Supp_3_Cog_top50_25up25down.tsv", sep = "\t")
lang50 = read.csv("results_paper/tables/Supp_6_Lang_top50_25up25down.tsv", sep = "\t")
mot50 = read.csv("results_paper/tables/Supp_9_Motor_top50_up25down25.tsv", sep = "\t")

genes_binary = read.csv("results_paper/tables/intersection_top50.tsv", sep = "\t")

# Add domain column
cog50 = cog50 %>% mutate(Domain = "Cognitive")
lang50 = lang50 %>% mutate(Domain = "Language")
mot50 = mot50 %>% mutate(Domain = "Motor")

# cog = cog %>% mutate(Domain = "Cognitive")
# lang = lang %>% mutate(Domain = "Language")
# mot = mot %>% mutate(Domain = "Motor")

# Combine into long format
genes = bind_rows(cog50, lang50, mot50)

# If you've already loaded MASS or another conflicting package, call select like this:
genes_binary = genes %>%
  dplyr::select(geneID, Domain) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Domain, values_from = present, values_fill = list(present = 0))

readr::write_tsv(genes_binary, "results_paper/tables/intersection_top50.tsv")

# Plot with ComplexUpset
p <- upset(genes_binary, intersect = c("Cognitive", "Language", "Motor"))
p

table1 <- p[[2]][["data"]]

readr::write_tsv(table1, "results_paper/tables/intersection_top50_upset_result.tsv")

genes_fixed <- genes %>%
  mutate(
    logFC = as.numeric(gsub(",", ".", logFC)),
    AveExpr = as.numeric(gsub(",", ".", AveExpr)),
    t = as.numeric(gsub(",", ".", t)),
    P.Value = as.numeric(gsub(",", ".", P.Value)),
    adj.P.Val = as.numeric(gsub(",", ".", adj.P.Val)),
    B = as.numeric(gsub(",", ".", B))
  )

genes_fixed <- genes_fixed %>%
  mutate(Regulation = ifelse(logFC >= 0, "Up", "Down"))

genes_top <- genes_fixed %>%
  group_by(Domain, Regulation) %>%
  arrange(Domain, Regulation, desc(abs(logFC))) %>%
  slice_head(n = 25) %>%
  ungroup()

# Create a combined 'Domain_Regulation' column
genes_top <- genes_top %>%
  mutate(Domain_Regulation = paste(Domain, Regulation, sep = "_"))

genes_binary <- genes_top %>%
  dplyr::select(geneID, Domain_Regulation) %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = Domain_Regulation, values_from = present, values_fill = list(present = 0))

p2 <- upset(genes_binary, intersect = c("Cognitive_Up", "Language_Up", "Motor_Up",
                                  "Cognitive_Down", "Language_Down", "Motor_Down"))
p2

table2 <- p2[[2]][["data"]]

readr::write_tsv(table2, "results_paper/tables/intersection_top50_up_down_result.tsv")

library(patchwork)

combined <- p | p2 +
  plot_annotation(
    title = "Top 25 up and down regulated genes intersection for each domain",
    tag_levels = 'a',
  )


combined
ggsave("results_paper/intersection.png", combined, width = 300, height = 200, units = "mm", dpi = 300)


genes_agg <- genes_binary %>%
  group_by(Motor_Up, Language_Up, Cognitive_Up, Motor_Down, Language_Down, Cognitive_Down) %>%
  summarize(gene_list = paste(geneID, collapse=", "), .groups="drop")

# Join that back to the original data so each row has a gene_list column
genes_binary_annot <- left_join(
  genes_binary,
  genes_agg,
  by = c("Motor_Up","Language_Up","Cognitive_Up","Motor_Down","Language_Down","Cognitive_Down")
)
# Plot with gene_list as the label in the bar
p3 <- upset(
  genes_binary,
  intersect = c("Cognitive", "Language", "Motor"),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,
        color='grey9',
        fill='grey80'
      )
      + ggfittext::geom_bar_text(
        mapping=aes(label=geneID),
        min.size=0,
        position='stack',
        contrast=FALSE,
        vjust=1.1,
      )
    )
  )
)

p3

p4 <- upset(
  genes_binary_annot,
  intersect = c("Motor_Up","Language_Up","Cognitive_Up","Motor_Down","Language_Down","Cognitive_Down"),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,
        color='grey9',
        fill='grey80'
      )
      + ggfittext::geom_bar_text(
        mapping=aes(label=geneID),
        min.size=0,
        position='stack',
        contrast=FALSE,
        vjust=1.1,
      )
    )
  )
)
p4

library(patchwork)

combined <- p3 | p4 +
  plot_annotation(
    title = "Top 50 genes intersection for each domain",
    tag_levels = 'a',
  )


combined
ggsave("results_paper/Fig6_intersection_genes.png", p3, width = 500, height = 500, units = "mm", dpi = 400)
ggsave("results_paper/Fig5_intersection_genes.png", p4, width = 500, height = 500, units = "mm", dpi = 400)





