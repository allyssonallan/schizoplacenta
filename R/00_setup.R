# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 00_setup.R: packages, paths, metadata, counts, gene sets
# Source via run_all.R from the project root (folder that contains 0_TableCode.tsv).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(edgeR); library(limma)
  library(ggplot2); library(ggrepel); library(patchwork); library(cowplot)
  library(gridExtra); library(ComplexHeatmap); library(circlize); library(RColorBrewer)
  library(stringr); library(scales); library(grid)
})
set.seed(20260928)
ROOT <- Sys.getenv("SCHIZO_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "results_revision"); TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures"); RDS <- file.path(OUT, "rds")
for (d in c(TAB, FIG, RDS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
wt <- function(x, f) fwrite(as.data.table(x), file.path(TAB, f), sep = "\t", na = "NA")
DOM <- c(cog = "Cognitive", lang = "Language", mot = "Motor")
# ------------------------------------------------------------------ metadata
visits <- fread(file.path(ROOT, "0_TableCode.tsv"))
first  <- visits[!duplicated(mascgenRN)]
means  <- visits[, .(cog = mean(QuantCRN, na.rm = TRUE), lang = mean(QuantLRN, na.rm = TRUE), mot = mean(QuantQRN, na.rm = TRUE)), by = mascgenRN]
close24 <- rbindlist(lapply(c(cog = "QuantCRN", lang = "QuantLRN", mot = "QuantQRN"), function(sc)
  visits[!is.na(get(sc))][order(abs(ageRN - 24), -ageRN), .SD[1], by = mascgenRN][, .(mascgenRN, v = get(sc))]), idcol = "dom")
close24 <- dcast(close24, mascgenRN ~ dom, value.var = "v"); setnames(close24, c("cog","lang","mot"), c("cog24","lang24","mot24"))
part <- merge(first[, .(participant = mascgenRN, fdir = masctranscRN1, f = masctranscRN, mdir = masctransM1, m = masctransM, sex, ageM, RINrn, RINm)], means, by.x = "participant", by.y = "mascgenRN", sort = FALSE)
part <- merge(part, close24, by.x = "participant", by.y = "mascgenRN", sort = FALSE)
meta <- rbind(part[, .(sample = m, dir = mdir, participant, face = "maternal", sex, ageM, rin = RINm, cog, lang, mot, cog24, lang24, mot24)],
              part[, .(sample = f, dir = fdir, participant, face = "fetal",    sex, ageM, rin = RINrn, cog, lang, mot, cog24, lang24, mot24)])
meta[, face := factor(face, levels = c("maternal", "fetal"))]
meta[, sex := factor(sex, levels = c("F", "M"))]
meta[, participant := factor(participant)]
meta[, age_z := as.numeric(scale(ageM))]
meta[, rin_z := as.numeric(scale(rin))]
for (d in c("cog", "lang", "mot")) {
  meta[, (paste0(d, "_grp")) := factor(ifelse(get(d) < 90, "DD", "TD"), levels = c("TD", "DD"))]
  meta[, (paste0(d, "24_grp")) := factor(ifelse(get(paste0(d, "24")) < 90, "DD", "TD"), levels = c("TD", "DD"))]
}
wt(meta, "S_sample_metadata.tsv")
wt(part[, .(participant, sex, maternal_age = ageM, cognitive_mean = cog, language_mean = lang, motor_mean = mot,
            cognitive_group = ifelse(cog < 90, "DD", "TD"), language_group = ifelse(lang < 90, "DD", "TD"), motor_group = ifelse(mot < 90, "DD", "TD"))],
   "S_participant_outcomes.tsv")

# ------------------------------------------------------------------ counts
cache <- file.path(ROOT, "results_revision", "rds", "gene_counts_lengthScaledTPM.rds")
if (file.exists(cache)) {
  counts <- readRDS(cache)
} else {
  suppressPackageStartupMessages({ library(tximport); library(ensembldb); library(EnsDb.Hsapiens.v86) })
  files <- file.path(ROOT, meta$dir, "abundance.h5"); names(files) <- meta$sample
  tx <- as.data.frame(transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name")))[, c("tx_id", "gene_name")]
  counts <- tximport(files, type = "kallisto", tx2gene = tx, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)$counts
  saveRDS(counts, cache)
}
counts <- counts[, meta$sample]
dge0 <- DGEList(counts); keep <- rowSums(cpm(dge0) > 1) >= 5
dge  <- calcNormFactors(dge0[keep, , keep.lib.sizes = FALSE], method = "TMM")
logcpm <- cpm(dge, log = TRUE)
cat("Genes retained:", nrow(dge), " samples:", ncol(dge), "\n")

# ------------------------------------------------------------------ gene sets
read_gmt <- function(f) { l <- strsplit(readLines(f), "\t"); setNames(lapply(l, function(x) x[-(1:2)]), sapply(l, `[`, 1)) }
gs <- list()
gs$C2CP  <- read_gmt(file.path(ROOT, "c2.cp.v7.1.symbols.gmt"))
gs$C2ALL <- read_gmt(file.path(ROOT, "c2.all.v7.1.symbols.gmt"))
for (cc in c("C7", "C5")) {  # optional exports from msigdbr (see response letter)
  f <- file.path(ROOT, paste0("msigdb_", cc, ".tsv"))
  if (file.exists(f)) {
    x <- fread(f)
    if (cc == "C5") { gs$HP <- split(x[grepl("HPO", gs_subcat)]$gene_symbol, x[grepl("HPO", gs_subcat)]$gs_name); c5 <- x } else gs[[cc]] <- split(x$gene_symbol, x$gs_name)
  }
}
# over-representation sources (g:Profiler-like, Figs. S4/S7/S10): GO BP/CC/MF and HPO from the C5 export, KEGG/Reactome/WikiPathways/PID/BioCarta from C2:CP
ora_sets <- list()
if (exists("c5")) for (sc in c("GO:BP", "GO:CC", "GO:MF", "HPO")) { y <- c5[gs_subcat == sc]; ora_sets[[sub("HPO", "HP", sc)]] <- split(y$gene_symbol, y$gs_name) }
for (src in c("KEGG", "REACTOME", "PID", "BIOCARTA")) ora_sets[[sub("REACTOME", "REAC", src)]] <- gs$C2CP[grepl(paste0("^", src, "_"), names(gs$C2CP))]
extf <- file.path("data", "external_gene_sets.tsv")
ext <- fread(extf)
extsets <- split(ext$gene, ext$gene_set)
# lists used in the original manuscript (Supplementary Tables 12-13)
extsets[["SCZ GWAS placenta-implicated (Ursini 2023 S1; 56)"]] <- c("CD14","CENPM","CRTC2","DDX28","ESAM","GALNT10","GLTP","NRGN","NUDT1","RNF38","RP5_1148A213","RPS17","SLC25A3","TMBIM6","DGCR8","TAC3","WNT7","IL-2R","LINC02050","WDR6","XXcos-LUCA16.1","RP4-555D20.2","RP4-555D20.1","LINC02484","RP3-407E4.3","RP5-1148A21.3","RP4-778K6.3","MIR4688","VSIG2","ESAM-AS1","TMTC1","RILPL2","GPR135","CTD-2574D22.3","CIAO2B","RP11-818O24.3","MIR138-1","CUL7","ATP5MGP8","ATP5F1B","P2RX4","CYP19A1","ATP13A1","IRF3","ATF4","HSPE1","STAB1","BAP1","FXR1","GAB1","TRIM8","TRIM8-DT","EIF5","FURIN","MAPK3","SREBP2")
extsets[["PlacPRS1 schizophrenia genes (Ursini 2021; 66)"]] <- c("ACO2","ALKBH5","ARHGAP1","ARL6IP4","ARPC3","ATF4","ATP2A2","ATPAF2","BZW1","CCDC25","CLK1","CLUH","CNOT1","CNOT7","CREB3L2","CYB5R3","DMTF1","DPP4","EFHD1","ENO1","ERRFI1","ESAM","FLII","FXR1","H2BC12","HSPA8","HSPD1","HSPE1","IFIH1","INPP5D","LARGE2","MORF4L1","NAGA","NCK1","NDRG4","NDUFA6","PAFAH1B1","PAPPA2","PARK7","PBK","PCCB","PHF5A","PKNOX2","PMM1","PSMA4","RFTN2","RORA","RPS19BP1","SBNO1","SCARA3","SEPTIN10","SERPINE2","SF3B1","SIL1","SMDT1","SNX3","SREBF1","SREBF2","TCTN1","THEM6","THOC7","TXNL1","VPS26B","WDR61","XRCC6","ZNF682")
ext_manifest <- data.table(gene_set = names(extsets), genes_listed = lengths(extsets), genes_expressed = sapply(extsets, function(g) sum(g %in% rownames(dge))))
wt(ext_manifest, "S_external_gene_set_manifest.tsv")
extsets_ok <- extsets[ext_manifest$genes_expressed >= 5]

