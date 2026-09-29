# =============================================================================
# schizoplacenta — revised analysis (covariate-adjusted, pregnancy-paired)
# 01_differential_expression.R: primary model, sensitivity models, sex model
# Source via run_all.R from the project root (folder that contains 0_TableCode.tsv).
# =============================================================================
fit_model <- function(dom, outcome = c("group", "continuous", "group24"), adjust = TRUE, rin = FALSE, drop = NULL) {
  outcome <- match.arg(outcome)
  m <- copy(meta); d <- dge
  if (!is.null(drop)) { k <- !(m$participant %in% drop); m <- m[k]; d <- calcNormFactors(d[, k, keep.lib.sizes = FALSE]); m[, participant := droplevels(participant)] }
  if (rin) { k <- !is.na(m$rin); m <- m[k]; d <- calcNormFactors(d[, k, keep.lib.sizes = FALSE]); m[, participant := droplevels(participant)]; m[, rin_z := as.numeric(scale(rin))] }
  m[, x := switch(outcome, group = get(paste0(dom, "_grp")), group24 = get(paste0(dom, "24_grp")), continuous = as.numeric(scale(get(dom))))]
  rhs <- "face * x"; if (adjust) rhs <- paste(rhs, "+ age_z + sex"); if (rin) rhs <- paste(rhs, "+ rin_z")
  design <- model.matrix(as.formula(paste("~", rhs)), m); colnames(design) <- make.names(colnames(design))
  stopifnot(qr(design)$rank == ncol(design))
  v0 <- voom(d, design); dc <- duplicateCorrelation(v0, design, block = m$participant)
  v  <- voom(d, design, block = m$participant, correlation = dc$consensus.correlation)
  fit <- lmFit(v, design, block = m$participant, correlation = dc$consensus.correlation)
  xn <- if (outcome == "continuous") "x" else "xDD"; ix <- paste0("facefetal.", xn)
  cm <- makeContrasts(contrasts = c(paste0(xn, "+", ix), xn, ix), levels = design); colnames(cm) <- c("fetal", "maternal", "interaction")
  fit2 <- eBayes(contrasts.fit(fit, cm), robust = TRUE)
  tabs <- rbindlist(lapply(colnames(cm), function(cc) as.data.table(topTable(fit2, coef = cc, number = Inf, sort.by = "P"), keep.rownames = "gene")[, contrast := cc]))
  tabs[, domain := dom]
  list(tabs = tabs, v = v, design = design, cm = cm, meta = m, corr = dc$consensus.correlation, fit = fit2)
}
summ <- function(t) t[, .(nominal_P05 = sum(P.Value < 0.05), nominal_P05_absLFC1 = sum(P.Value < 0.05 & abs(logFC) >= 1),
                          FDR05 = sum(adj.P.Val < 0.05), FDR05_absLFC1 = sum(adj.P.Val < 0.05 & abs(logFC) >= 1)), by = .(domain, contrast)]

primary <- lapply(names(DOM), function(d) fit_model(d)); names(primary) <- names(DOM)
prim_tabs <- rbindlist(lapply(primary, `[[`, "tabs"))
for (d in names(DOM)) for (cc in c("fetal", "maternal", "interaction")) wt(prim_tabs[domain == d & contrast == cc], sprintf("DE_primary_%s_%s.tsv", DOM[d], cc))
prim_summary <- summ(prim_tabs)[, model := "primary adjusted paired (DD vs TD)"]
diag <- rbindlist(lapply(names(DOM), function(d) data.table(domain = DOM[d], n_DD = uniqueN(primary[[d]]$meta[x == "DD", participant]), n_TD = uniqueN(primary[[d]]$meta[x == "TD", participant]),
  design_columns = ncol(primary[[d]]$design), rank = qr(primary[[d]]$design)$rank, residual_df = nrow(primary[[d]]$design) - ncol(primary[[d]]$design), pregnancy_correlation = round(primary[[d]]$corr, 3))))
wt(diag, "S_model_diagnostics.tsv")

# sensitivity models
sens <- list()
for (d in names(DOM)) {
  sens[[paste0(d, "_continuous")]]  <- fit_model(d, "continuous")
  sens[[paste0(d, "_unadjusted")]]  <- fit_model(d, "group", adjust = FALSE)
  sens[[paste0(d, "_RIN")]]         <- fit_model(d, "group", rin = TRUE)
  sens[[paste0(d, "_closest24")]]   <- fit_model(d, "group24")
  sens[[paste0(d, "_LOO_RN01")]]    <- fit_model(d, "group", drop = "RN01")
  sens[[paste0(d, "_cont_LOO_RN01")]] <- fit_model(d, "continuous", drop = "RN01")
}
sens_tabs <- rbindlist(lapply(names(sens), function(n) sens[[n]]$tabs[, model := sub("^[a-z]+_", "", n)]))
for (n in names(sens)) { d <- sub("_.*", "", n); wt(sens[[n]]$tabs[contrast == "fetal"], sprintf("DE_sensitivity_%s_%s_fetal.tsv", DOM[d], sub("^[a-z]+_", "", n))) }
sens_summary <- rbindlist(list(prim_summary, sens_tabs[, c(summ(.SD)), by = model][, .(domain, contrast, nominal_P05, nominal_P05_absLFC1, FDR05, FDR05_absLFC1, model)]))
sens_summary[, domain := DOM[domain]]
wt(sens_summary, "S_DE_summary_all_models.tsv")

# stability of fetal-facing effects: adjusted vs unadjusted, with vs without RN01
stab <- rbindlist(lapply(names(DOM), function(d) {
  a <- prim_tabs[domain == d & contrast == "fetal", .(gene, adj = logFC, t_adj = t)]
  u <- sens[[paste0(d, "_unadjusted")]]$tabs[contrast == "fetal", .(gene, unadj = logFC)]
  l <- sens[[paste0(d, "_LOO_RN01")]]$tabs[contrast == "fetal", .(gene, loo = logFC, t_loo = t)]
  z <- merge(merge(a, u, by = "gene"), l, by = "gene")
  data.table(domain = DOM[d], rho_adjusted_vs_unadjusted = cor(z$adj, z$unadj, method = "spearman"), sign_agreement_adj_unadj = mean(sign(z$adj) == sign(z$unadj)),
             rho_with_vs_without_RN01 = cor(z$adj, z$loo, method = "spearman"), rho_t_with_vs_without_RN01 = cor(z$t_adj, z$t_loo, method = "spearman"))
}))
wt(stab, "S_effect_stability.tsv")

# genes highlighted in the previous version, across models
key <- c("NTRK2", "SCGN", "HLA-DRB5", "TNFSF14", "MEGF10", "WASH1", "GRIN2B", "DGCR8", "MAPK3", "FURIN", "ABT1", "BTNL9", "PKD1P5", "USP32P1", "CNTNAP3P2", "KCNQ3", "EGR2", "ASTN1", "KCNN2", "NLGN3", "TECTA", "CRHBP")
keytab <- rbind(prim_tabs[gene %in% key & contrast == "fetal", .(model = "primary", domain, gene, logFC, P.Value, adj.P.Val)],
                sens_tabs[gene %in% key & contrast == "fetal", .(model, domain, gene, logFC, P.Value, adj.P.Val)])
keytab[, domain := DOM[domain]]
wt(keytab[order(gene, domain, model)], "S_previously_highlighted_genes_all_models.tsv")
# raw fetal-facing expression of those genes by participant
rawkey <- as.data.table(t(logcpm[intersect(key, rownames(logcpm)), meta[face == "fetal", sample]]), keep.rownames = "sample")
rawkey <- merge(meta[face == "fetal", .(sample, participant, cog, lang, mot)], rawkey, by = "sample")[order(cog)]
wt(rawkey, "S_fetal_log2CPM_previously_highlighted_genes.tsv")

# ------------------------------------------------------------------ sex model
ms <- copy(meta); ds <- model.matrix(~ sex + face + age_z, ms)
vs0 <- voom(dge, ds); dcs <- duplicateCorrelation(vs0, ds, block = ms$participant)
vs <- voom(dge, ds, block = ms$participant, correlation = dcs$consensus.correlation)
fs <- eBayes(lmFit(vs, ds, block = ms$participant, correlation = dcs$consensus.correlation), robust = TRUE)
sex_tab <- as.data.table(topTable(fs, coef = "sexM", number = Inf, sort.by = "P"), keep.rownames = "gene")
wt(sex_tab, "DE_sex_adjusted.tsv")
sex_sum <- data.table(model = "sex ~ sex + face + maternal_age, pregnancy-blocked", FDR05 = sum(sex_tab$adj.P.Val < 0.05), FDR05_absLFC1 = sum(sex_tab$adj.P.Val < 0.05 & abs(sex_tab$logFC) >= 1),
                      top10 = paste(sex_tab$gene[1:10], collapse = ";"), pregnancy_correlation = round(dcs$consensus.correlation, 3))
wt(sex_sum, "S_sex_model_summary.tsv")

saveRDS(list(primary = primary, sens = sens, prim_tabs = prim_tabs, sens_tabs = sens_tabs, logcpm = logcpm, meta = meta), file.path(RDS, "de_models.rds"))
