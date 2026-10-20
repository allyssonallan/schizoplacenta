#!/usr/bin/env python3
"""schizoplacenta — pre-ranked GSEA (gseapy, same algorithm as GSEA/fgsea: weighted KS statistic, gene permutations)
Collections are processed in chunks of 400 sets to bound memory; nominal P values come from the permutations and the
FDR is Benjamini-Hochberg across the whole collection (as in fgsea). Each ranking x collection result is cached in
results_revision/rds/gsea_parts/ so an interrupted run resumes.
usage: python3 R/02b_gsea_prerank.py <project_root> [nperm]
Reads results_revision/rds/gsea_rankings.tsv (moderated t, DD - TD), the two C2 GMTs in the project root,
results_revision/rds/C7.gmt and external_sets.gmt; writes results_revision/tables/GSEA_prerank_all.tsv.
C7 is run for the primary model only (fetal- and maternal-facing), C2 all for the fetal-facing primary contrast only; C2:CP and the published sets for all models."""
import sys, os, warnings, pandas as pd, gseapy as gp
warnings.filterwarnings("ignore")
root = sys.argv[1]; nperm = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
rds = os.path.join(root, "results_revision", "rds"); tab = os.path.join(root, "results_revision", "tables")
rk = pd.read_csv(os.path.join(rds, "gsea_rankings.tsv"), sep="\t").set_index("gene")
colls = {"C2CP": os.path.join(root, "c2.cp.v7.1.symbols.gmt"), "C2ALL": os.path.join(root, "c2.all.v7.1.symbols.gmt"),
         "external": os.path.join(rds, "external_sets.gmt"), "C7": os.path.join(rds, "C7.gmt")}
colls = {k: v for k, v in colls.items() if os.path.exists(v)}
def read_gmt(f):
    sets = {}
    for line in open(f):
        p = line.rstrip("\n").split("\t"); sets[p[0]] = [g for g in p[2:] if g]
    return sets
parts = os.path.join(rds, "gsea_parts"); os.makedirs(parts, exist_ok=True)
out = []
for col in rk.columns:
    domain, model, contrast = col.split("|")
    rnk = rk[col].dropna().sort_values(ascending=False)
    for name, gmt in colls.items():
        if name in ("C2ALL", "C7") and model != "primary": continue
        if name == "C2ALL" and contrast != "fetal": continue  # C2 all: fetal-facing primary only (Datasets 5/8/11)
        pf = os.path.join(parts, f"{domain}_{model}_{contrast}_{name}.tsv")
        if os.path.exists(pf): out.append(pd.read_csv(pf, sep="\t")); continue
        sets = read_gmt(gmt); names = list(sets); chunks = [dict((k, sets[k]) for k in names[i:i + 400]) for i in range(0, len(names), 400)]
        res = []
        for ch in chunks:
            r = gp.prerank(rnk=rnk, gene_sets=ch, min_size=5, max_size=1000, permutation_num=nperm, seed=1, threads=2, outdir=None, verbose=False).res2d
            res.append(r)
        res = pd.concat(res).rename(columns={"Term": "gene_set", "NOM p-val": "PValue", "Lead_genes": "leading_edge"})
        # BH across the whole collection (P floor = 1/nperm)
        pv = res["PValue"].clip(lower=1.0 / nperm).values; n = len(pv); order = pv.argsort(); ranked = pv[order] * n / (pd.Series(range(1, n + 1)).values)
        fdr = pd.Series(ranked).iloc[::-1].cummin().iloc[::-1].clip(upper=1).values; out_fdr = pd.Series(index=order, data=fdr).sort_index().values
        res["FDR"] = out_fdr; res["PValue"] = pv
        res["NGenes"] = res["Tag %"].str.split("/").str[1].astype(int)
        res["Direction"] = ["Up" if v > 0 else "Down" for v in res["NES"]]
        res = res.assign(domain=domain, model=model, contrast=contrast, collection=name)[["domain", "model", "collection", "contrast", "gene_set", "NGenes", "ES", "NES", "Direction", "PValue", "FDR", "leading_edge"]]
        res.to_csv(pf, sep="\t", index=False)
        out.append(res); print(f"{domain} {model} {contrast} {name}: {int((res.FDR < 0.05).sum())} sets FDR<0.05", flush=True)
pd.concat(out).to_csv(os.path.join(tab, "GSEA_prerank_all.tsv"), sep="\t", index=False)
print("done")
