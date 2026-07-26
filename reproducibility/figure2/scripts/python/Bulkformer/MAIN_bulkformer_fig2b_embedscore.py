#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import torch
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import rankdata, spearmanr
import seaborn as sns

# ------------------------------------------------------------
# CONFIG —  paths
# ------------------------------------------------------------
BULKFORMER_ROOT = Path(os.environ.get("BULKFORMER_ROOT", "~/BulkFormer")).expanduser()
RESULTS_DIR = Path(os.environ.get("BULKFORMER_RESULTS_DIR", "reproducibility/figure2/results/bulkformer")).expanduser()
OUT_DIR = RESULTS_DIR

EXPR_CSV = RESULTS_DIR / "fig2b_expr.csv"
META_CSV = RESULTS_DIR / "fig2b_meta.csv"
GMT_FILE = RESULTS_DIR / "fig2b_gene_sets.gmt"
GENE_EMB_PATH  = BULKFORMER_ROOT / "data" / "esm2_feature_concat.pt"
GENE_INFO_PATH = BULKFORMER_ROOT / "data" / "bulkformer_gene_info.csv"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def read_gmt(path: Path):
    gs = {}
    with open(path, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name, _desc, *genes = parts
            gs[name] = list(dict.fromkeys(g.upper() for g in genes if g))
    return gs

def ensure_exists(p, label):
    if not p.exists():
        print(f"[ERROR] Missing {label}: {p}")
        sys.exit(1)

def l2_normalize(v):
    n = np.linalg.norm(v) + 1e-8
    return v / n

def cosine_sim(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------
for p, lab in [(BULKFORMER_ROOT, "BULKFORMER_ROOT"),
               (RESULTS_DIR, "RESULTS_DIR"),
               (EXPR_CSV, "EXPR_CSV"),
               (META_CSV, "META_CSV"),
               (GMT_FILE, "GMT_FILE"),
               (GENE_EMB_PATH, "GENE_EMB_PATH"),
               (GENE_INFO_PATH, "GENE_INFO_PATH")]:
    ensure_exists(p, lab)

expr_df = pd.read_csv(EXPR_CSV)
expr_df["gene_symbol"] = expr_df["gene_symbol"].astype(str).str.upper()
expr_df = expr_df.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")

meta_df = pd.read_csv(META_CSV)
meta_df["sample"] = meta_df["sample"].astype(str)
gene_sets = read_gmt(GMT_FILE)

# ------------------------------------------------------------
# Load BulkFormer embeddings + mapping
# ------------------------------------------------------------
gene_emb_t = torch.load(GENE_EMB_PATH, map_location="cpu")
if isinstance(gene_emb_t, dict) and "emb" in gene_emb_t:
    gene_emb_t = gene_emb_t["emb"]
gene_emb = gene_emb_t.detach().cpu().float().numpy()
embed_dim = gene_emb.shape[1]

gi = pd.read_csv(GENE_INFO_PATH)
symbol_col = next((c for c in gi.columns if c.lower() in
                   ("symbol","gene_symbol","hgnc_symbol","gene","gene_name")),
                   gi.columns[0])
gi[symbol_col] = gi[symbol_col].astype(str).str.upper()
gene_index_to_symbol = gi[symbol_col].tolist()
symbol_to_index = {s: i for i, s in enumerate(gene_index_to_symbol)}

# ------------------------------------------------------------
# Center + normalize embeddings
# ------------------------------------------------------------
mu = gene_emb.mean(axis=0, keepdims=True)
gene_emb_c = gene_emb - mu
gene_emb_c /= np.linalg.norm(gene_emb_c, axis=1, keepdims=True) + 1e-8

# ------------------------------------------------------------
# Scoring methods
# ------------------------------------------------------------
def cosine_score(sample_vec, sig_idx, gene_emb_c):
    vec = gene_emb_c[sig_idx, :].mean(axis=0)
    return cosine_sim(sample_vec, vec)

def spearman_score(sample_vec, sig_idx, gene_emb_c, align_direction=True):
    proj = gene_emb_c @ sample_vec
    mask = np.zeros(len(proj))
    mask[sig_idx] = 1.0
    r, _ = spearmanr(proj, mask)
    if not np.isfinite(r):
        r = 0.0
    return -r if align_direction else r  # flip sign for consistency

def gsea_projection_score(sample_vec, sig_idx, gene_emb_c):
    proj = gene_emb_c @ sample_vec
    proj = proj + np.random.normal(0, 1e-9, size=proj.shape)
    ranked_idx = np.argsort(proj)[::-1]
    hits = np.isin(ranked_idx, sig_idx)
    Nh = hits.sum()
    N = len(proj)
    if Nh == 0:
        return 0.0
    Phit = np.cumsum(hits / Nh)
    Pmiss = np.cumsum(~hits / (N - Nh))
    ES = np.max(Phit - Pmiss)
    ES_neg = np.min(Phit - Pmiss)
    return ES if abs(ES) > abs(ES_neg) else ES_neg



methods = {
    "cosine": cosine_score,
    "spearman": spearman_score,
    "gsea_projection": gsea_projection_score
}

# ------------------------------------------------------------
# Compute signed_ranked variant only
# ------------------------------------------------------------
rows_all, regs_all, cov_all = [], [], []

samples = [s for s in meta_df["sample"] if s in expr_df.columns]
common_symbols = [g for g in expr_df.index if g in symbol_to_index]
expr_sub = expr_df.loc[common_symbols, samples].copy()
idx_all = np.array([symbol_to_index[g] for g in common_symbols], int)
emb_sub = gene_emb_c[idx_all, :]

variant = "signed_ranked"
print(f"\n=== Variant: {variant} ===")

for method_name, func in methods.items():
    print(f"Running {method_name}...")
    sig_cov = []
    sig_idx_dict = {}
    for sig_name, genes in gene_sets.items():
        idx = [symbol_to_index[g] for g in genes if g in symbol_to_index]
        cov = len(idx) / max(1, len(genes))
        sig_cov.append((sig_name, len(genes), len(idx), cov))
        sig_idx_dict[sig_name] = idx
    cov_df = pd.DataFrame(sig_cov, columns=["signature","n_genes","n_mapped","coverage"])
    cov_df["method"] = method_name
    cov_df["variant"] = variant
    cov_all.append(cov_df)

    # compute scores
    rows = []
    for s in samples:
        v = expr_sub[s].values.astype(np.float32)
        m, sd = v.mean(), v.std(ddof=0) or 1.0
        z = (v - m) / sd
        ranks = z.argsort().argsort().astype(np.float32)
        w = (ranks + 1) / (len(ranks) + 1)
        w /= w.sum()
        sv = (emb_sub * w[:, None]).sum(axis=0)
        sv = l2_normalize(sv)

        for sig, idx in sig_idx_dict.items():
            if not idx:
                continue
            gv = gene_emb_c[idx, :].mean(axis=0)
            gv = -l2_normalize(gv)  # all AKT1 KO are DOWN, invert direction
            if method_name == "spearman":
                sc = func(sv, idx, gene_emb_c, align_direction=True)
            else:
                sc = func(sv, idx, gene_emb_c)
            rows.append((s, sig, method_name, variant, sc))
    df_scores = pd.DataFrame(rows, columns=["sample","signature","method","variant","score"])
    df_scores = df_scores.merge(meta_df, on="sample", how="left")
    rows_all.append(df_scores)

    # regression vs reference
    regs = []
    for sig in df_scores["signature"].unique():
        sub = df_scores[df_scores["signature"] == sig].dropna(subset=["score","reference"])
        if len(sub) < 3:
            continue
        x = sub["reference"].values.reshape(-1,1)
        y = sub["score"].values
        lr = LinearRegression().fit(x, y)
        yhat = lr.predict(x)
        r2 = r2_score(y, yhat)
        slope = lr.coef_[0]
        n = len(x)
        se = np.sqrt(np.sum((y - yhat)**2)/(n-2)) / np.sqrt(np.sum((x-x.mean())**2))
        t = slope / se
        p = 2*(1 - st.t.cdf(abs(t), df=n-2))
        regs.append((sig, method_name, variant, r2, slope, p))
    regs_df = pd.DataFrame(regs, columns=["signature","method","variant","R2","slope","p_value"])
    regs_all.append(regs_df)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------
scores_df = pd.concat(rows_all, ignore_index=True)
cov_df = pd.concat(cov_all, ignore_index=True)
reg_df = pd.concat(regs_all, ignore_index=True)

scores_df.to_csv(OUT_DIR / "bulkformer_all_scores.csv", index=False)
cov_df.to_csv(OUT_DIR / "bulkformer_all_coverage.csv", index=False)
reg_df.to_csv(OUT_DIR / "bulkformer_all_regression.csv", index=False)

print("\n✅ Finished all methods (signed_ranked only).")
print(f"Scores → {OUT_DIR/'bulkformer_all_scores.csv'}")
print(f"Regression → {OUT_DIR/'bulkformer_all_regression.csv'}")
print(f"Coverage → {OUT_DIR/'bulkformer_all_coverage.csv'}")

# ------------------------------------------------------------
# === VISUALIZATION: per-method R2 and scatter plots ==========
# ------------------------------------------------------------
from matplotlib.backends.backend_pdf import PdfPages

pdf_path = OUT_DIR / "bulkformer_method_comparison_plots.pdf"
print(f"\n📊 Generating comparison plots → {pdf_path}")

with PdfPages(pdf_path) as pdf:
    # --- 1️⃣ R2 barplots, one per method ---
    for method_name in reg_df["method"].unique():
        sub = reg_df[reg_df["method"] == method_name].copy()
        if sub.empty:
            continue
        plt.figure(figsize=(5, 4))
        sns.barplot(data=sub, x="signature", y="R2", color="steelblue")
        plt.title(f"R² per signature — {method_name}")
        plt.ylim(0, 1)  # fixed y-axis
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("R²")
        plt.xlabel("Signature")
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # --- 2️⃣ Scatter plots (score vs reference) ---
    for method_name in reg_df["method"].unique():
        sub_scores = scores_df[scores_df["method"] == method_name].copy()
        if sub_scores.empty:
            continue
        sigs = sub_scores["signature"].unique()
        for sig in sigs:
            sig_sub = sub_scores[sub_scores["signature"] == sig]
            if sig_sub.empty:
                continue
            plt.figure(figsize=(4, 4))

            # seaborn regression plot with 95% CI
            sns.regplot(
                data=sig_sub,
                x="reference",
                y="score",
                scatter_kws={"color": "steelblue", "alpha": 0.8, "s": 40},
                line_kws={"color": "black", "lw": 1.2},
                ci=95
            )

            # compute R² and p-value for title
            if len(sig_sub) >= 3:
                lr = LinearRegression().fit(sig_sub["reference"].values.reshape(-1, 1), sig_sub["score"].values)
                yhat = lr.predict(sig_sub["reference"].values.reshape(-1, 1))
                r2 = r2_score(sig_sub["score"].values, yhat)
                p = reg_df.loc[
                    (reg_df["method"] == method_name) & (reg_df["signature"] == sig),
                    "p_value"
                ].values
                p_str = f"p={p[0]:.3g}" if len(p) else ""
                plt.title(f"{sig}\n{method_name} — R²={r2:.2f} {p_str}")
            else:
                plt.title(f"{sig}\n{method_name}")

            plt.xlabel("pAKT reference")
            plt.ylabel("Enrichment score")
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
print(f"✅ All comparison plots with 95% CI saved in → {pdf_path}")


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# Assuming you already have: gene_emb_c, gene_sets, symbol_to_index, l2_normalize

# 1️⃣ Compute signature vectors
sig_vecs = {}
for sig_name, genes in gene_sets.items():
    idx = [symbol_to_index[g] for g in genes if g in symbol_to_index]
    if len(idx) == 0:
        continue
    gv = gene_emb_c[idx, :].mean(axis=0)
    gv = -l2_normalize(gv)  # AKT1 KO are downregulated
    sig_vecs[sig_name] = gv

sig_names = list(sig_vecs.keys())
mat = np.vstack([sig_vecs[s] for s in sig_names])

# 2️⃣ Compute pairwise cosine similarities
cosine_mat = np.dot(mat, mat.T)
# Ensure numerical stability
cosine_mat = np.clip(cosine_mat, -1.0, 1.0)

# 3️⃣ Compute angular distances (degrees)
angles = np.degrees(np.arccos(cosine_mat))
cosine_df = pd.DataFrame(cosine_mat, index=sig_names, columns=sig_names)
angle_df = pd.DataFrame(angles, index=sig_names, columns=sig_names)

# Save to file
cosine_df.to_csv(OUT_DIR / "bulkformer_signature_cosine_similarity.csv")
angle_df.to_csv(OUT_DIR / "bulkformer_signature_angles.csv")

# 4a️⃣ Heatmap of cosine similarities
plt.figure(figsize=(6, 5))
sns.heatmap(cosine_df, annot=True, fmt=".2f", cmap="vlag", vmin=-1, vmax=1)
plt.title("BulkFormer AKT1 signature similarity (cosine)")
plt.tight_layout()
plt.savefig(OUT_DIR / "bulkformer_signature_cosine_heatmap.pdf")
plt.close()

# 4b️⃣ PCA projection of signature directions
pca = PCA(n_components=2)
coords = pca.fit_transform(mat)
pca_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=sig_names)

plt.figure(figsize=(5, 5))
sns.scatterplot(x="PC1", y="PC2", data=pca_df, s=120)
for sig in sig_names:
    plt.text(pca_df.loc[sig, "PC1"], pca_df.loc[sig, "PC2"], sig, fontsize=9, ha="center", va="center", weight="bold")
plt.title("AKT1 signature directions in BulkFormer space")
plt.tight_layout()
plt.savefig(OUT_DIR / "bulkformer_signature_PCA.pdf")
plt.close()

print("✅ Saved context-drift analyses:")
print("  - Cosine matrix → bulkformer_signature_cosine_similarity.csv")
print("  - Angle matrix → bulkformer_signature_angles.csv")
print("  - Heatmap → bulkformer_signature_cosine_heatmap.pdf")
print("  - PCA plot → bulkformer_signature_PCA.pdf")
