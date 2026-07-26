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
