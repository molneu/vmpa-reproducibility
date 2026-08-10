# Supplementary Figure 2: foundation-model analyses

This folder contains the scripts used to compare classical enrichment scores with BulkFormer, Geneformer, and UCE-derived scores for the AKT1 analyses shown in Supplementary Figure 2.

The scripts do not include pretrained model weights, H5AD objects, or intermediate embedding outputs. These inputs are retained in the local foundation-model analysis archive and are passed to the runner as a separate data-export directory.

## Analyses

- `score_bulkformer_activity.py` generates rank-weighted BulkFormer sample embeddings and computes cosine and gene-membership Spearman scores. Scores are compared with measured pAKT(Ser473) using linear regression.
- `summarize_geneformer_activity.py` summarizes Geneformer cell-embedding cosine scores for vehicle and 10 uM ZSTK474-treated A172 and U87 cells.
- `summarize_geneformer_gene_level.py` summarizes Geneformer gene-embedding cosine and gene-membership Spearman scores for the same comparisons.
- `score_uce_activity.py` generates continuous AKT1 signature directions separately within A172 and U87 UCE embedding spaces. Each direction is estimated using all cells from the covariance between L2-normalized cell embeddings and mean expression of the mapped signature genes. Cell scores are projections onto the corresponding direction.
- `plot_final_foundation_model_pheatmaps.R` generates the BulkFormer, UCE, and Geneformer significance heatmaps using `pheatmap`.

Single-cell comparisons use vehicle and 10 uM ZSTK474-treated cells and one-sided Mann-Whitney tests in the expected direction of lower AKT1 KO-down signature scores after treatment. AUCell and GSVA statistics are supplied as classical enrichment baselines.

## Running the analysis

The runner expects the local data-export directory used for the revision analysis and a local BulkFormer checkout containing the pretrained gene embeddings.

```bash
BULKFORMER_PYTHON=/path/to/bulkformer/python \
SINGLECELL_PYTHON=/path/to/geneformer/python \
RSCRIPT=/path/to/Rscript \
bash run_supplementary_figure2_foundation_models.sh \
  /path/to/embedding_methods_primary_data_export \
  /path/to/BulkFormer \
  local_outputs
```

Python requirements include `anndata`, `numpy`, `pandas`, `scipy`, `scikit-learn`, and `torch`. R requirements are `pheatmap` and `grid`.

Generated tables and plots should be written to `local_outputs/`, which is excluded from version control.
