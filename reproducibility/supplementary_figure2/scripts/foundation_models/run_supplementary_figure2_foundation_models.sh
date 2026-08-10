#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
  printf 'Usage: %s <data_export_dir> <bulkformer_dir> <output_dir>\n' "$0" >&2
  exit 1
fi

DATA_EXPORT_DIR="$1"
BULKFORMER_DIR="$2"
OUTPUT_DIR="$3"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BULKFORMER_PYTHON="${BULKFORMER_PYTHON:-python3}"
SINGLECELL_PYTHON="${SINGLECELL_PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"

mkdir -p \
  "$OUTPUT_DIR/tables/bulkformer" \
  "$OUTPUT_DIR/tables/geneformer_cell" \
  "$OUTPUT_DIR/tables/geneformer_gene" \
  "$OUTPUT_DIR/tables/uce4" \
  "$OUTPUT_DIR/tables/uce33" \
  "$OUTPUT_DIR/plots"

"$BULKFORMER_PYTHON" "$SCRIPT_DIR/score_bulkformer_activity.py" \
  --expression "$DATA_EXPORT_DIR/inputs/bulk/fig2b_expr.csv" \
  --metadata "$DATA_EXPORT_DIR/inputs/bulk/fig2b_meta.csv" \
  --gmt "$DATA_EXPORT_DIR/inputs/bulk/fig2b_gene_sets.gmt" \
  --gene-embeddings "$BULKFORMER_DIR/data/esm2_feature_concat.pt" \
  --gene-info "$BULKFORMER_DIR/data/bulkformer_gene_info.csv" \
  --outdir "$OUTPUT_DIR/tables/bulkformer"

"$SINGLECELL_PYTHON" "$SCRIPT_DIR/summarize_geneformer_activity.py" \
  --a172-h5ad "$DATA_EXPORT_DIR/inputs/single_cell/a172_geneformer_input.h5ad" \
  --u87-h5ad "$DATA_EXPORT_DIR/inputs/single_cell/u87_geneformer_input.h5ad" \
  --cell-scores "$DATA_EXPORT_DIR/primary_outputs/geneformer/geneformer_cell_cosine_scores.csv" \
  --gene-scores "$DATA_EXPORT_DIR/primary_outputs/geneformer/geneformer_genelevel_all_scores.csv" \
  --treated-dose 10 \
  --outdir "$OUTPUT_DIR/tables/geneformer_cell"

"$SINGLECELL_PYTHON" "$SCRIPT_DIR/summarize_geneformer_gene_level.py" \
  --a172-h5ad "$DATA_EXPORT_DIR/inputs/single_cell/a172_geneformer_input.h5ad" \
  --u87-h5ad "$DATA_EXPORT_DIR/inputs/single_cell/u87_geneformer_input.h5ad" \
  --gene-scores "$DATA_EXPORT_DIR/primary_outputs/geneformer/geneformer_genelevel_all_scores.csv" \
  --treated-dose 10 \
  --outdir "$OUTPUT_DIR/tables/geneformer_gene"

"$SINGLECELL_PYTHON" "$SCRIPT_DIR/score_uce_activity.py" \
  --a172 "$DATA_EXPORT_DIR/inputs/uce4_embeddings/a172_cleana172_geneformer_input_uce_adata.h5ad" \
  --u87 "$DATA_EXPORT_DIR/inputs/uce4_embeddings/u87_cleanu87_geneformer_input_uce_adata.h5ad" \
  --gmt "$DATA_EXPORT_DIR/inputs/single_cell/fig2b_gene_sets.gmt" \
  --model-label UCE-4 \
  --treated-dose 10 \
  --outdir "$OUTPUT_DIR/tables/uce4"

"$SINGLECELL_PYTHON" "$SCRIPT_DIR/score_uce_activity.py" \
  --a172 "$DATA_EXPORT_DIR/inputs/uce33_embeddings/a172_geneformer_input_uce_adata.h5ad" \
  --u87 "$DATA_EXPORT_DIR/inputs/uce33_embeddings/u87_geneformer_input_uce_adata.h5ad" \
  --gmt "$DATA_EXPORT_DIR/inputs/single_cell/fig2b_gene_sets.gmt" \
  --model-label UCE-33 \
  --treated-dose 10 \
  --outdir "$OUTPUT_DIR/tables/uce33"

"$RSCRIPT" "$SCRIPT_DIR/plot_final_foundation_model_pheatmaps.R" \
  "$DATA_EXPORT_DIR/primary_outputs/activity_model_comparison/gsva_repro_bulk_stats.csv" \
  "$OUTPUT_DIR/tables/bulkformer/bulkformer_activity_stats.csv" \
  "$DATA_EXPORT_DIR/primary_outputs/activity_model_comparison/aucell_repro_singlecell_stats.csv" \
  "$OUTPUT_DIR/tables/uce4/uce_4_cell_stats.csv" \
  "$OUTPUT_DIR/tables/uce33/uce_33_cell_stats.csv" \
  "$OUTPUT_DIR/tables/geneformer_cell/geneformer_cell_10uM_stats.csv" \
  "$OUTPUT_DIR/tables/geneformer_gene/geneformer_gene_level_10uM_stats.csv" \
  "$OUTPUT_DIR/plots"
