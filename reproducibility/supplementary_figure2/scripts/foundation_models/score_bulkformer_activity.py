#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import linregress, spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--gmt", type=Path, required=True)
    parser.add_argument("--gene-embeddings", type=Path, required=True)
    parser.add_argument("--gene-info", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def l2_normalize(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    denominator = np.linalg.norm(values, axis=axis, keepdims=axis is not None) + 1e-8
    return values / denominator


def read_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    gene_column = "gene_symbol" if "gene_symbol" in frame.columns else frame.columns[0]
    frame = frame.rename(columns={gene_column: "gene_symbol"})
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
    frame = frame.drop_duplicates("gene_symbol").set_index("gene_symbol")
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def read_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name, _description, *genes = fields
            gene_sets[name] = list(dict.fromkeys(gene.strip().upper() for gene in genes if gene.strip()))
    return gene_sets


def load_gene_embeddings(embedding_path: Path, gene_info_path: Path) -> tuple[np.ndarray, dict[str, int]]:
    loaded = torch.load(embedding_path, map_location="cpu", weights_only=False)
    if isinstance(loaded, dict):
        if "emb" not in loaded:
            raise KeyError(f"No 'emb' matrix found in {embedding_path}")
        loaded = loaded["emb"]
    embeddings = loaded.detach().cpu().float().numpy()
    embeddings = embeddings - embeddings.mean(axis=0, keepdims=True)
    embeddings = l2_normalize(embeddings, axis=1)

    gene_info = pd.read_csv(gene_info_path)
    candidates = {"symbol", "gene_symbol", "hgnc_symbol", "gene", "gene_name"}
    symbol_column = next(
        (column for column in gene_info.columns if column.lower() in candidates),
        gene_info.columns[0],
    )
    symbols = gene_info[symbol_column].astype(str).str.upper().tolist()
    if len(symbols) != len(embeddings):
        raise ValueError("Gene information and embedding matrix have different row counts")
    return embeddings, {symbol: index for index, symbol in enumerate(symbols)}


def expression_weighted_embeddings(
    expression: pd.DataFrame,
    gene_embeddings: np.ndarray,
    symbol_to_index: dict[str, int],
) -> tuple[list[str], np.ndarray]:
    symbols = [symbol for symbol in expression.index if symbol in symbol_to_index]
    indices = np.array([symbol_to_index[symbol] for symbol in symbols], dtype=int)
    mapped_embeddings = gene_embeddings[indices]
    values = expression.loc[symbols].to_numpy(np.float32).T
    sample_embeddings = []
    for row in values:
        mean = row.mean()
        standard_deviation = row.std(ddof=0) or 1.0
        standardized = (row - mean) / standard_deviation
        ranks = standardized.argsort().argsort().astype(np.float32)
        weights = (ranks + 1.0) / (len(ranks) + 1.0)
        weights /= weights.sum()
        sample_embeddings.append(weights @ mapped_embeddings)
    return list(expression.columns), l2_normalize(np.vstack(sample_embeddings), axis=1)


def signature_indices(
    gene_sets: dict[str, list[str]],
    symbol_to_index: dict[str, int],
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    indices_by_signature: dict[str, np.ndarray] = {}
    coverage_rows = []
    for signature, genes in gene_sets.items():
        indices = np.array(
            [symbol_to_index[gene] for gene in genes if gene in symbol_to_index],
            dtype=int,
        )
        coverage_rows.append(
            {
                "signature": signature,
                "n_genes": len(genes),
                "n_mapped": len(indices),
                "coverage": len(indices) / max(1, len(genes)),
            }
        )
        if len(indices):
            indices_by_signature[signature] = indices
    return indices_by_signature, pd.DataFrame(coverage_rows)


def score_default_samples(
    sample_names: list[str],
    sample_embeddings: np.ndarray,
    indices_by_signature: dict[str, np.ndarray],
    gene_embeddings: np.ndarray,
) -> pd.DataFrame:
    rows = []
    for sample, sample_embedding in zip(sample_names, sample_embeddings, strict=True):
        gene_projections = gene_embeddings @ sample_embedding
        for signature, indices in indices_by_signature.items():
            signature_vector = gene_embeddings[indices].mean(axis=0)
            cosine_score = float(
                np.dot(sample_embedding, signature_vector)
                / (
                    np.linalg.norm(sample_embedding) * np.linalg.norm(signature_vector)
                    + 1e-8
                )
            )
            membership = np.zeros(len(gene_projections), dtype=float)
            membership[indices] = 1.0
            membership_spearman = spearmanr(gene_projections, membership).statistic
            if not np.isfinite(membership_spearman):
                membership_spearman = 0.0
            rows.append(
                (
                    "BulkFormer",
                    sample,
                    signature,
                    "cosine",
                    cosine_score,
                )
            )
            rows.append(
                (
                    "BulkFormer",
                    sample,
                    signature,
                    "gene_membership_spearman",
                    float(-membership_spearman),
                )
            )
    return pd.DataFrame(rows, columns=["model", "sample", "signature", "method", "score"])


def regression_statistics(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, signature, method), group in scores.groupby(
        ["model", "signature", "method"], sort=False
    ):
        group = group.dropna(subset=["score", "reference"])
        if len(group) < 3:
            continue
        fit = linregress(group["reference"].to_numpy(float), group["score"].to_numpy(float))
        two_sided_p = float(fit.pvalue)
        expected_p = two_sided_p / 2.0 if fit.slope > 0 else 1.0 - two_sided_p / 2.0
        expected_p = max(expected_p, np.finfo(float).tiny)
        rows.append(
            {
                "model": model,
                "signature": signature,
                "method": method,
                "n": len(group),
                "R2": float(fit.rvalue**2),
                "slope": float(fit.slope),
                "signed_R2": float(np.sign(fit.slope) * fit.rvalue**2),
                "p_value_two_sided": two_sided_p,
                "neglog10_p_two_sided": float(-np.log10(max(two_sided_p, np.finfo(float).tiny))),
                "expected_direction": "positive_slope_with_pAKT",
                "expected_direction_hit": bool(fit.slope > 0),
                "expected_one_tailed_p": expected_p,
                "expected_neglog10_p": float(-np.log10(expected_p)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    expression = read_expression(args.expression)
    metadata = pd.read_csv(args.metadata)
    metadata["sample"] = metadata["sample"].astype(str)
    gene_sets = read_gmt(args.gmt)
    gene_embeddings, symbol_to_index = load_gene_embeddings(args.gene_embeddings, args.gene_info)

    sample_names, sample_embeddings = expression_weighted_embeddings(
        expression, gene_embeddings, symbol_to_index
    )
    indices_by_signature, coverage = signature_indices(gene_sets, symbol_to_index)
    scores = score_default_samples(
        sample_names,
        sample_embeddings,
        indices_by_signature,
        gene_embeddings,
    )
    scores = scores.merge(metadata, on="sample", how="left", validate="many_to_one")
    if scores["reference"].isna().any():
        raise ValueError("Some expression samples could not be matched to pAKT metadata")
    statistics = regression_statistics(scores)

    args.outdir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.outdir / "bulkformer_activity_scores.csv", index=False)
    statistics.to_csv(args.outdir / "bulkformer_activity_stats.csv", index=False)
    coverage.to_csv(args.outdir / "bulkformer_signature_coverage.csv", index=False)
    print(statistics.to_string(index=False))


if __name__ == "__main__":
    main()
