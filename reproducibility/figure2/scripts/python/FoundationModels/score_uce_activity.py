#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score


def read_gmt(path: Path) -> dict[str, list[str]]:
    signatures = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                signatures[fields[0]] = list(
                    dict.fromkeys(gene.strip().upper() for gene in fields[2:] if gene.strip())
                )
    return signatures


def normalize_rows(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.maximum(norms, np.finfo(values.dtype).eps)


def normalize_vector(values: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(values)
    return values / max(norm, np.finfo(values.dtype).eps)


def embedding_key(obj: ad.AnnData) -> str:
    if "X_uce" in obj.obsm:
        return "X_uce"
    matches = [key for key in obj.obsm if "uce" in key.lower()]
    if len(matches) != 1:
        raise KeyError(f"Expected one UCE embedding in obsm; found {list(obj.obsm.keys())}")
    return matches[0]


def mean_signature_expression(obj: ad.AnnData, genes: list[str]) -> tuple[np.ndarray | None, int]:
    gene_index = {str(gene).upper(): index for index, gene in enumerate(obj.var_names)}
    indices = [gene_index[gene] for gene in genes if gene in gene_index]
    if not indices:
        return None, 0
    values = obj.X[:, indices]
    means = np.asarray(values.mean(axis=1)).ravel() if sparse.issparse(values) else values.mean(axis=1)
    return np.asarray(means, dtype=float), len(indices)


def build_continuous_axes(
    obj: ad.AnnData,
    signatures: dict[str, list[str]],
    source: str,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    embeddings = normalize_rows(np.asarray(obj.obsm[embedding_key(obj)], dtype=np.float32))
    centered_embeddings = embeddings - embeddings.mean(axis=0, keepdims=True)
    axes = {}
    rows = []
    for signature, genes in signatures.items():
        scores, mapped = mean_signature_expression(obj, genes)
        if scores is None:
            rows.append(
                {
                    "source": source,
                    "signature": signature,
                    "n_signature_genes": len(genes),
                    "n_mapped_genes": 0,
                    "status": "no_mapped_genes",
                }
            )
            continue
        centered_scores = scores - scores.mean()
        denominator = float(np.dot(centered_scores, centered_scores))
        if denominator <= 0:
            rows.append(
                {
                    "source": source,
                    "signature": signature,
                    "n_signature_genes": len(genes),
                    "n_mapped_genes": mapped,
                    "status": "constant_score",
                }
            )
            continue
        axis = centered_embeddings.T @ centered_scores / denominator
        axes[signature] = normalize_vector(axis)
        rows.append(
            {
                "source": source,
                "signature": signature,
                "n_signature_genes": len(genes),
                "n_mapped_genes": mapped,
                "n_cells": obj.n_obs,
                "axis_orientation": "continuous_higher_signature_score",
                "axis_definition": "all_observations_continuous",
                "status": "ok",
            }
        )
    return axes, pd.DataFrame(rows)


def score_cells(
    dataset: str,
    obj: ad.AnnData,
    axes: dict[str, np.ndarray],
    model: str,
) -> pd.DataFrame:
    embeddings = normalize_rows(np.asarray(obj.obsm[embedding_key(obj)], dtype=np.float32))
    treatment_column = next((column for column in ("treatment", "condition") if column in obj.obs), None)
    dose_column = next((column for column in ("dose", "DOSE") if column in obj.obs), None)
    treatment = (
        obj.obs[treatment_column].astype(str).str.lower().to_numpy()
        if treatment_column
        else np.repeat("unknown", obj.n_obs)
    )
    dose = pd.to_numeric(obj.obs[dose_column], errors="coerce").to_numpy() if dose_column else np.repeat(np.nan, obj.n_obs)
    rows = []
    for signature, axis in axes.items():
        scores = embeddings @ axis
        rows.extend(
            {
                "model": model,
                "dataset": dataset,
                "cell_id": str(obj.obs_names[index]),
                "treatment": treatment[index],
                "dose": dose[index],
                "signature": signature,
                "score": float(scores[index]),
            }
            for index in range(obj.n_obs)
        )
    return pd.DataFrame(rows)


def summarize(
    scores: pd.DataFrame,
    control: str,
    treated: str,
    treated_dose: float,
) -> pd.DataFrame:
    rows = []
    for keys, group in scores.groupby(
        ["model", "dataset", "signature"],
        sort=False,
        dropna=False,
    ):
        control_mask = group["treatment"].eq(control)
        treated_mask = group["treatment"].eq(treated) & np.isclose(group["dose"], treated_dose)
        control_scores = group.loc[control_mask, "score"].to_numpy(float)
        treated_scores = group.loc[treated_mask, "score"].to_numpy(float)
        if len(control_scores) < 2 or len(treated_scores) < 2:
            continue
        labels = np.concatenate([np.zeros(len(control_scores)), np.ones(len(treated_scores))])
        combined = np.concatenate([control_scores, treated_scores])
        auroc = float(roc_auc_score(labels, combined))
        p_expected = float(mannwhitneyu(treated_scores, control_scores, alternative="less").pvalue)
        p_two_sided = float(mannwhitneyu(treated_scores, control_scores, alternative="two-sided").pvalue)
        delta = float(treated_scores.mean() - control_scores.mean())
        rows.append(
            {
                "model": keys[0],
                "dataset": keys[1],
                "signature": keys[2],
                "n_control": len(control_scores),
                "n_treated": len(treated_scores),
                "treated_dose": treated_dose,
                "mean_control": float(control_scores.mean()),
                "mean_treated": float(treated_scores.mean()),
                "delta_mean": delta,
                "auroc": auroc,
                "cliffs_delta": 2 * auroc - 1,
                "expected_direction": "treated_lower",
                "expected_direction_hit": delta < 0,
                "expected_one_tailed_p": p_expected,
                "expected_neglog10_p": -np.log10(max(p_expected, np.finfo(float).tiny)),
                "mw_p_two_sided": p_two_sided,
                "axis_orientation": "continuous_higher_signature_score",
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a172", type=Path, required=True)
    parser.add_argument("--u87", type=Path, required=True)
    parser.add_argument("--gmt", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--model-label", required=True)
    parser.add_argument("--control", default="vehicle")
    parser.add_argument("--treated", default="zstk474")
    parser.add_argument("--treated-dose", type=float, default=10.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    signatures = read_gmt(args.gmt)
    a172 = ad.read_h5ad(args.a172)
    u87 = ad.read_h5ad(args.u87)

    a172_axes, a172_axis_info = build_continuous_axes(a172, signatures, "A172")
    u87_axes, u87_axis_info = build_continuous_axes(u87, signatures, "U87")
    all_scores = pd.concat(
        [
            score_cells(
                "A172", a172, a172_axes, args.model_label
            ),
            score_cells(
                "U87", u87, u87_axes, args.model_label
            ),
        ],
        ignore_index=True,
    )
    stats = summarize(all_scores, args.control.lower(), args.treated.lower(), args.treated_dose)
    axis_info = pd.concat([a172_axis_info, u87_axis_info], ignore_index=True)
    axis_info.insert(0, "model", args.model_label)

    prefix = args.model_label.lower().replace("-", "_").replace(" ", "_")
    all_scores.to_csv(args.outdir / f"{prefix}_cell_scores.csv", index=False)
    stats.to_csv(args.outdir / f"{prefix}_cell_stats.csv", index=False)
    axis_info.to_csv(args.outdir / f"{prefix}_axis_info.csv", index=False)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
