#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a172-h5ad", type=Path, required=True)
    parser.add_argument("--u87-h5ad", type=Path, required=True)
    parser.add_argument("--cell-scores", type=Path, required=True)
    parser.add_argument("--gene-scores", type=Path)
    parser.add_argument("--treated-dose", type=float, default=10.0)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def read_metadata(path: Path, dataset: str) -> pd.DataFrame:
    obj = ad.read_h5ad(path, backed="r")
    required = {"treatment", "dose"}
    missing = required.difference(obj.obs.columns)
    if missing:
        raise KeyError(f"{path} is missing metadata columns: {sorted(missing)}")
    metadata = obj.obs.loc[:, ["treatment", "dose"]].copy()
    metadata.insert(0, "cell_id", obj.obs_names.astype(str))
    metadata.insert(0, "dataset", dataset)
    metadata["treatment"] = metadata["treatment"].astype(str).str.lower()
    metadata["dose"] = pd.to_numeric(metadata["dose"], errors="coerce")
    obj.file.close()
    if metadata.duplicated(["dataset", "cell_id"]).any():
        raise ValueError(f"Duplicate cell IDs found in {path}")
    return metadata


def read_scores(
    path: Path,
    model: str,
    score_column: str,
    orientation_multiplier: float,
    method: str | None = None,
) -> pd.DataFrame:
    scores = pd.read_csv(path)
    if method is not None:
        if "method" not in scores.columns:
            raise KeyError(f"{path} has no method column")
        scores = scores.loc[scores["method"] == method].copy()
    required = {"dataset", "cell_id", "signature", score_column}
    missing = required.difference(scores.columns)
    if missing:
        raise KeyError(f"{path} is missing columns: {sorted(missing)}")
    scores = scores.loc[:, ["dataset", "cell_id", "signature", score_column]].copy()
    scores["dataset"] = scores["dataset"].astype(str)
    scores["cell_id"] = scores["cell_id"].astype(str)
    scores["model"] = model
    scores["activity_score"] = pd.to_numeric(scores[score_column], errors="raise") * orientation_multiplier
    return scores.drop(columns=score_column)


def verify_direct_duplicate(cell_scores: pd.DataFrame, gene_score_path: Path) -> None:
    gene_scores = read_scores(
        gene_score_path,
        model="Geneformer direct gene cosine",
        score_column="score",
        orientation_multiplier=-1.0,
        method="cosine",
    )
    keys = ["dataset", "cell_id", "signature"]
    merged = cell_scores.merge(gene_scores, on=keys, suffixes=("_cell", "_gene"), validate="one_to_one")
    difference = np.abs(merged["activity_score_cell"] - merged["activity_score_gene"])
    if len(merged) != len(cell_scores) or float(difference.max()) > 1e-7:
        raise ValueError("Direct cell and direct gene cosine scores are not identical")


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (model, dataset, signature), group in scores.groupby(
        ["model", "dataset", "signature"], sort=False
    ):
        control = group.loc[group["group"] == "vehicle", "activity_score"].to_numpy(float)
        treated = group.loc[group["group"] == "zstk474_10uM", "activity_score"].to_numpy(float)
        if len(control) < 2 or len(treated) < 2:
            raise ValueError(f"Insufficient cells for {model}, {dataset}, {signature}")
        less_test = mannwhitneyu(treated, control, alternative="less")
        two_sided_test = mannwhitneyu(treated, control, alternative="two-sided")
        auroc = float(
            (treated[:, None] > control[None, :]).mean()
            + 0.5 * (treated[:, None] == control[None, :]).mean()
        )
        p_expected = max(float(less_test.pvalue), np.finfo(float).tiny)
        rows.append(
            {
                "model": model,
                "dataset": dataset,
                "signature": signature,
                "n_control": len(control),
                "n_treated": len(treated),
                "treated_dose_uM": float(group.loc[group["group"] == "zstk474_10uM", "dose"].iloc[0]),
                "mean_control": float(control.mean()),
                "mean_treated": float(treated.mean()),
                "delta_mean": float(treated.mean() - control.mean()),
                "auroc": auroc,
                "cliffs_delta": 2.0 * auroc - 1.0,
                "expected_direction": "treated_lower",
                "expected_direction_hit": bool(auroc < 0.5),
                "expected_one_tailed_p": p_expected,
                "expected_neglog10_p": float(-np.log10(p_expected)),
                "mw_p_two_sided": float(two_sided_test.pvalue),
                "score_orientation": "AKT1_KO_down_aligned",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    metadata = pd.concat(
        [
            read_metadata(args.a172_h5ad, "A172"),
            read_metadata(args.u87_h5ad, "U87"),
        ],
        ignore_index=True,
    )

    scores = read_scores(
        args.cell_scores,
        model="Geneformer cell embedding",
        score_column="cosine_score",
        orientation_multiplier=-1.0,
    )
    if args.gene_scores is not None:
        verify_direct_duplicate(scores, args.gene_scores)
    scores = scores.merge(metadata, on=["dataset", "cell_id"], how="left", validate="many_to_one")
    if scores[["treatment", "dose"]].isna().any().any():
        raise ValueError("Some score rows could not be matched to H5AD metadata")

    control = (scores["treatment"] == "vehicle") & (scores["dose"] == 0)
    treated = (scores["treatment"] == "zstk474") & (scores["dose"] == args.treated_dose)
    scores = scores.loc[control | treated].copy()
    scores["group"] = np.where(control.loc[scores.index], "vehicle", "zstk474_10uM")

    stats = summarize(scores)
    args.outdir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.outdir / "geneformer_cell_10uM_scores.csv", index=False)
    stats.to_csv(args.outdir / "geneformer_cell_10uM_stats.csv", index=False)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
