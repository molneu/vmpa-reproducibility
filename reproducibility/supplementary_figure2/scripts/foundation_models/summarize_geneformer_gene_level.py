#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a172-h5ad", type=Path, required=True)
    parser.add_argument("--u87-h5ad", type=Path, required=True)
    parser.add_argument("--gene-scores", type=Path, required=True)
    parser.add_argument("--treated-dose", type=float, default=10.0)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def read_metadata(path: Path, dataset: str) -> pd.DataFrame:
    obj = ad.read_h5ad(path, backed="r")
    metadata = obj.obs.loc[:, ["treatment", "dose"]].copy()
    metadata.insert(0, "cell_id", obj.obs_names.astype(str))
    metadata.insert(0, "dataset", dataset)
    metadata["treatment"] = metadata["treatment"].astype(str).str.lower()
    metadata["dose"] = pd.to_numeric(metadata["dose"], errors="raise")
    obj.file.close()
    return metadata


def read_gene_scores(path: Path) -> pd.DataFrame:
    scores = pd.read_csv(path)
    scores = scores.loc[
        scores["method"].isin(["cosine", "spearman"]),
        ["dataset", "cell_id", "signature", "method", "score"],
    ].copy()
    scores["dataset"] = scores["dataset"].astype(str)
    scores["cell_id"] = scores["cell_id"].astype(str)
    raw_scores = pd.to_numeric(scores["score"], errors="raise")
    scores["activity_score"] = np.where(scores["method"].eq("spearman"), raw_scores, -raw_scores)
    scores["representation"] = "pretrained_gene_embeddings"
    return scores.drop(columns="score")


def summarize(scores: pd.DataFrame, treated_dose: float) -> pd.DataFrame:
    rows = []
    grouping = ["representation", "dataset", "signature", "method"]
    for keys, group in scores.groupby(grouping, sort=False):
        control = group.loc[
            group["treatment"].eq("vehicle") & group["dose"].eq(0),
            "activity_score",
        ].to_numpy(float)
        treated = group.loc[
            group["treatment"].eq("zstk474") & np.isclose(group["dose"], treated_dose),
            "activity_score",
        ].to_numpy(float)
        if len(control) < 2 or len(treated) < 2:
            raise ValueError(f"Insufficient cells for {keys}")
        labels = np.concatenate([np.zeros(len(control)), np.ones(len(treated))])
        values = np.concatenate([control, treated])
        auroc = float(roc_auc_score(labels, values))
        expected_p = max(
            float(mannwhitneyu(treated, control, alternative="less").pvalue),
            np.finfo(float).tiny,
        )
        two_sided_p = float(mannwhitneyu(treated, control, alternative="two-sided").pvalue)
        rows.append(
            {
                "representation": keys[0],
                "dataset": keys[1],
                "signature": keys[2],
                "method": keys[3],
                "n_control": len(control),
                "n_treated": len(treated),
                "treated_dose_uM": treated_dose,
                "mean_control": float(control.mean()),
                "mean_treated": float(treated.mean()),
                "delta_mean": float(treated.mean() - control.mean()),
                "auroc": auroc,
                "cliffs_delta": 2.0 * auroc - 1.0,
                "expected_direction": "treated_lower",
                "expected_direction_hit": bool(auroc < 0.5),
                "expected_one_tailed_p": expected_p,
                "expected_neglog10_p": float(-np.log10(expected_p)),
                "mw_p_two_sided": two_sided_p,
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
    scores = read_gene_scores(args.gene_scores)
    scores = scores.merge(metadata, on=["dataset", "cell_id"], how="left", validate="many_to_one")
    if scores[["treatment", "dose"]].isna().any().any():
        raise ValueError("Some score rows could not be matched to H5AD metadata")
    control = scores["treatment"].eq("vehicle") & scores["dose"].eq(0)
    treated = scores["treatment"].eq("zstk474") & np.isclose(scores["dose"], args.treated_dose)
    scores = scores.loc[control | treated].copy()
    statistics = summarize(scores, args.treated_dose)
    args.outdir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.outdir / "geneformer_gene_level_10uM_scores.csv", index=False)
    statistics.to_csv(args.outdir / "geneformer_gene_level_10uM_stats.csv", index=False)
    print(statistics.to_string(index=False))


if __name__ == "__main__":
    main()
