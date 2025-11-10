"""
MIT License

Cross-sample comparison utilities.
"""

from __future__ import annotations

import itertools
import math
from typing import Dict, Iterable, List, Tuple

import pandas as pd


def family_matrix(curated: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    """Build samples × families matrix with counts and bp percentages."""

    keep = curated[curated["decision"].isin(["keep", "split"])].copy()
    if keep.empty:
        return pd.DataFrame(columns=["family_curated"])
    merged = keep.merge(primary[["seq_id", "sample", "length"]], on="seq_id", how="left")
    families = sorted(merged["family_curated"].unique())
    samples = sorted(merged["sample"].dropna().unique())
    total_bp = merged.groupby("sample")["length"].sum().to_dict()

    matrix = pd.DataFrame({"family_curated": families})
    for sample in samples:
        counts = (
            merged[merged["sample"] == sample]
            .groupby("family_curated")["seq_id"]
            .count()
            .reindex(families)
            .fillna(0)
            .astype(int)
        )
        bp = (
            merged[merged["sample"] == sample]
            .groupby("family_curated")["length"]
            .sum()
            .reindex(families)
            .fillna(0.0)
        )
        bp_pct = bp / (total_bp.get(sample, 1) or 1) * 100.0
        matrix[f"{sample}.count"] = counts.values
        matrix[f"{sample}.bp_pct"] = bp_pct.round(3).values
    return matrix


def subfamily_matrix(assignments: Dict[str, str], primary: pd.DataFrame) -> pd.DataFrame:
    """Samples × subfamilies matrix."""

    if not assignments:
        return pd.DataFrame(columns=["family_subfamily"])
    membership = pd.DataFrame(
        [(seq_id, subfam) for seq_id, subfam in assignments.items()],
        columns=["seq_id", "subfamily"],
    )
    merged = membership.merge(primary[["seq_id", "sample", "length"]], on="seq_id", how="left")
    subfamilies = sorted(merged["subfamily"].unique())
    samples = sorted(merged["sample"].dropna().unique())
    matrix = pd.DataFrame({"family_subfamily": subfamilies})

    totals = merged.groupby("sample")["length"].sum().to_dict()
    for sample in samples:
        sample_df = merged[merged["sample"] == sample]
        counts = (
            sample_df.groupby("subfamily")["seq_id"].count().reindex(subfamilies).fillna(0).astype(int)
        )
        bp = sample_df.groupby("subfamily")["length"].sum().reindex(subfamilies).fillna(0.0)
        bp_pct = bp / (totals.get(sample, 1) or 1) * 100.0
        matrix[f"{sample}.count"] = counts.values
        matrix[f"{sample}.bp_pct"] = bp_pct.round(3).values
    return matrix


def pairwise_comparison(family_mat: pd.DataFrame, presence_threshold: int = 3) -> pd.DataFrame:
    """Pairwise statistics across samples."""

    sample_names = sorted(
        {col.split(".")[0] for col in family_mat.columns if col.endswith(".count")}
    )
    rows = []
    for a, b in itertools.combinations(sample_names, 2):
        counts_a = family_mat[f"{a}.count"]
        counts_b = family_mat[f"{b}.count"]
        presence_a = counts_a >= presence_threshold
        presence_b = counts_b >= presence_threshold
        union = (presence_a | presence_b).sum()
        intersection = (presence_a & presence_b).sum()
        jaccard = float(intersection / union) if union else 0.0
        pair_key = f"{a}~{b}"
        for idx, family in enumerate(family_mat["family_curated"]):
            val_a = counts_a.iloc[idx]
            val_b = counts_b.iloc[idx]
            log2fc = math.log2((val_a + 1) / (val_b + 1))
            diff = abs(val_a - val_b)
            if presence_a.iloc[idx] and presence_b.iloc[idx]:
                presence = "both"
            elif presence_a.iloc[idx]:
                presence = "only_A"
            elif presence_b.iloc[idx]:
                presence = "only_B"
            else:
                continue
            rows.append(
                {
                    "pair": pair_key,
                    "family": family,
                    "log2FC": round(log2fc, 3),
                    "abs_diff": diff,
                    "present_in": presence,
                    "jaccard_over_families": round(jaccard, 3),
                }
            )
    return pd.DataFrame(rows, columns=["pair", "family", "log2FC", "abs_diff", "present_in", "jaccard_over_families"])


__all__ = ["family_matrix", "subfamily_matrix", "pairwise_comparison"]
