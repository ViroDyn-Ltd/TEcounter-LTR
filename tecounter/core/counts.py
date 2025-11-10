"""
MIT License

Abundance calculations for LTR families and subfamilies.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import pandas as pd

from ..io.gff import GFFRecord


def primary_hits(normalized: pd.DataFrame) -> pd.DataFrame:
    """Select the top-scoring hit per seq_id."""

    if normalized.empty:
        return normalized.copy()
    work = normalized.copy()
    work["score_cov"] = work["score"].astype(float) * work["cov"].astype(float)
    work = work.sort_values("score_cov", ascending=False)
    dedup = work.drop_duplicates("seq_id", keep="first")
    keep_cols = [
        "seq_id",
        "sample",
        "species",
        "length",
        "score",
        "cov",
        "strand",
        "order",
        "superfamily",
        "family_raw",
    ]
    for col in keep_cols:
        if col not in dedup.columns:
            dedup[col] = ""
    return dedup[keep_cols]


def family_summary(
    curated: pd.DataFrame,
    primary: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarize copy counts and statistics per curated family.
    """

    if curated.empty:
        return pd.DataFrame(
            columns=[
                "family_curated",
                "n_copies",
                "total_bp",
                "median_len",
                "truncated_pct",
                "chimera_suspect_pct",
                "mean_score",
                "mean_cov",
            ]
        )

    keep = curated[curated["decision"].isin(["keep", "split"])].copy()
    if keep.empty:
        return family_summary(curated.head(0), primary, domains)

    merged = (
        keep.merge(primary, on="seq_id", how="left", suffixes=("", "_prim"))
        .fillna({"length": 0, "score": 0, "cov": 0, "truncated": 0, "chimera_suspect": 0})
    )
    grouped = merged.groupby("family_curated", sort=False)
    rows = []
    for family, group in grouped:
        lengths = group["length"].astype(float)
        rows.append(
            {
                "family_curated": family,
                "n_copies": int(len(group)),
                "total_bp": int(lengths.sum()),
                "median_len": float(lengths.median()) if len(lengths) else 0.0,
                "truncated_pct": round(group["truncated"].mean(), 3),
                "chimera_suspect_pct": round(group["chimera_suspect"].mean(), 3),
                "mean_score": round(group["score"].mean(), 3),
                "mean_cov": round(group["cov"].mean(), 3),
            }
        )
    return pd.DataFrame(rows)


def subfamily_summary(
    assignments: Dict[str, str],
    representatives: Dict[str, str],
    primary: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summaries per subfamily using membership assignments.
    """

    rows = []
    if not assignments:
        return pd.DataFrame(
            columns=[
                "family_curated",
                "subfamily_id",
                "n_copies",
                "total_bp",
                "median_len",
                "representative_id",
            ]
        )
    primary_lookup = primary.set_index("seq_id").to_dict(orient="index")
    membership_by_family: Dict[Tuple[str, str], List[str]] = defaultdict(list)
    for seq_id, subfam in assignments.items():
        family = subfam.split(".S")[0]
        membership_by_family[(family, subfam)].append(seq_id)

    for (family, subfam), seq_ids in membership_by_family.items():
        lengths = [primary_lookup.get(seq_id, {}).get("length", 0) for seq_id in seq_ids]
        lengths = [float(val) for val in lengths if val is not None]
        rows.append(
            {
                "family_curated": family,
                "subfamily_id": subfam,
                "n_copies": len(seq_ids),
                "total_bp": int(sum(lengths)),
                "median_len": float(pd.Series(lengths).median()) if lengths else 0.0,
                "representative_id": representatives.get(subfam, ""),
            }
        )
    return pd.DataFrame(rows)


def chromosome_windows(
    curated: pd.DataFrame,
    primary: pd.DataFrame,
    gff_index: Dict[str, GFFRecord],
    window: int,
) -> pd.DataFrame:
    """Aggregate counts per chromosome windows when coordinates are available."""

    if not gff_index:
        return pd.DataFrame(
            columns=["sample", "chr", "window_start", "window_end", "n_copies", "bp_covered", "families_csv"]
        )

    merged = curated.merge(primary[["seq_id", "sample", "length"]], on="seq_id", how="left")
    records: Dict[Tuple[str, str, int], Dict[str, object]] = {}

    for _, row in merged.iterrows():
        seq_id = row["seq_id"]
        record = gff_index.get(seq_id)
        if not record:
            continue
        chrom = record.seqid
        start = int(record.start)
        window_start = ((start - 1) // window) * window + 1
        key = (row["sample"], chrom, window_start)
        bucket = records.setdefault(
            key,
            {
                "sample": row["sample"],
                "chr": chrom,
                "window_start": window_start,
                "window_end": window_start + window - 1,
                "n_copies": 0,
                "bp_covered": 0,
                "families": Counter(),
            },
        )
        bucket["n_copies"] += 1
        bucket["bp_covered"] += int(row.get("length") or (record.end - record.start + 1))
        bucket["families"][row["family_curated"]] += 1

    rows = []
    for bucket in records.values():
        families_csv = ",".join(f"{name}:{count}" for name, count in bucket["families"].most_common())
        rows.append(
            {
                "sample": bucket["sample"],
                "chr": bucket["chr"],
                "window_start": bucket["window_start"],
                "window_end": bucket["window_end"],
                "n_copies": bucket["n_copies"],
                "bp_covered": bucket["bp_covered"],
                "families_csv": families_csv,
            }
        )
    return pd.DataFrame(rows)


__all__ = ["primary_hits", "family_summary", "subfamily_summary", "chromosome_windows"]
