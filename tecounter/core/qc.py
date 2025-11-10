"""
MIT License

Conflict resolution, chimera detection, and confidence scoring.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd

from ..io.fasta import FastaResource
from ..util.logging import get_logger

LOGGER = get_logger()


@dataclass
class CurationResult:
    curated: pd.DataFrame
    updated_domains: pd.DataFrame
    unknown_flags: pd.Series
    stats: Dict[str, int]


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _score_z(df: pd.DataFrame) -> pd.Series:
    scores = df["score"].astype(float)
    grouped = df.groupby(["sample", "superfamily"])["score"]
    mean = grouped.transform("mean")
    std = grouped.transform("std").replace(0, 1.0)
    return (scores - mean) / std


def _chimera_by_scores(top_row: pd.Series, second_row: pd.Series | None) -> bool:
    if second_row is None:
        return False
    if top_row["superfamily"] == second_row["superfamily"]:
        return False
    return top_row["score_z"] >= 0.5 and second_row["score_z"] >= 0.5


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def _kmer_change_point(seq: str, k: int = 6, window: int = 1000, stride: int = 500) -> bool:
    seq = seq.upper()
    if len(seq) < window * 2:
        return False
    from sklearn.metrics.pairwise import cosine_distances

    windows = []
    for start in range(0, len(seq) - window + 1, stride):
        fragment = seq[start : start + window]
        windows.append(_kmer_vector(fragment, k))
    if len(windows) < 2:
        return False
    matrix = np.vstack(windows)
    distances = cosine_distances(matrix[:-1], matrix[1:])
    return bool((distances.diagonal() > 0.35).any())


def _kmer_vector(seq: str, k: int) -> np.ndarray:
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    size = 4 ** k
    arr = np.zeros(size, dtype=np.float32)
    for idx in range(len(seq) - k + 1):
        kmer = seq[idx : idx + k]
        try:
            val = 0
            for char in kmer:
                val = (val << 2) | mapping[char]
        except KeyError:
            continue
        arr[val] += 1
    total = arr.sum()
    if total:
        arr /= total
    return arr


def _is_unknown(
    family_curated: str,
    family_raw: str,
    superfamily: str,
    confidence: float,
    patterns: list[str],
) -> bool:
    fam_lower = (family_raw or "").strip().lower()
    curated_lower = (family_curated or "").strip().lower()
    if not fam_lower:
        fam_lower = ""
    matches_pattern = bool(fam_lower == "" or any(pat in fam_lower for pat in patterns))
    ambiguous_superfamily = (
        "?" in (superfamily or "")
        or "/" in (superfamily or "")
        or (superfamily or "").strip() == ""
    )
    if matches_pattern:
        return True
    if ambiguous_superfamily and confidence < 0.5:
        return True
    if "unresolved" in curated_lower:
        return True
    return False


def resolve_conflicts(
    normalized: pd.DataFrame,
    domains: pd.DataFrame,
    unknown_cfg,
    filter_unknowns: bool,
    strict_chimera: str = "keep",
    fasta: FastaResource | None = None,
) -> CurationResult:
    """
    Resolve per-copy labels, apply chimera detection, and compute confidence.
    """

    curated_records = []
    unknown_flags = {}
    stats = {
        "unknown_discarded": 0,
        "unknown_tagged": 0,
        "unknown_kept": 0,
    }

    work = normalized.copy()
    if work.empty:
        empty = pd.DataFrame(
            columns=["seq_id", "order", "superfamily", "family_curated", "family_raw", "decision", "confidence", "rationale"]
        )
        merged_domains = domains.copy()
        merged_domains["chimera_suspect"] = 0
        return CurationResult(empty, merged_domains, pd.Series(dtype=int), stats)

    work["score_z"] = _score_z(work)
    work["score_cov"] = work["score"] * work["cov"].clip(lower=0.0)
    domain_lookup = domains.set_index("seq_id").to_dict(orient="index")

    for seq_id, group in work.groupby("seq_id"):
        rows = group.sort_values("score_cov", ascending=False)
        top = rows.iloc[0]
        second = rows.iloc[1] if len(rows) > 1 else None
        family_raw = top["family_raw"] or top["superfamily"]
        family_curated = family_raw
        rationale = "top_score"
        decision = "keep"

        if second is not None:
            diff = abs(top["score_cov"] - second["score_cov"]) / (top["score_cov"] + 1e-6)
            if top["superfamily"] == second["superfamily"] and diff < 0.05:
                family_curated = f"{top['superfamily']}:Unresolved"
                rationale = "ambiguous_family"

        chimera = 1 if _chimera_by_scores(top, second) else 0
        if chimera and strict_chimera == "drop":
            decision = "drop"
        elif chimera and strict_chimera == "split":
            decision = "split"

        if fasta and chimera == 0:
            seq = fasta.get_sequence(seq_id)
            if seq and _kmer_change_point(seq):
                chimera = 1
                if strict_chimera == "drop":
                    decision = "drop"
                elif strict_chimera == "split":
                    decision = "split"
                rationale = "kmer_changepoint"

        domain_info = domain_lookup.get(seq_id, {})
        intact = domain_info.get("intact_core", 0)
        truncated = domain_info.get("truncated", 0)
        domain_hits = len(
            [part for part in (domain_info.get("domains_csv") or "").split(",") if part]
        )

        confidence = _sigmoid(
            0.5 * float(top["score_z"])
            + 0.8 * float(top["cov"])
            + 0.4 * intact
            - 0.7 * chimera
        )
        confidence = max(0.0, min(1.0, confidence))

        patterns = unknown_cfg.patterns if unknown_cfg.enabled else []
        unknown = False
        if unknown_cfg.enabled:
            unknown = _is_unknown(
                family_curated,
                family_raw,
                top["superfamily"],
                confidence,
                patterns,
            )

        unknown_flag = 1 if unknown else 0

        min_evidence_trigger = unknown and domain_hits < unknown_cfg.min_evidence
        if min_evidence_trigger:
            decision = "drop"
            stats["unknown_discarded"] += 1
        elif unknown and filter_unknowns:
            if unknown_cfg.policy == "drop":
                stats["unknown_discarded"] += 1
                decision = "drop"
            elif unknown_cfg.policy == "keep":
                stats["unknown_kept"] += 1
            elif unknown_cfg.policy == "tag":
                family_curated = "LTR_Unknown"
                stats["unknown_tagged"] += 1
        elif unknown:
            stats["unknown_kept"] += 1

        curated_records.append(
            {
                "seq_id": seq_id,
                "order": top["order"],
                "superfamily": top["superfamily"],
                "family_curated": family_curated,
                "family_raw": family_raw,
                "decision": decision,
                "confidence": round(confidence, 3),
                "rationale": rationale,
                "intact_core": intact,
                "truncated": truncated,
                "chimera_suspect": chimera,
            }
        )
        unknown_flags[seq_id] = unknown_flag

    curated_df = pd.DataFrame.from_records(
        curated_records,
        columns=[
            "seq_id",
            "order",
            "superfamily",
            "family_curated",
            "family_raw",
            "decision",
            "confidence",
            "rationale",
            "intact_core",
            "truncated",
            "chimera_suspect",
        ],
    )

    merged_domains = domains.merge(
        curated_df[["seq_id", "chimera_suspect"]],
        on="seq_id",
        how="left",
        suffixes=("", "_curated"),
    )
    merged_domains["chimera_suspect"] = merged_domains["chimera_suspect"].fillna(
        merged_domains.get("chimera_suspect_curated", 0)
    )
    merged_domains = merged_domains.drop(columns=[col for col in merged_domains.columns if col.endswith("_curated")])

    return CurationResult(
        curated=curated_df,
        updated_domains=merged_domains,
        unknown_flags=pd.Series(unknown_flags, name="unknown_flag"),
        stats=stats,
    )


__all__ = ["CurationResult", "resolve_conflicts"]
