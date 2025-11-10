"""
MIT License

Representative sequence selection for families and subfamilies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import pandas as pd

from ..io.fasta import FastaResource


@dataclass
class Representative:
    name: str
    seq_id: str
    basis: str
    selection_reason: str
    sequence: str
    sample: str


def _candidate_score(
    seq_id: str,
    curated_lookup: Dict[str, Dict[str, object]],
    length_lookup: Dict[str, float],
) -> Tuple[int, float, float]:
    entry = curated_lookup.get(seq_id, {})
    intact = entry.get("intact_core", 0)
    chimera = entry.get("chimera_suspect", 0)
    confidence = entry.get("confidence", 0.0)
    length = length_lookup.get(seq_id, 0.0)
    # Reward intact, penalize chimera, prefer long sequences
    return (
        int(intact),
        float(length),
        float(confidence) - (0.5 if chimera else 0.0),
    )


def select_representatives(
    curated: pd.DataFrame,
    primary: pd.DataFrame,
    assignments: Dict[str, str],
    fasta: FastaResource | None,
    sample_name: str,
) -> Tuple[List[Representative], Dict[str, str]]:
    """
    Pick representative sequences for each family and subfamily.
    """

    if fasta is None or not fasta.available():
        return [], {}

    curated_lookup = curated.set_index("seq_id").to_dict(orient="index")
    length_lookup = primary.set_index("seq_id")["length"].to_dict()

    representatives: List[Representative] = []
    subfamily_members: Dict[str, List[str]] = {}

    families = curated[curated["decision"].isin(["keep", "split"])].groupby("family_curated")
    for family, group in families:
        candidates = group["seq_id"].tolist()
        if not candidates:
            continue
        best = max(candidates, key=lambda seq_id: _candidate_score(seq_id, curated_lookup, length_lookup))
        seq = fasta.get_sequence(best)
        if not seq:
            continue
        representatives.append(
            Representative(
                name=family,
                seq_id=best,
                basis="family",
                selection_reason="longest_intact",
                sequence=seq,
                sample=sample_name,
            )
        )

    # Subfamily representatives (excluding S0 fallback)
    for seq_id, subfam in assignments.items():
        if subfam.endswith(".S0"):
            continue
        subfamily_members.setdefault(subfam, []).append(seq_id)

    subfamily_index: Dict[str, str] = {}
    for subfam, seq_ids in subfamily_members.items():
        best = max(seq_ids, key=lambda seq_id: _candidate_score(seq_id, curated_lookup, length_lookup))
        seq = fasta.get_sequence(best)
        if not seq:
            continue
        representatives.append(
            Representative(
                name=subfam,
                seq_id=best,
                basis="subfamily",
                selection_reason="cluster_centroid",
                sequence=seq,
                sample=sample_name,
            )
        )
        subfamily_index[subfam] = best

    return representatives, subfamily_index


__all__ = ["Representative", "select_representatives"]
