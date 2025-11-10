"""
MIT License

Deterministic LTR subfamily clustering via k-mer features.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

from ..io.fasta import FastaResource


@dataclass
class SubfamilyAssignment:
    family: str
    subfamily_id: str
    members: List[str]
    method: str
    stability: float


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def _canonical_sequence(seq: str) -> str:
    rc = _reverse_complement(seq)
    return seq if seq <= rc else rc


def _kmer_vector(seq: str, k: int = 6) -> np.ndarray:
    alphabet = {"A": 0, "C": 1, "G": 2, "T": 3}
    size = 4 ** k
    arr = np.zeros(size, dtype=np.float32)
    seq = seq.upper()
    for idx in range(len(seq) - k + 1):
        val = 0
        skip = False
        for base in seq[idx : idx + k]:
            if base not in alphabet:
                skip = True
                break
            val = (val << 2) | alphabet[base]
        if skip:
            continue
        arr[val] += 1.0
    total = arr.sum()
    if total > 0:
        arr /= total
    return arr


def _single_link_partition(features: np.ndarray, threshold: float = 0.35) -> np.ndarray:
    if len(features) <= 1:
        return np.zeros(len(features), dtype=int)
    model = AgglomerativeClustering(
        n_clusters=None,
        linkage="single",
        metric="cosine",
        distance_threshold=threshold,
        compute_full_tree=True,
    )
    try:
        labels = model.fit_predict(features)
    except ValueError:
        labels = np.zeros(len(features), dtype=int)
    return labels


def _bootstrap_stability(
    features: np.ndarray,
    labels: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 50,
    sample_frac: float = 0.7,
) -> Dict[int, float]:
    n_features = features.shape[1]
    unique_labels = np.unique(labels)
    if n_features == 0 or len(unique_labels) == 0:
        return {label: 1.0 for label in unique_labels}
    stability = {int(label): 0.0 for label in unique_labels}
    indices = np.arange(features.shape[0])
    for _ in range(n_boot):
        k = max(1, int(n_features * sample_frac))
        cols = rng.choice(n_features, size=k, replace=False)
        subset = features[:, cols]
        boot_labels = _single_link_partition(subset)
        for label in unique_labels:
            orig_members = indices[labels == label]
            if orig_members.size == 0:
                continue
            best = 0.0
            for boot_label in np.unique(boot_labels):
                boot_members = indices[boot_labels == boot_label]
                if boot_members.size == 0:
                    continue
                intersection = np.intersect1d(orig_members, boot_members, assume_unique=True).size
                union = np.union1d(orig_members, boot_members).size
                if union == 0:
                    continue
                best = max(best, intersection / union)
            stability[int(label)] += best
    return {label: round(total / n_boot, 3) for label, total in stability.items()}


def cluster_families(
    curated: pd.DataFrame,
    fasta: FastaResource | None,
    k: int = 6,
    min_cluster_size: int = 15,
    method: str = "kmer",
) -> Tuple[List[SubfamilyAssignment], Dict[str, str]]:
    """
    Cluster LTR copies into subfamilies within each curated family.
    """

    if method == "none" or fasta is None or not fasta.available():
        return [], {}

    assignments: List[SubfamilyAssignment] = []
    membership: Dict[str, str] = {}
    rng = np.random.default_rng(1337)

    grouped: Dict[str, List[str]] = {}
    for entry in curated.to_dict(orient="records"):
        seq_id = entry["seq_id"]
        family = entry["family_curated"]
        decision = entry.get("decision", "keep")
        if decision not in {"keep", "split"}:
            continue
        if not family:
            continue
        grouped.setdefault(family, []).append(seq_id)

    for family, members in grouped.items():
        seq_vectors = []
        valid_ids = []
        for seq_id in members:
            sequence = fasta.get_sequence(seq_id)
            if not sequence:
                continue
            canonical = _canonical_sequence(sequence)
            seq_vectors.append(_kmer_vector(canonical, k=k))
            valid_ids.append(seq_id)
        if not valid_ids:
            continue
        if len(valid_ids) < min_cluster_size:
            subfamily_id = f"{family}.S1"
            assignments.append(
                SubfamilyAssignment(
                    family=family,
                    subfamily_id=subfamily_id,
                    members=valid_ids,
                    method=method,
                    stability=1.0,
                )
            )
            for seq_id in valid_ids:
                membership[seq_id] = subfamily_id
            continue

        features = np.vstack(seq_vectors)
        labels = _single_link_partition(features)
        stability = _bootstrap_stability(features, labels, rng)
        label_to_members: Dict[int, List[str]] = {}
        for label, seq_id in zip(labels, valid_ids):
            label_to_members.setdefault(int(label), []).append(seq_id)

        cluster_idx = 1
        for label, seq_ids in sorted(label_to_members.items(), key=lambda item: (-len(item[1]), item[0])):
            subfamily_id = f"{family}.S{cluster_idx}"
            cluster_idx += 1
            score = stability.get(label, 0.0)
            if score < 0.7:
                # Unstable clusters fall back to S0
                for seq_id in seq_ids:
                    membership[seq_id] = f"{family}.S0"
                continue
            assignments.append(
                SubfamilyAssignment(
                    family=family,
                    subfamily_id=subfamily_id,
                    members=seq_ids,
                    method=method,
                    stability=score,
                )
            )
            for seq_id in seq_ids:
                membership[seq_id] = subfamily_id

        # Ensure fallback cluster exists if some members unassigned
        for seq_id in valid_ids:
            membership.setdefault(seq_id, f"{family}.S0")

    return assignments, membership


__all__ = ["SubfamilyAssignment", "cluster_families"]
