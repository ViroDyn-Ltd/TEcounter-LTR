"""
MIT License

Domain summarization for LTR retrotransposons.
"""

from __future__ import annotations

import re
from typing import Iterable, List

import pandas as pd

ESSENTIAL_CORE = ["gag", "PR", "RT", "RH", "INT"]
OPTIONAL_CORE = {"PR"}

DOMAIN_KEYWORDS = [
    ("gag", "gag"),
    ("capsid", "gag"),
    ("matrix", "gag"),
    ("protease", "PR"),
    ("prot", "PR"),
    ("AP", "PR"),
    ("reverse_transcriptase", "RT"),
    ("rvt", "RT"),
    ("RT", "RT"),
    ("RNaseH", "RH"),
    ("rh", "RH"),
    ("integrase", "INT"),
    ("int", "INT"),
]


def _canonical_domain(value: str) -> str | None:
    token = (value or "").lower()
    for key, mapped in DOMAIN_KEYWORDS:
        if key.lower() in token:
            return mapped
    return None


def _ordered_subset(domains: List[str], targets: Iterable[str]) -> bool:
    start = 0
    for target in targets:
        found = False
        for idx in range(start, len(domains)):
            if domains[idx] == target:
                found = True
                start = idx + 1
                break
        if not found:
            return False
    return True


def _has_core(domains: List[str]) -> bool:
    if not domains:
        return False
    if _ordered_subset(domains, ESSENTIAL_CORE):
        return True
    reduced = [token for token in ESSENTIAL_CORE if token not in OPTIONAL_CORE]
    return _ordered_subset(domains, reduced)


def call_ltr_domains(df: pd.DataFrame, min_cov: float = 0.6) -> pd.DataFrame:
    """
    Summarize ordered domains per TE copy from normalized rows.
    """

    records = []
    grouped = df.groupby("seq_id", sort=False)
    splitter = re.compile(r"[,\s|]+")
    for seq_id, group in grouped:
        domains: List[str] = []
        for hmm in group["hmm_model"]:
            tokens = [token for token in splitter.split(str(hmm)) if token]
            for token in tokens:
                domain = _canonical_domain(token)
                if domain:
                    if not domains or domains[-1] != domain:
                        domains.append(domain)
        cov = float(group["cov"].max()) if "cov" in group.columns else 0.0
        domains_csv = ",".join(domains)
        intact = 1 if _has_core(domains) else 0
        truncated = 1 if (cov < min_cov or not intact) else 0
        records.append(
            {
                "seq_id": seq_id,
                "domains_csv": domains_csv,
                "intact_core": intact,
                "truncated": truncated,
                "chimera_suspect": 0,
            }
        )
    domain_df = pd.DataFrame.from_records(
        records,
        columns=["seq_id", "domains_csv", "intact_core", "truncated", "chimera_suspect"],
    )
    return domain_df


__all__ = ["call_ltr_domains", "ESSENTIAL_CORE"]
