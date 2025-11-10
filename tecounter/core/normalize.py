"""
MIT License

Normalization of TEsorter TSV outputs into TEcounter canonical schema.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from ..util.logging import get_logger

LOGGER = get_logger()


@dataclass
class UnknownConfig:
    enabled: bool
    policy: str
    patterns: List[str]
    min_evidence: int


@dataclass
class NormalizationResult:
    table: pd.DataFrame
    total_rows: int
    non_ltr_filtered: int


CANONICAL_COLUMNS = [
    "seq_id",
    "sample",
    "species",
    "order",
    "superfamily",
    "family_raw",
    "hmm_model",
    "score",
    "evalue",
    "cov",
    "qstart",
    "qend",
    "strand",
    "length",
    "source_file",
    "note",
    "unknown_flag",
]

STRING_COLUMNS = {
    "seq_id",
    "sample",
    "species",
    "order",
    "superfamily",
    "family_raw",
    "hmm_model",
    "evalue",
    "strand",
    "source_file",
    "note",
}

COLUMN_ALIASES: Dict[str, Tuple[str, ...]] = {
    "seq_id": ("seq_id", "SequenceID", "seq", "name", "element", "#TE"),
    "order": ("order", "classification", "class", "order_name", "Order"),
    "superfamily": ("superfamily", "subclass", "lineage", "Superfamily"),
    "family_raw": ("family", "family_raw", "clade", "familyName", "Clade"),
    "hmm_model": ("hmm", "model", "hmm_model", "Domains"),
    "score": ("score", "bit_score", "bitscore"),
    "evalue": ("evalue", "expect", "E-value"),
    "cov": ("cov", "coverage", "complete"),
    "qstart": ("qstart", "query_start", "start"),
    "qend": ("qend", "query_end", "end"),
    "strand": ("strand", "orientation"),
    "length": ("length", "qlen", "query_length"),
    "note": ("note",),
}

LTR_KEYWORDS = ("ltr", "ty1", "ty3", "copia", "gypsy", "bel", "pao")


def make_unknown_config(
    enabled: bool = True,
    policy: str = "drop",
    patterns: str | None = None,
    min_evidence: int = 1,
) -> UnknownConfig:
    pattern_list = (
        [pat.strip().lower() for pat in patterns.split(",") if pat.strip()]
        if patterns
        else ["unknown", "unk", "unclassified", "na", "?", "."]
    )
    return UnknownConfig(
        enabled=enabled,
        policy=policy,
        patterns=pattern_list,
        min_evidence=min_evidence,
    )


def _extract_taxonomy(value: str) -> Tuple[str, str, str]:
    order = ""
    superfamily = ""
    family = ""
    tokens = (value or "").replace("::", "/").split("/")
    if tokens:
        order = tokens[0].strip()
    if len(tokens) > 1:
        superfamily = tokens[1].strip()
    if len(tokens) > 2:
        family = tokens[2].strip()
    return order, superfamily, family


def _parse_order(row: pd.Series) -> Tuple[str, str]:
    order = str(row.get("order", "")).strip()
    superfamily = str(row.get("superfamily", "")).strip()
    if not order and row.get("classification"):
        order, superfamily, _ = _extract_taxonomy(str(row["classification"]))
    if not order and superfamily:
        if any(key in superfamily.lower() for key in LTR_KEYWORDS):
            order = "LTR"
    return order or "", superfamily or ""


def normalize_tesorter(
    path: str | Path,
    sample: str,
    species: str,
    ltrs_only: bool,
) -> NormalizationResult:
    """
    Normalize a TEsorter TSV file into canonical TEcounter schema.
    """

    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(f"TEsorter TSV not found: {path}")
    raw = pd.read_csv(input_path, sep="\t", dtype=str, low_memory=False).fillna("")
    total_rows = len(raw)
    if total_rows == 0:
        df = pd.DataFrame(columns=CANONICAL_COLUMNS)
        return NormalizationResult(df, 0, 0)

    rename_map: Dict[str, str] = {}
    lower_to_original = {col.lower(): col for col in raw.columns}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            original = lower_to_original.get(alias.lower())
            if original:
                rename_map[original] = canonical
                break

    df = raw.rename(columns=rename_map).copy()
    df["sample"] = sample
    df["species"] = species

    if "cov" in df.columns:
        df["cov"] = (
            df["cov"]
            .astype(str)
            .str.lower()
            .replace(
                {
                    "yes": "1.0",
                    "true": "1.0",
                    "complete": "1.0",
                    "no": "0.5",
                    "false": "0.0",
                    "partial": "0.5",
                }
            )
        )

    for column in ("score", "cov", "qstart", "qend", "length"):
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0)

    defaults = {
        "order": "",
        "superfamily": "",
        "family_raw": "",
        "hmm_model": "",
        "strand": "+",
        "length": 0,
        "qstart": 0,
        "qend": 0,
    }
    for col, value in defaults.items():
        if col not in df.columns:
            df[col] = value

    def derive_taxonomy(row: pd.Series) -> pd.Series:
        order, superfamily = _parse_order(row)
        family_raw = row.get("family_raw") or ""
        return pd.Series({"order": order, "superfamily": superfamily, "family_raw": family_raw})

    tax_df = df.apply(derive_taxonomy, axis=1)
    df["order"] = tax_df["order"]
    df["superfamily"] = tax_df["superfamily"]
    df["family_raw"] = tax_df["family_raw"]

    ltr_mask = df["order"].str.upper().str.contains("LTR")
    alt_ltr = df["superfamily"].str.lower().apply(lambda x: any(k in x for k in LTR_KEYWORDS))
    ltr_mask |= alt_ltr

    non_ltr_filtered = int((~ltr_mask).sum()) if ltrs_only else 0
    if ltrs_only:
        df = df[ltr_mask].copy()

    df["unknown_flag"] = 0
    df["note"] = ""
    computed_length = (df["qend"] - df["qstart"]).abs() + 1
    if "length" in df.columns:
        mask = (df["length"] <= 0) | df["length"].isna()
        df.loc[mask, "length"] = computed_length.loc[mask]
    else:
        df["length"] = computed_length
    df["length"] = df["length"].fillna(0).astype(int)

    missing_cols = [col for col in CANONICAL_COLUMNS if col not in df.columns]
    for col in missing_cols:
        df[col] = "" if col in STRING_COLUMNS else 0

    df = df[CANONICAL_COLUMNS]
    df["source_file"] = str(input_path)

    LOGGER.info("Normalized %s rows (%s non-LTR filtered)", len(df), non_ltr_filtered)
    return NormalizationResult(df, total_rows, non_ltr_filtered)


__all__ = [
    "UnknownConfig",
    "NormalizationResult",
    "normalize_tesorter",
    "make_unknown_config",
    "CANONICAL_COLUMNS",
]
