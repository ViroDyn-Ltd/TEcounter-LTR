"""
MIT License

TSV/CSV/JSONL helpers for TEcounter.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import json

import pandas as pd


def write_table(df: pd.DataFrame, path: str | Path, fmt: str = "tsv") -> None:
    """
    Persist a DataFrame in the requested serialization format.

    Parameters
    ----------
    df:
        DataFrame to serialize.
    path:
        Output file path.
    fmt:
        One of ``tsv``, ``csv`` or ``jsonl``.
    """
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "tsv":
        df.to_csv(out_path, sep="\t", index=False)
    elif fmt == "csv":
        df.to_csv(out_path, sep=",", index=False)
    elif fmt == "jsonl":
        with out_path.open("w", encoding="utf-8") as handle:
            for record in df.to_dict(orient="records"):
                handle.write(json.dumps(record) + "\n")
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def write_lines(lines: Iterable[str], path: str | Path) -> None:
    """Write plain-text lines joined by newlines."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


__all__ = ["write_table", "write_lines"]
