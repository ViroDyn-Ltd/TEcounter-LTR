"""
MIT License

Sample manifest handling.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd

REQUIRED_COLUMNS = ["sample_id", "species"]


def read_manifest(path: str | Path) -> pd.DataFrame:
    manifest_path = Path(path)
    df = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Manifest missing columns: {', '.join(missing)}")
    return df


__all__ = ["read_manifest", "REQUIRED_COLUMNS"]
