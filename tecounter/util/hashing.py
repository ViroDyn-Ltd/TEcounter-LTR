"""
MIT License

Deterministic hashing utilities used for representative identifiers.
"""

from __future__ import annotations

import hashlib


def stable_hash(text: str, length: int = 12) -> str:
    """Return a lowercase deterministic hash prefix."""
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return digest[:length]


__all__ = ["stable_hash"]
