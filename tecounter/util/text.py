"""
MIT License

Text helpers for TEcounter.
"""

from __future__ import annotations

import re
import textwrap
from typing import Iterable

_SLUG_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def slugify(value: str, default: str = "value") -> str:
    """Generate filesystem-safe slug."""
    cleaned = _SLUG_RE.sub("_", value.strip()) or default
    return re.sub(r"_+", "_", cleaned)


def wrap_block(text: str, width: int = 100, indent: str = "") -> str:
    """Wrap a block of text to the given width with optional indent."""
    return "\n".join(textwrap.fill(line, width=width, subsequent_indent=indent) for line in text.splitlines())


def comma_join(items: Iterable[str]) -> str:
    """Join iterable entries using commas while skipping empties."""
    filtered = [item for item in items if item]
    return ",".join(filtered)


__all__ = ["slugify", "wrap_block", "comma_join"]
