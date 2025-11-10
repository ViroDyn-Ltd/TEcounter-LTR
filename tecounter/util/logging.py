"""
MIT License

Lightweight logging helpers for TEcounter.
"""

from __future__ import annotations

import logging
import sys
from typing import Optional

_LOGGER: Optional[logging.Logger] = None


def get_logger(name: str = "tecounter") -> logging.Logger:
    """Return a process-wide logger configured for CLI use."""
    global _LOGGER
    if _LOGGER is None:
        logger = logging.getLogger(name)
        handler = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(
            "%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%dT%H:%M:%S"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
        _LOGGER = logger
    return _LOGGER


__all__ = ["get_logger"]
