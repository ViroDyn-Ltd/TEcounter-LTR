"""
MIT License

Minimal GFF3 parsing utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional


@dataclass
class GFFRecord:
    """Representation of a single GFF3 record."""

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: Dict[str, str]


def parse_attributes(field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for chunk in field.split(";"):
        if not chunk:
            continue
        if "=" in chunk:
            key, value = chunk.split("=", 1)
        else:
            key, value = chunk, ""
        out[key.strip()] = value.strip()
    return out


def read_gff3(path: str | Path) -> List[GFFRecord]:
    """Load a GFF3 file into memory."""
    records: List[GFFRecord] = []
    gff_path = Path(path)
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            records.append(
                GFFRecord(
                    seqid=seqid,
                    source=source,
                    type=ftype,
                    start=int(start),
                    end=int(end),
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes=parse_attributes(attrs),
                )
            )
    return records


def build_attribute_index(records: Iterable[GFFRecord], attribute: str = "ID") -> Dict[str, GFFRecord]:
    """Map attribute value -> record to enable coordinate lookups."""
    mapping: Dict[str, GFFRecord] = {}
    for record in records:
        attr_val = record.attributes.get(attribute)
        if attr_val:
            mapping[attr_val] = record
    return mapping


__all__ = ["GFFRecord", "read_gff3", "build_attribute_index"]
