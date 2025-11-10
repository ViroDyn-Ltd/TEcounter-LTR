"""
MIT License

FASTA reading utilities for TEcounter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional
import logging

from Bio import SeqIO

LOGGER = logging.getLogger(__name__)


@dataclass
class FastaRecord:
    """Minimal FASTA record snapshot."""

    seq_id: str
    sequence: str


class FastaResource:
    """Provide deterministic, cached access to FASTA sequences."""

    def __init__(self, path: str | Path | None) -> None:
        self.path = Path(path) if path else None
        self._index: Optional[Dict[str, SeqIO.SeqRecord]] = None
        self._records: Optional[Dict[str, SeqIO.SeqRecord]] = None
        self._lookup: Dict[str, str] = {}
        if self.path:
            try:
                self._index = SeqIO.index(str(self.path), "fasta")
            except ValueError:
                LOGGER.warning("FASTA %s has duplicate IDs; falling back to in-memory store", self.path)
                self._index = None
                self._records = {}
                for record in SeqIO.parse(str(self.path), "fasta"):
                    self._records.setdefault(record.id, record)
        self._build_lookup()

    def _build_lookup(self) -> None:
        store = self._index or self._records
        if not store:
            return
        try:
            keys = list(store.keys())
        except AttributeError:
            keys = list(store)
        for key in keys:
            norm = self._normalize_id(key)
            if norm not in self._lookup:
                self._lookup[norm] = key

    @staticmethod
    def _normalize_id(seq_id: str) -> str:
        token = seq_id.split()[0]
        if "#" in token:
            token = token.split("#")[0]
        return token

    def available(self) -> bool:
        return self._index is not None or self._records is not None

    def _store(self):
        return self._index or self._records

    def iter_lengths(self) -> Iterator[tuple[str, int]]:
        store = self._store()
        if not store:
            return
        if self._lookup:
            for alias, key in self._lookup.items():
                seq = store[key]
                yield alias, len(seq.seq)
        else:
            for seq_id in store:
                seq = store[seq_id]
                yield seq_id, len(seq.seq)

    def get_sequence(self, seq_id: str) -> Optional[str]:
        store = self._store()
        if not store:
            return None
        key = seq_id if seq_id in store else self._lookup.get(self._normalize_id(seq_id))
        if not key:
            return None
        record = store[key]
        return str(record.seq).upper()

    def get_length(self, seq_id: str) -> Optional[int]:
        sequence = self.get_sequence(seq_id)
        if sequence is None:
            return None
        return len(sequence)

    def iter_records(self, ids: Optional[Iterable[str]] = None) -> Iterator[FastaRecord]:
        store = self._store()
        if not store:
            return
        if ids is None:
            ids = list(self._lookup.keys()) if self._lookup else list(store.keys())
        for seq_id in ids:
            seq = self.get_sequence(seq_id)
            if seq is None:
                continue
            yield FastaRecord(seq_id=seq_id, sequence=seq)


def longest_orf_length(seq: str) -> int:
    """Return the length of the longest ORF (simplified)."""

    stop_codons = {"TAA", "TAG", "TGA"}
    seq = seq.upper()
    best = 0
    for frame in range(3):
        start = None
        for idx in range(frame, len(seq) - 2, 3):
            codon = seq[idx : idx + 3]
            if codon == "ATG" and start is None:
                start = idx
            elif codon in stop_codons and start is not None:
                best = max(best, idx + 3 - start)
                start = None
        if start is not None:
            best = max(best, len(seq) - start)
    return best


__all__ = ["FastaResource", "FastaRecord", "longest_orf_length"]
