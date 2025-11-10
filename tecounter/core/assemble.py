"""
MIT License

Assemble the publish-ready text bundle for TEcounter outputs.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

from ..io.tsv import write_lines


def _format_top_families(summary: Dict[str, int]) -> str:
    parts = [f"{family}:{count}" for family, count in summary.items()]
    return " | ".join(parts) if parts else "NA"


def write_all_txt(context: Dict[str, object], out_dir: str | Path) -> None:
    """Render ALL.txt with ordered sections."""

    out_path = Path(out_dir) / "ALL.txt"
    lines: List[str] = []
    lines.append("# TECOUNTER RUN")
    lines.append(f"version: {context['version']}")
    lines.append(f"date_utc: {context['date_utc']}")
    lines.append(f"ltrs_only: {str(context['ltrs_only']).lower()}")
    lines.append(f"unknown_policy: {context['unknown_policy']}")
    lines.append(f"unknown_patterns: {','.join(context['unknown_patterns'])}")
    lines.append(f"non_ltr_filtered: {context['non_ltr_filtered']}")
    lines.append(f"unknown_filtered: {context['unknown_filtered']}")
    lines.append(f"n_samples: {context['n_samples']}")
    lines.append("")

    for sample in context.get("sample_summaries", []):
        lines.append(f"## SAMPLE {sample['sample']} ({sample['species']})")
        lines.append(f"copies_total_LTR: {sample['copies_total']}")
        lines.append(f"bp_total_LTR: {sample['bp_total']}")
        lines.append(f"families_detected_LTR: {sample['families_detected']}")
        lines.append(
            "top_families_by_copy: "
            + _format_top_families(sample.get("top_families", {}))
        )
        lines.append("")

    lines.append("### DOMAIN_SUMMARY (top 10 LTR families by n_copies)")
    for entry in context.get("domain_summary", []):
        lines.append(
            f"family {entry['family']}  domains {entry['domains']}  intact_pct {entry['intact_pct']}  truncated_pct {entry['truncated_pct']}"
        )
    lines.append("")

    lines.append("### SUBFAMILIES")
    for entry in context.get("subfamily_summary", []):
        lines.append(
            f"family {entry['family']}  subfamily {entry['subfamily']}  n={entry['n']}  rep={entry['rep']}"
        )
    lines.append("")

    if context.get("chr_windows"):
        lines.append(f"### CHR WINDOWS ({context['window']} bp)")
        for entry in context["chr_windows"]:
            lines.append(
                f"chr {entry['chr']}  {entry['window_start']}-{entry['window_end']}  n={entry['n_copies']}  bp={entry['bp_covered']}  families {entry['families_csv']}"
            )
        lines.append("")

    lines.append("## CROSS-SAMPLE")
    cross = context.get("cross_sample", {})
    lines.append(
        f"matrix_families rows={cross.get('rows', 0)} cols={cross.get('cols', 0)} filled={cross.get('fill', '0%')}"
    )
    lines.append("top_shifted_families (log2FC max across pairs):")
    for entry in cross.get("top_shifted", []):
        lines.append(f"{entry['family']} {entry['value']} ({entry['pair']})")
    lines.append("pairwise Jaccard (presence/absence):")
    for entry in cross.get("jaccard", []):
        lines.append(f"{entry['pair']}: {entry['value']}")
    lines.append("")

    lines.append("## FILTERED SUMMARY")
    filtered = context.get("filtered_summary", {})
    lines.append(f"non_LTR_discarded: {filtered.get('non_ltr', 0)}")
    lines.append(f"unknown_LTR_discarded: {filtered.get('unknown_drop', 0)}")
    lines.append(f"unknown_LTR_kept_tagged: {filtered.get('unknown_tag', 0)}")

    write_lines(lines, out_path)


def write_run_metadata(path: str | Path, metadata: Dict[str, object]) -> None:
    """Write run.jsonl with a single NDJSON record."""

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(metadata)
    if "date_utc" not in metadata:
        metadata["date_utc"] = datetime.now(timezone.utc).isoformat()
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(metadata) + "\n")


def write_readme(path: str | Path) -> None:
    """Emit README.txt summarizing usage and schemas."""

    out_path = Path(path)
    content = """TEcounter (LTR-specialized)
================================

Quickstart
----------
Single sample:
    python -m tecounter --tesorter input.tsv --sample SAMPLE --species Species --out out_dir

Multi-sample folder mode:
    python -m tecounter folder --in runs/ --manifest manifest.tsv --out out_dir

Key behaviors
-------------
* Processes only LTR retrotransposons unless --ltrs-only false.
* Deterministic k-mer clustering (Agglomerative single-link) within families.
* Unknown filtering controlled via --filter-unknowns, --unknown-policy, and --unknown-patterns.

Schema overview
---------------
01_normalized.tsv      Canonical TEsorter rows (LTR only when requested).
02_domains.tsv         Per-copy domain calls and integrity flags.
03_curated.tsv         Final family assignments with confidence.
04_subfamilies.tsv     Subfamily membership after k-mer clustering.
05_family_summary.tsv  Abundance metrics per family.
06_subfamily_summary.tsv Summary metrics per subfamily.
07_chr_windows.tsv     Windowed chromosome counts (if GFF provided).
08_representatives.*   FASTA + index of representative sequences.
09_family_matrix.tsv   Samples × families counts/bp%.
10_subfamily_matrix.tsv Samples × subfamilies counts/bp%.
11_pairwise_comparison.tsv Pairwise log2FC and Jaccard summaries.
filtered.summary.txt   Totals for filtered non-LTR/unknown rows.
ALL.txt                Human-readable rollup of the run.

This bundle is text-only (TSV/CSV/TXT/FASTA)."""
    write_lines(content.splitlines(), out_path)


__all__ = ["write_all_txt", "write_run_metadata", "write_readme"]
