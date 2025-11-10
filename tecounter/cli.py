"""
MIT License

Command-line interface for TEcounter.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from .core import counts
from .core.assemble import write_all_txt, write_readme, write_run_metadata
from .core.cluster import SubfamilyAssignment
from .core.compare import family_matrix as build_family_matrix, pairwise_comparison, subfamily_matrix as build_subfamily_matrix
from .core.normalize import NormalizationResult, make_unknown_config
from .core.pipeline import SampleConfig, SampleResult, run_sample
from .core.qc import CurationResult
from .core.representatives import Representative
from .io.manifest import read_manifest
from .io.tsv import write_lines, write_table
from .util.logging import get_logger

LOGGER = get_logger()


def _bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean: {value}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="tecounter", description="LTR-specialized TE counter")
    add_single_args(parser)
    subparsers = parser.add_subparsers(dest="command")

    folder_parser = subparsers.add_parser("folder", help="Process multiple samples from a manifest")
    add_folder_args(folder_parser)

    dump_parser = subparsers.add_parser("dump", help="Emit a specific artifact from an output folder")
    dump_parser.set_defaults(handler=run_dump_command)
    dump_parser.add_argument("--out", required=True, help="Run output folder")
    dump_parser.add_argument(
        "--what",
        required=True,
        choices=[
            "normalized",
            "domains",
            "curated",
            "subfamilies",
            "family_summary",
            "subfamily_summary",
            "chr_windows",
            "representatives_index",
            "family_matrix",
            "subfamily_matrix",
            "pairwise",
            "all_txt",
        ],
    )
    return parser


def add_common_processing_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cluster", choices=["none", "kmer"], default="kmer")
    parser.add_argument("--k", type=int, default=6, help="k-mer size")
    parser.add_argument("--min-cov", type=float, default=0.60)
    parser.add_argument("--min-score", type=float, default=80.0)
    parser.add_argument("--strict-chimera", choices=["keep", "drop", "split"], default="keep")
    parser.add_argument("--window", type=int, default=100000)
    parser.add_argument("--emit", choices=["tsv", "csv", "jsonl"], default="tsv")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--ltrs-only", type=_bool, default=True)
    parser.add_argument("--filter-unknowns", type=_bool, default=True)
    parser.add_argument("--unknown-policy", choices=["drop", "keep", "tag"], default="drop")
    parser.add_argument("--unknown-patterns", default="unknown,unk,unclassified,NA,?,.")
    parser.add_argument("--min-evidence", type=int, default=1)


def add_single_args(parser: argparse.ArgumentParser) -> None:
    parser.set_defaults(handler=run_single_command)
    parser.add_argument("--tesorter", required=False, help="TEsorter TSV path")
    parser.add_argument("--fasta", help="FASTA file with TE copies")
    parser.add_argument("--gff", help="GFF3 annotations for coordinates")
    parser.add_argument("--sample", required=False, help="Sample identifier")
    parser.add_argument("--species", required=False, help="Species name")
    parser.add_argument("--out", required=False, default="out", help="Output directory")
    add_common_processing_args(parser)


def add_folder_args(parser: argparse.ArgumentParser) -> None:
    parser.set_defaults(handler=run_folder_command)
    parser.add_argument("--in", dest="input_dir", required=True, help="Input folder containing *.tesorter.tsv")
    parser.add_argument("--manifest", required=True, help="Sample manifest TSV")
    parser.add_argument("--fasta-dir", help="Folder containing per-sample FASTA files")
    parser.add_argument("--gff-dir", help="Folder containing per-sample GFF3 files")
    parser.add_argument("--out", required=True, help="Output directory for combined bundle")
    parser.add_argument("--per-sample-out", action="store_true", help="Also write per-sample bundles")
    add_common_processing_args(parser)


def dispatch(args: argparse.Namespace) -> None:
    random.seed(1337)
    np.random.seed(1337)
    if getattr(args, "handler", None) is None:
        args.handler = run_single_command
    args.handler(args)


def run_single_command(args: argparse.Namespace) -> None:
    if not args.tesorter or not args.sample or not args.species or not args.out:
        raise SystemExit("Single-sample mode requires --tesorter, --sample, --species, and --out")
    unknown_cfg = make_unknown_config(
        enabled=_bool(args.filter_unknowns),
        policy=args.unknown_policy,
        patterns=args.unknown_patterns,
        min_evidence=args.min_evidence,
    )
    sample_cfg = SampleConfig(
        tesorter=args.tesorter,
        sample=args.sample,
        species=args.species,
        fasta=args.fasta,
        gff=args.gff,
        ltrs_only=_bool(args.ltrs_only),
        filter_unknowns=_bool(args.filter_unknowns),
        unknown_cfg=unknown_cfg,
        min_cov=args.min_cov,
        min_score=args.min_score,
        cluster_method=args.cluster,
        k=args.k,
        window=args.window,
        strict_chimera=args.strict_chimera,
    )
    result = run_sample(sample_cfg)
    out_dir = Path(args.out)
    write_bundle(result, out_dir, emit=args.emit, component_results=[result])
    LOGGER.info("Bundle written to %s", out_dir)


def run_folder_command(args: argparse.Namespace) -> None:
    manifest_df = read_manifest(args.manifest)
    input_dir = Path(args.input_dir)
    tesorter_map = {}
    for path in input_dir.glob("*.tesorter.tsv"):
        name = path.name.replace(".tesorter.tsv", "")
        tesorter_map[name] = path
    results: List[SampleResult] = []
    unknown_cfg = make_unknown_config(
        enabled=_bool(args.filter_unknowns),
        policy=args.unknown_policy,
        patterns=args.unknown_patterns,
        min_evidence=args.min_evidence,
    )
    for _, row in manifest_df.iterrows():
        sample = row["sample_id"]
        species = row.get("species", "Unknown")
        tesorter_path = tesorter_map.get(sample)
        if not tesorter_path:
            raise SystemExit(f"Missing TEsorter TSV for sample {sample}")
        fasta_path = None
        if args.fasta_dir:
            fasta_candidate = Path(args.fasta_dir) / f"{sample}.fa"
            if not fasta_candidate.exists():
                fasta_candidate = Path(args.fasta_dir) / f"{sample}.fasta"
            if fasta_candidate.exists():
                fasta_path = str(fasta_candidate)
        gff_path = None
        if args.gff_dir:
            gff_candidate = Path(args.gff_dir) / f"{sample}.gff3"
            if gff_candidate.exists():
                gff_path = str(gff_candidate)
        sample_cfg = SampleConfig(
            tesorter=str(tesorter_path),
            sample=sample,
            species=species,
            fasta=fasta_path,
            gff=gff_path,
            ltrs_only=_bool(args.ltrs_only),
            filter_unknowns=_bool(args.filter_unknowns),
            unknown_cfg=unknown_cfg,
            min_cov=args.min_cov,
            min_score=args.min_score,
            cluster_method=args.cluster,
            k=args.k,
            window=args.window,
            strict_chimera=args.strict_chimera,
        )
        result = run_sample(sample_cfg)
        results.append(result)
        if args.per_sample_out:
            sample_out = Path(args.out) / sample
            write_bundle(result, sample_out, emit=args.emit, component_results=[result])
    if not results:
        raise SystemExit("Manifest produced zero samples")
    combined = combine_results(results, unknown_cfg, args)
    write_bundle(combined, Path(args.out), emit=args.emit, component_results=results)
    LOGGER.info("Combined bundle written to %s", args.out)


def combine_results(
    results: List[SampleResult],
    unknown_cfg,
    args: argparse.Namespace,
) -> SampleResult:
    normalization_table = pd.concat([res.normalization.table for res in results], ignore_index=True)
    total_rows = sum(res.normalization.total_rows for res in results)
    non_ltr_filtered = sum(res.normalization.non_ltr_filtered for res in results)
    normalization = NormalizationResult(normalization_table, total_rows, non_ltr_filtered)

    curated_df = pd.concat([res.curation.curated for res in results], ignore_index=True)
    updated_domains = pd.concat([res.domains for res in results], ignore_index=True)
    unknown_flags = pd.concat([res.curation.unknown_flags for res in results])
    stats = {
        "unknown_discarded": sum(res.curation.stats.get("unknown_discarded", 0) for res in results),
        "unknown_tagged": sum(res.curation.stats.get("unknown_tagged", 0) for res in results),
        "unknown_kept": sum(res.curation.stats.get("unknown_kept", 0) for res in results),
    }
    curation = CurationResult(
        curated=curated_df,
        updated_domains=updated_domains,
        unknown_flags=unknown_flags,
        stats=stats,
    )
    primary = pd.concat([res.primary for res in results], ignore_index=True)

    assignments: List[SubfamilyAssignment] = []
    assignment_map: Dict[str, str] = {}
    representatives: List[Representative] = []
    subfamily_rep_map: Dict[str, str] = {}
    for res in results:
        sample_suffix = res.config.sample
        for assignment in res.assignments:
            subfamily_id = f"{assignment.subfamily_id}@{sample_suffix}"
            assignments.append(
                SubfamilyAssignment(
                    family=assignment.family,
                    subfamily_id=subfamily_id,
                    members=assignment.members,
                    method=assignment.method,
                    stability=assignment.stability,
                )
            )
            for member in assignment.members:
                assignment_map[member] = subfamily_id
        for key, value in res.subfamily_rep_map.items():
            subfamily_rep_map[f"{key}@{sample_suffix}"] = value
        representatives.extend(res.representatives)

    family_summary_df = counts.family_summary(curated_df, primary)
    subfamily_summary_df = counts.subfamily_summary(assignment_map, subfamily_rep_map, primary)
    chr_windows_df = pd.concat([res.chr_windows for res in results], ignore_index=True)

    combined_cfg = SampleConfig(
        tesorter="combined",
        sample="Combined",
        species="multi_sample",
        fasta=None,
        gff=None,
        ltrs_only=_bool(args.ltrs_only),
        filter_unknowns=_bool(args.filter_unknowns),
        unknown_cfg=unknown_cfg,
        min_cov=args.min_cov,
        min_score=args.min_score,
        cluster_method=args.cluster,
        k=args.k,
        window=args.window,
        strict_chimera=args.strict_chimera,
    )
    return SampleResult(
        config=combined_cfg,
        normalization=normalization,
        domains=updated_domains,
        curation=curation,
        primary=primary,
        assignments=assignments,
        assignment_map=assignment_map,
        representatives=representatives,
        subfamily_rep_map=subfamily_rep_map,
        family_summary=family_summary_df,
        subfamily_summary=subfamily_summary_df,
        chr_windows=chr_windows_df,
        gff_index={},
    )


def write_bundle(
    result: SampleResult,
    out_dir: Path,
    emit: str,
    component_results: Sequence[SampleResult],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    tables_written: List[str] = []

    table_files = {
        "01_normalized.tsv": result.normalization.table,
        "02_domains.tsv": build_domain_table(result),
        "03_curated.tsv": result.curation.curated,
        "04_subfamilies.tsv": build_subfamily_table(result),
        "05_family_summary.tsv": result.family_summary,
        "06_subfamily_summary.tsv": result.subfamily_summary,
    }
    for name, df in table_files.items():
        write_table(df, out_dir / name, fmt=emit)
        tables_written.append(str(out_dir / name))

    if not result.chr_windows.empty:
        write_table(result.chr_windows, out_dir / "07_chr_windows.tsv", fmt=emit)
        tables_written.append(str(out_dir / "07_chr_windows.tsv"))

    fasta_written = write_representatives(out_dir, result.representatives, fmt=emit)
    tables_written.extend(fasta_written)

    fam_matrix = build_family_matrix(result.curation.curated, result.primary)
    subfam_matrix = build_subfamily_matrix(result.assignment_map, result.primary)
    pairwise = pairwise_comparison(fam_matrix)

    write_table(fam_matrix, out_dir / "09_family_matrix.tsv", fmt=emit)
    write_table(subfam_matrix, out_dir / "10_subfamily_matrix.tsv", fmt=emit)
    write_table(pairwise, out_dir / "11_pairwise_comparison.tsv", fmt=emit)
    tables_written.extend(
        [
            str(out_dir / "09_family_matrix.tsv"),
            str(out_dir / "10_subfamily_matrix.tsv"),
            str(out_dir / "11_pairwise_comparison.tsv"),
        ]
    )

    filtered_summary = out_dir / "filtered.summary.txt"
    write_filtered_summary(
        filtered_summary,
        result.normalization.non_ltr_filtered,
        result.curation.stats,
    )
    tables_written.append(str(filtered_summary))

    write_readme(out_dir / "README.txt")
    tables_written.append(str(out_dir / "README.txt"))

    context = build_all_context(result, component_results, fam_matrix, pairwise)
    write_all_txt(context, out_dir)
    tables_written.append(str(out_dir / "ALL.txt"))

    run_metadata = {
        "date_utc": context["date_utc"],
        "version": context["version"],
        "n_samples": context["n_samples"],
        "non_ltr_filtered": context["non_ltr_filtered"],
        "unknown_filtered": context["unknown_filtered"],
        "unknown_policy": context["unknown_policy"],
        "unknown_patterns": context["unknown_patterns"],
        "files": tables_written,
    }
    write_run_metadata(out_dir / "run.jsonl", run_metadata)


def build_domain_table(result: SampleResult) -> pd.DataFrame:
    curated = result.curation.curated[["seq_id", "family_curated"]]
    domain_df = result.domains.merge(curated, on="seq_id", how="left")
    sample_lookup = result.primary.set_index("seq_id")["sample"].to_dict()
    domain_df["sample"] = domain_df["seq_id"].map(sample_lookup).fillna("")
    columns = ["seq_id", "sample", "family_curated", "domains_csv", "intact_core", "truncated", "chimera_suspect"]
    return domain_df[columns]


def build_subfamily_table(result: SampleResult) -> pd.DataFrame:
    rows = []
    for assignment in result.assignments:
        rep = result.subfamily_rep_map.get(assignment.subfamily_id, "")
        rows.append(
            {
                "family_curated": assignment.family,
                "subfamily_id": assignment.subfamily_id,
                "n_members": len(assignment.members),
                "members_csv": ",".join(assignment.members),
                "method": assignment.method,
                "stability": assignment.stability,
                "representative_id": rep,
            }
        )
    columns = ["family_curated", "subfamily_id", "n_members", "members_csv", "method", "stability", "representative_id"]
    return pd.DataFrame(rows, columns=columns)


def write_representatives(out_dir: Path, representatives: Sequence[Representative], fmt: str) -> List[str]:
    fasta_path = out_dir / "08_representatives.fasta"
    index_path = out_dir / "08_representatives.index.tsv"
    if not representatives:
        fasta_path.write_text("", encoding="utf-8")
        empty_df = pd.DataFrame(columns=["name", "seq_id", "basis", "selection_reason", "sample"])
        write_table(empty_df, index_path, fmt=fmt)
        return [str(fasta_path), str(index_path)]

    with fasta_path.open("w", encoding="utf-8") as handle:
        for rep in representatives:
            header = f">{rep.name}|{rep.seq_id}|{rep.basis}|sample={rep.sample}"
            handle.write(header + "\n")
            wrapped = textwrap.wrap(rep.sequence, 80)
            for line in wrapped:
                handle.write(line + "\n")

    index_df = pd.DataFrame(
        [
            {
                "name": rep.name,
                "seq_id": rep.seq_id,
                "basis": rep.basis,
                "selection_reason": rep.selection_reason,
                "sample": rep.sample,
            }
            for rep in representatives
        ]
    )
    write_table(index_df, index_path, fmt=fmt)
    return [str(fasta_path), str(index_path)]


def write_filtered_summary(path: Path, non_ltr_filtered: int, stats: Dict[str, int]) -> None:
    lines = [
        f"non_LTR_filtered\t{non_ltr_filtered}",
        f"unknown_discarded\t{stats.get('unknown_discarded', 0)}",
        f"unknown_tagged\t{stats.get('unknown_tagged', 0)}",
        f"unknown_kept\t{stats.get('unknown_kept', 0)}",
    ]
    write_lines(lines, path)


def build_all_context(
    result: SampleResult,
    component_results: Sequence[SampleResult],
    family_mat: pd.DataFrame,
    pairwise_df: pd.DataFrame,
) -> Dict[str, object]:
    now = datetime.now(timezone.utc).isoformat()
    top_families = result.family_summary.sort_values("n_copies", ascending=False).head(10)
    dom_df = result.domains.merge(
        result.curation.curated[["seq_id", "family_curated"]],
        on="seq_id",
        how="left",
    )
    domain_block = []
    for _, entry in top_families.iterrows():
        family = entry["family_curated"]
        fam_domains = dom_df[dom_df["family_curated"] == family]
        if fam_domains.empty:
            domains_repr = ""
            intact_mean = 0.0
            trunc_mean = 0.0
        else:
            common = fam_domains["domains_csv"].mode()
            domains_repr = common.iat[0] if not common.empty else ""
            intact_mean = round(fam_domains["intact_core"].mean(), 3)
            trunc_mean = round(fam_domains["truncated"].mean(), 3)
        domain_block.append(
            {
                "family": family,
                "domains": domains_repr,
                "intact_pct": intact_mean,
                "truncated_pct": trunc_mean,
            }
        )
    subfamily_block = []
    for _, row in result.subfamily_summary.head(10).iterrows():
        subfamily_block.append(
            {
                "family": row["family_curated"],
                "subfamily": row["subfamily_id"],
                "n": row["n_copies"],
                "rep": row["representative_id"],
            }
        )

    chr_rows = result.chr_windows.to_dict(orient="records") if not result.chr_windows.empty else []

    sample_summaries = []
    for sample_res in component_results:
        keep = sample_res.curation.curated[sample_res.curation.curated["decision"].isin(["keep", "split"])]
        lengths = sample_res.primary.set_index("seq_id")["length"].to_dict()
        bp_total = int(sum(lengths.get(seq_id, 0) for seq_id in keep["seq_id"]))
        top = (
            sample_res.family_summary.sort_values("n_copies", ascending=False)
            .head(3)
            .set_index("family_curated")["n_copies"]
            .to_dict()
        )
        sample_summaries.append(
            {
                "sample": sample_res.config.sample,
                "species": sample_res.config.species,
                "copies_total": int(len(keep)),
                "bp_total": bp_total,
                "families_detected": int(keep["family_curated"].nunique()),
                "top_families": top,
            }
        )

    rows, cols = family_mat.shape
    sample_cols = [col for col in family_mat.columns if col.endswith(".count")]
    total_cells = len(family_mat) * max(1, len(sample_cols))
    non_zero = (
        family_mat[sample_cols] > 0
    ).sum().sum() if sample_cols else 0
    fill_pct = (non_zero / total_cells * 100.0) if total_cells else 0.0
    jaccard_entries = pairwise_df.drop_duplicates("pair")[["pair", "jaccard_over_families"]].to_dict(orient="records")
    cross = {
        "rows": rows,
        "cols": cols,
        "fill": f"{fill_pct:.1f}%",
        "top_shifted": [],
        "jaccard": [{"pair": entry["pair"], "value": entry["jaccard_over_families"]} for entry in jaccard_entries],
    }
    if not pairwise_df.empty:
        top_shifted = pairwise_df.reindex(pairwise_df["log2FC"].abs().sort_values(ascending=False).index).head(5)
        cross["top_shifted"] = [
            {"family": row["family"], "value": row["log2FC"], "pair": row["pair"]} for _, row in top_shifted.iterrows()
        ]

    return {
        "version": "0.1.0",
        "date_utc": now,
        "ltrs_only": result.config.ltrs_only,
        "unknown_policy": result.config.unknown_cfg.policy,
        "unknown_patterns": result.config.unknown_cfg.patterns,
        "non_ltr_filtered": result.normalization.non_ltr_filtered,
        "unknown_filtered": result.curation.stats.get("unknown_discarded", 0),
        "n_samples": len(component_results),
        "sample_summaries": sample_summaries,
        "domain_summary": domain_block,
        "subfamily_summary": subfamily_block,
        "chr_windows": chr_rows,
        "window": result.config.window,
        "cross_sample": cross,
        "filtered_summary": {
            "non_ltr": result.normalization.non_ltr_filtered,
            "unknown_drop": result.curation.stats.get("unknown_discarded", 0),
            "unknown_tag": result.curation.stats.get("unknown_tagged", 0),
        },
    }


def run_dump_command(args: argparse.Namespace) -> None:
    mapping = {
        "normalized": "01_normalized.tsv",
        "domains": "02_domains.tsv",
        "curated": "03_curated.tsv",
        "subfamilies": "04_subfamilies.tsv",
        "family_summary": "05_family_summary.tsv",
        "subfamily_summary": "06_subfamily_summary.tsv",
        "chr_windows": "07_chr_windows.tsv",
        "representatives_index": "08_representatives.index.tsv",
        "family_matrix": "09_family_matrix.tsv",
        "subfamily_matrix": "10_subfamily_matrix.tsv",
        "pairwise": "11_pairwise_comparison.tsv",
        "all_txt": "ALL.txt",
    }
    target = Path(args.out) / mapping[args.what]
    if not target.exists():
        raise SystemExit(f"Artifact not found: {target}")
    with target.open("r", encoding="utf-8") as handle:
        for line in handle:
            sys.stdout.write(line)


__all__ = ["build_parser", "dispatch"]
