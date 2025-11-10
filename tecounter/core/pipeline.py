"""
MIT License

High-level orchestration for TEcounter sample processing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from . import counts
from .cluster import SubfamilyAssignment, cluster_families
from .domains import call_ltr_domains
from .normalize import NormalizationResult, UnknownConfig, normalize_tesorter
from .qc import CurationResult, resolve_conflicts
from .representatives import Representative, select_representatives
from ..io.fasta import FastaResource
from ..io.gff import build_attribute_index, read_gff3
from ..util.logging import get_logger

LOGGER = get_logger()


@dataclass
class SampleConfig:
    tesorter: str
    sample: str
    species: str
    fasta: Optional[str]
    gff: Optional[str]
    ltrs_only: bool
    filter_unknowns: bool
    unknown_cfg: UnknownConfig
    min_cov: float
    min_score: float
    cluster_method: str
    k: int
    window: int
    strict_chimera: str


@dataclass
class SampleResult:
    config: SampleConfig
    normalization: NormalizationResult
    domains: pd.DataFrame
    curation: CurationResult
    primary: pd.DataFrame
    assignments: List[SubfamilyAssignment]
    assignment_map: Dict[str, str]
    representatives: List[Representative]
    subfamily_rep_map: Dict[str, str]
    family_summary: pd.DataFrame
    subfamily_summary: pd.DataFrame
    chr_windows: pd.DataFrame
    gff_index: Dict[str, object]


def run_sample(config: SampleConfig) -> SampleResult:
    """Execute the full pipeline for one sample and return structured data."""

    LOGGER.info("Processing sample %s", config.sample)
    normalization = normalize_tesorter(
        config.tesorter,
        sample=config.sample,
        species=config.species,
        ltrs_only=config.ltrs_only,
    )

    fasta_resource = FastaResource(config.fasta)

    if fasta_resource.available():
        seq_lengths = dict(fasta_resource.iter_lengths())
        if seq_lengths:
            mask = (normalization.table["length"] <= 1) | normalization.table["length"].isna()
            normalization.table.loc[mask, "length"] = (
                normalization.table.loc[mask, "seq_id"].map(seq_lengths).fillna(normalization.table.loc[mask, "length"])
            )

    domain_df = call_ltr_domains(normalization.table, min_cov=config.min_cov)
    curation = resolve_conflicts(
        normalization.table,
        domain_df,
        unknown_cfg=config.unknown_cfg,
        filter_unknowns=config.filter_unknowns,
        strict_chimera=config.strict_chimera,
        fasta=fasta_resource if fasta_resource.available() else None,
    )

    # Update normalized table with unknown flags
    if not curation.unknown_flags.empty:
        normalization.table["unknown_flag"] = normalization.table["seq_id"].map(curation.unknown_flags).fillna(0).astype(int)

    primary = counts.primary_hits(normalization.table)

    gff_index = {}
    if config.gff:
        gff_records = read_gff3(config.gff)
        gff_index = build_attribute_index(gff_records, attribute="ID")

    assignments, membership = cluster_families(
        curation.curated,
        fasta=fasta_resource if fasta_resource.available() else None,
        k=config.k,
        min_cluster_size=15,
        method=config.cluster_method,
    )

    representatives, subfamily_rep_map = select_representatives(
        curation.curated,
        primary,
        membership,
        fasta=fasta_resource if fasta_resource.available() else None,
        sample_name=config.sample,
    )

    family_summary_df = counts.family_summary(curation.curated, primary)
    subfamily_summary_df = counts.subfamily_summary(membership, subfamily_rep_map, primary)
    chr_windows_df = counts.chromosome_windows(
        curation.curated,
        primary,
        gff_index=gff_index,
        window=config.window,
    )

    return SampleResult(
        config=config,
        normalization=normalization,
        domains=curation.updated_domains,
        curation=curation,
        primary=primary,
        assignments=assignments,
        assignment_map=membership,
        representatives=representatives,
        subfamily_rep_map=subfamily_rep_map,
        family_summary=family_summary_df,
        subfamily_summary=subfamily_summary_df,
        chr_windows=chr_windows_df,
        gff_index=gff_index,
    )


__all__ = ["SampleConfig", "SampleResult", "run_sample"]
