# TEcounter (LTR-specialized)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17574578.svg)](https://doi.org/10.5281/zenodo.17574578)

TEcounter is a deterministic, text-only command-line workflow for summarising LTR retrotransposons from TEsorter outputs. It normalises mixed-order tables, applies LTR-specific domain and QC logic, clusters subfamilies via deterministic k-mer profiles, quantifies abundance (including optional chromosomal windows), and produces cross-sample comparisons. All outputs are plain text (TSV/CSV/JSONL selectable) plus FASTA representatives.

## Quickstart

```
python -m tecounter \
  --tesorter /path/to/SAMPLE_LTRs.cls.tsv \
  --sample SAMPLE_ID \
  --species SPECIES_NAME \
  --out out/SAMPLE_ID \
  --cluster kmer --k 6 \
  --ltrs-only true \
  --filter-unknowns true --unknown-policy drop
```

To process a directory of samples (requires manifest columns `sample_id`, `species`):

```
python -m tecounter folder \
  --in /path/to/tesorter_folder \
  --manifest /path/to/manifest.tsv \
  --out out/All \
  --per-sample-out \
  --filter-unknowns true --unknown-policy tag
```

 Dump a specific artifact from a completed run:

 ```
 python -m tecounter dump --out out/Athaliana --what family_matrix
 ```

## Features

- LTR-only processing by default with provenance counts for filtered non-LTR orders.
- Tolerant TSV ingestion (column aliasing) and streaming-friendly FASTA access via Biopython indices.
- LTR domain calls (`gag, PR, RT, RH, INT`) with intact/truncated logic and lightweight chimera detection.
- Unknown-handling controls: `--filter-unknowns`, `--unknown-policy`, `--unknown-patterns`, `--min-evidence`.
- Deterministic k-mer clustering (single-link, cosine distance) with bootstrap stability (50 replicates).
- Abundance summaries, chromosome windows (when GFF provided), cross-sample matrices, pairwise log2FC/Jaccard.
- Representative FASTA/index per family/subfamily (longest intact sequences, tagged with sample/basis).
- Publish-ready bundle containing REQUIRED text artifacts plus `ALL.txt`, `README.txt`, `run.jsonl` metadata.

## Project Layout

```
tecounter/           Python package (CLI + core modules)
pyproject.toml       Build metadata and dependencies
README.md            This guide
LICENSE              MIT license
```

Network-free execution is supported. All randomness is deterministically seeded (`random.seed(1337)` + NumPy RNGs).

## Command Reference

### Single-sample mode

```
python -m tecounter \
  --tesorter <path/to/sample.tesorter.tsv> \
  --sample <sample_id> \
  --species <Species_name> \
  --out <output_folder> \
  [--fasta sample.fa] [--gff sample.gff3] \
  [--cluster {kmer,none}] [--k 6] \
  [--ltrs-only true|false] \
  [--filter-unknowns true|false --unknown-policy {drop,keep,tag}] \
  [--unknown-patterns pattern1,pattern2,...] \
  [--min-evidence 1]
```

Mandatory arguments are the TEsorter TSV, the sample ID, species name, and the destination folder. FASTA/GFF are optional; omit them when you only need per-copy tabular summaries. Set `--cluster none` for speed when you do not need subfamily calls.

### Folder (multi-sample) mode

```
python -m tecounter folder \
  --in <directory_with_tesorter_tsvs> \
  --manifest <manifest.tsv> \
  --out <bundle_folder> \
  [--per-sample-out] \
  [other flags as in single-sample mode]
```

The manifest must contain at least `sample_id` and `species` columns. File discovery is name driven: a row with `sample_id=SampleA` expects `<in>/SampleA.tesorter.tsv`. When `--per-sample-out` is present, each sample also gets its own bundle under `<out>/<sample_id>/` in addition to the combined cross-sample bundle.

### Dump mode

```
python -m tecounter dump --out <bundle_dir> --what <artifact>
```

`--what` can be any of `normalized`, `domains`, `curated`, `subfamilies`, `family_summary`, `subfamily_summary`, `chr_windows`, `representatives_index`, `family_matrix`, `subfamily_matrix`, `pairwise`, or `all_txt`. The command streams the requested artifact to stdout, so you can redirect or pipe it as needed (e.g., `python -m tecounter dump --out out/Run --what family_matrix > fam.tsv`).
