# TEcounter CLI Reference

This guide explains every command exposed by `tecounter`, the expected inputs, and how to interpret the resulting outputs. All commands are executed via `python -m tecounter ...` unless you install the package (which also provides the `tecounter` console script).

## 1. Single-Sample Mode

### Synopsis
```
python -m tecounter \
  --tesorter /path/to/SAMPLE_LTRs.cls.tsv \
  --sample SAMPLE_ID \
  --species SPECIES_NAME \
  --out /path/to/output_dir \
  [--fasta /path/to/sample.fa] \
  [--gff /path/to/sample.gff3] \
  [--cluster {kmer,none}] \
  [--k 6] \
  [--threads 4] \
  [--ltrs-only true|false] \
  [--filter-unknowns true|false] \
  [--unknown-policy {drop,keep,tag}] \
  [--unknown-patterns pattern1,pattern2,...] \
  [--min-evidence 1] \
  [--strict-chimera {keep,drop,split}] \
  [--window 100000] \
  [--emit {tsv,csv,jsonl}]
```

### Required Parameters
- `--tesorter`: Path to a TEsorter classification table (TSV). Mixed orders are accepted; non-LTR rows are filtered unless `--ltrs-only false`.
- `--sample`: Identifier used throughout the bundle (e.g., `Athaliana`).
- `--species`: Descriptive label (e.g., `Arabidopsis_thaliana`).
- `--out`: Directory for the output bundle. The command creates the folder if it does not exist.

### Optional Inputs
- `--fasta`: FASTA file with TE sequences; enables sequence-aware tasks (length backfilling, ORF checks, representative FASTA). When omitted, tables still generate; FASTA-specific artifacts remain empty.
- `--gff`: GFF3 annotations mapping `seq_id -> coordinates`; required for `07_chr_windows.tsv`.

### Behavioral Flags
- `--cluster`: Subfamily detection strategy. `kmer` (default) runs deterministic clustering on each curated family. `none` skips clustering and emits header-only tables.
- `--k`: K-mer size for clustering (default 6).
- `--threads`: Reserved for future use; currently single-threaded but kept for CLI compatibility.
- `--ltrs-only`: When `true` (default) restricts processing to LTR orders/superfamilies.
- `--filter-unknowns`: Master toggle for unknown handling (default `true`).
- `--unknown-policy`: Action for rows flagged as unknown. `drop` excludes them entirely, `keep` keeps them but marks `unknown_flag=1`, and `tag` replaces `family_curated` with `LTR_Unknown` while retaining `family_raw`.
- `--unknown-patterns`: Comma-separated list matched against `family_raw` (case-insensitive). Defaults to `unknown,unk,unclassified,NA,?,.`.
- `--min-evidence`: Minimum number of domain hits before an `unknown` entry can survive the filter (default 1).
- `--strict-chimera`: Behavior for chimera suspects. `keep` (default) reports them but retains the row, `drop` excludes them from downstream tables, and `split` writes `decision=split` for manual review.
- `--window`: Window size in bp for chromosome summaries (default 100,000). Used only when GFF data is present.
- `--emit`: Serialization format (`tsv`, `csv`, `jsonl`) applied to every table.

### Outputs
Running the single-sample command produces a text bundle under `--out` with:
- `01_normalized.tsv` through `11_pairwise_comparison.tsv`
- `08_representatives.fasta` + index (if FASTA supplied)
- `ALL.txt`, `README.txt`, `filtered.summary.txt`, and `run.jsonl`

Each file name is fixed (see README for schemas).

## 2. Folder (Multi-Sample) Mode

### Synopsis
```
python -m tecounter folder \
  --in /path/to/tesorter_folder \
  --manifest /path/to/manifest.tsv \
  --out /path/to/output_dir \
  [--per-sample-out] \
  [same optional flags as single-sample]
```

### Required Inputs
- `--in`: Directory containing one TEsorter TSV per sample. File discovery is name-based: a manifest row with `sample_id=XYZ` expects `<in>/XYZ.tesorter.tsv` (or any suffix you provide, as long as the filenames match `sample_id`).
- `--manifest`: Tab-delimited table with at least `sample_id` and `species` columns. Optional columns (e.g., `group`, `assembly`) are preserved for bookkeeping but not consumed directly.
- `--out`: Destination for the combined bundle.

### Behavior
- TEcounter iterates over each manifest row, running the single-sample pipeline internally.
- When `--per-sample-out` is set, TEcounter writes each sample bundle under `<out>/<sample_id>/` in addition to creating the combined cross-sample bundle at `<out>/`.
- Combined artifacts include `09_family_matrix.tsv`, `10_subfamily_matrix.tsv`, and `11_pairwise_comparison.tsv` populated across all samples.

## 3. Dump Command

### Synopsis
```
python -m tecounter dump --out /path/to/bundle --what <artifact>
```

`--what` may be any of:
`normalized`, `domains`, `curated`, `subfamilies`, `family_summary`, `subfamily_summary`, `chr_windows`, `representatives_index`, `family_matrix`, `subfamily_matrix`, `pairwise`, `all_txt`.

This command streams the artifact to stdout, allowing shell filters or redirection, for example:
```
python -m tecounter dump --out out/Run --what family_summary | column -t | head
```

## 4. Flag Reference

| Flag | Description |
|------|-------------|
| `--tesorter` | Input TEsorter TSV. Accepts gzip-compressed files when the system `pandas` build supports it. |
| `--sample` | Sample identifier, propagated to every table/FASTA header. |
| `--species` | Human-readable label for `ALL.txt` and summary tables. |
| `--out` | Output directory; created if missing. |
| `--fasta` | Optional TE FASTA file. Required for `08_representatives.*` and for backfilling missing lengths. |
| `--gff` | Optional GFF3 with `ID` attributes matching `seq_id`. Enables `07_chr_windows.tsv`. |
| `--cluster` | Subfamily detection method. `kmer` runs deterministic clustering; `none` emits headers only. |
| `--k` | K-mer size used by clustering (default 6). |
| `--threads` | Reserved; currently single-threaded. |
| `--ltrs-only` | Restrict to LTR orders/superfamilies. Toggle to `false` to pass through all orders (still reported in provenance). |
| `--filter-unknowns` | Enable/disable unknown filtering logic. |
| `--unknown-policy` | Strategy for handling rows flagged as unknown (`drop`, `keep`, `tag`). |
| `--unknown-patterns` | Comma-separated patterns indicating unknown families. |
| `--min-evidence` | Require at least N domain hits for an unknown to survive filtering. |
| `--strict-chimera` | How to treat chimera suspects (`keep`, `drop`, `split`). |
| `--window` | Chromosome window size (bp) for `07_chr_windows.tsv`. |
| `--emit` | Serialization for tables (`tsv`, `csv`, `jsonl`). |
| `--strict-chimera` | Control chimera handling (keep/drop/split). |
| `--filter-unknowns false` | Bypass unknown filtering entirely. |
| `--ltrs-only false` | Include all TE orders. |

## 5. Example Workflow (Summary)

1. Run single sample (no FASTA/GFF required):
   ```
   python -m tecounter --tesorter /data/project/S1_LTRs.cls.tsv --sample S1 --species Species_one --out out/S1 --cluster none
   ```
2. Run folder mode with manifest:
   ```
   python -m tecounter folder --in /data/project/tesorter --manifest /data/project/manifest.tsv --out out/all --per-sample-out
   ```
3. Inspect family summary for the combined run:
   ```
   python -m tecounter dump --out out/all --what family_summary | column -t | head
   ```

Use these commands as templates and adapt paths/sample names to your dataset. Every run is deterministic when using the same inputs and flags.

## 6. Example Output (drMalDome5.1)

To illustrate the shape of the tables, we ran TEcounter on the `drMalDome5.1_LTRs.cls.tsv` dataset (sample ID `drMalDome5.1`, species `Malus_domestica`) with `--cluster none` and default filters. The resulting family summary is stored in `docs/examples/drMalDome5.1_family_summary.tsv`. The first few lines are reproduced below:

```
family_curated  n_copies  total_bp  median_len  truncated_pct  chimera_suspect_pct  mean_score  mean_cov
Retand          1550      1550      1.0         0.598          0.0                   0.0         0.701
Angela          996       996       1.0         1.0            0.0                   0.0         0.679
Tork            732       732       1.0         1.0            0.0                   0.0         0.924
Athila          1385      1385      1.0         0.505          0.0                   0.0         0.747
Tekay           636       636       1.0         0.553          0.0                   0.0         0.723
Ale             2093      2093      1.0         1.0            0.0                   0.0         0.894
```

The TSV in `docs/examples/` contains the complete table (23 families). Use it as a reference for column names and value ranges when integrating TEcounter outputs into downstream pipelines.
