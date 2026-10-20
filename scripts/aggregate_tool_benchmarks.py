#!/usr/bin/env python3
"""
Concatenate per-job Snakemake `benchmark:` TSVs (wall time + peak memory, one file per
tool rule invocation in tools.smk) into a single summary table for the HTML report's
Resource Usage tab.

The benchmark TSVs themselves carry no condition info (just Snakemake's fixed
s/h:m:s/max_rss/... columns — see snakemake.benchmark.BenchmarkRecord.get_header()), so
tool/data_type/cancer_type/relevance/num_muts/replicate are recovered from each file's
path, which follows the fixed layout written in tools.smk:
    {output_dir}/benchmarks/{data_type}/{cancer_type}/{relevance}/{num_muts}/rep_{rep}/{tool}_benchmark.tsv

Usage:
    python aggregate_tool_benchmarks.py \\
        --input-files file1.tsv file2.tsv ... \\
        --output all_tool_benchmarks.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_condition(path: str) -> dict:
    """Recover tool/data_type/cancer_type/relevance/num_muts/replicate from a benchmark
    file's path (see module docstring for the fixed layout). Returns {} if the path
    doesn't match (so the row can be skipped with a warning rather than crash)."""
    p = Path(path)
    rep_dir, num_muts, relevance, cancer_type, data_type = (
        p.parent.name, p.parent.parent.name, p.parent.parent.parent.name,
        p.parent.parent.parent.parent.name, p.parent.parent.parent.parent.parent.name,
    )
    if not rep_dir.startswith("rep_") or not p.name.endswith("_benchmark.tsv"):
        return {}
    return dict(
        tool=p.name[: -len("_benchmark.tsv")],
        data_type=data_type,
        cancer_type=cancer_type,
        relevance=relevance,
        num_muts=num_muts,
        replicate=rep_dir[len("rep_"):],
    )


def aggregate_benchmarks(input_files: list, output_file: str) -> None:
    print(f"Concatenating {len(input_files)} tool benchmark files...")
    rows = []
    for f in input_files:
        cond = parse_condition(f)
        if not cond:
            print(f"  WARNING: could not parse condition from path {f}", file=sys.stderr)
            continue
        try:
            record = pd.read_csv(f, sep='\t').iloc[0].to_dict()
        except Exception as e:
            print(f"  WARNING: could not read {f}: {e}", file=sys.stderr)
            continue
        rows.append({**cond, **record})

    if not rows:
        print("Error: no valid benchmark files found", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows)
    df['num_muts'] = pd.to_numeric(df['num_muts'], errors='coerce')
    df['replicate'] = pd.to_numeric(df['replicate'], errors='coerce')
    # Snakemake reports memory in MiB under the "max_rss" column name.
    df = df.rename(columns={'s': 'seconds', 'max_rss': 'max_rss_mb'})
    print(f"  Combined {len(rows)} files → {len(df)} rows, {df['tool'].nunique()} tools")

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep='\t', index=False)
    print(f"Tool benchmarks written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-files', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    aggregate_benchmarks(args.input_files, args.output)


if __name__ == '__main__':
    main()
