#!/usr/bin/env python3
"""
Aggregate per-sample evaluation metrics from multiple per-run TSV files.

All metadata (cancer_type, data_type, tool, penalty, replicate, num_muts) is already
present as columns in each input file written by evaluate_tool_performance.py,
so this script is a plain concatenation.

Usage:
    python aggregate_per_sample_metrics.py \\
        --input-files file1.tsv file2.tsv ... \\
        --output combined_per_sample.tsv
"""

import argparse
import pandas as pd
from pathlib import Path
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description='Aggregate per-sample evaluation metrics')
    parser.add_argument('--input-files', nargs='+', required=True,
                        help='Per-sample metric TSV files to aggregate')
    parser.add_argument('--output', required=True,
                        help='Output file path for combined metrics')
    args = parser.parse_args()

    print(f"Aggregating {len(args.input_files)} per-sample metric files...", file=sys.stderr)

    dfs = []
    for f in args.input_files:
        try:
            dfs.append(pd.read_csv(f, sep='\t'))
        except Exception as e:
            print(f"  WARNING: Could not read {f}: {e}", file=sys.stderr)

    if not dfs:
        print("Error: No valid input files found", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(dfs, ignore_index=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, sep='\t', index=False)

    print(f"Combined {len(dfs)} files → {len(combined)} rows", file=sys.stderr)
    print(f"Tools: {sorted(combined['tool'].unique())}", file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
