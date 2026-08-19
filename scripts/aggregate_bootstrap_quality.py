"""Aggregate per-file bootstrap quality TSVs into a single table.

Filenames are expected to follow the pattern:
    {data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_bootstrap_quality.tsv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

_PAT = re.compile(r'^([^_]+)_([^_]+)_(\d+)_([^_]+)_(.+)_w(\d+)$')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-files", nargs="*", default=[], help="Bootstrap quality TSV files to aggregate")
    p.add_argument("--output", required=True, help="Output TSV path")
    return p.parse_args()


def main():
    args = parse_args()

    if not args.input_files:
        Path(args.output).touch()
        print("No input files provided — writing empty output")
        return

    all_results = []
    for filepath in args.input_files:
        df = pd.read_csv(filepath, sep='\t')
        stem = Path(filepath).stem.replace('_bootstrap_quality', '')
        m = _PAT.match(stem)
        if m:
            df['data_type']   = m.group(1)
            df['cancer_type'] = m.group(2)
            df['num_muts']    = m.group(3)
            df['relevance']   = m.group(4)
            df['tool']        = m.group(5)
            df['realization'] = int(m.group(6))
        all_results.append(df)

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(args.output, sep='\t', index=False, float_format='%.6f')

    print(f"Aggregated {len(args.input_files)} bootstrap quality files")
    print(f"Total samples: {len(combined)}")
    if {'z_L1', 'z_L2', 'z_cos'}.issubset(combined.columns):
        print(
            f"Mean z-scores: z_L1={combined['z_L1'].mean():.2f}, "
            f"z_L2={combined['z_L2'].mean():.2f}, z_cos={combined['z_cos'].mean():.2f}"
        )


if __name__ == "__main__":
    main()
