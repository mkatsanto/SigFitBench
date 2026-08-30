#!/usr/bin/env python3
"""
Concatenate per-sample evaluation TSVs and derive per-run summary statistics.

Reads all individual per-sample TSVs produced by evaluate_tool_performance,
concatenates them into all_per_sample_metrics.tsv, then groups by
(cancer_type, data_type, tool, replicate, num_muts) to compute
aggregate statistics across samples in each run.

Usage:
    python compute_summary_metrics.py \\
        --input-files file1.tsv file2.tsv ... \\
        --output-per-sample all_per_sample_metrics.tsv \\
        --output all_evaluation_metrics.tsv
"""

import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path


def compute_summary(input_files: list, per_sample_output: str, output_file: str) -> None:
    print(f"Concatenating {len(input_files)} per-sample metric files...")
    dfs = []
    for f in input_files:
        try:
            dfs.append(pd.read_csv(f, sep='\t'))
        except Exception as e:
            print(f"  WARNING: could not read {f}: {e}", file=sys.stderr)

    if not dfs:
        print("Error: no valid input files found", file=sys.stderr)
        sys.exit(1)

    df = pd.concat(dfs, ignore_index=True)
    print(f"  Combined {len(dfs)} files → {len(df)} rows")

    per_sample_path = Path(per_sample_output)
    per_sample_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(per_sample_path, sep='\t', index=False)
    print(f"Per-sample metrics written to {per_sample_output}")

    group_cols = ['cancer_type', 'data_type', 'tool', 'replicate', 'num_muts']
    agg = df.groupby(group_cols, dropna=False).agg(
        num_samples=('sample', 'count'),
        mean_num_muts=('num_muts', 'first'),
        TAE=('TAE', 'mean'),
        TAE_std=('TAE', 'std'),
        nRMSE=('nRMSE', lambda x: float(np.sqrt(np.mean(np.array(x) ** 2)))),
        weight_tot=('weight_tot', 'first'),
        n_eff=('n_eff', 'mean'),
        MAE_active=('MAE_active', 'mean'),
        wtot_FP=('wtot_FP', 'mean'),
        wtot_FP_std=('wtot_FP', 'std'),
        num_FP_sigs=('num_FP', 'mean'),
        wtot_FN=('wtot_FN', 'mean'),
        num_FN_sigs=('num_FN', 'mean'),
        P_mean=('P', 'mean'), P_std=('P', 'std'),
        R_mean=('R', 'mean'), R_std=('R', 'std'),
        S_mean=('S', 'mean'), S_std=('S', 'std'),
        F1_mean=('F1', 'mean'), F1_std=('F1', 'std'),
        MCC_mean=('MCC', 'mean'), MCC_std=('MCC', 'std'),
    ).reset_index()

    n = agg['num_samples'].clip(lower=1)
    for m in ['P', 'R', 'S', 'F1', 'MCC']:
        agg[f'{m}_sem'] = agg[f'{m}_std'] / np.sqrt(n)

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(output_path, sep='\t', index=False, float_format='%.4f')

    print(f"Summary metrics written: {len(agg)} rows")
    print(f"Tools: {sorted(agg['tool'].unique())}")


def main() -> None:
    parser = argparse.ArgumentParser(description='Concatenate per-sample TSVs and compute per-run summary metrics')
    parser.add_argument('--input-files', required=True, nargs='+',
                        help='Individual per-sample metric TSV files from evaluate_tool_performance')
    parser.add_argument('--output-per-sample', required=True,
                        help='Output path for concatenated per-sample TSV')
    parser.add_argument('--output', required=True,
                        help='Output path for summary statistics TSV')
    args = parser.parse_args()

    compute_summary(args.input_files, args.output_per_sample, args.output)


if __name__ == '__main__':
    main()
