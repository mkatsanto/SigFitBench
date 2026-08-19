#!/usr/bin/env python3
"""
Aggregate all evaluation metrics into a single TSV with all wildcards as columns.
"""
import argparse
import pandas as pd
from pathlib import Path
import sys
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate all evaluation metrics into one table.")
    parser.add_argument('--input-files', nargs='+', required=True, help='List of all *_metrics.dat files')
    parser.add_argument('--output', required=True, help='Output TSV file')
    return parser.parse_args()

def parse_wildcards_from_filename(filename):
    # Example patterns:
    # 1. Regular: synthetic_Stomach-AdenoCA_10000_standard_SPA_w0_metrics.dat
    # 2. Penalty: synthetic_Stomach-AdenoCA_10000_standard_extended_penalty_0.1_w0_metrics.dat
    
    stem = Path(filename).stem
    
    # Try penalty pattern first: data_type_cancer_type_num_muts_relevance_tool_penalty_value_wRep_metrics
    m = re.match(
        r'(?P<data_type>[^_]+)_(?P<cancer_type>[^_]+)_(?P<num_muts>\d+)_(?P<relevance>[^_]+)_(?P<tool>.+?)_penalty_(?P<penalty>[\d.]+)_w(?P<rep>\d+)_metrics',
        stem
    )
    if m:
        wildcards = m.groupdict()
        print(f"DEBUG: Parsed penalty method {filename} -> tool={wildcards['tool']}, penalty={wildcards['penalty']}, rep={wildcards['rep']}", file=sys.stderr)
        return wildcards
    
    # Try regular pattern: data_type_cancer_type_num_muts_relevance_tool_wRep_metrics
    m = re.match(
        r'(?P<data_type>[^_]+)_(?P<cancer_type>[^_]+)_(?P<num_muts>\d+)_(?P<relevance>[^_]+)_(?P<tool>.+?)_w(?P<rep>\d+)_metrics',
        stem
    )
    if m:
        wildcards = m.groupdict()
        wildcards['penalty'] = None  # No penalty for regular methods
        print(f"DEBUG: Parsed regular method {filename} -> tool={wildcards['tool']}, rep={wildcards['rep']}", file=sys.stderr)
        return wildcards
    
    print(f"WARNING: Could not parse filename: {filename}", file=sys.stderr)
    return {}

def select_best_penalty(df):
    """
    For penalty methods, select the penalty value with the lowest average TAE.
    Returns a dataframe with only the best penalty value for each method.
    """
    # Separate penalty methods from regular methods
    has_penalty = df['penalty'].notna()
    regular_methods = df[~has_penalty].copy()
    penalty_methods = df[has_penalty].copy()
    
    if penalty_methods.empty:
        return df
    
    print(f"DEBUG: Found {len(penalty_methods)} rows with penalty values", file=sys.stderr)
    
    # Group by all wildcards except penalty and rep
    # Calculate average TAE for each (tool, cancer_type, num_muts, relevance, penalty) combination
    groupby_cols = ['data_type', 'cancer_type', 'num_muts', 'relevance', 'tool', 'penalty']
    
    # Check if TAE column exists
    if 'TAE' not in penalty_methods.columns:
        print("WARNING: TAE column not found, cannot select best penalty. Using all penalties.", file=sys.stderr)
        return df
    
    # Calculate average TAE for each penalty value
    penalty_tae = penalty_methods.groupby(groupby_cols)['TAE'].mean().reset_index()
    penalty_tae.rename(columns={'TAE': 'avg_TAE'}, inplace=True)
    
    print(f"DEBUG: Penalty TAE summary:\n{penalty_tae}", file=sys.stderr)
    
    # For each tool/cancer_type/num_muts/relevance combination, find the penalty with lowest avg TAE
    idx_cols = ['data_type', 'cancer_type', 'num_muts', 'relevance', 'tool']
    best_penalties = penalty_tae.loc[penalty_tae.groupby(idx_cols)['avg_TAE'].idxmin()]
    
    print(f"DEBUG: Best penalties selected:\n{best_penalties[idx_cols + ['penalty', 'avg_TAE']]}", file=sys.stderr)
    
    # Filter penalty_methods to keep only rows with the best penalty for each combination
    penalty_methods_filtered = []
    for _, row in best_penalties.iterrows():
        mask = True
        for col in idx_cols + ['penalty']:
            mask &= (penalty_methods[col] == row[col])
        filtered = penalty_methods[mask].copy()
        # Add the optimal penalty value as metadata (for hover text in reports)
        filtered['optimal_penalty'] = row['penalty']
        penalty_methods_filtered.append(filtered)
    
    if penalty_methods_filtered:
        penalty_methods_best = pd.concat(penalty_methods_filtered, ignore_index=True)
        print(f"DEBUG: Selected {len(penalty_methods_best)} rows with best penalties", file=sys.stderr)
    else:
        penalty_methods_best = pd.DataFrame()
    
    # Combine regular methods with best penalty methods
    result = pd.concat([regular_methods, penalty_methods_best], ignore_index=True)
    
    return result

def main():
    args = parse_args()
    all_dfs = []
    wildcard_cols = ['data_type', 'cancer_type', 'num_muts', 'relevance', 'tool']
    for f in args.input_files:
        try:
            df = pd.read_csv(f, sep='\t')
            wildcards = parse_wildcards_from_filename(f)
            # Ensure all wildcards are present as columns
            for col in wildcard_cols:
                value = wildcards.get(col, None)
                df[col] = value
            # Add penalty and rep columns
            df['penalty'] = wildcards.get('penalty', None)
            df['rep'] = wildcards.get('rep', None)
            # Preserve num_samples if it exists, otherwise try to add it
            if 'num_samples' not in df.columns:
                if 'sample' in df.columns:
                    num_samples = df['sample'].nunique()
                else:
                    num_samples = None
                df['num_samples'] = num_samples
            all_dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)
    if not all_dfs:
        print("No input files loaded!", file=sys.stderr)
        sys.exit(1)
    
    # Concatenate all data
    result = pd.concat(all_dfs, ignore_index=True)
    
    # Select best penalty for penalty methods
    result = select_best_penalty(result)
    
    # Reorder columns: wildcards and metadata first
    metadata_cols = wildcard_cols + ['num_samples', 'penalty', 'optimal_penalty', 'rep']
    cols = [c for c in metadata_cols if c in result.columns] + [c for c in result.columns if c not in metadata_cols]
    result = result[cols]
    
    result.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(result)} rows to {args.output}")

if __name__ == '__main__':
    main()
