#!/usr/bin/env python3
"""
Aggregate evaluation metrics across all tools for a specific dataset configuration.

This script combines metrics from all tools for one specific simulation condition
(e.g., synthetic data with 1000 mutations using standard signatures).

Usage:
    python aggregate_evaluation.py \\
        --input-files <metric_file1> <metric_file2> ... \\
        --output <summary_file>
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Aggregate evaluation metrics across tools"
    )
    parser.add_argument(
        "--input-files",
        nargs="+",
        required=True,
        help="Paths to individual tool metric files"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output summary file"
    )
    parser.add_argument(
        "--cancer-type",
        required=False,
        help="Cancer type being aggregated (optional, for logging/validation)"
    )
    return parser.parse_args()


def extract_tool_name(filepath):
    """Extract tool name from metric file path.
    
    Expected format: {data_type}_{num_muts}_{relevance}_{tool}_metrics.dat
    """
    filename = Path(filepath).stem
    # Remove _metrics suffix
    if filename.endswith('_metrics'):
        filename = filename[:-8]
    
    # Split by underscore and take the last part (tool name)
    parts = filename.split('_')
    if len(parts) >= 4:
        # Format is: data_type_num_relevance_tool
        return parts[-1]
    else:
        return filename


def extract_metadata(filepath):
    """Extract metadata from metric file path.
    
    Expected format: {data_type}_{identifier}_{relevance}_{tool}_metrics.dat
    where:
    - data_type can be 'synthetic', 'realsim', or 'samples'
    - identifier can be num_muts (e.g., "1000000") or sample_id
    - relevance is 'standard' or 'relevant'
    - tool is the tool name (e.g., 'MuSiCal', 'SPA')
    
    Returns:
        dict with 'data_type', 'identifier', 'relevance', 'tool', 'cancer_type' (if present)
    """
    filename = Path(filepath).stem
    # Remove _metrics suffix
    if filename.endswith('_metrics'):
        filename = filename[:-8]
    parts = filename.split('_')
    # Try to extract cancer_type if present (for synthetic: {data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_metrics.dat)
    cancer_type = None
    if len(parts) >= 5:
        # Format: {data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}
        data_type = parts[0]
        cancer_type = parts[1]
        identifier = parts[2]
        relevance = parts[3]
        tool = '_'.join(parts[4:])
    elif len(parts) >= 4:
        # Format: {data_type}_{identifier}_{relevance}_{tool}
        data_type = parts[0]
        identifier = parts[1]
        relevance = parts[2]
        tool = '_'.join(parts[3:])
    else:
        raise ValueError(f"Could not parse filename: {filepath}")
    # Apply prefix rules based on user requirements:
    if relevance == 'relevant':
        data_type = f'relevant_{data_type}'
    if data_type.endswith('samples'):
        data_type = data_type.replace('samples', 'real')
    elif data_type.endswith('relevant_samples'):
        data_type = data_type.replace('samples', 'real')
    meta = {
        'data_type': data_type,
        'identifier': identifier,
        'relevance': relevance,
        'tool': tool
    }
    if cancer_type is not None:
        meta['cancer_type'] = cancer_type
    return meta


def aggregate_metrics(input_files):
    """Aggregate metrics from multiple tool evaluation files.
    
    Args:
        input_files: List of paths to metric files
        
    Returns:
        DataFrame with aggregated metrics (one row per tool)
    """
    all_metrics = []
    
    for filepath in input_files:
        try:
            # Read the metrics file
            df = pd.read_csv(filepath, sep='\t')
            # Extract metadata (includes tool name and cancer_type if present)
            metadata = extract_metadata(filepath)
            # Add metadata as columns
            for key, value in metadata.items():
                df[key] = value
            all_metrics.append(df)
        except Exception as e:
            print(f"Warning: Failed to read {filepath}: {e}", file=sys.stderr)
            continue
    
    if not all_metrics:
        raise ValueError("No valid metric files were loaded")
    
    # Combine all metrics
    combined = pd.concat(all_metrics, ignore_index=True)
    
    # Convert identifier to numeric if it represents mutation counts
    # (will fail silently for sample IDs, which is fine)
    if 'identifier' in combined.columns:
        combined['num_muts'] = pd.to_numeric(combined['identifier'], errors='coerce')
    
    return combined


def main():
    """Main entry point."""
    args = parse_args()
    
    print(f"Aggregating metrics from {len(args.input_files)} files...")
    
    # Aggregate metrics
    summary = aggregate_metrics(args.input_files)
    # If multiple cancer types are present, keep the column for downstream use
    if 'cancer_type' in summary.columns:
        n_types = summary['cancer_type'].nunique()
        print(f"Cancer types in aggregation: {summary['cancer_type'].unique().tolist()}")
    else:
        n_types = 0
    # Sort by tool name and cancer_type for consistent output
    sort_cols = ['tool']
    if 'cancer_type' in summary.columns:
        sort_cols = ['cancer_type'] + sort_cols
    summary = summary.sort_values(sort_cols).reset_index(drop=True)
    # Save summary
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, sep='\t', index=False)
    print(f"Summary saved to {output_path}")
    print(f"Aggregated metrics for {len(summary)} tool(s) across {n_types} cancer type(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
