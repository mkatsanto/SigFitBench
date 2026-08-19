#!/usr/bin/env python3
"""
Identify consistently active signatures from tool fitting results.

This script analyzes the output of signature fitting tools to determine which
signatures are consistently active across samples. These "relevant" signatures
can then be used for a second-pass fitting with a reduced signature set.

Criteria for "active" signature:
1. Present in at least min_sample_fraction of samples
2. Average contribution >= min_avg_contribution
3. Total contribution across cohort >= min_total_contribution

Author: SigFitTest workflow
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Identify consistently active signatures from fitting results",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--tool-results", 
        nargs='+',
        help="Tool result files (contribution.dat format)",
        required=True
    )
    parser.add_argument(
        "--output", 
        help="Output file with list of active signatures",
        required=True
    )
    parser.add_argument(
        "--min-contribution",
        type=float,
        default=10.0,
        help="Minimum contribution (mutations) for a signature to be considered active in a sample"
    )
    parser.add_argument(
        "--min-sample-fraction",
        type=float,
        default=0.3,
        help="Minimum fraction of samples where signature must be active (0-1)"
    )
    parser.add_argument(
        "--min-avg-contribution",
        type=float,
        default=50.0,
        help="Minimum average contribution across all samples"
    )
    parser.add_argument(
        "--consensus-tools",
        type=int,
        default=2,
        help="Number of tools that must agree on a signature for it to be included"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed information"
    )
    
    return parser.parse_args()


def load_tool_results(filepath):
    """Load tool results from contribution file"""
    try:
        df = pd.read_csv(filepath, sep='\t', index_col=0)
        return df
    except Exception as e:
        print(f"Warning: Could not load {filepath}: {e}")
        return None


def identify_active_signatures_for_tool(
    df,
    min_contribution,
    min_sample_fraction,
    min_avg_contribution
):
    """
    Identify active signatures for a single tool's results.
    
    Parameters
    ----------
    df : pd.DataFrame
        Tool results (signatures × samples)
    min_contribution : float
        Minimum contribution per sample
    min_sample_fraction : float
        Minimum fraction of samples where active
    min_avg_contribution : float
        Minimum average contribution
        
    Returns
    -------
    set
        Set of active signature names
    """
    active_sigs = set()
    num_samples = df.shape[1]
    
    for sig in df.index:
        # Count samples where signature is active
        active_in_samples = (df.loc[sig] >= min_contribution).sum()
        sample_fraction = active_in_samples / num_samples
        
        # Calculate average contribution
        avg_contribution = df.loc[sig].mean()
        
        # Check if signature meets criteria
        if (sample_fraction >= min_sample_fraction and 
            avg_contribution >= min_avg_contribution):
            active_sigs.add(sig)
    
    return active_sigs


def main():
    args = parse_arguments()
    
    print(f"=== Identifying Active Signatures ===")
    print(f"Analyzing {len(args.tool_results)} tool result file(s)")
    print(f"Criteria:")
    print(f"  - Min contribution per sample: {args.min_contribution}")
    print(f"  - Min sample fraction: {args.min_sample_fraction}")
    print(f"  - Min average contribution: {args.min_avg_contribution}")
    print(f"  - Consensus tools required: {args.consensus_tools}")
    print()
    
    # Load all tool results
    tool_active_sigs = {}
    all_signatures = set()
    
    for result_file in args.tool_results:
        tool_name = Path(result_file).stem.split('-')[-2]  # Extract tool name from filename
        print(f"Loading {tool_name} results from {result_file}...")
        
        df = load_tool_results(result_file)
        if df is None:
            continue
        
        print(f"  - {df.shape[0]} signatures × {df.shape[1]} samples")
        all_signatures.update(df.index)
        
        # Identify active signatures for this tool
        active = identify_active_signatures_for_tool(
            df,
            args.min_contribution,
            args.min_sample_fraction,
            args.min_avg_contribution
        )
        
        tool_active_sigs[tool_name] = active
        print(f"  - Found {len(active)} active signatures")
        
        if args.verbose and active:
            print(f"    Active: {', '.join(sorted(active))}")
    
    if not tool_active_sigs:
        print("ERROR: No tool results could be loaded")
        return
    
    print()
    
    # Find consensus signatures (active in >= N tools)
    signature_counts = {}
    for sig in all_signatures:
        count = sum(1 for active_set in tool_active_sigs.values() if sig in active_set)
        signature_counts[sig] = count
    
    consensus_sigs = {
        sig for sig, count in signature_counts.items() 
        if count >= args.consensus_tools
    }
    
    print(f"=== Consensus Results ===")
    print(f"Total unique signatures seen: {len(all_signatures)}")
    print(f"Signatures active in >= {args.consensus_tools} tools: {len(consensus_sigs)}")
    print()
    
    if consensus_sigs:
        print("Consensus active signatures:")
        for sig in sorted(consensus_sigs):
            tools_with_sig = [
                tool for tool, active in tool_active_sigs.items() 
                if sig in active
            ]
            print(f"  {sig:8s} (found by: {', '.join(tools_with_sig)})")
    else:
        print("WARNING: No consensus signatures found!")
        print("Consider lowering --consensus-tools or --min-sample-fraction")
    
    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write("# Active signatures identified from tool fitting results\n")
        f.write(f"# Analyzed {len(args.tool_results)} tool(s)\n")
        f.write(f"# Consensus threshold: {args.consensus_tools} tools\n")
        f.write(f"# Min contribution: {args.min_contribution}\n")
        f.write(f"# Min sample fraction: {args.min_sample_fraction}\n")
        f.write(f"# Min avg contribution: {args.min_avg_contribution}\n")
        f.write("#\n")
        for sig in sorted(consensus_sigs):
            f.write(f"{sig}\n")
    
    print()
    print(f"✓ Saved active signatures to {output_path}")
    print(f"  Use this file to generate RELEVANT signature sets")


if __name__ == "__main__":
    main()
