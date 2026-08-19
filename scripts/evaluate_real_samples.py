#!/usr/bin/env python3
"""
Evaluate signature fitting results for real samples (without ground truth).

For real patient samples, we don't have ground truth signature weights,
but we can still calculate useful descriptive and consistency metrics:

Metrics calculated:
- num_active_sigs: Number of signatures with non-zero weights
- mean_activity: Mean signature weight across active signatures
- std_activity: Standard deviation of signature weights
- max_activity: Maximum signature weight
- entropy: Shannon entropy of signature distribution (diversity measure)
- sparsity: Fraction of signatures that are inactive (zeros)
- gini_coefficient: Inequality measure (0=perfect equality, 1=one sig dominates)

These metrics help characterize:
1. Signature diversity vs sparsity
2. Whether tools agree on signature assignments
3. Consistency of results across similar samples

Usage:
    python evaluate_real_samples.py \\
        --tool-results <path> \\
        --output <path> \\
        --tool <tool_name> \\
        --sample-id <identifier>
"""

import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import entropy


def calculate_gini_coefficient(weights):
    """Calculate Gini coefficient for signature weight distribution.
    
    Gini coefficient measures inequality in the distribution:
    - 0 = perfect equality (all signatures have equal weight)
    - 1 = perfect inequality (one signature has all weight)
    
    Args:
        weights: Array of signature weights (non-negative)
        
    Returns:
        float: Gini coefficient between 0 and 1
    """
    if len(weights) == 0 or np.sum(weights) == 0:
        return 0.0
    
    # Sort weights
    sorted_weights = np.sort(weights)
    n = len(sorted_weights)
    
    # Calculate Gini coefficient
    # Formula: G = (2 * sum(i * x_i)) / (n * sum(x_i)) - (n+1)/n
    cumsum = np.cumsum(sorted_weights)
    gini = (2 * np.sum((np.arange(1, n+1)) * sorted_weights)) / (n * np.sum(sorted_weights)) - (n + 1) / n
    
    return gini


def evaluate_real_sample(tool_results_df, tool_name, sample_id):
    """Calculate descriptive metrics for real sample results.
    
    Args:
        tool_results_df: DataFrame with signature weights (samples x signatures)
        tool_name: Name of the tool
        sample_id: Sample identifier
        
    Returns:
        dict: Evaluation metrics
    """
    # Calculate metrics for each sample
    sample_metrics = []
    
    for idx, row in tool_results_df.iterrows():
        weights = row.values
        
        # Filter out non-numeric columns if any
        if isinstance(weights[0], str):
            # Skip sample name column if present
            weights = pd.to_numeric(row[1:], errors='coerce')
            weights = weights[~np.isnan(weights)]
        
        weights = np.array(weights, dtype=float)
        total_weight = np.sum(weights)
        
        # Normalize weights
        if total_weight > 0:
            norm_weights = weights / total_weight
        else:
            norm_weights = weights
        
        # Calculate metrics
        metrics = {
            'num_active_sigs': np.sum(weights > 0),
            'mean_activity': np.mean(weights[weights > 0]) if np.any(weights > 0) else 0.0,
            'std_activity': np.std(weights[weights > 0]) if np.any(weights > 0) else 0.0,
            'max_activity': np.max(weights),
            'total_weight': total_weight,
            'entropy': entropy(norm_weights[norm_weights > 0]) if np.any(norm_weights > 0) else 0.0,
            'sparsity': np.sum(weights == 0) / len(weights),
            'gini_coefficient': calculate_gini_coefficient(weights)
        }
        
        sample_metrics.append(metrics)
    
    # Average across all samples
    metrics_df = pd.DataFrame(sample_metrics)
    avg_metrics = {
        'tool': tool_name,
        'sample_id': sample_id,
        'num_samples': len(tool_results_df),
        'num_active_sigs_mean': metrics_df['num_active_sigs'].mean(),
        'num_active_sigs_std': metrics_df['num_active_sigs'].std(),
        'mean_activity_mean': metrics_df['mean_activity'].mean(),
        'mean_activity_std': metrics_df['mean_activity'].std(),
        'max_activity_mean': metrics_df['max_activity'].mean(),
        'max_activity_std': metrics_df['max_activity'].std(),
        'total_weight_mean': metrics_df['total_weight'].mean(),
        'total_weight_std': metrics_df['total_weight'].std(),
        'entropy_mean': metrics_df['entropy'].mean(),
        'entropy_std': metrics_df['entropy'].std(),
        'sparsity_mean': metrics_df['sparsity'].mean(),
        'sparsity_std': metrics_df['sparsity'].std(),
        'gini_mean': metrics_df['gini_coefficient'].mean(),
        'gini_std': metrics_df['gini_coefficient'].std()
    }
    
    return avg_metrics


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate signature fitting results for real samples (no ground truth)'
    )
    parser.add_argument('--tool-results', required=True,
                        help='Tool results file (contribution matrix)')
    parser.add_argument('--output', required=True,
                        help='Output metrics file')
    parser.add_argument('--tool', required=True,
                        help='Tool name')
    parser.add_argument('--sample-id', required=True,
                        help='Sample identifier')
    
    args = parser.parse_args()
    
    try:
        # Load tool results
        results_df = pd.read_csv(args.tool_results, sep='\t', index_col=0)
        
        # Calculate metrics
        metrics = evaluate_real_sample(results_df, args.tool, args.sample_id)
        
        # Save metrics
        metrics_df = pd.DataFrame([metrics])
        
        # Create output directory if needed
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        
        # Save to file
        metrics_df.to_csv(args.output, sep='\t', index=False)
        
        print(f"Successfully evaluated {args.tool} for sample {args.sample_id}")
        print(f"Average active signatures: {metrics['num_active_sigs_mean']:.2f}")
        print(f"Average entropy: {metrics['entropy_mean']:.3f}")
        print(f"Average sparsity: {metrics['sparsity_mean']:.3f}")
        print(f"Average Gini coefficient: {metrics['gini_mean']:.3f}")
        
    except Exception as e:
        print(f"Error evaluating real sample: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
