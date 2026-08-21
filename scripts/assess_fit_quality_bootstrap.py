#!/usr/bin/env python3
"""
Assess signature fitting quality using bootstrap resampling.

This script implements the bootstrap-based quality assessment from the original
SigFitTest fit_and_assess() function. For each fitted sample, it:
1. Generates synthetic bootstrap samples using the fitted signature activities
2. Refits all bootstrap samples with the same tool
3. Calculates reconstruction metrics (L1, L2, cosine) for original and bootstrap fits
4. Computes z-scores to assess if the original fit quality is unusual

The z-score indicates whether the fit is unusually good/bad compared to what
we'd expect if the estimated signature activities were correct.

Usage:
    python assess_fit_quality_bootstrap.py \\
        --tool-results <fitted_weights.dat> \\
        --input-catalog <counts.dat> \\
        --reference-signatures <cosmic.txt> \\
        --output <assessment.tsv> \\
        --tool <tool_name> \\
        --num-bootstrap 100 \\
        --random-seed 42

Output format (TSV):
    sample  num_muts  L1  L2  cos  L1_bootstrap  L2_bootstrap  cos_bootstrap  z_L1  z_L2  z_cos
"""

import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial.distance import cosine as cosine_distance


def load_tool_results(results_file):
    """Load tool results (fitted signature weights)."""
    try:
        df = pd.read_csv(results_file, sep=',')
        if 'signature' in df.columns:
            df = df.set_index('signature', drop=True)
        elif 'Unnamed: 0' in df.columns:
            df = df.rename(columns={'Unnamed: 0': 'signature'})
            df = df.set_index('signature', drop=True)
        elif df.columns[0] == '0' or (len(df) > 0 and str(df.iloc[0, 0]).startswith('SBS')):
            df = df.set_index(df.columns[0], drop=True)
        return df
    except Exception:
        try:
            df = pd.read_csv(results_file, sep='\t', index_col=0)
            return df
        except Exception as e:
            print(f"Error loading tool results from {results_file}: {e}", file=sys.stderr)
            sys.exit(1)


def load_catalog(catalog_file):
    """Load mutation catalog (96 mutation types × samples)."""
    try:
        df = pd.read_csv(catalog_file, sep='\t', index_col=0)
        return df
    except Exception:
        try:
            df = pd.read_csv(catalog_file, sep=',', index_col=0)
            return df
        except Exception as e:
            print(f"Error loading catalog from {catalog_file}: {e}", file=sys.stderr)
            sys.exit(1)


def load_reference_signatures(ref_file):
    """Load reference signature catalog (COSMIC signatures)."""
    try:
        df = pd.read_csv(ref_file, sep='\t', index_col='Type')
        return df
    except Exception:
        try:
            df = pd.read_csv(ref_file, sep=',')
            if 'Type' in df.columns:
                df = df.set_index('Type', drop=True)
            elif df.columns[0] in ['sub', 'tri']:
                df = df.set_index(df.columns[0], drop=True)
            return df
        except Exception as e:
            print(f"Error loading reference signatures from {ref_file}: {e}", file=sys.stderr)
            sys.exit(1)


def generate_bootstrap_samples(signature_catalog, fitted_weights, num_mutations, num_bootstrap, rng):
    """
    Generate bootstrap samples from fitted signature activities.
    
    Args:
        signature_catalog: DataFrame (96 mutation types × signatures)
        fitted_weights: Series (signature activities for one sample)
        num_mutations: int (total mutations to generate)
        num_bootstrap: int (number of bootstrap samples)
        rng: numpy random generator
        
    Returns:
        DataFrame (96 mutation types × num_bootstrap samples)
    """
    # Normalize fitted weights to probabilities
    if fitted_weights.sum() == 0:
        print(f"  Warning: Fitted weights sum to zero, cannot generate bootstrap", file=sys.stderr)
        return None
    
    signature_probs = fitted_weights / fitted_weights.sum()
    
    # Ensure signature catalog has all fitted signatures
    # Filter out unknown signatures and signatures not in reference catalog
    valid_sigs = []
    for sig in fitted_weights.index:
        if sig.startswith('Unknown_'):
            # Skip unknown signatures for bootstrap
            continue
        if sig not in signature_catalog.columns:
            print(f"  Warning: Signature {sig} not in reference catalog, skipping", file=sys.stderr)
            continue
        valid_sigs.append(sig)
    
    if len(valid_sigs) == 0:
        print(f"  Warning: No valid signatures found for bootstrap (all unknown or missing)", file=sys.stderr)
        return None
    
    # Keep only valid signatures and renormalize
    signature_probs = fitted_weights[valid_sigs]
    signature_probs = signature_probs / signature_probs.sum()
    
    # Calculate expected profile from fitted activities
    expected_profile = signature_catalog[signature_probs.index].dot(signature_probs)
    expected_profile = expected_profile / expected_profile.sum()
    
    # Verify expected profile is valid
    if np.any(np.isnan(expected_profile)) or np.any(expected_profile < 0) or np.any(expected_profile > 1):
        print(f"  Warning: Invalid expected profile (NaN or out of bounds), cannot generate bootstrap", file=sys.stderr)
        return None
    if not np.isclose(expected_profile.sum(), 1.0, atol=1e-6):
        print(f"  Warning: Expected profile doesn't sum to 1 ({expected_profile.sum()}), cannot generate bootstrap", file=sys.stderr)
        return None
    
    # Generate bootstrap samples using multinomial distribution
    bootstrap_samples = []
    for _ in range(num_bootstrap):
        counts = rng.multinomial(num_mutations, expected_profile)
        bootstrap_samples.append(counts)
    
    bootstrap_df = pd.DataFrame(
        bootstrap_samples,
        columns=signature_catalog.index
    ).T
    bootstrap_df.columns = [f'bootstrap_{i}' for i in range(num_bootstrap)]
    
    return bootstrap_df


def calculate_reconstruction_metrics(observed_counts, signature_catalog, fitted_weights):
    """
    Calculate reconstruction quality metrics (L1, L2, cosine similarity).
    
    Args:
        observed_counts: Series (96 mutation counts)
        signature_catalog: DataFrame (96 mutation types × signatures)
        fitted_weights: Series (signature activities)
        
    Returns:
        dict with L1, L2, cos metrics
    """
    num_muts = observed_counts.sum()
    if num_muts == 0:
        return {'L1': np.nan, 'L2': np.nan, 'cos': np.nan}
    
    # Normalize observed counts
    observed_profile = observed_counts / num_muts
    
    # Reconstruct profile from fitted weights
    fitted_weights_normalized = fitted_weights / fitted_weights.sum() if fitted_weights.sum() > 0 else fitted_weights
    
    # Ensure all signatures in fitted_weights are in signature_catalog
    valid_sigs = [sig for sig in fitted_weights_normalized.index if sig in signature_catalog.columns]
    if len(valid_sigs) == 0:
        return {'L1': np.nan, 'L2': np.nan, 'cos': np.nan}
    
    reconstructed_profile = signature_catalog[valid_sigs].dot(fitted_weights_normalized[valid_sigs])
    
    # Calculate metrics
    L1 = (observed_profile - reconstructed_profile).abs().sum()
    L2 = np.sqrt(np.power(observed_profile - reconstructed_profile, 2).sum())
    
    # Cosine similarity (1 - cosine distance)
    if reconstructed_profile.sum() > 0:
        cos = 1 - cosine_distance(observed_profile, reconstructed_profile)
    else:
        cos = 0
    
    return {'L1': L1, 'L2': L2, 'cos': cos}


def fit_bootstrap_samples(bootstrap_catalog, signature_catalog, tool, temp_dir):
    """
    Fit bootstrap samples using the specified tool.
    
    This is a simplified fitting that uses basic NNLS (non-negative least squares)
    for compatibility across all tools without needing to call external R/Python scripts.
    
    Args:
        bootstrap_catalog: DataFrame (96 types × num_bootstrap samples)
        signature_catalog: DataFrame (96 types × signatures)
        tool: str (tool name, for logging)
        temp_dir: Path (for temporary files if needed)
        
    Returns:
        DataFrame (signatures × num_bootstrap samples) with fitted weights
    """
    from scipy.optimize import nnls
    
    print(f"  Fitting {bootstrap_catalog.shape[1]} bootstrap samples with NNLS...", file=sys.stderr)
    
    # Align mutation types
    common_types = bootstrap_catalog.index.intersection(signature_catalog.index)
    bootstrap_catalog = bootstrap_catalog.loc[common_types]
    signature_catalog = signature_catalog.loc[common_types]
    
    # Fit each bootstrap sample
    fitted_weights_list = []
    
    for sample_col in bootstrap_catalog.columns:
        counts = bootstrap_catalog[sample_col].values
        
        # NNLS fitting: minimize ||counts - W @ h||^2 subject to h >= 0
        # W = signature_catalog (96 × num_sigs), h = weights (num_sigs,)
        weights, residual = nnls(signature_catalog.values, counts)
        
        fitted_weights_list.append(weights)
    
    fitted_df = pd.DataFrame(
        np.array(fitted_weights_list).T,
        index=signature_catalog.columns,
        columns=bootstrap_catalog.columns
    )
    
    return fitted_df


def assess_sample_quality(sample_name, observed_counts, fitted_weights, signature_catalog,
                          num_bootstrap, rng, tool, temp_dir):
    """
    Assess fit quality for a single sample using bootstrap.
    
    Returns:
        dict with assessment metrics
    """
    num_muts = observed_counts.sum()
    
    print(f"  Processing sample: {sample_name} ({num_muts} mutations)", file=sys.stderr)
    
    # Calculate original reconstruction metrics
    original_metrics = calculate_reconstruction_metrics(observed_counts, signature_catalog, fitted_weights)
    
    print(f"    Original: L1={original_metrics['L1']:.4f}, L2={original_metrics['L2']:.4f}, cos={original_metrics['cos']:.4f}",
          file=sys.stderr)
    
    # Generate bootstrap samples
    bootstrap_catalog = generate_bootstrap_samples(
        signature_catalog, fitted_weights, int(num_muts), num_bootstrap, rng
    )
    
    if bootstrap_catalog is None:
        return {
            'sample': sample_name,
            'num_muts': num_muts,
            **original_metrics,
            'L1_bootstrap_mean': np.nan,
            'L1_bootstrap_std': np.nan,
            'L2_bootstrap_mean': np.nan,
            'L2_bootstrap_std': np.nan,
            'cos_bootstrap_mean': np.nan,
            'cos_bootstrap_std': np.nan,
            'z_L1': np.nan,
            'z_L2': np.nan,
            'z_cos': np.nan
        }
    
    # Fit bootstrap samples
    bootstrap_fitted = fit_bootstrap_samples(bootstrap_catalog, signature_catalog, tool, temp_dir)
    
    # Calculate reconstruction metrics for each bootstrap sample
    bootstrap_L1 = []
    bootstrap_L2 = []
    bootstrap_cos = []
    
    for boot_col in bootstrap_catalog.columns:
        boot_metrics = calculate_reconstruction_metrics(
            bootstrap_catalog[boot_col],
            signature_catalog,
            bootstrap_fitted[boot_col]
        )
        bootstrap_L1.append(boot_metrics['L1'])
        bootstrap_L2.append(boot_metrics['L2'])
        bootstrap_cos.append(boot_metrics['cos'])
    
    bootstrap_L1 = np.array(bootstrap_L1)
    bootstrap_L2 = np.array(bootstrap_L2)
    bootstrap_cos = np.array(bootstrap_cos)
    
    # Calculate z-scores
    # z = (original - bootstrap_mean) / bootstrap_std
    # Positive z for L1/L2 means original fit is worse than expected
    # Negative z for cos means original fit is worse than expected
    z_L1 = (original_metrics['L1'] - bootstrap_L1.mean()) / bootstrap_L1.std() if bootstrap_L1.std() > 0 else 0
    z_L2 = (original_metrics['L2'] - bootstrap_L2.mean()) / bootstrap_L2.std() if bootstrap_L2.std() > 0 else 0
    # For cosine: we want higher is better, so flip the sign
    z_cos = (bootstrap_cos.mean() - original_metrics['cos']) / bootstrap_cos.std() if bootstrap_cos.std() > 0 else 0
    
    print(f"    Bootstrap: L1={bootstrap_L1.mean():.4f}±{bootstrap_L1.std():.4f}, "
          f"L2={bootstrap_L2.mean():.4f}±{bootstrap_L2.std():.4f}, "
          f"cos={bootstrap_cos.mean():.4f}±{bootstrap_cos.std():.4f}", file=sys.stderr)
    print(f"    Z-scores: z_L1={z_L1:.2f}, z_L2={z_L2:.2f}, z_cos={z_cos:.2f}", file=sys.stderr)
    
    return {
        'sample': sample_name,
        'num_muts': num_muts,
        'L1': original_metrics['L1'],
        'L2': original_metrics['L2'],
        'cos': original_metrics['cos'],
        'L1_bootstrap_mean': bootstrap_L1.mean(),
        'L1_bootstrap_std': bootstrap_L1.std(),
        'L2_bootstrap_mean': bootstrap_L2.mean(),
        'L2_bootstrap_std': bootstrap_L2.std(),
        'cos_bootstrap_mean': bootstrap_cos.mean(),
        'cos_bootstrap_std': bootstrap_cos.std(),
        'z_L1': z_L1,
        'z_L2': z_L2,
        'z_cos': z_cos
    }


def main():
    parser = argparse.ArgumentParser(
        description='Assess signature fitting quality using bootstrap resampling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--tool-results', required=True,
                        help='Path to fitted signature weights file')
    parser.add_argument('--input-catalog', required=True,
                        help='Path to input mutation catalog')
    parser.add_argument('--reference-signatures', required=True,
                        help='Path to reference signature catalog (COSMIC)')
    parser.add_argument('--output', required=True,
                        help='Output file for quality assessment results')
    parser.add_argument('--tool', required=True,
                        help='Tool name (for logging)')
    parser.add_argument('--num-bootstrap', type=int, default=100,
                        help='Number of bootstrap samples per original sample (default: 100)')
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--temp-dir', default=None,
                        help='Temporary directory for intermediate files')
    
    args = parser.parse_args()
    
    # Set up random number generator
    rng = np.random.default_rng(args.random_seed)
    
    # Create temp directory
    temp_dir = Path(args.temp_dir) if args.temp_dir else Path(args.output).parent / 'tmp'
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Bootstrap quality assessment for {args.tool}", file=sys.stderr)
    print(f"  Tool results: {args.tool_results}", file=sys.stderr)
    print(f"  Input catalog: {args.input_catalog}", file=sys.stderr)
    print(f"  Reference signatures: {args.reference_signatures}", file=sys.stderr)
    print(f"  Number of bootstrap samples: {args.num_bootstrap}", file=sys.stderr)
    print(f"  Random seed: {args.random_seed}", file=sys.stderr)
    
    # Load data
    print("Loading data...", file=sys.stderr)
    fitted_weights = load_tool_results(args.tool_results)
    input_catalog = load_catalog(args.input_catalog)
    signature_catalog = load_reference_signatures(args.reference_signatures)
    
    print(f"  Fitted weights: {fitted_weights.shape[0]} signatures × {fitted_weights.shape[1]} samples", file=sys.stderr)
    print(f"  Input catalog: {input_catalog.shape[0]} mutation types × {input_catalog.shape[1]} samples", file=sys.stderr)
    print(f"  Reference signatures: {signature_catalog.shape[0]} mutation types × {signature_catalog.shape[1]} signatures",
          file=sys.stderr)
    
    # Ensure samples match
    common_samples = fitted_weights.columns.intersection(input_catalog.columns)
    if len(common_samples) == 0:
        print("ERROR: No common samples found between fitted weights and input catalog", file=sys.stderr)
        sys.exit(1)
    
    if len(common_samples) < fitted_weights.shape[1]:
        print(f"  Warning: Only {len(common_samples)}/{fitted_weights.shape[1]} samples found in both files",
              file=sys.stderr)
    
    fitted_weights = fitted_weights[common_samples]
    input_catalog = input_catalog[common_samples]
    
    # Assess quality for each sample
    print(f"\nAssessing fit quality for {len(common_samples)} samples...", file=sys.stderr)
    results = []
    
    for i, sample in enumerate(common_samples):
        print(f"\n[{i+1}/{len(common_samples)}]", file=sys.stderr)
        
        sample_result = assess_sample_quality(
            sample_name=sample,
            observed_counts=input_catalog[sample],
            fitted_weights=fitted_weights[sample],
            signature_catalog=signature_catalog,
            num_bootstrap=args.num_bootstrap,
            rng=rng,
            tool=args.tool,
            temp_dir=temp_dir
        )
        results.append(sample_result)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    
    print(f"\nResults saved to {args.output}", file=sys.stderr)
    print(f"  Processed {len(results)} samples", file=sys.stderr)
    
    # Summary statistics
    print("\nSummary statistics:", file=sys.stderr)
    print(f"  Mean z_L1: {results_df['z_L1'].mean():.2f} ± {results_df['z_L1'].std():.2f}", file=sys.stderr)
    print(f"  Mean z_L2: {results_df['z_L2'].mean():.2f} ± {results_df['z_L2'].std():.2f}", file=sys.stderr)
    print(f"  Mean z_cos: {results_df['z_cos'].mean():.2f} ± {results_df['z_cos'].std():.2f}", file=sys.stderr)
    
    # Flag unusual samples (|z| > 2)
    unusual_samples = results_df[(results_df['z_L1'].abs() > 2) | 
                                  (results_df['z_L2'].abs() > 2) | 
                                  (results_df['z_cos'].abs() > 2)]
    if len(unusual_samples) > 0:
        print(f"\n  WARNING: {len(unusual_samples)} samples with unusual fit quality (|z| > 2):", file=sys.stderr)
        for _, row in unusual_samples.iterrows():
            print(f"    {row['sample']}: z_L1={row['z_L1']:.2f}, z_L2={row['z_L2']:.2f}, z_cos={row['z_cos']:.2f}",
                  file=sys.stderr)


if __name__ == '__main__':
    main()
