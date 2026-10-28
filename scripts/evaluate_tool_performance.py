#!/usr/bin/env python3
"""
Evaluate signature fitting performance for a specific tool and dataset.

This script corresponds to the eval_results() function from the original SigFitTest
(MS_run_and_evaluate.py). It calculates comprehensive evaluation metrics by comparing
fitted signature weights against ground truth.

Metrics calculated:
- TAE: Total Absolute Error (mean fitting error across samples)
- MAE_active: Mean Absolute Error for active signatures only
- Precision, Recall, F1, MCC: Classification metrics for signature detection
- FP/FN statistics: False positive and false negative signature counts and weights
- n_eff: Effective number of signatures per sample
- weight_tot: Total assigned signature weight (normalized)

Usage (single tool):
    python evaluate_tool_performance.py \\
        --tool-results <path> \\
        --true-weights <path> \\
        --reference-signatures <path> \\
        --output <path> \\
        --tool <tool_name> \\
        --num-muts <mutation_count> \\
        --min-contribution <threshold>

Usage (multiple tools with plots):
    python evaluate_tool_performance.py \\
        --tool-results <path1> <path2> ... \\
        --tools <tool1> <tool2> ... \\
        --true-weights <path> \\
        --reference-signatures <path1> <path2> ... \\
        --output <path> \\
        --num-muts <mutation_count> \\
        --plots-dir <directory> \\
        --min-contribution <threshold>
"""

import argparse
import sys
from typing import Tuple, List, Optional, Dict, Any
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================================
# Constants
# ============================================================================

# Minimum mutation contribution to consider a signature active (default)
DEFAULT_MIN_CONTRIBUTION = 10.0

# Threshold for considering signature weights equivalent in comparisons
WEIGHT_EPSILON = 1e-10


# ============================================================================
# ENP helpers
# ============================================================================

def calculate_enp_from_channels(channel_proportions: np.ndarray) -> float:
    """Return ENP = 1/sum(p_i^2) for a signature's 96-channel profile.

    High ENP → flat/diffuse spectrum (harder to attribute).
    Low ENP  → peaked/concentrated spectrum (easier to attribute).
    """
    p = np.asarray(channel_proportions, dtype=float)
    p = p[p > 0]
    if len(p) == 0:
        return np.nan
    p = p / p.sum()
    return 1.0 / float(np.sum(p ** 2))


def load_cosmic_profiles(cosmic_file: str) -> Optional[pd.DataFrame]:
    """Load the full COSMIC signature profile matrix (channels × signatures)."""
    try:
        return pd.read_csv(cosmic_file, sep='\t', index_col=0)
    except Exception as e:
        print(f"Warning: could not load COSMIC profiles from {cosmic_file}: {e}", file=sys.stderr)
        return None


# ============================================================================
# Helper Functions
# ============================================================================

def align_signature_dataframes(
    df_fitted: pd.DataFrame,
    df_true: pd.DataFrame,
    ref_signatures: List[str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Align fitted and true signature DataFrames to have matching signature sets.
    
    Ensures both dataframes contain:
    - All reference signatures
    - All signatures present in either dataframe
    - Signatures in the same sorted order
    
    Missing signatures are added with zero values.
    
    Args:
        df_fitted: DataFrame with fitted signature weights (signatures × samples)
        df_true: DataFrame with true signature weights (signatures × samples)
        ref_signatures: List of reference signatures that should be included
        
    Returns:
        Tuple of (aligned_fitted, aligned_true) DataFrames with matching indices
    """
    # Ensure all reference signatures are included in both dataframes
    for sig in ref_signatures:
        if sig not in df_true.index:
            df_true.loc[sig] = 0
        if sig not in df_fitted.index:
            df_fitted.loc[sig] = 0
    
    # Ensure all signatures from results are in true weights and vice versa
    for sig in df_fitted.index:
        if sig not in df_true.index:
            df_true.loc[sig] = 0
    
    for sig in df_true.index:
        if sig not in df_fitted.index:
            df_fitted.loc[sig] = 0
    
    # Sort indices to match
    df_fitted = df_fitted.sort_index()
    df_true = df_true.sort_index()
    
    return df_fitted, df_true


def tool_failure(results_file: str) -> Optional[str]:
    """The reason if the tool recorded a failure in its result file ("# TOOL FAILED: ..."),
    otherwise None."""
    with open(results_file) as f:
        first_line = f.readline().strip()
    if first_line.startswith("# TOOL FAILED"):
        return first_line.split(":", 1)[-1].strip() or "see the tool's log"
    return None


def load_tool_results(results_file: str) -> pd.DataFrame:
    """Load tool results (fitted signature weights).
    Expected format: CSV or TSV with signatures as rows, samples as columns.
    First column should be signature names (may be named 'signature', 'Unnamed: 0', etc.)
    """
    try:
        # Try comma-separated first (most common)
        df = pd.read_csv(results_file, sep=',')
        # Handle various index column names
        if 'signature' in df.columns:
            df = df.set_index('signature', drop=True)
        elif 'Unnamed: 0' in df.columns:
            df = df.rename(columns={'Unnamed: 0': 'signature'})
            df = df.set_index('signature', drop=True)
        elif df.columns[0] == '0' or df.iloc[0, 0] in ['SBS1', 'SBS2', 'SBS3']:  # First column looks like signatures
            df = df.set_index(df.columns[0], drop=True)
        return df
    except Exception:
        # Try tab-separated
        try:
            df = pd.read_csv(results_file, sep='\t', index_col=0)
            return df
        except Exception as e:
            print(f"Error loading tool results from {results_file}: {e}", file=sys.stderr)
            sys.exit(1)


def load_true_weights(weights_file):
    """Load ground truth signature weights.
    
    Expected format: TSV with signatures as rows, samples as columns.
    """
    try:
        df = pd.read_csv(weights_file, sep='\t', index_col=0)
        return df
    except Exception as e:
        print(f"Error loading true weights from {weights_file}: {e}", file=sys.stderr)
        sys.exit(1)


def load_reference_signatures(ref_file: str) -> List[str]:
    """Load reference signature catalog (e.g., COSMIC signatures).
    
    Returns list of signature names that were available for fitting.
    """
    try:
        # Try TSV first (COSMIC format)
        df = pd.read_csv(ref_file, sep='\t', index_col='Type')
        return df.columns.tolist()
    except Exception:
        try:
            # Try CSV
            df = pd.read_csv(ref_file, sep=',')
            if 'Type' in df.columns:
                df = df.set_index('Type', drop=True)
            elif df.columns[0] in ['sub', 'tri']:  # mmsig format
                return [col for col in df.columns if col not in ['sub', 'tri', 'Type', 'SubType']]
            return df.columns.tolist()
        except Exception as e:
            print(f"Error loading reference signatures from {ref_file}: {e}", file=sys.stderr)
            sys.exit(1)


def evaluate_performance(
    df_fitted: pd.DataFrame,
    df_true: pd.DataFrame,
    ref_signatures: List[str],
    muts: pd.Series,
    min_contribution: float = DEFAULT_MIN_CONTRIBUTION,
    out_of_reference_sigs: Optional[List[str]] = None,
    cosmic_profiles: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    """Calculate evaluation metrics comparing fitted weights to ground truth.
    
    This implements the core evaluation logic from SigFitTest's eval_results() function.
    
    Args:
        df_fitted: DataFrame with fitted signature weights (signatures × samples)
        df_true: DataFrame with true signature weights (signatures × samples)
        ref_signatures: List of reference signatures available for fitting
        muts: Series with number of mutations per sample
        min_contribution: Predicted signatures with fewer mutations count as not found
            (TAE and other errors still use the full weights)
        out_of_reference_sigs: List of out-of-reference signature names to exclude from ground truth
        
    Returns:
        Dictionary of evaluation metrics
    """
    # Align signature sets using helper function
    df_fitted, df_true = align_signature_dataframes(df_fitted, df_true, ref_signatures)
    
    # Convert true weights from proportions to absolute counts
    # Following SigFitTest: true_res = (contribs.T * muts).T
    # This multiplies each column by the corresponding mutation count
    df_true = df_true * muts
    
    # ========================================================================
    # EVALUATION APPROACH: Keep all ground truth, exclude predicted non-fitting-catalog
    # ========================================================================
    # KEEP ALL signatures in ground truth (simulation catalog), including:
    # - Known signatures (also in fitting catalog): counted as TP or FN
    # - Unknown signatures (in simulation but NOT in fitting catalog): counted as FN if missed
    #   * Unknowns are determined by: (signature in truth) AND (NOT in fitting catalog)
    #   * e.g., SBS86 in simulation catalog but not in fitting catalog = unknown
    #   * Methods are penalized for missing this signal (false negative)
    # 
    # Unknown signatures: in the simulation but not in the fitting reference.
    # They stay in the evaluation, so a tool that misses their mutations is penalised.
    unknown_in_truth_mask = ~df_true.index.isin(ref_signatures)
    if unknown_in_truth_mask.any():
        unknown_sigs = df_true[unknown_in_truth_mask].index.tolist()
        unknown_weight = df_true[unknown_in_truth_mask].sum().sum()
        print(f"  Ground truth contains {unknown_in_truth_mask.sum()} UNKNOWN signature(s) (in simulation but not in fitting catalog): {unknown_sigs}", file=sys.stderr)
        print(f"  Total unknown weight: {unknown_weight:.1f} mutations", file=sys.stderr)
        print(f"  KEEPING unknowns in evaluation - methods penalized if they miss this signal", file=sys.stderr)
    
    # The out_of_reference_sigs list is not needed: unknown signatures are found above by
    # comparing the simulation with the fitting reference
    if out_of_reference_sigs and len(out_of_reference_sigs) > 0:
        print(f"  Note: out_of_reference_sigs parameter provided but not needed", file=sys.stderr)
        print(f"        Unknowns automatically identified by comparing simulation vs fitting catalog", file=sys.stderr)
    
    # Signatures a tool creates itself (e.g. Unknown_X, SigA) are not in the fitting
    # reference and cannot be matched to the truth, so they are left out of the evaluation
    non_fitting_mask = ~df_fitted.index.isin(ref_signatures)
    if non_fitting_mask.any():
        non_fitting_sigs = df_fitted[non_fitting_mask].index.tolist()
        non_fitting_weight = df_fitted[non_fitting_mask].sum().sum()
        print(f"  Predictions contain {non_fitting_mask.sum()} non-fitting-catalog signature(s): {non_fitting_sigs}", file=sys.stderr)
        print(f"  Total non-fitting-catalog weight: {non_fitting_weight:.1f} mutations", file=sys.stderr)
        print(f"  EXCLUDING non-fitting-catalog signatures from evaluation (internally created)", file=sys.stderr)
        df_fitted = df_fitted[~non_fitting_mask]
    else:
        print(f"  All {len(df_fitted)} prediction signatures are in fitting catalog", file=sys.stderr)
    
    # Re-align DataFrames after removing non-fitting-catalog predictions
    df_fitted, df_true = align_signature_dataframes(df_fitted, df_true, ref_signatures)
    
    # Calculate total absolute error for each sample
    err = df_fitted.sub(df_true, axis=0).abs().sum() / (2 * muts)
    TAE = err.mean()
    TAE_std = err.std()
    nRMSE = np.sqrt(np.power(err, 2).mean())
    
    # Initialize counters for per-sample metrics
    wtot_FP, wtot_FP_squared, num_FP_sigs, wtot_FN, num_FN_sigs, MAE_active = 0, 0, 0, 0, 0, 0
    P_vals, R_vals, S_vals, F_vals, MCC_vals = [], [], [], [], []
    
    # Calculate metrics for each sample
    per_sample_metrics = []  # Store per-sample metrics for detailed analysis
    
    for sample in df_true.columns:
        wtot_FP_one_sample = 0
        num_TP, num_TN, num_FP, num_FN = 0, 0, 0, 0
        sample_MAE_active = 0
        sample_wtot_FN = 0
        sample_num_FN = 0
        
        # Count known vs unknown signatures in this sample's ground truth.
        # "Unknown" = active in truth AND not in the fitting reference catalog.
        num_active_sigs = 0
        num_unknown_active = 0
        unknown_weight = 0.0  # mutation count on unknown signatures (for theoretical_min_TAE)
        for sig in df_true.index:
            true_weight = df_true.loc[sig, sample]
            if true_weight > 0:
                num_active_sigs += 1
                if sig not in ref_signatures:
                    num_unknown_active += 1
                    unknown_weight += true_weight

        num_known_active = num_active_sigs - num_unknown_active

        # ENP from 96-channel profiles: highest (flattest) known and unknown signature
        highest_enp_known = np.nan
        highest_enp_unknown = np.nan
        if cosmic_profiles is not None:
            known_enps, unknown_enps = [], []
            for sig in df_true.index:
                if df_true.loc[sig, sample] > 0 and sig in cosmic_profiles.columns:
                    enp = calculate_enp_from_channels(cosmic_profiles[sig].values)
                    if not np.isnan(enp):
                        if sig in ref_signatures:
                            known_enps.append(enp)
                        else:
                            unknown_enps.append(enp)
            if known_enps:
                highest_enp_known = max(known_enps)
            if unknown_enps:
                highest_enp_unknown = max(unknown_enps)

        for sig in df_true.index:
            true_weight = df_true.loc[sig, sample]
            fitted_weight = df_fitted.loc[sig, sample]
            
            if true_weight > 0:  # Active signature in ground truth
                # Calculate error for active signatures
                sample_MAE_active += np.abs(fitted_weight - true_weight) / muts[sample]
                MAE_active += np.abs(fitted_weight - true_weight) / muts[sample]
                
                if fitted_weight < min_contribution:  # False negative
                    sample_wtot_FN += true_weight / muts[sample]
                    sample_num_FN += 1
                    wtot_FN += true_weight / muts[sample]
                    num_FN_sigs += 1
                    num_FN += 1
                else:  # True positive
                    num_TP += 1
            else:  # Inactive signature in ground truth
                if fitted_weight >= min_contribution:  # False positive
                    wtot_FP += fitted_weight / muts[sample]
                    wtot_FP_one_sample += fitted_weight
                    num_FP_sigs += 1
                    num_FP += 1
                else:  # True negative
                    num_TN += 1
        
        # Normalize FP weight for this sample
        wtot_FP_one_sample /= muts[sample]
        wtot_FP_squared += wtot_FP_one_sample * wtot_FP_one_sample
        
        # Calculate TAE for this sample
        sample_TAE = err[sample]
        
        # Calculate nRMSE for this sample (RMSE of signature weight errors)
        # nRMSE = sqrt(mean((w_fitted - w_true)^2))
        sig_errors = (df_fitted[sample] - df_true[sample]) / muts[sample]
        sample_nRMSE = np.sqrt(np.mean(sig_errors ** 2))
        
        # Calculate classification metrics for this sample
        if num_TP + num_FP > 0:
            Psample = num_TP / (num_TP + num_FP)  # Precision
            Rsample = num_TP / (num_TP + num_FN) if (num_TP + num_FN) > 0 else 0  # Recall
        else:
            Psample = 0
            Rsample = 0
        
        if num_TN + num_FP > 0:
            Ssample = num_TN / (num_TN + num_FP)  # Specificity
        else:
            Ssample = 0
        
        if Psample + Rsample > 0:
            Fsample = 2 * Psample * Rsample / (Psample + Rsample)  # F1 score
        else:
            Fsample = 0
        
        # Matthews correlation coefficient
        if (Psample + Rsample > 0 and num_TN + num_FP > 0 and num_TN + num_FN > 0):
            denom = np.sqrt((num_TP + num_FP) * (num_TP + num_FN) * (num_TN + num_FP) * (num_TN + num_FN))
            MCCsample = (num_TP * num_TN - num_FP * num_FN) / denom
        else:
            MCCsample = 0
        
        P_vals.append(Psample)
        R_vals.append(Rsample)
        S_vals.append(Ssample)
        F_vals.append(Fsample)
        MCC_vals.append(MCCsample)
        
        # Calculate theoretical maximum performance for this sample
        # Best case: detect all known sigs, miss unknowns, no FP
        theoretical_max_recall = num_known_active / num_active_sigs if num_active_sigs > 0 else 1.0
        theoretical_max_precision = 1.0  # No FP in best case
        theoretical_max_f1 = 2 * theoretical_max_precision * theoretical_max_recall / (theoretical_max_precision + theoretical_max_recall) if (theoretical_max_precision + theoretical_max_recall) > 0 else 0.0
        
        # MCC calculation for best case
        # TP=num_known_active, FP=0, FN=num_unknown_active, TN=(total_ref_sigs - num_active_sigs)
        total_ref_sigs = len(df_true.index)  # Total signatures in reference
        TN_theoretical = total_ref_sigs - num_active_sigs
        TP_theoretical = num_known_active
        FP_theoretical = 0
        FN_theoretical = num_unknown_active
        
        numerator = TP_theoretical * TN_theoretical - FP_theoretical * FN_theoretical
        denominator = np.sqrt((TP_theoretical + FP_theoretical) * (TP_theoretical + FN_theoretical) * (TN_theoretical + FP_theoretical) * (TN_theoretical + FN_theoretical))
        theoretical_max_mcc = numerator / denominator if denominator > 0 else 0.0

        # Best case for TAE: known signatures fit perfectly (zero error), unknown
        # signatures' mutations can't be placed anywhere (no catalog entry for them) and
        # aren't misattributed onto known signatures either -- so the only irreducible
        # error is the unknown signatures' own true weight. Same normalization as TAE
        # itself (err.../(2*muts)). Zero when there are no unknowns (e.g. a custom
        # dataset's own catalog, with nothing excluded from the fitting reference).
        theoretical_min_tae = unknown_weight / (2 * muts[sample]) if muts[sample] > 0 else 0.0

        # Store per-sample metrics
        per_sample_metrics.append({
            'sample': sample,
            'num_muts': muts[sample],
            'TAE': sample_TAE,
            'nRMSE': sample_nRMSE,
            'P': Psample,
            'R': Rsample,
            'S': Ssample,
            'F1': Fsample,
            'MCC': MCCsample,
            'MAE_active': sample_MAE_active / 2,
            'theoretical_max_Recall': theoretical_max_recall,
            'theoretical_max_Precision': theoretical_max_precision,
            'theoretical_max_F1': theoretical_max_f1,
            'theoretical_max_MCC': theoretical_max_mcc,
            'theoretical_min_FN': num_unknown_active,
            'theoretical_min_TAE': theoretical_min_tae,
            'wtot_FP': wtot_FP_one_sample,
            'num_FP': num_FP,
            'wtot_FN': sample_wtot_FN,
            'num_FN': sample_num_FN,
            'num_TP': num_TP,
            'num_TN': num_TN,
            'highest_ENP_known': highest_enp_known,
            'highest_ENP_unknown': highest_enp_unknown,
        })
    
    # Normalize metrics across samples
    num_samples = df_fitted.shape[1]
    MAE_active /= (2 * num_samples)
    wtot_FP /= num_samples
    wtot_FP_squared /= num_samples
    wtot_FP_squared -= wtot_FP * wtot_FP
    num_FP_sigs /= num_samples
    wtot_FN /= num_samples
    num_FN_sigs /= num_samples
    
    # Calculate additional summary metrics
    weight_tot = (df_fitted.sum() / muts).mean()  # Total assigned weight (normalized)
    
    # Effective number of signatures per sample
    w_norm = df_fitted / df_fitted.sum()
    n_eff = np.ma.masked_invalid(1 / np.power(w_norm, 2).sum()).mean()
    if np.ma.is_masked(n_eff):
        n_eff = np.nan
    
    summary_metrics = {
        'num_samples': num_samples,
        'mean_num_muts': muts.mean(),
        'TAE': TAE,
        'TAE_std': TAE_std,
        'nRMSE': nRMSE,
        'weight_tot': weight_tot,
        'n_eff': n_eff,
        'MAE_active': MAE_active,
        'wtot_FP': wtot_FP,
        'wtot_FP_std': np.sqrt(wtot_FP_squared),
        'num_FP_sigs': num_FP_sigs,
        'wtot_FN': wtot_FN,
        'num_FN_sigs': num_FN_sigs,
        'P_mean': np.mean(P_vals),
        'P_std': np.std(P_vals),
        'P_sem': np.std(P_vals) / np.sqrt(num_samples) if num_samples > 0 else 0,
        'R_mean': np.mean(R_vals),
        'R_std': np.std(R_vals),
        'R_sem': np.std(R_vals) / np.sqrt(num_samples) if num_samples > 0 else 0,
        'S_mean': np.mean(S_vals),
        'S_std': np.std(S_vals),
        'S_sem': np.std(S_vals) / np.sqrt(num_samples) if num_samples > 0 else 0,
        'F1_mean': np.mean(F_vals),
        'F1_std': np.std(F_vals),
        'F1_sem': np.std(F_vals) / np.sqrt(num_samples) if num_samples > 0 else 0,
        'MCC_mean': np.mean(MCC_vals),
        'MCC_std': np.std(MCC_vals),
        'MCC_sem': np.std(MCC_vals) / np.sqrt(num_samples) if num_samples > 0 else 0
    }
    
    return summary_metrics, per_sample_metrics


def create_comparison_plots(all_metrics, plots_dir, data_type, num_muts, relevance):
    """Create comparison plots across all tools.
    
    Args:
        all_metrics: List of dicts, one per tool with metrics
        plots_dir: Directory to save plots
        data_type: Type of data (synthetic, realsim)
        num_muts: Number of mutations
        relevance: Signature set relevance (standard, relevant)
    """
    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)
    
    # Set style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(all_metrics)
    df = df.sort_values('tool')
    
    # Create filename prefix
    prefix = f"{data_type}_{num_muts}_{relevance}"
    
    # 1. Bar chart of key metrics
    create_metrics_bar_chart(df, plots_path, prefix)
    
    # 2. Precision-Recall scatter
    create_precision_recall_plot(df, plots_path, prefix)
    
    # 3. Error metrics comparison
    create_error_comparison(df, plots_path, prefix)
    
    # 4. False positive/negative analysis
    create_fp_fn_plot(df, plots_path, prefix)
    
    print(f"  Plots saved to {plots_path}", file=sys.stderr)


def create_metrics_bar_chart(df, plots_path, prefix):
    """Create grouped bar chart of key metrics."""
    metrics = ['TAE', 'F1_mean', 'P_mean', 'R_mean', 'MCC_mean']
    metric_labels = ['TAE\n(lower better)', 'F1 Score', 'Precision', 'Recall', 'MCC']
    
    fig, axes = plt.subplots(1, len(metrics), figsize=(16, 4))
    
    colors = sns.color_palette("husl", len(df))
    
    for idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[idx]
        
        # For TAE, invert so lower is better visually
        if metric == 'TAE':
            values = df[metric].values
            bars = ax.bar(range(len(df)), values, color=colors)
            ax.set_ylabel(label)
            # Add horizontal line at 0.05 (typical good performance threshold)
            ax.axhline(y=0.05, color='red', linestyle='--', alpha=0.5, linewidth=1)
        else:
            values = df[metric].values
            bars = ax.bar(range(len(df)), values, color=colors)
            ax.set_ylabel(label)
            ax.set_ylim([0, 1.05])
            # Add horizontal line at 0.9 (typical good performance threshold)
            ax.axhline(y=0.9, color='green', linestyle='--', alpha=0.5, linewidth=1)
        
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df['tool'].values, rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}', ha='center', va='bottom', fontsize=8)
    
    plt.suptitle(f'Tool Performance Comparison - {prefix}', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(plots_path / f'{prefix}_metrics_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(plots_path / f'{prefix}_metrics_comparison.pdf', bbox_inches='tight')
    plt.close()
    print(f"    Created: {prefix}_metrics_comparison.png/pdf", file=sys.stderr)


def create_precision_recall_plot(df, plots_path, prefix):
    """Create precision-recall scatter plot."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Scatter plot with tool names
    colors = sns.color_palette("husl", len(df))
    for idx, row in df.iterrows():
        ax.scatter(row['R_mean'], row['P_mean'], s=200, alpha=0.7, 
                  c=[colors[idx]], edgecolors='black', linewidth=1.5)
        ax.text(row['R_mean'], row['P_mean'], row['tool'], 
               ha='center', va='center', fontsize=9, fontweight='bold')
    
    # Add F1 iso-contours
    recall_vals = np.linspace(0, 1, 100)
    for f1 in [0.5, 0.7, 0.9]:
        precision_vals = (f1 * recall_vals) / (2 * recall_vals - f1 + 1e-10)
        precision_vals = np.clip(precision_vals, 0, 1)
        ax.plot(recall_vals, precision_vals, '--', alpha=0.3, 
               label=f'F1={f1}', color='gray')
    
    ax.set_xlabel('Recall (Sensitivity)', fontsize=12)
    ax.set_ylabel('Precision', fontsize=12)
    ax.set_title(f'Precision-Recall Trade-off - {prefix}', fontsize=14)
    ax.set_xlim([0, 1.05])
    ax.set_ylim([0, 1.05])
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower left')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.2)  # Diagonal reference
    
    plt.tight_layout()
    plt.savefig(plots_path / f'{prefix}_precision_recall.png', dpi=300, bbox_inches='tight')
    plt.savefig(plots_path / f'{prefix}_precision_recall.pdf', bbox_inches='tight')
    plt.close()
    print(f"    Created: {prefix}_precision_recall.png/pdf", file=sys.stderr)


def create_error_comparison(df, plots_path, prefix):
    """Create comparison of error metrics."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    colors = sns.color_palette("husl", len(df))
    
    # TAE with error bars
    ax = axes[0]
    x_pos = range(len(df))
    ax.bar(x_pos, df['TAE'], yerr=df['TAE_std'], capsize=5, color=colors, alpha=0.7)
    ax.set_ylabel('TAE (Total Absolute Error)')
    ax.set_xlabel('Tool')
    ax.set_title('Fitting Error')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['tool'].values, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # MAE active signatures
    ax = axes[1]
    ax.bar(x_pos, df['MAE_active'], color=colors, alpha=0.7)
    ax.set_ylabel('MAE (Active Signatures Only)')
    ax.set_xlabel('Tool')
    ax.set_title('Error on Active Signatures')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['tool'].values, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # nRMSE
    ax = axes[2]
    ax.bar(x_pos, df['nRMSE'], color=colors, alpha=0.7)
    ax.set_ylabel('nRMSE (Normalized RMSE)')
    ax.set_xlabel('Tool')
    ax.set_title('Normalized Root Mean Square Error')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['tool'].values, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle(f'Error Metrics Comparison - {prefix}', fontsize=14, y=1.00)
    plt.tight_layout()
    plt.savefig(plots_path / f'{prefix}_error_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig(plots_path / f'{prefix}_error_metrics.pdf', bbox_inches='tight')
    plt.close()
    print(f"    Created: {prefix}_error_metrics.png/pdf", file=sys.stderr)


def create_fp_fn_plot(df, plots_path, prefix):
    """Create false positive/negative analysis plot."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    colors = sns.color_palette("husl", len(df))
    x_pos = range(len(df))
    
    # False positives
    ax = axes[0]
    ax.bar(x_pos, df['num_FP_sigs'], alpha=0.7, color=colors, label='Number of FP sigs')
    ax.set_ylabel('Count / Weight', color='C0')
    ax.set_xlabel('Tool')
    ax.set_title('False Positive Analysis')
    ax.tick_params(axis='y', labelcolor='C0')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['tool'].values, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add FP weight on secondary axis
    ax2 = ax.twinx()
    ax2.plot(x_pos, df['wtot_FP'], 'ro-', linewidth=2, markersize=8, label='Total FP weight')
    ax2.set_ylabel('FP Weight', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    # False negatives
    ax = axes[1]
    ax.bar(x_pos, df['num_FN_sigs'], alpha=0.7, color=colors, label='Number of FN sigs')
    ax.set_ylabel('Count / Weight', color='C0')
    ax.set_xlabel('Tool')
    ax.set_title('False Negative Analysis')
    ax.tick_params(axis='y', labelcolor='C0')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['tool'].values, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add FN weight on secondary axis
    ax2 = ax.twinx()
    ax2.plot(x_pos, df['wtot_FN'], 'ro-', linewidth=2, markersize=8, label='Total FN weight')
    ax2.set_ylabel('FN Weight', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    plt.suptitle(f'False Positive/Negative Analysis - {prefix}', fontsize=14, y=1.00)
    plt.tight_layout()
    plt.savefig(plots_path / f'{prefix}_fp_fn_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(plots_path / f'{prefix}_fp_fn_analysis.pdf', bbox_inches='tight')
    plt.close()
    print(f"    Created: {prefix}_fp_fn_analysis.png/pdf", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate signature fitting performance against ground truth',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--tool-results', required=True, nargs='+',
                        help='Path(s) to tool results file(s) (fitted signature weights)')
    parser.add_argument('--true-weights', required=True,
                        help='Path to ground truth weights file')
    parser.add_argument('--reference-signatures', required=True, nargs='+',
                        help='Path(s) to reference signature catalog(s) (COSMIC)')
    parser.add_argument('--output', required=True,
                        help='Path to output metrics file')
    parser.add_argument('--tool', 
                        help='Name of the tool being evaluated (single-tool mode)')
    parser.add_argument('--tools', nargs='+',
                        help='Names of tools being evaluated (multi-tool mode)')
    parser.add_argument('--num-muts', required=True, type=int,
                        help='Number of mutations in this simulation')
    parser.add_argument('--min-contribution', type=float, default=10.0,
                        help='Minimum mutation count to consider a signature active')
    parser.add_argument('--out-of-reference-sigs', nargs='*', default=[],
                        help='List of out-of-reference signature names to exclude from ground truth evaluation')
    parser.add_argument('--realization', type=str, default='0',
                        help='Realization/replicate number (used to track which simulation run this is from)')
    parser.add_argument('--cancer-type', type=str, default=None,
                        help='Cancer type for this evaluation (e.g., Breast, Lung)')
    parser.add_argument('--data-type', type=str, default=None,
                        help='Data type: synthetic or realsim (used in output and plot titles)')
    parser.add_argument('--plots-dir',
                        help='Directory to save comparison plots (multi-tool mode)')
    parser.add_argument('--relevance', default='standard',
                        help='Relevance mode for plot titles (standard, relevant)')
    parser.add_argument('--cosmic-file', default=None,
                        help='Path to full COSMIC signature file for ENP computation')
    args = parser.parse_args()
    
    # Determine mode: single tool or multiple tools
    if args.tool:
        # Single-tool mode (backward compatibility)
        tools = [args.tool]
        if len(args.tool_results) != 1:
            print("Error: Single-tool mode requires exactly one --tool-results file", file=sys.stderr)
            sys.exit(1)
        if len(args.reference_signatures) != 1:
            print("Error: Single-tool mode requires exactly one --reference-signatures file", file=sys.stderr)
            sys.exit(1)
    elif args.tools:
        # Multi-tool mode
        tools = args.tools
        if len(args.tool_results) != len(tools):
            print(f"Error: Number of --tool-results ({len(args.tool_results)}) must match number of --tools ({len(tools)})", file=sys.stderr)
            sys.exit(1)
        if len(args.reference_signatures) != len(tools):
            print(f"Error: Number of --reference-signatures ({len(args.reference_signatures)}) must match number of --tools ({len(tools)})", file=sys.stderr)
            sys.exit(1)
    else:
        print("Error: Must specify either --tool (single-tool mode) or --tools (multi-tool mode)", file=sys.stderr)
        sys.exit(1)
    
    # Load COSMIC profiles for ENP computation (optional)
    cosmic_profiles = load_cosmic_profiles(args.cosmic_file) if args.cosmic_file else None

    # Load true weights (same for all tools)
    print(f"Loading true weights from {args.true_weights}...", file=sys.stderr)
    df_true = load_true_weights(args.true_weights)
    
    # Get mutation counts per sample
    num_samples = df_true.shape[1]
    muts = pd.Series([args.num_muts] * num_samples, index=df_true.columns)
    
    # Evaluate each tool
    all_metrics = []
    all_per_sample_metrics = []
    
    for idx, tool_name in enumerate(tools):
        print(f"\nEvaluating {tool_name}...", file=sys.stderr)

        # A failed tool run has no fit to score: it adds no rows (the report lists it as failed)
        reason = tool_failure(args.tool_results[idx])
        if reason:
            print(f"  {tool_name} failed in this run, nothing to evaluate: {reason}", file=sys.stderr)
            continue

        # Load tool-specific data
        print(f"  Loading results from {args.tool_results[idx]}...", file=sys.stderr)
        df_fitted = load_tool_results(args.tool_results[idx])
        
        print(f"  Loading reference signatures from {args.reference_signatures[idx]}...", file=sys.stderr)
        ref_sigs = load_reference_signatures(args.reference_signatures[idx])
        
        # Evaluate
        metrics, per_sample = evaluate_performance(df_fitted, df_true, ref_sigs, muts, args.min_contribution, args.out_of_reference_sigs, cosmic_profiles)
        metrics['tool'] = tool_name
        all_metrics.append(metrics)
        
        # Attach tool name and per-run constants to each per-sample row
        for sample_metric in per_sample:
            sample_metric['tool'] = tool_name
            sample_metric['weight_tot'] = metrics['weight_tot']
            sample_metric['n_eff'] = metrics['n_eff']
        all_per_sample_metrics.extend(per_sample)
        
        # Print summary
        print(f"  TAE: {metrics['TAE']:.4f}", file=sys.stderr)
        print(f"  Precision: {metrics['P_mean']:.4f}", file=sys.stderr)
        print(f"  Recall: {metrics['R_mean']:.4f}", file=sys.stderr)
        print(f"  F1: {metrics['F1_mean']:.4f}", file=sys.stderr)
        print(f"  MCC: {metrics['MCC_mean']:.4f}", file=sys.stderr)
    
    # Write per-sample metrics (single output — all metadata already present as columns)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cancer_type_str = args.cancer_type if args.cancer_type is not None else 'None'
    data_type_str = args.data_type if args.data_type is not None else 'None'

    def _fmt(v):
        return f"{v:.4f}" if (v is not None and not (isinstance(v, float) and np.isnan(v))) else "nan"

    with open(output_path, 'w') as f:
        f.write('cancer_type\tdata_type\ttool\treplicate\tsample\tnum_muts\t')
        f.write('TAE\tnRMSE\tweight_tot\tn_eff\t')
        f.write('P\tR\tS\tF1\tMCC\t')
        f.write('MAE_active\twtot_FP\tnum_FP\twtot_FN\tnum_FN\tnum_TP\tnum_TN\t')
        f.write('theoretical_max_Recall\ttheoretical_max_Precision\ttheoretical_max_F1\ttheoretical_max_MCC\ttheoretical_min_FN\ttheoretical_min_TAE\t')
        f.write('highest_ENP_known\thighest_ENP_unknown\n')

        for sm in all_per_sample_metrics:
            f.write(f"{cancer_type_str}\t{data_type_str}\t{sm['tool']}\t{args.realization}\t{sm['sample']}\t{sm['num_muts']}\t")
            f.write(f"{sm['TAE']:.4f}\t{sm['nRMSE']:.4f}\t{sm['weight_tot']:.4f}\t{sm['n_eff']:.4f}\t")
            f.write(f"{sm['P']:.4f}\t{sm['R']:.4f}\t{sm['S']:.4f}\t{sm['F1']:.4f}\t{sm['MCC']:.4f}\t")
            f.write(f"{sm['MAE_active']:.4f}\t{sm['wtot_FP']:.4f}\t{sm['num_FP']}\t{sm['wtot_FN']:.4f}\t{sm['num_FN']}\t")
            f.write(f"{sm['num_TP']}\t{sm['num_TN']}\t")
            f.write(f"{sm['theoretical_max_Recall']:.4f}\t{sm['theoretical_max_Precision']:.4f}\t")
            f.write(f"{sm['theoretical_max_F1']:.4f}\t{sm['theoretical_max_MCC']:.4f}\t{sm['theoretical_min_FN']}\t{sm['theoretical_min_TAE']:.4f}\t")
            f.write(f"{_fmt(sm['highest_ENP_known'])}\t{_fmt(sm['highest_ENP_unknown'])}\n")

    print(f"\nEvaluation complete. Per-sample metrics written to {args.output}", file=sys.stderr)
    
    # Create comparison plots if requested (multi-tool mode with 2+ tools)
    if args.plots_dir and len(all_metrics) >= 2:
        print(f"\nCreating comparison plots...", file=sys.stderr)
        create_comparison_plots(
            all_metrics, 
            args.plots_dir,
            args.data_type or 'dataset',
            args.num_muts,
            args.relevance
        )
    elif args.plots_dir and len(all_metrics) < 2:
        print(f"\nSkipping plots: Need at least 2 tools for comparison (got {len(all_metrics)})", file=sys.stderr)


if __name__ == '__main__':
    main()
