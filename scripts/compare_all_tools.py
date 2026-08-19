#!/usr/bin/env python3
"""
Create comprehensive comparison across all tools, datasets, and simulation conditions.

This script combines all evaluation summaries (across mutation counts, data types,
and signature sets) to create a complete performance comparison as in the SigFitTest paper.

Usage:
    python compare_all_tools.py \\
        --input-files <summary1> <summary2> ... \\
        --output <comparison_file> \\
        --figures-dir <figures_directory>
"""

import argparse
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare all tools across all simulation conditions"
    )
    parser.add_argument(
        "--input-files",
        nargs="+",
        required=True,
        help="Paths to summary files from aggregate_evaluation"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output comparison file"
    )
    parser.add_argument(
        "--figures-dir",
        required=True,
        help="Directory to save comparison figures"
    )
    return parser.parse_args()


def parse_filename(filepath):
    """Extract metadata from summary filename.
    
    Expected format: {data_type}_{num_muts}_{relevance}_summary.dat
    
    Returns:
        dict with 'data_type', 'num_muts', 'relevance'
    """
    filename = Path(filepath).stem
    
    # Remove _summary suffix
    if filename.endswith('_summary'):
        filename = filename[:-8]
    
    parts = filename.split('_')
    
    if len(parts) >= 3:
        return {
            'data_type': parts[0],
            'num_muts': parts[1],
            'relevance': parts[2]
        }
    else:
        raise ValueError(f"Could not parse filename: {filepath}")


def load_all_summaries(input_files):
    """Load all summary files and add metadata columns.
    
    Args:
        input_files: List of paths to summary files
        
    Returns:
        DataFrame with all summaries combined
    """
    all_data = []
    
    for filepath in input_files:
        try:
            # Parse filename to get metadata
            metadata = parse_filename(filepath)
            
            # Read the summary
            df = pd.read_csv(filepath, sep='\t')
            
            # Add metadata columns
            for key, value in metadata.items():
                df[key] = value
            
            all_data.append(df)
            
        except Exception as e:
            print(f"Warning: Failed to read {filepath}: {e}", file=sys.stderr)
            continue
    
    if not all_data:
        raise ValueError("No valid summary files were loaded")
    
    # Combine all data
    combined = pd.concat(all_data, ignore_index=True)
    
    # Convert num_muts to integer for proper sorting
    combined['num_muts'] = pd.to_numeric(combined['num_muts'], errors='coerce')
    
    return combined


def compute_rankings(df):
    """Compute tool rankings for each metric and condition.
    
    Args:
        df: DataFrame with all evaluation data
        
    Returns:
        DataFrame with rankings added
    """
    # Metrics where higher is better
    higher_better = ['Precision', 'Recall', 'F1', 'MCC', 'Specificity']
    # Metrics where lower is better
    lower_better = ['TAE', 'TAE_std', 'nRMSE', 'MAE_active', 'wtot_FP', 'num_FP_sigs', 'wtot_FN', 'num_FN_sigs']
    
    df_ranked = df.copy()
    
    # For each combination of data_type, num_muts, relevance
    for (data_type, num_muts, relevance), group in df.groupby(['data_type', 'num_muts', 'relevance']):
        # Get the index of this group
        idx = group.index
        
        # Rank each metric
        for metric in higher_better:
            if metric in df.columns:
                # Higher is better: rank 1 = highest value
                df_ranked.loc[idx, f'{metric}_rank'] = group[metric].rank(ascending=False, method='min')
        
        for metric in lower_better:
            if metric in df.columns:
                # Lower is better: rank 1 = lowest value
                df_ranked.loc[idx, f'{metric}_rank'] = group[metric].rank(ascending=True, method='min')
    
    return df_ranked


def create_summary_table(df):
    """Create summary statistics across all conditions.
    
    Args:
        df: DataFrame with rankings
        
    Returns:
        DataFrame with summary statistics per tool
    """
    # Identify ranking columns
    rank_cols = [col for col in df.columns if col.endswith('_rank')]
    
    if not rank_cols:
        print("Warning: No ranking columns found", file=sys.stderr)
        return None
    
    # Group by tool and compute mean rank for each metric
    summary = df.groupby('tool')[rank_cols].agg(['mean', 'median', 'std'])
    
    # Also compute overall mean rank across all metrics
    df['overall_rank'] = df[rank_cols].mean(axis=1)
    overall_summary = df.groupby('tool')['overall_rank'].agg(['mean', 'median', 'std'])
    overall_summary.columns = ['overall_mean_rank', 'overall_median_rank', 'overall_std_rank']
    
    # Combine
    full_summary = pd.concat([summary, overall_summary], axis=1)
    
    # Sort by overall mean rank
    full_summary = full_summary.sort_values('overall_mean_rank')
    
    return full_summary


def save_figures(df, figures_dir):
    """Create and save comparison figures reproducing SigFitTest publication plots.
    
    Args:
        df: DataFrame with all evaluation data
        figures_dir: Directory to save figures
    """
    figures_path = Path(figures_dir)
    figures_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating figures in: {figures_path}")
    
    # Set publication-quality style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
    
    # 1. Main comparison figure: TAE and F1 across mutation counts and data types
    # Similar to Figure 2 in the publication
    create_main_comparison_figure(df, figures_path)
    
    # 2. Performance metrics heatmap
    create_performance_heatmap(df, figures_path)
    
    # 3. Box plots for key metrics
    create_metric_boxplots(df, figures_path)
    
    # 4. Ranking visualization
    create_ranking_plot(df, figures_path)
    
    # 5. Performance vs mutation count trends
    create_mutation_trend_plots(df, figures_path)
    
    print(f"Figures saved to {figures_path}")


def create_main_comparison_figure(df, figures_path):
    """Create main comparison figure showing TAE and F1 across conditions.
    
    Similar to update_results-set6.png from the publication.
    """
    # Filter for main metrics
    metrics_to_plot = ['TAE', 'F1', 'Precision', 'Recall']
    
    # Get unique tools, mutation counts, data types
    tools = sorted(df['tool'].unique())
    num_muts_list = sorted(df['num_muts'].unique())
    data_types = sorted(df['data_type'].unique())
    
    # Create subplots: 2 rows (TAE, F1) x number of data types
    n_data_types = len(data_types)
    fig, axes = plt.subplots(2, n_data_types, figsize=(5 * n_data_types, 10))
    
    if n_data_types == 1:
        axes = axes.reshape(-1, 1)
    
    for col_idx, data_type in enumerate(data_types):
        df_subset = df[df['data_type'] == data_type]
        
        # TAE plot (row 0)
        ax_tae = axes[0, col_idx]
        for tool in tools:
            tool_data = df_subset[df_subset['tool'] == tool]
            if not tool_data.empty:
                ax_tae.plot(tool_data['num_muts'], tool_data['TAE'], 
                           marker='o', label=tool, linewidth=2)
        
        ax_tae.set_xscale('log')
        ax_tae.set_xlabel('Number of mutations')
        ax_tae.set_ylabel('TAE (Total Absolute Error)')
        ax_tae.set_title(f'{data_type} - TAE')
        ax_tae.grid(True, alpha=0.3)
        if col_idx == n_data_types - 1:
            ax_tae.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # F1 plot (row 1)
        ax_f1 = axes[1, col_idx]
        for tool in tools:
            tool_data = df_subset[df_subset['tool'] == tool]
            if not tool_data.empty:
                ax_f1.plot(tool_data['num_muts'], tool_data['F1'], 
                          marker='o', label=tool, linewidth=2)
        
        ax_f1.set_xscale('log')
        ax_f1.set_xlabel('Number of mutations')
        ax_f1.set_ylabel('F1 Score')
        ax_f1.set_title(f'{data_type} - F1 Score')
        ax_f1.set_ylim([0, 1.05])
        ax_f1.grid(True, alpha=0.3)
        if col_idx == n_data_types - 1:
            ax_f1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(figures_path / 'main_comparison_TAE_F1.png', dpi=300, bbox_inches='tight')
    plt.savefig(figures_path / 'main_comparison_TAE_F1.pdf', bbox_inches='tight')
    plt.close()
    print("  Created: main_comparison_TAE_F1.png/pdf")


def create_performance_heatmap(df, figures_path):
    """Create heatmap of tool performance across conditions."""
    # Pivot data for heatmap: tools x conditions
    # Use TAE as the main metric
    
    # Create condition labels
    df['condition'] = df['data_type'] + '_' + df['num_muts'].astype(str) + '_' + df['relevance']
    
    # Pivot
    heatmap_data = df.pivot_table(
        index='tool',
        columns='condition',
        values='TAE',
        aggfunc='mean'
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(12, len(heatmap_data.columns) * 0.8), 
                                     max(6, len(heatmap_data.index) * 0.5)))
    
    # Create heatmap (use reversed colormap so darker = better/lower TAE)
    sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='RdYlGn_r', 
                cbar_kws={'label': 'TAE'}, ax=ax, linewidths=0.5)
    
    ax.set_xlabel('Condition (data_type_mutations_relevance)')
    ax.set_ylabel('Tool')
    ax.set_title('Tool Performance Heatmap (TAE across conditions)')
    
    plt.tight_layout()
    plt.savefig(figures_path / 'performance_heatmap_TAE.png', dpi=300, bbox_inches='tight')
    plt.savefig(figures_path / 'performance_heatmap_TAE.pdf', bbox_inches='tight')
    plt.close()
    print("  Created: performance_heatmap_TAE.png/pdf")
    
    # Also create F1 heatmap
    heatmap_data_f1 = df.pivot_table(
        index='tool',
        columns='condition',
        values='F1',
        aggfunc='mean'
    )
    
    fig, ax = plt.subplots(figsize=(max(12, len(heatmap_data_f1.columns) * 0.8), 
                                     max(6, len(heatmap_data_f1.index) * 0.5)))
    
    sns.heatmap(heatmap_data_f1, annot=True, fmt='.3f', cmap='RdYlGn', 
                cbar_kws={'label': 'F1 Score'}, ax=ax, linewidths=0.5)
    
    ax.set_xlabel('Condition (data_type_mutations_relevance)')
    ax.set_ylabel('Tool')
    ax.set_title('Tool Performance Heatmap (F1 Score across conditions)')
    
    plt.tight_layout()
    plt.savefig(figures_path / 'performance_heatmap_F1.png', dpi=300, bbox_inches='tight')
    plt.savefig(figures_path / 'performance_heatmap_F1.pdf', bbox_inches='tight')
    plt.close()
    print("  Created: performance_heatmap_F1.png/pdf")


def create_metric_boxplots(df, figures_path):
    """Create box plots comparing tools across key metrics."""
    metrics = ['TAE', 'F1', 'Precision', 'Recall', 'MCC']
    
    # Create subplot for each metric
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        
        # Box plot
        tools = sorted(df['tool'].unique())
        data_by_tool = [df[df['tool'] == tool][metric].dropna() for tool in tools]
        
        bp = ax.boxplot(data_by_tool, labels=tools, patch_artist=True)
        
        # Color boxes
        colors = sns.color_palette("husl", len(tools))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        
        ax.set_ylabel(metric)
        ax.set_xlabel('Tool')
        ax.set_title(f'{metric} Distribution Across All Conditions')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Remove empty subplot
    fig.delaxes(axes[-1])
    
    plt.tight_layout()
    plt.savefig(figures_path / 'metric_boxplots.png', dpi=300, bbox_inches='tight')
    plt.savefig(figures_path / 'metric_boxplots.pdf', bbox_inches='tight')
    plt.close()
    print("  Created: metric_boxplots.png/pdf")


def create_ranking_plot(df, figures_path):
    """Create visualization of tool rankings across metrics."""
    # Get ranking columns
    rank_cols = [col for col in df.columns if col.endswith('_rank')]
    
    if not rank_cols:
        print("  Skipping ranking plot: no ranking columns found")
        return
    
    # Calculate mean rank per tool for each metric
    ranking_data = df.groupby('tool')[rank_cols].mean()
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, len(ranking_data.index) * 0.4)))
    
    # Create horizontal bar plot
    y_pos = np.arange(len(ranking_data.index))
    mean_ranks = ranking_data.mean(axis=1).sort_values()
    
    bars = ax.barh(y_pos, mean_ranks.values)
    
    # Color bars by rank
    colors = plt.cm.RdYlGn_r(mean_ranks.values / mean_ranks.max())
    for bar, color in zip(bars, colors):
        bar.set_color(color)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(mean_ranks.index)
    ax.set_xlabel('Mean Rank (lower is better)')
    ax.set_title('Overall Tool Rankings Across All Metrics and Conditions')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(figures_path / 'tool_rankings.png', dpi=300, bbox_inches='tight')
    plt.savefig(figures_path / 'tool_rankings.pdf', bbox_inches='tight')
    plt.close()
    print("  Created: tool_rankings.png/pdf")


def create_mutation_trend_plots(df, figures_path):
    """Create plots showing how performance changes with mutation count."""
    metrics = ['TAE', 'F1', 'Precision', 'Recall']
    
    # Get unique data types
    data_types = sorted(df['data_type'].unique())
    
    for data_type in data_types:
        df_subset = df[df['data_type'] == data_type]
        
        # Create subplot for each metric
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
        
        for idx, metric in enumerate(metrics):
            ax = axes[idx]
            
            # Plot each tool
            tools = sorted(df_subset['tool'].unique())
            for tool in tools:
                tool_data = df_subset[df_subset['tool'] == tool]
                if not tool_data.empty:
                    # Group by num_muts and calculate mean and std
                    grouped = tool_data.groupby('num_muts')[metric].agg(['mean', 'std'])
                    
                    ax.plot(grouped.index, grouped['mean'], marker='o', 
                           label=tool, linewidth=2)
                    
                    # Add error bars if std is available
                    if not grouped['std'].isna().all():
                        ax.fill_between(grouped.index, 
                                       grouped['mean'] - grouped['std'],
                                       grouped['mean'] + grouped['std'],
                                       alpha=0.2)
            
            ax.set_xscale('log')
            ax.set_xlabel('Number of mutations')
            ax.set_ylabel(metric)
            ax.set_title(f'{metric} vs Mutation Count')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)
        
        plt.suptitle(f'Performance Trends - {data_type}', fontsize=14, y=1.00)
        plt.tight_layout()
        
        safe_data_type = data_type.replace('/', '_')
        plt.savefig(figures_path / f'mutation_trends_{safe_data_type}.png', 
                   dpi=300, bbox_inches='tight')
        plt.savefig(figures_path / f'mutation_trends_{safe_data_type}.pdf', 
                   bbox_inches='tight')
        plt.close()
        print(f"  Created: mutation_trends_{safe_data_type}.png/pdf")


def main():
    """Main entry point."""
    args = parse_args()
    
    print(f"Comparing results from {len(args.input_files)} summary files...")
    
    # Load all summaries
    all_data = load_all_summaries(args.input_files)
    print(f"Loaded data for {all_data['tool'].nunique()} tools across {len(all_data)} conditions")
    
    # Compute rankings
    ranked_data = compute_rankings(all_data)
    
    # Create summary table
    summary_table = create_summary_table(ranked_data)
    
    # Save detailed comparison (all data with rankings)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ranked_data.to_csv(output_path, sep='\t', index=False)
    print(f"Detailed comparison saved to {output_path}")
    
    # Save summary table
    if summary_table is not None:
        summary_path = output_path.parent / "tool_ranking_summary.dat"
        summary_table.to_csv(summary_path, sep='\t')
        print(f"Ranking summary saved to {summary_path}")
        
        # Print top tools
        print("\nTop 5 tools by overall mean rank:")
        print(summary_table.head()[['overall_mean_rank', 'overall_median_rank']].to_string())
    
    # Create figures
    save_figures(ranked_data, args.figures_dir)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
