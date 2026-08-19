#!/usr/bin/env python3
"""
Generate a clean HTML benchmark report for mutational signature fitting tools.

This script creates an interactive HTML report with:
- Summary statistics
- Interactive Plotly visualizations  
- Per-cancer-type analysis
- Aggregate performance metrics

Author: SigFitTest Workflow
"""

import argparse
import os
import re
from pathlib import Path
from datetime import datetime
import yaml

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ============================================================================
# Configuration Classes - Centralize all tool and column configurations
# ============================================================================

class ReportMetadata:
    """Centralize all tool colors, symbols, column mappings, and display order."""
    
    def __init__(self):
        # Tool colors — same-family methods share similar hues, families are
        # maximally distinct from each other.
        self.tool_colors = {
            # MuSiCal family — blues (dark → light)
            'MuSiCal':                 '#0C5DA5',
            'MuSiCal-FN':              '#4898CC',
            'MuSiCal-BI':              '#87C0E0',
            # SPA family — amber/orange (distinct warm tone)
            'SPA':                     '#E69F00',
            'SPA_low':                 '#B36B00',
            'FSPA':                    '#FFCC66',
            # SigID — teal (auto-penalty unknown-signature discovery)
            'SigID':                   '#007A73',
            # SigTrim — red (unassigned-mutation trimming)
            'SigTrim':                 '#CC3311',
            # external tools — purples/neutrals/pinks, kept distinct from above
            'deconstructSigs':         '#8B4FCF',
            'MutationalPatterns':      '#5E5E5E',
            'sigfit':                  '#B8860B',
            'yapsa':                   '#6A5ACD',
            'sigspack':                '#A0522D',
            'signature_tools':         '#9932CC',
            'siglasso':                '#556B2F',
            'mmsig':                   '#C2185B',
            'sigminer':                '#455A64',
        }

        # Base method extraction (for grouping tools by approach)
        self.base_method_map = {
            'MuSiCal': 'MuSiCal', 'MuSiCal-FN': 'MuSiCal', 'MuSiCal-BI': 'MuSiCal',
            'SPA': 'SPA', 'SPA_low': 'SPA', 'FSPA': 'SPA',
            'SigID': 'SigID',
            'SigTrim': 'SigTrim',
        }

        # Symbols for different approaches
        self.method_symbols = {
            'MuSiCal': 'circle',
            'SPA': 'square',
            'SigID': 'diamond',
            'SigTrim': 'triangle-up',
        }

        # Display order for consistent presentation
        self.tool_display_order = [
            'SPA', 'SPA_low', 'FSPA',
            'MuSiCal', 'MuSiCal-FN', 'MuSiCal-BI',
            'SigID', 'SigTrim',
            'deconstructSigs', 'MutationalPatterns',
            'sigfit', 'yapsa', 'sigspack',
            'signature_tools', 'siglasso',
            'mmsig', 'sigminer'
        ]
        
        # Column name mappings (for flexibility in input data)
        self.column_aliases = {
            'precision': ['precision', 'P'],
            'recall': ['recall', 'R'],
            'f1_score': ['f1_score', 'F1'],
            'mcc': ['mcc', 'MCC'],
            'tae': ['tae', 'TAE'],
            'nrmse': ['nrmse', 'nRMSE'],
            'fp': ['fp', 'FP'],
            'fn': ['fn', 'FN']
        }
    
    def _find_col(self, df, aliases):
        """Find the first matching column name from a list of aliases."""
        for alias in aliases:
            if alias in df.columns:
                return alias
        return None
    
    def get_base_method(self, tool):
        """Get the base method for a tool."""
        return self.base_method_map.get(tool, tool)
    
    def get_ordered_tools(self, available_tools):
        """Get tools in display order, filtering to only those available."""
        ordered = [t for t in self.tool_display_order if t in available_tools]
        remaining = sorted(set(available_tools) - set(ordered))
        return ordered + remaining


class CancerData:
    """Pre-filtered dataset for a specific cancer type with computed ranges."""
    
    def __init__(self, df, df_samples, cancer_type, metadata):
        self.cancer_type = cancer_type
        self.metadata = metadata
        
        # Filter dataframes
        self.df = df[df['cancer_type'] == cancer_type].copy()
        self.df_samples = None
        if df_samples is not None:
            self.df_samples = df_samples[df_samples['cancer_type'] == cancer_type].copy()
        
        # Extract key info
        self.tools = metadata.get_ordered_tools(set(self.df['tool'].unique()))
        self.mutation_burdens = sorted(self.df['num_muts'].unique())
        
        # Pre-compute ranges for consistent axis scaling
        self.ranges = self._compute_ranges()
    
    def _compute_ranges(self):
        """Compute value ranges for consistent plotting."""
        ranges = {}
        
        # Precision and Recall (0-1.1 with padding)
        if 'precision' in self.df.columns:
            ranges['precision'] = (0, 1.1)
        if 'recall' in self.df.columns:
            ranges['recall'] = (0, 1.1)
        
        # Other metrics - compute from data
        for col in ['f1_score', 'tae', 'mcc', 'nrmse']:
            if col in self.df.columns:
                vals = self.df[col].dropna()
                if len(vals) > 0:
                    vmin, vmax = vals.min(), vals.max()
                    padding = (vmax - vmin) * 0.1
                    ranges[col] = (vmin - padding, vmax + padding)
        
        return ranges


# ============================================================================
# Helper Functions for Penalty Display
# ============================================================================

def get_hover_text_with_penalty(df, tool, text_items=None):
    """
    Generate hover text that includes penalty information for penalty methods.
    
    Args:
        df: DataFrame containing the data
        tool: Tool name
        text_items: List of additional text items to include (e.g., ["Mutations: 10000"])
    
    Returns:
        List of hover text strings, one per data point
    """
    tool_data = df[df['tool'] == tool]
    
    hover_texts = []
    for idx, row in tool_data.iterrows():
        items = []
        if text_items:
            items.extend(text_items)
        
        # Add penalty information if available
        if 'optimal_penalty' in row and pd.notna(row['optimal_penalty']):
            items.append(f"Penalty: {row['optimal_penalty']}")
        
        hover_texts.append('<br>'.join(items) if items else '')
    
    return hover_texts


# ============================================================================
# HTML Template with CSS Styling
# ============================================================================

def get_html_template():
    """Return the HTML template with embedded CSS."""
    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        /* Base Styles */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f7fa;
        }}
        
        /* Header */
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 2rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .header h1 {{
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
        }}
        
        /* Navigation */
        .nav {{
            background: white;
            padding: 1rem 2rem;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            position: sticky;
            top: 0;
            z-index: 100;
        }}
        
        .nav a {{
            color: #667eea;
            text-decoration: none;
            margin-right: 2rem;
            font-weight: 500;
            transition: all 0.3s;
        }}
        
        .nav a:hover {{
            color: #764ba2;
        }}
        
        /* Main Container */
        .container {{
            max-width: 100%;
            margin: 2rem auto;
            padding: 0 2rem;
        }}
        
        /* Sections */
        .section {{
            background: white;
            padding: 2rem;
            margin-bottom: 2rem;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        
        .section h2 {{
            color: #667eea;
            margin-bottom: 1.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid #667eea;
        }}
        
        /* Tables */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
        }}
        
        th {{
            background: #667eea;
            color: white;
            padding: 0.75rem;
            text-align: left;
            font-weight: 600;
        }}
        
        td {{
            padding: 0.75rem;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        tr:hover {{
            background: #f9f9f9;
        }}
        
        /* Tabs */
        .tabs {{
            display: flex;
            border-bottom: 2px solid #e0e0e0;
            margin-bottom: 1.5rem;
        }}
        
        .tab {{
            padding: 0.75rem 1.5rem;
            cursor: pointer;
            border: none;
            background: none;
            font-size: 1rem;
            color: #666;
            transition: all 0.3s;
        }}
        
        .tab:hover {{
            background: #f0f0f0;
        }}
        
        .tab.active {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            font-weight: 600;
        }}
        
        .tab-content {{
            display: none;
        }}
        
        .tab-content.active {{
            display: block;
        }}
        
        /* Plots */
        .plot-container {{
            margin: 2rem 0;
            width: 100%;
            overflow-x: visible;
        }}
        
        /* Force Plotly divs to fill container */
        .plot-container > div,
        .plot-box > div,
        .js-plotly-plot {{
            width: 100% !important;
        }}

        /* Flex layout for plot + controls */
        .plot-flex {{
            display: flex;
            gap: 1rem;
            align-items: flex-start;
            width: 100%;
        }}

        .plot-box {{
            flex: 1 1 auto;
            min-width: 0; /* allow flex child to shrink so plotly's div can be responsive */
            width: 100%;
        }}

        .plot-controls {{
            flex: 0 0 200px;
            width: 200px;
            max-height: 600px;
            overflow-y: auto;
            overflow-x: hidden;
            background: #fbfbfb;
            border: 1px solid #ebebeb;
            padding: 0.75rem;
            border-radius: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04);
            font-size: 13px;
        }}

        .plot-controls .controls-header {{ 
            margin-bottom: 0.5rem; 
            position: sticky;
            top: 0;
            background: #fbfbfb;
            padding-bottom: 0.25rem;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .plot-controls .controls-section {{ 
            margin-bottom: 0.75rem; 
        }}
        
        .plot-controls label {{
            display: block;
            margin: 0.25rem 0;
            cursor: pointer;
        }}
        
        .plot-controls label:hover {{
            background: #f0f0f0;
            padding-left: 2px;
        }}
        
        /* Info Box */
        .info-box {{
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 1rem;
            margin: 1rem 0;
            border-radius: 4px;
        }}
        
        /* Footer */
        .footer {{
            text-align: center;
            padding: 2rem;
            color: #666;
            font-size: 0.9rem;
        }}
    </style>
    <script>
        function switchTab(tabName, event) {{
            // Hide all tab contents
            document.querySelectorAll('.tab-content').forEach(content => {{
                content.classList.remove('active');
            }});
            
            // Remove active from all tabs
            document.querySelectorAll('.tab').forEach(tab => {{
                tab.classList.remove('active');
            }});
            
            // Show selected tab content
            document.getElementById(tabName).classList.add('active');
            if (event && event.target) {{
                event.target.classList.add('active');
            }}
            
            // Force resize all Plotly plots in the newly visible tab
            setTimeout(function() {{
                var activeTab = document.getElementById(tabName);
                if (activeTab) {{
                    var plots = activeTab.querySelectorAll('.js-plotly-plot');
                    plots.forEach(function(plot) {{
                        Plotly.Plots.resize(plot);
                    }});
                }}
            }}, 50);
        }}
        
        // Auto-resize all plots when page loads and when window resizes
        window.addEventListener('load', function() {{
            setTimeout(function() {{
                var allPlots = document.querySelectorAll('.js-plotly-plot');
                allPlots.forEach(function(plot) {{
                    Plotly.Plots.resize(plot);
                }});
            }}, 200);
        }});
        
        window.addEventListener('resize', function() {{
            var allPlots = document.querySelectorAll('.js-plotly-plot');
            allPlots.forEach(function(plot) {{
                Plotly.Plots.resize(plot);
            }});
        }});
    </script>
</head>
<body>
    <div class="header">
        <h1>{title}</h1>
        <p>{subtitle}</p>
    </div>
    
    <div class="nav">
        <a href="#summary">Summary</a>
        <a href="#experiment">Experiment Details</a>
        <a href="#visualizations">Visualizations</a>
    </div>
    
    <div class="container">
        {content}
    </div>
    
    <div class="footer">
        Generated by SigFitTest Workflow | {timestamp}
    </div>
</body>
</html>
"""


# ============================================================================
# Data Loading
# ============================================================================

def load_evaluation_data(tsv_files, return_raw=False):
    """
    Load evaluation metrics from TSV file(s).
    
    Args:
        tsv_files: List of paths to TSV files with evaluation metrics
        return_raw: If True, return tuple (aggregated_df, raw_df) where raw_df contains
                   individual realization data before aggregation
        
    Returns:
        pandas.DataFrame with combined evaluation data (aggregated across realizations if present)
        OR tuple (aggregated_df, raw_df) if return_raw=True
        
    Note: The "raw" data here means per-realization summary metrics (already aggregated per sample).
          For true per-sample data, use load_per_sample_data().
    """
    print(f"Loading {len(tsv_files)} file(s)...")
    
    dataframes = []
    for tsv_file in tsv_files:
        print(f"  Reading: {tsv_file}")
        df = pd.read_csv(tsv_file, sep='\t')
        print(f"    Shape: {df.shape}")
        dataframes.append(df)
    
    # Combine all dataframes
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Clean tool names: remove realization suffix (_w0, _w1, etc.)
    if 'tool' in combined_df.columns:
        combined_df['tool'] = combined_df['tool'].str.replace(r'_w\d+$', '', regex=True)
        print(f"  Cleaned tool names (removed realization suffixes)")
    
    # Check if realization column exists
    has_realizations = 'realization' in combined_df.columns
    
    if not has_realizations:
        print("  No realization column found - data is already aggregated")
        if return_raw:
            return combined_df, combined_df.copy()
        else:
            return combined_df
    
    # Store raw per-realization data (for stability plots)
    raw_df = combined_df.copy()
    num_realizations = raw_df['realization'].nunique()
    print(f"  Preserved raw per-realization data: {len(raw_df)} rows across {num_realizations} realizations")
    
    # Create aggregated view (averaging across realizations)
    print(f"  Creating aggregated view (averaging across {num_realizations} realizations)...")
    
    # Group by all columns except realization and metrics
    group_cols = ['tool', 'data_type', 'relevance', 'cancer_type', 'num_muts']
    
    # Check if data already has _mean/_std columns (from per-sample aggregation)
    has_mean_std_cols = any(col.endswith('_mean') for col in combined_df.columns)
    
    if has_mean_std_cols:
        # Data already has _mean and _std columns from per-sample aggregation
        # Average these across realizations
        print("    Data has _mean/_std columns, averaging across realizations...")
        agg_dict = {}
        
        # Keep these columns as-is (same across realizations)
        if 'num_samples' in combined_df.columns:
            agg_dict['num_samples'] = 'first'
        if 'mean_num_muts' in combined_df.columns:
            agg_dict['mean_num_muts'] = 'first'
        
        # Get all numeric columns
        numeric_cols = combined_df.select_dtypes(include=[np.number]).columns.tolist()
        
        # Process all numeric columns - just average them
        for col in numeric_cols:
            if col in group_cols + ['realization', 'num_samples', 'mean_num_muts']:
                continue
            agg_dict[col] = 'mean'
        
        combined_df = combined_df.groupby(group_cols, as_index=False).agg(agg_dict)
    else:
        # Data doesn't have _mean/_std yet, create them by aggregating
        print("    Creating _mean/_std columns from realizations...")
        numeric_cols = combined_df.select_dtypes(include=[np.number]).columns.tolist()
        metric_cols = [col for col in numeric_cols if col not in group_cols + ['realization']]
        
        # Build aggregation dict
        agg_dict = {}
        if 'num_samples' in combined_df.columns:
            agg_dict['num_samples'] = 'first'
        if 'mean_num_muts' in combined_df.columns:
            agg_dict['mean_num_muts'] = 'first'
        
        for col in metric_cols:
            agg_dict[col] = ['mean', 'std']
        
        combined_df = combined_df.groupby(group_cols, as_index=False).agg(agg_dict)
        
        # Flatten column names
        combined_df.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col for col in combined_df.columns]
    
    print(f"  Aggregated shape: {combined_df.shape}")
    
    if return_raw:
        return combined_df, raw_df
    else:
        return combined_df


def reaggregate_from_per_sample(df_ps):
    """Re-aggregate per-sample data (after per-sample oracle penalty selection) into summary format.

    Each sample has already been assigned its best-TAE penalty independently.
    Grouping here drops penalty so the resulting df has one row per
    (cancer_type, data_type, tool, num_muts, replicate).
    """
    group_cols = [c for c in ['cancer_type', 'data_type', 'tool', 'num_muts', 'replicate']
                  if c in df_ps.columns]

    agg_dict = dict(
        num_samples=('sample', 'count'),
        mean_num_muts=('num_muts', 'first'),
        TAE=('TAE', 'mean'),
        TAE_std=('TAE', 'std'),
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
    )
    if 'nRMSE' in df_ps.columns:
        agg_dict['nRMSE'] = ('nRMSE', lambda x: float(np.sqrt(np.mean(np.array(x) ** 2))))
    for col in ['theoretical_min_TAE', 'theoretical_max_Recall',
                'theoretical_max_Precision', 'theoretical_max_F1']:
        if col in df_ps.columns:
            agg_dict[col] = (col, 'mean')

    agg = df_ps.groupby(group_cols).agg(**agg_dict).reset_index()

    n = agg['num_samples'].clip(lower=1)
    for m in ['P', 'R', 'S', 'F1', 'MCC']:
        agg[f'{m}_sem'] = agg[f'{m}_std'] / np.sqrt(n)

    return agg


def load_per_sample_data(tsv_files):
    """
    Load per-sample metrics from TSV file(s).
    Does NOT deduplicate because we want all individual samples.
    
    Args:
        tsv_files: List of paths to TSV files with per-sample metrics
        
    Returns:
        pandas.DataFrame with combined per-sample data (all samples preserved)
    """
    print(f"Loading {len(tsv_files)} per-sample file(s)...")
    
    dataframes = []
    for tsv_file in tsv_files:
        print(f"  Reading: {tsv_file}")
        df = pd.read_csv(tsv_file, sep='\t')
        print(f"    Shape: {df.shape}")
        dataframes.append(df)
    
    # Combine all dataframes WITHOUT deduplication
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Clean tool names: remove realization suffix (_w0, _w1, etc.)
    if 'tool' in combined_df.columns:
        combined_df['tool'] = combined_df['tool'].str.replace(r'_w\d+$', '', regex=True)
        print(f"  Cleaned tool names (removed realization suffixes)")
    
    print(f"Combined per-sample data shape: {combined_df.shape}")
    print(f"  Total samples: {len(combined_df)}")
    print(f"  Tools: {sorted(combined_df['tool'].unique())}")
    print(f"  Cancer types: {sorted(combined_df['cancer_type'].unique())}")
    
    return combined_df


# ============================================================================
# HTML Section Generators
# ============================================================================

def create_summary_section(df, df_raw=None, unknown_sig_config=None):
    """Create summary statistics section with experimental parameters.
    
    Args:
        df: DataFrame with evaluation metrics
        df_raw: Optional raw per-realization data
        unknown_sig_config: Optional dict with unknown signature configuration
                          {'num_unknowns': int, 'weights': list}
    """
    n_tools = df['tool'].nunique()
    n_cancer_types = df['cancer_type'].nunique()
    n_conditions = len(df)
    num_samples = int(df['num_samples'].iloc[0]) if 'num_samples' in df.columns else 'N/A'
    
    # Get number of realizations
    num_realizations = 'N/A'
    if df_raw is not None and 'realization' in df_raw.columns:
        num_realizations = int(df_raw['realization'].nunique())
    
    # Get unique tools in custom order
    tool_order = get_tool_display_order()
    available_tools = set(df['tool'].unique())
    tools_list = [t for t in tool_order if t in available_tools]
    # Add any tools not in our predefined order
    remaining_tools = sorted(available_tools - set(tools_list))
    tools_list.extend(remaining_tools)
    
    # Capitalize tool names for display: musical -> MuSiCal, spa -> SPA, etc.
    tool_display_names = []
    for tool in tools_list:
        if tool.lower() == 'musical':
            tool_display_names.append('MuSiCal')
        elif tool.lower() == 'spa':
            tool_display_names.append(tool.upper())
        else:
            tool_display_names.append(tool)
    tools = ', '.join(tool_display_names)
    
    cancer_types = ', '.join(sorted(df['cancer_type'].unique()))
    
    # Get mutation counts tested
    mutation_counts = sorted(df['num_muts'].unique())
    mut_counts_str = ', '.join([f"{int(m):,}" for m in mutation_counts])
    
    # Get data type and relevance info
    data_type = df['data_type'].iloc[0] if 'data_type' in df.columns else 'N/A'
    relevance = df['relevance'].iloc[0] if 'relevance' in df.columns else 'N/A'
    
    # Format unknown signature information
    if unknown_sig_config and unknown_sig_config.get('num_unknowns', 0) > 0:
        num_unknowns = unknown_sig_config['num_unknowns']
        weights = unknown_sig_config.get('weights', [])
        if weights:
            weights_pct = [f"{w*100:.1f}%" for w in weights]
            unknown_weights_str = ', '.join(weights_pct)
        else:
            unknown_weights_str = 'N/A'
    else:
        num_unknowns = 0
        unknown_weights_str = 0
    
    html = f'''
    <section id="summary" class="section">
        <h2>Summary</h2>
        
        <div class="info-box" style="background: #f5f5f5; border-left-color: #666; margin-bottom: 1.5rem;">
            <h3 style="margin-bottom: 0.5rem;">Background</h3>
            <p style="margin-bottom: 0.75rem;">
                Somatic mutations accumulate in cancer cells through distinct biological processes, each leaving a characteristic pattern
                across the 96 single-base substitution (SBS) contexts — a <strong>mutational signature</strong>. Decomposing a tumour's
                mutation catalogue into its contributing signatures (<em>signature fitting</em>) helps identify the active mutational
                processes, which in turn can inform treatment decisions, aetiology studies, and biomarker discovery.
            </p>
            <p style="margin-bottom: 0.75rem;">
                Many tools exist for signature fitting, but their relative performance is rarely evaluated under fair, controlled conditions.
                Key questions include: <em>How well does a tool recover the correct signatures?</em> <em>Does it over-call rare signatures
                or miss low-contribution ones?</em> <em>How accurately does it estimate the contribution of each signature?</em>
                <em>How does performance degrade at low mutation burdens?</em>
            </p>
            <p style="margin-bottom: 0.75rem;">
                <strong>SigFitTest</strong> addresses this by generating synthetic tumour samples with a <strong>known ground truth</strong>:
                signature identities and their exact contributions are fixed before fitting, so every tool's output can be compared against
                the true answer using a range of metrics (F1, Precision, Recall, TAE, nRMSE, MCC). The benchmark spans multiple
                cancer types (each with its own COSMIC tissue-specific signature landscape), multiple mutation burdens, and optionally
                includes <strong>out-of-reference (unknown) signatures</strong> — signals present in the data but absent from the tool's
                reference catalog — to test robustness to novel or rare mutational processes.
            </p>
            <p style="margin: 0; font-size: 0.9em; color: #555;">
                Based on the methodology described in:
                <span style="display: block; margin: 0.4rem 0 0 1rem;">
                    <strong>Matúš Medo, Charlotte K Y Ng &amp; Michaela Medová</strong> —
                    <em>"A comprehensive comparison of tools for fitting mutational signatures"</em>
                </span>
            </p>
        </div>
        
        <details>
        <summary style="font-weight: 600; color: #374151; padding: 0.5rem 0; cursor: pointer; user-select: none; list-style: none; display: flex; align-items: center; gap: 0.4rem;">
            <span>▶</span> <span>Benchmark Overview (click to expand)</span>
        </summary>
        <table>
            <tr>
                <th style="width: 22%;">Parameter</th>
                <th style="width: 28%;">Value</th>
                <th style="width: 50%; color: #555; font-weight: normal; font-style: italic;">What it means</th>
            </tr>
            <tr>
                <td><strong>Tools Compared</strong></td>
                <td>{n_tools} ({tools})</td>
                <td style="color: #555; font-size: 0.9em;">The signature fitting algorithms being benchmarked head-to-head.</td>
            </tr>
            <tr>
                <td><strong>Cancer Types</strong></td>
                <td>{n_cancer_types} ({cancer_types})</td>
                <td style="color: #555; font-size: 0.9em;">Each cancer type uses its own COSMIC tissue-specific signature catalog, so the set of expected signatures differs per type.</td>
            </tr>
            <tr>
                <td><strong>Mutation Burdens</strong></td>
                <td>{mut_counts_str}</td>
                <td style="color: #555; font-size: 0.9em;">Number of somatic mutations simulated per sample. Low burdens are harder — fewer mutations = noisier signature estimates.</td>
            </tr>
            <tr>
                <td><strong>Samples per Condition</strong></td>
                <td>{num_samples}</td>
                <td style="color: #555; font-size: 0.9em;">Number of independent synthetic samples generated for each (cancer type × mutation burden) combination.</td>
            </tr>
            <tr>
                <td><strong>Realizations</strong></td>
                <td>{num_realizations}</td>
                <td style="color: #555; font-size: 0.9em;">Independent repetitions of the full simulation with different random seeds. Used to assess algorithm stability (variance across runs).</td>
            </tr>
            <tr>
                <td><strong>Unknown Signatures</strong></td>
                <td>{num_unknowns}</td>
                <td style="color: #555; font-size: 0.9em;">Signatures present in the simulated data but deliberately withheld from the tool's reference catalog. Tests how tools handle signals they cannot explain.</td>
            </tr>
            <tr>
                <td><strong>Unknown Sig. Weights</strong></td>
                <td>{unknown_weights_str}</td>
                <td style="color: #555; font-size: 0.9em;">Fraction of each sample's mutations contributed by the unknown signature(s). Larger weight = stronger out-of-reference signal.</td>
            </tr>
            <tr>
                <td><strong>Data Type</strong></td>
                <td>{data_type}</td>
                <td style="color: #555; font-size: 0.9em;"><em>synthetic</em>: mutations drawn from COSMIC profiles with known ground-truth weights.</td>
            </tr>
            <tr>
                <td><strong>Signature Relevance</strong></td>
                <td>{relevance}</td>
                <td style="color: #555; font-size: 0.9em;"><em>standard</em>: tools see all COSMIC signatures as candidates. <em>relevant</em>: reference is pre-filtered to tissue-specific signatures only.</td>
            </tr>
            <tr>
                <td><strong>Total Conditions</strong></td>
                <td>{n_conditions}</td>
                <td style="color: #555; font-size: 0.9em;">tools × cancer types × mutation burdens — each condition is one row in the results tables.</td>
            </tr>
        </table>
        </details>

        <div style="margin-top: 1rem; padding: 0.75rem; background: #f8f9fa; border-radius: 4px; font-size: 0.9em;">
            <strong>Experimental Design:</strong> Synthetic data with configurable unknown signatures excluded from reference catalog.
            Each sample has randomly selected unknown signatures from a predefined pool, mimicking realistic scenarios with rare or novel mutational processes.
            The number and contribution of unknown signatures varies based on configuration settings.
        </div>

        <details style="margin-top: 1.5rem;">
            <summary style="font-weight: 600; color: #374151; padding: 0.5rem 0; cursor: pointer; user-select: none; list-style: none; display: flex; align-items: center; gap: 0.4rem;">
                <span>▶</span> <span>Metric Reference (click to expand)</span>
            </summary>
            <div style="padding: 1rem; background: #f8f9fa; border-radius: 4px; margin-top: 0.5rem; font-size: 0.9em; border: 1px solid #e5e7eb;">
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 0.6rem; color: #374151; line-height: 1.5;">
                    <div><strong>P — Precision</strong>: fraction of reported signatures that are correct (penalizes over-fitting / false positives). Range 0–1, higher is better.</div>
                    <div><strong>R — Recall</strong>: fraction of true signatures that were found (penalizes missing signatures / false negatives). Range 0–1, higher is better.</div>
                    <div><strong>F1</strong>: harmonic mean of P and R — balances both. Range 0–1, higher is better. 0.8+ is excellent.</div>
                    <div><strong>TAE — Total Absolute Error</strong>: sum of |estimated − true| contributions across all signatures. Lower is better.</div>
                    <div><strong>nRMSE — Normalized RMSE</strong>: RMS of contribution errors normalized by mutation count. Lower is better, more sensitive to large outliers than TAE.</div>
                    <div><strong>MCC — Matthews Correlation Coefficient</strong>: binary classification quality (signature present/absent), accounting for class imbalance. Range −1 to +1, higher is better.</div>
                    <div><strong>FP / FN</strong>: number of False Positive (over-called) / False Negative (missed) signature calls per sample.</div>
                </div>
            </div>
        </details>
    </section>
    '''
    return html


def _derive_cosmic_label(base_file):
    """Build a human-readable COSMIC catalog label from the configured base_file path.

    e.g. "input/COSMIC_v3.3.1_SBS_GRCh38.txt" -> "COSMIC v3.3.1 SBS signatures
    (single base substitutions, GRCh38)". Falls back to a generic label if the
    filename doesn't match the expected COSMIC naming convention.
    """
    if base_file:
        m = re.match(r"COSMIC_(v[\d.]+)_SBS_(GRCh\d+)", os.path.basename(base_file))
        if m:
            version, genome = m.groups()
            return f"COSMIC {version} SBS signatures (single base substitutions, {genome})"
    return "COSMIC SBS signatures"


def create_experiment_section(df, df_raw=None, unknown_sig_config=None, cosmic_version_label=None):
    """Create experiment details section with comprehensive explanations.

    Args:
        df: DataFrame with evaluation metrics
        df_raw: Optional raw per-realization data
        unknown_sig_config: Optional dict with unknown signature configuration
                          {'num_unknowns': int, 'weights': list}
        cosmic_version_label: Optional human-readable COSMIC catalog label derived
                          from config's cosmic_signatures.base_file
    """
    cosmic_version_label = cosmic_version_label or "COSMIC SBS signatures"
    num_samples = int(df['num_samples'].iloc[0]) if 'num_samples' in df.columns else 'N/A'
    
    # Get number of realizations
    num_realizations = 'N/A'
    if df_raw is not None and 'realization' in df_raw.columns:
        num_realizations = int(df_raw['realization'].nunique())
    
    # Build unknown signatures description from config
    if unknown_sig_config and unknown_sig_config.get('num_unknowns', 0) > 0:
        num_unknowns = unknown_sig_config['num_unknowns']
        weights = unknown_sig_config.get('weights', [])
        
        unknown_sig_description = f"{num_unknowns} unknown signature(s) per sample, randomly selected from excluded COSMIC signatures"
        
        if weights:
            # Convert weights to percentages
            weights_pct = [f"{w*100:.1f}%" for w in weights]
            unknown_sig_description += f"<br><strong>Contribution weights:</strong> {', '.join(weights_pct)} of total mutations per unknown signature"
    else:
        # Fallback: Try to infer from dataframe or use generic description
        num_unknowns = 'N/A'
        if 'num_unknown_sigs' in df.columns:
            num_unknowns = int(df['num_unknown_sigs'].iloc[0])
            if num_unknowns > 0:
                unknown_sig_description = f"{num_unknowns} unknown signature(s) per sample, randomly selected from excluded COSMIC signatures"
            else:
                unknown_sig_description = "No unknown signatures (all signatures available in reference catalog)"
        else:
            unknown_sig_description = "Selected signatures excluded from reference catalog (randomly chosen per sample)"
    
    html = f'''
    <section id="experiment" class="section">
        <h2>Experimental Setup</h2>
        
        <details>
        <summary style="font-weight: 600; color: #374151; padding: 0.5rem 0; cursor: pointer; user-select: none; list-style: none; display: flex; align-items: center; gap: 0.4rem;">
            <span>▶</span> <span>Data Generation (click to expand)</span>
        </summary>
        <div style="padding: 1rem; background: #f8f9fa; border-radius: 4px; margin-top: 0.5rem; margin-bottom: 1rem;">
            <table style="width: 100%; margin-top: 0;">
                <tr>
                    <th style="width: 30%;">Parameter</th>
                    <th>Details</th>
                </tr>
                <tr>
                    <td><strong>Signature Catalog</strong></td>
                    <td>{cosmic_version_label}</td>
                </tr>
                <tr>
                    <td><strong>Unknown Signatures</strong></td>
                    <td>{unknown_sig_description}</td>
                </tr>
                <tr>
                    <td><strong>Known Signatures</strong></td>
                    <td>Cancer-type-specific COSMIC signatures from reference catalog</td>
                </tr>
                <tr>
                    <td><strong>Samples per Condition</strong></td>
                    <td>{num_samples}</td>
                </tr>
                <tr>
                    <td><strong>Independent Realizations</strong></td>
                    <td>{num_realizations} (for assessing algorithm stability)</td>
                </tr>
            </table>
        </div>
        </details>

        <details>
        <summary style="font-weight: 600; color: #374151; padding: 0.5rem 0; cursor: pointer; user-select: none; list-style: none; display: flex; align-items: center; gap: 0.4rem;">
            <span>▶</span> <span>Software Versions (click to expand)</span>
        </summary>
        <div style="padding: 1rem; background: #f8f9fa; border-radius: 4px; margin-top: 0.5rem; margin-bottom: 1rem;">
            <table style="width: 100%; margin-top: 0;">
                <tr>
                    <th style="width: 30%;">Tool</th>
                    <th>Version & Implementation</th>
                </tr>
                <tr>
                    <td><strong>MuSiCal</strong></td>
                    <td>Python implementation (custom, based on Medo et al. 2023)</td>
                </tr>
                <tr>
                    <td><strong>SPA</strong></td>
                    <td>Python 3.x - Signature Presence Analysis</td>
                </tr>
                <tr>
                    <td><strong>SigID</strong></td>
                    <td>SPA-style greedy add/remove with automatic per-sample L1 penalty selection for unknown-signature discovery</td>
                </tr>
                <tr>
                    <td><strong>SigTrim</strong></td>
                    <td>Iterative unassigned-mutation trimming followed by SPA-style greedy add/remove</td>
                </tr>
                <tr>
                    <td><strong>Python Environment</strong></td>
                    <td>Python 3.10+, NumPy 1.24+, Pandas 2.0+, SciPy 1.11+</td>
                </tr>
                <tr>
                    <td><strong>Workflow Manager</strong></td>
                    <td>Snakemake 7.32+ for reproducible execution</td>
                </tr>
            </table>
        </div>
        </details>
    </section>
    '''
    return html


def create_performance_metrics_section():
    """Create Performance Metrics section with comprehensive explanations."""
    html = '''
    <section id="metrics" class="section">
        <h2>Performance Metrics</h2>
        
        <div class="info-box" style="background: #e3f2fd; border-left-color: #2196f3; margin-bottom: 1.5rem;">
            <h4 style="margin-bottom: 0.5rem;">Acronym Reference</h4>
            <p style="margin-bottom: 0.3rem;"><strong>P</strong> = Precision</p>
            <p style="margin-bottom: 0.3rem;"><strong>R</strong> = Recall (Sensitivity)</p>
            <p style="margin-bottom: 0.3rem;"><strong>F1</strong> = F1 Score (harmonic mean of P and R)</p>
            <p style="margin-bottom: 0.3rem;"><strong>TAE</strong> = Total Absolute Error</p>
            <p style="margin-bottom: 0.3rem;"><strong>nRMSE</strong> = Normalized Root Mean Squared Error</p>
            <p style="margin-bottom: 0.3rem;"><strong>MCC</strong> = Matthews Correlation Coefficient</p>
            <p style="margin-bottom: 0.3rem;"><strong>FP</strong> = False Positive | <strong>FN</strong> = False Negative | <strong>TP</strong> = True Positive</p>
        </div>
        
        <table>
            <tr>
                <th style="width: 20%;">Metric</th>
                <th style="width: 50%;">Description</th>
                <th style="width: 30%;">Interpretation</th>
            </tr>
            <tr>
                <td><strong>F1 Score</strong></td>
                <td>Harmonic mean of precision and recall. Balances the trade-off between 
                    identifying correct signatures and avoiding false positives.
                    <br><br><strong>Formula:</strong> <code>F1 = 2&middot;(P&middot;R)/(P+R) = 2&middot;TP / (2&middot;TP + FP + FN)</code></td>
                <td>Higher is better<br>Range: 0-1<br>0.8+ is excellent</td>
            </tr>
            <tr>
                <td><strong>Precision (P)</strong></td>
                <td>Proportion of identified signatures that are truly present. 
                    Low precision means many false positives.
                    <br><br><strong>Formula:</strong> <code>P = TP / (TP + FP)</code></td>
                <td>Higher is better<br>Range: 0-1<br>Penalizes over-fitting</td>
            </tr>
            <tr>
                <td><strong>Recall (R)</strong></td>
                <td>Proportion of true signatures that were successfully identified. 
                    Low recall means missing real signatures.
                    <br><br><strong>Formula:</strong> <code>R = TP / (TP + FN)</code></td>
                <td>Higher is better<br>Range: 0-1<br>Penalizes under-fitting</td>
            </tr>
            <tr>
                <td><strong>TAE</strong></td>
                <td>Total Absolute Error - sum of absolute differences between estimated and 
                    true signature contributions.
                    <br><br><strong>Formula:</strong> <code>TAE = &sum;|w<sub>i</sub><sup>true</sup> &minus; w<sub>i</sub><sup>est</sup>|</code>
                    <br>where <code>w<sub>i</sub></code> is the contribution (weight) of signature <code>i</code></td>
                <td>Lower is better<br>Measures weight accuracy</td>
            </tr>
            <tr>
                <td><strong>nRMSE</strong></td>
                <td>Normalized Root Mean Square Error - RMS of errors in signature 
                    weight estimation, normalized by total mutation count.
                    <br><br><strong>Formula:</strong> <code>nRMSE = &radic;( (1/n) &sum;((w<sub>i</sub><sup>true</sup> &minus; w<sub>i</sub><sup>est</sup>) / N)&sup2; )</code>
                    <br>where <code>n</code> is the number of signatures and <code>N</code> is the total number of mutations</td>
                <td>Lower is better<br>Sensitive to outliers<br>Per-sample metric</td>
            </tr>
            <tr>
                <td><strong>MCC</strong></td>
                <td>Matthews Correlation Coefficient - quality of binary classification 
                    (signature present/absent), accounting for class imbalance.
                    <br><br><strong>Formula:</strong> <code>MCC = (TP&middot;TN &minus; FP&middot;FN) / &radic;((TP+FP)(TP+FN)(TN+FP)(TN+FN))</code></td>
                <td>Higher is better<br>Range: -1 to +1<br>0.5+ is good</td>
            </tr>
            <tr>
                <td><strong>False Positives (FP)</strong></td>
                <td>Number of signatures incorrectly identified as present when they are absent (over-fitting indicators).
                    <br><br><strong>Definition:</strong> Signatures reported by the tool but not in the ground truth</td>
                <td>Lower is better<br>Indicates over-fitting</td>
            </tr>
            <tr>
                <td><strong>False Negatives (FN)</strong></td>
                <td>Number of true signatures missed by the tool (under-fitting indicators).
                    <br><br><strong>Definition:</strong> Signatures in the ground truth but not reported by the tool</td>
                <td>Lower is better<br>Indicates under-fitting</td>
            </tr>
            <tr>
                <td><strong>True Positives (TP)</strong></td>
                <td>Number of signatures correctly identified as present.
                    <br><br><strong>Definition:</strong> Signatures both in ground truth and reported by the tool</td>
                <td>Higher is better<br>Correct detection</td>
            </tr>
            <tr>
                <td><strong>True Negatives (TN)</strong></td>
                <td>Number of signatures correctly identified as absent.
                    <br><br><strong>Definition:</strong> Signatures neither in ground truth nor reported by the tool</td>
                <td>Higher is better<br>Correct rejection</td>
            </tr>
        </table>
        
        <div class="info-box" style="margin-top: 1.5rem;">
            <h4>What Makes a Good Tool?</h4>
            <ul style="margin-left: 1.5rem; margin-top: 0.5rem;">
                <li><strong>High F1 Score:</strong> Correctly identifies signatures without too many false positives</li>
                <li><strong>Low TAE/nRMSE:</strong> Accurately estimates signature contributions</li>
                <li><strong>Robustness:</strong> Performs well across different mutation burdens and cancer types</li>
                <li><strong>Stability:</strong> Consistent results across independent realizations</li>
            </ul>
        </div>
    </section>
    '''
    return html


# ============================================================================
# Plotting Functions
# ============================================================================

def get_tool_colors():
    """Return consistent color scheme for tools (uses ReportMetadata)."""
    return ReportMetadata().tool_colors


def get_tool_base_method(tool):
    """Return the base method for a tool (uses ReportMetadata)."""
    return ReportMetadata().get_base_method(tool)


def get_method_symbols():
    """Return marker symbols for different base methods (uses ReportMetadata)."""
    return ReportMetadata().method_symbols


def get_tool_display_order():
    """Return preferred order for tool display (uses ReportMetadata)."""
    return ReportMetadata().tool_display_order


def create_precision_recall_plot(df, cancer_type, tools, colors, df_per_sample=None):
    """Create precision-recall scatter plot with marker size by mutation count."""
    fig = go.Figure()
    
    method_symbols = get_method_symbols()
    
    for tool in tools:
        tool_data = df[df['tool'] == tool]
        
        if len(tool_data) == 0:
            continue
        
        # Get base method to determine symbol
        base_method = get_tool_base_method(tool)
        symbol = method_symbols.get(base_method, 'circle')
        
        # Scale marker size by mutation count (log scale for visibility)
        # Reduced from 8+15 to 5+10 for smaller overall markers while keeping size variation
        marker_sizes = 5 + 10 * np.log10(tool_data['num_muts']) / np.log10(tool_data['num_muts'].max())
        
        # Build custom hover text with mutation count and penalty (if applicable)
        hover_text_items = []
        for _, row in tool_data.iterrows():
            items = [f"Mutations: {int(row['num_muts']):,}"]
            if 'optimal_penalty' in row and pd.notna(row['optimal_penalty']):
                items.append(f"Penalty: {row['optimal_penalty']}")
            hover_text_items.append('<br>'.join(items))
            
        fig.add_trace(go.Scatter(
            x=tool_data['R_mean'],
            y=tool_data['P_mean'],
            mode='markers',
            name=tool,
            marker=dict(
                size=marker_sizes,
                color=colors.get(tool, '#999'),
                symbol=symbol,
                opacity=0.6,  # Transparency so overlapping markers are visible
                line=dict(width=1, color='black')
            ),
            error_x=dict(array=tool_data['R_sem'], type='data'),
            error_y=dict(array=tool_data['P_sem'], type='data'),
            text=hover_text_items,
            hovertemplate='<b>%{fullData.name}</b><br>Recall: %{x:.3f}<br>Precision: %{y:.3f}<br>%{text}<extra></extra>'
        ))
    
    # Add F1 score reference lines
    for f1_val, f1_label in [(0.5, 'F1=0.5'), (0.7, 'F1=0.7'), (0.9, 'F1=0.9')]:
        recall = np.linspace(0.01, 1, 100)
        precision = f1_val * recall / (2 * recall - f1_val * recall)
        precision = np.clip(precision, 0, 1)
        
        fig.add_trace(go.Scatter(
            x=recall,
            y=precision,
            mode='lines',
            name=f1_label,
            line=dict(dash='dash', color='gray', width=1),
            opacity=0.5,
            hoverinfo='name'
        ))
    
    # Add theoretical threshold lines if unknown signatures present
    # Note: With unknown signatures present, theoretical max Recall is < 1.0
    # because unknowns will be counted as FN in the evaluation
    # Add theoretical max Recall and Precision lines if available in per-sample metrics
    if df_per_sample is not None:
        cancer_per_sample = df_per_sample[df_per_sample['cancer_type'] == cancer_type]
        if len(cancer_per_sample) > 0:
            # The per-sample evaluation step writes theoretical_max_Recall and theoretical_max_Precision
            # We will plot the average across samples and always show the line so the reader can see the
            # theoretical performance ceiling (even if it's near 1.0).
            if 'theoretical_max_Recall' in cancer_per_sample.columns:
                theoretical_recall = cancer_per_sample['theoretical_max_Recall'].mean()
                if not np.isnan(theoretical_recall):
                    fig.add_trace(go.Scatter(
                        x=[theoretical_recall, theoretical_recall],
                        y=[0, 1],
                        mode='lines',
                        name='Max Recall',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=True,
                        hovertemplate=f'Theoretical Max Recall: {theoretical_recall:.3f}<br>(Unknown sigs counted as FN)<extra></extra>'
                    ))

            if 'theoretical_max_Precision' in cancer_per_sample.columns:
                theoretical_precision = cancer_per_sample['theoretical_max_Precision'].mean()
                if not np.isnan(theoretical_precision):
                    fig.add_trace(go.Scatter(
                        x=[0, 1],
                        y=[theoretical_precision, theoretical_precision],
                        mode='lines',
                        name='Max Precision',
                        line=dict(color='blue', width=2, dash='dot'),
                        showlegend=True,
                        hovertemplate=f'Theoretical Max Precision: {theoretical_precision:.3f}<extra></extra>'
                    ))
    
    fig.update_layout(
        title=dict(
            text=f'{cancer_type}',
            font=dict(size=16)
        ),
        xaxis_title='R (Recall) ↑',
        yaxis_title='P (Precision) ↑',
        xaxis=dict(range=[0, 1.05]),
        yaxis=dict(range=[0, 1.05]),
        height=500,
        hovermode='closest',
        showlegend=True,
        legend=dict(itemdoubleclick='toggleothers'),
        autosize=True
    )
    
    return fig


def create_violin_plots(df, cancer_type, tools, colors):
    """Create scatter plots with error bars showing metric variability across conditions.
    
    Note: This shows aggregated data (mean ± std across 20 samples per condition).
    Error bars represent the variability across samples for each condition.
    """
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=['F1 Score (Aggregated)', 'TAE (Aggregated)', 
                       'MCC (Aggregated)', 'nRMSE (Aggregated)'],
        vertical_spacing=0.15,
        horizontal_spacing=0.12
    )
    
    # Metrics with their standard deviations
    metrics = [
        ('F1_mean', 'F1_std', 1, 1, 'F1 ↑'),
        ('TAE', 'TAE_std', 1, 2, 'TAE ↓'),
        ('MCC_mean', 'MCC_std', 2, 1, 'MCC ↑'),
        ('nRMSE', 'nRMSE_std', 2, 2, 'nRMSE ↓'),
    ]
    
    for tool in tools:
        tool_data = df[df['tool'] == tool].sort_values('num_muts')
        
        if len(tool_data) == 0:
            continue
        
        for metric_mean, metric_std, row, col, ylabel in metrics:
            # Create error bars if std is available
            if metric_std and metric_std in tool_data.columns:
                error_y = dict(
                    type='data',
                    array=tool_data[metric_std],
                    visible=True
                )
            else:
                error_y = None
            
            fig.add_trace(
                go.Scatter(
                    x=[tool] * len(tool_data),
                    y=tool_data[metric_mean],
                    mode='markers',
                    name=tool,
                    marker=dict(
                        size=10,
                        color=colors.get(tool, '#999'),
                        opacity=0.7
                    ),
                    error_y=error_y,
                    showlegend=(row == 1 and col == 1),
                    legendgroup=tool,
                    text=[f"{int(m):,} muts<br>n={int(n)} samples" for m, n in zip(tool_data['num_muts'], tool_data['num_samples'])],
                    hovertemplate='<b>%{fullData.name}</b><br>%{text}<br>Mean: %{y:.3f}<extra></extra>'
                ),
                row=row, col=col
            )
    
    # Update y-axes labels with direction indicators
    fig.update_yaxes(title_text='F1 ↑', row=1, col=1)
    fig.update_yaxes(title_text='TAE ↓', row=1, col=2)
    fig.update_yaxes(title_text='MCC ↑', row=2, col=1)
    fig.update_yaxes(title_text='nRMSE ↓', row=2, col=2)
    
    # Update x-axes
    for row in [1, 2]:
        for col in [1, 2]:
            fig.update_xaxes(title_text='Tool', row=row, col=col)
    
    fig.update_layout(
        height=700,
        title_text=f'Performance Variability Across Mutation Burdens<br><sub>Error bars show ± std dev across samples</sub>',
        showlegend=True,
        legend=dict(x=1.05, y=1, itemdoubleclick='toggleothers'),
        autosize=True
    )
    
    return fig


def create_per_sample_violin_plots(df_samples, cancer_type, tools, colors):
    """Create boxplots with overlaid sample markers for per-sample data.
    Shows one mutation burden at a time with a dropdown selector.
    X-axis shows methods (tools), each point is a sample.
    
    Args:
        df_samples: DataFrame with per-sample metrics
        cancer_type: Cancer type to filter
        tools: List of tool names
        colors: Dict mapping tool names to colors
    """
    # Filter for this cancer type
    cancer_samples = df_samples[df_samples['cancer_type'] == cancer_type]
    
    if len(cancer_samples) == 0:
        return None
    
    # Get unique mutation burdens sorted
    mut_burdens = sorted(cancer_samples['num_muts'].unique())
    
    # Use 3x2 layout for metrics (boxplots with sample markers)
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=['F1 ↑', 'TAE ↓', 'MCC ↑', 'nRMSE ↓', 'FP Sigs', 'FN Sigs'],
        vertical_spacing=0.18,  # Increased from 0.12 to prevent x-label overlap
        horizontal_spacing=0.12
    )
    
    # Metrics configuration
    metrics = [
        ('F1', 1, 1, 'F1 ↑'),
        ('TAE', 1, 2, 'TAE ↓'),
        ('MCC', 2, 1, 'MCC ↑'),
        ('nRMSE', 2, 2, 'nRMSE ↓'),
        ('num_FP', 3, 1, 'FP Count'),
        ('num_FN', 3, 2, 'FN Count'),
    ]
    metrics = [(m, r, c, l) for m, r, c, l in metrics if m in cancer_samples.columns]
    
    if not metrics:
        print(f"  No valid metrics found for {cancer_type}")
        return None
    
    # Track all traces for visibility control
    all_traces = []
    trace_burdens = []  # Track which burden each trace belongs to
    
    # Create traces for EACH burden separately
    for burden in mut_burdens:
        burden_data = cancer_samples[cancer_samples['num_muts'] == burden]
        
        # For each metric/subplot, calculate median performance and rank tools
        # This creates a separate ordering for each panel
        tool_orders = {}  # Key: (metric, row, col), Value: list of tools sorted by median
        
        for metric, row, col, ylabel in metrics:
            # Calculate median for each tool in this metric
            tool_medians = {}
            for tool in tools:
                tool_data = burden_data[burden_data['tool'] == tool]
                if len(tool_data) > 0:
                    grouped = tool_data.groupby('sample')[metric].mean()
                    if len(grouped) > 0:
                        tool_medians[tool] = grouped.median()
            
            # Sort tools by median performance
            # For metrics where lower is better (TAE, nRMSE, FP, FN), reverse sort
            reverse_sort = metric in ['TAE', 'nRMSE', 'num_FP', 'num_FN']
            sorted_tools = sorted(tool_medians.keys(), 
                                key=lambda t: tool_medians[t], 
                                reverse=not reverse_sort)  # Higher is better by default
            tool_orders[(metric, row, col)] = sorted_tools
        
        # Add boxplots with overlaid jittered points for each tool
        # X-axis = tools (ranked by median performance per metric), Y-axis = metric values
        for metric, row, col, ylabel in metrics:
            # Get the sorted tool order for this metric
            sorted_tools = tool_orders.get((metric, row, col), tools)
            
            for tool in sorted_tools:
                tool_data = burden_data[burden_data['tool'] == tool]
                
                if len(tool_data) == 0:
                    continue
                
                # Group by sample to get per-sample mean across realizations
                grouped = tool_data.groupby('sample')[metric].mean().reset_index()
                
                if len(grouped) == 0:
                    continue
                
                # X-axis: tool name (repeated for each sample)
                # Y-axis: metric values for each sample
                x_data = [tool] * len(grouped)
                y_data = grouped[metric].values
                text_data = grouped['sample'].values  # Sample names for hover
                
                # Add ONE boxplot trace per tool per burden per subplot
                trace = go.Box(
                    x=x_data,
                    y=y_data,
                    name=tool,
                    marker=dict(
                        color=colors.get(tool, '#999'),
                        opacity=0.9,
                        line=dict(color='black', width=0.5),
                        size=4
                    ),
                    line=dict(color=colors.get(tool, '#999')),
                    fillcolor=colors.get(tool, '#999'),
                    opacity=0.5,
                    boxmean=False,
                    showlegend=bool(row == 1 and col == 1 and burden == mut_burdens[0]),  # Show legend only once per tool
                    legendgroup=tool,
                    offsetgroup=tool,  # Separate boxes by tool
                    boxpoints='all',  # Show all points
                    pointpos=0,  # Center points over box
                    jitter=0.3,  # Add horizontal jitter to points
                    text=text_data,
                    hovertemplate=f'<b>{tool}</b><br>Sample: %{{text}}<br>{ylabel}: %{{y:.3f}}<extra></extra>',
                    visible=bool(burden == mut_burdens[0]),  # Only first burden visible initially
                    width=0.6  # Increased box width from default (0.5) to 0.6
                )
                fig.add_trace(trace, row=row, col=col)
                all_traces.append(trace)
                trace_burdens.append(burden)  # Track burden for this trace
            
            # Update x-axis category order for this subplot to match sorted tools
            fig.update_xaxes(categoryorder='array', categoryarray=sorted_tools, row=row, col=col)
            # Update x-axis category order for this subplot to match sorted tools
            fig.update_xaxes(categoryorder='array', categoryarray=sorted_tools, row=row, col=col)
        
        # Add theoretical threshold lines for this burden (if available)
        if 'theoretical_min_TAE' in burden_data.columns:
            theoretical_tae = burden_data['theoretical_min_TAE'].iloc[0]
            if not np.isnan(theoretical_tae) and theoretical_tae > 0:
                # Add theoretical limits for relevant metrics
                for metric, row, col, ylabel in metrics:
                    if metric == 'TAE':
                        # Get the sorted tool order for TAE metric
                        sorted_tools_tae = tool_orders.get((metric, row, col), tools)
                        # Add horizontal line for minimum TAE (spans across all tools in sorted order)
                        trace = go.Scatter(
                            x=sorted_tools_tae,
                            y=[theoretical_tae] * len(sorted_tools_tae),
                            name='Theoretical Best',
                            mode='lines',
                            line=dict(color='#ef4444', width=2, dash='dash'),
                            showlegend=bool(burden == mut_burdens[0]),
                            legendgroup='threshold',
                            hovertemplate=f'Theoretical Min TAE: {theoretical_tae:.4f}<br>(Due to unknown signatures)<extra></extra>',
                            visible=bool(burden == mut_burdens[0])
                        )
                        fig.add_trace(trace, row=row, col=col)
                        all_traces.append(trace)
                        trace_burdens.append(burden)  # Track burden for threshold trace
    
    # Create dropdown menu for mutation burden selection
    buttons = []
    
    for burden_idx, burden in enumerate(mut_burdens):
        # Calculate which traces should be visible for this burden
        # Use explicit burden tracking instead of rigid indexing
        visible = [tb == burden for tb in trace_burdens]
        
        buttons.append(
            dict(
                label=f"{int(burden):,} mutations",
                method='update',
                args=[
                    {'visible': visible},
                    {'title': {'text': f'Performance per Method - {int(burden):,} Mutations'}}
                ]
            )
        )
    
    # Update y-axes labels
    fig.update_yaxes(title_text='F1 ↑', row=1, col=1)
    fig.update_yaxes(title_text='TAE ↓', row=1, col=2)
    fig.update_yaxes(title_text='MCC ↑', row=2, col=1)
    fig.update_yaxes(title_text='nRMSE ↓', row=2, col=2)
    fig.update_yaxes(title_text='FP Count', row=3, col=1)
    fig.update_yaxes(title_text='FN Count', row=3, col=2)
    
    # Update x-axes labels
    for col in range(1, 3):
        for row in range(1, 4):
            fig.update_xaxes(
                title_text='Method' if row == 3 else '',  # Only show "Method" label on bottom row
                tickangle=-90,  # Changed from -45 to -90 for vertical labels
                tickfont=dict(size=10),
                row=row, col=col
            )
    
    fig.update_layout(
        height=1000,  # Increased from 900 to 1000 to accommodate spacing and vertical labels
        title_text=f'Performance per Method - {int(mut_burdens[0]):,} Mutations',
        showlegend=True,
        boxmode='group',  # Group boxes by tool
        scattermode='overlay',  # Allow scatter points to overlay on boxplots
        updatemenus=[
            dict(
                buttons=buttons,
                direction='down',
                pad={'r': 10, 't': 10},
                showactive=True,
                x=1.15,
                xanchor='left',
                y=1.0,
                yanchor='top',
                bgcolor='white',
                bordercolor='#ccc',
                borderwidth=1
            )
        ],
        legend=dict(
            orientation='v',
            x=1.02,
            xanchor='left',
            y=0.5,
            yanchor='middle',
            itemdoubleclick='toggleothers'
        ),
        margin=dict(l=60, r=250, t=100, b=120),  # Increased bottom margin from 80 to 120 for vertical labels
        autosize=True
    )
    
    return fig



def create_per_method_sample_plots(df_samples, cancer_type, tools, colors):
    """Create scatter plots with jittered markers for sample-by-sample data.
    Shows one mutation burden at a time with a dropdown selector.
    
    Args:
        df_samples: DataFrame with per-sample metrics
        cancer_type: Cancer type to filter
        tools: List of tool names
        colors: Dict mapping tool names to colors
    """
    # Filter for this cancer type
    cancer_samples = df_samples[df_samples['cancer_type'] == cancer_type]
    
    if len(cancer_samples) == 0:
        return None
    
    # Get unique mutation burdens sorted
    mut_burdens = sorted(cancer_samples['num_muts'].unique())
    
    # Use 3x2 layout for metrics
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=['F1 ↑', 'TAE ↓', 'MCC ↑', 'nRMSE ↓', 'FP Sigs', 'FN Sigs'],
        vertical_spacing=0.18,
        horizontal_spacing=0.12
    )
    
    # Metrics configuration
    metrics = [
        ('F1', 1, 1, 'F1 ↑'),
        ('TAE', 1, 2, 'TAE ↓'),
        ('MCC', 2, 1, 'MCC ↑'),
        ('nRMSE', 2, 2, 'nRMSE ↓'),
        ('num_FP', 3, 1, 'FP Count'),
        ('num_FN', 3, 2, 'FN Count'),
    ]
    metrics = [(m, r, c, l) for m, r, c, l in metrics if m in cancer_samples.columns]
    
    if not metrics:
        print(f"  No valid metrics found for {cancer_type}")
        return None
    
    # Track all traces for visibility control
    all_traces = []
    trace_burdens = []
    
    # Create traces for EACH burden separately
    for burden in mut_burdens:
        burden_data = cancer_samples[cancer_samples['num_muts'] == burden]
        
        # Get unique samples for this burden
        samples = sorted(burden_data['sample'].unique())
        sample_labels = samples
        
        # Add colored markers for each method (with jitter)
        for metric, row, col, ylabel in metrics:
            for tool in tools:
                tool_data = burden_data[burden_data['tool'] == tool]
                
                if len(tool_data) == 0:
                    continue
                
                # Collect all data points for this tool across all samples
                x_data = []
                y_data = []
                text_data = []
                
                for sample_idx, sample in enumerate(samples):
                    sample_tool_data = tool_data[tool_data['sample'] == sample]
                    if len(sample_tool_data) == 0:
                        continue
                    
                    # Add horizontal jitter to spread points within each sample
                    jitter_amount = 0.3
                    x_jittered = [sample_idx + np.random.uniform(-jitter_amount, jitter_amount) 
                                 for _ in range(len(sample_tool_data))]
                    
                    x_data.extend(x_jittered)
                    y_data.extend(sample_tool_data[metric].values)
                    text_data.extend([f"{tool} - {sample}" for _ in range(len(sample_tool_data))])
                
                if len(x_data) == 0:
                    continue
                
                # Add colored scatter points for this tool (no boxplots)
                trace = go.Scatter(
                    x=x_data,
                    y=y_data,
                    name=tool,
                    mode='markers',
                    marker=dict(
                        color=colors.get(tool, '#999'),
                        size=6,
                        opacity=0.7,
                        line=dict(color='black', width=0.5)
                    ),
                    showlegend=bool(row == 1 and col == 1 and burden == mut_burdens[0]),
                    legendgroup=tool,
                    text=text_data,
                    hovertemplate=f'<b>%{{text}}</b><br>{ylabel}: %{{y:.3f}}<extra></extra>',
                    visible=bool(burden == mut_burdens[0])
                )
                fig.add_trace(trace, row=row, col=col)
                all_traces.append(trace)
                trace_burdens.append(burden)
            
            # Add vertical separator lines between samples
            for sample_idx in range(len(samples) - 1):
                # Line position is between sample_idx and sample_idx+1
                x_line = sample_idx + 0.5
                
                # Get y-axis range for this metric to draw full-height line
                metric_data = burden_data[metric].dropna()
                if len(metric_data) > 0:
                    y_min = metric_data.min()
                    y_max = metric_data.max()
                    y_range = y_max - y_min
                    # Extend slightly beyond data range
                    y_min_plot = y_min - 0.05 * y_range
                    y_max_plot = y_max + 0.05 * y_range
                else:
                    y_min_plot = 0
                    y_max_plot = 1
                
                separator_trace = go.Scatter(
                    x=[x_line, x_line],
                    y=[y_min_plot, y_max_plot],
                    mode='lines',
                    line=dict(color='lightgray', width=1, dash='dot'),
                    showlegend=False,
                    hoverinfo='skip',
                    visible=bool(burden == mut_burdens[0])
                )
                fig.add_trace(separator_trace, row=row, col=col)
                all_traces.append(separator_trace)
                trace_burdens.append(burden)
            
            # Update x-axis to show sample names
            fig.update_xaxes(
                tickmode='array',
                tickvals=list(range(len(samples))),
                ticktext=samples,
                row=row, col=col
            )
        
        # Add theoretical threshold lines for this burden (if available)
        if 'theoretical_min_TAE' in burden_data.columns:
            theoretical_tae = burden_data['theoretical_min_TAE'].iloc[0]
            if not np.isnan(theoretical_tae) and theoretical_tae > 0:
                for metric, row, col, ylabel in metrics:
                    if metric == 'TAE':
                        trace = go.Scatter(
                            x=list(range(len(sample_labels))),
                            y=[theoretical_tae] * len(sample_labels),
                            name='Theoretical Best',
                            mode='lines',
                            line=dict(color='#ef4444', width=2, dash='dash'),
                            showlegend=bool(burden == mut_burdens[0]),
                            legendgroup='threshold',
                            hovertemplate=f'Theoretical Min TAE: {theoretical_tae:.4f}<br>(Due to unknown signatures)<extra></extra>',
                            visible=bool(burden == mut_burdens[0])
                        )
                        fig.add_trace(trace, row=row, col=col)
                        all_traces.append(trace)
                        trace_burdens.append(burden)
    
    # Create dropdown menu for mutation burden selection
    buttons = []
    
    for burden_idx, burden in enumerate(mut_burdens):
        visible = [tb == burden for tb in trace_burdens]
        
        buttons.append(
            dict(
                label=f"{int(burden):,} mutations",
                method='update',
                args=[
                    {'visible': visible},
                    {'title': {'text': f'Sample-by-Sample Performance - {int(burden):,} Mutations'}}
                ]
            )
        )
    
    # Update y-axes labels
    fig.update_yaxes(title_text='F1 ↑', row=1, col=1)
    fig.update_yaxes(title_text='TAE ↓', row=1, col=2)
    fig.update_yaxes(title_text='MCC ↑', row=2, col=1)
    fig.update_yaxes(title_text='nRMSE ↓', row=2, col=2)
    fig.update_yaxes(title_text='FP Count', row=3, col=1)
    fig.update_yaxes(title_text='FN Count', row=3, col=2)
    
    # Update x-axes labels
    for col in range(1, 3):
        for row in range(1, 4):
            fig.update_xaxes(
                title_text='Sample' if row == 3 else '',
                tickangle=-90,
                tickfont=dict(size=8),
                row=row, col=col
            )
    
    fig.update_layout(
        height=900,
        title_text=f'Sample-by-Sample Performance - {int(mut_burdens[0]):,} Mutations',
        showlegend=True,
        scattermode='overlay',
        updatemenus=[
            dict(
                buttons=buttons,
                direction='down',
                pad={'r': 10, 't': 10},
                showactive=True,
                x=1.15,
                xanchor='left',
                y=1.0,
                yanchor='top',
                bgcolor='white',
                bordercolor='#ccc',
                borderwidth=1
            )
        ],
        legend=dict(
            orientation='v',
            x=1.02,
            xanchor='left',
            y=0.5,
            yanchor='middle',
            itemdoubleclick='toggleothers'
        ),
        margin=dict(l=60, r=250, t=80, b=100),
        autosize=True
    )
    
    return fig


def create_optimization_stability_plots(df_samples, cancer_type, tools, colors):
    """
    Create box plots showing distribution of per-sample metrics.
    This shows optimization stability - tighter distributions indicate more stable results.
    
    Args:
        df_samples: DataFrame with per-sample metrics
        cancer_type: Cancer type to filter
        tools: List of tool names
        colors: Dict mapping tool names to colors
        
    Returns:
        Plotly figure or None if no data
    """
    # Filter for this cancer type
    cancer_samples = df_samples[df_samples['cancer_type'] == cancer_type]
    
    if len(cancer_samples) == 0:
        return None
    
    # Get unique mutation burdens
    mut_burdens = sorted(cancer_samples['num_muts'].unique())

    # Metrics we care about for stability (per-sample)
    metrics_info = [
        ('F1', 1, 1, 'F1 ↑'),
        ('TAE', 1, 2, 'TAE ↓'),
        ('MCC', 1, 3, 'MCC ↑'),
        ('P', 2, 1, 'Precision ↑'),
        ('R', 2, 2, 'Recall ↑'),
        ('num_FP', 2, 3, 'FP Count ↓'),
    ]

    # If we have multiple realizations, compute per-sample variance across realizations
    if 'realization' in cancer_samples.columns and cancer_samples['realization'].nunique() > 1:
        # Group by tool, num_muts, sample and compute variance for each metric
        var_df = cancer_samples.groupby(['tool', 'num_muts', 'sample']).agg(
            {
                'F1': 'var',
                'TAE': 'var',
                'MCC': 'var',
                'P': 'var',
                'R': 'var',
                'num_FP': 'var'
            }
        ).reset_index()

        # Rename variance columns to indicate they are variances
        var_df = var_df.rename(columns={
            'F1': 'F1_var', 'TAE': 'TAE_var', 'MCC': 'MCC_var',
            'P': 'P_var', 'R': 'R_var', 'num_FP': 'num_FP_var'
        })

        # Create 2x3 grid of violin plots showing variance distributions per sample
        fig = make_subplots(rows=2, cols=3,
                            subplot_titles=['F1', 'TAE', 'MCC', 'Precision', 'Recall', 'FP Count'],
                            vertical_spacing=0.15, horizontal_spacing=0.10)

        # Track which tools have been shown in legend
        shown_tools = set()

        for tool in tools:
            tool_color = colors.get(tool, '#999')
            
            for i, burden in enumerate(mut_burdens):
                subset = var_df[(var_df['tool'] == tool) & (var_df['num_muts'] == burden)]
                if len(subset) == 0:
                    continue
                
                show_legend = tool not in shown_tools

                # For each metric, add violin of the variance values across samples
                for metric, row, col, ylabel in metrics_info:
                    var_col = f"{metric}_var"
                    if var_col not in subset.columns:
                        continue
                    fig.add_trace(
                        go.Violin(
                            y=subset[var_col].fillna(0),
                            x=[i] * len(subset),
                            name=tool,
                            legendgroup=tool,
                            line_color=tool_color,
                            fillcolor=tool_color,
                            opacity=0.6,
                            box_visible=True,
                            meanline_visible=True,
                            showlegend=bool(show_legend and (row == 1 and col == 1)),
                            meta={'tool': tool, 'mut_burden': int(burden)}
                        ),
                        row=row, col=col
                    )
                
                if show_legend:
                    shown_tools.add(tool)
        
        # Update x-axes to show mutation burdens
        for row, col in [(1,1), (1,2), (1,3), (2,1), (2,2), (2,3)]:
            fig.update_xaxes(
                tickmode='array',
                tickvals=list(range(len(mut_burdens))),
                ticktext=[f'{int(m):,}' for m in mut_burdens],
                title_text='Mutation Burden',
                row=row, col=col
            )

        # Update y-axes labels for variance plots
        fig.update_yaxes(title_text='Variance', row=1, col=1)
        fig.update_yaxes(title_text='Variance', row=1, col=2)
        fig.update_yaxes(title_text='Variance', row=1, col=3)
        fig.update_yaxes(title_text='Variance', row=2, col=1)
        fig.update_yaxes(title_text='Variance', row=2, col=2)
        fig.update_yaxes(title_text='Variance', row=2, col=3)

        fig.update_layout(
            height=700,
            showlegend=True,
            autosize=True
        )

        return fig

    # Fallback: if no realizations, keep previous box-based view (distribution of per-sample metrics)
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=['F1 Score', 'TAE (Total Absolute Error)', 'MCC',
                       'Precision', 'Recall', 'False Positives'],
        vertical_spacing=0.15,
        horizontal_spacing=0.10
    )

    for tool in tools:
        tool_color = colors.get(tool, '#999')

        for burden in mut_burdens:
            subset = cancer_samples[(cancer_samples['tool'] == tool) & 
                                   (cancer_samples['num_muts'] == burden)]

            if len(subset) == 0:
                continue

            burden_label = f'{int(burden):,}'
            x_label = f'{tool}<br>{burden_label}'
            show_legend = bool(burden == mut_burdens[0])  # Only show legend for first burden

            for metric, row, col, ylabel in metrics_info:
                fig.add_trace(
                    go.Box(
                        y=subset[metric],
                        name=f'{tool} ({burden_label})',
                        legendgroup=tool,
                        marker_color=tool_color,
                        showlegend=show_legend and (row == 1 and col == 1),
                        boxmean='sd',  # Show mean and std dev
                        boxpoints='all',  # Show all individual points
                        jitter=0.3,  # Add jitter to see overlapping points
                        pointpos=-1.8,  # Position points to the left of box
                        x=[x_label] * len(subset),
                        text=[f"Sample: {s}" for s in subset['sample']],
                        hovertemplate='%{text}<br>Value: %{y:.3f}<extra></extra>',
                        marker=dict(
                            size=4,
                            opacity=0.6,
                            line=dict(width=0.5, color='white')
                        )
                    ),
                    row=row, col=col
                )

    # Update axes
    for row in [1, 2]:
        for col in [1, 2, 3]:
            fig.update_xaxes(title_text='Tool / Mutation Burden', row=row, col=col)

    fig.update_yaxes(title_text='F1 ↑', row=1, col=1)
    fig.update_yaxes(title_text='TAE ↓', row=1, col=2)
    fig.update_yaxes(title_text='MCC ↑', row=1, col=3)
    fig.update_yaxes(title_text='Precision ↑', row=2, col=1)
    fig.update_yaxes(title_text='Recall ↑', row=2, col=2)
    fig.update_yaxes(title_text='FP Count ↓', row=2, col=3)

    fig.update_layout(
        height=700,
        title_text=f'Optimization Stability: Per-Sample Metric Distributions ({cancer_type})',
        showlegend=True,
        autosize=True
    )

    return fig


def create_realization_stability_plot(df_per_sample, cancer_type, tools, colors):
    """
    Create violin plots showing variance of per-sample metrics across realizations.
    
    This computes the variance for each sample across different realizations to assess 
    algorithm stability. Lower variance = more stable/reproducible algorithm.
    
    Args:
        df_per_sample: DataFrame with per-sample, per-realization results
        cancer_type: Cancer type to filter
        tools: List of tool names
        colors: Dict mapping tools to colors
    """
    # Filter for the cancer type
    df = df_per_sample[df_per_sample['cancer_type'] == cancer_type].copy()
    
    # Clean tool names (remove _w0, _w1 suffixes if present)
    if 'tool' in df.columns:
        df['tool'] = df['tool'].str.replace(r'_w\d+$', '', regex=True)
    
    # Check if we have realization data
    if 'realization' not in df.columns:
        print(f"  No realization column found for {cancer_type}")
        return None
    
    if 'sample' not in df.columns:
        print(f"  No sample column found for {cancer_type} - cannot track individual samples")
        return None
    
    num_realizations = df['realization'].nunique()
    if num_realizations < 2:
        print(f"  Only {num_realizations} realization(s) found for {cancer_type}")
        return None
    
    print(f"  Computing per-sample variance across {num_realizations} realizations...")
    
    # Compute variance for each sample across realizations
    # Only include columns that exist in the dataframe
    agg_dict = {}
    for col in ['F1', 'TAE', 'MCC', 'nRMSE', 'num_FP', 'num_FN']:
        if col in df.columns:
            agg_dict[col] = 'var'
    
    if not agg_dict:
        print(f"  No variance metrics available for {cancer_type}")
        return None
    
    var_df = df.groupby(['tool', 'num_muts', 'sample']).agg(agg_dict).reset_index()
    
    if len(var_df) == 0:
        print(f"  No variance data computed for {cancer_type}")
        return None
    
    # Get unique mutation burdens sorted
    mut_burdens = sorted(var_df['num_muts'].unique())
    
    # Use 3x2 layout (3 rows, 2 columns) for better visibility
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=['F1 Variance ↓', 'TAE Variance ↓', 'MCC Variance ↓', 
                        'nRMSE Variance ↓', 'FP Variance', 'FN Variance'],
        vertical_spacing=0.12,
        horizontal_spacing=0.12
    )
    
    # Metrics to plot - update positions for 3x2 layout
    all_metrics = [
        ('F1', 1, 1, 'F1 Variance ↓'),
        ('TAE', 1, 2, 'TAE Variance ↓'),
        ('MCC', 2, 1, 'MCC Variance ↓'),
        ('nRMSE', 2, 2, 'nRMSE Variance ↓'),
        ('num_FP', 3, 1, 'FP Variance'),
        ('num_FN', 3, 2, 'FN Variance'),
    ]
    metrics = [(m, r, c, l) for m, r, c, l in all_metrics if m in var_df.columns]
    
    # Get all unique mutation burdens and format them consistently
    mut_burden_labels = [f"{int(m):,}" for m in sorted(mut_burdens)]
    
    # Create box plots - iterate through tools AND mutation burdens for proper filtering
    for tool in tools:
        for mut_burden in mut_burdens:
            tool_mut_data = var_df[(var_df['tool'] == tool) & (var_df['num_muts'] == mut_burden)]
            
            if len(tool_mut_data) == 0:
                continue
            
            # Box plots for each metric showing individual points
            for metric, row, col, ylabel in metrics:
                if metric not in tool_mut_data.columns:
                    continue
                
                # Convert mutation burden to categorical string
                x_label = f"{int(mut_burden):,}"
                
                fig.add_trace(
                    go.Box(
                        x=[x_label] * len(tool_mut_data),
                        y=tool_mut_data[metric],
                        name=tool,
                        marker=dict(
                            color=colors.get(tool, '#999'),
                            size=6,
                            opacity=0.6,
                            line=dict(width=0.5, color='white')
                        ),
                        boxpoints='all',  # Show all individual points
                        pointpos=0,  # Center points
                        jitter=0.3,  # Spread points horizontally
                        showlegend=bool(row == 1 and col == 1 and mut_burden == mut_burdens[0]),  # Only show legend once per tool
                        legendgroup=tool,
                        line=dict(color=colors.get(tool, '#999'), width=2),
                        fillcolor=colors.get(tool, '#999'),
                        opacity=0.5,
                        offsetgroup=tool,  # Critical for grouping boxes side-by-side
                        text=[f"Tool: {tool}<br>Sample: {s}<br>{int(mut_burden):,} muts<br>Variance: {v:.4f}" 
                              for s, v in zip(tool_mut_data['sample'], tool_mut_data[metric])],
                        hovertemplate='%{text}<extra></extra>',
                        visible=True,
                        meta={'tool': tool, 'metric': metric, 'mut_burden': str(int(mut_burden))}
                    ),
                    row=row, col=col
                )
    
    # Update y-axes labels with metric names and variance indication
    fig.update_yaxes(title_text='F1 Variance ↓', row=1, col=1)
    fig.update_yaxes(title_text='TAE Variance ↓', row=1, col=2)
    fig.update_yaxes(title_text='MCC Variance ↓', row=2, col=1)
    fig.update_yaxes(title_text='nRMSE Variance ↓', row=2, col=2)
    fig.update_yaxes(title_text='FP Variance', row=3, col=1)
    fig.update_yaxes(title_text='FN Variance', row=3, col=2)
    
    # Update x-axes to show mutation burden
    fig.update_xaxes(title_text='Mutation Burden', row=1, col=1, type='category')
    fig.update_xaxes(title_text='Mutation Burden', row=1, col=2, type='category')
    fig.update_xaxes(title_text='Mutation Burden', row=2, col=1, type='category')
    fig.update_xaxes(title_text='Mutation Burden', row=2, col=2, type='category')
    fig.update_xaxes(title_text='Mutation Burden', row=3, col=1, type='category')
    fig.update_xaxes(title_text='Mutation Burden', row=3, col=2, type='category')
    
    fig.update_layout(
        height=900,  # Increased height for 3 rows
        title_text=f'Algorithm Stability (Variance Across Realizations)<br><sub>All mutation burdens - Lower variance = more stable/reproducible</sub>',
        showlegend=True,
        boxmode='group',  # Group boxes side-by-side instead of overlapping
        legend=dict(
            orientation='v',
            x=1.02,
            xanchor='left',
            y=1,
            yanchor='top'
        ),
        margin=dict(l=60, r=120, t=100, b=60),
        autosize=True
    )

    return fig


def create_precision_recall_subplots(df, cancer_type, tools, colors, df_per_sample=None):
    """
    Create Precision-Recall plots (one per mutation burden).
    
    Args:
        df: DataFrame with evaluation metrics (aggregated by cancer type)
        cancer_type: Cancer type name (for display; if 'cancer_type' column not present, uses df as-is)
        tools: List of tool names
        colors: Dict mapping tool names to colors
        df_per_sample: DataFrame with per-sample data for computing variance and theoretical thresholds
    """
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    
    # Filter by cancer type if the column exists, otherwise use the dataframe as-is
    if 'cancer_type' in df.columns:
        cancer_df = df[df['cancer_type'] == cancer_type]
    else:
        cancer_df = df
    
    # Get unique mutation burdens
    mut_burdens = sorted(cancer_df['num_muts'].unique())
    n_burdens = len(mut_burdens)
    
    # Determine layout: up to 4 columns
    n_cols = min(n_burdens, 4)
    n_rows = (n_burdens + n_cols - 1) // n_cols  # Ceiling division
    
    # Create subplot titles
    subplot_titles = [f'{int(burden):,} muts' for burden in mut_burdens]
    
    # Create specs for subplots
    specs = []
    for row in range(n_rows):
        row_spec = []
        for col in range(n_cols):
            idx = row * n_cols + col
            if idx < n_burdens:
                row_spec.append({})
            else:
                row_spec.append(None)  # Empty subplot
        specs.append(row_spec)
    
    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        subplot_titles=subplot_titles,
        specs=specs,
        vertical_spacing=0.12,
        horizontal_spacing=0.10
    )
    
    # ===== Precision-Recall Plots per Mutation Burden =====
    method_symbols = get_method_symbols()
    
    for burden_idx, burden in enumerate(mut_burdens):
        # Calculate subplot position
        pr_row = (burden_idx // n_cols) + 1
        pr_col = (burden_idx % n_cols) + 1
        
        # Filter data for this burden
        burden_data = cancer_df[cancer_df['num_muts'] == burden]
        
        # Plot each tool
        for tool in tools:
            tool_data = burden_data[burden_data['tool'] == tool]
            
            if len(tool_data) == 0:
                continue
            
            # Get base method to determine symbol
            base_method = get_tool_base_method(tool)
            symbol = method_symbols.get(base_method, 'circle')
            
            # Build custom hover text
            hover_text_items = []
            for _, row in tool_data.iterrows():
                items = [f"Mutations: {int(row['num_muts']):,}"]
                if 'optimal_penalty' in row and pd.notna(row['optimal_penalty']):
                    items.append(f"Penalty: {row['optimal_penalty']}")
                hover_text_items.append('<br>'.join(items))
                
            fig.add_trace(go.Scatter(
                x=tool_data['R_mean'],
                y=tool_data['P_mean'],
                mode='markers',
                name=tool,
                marker=dict(
                    size=10,  # Fixed size since we're showing one burden per plot
                    color=colors.get(tool, '#999'),
                    symbol=symbol,
                    opacity=0.7,
                    line=dict(width=1, color='black')
                ),
                error_x=dict(array=tool_data['R_sem'], type='data'),
                error_y=dict(array=tool_data['P_sem'], type='data'),
                text=hover_text_items,
                hovertemplate='<b>%{fullData.name}</b><br>Recall: %{x:.3f}<br>Precision: %{y:.3f}<br>%{text}<extra></extra>',
                showlegend=(burden_idx == 0),  # Show legend only for first subplot
                legendgroup=tool
            ), row=pr_row, col=pr_col)
        
        # Add F1 score reference lines to each P-R plot
        for f1_val, f1_label in [(0.5, 'F1=0.5'), (0.7, 'F1=0.7'), (0.9, 'F1=0.9')]:
            recall = np.linspace(0.01, 1, 100)
            precision = f1_val * recall / (2 * recall - f1_val * recall)
            precision = np.clip(precision, 0, 1)
            
            fig.add_trace(go.Scatter(
                x=recall,
                y=precision,
                mode='lines',
                name=f1_label,
                line=dict(dash='dash', color='gray', width=0.8),
                opacity=0.4,
                hoverinfo='name',
                showlegend=(burden_idx == 0),  # Show legend only for first subplot
                legendgroup='f1_lines'
            ), row=pr_row, col=pr_col)
        
        # Add theoretical threshold lines if available (per burden)
        if df_per_sample is not None:
            # Filter per-sample data by cancer type if the column exists
            if 'cancer_type' in df_per_sample.columns:
                cancer_per_sample = df_per_sample[df_per_sample['cancer_type'] == cancer_type]
            else:
                cancer_per_sample = df_per_sample
            burden_per_sample = cancer_per_sample[cancer_per_sample['num_muts'] == burden]
            
            if len(burden_per_sample) > 0:
                if 'theoretical_max_Recall' in burden_per_sample.columns:
                    theoretical_recall = burden_per_sample['theoretical_max_Recall'].mean()
                    if not np.isnan(theoretical_recall) and theoretical_recall < 1.0:
                        fig.add_trace(go.Scatter(
                            x=[theoretical_recall, theoretical_recall],
                            y=[0, 1],
                            mode='lines',
                            name='Max Recall',
                            line=dict(color='red', width=1.5, dash='dash'),
                            showlegend=(burden_idx == 0),
                            legendgroup='thresholds',
                            hovertemplate=f'Theoretical Max Recall: {theoretical_recall:.3f}<br>(Unknown sigs counted as FN)<extra></extra>'
                        ), row=pr_row, col=pr_col)

                if 'theoretical_max_Precision' in burden_per_sample.columns:
                    theoretical_precision = burden_per_sample['theoretical_max_Precision'].mean()
                    if not np.isnan(theoretical_precision) and theoretical_precision < 1.0:
                        fig.add_trace(go.Scatter(
                            x=[0, 1],
                            y=[theoretical_precision, theoretical_precision],
                            mode='lines',
                            name='Max Precision',
                            line=dict(color='blue', width=1.5, dash='dot'),
                            showlegend=(burden_idx == 0),
                            legendgroup='thresholds',
                            hovertemplate=f'Theoretical Max Precision: {theoretical_precision:.3f}<extra></extra>'
                        ), row=pr_row, col=pr_col)
        
        # Update P-R plot axes
        fig.update_xaxes(title_text='Recall', range=[0, 1.05], row=pr_row, col=pr_col)
        fig.update_yaxes(title_text='Precision', range=[0, 1.05], row=pr_row, col=pr_col)
    
    # Calculate total height based on number of rows
    total_height = 400 * n_rows
    
    fig.update_layout(
        height=total_height,
        title_text=f'{cancer_type} - Precision-Recall Trade-off<br><sub>Error bars: ±1 SEM across samples</sub>',
        showlegend=True,
        legend=dict(x=1.05, y=1, itemdoubleclick='toggleothers'),
        hovermode='closest',
        autosize=True
    )
    
    return fig


def create_metrics_grid(df, cancer_type, tools, colors, df_per_sample=None):
    """
    Create 2x3 grid of performance metrics vs mutation count.
    Shows mean across all samples/realizations with variance as error bars.
    
    Args:
        df: DataFrame with evaluation metrics
        cancer_type: Cancer type name (for display; if 'cancer_type' column not present, uses df as-is)
        tools: List of tool names
        colors: Dict mapping tool names to colors
        df_per_sample: DataFrame with per-sample metrics for error bars
    """
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    
    # Filter by cancer type if the column exists, otherwise use the dataframe as-is
    if 'cancer_type' in df.columns:
        cancer_df = df[df['cancer_type'] == cancer_type]
    else:
        cancer_df = df
    
    # Create 2x3 subplot grid
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=['F1 ↑', 'TAE ↓', 'MCC ↑', 'nRMSE ↓', 'FP Sigs', 'FN Sigs'],
        vertical_spacing=0.15,
        horizontal_spacing=0.12
    )
    
    metric_mapping = [
        ('F1', ['F1', 'F1_mean'], 1, 1),
        ('TAE', ['TAE', 'TAE_mean'], 1, 2),
        ('MCC', ['MCC', 'MCC_mean'], 1, 3),
        ('nRMSE', ['nRMSE', 'nRMSE_mean'], 2, 1),
        ('num_FP', ['num_FP_sigs', 'num_FP'], 2, 2),
        ('num_FN', ['num_FN_sigs', 'num_FN'], 2, 3),
    ]
    
    # Use per-sample data if available for proper error bars
    use_per_sample = df_per_sample is not None and len(df_per_sample) > 0
    if use_per_sample:
        # Filter per-sample data by cancer type if the column exists
        if 'cancer_type' in df_per_sample.columns:
            cancer_per_sample = df_per_sample[df_per_sample['cancer_type'] == cancer_type]
        else:
            cancer_per_sample = df_per_sample
    
    for tool in tools:
        tool_data = cancer_df[cancer_df['tool'] == tool]
        
        if len(tool_data) == 0:
            continue
        
        for metric_name, possible_cols, row, col in metric_mapping:
            # Use per-sample data for aggregation if available
            if use_per_sample and metric_name in cancer_per_sample.columns:
                tool_per_sample = cancer_per_sample[cancer_per_sample['tool'] == tool]
                
                if len(tool_per_sample) == 0:
                    continue
                
                # Group by num_muts and compute statistics
                grouped = tool_per_sample.groupby('num_muts')[metric_name].agg(['mean', 'std', 'count']).reset_index()
                grouped['sem'] = grouped['std'] / np.sqrt(grouped['count'])
                
                error_y = dict(
                    array=grouped['sem'],
                    type='data',
                    visible=True,
                    thickness=1.5,
                    width=3
                )
                
                fig.add_trace(
                    go.Scatter(
                        x=grouped['num_muts'],
                        y=grouped['mean'],
                        mode='lines+markers',
                        name=tool,
                        line=dict(color=colors.get(tool, '#999'), width=2),
                        marker=dict(size=8, color=colors.get(tool, '#999')),
                        error_y=error_y,
                        showlegend=(row == 1 and col == 1),  # Show legend only once
                        legendgroup=tool,
                        text=[f"{int(m):,} muts<br>n={int(n)}<br>Std: {s:.3f}" 
                              for m, n, s in zip(grouped['num_muts'], grouped['count'], grouped['std'])],
                        hovertemplate='<b>%{fullData.name}</b><br>%{text}<br>Mean: %{y:.3f}<extra></extra>'
                    ),
                    row=row, col=col
                )
            else:
                # Fallback: use aggregated data
                metric_col = None
                for col_name in possible_cols:
                    if col_name in tool_data.columns:
                        metric_col = col_name
                        break
                
                if metric_col is None:
                    continue
                
                fig.add_trace(
                    go.Scatter(
                        x=tool_data['num_muts'],
                        y=tool_data[metric_col],
                        mode='lines+markers',
                        name=tool,
                        line=dict(color=colors.get(tool, '#999'), width=2),
                        marker=dict(size=8, color=colors.get(tool, '#999')),
                        showlegend=(row == 1 and col == 1),  # Show legend only once
                        legendgroup=tool,
                        hovertemplate='<b>%{fullData.name}</b><br>Mean: %{y:.3f}<extra></extra>'
                    ),
                    row=row, col=col
                )
    
    # Add theoretical threshold lines for metrics
    threshold_source = df_per_sample if df_per_sample is not None else cancer_df
    # Filter threshold source by cancer type if the column exists
    if 'cancer_type' in threshold_source.columns:
        cancer_threshold = threshold_source[threshold_source['cancer_type'] == cancer_type]
    else:
        cancer_threshold = threshold_source
    
    if len(cancer_threshold) > 0:
        all_muts = sorted(cancer_df['num_muts'].unique())
        
        # TAE threshold
        if 'theoretical_min_TAE' in cancer_threshold.columns:
            theoretical_tae = cancer_threshold['theoretical_min_TAE'].iloc[0]
            if not np.isnan(theoretical_tae) and theoretical_tae > 0:
                fig.add_trace(
                    go.Scatter(
                        x=all_muts,
                        y=[theoretical_tae] * len(all_muts),
                        mode='lines',
                        name='Theoretical Limit',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=True,
                        legendgroup='threshold',
                        hovertemplate=f'Theoretical Min TAE: {theoretical_tae:.3f}<extra></extra>'
                    ),
                    row=1, col=2  # TAE subplot
                )
        
        # F1 threshold
        if 'theoretical_max_F1' in cancer_threshold.columns:
            avg_f1 = cancer_threshold['theoretical_max_F1'].mean()
            if not np.isnan(avg_f1) and avg_f1 < 1.0:
                fig.add_trace(
                    go.Scatter(
                        x=all_muts,
                        y=[avg_f1] * len(all_muts),
                        mode='lines',
                        name='Theoretical Limit',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=False,
                        legendgroup='threshold',
                        hovertemplate=f'Theoretical Max F1: {avg_f1:.3f}<extra></extra>'
                    ),
                    row=1, col=1  # F1 subplot
                )
        
        # MCC threshold
        if 'theoretical_max_MCC' in cancer_threshold.columns:
            avg_mcc = cancer_threshold['theoretical_max_MCC'].mean()
            if not np.isnan(avg_mcc) and avg_mcc < 1.0:
                fig.add_trace(
                    go.Scatter(
                        x=all_muts,
                        y=[avg_mcc] * len(all_muts),
                        mode='lines',
                        name='Theoretical Limit',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=False,
                        legendgroup='threshold',
                        hovertemplate=f'Theoretical Max MCC: {avg_mcc:.3f}<extra></extra>'
                    ),
                    row=1, col=3  # MCC subplot
                )
    
    # Update x-axes for metrics to log scale
    for col in [1, 2, 3]:
        fig.update_xaxes(type='log', title_text='Mutations', row=1, col=col)
        fig.update_xaxes(type='log', title_text='Mutations', row=2, col=col)
    
    fig.update_layout(
        height=700,
        title_text=f'{cancer_type} - Performance Metrics<br><sub>Mean ± SEM across all samples vs mutation burden</sub>',
        showlegend=True,
        legend=dict(x=1.05, y=1, itemdoubleclick='toggleothers'),
        hovermode='closest',
        autosize=True
    )
    
    return fig


def create_old_metrics_grid_UNUSED(df, cancer_type, tools, colors, df_per_sample=None):
    """
    Create 2x3 grid of performance metrics vs mutation count.
    Shows mean across all samples/realizations with variance as error bars.
    
    Args:
        df: DataFrame with evaluation metrics (aggregated by cancer type)
        cancer_type: Cancer type name (for display)
        tools: List of tool names
        colors: Dict mapping tool names to colors
        df_per_sample: DataFrame with per-sample data for computing variance and theoretical thresholds
    """
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=['F1 ↑', 'TAE ↓', 'MCC ↑',
                       'nRMSE ↓', 'FP Sigs', 'FN Sigs'],
        vertical_spacing=0.15,
        horizontal_spacing=0.1
    )
    
    # Map metric display names to potential column names in the data
    metric_mapping = [
        ('F1', ['F1', 'F1_mean'], 1, 1),
        ('TAE', ['TAE', 'TAE_mean'], 1, 2),
        ('MCC', ['MCC', 'MCC_mean'], 1, 3),
        ('nRMSE', ['nRMSE', 'nRMSE_mean'], 2, 1),
        ('num_FP', ['num_FP_sigs', 'num_FP'], 2, 2),
        ('num_FN', ['num_FN_sigs', 'num_FN'], 2, 3),
    ]
    
    # Use per-sample data if available for proper error bars
    use_per_sample = df_per_sample is not None and len(df_per_sample) > 0
    if use_per_sample:
        cancer_per_sample = df_per_sample[df_per_sample['cancer_type'] == cancer_type]
    
    for tool in tools:
        tool_data = df[df['tool'] == tool]
        
        if len(tool_data) == 0:
            continue
        
        for metric_name, possible_cols, row, col in metric_mapping:
            # Use per-sample data for aggregation if available
            if use_per_sample and metric_name in cancer_per_sample.columns:
                tool_per_sample = cancer_per_sample[cancer_per_sample['tool'] == tool]
                
                if len(tool_per_sample) == 0:
                    continue
                
                # Group by num_muts and compute statistics across all samples×realizations
                grouped = tool_per_sample.groupby('num_muts')[metric_name].agg(['mean', 'std', 'count']).reset_index()
                
                # Use SEM for error bars
                grouped['sem'] = grouped['std'] / np.sqrt(grouped['count'])
                
                error_y = dict(
                    array=grouped['sem'],
                    type='data',
                    visible=True,
                    thickness=1.5,
                    width=3
                )
                
                fig.add_trace(
                    go.Scatter(
                        x=grouped['num_muts'],
                        y=grouped['mean'],
                        mode='lines+markers',
                        name=tool,
                        line=dict(color=colors.get(tool, '#999'), width=2),
                        marker=dict(size=8, color=colors.get(tool, '#999')),
                        error_y=error_y,
                        showlegend=(row == 1 and col == 1),
                        legendgroup=tool,
                        text=[f"{int(m):,} muts<br>n={int(n)}<br>Std: {s:.3f}" 
                              for m, n, s in zip(grouped['num_muts'], grouped['count'], grouped['std'])],
                        hovertemplate='<b>%{fullData.name}</b><br>%{text}<br>Mean: %{y:.3f}<extra></extra>'
                    ),
                    row=row, col=col
                )
            else:
                # Fallback: use aggregated data (no error bars)
                metric_col = None
                for col_name in possible_cols:
                    if col_name in tool_data.columns:
                        metric_col = col_name
                        break
                
                if metric_col is None:
                    continue
                
                fig.add_trace(
                    go.Scatter(
                        x=tool_data['num_muts'],
                        y=tool_data[metric_col],
                        mode='lines+markers',
                        name=tool,
                        line=dict(color=colors.get(tool, '#999'), width=2),
                        marker=dict(size=8, color=colors.get(tool, '#999')),
                        showlegend=(row == 1 and col == 1),
                        legendgroup=tool,
                        hovertemplate='<b>%{fullData.name}</b><br>Mean: %{y:.3f}<extra></extra>'
                    ),
                    row=row, col=col
                )
    
    # Add theoretical threshold lines if available
    # Try per-sample data first, then fall back to aggregated data
    threshold_source = df_per_sample if df_per_sample is not None else df
    cancer_threshold = threshold_source[threshold_source['cancer_type'] == cancer_type]
    
    if len(cancer_threshold) > 0:
        # Get mutation burdens for x-axis
        all_muts = sorted(df['num_muts'].unique())
        
        # TAE threshold (global value - same for all samples)
        if 'theoretical_min_TAE' in cancer_threshold.columns:
            theoretical_tae = cancer_threshold['theoretical_min_TAE'].iloc[0]
            if not np.isnan(theoretical_tae) and theoretical_tae > 0:
                fig.add_trace(
                    go.Scatter(
                        x=all_muts,
                        y=[theoretical_tae] * len(all_muts),
                        mode='lines',
                        name='Theoretical Limit',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=True,
                        legendgroup='threshold',
                        hovertemplate=f'Theoretical Min TAE: {theoretical_tae:.3f}<extra></extra>'
                    ),
                    row=1, col=2  # TAE subplot
                )
        
        # F1 threshold
        if 'theoretical_max_F1' in cancer_threshold.columns:
            # Average across all samples if using per-sample data
            avg_f1 = cancer_threshold['theoretical_max_F1'].mean()
            if not np.isnan(avg_f1) and avg_f1 < 1.0:
                fig.add_trace(
                    go.Scatter(
                        x=all_muts,
                        y=[avg_f1] * len(all_muts),
                        mode='lines',
                        name='Theoretical Limit',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=False,
                        legendgroup='threshold',
                        hovertemplate=f'Theoretical Max F1: {avg_f1:.3f}<extra></extra>'
                    ),
                    row=1, col=1  # F1 subplot
                )
        
        # MCC threshold
        if 'theoretical_max_MCC' in cancer_threshold.columns:
            # Average across all samples if using per-sample data
            avg_mcc = cancer_threshold['theoretical_max_MCC'].mean()
            if not np.isnan(avg_mcc) and avg_mcc < 1.0:
                fig.add_trace(
                    go.Scatter(
                        x=all_muts,
                        y=[avg_mcc] * len(all_muts),
                        mode='lines',
                        name='Theoretical Limit',
                        line=dict(color='red', width=2, dash='dash'),
                        showlegend=False,
                        legendgroup='threshold',
                        hovertemplate=f'Theoretical Max MCC: {avg_mcc:.3f}<extra></extra>'
                    ),
                    row=1, col=3  # MCC subplot
                )
        
        # Note: No theoretical thresholds for nRMSE, FP, or FN
        # - nRMSE: theoretical minimum is 0 (perfect reconstruction)
        # - FP: theoretical minimum is 0 (no false positives)
        # - FN: theoretical minimum is 0 (no false negatives)
        # These are trivial limits and don't provide useful context
    
    # Update x-axes to log scale
    for col in range(1, 4):
        fig.update_xaxes(type='log', title_text='Mutations', row=1, col=col)
        fig.update_xaxes(type='log', title_text='Mutations', row=2, col=col)
    
    fig.update_layout(
        height=700,
        title_text=f'Performance Metrics vs Mutation Burden<br><sub>Mean ± SEM across all samples and realizations</sub>',
        showlegend=True,
        legend=dict(x=1.05, y=1, itemdoubleclick='toggleothers'),
        autosize=True
    )
    
    return fig


def create_aggregate_plot(df, tools, colors):
    """Create aggregate precision-recall plot across all cancer types."""
    # Group by tool and num_muts, calculate mean, std, and count for SEM
    agg_df = df.groupby(['tool', 'num_muts']).agg({
        'P_mean': ['mean', 'std', 'count'],
        'R_mean': ['mean', 'std', 'count'],
        'F1_mean': ['mean', 'std', 'count']
    }).reset_index()
    
    agg_df.columns = ['tool', 'num_muts', 'P_mean', 'P_std', 'P_count', 'R_mean', 'R_std', 'R_count', 'F1_mean', 'F1_std', 'F1_count']
    
    # Calculate SEM
    agg_df['P_sem'] = agg_df['P_std'] / np.sqrt(agg_df['P_count'])
    agg_df['R_sem'] = agg_df['R_std'] / np.sqrt(agg_df['R_count'])
    
    fig = go.Figure()
    
    method_symbols = get_method_symbols()
    
    for tool in tools:
        tool_data = agg_df[agg_df['tool'] == tool]
        
        if len(tool_data) == 0:
            continue
        
        # Get base method to determine symbol
        base_method = get_tool_base_method(tool)
        symbol = method_symbols.get(base_method, 'circle')
        
        fig.add_trace(go.Scatter(
            x=tool_data['R_mean'],
            y=tool_data['P_mean'],
            mode='markers',
            name=tool,
            marker=dict(
                size=12,
                color=colors.get(tool, '#999'),
                symbol=symbol,
                line=dict(width=1, color='black')
            ),
            error_x=dict(array=tool_data['R_sem'], type='data'),
            error_y=dict(array=tool_data['P_sem'], type='data'),
            text=[f"Mutations: {m}" for m in tool_data['num_muts']],
            hovertemplate='<b>%{fullData.name}</b><br>Recall: %{x:.3f}<br>Precision: %{y:.3f}<br>%{text}<extra></extra>'
        ))
    
    fig.update_layout(
        title=dict(
            text='Aggregate Precision-Recall (Across All Cancer Types)<br><sub>Error bars: ±1 SEM (Standard Error of Mean)</sub>',
            font=dict(size=16)
        ),
        xaxis_title='R (Recall) ↑',
        yaxis_title='P (Precision) ↑',
        xaxis=dict(range=[0, 1.05]),
        yaxis=dict(range=[0, 1.05]),
        height=500,
        hovermode='closest',
        autosize=True
    )
    
    return fig


def create_bootstrap_zscore_histograms(df_bootstrap, cancer_type, tools, colors):
    """
    Create barplots showing distribution of z-scores for each mutation burden.
    Values >3 are binned into a ">3" category.
    
    Args:
        df_bootstrap: DataFrame with bootstrap quality assessment results
        cancer_type: Cancer type to filter
        tools: List of tool names
        colors: Dict mapping tool names to colors
        
    Returns:
        Dict mapping mutation burden to Plotly figure, or None if no data
    """
    cancer_data = df_bootstrap[df_bootstrap['cancer_type'] == cancer_type]
    
    if len(cancer_data) == 0:
        return None
    
    burdens = sorted(cancer_data['num_muts'].unique())
    
    # Create a separate figure for each mutation burden
    figures = {}
    
    for burden in burdens:
        burden_data = cancer_data[cancer_data['num_muts'] == burden]
        
        # Create subplots for the three z-score metrics
        fig = make_subplots(
            rows=1, cols=3,
            subplot_titles=('Z-score L1', 'Z-score L2', 'Z-score Cosine Similarity'),
            horizontal_spacing=0.12
        )
        
        # Define bins: <-3, -3 to -2, -2 to -1, -1 to 0, 0 to 1, 1 to 2, 2 to 3, >3
        bin_edges = [-np.inf, -3, -2, -1, 0, 1, 2, 3, np.inf]
        bin_labels = ['<-3', '[-3,-2)', '[-2,-1)', '[-1,0)', '[0,1)', '[1,2)', '[2,3)', '>3']
        
        for col_idx, metric in enumerate(['z_L1', 'z_L2', 'z_cos'], start=1):
            # Calculate histogram for each tool
            for tool in tools:
                tool_data = burden_data[burden_data['tool'] == tool]
                
                if len(tool_data) == 0:
                    continue
                
                tool_color = colors.get(tool, '#999')
                
                # Bin the z-scores
                z_values = tool_data[metric].values
                counts, _ = np.histogram(z_values, bins=bin_edges)
                
                # Create bar plot
                fig.add_trace(
                    go.Bar(
                        x=bin_labels,
                        y=counts,
                        name=tool,
                        marker_color=tool_color,
                        showlegend=(col_idx == 1),  # Only show legend in first subplot
                        legendgroup=tool,
                        hovertemplate=f'{tool}<br>Bin: %{{x}}<br>Count: %{{y}}<extra></extra>'
                    ),
                    row=1, col=col_idx
                )
        
        # Update axes
        for col in [1, 2, 3]:
            fig.update_xaxes(title_text="Z-score Bin", row=1, col=col)
        
        fig.update_yaxes(title_text="Number of Samples", row=1, col=1)
        fig.update_yaxes(title_text="Number of Samples", row=1, col=2)
        fig.update_yaxes(title_text="Number of Samples", row=1, col=3)
        
        burden_label = f'{int(burden):,}'
        fig.update_layout(
            height=450,
            title_text=f'Bootstrap Z-score Distribution - {burden_label} mutations ({cancer_type})<br>'
                      f'<sub>Bars show how many samples fall in each z-score range. Lower |z| indicates better fit.</sub>',
            showlegend=True,
            barmode='group',  # Group bars by tool
            autosize=True
        )
        
        figures[burden] = fig
    
    return figures


def create_bootstrap_summary_by_burden(df_bootstrap, cancer_type, tools, colors):
    """
    Create separate bar plots for each mutation burden showing z-score distributions.
    
    Args:
        df_bootstrap: DataFrame with bootstrap quality assessment results
        cancer_type: Cancer type to filter
        tools: List of tool names
        colors: Dict mapping tool names to colors
        
    Returns:
        Dict mapping mutation burden to Plotly figure, or None if no data
    """
    cancer_data = df_bootstrap[df_bootstrap['cancer_type'] == cancer_type]
    
    if len(cancer_data) == 0:
        return None
    
    burdens = sorted(cancer_data['num_muts'].unique())
    
    # Create a separate figure for each mutation burden
    figures = {}
    
    for burden in burdens:
        burden_data = cancer_data[cancer_data['num_muts'] == burden]
        
        # Create subplots for the three z-score metrics
        fig = make_subplots(
            rows=1, cols=3,
            subplot_titles=('Z-score L1', 'Z-score L2', 'Z-score Cosine Similarity'),
            horizontal_spacing=0.12
        )
        
        # Calculate mean z-scores for each tool
        for tool in tools:
            tool_data = burden_data[burden_data['tool'] == tool]
            
            if len(tool_data) == 0:
                continue
            
            tool_color = colors.get(tool, '#999')
            
            # Calculate mean and std for each metric
            z_L1_mean = tool_data['z_L1'].mean()
            z_L1_std = tool_data['z_L1'].std()
            z_L2_mean = tool_data['z_L2'].mean()
            z_L2_std = tool_data['z_L2'].std()
            z_cos_mean = tool_data['z_cos'].mean()
            z_cos_std = tool_data['z_cos'].std()
            
            # Z-score L1 bar
            fig.add_trace(
                go.Bar(
                    x=[tool],
                    y=[z_L1_mean],
                    error_y=dict(type='data', array=[z_L1_std]),
                    name=tool,
                    marker_color=tool_color,
                    showlegend=True,
                    legendgroup=tool,
                    hovertemplate=f'{tool}<br>Mean: %{{y:.2f}}<br>Std: {z_L1_std:.2f}<extra></extra>'
                ),
                row=1, col=1
            )
            
            # Z-score L2 bar
            fig.add_trace(
                go.Bar(
                    x=[tool],
                    y=[z_L2_mean],
                    error_y=dict(type='data', array=[z_L2_std]),
                    name=tool,
                    marker_color=tool_color,
                    showlegend=False,
                    legendgroup=tool,
                    hovertemplate=f'{tool}<br>Mean: %{{y:.2f}}<br>Std: {z_L2_std:.2f}<extra></extra>'
                ),
                row=1, col=2
            )
            
            # Z-score cos bar
            fig.add_trace(
                go.Bar(
                    x=[tool],
                    y=[z_cos_mean],
                    error_y=dict(type='data', array=[z_cos_std]),
                    name=tool,
                    marker_color=tool_color,
                    showlegend=False,
                    legendgroup=tool,
                    hovertemplate=f'{tool}<br>Mean: %{{y:.2f}}<br>Std: {z_cos_std:.2f}<extra></extra>'
                ),
                row=1, col=3
            )
        
        # Add reference lines at z = ±2 and ±3
        for col in [1, 2, 3]:
            fig.add_hline(y=2, line_dash="dash", line_color="red", opacity=0.5, row=1, col=col, 
                         annotation_text="95% CI" if col == 3 else None, annotation_position="right")
            fig.add_hline(y=-2, line_dash="dash", line_color="red", opacity=0.5, row=1, col=col)
            fig.add_hline(y=3, line_dash="dash", line_color="darkred", opacity=0.3, row=1, col=col,
                         annotation_text="99.7% CI" if col == 3 else None, annotation_position="right")
            fig.add_hline(y=-3, line_dash="dash", line_color="darkred", opacity=0.3, row=1, col=col)
            fig.add_hline(y=0, line_dash="solid", line_color="gray", opacity=0.5, row=1, col=col)
        
        fig.update_xaxes(title_text="Tool", row=1, col=1)
        fig.update_xaxes(title_text="Tool", row=1, col=2)
        fig.update_xaxes(title_text="Tool", row=1, col=3)
        
        fig.update_yaxes(title_text="Z-score (lower = better fit)", row=1, col=1)
        fig.update_yaxes(title_text="Z-score (lower = better fit)", row=1, col=2)
        fig.update_yaxes(title_text="Z-score (higher = better fit)", row=1, col=3)
        
        burden_label = f'{int(burden):,}'
        fig.update_layout(
            height=450,
            title_text=f'Bootstrap Quality Assessment - {burden_label} mutations<br>'
                      f'<sub>Mean Z-scores with error bars showing std dev across bootstrap replicates</sub>',
            showlegend=True,
            autosize=True
        )
        
        figures[burden] = fig
    
    return figures


def create_bootstrap_zscore_distribution(df_bootstrap, tools, colors):
    """
    Create histogram/violin plots showing distribution of z-scores for bootstrap quality assessment.
    
    Args:
        df_bootstrap: DataFrame with bootstrap quality assessment results
        tools: List of tool names
        colors: Dict mapping tool names to colors
        
    Returns:
        Plotly figure
    """
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('Z-score L1 (Lower is Better Fit)', 
                       'Z-score L2 (Lower is Better Fit)', 
                       'Z-score Cosine (Lower is Better Fit)'),
        horizontal_spacing=0.1
    )
    
    for tool in tools:
        tool_data = df_bootstrap[df_bootstrap['tool'] == tool]
        
        if len(tool_data) == 0:
            continue
        
        # Z-score L1
        fig.add_trace(
            go.Violin(
                y=tool_data['z_L1'],
                name=tool,
                legendgroup=tool,
                scalegroup=tool,
                line_color=colors.get(tool, '#999'),
                showlegend=True,
                box_visible=True,
                meanline_visible=True
            ),
            row=1, col=1
        )
        
        # Z-score L2
        fig.add_trace(
            go.Violin(
                y=tool_data['z_L2'],
                name=tool,
                legendgroup=tool,
                scalegroup=tool,
                line_color=colors.get(tool, '#999'),
                showlegend=False,
                box_visible=True,
                meanline_visible=True
            ),
            row=1, col=2
        )
        
        # Z-score cosine
        fig.add_trace(
            go.Violin(
                y=tool_data['z_cos'],
                name=tool,
                legendgroup=tool,
                scalegroup=tool,
                line_color=colors.get(tool, '#999'),
                showlegend=False,
                box_visible=True,
                meanline_visible=True
            ),
            row=1, col=3
        )
    
    # Add reference lines at z = ±2 (95% confidence) and z = ±3 (99.7% confidence)
    for col in [1, 2, 3]:
        fig.add_hline(y=2, line_dash="dash", line_color="red", opacity=0.3, row=1, col=col)
        fig.add_hline(y=-2, line_dash="dash", line_color="red", opacity=0.3, row=1, col=col)
        fig.add_hline(y=3, line_dash="dash", line_color="darkred", opacity=0.2, row=1, col=col)
        fig.add_hline(y=-3, line_dash="dash", line_color="darkred", opacity=0.2, row=1, col=col)
        fig.add_hline(y=0, line_dash="solid", line_color="gray", opacity=0.5, row=1, col=col)
    
    fig.update_yaxes(title_text="Z-score", row=1, col=1)
    fig.update_yaxes(title_text="Z-score", row=1, col=2)
    fig.update_yaxes(title_text="Z-score", row=1, col=3)
    
    fig.update_layout(
        height=500,
        showlegend=True,
        title_text='Bootstrap Quality Assessment: Z-score Distributions',
        hovermode='closest',
        autosize=True
    )
    
    return fig


def create_bootstrap_scatter(df_bootstrap, tools, colors):
    """
    Create scatter plot comparing original fit quality vs bootstrap mean.
    
    Args:
        df_bootstrap: DataFrame with bootstrap quality assessment results
        tools: List of tool names
        colors: Dict mapping tool names to colors
        
    Returns:
        Plotly figure
    """
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('L1 Distance', 'L2 Distance', 'Cosine Similarity'),
        horizontal_spacing=0.1
    )
    
    for tool in tools:
        tool_data = df_bootstrap[df_bootstrap['tool'] == tool]
        
        if len(tool_data) == 0:
            continue
        
        # L1: Original vs Bootstrap
        fig.add_trace(
            go.Scatter(
                x=tool_data['L1_bootstrap_mean'],
                y=tool_data['L1'],
                mode='markers',
                name=tool,
                legendgroup=tool,
                marker=dict(
                    size=6,
                    color=colors.get(tool, '#999'),
                    opacity=0.6,
                    line=dict(width=0.5, color='black')
                ),
                text=[f"Sample: {s}<br>Z-score: {z:.2f}" 
                      for s, z in zip(tool_data['sample'], tool_data['z_L1'])],
                hovertemplate='<b>%{fullData.name}</b><br>Bootstrap: %{x:.4f}<br>Original: %{y:.4f}<br>%{text}<extra></extra>',
                showlegend=True
            ),
            row=1, col=1
        )
        
        # L2: Original vs Bootstrap
        fig.add_trace(
            go.Scatter(
                x=tool_data['L2_bootstrap_mean'],
                y=tool_data['L2'],
                mode='markers',
                name=tool,
                legendgroup=tool,
                marker=dict(
                    size=6,
                    color=colors.get(tool, '#999'),
                    opacity=0.6,
                    line=dict(width=0.5, color='black')
                ),
                text=[f"Sample: {s}<br>Z-score: {z:.2f}" 
                      for s, z in zip(tool_data['sample'], tool_data['z_L2'])],
                hovertemplate='<b>%{fullData.name}</b><br>Bootstrap: %{x:.4f}<br>Original: %{y:.4f}<br>%{text}<extra></extra>',
                showlegend=False
            ),
            row=1, col=2
        )
        
        # Cosine: Original vs Bootstrap
        fig.add_trace(
            go.Scatter(
                x=tool_data['cos_bootstrap_mean'],
                y=tool_data['cos'],
                mode='markers',
                name=tool,
                legendgroup=tool,
                marker=dict(
                    size=6,
                    color=colors.get(tool, '#999'),
                    opacity=0.6,
                    line=dict(width=0.5, color='black')
                ),
                text=[f"Sample: {s}<br>Z-score: {z:.2f}" 
                      for s, z in zip(tool_data['sample'], tool_data['z_cos'])],
                hovertemplate='<b>%{fullData.name}</b><br>Bootstrap: %{x:.4f}<br>Original: %{y:.4f}<br>%{text}<extra></extra>',
                showlegend=False
            ),
            row=1, col=3
        )
    
    # Add diagonal reference line (y = x) to all subplots
    for col in [1, 2, 3]:
        fig.add_trace(
            go.Scatter(
                x=[0, 1],
                y=[0, 1],
                mode='lines',
                line=dict(color='gray', dash='dash'),
                showlegend=False,
                hoverinfo='skip'
            ),
            row=1, col=col
        )
    
    fig.update_xaxes(title_text="Bootstrap Mean", row=1, col=1)
    fig.update_xaxes(title_text="Bootstrap Mean", row=1, col=2)
    fig.update_xaxes(title_text="Bootstrap Mean", row=1, col=3)
    
    fig.update_yaxes(title_text="Original Fit", row=1, col=1)
    fig.update_yaxes(title_text="Original Fit", row=1, col=2)
    fig.update_yaxes(title_text="Original Fit", row=1, col=3)
    
    fig.update_layout(
        height=500,
        showlegend=True,
        title_text='Bootstrap Quality: Original Fit vs Bootstrap Mean (Points on diagonal = expected)',
        hovermode='closest',
        autosize=True
    )
    
    return fig


def create_enp_vs_f1_plot(df_per_sample, cancer_type, tools, colors):
    """Create scatter plot of ENP vs F1 score for unknown signature scenarios.
    
    Args:
        df_per_sample: DataFrame with per-sample metrics including F1 and ENP
        cancer_type: Cancer type to plot
        tools: List of tools to include
        colors: Dictionary mapping tools to colors
    
    Returns:
        plotly.graph_objects.Figure or None if no data
    """
    if df_per_sample is None:
        return None
    
    # Check if ENP column exists
    if 'ENP' not in df_per_sample.columns:
        print(f"  Warning: ENP column not found in per-sample data")
        return None
    
    # Filter for this cancer type
    cancer_samples = df_per_sample[df_per_sample['cancer_type'] == cancer_type].copy()
    
    if len(cancer_samples) == 0:
        return None
    
    # Remove rows with NaN ENP or F1
    cancer_samples = cancer_samples.dropna(subset=['ENP', 'F1'])
    
    if len(cancer_samples) == 0:
        print(f"    No valid ENP/F1 data for {cancer_type}")
        return None
    
    print(f"    Plotting {len(cancer_samples)} samples with valid ENP and F1 for {cancer_type}")
    
    # Get unique mutation burdens for grouping
    mutation_burdens = sorted(cancer_samples['num_muts_condition'].unique())
    min_muts = min(mutation_burdens)
    max_muts = max(mutation_burdens)
    
    # Helper function to darken color based on mutation burden
    def darken_color(hex_color, mutation_burden):
        """Darken a hex color based on mutation burden (higher = darker)"""
        # Convert hex to RGB
        hex_color = hex_color.lstrip('#')
        r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        
        # Calculate darkness factor (0.6 for min, 1.0 for max)
        if max_muts == min_muts:
            factor = 1.0
        else:
            factor = 0.6 + 0.4 * (mutation_burden - min_muts) / (max_muts - min_muts)
        
        # Apply darkening
        r = int(r * factor)
        g = int(g * factor)
        b = int(b * factor)
        
        return f'#{r:02x}{g:02x}{b:02x}'
    
    # Create scatter plot
    fig = go.Figure()
    
    method_symbols = get_method_symbols()
    
    # First pass: add traces for each tool (grouped by tool, colored by mutation burden)
    for tool in tools:
        tool_data = cancer_samples[cancer_samples['tool'] == tool]
        if len(tool_data) == 0:
            continue
        
        # Get base method to determine symbol
        base_method = get_tool_base_method(tool)
        symbol = method_symbols.get(base_method, 'circle')
        
        # Group by mutation burden for coloring and filtering
        for idx, num_muts in enumerate(mutation_burdens):
            muts_data = tool_data[tool_data['num_muts_condition'] == num_muts]
            if len(muts_data) == 0:
                continue
            
            # Add slight jitter to show overlapping points
            jitter_amount = 0.005  # Small jitter on F1 axis
            jittered_f1 = muts_data['F1'].values + np.random.uniform(-jitter_amount, jitter_amount, len(muts_data))
            
            # Get darkened color based on mutation burden
            point_color = darken_color(colors.get(tool, '#999999'), num_muts)
            
            fig.add_trace(go.Scatter(
                x=muts_data['ENP'],
                y=jittered_f1,
                mode='markers',
                name=tool,
                legendgroup=tool,
                legendgrouptitle_text='Tools' if tool == tools[0] else None,
                showlegend=(idx == 0),  # Only show legend entry once per tool
                marker=dict(
                    color=point_color,
                    size=10,
                    symbol=symbol,
                    opacity=0.7,
                    line=dict(width=1, color='white')
                ),
                customdata=[[num_muts, 'tool']] * len(muts_data),  # Mark as tool trace
                text=[f"Tool: {tool}<br>Sample: {row['sample']}<br>Muts: {num_muts}<br>Realization: {row.get('realization', 'N/A')}<br>highest_ENP: {row.get('highest_ENP', row.get('ENP', 'N/A')):.2f}<br>F1: {row['F1']:.3f}" 
                      for _, row in muts_data.iterrows()],
                hovertemplate='%{text}<extra></extra>'
            ))
    
    # Second pass: add invisible traces for mutation burden filtering
    for num_muts in mutation_burdens:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            name=f"{num_muts} muts",
            legendgroup='mutations',
            legendgrouptitle_text='Mutation Burden' if num_muts == mutation_burdens[0] else None,
            showlegend=True,
            marker=dict(
                size=10,
                color='rgba(128,128,128,0.5)',
                line=dict(width=1, color='gray')
            ),
            customdata=[[num_muts, 'filter']],  # Mark as filter trace
            visible=True
        ))
    
    fig.update_layout(
        title=dict(
            text=f'Effective Number of Parties vs F1 Score<br><sub>{cancer_type}</sub>',
            font=dict(size=16)
        ),
        xaxis=dict(
            title='ENP (Effective Number of Parties)',
            gridcolor='#eee',
            showline=True,
            linewidth=1,
            linecolor='#ddd'
        ),
        yaxis=dict(
            title='F1 Score',
            gridcolor='#eee',
            showline=True,
            linewidth=1,
            linecolor='#ddd',
            range=[0, 1.05]
        ),
        plot_bgcolor='white',
        hovermode='closest',
        showlegend=True,
        legend=dict(
            x=1.02,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='#ddd',
            borderwidth=1
        ),
        height=500,
        margin=dict(l=80, r=150, t=80, b=60)
    )
    
    return fig


def create_enp_vs_tae_plot(df_per_sample, cancer_type, tools, colors):
    """Create scatter plot of ENP vs TAE for unknown signature scenarios.
    
    Args:
        df_per_sample: DataFrame with per-sample metrics including TAE and ENP
        cancer_type: Cancer type to plot
        tools: List of tools to include
        colors: Dictionary mapping tools to colors
    
    Returns:
        plotly.graph_objects.Figure or None if no data
    """
    if df_per_sample is None:
        return None
    
    # Check if ENP and TAE columns exist
    if 'ENP' not in df_per_sample.columns or 'TAE' not in df_per_sample.columns:
        print(f"  Warning: ENP or TAE column not found in per-sample data")
        return None
    
    # Filter for this cancer type
    cancer_samples = df_per_sample[df_per_sample['cancer_type'] == cancer_type].copy()
    
    if len(cancer_samples) == 0:
        return None
    
    # Remove rows with NaN ENP or TAE
    cancer_samples = cancer_samples.dropna(subset=['ENP', 'TAE'])
    
    if len(cancer_samples) == 0:
        print(f"    No valid ENP/TAE data for {cancer_type}")
        return None
    
    print(f"    Plotting {len(cancer_samples)} samples with valid ENP and TAE for {cancer_type}")
    
    # Get unique mutation burdens for grouping
    mutation_burdens = sorted(cancer_samples['num_muts_condition'].unique())
    min_muts = min(mutation_burdens)
    max_muts = max(mutation_burdens)
    
    # Helper function to darken color based on mutation burden
    def darken_color(hex_color, mutation_burden):
        """Darken a hex color based on mutation burden (higher = darker)"""
        # Convert hex to RGB
        hex_color = hex_color.lstrip('#')
        r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        
        # Calculate darkness factor (0.6 for min, 1.0 for max)
        if max_muts == min_muts:
            factor = 1.0
        else:
            factor = 0.6 + 0.4 * (mutation_burden - min_muts) / (max_muts - min_muts)
        
        # Apply darkening
        r = int(r * factor)
        g = int(g * factor)
        b = int(b * factor)
        
        return f'#{r:02x}{g:02x}{b:02x}'
    
    # Create scatter plot
    fig = go.Figure()
    
    method_symbols = get_method_symbols()
    
    # First pass: add traces for each tool (grouped by tool, colored by mutation burden)
    for tool in tools:
        tool_data = cancer_samples[cancer_samples['tool'] == tool]
        if len(tool_data) == 0:
            continue
        
        # Get base method to determine symbol
        base_method = get_tool_base_method(tool)
        symbol = method_symbols.get(base_method, 'circle')
        
        # Group by mutation burden for coloring and filtering
        for idx, num_muts in enumerate(mutation_burdens):
            muts_data = tool_data[tool_data['num_muts_condition'] == num_muts]
            if len(muts_data) == 0:
                continue
            
            # Add slight jitter to show overlapping points
            jitter_amount = 0.002  # Small jitter on TAE axis
            jittered_tae = muts_data['TAE'].values + np.random.uniform(-jitter_amount, jitter_amount, len(muts_data))
            
            # Get darkened color based on mutation burden
            point_color = darken_color(colors.get(tool, '#999999'), num_muts)
            
            fig.add_trace(go.Scatter(
                x=muts_data['ENP'],
                y=jittered_tae,
                mode='markers',
                name=tool,
                legendgroup=tool,
                legendgrouptitle_text='Tools' if tool == tools[0] else None,
                showlegend=(idx == 0),  # Only show legend entry once per tool
                marker=dict(
                    color=point_color,
                    size=10,
                    symbol=symbol,
                    opacity=0.7,
                    line=dict(width=1, color='white')
                ),
                customdata=[[num_muts, 'tool']] * len(muts_data),  # Mark as tool trace
                text=[f"Tool: {tool}<br>Sample: {row['sample']}<br>Muts: {num_muts}<br>Realization: {row.get('realization', 'N/A')}<br>highest_ENP: {row.get('highest_ENP', row.get('ENP', 'N/A')):.2f}<br>TAE: {row['TAE']:.4f}" 
                      for _, row in muts_data.iterrows()],
                hovertemplate='%{text}<extra></extra>'
            ))
    
    # Second pass: add invisible traces for mutation burden filtering
    for num_muts in mutation_burdens:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            name=f"{num_muts} muts",
            legendgroup='mutations',
            legendgrouptitle_text='Mutation Burden' if num_muts == mutation_burdens[0] else None,
            showlegend=True,
            marker=dict(
                size=10,
                color='rgba(128,128,128,0.5)',
                line=dict(width=1, color='gray')
            ),
            customdata=[[num_muts, 'filter']],  # Mark as filter trace
            visible=True
        ))
    
    fig.update_layout(
        title=dict(
            text=f'Effective Number of Parties vs Total Absolute Error<br><sub>{cancer_type}</sub>',
            font=dict(size=16)
        ),
        xaxis=dict(
            title='ENP (Effective Number of Parties)',
            gridcolor='#eee',
            showline=True,
            linewidth=1,
            linecolor='#ddd'
        ),
        yaxis=dict(
            title='TAE (Total Absolute Error)',
            gridcolor='#eee',
            showline=True,
            linewidth=1,
            linecolor='#ddd'
        ),
        plot_bgcolor='white',
        hovermode='closest',
        showlegend=True,
        legend=dict(
            x=1.02,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='#ddd',
            borderwidth=1
        ),
        height=500,
        margin=dict(l=80, r=150, t=80, b=60)
    )
    
    return fig


def create_interactive_group_comparison_html(df_per_sample, cancer_type, tools, colors):
    """
    Create HTML/JavaScript for interactive group comparison tool.
    Allows filtering by: mutation burden, ENP of flattest known signature,
    ENP of flattest unknown signature.
    """
    import json

    if df_per_sample is None:
        return ""

    cancer_samples = df_per_sample[df_per_sample['cancer_type'] == cancer_type].copy()
    if len(cancer_samples) == 0:
        return ""

    available_metrics = [m for m in ['F1', 'TAE', 'MCC', 'nRMSE', 'num_FP', 'num_FN', 'P', 'R']
                         if m in cancer_samples.columns]
    if not available_metrics:
        return ""

    mutation_burdens = sorted(cancer_samples['num_muts'].unique())

    has_enp_known    = 'highest_ENP_known'    in cancer_samples.columns
    has_enp_unknown  = 'highest_ENP_unknown'  in cancer_samples.columns

    def _range(col):
        valid = cancer_samples[col].dropna()
        if len(valid) == 0:
            return (0.0, 0.0)
        return (float(valid.min()), float(valid.max()))

    enp_known_range   = _range('highest_ENP_known')   if has_enp_known   else (0.0, 0.0)
    enp_unknown_range = _range('highest_ENP_unknown') if has_enp_unknown else (0.0, 0.0)

    # Performance metric filters
    perf_filter_cols = [('F1', 'F1 Score'), ('TAE', 'TAE'), ('MCC', 'MCC'),
                        ('P', 'Precision'), ('R', 'Recall'), ('nRMSE', 'nRMSE')]
    perf_filters = [(col, label, _range(col))
                    for col, label in perf_filter_cols if col in cancer_samples.columns]

    plot_id = cancer_type.replace('-', '_') + '_interactive'

    columns_to_include = ['sample', 'tool', 'num_muts'] + available_metrics
    for col in ['highest_ENP_known', 'highest_ENP_unknown']:
        if col in cancer_samples.columns:
            columns_to_include.append(col)

    data_dict = {
        'samples': cancer_samples[columns_to_include].where(
            cancer_samples[columns_to_include].notna(), other=None
        ).to_dict('records'),
        'tools': tools,
        'colors': colors,
        'metrics': available_metrics,
        'mutation_burdens': [int(m) for m in mutation_burdens],
        'has_enp_known':   has_enp_known,
        'has_enp_unknown': has_enp_unknown,
        'enp_known_range':   list(enp_known_range),
        'enp_unknown_range': list(enp_unknown_range),
        'perf_filters': {col: {'label': label, 'range': list(rng)}
                         for col, label, rng in perf_filters},
    }
    data_json = json.dumps(data_dict)

    def build_filter_html(group_letter, group_color):
        html_f = f"""
        <div style="background: white; padding: 1rem; border-radius: 4px; border: 2px solid {group_color};">
            <h5 style="margin-top: 0; color: {group_color};">Group {group_letter}</h5>

            <div style="margin-bottom: 0.8rem;">
                <label style="display: block; margin-bottom: 0.3rem; font-weight: 500;">Mutation Burden:</label>
                <select id="{plot_id}_group{group_letter}_burden"
                        style="width: 100%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;">
                    {''.join(f'<option value="{int(m)}">{int(m):,} mutations</option>' for m in mutation_burdens)}
                </select>
            </div>
        """

        if has_enp_known:
            lo, hi = enp_known_range
            html_f += f"""
            <div style="margin-bottom: 0.8rem;">
                <label style="display: block; margin-bottom: 0.3rem; font-weight: 500;">
                    ENP Known — flattest known signature ({lo:.1f} – {hi:.1f}):
                </label>
                <div style="display: flex; gap: 0.5rem; align-items: center;">
                    <input type="number" id="{plot_id}_group{group_letter}_enp_known_min" placeholder="Min"
                           style="width: 45%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;"
                           step="0.1" min="{lo:.1f}" max="{hi:.1f}">
                    <span style="font-size: 0.9em;">to</span>
                    <input type="number" id="{plot_id}_group{group_letter}_enp_known_max" placeholder="Max"
                           style="width: 45%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;"
                           step="0.1" min="{lo:.1f}" max="{hi:.1f}">
                </div>
            </div>
            """

        if has_enp_unknown:
            lo, hi = enp_unknown_range
            html_f += f"""
            <div style="margin-bottom: 0.8rem;">
                <label style="display: block; margin-bottom: 0.3rem; font-weight: 500;">
                    ENP Unknown — flattest unknown signature ({lo:.1f} – {hi:.1f}):
                </label>
                <div style="display: flex; gap: 0.5rem; align-items: center;">
                    <input type="number" id="{plot_id}_group{group_letter}_enp_unknown_min" placeholder="Min"
                           style="width: 45%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;"
                           step="0.1" min="{lo:.1f}" max="{hi:.1f}">
                    <span style="font-size: 0.9em;">to</span>
                    <input type="number" id="{plot_id}_group{group_letter}_enp_unknown_max" placeholder="Max"
                           style="width: 45%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;"
                           step="0.1" min="{lo:.1f}" max="{hi:.1f}">
                </div>
            </div>
            """

        for col, label, (lo, hi) in perf_filters:
            step = "0.001" if hi <= 1.0 else "0.01"
            html_f += f"""
            <div style="margin-bottom: 0.8rem;">
                <label style="display: block; margin-bottom: 0.3rem; font-weight: 500;">
                    {label} ({lo:.3f} – {hi:.3f}):
                </label>
                <div style="display: flex; gap: 0.5rem; align-items: center;">
                    <input type="number" id="{plot_id}_group{group_letter}_perf_{col}_min" placeholder="Min"
                           style="width: 45%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;"
                           step="{step}" min="{lo:.3f}" max="{hi:.3f}">
                    <span style="font-size: 0.9em;">to</span>
                    <input type="number" id="{plot_id}_group{group_letter}_perf_{col}_max" placeholder="Max"
                           style="width: 45%; padding: 0.4rem; border: 1px solid #ddd; border-radius: 4px;"
                           step="{step}" min="{lo:.3f}" max="{hi:.3f}">
                </div>
            </div>
            """

        html_f += f"""
            <div style="font-size: 0.9em; color: #666; margin-top: 0.8rem; padding: 0.5rem; background: #f9f9f9; border-radius: 4px;">
                <span id="{plot_id}_group{group_letter}_count" style="font-weight: 600;">0 samples selected</span>
            </div>
        </div>
        """
        return html_f

    html = f"""
<div style="margin: 2rem 0; padding: 1.5rem; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
    <h4 style="margin-top: 0; color: white; text-shadow: 1px 1px 2px rgba(0,0,0,0.2);">
        Interactive Signature Complexity Analysis
    </h4>
    <p style="color: rgba(255,255,255,0.95); margin-bottom: 1.5rem; font-size: 0.95em;">
        Compare tool performance across two sample groups defined by mutation burden and signature complexity
        (ENP = Effective Number of Parties — higher means a flatter, harder-to-attribute spectrum).
    </p>

    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1.5rem; margin-bottom: 1.5rem;">
        {build_filter_html('A', '#4CAF50')}
        {build_filter_html('B', '#FF9800')}
    </div>

    <div style="background: white; padding: 1rem; border-radius: 4px; margin-bottom: 1rem; display: flex; align-items: center; gap: 1rem; flex-wrap: wrap;">
        <div style="flex: 1; min-width: 200px;">
            <label style="display: block; margin-bottom: 0.3rem; font-weight: 500; color: #333;">Metric to Compare:</label>
            <select id="{plot_id}_metric" style="width: 100%; padding: 0.5rem; border: 1px solid #ddd; border-radius: 4px; font-size: 1em;">
                {''.join(f'<option value="{m}">{m}</option>' for m in available_metrics)}
            </select>
        </div>
        <div style="display: flex; align-items: flex-end;">
            <button onclick="updateComparison_{plot_id}()"
                    style="padding: 0.7rem 2rem; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border: none; border-radius: 4px; cursor: pointer; font-weight: 600; font-size: 1em; box-shadow: 0 2px 4px rgba(0,0,0,0.2);">
                Update Comparison
            </button>
        </div>
    </div>

    <div style="background: white; padding: 1rem; border-radius: 4px;">
        <div id="{plot_id}_plot" style="width: 100%; height: 500px;"></div>
    </div>

    <div style="margin-top: 1rem; padding: 1rem; background: rgba(255,255,255,0.1); border-radius: 4px; color: rgba(255,255,255,0.9); font-size: 0.85em;">
        <strong>Tip:</strong> Compare samples where the flattest known signature has high vs low ENP to see
        whether diffuse reference signatures hurt fitting accuracy. Use ENP Unknown to assess how the
        difficulty of the hidden signal affects tool performance.
    </div>
</div>

<script>
(function() {{
    var plotData = {data_json};
    var plotId = '{plot_id}';

    function getGroupFilters(group) {{
        var filters = {{}};
        filters.burden = parseInt(document.getElementById(plotId + '_group' + group + '_burden').value);

        if (plotData.has_enp_known) {{
            var el_min = document.getElementById(plotId + '_group' + group + '_enp_known_min');
            var el_max = document.getElementById(plotId + '_group' + group + '_enp_known_max');
            filters.enp_known_min = el_min && el_min.value !== '' ? parseFloat(el_min.value) : null;
            filters.enp_known_max = el_max && el_max.value !== '' ? parseFloat(el_max.value) : null;
        }}

        if (plotData.has_enp_unknown) {{
            var el_min = document.getElementById(plotId + '_group' + group + '_enp_unknown_min');
            var el_max = document.getElementById(plotId + '_group' + group + '_enp_unknown_max');
            filters.enp_unknown_min = el_min && el_min.value !== '' ? parseFloat(el_min.value) : null;
            filters.enp_unknown_max = el_max && el_max.value !== '' ? parseFloat(el_max.value) : null;
        }}

        filters.perf = {{}};
        Object.keys(plotData.perf_filters).forEach(function(col) {{
            var el_min = document.getElementById(plotId + '_group' + group + '_perf_' + col + '_min');
            var el_max = document.getElementById(plotId + '_group' + group + '_perf_' + col + '_max');
            filters.perf[col] = {{
                min: el_min && el_min.value !== '' ? parseFloat(el_min.value) : null,
                max: el_max && el_max.value !== '' ? parseFloat(el_max.value) : null
            }};
        }});

        return filters;
    }}

    function filterSamples(filters) {{
        return plotData.samples.filter(function(s) {{
            if (s.num_muts !== filters.burden) return false;

            if (filters.enp_known_min !== null && !isNaN(filters.enp_known_min)) {{
                if (s.highest_ENP_known === null || s.highest_ENP_known < filters.enp_known_min) return false;
            }}
            if (filters.enp_known_max !== null && !isNaN(filters.enp_known_max)) {{
                if (s.highest_ENP_known === null || s.highest_ENP_known > filters.enp_known_max) return false;
            }}

            if (filters.enp_unknown_min !== null && !isNaN(filters.enp_unknown_min)) {{
                if (s.highest_ENP_unknown === null || s.highest_ENP_unknown < filters.enp_unknown_min) return false;
            }}
            if (filters.enp_unknown_max !== null && !isNaN(filters.enp_unknown_max)) {{
                if (s.highest_ENP_unknown === null || s.highest_ENP_unknown > filters.enp_unknown_max) return false;
            }}

            var perf = filters.perf || {{}};
            for (var col in perf) {{
                var pf = perf[col];
                var val = s[col];
                if (pf.min !== null && !isNaN(pf.min) && (val === null || val === undefined || val < pf.min)) return false;
                if (pf.max !== null && !isNaN(pf.max) && (val === null || val === undefined || val > pf.max)) return false;
            }}

            return true;
        }});
    }}

    function updateComparison_{plot_id}() {{
        var filtersA = getGroupFilters('A');
        var filtersB = getGroupFilters('B');
        var metric = document.getElementById(plotId + '_metric').value;

        var groupA_samples = filterSamples(filtersA);
        var groupB_samples = filterSamples(filtersB);

        var uniqueSamplesA = new Set(groupA_samples.map(function(s) {{ return s.sample; }})).size;
        var uniqueSamplesB = new Set(groupB_samples.map(function(s) {{ return s.sample; }})).size;
        document.getElementById(plotId + '_groupA_count').textContent =
            uniqueSamplesA + ' samples (' + groupA_samples.length + ' data points)';
        document.getElementById(plotId + '_groupB_count').textContent =
            uniqueSamplesB + ' samples (' + groupB_samples.length + ' data points)';

        var traces = [];
        plotData.tools.forEach(function(tool) {{
            var dataA = groupA_samples.filter(function(s) {{ return s.tool === tool; }})
                .map(function(s) {{ return s[metric]; }})
                .filter(function(v) {{ return v !== null && v !== undefined && !isNaN(v); }});
            if (dataA.length > 0) {{
                traces.push({{
                    y: dataA, name: tool, type: 'box',
                    marker: {{ color: plotData.colors[tool] || '#888' }},
                    boxmean: 'sd', offsetgroup: tool, legendgroup: tool,
                    showlegend: true, xaxis: 'x1'
                }});
            }}

            var dataB = groupB_samples.filter(function(s) {{ return s.tool === tool; }})
                .map(function(s) {{ return s[metric]; }})
                .filter(function(v) {{ return v !== null && v !== undefined && !isNaN(v); }});
            if (dataB.length > 0) {{
                traces.push({{
                    y: dataB, name: tool, type: 'box',
                    marker: {{ color: plotData.colors[tool] || '#888' }},
                    boxmean: 'sd', offsetgroup: tool, legendgroup: tool,
                    showlegend: false, xaxis: 'x2'
                }});
            }}
        }});

        var layout = {{
            title: {{ text: 'Group Comparison: ' + metric, font: {{ size: 16 }} }},
            xaxis:  {{ domain: [0, 0.45], title: {{ text: 'Group A', font: {{ color: '#4CAF50' }} }}, showticklabels: false }},
            xaxis2: {{ domain: [0.55, 1], title: {{ text: 'Group B', font: {{ color: '#FF9800' }} }}, showticklabels: false }},
            yaxis: {{ title: {{ text: metric }}, zeroline: true }},
            boxmode: 'group', showlegend: true,
            legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: -0.15 }},
            hovermode: 'closest',
            margin: {{ t: 60, b: 80, l: 60, r: 40 }}
        }};

        Plotly.newPlot(plotId + '_plot', traces, layout, {{responsive: true}});
    }}

    var allInputs = [
        plotId + '_groupA_burden', plotId + '_groupB_burden',
        plotId + '_groupA_enp_known_min',    plotId + '_groupA_enp_known_max',
        plotId + '_groupB_enp_known_min',    plotId + '_groupB_enp_known_max',
        plotId + '_groupA_enp_unknown_min',  plotId + '_groupA_enp_unknown_max',
        plotId + '_groupB_enp_unknown_min',  plotId + '_groupB_enp_unknown_max',
        plotId + '_metric'
    ];
    Object.keys(plotData.perf_filters).forEach(function(col) {{
        ['A', 'B'].forEach(function(g) {{
            allInputs.push(plotId + '_group' + g + '_perf_' + col + '_min');
            allInputs.push(plotId + '_group' + g + '_perf_' + col + '_max');
        }});
    }});
    allInputs.forEach(function(id) {{
        var elem = document.getElementById(id);
        if (elem) elem.addEventListener('change', updateComparison_{plot_id});
    }});

    setTimeout(updateComparison_{plot_id}, 100);
}})();
</script>
"""
    
    return html
def create_visualizations_section(df, df_per_sample=None, df_bootstrap=None, df_raw=None,
                                   results_dir=None, id_prefix='', include_plotly_cdn=True):
    """Create main visualizations section with tabs.

    Args:
        df: DataFrame with aggregated metrics
        df_per_sample: Optional DataFrame with per-sample metrics for violin plots
        df_bootstrap: Optional DataFrame with bootstrap quality metrics
        df_raw: Optional DataFrame with individual realization data (before aggregation)
        results_dir: Optional base results directory for reading contribution files
        id_prefix: String prefix for all HTML element IDs (needed when called multiple times)
        include_plotly_cdn: Whether to include plotly.js CDN link in the first figure
    """
    # Get tool order and filter to only tools present in data
    tool_order = get_tool_display_order()
    available_tools = set(df['tool'].unique())
    tools = [t for t in tool_order if t in available_tools]
    # Add any tools not in our predefined order (sorted alphabetically at the end)
    remaining_tools = sorted(available_tools - set(tools))
    tools.extend(remaining_tools)

    cancer_types = sorted(df['cancer_type'].unique())
    colors = get_tool_colors()

    html = '<section id="visualizations" class="section">\n'
    html += '<h2>📈 Visualizations</h2>\n'

    # Create tabs
    html += '<div class="tabs">\n'
    for i, cancer_type in enumerate(cancer_types):
        tab_id = f'{id_prefix}{cancer_type.replace("-", "_")}'
        active = 'active' if i == 0 else ''
        html += f'<button class="tab {active}" onclick="switchTab(\'{tab_id}-tab\', event)">{cancer_type}</button>\n'
    html += f'<button class="tab" onclick="switchTab(\'{id_prefix}aggregate-tab\', event)">Aggregate Analysis</button>\n'
    html += '</div>\n'

    # Per-cancer-type tabs
    for i, cancer_type in enumerate(cancer_types):
        tab_id = f'{id_prefix}{cancer_type.replace("-", "_")}'
        active = 'active' if i == 0 else ''
        cancer_df = df[df['cancer_type'] == cancer_type]

        html += f'<div id="{tab_id}-tab" class="tab-content {active}">\n'
        html += f'<h3>{cancer_type}</h3>\n'
        
        # SECTION 1: Tool Comparison
        html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #2c5aa0; border-radius: 4px;">\n'
        html += '<h4 style="margin-top: 0; color: #2c5aa0;">🔍 Tool Comparison — which tools identify and quantify signatures most accurately?</h4>\n'
        html += '<p style="color: #555; font-size: 0.9em; margin: 0 0 1.2rem 0;">Each point is the mean across all simulated samples at a given mutation burden. Higher F1 / lower TAE = better. Expand <em>Metric Reference</em> in the Summary section for definitions.</p>\n'

        # Precision-Recall Trade-off
        fig_pr = create_precision_recall_subplots(cancer_df, cancer_type, tools, colors, df_per_sample)
        html += '<div class="plot-container">\n'
        html += fig_pr.to_html(include_plotlyjs=True if (i == 0 and include_plotly_cdn) else False, div_id=f'pr_{tab_id}', config={'responsive': True})
        html += '</div>\n'
        html += '<p style="color: #666; margin: 1rem 0;"><strong>Precision–Recall Trade-off:</strong> <em>Precision</em> = how many reported signatures were actually present; <em>Recall</em> = how many true signatures were found. Error bars = ±1 SEM. Diagonal lines = F1 contours. A good tool sits in the top-right corner.</p>\n'

        # Performance Metrics
        fig_metrics = create_metrics_grid(cancer_df, cancer_type, tools, colors, df_per_sample)
        html += '<div class="plot-container" style="margin-top: 2rem;">\n'
        html += fig_metrics.to_html(include_plotlyjs=False, div_id=f'metrics_{tab_id}', config={'responsive': True})
        html += '</div>\n'
        html += '<p style="color: #666; margin: 1rem 0;"><strong>Performance vs Mutation Burden:</strong> Mean ± SEM across all samples (x-axis log scale). F1 / P / R / MCC: higher is better. TAE / nRMSE / FP / FN: lower is better.</p>\n'
        html += '</div>\n'  # Close tool comparison section
        
        # SECTION 2: Algorithm Stability
        if df_per_sample is not None:
            print(f"  Checking stability plot for {cancer_type}...")
            print(f"    df_per_sample shape: {df_per_sample.shape}")
            print(f"    Has realization column: {'realization' in df_per_sample.columns}")
            if 'realization' in df_per_sample.columns:
                print(f"    Realizations: {df_per_sample['realization'].unique()}")
            fig_stability = create_realization_stability_plot(df_per_sample, cancer_type, tools, colors)
            if fig_stability is not None:
                html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #d97706; border-radius: 4px;">\n'
                html += '<h4 style="margin-top: 0; color: #d97706;">Algorithm Stability — does the tool give consistent results across independent repeats?</h4>\n'
                html += '<p style="color: #666; margin-bottom: 1rem;">Each realization uses the same signature mix but a different random draw of mutations. The variance shown here measures how much a tool\'s result changes across these repeats — a stable tool gives the same answer regardless of the random seed. Lower variance = more reproducible.</p>\n'
                
                # Render the figure HTML fragment (no plotly.js) so we can embed custom controls
                fig_stability_html = fig_stability.to_html(include_plotlyjs=False, full_html=False, div_id=f'stability_{tab_id}', config={'responsive': True})

                # Build controls: tool checkboxes and mutation-burden checkboxes
                burdens_stability = sorted(df_per_sample[df_per_sample['cancer_type'] == cancer_type]['num_muts'].unique())

                controls_stability = []
                controls_stability.append('<div class="controls-header"><h4>Filters</h4></div>')
                # Tools
                controls_stability.append('<div class="controls-section"><strong>Tools</strong><br/>')
                for tool in tools:
                    controls_stability.append(f'<label><input type="checkbox" class="tool-checkbox-stability-{tab_id}" value="{tool}" checked> {tool}</label><br/>')
                controls_stability.append('</div>')

                # Mutation burdens
                controls_stability.append('<div class="controls-section"><strong>Mutation burden</strong><br/>')
                for mb in burdens_stability:
                    controls_stability.append(f'<label><input type="checkbox" class="burden-checkbox-stability-{tab_id}" value="{int(mb)}" checked> {int(mb):,} muts</label><br/>')
                controls_stability.append('</div>')

                controls_stability_html = '\n'.join(controls_stability)

                # JavaScript to hook checkboxes and toggle trace visibility based on trace.meta
                script_stability = f"""
<script>
(function(){{
  // Wait until plot div exists
  var gd = document.getElementById('stability_{tab_id}');
  if(!gd) return;

  function getSelections(){{
    var tools = Array.from(document.querySelectorAll('.tool-checkbox-stability-{tab_id}:checked')).map(function(n){{return n.value;}});
    var burdens = Array.from(document.querySelectorAll('.burden-checkbox-stability-{tab_id}:checked')).map(function(n){{return n.value;}});
    return {{tools: tools, burdens: burdens}};
  }}

  function updateVisibility(){{
    var sel = getSelections();
    var vis = gd.data.map(function(trace){{
      var meta = trace.meta || {{}};
      var t = meta.tool || trace.name || '';
      var b = String(meta.mut_burden || meta.burden || '');
      var toolMatch = (sel.tools.length === 0) || (sel.tools.indexOf(t) !== -1);
      var burdenMatch = (sel.burdens.length === 0) || (sel.burdens.indexOf(b) !== -1);
      return (toolMatch && burdenMatch);
    }});
    // Apply visibility array
    Plotly.restyle(gd, {{'visible': vis}});
  }}

  // Attach listeners
  Array.from(document.querySelectorAll('.tool-checkbox-stability-{tab_id}')).forEach(function(cb){{cb.addEventListener('change', updateVisibility);}});
  Array.from(document.querySelectorAll('.burden-checkbox-stability-{tab_id}')).forEach(function(cb){{cb.addEventListener('change', updateVisibility);}});

  // Initial apply
  setTimeout(updateVisibility, 100);
}})();
</script>
"""

                # Put controls to the right of the plot using a flex layout
                html += '<div class="plot-container">\n'
                html += '<div class="plot-flex">\n'
                html += '<div class="plot-box">\n'
                html += fig_stability_html
                html += '</div>\n'
                html += f'<div class="plot-controls" id="controls_stability_{tab_id}">\n{controls_stability_html}\n</div>\n'
                html += '</div>\n'
                html += script_stability
                html += '</div>\n'
                html += '</div>\n'  # Close stability section
            else:
                print(f"    Stability plot returned None for {cancer_type}")

        # SECTION 3: Performance per Method (formerly Per-Sample Variation)
        if df_per_sample is not None:
            fig_per_sample = create_per_sample_violin_plots(df_per_sample, cancer_type, tools, colors)
            if fig_per_sample is not None:
                html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #059669; border-radius: 4px;">\n'
                html += '<h4 style="margin-top: 0; color: #059669;">Performance Distribution per Tool</h4>\n'
                html += '<p style="color: #666; margin-bottom: 1rem;">Each violin/strip shows the spread of a metric across all simulated samples for one tool. A narrow, high-placed distribution = reliable and accurate. Use the dropdown to switch mutation burden.</p>\n'
                
                # Render the figure HTML with built-in dropdown (no custom controls needed)
                fig_html = fig_per_sample.to_html(include_plotlyjs=False, full_html=False, div_id=f'per_sample_{tab_id}', config={'responsive': True})

                # Simple layout without custom controls (dropdown is built into the plot)
                html += '<div class="plot-container">\n'
                html += fig_html
                html += '</div>\n'
                html += '</div>\n\n'
        
        # SECTION 3.5: Sample-by-Sample Performance Comparison
        if df_per_sample is not None:
            fig_per_method = create_per_method_sample_plots(df_per_sample, cancer_type, tools, colors)
            if fig_per_method is not None:
                html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #059669; border-radius: 4px;">\n'
                html += '<h4 style="margin-top: 0; color: #059669;">Sample-by-Sample Comparison — which samples are hard for all tools?</h4>\n'
                html += '<p style="color: #666; margin-bottom: 1rem;">Each column is one simulated sample. Grey boxplot = spread of all tools on that sample; colored markers = individual tool results. Samples where the boxplot is low or tools spread far apart are the hardest cases. Use the dropdown to switch mutation burden.</p>\n'
                
                # Render the figure HTML with built-in dropdown (no custom controls needed)
                fig_html = fig_per_method.to_html(include_plotlyjs=False, full_html=False, div_id=f'per_method_{tab_id}', config={'responsive': True})

                # Simple layout without custom controls (dropdown is built into the plot)
                html += '<div class="plot-container">\n'
                html += fig_html
                html += '</div>\n'
                html += '</div>\n\n'
        
        # SECTION 4: Interactive Feature vs Metric Analysis
        if df_per_sample is not None:
            html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #7c3aed; border-radius: 4px;">\n'
            html += '<h4 style="margin-top: 0; color: #7c3aed;">📊 Does Signature Complexity Affect Tool Performance?</h4>\n'
            html += '<p style="color: #666; margin-bottom: 1rem;">Some samples are harder than others because their active signatures have a more spread-out (diffuse) mutational spectrum — measured here by ENP (Effective Number of Parties, range 1–96). ENP = 1 means all mutations fall in a single channel; ENP = 96 means they are perfectly uniform across all channels. Filter samples by ENP range or other metrics, then compare how each tool performs on the selected subset.</p>\n'
            
            # Add explanations of the features
            html += '<div style="background: white; padding: 1rem; border-radius: 6px; margin-bottom: 1.5rem; border: 1px solid #e5e7eb;">\n'
            html += '<h5 style="margin-top: 0; color: #374151; font-size: 14px; font-weight: 600;">📐 Sample Features Explained:</h5>\n'
            html += '<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 0.75rem; font-size: 13px; color: #4b5563;">\n'
            html += '<div><strong>highest_ENP_known</strong>: Complexity of the most diffuse <em>known</em> (in-reference) signature active in this sample. Scale 1–96; higher = more uniform spectrum = harder to attribute correctly.</div>\n'
            html += '<div><strong>highest_ENP_unknown</strong>: Complexity of the most diffuse <em>unknown</em> (out-of-reference) signature active in this sample. Higher = the hidden signal is more spread out and harder to distinguish from background noise.</div>\n'
            html += '</div>\n'
            html += '</div>\n'
            
            # Add feature-metric scatter widget
            from html_templates import create_feature_metric_scatter_widget
            
            # Features are now included directly in per_sample_metrics file
            if results_dir:
                per_sample_tsv = Path(results_dir) / 'all_per_sample_metrics.tsv'
            else:
                per_sample_tsv = Path('all_per_sample_metrics.tsv')
            
            # Check if the file has feature columns (check for new names)
            has_features = False
            if per_sample_tsv.exists():
                try:
                    sample_df = pd.read_csv(per_sample_tsv, sep='\t', nrows=1)
                    feature_cols = ['highest_ENP_known', 'highest_ENP_unknown']
                    has_features = any(col in sample_df.columns for col in feature_cols)
                except:
                    pass
            
            if has_features:
                html += '<div style="margin-top: 1rem;">\n'
                # Pass the dataframe directly to avoid CORS issues with file:// protocol
                widget_html = create_feature_metric_scatter_widget(
                    per_sample_tsv.name,  # Kept for compatibility
                    per_sample_tsv.name,  # Kept for compatibility
                    widget_id=f'feature_metric_{tab_id}',
                    data_df=df_per_sample  # Pass actual data to embed in HTML
                )
                html += widget_html
                html += '</div>\n'
            else:
                html += '<div style="padding: 1rem; background: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; margin: 1rem 0;">\n'
                html += '<p style="margin: 0; color: #856404;"><strong>Note:</strong> Sample features not found. Re-run the workflow to compute features (ENP, condition number, KL divergence, etc.).</p>\n'
                html += '</div>\n'
            
            html += '</div>\n'  # Close complexity analysis section
        
        # Bootstrap quality assessment for this cancer type (DISABLED)
        # if df_bootstrap is not None and len(df_bootstrap) > 0:
        #     html += '<h4 style="margin-top: 2rem;">Bootstrap Quality Assessment</h4>\n'
        #     (bootstrap code commented out)
        
        html += '</div>\n'
    
    # Aggregate tab
    html += f'<div id="{id_prefix}aggregate-tab" class="tab-content">\n'
    html += '<h3>Aggregate Analysis</h3>\n'
    html += '<p>Performance averaged across all cancer types.</p>\n'
    
    # Precision-Recall plot
    fig_agg = create_aggregate_plot(df, tools, colors)
    html += '<div class="plot-container">\n'
    html += fig_agg.to_html(include_plotlyjs=False, div_id='aggregate_plot', config={'responsive': True})
    html += '</div>\n'
    
    # Aggregate metrics grid (pooled across cancer types)
    # Group by tool and num_muts, calculate mean and std across cancer types/realizations
    agg_metrics_df = df.groupby(['tool', 'num_muts']).agg({
        'F1_mean': 'mean',
        'F1_std': 'mean',  # Average the within-sample std
        'TAE': 'mean',
        'TAE_std': 'mean',
        'MCC_mean': 'mean',
        'MCC_std': 'mean',
        'nRMSE': ['mean', 'std'],  # Calculate std across realizations
        'num_FP_sigs': ['mean', 'std'],  # Calculate std across realizations
        'num_FN_sigs': ['mean', 'std'],  # Calculate std across realizations
        'num_samples': 'first'
    }).reset_index()
    
    # Flatten multi-level columns and rename
    agg_metrics_df.columns = ['tool', 'num_muts', 'F1_mean', 'F1_std', 'TAE', 'TAE_std', 
                               'MCC_mean', 'MCC_std', 'nRMSE', 'nRMSE_std', 
                               'num_FP_sigs', 'num_FP_sigs_std', 
                               'num_FN_sigs', 'num_FN_sigs_std', 'num_samples']
    
    fig_agg_metrics = create_metrics_grid(agg_metrics_df, 'All Cancer Types (Pooled)', tools, colors, df_per_sample)
    html += '<div class="plot-container">\n'
    html += fig_agg_metrics.to_html(include_plotlyjs=False, div_id='aggregate_metrics_grid', config={'responsive': True})
    html += '</div>\n'
    
    html += '</div>\n'
    
    # Bootstrap Quality Assessment tab (if data available)
    if df_bootstrap is not None and len(df_bootstrap) > 0:
        html += '<div id="bootstrap-tab" class="tab-content">\n'
        html += '<h3>Bootstrap Quality Assessment</h3>\n'
        html += '''
        <div class="info-box" style="background: #fff3cd; border-left-color: #ffc107; margin-bottom: 1.5rem;">
            <h4>📊 Understanding Bootstrap Quality Assessment</h4>
            <p style="margin-bottom: 0.5rem;">
                This analysis validates fitting quality by comparing the original fit to bootstrap samples 
                generated from the fitted signature activities. If the fitted activities are correct, the 
                original fit quality should be similar to the bootstrap distribution.
            </p>
            
            <h5 style="margin-top: 1rem; margin-bottom: 0.5rem;">Reconstruction Quality Metrics:</h5>
            <ul style="margin-left: 1.5rem; margin-bottom: 0.5rem;">
                <li><strong>L1 Distance (Manhattan):</strong> Sum of absolute differences between observed and reconstructed mutation profiles. 
                    Measures total deviation across all 96 mutation types. Range: [0, 2], where 0 = perfect match, 2 = completely different. 
                    <em>Lower is better.</em></li>
                <li><strong>L2 Distance (Euclidean):</strong> Square root of sum of squared differences. 
                    Penalizes large deviations more heavily than L1. Range: [0, √2], where 0 = perfect match. 
                    <em>Lower is better.</em></li>
                <li><strong>Cosine Similarity:</strong> Measures angular similarity between mutation profiles (ignores magnitude). 
                    Range: [0, 1], where 1 = perfect match, 0 = completely different directions. 
                    <em>Higher is better.</em></li>
            </ul>
            
            <h5 style="margin-top: 1rem; margin-bottom: 0.5rem;">Z-score Interpretation:</h5>
            <ul style="margin-left: 1.5rem; margin-bottom: 0.5rem;">
                <li><strong>Z-scores near 0:</strong> Fit quality matches expectations ✓ (good)</li>
                <li><strong>|Z| > 2:</strong> Unusual fit (outside 95% confidence)</li>
                <li><strong>|Z| > 3:</strong> Highly unusual (outside 99.7% confidence)</li>
                <li><strong>Positive Z (L1/L2):</strong> Original fit is worse than expected (higher distance = worse fit) → potential underfitting/convergence issues</li>
                <li><strong>Negative Z (L1/L2):</strong> Original fit is better than expected (lower distance) → potential overfitting to noise</li>
                <li><strong>Positive Z (cosine):</strong> Original fit is worse than expected (lower similarity) → poor reconstruction</li>
                <li><strong>Negative Z (cosine):</strong> Original fit is better than expected (higher similarity) → unusually good fit</li>
            </ul>
        </div>
        '''
        
        # Bootstrap plots (DISABLED)
        # fig_zscore = create_bootstrap_zscore_distribution(df_bootstrap, tools, colors)
        # fig_scatter = create_bootstrap_scatter(df_bootstrap, tools, colors)
        # (bootstrap code commented out)
        
        # Bootstrap summary statistics table (DISABLED)
        # html += '<h4>Summary Statistics</h4>\n'
        # (table code commented out)
    
    html += '</section>\n'
    
    return html


# ============================================================================
# Main Report Generator
# ============================================================================

def generate_report(df, output_file, title, subtitle, df_per_sample=None, df_bootstrap=None, df_raw=None, results_dir=None, standard_penalties=None, unknown_sig_config=None, cosmic_version_label=None):
    """
    Generate complete HTML report.

    Args:
        df: pandas.DataFrame with aggregated evaluation metrics
        output_file: Path to output HTML file
        title: Report title
        subtitle: Report subtitle
        df_per_sample: Optional pandas.DataFrame with per-sample metrics
        df_bootstrap: Optional pandas.DataFrame with bootstrap quality assessment results
        df_raw: Optional pandas.DataFrame with raw per-realization data for stability plots
        results_dir: Optional base results directory for reading contribution files
        standard_penalties: Optional dict of standard penalty values per tool
        unknown_sig_config: Optional dict with unknown signature configuration {'num_unknowns': int, 'weights': list}
        cosmic_version_label: Optional human-readable COSMIC catalog label derived from
                          config's cosmic_signatures.base_file
    """
    print("\nGenerating report sections...")

    # Generate all sections
    summary_html = create_summary_section(df, df_raw, unknown_sig_config)
    experiment_html = create_experiment_section(df, df_raw, unknown_sig_config, cosmic_version_label)

    # Build visualizations section — one sub-section per data_type when multiple exist
    data_types = sorted(df['data_type'].unique()) if 'data_type' in df.columns else ['synthetic']
    if len(data_types) <= 1:
        visualizations_html = create_visualizations_section(
            df, df_per_sample, df_bootstrap, df_raw, results_dir
        )
    else:
        # Toggle buttons
        vis_html = '<section id="visualizations" class="section">\n'
        vis_html += '<h2>\U0001f4c8 Visualizations</h2>\n'
        vis_html += '<div style="margin-bottom:1rem;">\n'
        vis_html += '<span style="font-weight:600;margin-right:0.5rem;">Data type:</span>\n'
        for j, dt in enumerate(data_types):
            active_style = 'background:#2c5aa0;color:#fff;' if j == 0 else 'background:#e5e7eb;color:#374151;'
            vis_html += (
                f'<button id="dt-btn-{dt}" onclick="switchDT(\'{dt}\')" '
                f'style="margin-right:0.4rem;padding:0.3rem 0.9rem;border:none;border-radius:4px;'
                f'cursor:pointer;font-size:14px;{active_style}">{dt}</button>\n'
            )
        vis_html += '</div>\n'
        vis_html += '<script>\n'
        vis_html += 'function switchDT(dt) {\n'
        vis_html += '  document.querySelectorAll(".dt-section").forEach(s => s.style.display = "none");\n'
        vis_html += '  document.querySelectorAll("[id^=\'dt-btn-\']").forEach(b => { b.style.background="#e5e7eb"; b.style.color="#374151"; });\n'
        vis_html += '  var sec = document.getElementById("dt-section-" + dt);\n'
        vis_html += '  if (sec) sec.style.display = "block";\n'
        vis_html += '  var btn = document.getElementById("dt-btn-" + dt);\n'
        vis_html += '  if (btn) { btn.style.background="#2c5aa0"; btn.style.color="#fff"; }\n'
        vis_html += '}\n'
        vis_html += '</script>\n'

        for dt_idx, dt in enumerate(data_types):
            dt_df = df[df['data_type'] == dt]
            dt_per_sample = None
            if df_per_sample is not None:
                dt_per_sample = (
                    df_per_sample[df_per_sample['data_type'] == dt]
                    if 'data_type' in df_per_sample.columns
                    else df_per_sample
                )
            display = 'block' if dt_idx == 0 else 'none'
            vis_html += f'<div id="dt-section-{dt}" class="dt-section" style="display:{display};">\n'
            vis_html += f'<h3 style="color:#6b7280;font-size:0.9rem;text-transform:uppercase;letter-spacing:0.05em;margin:0 0 1rem 0;">{dt}</h3>\n'
            # Strip the outer <section> wrapper that create_visualizations_section adds
            inner = create_visualizations_section(
                dt_df, dt_per_sample, df_bootstrap, df_raw, results_dir,
                id_prefix=f'{dt}_',
                include_plotly_cdn=(dt_idx == 0),
            )
            # Remove opening <section ...> and closing </section> tags so they don't nest
            inner = inner.replace('<section id="visualizations" class="section">\n', '', 1)
            inner = inner.replace('<h2>\U0001f4c8 Visualizations</h2>\n', '', 1)
            if inner.endswith('</section>\n'):
                inner = inner[:-len('</section>\n')]
            vis_html += inner
            vis_html += '</div>\n'  # close dt-section

        vis_html += '</section>\n'
        visualizations_html = vis_html

    content = summary_html + experiment_html + visualizations_html
    
    # Fill template
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    html = get_html_template().format(
        title=title,
        subtitle=subtitle,
        content=content,
        timestamp=timestamp
    )
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"\n✅ Report generated: {output_file}")
    print(f"   Open in browser: file://{Path(output_file).absolute()}")
    
    if df_bootstrap is not None and len(df_bootstrap) > 0:
        print(f"   ✓ Bootstrap quality assessment included ({len(df_bootstrap)} samples)")



# ============================================================================
# Command Line Interface
# ============================================================================

def main(summaries, per_sample, bootstrap_quality, output, title, subtitle,
         results_dir, config_file):
    """Main entry point."""
    # Load standard penalties and unknown signature info from config if provided
    standard_penalties = {}
    unknown_sig_config = None
    cosmic_version_label = None

    if config_file:
        print(f"Attempting to load config file: {config_file}")
        try:
            if not os.path.exists(config_file):
                print(f"ERROR: Config file does not exist: {config_file}")
            else:
                print(f"Config file exists, loading...")
            with open(config_file, 'r') as f:
                config_data = yaml.safe_load(f)
                standard_penalties = config_data.get('standard_penalties', {})
                cosmic_version_label = _derive_cosmic_label(
                    config_data.get('cosmic_signatures', {}).get('base_file')
                )
                out_of_ref_weights = config_data.get('out_of_reference_weights', [])
                if out_of_ref_weights:
                    unknown_sig_config = {
                        'num_unknowns': len(out_of_ref_weights),
                        'weights': out_of_ref_weights
                    }
                    print(f"  - With fixed weights: {out_of_ref_weights}")
                else:
                    print("\nNo unknown signatures configured (out_of_reference_sigs is empty)")
        except Exception as e:
            print(f"\n❌ ERROR: Could not load config file {config_file}: {e}")
    else:
        print("No config file provided")

    # Infer results_dir from output path if not provided
    if not results_dir and output:
        results_dir = str(Path(output).parent)

    # Load summary data only for df_raw (used in realization stability plots)
    _, df_raw = load_evaluation_data(summaries, return_raw=True)

    # Load per-sample metrics and apply per-sample oracle penalty selection:
    # for each (cancer_type, data_type, tool, num_muts, replicate, sample) keep
    # the row with the lowest individual-sample TAE across all penalty values.
    # Non-penalty tools have a single penalty value and are kept as-is.
    df_per_sample = None
    if per_sample:
        print(f"\nLoading per-sample metrics from {len(per_sample)} file(s)...")
        df_per_sample = load_per_sample_data(per_sample)
        if df_per_sample is not None and 'penalty' in df_per_sample.columns and 'TAE' in df_per_sample.columns:
            pen_tools_ps = set(
                df_per_sample.groupby('tool')['penalty']
                .nunique(dropna=True)
                .pipe(lambda s: s[s > 1].index)
            )
            if pen_tools_ps:
                pen_ps = df_per_sample[df_per_sample['tool'].isin(pen_tools_ps)]
                # Find single best penalty per tool: lowest mean TAE across all samples
                best_penalty_per_tool = (
                    pen_ps.groupby(['tool', 'penalty'])['TAE']
                    .mean()
                    .groupby('tool')
                    .idxmin()
                    .apply(lambda x: x[1])
                )
                print(f"  Best penalty per tool: {best_penalty_per_tool.to_dict()}")
                mask = pen_ps.apply(lambda r: r['penalty'] == best_penalty_per_tool[r['tool']], axis=1)
                pen_ps = pen_ps[mask]
                no_pen_ps = df_per_sample[~df_per_sample['tool'].isin(pen_tools_ps)]
                df_per_sample = pd.concat([pen_ps, no_pen_ps], ignore_index=True)
                print(f"  Best-overall penalty selection for {sorted(pen_tools_ps)}. Shape: {df_per_sample.shape}")

    # Re-aggregate the best-overall-penalty-selected per-sample data into the summary df.
    if df_per_sample is not None and not df_per_sample.empty:
        df = reaggregate_from_per_sample(df_per_sample)
        print(f"  Summary df re-aggregated from best-overall-penalty data. Shape: {df.shape}")
    else:
        # Fallback: use pre-aggregated summary with per-replicate lowest-TAE selection
        df, _ = load_evaluation_data(summaries, return_raw=True)
        if 'penalty' in df.columns and 'TAE' in df.columns:
            penalty_tools = set(
                df.groupby('tool')['penalty'].nunique(dropna=True)
                .pipe(lambda s: s[s > 1].index)
            )
            if penalty_tools:
                pen_df = df[df['tool'].isin(penalty_tools)]
                # Find single best penalty per tool: lowest mean TAE across all groups
                best_penalty_per_tool = (
                    pen_df.groupby(['tool', 'penalty'])['TAE']
                    .mean()
                    .groupby('tool')
                    .idxmin()
                    .apply(lambda x: x[1])
                )
                print(f"  Best penalty per tool (fallback): {best_penalty_per_tool.to_dict()}")
                mask = pen_df.apply(lambda r: r['penalty'] == best_penalty_per_tool[r['tool']], axis=1)
                pen_df = pen_df[mask]
                no_pen_df = df[~df['tool'].isin(penalty_tools)]
                df = pd.concat([pen_df, no_pen_df], ignore_index=True)
        print(f"  WARNING: Using fallback per-replicate penalty selection. Shape: {df.shape}")

    # Load bootstrap quality assessment if provided
    df_bootstrap = None
    if bootstrap_quality:
        print(f"\nLoading bootstrap quality assessment from {len(bootstrap_quality)} file(s)...")
        df_bootstrap = load_evaluation_data(bootstrap_quality)

    generate_report(df, output, title, subtitle, df_per_sample, df_bootstrap, df_raw, results_dir, standard_penalties, unknown_sig_config, cosmic_version_label)


if 'snakemake' in dir():
    import sys
    sys.stdout = open(snakemake.log.stdout, 'w')
    sys.stderr = open(snakemake.log.stderr, 'w')
    main(
        summaries=[snakemake.input.summaries],
        per_sample=[snakemake.input.per_sample_metrics],
        bootstrap_quality=None,
        output=snakemake.output.report,
        title=snakemake.params.title,
        subtitle=snakemake.params.subtitle,
        results_dir=snakemake.params.results_dir,
        config_file=getattr(snakemake.input, 'config_file', None) or None,
    )
elif __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate HTML benchmark report for signature fitting tools',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--summaries', nargs='+', required=True)
    parser.add_argument('--per-sample', nargs='+')
    parser.add_argument('--bootstrap-quality', nargs='+')
    parser.add_argument('--output', required=True)
    parser.add_argument('--title', default='SigFitTest Benchmark Report')
    parser.add_argument('--subtitle', default='Mutational Signature Fitting Tool Comparison')
    parser.add_argument('--results-dir')
    parser.add_argument('--config')
    args = parser.parse_args()
    main(
        summaries=args.summaries,
        per_sample=args.per_sample,
        bootstrap_quality=args.bootstrap_quality,
        output=args.output,
        title=args.title,
        subtitle=args.subtitle,
        results_dir=args.results_dir,
        config_file=args.config,
    )
