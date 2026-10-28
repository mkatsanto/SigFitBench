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
import html
import json
import os
import re
from pathlib import Path
from datetime import datetime
import yaml

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DEFAULT_MIN_CONTRIBUTION = 10.0  # fallback; must match evaluate_tool_performance.py's default

# Names shown in the report for each data type
DATA_TYPE_LABELS = {'synthetic': 'Realistic simulations', 'custom': 'Synthetic'}
DATA_TYPE_NOTES = {
    'synthetic': 'Signature weights taken from real tumours (PCAWG); mutations drawn from them.',
    'custom': 'Signature weights from your own table; mutations drawn from them.',
}


# ============================================================================
# Configuration Classes - Centralize all tool and column configurations
# ============================================================================

class ReportMetadata:
    """Centralize all tool colors, symbols, column mappings, and display order."""
    
    def __init__(self):
        # Tool colors: each tool a clearly different hue; our own methods are the darkest
        self.tool_colors = {
            # our tools
            'SigTrim':                 '#1B3F8B',  # navy
            'SigID':                   '#E31A1C',  # red
            # other tools
            'SPA':                     '#1F77B4',  # blue
            'SPA_low':                 '#17BECF',  # cyan
            'FSPA':                    '#9467BD',  # purple
            'MuSiCal':                 '#FF7F0E',  # orange
            'MuSiCal-FN':              '#FFBB78',  # light orange
            'MuSiCal-BI':              '#C49C94',  # tan
            'mSigAct':                 '#2CA02C',  # green
            'StarSignDNA':             '#E377C2',  # pink
            'siglasso':                '#8C564B',  # brown
            'deconstructSigs':         '#BCBD22',  # olive
            'MutationalPatterns':      '#7F7F7F',  # grey
            'sigfit':                  '#E7BA52',  # gold
            'yapsa':                   '#AEC7E8',  # light blue
            'sigspack':                '#FF9896',  # salmon
            'signature_tools':         '#C5B0D5',  # lavender
            'mmsig':                   '#7B4173',  # plum
            'sigminer':                '#393B79',  # indigo
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
        :root {{
            --bg: #f6f7fb;
            --surface: #ffffff;
            --border: #e5e7eb;
            --border-soft: #eef0f4;
            --text: #1f2430;
            --text-muted: #667085;
            --text-faint: #9aa1b1;
            --primary: #4338ca;
            --primary-dark: #33279c;
            --primary-light: #eef1fd;
            --primary-lighter: #f5f6fe;
            --radius: 10px;
            --shadow-sm: 0 1px 2px rgba(16, 24, 40, 0.05);
            --shadow-md: 0 4px 14px rgba(16, 24, 40, 0.06);
        }}

        /* Base Styles */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        html {{
            scroll-behavior: smooth;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            font-size: 15px;
            line-height: 1.65;
            color: var(--text);
            background: var(--bg);
            -webkit-font-smoothing: antialiased;
        }}

        h1, h2, h3, h4, h5 {{
            line-height: 1.3;
            font-weight: 650;
            letter-spacing: -0.01em;
        }}

        a {{ color: var(--primary); }}

        /* Header */
        .header {{
            background: linear-gradient(120deg, #3730a3 0%, #4f46e5 55%, #4338ca 100%);
            color: white;
            padding: 2.75rem 2.5rem 2.25rem;
            box-shadow: var(--shadow-md);
        }}

        .header-inner {{
            max-width: 1400px;
            margin: 0 auto;
        }}

        .header-kicker {{
            display: inline-block;
            font-size: 0.72rem;
            font-weight: 700;
            letter-spacing: 0.12em;
            text-transform: uppercase;
            color: rgba(255,255,255,0.75);
            background: rgba(255,255,255,0.12);
            padding: 0.25rem 0.65rem;
            border-radius: 999px;
            margin-bottom: 0.9rem;
        }}

        .header h1 {{
            font-size: 2.1rem;
            font-weight: 700;
            margin-bottom: 0.4rem;
        }}

        .header p {{
            color: rgba(255,255,255,0.85);
            font-size: 1rem;
            font-weight: 400;
            max-width: 60ch;
        }}

        /* Navigation */
        .nav {{
            background: rgba(255,255,255,0.92);
            backdrop-filter: saturate(180%) blur(8px);
            border-bottom: 1px solid var(--border);
            padding: 0.7rem 2.5rem;
            position: sticky;
            top: 0;
            z-index: 100;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }}

        .nav-brand {{
            font-weight: 700;
            font-size: 0.85rem;
            color: var(--text);
            letter-spacing: -0.01em;
        }}

        .nav-links a {{
            color: var(--text-muted);
            text-decoration: none;
            margin-left: 1.75rem;
            font-weight: 600;
            font-size: 0.88rem;
            padding-bottom: 0.2rem;
            border-bottom: 2px solid transparent;
            transition: color 0.2s, border-color 0.2s;
        }}

        .nav-links a:hover {{
            color: var(--primary);
            border-bottom-color: var(--primary-light);
        }}

        /* Main Container */
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 2rem 2.5rem 3rem;
        }}

        /* Intro: quick-stats strip + orientation guide */
        .intro {{
            margin-bottom: 1.75rem;
        }}

        .stat-strip {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.75rem;
            margin-bottom: 1.25rem;
        }}

        .stat-badge {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: var(--radius);
            padding: 0.7rem 1.1rem;
            box-shadow: var(--shadow-sm);
            display: flex;
            flex-direction: column;
            gap: 0.15rem;
            min-width: 128px;
        }}

        .stat-badge-value {{
            font-size: 1.25rem;
            font-weight: 700;
            color: var(--primary);
            letter-spacing: -0.01em;
        }}

        .stat-badge-label {{
            font-size: 0.72rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-faint);
        }}

        .intro-guide {{
            background: var(--primary-lighter);
            border: 1px solid var(--primary-light);
            border-radius: var(--radius);
            padding: 1.1rem 1.4rem;
            font-size: 0.92rem;
            color: var(--text);
        }}

        .intro-guide p {{
            margin-bottom: 0.6rem;
        }}

        .intro-guide ul {{
            margin: 0;
            padding-left: 1.2rem;
            display: grid;
            gap: 0.35rem;
        }}

        .intro-guide li {{ color: var(--text); }}
        .intro-guide a {{ font-weight: 600; text-decoration: none; }}
        .intro-guide a:hover {{ text-decoration: underline; }}

        /* Sections */
        .section {{
            background: var(--surface);
            padding: 1.85rem 2rem;
            margin-bottom: 1.5rem;
            border-radius: var(--radius);
            border: 1px solid var(--border);
            box-shadow: var(--shadow-sm);
            scroll-margin-top: 4.5rem;
        }}

        .section h2 {{
            font-size: 1.35rem;
            color: var(--text);
            margin-bottom: 1.25rem;
            padding-bottom: 0.6rem;
            border-bottom: 1px solid var(--border-soft);
        }}

        .section h3 {{
            font-size: 1.05rem;
            color: var(--text);
            margin: 1.5rem 0 0.75rem;
        }}

        .section h3:first-child {{ margin-top: 0; }}

        /* Tables */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
            font-size: 0.9rem;
        }}

        th {{
            background: var(--primary-light);
            color: var(--primary-dark);
            padding: 0.65rem 0.85rem;
            text-align: left;
            font-weight: 700;
            font-size: 0.8rem;
            text-transform: uppercase;
            letter-spacing: 0.02em;
            border-bottom: 1px solid var(--border);
        }}

        td {{
            padding: 0.65rem 0.85rem;
            border-bottom: 1px solid var(--border-soft);
        }}

        tr:hover td {{
            background: var(--primary-lighter);
        }}

        /* Tabs */
        .tabs {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.25rem;
            border-bottom: 1px solid var(--border);
            margin-bottom: 1.5rem;
        }}

        .tab {{
            padding: 0.65rem 1.25rem;
            cursor: pointer;
            border: none;
            background: none;
            font-size: 0.92rem;
            font-weight: 600;
            color: var(--text-muted);
            border-bottom: 2.5px solid transparent;
            border-radius: 6px 6px 0 0;
            transition: color 0.15s, background 0.15s, border-color 0.15s;
        }}

        .tab:hover {{
            background: var(--primary-lighter);
            color: var(--primary);
        }}

        .tab.active {{
            color: var(--primary);
            border-bottom-color: var(--primary);
            background: var(--primary-lighter);
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
            background: var(--bg);
            border: 1px solid var(--border);
            padding: 0.75rem;
            border-radius: 8px;
            font-size: 13px;
        }}

        .plot-controls .controls-header {{
            margin-bottom: 0.5rem;
            position: sticky;
            top: 0;
            background: var(--bg);
            padding-bottom: 0.25rem;
            border-bottom: 1px solid var(--border);
        }}

        .plot-controls .controls-section {{
            margin-bottom: 0.75rem;
        }}

        .plot-controls label {{
            display: block;
            margin: 0.25rem 0;
            cursor: pointer;
            border-radius: 4px;
        }}

        .plot-controls label:hover {{
            background: var(--primary-light);
            padding-left: 2px;
        }}

        /* Info Box */
        .info-box {{
            background: var(--primary-lighter);
            border-left: 3px solid var(--primary);
            padding: 1.1rem 1.3rem;
            margin: 1rem 0;
            border-radius: 0 8px 8px 0;
        }}

        /* Collapsible details panels (Benchmark Overview, Metric Reference, etc.) */
        details {{
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 0 1rem;
            margin: 0.9rem 0;
            background: var(--surface);
            transition: background 0.15s;
        }}

        details[open] {{
            background: var(--surface);
        }}

        details summary {{
            list-style: none;
        }}

        details summary::-webkit-details-marker {{
            display: none;
        }}

        details summary:hover {{
            color: var(--primary) !important;
        }}

        details table, details > div {{
            margin-left: 0;
            margin-right: 0;
        }}

        /* Footer */
        .footer {{
            text-align: center;
            padding: 2.5rem 1rem;
            color: var(--text-faint);
            font-size: 0.85rem;
            border-top: 1px solid var(--border);
            margin-top: 1rem;
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
        <div class="header-inner">
            <span class="header-kicker">Benchmark Report</span>
            <h1>{title}</h1>
            <p>{subtitle}</p>
        </div>
    </div>

    <div class="nav">
        <span class="nav-brand">SigFitBench</span>
        <span class="nav-links">
            <a href="#visualizations">Results</a>
            <a href="#appendix">Appendix</a>
        </span>
    </div>

    <div class="container">
        {content}
    </div>

    <div class="footer">
        Generated by SigFitBench | {timestamp}
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


def mean_sem(df, group_cols, value_cols):
    """Group df and return mean, std, count and SEM (std / sqrt(count)) per value column.

    Output columns: <col>_mean, <col>_std, <col>_count, <col>_sem.
    """
    out = df.groupby(group_cols)[value_cols].agg(['mean', 'std', 'count'])
    out.columns = [f'{c}_{stat}' for c, stat in out.columns]
    out = out.reset_index()
    for c in value_cols:
        out[f'{c}_sem'] = out[f'{c}_std'] / np.sqrt(out[f'{c}_count'])
    return out


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

def create_report_intro(df, unknown_sig_config=None, tool_overrides=None):
    """Compact orientation block shown right under the header, before the results tabs:
    a row of quick-stat badges plus a short guide to what's in the report and where to
    find it. Keeps readers from having to scroll past a page of background/methodology
    text (that content still exists — see create_appendix_section — just moved to the
    end of the report as reference material).

    Args:
        df: DataFrame with evaluation metrics
        unknown_sig_config: Optional dict with unknown signature configuration
                          {'num_unknowns': int, 'weights': list}
        tool_overrides: Optional dict {tool_name: {param: value}} from the config's
                        `tool_overrides` section (parameters set explicitly for a tool)
    """
    n_tools = df['tool'].nunique()
    cancer_types_list = sorted(df['cancer_type'].unique())
    n_cancer_types = len(cancer_types_list)
    num_samples = int(df['num_samples'].iloc[0]) if 'num_samples' in df.columns else None

    mutation_counts = sorted(df['num_muts'].unique())
    if mutation_counts:
        burden_range = f"{int(mutation_counts[0]):,} – {int(mutation_counts[-1]):,}"
    else:
        burden_range = 'N/A'

    badges = [
        ('Tools Compared', str(n_tools)),
        ('Cancer Types', str(n_cancer_types)),
        ('Mutation Burdens', burden_range),
        ('Samples / Condition', str(num_samples) if num_samples is not None else 'N/A'),
    ]

    if unknown_sig_config and unknown_sig_config.get('num_unknowns', 0) > 0:
        weights = unknown_sig_config.get('weights', [])
        weight_label = ', '.join(f"{w*100:.0f}%" for w in weights) if weights else \
            f"{unknown_sig_config['num_unknowns']} sig(s)"
        badges.append(('Unknown Sig. Weight', weight_label))

    badges_html = ''.join(
        f'<div class="stat-badge"><span class="stat-badge-value">{value}</span>'
        f'<span class="stat-badge-label">{label}</span></div>'
        for label, value in badges
    )

    overrides_html = ''
    if tool_overrides:
        overrides = {k.lower(): v for k, v in tool_overrides.items() if v}
        items = []
        for tool in sorted(df['tool'].unique()):
            params = overrides.get(str(tool).lower())
            if params:
                param_str = ', '.join(f"<code>{k}={v}</code>" for k, v in params.items())
                items.append(f"<strong>{tool}</strong> ({param_str})")
        if items:
            overrides_html = (
                '<p>Tools run with explicitly specified (non-default) parameters: '
                + '; '.join(items) + '. All other tools use their default settings.</p>'
            )

    cancer_types_str = ', '.join(cancer_types_list)
    plural = 's' if n_cancer_types != 1 else ''

    html = f'''
    <div class="intro">
        <div class="stat-strip">{badges_html}</div>
        <div class="intro-guide">
            <p>
                This report compares <strong>{n_tools} signature-fitting tools</strong> on synthetic tumour samples with a
                <strong>known ground truth</strong>, across <strong>{n_cancer_types} cancer type{plural}</strong>
                ({cancer_types_str}) and mutation burdens from {burden_range}.
            </p>
            {overrides_html}
            <ul>
                <li>Each tab below is one cancer type's results: precision/recall, error metrics, tool agreement, and per-sample signature activity.</li>
                <li>The <strong>Aggregate Analysis</strong> tab pools results across all cancer types.</li>
                <li>Full background, parameter tables, and metric definitions are in the <a href="#appendix">Appendix</a> at the end of this report.</li>
            </ul>
        </div>
    </div>
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


def create_resource_usage_section(df_benchmarks, colors):
    """Wall time and peak memory per tool, from Snakemake's `benchmark:` directive
    (see tools.smk / aggregate_tool_benchmarks.py). One row of panels per data type; within
    it, cancer types, relevance modes and replicates are pooled — for each tool and mutation
    burden, points are the mean over every job run at that burden, with error bars at +-1 std
    (the spread between jobs, not the uncertainty of the mean). Each job fits a whole
    sample matrix in one process, so there is no per-sample timing to show spread over.
    Returns '' if no benchmark data is available.
    """
    if df_benchmarks is None or len(df_benchmarks) == 0:
        return ''
    if 'seconds' not in df_benchmarks.columns or 'max_rss_mb' not in df_benchmarks.columns:
        return ''

    tools = get_tool_display_order()
    tools = [t for t in tools if t in set(df_benchmarks['tool'].unique())] + \
            sorted(set(df_benchmarks['tool'].unique()) - set(tools))
    data_types = sorted(df_benchmarks['data_type'].unique())

    agg = mean_sem(df_benchmarks, ['data_type', 'tool', 'num_muts'], ['seconds', 'max_rss_mb'])

    fig = make_subplots(rows=len(data_types), cols=2, horizontal_spacing=0.12,
                        vertical_spacing=0.12 if len(data_types) > 1 else 0.0,
                        subplot_titles=[f'{title} — {dt}' if len(data_types) > 1 else title
                                        for dt in data_types
                                        for title in ('Wall time', 'Peak memory')])
    method_symbols = get_method_symbols()
    for row, data_type in enumerate(data_types, start=1):
        for tool in tools:
            tool_data = agg[(agg['data_type'] == data_type) & (agg['tool'] == tool)].sort_values('num_muts')
            if tool_data.empty:
                continue
            symbol = method_symbols.get(get_tool_base_method(tool), 'circle')
            color = colors.get(tool, '#999')
            for col, metric, unit in [(1, 'seconds', 's'), (2, 'max_rss_mb', 'MB')]:
                fig.add_trace(go.Scatter(
                    x=tool_data['num_muts'], y=tool_data[f'{metric}_mean'], mode='lines+markers',
                    name=tool, legendgroup=tool, showlegend=(row == 1 and col == 1),
                    line=dict(color=color, width=2), marker=dict(color=color, size=9, symbol=symbol),
                    error_y=dict(array=tool_data[f'{metric}_std'].fillna(0), type='data'),
                    customdata=np.stack([tool_data[f'{metric}_std'].fillna(0),
                                         tool_data[f'{metric}_count']], axis=-1),
                    hovertemplate=(f'<b>{tool}</b> ({data_type})<br>%{{x:,}} mutations<br>'
                                   f'%{{y:.3g}} ± %{{customdata[0]:.2g}} {unit} (std)<br>'
                                   f'n = %{{customdata[1]}} jobs<extra></extra>'),
                ), row=row, col=col)
        fig.update_xaxes(type='log', title_text='Mutations', row=row, col=1)
        fig.update_xaxes(type='log', title_text='Mutations', row=row, col=2)
        fig.update_yaxes(title_text='Seconds', row=row, col=1)
        fig.update_yaxes(title_text='MB (max RSS)', row=row, col=2)

    fig.update_layout(
        title_text='Resource Usage by Mutation Burden<br><sub>Mean ± std across jobs '
                   '(cancer types / replicates); hover for n</sub>',
        height=480 * len(data_types), hovermode='closest', autosize=True,
    )

    html = '<section id="resources" class="section">\n'
    html += '<h2>⚡ Resource Usage</h2>\n'
    html += (
        '<p style="color:#666; margin-bottom:1rem;">Wall time and peak memory (max RSS) per '
        'tool, recorded by Snakemake for every fitting job actually run, shown separately per '
        'data type and pooled across cancer types and replicates. Each job fits one whole '
        'sample matrix, so error bars (±1 std) show the spread between jobs, not between '
        'individual samples. Reflects this benchmark\'s hardware and job scheduling '
        '(local vs. cluster), not each tool\'s theoretical best case.</p>\n'
    )
    html += '<div class="plot-container">\n'
    html += fig.to_html(include_plotlyjs=False, full_html=False, div_id='resource_usage_plot',
                        config={'responsive': True})
    html += '</div>\n'
    html += '</section>\n'
    return html


def create_appendix_section(df, df_raw=None, unknown_sig_config=None, cosmic_version_label=None,
                             min_contribution=None):
    """Background, full parameter tables, metric reference, and experimental setup —
    reference material shown at the end of the report (see create_report_intro for the
    compact quick-orientation version shown at the top instead).

    Args:
        df: DataFrame with evaluation metrics
        df_raw: Optional raw per-realization data
        unknown_sig_config: Optional dict with unknown signature configuration
                          {'num_unknowns': int, 'weights': list}
        cosmic_version_label: Optional human-readable COSMIC catalog label derived
                          from config's cosmic_signatures.base_file
        min_contribution: Optional float — config's evaluation.min_contribution, the
                          mutation-count threshold a signature's estimated contribution
                          must clear to count as "active" for P/R/F1/MCC/FP/FN and for
                          the Signature Activity widget's predicted markers.
    """
    n_tools = df['tool'].nunique()
    n_cancer_types = df['cancer_type'].nunique()
    n_conditions = len(df)
    num_samples = int(df['num_samples'].iloc[0]) if 'num_samples' in df.columns else 'N/A'

    num_realizations = 'N/A'
    if df_raw is not None and 'realization' in df_raw.columns:
        num_realizations = int(df_raw['realization'].nunique())

    # Get unique tools in custom order, with display-friendly capitalization
    tool_order = get_tool_display_order()
    available_tools = set(df['tool'].unique())
    tools_list = [t for t in tool_order if t in available_tools]
    remaining_tools = sorted(available_tools - set(tools_list))
    tools_list.extend(remaining_tools)
    tool_display_names = []
    for tool in tools_list:
        if tool.lower() == 'musical':
            tool_display_names.append('MuSiCal')
        elif tool.lower() == 'spa':
            tool_display_names.append(tool.upper())
        else:
            tool_display_names.append(tool)
    tools_str = ', '.join(tool_display_names)

    cancer_types_str = ', '.join(sorted(df['cancer_type'].unique()))

    mutation_counts = sorted(df['num_muts'].unique())
    mut_counts_str = ', '.join([f"{int(m):,}" for m in mutation_counts])

    data_type = df['data_type'].iloc[0] if 'data_type' in df.columns else 'N/A'
    data_type_note = DATA_TYPE_NOTES.get(data_type, '')
    data_type = DATA_TYPE_LABELS.get(data_type, data_type)
    relevance_val = df['relevance'].iloc[0] if 'relevance' in df.columns else 'N/A'

    if unknown_sig_config and unknown_sig_config.get('num_unknowns', 0) > 0:
        num_unknowns = unknown_sig_config['num_unknowns']
        weights = unknown_sig_config.get('weights', [])
        unknown_weights_str = ', '.join(f"{w*100:.1f}%" for w in weights) if weights else 'N/A'
    else:
        num_unknowns = 0
        unknown_weights_str = 0

    cosmic_version_label = cosmic_version_label or "COSMIC SBS signatures"
    min_contribution_str = f"{min_contribution:g} mutations" if min_contribution is not None else "N/A"

    # Build unknown signatures description from config (for the Data Generation table below)
    if unknown_sig_config and unknown_sig_config.get('num_unknowns', 0) > 0:
        weights = unknown_sig_config.get('weights', [])
        unknown_sig_description = f"{num_unknowns} unknown signature(s) per sample, randomly selected from excluded COSMIC signatures"
        if weights:
            weights_pct = [f"{w*100:.1f}%" for w in weights]
            unknown_sig_description += f"<br><strong>Contribution weights:</strong> {', '.join(weights_pct)} of total mutations per unknown signature"
    else:
        if 'num_unknown_sigs' in df.columns:
            inferred_unknowns = int(df['num_unknown_sigs'].iloc[0])
            if inferred_unknowns > 0:
                unknown_sig_description = f"{inferred_unknowns} unknown signature(s) per sample, randomly selected from excluded COSMIC signatures"
            else:
                unknown_sig_description = "No unknown signatures (all signatures available in reference catalog)"
        else:
            unknown_sig_description = "Selected signatures excluded from reference catalog (randomly chosen per sample)"

    html = f'''
    <section id="appendix" class="section">
        <h2>Appendix: Background &amp; Reference</h2>

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
                <td>{n_tools} ({tools_str})</td>
                <td style="color: #555; font-size: 0.9em;">The signature fitting algorithms being benchmarked head-to-head.</td>
            </tr>
            <tr>
                <td><strong>Cancer Types</strong></td>
                <td>{n_cancer_types} ({cancer_types_str})</td>
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
                <td style="color: #555; font-size: 0.9em;">{data_type_note}</td>
            </tr>
            <tr>
                <td><strong>Signature Relevance</strong></td>
                <td>{relevance_val}</td>
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
                <p style="margin: 0 0 0.75rem 0; padding: 0.6rem 0.8rem; background: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; color: #856404;">
                    <strong>Active-signature threshold:</strong> a signature counts as "present" in a sample only if its
                    estimated contribution is &ge; <strong>{min_contribution_str}</strong> (config <code>evaluation.min_contribution</code>).
                    Everything below that line is treated as absent — it's what determines P / R / F1 / MCC / FP / FN
                    below, and predicted signatures below this threshold are hidden from the Signature Activity widget's markers.
                </p>
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

        <h3 style="margin-top: 2rem;">Experimental Setup</h3>

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
                <tr>
                    <td><strong>Active-Signature Threshold</strong></td>
                    <td>{min_contribution_str} (<code>evaluation.min_contribution</code>) — see Metric Reference above</td>
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

# Candidate colors for tools missing from ReportMetadata.tool_colors (Plotly Dark24 + Alphabet)
_FALLBACK_COLORS = [
    '#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86',
    '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE',
    '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038',
    '#AA0DFE', '#3283FE', '#85660D', '#782AB6', '#565656', '#1C8356', '#16FF32', '#F7E1A0',
    '#E2E2E2', '#1CBE4F', '#C4451C', '#DEA0FD', '#FE00FA', '#325A9B', '#FEAF16', '#F8A19F',
    '#90AD1C', '#F6222E', '#1CFFCE', '#2ED9FF', '#B10DA1', '#C075A6', '#FC1CBF', '#B00068',
]


def _hex_to_rgb(hex_color):
    h = hex_color.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def get_tool_colors(tools=None):
    """Return the tool -> hex color map, shared by all reports.

    Tools in ReportMetadata.tool_colors keep their fixed color. Any other tool in `tools`
    is assigned (in alphabetical order, so deterministic) the fallback color that is
    farthest in RGB space from every color already in use, so new tools stay
    distinguishable from existing ones and from each other.
    """
    colors = dict(ReportMetadata().tool_colors)
    used = [_hex_to_rgb(c) for c in colors.values()]
    for tool in sorted(set(tools or []) - set(colors)):
        def min_dist(candidate):
            r = _hex_to_rgb(candidate)
            return min(sum((a - b) ** 2 for a, b in zip(r, u)) for u in used)
        best = max(_FALLBACK_COLORS, key=min_dist)
        colors[tool] = best
        used.append(_hex_to_rgb(best))
    return colors


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


def _load_sample_signature_tables(df_samples, cancer_type, results_dir):
    """Per mutation burden, one row per simulated sample with its mean TAE / F1 / ENP across
    tools and its ground-truth signature set (frozenset of signatures with weight > 0, from
    results_dir/generated_data/{data_type}/weights_{cancer}-w_{rep}-{burden}.dat — "synthetic"
    and "custom" both have ground truth and use the same layout).

    Returns a list of (burden, DataFrame[replicate, sample, TAE, F1?, highest_ENP_known?, sigs])
    sorted by burden (hardest sample first), or None if metrics or ground-truth weights are
    unavailable.
    """
    if not results_dir or 'TAE' not in df_samples.columns:
        return None
    cancer_samples = df_samples[df_samples['cancer_type'] == cancer_type]
    if cancer_samples.empty:
        return None
    has_rep = 'replicate' in cancer_samples.columns
    metric_cols = [c for c in ['TAE', 'F1', 'highest_ENP_known'] if c in cancer_samples.columns]
    data_type = cancer_samples['data_type'].iloc[0] if 'data_type' in cancer_samples.columns else 'synthetic'
    weights_dir = Path(results_dir) / 'generated_data' / data_type

    tables = []
    for burden in sorted(cancer_samples['num_muts'].unique()):
        rows = cancer_samples[cancer_samples['num_muts'] == burden]
        group_cols = ['replicate', 'sample'] if has_rep else ['sample']
        tae = rows.groupby(group_cols)[metric_cols].mean().reset_index()
        if not has_rep:
            tae['replicate'] = 0
        presence = {}
        for rep in tae['replicate'].unique():
            path = weights_dir / f'weights_{cancer_type}-w_{int(rep)}-{int(burden)}.dat'
            if not path.exists():
                continue
            w = pd.read_csv(path, sep='\t', index_col=0) > 0
            for sample in w.columns:
                presence[(rep, sample)] = frozenset(w.index[w[sample]])
        tae['sigs'] = [presence.get((r, sm)) for r, sm in zip(tae['replicate'], tae['sample'])]
        tae = tae.dropna(subset=['sigs'])
        if len(tae) < 2:
            continue
        tables.append((burden, tae.sort_values('TAE', ascending=False).reset_index(drop=True)))
    return tables or None


def _bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values (q-values), same order as the input."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for rank, idx in zip(range(n, 0, -1), order[::-1]):
        running = min(running, p[idx] * n / rank)
        q[idx] = running
    return q


def create_hard_samples_row_plot(df_samples, cancer_type, results_dir, out_of_reference_sigs=None,
                                  top_n=15, top_k=12, min_group_n=5):
    """One row of panels, sharing a single mutation-burden dropdown (menu sits outside the
    plot area, to the right, so it never overlaps the subplot titles):
      1. The `top_n` hardest samples — ranked by mean TAE across all tools (pooling tools,
         not per-tool). TAE is used because it measures how much of the correct signal was
         misattributed, a more direct "how wrong" measure than F1.
      2. Which ground-truth signatures are over-represented among those hardest samples
         (orange) compared with their frequency across all samples (grey). Unknown
         (out-of-reference) signatures are outlined in red instead of labelled, to keep
         the axis readable.
      3. The effect on mean TAE of each signature being present vs. absent across all
         samples (Mann-Whitney U test, Benjamini-Hochberg corrected; q < 0.05 marked '*').
         Bars for unknown signatures are red, known signatures blue.
    No in-plot legend is shown (color meaning is fixed across panels and described in the
    caller's caption); this keeps the row compact enough that the burden dropdown never
    overlaps a subplot.
      4. (only if highest_ENP_known is available) Difficulty (mean TAE) vs. ENP — the
         spectral flatness of the sample's most diffuse known signature — hardest samples
         highlighted in red, to show whether hard samples share a diffuse dominant signature.
    All panels use each sample's simulated ground-truth signature set (present = weight > 0,
    from results_dir/generated_data/{data_type}/weights_{cancer}-w_{rep}-{burden}.dat).
    Returns None if per-sample TAE or the ground-truth weights are unavailable.
    """
    from scipy.stats import mannwhitneyu

    tables = _load_sample_signature_tables(df_samples, cancer_type, results_dir)
    if tables is None:
        return None
    unknown = set(out_of_reference_sigs or [])
    has_enp = all('highest_ENP_known' in tae.columns for _, tae in tables)
    n_cols = 4 if has_enp else 3

    def test_features(values, sigs_list, names):
        results = []
        for name in names:
            mask = np.array([name in x for x in sigs_list])
            n1, n0 = int(mask.sum()), int((~mask).sum())
            if n1 < min_group_n or n0 < min_group_n:
                continue
            g1, g0 = values[mask], values[~mask]
            try:
                p = float(mannwhitneyu(g1, g0, alternative='two-sided').pvalue)
            except ValueError:
                p = 1.0
            results.append(dict(name=name, effect=float(g1.mean() - g0.mean()),
                                m1=float(g1.mean()), m0=float(g0.mean()), n1=n1, n0=n0, p=p))
        if results:
            for r, q in zip(results, _bh_adjust([r['p'] for r in results])):
                r['q'] = float(q)
        return sorted(results, key=lambda r: abs(r['effect']), reverse=True)[:top_k][::-1]

    titles = [f'Hardest {top_n} samples', 'Over-represented signatures', 'Effect on TAE']
    if has_enp:
        titles.append('Difficulty vs. ENP')
    fig = make_subplots(rows=1, cols=n_cols, horizontal_spacing=0.06, subplot_titles=titles)
    fig.update_annotations(font_size=12)

    trace_burdens = []
    for i, (burden, tae) in enumerate(tables):
        visible = (i == 0)
        n_all = len(tae)
        hardest = tae.head(min(top_n, n_all))

        # Panel 1: ranked bar of the hardest samples (reversed so #1 renders at the top).
        hardest_display = hardest.iloc[::-1]
        tae_max = float(tae['TAE'].max())
        fig.add_trace(go.Bar(
            x=hardest_display['TAE'], y=[str(s) for s in hardest_display['sample']], orientation='h',
            marker=dict(color=hardest_display['TAE'], colorscale='RdYlGn_r', cmin=0, cmax=tae_max,
                       showscale=False, line=dict(color='black', width=0.5)),
            text=[f"{t:.2f}" for t in hardest_display['TAE']], textposition='outside',
            hovertemplate='Sample %{y}<br>Mean TAE: %{x:.3f}<extra></extra>',
            showlegend=False, visible=visible), row=1, col=1)
        trace_burdens.append(burden)

        # Panel 2: signature frequency, hardest samples vs. all samples. Unknown signatures
        # get a red outline instead of a "(unknown)" label suffix, to keep names short.
        all_sigs = sorted(set().union(*tae['sigs']))
        pct_all = {g: 100 * sum(g in x for x in tae['sigs']) / n_all for g in all_sigs}
        pct_hard = {g: 100 * sum(g in x for x in hardest['sigs']) / len(hardest) for g in all_sigs}
        top = sorted(all_sigs, key=lambda g: pct_hard[g] - pct_all[g], reverse=True)[:top_k][::-1]
        line_colors = ['#dc2626' if g in unknown else 'rgba(0,0,0,0)' for g in top]
        line_widths = [2 if g in unknown else 0 for g in top]
        fig.add_trace(go.Bar(
            y=top, x=[pct_all[g] for g in top], orientation='h', name='All samples',
            marker=dict(color='#9ca3af', line=dict(color=line_colors, width=line_widths)),
            showlegend=False, visible=visible,
            hovertemplate='%{y}<br>%{x:.0f}% of all samples<extra></extra>'), row=1, col=2)
        fig.add_trace(go.Bar(
            y=top, x=[pct_hard[g] for g in top], orientation='h', name=f'Hardest {top_n}',
            marker=dict(color='#f97316', line=dict(color=line_colors, width=line_widths)),
            showlegend=False, visible=visible,
            hovertemplate='%{y}<br>%{x:.0f}% of hardest samples<extra></extra>'), row=1, col=2)
        trace_burdens += [burden, burden]

        # Panel 3: present-vs-absent effect on TAE, per signature. Red = unknown, blue = known.
        values = tae['TAE'].to_numpy(dtype=float)
        res = test_features(values, list(tae['sigs']), all_sigs)
        bar_colors = ['#dc2626' if r['name'] in unknown else '#2563eb' for r in res]
        fig.add_trace(go.Bar(
            y=[r['name'] for r in res], x=[r['effect'] for r in res], orientation='h',
            marker_color=bar_colors, text=['*' if r['q'] < 0.05 else '' for r in res], textposition='outside',
            customdata=[[r['m1'], r['m0'], r['n1'], r['n0'], r['p'], r['q']] for r in res],
            hovertemplate=('%{y}<br>Effect (present − absent): %{x:+.3f}'
                           '<br>Mean TAE present: %{customdata[0]:.3f} (n=%{customdata[2]})'
                           '<br>Mean TAE absent: %{customdata[1]:.3f} (n=%{customdata[3]})'
                           '<br>p=%{customdata[4]:.3g}, q=%{customdata[5]:.3g}<extra></extra>'),
            showlegend=False, visible=visible), row=1, col=3)
        trace_burdens.append(burden)

        # Panel 4 (optional): difficulty vs. ENP of the sample's most diffuse known signature.
        if has_enp:
            hard_ids = set(hardest['sample'])
            is_hard = tae['sample'].isin(hard_ids)
            other, hard = tae[~is_hard], tae[is_hard]
            fig.add_trace(go.Scatter(
                x=other['highest_ENP_known'], y=other['TAE'], mode='markers',
                marker=dict(color='#9ca3af', size=6, opacity=0.6), name='Other samples',
                showlegend=False, visible=visible,
                hovertemplate='Sample %{text}<br>ENP: %{x:.1f}<br>Mean TAE: %{y:.3f}<extra></extra>',
                text=[str(s) for s in other['sample']]), row=1, col=4)
            fig.add_trace(go.Scatter(
                x=hard['highest_ENP_known'], y=hard['TAE'], mode='markers',
                marker=dict(color='#dc2626', size=8, line=dict(color='black', width=1)),
                name=f'Hardest {top_n}', showlegend=False, visible=visible,
                hovertemplate='Sample %{text}<br>ENP: %{x:.1f}<br>Mean TAE: %{y:.3f}<extra></extra>',
                text=[str(s) for s in hard['sample']]), row=1, col=4)
            trace_burdens += [burden, burden]

    buttons = [dict(
        label=f"{int(b):,} mutations", method='update',
        args=[{'visible': [tb == b for tb in trace_burdens]},
              {'title': {'text': f'Hardest Samples — {int(b):,} Mutations'}}])
        for b, _ in tables]
    fig.update_layout(
        title_text=f'Hardest Samples — {int(tables[0][0]):,} Mutations',
        barmode='group', height=max(420, 28 * max(top_n, top_k)),
        margin=dict(l=110, r=140, t=80, b=60),
        # Menu placed outside the plot area (x > 1) so it can never overlap a subplot title.
        updatemenus=[dict(buttons=buttons, direction='down', pad={'r': 5, 't': 5}, showactive=True,
                          x=1.02, xanchor='left', y=1.0, yanchor='top',
                          bgcolor='white', bordercolor='#ccc', borderwidth=1)],
        autosize=True,
    )
    fig.update_xaxes(title_text='Mean TAE', row=1, col=1)
    fig.update_xaxes(title_text='% of samples', row=1, col=2)
    fig.update_xaxes(title_text='Δ mean TAE (present − absent)', row=1, col=3, zeroline=True, zerolinecolor='black')
    if has_enp:
        fig.update_xaxes(title_text='ENP (known)', row=1, col=4)
        fig.update_yaxes(title_text='Mean TAE', row=1, col=4)
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
            denom = 2 * recall - f1_val
            with np.errstate(divide='ignore', invalid='ignore'):
                precision = np.where(denom > 0, f1_val * recall / denom, np.nan)
            precision = np.where((precision >= 0) & (precision <= 1), precision, np.nan)
            
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


def create_metrics_grid(df, cancer_type, tools, colors, df_per_sample=None, pooled=False):
    """
    Create 3x3 grid (9 panels) of performance metrics vs mutation count, including the
    % of mutations left unassigned (1 - weight_tot). Shows mean across all
    samples/realizations with SEM as error bars.

    Args:
        df: DataFrame with evaluation metrics
        cancer_type: Cancer type name (for display; if 'cancer_type' column not present, uses df as-is)
        tools: List of tool names
        colors: Dict mapping tool names to colors
        df_per_sample: DataFrame with per-sample metrics for error bars
        pooled: If True, this is a cross-cancer-type aggregate view -- df and df_per_sample are
            used as-is (no cancer_type filtering), since `cancer_type` is just a display label
            here (e.g. "All Cancer Types (Pooled)") and won't match any real row. Without this,
            filtering df_per_sample by that label matches zero rows and silently drops every
            tool/metric from the grid.
    """
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    # Filter by cancer type if the column exists, otherwise use the dataframe as-is
    if pooled or 'cancer_type' not in df.columns:
        cancer_df = df
    else:
        cancer_df = df[df['cancer_type'] == cancer_type]

    # Create 3x3 subplot grid
    fig = make_subplots(
        rows=3, cols=3,
        subplot_titles=['F1 ↑', 'TAE ↓', 'MCC ↑', 'nRMSE ↓', 'FP Sigs', 'FN Sigs',
                        '% Unassigned Mutations', 'Precision ↑', 'Recall ↑'],
        vertical_spacing=0.1,
        horizontal_spacing=0.12
    )

    metric_mapping = [
        ('F1', ['F1', 'F1_mean'], 1, 1),
        ('TAE', ['TAE', 'TAE_mean'], 1, 2),
        ('MCC', ['MCC', 'MCC_mean'], 1, 3),
        ('nRMSE', ['nRMSE', 'nRMSE_mean'], 2, 1),
        ('num_FP', ['num_FP_sigs', 'num_FP'], 2, 2),
        ('num_FN', ['num_FN_sigs', 'num_FN'], 2, 3),
        ('pct_unassigned', ['pct_unassigned'], 3, 1),
        ('P', ['P', 'P_mean'], 3, 2),
        ('R', ['R', 'R_mean'], 3, 3),
    ]

    # Use per-sample data if available for proper error bars
    use_per_sample = df_per_sample is not None and len(df_per_sample) > 0
    if use_per_sample:
        # Filter per-sample data by cancer type if the column exists
        if pooled or 'cancer_type' not in df_per_sample.columns:
            cancer_per_sample = df_per_sample
        else:
            cancer_per_sample = df_per_sample[df_per_sample['cancer_type'] == cancer_type]
        if 'weight_tot' in cancer_per_sample.columns:
            # weight_tot = assigned / observed mutations, so unassigned % = (1 - weight_tot) * 100
            # (negative when a tool assigns more mutations than the sample has)
            cancer_per_sample = cancer_per_sample.assign(
                pct_unassigned=(1 - cancer_per_sample['weight_tot']) * 100)

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
        fig.update_xaxes(type='log', title_text='Mutations', row=3, col=col)
    fig.update_yaxes(title_text='% of mutations', row=3, col=1)
    
    fig.update_layout(
        height=1000,
        title_text=f'{cancer_type} - Performance Metrics<br><sub>Mean ± SEM across all samples vs mutation burden</sub>',
        showlegend=True,
        legend=dict(x=1.05, y=1, itemdoubleclick='toggleothers'),
        hovermode='closest',
        autosize=True
    )
    
    return fig


def create_aggregate_plot(df, tools, colors, df_per_sample=None):
    """Create aggregate precision-recall plot across all cancer types.

    Error bars are +-1 SEM. When per-sample data is available, SEM is computed over the
    individual samples (std / sqrt(n samples)); otherwise it falls back to the SEM of the
    per-group means in `df` (i.e. across cancer types/replicates, not samples).
    """
    if df_per_sample is not None and len(df_per_sample) > 0:
        agg_df = mean_sem(df_per_sample, ['tool', 'num_muts'], ['P', 'R', 'F1'])
    else:
        group_means = df.rename(columns={'P_mean': 'P', 'R_mean': 'R', 'F1_mean': 'F1'})
        agg_df = mean_sem(group_means, ['tool', 'num_muts'], ['P', 'R', 'F1'])

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

    # F1 reference curves (same grey dashed contours as the per-cancer-type PR plot).
    for f1_val, f1_label in [(0.5, 'F1=0.5'), (0.7, 'F1=0.7'), (0.9, 'F1=0.9')]:
        recall = np.linspace(0.01, 1, 100)
        denom = 2 * recall - f1_val
        precision = np.where(denom > 0, f1_val * recall / denom, np.nan)
        precision = np.clip(precision, 0, 1)
        fig.add_trace(go.Scatter(
            x=recall, y=precision, mode='lines', name=f1_label,
            line=dict(dash='dash', color='gray', width=1), opacity=0.5, hoverinfo='name'
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


def create_feature_metric_scatter_widget(
        features_tsv_path: str,
        per_sample_metrics_tsv: str,
        widget_id: str = 'feature_metric_widget',
        data_df=None,
        data_json: str = None,
        colors: dict = None
) -> str:
    """Return HTML + JS for interactive boxplot: feature binned low/high vs metric distribution.

    Args:
        features_tsv_path: path to TSV file with scalar features per sample (DEPRECATED if data_df provided)
        per_sample_metrics_tsv: path to TSV with per-sample metrics (DEPRECATED if data_df provided)
        widget_id: unique DOM id prefix for this widget
        data_df: pandas DataFrame with all data (preferred over file paths to avoid CORS issues)
        data_json: pre-serialized {header, rows} JSON string; if given, skips serializing data_df again
            (callers embedding this widget once per tab should serialize the shared data once and reuse it)
        colors: tool -> hex color mapping, e.g. from get_tool_colors(); falls back to gray if omitted

    The function produces a small interactive UI where the user selects:
        - feature (to bin into low/high)
        - metric (Y-axis for boxplot)
        - mutation burden filter

    Shows separate boxplots for each tool, split by low vs high feature values.
    Output: HTML string to embed in report.
    """
    html = f"""
<div style="background: white; padding: 1.5rem; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); margin: 1rem 0;">
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1rem; margin-bottom: 1.5rem;">
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">
                📊 Feature (binned Low/High)
            </label>
            <select id="{widget_id}_feature" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;">
                <option value="">Select feature...</option>
            </select>
        </div>
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">
                📈 Metric (distribution)
            </label>
            <select id="{widget_id}_metric" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;">
                <option value="">Select metric...</option>
            </select>
        </div>
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">
                🧬 Mutation Burden
            </label>
            <select id="{widget_id}_muts" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;">
                <option value="">All</option>
            </select>
        </div>
    </div>

    <div style="text-align: center; margin-bottom: 1rem;">
        <button id="{widget_id}_update" style="
            background: linear-gradient(135deg, #7c3aed 0%, #5b21b6 100%);
            color: white;
            border: none;
            padding: 0.75rem 2rem;
            border-radius: 8px;
            font-size: 16px;
            font-weight: 600;
            cursor: pointer;
            box-shadow: 0 4px 6px rgba(124, 58, 237, 0.3);
            transition: all 0.3s ease;
        " onmouseover="this.style.transform='translateY(-2px)'; this.style.boxShadow='0 6px 12px rgba(124, 58, 237, 0.4)';"
           onmouseout="this.style.transform='translateY(0)'; this.style.boxShadow='0 4px 6px rgba(124, 58, 237, 0.3)';">
            🎯 Update Plot
        </button>
    </div>

    <div id="{widget_id}_plot" style="width: 100%; height: 600px; border: 2px solid #f3f4f6; border-radius: 8px; background: #fafafa;"></div>

    <div id="{widget_id}_stats" style="margin-top: 1rem; padding: 1rem; background: #f9fafb; border-radius: 6px; display: none;">
        <div style="display: flex; justify-content: center; align-items: center; gap: 2rem; font-size: 14px;">
            <div>
                <span style="color: #6b7280;">Median threshold:</span>
                <strong id="{widget_id}_threshold" style="color: #111827; margin-left: 0.5rem;">-</strong>
            </div>
            <div>
                <span style="color: #6b7280;">N (Low):</span>
                <strong id="{widget_id}_n_low" style="color: #111827; margin-left: 0.5rem;">-</strong>
            </div>
            <div>
                <span style="color: #6b7280;">N (High):</span>
                <strong id="{widget_id}_n_high" style="color: #111827; margin-left: 0.5rem;">-</strong>
            </div>
        </div>
    </div>
</div>

<script>
(function() {{
    const widgetId = '{widget_id}';

    // Feature columns (complexity features)
    const featureCols = ['highest_ENP_known', 'highest_ENP_unknown'];

    // Metric columns (performance metrics)
    const metricCols = ['TAE', 'F1', 'MCC', 'P', 'R', 'theoretical_min_TAE'];

    // Metadata columns
    const metadataCols = ['data_type', 'cancer_type', 'num_muts_condition', 'relevance',
                          'tool', 'realization', 'sample', 'num_muts'];

    // Embedded data (avoids CORS issues with file:// protocol)
    const embeddedData = __DATA__;

    const features = embeddedData;
    const metrics = embeddedData;  // Same data source

    // Build feature list (only feature columns)
    const featureSelect = document.getElementById(widgetId + '_feature');

    features.header.forEach(h => {{
        if (featureCols.includes(h)) {{
            const opt = document.createElement('option');
            opt.value = h;
            opt.text = h;
            featureSelect.appendChild(opt);
        }}
    }});

    // Build metric list (only metric columns)
    const metricSelect = document.getElementById(widgetId + '_metric');

    metrics.header.forEach(h => {{
        if (metricCols.includes(h)) {{
            const opt = document.createElement('option');
            opt.value = h;
                opt.text = h;
                metricSelect.appendChild(opt);
            }}
        }});

        // Build mutation burdens list
        const toolIdx = metrics.header.indexOf('tool');
        const mutsIdx = metrics.header.indexOf('num_muts');
        const muts = new Set();
        if (mutsIdx >= 0) {{
            metrics.rows.forEach(r => {{
                if (r[mutsIdx]) muts.add(r[mutsIdx]);
            }});
        }}
        const mutsSelect = document.getElementById(widgetId + '_muts');
        Array.from(muts).sort((a,b)=>+a-+b).forEach(m => {{
            const opt = document.createElement('option');
            opt.value = m;
            opt.text = m;
            mutsSelect.appendChild(opt);
        }});

        // Build data structure for plotting
        const sampleIdx = metrics.header.indexOf('sample');
        const dataByToolMuts = new Map(); // tool_muts -> array of {{feature_vals, metric_vals, sample}}

        metrics.rows.forEach((row, idx) => {{
            const sample = row[sampleIdx];
            const tool = toolIdx >= 0 ? row[toolIdx] : '';
            const numMuts = mutsIdx >= 0 ? row[mutsIdx] : '';
            const key = tool + '_' + numMuts;

            // Get feature values for this row
            const featureVals = {{}};
            features.header.forEach((h, i) => {{
                if (featureCols.includes(h)) {{
                    featureVals[h] = row[i];  // Same row if same file
                }}
            }});

            // Get metric values for this row
            const metricVals = {{}};
            metrics.header.forEach((h, i) => {{
                if (metricCols.includes(h)) {{
                    metricVals[h] = row[i];
                }}
            }});

            if (!dataByToolMuts.has(key)) {{
                dataByToolMuts.set(key, []);
            }}
            dataByToolMuts.get(key).push({{
                sample,
                tool,
                numMuts,
                features: featureVals,
                metrics: metricVals
            }});
        }});

        function updatePlot() {{
            const feat = featureSelect.value;
            const metric = metricSelect.value;

            if (!feat || !metric) {{
                document.getElementById('{widget_id}_plot').innerHTML =
                    '<div style="display:flex; align-items:center; justify-content:center; height:100%; color:#9ca3af; font-size:18px;">📊 Select a feature and metric to visualize</div>';
                return;
            }}

            const selectedMuts = mutsSelect.value;

            // Collect all data points with feature and metric values
            const allData = [];

            dataByToolMuts.forEach((dataPoints, key) => {{
                dataPoints.forEach(point => {{
                    const numMuts = point.numMuts;

                    // Filter by selected mutation burden
                    if (selectedMuts && numMuts !== selectedMuts) return;

                    // Get feature and metric values
                    const featureVal = parseFloat(point.features[feat]);
                    const metricVal = parseFloat(point.metrics[metric]);

                    // Skip if invalid
                    if (isNaN(featureVal) || isNaN(metricVal)) return;

                    allData.push({{
                        tool: point.tool,
                        feature: featureVal,
                        metric: metricVal,
                        sample: point.sample,
                        numMuts: numMuts
                    }});
                }});
            }});

            if (allData.length === 0) {{
                document.getElementById('{widget_id}_plot').innerHTML =
                    '<div style="display:flex; align-items:center; justify-content:center; height:100%; color:#ef4444; font-size:16px;">❌ No data available for selected options</div>';
                return;
            }}

            // Calculate median threshold for binning (proper median calculation)
            const featureValues = allData.map(d => d.feature);
            featureValues.sort((a,b) => a - b);
            let medianThreshold;
            if (featureValues.length % 2 === 0) {{
                // Even number: average of two middle values
                medianThreshold = (featureValues[featureValues.length/2 - 1] + featureValues[featureValues.length/2]) / 2;
            }} else {{
                // Odd number: middle value
                medianThreshold = featureValues[Math.floor(featureValues.length / 2)];
            }}

            // Bin data into low/high by feature value (using < for high to handle ties at median)
            const dataLow = allData.filter(d => d.feature < medianThreshold);
            const dataHigh = allData.filter(d => d.feature >= medianThreshold);

            // If all values are the same (or very skewed), split evenly to avoid empty groups
            if (dataLow.length === 0 || dataHigh.length === 0) {{
                allData.sort((a, b) => a.feature - b.feature);
                const midPoint = Math.floor(allData.length / 2);
                dataLow.length = 0;
                dataHigh.length = 0;
                dataLow.push(...allData.slice(0, midPoint));
                dataHigh.push(...allData.slice(midPoint));
            }}

            // Group by tool for each bin
            const tools = [...new Set(allData.map(d => d.tool))].sort();

            // Tool colors - embedded from ReportMetadata.tool_colors (see get_tool_colors())
            const toolColors = __TOOL_COLORS__;

            // Create boxplot traces - separate trace for each tool, grouped by Low/High
            const traces = [];

            // For each tool, create a Low trace and a High trace
            tools.forEach(tool => {{
                const toolDataLow = dataLow.filter(d => d.tool === tool);
                const toolDataHigh = dataHigh.filter(d => d.tool === tool);

                // Low trace for this tool
                if (toolDataLow.length > 0) {{
                    traces.push({{
                        y: toolDataLow.map(d => d.metric),
                        x: Array(toolDataLow.length).fill(tool),
                        name: 'Low',
                        type: 'box',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: 0,
                        marker: {{
                            size: 4,
                            opacity: 0.5,
                            color: toolColors[tool] || '#6b7280'
                        }},
                        line: {{
                            color: toolColors[tool] || '#6b7280',
                            width: 2
                        }},
                        fillcolor: `rgba${{(toolColors[tool] || '#6b7280').replace('#','').match(/.{{2}}/g).map(x=>parseInt(x,16)).join(',')}},.3)`,
                        text: toolDataLow.map(d => `Sample: ${{d.sample}}<br>${{feat}}: ${{d.feature.toFixed(3)}}<br>${{metric}}: ${{d.metric.toFixed(3)}}`),
                        hovertemplate: `<b>${{tool}} - Low</b><br>%{{text}}<extra></extra>`,
                        offsetgroup: 'Low',
                        legendgroup: 'Low',
                        showlegend: tool === tools[0]  // Only show legend for first tool
                    }});
                }}

                // High trace for this tool
                if (toolDataHigh.length > 0) {{
                    traces.push({{
                        y: toolDataHigh.map(d => d.metric),
                        x: Array(toolDataHigh.length).fill(tool),
                        name: 'High',
                        type: 'box',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: 0,
                        marker: {{
                            size: 4,
                            opacity: 0.6,
                            color: toolColors[tool] || '#6b7280'
                        }},
                        line: {{
                            color: toolColors[tool] || '#6b7280',
                            width: 2
                        }},
                        fillcolor: `rgba${{(toolColors[tool] || '#6b7280').replace('#','').match(/.{{2}}/g).map(x=>parseInt(x,16)).join(',')}},.5)`,
                        text: toolDataHigh.map(d => `Sample: ${{d.sample}}<br>${{feat}}: ${{d.feature.toFixed(3)}}<br>${{metric}}: ${{d.metric.toFixed(3)}}`),
                        hovertemplate: `<b>${{tool}} - High</b><br>%{{text}}<extra></extra>`,
                        offsetgroup: 'High',
                        legendgroup: 'High',
                        showlegend: tool === tools[0]  // Only show legend for first tool
                    }});
                }}
            }});

            // Get theoretical minimum TAE (if available and metric is TAE)
            let theoreticalMinTAE = null;
            if (metric === 'TAE' && allData.length > 0 && allData[0].metrics && 'theoretical_min_TAE' in allData[0].metrics) {{
                theoreticalMinTAE = allData[0].metrics['theoretical_min_TAE'];
            }}

            // Add theoretical maximum line if applicable
            if (theoreticalMinTAE !== null && theoreticalMinTAE > 0) {{
                traces.push({{
                    x: tools,
                    y: Array(tools.length).fill(theoreticalMinTAE),
                    name: 'Theoretical Best (Unknown Sigs)',
                    type: 'scatter',
                    mode: 'lines',
                    line: {{
                        color: '#ef4444',
                        width: 2,
                        dash: 'dash'
                    }},
                    showlegend: true,
                    hovertemplate: `Theoretical Minimum TAE: ${{theoreticalMinTAE.toFixed(4)}}<br>(Due to unknown signatures)<extra></extra>`
                }});
            }}

            const layout = {{
                title: {{
                    text: `<b>${{metric}}</b> distribution by <b>${{feat}}</b> (Low vs High)`,
                    font: {{ size: 16, color: '#111827' }}
                }},
                xaxis: {{
                    title: {{ text: 'Tool', font: {{ size: 13, color: '#374151' }} }},
                    gridcolor: '#e5e7eb',
                    categoryorder: 'array',
                    categoryarray: tools
                }},
                yaxis: {{
                    title: {{ text: metric, font: {{ size: 14, color: '#374151' }} }},
                    gridcolor: '#e5e7eb',
                    zeroline: true,
                    zerolinecolor: '#9ca3af',
                    zerolinewidth: 1
                }},
                showlegend: true,
                legend: {{
                    x: 1.02,
                    y: 1,
                    bgcolor: 'rgba(255,255,255,0.9)',
                    bordercolor: '#e5e7eb',
                    borderwidth: 1,
                    title: {{ text: 'Group', font: {{ size: 12, color: '#374151' }} }}
                }},
                plot_bgcolor: 'white',
                paper_bgcolor: '#fafafa',
                margin: {{ l: 70, r: 120, t: 70, b: 90 }},
                boxmode: 'group',
                boxgap: 0.2,
                boxgroupgap: 0.1
            }};

            Plotly.newPlot('{widget_id}_plot', traces, layout, {{responsive:true, displayModeBar: true}});

            // Update stats
            document.getElementById('{widget_id}_stats').style.display = 'block';
            document.getElementById('{widget_id}_threshold').textContent = medianThreshold.toFixed(3);
            document.getElementById('{widget_id}_n_low').textContent = dataLow.length;
            document.getElementById('{widget_id}_n_high').textContent = dataHigh.length;
        }}

        document.getElementById('{widget_id}_update').addEventListener('click', updatePlot);

        // Pre-select first options and auto-plot
        if (featureSelect.options.length > 1) featureSelect.selectedIndex = 1;
        if (metricSelect.options.length > 1) metricSelect.selectedIndex = 1;
        if (mutsSelect.options.length) mutsSelect.selectedIndex = 0;

        // Auto-update on load
        setTimeout(updatePlot, 100);
}})();
</script>
"""
    # Prepare data for embedding (avoid CORS issues with file:// protocol)
    import json

    if data_json is None:
        if data_df is not None:
            # Use provided dataframe
            df = data_df
        else:
            # Fallback to reading file (but this won't work with file:// protocol)
            df = pd.read_csv(per_sample_metrics_tsv, sep='\t')

        # Convert dataframe to {{header: [...], rows: [...]}} format for JavaScript
        data_dict = {
            'header': df.columns.tolist(),
            'rows': df.astype(str).values.tolist()  # Convert all to strings for JSON
        }
        data_json = json.dumps(data_dict)

    # Embed data and tool colors directly in JavaScript to avoid fetch/CORS issues
    html = html.replace('__DATA__', data_json)
    html = html.replace('__TOOL_COLORS__', json.dumps(colors or {}))

    return html


def build_signature_activity_data(results_dir, df, tools, colors, out_of_reference_sigs=None,
                                   relevance='standard', selected_rep=None, min_contribution=None):
    """Load per-sample ground-truth vs. predicted signature activities for the
    Signature Activity widget (see create_signature_activity_widget below).

    Reads two things straight off disk (these aren't in any aggregated dataframe
    already loaded by the report). Both "synthetic" and "custom" data types have
    ground truth and use the same directory layout (just "custom" uses the dataset
    name where "synthetic" uses a cancer type), so both are included here:
      - ground truth: results_dir/generated_data/{data_type}/weights_{cancer}-w_{rep}-{burden}.dat
        (signature fractions per sample; multiplied by burden to get mutation counts)
      - predictions: results_dir/{data_type}/{cancer}/{relevance}/{burden}/rep_{rep}/{tool}_contribution.dat

    Predicted signatures below `min_contribution` mutations are dropped (not just those
    at exactly 0): the evaluation metrics (P/R/F1/MCC/FP/FN, via evaluate_tool_performance.py
    --min-contribution) already treat them as "not called," so showing them as markers here
    would contradict the metrics computed from the same threshold. Ground truth is left at
    the raw >0 cut — it's what was actually simulated, not a detection call.

    Returns a JSON string ready to embed client-side, or None if no ground-truth data
    (synthetic or custom) is available for this run (e.g. real-sample-only reports).
    """
    min_contribution = DEFAULT_MIN_CONTRIBUTION if min_contribution is None else min_contribution
    if not results_dir or df is None or df.empty:
        return None

    gt_df = df[df['data_type'].isin(['synthetic', 'custom'])] if 'data_type' in df.columns else df
    if gt_df.empty or 'cancer_type' not in gt_df.columns or 'num_muts' not in gt_df.columns:
        return None

    has_data_type_col = 'data_type' in gt_df.columns
    cancer_data_types = sorted(gt_df[['cancer_type', 'data_type']].drop_duplicates().itertuples(index=False)) \
        if has_data_type_col else [(c, 'synthetic') for c in sorted(gt_df['cancer_type'].unique())]
    burdens = sorted(int(b) for b in gt_df['num_muts'].unique())
    if selected_rep is None:
        selected_rep = int(gt_df['replicate'].min()) if 'replicate' in gt_df.columns else 0

    unknown_sigs = set(out_of_reference_sigs or [])
    results_root = Path(results_dir)
    sorted_tools = sorted(tools)

    marker_options = ['circle', 'square', 'diamond', 'cross', 'x', 'triangle-up',
                       'triangle-down', 'pentagon', 'hexagon', 'star']
    tool_markers = {t: marker_options[i % len(marker_options)] for i, t in enumerate(sorted_tools)}

    samples = None
    truth_data = {}
    pred_data = {}

    for cancer_type, data_type in cancer_data_types:
        for burden in burdens:
            weights_path = (results_root / 'generated_data' / data_type /
                             f'weights_{cancer_type}-w_{selected_rep}-{burden}.dat')
            if not weights_path.exists():
                continue
            truth_df = pd.read_csv(weights_path, sep='\t', index_col=0)
            if samples is None:
                samples = sorted(truth_df.columns.tolist())

            for sample in truth_df.columns:
                counts = truth_df[sample] * burden
                active = counts[counts > 0]
                if not active.empty:
                    truth_data[f"{cancer_type}|{burden}|{sample}"] = {
                        sig: round(float(v), 2) for sig, v in active.items()
                    }

            base_dir = results_root / data_type / cancer_type / relevance / str(burden) / f'rep_{selected_rep}'
            for tool in sorted_tools:
                contrib_path = base_dir / f'{tool}_contribution.dat'
                if not contrib_path.exists():
                    continue
                pred_df = pd.read_csv(contrib_path, index_col=0)
                for sample in pred_df.columns:
                    active = pred_df[sample][pred_df[sample] >= min_contribution]
                    if active.empty:
                        continue
                    pred_data[f"{cancer_type}|{burden}|{sample}|{tool}"] = {
                        sig: round(float(v), 2) for sig, v in active.items()
                    }

    if not truth_data:
        return None

    return json.dumps({
        'cancer_types': sorted({c for c, _ in cancer_data_types}),
        'burdens': burdens,
        'samples': samples or [],
        'tools': sorted_tools,
        'tool_colors': colors or {},
        'tool_markers': tool_markers,
        'unknown_sigs': sorted(unknown_sigs),
        'selected_rep': selected_rep,
        'truth': truth_data,
        'predictions': pred_data,
    })


def create_signature_activity_widget(data_json, cancer_type, widget_id='sig_activity'):
    """Return HTML + JS for the interactive ground-truth-vs-predicted signature
    activity plot: pick a mutation burden / sample (cancer type is fixed to the
    enclosing tab) and see every tool's predicted signature activities plotted
    against the ground truth for that one sample. Analogous to the "Signature
    Activity Comparison" plot in the SigFitTest_workflow penalty analysis report,
    minus the penalty selector (SigFitBench tools each produce a single
    contribution file, not a penalty sweep) and minus the cancer-type selector
    (one widget instance is embedded per cancer-type tab instead).
    """
    html = f"""
<div style="background: white; padding: 1.5rem; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); margin: 1rem 0;">
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin-bottom: 1rem;">
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">Mutation Burden</label>
            <select id="{widget_id}_burden" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;"></select>
        </div>
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">Sample</label>
            <select id="{widget_id}_sample" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;"></select>
        </div>
    </div>
    <div id="{widget_id}_plot" style="width: 100%; height: 650px; border: 2px solid #f3f4f6; border-radius: 8px; background: #fafafa;"></div>
    <p style="color: #6b7280; font-size: 13px; margin-top: 0.75rem; font-style: italic; text-align: center;">
        Horizontal lines: ground-truth activity (red = reference signature, green = out-of-reference/unknown).
        Markers: each tool's predicted activity for that signature. Y-axis is log-scale mutation counts.
        Click a tool in the legend to show/hide it, same as the other plots in this report &mdash;
        the x-axis updates to drop any signature only that tool predicted.
    </p>
</div>

<script>
(function() {{
    const widgetId = '{widget_id}';
    const data = __SIG_ACTIVITY_DATA__;
    const cancer = {cancer_type!r};
    const plotDiv = widgetId + '_plot';

    if (!data || !data.cancer_types || !data.cancer_types.includes(cancer)) {{
        document.getElementById(plotDiv).innerHTML =
            '<div style="display:flex; align-items:center; justify-content:center; height:100%; color:#9ca3af; font-size:16px;">No per-sample ground-truth data available for this cancer type</div>';
        return;
    }}

    const burdenSelect = document.getElementById(widgetId + '_burden');
    data.burdens.forEach(b => {{
        const opt = document.createElement('option');
        opt.value = b; opt.text = b.toLocaleString();
        burdenSelect.appendChild(opt);
    }});

    const sampleSelect = document.getElementById(widgetId + '_sample');
    data.samples.forEach(s => {{
        const opt = document.createElement('option');
        opt.value = s; opt.text = s;
        sampleSelect.appendChild(opt);
    }});

    const unknownSigs = new Set(data.unknown_sigs || []);
    const tools = data.tools || [];

    // Which tools are currently shown. Starts as "all of them"; updated by reading Plotly's
    // own trace.visible state after a legend click (see attachRestyleListener below) — the
    // legend itself behaves exactly like every other plot's legend (plain Plotly default
    // click-to-toggle, nothing intercepted or prevented).
    let hiddenTools = new Set();
    let syncingAxis = false;  // true while we're applying our own axis-rebuild react() call

    function horizontalLineTrace(xPos, truthRaw, sigs, color, name) {{
        const x = [], y = [], text = [];
        sigs.forEach(sig => {{
            for (let k = 0; k < 20; k++) {{
                x.push(xPos[sig] - 0.35 + k * 0.035);
                y.push(truthRaw[sig]);
                text.push(`${{sig}}: ${{truthRaw[sig].toFixed(1)}}`);
            }}
            x.push(null); y.push(null); text.push(null);
        }});
        if (!x.length) return null;
        return {{
            x, y, mode: 'lines', type: 'scatter',
            line: {{ color, width: 3 }}, name, text,
            hovertemplate: '%{{text}}<extra></extra>',
        }};
    }}

    function buildFigure() {{
        const burden = burdenSelect.value;
        const sample = sampleSelect.value;
        const visibleTools = tools.filter(t => !hiddenTools.has(t));

        const truthRaw = data.truth[`${{cancer}}|${{burden}}|${{sample}}`] || {{}};
        const truthSigs = new Set(Object.keys(truthRaw));

        // Load every tool's predictions (not just visible ones) so a re-shown tool's trace has
        // real data immediately; only *visible* tools' predictions decide which signature
        // columns exist on the axis.
        const toolPreds = {{}};
        const predOnlySigsSet = new Set();
        tools.forEach(tool => {{
            const preds = data.predictions[`${{cancer}}|${{burden}}|${{sample}}|${{tool}}`] || {{}};
            toolPreds[tool] = preds;
            if (!hiddenTools.has(tool)) {{
                Object.keys(preds).forEach(sig => {{ if (!truthSigs.has(sig)) predOnlySigsSet.add(sig); }});
            }}
        }});

        const truthKnown = [...truthSigs].filter(s => !unknownSigs.has(s)).sort();
        const truthUnknown = [...truthSigs].filter(s => unknownSigs.has(s)).sort();
        const predOnly = [...predOnlySigsSet].sort();
        const orderedSigs = [...truthKnown, ...truthUnknown, ...predOnly];

        const titleBase = `${{cancer}} | ${{Number(burden).toLocaleString()}} mutations | ${{sample}}`;

        if (orderedSigs.length === 0) {{
            return {{
                traces: [],
                layout: {{
                    title: {{ text: `No data: ${{titleBase}}` }},
                    annotations: [{{
                        text: 'No Data Available<br>This combination may not have been run or all signatures have zero activity.',
                        xref: 'paper', yref: 'paper', x: 0.5, y: 0.5, xanchor: 'center', yanchor: 'middle',
                        showarrow: false, font: {{ size: 16, color: '#666' }}
                    }}],
                    xaxis: {{ visible: false }}, yaxis: {{ visible: false }},
                    height: 650, margin: {{ t: 80, b: 60, l: 60, r: 60 }}
                }}
            }};
        }}

        const xPos = {{}};
        orderedSigs.forEach((sig, i) => {{ xPos[sig] = i; }});

        // Fixed offset basis (over the full tool list, not just visible ones) so a tool's
        // marker position doesn't jump around every time some other tool is toggled.
        const n = tools.length;
        const toolOffset = {{}};
        tools.forEach((t, i) => {{ toolOffset[t] = n > 1 ? -0.15 + (i / (n - 1)) * 0.3 : 0; }});

        const traces = [];
        const truthKnownTrace = horizontalLineTrace(xPos, truthRaw, truthKnown, '#C1121F', 'Truth (Reference)');
        if (truthKnownTrace) traces.push(truthKnownTrace);
        const truthUnknownTrace = horizontalLineTrace(xPos, truthRaw, truthUnknown, '#2DA44E', 'Truth (Unknown)');
        if (truthUnknownTrace) traces.push(truthUnknownTrace);

        // Every tool gets a trace with its real data, whether currently shown or not — a
        // trace forced to completely empty x/y combined with visible:'legendonly' can make
        // Plotly drop the legend entry outright instead of just dimming it, which is what
        // made hidden tools un-clickable before. Hiding is done purely via `visible` below.
        tools.forEach(tool => {{
            const color = (data.tool_colors || {{}})[tool] || '#6b7280';
            const sym = (data.tool_markers || {{}})[tool] || 'circle';
            const marker = {{ symbol: sym, size: 9, color, line: {{ color: 'black', width: 1 }}, opacity: 0.9 }};
            const preds = toolPreds[tool];
            const px = [], py = [], ptext = [];
            orderedSigs.forEach(sig => {{
                if (!preds[sig] || preds[sig] <= 0) return;
                px.push(xPos[sig] + toolOffset[tool]);
                py.push(preds[sig]);
                ptext.push(`${{tool}}: ${{sig}} = ${{preds[sig].toFixed(1)}}`);
            }});
            traces.push({{
                x: px, y: py, mode: 'markers', type: 'scatter', marker,
                name: tool, text: ptext,
                hovertemplate: '%{{text}}<extra></extra>',
                visible: hiddenTools.has(tool) ? 'legendonly' : true,
            }});
        }});

        return {{
            traces,
            layout: {{
                title: {{ text: `${{titleBase}}<br><sub>Rep ${{data.selected_rep}} | Line = ground truth, markers = predicted</sub>` }},
                xaxis: {{
                    tickmode: 'array',
                    tickvals: orderedSigs.map((_, i) => i),
                    ticktext: orderedSigs,
                    tickangle: -45,
                    tickfont: {{ size: 10 }},
                    title: 'Signature',
                }},
                yaxis: {{ type: 'log', title: 'Mutations' }},
                legend: {{ orientation: 'v', yanchor: 'top', y: 1, xanchor: 'left', x: 1.02 }},
                height: 650,
                margin: {{ t: 80, b: 100, l: 60, r: 160 }},
                hovermode: 'closest',
            }}
        }};
    }}

    function attachRestyleListener(gd) {{
        if (gd.__sigActivityRestyleBound) return;
        gd.__sigActivityRestyleBound = true;
        // Fires after Plotly's own (un-intercepted) legend click has already applied its
        // default visibility toggle. We read that resulting state back and rebuild the axis
        // to drop/restore signatures the now-hidden/shown tool exclusively predicts. Using
        // Plotly.react (not newPlot) here, after the click has fully finished, avoids the
        // teardown/race that made legend entries disappear when the rebuild ran synchronously
        // inside the click handler itself.
        gd.on('plotly_restyle', () => {{
            if (syncingAxis) return;  // ignore restyle events caused by our own render() below
            hiddenTools = new Set(tools.filter(tool => {{
                const trace = (gd.data || []).find(t => t.name === tool);
                return trace && trace.visible === 'legendonly';
            }}));
            render();
        }});
    }}

    function render() {{
        const {{ traces, layout }} = buildFigure();
        syncingAxis = true;
        Plotly.react(plotDiv, traces, layout, {{responsive: true}}).then(gd => {{
            syncingAxis = false;
            attachRestyleListener(gd);
        }});
    }}

    burdenSelect.addEventListener('change', render);
    sampleSelect.addEventListener('change', render);

    render();
}})();
</script>
"""
    return html.replace('__SIG_ACTIVITY_DATA__', data_json)


def create_visualizations_section(df, df_per_sample=None, df_bootstrap=None, df_raw=None,
                                   results_dir=None, id_prefix='', include_plotly_cdn=True,
                                   out_of_reference_sigs=None, relevance='standard',
                                   unknown_sig_config=None, min_contribution=None):
    """Create main visualizations section with tabs.

    Args:
        df: DataFrame with aggregated metrics
        df_per_sample: Optional DataFrame with per-sample metrics for violin plots
        df_bootstrap: Optional DataFrame with bootstrap quality metrics
        df_raw: Optional DataFrame with individual realization data (before aggregation)
        results_dir: Optional base results directory for reading contribution files
        id_prefix: String prefix for all HTML element IDs (needed when called multiple times)
        include_plotly_cdn: Whether to include plotly.js CDN link in the first figure
        unknown_sig_config: Optional dict {'num_unknowns': int, 'weights': list} with the
            simulated unknown-signature configuration (currently unused here; kept for a
            uniform call signature with create_report_intro/create_appendix_section).
    """

    # Get tool order and filter to only tools present in data
    tool_order = get_tool_display_order()
    available_tools = set(df['tool'].unique())
    tools = [t for t in tool_order if t in available_tools]
    # Add any tools not in our predefined order (sorted alphabetically at the end)
    remaining_tools = sorted(available_tools - set(tools))
    tools.extend(remaining_tools)

    cancer_types = sorted(df['cancer_type'].unique())
    colors = get_tool_colors(tools)

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

    # Precompute the signature-activity widget's data once (covers all cancer types);
    # each per-cancer-type tab below embeds its own widget instance scoped to that cancer type.
    sig_activity_json = build_signature_activity_data(
        results_dir, df, tools, colors, out_of_reference_sigs, relevance, min_contribution=min_contribution
    )

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
        
        # SECTION 3.5: Sample-by-Sample Performance Comparison — the per-sample jittered
        # scatter (create_per_method_sample_plots) becomes unreadable past ~10 samples per
        # burden (overlapping x-axis labels, illegible jitter clusters), so cancer types
        # with more samples get a ranked "hardest samples" summary instead.
        if df_per_sample is not None:
            cancer_sample_count = (
                df_per_sample.loc[df_per_sample['cancer_type'] == cancer_type, 'sample'].nunique()
                if 'sample' in df_per_sample.columns else 0
            )
            if cancer_sample_count > 10:
                fig_per_method = create_hard_samples_row_plot(
                    df_per_sample, cancer_type, results_dir, out_of_reference_sigs)
                section_title = 'Hardest Samples — which samples are hard for all tools?'
                section_desc = (
                    f'With {cancer_sample_count} samples per condition, an individual per-sample plot is too '
                    'crowded to read, so samples are ranked by mean TAE across all tools instead (higher TAE '
                    '= more of the true signal misattributed, i.e. harder for everyone). Left to right: the '
                    'hardest samples; how often each ground-truth signature occurs in the hardest samples '
                    '(orange) vs. all samples (grey), unknown signatures outlined in red; effect of each '
                    'signature\'s presence on mean TAE (red = unknown, blue = known, * = significant q &lt; 0.05); '
                    'and, where available, difficulty vs. ENP (spectral flatness) of the sample\'s most diffuse '
                    'known signature, hardest samples in red. Use the dropdown to switch mutation burden.'
                )
                div_id = f'hard_samples_{tab_id}'
            else:
                fig_per_method = create_per_method_sample_plots(df_per_sample, cancer_type, tools, colors)
                section_title = 'Sample-by-Sample Comparison — which samples are hard for all tools?'
                section_desc = (
                    'Each column is one simulated sample. Grey boxplot = spread of all tools on that '
                    'sample; colored markers = individual tool results. Samples where the boxplot is low '
                    'or tools spread far apart are the hardest cases. Use the dropdown to switch mutation burden.'
                )
                div_id = f'per_method_{tab_id}'

            if fig_per_method is not None:
                html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #059669; border-radius: 4px;">\n'
                html += f'<h4 style="margin-top: 0; color: #059669;">{section_title}</h4>\n'
                html += f'<p style="color: #666; margin-bottom: 1rem;">{section_desc}</p>\n'

                # Render the figure HTML with built-in dropdown (no custom controls needed)
                fig_html = fig_per_method.to_html(include_plotlyjs=False, full_html=False, div_id=div_id, config={'responsive': True})

                # Simple layout without custom controls (dropdown is built into the plot)
                html += '<div class="plot-container">\n'
                html += fig_html
                html += '</div>\n'
                html += '</div>\n\n'
        
        # SECTION 5: Signature Activity per Sample (ground truth vs. predicted)
        if sig_activity_json is not None:
            html += '<div style="margin: 2rem 0; padding: 1.5rem; background: #f8f9fa; border-left: 4px solid #0891b2; border-radius: 4px;">\n'
            html += f'<h4 style="margin-top: 0; color: #0891b2;" id="{tab_id}-signature-activity">🧬 Signature Activity per Sample</h4>\n'
            html += (
                '<p style="color:#666; margin-bottom: 1rem;">'
                'Ground-truth vs. predicted signature activities for one sample at a time '
                '&mdash; pick a mutation burden and sample below.</p>\n'
            )
            html += create_signature_activity_widget(
                sig_activity_json, cancer_type, widget_id=f'sig_activity_{tab_id}'
            )
            html += '</div>\n'  # Close signature activity section

        html += '</div>\n'
    
    # Aggregate tab
    html += f'<div id="{id_prefix}aggregate-tab" class="tab-content">\n'
    html += '<h3>Aggregate Analysis</h3>\n'
    html += '<p>Performance averaged across all cancer types.</p>\n'
    
    # Precision-Recall plot
    fig_agg = create_aggregate_plot(df, tools, colors, df_per_sample)
    html += '<div class="plot-container">\n'
    html += fig_agg.to_html(include_plotlyjs=False, div_id=f'{id_prefix}aggregate_plot', config={'responsive': True})
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
    
    fig_agg_metrics = create_metrics_grid(agg_metrics_df, 'All Cancer Types (Pooled)', tools, colors, df_per_sample, pooled=True)
    html += '<div class="plot-container">\n'
    html += fig_agg_metrics.to_html(include_plotlyjs=False, div_id=f'{id_prefix}aggregate_metrics_grid', config={'responsive': True})
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
    
    html += '</section>\n'
    
    return html


# ============================================================================
# Main Report Generator
# ============================================================================

def find_failed_tool_runs(results_dir, data_type):
    """Tool runs that recorded a failure in their result file ("# TOOL FAILED: <reason>").

    Looks in results_dir/{data_type}/{dataset}/{relevance}/{num_muts}/rep_{rep}/ (the RUN_DIR
    layout in workflow/rules/common.smk).
    """
    failed = []
    if not results_dir or not data_type:
        return failed
    for path in sorted((Path(results_dir) / data_type).glob("*/*/*/rep_*/*_contribution.dat")):
        with open(path) as f:
            first_line = f.readline().strip()
        if first_line.startswith("# TOOL FAILED"):
            dataset, relevance, num_muts, rep = path.parts[-5:-1]
            failed.append({
                "tool": path.name[:-len("_contribution.dat")],
                "run": f"{dataset}, {num_muts} mutations, {relevance}, {rep.replace('_', ' ')}",
                "reason": first_line.split(":", 1)[-1].strip() or "see the tool's log",
            })
    return failed


def create_failed_tools_section(failed):
    """Section listing the tools that failed, in which runs, and why. Empty if none failed."""
    if not failed:
        return ""
    by_tool = {}
    for row in failed:
        by_tool.setdefault(row["tool"], []).append(row)
    out = '<section id="failed-tools" class="section">\n'
    out += '<h2>\u26a0\ufe0f Failed tools</h2>\n'
    out += ('<p>These tools did not produce a result in the runs listed, so they are left out of '
            'the tables and plots for those runs. The full error is in each run\'s tool log.</p>\n')
    out += ('<table style="width:100%;border-collapse:collapse;font-size:14px;">\n'
            '<tr style="text-align:left;border-bottom:2px solid #d1d5db;">'
            '<th style="padding:6px;">Tool</th><th style="padding:6px;">Failed runs</th>'
            '<th style="padding:6px;">Reason</th><th style="padding:6px;">Runs</th></tr>\n')
    for tool, rows in sorted(by_tool.items()):
        reasons = sorted({r["reason"] for r in rows})
        runs = [r["run"] for r in rows]
        shown = ", ".join(runs[:6]) + (f" and {len(runs) - 6} more" if len(runs) > 6 else "")
        out += ('<tr style="border-bottom:1px solid #e5e7eb;vertical-align:top;">'
                f'<td style="padding:6px;font-weight:600;">{html.escape(tool)}</td>'
                f'<td style="padding:6px;">{len(rows)}</td>'
                f'<td style="padding:6px;">{"<br>".join(html.escape(r) for r in reasons)}</td>'
                f'<td style="padding:6px;color:#6b7280;">{html.escape(shown)}</td></tr>\n')
    out += '</table>\n</section>\n'
    return out


def generate_report(df, output_file, title, subtitle, df_per_sample=None, df_bootstrap=None, df_raw=None, results_dir=None, unknown_sig_config=None, cosmic_version_label=None, out_of_reference_sigs=None, relevance='standard', tool_overrides=None, df_benchmarks=None, min_contribution=None, failed_runs=None):
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
        unknown_sig_config: Optional dict with unknown signature configuration {'num_unknowns': int, 'weights': list}
        cosmic_version_label: Optional human-readable COSMIC catalog label derived from
                          config's cosmic_signatures.base_file
        out_of_reference_sigs: Optional list of signature names simulated but excluded from the
                          fitting catalog, used to label truth signatures as reference vs. unknown
                          in the Signature Activity widget
        relevance: Relevance mode subdirectory name used under results_dir/synthetic/{cancer}/{relevance}/...
        df_benchmarks: Optional pandas.DataFrame with per-tool wall time / peak memory
                          (all_tool_benchmarks.tsv — see aggregate_tool_benchmarks.py)
        min_contribution: Optional float — config's evaluation.min_contribution (the
                          mutation-count threshold for a signature to count as "active";
                          see create_appendix_section / build_signature_activity_data)
    """
    print("\nGenerating report sections...")

    # Generate all sections
    intro_html = create_report_intro(df, unknown_sig_config, tool_overrides)
    tool_order = get_tool_display_order()
    available_tools = set(df['tool'].unique())
    all_tools = [t for t in tool_order if t in available_tools]
    all_tools += sorted(available_tools - set(all_tools))
    resources_html = create_resource_usage_section(df_benchmarks, get_tool_colors(all_tools))
    appendix_html = create_appendix_section(df, df_raw, unknown_sig_config, cosmic_version_label, min_contribution)

    # Build visualizations section — one sub-section per data_type when multiple exist
    data_types = sorted(df['data_type'].unique()) if 'data_type' in df.columns else ['synthetic']
    if len(data_types) <= 1:
        visualizations_html = create_visualizations_section(
            df, df_per_sample, df_bootstrap, df_raw, results_dir,
            out_of_reference_sigs=out_of_reference_sigs, relevance=relevance,
            unknown_sig_config=unknown_sig_config, min_contribution=min_contribution
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
                out_of_reference_sigs=out_of_reference_sigs,
                relevance=relevance,
                unknown_sig_config=unknown_sig_config,
                min_contribution=min_contribution,
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

    content = (intro_html + create_failed_tools_section(failed_runs) + visualizations_html
               + resources_html + appendix_html)
    
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

def keep_data_type(df, data_type):
    """Keep only the rows of one data type (all rows if data_type is None)."""
    if df is None or data_type is None or 'data_type' not in df.columns:
        return df
    return df[df['data_type'] == data_type].reset_index(drop=True)


def main(summaries, per_sample, bootstrap_quality, output, title, subtitle,
         results_dir, config_file, tool_benchmarks=None, data_type=None):
    """Main entry point. With data_type, the report shows only that data type."""
    # Load unknown signature info from config if provided
    unknown_sig_config = None
    tool_overrides = None
    cosmic_version_label = None
    out_of_reference_sigs = []
    relevance = 'standard'
    min_contribution = None

    if config_file:
        print(f"Attempting to load config file: {config_file}")
        try:
            if not os.path.exists(config_file):
                print(f"ERROR: Config file does not exist: {config_file}")
            else:
                print(f"Config file exists, loading...")
            with open(config_file, 'r') as f:
                config_data = yaml.safe_load(f)
                cosmic_version_label = _derive_cosmic_label(
                    config_data.get('cosmic_signatures', {}).get('base_file')
                )
                out_of_reference_sigs = config_data.get('out_of_reference_sigs', []) or []
                relevance = (config_data.get('relevance_modes') or ['standard'])[0]
                tool_overrides = config_data.get('tool_overrides') or {}
                min_contribution = config_data.get('evaluation', {}).get('min_contribution')
                out_of_ref_weights = config_data.get('out_of_reference_weights', [])
                if out_of_ref_weights:
                    unknown_sig_config = {
                        'num_unknowns': len(out_of_ref_weights),
                        'weights': out_of_ref_weights
                    }
                    print(f"  - With fixed weights: {out_of_ref_weights}")
                else:
                    print("\nNo unknown signatures configured (out_of_reference_sigs is empty)")
                if data_type == 'custom':
                    # Custom datasets have their own catalogue and hidden signatures
                    datasets = (config_data.get('custom_datasets') or {}).values()
                    catalogues = {ds.get('signature_file') or
                                  config_data.get('cosmic_signatures', {}).get('base_file')
                                  for ds in datasets}
                    if len(catalogues) == 1:
                        cosmic_version_label = _derive_cosmic_label(catalogues.pop())
                    out_of_reference_sigs = sorted({s for ds in datasets
                                                    for s in ds.get('out_of_reference_sigs') or []})
                    unknown_sig_config = None
        except Exception as e:
            print(f"\n❌ ERROR: Could not load config file {config_file}: {e}")
    else:
        print("No config file provided")

    # Infer results_dir from output path if not provided
    if not results_dir and output:
        results_dir = str(Path(output).parent)

    # Load summary data only for df_raw (used in realization stability plots)
    _, df_raw = load_evaluation_data(summaries, return_raw=True)
    df_raw = keep_data_type(df_raw, data_type)

    df_per_sample = None
    if per_sample:
        print(f"\nLoading per-sample metrics from {len(per_sample)} file(s)...")
        df_per_sample = keep_data_type(load_per_sample_data(per_sample), data_type)

    if df_per_sample is not None and not df_per_sample.empty:
        df = reaggregate_from_per_sample(df_per_sample)
        print(f"  Summary df re-aggregated from per-sample data. Shape: {df.shape}")
    else:
        df, _ = load_evaluation_data(summaries, return_raw=True)
        df = keep_data_type(df, data_type)

    # Load bootstrap quality assessment if provided
    df_bootstrap = None
    if bootstrap_quality:
        print(f"\nLoading bootstrap quality assessment from {len(bootstrap_quality)} file(s)...")
        df_bootstrap = load_evaluation_data(bootstrap_quality)

    # Load per-tool wall time / peak memory if provided (all_tool_benchmarks.tsv, written
    # by aggregate_tool_benchmarks.py from Snakemake's benchmark: TSVs — see tools.smk)
    df_benchmarks = None
    if tool_benchmarks:
        print(f"\nLoading tool resource benchmarks from {len(tool_benchmarks)} file(s)...")
        df_benchmarks = keep_data_type(load_evaluation_data(tool_benchmarks), data_type)

    failed_runs = find_failed_tool_runs(results_dir, data_type)
    if failed_runs:
        print(f"\n{len(failed_runs)} failed tool run(s) will be listed in the report")

    generate_report(df, output, title, subtitle, df_per_sample, df_bootstrap, df_raw, results_dir, unknown_sig_config, cosmic_version_label, out_of_reference_sigs, relevance, tool_overrides, df_benchmarks, min_contribution, failed_runs)


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
        tool_benchmarks=[snakemake.input.tool_benchmarks] if getattr(snakemake.input, 'tool_benchmarks', None) else None,
        data_type=snakemake.params.data_type,
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
    parser.add_argument('--title', default='SigFitBench Benchmark Report')
    parser.add_argument('--subtitle', default='Mutational Signature Fitting Tool Comparison')
    parser.add_argument('--results-dir')
    parser.add_argument('--config')
    parser.add_argument('--tool-benchmarks', nargs='+')
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
        tool_benchmarks=args.tool_benchmarks,
    )
