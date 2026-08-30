#!/usr/bin/env python3
"""
Generate an HTML report for no-ground-truth data (realsim / external): cross-tool
agreement instead of TAE/Precision/Recall against a known answer.

This is deliberately a separate report from create_html_report.py (the
ground-truth/simulation report) — it reads a completely different input
(all_tool_agreement.tsv, produced by evaluate_tool_agreement.py +
aggregate_tool_agreement.py) and never touches the simulation subworkflow's
data or code path. It only imports get_html_template()/get_tool_colors() from
create_html_report.py for a visually consistent header/nav/tabs/table style.

Author: SigFitTest Workflow
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from create_html_report import get_html_template, get_tool_colors


# ============================================================================
# Metric display config: (column prefix, display title, colorscale, range)
# ============================================================================

METRICS = [
    ("cosine_similarity", "Cosine Similarity", "Blues", (0, 1)),
    ("pearson_correlation", "Pearson Correlation", "RdBu", (-1, 1)),
    ("jaccard_active", "Jaccard (Active Signatures)", "Greens", (0, 1)),
    ("tae", "TAE (tool vs. tool)", "Reds", None),  # unbounded, no known max
]


def build_matrix(df_cell: pd.DataFrame, tools: list, value_col: str) -> np.ndarray:
    """Build a symmetric tools x tools matrix from a long tool_1/tool_2 table.
    Diagonal is 1.0 for similarity-style metrics (cosine/pearson/jaccard) and
    0.0 for tae (a tool agrees perfectly with itself)."""
    n = len(tools)
    mat = np.full((n, n), np.nan)
    idx = {t: i for i, t in enumerate(tools)}
    diag_value = 0.0 if value_col == "mean_tae" else 1.0
    np.fill_diagonal(mat, diag_value)
    for _, row in df_cell.iterrows():
        i, j = idx.get(row["tool_1"]), idx.get(row["tool_2"])
        if i is None or j is None:
            continue
        mat[i, j] = row[value_col]
        mat[j, i] = row[value_col]
    return mat


def create_agreement_heatmaps(df_cell: pd.DataFrame, tools: list, colors: dict, cell_label: str) -> go.Figure:
    """2x2 grid of heatmaps (cosine, pearson, jaccard, tae) for one (cancer_type,
    relevance, num_muts) cell, tools x tools."""
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[title for _, title, _, _ in METRICS],
        horizontal_spacing=0.15, vertical_spacing=0.2,
    )
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
    for (col_prefix, title, colorscale, zrange), (row, col) in zip(METRICS, positions):
        value_col = f"mean_{col_prefix}"
        if value_col not in df_cell.columns:
            continue
        mat = build_matrix(df_cell, tools, value_col)
        heat_kwargs = dict(
            z=mat, x=tools, y=tools, colorscale=colorscale,
            text=np.round(mat, 3), texttemplate="%{text}",
            hovertemplate=f"%{{y}} vs %{{x}}<br>{title}: %{{z:.3f}}<extra></extra>",
            showscale=False,
        )
        if zrange:
            heat_kwargs["zmin"], heat_kwargs["zmax"] = zrange
        fig.add_trace(go.Heatmap(**heat_kwargs), row=row, col=col)

    fig.update_layout(
        title=f"Cross-Tool Agreement — {cell_label}",
        height=780,
        margin=dict(t=90, b=40, l=60, r=40),
    )
    return fig


def create_summary_table(df_cell: pd.DataFrame) -> str:
    """Small HTML table of every tool-pair's mean metrics + sample count for this cell."""
    cols = [c for c in [
        "tool_1", "tool_2", "n_samples",
        "mean_cosine_similarity", "mean_pearson_correlation",
        "mean_jaccard_active", "mean_tae",
    ] if c in df_cell.columns]
    html = "<table>\n<tr>" + "".join(f"<th>{c.replace('mean_', '').replace('_', ' ').title()}</th>" for c in cols) + "</tr>\n"
    for _, row in df_cell[cols].iterrows():
        html += "<tr>" + "".join(
            f"<td>{row[c]:.3f}</td>" if isinstance(row[c], (float, np.floating)) else f"<td>{row[c]}</td>"
            for c in cols
        ) + "</tr>\n"
    html += "</table>\n"
    return html


_SYMBOL_CYCLE = [
    "circle", "square", "diamond", "triangle-up", "triangle-down",
    "star", "cross", "x", "pentagon", "hexagon2", "triangle-left", "triangle-right",
]


def create_burden_vs_unknown_section(df_activity: pd.DataFrame, tools: list, include_plotly_cdn: bool) -> str:
    """Aggregate scatter across every cancer type/dataset: one point per (sample,
    tool). X-axis is the sample's true mutation burden (from the input catalog,
    the same for every tool — not any tool's own assigned total). Y-axis is the
    percentage of that burden the tool left unassigned ("unknown"). Marker shape
    encodes tool, color encodes cancer type/dataset.

    Both tools and cancer types are deselectable one-by-one from the legend:
    - Each tool is its own real trace, so clicking its legend entry natively
      hides/shows all of that tool's points (standard Plotly behavior, no JS).
    - Cancer types share one legend with decorative (data-less) entries, since
      a single Plotly legend click can only toggle one trace's own visibility,
      not filter points out of a *different* trace by an attribute value.
      A small script (mirroring the legend-click fix used for the Signature
      Activity widget in create_html_report.py: let the default click happen,
      then react to it) listens for the resulting plotly_restyle event, reads
      which cancer entries are currently active, and calls Plotly.restyle to
      rebuild each tool trace's points restricted to the active cancer types —
      still driven entirely by ordinary legend clicks, just resolved in JS
      because two independent legend filters over one dataset isn't natively
      expressible as trace visibility alone.
    """
    div_id = "burden_unknown_plot"
    plot_df = df_activity.dropna(subset=["pct_unknown", "total_mutations"]).copy()
    cancer_types = sorted(plot_df["cancer_type"].unique())
    tool_names = sorted(tools)
    tool_symbols = {t: _SYMBOL_CYCLE[i % len(_SYMBOL_CYCLE)] for i, t in enumerate(tool_names)}
    palette = px.colors.qualitative.Dark24 if len(cancer_types) > 10 else px.colors.qualitative.Plotly
    cancer_colors = {c: palette[i % len(palette)] for i, c in enumerate(cancer_types)}

    fig = go.Figure()
    if plot_df.empty:
        fig.update_layout(title="Mutation Burden vs. Percentage Unassigned — No Data", height=300)
        return fig.to_html(include_plotlyjs=("cdn" if include_plotly_cdn else False), full_html=False, div_id=div_id)

    for i, tool in enumerate(tool_names):
        tsub = plot_df[plot_df["tool"] == tool]
        fig.add_trace(go.Scatter(
            x=tsub["total_mutations"], y=tsub["pct_unknown"], mode="markers", name=tool,
            marker=dict(
                symbol=tool_symbols[tool],
                color=[cancer_colors.get(c, "#888888") for c in tsub["cancer_type"]],
                size=9, opacity=0.75, line=dict(color="black", width=0.5),
            ),
            text=[f"{s} — {tool}, {c}" for s, c in zip(tsub["sample"], tsub["cancer_type"])],
            hovertemplate="%{text}<br>Mutations: %{x:,.0f}<br>Unassigned: %{y:.1f}%<extra></extra>",
            legendgrouptitle_text="Method" if i == 0 else None,
        ))
    for i, cancer_type in enumerate(cancer_types):
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(symbol="circle", color=cancer_colors[cancer_type], size=10),
            name=cancer_type, legendgroup="cancer",
            legendgrouptitle_text="Cancer type / dataset" if i == 0 else None,
        ))
    fig.update_layout(
        title="Mutation Burden vs. Percentage of Mutations Unassigned — All Cancer Types",
        xaxis=dict(title="Mutations per sample (original, from input catalog)", type="log"),
        yaxis=dict(title="Mutations left unassigned by tool (%)"),
        height=680,
        legend=dict(orientation="v", x=1.02, y=1, groupclick="toggleitem"),
        margin=dict(t=60, b=60, l=60, r=220),
        hovermode="closest",
    )

    fig_html = fig.to_html(include_plotlyjs=("cdn" if include_plotly_cdn else False), full_html=False, div_id=div_id)

    full_data = [
        {"tool": t, "cancer": c, "x": float(x), "y": float(y), "text": f"{s} — {t}, {c}", "color": cancer_colors[c]}
        for t, c, s, x, y in zip(
            plot_df["tool"], plot_df["cancer_type"], plot_df["sample"],
            plot_df["total_mutations"], plot_df["pct_unknown"],
        )
    ]
    script = f'''
<script>
(function() {{
    var divId = "{div_id}";
    var fullData = {json.dumps(full_data)};
    var toolNames = {json.dumps(tool_names)};
    var nToolTraces = {len(tool_names)};
    var restyleGuard = false;
    var gd = document.getElementById(divId);
    if (!gd) return;
    gd.on("plotly_restyle", function() {{
        if (restyleGuard) return;
        var activeCancers = {{}};
        var traces = gd.data;
        for (var i = nToolTraces; i < traces.length; i++) {{
            var tr = traces[i];
            if (tr.visible === undefined || tr.visible === true) {{
                activeCancers[tr.name] = true;
            }}
        }}
        var xs = [], ys = [], texts = [], colors = [], idxs = [];
        for (var ti = 0; ti < nToolTraces; ti++) {{
            var toolName = toolNames[ti];
            var fx = [], fy = [], ftext = [], fcolor = [];
            for (var j = 0; j < fullData.length; j++) {{
                var d = fullData[j];
                if (d.tool === toolName && activeCancers[d.cancer]) {{
                    fx.push(d.x); fy.push(d.y); ftext.push(d.text); fcolor.push(d.color);
                }}
            }}
            xs.push(fx); ys.push(fy); texts.push(ftext); colors.push(fcolor); idxs.push(ti);
        }}
        restyleGuard = true;
        Plotly.restyle(gd, {{x: xs, y: ys, text: texts, "marker.color": colors}}, idxs).then(function() {{
            restyleGuard = false;
        }});
    }});
}})();
</script>
'''
    return fig_html + script


def create_signature_activity_boxplot(
    df_sw_cell: pd.DataFrame, active_signatures: list, colors: dict, cell_label: str,
    unknown_label: str = "Unknown_0",
) -> go.Figure:
    """Per-signature, per-tool boxplot of fitted activity across every sample in
    this cell, for every signature any tool called active. One box trace per
    tool: click a tool's legend entry to hide/show its boxes — the same native
    Plotly legend toggle used everywhere else in this report, rather than a
    separate tool picker widget.

    Activity is each signature's weight as a percentage of that tool's total
    assigned weight for that sample (not the raw mutation count), so boxes are
    comparable across samples regardless of how many mutations each one has.

    X-axis signatures are ordered by overall median activity, descending, so
    the most consensus-active signatures come first — except unknown_label
    (SigID's catch-all "Unknown_0" placeholder for unclassified signal), which
    is always pinned last since it isn't a real signature to rank alongside
    the others.
    """
    totals = (
        df_sw_cell.groupby(["tool", "sample"])["weight"].sum()
        .rename("tool_sample_total").reset_index()
    )
    sub = df_sw_cell[df_sw_cell["signature"].isin(active_signatures)].merge(
        totals, on=["tool", "sample"], how="left"
    )
    sub["activity"] = np.where(
        sub["tool_sample_total"] > 0, sub["weight"] / sub["tool_sample_total"] * 100.0, np.nan
    )

    order = sub.groupby("signature")["activity"].median().sort_values(ascending=False).index.tolist()
    if unknown_label in order:
        order.remove(unknown_label)
        order.append(unknown_label)

    fig = go.Figure()
    for tool in sorted(sub["tool"].unique()):
        tsub = sub[sub["tool"] == tool]
        fig.add_trace(go.Box(
            x=tsub["signature"], y=tsub["activity"], name=tool,
            marker_color=colors.get(tool, "#6b7280"),
            boxpoints="outliers",
        ))
    fig.update_layout(
        title=f"Predicted Signature Activity by Tool — {cell_label}",
        xaxis=dict(title="Signature", categoryorder="array", categoryarray=order),
        yaxis=dict(title="Activity (% of sample's assigned mutations)"),
        boxmode="group",
        height=600,
        legend=dict(orientation="h", y=-0.22, x=0.5, xanchor="center"),
        margin=dict(t=60, b=140, l=60, r=40),
    )
    return fig


def create_cell_section(df, df_sig, df_sw, cancer_type, relevance, tools, colors, id_prefix, include_plotly_cdn):
    df_cell = df[(df["cancer_type"] == cancer_type) & (df["relevance"] == relevance)]
    if df_cell.empty:
        return ""
    tab_id = f"{id_prefix}{cancer_type.replace('-', '_')}"
    html = f'<div id="{tab_id}-tab" class="tab-content">\n'
    html += f"<h3>{cancer_type} ({relevance})</h3>\n"
    html += (
        '<p style="color:#666;margin-bottom:1rem;">Pairwise agreement between every selected '
        "tool's fitted signature weights for this cancer type, averaged across all its samples "
        "— there's no ground truth for real samples, so tools are compared to each other instead "
        "of to a known answer. Cosine/Pearson/Jaccard: higher = more agreement. TAE: lower = more "
        "agreement (0 on the diagonal, a tool trivially agrees with itself).</p>\n"
    )
    fig = create_agreement_heatmaps(df_cell, tools, colors, f"{cancer_type} ({relevance})")
    html += '<div class="plot-container">\n'
    html += fig.to_html(include_plotlyjs=("cdn" if include_plotly_cdn else False), full_html=False, div_id=f"agree_{tab_id}")
    html += "</div>\n"
    html += "<h4>Summary</h4>\n"
    html += create_summary_table(df_cell)

    df_sig_cell = df_sig[(df_sig["cancer_type"] == cancer_type) & (df_sig["relevance"] == relevance)] if df_sig is not None else pd.DataFrame()
    df_sw_cell = df_sw[(df_sw["cancer_type"] == cancer_type) & (df_sw["relevance"] == relevance)] if df_sw is not None else pd.DataFrame()
    active_sigs = df_sig_cell[df_sig_cell["n_samples_active"] > 0]["signature"].unique().tolist() if not df_sig_cell.empty else []
    if active_sigs and not df_sw_cell.empty:
        html += '<h4 style="margin-top:2rem;">Which Signatures Are Novel Across Methods?</h4>\n'
        html += (
            '<p style="color:#666;margin-bottom:1rem;">The paper this benchmark is based on '
            "(<a href=\"https://www.nature.com/articles/s41467-024-53711-6\" target=\"_blank\">Medo et al. 2024, "
            "<em>Nat Commun</em></a>) found that tool disagreement on real PCAWG samples concentrates on specific "
            "signatures rather than being spread evenly. Boxes below show each tool's fitted activity (weight) "
            "distribution across all samples, for every signature any tool called active — a signature whose "
            "boxes are non-zero for only one or two tools is that tool's own call, not something every other "
            "tool independently missed. Click a tool's name in the legend to hide/show its boxes.</p>\n"
        )
        box_fig = create_signature_activity_boxplot(df_sw_cell, active_sigs, colors, f"{cancer_type} ({relevance})")
        html += '<div class="plot-container">\n'
        html += box_fig.to_html(include_plotlyjs=False, full_html=False, div_id=f"activity_box_{tab_id}")
        html += "</div>\n"

    html += "</div>\n"
    return html


def create_report_intro(df: pd.DataFrame) -> str:
    tools = sorted(set(df["tool_1"]) | set(df["tool_2"]))
    cancer_types = sorted(df["cancer_type"].unique())
    badges = [
        ("Tools Compared", str(len(tools))),
        ("Cancer Types / Datasets", str(len(cancer_types))),
        ("Tool Pairs", str(len(df[["tool_1", "tool_2"]].drop_duplicates()))),
        ("Total Samples (max)", str(int(df["n_samples"].max())) if "n_samples" in df.columns and len(df) else "N/A"),
    ]
    badges_html = "".join(
        f'<div class="stat-badge"><span class="stat-badge-value">{v}</span>'
        f'<span class="stat-badge-label">{l}</span></div>'
        for l, v in badges
    )
    cancer_types_str = ", ".join(cancer_types)
    return f'''
    <div class="intro">
        <div class="stat-strip">{badges_html}</div>
        <div class="intro-guide">
            <p>
                This report compares <strong>{len(tools)} tools'</strong> fitted signature weights against
                <strong>each other</strong> (not against a known answer) across <strong>{len(cancer_types)}</strong>
                cancer type(s)/dataset(s) ({cancer_types_str}) — because real samples have no ground truth to
                measure against directly.
            </p>
            <ul>
                <li>Each tab below is one cancer type / external dataset: heatmaps of pairwise tool agreement.</li>
                <li><strong>Cosine similarity</strong> and <strong>Pearson correlation</strong> measure shape agreement
                    of the fitted signature-weight vectors; <strong>Jaccard</strong> measures overlap of which
                    signatures each tool called "active"; <strong>TAE</strong> measures the raw magnitude of
                    disagreement (same formula used against ground truth elsewhere, with one tool standing in
                    for "truth").</li>
                <li>High agreement across tools is reassuring but not proof of correctness — all tools could
                    still share the same systematic error.</li>
            </ul>
        </div>
    </div>
    '''


def generate_report(df: pd.DataFrame, df_sig: pd.DataFrame, df_activity: pd.DataFrame, df_sw: pd.DataFrame, output_file: str, title: str, subtitle: str) -> None:
    colors = get_tool_colors()

    if df.empty:
        content = (
            '<div class="section"><h2>No Data</h2>'
            "<p>all_tool_agreement.tsv has no rows — check that realsim/external tools actually "
            "produced contribution files, and that get_agreement_tools() selected more than one tool "
            "(agreement needs at least 2 tools to compare).</p></div>"
        )
    else:
        tools = sorted(set(df["tool_1"]) | set(df["tool_2"]))
        cancer_types = sorted(df["cancer_type"].unique())
        relevances = sorted(df["relevance"].unique())

        intro_html = create_report_intro(df)

        has_activity = df_activity is not None and not df_activity.empty

        tabs_html = '<section id="visualizations" class="section">\n<h2>📈 Cross-Tool Agreement</h2>\n'
        tabs_html += '<div class="tabs">\n'
        if has_activity:
            tabs_html += '<button class="tab active" onclick="switchTab(\'burden_unknown-tab\', event)">All Cancers — Burden vs. Unknown</button>\n'
        first_cell_tab = not has_activity
        for cancer_type in cancer_types:
            for relevance in relevances:
                if df[(df["cancer_type"] == cancer_type) & (df["relevance"] == relevance)].empty:
                    continue
                tab_id = f"{cancer_type.replace('-', '_')}"
                active = "active" if first_cell_tab else ""
                tabs_html += f'<button class="tab {active}" onclick="switchTab(\'{tab_id}-tab\', event)">{cancer_type}</button>\n'
                first_cell_tab = False
        tabs_html += "</div>\n"

        if has_activity:
            tabs_html += '<div id="burden_unknown-tab" class="tab-content active">\n'
            tabs_html += "<h3>Mutation Burden vs. Percentage of Mutations Unassigned</h3>\n"
            tabs_html += (
                '<p style="color:#666;margin-bottom:1rem;">Every point is one sample fit by one tool, pooled '
                "across all cancer types/datasets. X-axis: the sample's true mutation burden, taken from the "
                "input catalog rather than any tool's own fit (log scale). Y-axis: the percentage of that "
                "burden the tool's fit left unassigned — <code>(true total − tool's assigned total) / true "
                "total</code>. Marker shape identifies the method, color identifies the cancer type/dataset, "
                "so both cross-tool and cross-cancer patterns in the unassigned fraction are visible together. "
                "A tool that leaves a growing share unassigned as burden increases may be under-fitting "
                "higher-mutation-count samples; one whose unassigned share differs sharply by cancer type may "
                "be sensitive to catalog composition for that tissue.</p>\n"
            )
            tabs_html += '<div class="plot-container">\n'
            tabs_html += create_burden_vs_unknown_section(df_activity, tools, include_plotly_cdn=True)
            tabs_html += "</div>\n</div>\n"

        first = not has_activity
        for cancer_type in cancer_types:
            for relevance in relevances:
                section = create_cell_section(df, df_sig, df_sw, cancer_type, relevance, tools, colors, "", include_plotly_cdn=first)
                if section:
                    section = section.replace('class="tab-content"', 'class="tab-content active"', 1) if first else section
                    tabs_html += section
                    first = False
        tabs_html += "</section>\n"

        content = intro_html + tabs_html

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html = get_html_template().format(title=title, subtitle=subtitle, content=content, timestamp=timestamp)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        f.write(html)
    print(f"\n✅ Agreement report generated: {output_file}")


def _load_tsv_or_empty(path):
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def main(tool_agreement_file, signature_calls_file, sample_activity_file, signature_weights_file, output, title, subtitle):
    df = _load_tsv_or_empty(tool_agreement_file)
    df_sig = _load_tsv_or_empty(signature_calls_file)
    df_activity = _load_tsv_or_empty(sample_activity_file)
    df_sw = _load_tsv_or_empty(signature_weights_file)
    generate_report(df, df_sig, df_activity, df_sw, output, title, subtitle)


if "snakemake" in dir():
    import sys
    sys.stdout = open(snakemake.log.stdout, "w")
    sys.stderr = open(snakemake.log.stderr, "w")
    main(
        tool_agreement_file=snakemake.input.tool_agreement,
        signature_calls_file=snakemake.input.signature_calls,
        sample_activity_file=snakemake.input.sample_activity,
        signature_weights_file=snakemake.input.signature_weights,
        output=snakemake.output[0],
        title=snakemake.params.title,
        subtitle=snakemake.params.subtitle,
    )
elif __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate cross-tool agreement HTML report")
    parser.add_argument("--tool-agreement", required=True)
    parser.add_argument("--signature-calls", required=True)
    parser.add_argument("--sample-activity", required=True)
    parser.add_argument("--signature-weights", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--title", default="SigFitBench Cross-Tool Agreement Report")
    parser.add_argument("--subtitle", default="No-ground-truth comparison (realsim / external)")
    args = parser.parse_args()
    main(args.tool_agreement, args.signature_calls, args.sample_activity, args.signature_weights, args.output, args.title, args.subtitle)
