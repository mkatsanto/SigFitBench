#!/usr/bin/env python3
"""
HTML report for real samples, where the true signatures are not known, so the
tools are compared with each other.

  Overview   main findings in a few sentences, agreement vs. number of mutations,
             signatures found per sample by each tool, agreement per group and tool pair.
  By group   one group at a time (cancer type, plus relevance / mutation number when
             the run has more than one): agreement per tool pair, which signatures
             each tool finds, and the samples the tools disagree on most.

Agreement between two tools for a sample is the cosine similarity of their
signature weights (1 = same proportions, 0 = nothing in common).
The page style comes from create_html_report.py so both reports look alike.
"""

import argparse
import html
import os
from datetime import datetime
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from create_html_report import PLOT_CONFIG, get_html_template, get_tool_colors

MIN_PCT_SAMPLES = 5  # signatures found in fewer samples than this (%) by every tool are not drawn
N_WORST_SAMPLES = 5  # samples listed in "samples the tools disagree on most"

PAGE_EXTRAS = """
<style>
  .note { color:#555; margin:0.25rem 0 1rem 0; max-width:60rem; }
  .findings { background:#f3f6ff; border:1px solid #dbe3ff; border-radius:10px; padding:1rem 1.25rem; margin:1rem 0 2rem; }
  .group { display:none; } .group.shown { display:block; }
  select { font-size:1rem; padding:0.4rem 0.6rem; border-radius:6px; max-width:100%; margin-bottom:1rem; }
  details summary { cursor:pointer; color:#1d4ed8; margin-top:1rem; }
</style>
<script>
function showGroup(id) {
  document.querySelectorAll('.group').forEach(el => el.classList.toggle('shown', el.id === id));
  setTimeout(() => document.querySelectorAll('#' + id + ' .js-plotly-plot').forEach(p => Plotly.Plots.resize(p)), 30);
}
document.addEventListener('DOMContentLoaded', () => {
  const sel = document.getElementById('group-select');
  if (sel) document.querySelectorAll('.tab').forEach(t => t.addEventListener('click', () => showGroup(sel.value)));
  if (sel) showGroup(sel.value);
});
</script>
"""


# ---------------------------------------------------------------- calculations

def add_group(df: pd.DataFrame, extra_cols: list) -> pd.DataFrame:
    """Label each row with its group: cancer type, plus any setting that varies in this run."""
    if df.empty:
        return df
    label = df["cancer_type"].astype(str)
    for c in extra_cols:
        label = label + " · " + df[c].astype(str)
    return df.assign(group=label)


def per_sample_agreement(df_sw: pd.DataFrame, tools: list) -> pd.DataFrame:
    """Cosine similarity of signature weights for every sample and tool pair, plus their mean."""
    if df_sw.empty:
        return pd.DataFrame(columns=["group", "sample", "mean_agreement"])
    wide = df_sw.pivot_table(index=["group", "sample", "signature"], columns="tool", values="weight", aggfunc="sum")
    pairs = {}
    for t1, t2 in combinations([t for t in tools if t in wide.columns], 2):
        both = wide[[t1, t2]].dropna()
        by = [both.index.get_level_values(0), both.index.get_level_values(1)]
        dot = (both[t1] * both[t2]).groupby(by).sum()
        norm = np.sqrt((both[t1] ** 2).groupby(by).sum() * (both[t2] ** 2).groupby(by).sum())
        pairs[f"{t1} vs {t2}"] = dot / norm.replace(0, np.nan)
    out = pd.DataFrame(pairs)
    out.index.names = ["group", "sample"]
    return out.assign(mean_agreement=out.mean(axis=1)).reset_index()


# ---------------------------------------------------------------- figures

def to_html(fig, div_id, first=False):
    return '<div class="plot-container">' + fig.to_html(
        include_plotlyjs=first, full_html=False, div_id=div_id, config=PLOT_CONFIG) + "</div>\n"


def fig_agreement_vs_mutations(ps: pd.DataFrame) -> go.Figure:
    d = ps.dropna(subset=["mean_agreement", "mutations"]).query("mutations > 0")
    fig = go.Figure(go.Scatter(
        x=d["mutations"], y=d["mean_agreement"], mode="markers", showlegend=False,
        marker=dict(color="rgba(107,114,128,0.35)", size=6),
        text=[f"{html.escape(str(s))}<br>{g}" for s, g in zip(d["sample"], d["group"])],
        hovertemplate="%{text}<br>Mutations: %{x:,.0f}<br>Agreement: %{y:.2f}<extra></extra>",
    ))
    if len(d) >= 20:  # median line over log-spaced bins with at least 10 samples
        bins = pd.cut(np.log10(d["mutations"]), 12)
        med = d.groupby(bins, observed=True).agg(x=("mutations", "median"), y=("mean_agreement", "median"),
                                                 n=("mutations", "size")).query("n >= 10")
        fig.add_trace(go.Scatter(x=med["x"], y=med["y"], mode="lines+markers", showlegend=False,
                                 line=dict(color="#1d4ed8", width=3), hoverinfo="skip"))
    fig.update_layout(xaxis=dict(title="Mutations in sample", type="log"),
                      yaxis=dict(title="Agreement between tools", range=[0, 1.02]),
                      height=480, margin=dict(t=30, b=60, l=70, r=30))
    return fig


def fig_signatures_per_sample(act: pd.DataFrame, colors: dict) -> go.Figure:
    order = act.groupby("tool")["n_active_signatures"].median().sort_values(ascending=False).index
    fig = go.Figure([go.Box(y=act.loc[act["tool"] == t, "n_active_signatures"], name=t, boxpoints=False,
                            marker_color=colors.get(t, "#6b7280"), showlegend=False) for t in order])
    fig.update_layout(yaxis=dict(title="Signatures found per sample", rangemode="tozero"),
                      height=420, margin=dict(t=30, b=50, l=70, r=30))
    return fig


def fig_group_by_pair(df: pd.DataFrame, n_samples: pd.Series) -> go.Figure:
    """Mean agreement per group (rows, lowest first) and tool pair (columns)."""
    mat = df.assign(pair=df["tool_1"] + " vs " + df["tool_2"]).pivot_table(
        index="group", columns="pair", values="mean_cosine_similarity")
    mat = mat[mat.mean().sort_values(ascending=False).index]
    mat = mat.loc[mat.mean(axis=1).sort_values().index]
    fig = go.Figure(go.Heatmap(
        z=mat.values, x=list(mat.columns), y=[f"{g} ({n_samples[g]})" for g in mat.index],
        colorscale="Blues", zmin=min(0.5, np.nanmin(mat.values)), zmax=1, texttemplate="%{z:.2f}",
        hovertemplate="%{y}<br>%{x}<br>Agreement: %{z:.2f}<extra></extra>", colorbar=dict(title="Agreement")))
    fig.update_layout(height=max(400, 24 * len(mat) + 160), yaxis=dict(autorange="reversed"),
                      xaxis=dict(side="top", tickangle=-30), margin=dict(t=120, b=20, l=220, r=30))
    return fig


def fig_pair_heatmaps(cell: pd.DataFrame, tools: list) -> go.Figure:
    panels = [("mean_cosine_similarity", "Similarity of signature weights", "Blues"),
              ("mean_jaccard_active", "Overlap of signatures found", "Greens")]
    fig = make_subplots(rows=1, cols=2, subplot_titles=[p[1] for p in panels], horizontal_spacing=0.18)
    for k, (col, title, scale) in enumerate(panels, start=1):
        mat = pd.DataFrame(np.eye(len(tools)), index=tools, columns=tools)
        for _, r in cell.iterrows():
            mat.loc[r["tool_1"], r["tool_2"]] = mat.loc[r["tool_2"], r["tool_1"]] = r[col]
        fig.add_trace(go.Heatmap(z=mat.values, x=tools, y=tools, colorscale=scale, zmin=0, zmax=1,
                                 texttemplate="%{z:.2f}", showscale=False,
                                 hovertemplate=f"%{{y}} vs %{{x}}<br>{title}: %{{z:.2f}}<extra></extra>"),
                      row=1, col=k)
        fig.update_yaxes(autorange="reversed", row=1, col=k)
    fig.update_layout(height=60 * len(tools) + 200, margin=dict(t=50, b=40, l=90, r=30))
    return fig


def fig_signature_calls(sig: pd.DataFrame, sw: pd.DataFrame, tools: list):
    """% of samples in which each tool finds each signature. Returns (figure or None, signatures left out)."""
    mat = sig.pivot_table(index="signature", columns="tool", values="frac_samples_active").reindex(columns=tools) * 100
    keep = (mat.fillna(0) >= MIN_PCT_SAMPLES).any(axis=1)
    n_left_out = int(((mat.fillna(0) > 0).any(axis=1) & ~keep).sum())
    mat = mat[keep]
    if mat.empty:
        return None, n_left_out
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    share = sw.assign(share=sw["weight"] / sw.groupby(["tool", "sample"])["weight"].transform("sum").replace(0, np.nan) * 100)
    share = share.pivot_table(index="signature", columns="tool", values="share").reindex(index=mat.index, columns=tools)
    fig = go.Figure(go.Heatmap(
        z=mat.values, x=tools, y=list(mat.index), customdata=share.values, colorscale="Purples", zmin=0, zmax=100,
        texttemplate="%{z:.0f}", colorbar=dict(title="% samples"),
        hovertemplate="%{y} — %{x}<br>Found in %{z:.0f}% of samples"
                      "<br>Average share of mutations: %{customdata:.1f}%<extra></extra>"))
    fig.update_layout(height=max(300, 22 * len(mat) + 130), yaxis=dict(autorange="reversed"),
                      xaxis=dict(side="top"), margin=dict(t=60, b=20, l=110, r=30))
    return fig, n_left_out


# ---------------------------------------------------------------- text and tables

def table(df: pd.DataFrame) -> str:
    fmt = lambda v: "" if pd.isna(v) else (f"{v:,.0f}" if v.is_integer() else f"{v:.3f}") if isinstance(
        v, float) else html.escape(str(v))
    head = "".join(f"<th>{c}</th>" for c in df.columns)
    rows = "".join("<tr>" + "".join(f"<td>{fmt(v)}</td>" for v in r) + "</tr>" for r in df.itertuples(index=False))
    return f"<table><tr>{head}</tr>{rows}</table>"


def key_findings(df, act, ps, tools) -> str:
    items = []
    pairs = df.groupby(["tool_1", "tool_2"])["mean_cosine_similarity"].mean().sort_values()
    if len(pairs) >= 2:
        items.append(f"<b>{' and '.join(pairs.index[-1])}</b> agree most (average {pairs.iloc[-1]:.2f}); "
                     f"<b>{' and '.join(pairs.index[0])}</b> agree least ({pairs.iloc[0]:.2f}).")
    if not act.empty:
        med = act.groupby("tool")["n_active_signatures"].median().sort_values()
        if med.iloc[-1] > med.iloc[0]:
            items.append(f"Signatures found per sample (median): " + ", ".join(
                f"<b>{t}</b> {v:.0f}" for t, v in med[::-1].items()) + ".")
    d = ps.dropna(subset=["mean_agreement", "mutations"])
    if len(d) >= 40:
        lo, hi = d["mutations"].quantile([0.25, 0.75])
        a_lo = d.loc[d["mutations"] <= lo, "mean_agreement"].median()
        a_hi = d.loc[d["mutations"] >= hi, "mean_agreement"].median()
        items.append(f"Median agreement is {a_lo:.2f} in the quarter of samples with the fewest mutations "
                     f"(&le;{lo:,.0f}) and {a_hi:.2f} in the quarter with the most (&ge;{hi:,.0f}).")
    groups = df.groupby("group")["mean_cosine_similarity"].mean().sort_values()
    if len(groups) >= 6:
        items.append("Lowest agreement: " + ", ".join(f"{g} ({v:.2f})" for g, v in groups.head(3).items())
                     + ". Highest: " + ", ".join(f"{g} ({v:.2f})" for g, v in groups.tail(3)[::-1].items()) + ".")
    if not act.empty:
        all_groups = set(act["group"])
        for t in tools:
            absent = sorted(all_groups - set(act.loc[act["tool"] == t, "group"]))
            if absent:
                items.append(f"<b>{t}</b> has no results for: {', '.join(absent)}.")
        if "pct_unknown" in act.columns and act["pct_unknown"].notna().any():
            unk = act.groupby("tool")["pct_unknown"].agg(["mean", "max"]).query("max > 0")
            items.append("Mutations left unassigned: " + (", ".join(
                f"<b>{t}</b> {r['mean']:.1f}% on average (at most {r['max']:.1f}%)" for t, r in unk.iterrows())
                if len(unk) else "none") + (f"; the other tools assign every mutation" if 0 < len(unk) < len(tools) else "") + ".")
    return '<div class="findings"><h3>Main findings</h3><ul>' + "".join(f"<li>{i}</li>" for i in items) + "</ul></div>"


def group_section(g, gid, cell, sig, sw, ps, act, rank, n_groups) -> str:
    tools = sorted(set(cell["tool_1"]) | set(cell["tool_2"]))  # only tools with results here
    h = (f'<div id="{gid}" class="group"><h3>{html.escape(g)}</h3>'
         f'<p><b>{int(cell["n_samples"].max()):,}</b> samples · agreement <b>{cell["mean_cosine_similarity"].mean():.2f}</b> '
         f"(rank {rank} of {n_groups}, 1 = highest)</p>"
         "<h4>Agreement between each pair of tools</h4>"
         '<p class="note">Left: how similar the signature weights are (1 = same proportions). '
         "Right: share of signatures found by either tool that both tools found. Averages over all samples.</p>"
         + to_html(fig_pair_heatmaps(cell, tools), f"{gid}_pairs"))
    if not sig.empty and not sw.empty:
        fig, n_left_out = fig_signature_calls(sig, sw, tools)
        if fig is not None:
            h += ("<h4>Which signatures each tool finds</h4>"
                  '<p class="note">Percentage of samples in which each tool finds each signature; a row filled for one '
                  "tool only is a signature the others do not find. Hover for the average share of mutations."
                  + (f" {n_left_out} signatures found in under {MIN_PCT_SAMPLES}% of samples by every tool are not shown."
                     if n_left_out else "") + "</p>" + to_html(fig, f"{gid}_sigs"))
    worst = ps.dropna(subset=["mean_agreement"]).nsmallest(N_WORST_SAMPLES, "mean_agreement")
    if not worst.empty:
        counts = act.pivot_table(index="sample", columns="tool", values="n_active_signatures").reindex(columns=tools)
        worst = worst[["sample"] + (["mutations"] if "mutations" in worst else []) + ["mean_agreement"]].merge(
            counts, left_on="sample", right_index=True, how="left")
        worst.columns = ["Sample"] + (["Mutations"] if "mutations" in ps else []) + ["Agreement"] + tools
        h += ("<h4>Samples the tools disagree on most</h4>" + table(worst)
              + '<p class="note">Tool columns: number of signatures each tool found in the sample.</p>')
    cols = {"tool_1": "Tool 1", "tool_2": "Tool 2", "n_samples": "Samples",
            "mean_cosine_similarity": "Similarity (cosine)", "mean_pearson_correlation": "Pearson correlation",
            "mean_jaccard_active": "Overlap (Jaccard)", "mean_tae": "Total absolute error"}
    cols = {k: v for k, v in cols.items() if k in cell}
    h += f"<details><summary>All numbers for each tool pair</summary>{table(cell[list(cols)].rename(columns=cols))}</details>"
    return h + "</div>\n"


# ---------------------------------------------------------------- page

def generate_report(df, df_sig, df_act, df_sw, output_file, title, subtitle):
    if df.empty:
        content = ('<div class="section"><h2>No data</h2><p>summary/tool_agreement.tsv has no rows. '
                   "Agreement needs results from at least two tools.</p></div>")
    else:
        extra = [c for c in ["relevance", "num_muts"] if c in df and df[c].nunique() > 1]
        df, df_sig, df_act, df_sw = (add_group(x, extra) for x in (df, df_sig, df_act, df_sw))
        tools = sorted(set(df["tool_1"]) | set(df["tool_2"]))
        groups = sorted(df["group"].unique())
        n_samples = df.groupby("group")["n_samples"].max()

        ps = per_sample_agreement(df_sw, tools)
        if "total_mutations" in df_act:
            muts = df_act.groupby(["group", "sample"])["total_mutations"].max().rename("mutations")
            ps = ps.merge(muts.reset_index(), on=["group", "sample"], how="left")
        n_total = len(df_act[["group", "sample"]].drop_duplicates()) if not df_act.empty else int(n_samples.sum())
        median_agree = ps["mean_agreement"].median() if len(ps) else df["mean_cosine_similarity"].mean()

        badges = "".join(f'<div class="stat-badge"><span class="stat-badge-value">{v}</span>'
                         f'<span class="stat-badge-label">{k}</span></div>' for k, v in
                         [("Samples", f"{n_total:,}"), ("Groups", len(groups)), ("Tools", len(tools)),
                          ("Median agreement", f"{median_agree:.2f}")])
        intro = (f'<div class="intro"><div class="stat-strip">{badges}</div><div class="intro-guide"><p>'
                 f"Real samples have no known answer, so this report compares the tools ({', '.join(tools)}) with "
                 "each other. <b>Agreement</b> is how similar two tools' signature weights are for a sample "
                 "(1 = same proportions, 0 = nothing in common). Tools that agree can still share the same mistake."
                 "</p></div></div>")

        ov = '<div id="overview-tab" class="tab-content active">' + key_findings(df, df_act, ps, tools)
        plots = []
        if "mutations" in ps and ps["mutations"].notna().any():
            plots.append(("Do tools agree more on samples with more mutations?",
                          "Each point is one sample: the average agreement over all tool pairs. The line is the median.",
                          fig_agreement_vs_mutations(ps)))
        if "n_active_signatures" in df_act:
            plots.append(("How many signatures does each tool find?",
                          "Signatures found per sample, all groups together.",
                          fig_signatures_per_sample(df_act, get_tool_colors(tools))))
        plots.append(("Agreement by group", "Average agreement for each tool pair, lowest at the top; "
                      "number of samples in brackets.", fig_group_by_pair(df, n_samples)))
        for i, (head, note, fig) in enumerate(plots):
            ov += f'<h3>{head}</h3><p class="note">{note}</p>' + to_html(fig, f"overview_{i}", first=(i == 0))
        ov += "</div>\n"

        rank = df.groupby("group")["mean_cosine_similarity"].mean().rank(ascending=False, method="min").astype(int)
        pick = lambda x, g: x[x["group"] == g] if not x.empty else x
        options, sections = "", ""
        for i, g in enumerate(groups):
            cell = df[df["group"] == g]
            options += (f'<option value="g{i}">{html.escape(g)} — {n_samples[g]} samples, '
                        f'agreement {cell["mean_cosine_similarity"].mean():.2f}</option>')
            sections += group_section(g, f"g{i}", cell, pick(df_sig, g), pick(df_sw, g), pick(ps, g),
                                      pick(df_act, g), rank[g], len(groups))
        by = ('<div id="bygroup-tab" class="tab-content">'
              f'<select id="group-select" onchange="showGroup(this.value)">{options}</select>{sections}</div>')

        content = (PAGE_EXTRAS + intro + '<section class="section"><div class="tabs">'
                   '<button class="tab active" onclick="switchTab(\'overview-tab\', event)">Overview</button>'
                   '<button class="tab" onclick="switchTab(\'bygroup-tab\', event)">By cancer type</button>'
                   f"</div>{ov}{by}</section>")

    page = get_html_template().format(title=title, subtitle=subtitle, content=content,
                                      timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(output_file).write_text(page)
    print(f"Agreement report written: {output_file}")


def read_tsv(path):
    return pd.read_csv(path, sep="\t") if path and os.path.exists(path) and os.path.getsize(path) > 0 else pd.DataFrame()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="HTML report comparing tools on real samples")
    for arg in ["--tool-agreement", "--signature-calls", "--sample-activity", "--signature-weights", "--output"]:
        p.add_argument(arg, required=True)
    p.add_argument("--title", default="SigFitBench: Real samples")
    p.add_argument("--subtitle", default="Agreement between tools (no known answer)")
    a = p.parse_args()
    generate_report(read_tsv(a.tool_agreement), read_tsv(a.signature_calls), read_tsv(a.sample_activity),
                    read_tsv(a.signature_weights), a.output, a.title, a.subtitle)
