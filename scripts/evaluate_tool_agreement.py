#!/usr/bin/env python3
"""
Cross-tool agreement for one (data_type, cancer_type, relevance, num_muts) cell.

Used where no ground truth exists (real data): instead of comparing each tool's fit to a known answer, tools are
compared pairwise to each other, per sample, on their fitted signature-weight
vectors (aligned on the union of signature names, missing filled with 0):

  - cosine_similarity   dot(a,b) / (||a|| ||b||)          — shape agreement
  - pearson_correlation across signature weights            — linear agreement
  - tae                 sum(|a-b|) / (2*total_mutations)  — same TAE definition
                         used against ground truth elsewhere in this workflow,
                         just with one tool standing in for "truth"
  - jaccard_active      overlap of {sig: weight > min_contribution} sets

Per-sample values are averaged per tool pair into one summary row.
"""

import argparse
import re
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from evaluate_tool_performance import load_tool_results


def load_original_catalog(path: str) -> pd.Series:
    """Load the input mutation-count catalog (contexts x samples) used to fit
    every tool for this cell, and return each sample's true mutation total —
    the sum of its column, independent of any tool's fit. This is the
    denominator for "percentage of mutations left unassigned": a tool's own
    contribution-weight sum can fall short of this if it doesn't explain 100%
    of the sample's mutations.

    Separator varies by which format was generated (tab for most tools' plain
    to_csv(), comma for sigfit/sigLASSO), so it's detected from the first line
    rather than assumed, same reasoning as load_contribution above.
    """
    with open(path) as f:
        first_line = f.readline()
    sep = "," if first_line.count(",") > first_line.count("\t") else "\t"
    catalog = pd.read_csv(path, sep=sep, index_col=0)
    catalog = catalog.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return catalog.sum(axis=0)


def load_contribution(path: str) -> pd.DataFrame:
    """Load a {tool}_contribution.dat file: signatures x samples.

    Reuses load_tool_results (evaluate_tool_performance.py) rather than
    hardcoding sep="\\t" — contribution.dat format actually varies by tool
    (comma-separated for most Python tools' plain to_csv() output, tab-separated
    for R-script tools), which load_tool_results already auto-detects correctly
    for the ground-truth evaluation path. A hardcoded tab separator here silently
    parsed comma-separated files as a single column, producing an empty (0-sample)
    DataFrame for those tools and zero rows of agreement output — no error, just
    nothing to compare.
    """
    df = load_tool_results(path)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # Drop SigID's internal "Unknown_i" placeholder rows so pct_unknown here
    # matches evaluate_tool_performance.py, which also excludes them.
    keep_rows = [sig for sig in df.index if not str(sig).startswith("Unknown_")]
    return df.loc[keep_rows]


def cosine_similarity(a: pd.Series, b: pd.Series) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def pearson_correlation(a: pd.Series, b: pd.Series) -> float:
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(pearsonr(a, b)[0])


def tae(a: pd.Series, b: pd.Series) -> float:
    """Same definition as recalculate_TAE in evaluate_tool_performance.py, with
    tool_a's total standing in for the mutation count normaliser."""
    total = a.sum()
    if total == 0:
        return np.nan
    return float((a - b).abs().sum() / (2 * total))


def jaccard(active_a: set, active_b: set) -> float:
    union = active_a | active_b
    if not union:
        return np.nan
    return len(active_a & active_b) / len(union)


def pairwise_agreement(contributions: dict, min_contribution: float) -> pd.DataFrame:
    """Return a per-sample, per-tool-pair long table of agreement metrics."""
    rows = []
    for tool_a, tool_b in combinations(sorted(contributions), 2):
        df_a, df_b = contributions[tool_a], contributions[tool_b]
        common_samples = df_a.columns.intersection(df_b.columns)
        all_sigs = df_a.index.union(df_b.index)
        for sample in common_samples:
            a = df_a[sample].reindex(all_sigs, fill_value=0.0)
            b = df_b[sample].reindex(all_sigs, fill_value=0.0)
            active_a = set(a.index[a > min_contribution])
            active_b = set(b.index[b > min_contribution])
            rows.append({
                "tool_1": tool_a,
                "tool_2": tool_b,
                "sample": sample,
                "cosine_similarity": cosine_similarity(a, b),
                "pearson_correlation": pearson_correlation(a, b),
                "tae": tae(a, b),
                "jaccard_active": jaccard(active_a, active_b),
                "n_active_1": len(active_a),
                "n_active_2": len(active_b),
            })
    return pd.DataFrame(rows)


def summarise(per_sample: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the per-sample long table to one row per tool pair."""
    metric_cols = ["cosine_similarity", "pearson_correlation", "tae", "jaccard_active"]
    if per_sample.empty:
        cols = ["tool_1", "tool_2", "n_samples"]
        for m in metric_cols:
            cols += [f"mean_{m}", f"std_{m}"]
        return pd.DataFrame(columns=cols)

    grouped = per_sample.groupby(["tool_1", "tool_2"])
    agg = {"sample": "count"}
    for m in metric_cols:
        agg[m] = ["mean", "std"]
    summary = grouped.agg(agg)
    summary.columns = ["n_samples"] + [f"{stat}_{m}" for m in metric_cols for stat in ["mean", "std"]]
    return summary.reset_index()


def signature_call_consensus(contributions: dict, min_contribution: float) -> pd.DataFrame:
    """Per (tool, signature): how often that tool calls that signature active.

    Used to spot "novel"/tool-specific signature calls when there's no ground
    truth to check against directly — a signature lit up by only one tool across
    most samples is much more likely to be that tool's own artifact than a real
    signal every other tool independently missed. This is the same phenomenon
    the SigFitTest paper (Medo et al. 2024, Nat Commun) found driving tool
    disagreement on real PCAWG samples: it concentrates on specific signatures
    rather than being spread evenly, especially ones absent from the fitting
    reference catalog.
    """
    n_tools = len(contributions)
    all_sigs = sorted(set().union(*(df.index for df in contributions.values()))) if contributions else []
    rows = []
    for sig in all_sigs:
        n_samples_by_tool = {}
        for tool, df in contributions.items():
            if sig not in df.index:
                n_samples_by_tool[tool] = 0
                continue
            n_samples_by_tool[tool] = int((df.loc[sig] > min_contribution).sum())
        n_tools_calling = sum(1 for n in n_samples_by_tool.values() if n > 0)
        for tool, n_active in n_samples_by_tool.items():
            n_total = contributions[tool].shape[1]
            rows.append({
                "signature": sig,
                "tool": tool,
                "n_samples_active": n_active,
                "n_samples_total": n_total,
                "frac_samples_active": (n_active / n_total) if n_total else np.nan,
                "n_tools_calling_this_signature": n_tools_calling,
                "n_tools_total": n_tools,
            })
    return pd.DataFrame(rows)


def signature_weight_table(contributions: dict) -> pd.DataFrame:
    """Long-format (tool, signature, sample) -> fitted weight, for every signature
    each tool's fit includes (not just ones it called active). Used to draw
    per-signature, per-tool boxplots of predicted activity — the raw distribution
    behind the call-consensus numbers in signature_call_consensus, so a signature
    can be inspected sample-by-sample rather than only via its active/inactive
    summary.
    """
    rows = []
    for tool, df in contributions.items():
        for sig in df.index:
            for sample in df.columns:
                rows.append({
                    "tool": tool,
                    "signature": sig,
                    "sample": sample,
                    "weight": float(df.loc[sig, sample]),
                })
    return pd.DataFrame(rows)


def sample_activity_summary(
    contributions: dict, min_contribution: float, original_totals: pd.Series = None
) -> pd.DataFrame:
    """Per (sample, tool): true mutation burden, number of signatures called
    active, and the percentage of that burden the tool left unassigned.

    total_mutations is the sample's true mutation count from the input catalog
    (original_totals, same for every tool) when available, falling back to the
    tool's own assigned-weight sum otherwise. pct_unknown is
    (true_total - assigned_total) / true_total * 100, clipped at 0 — the
    fraction of the sample's mutations the tool's fit doesn't account for.
    """
    rows = []
    for tool, df in contributions.items():
        for sample in df.columns:
            col = df[sample]
            assigned_total = float(col.sum())
            true_total = (
                float(original_totals.get(sample, np.nan))
                if original_totals is not None else np.nan
            )
            if np.isnan(true_total) or true_total <= 0:
                true_total = assigned_total
                pct_unknown = np.nan
            else:
                pct_unknown = max(0.0, (true_total - assigned_total) / true_total * 100.0)
            rows.append({
                "sample": sample,
                "tool": tool,
                "total_mutations": true_total,
                "assigned_mutations": assigned_total,
                "n_active_signatures": int((col > min_contribution).sum()),
                "pct_unknown": pct_unknown,
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Compute cross-tool agreement (no ground truth needed)")
    parser.add_argument("--tool-results", nargs="+", required=True,
                        help="Paths to {tool}_contribution.dat files, one per tool")
    parser.add_argument("--tools", nargs="+", required=True,
                        help="Tool name for each --tool-results path, same order")
    parser.add_argument("--min-contribution", type=float, default=10.0,
                        help="Weight threshold above which a signature counts as 'active'")
    parser.add_argument("--data-type", required=True)
    parser.add_argument("--cancer-type", required=True)
    parser.add_argument("--relevance", required=True)
    parser.add_argument("--num-muts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--output-signature-calls", required=True,
                        help="Path for the per-signature, per-tool call-consensus table "
                             "(how often each tool calls each signature active)")
    parser.add_argument("--output-sample-activity", required=True,
                        help="Path for the per-sample, per-tool (mutation burden, "
                             "active signature count) table")
    parser.add_argument("--output-signature-weights", required=True,
                        help="Path for the long-format (tool, signature, sample) -> "
                             "fitted weight table, used to draw per-signature, "
                             "per-tool boxplots of predicted activity")
    parser.add_argument("--original-catalog", default=None,
                        help="Path to the input mutation-count catalog (contexts x "
                             "samples) used to fit every tool in this cell — gives "
                             "each sample's true mutation total, for computing what "
                             "percentage of it a tool left unassigned. Optional: "
                             "falls back to each tool's own assigned total if omitted.")
    args = parser.parse_args()

    if len(args.tool_results) != len(args.tools):
        raise ValueError("--tool-results and --tools must have the same length")

    contributions = {
        tool: load_contribution(path)
        for tool, path in zip(args.tools, args.tool_results)
    }

    per_sample = pairwise_agreement(contributions, args.min_contribution)
    summary = summarise(per_sample)
    summary.insert(0, "num_muts", args.num_muts)
    summary.insert(0, "relevance", args.relevance)
    summary.insert(0, "cancer_type", args.cancer_type)
    summary.insert(0, "data_type", args.data_type)

    summary.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output} ({len(summary)} tool pairs, {len(contributions)} tools)")

    sig_calls = signature_call_consensus(contributions, args.min_contribution)
    sig_calls.insert(0, "num_muts", args.num_muts)
    sig_calls.insert(0, "relevance", args.relevance)
    sig_calls.insert(0, "cancer_type", args.cancer_type)
    sig_calls.insert(0, "data_type", args.data_type)
    sig_calls.to_csv(args.output_signature_calls, sep="\t", index=False)
    print(f"Wrote {args.output_signature_calls} ({len(sig_calls)} (tool, signature) rows)")

    sig_weights = signature_weight_table(contributions)
    sig_weights.insert(0, "num_muts", args.num_muts)
    sig_weights.insert(0, "relevance", args.relevance)
    sig_weights.insert(0, "cancer_type", args.cancer_type)
    sig_weights.insert(0, "data_type", args.data_type)
    sig_weights.to_csv(args.output_signature_weights, sep="\t", index=False)
    print(f"Wrote {args.output_signature_weights} ({len(sig_weights)} (tool, signature, sample) rows)")

    original_totals = load_original_catalog(args.original_catalog) if args.original_catalog else None
    sample_activity = sample_activity_summary(contributions, args.min_contribution, original_totals)
    sample_activity.insert(0, "num_muts", args.num_muts)
    sample_activity.insert(0, "relevance", args.relevance)
    sample_activity.insert(0, "cancer_type", args.cancer_type)
    sample_activity.insert(0, "data_type", args.data_type)
    sample_activity.to_csv(args.output_sample_activity, sep="\t", index=False)
    print(f"Wrote {args.output_sample_activity} ({len(sample_activity)} (sample, tool) rows)")


if __name__ == "__main__":
    main()
