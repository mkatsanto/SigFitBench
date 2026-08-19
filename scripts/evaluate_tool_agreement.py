#!/usr/bin/env python3
"""
Cross-tool agreement for one (data_type, cancer_type, relevance, num_muts) cell.

Used where no ground truth exists (realsim, external without a ground_truth
file): instead of comparing each tool's fit to a known answer, tools are
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
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def load_contribution(path: str) -> pd.DataFrame:
    """Load a {tool}_contribution.dat file: signatures x samples."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


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


if __name__ == "__main__":
    main()
