#!/usr/bin/env python3
"""
Prepare a Real Mutational Catalog for SigFitTest Workflow

Subsets/reformats a catalog of real samples (no ground truth) into the counts file
every tool reads, in the same layout as the simulated samples from
MS_create_synthetic.py, so all downstream tool rules work without modification. Evaluated via cross-tool agreement, since there's no ground
truth.

Inputs
------
--input-catalog PATH      User catalog (TSV/CSV, contexts × samples).
--dataset-name NAME       Dataset name, for the log.
--counts-file PATH        File to write the counts to (folders created if absent).
--order-like PATH         Signature table (e.g. the COSMIC base_file) whose row order the
                          counts follow, the same order as the simulated samples and the
                          tools' reference signatures.
--matrix-order PATH       Instead of --order-like: file with a 96-context order (one
                          context per line, no header).

Author: SigFitTest workflow
"""

import argparse
import csv
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def sanitize_sample_name(name: str) -> str:
    """Replace any character outside [A-Za-z0-9._] with '.'.

    R's read.csv/read.table default to check.names=TRUE, which silently mangles
    column names into valid R identifiers this same way (make.names()) — e.g.
    "Tissue::SampleID" becomes "Tissue..SampleID". Doing it here, once, before
    any tool sees the data, means every tool (R or Python) gets identical sample
    names, instead of R-based tools silently diverging from Python-based tools
    reading the same source file (which broke cross-tool agreement matching).
    """
    return re.sub(r"[^A-Za-z0-9._]", ".", name)


def load_matrix_order(matrix_order_path: str) -> list:
    """Load standard 96-context order from file (one context per line, no header).

    input/mut_matrix_order_*.dat wraps each context in literal CSV double-quotes
    (e.g. `"A[C>A]A"`), matching how every other consumer of this file reads it
    (pd.read_csv, which strips CSV quoting automatically). A plain open()+strip()
    here would keep the quote characters and never match the bare "A[C>A]A"
    context names used elsewhere (e.g. a catalog's row index), silently zeroing
    every context after reindex_to_standard_order — so this uses csv.reader too,
    for the same quote-stripping behavior.
    """
    contexts = []
    with open(matrix_order_path, "r", newline="") as fh:
        for row in csv.reader(fh):
            if row and row[0].strip():
                contexts.append(row[0].strip())
    return contexts


def load_catalog(input_catalog: str) -> pd.DataFrame:
    """Load user catalog with flexible separator detection.

    Returns a DataFrame with contexts as rows and samples as columns.
    """
    # comment='#' skips leading comment lines (e.g. "# PCAWG mutational profiles...")
    # before delimiter sniffing runs — otherwise sep=None can latch onto punctuation
    # inside a comment line (e.g. a comma in "50,000") and misdetect the real
    # delimiter for the whole file.
    catalog = pd.read_csv(input_catalog, sep=None, index_col=0, engine="python", comment="#")
    catalog.index = catalog.index.astype(str).str.strip()
    catalog.columns = catalog.columns.astype(str).str.strip()
    catalog = catalog.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return catalog


def select_tissue(catalog: pd.DataFrame, tissue_prefix: str) -> pd.DataFrame:
    """Subset to columns whose 'Tissue::SampleID' prefix matches tissue_prefix.

    Used when a catalog's columns follow that naming (see get_real_cancer_types
    in workflow/rules/common.smk).
    """
    prefix = f"{tissue_prefix}::"
    selected = [c for c in catalog.columns if c.startswith(prefix)]
    if not selected:
        available = sorted({c.split("::")[0] for c in catalog.columns if "::" in c})
        raise ValueError(f"No samples found for tissue_prefix={tissue_prefix!r}. Available: {available}")
    return catalog[selected]


def reindex_to_standard_order(catalog: pd.DataFrame, context_order: list) -> pd.DataFrame:
    """Reindex catalog rows to the standard context order.

    Missing contexts are filled with 0; unknown contexts trigger a warning.
    """
    unknown = set(catalog.index) - set(context_order)
    if unknown:
        warnings.warn(
            f"The following contexts in the input catalog are NOT in the standard "
            f"matrix order and will be dropped: {sorted(unknown)}"
        )
    reindexed = catalog.reindex(index=context_order, fill_value=0)
    return reindexed


def save_counts(catalog: pd.DataFrame, counts_file: str) -> None:
    """Write the catalog (mutation types x samples, first column "Type"), the same layout
    the simulated samples have, so every tool reads it the same way."""
    Path(counts_file).parent.mkdir(parents=True, exist_ok=True)
    catalog = catalog.copy()
    catalog.index.name = "Type"
    catalog.to_csv(counts_file, sep="\t")
    print(f"  Wrote {counts_file}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a real mutational catalog for the SigFitTest workflow.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--input-catalog", required=True,
        help="Path to the user's catalog (contexts × samples TSV/CSV).",
    )
    parser.add_argument(
        "--dataset-name", required=True,
        help="Dataset name, for the log.",
    )
    parser.add_argument(
        "--counts-file", required=True,
        help="File to write the counts to.",
    )
    order = parser.add_mutually_exclusive_group(required=True)
    order.add_argument(
        "--order-like",
        help="Signature table (e.g. the COSMIC base_file) whose row order the counts follow.",
    )
    order.add_argument(
        "--matrix-order",
        help="File with a 96-context order (one context per line).",
    )
    parser.add_argument(
        "--tissue-prefix", default=None,
        help="If given, subset to columns named 'tissue_prefix::SampleID' (auto-detected "
             "per-tissue splitting for catalogs using that convention — see "
             "get_real_cancer_types in workflow/rules/common.smk). Omit to use every "
             "column in the catalog as one flat dataset.",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load catalog
    # ------------------------------------------------------------------
    print(f"Loading catalog from {args.input_catalog} ...", file=sys.stderr)
    catalog = load_catalog(args.input_catalog)
    print(
        f"  {catalog.shape[0]} contexts × {catalog.shape[1]} samples",
        file=sys.stderr,
    )
    col_sums = catalog.sum(axis=0)
    print(
        f"  Mutations per sample: min={col_sums.min():.0f}, "
        f"median={col_sums.median():.0f}, max={col_sums.max():.0f}",
        file=sys.stderr,
    )

    # ------------------------------------------------------------------
    # Subset to one tissue if the catalog uses "Tissue::SampleID" columns
    # ------------------------------------------------------------------
    if args.tissue_prefix:
        catalog = select_tissue(catalog, args.tissue_prefix)
        print(
            f"  Selected {catalog.shape[1]} samples for tissue_prefix={args.tissue_prefix!r}",
            file=sys.stderr,
        )

    # ------------------------------------------------------------------
    # Reindex to standard context order
    # ------------------------------------------------------------------
    if args.order_like:
        # Same order as the simulated samples and the reference signatures, so every tool
        # gets real samples in the order it gets in the benchmark
        print(f"Using the mutation order of {args.order_like} ...", file=sys.stderr)
        context_order = pd.read_csv(args.order_like, sep="\t", index_col=0).index.astype(str).str.strip().tolist()
    else:
        print(f"Loading matrix order from {args.matrix_order} ...", file=sys.stderr)
        context_order = load_matrix_order(args.matrix_order)
    catalog = reindex_to_standard_order(catalog, context_order)

    # ------------------------------------------------------------------
    # Sanitize sample names (after tissue selection, which needs the raw
    # "Tissue::SampleID" columns; before saving, so every tool gets the same
    # already-safe names)
    # ------------------------------------------------------------------
    catalog.columns = [sanitize_sample_name(c) for c in catalog.columns]

    # ------------------------------------------------------------------
    # Save the counts file
    # ------------------------------------------------------------------
    save_counts(catalog, args.counts_file)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(
        f"\nReal catalog preparation complete: "
        f"dataset={args.dataset_name}, "
        f"samples={catalog.shape[1]}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
