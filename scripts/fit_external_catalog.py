#!/usr/bin/env python3
"""
Prepare External Mutational Catalog for SigFitTest Workflow

Subsets/reformats a user-supplied mutation-channel catalog (no ground truth
required) into the standard per-tool formats used by the workflow. Outputs are
drop-in replacements for the files produced by MS_create_synthetic.py /
prepare_real_samples.py, so all downstream tool rules work without modification.
Evaluated via cross-tool agreement, like realsim, since there's no ground truth.

Inputs
------
--input-catalog PATH      User catalog (TSV/CSV, contexts × samples).
--dataset-name NAME       Dataset name; used as the cancer_type wildcard.
--output-dir DIR          Directory to write output files into (created if absent).
--output-formats FMT...   List of format names required (e.g. MutationalPatterns
                          deconstructSigs sigfit …).
--matrix-order PATH       File with the standard 96-context order (one context per
                          line, no header).
--num-muts-label INT      Integer label used in output filenames (= wildcards.num_muts
                          from the calling Snakemake rule).

Outputs (written to --output-dir)
----------------------------------
counts_{dataset_name}-w_0-{num_muts_label}_{format}.dat   for each format

Author: SigFitTest workflow
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Format constants — same groupings as MS_create_synthetic.py
# ---------------------------------------------------------------------------

# Formats where the index has NO special name (plain CSV/TSV)
NO_INDEX_NAME_FORMATS = {"deconstructSigs", "sigfit", "sigLASSO"}

# Formats where index.name must be set to 'Type'
TYPE_INDEX_NAME_FORMATS = {
    "MutationalPatterns", "mmsig", "YAPSA", "sigminer",
    "signature_tools", "SigsPack"
}


def load_matrix_order(matrix_order_path: str) -> list:
    """Load standard 96-context order from file (one context per line, no header)."""
    contexts = []
    with open(matrix_order_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line:
                contexts.append(line)
    return contexts


def load_catalog(input_catalog: str) -> pd.DataFrame:
    """Load user catalog with flexible separator detection.

    Returns a DataFrame with contexts as rows and samples as columns.
    """
    catalog = pd.read_csv(input_catalog, sep=None, index_col=0, engine="python")
    catalog.index = catalog.index.astype(str).str.strip()
    catalog.columns = catalog.columns.astype(str).str.strip()
    catalog = catalog.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return catalog


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


def save_catalog_formats(
    catalog: pd.DataFrame,
    output_dir: Path,
    dataset_name: str,
    num_muts_label: str,
    output_formats: list,
) -> None:
    """Save catalog in each requested tool format (same logic as MS_create_synthetic.py)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"counts_{dataset_name}-w_0-{num_muts_label}"

    for fmt in output_formats:
        out_path = output_dir / f"{prefix}_{fmt}.dat"
        if fmt in NO_INDEX_NAME_FORMATS:
            catalog.to_csv(out_path, sep="\t")
        elif fmt in TYPE_INDEX_NAME_FORMATS:
            cat_copy = catalog.copy()
            cat_copy.index.name = "Type"
            cat_copy.to_csv(out_path, sep="\t")
        else:
            warnings.warn(f"Unknown format '{fmt}'; writing as plain TSV.")
            catalog.to_csv(out_path, sep="\t")
        print(f"  Wrote {out_path}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare an external mutational catalog for the SigFitTest workflow.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--input-catalog", required=True,
        help="Path to the user's catalog (contexts × samples TSV/CSV).",
    )
    parser.add_argument(
        "--dataset-name", required=True,
        help="Dataset name; used as cancer_type wildcard in filenames.",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory to write output files into.",
    )
    parser.add_argument(
        "--output-formats", nargs="+", required=True,
        help="List of format names to generate (e.g. MutationalPatterns deconstructSigs).",
    )
    parser.add_argument(
        "--matrix-order", required=True,
        help="File with standard 96-context order (one context per line).",
    )
    parser.add_argument(
        "--num-muts-label", required=True,
        help="Integer label used in output filenames (= wildcards.num_muts).",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)

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
    # Reindex to standard context order
    # ------------------------------------------------------------------
    print(f"Loading matrix order from {args.matrix_order} ...", file=sys.stderr)
    context_order = load_matrix_order(args.matrix_order)
    catalog = reindex_to_standard_order(catalog, context_order)

    # ------------------------------------------------------------------
    # Save per-format count files
    # ------------------------------------------------------------------
    print("Saving catalog in tool formats ...", file=sys.stderr)
    save_catalog_formats(
        catalog, output_dir, args.dataset_name,
        args.num_muts_label, args.output_formats
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(
        f"\nExternal catalog preparation complete: "
        f"dataset={args.dataset_name}, "
        f"num_muts_label={args.num_muts_label}, "
        f"samples={catalog.shape[1]}, "
        f"formats={args.output_formats}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
