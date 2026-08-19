#!/usr/bin/env python3
"""Concatenate per-cell cross-tool agreement TSVs into one table."""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Aggregate cross-tool agreement TSVs")
    parser.add_argument("--input-files", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    frames = [pd.read_csv(f, sep="\t") for f in args.input_files]
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    combined.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output} ({len(combined)} rows from {len(frames)} files)")


if __name__ == "__main__":
    main()
