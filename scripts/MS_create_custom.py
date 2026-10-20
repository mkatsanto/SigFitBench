#!/usr/bin/env python3
"""
Simulate mutation counts from your own table of signature weights.

Each row of the table is one sample and keeps its name. Its weights give the share of
mutations from each signature, and the sample gets num_muts mutations drawn from that mix.
Repeated runs use the same weights and differ only by chance.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from MS_create_synthetic import (
    check_signatures_have_profiles,
    generate_synthetic_catalog_from_contributions,
    load_signature_data,
    save_synthetic_catalogs,
)


def table_weights(table: pd.DataFrame, table_file: str) -> pd.DataFrame:
    """Make each sample's weights add up to 1. Stops if a weight is negative or a sample has none."""
    table = table.astype(float)
    if (table < 0).any().any():
        raise ValueError(f"{table_file} contains negative weights")
    totals = table.sum(axis=1)
    if (totals == 0).any():
        raise ValueError(f"{table_file} has samples with all weights 0: "
                         f"{', '.join(map(str, totals.index[totals == 0]))}")
    return table.div(totals, axis=0)


def create_custom_catalog(
    weights_file: str,
    num_muts: int,
    output_dir: str,
    input_sigs: str,
    dataset: str,
    realization: int,
    random_seed: int = 0,
    output_formats: list = None,
) -> None:
    print(f'=== Creating custom mutational catalog: {dataset} ===')
    signature_profiles, table = load_signature_data(input_sigs, weights_file)
    check_signatures_have_profiles(signature_profiles, table, input_sigs, weights_file)
    weights = table_weights(table, weights_file)
    print(f'  Samples: {weights.shape[0]} (rows of {weights_file}) | Mutations: {num_muts}')

    # Each repeat uses its own seed, so it draws different mutations from the same weights
    rng = np.random.default_rng(random_seed + realization)
    catalog = generate_synthetic_catalog_from_contributions(
        rng,
        pd.Series(num_muts, index=range(weights.shape[0])),
        weights,
        signature_profiles,
    )

    prefix = Path(output_dir) / f'counts_{dataset}-w_{realization}-{num_muts}'
    save_synthetic_catalogs(catalog, str(prefix), output_formats)
    weights_path = Path(output_dir) / f'weights_{dataset}-w_{realization}-{num_muts}.dat'
    weights.T.to_csv(weights_path, sep='\t')
    print(f'Saved {prefix}_*.dat and {weights_path}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--weights_file', required=True,
                        help='Table of signature weights (rows = samples, columns = signatures)')
    parser.add_argument('--num_muts', type=int, required=True, help='Mutations per sample')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--input_sigs', required=True, help='COSMIC signature profiles file')
    parser.add_argument('--dataset', required=True, help='Dataset name used in output filenames')
    parser.add_argument('--realization', type=int, required=True, help='Repeat number')
    parser.add_argument('--random_seed', type=int, default=0, help='Random seed')
    parser.add_argument('--output_formats', nargs='+',
                        default=['deconstructSigs', 'MutationalPatterns'], help='Output formats')
    args = parser.parse_args()

    create_custom_catalog(
        weights_file=args.weights_file,
        num_muts=args.num_muts,
        output_dir=args.output_dir,
        input_sigs=args.input_sigs,
        dataset=args.dataset,
        realization=args.realization,
        random_seed=args.random_seed,
        output_formats=args.output_formats,
    )


if __name__ == '__main__':
    main()
