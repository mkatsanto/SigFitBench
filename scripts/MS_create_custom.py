#!/usr/bin/env python3
"""
Simulate mutation counts from your own table of signature weights.

Each row of the table is one sample and keeps its name. Its weights give the share of
mutations from each signature, and the sample gets num_muts mutations drawn from that mix.
Repeated runs use the same weights and differ only by chance.
"""

import argparse

import numpy as np
import pandas as pd

from MS_create_synthetic import (
    check_signatures_have_profiles,
    generate_synthetic_catalog_from_contributions,
    load_signature_data,
    save_run_files,
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
    counts_file: str,
    known_answer_file: str,
    input_sigs: str,
    dataset: str,
    realization: int,
    random_seed: int = 0,
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

    save_run_files(catalog, weights, counts_file, known_answer_file)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--weights_file', required=True,
                        help='Table of signature weights (rows = samples, columns = signatures)')
    parser.add_argument('--num_muts', type=int, required=True, help='Mutations per sample')
    parser.add_argument('--counts_file', required=True, help='File for the mutation counts')
    parser.add_argument('--known_answer_file', required=True,
                        help='File for the true signature weights of each sample')
    parser.add_argument('--input_sigs', required=True, help='COSMIC signature profiles file')
    parser.add_argument('--dataset', required=True, help='Dataset name, for the log')
    parser.add_argument('--realization', type=int, required=True, help='Repeat number')
    parser.add_argument('--random_seed', type=int, default=0, help='Random seed')
    args = parser.parse_args()

    create_custom_catalog(
        weights_file=args.weights_file,
        num_muts=args.num_muts,
        counts_file=args.counts_file,
        known_answer_file=args.known_answer_file,
        input_sigs=args.input_sigs,
        dataset=args.dataset,
        realization=args.realization,
        random_seed=args.random_seed,
    )


if __name__ == '__main__':
    main()
