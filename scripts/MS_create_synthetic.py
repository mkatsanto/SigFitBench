#!/usr/bin/env python3
"""
Generate synthetic mutational catalogs by sampling signature weights from empirical
data and drawing mutation counts via multinomial sampling over trinucleotide contexts.
Supports multiple tool-compatible output formats.
"""

import argparse
import sys
from pathlib import Path
from typing import Union, Tuple, Optional

import numpy as np
import pandas as pd


def validate_signature_contributions(
    contribs: Union[pd.Series, pd.DataFrame],
    sample_name: str,
    epsilon: float = 1e-6,
) -> None:
    """Exit with an error if contributions don't sum to 1 or contain negatives."""
    contrib_sum = contribs.sum() if isinstance(contribs, pd.Series) else contribs.sum().sum()
    if abs(contrib_sum - 1) > epsilon:
        print(f'ERROR: contributions must sum to 1 for {sample_name} '
              f'(off by {contrib_sum - 1:.2e})')
        sys.exit(1)

    min_val = contribs.min() if isinstance(contribs, pd.Series) else contribs.min().min()
    if min_val < 0:
        print(f'ERROR: negative contribution in {sample_name}: {min_val:.2e}')
        sys.exit(1)


def calculate_context_probabilities(
    contribs: pd.Series,
    input_sigs: pd.DataFrame,
) -> pd.Series:
    """Return normalised trinucleotide context probabilities for one sample."""
    context_weights = pd.Series(0.0, index=input_sigs.index)
    for sig_name in input_sigs.columns:
        if sig_name in contribs.index:
            context_weights += contribs[sig_name] * input_sigs[sig_name]

    total = context_weights.sum()
    if total > 0:
        return context_weights / total
    # Fallback: uniform when all weights are zero
    return pd.Series(1.0 / len(context_weights), index=context_weights.index)


def generate_synthetic_mutations(
    rng: np.random.Generator,
    sample_name: str,
    num_muts: int,
    contribs: pd.Series,
    input_sigs: pd.DataFrame,
    epsilon: float = 1e-6,
) -> pd.DataFrame:
    """Sample num_muts mutations for one sample and return per-context counts."""
    validate_signature_contributions(contribs, sample_name, epsilon)
    context_probs = calculate_context_probabilities(contribs, input_sigs)

    counts = np.zeros(len(context_probs), dtype=int)
    np.add.at(counts, rng.choice(len(context_probs), p=context_probs.values, size=num_muts), 1)

    return pd.DataFrame(counts, index=context_probs.index, columns=[sample_name])


def generate_synthetic_catalog_from_contributions(
    rng: np.random.Generator,
    mutations_per_sample: pd.Series,
    signature_contribs: Union[pd.Series, pd.DataFrame],
    signature_profiles: pd.DataFrame,
    epsilon: float = 1e-6,
) -> pd.DataFrame:
    """Generate a full mutational catalog (contexts × samples) from signature contributions."""
    catalogs = []

    if signature_contribs.ndim == 1:
        for i, n in enumerate(mutations_per_sample):
            catalogs.append(generate_synthetic_mutations(
                rng, f'S{i}', n, signature_contribs, signature_profiles, epsilon
            ))
    elif signature_contribs.ndim == 2:
        for i in range(signature_contribs.shape[0]):
            catalogs.append(generate_synthetic_mutations(
                rng, signature_contribs.index[i], mutations_per_sample.iloc[i],
                signature_contribs.iloc[i], signature_profiles, epsilon
            ))
    else:
        raise ValueError(f"Invalid signature contributions shape: {signature_contribs.shape}")

    return pd.concat(catalogs, axis=1)


def save_synthetic_catalogs(
    counts: pd.DataFrame,
    output_prefix: str,
    formats: Optional[list] = None,
) -> None:
    """Save the catalog in one or more tool-compatible formats.

    Formats without index name: deconstructSigs, sigfit, sigLASSO
    Formats with index.name='Type': MutationalPatterns, mmsig, YAPSA, sigminer,
                                     signature_tools, SigsPack
    """
    if formats is None:
        formats = ['deconstructSigs', 'MutationalPatterns']

    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

    no_index_name = {'deconstructSigs', 'sigfit', 'sigLASSO'}
    type_index_name = {'MutationalPatterns', 'mmsig', 'YAPSA', 'sigminer', 'signature_tools', 'SigsPack'}

    for fmt in formats:
        out = f"{output_prefix}_{fmt}.dat"
        if fmt in no_index_name:
            counts.to_csv(out, sep='\t')
        elif fmt in type_index_name:
            c = counts.copy()
            c.index.name = 'Type'
            c.to_csv(out, sep='\t')
        else:
            print(f"Warning: unknown output format '{fmt}', skipping")


def add_out_of_reference_signatures(
    contribs: pd.DataFrame,
    out_of_reference_weights: list,
    out_sigs: list,
    random_seed: int,
) -> pd.DataFrame:
    """Scale down known-signature weights and inject randomly chosen unknown signatures.

    Each sample i gets a sub-RNG seeded with [random_seed, 1, i] so the unknown
    signature assigned to sample i is determined solely by random_seed and i,
    independent of num_samples.
    """
    tot = sum(out_of_reference_weights)
    result = contribs.copy()
    result *= (1 - tot)  # shrink known signatures to make room

    for i, idx in enumerate(result.index):
        sample_rng = np.random.default_rng([random_seed, 1, i])
        chosen = sample_rng.choice(out_sigs, size=len(out_of_reference_weights), replace=False)
        for sig, weight in zip(chosen, out_of_reference_weights):
            if sig not in result.columns:
                result[sig] = 0.0
            result.loc[idx, sig] = weight

    return result


def normalize_empirical_contributions(
    empirical_contribs: pd.DataFrame,
    num_muts: int,
    min_mutations: int = 10,
) -> pd.DataFrame:
    """Rescale empirical counts to num_muts, drop signatures below min_mutations, renormalise.

    Mirrors the original SigFitTest generate_weights_empirical approach.
    """
    contribs = empirical_contribs.copy().astype(float)  # copy avoids dtype-coercion side effects

    for ix in contribs.index:
        total = contribs.loc[ix].sum()
        if total > 0:
            contribs.loc[ix] *= num_muts / total

    contribs[contribs < min_mutations] = 0

    zero_samples = (contribs.sum(axis=1) == 0).sum()
    if zero_samples > 0:
        raise ValueError(
            f"For num_muts={num_muts}, {zero_samples} sample(s) have all signature "
            f"contributions below {min_mutations} mutations. "
            f"Try increasing num_muts or decreasing min_mutations."
        )

    for ix in contribs.index:
        total = contribs.loc[ix].sum()
        if total > 0:
            contribs.loc[ix] /= total

    return contribs


def load_signature_data(
    signature_file: str,
    contributions_file: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load COSMIC signature profiles and empirical contributions; raise if either is missing."""
    sig_path = Path(signature_file)
    contrib_path = Path(contributions_file)

    if not sig_path.exists():
        raise FileNotFoundError(f"Signature profiles file not found: {sig_path}")
    if not contrib_path.exists():
        raise FileNotFoundError(f"Signature contributions file not found: {contrib_path}")

    return pd.read_csv(sig_path, sep='\t', index_col=0), \
           pd.read_csv(contrib_path, sep='\t', index_col=0)


def check_signatures_have_profiles(
    signature_profiles: pd.DataFrame,
    contribs: pd.DataFrame,
    signature_file: str,
    contributions_file: str,
) -> None:
    """Stop if a signature in the weights is missing from the signature catalogue,
    e.g. SBS40 (COSMIC v3.3) when the catalogue is v3.4, which splits it into SBS40a/b/c."""
    used = [s for s in contribs.columns if (contribs[s] > 0).any()]
    missing = sorted(set(used) - set(signature_profiles.columns))
    if missing:
        print(f"WARNING: COSMIC version mismatch. {contributions_file} uses signatures that are "
              f"not in {signature_file}: {', '.join(missing)}", file=sys.stderr)
        print("  Use a catalogue with the same COSMIC version (cosmic_signatures.base_file, "
              "or signature_file for a custom dataset).", file=sys.stderr)
        sys.exit(1)


def create_synthetic_catalog(
    signature_contributions_file: str,
    num_samples: int,
    num_muts: int,
    output_dir: str,
    input_sigs: str,
    epsilon: float = 1e-6,
    random_seed: int = 0,
    output_formats: Optional[list] = None,
    save_weights: bool = True,
    out_of_reference_weights: Optional[list] = None,
    out_sigs: Optional[list] = None,
    cancer_type: Optional[str] = None,
    realization: Optional[int] = None,
) -> None:
    """Create synthetic mutational catalogs and optionally save ground-truth weights."""
    mutation_seed = random_seed + (realization if realization is not None else 0)
    print(f'=== Creating synthetic mutational catalog ===')
    print(f'  Samples: {num_samples} | Mutations: {num_muts}')
    print(f'  Weights seed: {random_seed} (shared across all realisations)')
    print(f'  Mutation seed: {mutation_seed} (random_seed + realization)')
    print(f'  Output: {output_dir}')

    if out_of_reference_weights is None:
        out_of_reference_weights = []
    tot_out_of_reference = sum(out_of_reference_weights)

    if tot_out_of_reference > 0:
        if out_sigs is None or len(out_sigs) == 0:
            raise ValueError("out_sigs must be specified when using out-of-reference weights")
        if not all(0 <= w <= 1 for w in out_of_reference_weights):
            raise ValueError("All out-of-reference weights must be between 0 and 1")
        if tot_out_of_reference >= 1:
            raise ValueError(f"out-of-reference weights sum to {tot_out_of_reference:.3f}, must be < 1")
        print(f'  Out-of-reference: {out_of_reference_weights} → {out_sigs}')

    mutation_rng = np.random.default_rng(mutation_seed)

    print('Loading signature profiles and empirical contributions...')
    signature_profiles, empirical_contribs = load_signature_data(input_sigs, signature_contributions_file)
    print(f'  {signature_profiles.shape[1]} signatures | {empirical_contribs.shape[0]} empirical samples')
    check_signatures_have_profiles(signature_profiles, empirical_contribs, input_sigs,
                                   signature_contributions_file)

    # Each sample i uses sub-RNGs seeded with [random_seed, 0, i] for index selection
    # and [random_seed, 1, i] for unknown-signature assignment, so sample i's weights
    # are fully determined by random_seed and i, independent of num_samples.
    # mutation_rng is offset by realization so each rep produces different mutation counts.
    n_empirical = empirical_contribs.shape[0]
    sample_indices = np.array([
        int(np.random.default_rng([random_seed, 0, i]).integers(0, n_empirical))
        for i in range(num_samples)
    ])
    selected_contribs = empirical_contribs.iloc[sample_indices].reset_index(drop=True)

    normalized_contribs = normalize_empirical_contributions(selected_contribs, num_muts)
    normalized_contribs.index = [f'S{i}' for i in range(num_samples)]

    if tot_out_of_reference > 0:
        normalized_contribs = add_out_of_reference_signatures(
            normalized_contribs, out_of_reference_weights, out_sigs, random_seed
        )
        print(f'  Active signatures after adding unknowns: {(normalized_contribs.sum(axis=0) > 0).sum()}')

    print('Generating synthetic mutations...')
    synthetic_catalog = generate_synthetic_catalog_from_contributions(
        mutation_rng,
        pd.Series(num_muts, index=range(num_samples)),
        normalized_contribs,
        signature_profiles,
        epsilon,
    )

    # Infer cancer type from filename if not supplied.
    # Expected filename format: signature_contributions_{WGS|WES}_{cancer_type}.dat
    if cancer_type is None:
        parts = Path(signature_contributions_file).stem.split('_')
        if len(parts) >= 4 and parts[0] == 'signature' and parts[1] == 'contributions' and parts[2] in ['WGS', 'WES']:
            cancer_type = '_'.join(parts[3:])
        else:
            print(f"Warning: could not extract cancer type from {signature_contributions_file}")

    def _build_path(prefix: str) -> Path:
        if cancer_type:
            if realization is not None:
                return Path(output_dir) / f'{prefix}_{cancer_type}-w_{realization}-{num_muts}'
            return Path(output_dir) / f'{prefix}_{cancer_type}_{num_muts}'
        if realization is not None:
            return Path(output_dir) / f'{prefix}-w_{realization}-{num_muts}'
        return Path(output_dir) / f'{prefix}_{num_muts}'

    output_prefix = _build_path('counts')
    print('Saving synthetic catalog...')
    save_synthetic_catalogs(synthetic_catalog, str(output_prefix), output_formats)

    if save_weights:
        weights_path = _build_path('weights').with_suffix('.dat')
        normalized_contribs.T.to_csv(weights_path, sep='\t')
        print(f'Saved weights → {weights_path}')

    print('✓ Done.')
    for fmt in (output_formats or ['deconstructSigs', 'MutationalPatterns']):
        print(f'  {output_prefix}_{fmt}.dat')
    if save_weights:
        print(f'  {weights_path}')


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Create synthetic mutational catalogs from empirical signature contributions.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python MS_create_synthetic.py \\
    --signature_contributions_file contributions.dat \\
    --num_samples 100 --num_muts 1000

  python MS_create_synthetic.py \\
    --signature_contributions_file contributions.dat \\
    --num_samples 50 --num_muts 5000 --random_seed 42 --realization 2
        """,
    )

    parser.add_argument('--signature_contributions_file', required=True,
                        help='Empirical signature contributions file')
    parser.add_argument('--num_samples', type=int, required=True,
                        help='Number of synthetic samples to generate')
    parser.add_argument('--num_muts', type=int, required=True,
                        help='Mutations per sample')

    parser.add_argument('--output_dir', default='generated_data',
                        help='Output directory (default: generated_data)')
    parser.add_argument('--input_sigs', default='input/COSMIC_v3_SBS_GRCh38.txt',
                        help='COSMIC signature profiles file')
    parser.add_argument('--epsilon', type=float, default=1e-6,
                        help='Numerical tolerance for validation (default: 1e-6)')
    parser.add_argument('--random_seed', type=int, default=0,
                        help='Base seed: used as-is for weight sampling; offset by '
                             '--realization for mutation sampling (default: 0)')
    parser.add_argument('--output_formats', nargs='+',
                        default=['deconstructSigs', 'MutationalPatterns'],
                        help='Output formats (default: deconstructSigs MutationalPatterns)')
    parser.add_argument('--save_weights', action='store_true', default=True,
                        help='Save ground-truth signature weights (default: True)')
    parser.add_argument('--wgs_or_wes', default='WGS',
                        help='WGS or WES (informational, default: WGS)')
    parser.add_argument('--out_of_reference_weights', type=float, nargs='*', default=None,
                        help='Weights for out-of-reference (unknown) signatures')
    parser.add_argument('--out_sigs', type=str, nargs='*', default=None,
                        help='Available out-of-reference signature names')
    parser.add_argument('--cancer_type', default=None,
                        help='Cancer type for output filenames (inferred from contributions '
                             'filename if omitted)')
    parser.add_argument('--realization', type=int, default=None,
                        help='Realisation index; offsets the mutation seed so each rep '
                             'produces different counts with identical weights')

    args = parser.parse_args()

    create_synthetic_catalog(
        signature_contributions_file=args.signature_contributions_file,
        num_samples=args.num_samples,
        num_muts=args.num_muts,
        output_dir=args.output_dir,
        input_sigs=args.input_sigs,
        epsilon=args.epsilon,
        random_seed=args.random_seed,
        output_formats=args.output_formats,
        save_weights=args.save_weights,
        out_of_reference_weights=args.out_of_reference_weights,
        out_sigs=args.out_sigs,
        cancer_type=args.cancer_type,
        realization=args.realization,
    )


if __name__ == '__main__':
    main()
