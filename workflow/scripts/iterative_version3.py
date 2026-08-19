#!/usr/bin/env python3
"""
Iterative signature fitting method (Version 3): z-score based unassigned mutation
removal (Phase 1, as in Version 1), followed by an SPA-style coarse pre-filter and
greedy add/remove of known signatures (Phase 2, adapted from
extended_removal_addition.py) instead of a single backward-elimination pass.

Usage:
    python iterative_version3.py -i <input_data> -w <cosmic_signatures> -m <matrix_order> -o <output_file>
"""

import argparse
import numpy as np
import pandas as pd
import scipy as sp
import os
from pathlib import Path


def recalculate_TAE(df_pred, df_true):
    """
    Calculate the Total Absolute Error (TAE) between predicted and true values.
    """
    muts_true = df_true.sum()
    df_pred1 = df_pred.reindex(index=df_true.index, fill_value=0)
    df_true1 = df_pred1.reindex(index=df_true.index, fill_value=0)
    TAE = df_true1.sub(df_pred1, axis=0).abs().sum() / (2 * (muts_true))
    TAE = pd.DataFrame(TAE)
    return TAE


def cosine_similarity(vec1, vec2):
    """Calculate cosine similarity between two vectors."""
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    if norm1 == 0 or norm2 == 0:
        return 0.0
    return np.dot(vec1, vec2) / (norm1 * norm2)


def normalized_l2_distance(vec1, vec2):
    """Calculate normalized L2 distance between two vectors."""
    norm1 = np.linalg.norm(vec1)
    if norm1 == 0:
        return 0.0
    return np.linalg.norm(vec1 - vec2) / norm1


def _fit_score(x, W, metric="l2"):
    """NNLS-fit (score, h) for W against x; higher score is always better.
    Adapted from extended_removal_addition.py's _fit_score."""
    h, _ = sp.optimize.nnls(W, x)
    pred = W @ h
    if metric == "cosine":
        return cosine_similarity(x, pred), h
    if metric == "l2":
        return -normalized_l2_distance(x, pred), h
    raise ValueError(f"Unknown metric {metric!r}; expected 'cosine' or 'l2'")


def _backward_eliminate(x, W_curr, active_indices, thresh_backward, max_iter, background=frozenset(), metric="l2"):
    """Drop the least-damaging signature each round until every drop exceeds thresh_backward.
    Indices in `background` are never removed. Adapted from extended_removal_addition.py."""
    for _ in range(max_iter):
        if len(active_indices) <= 1:
            break

        score_current, _ = _fit_score(x, W_curr, metric)

        best_drop, best_keep = np.inf, None
        for position, signature in enumerate(active_indices):
            if signature in background:
                continue
            keep = [j for j in range(len(active_indices)) if j != position]
            drop = score_current - _fit_score(x, W_curr[:, keep], metric)[0]
            if drop < best_drop:
                best_drop, best_keep = drop, keep

        if best_keep is None or best_drop >= thresh_backward:
            break

        active_indices = [active_indices[j] for j in best_keep]
        W_curr = W_curr[:, best_keep]

    return W_curr, active_indices


def _initial_prune(x, W_known, background, initial_remove_penalty, max_iter, metric):
    """SPA-style coarse pre-filter: one lenient backward-elimination pass (looser than
    thresh_backward) to shrink the pool before the stricter add/remove loop.
    Adapted from extended_removal_addition.py."""
    if initial_remove_penalty <= 0:
        return list(range(W_known.shape[1]))
    _, active_indices = _backward_eliminate(
        x, W_known, list(range(W_known.shape[1])), initial_remove_penalty, max_iter,
        background=background, metric=metric,
    )
    return active_indices


def _add_remove_known_signatures(
    x, W_known, background, thresh_forward, thresh_backward, metric, max_iter, initial_active=None,
):
    """SPA-style greedy layered selection: each layer adds whichever candidate scores best
    after also being backward-eliminated; stops once a layer fails to improve by more than
    thresh_forward. Adapted from extended_removal_addition.py."""
    n_known = W_known.shape[1]
    active_indices = sorted(initial_active if initial_active is not None else background)
    current_score = (
        _fit_score(x, W_known[:, active_indices], metric)[0] if active_indices else -np.inf
    )

    for _ in range(max_iter):
        candidates = [i for i in range(n_known) if i not in active_indices]
        if not candidates:
            break

        layer_best_score, layer_best_active = -np.inf, None
        for candidate in candidates:
            trial_active = sorted(set(active_indices) | {candidate})
            W_pruned, pruned_active = _backward_eliminate(
                x, W_known[:, trial_active], trial_active, thresh_backward, max_iter,
                background=background, metric=metric,
            )
            score = _fit_score(x, W_pruned, metric)[0]
            if score > layer_best_score:
                layer_best_score, layer_best_active = score, pruned_active

        if layer_best_active is None or (layer_best_score - current_score) < thresh_forward:
            break

        active_indices, current_score = layer_best_active, layer_best_score

    return W_known[:, active_indices], active_indices


def leave_out_unassigned(x, W, true_weights, sig_names, sample_name=None, zscore_threshold=3,
                         thresh_backward=0.01, thresh_forward=0.01, initial_remove_penalty=0.05,
                         background_signatures=("SBS1", "SBS5"), max_iter=3000,
                         prune_signatures_flag=True, output_dir=None):
    """
    Iteratively remove significant unassigned mutations, then apply an SPA-style coarse
    pre-filter followed by greedy add/remove of known signatures (L2 distance based).

    Parameters:
    -----------
    x : np.ndarray
        Observed mutation counts (channels)
    W : np.ndarray
        Signature matrix (channels x signatures)
    true_weights : np.ndarray
        True weights for TAE calculation
    sig_names : list
        Signature names matching the columns of W, used to protect background_signatures
        from removal.
    sample_name : str
        Name of the sample (for logging)
    zscore_threshold : float
        Z-score threshold for unassigned mutation removal
    thresh_backward : float
        L2 distance increase threshold for backward elimination (default: 0.01)
    thresh_forward : float
        Minimum score improvement required for a greedy add layer to continue (default: 0.01)
    initial_remove_penalty : float
        Threshold for the coarse pre-filter pass before the add/remove loop (default: 0.05)
    background_signatures : tuple
        Signature names that are never removed (default: SBS1, SBS5)
    max_iter : int
        Maximum iterations for unassigned removal
    prune_signatures_flag : bool
        Whether to apply the SPA-style add/remove step
    output_dir : str or Path
        Directory to save convergence file

    Returns:
    --------
    np.ndarray: Final weights
    """
    x = np.array([0 if v is None else v for v in x], dtype=np.float64)
    n_sigs = W.shape[1]
    background = frozenset(
        sig_names.index(s) for s in background_signatures if s in sig_names
    )

    convergence_data = []
    unassigned_list = []
    x_fit = x.copy()

    print("\n=== PHASE 1: Iterative Unassigned Removal ===")
    for i in range(max_iter):
        weights, _ = sp.optimize.nnls(W, x_fit)

        predicted_counts = W @ weights
        unassigned = np.maximum(x_fit - predicted_counts, 0)
        assigned = np.minimum(predicted_counts, x_fit)
        overfitted = np.maximum(predicted_counts - x_fit, 0)

        overfitted_sum = np.sum(overfitted)
        unassigned_sum = np.sum(unassigned)
        assigned_sum = np.sum(assigned)
        x_fit_sum = np.sum(x_fit)
        predicted_sum = np.sum(predicted_counts)
        total_unassigned = np.sum(unassigned_list)

        TAE = recalculate_TAE(pd.DataFrame(weights), pd.DataFrame(true_weights))
        tae_value = TAE.values[0][0]

        n_active = np.sum(weights > 5)
        l2_dist = normalized_l2_distance(x_fit, predicted_counts)
        cos_sim = cosine_similarity(x_fit, predicted_counts)

        convergence_data.append({
            'iteration': i + 1,
            'phase': 'removal',
            'x_fit_sum': x_fit_sum,
            'predicted_sum': predicted_sum,
            'assigned': assigned_sum,
            'unassigned_current': unassigned_sum,
            'unassigned_total': total_unassigned,
            'overfitted': overfitted_sum,
            'TAE': tae_value,
            'l2_distance': l2_dist,
            'cosine_similarity': cos_sim,
            'n_active_signatures': n_active,
            'zscore_threshold': zscore_threshold,
            'thresh_backward': thresh_backward
        })

        if np.sum(unassigned) > 0:
            mean = np.mean(unassigned)
            std = np.std(unassigned)

            if std == 0:
                print("Standard deviation is zero. No significant unassigned mutations to remove.")
                break

            threshold = mean + zscore_threshold * std
            removed = np.zeros_like(unassigned)
            above = unassigned > threshold
            removed[above] = unassigned[above] - threshold

            if removed[above].sum() == 0:
                print("No significant unassigned mutations left to remove.")
                break

            x_fit -= removed
            x_fit = np.maximum(x_fit, 0)
            unassigned_list.append(np.sum(removed))
        else:
            unassigned_list.append(0)
            print("No unassigned mutations left.")
            break

        print(f"Iter {i+1}: x_fit_sum={x_fit_sum:.2f}, predicted_sum={predicted_sum:.2f}, "
              f"assigned={assigned_sum:.2f}, unassigned_total={total_unassigned:.2f}, "
              f"overfitted={overfitted_sum:.2f}, TAE={tae_value:.4f}, l2_dist={l2_dist:.4f}, "
              f"cos_sim={cos_sim:.4f}, active_sigs={n_active}")

        if overfitted_sum < 100:
            print("Overfitting below threshold. Stopping removal phase.")
            break

    print("\n=== PHASE 2: SPA-Style Add/Remove of Known Signatures (L2 Version) ===")
    if prune_signatures_flag:
        initial_active = _initial_prune(
            x_fit, W, background, initial_remove_penalty, max_iter=100, metric="l2"
        )
        W_final, active_indices = _add_remove_known_signatures(
            x_fit, W, background, thresh_forward, thresh_backward, "l2", max_iter=100,
            initial_active=initial_active,
        )
        h_active, _ = sp.optimize.nnls(W_final, x_fit)
        weights = np.zeros(n_sigs)
        weights[active_indices] = h_active

        predicted_counts = W @ weights
        unassigned = np.maximum(x_fit - predicted_counts, 0)
        assigned = np.minimum(predicted_counts, x_fit)
        overfitted = np.maximum(predicted_counts - x_fit, 0)

        overfitted_sum = np.sum(overfitted)
        unassigned_sum = np.sum(unassigned)
        assigned_sum = np.sum(assigned)
        predicted_sum = np.sum(predicted_counts)
        total_unassigned = np.sum(unassigned_list)

        TAE = recalculate_TAE(pd.DataFrame(weights), pd.DataFrame(true_weights))
        tae_value = TAE.values[0][0]
        n_active = np.sum(weights > 5)
        l2_dist = normalized_l2_distance(x_fit, predicted_counts)
        cos_sim = cosine_similarity(x_fit, predicted_counts)

        convergence_data.append({
            'iteration': len(convergence_data) + 1,
            'phase': 'add_remove',
            'x_fit_sum': x_fit_sum,
            'predicted_sum': predicted_sum,
            'assigned': assigned_sum,
            'unassigned_current': unassigned_sum,
            'unassigned_total': total_unassigned,
            'overfitted': overfitted_sum,
            'TAE': tae_value,
            'l2_distance': l2_dist,
            'cosine_similarity': cos_sim,
            'n_active_signatures': n_active,
            'zscore_threshold': zscore_threshold,
            'thresh_backward': thresh_backward
        })

        print(f"After SPA-style add/remove (L2): predicted_sum={predicted_sum:.2f}, assigned={assigned_sum:.2f}, "
              f"unassigned_total={total_unassigned:.2f}, overfitted={overfitted_sum:.2f}, "
              f"TAE={tae_value:.4f}, l2_dist={l2_dist:.4f}, cos_sim={cos_sim:.4f}, active_sigs={n_active}")
    else:
        print("\n=== PHASE 2: Final Fit (No Add/Remove) ===")
        weights, _ = sp.optimize.nnls(W, x_fit)

        TAE = recalculate_TAE(pd.DataFrame(weights), pd.DataFrame(true_weights))
        tae_value = TAE.values[0][0]
        n_active = np.sum(weights > 5)

        print(f"Final fit: TAE={tae_value:.4f}, active_sigs={n_active}")

    print(f"\nFinal: {n_active} active signatures, TAE={tae_value:.4f}")

    convergence_df = pd.DataFrame(convergence_data)

    if output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        sample_str = f"_{sample_name}" if sample_name else ""
        filepath = output_path / f"convergence{sample_str}_sigprofiler_style_l2.csv"

        convergence_df.to_csv(filepath, index=False)
        print(f"Convergence data saved to: {filepath}")

    return weights


def main():
    parser = argparse.ArgumentParser(description='Iterative signature fitting method (Version 3 with SPA-style add/remove)')
    parser.add_argument('-i', '--input', required=True, help='Input mutation matrix')
    parser.add_argument('-w', '--weights', required=True, help='COSMIC signature matrix')
    parser.add_argument('-m', '--matrix-order', required=True, help='Matrix order file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('--zscore-threshold', type=float, default=1.0,
                        help='Z-score threshold for unassigned mutation removal')
    parser.add_argument('--thresh-backward', type=float, default=0.01,
                        help='L2 distance increase threshold for backward elimination')
    parser.add_argument('--thresh-forward', type=float, default=0.01,
                        help='Minimum score improvement required for a greedy add layer to continue')
    parser.add_argument('--initial-remove-penalty', type=float, default=0.05,
                        help='Threshold for the coarse pre-filter pass before the add/remove loop')
    parser.add_argument('--max-iter', type=int, default=5000,
                        help='Maximum iterations')
    args = parser.parse_args()

    X = pd.read_csv(args.input, index_col=0, sep='\t')
    W = pd.read_csv(args.weights, sep='\t', index_col=0)
    MP_index = pd.read_csv(args.matrix_order, header=None).squeeze()

    W = W.reindex(MP_index)
    X = X.reindex(MP_index)

    sig_names = W.columns.tolist()

    outdir = os.path.dirname(args.output) or '.'
    os.makedirs(outdir, exist_ok=True)

    H_results = []
    for sample_name in X.columns:
        print(f"Processing sample: {sample_name}")

        h = leave_out_unassigned(
            x=X[sample_name].values,
            W=W.values,
            true_weights=np.zeros(W.shape[1]),
            sig_names=sig_names,
            sample_name=sample_name,
            zscore_threshold=args.zscore_threshold,
            thresh_backward=args.thresh_backward,
            thresh_forward=args.thresh_forward,
            initial_remove_penalty=args.initial_remove_penalty,
            max_iter=args.max_iter,
            prune_signatures_flag=True,
            output_dir=outdir
        )
        H_results.append(h)

    H_results = np.array(H_results).T
    all_sig_names = W.columns.tolist()
    H = pd.DataFrame(H_results, index=all_sig_names, columns=X.columns)

    H.to_csv(args.output)
    print(f"Results saved to {args.output}")


if __name__ == '__main__':
    main()
