#!/usr/bin/env python3
"""
Extended signature fitting method with L1 penalty for novel signature discovery.

Usage:
    python extended.py -i <input_data> -w <cosmic_signatures> -m <matrix_order> -o <output_file>
"""

import argparse
import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import minimize, curve_fit
import os


def recalculate_TAE(df_pred, df_true, known_sig_names=None) -> float:
    """Total Absolute Error between predicted and true values; nan if df_true is None.

    known_sig_names, if given, restricts df_pred to those signatures first, so dummy
    Unknown_i predictions can't coincidentally name-match unrelated ones in df_true.
    """
    if df_true is None:
        return float("nan")
    if known_sig_names is not None:
        df_pred = df_pred.reindex(known_sig_names, fill_value=0)
    domain = df_pred.index.union(df_true.index)
    df_pred = df_pred.reindex(domain, fill_value=0)
    df_true = df_true.reindex(domain, fill_value=0)
    tae = df_true.sub(df_pred, axis=0).abs().sum() / (2 * df_true.sum())
    return float(tae.iloc[0])


def create_enhanced_unknown_signatures(W, method='identity', n_unknown=None):
    """Create unknown signature columns instead of an identity matrix."""
    n_mut_types = W.shape[0]
    if n_unknown is None:
        n_unknown = min(n_mut_types, 20)

    if method == 'identity':
        lowest_channels = np.argsort(W.max(axis=1))[:n_unknown]
        unknown_sigs = np.zeros((n_mut_types, n_unknown))
        for col, channel in enumerate(lowest_channels):
            unknown_sigs[channel, col] = 1.0
        return unknown_sigs

    # method == 'orthogonal'
    U, s, _ = np.linalg.svd(W, full_matrices=True)
    null_space = U[:, len(s):]
    n_to_use = min(n_unknown, null_space.shape[1])
    if n_to_use == 0:
        unknown_sigs = np.random.random((n_mut_types, n_unknown))
        return unknown_sigs / unknown_sigs.sum(axis=0)

    candidates = np.abs(null_space[:, :n_to_use]).T
    valid = [col / col.sum() for col in candidates if col.sum() > 0]
    return np.column_stack(valid) if valid else np.ones((n_mut_types, 1)) / n_mut_types


def advanced_signature_filtering(
    W_extended, original_n_sigs, correlation_threshold=0.90, redundancy_threshold=0.99
):
    """Remove dummy signatures that are highly correlated with known ones or with each other."""
    correlation_matrix = np.corrcoef(W_extended.T)
    n_sigs = W_extended.shape[1]
    to_remove = set()

    for cosmic_idx in range(original_n_sigs):
        for dummy_idx in range(original_n_sigs, n_sigs):
            if dummy_idx not in to_remove and abs(correlation_matrix[cosmic_idx, dummy_idx]) > correlation_threshold:
                to_remove.add(dummy_idx)

    dummy_indices = [i for i in range(original_n_sigs, n_sigs) if i not in to_remove]
    for i, idx1 in enumerate(dummy_indices):
        for idx2 in dummy_indices[i + 1:]:
            if idx2 not in to_remove and abs(correlation_matrix[idx1, idx2]) > redundancy_threshold:
                to_remove.add(idx2)

    kept_indices = [i for i in range(n_sigs) if i not in to_remove]
    return W_extended[:, kept_indices]


def create_diverse_initialization(W_extended, x, start_idx, original_n_sigs):
    """Return an initialisation vector for SLSQP using one of four strategies."""
    n_sigs = W_extended.shape[1]
    n_unknown = n_sigs - original_n_sigs
    W_known = W_extended[:, :original_n_sigs]
    total_muts = np.sum(x)
    h_init = np.zeros(n_sigs)

    if start_idx == 0:
        h_init[:original_n_sigs], _ = sp.optimize.nnls(W_known, x)

    elif start_idx == 1:
        try:
            h_nnls, _ = sp.optimize.nnls(W_known, x)
            h_init[:original_n_sigs] = h_nnls * np.random.uniform(0.8, 1.2, original_n_sigs)
            if n_unknown > 0:
                residual_strength = np.sum(np.abs(x - W_known @ h_nnls))
                h_init[original_n_sigs:] = residual_strength * 0.1 / n_unknown
        except Exception:
            h_init[:original_n_sigs] = total_muts * 0.95 / original_n_sigs
            h_init[original_n_sigs:] = total_muts * 0.05 / n_unknown

    elif start_idx == 2:
        np.random.seed(42 + start_idx)
        n_active = min(10, original_n_sigs)
        active_indices = np.argsort(np.var(W_known, axis=0))[-n_active:]
        h_init[active_indices] = np.random.uniform(0.01, 0.1, n_active) * total_muts / n_active
        if n_unknown > 0:
            h_init[original_n_sigs:] = 1e-3

    else:  # start_idx >= 3: NNLS solution perturbed with Gaussian noise
        try:
            h_nnls, _ = sp.optimize.nnls(W_known, x)
            noise = np.random.normal(0, 0.1 * h_nnls, original_n_sigs)
            h_init[:original_n_sigs] = np.maximum(h_nnls + noise, 1e-6)
            if n_unknown > 0:
                h_init[original_n_sigs:] = np.random.uniform(1e-4, 1e-3, n_unknown) * total_muts
        except Exception:
            h_init = np.random.uniform(0.001, 0.05, n_sigs) * total_muts / n_sigs

    current_sum = np.sum(h_init)
    h_init = h_init * (total_muts / current_sum) if current_sum > 0 else np.full(n_sigs, total_muts / n_sigs)
    return np.maximum(h_init, 1e-10)


def solution_based_initialization(true_weights, W_extended, x, original_n_sigs):
    """Create initialization based on known weights.
    Used only in simulations to study the effect of initialization on objective convergence."""
    n_sigs = W_extended.shape[1]
    total_muts = np.sum(x)
    h_init = np.zeros(n_sigs)

    try:
        h_init[:original_n_sigs] = true_weights
    except Exception:
        h_init[:original_n_sigs] = 1e-3

    total = np.sum(h_init)
    if total > 0:
        return h_init * (total_muts / total)
    return np.full(n_sigs, total_muts / n_sigs)


def _cosine_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """Cosine similarity between two vectors; returns 0.0 if either is zero."""
    n1, n2 = np.linalg.norm(vec1), np.linalg.norm(vec2)
    if n1 == 0.0 or n2 == 0.0:
        return 0.0
    return float(np.dot(vec1, vec2) / (n1 * n2))


def _sigmoid(logp, L, U, x0, k):
    """Logistic curve in log10(penalty): L = lower asymptote, U = upper asymptote,
    x0 = inflection point (log10 scale), k = steepness."""
    return L + (U - L) / (1 + np.exp(-k * (logp - x0)))


def _sigmoid_d1(logp, L, U, x0, k):
    """dy/dlogp for _sigmoid."""
    e = np.exp(-k * (logp - x0))
    return (U - L) * k * e / (1 + e) ** 2


def _sigmoid_d2(logp, L, U, x0, k):
    """d2y/dlogp2 for _sigmoid."""
    e = np.exp(-k * (logp - x0))
    return (U - L) * k ** 2 * e * (e - 1) / (1 + e) ** 3


def _robust_sigmoid_fit(logp: np.ndarray, y: np.ndarray, n_iter: int = 4, mad_thresh: float = 3.0):
    """Fit a sigmoid to (logp, y), iteratively dropping points whose residual from the
    current fit is a large outlier (> mad_thresh scaled MAD), so a few noisy grid points
    don't distort where the curve is judged to bend. The inflection x0 is bounded to the
    observed penalty range so a non-monotonic curve can't push the fit into a wild
    extrapolation. Returns (params, inlier_mask); (None, all-True mask) on failure.
    """
    mask = np.ones_like(logp, dtype=bool)
    params = None
    y_range = y.max() - y.min() + 1e-9
    for _ in range(n_iter):
        if mask.sum() < 4:
            break
        lo, hi = logp[mask].min(), logp[mask].max()
        p0 = [y[mask].min(), y[mask].max(), np.median(logp[mask]), 1.0]
        bounds = (
            [y.min() - 0.5 * y_range, y.min(), lo, 1e-3],
            [y.max(), y.max() + 0.5 * y_range, hi, 50.0],
        )
        try:
            new_params, _ = curve_fit(_sigmoid, logp[mask], y[mask], p0=p0, bounds=bounds, maxfev=20000)
        except (RuntimeError, ValueError):
            break
        params = new_params
        resid = y - _sigmoid(logp, *params)
        mad = np.median(np.abs(resid[mask] - np.median(resid[mask]))) + 1e-12
        new_mask = np.abs(resid) <= mad_thresh * 1.4826 * mad
        if new_mask.sum() == mask.sum():
            mask = new_mask
            break
        mask = new_mask
    return params, mask


def _knee_penalty(grid: np.ndarray, l2_errors: list, k_min_reliable: float = 1e-2, rise_min_frac: float = 0.1):
    """Locate the penalty at the *low* knee — the max-curvature point on the rising part
    of the l2_error-vs-log10(penalty) curve, before its inflection — from a robustly-fit
    sigmoid's first and second derivatives:
        curvature(logp) = |f''(logp)| / (1 + f'(logp)^2)^1.5
    maximised over logp <= x0 (itself bounded to the observed grid, so the knee can never
    be extrapolated outside the data). Returns None if fewer than 4 finite points are
    available, or the fit stays flat/degenerate (no clean rising sigmoid) — callers should
    fall back to another selection rule in that case.
    """
    logp = np.log10(grid)
    y = np.asarray(l2_errors, dtype=float)
    valid = ~np.isnan(y)
    if valid.sum() < 4:
        return None
    logp, y = logp[valid], y[valid]

    params, _ = _robust_sigmoid_fit(logp, y)
    if params is None:
        return None
    L, U, x0, k = params
    if k <= k_min_reliable or (U - L) < rise_min_frac * (y.max() - y.min() + 1e-9):
        return None

    logp_fine = np.linspace(logp.min(), x0, 400)
    curvature = np.abs(_sigmoid_d2(logp_fine, *params)) / (1 + _sigmoid_d1(logp_fine, *params) ** 2) ** 1.5
    return float(10 ** logp_fine[np.argmax(curvature)])


# Signature groups that SigProfilerAssignment always keeps together: if any member
# of a group is protected/selected, the rest of the group is protected/selected too.
_CONNECTED_SIG_GROUPS = [
    ["SBS2", "SBS13"],
    ["SBS7a", "SBS7b", "SBS7c", "SBS7d"],
    ["SBS10a", "SBS10b"],
    ["SBS17a", "SBS17b"],
]


def _expand_connected_sigs(names: set, all_names: list) -> set:
    """Add the rest of any connected-signature group `names` touches."""
    all_names = set(all_names)
    names = set(names) & all_names
    for group in _CONNECTED_SIG_GROUPS:
        group_present = set(group) & all_names
        if names & group_present:
            names |= group_present
    return names


def _fit_score(x: np.ndarray, W: np.ndarray, metric: str) -> tuple:
    """NNLS-fit (score, h) for W against x; higher score is always better."""
    h, _ = sp.optimize.nnls(W, x)
    pred = W @ h
    if metric == "cosine":
        return _cosine_similarity(x, pred), h
    if metric == "l2":
        # negate: lower relative L2 error = better = higher score
        return -np.linalg.norm(x - pred) / (np.linalg.norm(x) + 1e-30), h
    raise ValueError(f"Unknown metric {metric!r}; expected 'cosine' or 'l2'")


def _unknown_local(active_indices, dummy_indices):
    """Positions within active_indices that are unknown signatures."""
    return [pos for pos, sig in enumerate(active_indices) if sig in dummy_indices]


def _penalised_objective(
    weights: np.ndarray,
    W_curr: np.ndarray,
    x: np.ndarray,
    unknown_local_positions: list,
    unknown_penalty: float,
) -> float:
    """SSE plus L1 penalty on the unknown-signature weights."""
    predicted = W_curr @ np.maximum(weights, 0)
    sse = 0.5 * np.sum((predicted - x) ** 2)
    l1 = (unknown_penalty * np.sum(np.abs(weights[unknown_local_positions]))
          if unknown_local_positions else 0.0)
    return sse + l1


def _penalised_jacobian(
    weights: np.ndarray,
    W_curr: np.ndarray,
    x: np.ndarray,
    unknown_local_positions: list,
    unknown_penalty: float,
) -> np.ndarray:
    """Gradient of _penalised_objective with respect to weights."""
    w = np.maximum(weights, 0)
    grad = W_curr.T @ (W_curr @ w - x)
    if unknown_local_positions:
        grad[unknown_local_positions] += unknown_penalty * np.sign(w[unknown_local_positions])
    return grad


def _backward_eliminate(
    x, W_curr, active_indices, thresh_backward, max_iter, background=frozenset(), metric="cosine"
) -> tuple:
    """Drop the least-damaging signature each round until every drop exceeds thresh_backward.
    Indices in `background` are never removed."""
    for _ in range(max_iter):
        if len(active_indices) <= 1:
            break

        score_current, _ = _fit_score(x, W_curr, metric)

        # try dropping each removable signature one at a time; keep whichever
        # drop damages the fit the least
        best_drop, best_keep = np.inf, None
        for position, signature in enumerate(active_indices):
            if signature in background:
                continue
            keep = [j for j in range(len(active_indices)) if j != position]
            drop = score_current - _fit_score(x, W_curr[:, keep], metric)[0]
            if drop < best_drop:
                best_drop, best_keep = drop, keep

        if best_keep is None or best_drop >= thresh_backward:
            break  # all remaining signatures are significant, or none removable

        active_indices = [active_indices[j] for j in best_keep]
        W_curr = W_curr[:, best_keep]

    return W_curr, active_indices


def _initial_prune(x, W_known, background, initial_remove_penalty, max_iter, metric) -> list:
    """SPA-style coarse pre-filter: fit the full catalog, then one lenient backward-elimination
    pass (looser than thresh_backward) to shrink the pool before the stricter add/remove loop."""
    if initial_remove_penalty <= 0:
        return list(range(W_known.shape[1]))
    _, active_indices = _backward_eliminate(
        x, W_known, list(range(W_known.shape[1])), initial_remove_penalty, max_iter,
        background=background, metric=metric,
    )
    return active_indices


def _add_remove_known_signatures(
    x, W_known, background, thresh_forward, thresh_backward, metric, max_iter,
    initial_active=None,
) -> tuple:
    """SPA-style greedy layered selection: each layer adds whichever candidate scores best
    after also being backward-eliminated; stops once a layer fails to improve by more than
    thresh_forward. initial_active, if given, warm-starts the active set (e.g. survivors
    of a coarse pre-filter) instead of starting from just `background`."""
    n_known = W_known.shape[1]
    active_indices = sorted(initial_active if initial_active is not None else background)
    current_score = (
        _fit_score(x, W_known[:, active_indices], metric)[0] if active_indices else -np.inf
    )

    for _ in range(max_iter):
        candidates = [i for i in range(n_known) if i not in active_indices]
        if not candidates:
            break

        # try adding each remaining candidate one at a time, backward-eliminate the
        # trial set, and keep whichever candidate scores best after that pruning
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


def _add_unknown_signatures(x, W_curr, W_extended, active_indices, channel_number):
    """Phase 2: append one-hot unknowns for the highest-residual channels."""
    h_known, _ = sp.optimize.nnls(W_curr, x)
    pos_residuals = np.maximum(x - W_curr @ h_known, 0)
    if pos_residuals.sum() == 0:
        return W_extended, active_indices, set(), W_curr

    top_channels = np.argsort(pos_residuals)[-channel_number:]
    n_added = len(top_channels)
    unknown_sigs = np.zeros((W_extended.shape[0], n_added))
    unknown_sigs[top_channels, np.arange(n_added)] = 1.0

    W_extended = np.column_stack([W_extended, unknown_sigs])
    new_unknown_cols = list(range(W_extended.shape[1] - n_added, W_extended.shape[1]))
    active_indices = active_indices + new_unknown_cols
    return W_extended, active_indices, set(new_unknown_cols), W_extended[:, active_indices]


def _combine_unknowns(
    x, W_curr, W_known, active_indices, unknown_local, unknown_penalty,
    background, thresh_backward, metric, max_iter,
):
    """Phase 3: merge multiple unknown signatures into one profile and re-prune."""
    n_known_sigs = W_known.shape[1]
    h_temp, _ = sp.optimize.nnls(W_curr, x)
    result_combine = minimize(
        fun=lambda w: _penalised_objective(w, W_curr, x, unknown_local, unknown_penalty),
        jac=lambda w: _penalised_jacobian(w, W_curr, x, unknown_local, unknown_penalty),
        x0=h_temp,
        bounds=[(0, None)] * len(active_indices),
        method="SLSQP",
        options={"maxiter": 10_000, "ftol": 1e-12, "disp": False},
    )
    unknown_weights = np.maximum(result_combine.x[unknown_local], 0)
    total_unknown = unknown_weights.sum()

    # Merge into a single weighted-average profile
    unknown_profiles = W_curr[:, unknown_local]
    combined_profile = (
        unknown_profiles @ (unknown_weights / total_unknown)
        if total_unknown > 0
        else unknown_profiles.mean(axis=1)
    )

    # Rebuild: full known set + combined unknown as the last column
    W_curr = np.column_stack([W_known, combined_profile.reshape(-1, 1)])
    active_indices = list(range(n_known_sigs)) + [n_known_sigs]
    dummy_indices = {n_known_sigs}

    if thresh_backward > 0 and len(active_indices) > 1:
        W_curr, active_indices = _backward_eliminate(
            x, W_curr, active_indices, thresh_backward, max_iter,
            background=background, metric=metric,
        )

    # If backward elimination removed all knowns, re-add the best-fit one so
    # Phase 4's max_unknown_fraction constraint has known sigs to work with
    if all(i in dummy_indices for i in active_indices):
        best_k, _ = max(
            (
                (k, _fit_score(x, np.column_stack([W_known[:, k], combined_profile]), metric)[0])
                for k in range(n_known_sigs)
            ),
            key=lambda pair: pair[1],
        )
        W_full = np.column_stack([W_known, combined_profile.reshape(-1, 1)])
        active_indices = [best_k, n_known_sigs]
        W_curr = W_full[:, active_indices]

    return W_curr, active_indices, dummy_indices


def _final_fit(
    x, W_curr, active_indices, unknown_local, unknown_penalty,
    max_unknown_fraction, n_known_active, sample_name,
):
    """Phase 4: multi-start SLSQP fit with an L1 penalty on unknown weights."""
    known_local = [i for i in range(len(active_indices)) if i not in set(unknown_local)]

    constraints = []
    if unknown_local and known_local and 0.0 < max_unknown_fraction < 1.0:
        f = max_unknown_fraction
        constraints = [{"type": "ineq", "fun": lambda w, kl=known_local, ul=unknown_local, f=f:
                         f * w[kl].sum() - (1.0 - f) * w[ul].sum()}]

    best_result, best_obj, best_start_idx = None, np.inf, -1
    for start_idx in range(4):
        h_init = create_diverse_initialization(W_curr, x, start_idx, n_known_active)
        res = minimize(
            fun=lambda w: _penalised_objective(w, W_curr, x, unknown_local, unknown_penalty),
            jac=lambda w: _penalised_jacobian(w, W_curr, x, unknown_local, unknown_penalty),
            x0=h_init,
            bounds=[(0, None)] * len(active_indices),
            constraints=constraints,
            method="SLSQP",
            options={"maxiter": 10_000, "ftol": 1e-12, "disp": False},
        )
        obj = _penalised_objective(res.x, W_curr, x, unknown_local, unknown_penalty)
        print(f"[phase4 init] sample={sample_name}\tstart_idx={start_idx}\tobj={obj:.6g}\tconverged={res.success}", flush=True)
        if obj < best_obj:
            best_obj, best_result, best_start_idx = obj, res, start_idx

    print(f"[phase4 init] sample={sample_name}\tbest_start_idx={best_start_idx}\tbest_obj={best_obj:.6g}", flush=True)
    return np.maximum(best_result.x, 0)


def _prefit_known_and_unknown(
    x, W_known, sig_names, background_signatures, thresh_forward, thresh_backward,
    initial_remove_penalty, metric, max_iter, allow_unknown, channel_number,
) -> tuple:
    """Phases 0-2: SPA-style coarse pre-filter, greedy known-signature add/remove, and
    one-hot unknown placement. None of these depend on unknown_penalty, so callers that
    sweep a penalty grid (find_optimal_penalty) can run this once per sample and reuse
    the result for every grid point instead of repeating the expensive greedy
    add/remove search (the dominant cost) from scratch at each penalty value.
    """
    n_known_sigs = W_known.shape[1]
    W_extended = W_known
    dummy_indices = set()  # global column indices of unknown sigs

    protected_names = _expand_connected_sigs(set(background_signatures), sig_names)
    background = frozenset(sig_names.index(s) for s in protected_names)

    # Phase 0: SPA-style coarse pre-filter of the full catalog
    initial_active = _initial_prune(
        x, W_known, background, initial_remove_penalty, max_iter, metric
    )

    # Phase 1: greedy add/remove of known signatures, warm-started from the pre-filter
    W_curr, active_indices = _add_remove_known_signatures(
        x, W_known, background, thresh_forward, thresh_backward, metric, max_iter,
        initial_active=initial_active,
    )

    # Phase 2: one-hot unknowns for high-residual channels
    if allow_unknown:
        W_extended, active_indices, dummy_indices, W_curr = _add_unknown_signatures(
            x, W_curr, W_extended, active_indices, channel_number
        )

    return W_curr, active_indices, dummy_indices, W_extended, background, n_known_sigs


def enrich_signature_references(
    x: np.ndarray,
    W: "pd.DataFrame | np.ndarray",
    sample_name: str,
    outdir_prefix: str,
    channel_number: int,
    thresh_backward: float = 0.01,
    thresh_forward: float = 0.01,
    initial_remove_penalty: float = 0.05,
    metric: str = "l2",
    background_signatures: tuple = ("SBS1", "SBS5"),
    max_iter: int = 1000,
    allow_unknown: bool = False,
    unknown_penalty: float = 2.0,
    max_unknown_fraction: float = 0.8,
    true_weights: "np.ndarray | None" = None,
    save_convergence: bool = True,
    _prefit: tuple = None,
) -> pd.DataFrame:
    """
    Fit mutational signatures with optional unknown-signature discovery.

    Phases: (0) coarse pre-filter of the full catalog via initial_remove_penalty;
    (1) SPA-style greedy add/remove of known signatures, protecting background_signatures;
    (2) add one-hot unknowns for high-residual channels; (3) combine unknowns into one
    profile and re-prune; (4) final SLSQP fit with an L1 penalty on unknown weights.
    `metric` ("cosine" or "l2") drives phases 0, 1 & 3.

    _prefit, if given, is the tuple returned by _prefit_known_and_unknown() -- phases
    0-2 rerun with identical results whenever only unknown_penalty changes, so a caller
    sweeping a penalty grid (find_optimal_penalty) can compute it once and pass it here
    to skip straight to phases 3-4.

    Returns (weights DataFrame, predicted x) for sample_name.
    """
    if isinstance(W, pd.DataFrame):
        sig_names = W.columns.tolist()
        W_known = W.values.copy()
    else:
        sig_names = [f"Signature_{i}" for i in range(W.shape[1])]
        W_known = W.copy()

    if _prefit is not None:
        W_curr, active_indices, dummy_indices, W_extended, background, n_known_sigs = _prefit
    else:
        W_curr, active_indices, dummy_indices, W_extended, background, n_known_sigs = _prefit_known_and_unknown(
            x, W_known, sig_names, background_signatures, thresh_forward, thresh_backward,
            initial_remove_penalty, metric, max_iter, allow_unknown, channel_number,
        )

    # Phase 3: combine multiple unknowns into one profile and re-prune
    unknown_local = _unknown_local(active_indices, dummy_indices)
    if allow_unknown and len(unknown_local) > 1:
        W_curr, active_indices, dummy_indices = _combine_unknowns(
            x, W_curr, W_known, active_indices, unknown_local, unknown_penalty,
            background, thresh_backward, metric, max_iter,
        )
        unknown_local = _unknown_local(active_indices, dummy_indices)

    # Phase 4: final penalised SLSQP fit
    known_local = [i for i in range(len(active_indices)) if i not in set(unknown_local)]
    weights_final = _final_fit(
        x, W_curr, active_indices, unknown_local, unknown_penalty,
        max_unknown_fraction, len(known_local), sample_name,
    )

    n_unknown_total = W_extended.shape[1] - n_known_sigs
    all_names = sig_names + [f"Unknown_{i}" for i in range(n_unknown_total)]
    h_out = np.zeros(len(all_names))
    h_out[active_indices] = weights_final

    h_df = pd.DataFrame(h_out, index=all_names, columns=[sample_name])
    return h_df, W_curr @ weights_final


def find_optimal_penalty(
    x: np.ndarray,
    W: "pd.DataFrame | np.ndarray",
    sample_name: str,
    outdir_prefix: str,
    channel_number: int = 96,
    thresh_backward: float = 0.01,
    thresh_forward: float = 0.01,
    metric: str = "l2",
    background_signatures: tuple = ("SBS1", "SBS5"),
    max_iter: int = 1000,
    min_penalty: float = 0.1,
    max_penalty: float = 10,
    true_weights: "np.ndarray | None" = None,
    n_grid: int = 20,
    acceptance_frac: float = 0.005,
    max_unknown_fraction: float = 0.8,
    true_weights_df: "pd.DataFrame | None" = None,
    selection: str = "knee",
) -> pd.DataFrame:
    """
    Find the L1 penalty where unknowns genuinely improve the fit over knowns alone.

    selection="knee" (default): picks the *low* knee of the l2_error-vs-log10(penalty)
    curve — the max-curvature point on the rising part of the curve, before its inflection
    — via a robustly-fit sigmoid (see _knee_penalty). Falls back to the threshold rule
    below if the curve can't be fit reliably (e.g. too few finite points, or no clean
    rising sigmoid).

    selection="threshold": selects the highest λ where improvement(λ) = l2_error_known -
    l2_error(λ) exceeds acceptance_frac × l2_error_known. If no λ clears this threshold,
    returns a known-only fit.

    true_weights_df, if given, is the ground-truth per-signature weight column for this sample
    (index = signature name, including any out-of-reference signatures) — used ONLY to log a
    diagnostic TAE alongside the L2 error at each λ. It never influences the fit itself or the
    selection of optimal_penalty.
    """
    if selection not in ("knee", "threshold"):
        raise ValueError(f"Unknown selection {selection!r}; expected 'knee' or 'threshold'")
    W_arr = W.values if isinstance(W, pd.DataFrame) else W
    norm_x = float(np.linalg.norm(x))
    known_sig_names = W.columns.tolist() if isinstance(W, pd.DataFrame) else [f"Signature_{i}" for i in range(W_arr.shape[1])]
    true_counts_df = true_weights_df * np.sum(x) if true_weights_df is not None else None  # proportions -> counts
    fit_kwargs = {
        "x": x, "W": W, "sample_name": sample_name, "outdir_prefix": outdir_prefix,
        "channel_number": channel_number, "thresh_backward": thresh_backward,
        "thresh_forward": thresh_forward, "metric": metric,
        "background_signatures": background_signatures, "max_iter": max_iter,
        "allow_unknown": True, "max_unknown_fraction": max_unknown_fraction,
        "true_weights": true_weights,
    }

    h_known, _ = sp.optimize.nnls(W_arr, x)
    l2_error_known = float(np.linalg.norm(x - W_arr @ h_known)) / (norm_x + 1e-30)
    improvement_threshold = acceptance_frac * l2_error_known
    TAE_known = recalculate_TAE(pd.DataFrame(h_known, index=known_sig_names, columns=[sample_name]), true_counts_df, known_sig_names)
    print(f"[penalty search] sample={sample_name}\tl2_error_known={l2_error_known:.4f}\tTAE_known={TAE_known:.4f}\tthreshold={improvement_threshold:.4f}", flush=True)
    print(f"[penalty search] sample={sample_name}\tpenalty\tl2_error\timprovement\tTAE\tunknown_frac", flush=True)

    # Phases 0-2 (greedy known-signature selection + one-hot unknown placement) don't
    # depend on unknown_penalty, so compute them once here instead of once per grid point.
    prefit = _prefit_known_and_unknown(
        x, W_arr, known_sig_names, background_signatures, thresh_forward, thresh_backward,
        0.05, metric, max_iter, True, channel_number,
    )

    grid = np.logspace(np.log10(min_penalty), np.log10(max_penalty), n_grid)
    l2_errors = []
    for lam in grid:
        try:
            h_lam, x_pred = enrich_signature_references(**fit_kwargs, unknown_penalty=lam, save_convergence=False, _prefit=prefit)
            l2_error = float(np.linalg.norm(x - x_pred)) / (norm_x + 1e-30)
            tae = recalculate_TAE(h_lam, true_counts_df, known_sig_names)
            unknown_mask = h_lam.index.str.startswith("Unknown_")
            unknown_frac = float(h_lam.loc[unknown_mask, sample_name].sum() / (np.sum(x) + 1e-30))
        except Exception:
            l2_error, tae, unknown_frac = np.nan, float("nan"), np.nan
        print(f"[penalty search] sample={sample_name}\t{lam}\t{l2_error}\t{l2_error_known - l2_error}\t{tae}\t{unknown_frac}", flush=True)
        l2_errors.append(l2_error)

    improvements = l2_error_known - np.array(l2_errors)
    useful = improvements >= improvement_threshold

    optimal_penalty = None
    if selection == "knee":
        optimal_penalty = _knee_penalty(grid, l2_errors)
        if optimal_penalty is not None:
            print(f"  [penalty search] knee λ={optimal_penalty:.4g}  (curvature-based)")
        else:
            print("  [penalty search] knee fit unreliable → falling back to threshold rule")

    if optimal_penalty is None:
        if not useful.any():
            print(f"  [penalty search] unknowns do not improve fit (max={np.nanmax(improvements):.4f} < {improvement_threshold:.4f}) → known-only")
            return pd.DataFrame(h_known, index=known_sig_names, columns=[sample_name])
        optimal_penalty = float(grid[useful][-1])
        print(f"  [penalty search] optimal λ={optimal_penalty:.4g}  (threshold-based)  improvement={improvements[useful][-1]:.4f}  threshold={improvement_threshold:.4f}")

    h, _ = enrich_signature_references(**fit_kwargs, unknown_penalty=optimal_penalty, save_convergence=True, _prefit=prefit)
    return h


def main():
    parser = argparse.ArgumentParser(description='Extended signature fitting method')
    parser.add_argument('-i', '--input', required=True, help='Input mutation matrix')
    parser.add_argument('-w', '--weights', required=True, help='COSMIC signature matrix')
    parser.add_argument('-m', '--matrix-order', required=True, help='Matrix order file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('--unknown-penalty', type=str, default='0.001',
                        help='L1 penalty for unknown signatures (single value or comma-separated list). '
                             'Ignored when --auto-penalty is set.')
    parser.add_argument('--thresh-backward', type=float, default=0.001,
                        help='Cosine similarity threshold for backward elimination.')
    parser.add_argument('--max-iter', type=int, default=1000,
                        help='Maximum iterations.')
    parser.add_argument('--auto-penalty', action='store_true',
                        help='Automatically find the optimal L1 penalty per sample via binary '
                             'search on the reconstruction residual. Overrides --unknown-penalty.')
    parser.add_argument('--min-penalty', type=float, default=0.001,
                        help='Lower bound of the penalty search range (used with --auto-penalty).')
    parser.add_argument('--max-penalty', type=float, default=10000.0,
                        help='Upper bound of the penalty search range (used with --auto-penalty).')
    parser.add_argument('--true-weights', type=str, default=None,
                        help='Optional ground-truth weights file (signatures x samples, tab-separated, '
                             'index_col=0) used ONLY to log a diagnostic TAE per penalty during '
                             '--auto-penalty. Never influences the fit itself. For simulated-data testing only.')
    args = parser.parse_args()

    X = pd.read_csv(args.input, index_col=0, sep='\t')
    W = pd.read_csv(args.weights, sep='\t', index_col=0)
    MP_index = pd.read_csv(args.matrix_order, header=None).squeeze()

    W = W.reindex(MP_index)
    X = X.reindex(MP_index)

    true_weights_all = None
    if args.true_weights:
        true_weights_all = pd.read_csv(args.true_weights, sep='\t', index_col=0)

    outdir = os.path.dirname(args.output) or '.'
    os.makedirs(outdir, exist_ok=True)

    if args.auto_penalty:
        # ── Automatic per-sample penalty search ───────────────────────────────
        H_results = []
        for sample_name in X.columns:
            print(f"[auto-penalty] sample: {sample_name}")
            true_weights_df = (
                true_weights_all[[sample_name]] if true_weights_all is not None else None
            )
            h = find_optimal_penalty(
                x=X[sample_name].values,
                W=W,
                sample_name=sample_name,
                outdir_prefix=outdir,
                channel_number=W.shape[0],
                thresh_backward=args.thresh_backward,
                max_iter=args.max_iter,
                min_penalty=args.min_penalty,
                max_penalty=args.max_penalty,
                true_weights=None,
                true_weights_df=true_weights_df,
            )
            H_results.append(h)

        H = pd.concat(H_results, axis=1)
        H.to_csv(args.output)

    else:
        # ── Fixed-penalty sweep (original behaviour) ───────────────────────────
        penalty_str = args.unknown_penalty.strip('[]').replace(' ', '')
        penalties = [float(p) for p in penalty_str.split(',')]

        for penalty in penalties:
            H_results = []
            for sample_name in X.columns:
                h, _ = enrich_signature_references(
                    x=X[sample_name].values,
                    W=W,
                    sample_name=sample_name,
                    outdir_prefix=outdir,
                    channel_number=W.shape[0],
                    thresh_backward=args.thresh_backward,
                    max_iter=args.max_iter,
                    allow_unknown=True,
                    unknown_penalty=penalty,
                    true_weights=None,
                )
                H_results.append(h)

            H = pd.concat(H_results, axis=1)
            base_name = args.output.replace('_contribution.dat', '')
            output_file = f"{base_name}_contribution_penalty_{penalty}.dat"
            H.to_csv(output_file)


if __name__ == '__main__':
    main()
