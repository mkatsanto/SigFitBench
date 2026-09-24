#!/usr/bin/env python3
"""
Extended signature fitting method (DM variant): a one-hot dummy signature is added
for every mutation-type channel, then an L2-penalized SLSQP fit is used together with
SPA-style epsilon-based backward elimination (dummy signatures penalized more heavily
than real ones) to prune down to a parsimonious signature set, with SPA's connected
signature groups restored afterwards.

Usage:
    python extended_dm.py -i <input_data> -w <cosmic_signatures> -m <matrix_order> -o <output_file>
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize


# L1-style penalty on dummy-signature weights (fit against mutation-fraction scale,
# not raw counts). Kept fixed rather than exposed as a CLI/config parameter.
PENALTY_MAG = 0.1

# custom arguments for minimize()
MINIMIZE_OPTIONS = {"ftol": 1e-14, "maxiter": 10000}

# Signature groups that SigProfilerAssignment always keeps together. If only some
# members of a group survive pruning, the rest are added back.
CONNECTED_SIG_GROUPS = [
    ["SBS2", "SBS13"],
    ["SBS7a", "SBS7b", "SBS7c", "SBS7d"],
    ["SBS10a", "SBS10b"],
    ["SBS17a", "SBS17b"],
]


def _add_which_connected(active_sigs):
    """Names from CONNECTED_SIG_GROUPS whose group-mates are active but they are not."""
    add_which = set()
    for sig in active_sigs:
        for group in CONNECTED_SIG_GROUPS:
            if sig in group:
                add_which |= {s for s in group if s not in active_sigs}
    return add_which


def _l2_penalized(x, ref, num_true_sigs, catalog):
    """L2 reconstruction error plus a penalty on dummy-signature weights."""
    return (
        np.sqrt(np.power(np.matmul(ref, x) - catalog, 2).sum()) + PENALTY_MAG * x[num_true_sigs:].sum()
    )


def _jacobian_l2_penalized(x, ref, num_true_sigs, catalog):
    """Gradient of _l2_penalized with respect to x."""
    resid = np.matmul(ref, x) - catalog
    l2_norm = np.sqrt(np.power(resid, 2).sum())
    ref_values = ref.values if hasattr(ref, "values") else ref
    grad = np.matmul(ref_values.T, resid) / l2_norm
    grad = grad + PENALTY_MAG * np.concatenate(
        [np.zeros(num_true_sigs), np.ones(x.size - num_true_sigs)]
    )
    return grad


def _epsilon_spa(x, ref, catalog):
    """Relative L2 reconstruction error (SPA's epsilon)."""
    return np.sqrt(np.power(np.matmul(ref, x) - catalog, 2).sum()) / np.sqrt(np.power(catalog, 2).sum())


def _fit(ref, catalog):
    """SLSQP fit of _l2_penalized from a random initialisation on the simplex."""
    num_true_sigs = ref.columns.str.startswith("SBS").sum()
    bounds = [(0, None) for _ in range(ref.shape[1])]
    w_init = np.random.rand(ref.shape[1])
    w_init /= w_init.sum()
    res = minimize(
        _l2_penalized, w_init, args=(ref, num_true_sigs, catalog),
        jac=_jacobian_l2_penalized, bounds=bounds, method="SLSQP", options=MINIMIZE_OPTIONS
    )
    if res.success == False: print('NNLS failed to converge')
    return res


def _iteratively_remove_signatures(ref_reduced, eps_now, catalog, thresh_backward, dummy_sig_multiplier):
    """Drop the least-damaging signature each round (dummy sigs penalized
    dummy_sig_multiplier x more heavily) until every drop would increase epsilon by
    more than thresh_backward."""
    while True:
        smallest_increase, new_epsilon, chosen_sig = np.inf, None, None
        for sig in ref_reduced.columns:
            ref_test = ref_reduced.drop(columns=[sig])
            res = _fit(ref_test, catalog)
            eps_candidate = _epsilon_spa(res.x, ref_test, catalog)
            score_worsening = eps_candidate - eps_now
            if not sig.startswith("SBS"):
                score_worsening *= dummy_sig_multiplier
            if score_worsening < smallest_increase:
                smallest_increase, new_epsilon, chosen_sig = score_worsening, eps_candidate, sig

        if smallest_increase < thresh_backward:
            eps_now = new_epsilon
            ref_reduced = ref_reduced.drop(columns=[chosen_sig])
        else:
            break
    return ref_reduced, eps_now


def fit_sample(
    x, W, sig_names, sample_name,
    init_sig_rem_thr=0.01, thresh_backward=0.01, use_connected_sigs=True,
):
    """
    Fit one sample: extend the reference with a one-hot dummy signature per channel,
    run an L2-penalized SLSQP fit, coarsely drop signatures whose removal barely
    changes epsilon, iteratively backward-eliminate the rest (dummy signatures
    penalized dummy_sig_multiplier x more heavily than real ones, where the
    multiplier is set from the sample's total mutation count), then restore any
    connected-signature group (as defined by SigProfilerAssignment) left incomplete
    by the pruning.

    Returns a pandas Series (mutation-count scale) indexed by sig_names plus one
    dummy column per channel ("E0".."E{n_channels-1}").
    """
    ref_sigs = pd.DataFrame(np.asarray(W), columns=sig_names)
    n_channels = ref_sigs.shape[0]
    num_muts = float(np.sum(x))
    catalog = np.asarray(x, dtype=np.float64) / num_muts
    dummy_sig_multiplier = 2 if num_muts < 20000 else 4
    all_names = sig_names + [f"E{n}" for n in range(n_channels)]

    # construct the extended reference catalog with all reference signatures and dummies
    ref_extended = ref_sigs.copy()
    for n in range(n_channels):
        col = f"E{n}"
        ref_extended[col] = 0.0
        ref_extended.iloc[n, ref_extended.columns.get_loc(col)] = 1.0

    # initial fast removal of signature without optimization
    res = _fit(ref_extended, catalog)
    eps_a = _epsilon_spa(res.x, ref_extended, catalog)
    orig_sol = res.x.copy()
    drop_init = []
    for i, sig in enumerate(ref_extended.columns):
        test_sol = orig_sol.copy()
        test_sol[i] = 0.0
        eps_test = _epsilon_spa(test_sol, ref_extended, catalog)
        increase = eps_test - eps_a
        if not sig.startswith("SBS"):
            increase *= dummy_sig_multiplier
        if increase < init_sig_rem_thr:
            drop_init.append(sig)

    # iteratively remove signatures one by one
    sigs_keep = [s for s in ref_extended.columns if s not in drop_init]
    ref_reduced = ref_extended[sigs_keep]
    res = _fit(ref_reduced, catalog)
    eps_now = _epsilon_spa(res.x, ref_reduced, catalog)
    ref_reduced, eps_now = _iteratively_remove_signatures(
        ref_reduced, eps_now, catalog, thresh_backward, dummy_sig_multiplier
    )

    # add connected signatures
    if use_connected_sigs:
        add_which = _add_which_connected(set(ref_reduced.columns))
        for sig in add_which:
            ref_reduced[sig] = ref_sigs[sig]
        # order the columns as in the reference signatures (dummies last)
        col_order = []
        for sig in ref_sigs.columns:
            if sig in ref_reduced.columns: col_order.append(sig)
        for sig in ref_reduced.columns:
            if sig.startswith('E') and sig not in ref_sigs.columns:
                col_order.append(sig)
        ref_reduced = ref_reduced[col_order]

    # final optimization with the selected signatures used as a reference
    res = _fit(ref_reduced, catalog)


    # add normalization by hand and construct the resulting object
    weights = res.x / res.x.sum() if res.x.sum() > 0 else res.x
    h_out = pd.Series(0.0, index=all_names, name=sample_name)
    h_out[ref_reduced.columns] = weights * num_muts
    return h_out


def main():
    parser = argparse.ArgumentParser(description="Extended signature fitting method (DM variant)")
    parser.add_argument('-i', '--input', required=True, help='Input mutation matrix')
    parser.add_argument('-w', '--weights', required=True, help='COSMIC signature matrix')
    parser.add_argument('-m', '--matrix-order', required=True, help='Matrix order file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('--init-sig-rem-thr', type=float, default=0.01,
                        help='Epsilon-increase threshold for the coarse initial signature drop')
    parser.add_argument('--thresh-backward', type=float, default=0.01,
                        help='Epsilon-increase threshold for iterative backward elimination')
    parser.add_argument('--seed', type=int, default=0, help='Random seed for reproducibility')
    args = parser.parse_args()

    np.random.seed(args.seed)

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
        h = fit_sample(
            x=X[sample_name].values,
            W=W.values,
            sig_names=sig_names,
            sample_name=sample_name,
            init_sig_rem_thr=args.init_sig_rem_thr,
            thresh_backward=args.thresh_backward,
        )
        H_results.append(h)

    H = pd.concat(H_results, axis=1)
    H.to_csv(args.output)
    print(f"Results saved to {args.output}")


if __name__ == '__main__':
    main()
