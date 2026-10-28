#!/usr/bin/env python3
import os
import sys
import traceback


def record_failure(exc_type, exc, tb):
    """If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and
    finish, so the rest of the workflow and the HTML report still run. The error stays in
    the log."""
    traceback.print_exception(exc_type, exc, tb)
    reason = " ".join(str(exc).split()) or exc_type.__name__
    for path in snakemake.output:
        with open(path, "w") as f:
            f.write(f"# TOOL FAILED: {reason}\n")
    print(f"TOOL FAILED ({snakemake.rule}): {reason}", file=sys.__stderr__)
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


sys.excepthook = record_failure

import sys
import musical
import pandas as pd

sys.stdout = open(snakemake.log.stdout, "w")  # noqa: F821
sys.stderr = open(snakemake.log.stderr, "w")  # noqa: F821

X = pd.read_csv(snakemake.input.data, index_col=0, sep='\t')
W = pd.read_csv(snakemake.input.weights, sep='\t', index_col=0)
MP_index = pd.read_csv(snakemake.input.matrix_order, header=None).squeeze()

W = W.reindex(MP_index)
X = X.reindex(MP_index, axis=0)

# Fit one sample at a time
H = pd.concat(
    [musical.refit.refit(X[[s]], W, method='likelihood_bidirectional', thresh=0.001)[0]
     for s in X.columns],
    axis=1,
)
H.to_csv(snakemake.output.result)
