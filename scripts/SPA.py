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

sys.stdout = open(snakemake.log.stdout, "w")  
sys.stderr = open(snakemake.log.stderr, "w") 

from SigProfilerAssignment import Analyzer as Analyze
import pandas as pd

# Fit one sample at a time: each sample gets its own input file and output folder
counts = pd.read_csv(snakemake.input.data, sep='\t', index_col=0)
fits = []
for s in counts.columns:
    os.makedirs(s, exist_ok=True)
    sample_file = os.path.join(s, "counts.txt")
    counts[[s]].to_csv(sample_file, sep='\t')
    Analyze.cosmic_fit(
        samples=sample_file,
        input_type='matrix',
        output=s,
        signature_database=snakemake.input.signature_db,
        make_plots=False,
        sample_reconstruction_plots=False,
        initial_remove_penalty=snakemake.params.initial_remove_penalty,
        cpu=snakemake.threads,
    )
    fits.append(pd.read_csv(
        os.path.join(s, "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"),
        sep='\t', index_col='Samples'
    ).T)

df = pd.concat(fits, axis=1)
df.to_csv(snakemake.output.result)
