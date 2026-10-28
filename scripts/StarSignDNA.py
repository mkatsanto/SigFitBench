#!/usr/bin/env python3
"""
StarSignDNA signature refitting via its `starsigndna refit` command.

StarSignDNA fits exposures with a Poisson likelihood and an L1 (sparsity) penalty, whose
strength depends on the data type (genome or exome). Each sample is fitted on its own and its
exposures are the median over StarSignDNA's bootstrap runs.
Output: mutations per signature (signatures x samples).
"""
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

import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.stdout = open(snakemake.log.stdout, "w")  # noqa: F821
sys.stderr = open(snakemake.log.stderr, "w")  # noqa: F821

N_ITERATIONS = snakemake.params.n_iterations  # noqa: F821
DATA_TYPE = snakemake.params.data_type  # noqa: F821

counts = pd.read_csv(snakemake.input.data, sep="\t", index_col=0)  # noqa: F821
signatures = pd.read_csv(snakemake.input.signature_db, sep="\t", index_col=0)  # noqa: F821

# Fit one sample at a time. StarSignDNA input: samples x contexts and signatures x contexts,
# tab-separated
work = Path("starsigndna_work")
work.mkdir(exist_ok=True)
signature_file = work / "signatures.txt"
signatures.loc[counts.index].T.to_csv(signature_file, sep="\t", index_label=False)
out_dir = work / "output"
out_dir.mkdir(exist_ok=True)  # StarSignDNA only creates it for VCF input

result = pd.DataFrame(0.0, index=signatures.columns, columns=counts.columns)
for i, s in enumerate(counts.columns):
    matrix_file = work / f"catalogue_{i}.txt"
    counts[[s]].T.to_csv(matrix_file, sep="\t", index_label=False)
    cmd = ["starsigndna", "refit", str(matrix_file), str(signature_file),
           "--output-folder", str(out_dir), "--n-iterations", str(N_ITERATIONS),
           "--data-type", DATA_TYPE]
    print("running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=sys.stdout, stderr=sys.stderr)

    # StarSignDNA names its output files <matrix>_<signatures>_<iterations>; for one sample
    # it reports the median exposures over its bootstrap runs, as shares of the sample
    run_name = f"{matrix_file.stem}_{signature_file.stem}_{N_ITERATIONS}"
    median = pd.read_csv(out_dir / f"StarSign_exposure_median_{run_name}.txt", sep="\t", index_col=0)
    shares = median.iloc[:, 0] if median.shape[1] == 1 else median.iloc[0]
    result.loc[shares.index, s] = shares * counts[s].sum()

result.round(4).to_csv(snakemake.output.result)  # noqa: F821
