#!/usr/bin/env python3
"""
StarSignDNA signature refitting via its `starsigndna refit` command.

StarSignDNA fits exposures with a Poisson likelihood and an L1 (sparsity) penalty, whose
strength depends on the data type (genome or exome). With several samples it fits them together
and keeps only signatures whose average exposure across the samples passes its threshold.
Output: mutations per signature (signatures x samples).
"""
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

# StarSignDNA input: samples x contexts and signatures x contexts, tab-separated
work = Path("starsigndna_work")
work.mkdir(exist_ok=True)
matrix_file, signature_file = work / "catalogue.txt", work / "signatures.txt"
counts.T.to_csv(matrix_file, sep="\t", index_label=False)
signatures.loc[counts.index].T.to_csv(signature_file, sep="\t", index_label=False)

out_dir = work / "output"
out_dir.mkdir(exist_ok=True)  # StarSignDNA only creates it for VCF input
cmd =["starsigndna", "refit", str(matrix_file), str(signature_file),
       "--output-folder", str(out_dir), "--n-iterations", str(N_ITERATIONS),
       "--data-type", DATA_TYPE]
print("running:", " ".join(cmd), flush=True)
subprocess.run(cmd, check=True, stdout=sys.stdout, stderr=sys.stderr)

# StarSignDNA names its output files <matrix>_<signatures>_<iterations>
run_name = f"{matrix_file.stem}_{signature_file.stem}_{N_ITERATIONS}"
if counts.shape[1] > 1:
    proportions = pd.read_csv(out_dir / f"{run_name}_threshold.txt", sep="\t", index_col=0)
else:  # single sample: bootstrap median
    median = pd.read_csv(out_dir / f"StarSign_exposure_median_{run_name}.txt", sep="\t", index_col=0)
    proportions = median.T.set_axis(counts.columns, axis=0)

# Proportions -> mutation counts (samples without mutations stay 0)
result = pd.DataFrame(0.0, index=signatures.columns, columns=counts.columns)
fitted = proportions.mul(counts.sum(axis=0).reindex(proportions.index), axis=0).T
result.loc[fitted.index, fitted.columns] = fitted
result.round(4).to_csv(snakemake.output.result)  # noqa: F821
