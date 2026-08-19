#!/usr/bin/env python3
import sys

sys.stdout = open(snakemake.log.stdout, "w")  # noqa: F821
sys.stderr = open(snakemake.log.stderr, "w")  # noqa: F821

from SigProfilerAssignment import Analyzer as Analyze
import pandas as pd

Analyze.cosmic_fit(
    samples=snakemake.input.data,
    input_type='matrix',
    output=".",
    signature_database=snakemake.input.signature_db,
    make_plots=False,
    sample_reconstruction_plots=False,
    initial_remove_penalty=snakemake.params.initial_remove_penalty
)

df = pd.read_csv(
    "Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep='\t', index_col='Samples'
).T
df.to_csv(snakemake.output.result)
