#!/usr/bin/env python3
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

H, model = musical.refit.refit(X, W, method='likelihood_bidirectional', thresh=0.001)
H.to_csv(snakemake.output.result)
