# Results and metrics

## Where results go

```
{output_dir}/
├── report_synth_realistic.html            ← open this (one report per data type)
├── summary/                               one table each, all runs together
│   ├── evaluation_metrics.tsv             one row per tool, tissue, mutation count, repeat
│   ├── per_sample_metrics.tsv             one row per tool and sample
│   └── tool_benchmarks.tsv                run time (seconds) and peak memory (max_rss_mb)
├── runs/{data_type}/{tissue}/{N}_mutations/
│   ├── samples/counts_rep{r}.dat          mutation counts the tools fit
│   ├── samples/known_answer_rep{r}.dat    true signature mix of each sample
│   ├── fits/{tool}_rep{r}.dat             what each tool found: mutations per signature and sample
│   ├── scores/{tool}_rep{r}.tsv           scores of each sample
│   └── benchmarks/{tool}_rep{r}.tsv       run time and memory
└── internal/                              files the workflow needs, rarely opened
```

Everything about one run is in one folder: to see why a tool scored badly for Head-SCC at
1,000 mutations, compare `samples/known_answer_rep0.dat` with `fits/<tool>_rep0.dat` in
`runs/synth_realistic/Head-SCC/1000_mutations/`. `{r}` is the repeat (`num_realizations`),
starting at 0. For `synth_custom` data, the dataset name takes the place of the tissue.
For `real` data, the folder is `runs/real/{tissue}/all_mutations/`, since real samples
keep their own mutation counts, and the repeat is always 0.

Logs follow the same folders under `{log_dir}`, e.g.
`{log_dir}/runs/synth_realistic/Head-SCC/1000_mutations/SPA_rep0.stderr.log`.

If a tool fails, its `fits/{tool}_rep{r}.dat` holds one line starting with
`# TOOL FAILED:` and the reason. The rest of the run continues. The full error is in the
tool's `.stderr.log` in the same run folder under `{log_dir}`.

## Scores against the known answer

For `synth_realistic` and `synth_custom` data.

| Score | Range | Better | Meaning |
|---|---|---|---|
| **TAE** | 0–1 | lower | Total absolute error: half the summed difference between found and true weights. 0.1 means 10% of mutations went to the wrong signatures. |
| **Precision** | 0–1 | higher | Share of reported signatures that are really there. Low = false signatures. |
| **Recall** | 0–1 | higher | Share of real signatures that were found. Low = missed signatures. |
| **F1** | 0–1 | higher | Balance of precision and recall. |
| **MCC** | −1–1 | higher | Matthews correlation: agreement on which signatures are active, fair even when most signatures are absent. |

A signature counts as found when a tool gives it at least `evaluation.min_contribution`
mutations (default 10).

`summary/evaluation_metrics.tsv` gives each score as a mean over samples (`P_mean`, `R_mean`,
`F1_mean`, `MCC_mean`, ...) with its spread (`_std`, `_sem`). It also has:

| Column | Meaning |
|---|---|
| `wtot_FP` | Share of mutations given to signatures that are not there |
| `wtot_FN` | Share of mutations of real signatures that were missed |
| `num_FP_sigs`, `num_FN_sigs` | Average number of false and missed signatures per sample |

## The report

One report per data type: `report_synth_realistic.html` for `synth_realistic`,
`report_synth_custom.html` for `synth_custom`. Open it in a browser.

One tab per tissue, plus an "Aggregate Analysis" tab for all tissues together. Each tab
has:

- **Tool comparison:** scores against mutations per sample, and precision against recall
- **Performance distribution per tool:** spread of scores over samples
- **Sample-by-sample comparison:** samples that are hard for every tool (with more than 10 samples per run)
- **Signature complexity vs. performance:** only for runs with unknown signatures
- **Signature activity per sample:** known vs. found signature weights, one sample at a time

Methods, settings and score definitions are in the appendix at the end of the report.

## Comparing tools without a known answer

For `real` data, the report is `report_real.html`.
Tools are compared with each other, per tissue:

| Measure | Meaning |
|---|---|
| Cosine similarity | How similar two tools' signature weights are (1 = same) |
| Pearson correlation | Same idea, as a correlation |
| Jaccard (active signatures) | Overlap of the signatures two tools report as active (1 = same set) |
| TAE (tool vs. tool) | Share of mutations the two tools assign differently |
| Unassigned mutations | Share of each sample's mutations a tool leaves out |

The same data is in `summary/`, and for each tissue in its run folder (`runs/real/{tissue}/all_mutations/agreement/`):

| File | Contents |
|---|---|
| `tool_agreement.tsv` | The measures above for each pair of tools |
| `signature_calls.tsv` | For each signature and tool: in how many samples it is active, and how many tools report it |
| `signature_weights.tsv` | Weights per signature, sample and tool |
| `sample_activity.tsv` | Active signatures per sample and tool |

A signature reported by most tools is a safer call than one only a single tool reports.
