# Results and metrics

## Where results go

```
{output_dir}/
├── {data_type}/{tissue}/standard/{num_muts}/rep_{n}/
│   └── {tool}_contribution.dat          # mutations per signature (rows) and sample (columns)
├── generated_data/{data_type}/
│   ├── counts_*.dat                     # simulated catalogs given to the tools
│   └── weights_*.dat                    # the known signature mix of each sample
├── evaluation/                          # scores per tool and run
├── benchmarks/                          # run time and memory per job
├── all_evaluation_metrics.tsv           # one row per tool, tissue, mutation count, realization
├── all_per_sample_metrics.tsv           # one row per tool and sample
├── all_tool_benchmarks.tsv              # run time (seconds) and peak memory (max_rss_mb)
└── sigfitbench_report_*.html            # reports
```

If a tool fails, its `{tool}_contribution.dat` holds one line starting with
`# TOOL FAILED:` and the reason. The rest of the run continues. The full error is in
`{log_dir}/{data_type}/{tissue}/standard/{num_muts}/rep_{n}/{tool}_stderr.log`, with the
tool name in lower case (e.g. `musical_stderr.log`).

## Scores against the known answer

For `realsim` and `custom` data.

| Score | Range | Better | Meaning |
|---|---|---|---|
| **TAE** | 0–1 | lower | Total absolute error: half the summed difference between found and true weights. 0.1 means 10% of mutations went to the wrong signatures. |
| **Precision** | 0–1 | higher | Share of reported signatures that are really there. Low = false signatures. |
| **Recall** | 0–1 | higher | Share of real signatures that were found. Low = missed signatures. |
| **F1** | 0–1 | higher | Balance of precision and recall. |
| **MCC** | −1–1 | higher | Matthews correlation: agreement on which signatures are active, fair even when most signatures are absent. |

A signature counts as found when a tool gives it at least `evaluation.min_contribution`
mutations (default 10).

`all_evaluation_metrics.tsv` gives each score as a mean over samples (`P_mean`, `R_mean`,
`F1_mean`, `MCC_mean`, ...) with its spread (`_std`, `_sem`). It also has:

| Column | Meaning |
|---|---|
| `wtot_FP` | Share of mutations given to signatures that are not there |
| `wtot_FN` | Share of mutations of real signatures that were missed |
| `num_FP_sigs`, `num_FN_sigs` | Average number of false and missed signatures per sample |

## The report

One report per data type: `sigfitbench_report_realistic_simulations.html` for
`realsim`, `sigfitbench_report_synthetic.html` for `custom`. Open it in a browser.

One tab per tissue, plus an "Aggregate Analysis" tab for all tissues together. Each tab
has:

- **Tool comparison:** scores against mutations per sample, and precision against recall
- **Performance distribution per tool:** spread of scores over samples
- **Sample-by-sample comparison:** samples that are hard for every tool (with more than 10 samples per run)
- **Signature complexity vs. performance:** only for runs with hidden signatures
- **Signature activity per sample:** known vs. found signature weights, one sample at a time

Methods, settings and score definitions are in the appendix at the end of the report.

## Comparing tools without a known answer

For `real` data, the report is `sigfitbench_report_real_samples.html`.
Tools are compared with each other, per tissue:

| Measure | Meaning |
|---|---|
| Cosine similarity | How similar two tools' signature weights are (1 = same) |
| Pearson correlation | Same idea, as a correlation |
| Jaccard (active signatures) | Overlap of the signatures two tools report as active (1 = same set) |
| TAE (tool vs. tool) | Share of mutations the two tools assign differently |
| Unassigned mutations | Share of each sample's mutations a tool leaves out |

The same data is in:

| File | Contents |
|---|---|
| `all_tool_agreement.tsv` | The measures above for each pair of tools |
| `all_signature_calls.tsv` | For each signature and tool: in how many samples it is active, and how many tools report it |
| `all_signature_weights.tsv` | Weights per signature, sample and tool |
| `all_sample_activity.tsv` | Active signatures per sample and tool |

A signature reported by most tools is a safer call than one only a single tool reports.
