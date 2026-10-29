# Guide

This guide walks through the questions SigFitBench can answer. Each one has a
ready-to-run config in
[`config/scenarios/`](https://github.com/mkatsanto/SigFitBench/tree/main/config/scenarios).
Run any of them the same way:

```bash
conda activate sigfitbench
snakemake --profile profiles/local-conda \
  --snakefile workflow/Snakefile --configfile config/scenarios/<scenario>.yaml
```

Add `-n` first to see what would run without running it.

To make your own, copy the closest scenario, give it a new `output_dir` and `log_dir`,
and change the settings described below. Every scenario has the same
"Shared settings" block at the end; you rarely need to change it.

| Scenario | Known answer? | Time on a laptop |
|---|---|---|
| [Quick test](#quick-test) | Yes | minutes |
| [Unknown signatures](#unknown-signatures) | Yes | about an hour |
| [Mutation burden](#mutation-burden) | Yes | about an hour |
| [Custom mixtures](#custom-mixtures) | Yes | minutes |
| [Own catalog](#own-catalog) | No | minutes |
| [Real PCAWG](#real-pcawg) | No | minutes |
| [Full benchmark](#full-benchmark) | Yes | days: use a cluster |

Times are after the tool environments exist. The first run of any scenario also builds
them (30–60 minutes).

---

## Quick test

**Config:** [`quick_test.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/quick_test.yaml)

**Question:** does the installation work?

Three fast tools (`MutationalPatterns`, `deconstructSigs`, `SPA`), one tissue
(`Head-SCC`), three samples, 1,000 and 10,000 mutations. No unknown signatures.

**Check:** `results_quick_test/sigfitbench_report_realistic_simulations.html` opens and
shows three tools. If a tool is missing, see its log in `logs_quick_test/`.

**Next:** add the tools you want to compare to `tools_to_run` and run it again. Only
the new tools run; finished results are reused.

---

## Unknown signatures

**Config:** [`unknown_signatures.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/unknown_signatures.yaml)

**Question:** what happens when a tumour carries a signature the tools don't know about?

COSMIC adds signatures with nearly every release, and a tumour can carry a process that
is not in any catalogue. A tool given only known signatures has to put those mutations
somewhere. Most give them to whichever known signatures look closest, which adds false
signatures and distorts the weights of the real ones.

This scenario simulates that. The signatures in `out_of_reference_sigs` are removed from
the reference the tools fit against. Each sample then gets one of them, picked at random,
carrying 20% of its mutations.

**Settings that matter:**

```yaml
out_of_reference_sigs: ["SBS10c", "SBS10d", "SBS86", ...]   # signatures the tools never see
out_of_reference_weights: [0.2]   # one hidden signature per sample, 20% of its mutations
```

Each value in `out_of_reference_weights` adds **one** hidden signature to every sample:

| `out_of_reference_weights` | Hidden signatures per sample | Share of mutations hidden |
|---|---|---|
| `[0]` | none | 0% (baseline) |
| `[0.2]` | 1 | 20% |
| `[0.1, 0.1]` | 2 | 20% |
| `[0.1, 0.2]` | 2 | 30% |

**Comparing levels:** this is not a list of levels to compare. To compare 0%, 10% and
20%, make one copy of the config per level, each with its own `output_dir`, and run them
one after the other.

**Read the results:** compare precision with the 0% baseline. A tool whose precision
drops a lot is filling the gap with signatures that are not there. In the report, the
"Signature complexity vs. performance" section appears only for runs with hidden
signatures.

---

## Mutation burden

**Config:** [`mutation_burden.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/mutation_burden.yaml)

**Question:** how does accuracy change from samples with few mutations (exomes, panels)
to samples with many (whole genomes), and how much of the error is chance?

One tissue, 100 to 50,000 mutations per sample, 20 signature mixes, each drawn three
times.

**Settings that matter:**

```yaml
num_muts: [100, 500, 2000, 10000, 50000]   # one run per value
num_samples: 20          # 20 different signature mixes
num_realizations: 3      # each mix drawn 3 times with fresh random counts
```

`num_samples` and `num_realizations` test different things:

- More **samples** means more *different* tumours: a wider range of signature mixes.
- More **realizations** means the *same* tumours drawn again. The known answer stays the
  same; only chance changes. If a tool's score varies a lot between realizations, its
  errors at that burden are mostly bad luck, not a weakness of the method.

**Read the results:** the "Tool comparison" plots show each score against mutations per
sample. Look for where each tool's error levels off: below that burden, the data cannot
support a reliable fit.

---

## Custom mixtures

**Config:** [`custom_mixtures.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/custom_mixtures.yaml)
· **Table:** [`config/scenarios/custom_mixtures.tsv`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/custom_mixtures.tsv)

**Question:** how do the tools handle signature mixes I choose?

Use this to test specific cases: two similar signatures together, one signature at a
very small share, a mix typical of your cohort. Every row of your table becomes one
sample.

**The table** is tab-separated: one row per sample, one column per signature, weights
as shares. Rows are rescaled to add up to 1. Signatures you don't list get 0.

```
sample_id      SBS1   SBS2   SBS5   SBS13  SBS18  SBS40  SBS90
clock_only     0.5    0      0.5    0      0      0      0
apobec_high    0.2    0.3    0.2    0.3    0      0      0
mixed          0.25   0.1    0.25   0.1    0.1    0.2    0
with_unknown   0.3    0      0.3    0      0.1    0.1    0.2
```

**Settings that matter:**

```yaml
data_types: ["custom"]
custom_datasets:
  ExampleMix:                                    # name used in result folders and the report
    combinations_file: "config/scenarios/custom_mixtures.tsv"
    out_of_reference_sigs: ["SBS90"]             # in the table, hidden from the tools
num_muts: [1000, 10000]                          # every row gets this many mutations
num_realizations: 5                              # 5 random draws of each row
```

- `num_samples` is not used: the table sets the samples.
- Hidden signatures for a custom table go **inside** its `custom_datasets` entry. The
  top-level `out_of_reference_sigs` applies only to `realsim`.
- You can list several tables under `custom_datasets`; each gets its own tab in the report.
- `signature_file` (optional) uses another COSMIC version for this table. It is needed
  when the table uses signatures that are only in newer versions (e.g. `SBS40a`,
  `SBS96`, which are in v3.4 but not v3.3.1).

**A larger example:** [`config/config_custom.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/config_custom.yaml)
with the table [`config/custom.tsv`](https://github.com/mkatsanto/SigFitBench/blob/main/config/custom.tsv).
It has 100 samples. Each one has SBS3 and SBS5 plus three other signatures, every
signature at 20%. It is fitted against COSMIC v3.4 and runs together with a `realsim`
run on all 12 tissues.

**Read the results:** `results_custom_mixtures/sigfitbench_report_synthetic.html`. The
"Signature activity per sample" section compares the known mix of each row with what
each tool found.

---

## Own catalog

**Config:** [`own_catalog.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/own_catalog.yaml)
· **Catalog:** [`input/examples/example_catalog.tsv`](https://github.com/mkatsanto/SigFitBench/blob/main/input/examples/example_catalog.tsv)

**Question:** do the tools agree on my samples?

Real samples have no known answer, so the tools are compared with each other. Where
most tools agree, you can be more confident in a signature; where they disagree, the
result depends on the tool you pick.

**The catalog:** mutation counts with the 96 single-base substitution types as rows
(`A[C>A]A`, `A[C>A]C`, ...) and samples as columns. Tab- or comma-separated. Lines
starting with `#` are ignored.

```
Type      Breast::sample_1  Breast::sample_2  Head-SCC::sample_1  ...
A[C>A]A   102               66                214
A[C>A]C   73                54                103
```

Columns named `Tissue::Sample` are split into one group per tissue. Columns without
`::` are fitted together as one group named after the dataset.

**Settings that matter:**

```yaml
data_types: ["real"]
real_datasets:
  Example:
    catalog: "input/examples/example_catalog.tsv"
num_muts: []                 # not used: the catalog sets the counts
cosmic_signatures:
  base_file: "input/COSMIC_v3.3.1_SBS_GRCh38.txt"   # match your data's genome build
```

`input/` has COSMIC files for GRCh37 and GRCh38. Use the one matching the genome your
mutations were called on.

**Read the results:** `results_own_catalog/sigfitbench_report_real_samples.html` shows,
per tissue, how similar each pair of tools' results is, which signatures each tool
reports, and what share of mutations each tool leaves unassigned. See
[Results and metrics](results.md#comparing-tools-without-a-known-answer).

---

## Real PCAWG

**Config:** [`real_pcawg.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/real_pcawg.yaml)

**Question:** do the tools agree on real tumour genomes?

The same data type as [Own catalog](#own-catalog) (`real`), but with whole genomes from
PCAWG, split by tumour type, each at its own mutation count.

**Before you run it:** the PCAWG catalog is not in this repository. Put it at
`real_mutational_catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat`, or change
`catalog` to where yours is. Column names must use the `Tissue::Sample` form.

**Settings that matter:**

```yaml
data_types: ["real"]
real_datasets:
  PCAWG:
    catalog: "real_mutational_catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat"
    cancer_types: ["Breast-AdenoCA", "Skin-Melanoma"]   # leave out to fit every tumour type
```

PCAWG tumour-type names (`Breast-AdenoCA`) differ from the tissue names used for
simulation (`Breast`), so they have their own `cancer_types` list inside the
`real_datasets` entry, separate from the top-level `cancer_types`.

---

## Full benchmark

**Config:** [`full_benchmark.yaml`](https://github.com/mkatsanto/SigFitBench/blob/main/config/scenarios/full_benchmark.yaml)

**Question:** which tool should I use, for which tissue and mutation burden?

All 14 tools, all 12 tissues, 100 to 200,000 mutations per sample, 100 samples, one
hidden signature carrying 10% of mutations. This is the scale of the published
comparison. It also lists every tool setting with its default value.

Run it on a cluster:

```bash
snakemake --profile profiles/slurm-generic \
  --snakefile workflow/Snakefile --configfile config/scenarios/full_benchmark.yaml
```

Fill in your partition and account in `profiles/slurm-generic/config.yaml` first.
`mmsig`, `mSigAct` and `StarSignDNA` take most of the time; leave them out of
`tools_to_run` for a faster first pass.

---

## Combining scenarios

One config can hold several data types. For example, to score tools on simulated data
and then see whether they agree on your own samples, in a single run:

```yaml
data_types: ["realsim", "real"]
cancer_types: ["Breast"]
num_muts: [2000, 20000]
num_samples: 20
real_datasets:
  MyCohort:
    catalog: "input/my_cohort.tsv"
```

Each data type gets its own report in the same `output_dir`.
