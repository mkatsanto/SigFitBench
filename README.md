# SigFitBench

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.20-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Snakemake workflow that compares mutational signature fitting tools. It runs 14
published tools on the same mutation catalogs, scores every result against the same
known answer, and writes one interactive HTML report.

Many fitting tools exist, and they disagree: some report signatures that are not
there, others miss ones that are, and how well each does depends on the number of
mutations, the tumour type and how similar the active signatures are. SigFitBench
follows the approach of [Medo et al., 2024](https://doi.org/10.1038/s41467-024-53711-6)
([SigFitTest](https://github.com/8medom/SigFitTest)):

1. **Simulate** samples with a known signature mix, based on COSMIC signature activity in real tumours
2. **Fit** them with every selected tool, using the same input and the same reference signatures
3. **Score** each tool (error, precision, recall, F1, MCC)
4. **Report** where each tool works well and where it does not

📖 **Docs:** [Guide](docs/guide.md) · [Configuration reference](docs/configuration.md) · [Results and metrics](docs/results.md) · [Add your own tool](docs/add_a_tool.md)

---

## Contents

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [Scenarios](#scenarios)
4. [Fitting tools](#fitting-tools)
5. [Configuration](#configuration)
6. [Results](#results)
7. [Running on a cluster](#running-on-a-cluster)
8. [Troubleshooting](#troubleshooting)
9. [Citation](#citation)

---

## Installation

**Needs:** Linux or macOS, conda ([Miniconda](https://docs.anaconda.com/miniconda/) is enough), about 10 GB of disk and 8 GB of RAM.

```bash
git clone https://github.com/mkatsanto/SigFitBench.git
cd SigFitBench
conda env create -f install/environment.yml
conda activate sigfitbench
```

Each tool gets its own conda environment, created automatically the first time it runs
(`--use-conda`). This first run takes 30–60 minutes; later runs reuse the environments.

---

## Quick start

```bash
conda activate sigfitbench

# Preview what will run
snakemake --snakefile workflow/Snakefile --configfile config/scenarios/quick_test.yaml -n

# Run it
snakemake --profile profiles/local-conda \
  --snakefile workflow/Snakefile --configfile config/scenarios/quick_test.yaml
```

Then open `results_quick_test/report_synth_realistic.html` in a browser.

The `local-conda` profile turns on `--use-conda`, `--rerun-incomplete` and `--keep-going`,
and uses 12 cores. Change the core count with `--cores N`.

There is no default config: always pass `--configfile`.

---

## Scenarios

Each file in [`config/scenarios/`](config/scenarios/) is a complete, working config for one
question. Copy the one closest to yours and change it. The
[guide](docs/guide.md) explains each one: what it tests, which settings matter, and how
to read the results.

| Scenario | Config | Question it answers | Data type |
|---|---|---|---|
| Quick test | [`quick_test.yaml`](config/scenarios/quick_test.yaml) | Does the installation work? | `synth_realistic` |
| Mutation burden | [`mutation_burden.yaml`](config/scenarios/mutation_burden.yaml) | How does accuracy change from 100 to 50,000 mutations, and how much is chance? | `synth_realistic` |
| Custom mixtures | [`custom_mixtures.yaml`](config/scenarios/custom_mixtures.yaml) | How do tools handle signature mixes I choose? | `synth_custom` |
| Own catalog | [`own_catalog.yaml`](config/scenarios/own_catalog.yaml) | Do the tools agree on my samples? | `real` |
| Real PCAWG | [`real_pcawg.yaml`](config/scenarios/real_pcawg.yaml) | Do the tools agree on real tumours? | `real` (PCAWG file not included) |
| Full benchmark | [`full_benchmark.yaml`](config/scenarios/full_benchmark.yaml) | All tools, all tissues: for a cluster | `synth_realistic` |
| Unknown signatures | [`unknown_signatures.yaml`](config/scenarios/unknown_signatures.yaml) | What happens when tumours carry signatures the tools don't know? | `synth_realistic` |

`synth_realistic` and `synth_custom` are simulated, so they have a known answer and tools are scored
against it. `real` data has none, so tools are compared with each other instead. See
[Data types](#configuration) below.

---

## Fitting tools

Put any of these names in `tools_to_run`, spelled exactly as shown (the name is used in
the output file names). Most names are the package name; `SPA` is short for
SigProfilerAssignment. Always list the tools you want: `"all"` also selects methods
that are still in development and not yet released.

| `tools_to_run` | Package | How it fits | Settings (`tool_overrides`) |
|---|---|---|---|
| `MutationalPatterns` | [MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns/) | Non-negative least squares; with `strict`, removes signatures one by one while the fit stays good | `mutationalpatterns.strict` (**required**) |
| `deconstructSigs` | [deconstructSigs](https://github.com/raerose01/deconstructSigs) | Repeated regression, adding the signature that improves the fit most each step | — |
| `SPA` | [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment) | Non-negative least squares with forward and backward selection of signatures | — |
| `SPA_low` | SigProfilerAssignment | As `SPA`, with a lower removal penalty (0.02 instead of 0.05), so it keeps more signatures | — |
| `MuSiCal` | [MuSiCal](https://github.com/parklab/MuSiCal) | Likelihood-based refitting that adds and removes signatures | — |
| `signature_tools` | [signature.tools.lib](https://github.com/Nik-Zainal-Group/signature.tools.lib) | Fit repeated on bootstrap resamples; signatures with uncertain exposure are dropped | `method`, `exposureFilterType`, `nboot`, `nparallel` |
| `sigminer` | [sigminer](https://github.com/ShixiangWang/sigminer) | Least squares, quadratic programming or simulated annealing | `sigminer.method` (**required**) |
| `sigfit` | [sigfit](https://github.com/kgori/sigfit) | Bayesian fit sampled with MCMC (Stan) | — |
| `mmsig` | [mmsig](https://github.com/evenrus/mmsig) | Expectation–maximisation, removing signatures that do not change the profile | — |
| `siglasso` | [sigLASSO](https://github.com/gersteinlab/siglasso) | L1-penalised fit that also models the sampling of mutations | — |
| `sigspack` | [SigsPack](https://bioconductor.org/packages/SigsPack/) | Non-negative least squares (quadratic programming) | — |
| `yapsa` | [YAPSA](https://bioconductor.org/packages/YAPSA/) | Non-negative least squares (YAPSA's `LCD`) | — |
| `mSigAct` | [mSigAct](https://github.com/steverozen/mSigAct) | Maximum likelihood, then removes signatures while a likelihood-ratio test allows it | `max_level`, `p_thresh`, `likelihood`, `threads`, `extra` |
| `StarSignDNA` | [StarSignDNA](https://github.com/uio-bmi/StarSignDNA) | Poisson likelihood with an L1 penalty; reports the median over bootstrap runs | `n_iterations` |

Settings go under `tool_overrides`, keyed by the tool name in lower case:

```yaml
tool_overrides:
  mutationalpatterns:
    strict: true          # required whenever MutationalPatterns runs
  sigminer:
    method: "NNLS"        # required whenever sigminer runs: "NNLS", "QP" or "SA"
  msigact:
    max_level: 2
```

[`full_benchmark.yaml`](config/scenarios/full_benchmark.yaml) lists every setting with its
default.

**Run time.** Measured on a laptop, one job (two samples, 1,000–500,000 mutations each):

| Speed | Tools | Time per job |
|---|---|---|
| Fast | `yapsa`, `sigspack`, `MuSiCal`, `SPA`, `SPA_low`, `sigminer`, `siglasso`, `MutationalPatterns` | under 30 s |
| Moderate | `signature_tools`, `sigfit`, `deconstructSigs` | 30–75 s |
| Slow | `mmsig`, `mSigAct` | 25–40 min |
| Very slow | `StarSignDNA` | over 2 h |

Time grows with the number of samples. Run the slow tools on a cluster.

> **Install problems.** Some tool environments are harder to build than others,
> especially on Apple Silicon (`mmsig`, `sigfit`). See the tool's `workflow/envs/*.yaml`
> and `*.post-deploy.sh` for how it is installed, and
> [open an issue](https://github.com/mkatsanto/SigFitBench/issues) if it fails.

---

## Configuration

The most used settings are below. The
[configuration reference](docs/configuration.md) lists all of them.

```yaml
tools_to_run: ["MutationalPatterns", "SPA", "MuSiCal"]

data_types: ["synth_realistic"]          # "synth_realistic", "synth_custom", "real"
relevance_modes: ["standard"]

cancer_types: ["Head-SCC"]       # tissue profiles in cosmic_tissue_data/
num_muts: [1000, 10000]          # mutations per sample; one run per value
num_samples: 30                  # different signature mixes per tissue
num_realizations: 1              # repeats of the dataset with new random mutations

out_of_reference_sigs: ["SBS86", "SBS90", "SBS92"]   # unknown signatures: not given to the tools
out_of_reference_weights: [0.1]  # one unknown signature per value, with that share of mutations

output_dir: "results_my_run"
log_dir: "logs_my_run"
random_seed: 0
```

**Data types.** Where the mutation catalogs come from:

| `data_types` | Source | Known answer? |
|---|---|---|
| `synth_realistic` | Realistic simulations: each sample copies the signature mix of a real tumour (PCAWG for WGS, TCGA for WES) of a `cancer_types` tissue, then mutations are drawn at random from that mix | Yes |
| `synth_custom` | Simulated from your own table of signature mixes (`custom_datasets`) | Yes |
| `real` | Real samples: PCAWG genomes or your own catalog (`real_datasets`) | No: tools compared with each other |

Each data type gets its own report. The [guide](docs/guide.md#where-the-samples-come-from)
explains step by step how `synth_realistic` samples are made.

**`num_samples` vs. `num_realizations`.** `num_samples` sets how many *different* signature
mixes are used. `num_realizations` repeats the whole dataset: each repeat keeps the same
mixes (the same known answer) and draws new random mutations, and is scored on its own.
Use more realizations to tell a tool's real errors from bad luck, which matters most when
samples have few mutations.

**Reference signatures.** All tools fit against the same list: the COSMIC signatures in
`cosmic_signatures.base_file` (v3.3.1 in the scenarios), minus any unknown signatures.
COSMIC files up to v3.6 are in `input/`; see the
[configuration reference](docs/configuration.md#files-and-folders) for which version
fits which data type.

**Unknown signatures.** To test how tools cope with signatures they don't know, the
signatures in `out_of_reference_sigs` are taken out of the tools' list but still used to
simulate samples. Each value in `out_of_reference_weights` adds one of them to every
sample, with that share of its mutations: `[0.1, 0.1]` adds two, 10% each. One run tests
one level; to compare 0%, 10% and 20%, run three copies of the config with their own
`output_dir`. The [guide](docs/guide.md#unknown-signatures) explains this in detail.

---

## Results

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

| Score | Range | Better | Meaning |
|---|---|---|---|
| **TAE** | 0–1 | lower | Total absolute error: share of mutations assigned to the wrong signatures |
| **Precision** | 0–1 | higher | Share of reported signatures that are really there |
| **Recall** | 0–1 | higher | Share of real signatures that were found |
| **F1** | 0–1 | higher | Balance of precision and recall |
| **MCC** | −1–1 | higher | Matthews correlation: overall agreement on which signatures are active |

A signature counts as found when a tool gives it at least `evaluation.min_contribution`
mutations (default 10). See [Results and metrics](docs/results.md) for how to read the
report.

---

## Running on a cluster

[`profiles/slurm-generic`](profiles/slurm-generic/config.yaml) uses Snakemake's SLURM
executor (installed with the base environment). Fill in your partition and account under
`default-resources`, then:

```bash
snakemake --profile profiles/slurm-generic \
  --snakefile workflow/Snakefile --configfile config/scenarios/full_benchmark.yaml
```

- Time limits per tool are under `set-resources` in the profile.
- Up to 50 jobs are submitted at once; change with `--jobs N`.
- SLURM job output goes to `.snakemake/slurm_logs/`; tool logs go to `log_dir`.

[`profiles/slurm-conda`](profiles/slurm-conda/config.yaml) is the same profile filled in for
UBELIX (University of Bern), as an example.

**Workflow diagram.** `profiles/graphs` draws the job graph:

```bash
snakemake --profile profiles/graphs --configfile config/scenarios/quick_test.yaml \
  --rulegraph | dot -Tpdf > workflow_rules.pdf
```

---

## Troubleshooting

**A tool failed.** The run continues; the tool's result file
(`runs/.../fits/{tool}_rep{r}.dat`) starts with `# TOOL FAILED:` and the last error line.
The full log is in the same run folder under `{log_dir}`, e.g.
`{log_dir}/runs/synth_realistic/Head-SCC/1000_mutations/MuSiCal_rep0.stderr.log`.

**A tool environment is broken.** Delete it and rebuild only the environments:

```bash
rm -rf .snakemake/conda/<hash>
snakemake --snakefile workflow/Snakefile --configfile config/scenarios/quick_test.yaml \
  --use-conda --conda-create-envs-only --cores 2
```

**The run was interrupted.** Unlock the folder, then start again:

```bash
snakemake --snakefile workflow/Snakefile --configfile <your config> --unlock
```

**Out of memory.** Use fewer cores (`--cores 2`) or fewer samples (`num_samples`).

**`KeyError: 'strict'` or `KeyError: 'method'`.** Add `mutationalpatterns.strict` or
`sigminer.method` under `tool_overrides`.

Questions: [GitHub issues](https://github.com/mkatsanto/SigFitBench/issues) ·
maria.katsantoni@unibe.ch, matus.medo@unibe.ch, qixuan.wang@unibe.ch

---

## Citation

If you use SigFitBench, please cite:

> Medo, M., Ng, C. K. Y., Medová, M. A comprehensive comparison of tools for fitting
> mutational signatures. *Nat Commun* 15, 9467 (2024).
> [doi:10.1038/s41467-024-53711-6](https://doi.org/10.1038/s41467-024-53711-6)

> Mölder, F. et al. Sustainable data analysis with Snakemake. *F1000Res* 10, 33 (2021).
> [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)

Please also cite the tools whose results you use:

| Tool | Citation |
|---|---|
| `MutationalPatterns` | Blokzijl, F. et al. MutationalPatterns: comprehensive genome-wide analysis of mutational processes. *Genome Med* 10, 33 (2018). [doi:10.1186/s13073-018-0539-0](https://doi.org/10.1186/s13073-018-0539-0) |
| `deconstructSigs` | Rosenthal, R. et al. deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. *Genome Biol* 17, 31 (2016). [doi:10.1186/s13059-016-0893-4](https://doi.org/10.1186/s13059-016-0893-4) |
| `SPA`, `SPA_low` | Díaz-Gay, M. et al. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. *Bioinformatics* 39, btad756 (2023). [doi:10.1093/bioinformatics/btad756](https://doi.org/10.1093/bioinformatics/btad756) |
| `MuSiCal` | Jin, H. et al. Accurate and sensitive mutational signature analysis with MuSiCal. *Nat Genet* 56, 541–552 (2024). [doi:10.1038/s41588-024-01659-0](https://doi.org/10.1038/s41588-024-01659-0) |
| `signature_tools` | Degasperi, A. et al. A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies. *Nat Cancer* 1, 249–263 (2020). [doi:10.1038/s43018-020-0027-5](https://doi.org/10.1038/s43018-020-0027-5) |
| `sigminer` | Wang, S. et al. Copy number signature analysis tool and its application in prostate cancer reveals distinct mutational processes and clinical outcomes. *PLoS Genet* 17, e1009557 (2021). [doi:10.1371/journal.pgen.1009557](https://doi.org/10.1371/journal.pgen.1009557) |
| `sigfit` | Gori, K., Baez-Ortega, A. sigfit: flexible Bayesian inference of mutational signatures. *bioRxiv* 372896 (2018). [doi:10.1101/372896](https://doi.org/10.1101/372896) |
| `mmsig` | Rustad, E. H. et al. mmsig: a fitting approach to accurately identify somatic mutational signatures in hematological malignancies. *Commun Biol* 4, 424 (2021). [doi:10.1038/s42003-021-01938-0](https://doi.org/10.1038/s42003-021-01938-0) |
| `siglasso` | Li, S., Crawford, F. W., Gerstein, M. B. Using sigLASSO to optimize cancer mutation signatures jointly with sampling likelihood. *Nat Commun* 11, 3575 (2020). [doi:10.1038/s41467-020-17388-x](https://doi.org/10.1038/s41467-020-17388-x) |
| `sigspack` | Schumann, F. et al. SigsPack, a package for cancer mutational signatures. *BMC Bioinformatics* 20, 450 (2019). [doi:10.1186/s12859-019-3043-7](https://doi.org/10.1186/s12859-019-3043-7) |
| `yapsa` | Hübschmann, D. et al. Analysis of mutational signatures with yet another package for signature analysis. *Genes Chromosomes Cancer* 60, 314–331 (2021). [doi:10.1002/gcc.22918](https://doi.org/10.1002/gcc.22918) |
| `mSigAct` | Ng, A. W. T. et al. Aristolochic acids and their derivatives are widely implicated in liver cancers in Taiwan and throughout Asia. *Sci Transl Med* 9, eaan6446 (2017). [doi:10.1126/scitranslmed.aan6446](https://doi.org/10.1126/scitranslmed.aan6446) |
| `StarSignDNA` | See [github.com/uio-bmi/StarSignDNA](https://github.com/uio-bmi/StarSignDNA) for the current citation. |

---

## License

MIT, see [LICENSE](LICENSE).
