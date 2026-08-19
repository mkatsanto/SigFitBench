# SigFitTest Workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake-based workflow for comprehensive benchmarking and evaluation of mutational signature fitting tools — synthetic and real data, a dozen third-party tools, and two in-house methods aimed at the unknown-signature problem, all evaluated side by side.

## Introduction

### The Problem

Cancer genomes accumulate mutations through various mutational processes (e.g., UV exposure, smoking, DNA repair defects). Each process leaves a characteristic "mutational signature" — a pattern of mutations across different sequence contexts. Identifying which signatures are active in a tumor helps understand its etiology and guide treatment.

**However**, dozens of signature fitting tools exist and their accuracy varies significantly:
- Some tools detect false signatures (false positives)
- Others miss active signatures (false negatives)
- Performance depends on mutation burden, cancer type, and signature similarity
- No consensus on which tool to use for which scenario

### The Solution

This workflow implements the **SigFitTest** methodology ([Medo et al., 2024, *Nat Commun*](https://doi.org/10.1038/s41467-024-53711-6)) — a systematic benchmarking framework that:

1. **Generates realistic synthetic data** with known ground truth from COSMIC tissue distributions
2. **Runs signature fitting tools** in parallel with identical inputs
3. **Evaluates performance** using comprehensive metrics (TAE, Precision, Recall, F1, MCC)
4. **Identifies best tools** for different scenarios (mutation counts, cancer types, signature sets)
5. **Tests robustness** when some signatures are unknown to tools

**Original Repository**: [SigFitTest](https://github.com/8medom/SigFitTest)

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Configuration](#configuration)
4. [Method Descriptions](#method-descriptions)
5. [Understanding Results](#understanding-results)
6. [Advanced Usage](#advanced-usage)
7. [Troubleshooting](#troubleshooting)
8. [Performance Expectations](#performance-expectations)
9. [Citation](#citation)

---

## Installation

**Requirements**: Linux/macOS, 10 GB disk space, 8 GB RAM

### 1. Install Conda (if needed)

```bash
# macOS (Apple Silicon)
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh

# macOS (Intel)
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Linux
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Restart terminal after installation.

### 2. Create Environment

```bash
git clone https://github.com/mkatsanto/sigfittest_workflow.git
cd sigfittest_workflow

# Create base environment (10-20 min)
conda env create -f install/environment.yml

# Activate
conda activate sigfittest
```

Tool-specific conda environments are created automatically on first run via the `--use-conda` flag.

---

## Quick Start

### Local Execution

**Option A — plain flags:**

```bash
conda activate sigfittest
cd sigfittest_workflow

# Run with a specific config file
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --cores 8 --use-conda --rerun-incomplete
```

Replace `config_test.yaml` with any config from the `config/` directory (see [Configuration](#configuration) below).

**Useful flags**:
- `--cores N` — number of parallel jobs
- `--use-conda` — auto-create tool-specific conda environments (required on first run)
- `--rerun-incomplete` — resume after failures
- `--dry-run` / `-n` — preview jobs without executing
- `-p` — print shell commands as they run

**Option B — `local-conda` profile:**

`profiles/local-conda` bundles the same flags (`use-conda: true`, `rerun-incomplete: true`,
`cores: 3`) into a single reusable profile, so you don't have to repeat them every time:

```bash
conda activate sigfittest
cd sigfittest_workflow

snakemake --profile profiles/local-conda \
  --snakefile workflow/Snakefile --configfile config/config_test.yaml
```

Override `cores: 3` from the profile with `--cores N` on the command line if you want more
parallelism than the default.

**First run**: conda environments are created automatically (30–60 min). Subsequent runs skip this step.

### SLURM Cluster Execution (UBELIX)

`profiles/slurm-conda` is a ready-to-use Snakemake profile for UBELIX (University of Bern HPC). It
submits jobs via `sbatch` (`profiles/slurm-submit.py`), polls status via `profiles/slurm-status.py`,
and already sets `use-conda: true` and `rerun-incomplete: true`, so those flags don't need to be
repeated on the command line.

```bash
conda activate sigfittest
cd sigfittest_workflow

snakemake --profile profiles/slurm-conda \
  --snakefile workflow/Snakefile --configfile config/config_test.yaml
```

- **Account/QoS/partition**: pre-configured in `profiles/slurm-config.json`'s `__default__` block
  (`account: gratis`, `qos: job_gratis`, `partition: epyc2`) — edit these if your UBELIX allocation
  differs.
- **Per-rule time limits**: also in `profiles/slurm-config.json`, keyed by rule name (e.g. `extended`,
  `sigfit`, `assess_fit_quality_bootstrap`) and tuned to each rule's expected cost; add/adjust an entry
  there instead of passing `--time` on the command line.
- **Parallel jobs**: the profile submits up to `jobs: 10` at once by default — override with `--jobs N`.
- **Per-job SLURM stdout/stderr**: every config file must set `cluster_log_dir` (e.g.
  `cluster_log_dir: "logs_test/slurm"`), which every rule passes through as
  `params.cluster_log_path` for `slurm-config.json`'s `output` template
  (`{params.cluster_log_path}/{rule}.{wildcards}-%j-%N.out`). This is separate from each rule's own
  `log:` files under `log_dir`, which already work the same locally or on SLURM.

---

## Configuration

All config files live in `config/`. There is no default `config.yaml` — you must always pass `--configfile`.

### Available Configs

| File | Purpose |
|------|---------|
| `config_unknown_0.yaml` | Benchmark run: all 16 tools, 3 tissues, 0% unknown-signature fraction |
| `config_unknown_10.yaml` | Same, at 10% unknown-signature fraction |
| `config_unknown_20.yaml` | Same, at 20% unknown-signature fraction |
| `config_unknown_30.yaml` | Same, at 30% unknown-signature fraction |
| `config_mysigmix.yaml` | Example `"custom"` data type run, drawing signature mixes from a user-supplied combinations file |

All configs share the same structure (see [Key Parameters](#key-parameters) below) — copy one and
adjust `tools_to_run`, `num_muts`, `cancer_types`, `num_samples`, and `out_of_reference_weights` for
your own run. Runtime scales with `num_samples`, `num_muts`, `num_realizations`, and how many
tools/penalty values are selected — see [Performance Expectations](#performance-expectations).

### Key Parameters

```yaml
# Tools to run
tools_to_run: ["SPA", "MuSiCal", "SigID", "SigTrim"]
# Use "all" to run every tool below, or a list of specific names — casing must match the
# names in the Available Tools table below (it's used verbatim in output filenames).

# Simulation parameters
num_muts: [100000, 1000000]       # Mutations per sample
num_samples: 30                    # Samples per cohort
num_realizations: 1                # Independent runs
cancer_types: ["Head-SCC"]         # Cancer types (must match cosmic_tissue_data/)
data_types: ["synthetic"]          # "synthetic", "realsim", "external", "custom" — see below
relevance_modes: ["standard"]      # "standard", "relevant", or both

# Unknown signature testing (set weights to [] to disable)
out_of_reference_sigs: ["SBS86", "SBS90", "SBS92"]
out_of_reference_weights: [0.1, 0.2]   # Fraction of mutations from unknown sigs

# Output
output_dir: "results_unknown_10"
log_dir: "logs_unknown_10"
cluster_log_dir: "logs_unknown_10/slurm"     # SLURM job stdout/stderr — required by every config, see SLURM section above

# Reproducibility
random_seed: 0
```

### Available Tools

`tools_to_run` accepts `"all"` or any of the following names, grouped by the input format they need
(see `tool_formats` in each config — this is how the workflow decides which COSMIC-signature/catalog
file variant to generate for each tool). The exact casing shown here is required — it's used verbatim
in each tool's `{tool}_contribution.dat` output filename:

| Group | Tools |
|---|---|
| `tab_with_index` | `deconstructSigs`, `SPA`, `SPA_low`, `FSPA` |
| `csv_no_index` | `MutationalPatterns`, `MuSiCal`, `SigID`, `SigTrim` |
| `csv_with_index` | `sigfit`, `mmsig`, `sigminer`, `sigspack`, `yapsa` |
| `csv_quoted` | `siglasso` |
| `signature_tools` | `signature_tools` |

`SigID` and `SigTrim` are this workflow's own methods for the unknown-signature problem (see
[Method Descriptions](#method-descriptions)); the rest are established third-party tools. Neither
sweeps a penalty grid — both always emit a single `{tool}_contribution.dat` per condition.

Add per-tool overrides (thresholds, penalties, `metric: cosine`/`l2`, etc.) under `tool_overrides` —
see the blocks already in `config_unknown_0.yaml` for every tool's available knobs.

### Data Types

`data_types` controls where each condition's mutation catalog comes from. Combine any of these in one
run (each activates a different Snakemake rule):

| `data_types` value | Rule | Source |
|---|---|---|
| `"synthetic"` | `ms_create_synthetic_catalog` | Samples signature weights from empirical COSMIC tissue distributions (`cosmic_tissue_data/`), then draws mutations via multinomial sampling. |
| `"realsim"` | `create_realsim_catalog` | Real PCAWG samples for selected cancer types, at their native mutation count — no downsampling, no ground truth; `rep` is always 0. |
| `"external"` | `prepare_external_catalog` | Brings your own catalog + ground truth via `external_datasets` (see the commented example at the bottom of any config); `rep` is always 0. |
| `"custom"` | `create_custom_catalog` | Draws from a user-supplied `combinations_file` under `custom_datasets` (see `config_mysigmix.yaml` for a working example) — each sample randomly picks one row (signature mix) from the table. |

For `"realsim"`, cancer types come from `simulation.real_samples.cancer_types` (a separate list from
the top-level `cancer_types`, since real PCAWG labels don't match the synthetic ones, e.g.
`Breast-AdenoCA` vs. `Breast`). Since there's no ground truth, `SigID` is skipped for this mode (its
auto-penalty search requires ground truth as an input), and results are compared via cross-tool
agreement rather than TAE/Precision/Recall:

```yaml
data_types: ["realsim"]
simulation:
  real_samples:
    input_real_catalog: "real_mutational_catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat"
    cancer_types: ["Breast-AdenoCA", "Skin-Melanoma"]   # see config comments for the full list available
```

For `"custom"`, each key under `custom_datasets` becomes a `cancer_type` in the output paths:

```yaml
data_types: ["custom"]
custom_datasets:
  MySigMix:
    combinations_file: "input/my_combinations.dat"   # TSV: rows=combinations, cols=signatures
    out_of_reference_sigs: ["SBS90"]                  # optional per-dataset override
    out_of_reference_weights: [0.1]                   # optional per-dataset override
```

---

## Method Descriptions

Every fitting tool answers the same question — *which COSMIC signatures explain this
sample's mutation catalog, and how much does each contribute?* This workflow's own
methods (`workflow/scripts/`) go one step further and also try to answer *what if the
sample contains a signature that isn't in the reference at all?* That's the
**unknown-signature problem**, and it's the reason these custom methods exist alongside
the third-party tools above.

Two general approaches tackle it, described below at the level of "how it works" rather
than a fixed recipe — the codebase currently carries several experimental variants of
each (different metrics, protection rules, penalty-selection strategies) while the best
configuration is still being benchmarked. Once that settles, the two winning versions
will replace this placeholder with full step-by-step writeups.

### What is LOO (Leave-One-Out)?

Used throughout all methods for signature removal decisions. For a set of N active signatures:

```
current_fit = metric(x, W_active @ NNLS(W_active, x))

for each sig i in active set:
    fit NNLS on (active set minus sig i)
    cost_i = current_fit - metric(x, W_minus_i @ NNLS(W_minus_i, x))

remove the sig with the lowest cost_i if cost_i < threshold
```

The cost of removing sig i is how much the fit degrades without it. A near-zero cost means another sig is covering for it — strong signal of a false positive or correlated duplicate. LOO is repeated iteratively (one removal per round) because removing one sig changes the NNLS weights of all others.

Metrics used:
- **cosine**: `dot(x, Wh) / (‖x‖ ‖Wh‖)` — higher is better; drop on removal is the cost.
- **L2**: `‖x − Wh‖ / ‖x‖` — lower is better; rise on removal is the cost.

---

### SPA (`cosmic_fit`) — Reference Algorithm

**No unknowns.** The reference against which all extended methods are compared.

1. **NNLS on all signatures** → sparse active set (typically 10–20 non-zero weights). This is the starting point — not all 78.
2. **Initial backward elimination** (`initial_remove_penalty = 0.05`): LOO cosine loop on the NNLS-sparse set. Remove the sig with the smallest LOO drop if drop < 0.05. Repeat until convergence.
3. **Add-remove loop** (`add_penalty = 0.05`, `remove_penalty = 0.01`): for each sig not in active set, tentatively add it, run LOO removal at 0.01, score the result. Accept the best candidate only if its improvement > 0.05. Repeat until no candidate clears the bar.
4. **SBS1/SBS5 protection**: always present, protected from removal at every step.
5. **Connected family protection** (APOBEC SBS2/13, UV SBS7a–d, POLE SBS10a/b, SBS17a/b): when any family member is active, all are protected.
6. **Final output**: NNLS weights on the active known set. No SLSQP, no penalty sweep.

**Key design choices:**
- Starts from NNLS-sparse set (not all sigs simultaneously) → LOO costs are honest from the start.
- `add_penalty (0.05) >> remove_penalty (0.01)` → hysteresis gap. Prevents oscillation at the boundary.

---

### `iterative` — lightweight pruning approach

**Objective**: standard NNLS, no unknown-signature penalty.

**Steps:**

1. **Channel-level pruning** (z-score based, `zscore_threshold`): fit NNLS, compute per-channel residual, and drop channels whose residual z-score exceeds the threshold — they likely carry out-of-reference signal that would distort the known-signature fit. Repeat until no anomalous channels remain.
2. **Backward elimination** (`thresh_backward`): LOO loop on the pruned channel set until convergence.
3. **Final NNLS** on surviving signatures over the original channels.

No SBS1/SBS5 protection, no add step, no explicit unknown-signature output — it simply
keeps unexplained signal from distorting the known-signature fit rather than modelling
it directly. `iterative_version1` is identical except it uses **L2 distance** instead of
cosine for the backward-elimination step.

---

### `extended` — residual-driven unknown detection

**Objective**: `min_{h≥0} 0.5·‖Wh−x‖² + λ·Σ|h_k|` (k ∈ unknown sigs), minimised with SLSQP.

**Steps:**

1. **Backward elimination** on known signatures (LOO, cosine or L2 depending on variant).
2. **Unknown spike construction**: take the NNLS residual, find the channels with the
   largest unexplained signal, and build one identity "spike" signature per channel.
3. **Combine spikes** into a single unknown profile via SLSQP-weighted averaging.
4. **Second backward elimination** on known signatures + the combined unknown.
5. **Final SLSQP with an L1 penalty** on the unknown weight, so it only survives in the
   fit if it earns its place.

This is the family under active development. The codebase currently has several
experimental variants exploring: different LOO metrics (cosine vs. L2), SBS1/SBS5 and
connected-family protection, greedy add/remove search instead of pure pruning, automatic
vs. swept penalty selection, and bias corrections (e.g. a `dummy_sig_multiplier`) for the
unfair LOO advantage unknown signatures get from being built directly out of the
residual. Two of these variants will be promoted here with full writeups once
benchmarking confirms which perform best.

---

### At a Glance

| Method | Unknowns modelled | Status |
|---|---|---|
| SPA | ✗ | Reference algorithm — all custom methods are benchmarked against it |
| `iterative` (+ `_version1`) | ✗ | Stable — channel pruning + backward elimination (cosine / L2) |
| `extended` | ✓ (L1 penalty) | Under active development — several variants in evaluation |

---

## Understanding Results

### Output Structure

```
{output_dir}/
├── synthetic/                              # Tool fitting results, organised by condition
│   └── {cancer_type}/
│       └── {relevance}/
│           └── {num_muts}/
│               └── rep_{rep}/
│                   ├── {tool}_contribution.dat          # Fitted signature weights
│                   └── {tool}_contribution_penalty_{p}.dat  # Penalty sweep tools
│
├── evaluation/                             # Per-run evaluation files
│   └── *_per_sample.tsv
│
├── generated_data/                         # Synthetic catalogs and ground truth
│
├── all_evaluation_metrics.tsv              # Aggregated metrics across all tools and conditions
├── all_per_sample_metrics.tsv             # Per-sample metrics with features
├── sigfittest_report.html                 # Interactive HTML report (main output)
├── penalty_analysis_report.html           # TAE vs penalty sweep report
└── workflow_complete.txt                  # Written when workflow finishes successfully
```

### Key Metrics

| Metric | Range | Best | Description |
|--------|-------|------|-------------|
| **TAE** | [0, 1] | Lower | Total Absolute Error: sum(\|fitted − true\|) / (2 × muts) |
| **Precision** | [0, 1] | Higher | TP / (TP + FP) — few false positives |
| **Recall** | [0, 1] | Higher | TP / (TP + FN) — few missed signatures |
| **F1** | [0, 1] | Higher | Harmonic mean of Precision and Recall |
| **MCC** | [-1, 1] | Higher | Matthews Correlation — robust for imbalanced data |

### Interactive HTML Report

Open the report in a browser after the workflow completes:

```bash
open {output_dir}/sigfittest_report.html        # macOS
xdg-open {output_dir}/sigfittest_report.html    # Linux
```

The report contains:
- **Tool comparison**: precision-recall trade-offs, performance vs mutation burden
- **Per-sample analysis**: scatter plots per metric
- **Algorithm stability**: box plots of metric distributions
- **Penalty analysis**: TAE vs penalty value for penalty sweep tools

---

## Advanced Usage

### Dry Run (Preview)

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml --dry-run
```

### Unlock After Interruption

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml --unlock
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --cores 8 --use-conda --rerun-incomplete
```

### Workflow Visualisation

![Rule graph](docs/rulegraph.svg)

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --dag | dot -Tpdf > workflow_dag.pdf

snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --rulegraph | dot -Tpdf > workflow_rules.pdf
```

### Two-Pass Relevant Signature Detection

Auto-detect active signatures from a standard run, then refit with the reduced set:

```yaml
relevance_modes: ["standard", "relevant"]
relevant_signatures: "auto"
relevant_detection:
  min_contribution: 10.0       # Min mutations per sample
  min_sample_fraction: 0.3     # Active in ≥30% of samples
  consensus_tools: 2           # Detected by ≥2 tools
```

### Unknown Signature Testing

Simulate scenarios where tools lack some active signatures:

```yaml
out_of_reference_sigs: ["SBS86", "SBS90", "SBS92"]
out_of_reference_weights: [0.1, 0.2]   # 10% and 20% of mutations from unknowns
```

### Bootstrap Fit-Quality Assessment (under development)

Independently of ground-truth comparison, validate each tool's fit by resampling mutations from its
own reconstructed catalog and refitting `num_bootstrap` times, then comparing reconstruction-metric
distributions against the original fit (`assess_fit_quality_bootstrap` rule → `*_bootstrap_quality.tsv`,
aggregated into `all_bootstrap_quality.tsv`). This feature is still under development:

```yaml
evaluation:
  run_bootstrap_quality: true   # computationally intensive — off by default
  num_bootstrap: 100            # 50 for testing, 100 standard, 1000 for publication
```

Only applies to non-penalty-sweep tools (tools in `penalty_tools` are skipped, since they don't have
a single `{tool}_contribution.dat`).

---

## Troubleshooting

### Conda environment creation fails

```bash
conda update -n base conda
conda env create -f install/environment.yml
```

### Tool execution fails

```bash
# Check tool-specific logs
less logs/{condition}_{tool}.stderr.log

# Recreate a broken conda env
rm -rf .snakemake/conda/<hash>
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --cores 2 --use-conda --conda-create-envs-only
```

### Out of memory

```bash
# Reduce parallel jobs
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --cores 2 --use-conda

# Or reduce samples in the config
num_samples: 4
```

### Missing input files

```bash
# Verify COSMIC signatures are present
ls input/COSMIC_v3.3.1_SBS_GRCh38.txt

# Verify tissue data for your cancer type
ls cosmic_tissue_data/signature_contributions_WGS_*.dat
```

### Getting Help

- **Issues**: [GitHub Issues](https://github.com/mkatsanto/sigfittest_workflow/issues)
- **Original Paper**: [Medo et al., 2024](https://doi.org/10.1038/s41467-024-53711-6)

---

## Performance Expectations

Based on benchmarking with 8 samples, 1,000,000 mutations, 1 realization:

| Speed | Tools | Approx. Runtime |
|-------|-------|-----------------|
| Very fast | MutationalPatterns, deconstructSigs, sigminer_NNLS | < 2 min |
| Fast | SPA, sigminer_QP | 2–10 min |
| Moderate | MuSiCal, sigLASSO, YAPSA | 10–30 min |
| Slow | SigsPack, signature_tools, sigminer_SA | 30–120 min |
| Very slow | mmsig, sigfit (MCMC) | > 2 hours |

Runtime scales linearly with `num_samples` and `num_realizations`. A full analysis across many tools and mutation counts should be run on a cluster.

---

## Citation

If you use this workflow, please cite:

> Matúš Medo, Charlotte K. Y. Ng, Michaela Medová (2024).
> "A comprehensive comparison of tools for fitting mutational signatures."
> *Nature Communications* 15, 9467.
> https://doi.org/10.1038/s41467-024-53711-6

**Snakemake:**
> Mölder, F. et al. Sustainable data analysis with Snakemake. *F1000Res* 10, 33 (2021).
> https://doi.org/10.12688/f1000research.29032.2

### Third-Party Signature Fitting Tools

This workflow benchmarks the following tools — please also cite whichever ones you use in your results:

| Tool (`tools_to_run` name) | Citation |
|---|---|
| `SPA` (SigProfilerAssignment) | Díaz-Gay, M. et al. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. *Bioinformatics* 39, btad756 (2023). [doi:10.1093/bioinformatics/btad756](https://doi.org/10.1093/bioinformatics/btad756) |
| `mutationalpatterns` | Blokzijl, F., Janssen, R., van Boxtel, R. et al. MutationalPatterns: comprehensive genome-wide analysis of mutational processes. *Genome Med* 10, 33 (2018). [doi:10.1186/s13073-018-0539-0](https://doi.org/10.1186/s13073-018-0539-0) |
| `deconstructsigs` | Rosenthal, R., McGranahan, N., Herrero, J., Taylor, B.S., Swanton, C. deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. *Genome Biol* 17, 31 (2016). [doi:10.1186/s13059-016-0893-4](https://doi.org/10.1186/s13059-016-0893-4) |
| `MuSiCal` | Jin, H., Gulhan, D.C. et al. Accurate and sensitive mutational signature analysis with MuSiCal. *Nat Genet* 56, 541–552 (2024). [doi:10.1038/s41588-024-01659-0](https://doi.org/10.1038/s41588-024-01659-0) |
| `sigfit` | Gori, K., Baez-Ortega, A. sigfit: flexible Bayesian inference of mutational signatures. *bioRxiv* 372896 (2018, preprint). [doi:10.1101/372896](https://doi.org/10.1101/372896) |
| `mmsig` | Rustad, E.H., Nadeu, F., Angelopoulos, N. et al. mmsig: a fitting approach to accurately identify somatic mutational signatures in hematological malignancies. *Commun Biol* 4, 424 (2021). [doi:10.1038/s42003-021-01938-0](https://doi.org/10.1038/s42003-021-01938-0) |
| `siglasso` | Li, S., Crawford, F.W., Gerstein, M.B. Using sigLASSO to optimize cancer mutation signatures jointly with sampling likelihood. *Nat Commun* 11, 3575 (2020). [doi:10.1038/s41467-020-17388-x](https://doi.org/10.1038/s41467-020-17388-x) |
| `sigminer` | Wang, S., Li, H., Song, M. et al. Copy number signature analysis tool and its application in prostate cancer reveals distinct mutational processes and clinical outcomes. *PLoS Genet* 17, e1009557 (2021). [doi:10.1371/journal.pgen.1009557](https://doi.org/10.1371/journal.pgen.1009557) |
| `signature_tools` (signature.tools.lib) | Degasperi, A. et al. A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies. *Nat Cancer* 1, 249–263 (2020). [doi:10.1038/s43018-020-0027-5](https://doi.org/10.1038/s43018-020-0027-5) |
| `sigspack` | Schumann, F., Blanc, E., Messerschmidt, C. et al. SigsPack, a package for cancer mutational signatures. *BMC Bioinformatics* 20, 450 (2019). [doi:10.1186/s12859-019-3043-7](https://doi.org/10.1186/s12859-019-3043-7) |
| `yapsa` | Hübschmann, D., Jopp-Saile, L. et al. Analysis of mutational signatures with yet another package for signature analysis. *Genes Chromosomes Cancer* 60, 314–331 (2021). [doi:10.1002/gcc.22918](https://doi.org/10.1002/gcc.22918) |

The `iterative`/`iterative_version1` and `extended` methods are original to this workflow (see [Method Descriptions](#method-descriptions)) and don't need a separate citation beyond the SigFitTest paper above.

---

## License

This workflow is distributed under the same license as the original [SigFitTest](https://github.com/8medom/SigFitTest) repository.
