# SigFitBench

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.20-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A reproducible Snakemake workflow for benchmarking mutational signature fitting tools.
It runs a dozen established third-party tools alongside two purpose-built methods for the
unknown-signature problem — `SigID` and `SigTrim` (not yet public, see note below) — against
synthetic and real mutation catalogs, scores them all with the same ground truth and metrics,
and produces a single interactive HTML report comparing the results.

> **Note**: `FSPA` (`github.com/sirenfrappe/FSPA`), `SigTrim`, and `SigID` (both
> `github.com/mkatsanto/...`) are packaged in separate repos that aren't public yet — they'll be
> released soon. Until then they can't be invoked from `tools_to_run`; use the other tools listed
> under [Available Tools](#available-tools) in the meantime.

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

**Methodology paper's original repository**: [SigFitTest](https://github.com/8medom/SigFitTest)

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
git clone https://github.com/mkatsanto/SigFitBench.git
cd SigFitBench

# Create base environment (a few minutes)
conda env create -f install/environment.yml

# Activate
conda activate sigfitbench
```

Tool-specific conda environments are created automatically on first run via the `--use-conda` flag.

---

## Quick Start

### Local Execution

**Option A — plain flags:**

```bash
conda activate sigfitbench
cd SigFitBench

# Quick smoke test (3 fast tools, small data, finishes in minutes)
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
`cores: 12`) into a single reusable profile, so you don't have to repeat them every time:

```bash
conda activate sigfitbench
cd SigFitBench

snakemake --profile profiles/local-conda \
  --snakefile workflow/Snakefile --configfile config/config_test.yaml
```

Override the profile's `cores` setting with `--cores N` on the command line, or edit
`profiles/local-conda/config.yaml` directly, to match your machine's parallelism.

**First run**: conda environments are created automatically (30–60 min). Subsequent runs skip this step.

### SLURM Cluster Execution

`profiles/slurm-generic` is a cluster-agnostic Snakemake profile built on Snakemake 8's native
`executor: slurm` plugin (already pinned in `install/environment.yml` as
`snakemake-executor-plugin-slurm`) — no separate submission scripts needed. It already sets
`use-conda: true` and `rerun-incomplete: true`, so those flags don't need to be repeated on the
command line. Copy it and fill in the commented-out `slurm_partition`/`slurm_account`/`qos` fields
under `default-resources` for your own cluster (check `sinfo`/`sacctmgr show associations` or your
cluster's docs for the right values).

```bash
conda activate sigfitbench
cd SigFitBench

snakemake --profile profiles/slurm-generic \
  --snakefile workflow/Snakefile --configfile config/config_test.yaml
```

- **Account/QoS/partition**: set under `default-resources` in your profile's `config.yaml`.
- **Per-rule time limits**: `set-resources` in the same file, keyed by rule name (e.g. `sigminer`,
  `sigfit`), each a `runtime:` in minutes tuned to that rule's expected cost — cluster-agnostic as
  shipped; add/adjust an entry there instead of passing `--time` on the command line.
- **Parallel jobs**: the profile submits up to `jobs: 50` at once by default — override with `--jobs N`.
- **Per-job SLURM stdout/stderr**: handled automatically by the executor plugin under
  `.snakemake/slurm_logs/`. This is separate from each rule's own `log:` files under `log_dir`,
  which already work the same locally or on SLURM.

`profiles/slurm-conda` is the same profile already filled in for UBELIX (University of Bern HPC) —
use it as a worked example if you're on UBELIX, or as a reference for what real values look like
when adapting `slurm-generic` to your own cluster.

---

## Configuration

All config files live in `config/`. There is no default `config.yaml` — you must always pass `--configfile`.

### Available Configs

| File | Purpose |
|------|---------|
| `config_test.yaml` | Quick smoke test: `signature_tools`, `MutationalPatterns`, `SPA`, `MuSiCal`, 2 cancer types, 5 samples, ~10% unknown-signature fraction — confirm the install works before scaling up |
| `config_mysigmix.yaml` | Example `"custom"` data type run, drawing signature mixes from a user-supplied combinations file |

All configs share the same structure (see [Key Parameters](#key-parameters) below) — copy one and
adjust `tools_to_run`, `num_muts`, `cancer_types`, `num_samples`, and `out_of_reference_weights` for
your own run. Runtime scales with `num_samples`, `num_muts`, `num_realizations`, and how many
tools are selected — see [Performance Expectations](#performance-expectations).

### Key Parameters

```yaml
# Tools to run
tools_to_run: ["signature_tools", "MutationalPatterns", "SPA", "MuSiCal"]
# Use "all" to run every tool below, or a list of specific names — casing must match the
# names in the Available Tools table below (it's used verbatim in output filenames).

# Simulation parameters
num_muts: [100000, 1000000]       # Mutations per sample
num_samples: 30                    # Samples per cohort — each gets its own true signature mixture
num_realizations: 1                # Repeat draws of the *same* mixtures with fresh sampling noise
cancer_types: ["Head-SCC"]         # Cancer types (must match cosmic_tissue_data/)
data_types: ["synthetic"]          # "synthetic", "realsim", "external", "custom" — see below
relevance_modes: ["standard"]      # "standard", "relevant", or both

# Unknown signature testing (set weights to [] to disable)
out_of_reference_sigs: ["SBS86", "SBS90", "SBS92"]
out_of_reference_weights: [0.1, 0.2]   # Fraction of mutations from unknown sigs

# Output
output_dir: "results_test"
log_dir: "logs_test"
cluster_log_dir: "logs_test/slurm"     # SLURM job stdout/stderr — required by every config, see SLURM section above

# Reproducibility
random_seed: 0
```

**`num_samples` vs. `num_realizations`** — these vary genuinely different things, not just "more data" two ways:
- `num_samples` draws that many *distinct* ground-truth signature mixtures (which signatures are
  active, at what proportions) — more samples means testing a tool across a wider diversity of
  possible tumors.
- `num_realizations` reruns the *same* set of mixtures from `num_samples`, redrawing only the
  actual mutation counts sampled from them — every rep has identical ground truth, differing only
  in random sampling noise.

A single realization can't tell you whether a tool's error came from a real weakness in fitting
that particular signature mixture, or just an unlucky draw — this matters most at low `num_muts`,
where multinomial sampling noise is proportionally large relative to the signal. Multiple
realizations of the same mixtures let the benchmark average that noise out and separate a tool's
systematic bias from unavoidable sampling variance.

### Available Tools

`tools_to_run` accepts `"all"` or any of the following names. The exact casing shown here is
required — it's used verbatim in each tool's `{tool}_contribution.dat` output filename:

`deconstructSigs`, `SPA`, `SPA_low`, `FSPA`†, `MutationalPatterns`, `MuSiCal`, `SigID`†, `SigTrim`†,
`sigfit`, `mmsig`, `sigminer`, `sigspack`, `yapsa`, `siglasso`, `signature_tools`

† Not runnable yet — repo not public, see the note near the top of this README.

`SigID` and `SigTrim` are this workflow's own methods for the unknown-signature problem — `SigID`
models unknown signal via an L1-penalised search, `SigTrim` trims anomalous residual out before
fitting (see [Method Descriptions](#method-descriptions) for detail); the rest are established
third-party tools.

> **Known install issues**: some tools' conda environments are more fragile than others,
> particularly on Apple Silicon (`mmsig`, `sigfit`) — see each tool's `workflow/envs/*.yaml` and
> `*.post-deploy.sh` for how it's installed. If env creation fails for a tool you need, please
> [open an issue](https://github.com/mkatsanto/SigFitBench/issues).

Add per-tool overrides (thresholds, penalties, `metric: cosine`/`l2`, etc.) under `tool_overrides` —
see the blocks already in `config_test.yaml` for every tool's available knobs.

### Data Types

`data_types` controls where each condition's mutation catalog comes from. Combine any of these in one
run (each activates a different Snakemake rule):

| `data_types` value | Rule | Source |
|---|---|---|
| `"synthetic"` | `ms_create_synthetic_catalog` | Samples signature weights from empirical COSMIC tissue distributions (`cosmic_tissue_data/`), then draws mutations via multinomial sampling. |
| `"realsim"` | `create_realsim_catalog` | Real PCAWG samples for selected cancer types, at their native mutation count — no downsampling, no ground truth; `rep` is always 0. |
| `"external"` | `prepare_external_catalog` | Brings your own catalog, no ground truth (see the commented example in `config_mysigmix.yaml`) — every column is fit as its own sample, no per-tissue splitting like `realsim`; `num_muts` and `rep` are always `"none"`/`0`. |
| `"custom"` | `create_custom_catalog` | Draws from a user-supplied `combinations_file` under `custom_datasets` (see `config_mysigmix.yaml` for a working example) — each sample randomly picks one row (signature mix) from the table. |

For `"realsim"`, cancer types come from `simulation.real_samples.cancer_types` (a separate list from
the top-level `cancer_types`, since real PCAWG labels don't match the synthetic ones, e.g.
`Breast-AdenoCA` vs. `Breast`). Since there's no ground truth, `SigID` is skipped for this mode even
once its repo is public (its auto-penalty search requires ground truth as an input), and results are
compared via cross-tool agreement rather than TAE/Precision/Recall:

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
sample's mutation catalog, and how much does each contribute?* That question assumes the
reference catalog is complete. It isn't: COSMIC adds new signatures with nearly every
release (e.g. v3.4 added SBS96–99 on top of v3.3.1), and any given tumor can carry a
mutational process that's specific to that patient, study, or not yet catalogued at all.
When a tool is only allowed to explain a sample using its known reference, any truly
unmodelled signal has nowhere to go — it gets silently folded into whichever known
signature(s) look closest, inflating false positives and distorting the true known
signatures' weights. That's the **unknown-signature problem**, and it's what
`out_of_reference_sigs`/`out_of_reference_weights` are built to stress-test (see
[Unknown Signature Testing](#unknown-signature-testing)): specific COSMIC signatures are
included when generating the synthetic catalog but withheld from the reference every tool
fits against, so every tool is forced to cope with mutation burden it can't attribute to
anything it knows about.

The key distinction isn't "does this tool explicitly discover what an unknown signature
looks like" (only `SigID` does that) — it's "does this tool's algorithm have *any*
confidence-thresholding or residual-trimming mechanism that lets it refuse to
force-explain 100% of the signal," versus always fully absorbing unknown signal into
whichever known signatures look closest. A few tools, packaged separately at
[github.com/mkatsanto/SigID](https://github.com/mkatsanto/SigID) and
[github.com/mkatsanto/SigTrim](https://github.com/mkatsanto/SigTrim), share that kind of
thinking:

| Method | Approach |
|---|---|
| `SigID` | Models unknown signal directly: an L1-penalised search lets one or more data-driven "unknown" signatures join the active set when they genuinely improve reconstruction. |
| `SigTrim` | Detects mutation-context channels carrying anomalous residual and trims them out before fitting, so unknown signal can't distort the known-signature weights it reports. |
| `signature_tools` | `signature.tools.lib::Fit()`'s confidence-based exposure filtering (`giniScaledThreshold`) also refuses to force-explain low-confidence signal — though `scripts/signature_tools.R` discards that "unassigned" row before writing output, so it isn't currently surfaced as its own field. |
| `sigfit` | Its package supports `fit_extract_signatures()`, which fits known signatures and extracts novel ones simultaneously — but `scripts/sigfit.R` calls plain `fit_signatures()` instead, so this workflow doesn't currently exercise it. |
| `deconstructSigs` | `whichSignatures()` returns an explicit `unknown` field, but `scripts/deconstructSigs.R` only reads the per-signature weights, not that field. |

Everything else (`SPA`, `MutationalPatterns`, `MuSiCal` ) fits purely against the
known reference with no such mechanism, so any out-of-reference mutation burden gets
silently absorbed rather than flagged.

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
│                   └── {tool}_contribution.dat          # Fitted signature weights
│
├── evaluation/                             # Per-run evaluation files
│   └── *_per_sample.tsv
│
├── generated_data/                         # Synthetic catalogs and ground truth
│
├── all_evaluation_metrics.tsv              # Aggregated metrics across all tools and conditions
├── all_per_sample_metrics.tsv              # Per-sample metrics with features
└── sigfitbench_report.html                 # Interactive HTML report (main output)
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
open {output_dir}/sigfitbench_report.html        # macOS
xdg-open {output_dir}/sigfitbench_report.html    # Linux
```

A quick-stats strip and short orientation guide sit at the top; full methodology, parameter
tables, and metric definitions live in the Appendix at the end. One tab per cancer type
contains:

- **Tool comparison**: precision-recall trade-offs and performance metrics vs. mutation burden
- **Performance distribution per tool**: violin/strip plots of the metric spread across samples
- **Sample-by-sample comparison**: which samples are hard for every tool at once (a ranked
  "hardest samples" summary, related to signature spectral complexity, once a cancer type has
  more than 10 samples per condition)
- **Signature complexity vs. performance**: shown only when the run includes unknown
  (out-of-reference) signatures
- **Signature activity per sample**: ground-truth vs. predicted per-signature activity for one
  sample at a time

An **Aggregate Analysis** tab pools these views across all cancer types.

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

Generate a diagram of the workflow's rule graph or full DAG with `graphviz`'s `dot`:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --dag | dot -Tpdf > workflow_dag.pdf

snakemake --snakefile workflow/Snakefile --configfile config/config_test.yaml \
  --rulegraph | dot -Tpdf > workflow_rules.pdf
```

`profiles/graphs` bundles the flags these need (`dryrun: true`, `printshellcmds: true`, and the
`snakefile:` path) into a reusable profile, so you don't have to pass `--snakefile`/`--dry-run`
yourself:

```bash
snakemake --profile profiles/graphs --configfile config/config_test.yaml \
  --dag | dot -Tpdf > workflow_dag.pdf
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

- **Issues**: [GitHub Issues](https://github.com/mkatsanto/SigFitBench/issues)
- **Original Paper**: [Medo et al., 2024](https://doi.org/10.1038/s41467-024-53711-6)
- **Contact**: maria.katsantoni@unibe.ch, matus.medo@unibe.ch, qixuan.wang@unibe.ch

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

`SigID` and `SigTrim` as well as `FSPA`  are original to this workflow (see [Method Descriptions](#method-descriptions)) and they will be published in the future.

---

## License

MIT License — see [LICENSE](LICENSE) for the full text.
